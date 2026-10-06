/**
 * @file zcash_action.bench.cpp
 * @brief Benchmarks the Orchard Action circuit in every barretenberg configuration and prints one JSON object per
 * (configuration, number of Actions).
 *
 * Usage: zcash_action_bench <system> <comma-separated action counts> [repetitions]
 *   system: halo2-honk-pasta | halo2-honk-bn254 | ultra-zk-bn254 | ultra-zk-pasta
 *
 * Timings: `synthesize_s` is witness synthesis (halo2 table / circuit construction from the Action witnesses),
 * `prove_s` the proof construction from the synthesized witness, `prove_total_s` their sum (comparable to halo2's
 * create_proof, which synthesizes the circuit internally), `verify_s` the verification of one proof, `keygen_s` the
 * proving + verification key generation. Medians over the repetitions are reported. The thread count is the
 * HARDWARE_CONCURRENCY environment variable (default: all cores).
 */
#include "barretenberg/common/bb_bench.hpp"
#include "barretenberg/common/thread.hpp"
#include "barretenberg/flavor/ultra_zk_flavor.hpp"
#include "barretenberg/numeric/random/engine.hpp"
#include "barretenberg/srs/global_crs.hpp"
#include "barretenberg/ultra_honk/prover_instance.hpp"
#include "barretenberg/ultra_honk/ultra_prover.hpp"
#include "barretenberg/ultra_honk/ultra_verifier.hpp"
#include "barretenberg/zcash/halo2/action_circuit.hpp"
#include "barretenberg/zcash/honk/orchard_honk.hpp"
#include "barretenberg/zcash/ultra/action_circuit_ultra.hpp"
#include "barretenberg/zcash/ultra_pasta/action_circuit_ultra_pasta.hpp"
#include "barretenberg/zcash/ultra_pasta/ultra_pasta_honk.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace bb;
using namespace bb::zcash;

namespace {

double seconds_since(std::chrono::steady_clock::time_point start)
{
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

double median(std::vector<double> xs)
{
    std::sort(xs.begin(), xs.end());
    return xs[xs.size() / 2];
}

std::string json_array(const std::vector<double>& xs)
{
    std::ostringstream out;
    out << "[";
    for (size_t i = 0; i < xs.size(); ++i) {
        out << (i ? ", " : "") << xs[i];
    }
    out << "]";
    return out.str();
}

template <typename Cycle> struct Bundle {
    std::vector<typename Orchard<Cycle>::ActionWitness> witnesses;
    std::vector<typename Cycle::FF> public_inputs;
};

template <typename Cycle> Bundle<Cycle> random_bundle(size_t num_actions)
{
    auto& engine = numeric::get_debug_randomness();
    Bundle<Cycle> bundle;
    for (size_t i = 0; i < num_actions; ++i) {
        bundle.witnesses.push_back(Orchard<Cycle>::random_witness(engine));
        const auto pi = Orchard<Cycle>::evaluate(bundle.witnesses.back());
        BB_ASSERT(pi.has_value());
        for (const auto& f : pi->to_field_elements()) {
            bundle.public_inputs.push_back(f);
        }
    }
    return bundle;
}

struct Result {
    std::string system;
    size_t actions = 0;
    size_t rows = 0;  // used rows / gates
    size_t log_n = 0; // log2 of the proven polynomial size
    size_t proof_bytes = 0;
    // With 32-byte compressed curve points, as halo2 serializes them, and without the public inputs (which halo2 proofs
    // do not contain).
    size_t proof_bytes_compressed = 0;
    double keygen_s = 0;
    std::vector<double> synthesize_s;
    std::vector<double> prove_s;
    std::vector<double> verify_s;
    bool verified = false;

    void print() const
    {
        std::vector<double> total(prove_s.size());
        for (size_t i = 0; i < total.size(); ++i) {
            total[i] = synthesize_s[i] + prove_s[i];
        }
        const char* threads = std::getenv("HARDWARE_CONCURRENCY");
        printf("{\"system\": \"%s\", \"actions\": %zu, \"threads\": \"%s\", \"rows\": %zu, \"log_n\": %zu, "
               "\"proof_bytes\": %zu, \"proof_bytes_compressed\": %zu, \"keygen_s\": %.6f, \"synthesize_s\": %.6f, "
               "\"prove_s\": %.6f, \"prove_total_s\": %.6f, \"verify_s\": %.6f, \"prove_total_all_s\": %s, "
               "\"verified\": %s}\n",
               system.c_str(),
               actions,
               threads != nullptr ? threads : "default",
               rows,
               log_n,
               proof_bytes,
               proof_bytes_compressed,
               keygen_s,
               median(synthesize_s),
               median(prove_s),
               median(total),
               median(verify_s),
               json_array(total).c_str(),
               verified ? "true" : "false");
        fflush(stdout);
    }
};

template <typename Cycle> Result bench_halo2_honk(size_t num_actions, size_t reps)
{
    using Circuit = halo2::ActionCircuit<Cycle>;
    using Trace = halo2::AnchoredTrace<Cycle>;
    using Flavor = OrchardFlavor_<Cycle>;
    constexpr bool IS_PASTA = Flavor::IS_PASTA;
    Result result{ .system = IS_PASTA ? "halo2-honk-pasta" : "halo2-honk-bn254", .actions = num_actions };
    const auto bundle = random_bundle<Cycle>(num_actions);

    auto synthesize = [&]() {
        auto table = Circuit::build(bundle.witnesses, bundle.public_inputs);
        result.rows = table.num_rows;
        return Trace::build(table, Flavor::TRACE_OFFSET);
    };
    // Warm-up: generator tables and the CRS are derived once per process.
    auto trace = synthesize();
    auto start = std::chrono::steady_clock::now();
    OrchardProvingKey_<Cycle> pk(trace);
    result.keygen_s = seconds_since(start);
    result.log_n = pk.log_circuit_size;
    (void)orchard_prove(pk, trace);

    typename Flavor::Proof proof;
    for (size_t r = 0; r < reps; ++r) {
        start = std::chrono::steady_clock::now();
        trace = synthesize();
        result.synthesize_s.push_back(seconds_since(start));
        start = std::chrono::steady_clock::now();
        proof = orchard_prove(pk, trace);
        result.prove_s.push_back(seconds_since(start));
        start = std::chrono::steady_clock::now();
        result.verified = orchard_verify(*pk.vk, bundle.public_inputs, proof);
        result.verify_s.push_back(seconds_since(start));
    }
    result.proof_bytes = proof.size() * 32;
    // Curve points in the proof: 17 witness commitments, 3 Libra, log_n - 1 Gemini folds, the Shplonk quotient, then
    // the halo2 IPA's S and 2 log_n L/R (Pasta) or the KZG quotient (BN254). A Pasta point serializes to 64 bytes and
    // a BN254 point to 128 bytes (4 field elements); both compress to 32 bytes.
    const size_t log_n = result.log_n;
    const size_t num_points = 17 + 3 + (log_n - 1) + 1 + (IS_PASTA ? 1 + (2 * log_n) : 1);
    result.proof_bytes_compressed = result.proof_bytes - ((IS_PASTA ? 32 : 96) * num_points);
    return result;
}

Result bench_ultra_zk_bn254(size_t num_actions, size_t reps)
{
    using Flavor = UltraZKFlavor;
    using Builder = ultra::Builder;
    Result result{ .system = "ultra-zk-bn254", .actions = num_actions };
    const auto bundle = random_bundle<Bn254Cycle>(num_actions);
    auto synthesize = [&]() {
        Builder builder;
        ultra::ActionCircuitUltra::build(builder, bundle.witnesses, bundle.public_inputs);
        DefaultIO::add_default(builder);
        return builder;
    };

    // Key generation: circuit construction, trace population and the commitments to the precomputed polynomials.
    auto start = std::chrono::steady_clock::now();
    std::shared_ptr<Flavor::VerificationKey> vk;
    {
        auto builder = synthesize();
        result.rows = builder.get_num_finalized_gates_inefficient();
        auto instance = std::make_shared<ProverInstance_<Flavor>>(builder);
        vk = std::make_shared<Flavor::VerificationKey>(instance->get_precomputed());
    }
    result.keygen_s = seconds_since(start);
    result.log_n = vk->log_circuit_size;
    auto vk_and_hash = std::make_shared<Flavor::VKAndHash>(vk);

    HonkProof proof;
    for (size_t r = 0; r < reps + 1; ++r) {
        start = std::chrono::steady_clock::now();
        auto builder = synthesize();
        auto instance = std::make_shared<ProverInstance_<Flavor>>(builder);
        const double synthesize_s = seconds_since(start);
        start = std::chrono::steady_clock::now();
        UltraProver_<Flavor> prover(instance, vk);
        proof = prover.construct_proof();
        const double prove_s = seconds_since(start);
        start = std::chrono::steady_clock::now();
        UltraVerifier_<Flavor, DefaultIO> verifier(vk_and_hash);
        bool ok = verifier.verify_proof(proof).result;
        for (size_t i = 0; i < bundle.public_inputs.size(); ++i) {
            ok = ok && proof[i] == bundle.public_inputs[i];
        }
        const double verify_s = seconds_since(start);
        if (r == 0) {
            continue; // warm-up
        }
        result.synthesize_s.push_back(synthesize_s);
        result.prove_s.push_back(prove_s);
        result.verify_s.push_back(verify_s);
        result.verified = ok;
    }
    result.proof_bytes = proof.size() * 32;
    // Curve points: 9 witness commitments (incl. the Gemini masking polynomial), 3 Libra, CONST_PROOF_SIZE_LOG_N - 1
    // Gemini folds (the proof is padded to a constant size), the Shplonk quotient and the KZG quotient. A BN254 point
    // takes 4 field elements (128 bytes) here and 32 bytes compressed. The proof leads with the public inputs: the
    // Action's and the default IO (pairing points).
    const size_t num_points = 9 + 3 + (CONST_PROOF_SIZE_LOG_N - 1) + 2;
    const size_t num_public_inputs = bundle.public_inputs.size() + DefaultIO::PUBLIC_INPUTS_SIZE;
    result.proof_bytes_compressed = result.proof_bytes - (num_points * 96) - (num_public_inputs * 32);
    return result;
}

Result bench_ultra_zk_pasta(size_t num_actions, size_t reps)
{
    using Flavor = UltraPastaZKFlavor;
    using Builder = ultra_pasta::Builder;
    Result result{ .system = "ultra-zk-pasta", .actions = num_actions };
    const auto bundle = random_bundle<PastaCycle>(num_actions);
    auto synthesize = [&]() {
        Builder builder;
        ultra_pasta::ActionCircuitUltraPasta::build(builder, bundle.witnesses, bundle.public_inputs);
        return builder;
    };

    // Key generation: circuit construction, trace population and the commitments to the precomputed polynomials.
    auto start = std::chrono::steady_clock::now();
    std::shared_ptr<Flavor::VerificationKey> vk;
    {
        auto builder = synthesize();
        result.rows = builder.get_num_finalized_gates_inefficient();
        auto instance = std::make_shared<UltraPastaProverInstance>(builder);
        vk = std::make_shared<Flavor::VerificationKey>(instance->get_precomputed());
    }
    result.keygen_s = seconds_since(start);
    result.log_n = vk->log_circuit_size;

    Flavor::Proof proof;
    for (size_t r = 0; r < reps + 1; ++r) {
        start = std::chrono::steady_clock::now();
        auto builder = synthesize();
        auto instance = std::make_shared<UltraPastaProverInstance>(builder);
        const double synthesize_s = seconds_since(start);
        start = std::chrono::steady_clock::now();
        proof = ultra_pasta_prove(instance, vk);
        const double prove_s = seconds_since(start);
        start = std::chrono::steady_clock::now();
        std::vector<PastaCycle::FF> public_inputs;
        const bool ok = ultra_pasta_verify(vk, proof, &public_inputs) && public_inputs == bundle.public_inputs;
        const double verify_s = seconds_since(start);
        if (r == 0) {
            continue; // warm-up
        }
        result.synthesize_s.push_back(synthesize_s);
        result.prove_s.push_back(prove_s);
        result.verify_s.push_back(verify_s);
        result.verified = ok;
    }
    result.proof_bytes = proof.size() * 32;
    // Curve points: 9 witness commitments (incl. the Gemini masking polynomial), 3 Libra, log_n - 1 Gemini folds, the
    // Shplonk quotient, then the halo2 IPA's S and 2 log_n L/R. A Vesta point serializes to 64 bytes and compresses to
    // 32. The proof leads with the public inputs.
    const size_t num_points = 9 + 3 + (result.log_n - 1) + 1 + 1 + (2 * result.log_n);
    result.proof_bytes_compressed = result.proof_bytes - (num_points * 32) - (bundle.public_inputs.size() * 32);
    return result;
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 3) {
        fprintf(stderr, "usage: %s <system> <action counts> [repetitions]\n", argv[0]);
        return 1;
    }
    const std::string system = argv[1];
    std::vector<size_t> counts;
    {
        std::stringstream ss(argv[2]);
        std::string item;
        while (std::getline(ss, item, ',')) {
            counts.push_back(std::stoul(item));
        }
    }
    const size_t reps = argc > 3 ? std::stoul(argv[3]) : 3;
    // ZCASH_BENCH_PROFILE=1 prints barretenberg's hierarchical timers to stderr.
    const bool profile = std::getenv("ZCASH_BENCH_PROFILE") != nullptr;
    bb::detail::use_bb_bench = profile;
    if (system == "ultra-zk-bn254") {
        srs::init_bn254_file_crs_factory(srs::bb_crs_path());
    }
    for (const size_t n : counts) {
        if (system == "halo2-honk-pasta") {
            bench_halo2_honk<PastaCycle>(n, reps).print();
        } else if (system == "halo2-honk-bn254") {
            bench_halo2_honk<Bn254Cycle>(n, reps).print();
        } else if (system == "ultra-zk-bn254") {
            bench_ultra_zk_bn254(n, reps).print();
        } else if (system == "ultra-zk-pasta") {
            bench_ultra_zk_pasta(n, reps).print();
        } else {
            fprintf(stderr, "unknown system %s\n", system.c_str());
            return 1;
        }
    }
    if (profile) {
        bb::detail::GLOBAL_BENCH_STATS.print_aggregate_counts_hierarchical(std::cerr);
    }
    return 0;
}
