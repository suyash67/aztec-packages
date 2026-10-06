/**
 * @file zcash_action.bench.cpp
 * @brief Benchmarks the Orchard Action circuit in every barretenberg configuration and prints one JSON object per
 * (configuration, number of Actions).
 *
 * Usage: zcash_action_bench <system> <comma-separated action counts> [repetitions]
 *   system: halo2-honk-pasta
 *
 * Timings: `synthesize_s` is witness synthesis (halo2 table / circuit construction from the Action witnesses),
 * `prove_s` the proof construction from the synthesized witness, `prove_total_s` their sum (comparable to halo2's
 * create_proof, which synthesizes the circuit internally), `verify_s` the verification of one proof, `keygen_s` the
 * proving + verification key generation. Medians over the repetitions are reported. The thread count is the
 * HARDWARE_CONCURRENCY environment variable (default: all cores).
 */
#include "barretenberg/common/thread.hpp"
#include "barretenberg/numeric/random/engine.hpp"
#include "barretenberg/zcash/halo2/action_circuit.hpp"
#include "barretenberg/zcash/honk/orchard_honk.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <functional>
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
    size_t proof_bytes_compressed = 0; // with 32-byte compressed curve points, as halo2 serializes them
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

Result bench_halo2_honk_pasta(size_t num_actions, size_t reps)
{
    using Circuit = halo2::ActionCircuit<PastaCycle>;
    using Trace = halo2::AnchoredTrace<PastaCycle>;
    Result result{ .system = "halo2-honk-pasta", .actions = num_actions };
    const auto bundle = random_bundle<PastaCycle>(num_actions);

    auto synthesize = [&]() {
        auto table = Circuit::build(bundle.witnesses, bundle.public_inputs);
        result.rows = table.num_rows;
        return Trace::build(table, OrchardFlavor::TRACE_OFFSET);
    };
    // Warm-up: generator tables and the CRS are derived once per process.
    auto trace = synthesize();
    auto start = std::chrono::steady_clock::now();
    OrchardProvingKey pk(trace);
    result.keygen_s = seconds_since(start);
    result.log_n = pk.log_circuit_size;
    (void)orchard_prove(pk, trace);

    OrchardFlavor::Proof proof;
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
    // Curve points in the proof: 17 witness commitments, 3 Libra, log_n - 1 Gemini folds, Shplonk Q, IPA S and
    // 2 log_n L/R. Each serializes to 64 bytes here and to 32 bytes compressed.
    const size_t log_n = result.log_n;
    const size_t num_points = 17 + 3 + (log_n - 1) + 1 + 1 + (2 * log_n);
    result.proof_bytes_compressed = result.proof_bytes - (32 * num_points);
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
    for (const size_t n : counts) {
        if (system == "halo2-honk-pasta") {
            bench_halo2_honk_pasta(n, reps).print();
        } else {
            fprintf(stderr, "unknown system %s\n", system.c_str());
            return 1;
        }
    }
    return 0;
}
