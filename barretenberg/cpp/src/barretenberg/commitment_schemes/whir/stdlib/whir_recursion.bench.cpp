#include "barretenberg/commitment_schemes/whir/stdlib/recursion_harness.hpp"

#include "barretenberg/common/timer.hpp"
#include "barretenberg/special_public_inputs/special_public_inputs.hpp"
#include "barretenberg/srs/global_crs.hpp"
#include "barretenberg/ultra_honk/ultra_prover.hpp"
#include "barretenberg/ultra_honk/ultra_verifier.hpp"

#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

// Gate-count sweep for the in-circuit WHIR verifier: what each optimization is worth, how the cost
// moves with the schedule, and how many proofs an aggregator of a given size holds.
//
// Every row builds a real WHIR proof and a real UltraHonk circuit that verifies it, so the numbers
// are measured rather than modelled. Usage:
//
//   whir_recursion_bench [levers|schedule|columns|aggregate|all]
//
// The default runs everything.

namespace {

using namespace bb;
using namespace bb::whir;
using Harness = whir::recursion::RecursionHarness;
using GateReport = whir::recursion::GateReport;

/** @brief The layout a transparent UltraHonk proof opens with: 24 precomputed + 8 witness columns
 * in four commitment phases, 29 claimed unshifted and 5 also claimed shifted. Mirrors the generated
 * `params.nr` of the Noir verifier so the two are measuring the same statement. */
const std::vector<size_t> HONK_GROUPS = { 21, 3, 3, 2 };
const std::vector<WhirColumnRef> HONK_SHIFTED = { { 1, 0 }, { 1, 1 }, { 1, 2 }, { 2, 2 }, { 3, 1 } };

/** @brief ProveKit's shape: Spartan reduces the statement to a single committed vector. */
const std::vector<size_t> PROVEKIT_GROUPS = { 1 };

struct Row {
    std::string label;
    size_t gates;
    size_t queries;
    size_t proof_fields;
};

size_t measure(const WhirConfig& config,
               const std::vector<size_t>& groups,
               const std::vector<WhirColumnRef>& shifted,
               GateReport* report = nullptr,
               size_t* proof_fields = nullptr)
{
    const Harness::Instance instance = Harness::make_instance(config, groups, shifted);
    if (!Harness::verify_natively(instance)) {
        throw_or_abort("whir_recursion_bench: the proof under measurement does not verify");
    }
    if (proof_fields != nullptr) {
        *proof_fields = instance.proof.size();
    }
    Harness::Builder builder = Harness::build_circuit(instance, report);
    return builder.get_num_finalized_gates_inefficient();
}

void print_table(const std::string& title, const std::vector<Row>& rows, size_t baseline)
{
    std::cout << "\n" << title << "\n";
    std::cout << std::string(title.size(), '=') << "\n";
    std::cout << std::left << std::setw(46) << "configuration" << std::right << std::setw(10) << "gates"
              << std::setw(11) << "vs base" << std::setw(9) << "round-0" << std::setw(12) << "proof (Fr)"
              << "\n";
    for (const Row& row : rows) {
        std::cout << std::left << std::setw(46) << row.label << std::right << std::setw(10) << row.gates;
        if (baseline > 0) {
            const double ratio = 100.0 * static_cast<double>(row.gates) / static_cast<double>(baseline);
            std::cout << std::setw(10) << std::fixed << std::setprecision(0) << ratio << "%";
        } else {
            std::cout << std::setw(11) << "-";
        }
        std::cout << std::setw(9) << row.queries << std::setw(12) << row.proof_fields << "\n";
    }
}

/**
 * @brief What each recursion-profile setting is worth, one at a time and then together.
 * @details The base row is the protocol as it stands with only the change a circuit cannot do
 * without — per-query authentication paths, since the batched walk's shape depends on the query
 * values. Everything after it is optional.
 */
void levers()
{
    // The parameters the Noir verifier's `whir-r4-repaired-k4-k0_2` family was generated at.
    const size_t m = 10;
    const size_t lambda = 100;
    const size_t rate = 4;
    const size_t k = 4;
    const auto soundness = WhirSoundness::REPAIRED_LIST;

    auto base = [&](size_t k0, size_t pow_bits) {
        return Harness::config(m, lambda, rate, k, k0, pow_bits, soundness);
    };

    std::vector<Row> rows;
    auto add = [&](const std::string& label, WhirConfig config) {
        size_t proof_fields = 0;
        const size_t gates = measure(config, HONK_GROUPS, HONK_SHIFTED, nullptr, &proof_fields);
        rows.push_back({ label, gates, config.rounds.empty() ? 0 : config.rounds[0].num_queries, proof_fields });
    };

    WhirConfig plain = base(2, 0);
    plain.merkle_cap_levels = 0;
    add("per-query paths only (the base circuit)", plain);

    WhirConfig capped = base(2, 0);
    add("+ Merkle cap", capped);

    WhirConfig ground = base(2, 20);
    ground.merkle_cap_levels = 0;
    add("+ Poseidon2 grinding, 20 bits", ground);

    WhirConfig narrow = base(1, 0);
    narrow.merkle_cap_levels = 0;
    add("+ k0 = 1 (narrow first fold)", narrow);

    add("cap + grinding", base(2, 20));
    add("cap + grinding + k0 = 1 (recommended)", base(1, 20));
    add("cap + grinding(24) + k0 = 1", base(1, 24));

    print_table("Optimization levers - transparent-Honk layout, m=10, lambda=100, rate 2^-4, k=4, repaired list",
                rows,
                rows.front().gates);

    GateReport report;
    const WhirConfig best = base(1, 20);
    const size_t gates = measure(best, HONK_GROUPS, HONK_SHIFTED, &report);
    std::cout << "\nWhere the constraints go in the best configuration (" << gates << " gates)\n";
    for (const std::string& phase : report.phases()) {
        const size_t phase_gates = report.gates(phase);
        std::cout << "  " << std::left << std::setw(30) << phase << std::right << std::setw(9) << phase_gates
                  << std::setw(8) << std::fixed << std::setprecision(1)
                  << 100.0 * static_cast<double>(phase_gates) / static_cast<double>(gates) << "%\n";
    }
}

/** @brief The rate/soundness dial: a longer codeword needs fewer queries but a deeper tree. */
void schedule()
{
    const size_t m = 10;
    const size_t lambda = 100;
    std::vector<Row> rows;
    for (const size_t rate : std::vector<size_t>{ 2, 3, 4, 6 }) {
        for (const size_t pow_bits : std::vector<size_t>{ size_t(0), size_t(20) }) {
            const WhirConfig config =
                Harness::config(m, lambda, rate, /*k=*/4, /*k0=*/1, pow_bits, WhirSoundness::REPAIRED_LIST);
            size_t proof_fields = 0;
            const size_t gates = measure(config, HONK_GROUPS, HONK_SHIFTED, nullptr, &proof_fields);
            rows.push_back({ "rate 2^-" + std::to_string(rate) + ", grind " + std::to_string(pow_bits) + " bits",
                             gates,
                             config.rounds.empty() ? 0 : config.rounds[0].num_queries,
                             proof_fields });
        }
    }
    print_table("Rate and grinding - transparent-Honk layout, m=10, lambda=100, k=4, k0=1, repaired list",
                rows,
                rows.front().gates);

    rows.clear();
    for (const size_t vars : std::vector<size_t>{ 8, 10, 12, 14 }) {
        const WhirConfig config = Harness::config(
            vars, lambda, /*log_inv_rate=*/4, /*k=*/4, /*k0=*/1, /*pow_bits=*/20, WhirSoundness::REPAIRED_LIST);
        size_t proof_fields = 0;
        const size_t gates = measure(config, HONK_GROUPS, HONK_SHIFTED, nullptr, &proof_fields);
        rows.push_back({ "inner circuit 2^" + std::to_string(vars),
                         gates,
                         config.rounds.empty() ? 0 : config.rounds[0].num_queries,
                         proof_fields });
    }
    print_table(
        "Inner circuit size - transparent-Honk layout, lambda=100, rate 2^-4, k=4, k0=1, 20-bit grind", rows, 0);
}

/** @brief Column count is the round-0 leaf width, and round 0 is where the cost is. */
void columns()
{
    const size_t lambda = 100;
    std::vector<Row> rows;
    struct Shape {
        std::string label;
        std::vector<size_t> groups;
        std::vector<WhirColumnRef> shifted;
    };
    const std::vector<Shape> shapes = {
        { "ProveKit: 1 column, 1 group", PROVEKIT_GROUPS, {} },
        { "4 columns, 1 group", { 4 }, {} },
        { "8 columns, 2 groups", { 4, 4 }, { { 1, 0 } } },
        { "transparent Honk: 29 columns, 4 groups", HONK_GROUPS, HONK_SHIFTED },
    };
    for (const Shape& shape : shapes) {
        const WhirConfig config = Harness::config(
            12, lambda, /*log_inv_rate=*/4, /*k=*/4, /*k0=*/1, /*pow_bits=*/20, WhirSoundness::REPAIRED_LIST);
        size_t proof_fields = 0;
        const size_t gates = measure(config, shape.groups, shape.shifted, nullptr, &proof_fields);
        rows.push_back({ shape.label, gates, config.rounds.empty() ? 0 : config.rounds[0].num_queries, proof_fields });
    }
    print_table("Committed column count - m=12, lambda=100, rate 2^-4, k=4, k0=1, 20-bit grind", rows, 0);
}

/** @brief How many WHIR proofs fit in one outer UltraHonk circuit. */
void aggregate()
{
    const WhirConfig config = Harness::config(12,
                                              /*security_bits=*/100,
                                              /*log_inv_rate=*/4,
                                              /*k=*/4,
                                              /*k0=*/1,
                                              /*pow_bits=*/20,
                                              WhirSoundness::REPAIRED_LIST);

    std::cout << "\nAggregation - ProveKit shape (1 column), m=12, lambda=100, rate 2^-4, k0=1, 20-bit grind\n";
    std::cout << "==========================================================================================\n";
    std::cout << std::left << std::setw(12) << "proofs" << std::right << std::setw(12) << "gates" << std::setw(14)
              << "gates/proof" << std::setw(14) << "outer 2^k" << "\n";

    std::vector<Harness::Instance> instances;
    for (size_t i = 0; i < 4; ++i) {
        instances.push_back(Harness::make_instance(config, PROVEKIT_GROUPS, {}));
    }
    for (const size_t count : std::vector<size_t>{ 1, 2, 4 }) {
        Harness::Builder builder;
        for (size_t i = 0; i < count; ++i) {
            Harness::verify_in_circuit(builder, instances[i]);
        }
        const size_t gates = builder.get_num_finalized_gates_inefficient();
        size_t log_n = 1;
        while ((size_t(1) << log_n) < gates) {
            ++log_n;
        }
        std::cout << std::left << std::setw(12) << count << std::right << std::setw(12) << gates << std::setw(14)
                  << gates / count << std::setw(14) << ("2^" + std::to_string(log_n)) << "\n";
    }
    std::cout << "\nThe per-proof figure is what matters: an outer circuit of 2^k rows holds 2^k / (gates per"
                 " proof) of them.\n";
}

/**
 * @brief The whole point: N WHIR proofs verified in one circuit, proved once with UltraHonk + KZG.
 * @details The output is a 13 KB pairing-based proof that Ethereum checks for about 630k gas with
 * the `--optimized` Solidity verifier, whatever N is. That is the aggregation story — WHIR's own
 * proof never touches the chain.
 */
void outer(size_t count)
{
    const WhirConfig config = Harness::config(12,
                                              /*security_bits=*/100,
                                              /*log_inv_rate=*/4,
                                              /*k=*/4,
                                              /*k0=*/1,
                                              /*pow_bits=*/20,
                                              WhirSoundness::REPAIRED_LIST);

    std::cout << "\nOuter proof - " << count
              << " ProveKit-shape WHIR proofs aggregated into one UltraHonk + KZG proof\n";

    Harness::Builder builder;
    for (size_t i = 0; i < count; ++i) {
        const Harness::Instance instance = Harness::make_instance(config, PROVEKIT_GROUPS, {});
        if (!Harness::verify_natively(instance)) {
            throw_or_abort("whir_recursion_bench: inner proof does not verify");
        }
        Harness::verify_in_circuit(builder, instance);
    }
    // UltraHonk circuits carry a pairing-point accumulator on their public inputs.
    DefaultIO::add_default(builder);

    auto start = std::chrono::steady_clock::now();
    auto prover_instance = std::make_shared<ProverInstance_<UltraFlavor>>(builder);
    auto verification_key = std::make_shared<UltraFlavor::VerificationKey>(prover_instance->get_precomputed());
    const auto key_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count();

    UltraProver_<UltraFlavor> prover(prover_instance, verification_key);
    start = std::chrono::steady_clock::now();
    const HonkProof proof = prover.construct_proof();
    const auto prove_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count();

    auto vk_and_hash = std::make_shared<UltraFlavor::VKAndHash>(verification_key);
    start = std::chrono::steady_clock::now();
    UltraVerifier_<UltraFlavor, DefaultIO> verifier(vk_and_hash);
    const bool verified = verifier.verify_proof(proof).result;
    const auto verify_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count();

    std::cout << "  outer circuit          2^" << prover_instance->log_dyadic_size() << " rows ("
              << prover_instance->dyadic_size() << ")\n";
    std::cout << "  key generation         " << key_ms << " ms\n";
    std::cout << "  prove                  " << prove_ms << " ms\n";
    std::cout << "  verify                 " << verify_ms << " ms\n";
    std::cout << "  proof                  " << proof.size() << " field elements (" << proof.size() * sizeof(fr)
              << " bytes)\n";
    std::cout << "  verified               " << (verified ? "yes" : "NO") << "\n";
}

} // namespace

int main(int argc, char** argv)
{
    const std::string mode = argc > 1 ? argv[1] : "all";
    if (mode == "outer") {
        srs::init_file_crs_factory(srs::bb_crs_path());
        outer(argc > 2 ? static_cast<size_t>(std::stoul(argv[2])) : 1);
        std::cout << std::endl;
        return 0;
    }
    if (mode == "levers" || mode == "all") {
        levers();
    }
    if (mode == "schedule" || mode == "all") {
        schedule();
    }
    if (mode == "columns" || mode == "all") {
        columns();
    }
    if (mode == "aggregate" || mode == "all") {
        aggregate();
    }
    std::cout << std::endl;
    return 0;
}
