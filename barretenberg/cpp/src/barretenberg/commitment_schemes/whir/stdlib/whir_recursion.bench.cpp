#include "barretenberg/commitment_schemes/whir/stdlib/recursion_harness.hpp"

#include "barretenberg/commitment_schemes/whir/stdlib/transparent_honk_recursive_verifier.hpp"
#include "barretenberg/stdlib_circuit_builders/mock_circuits.hpp"
#include "barretenberg/ultra_honk/prover_instance.hpp"

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

/**
 * @brief The whole pipeline for a *complete* WHIR-Honk proof: inner, verifier circuit, outer.
 * @details Every column is measured on the same run — the inner proof is produced, verified
 * natively, verified again inside a circuit, and that circuit is proved with UltraHonk + KZG.
 */
void full(const std::vector<size_t>& inner_gate_counts, size_t pow_bits)
{
    using Recursive = whir::recursion::TransparentHonkRecursiveVerifier<UltraCircuitBuilder>;
    using InnerHonk = Recursive::NativeHonk;
    using FF = stdlib::field_t<UltraCircuitBuilder>;

    const auto build_inner = [](size_t num_gates) {
        UltraCircuitBuilder builder;
        MockCircuits::add_arithmetic_gates_with_public_inputs(builder, 16);
        MockCircuits::add_arithmetic_gates(builder, num_gates);
        MockCircuits::add_lookup_gates(builder, 2);
        const size_t rom_id = builder.create_ROM_array(4);
        for (size_t i = 0; i < 4; ++i) {
            builder.set_ROM_element(rom_id, i, builder.add_variable(fr(3 * i + 1)));
        }
        builder.read_ROM_array(rom_id, builder.add_variable(fr(2)));
        return builder;
    };
    const auto elapsed_ms = [](auto start) {
        return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count();
    };

    std::cout << "\nWhole WHIR-Honk proofs, verified recursively and settled with UltraHonk + KZG\n";
    std::cout << "lambda=100, rate 2^-4, k=4, k0=1, " << pow_bits
              << "-bit grind; UltraProveKitFlavor (Poseidon2 Merkle)\n";
    std::cout << std::string(100, '=') << "\n";
    std::cout << std::left << std::setw(9) << "inner" << std::right << std::setw(11) << "in.prove" << std::setw(11)
              << "in.verify" << std::setw(12) << "in.proof" << std::setw(12) << "rec.gates" << std::setw(10) << "outer"
              << std::setw(11) << "out.prove" << std::setw(11) << "out.verify" << std::setw(11) << "out.proof" << "\n";

    for (const size_t num_gates : inner_gate_counts) {
        UltraCircuitBuilder sizing = build_inner(num_gates);
        const size_t log_n = ProverInstance_<UltraFlavor>(sizing).log_dyadic_size();

        WhirConfig config = WhirConfig::create(log_n,
                                               /*security_bits=*/100,
                                               /*log_inv_rate=*/4,
                                               /*folding_factor_bits=*/4,
                                               /*final_poly_bits=*/4,
                                               WhirSoundness::REPAIRED_LIST,
                                               /*zk=*/false,
                                               /*max_stack_bits=*/0,
                                               /*initial_folding_factor_bits=*/1,
                                               pow_bits);
        config.enable_recursion_profile();

        UltraCircuitBuilder inner_builder = build_inner(num_gates);
        auto pk = InnerHonk::create_proving_key(inner_builder, config);
        auto start = std::chrono::steady_clock::now();
        const HonkProof inner_proof = InnerHonk::prove(pk);
        const auto inner_prove_ms = elapsed_ms(start);

        start = std::chrono::steady_clock::now();
        const bool inner_ok = InnerHonk::verify(pk.vk, config, inner_proof);
        const auto inner_verify_ms = elapsed_ms(start);
        if (!inner_ok) {
            throw_or_abort("whir_recursion_bench: the inner proof does not verify");
        }

        UltraCircuitBuilder outer;
        typename StdlibTranscript<UltraCircuitBuilder>::Proof stdlib_proof;
        stdlib_proof.reserve(inner_proof.size());
        for (const fr& element : inner_proof) {
            stdlib_proof.push_back(FF::from_witness(&outer, element));
        }
        const std::vector<FF> inner_public_inputs = Recursive::verify(outer, pk.vk, config, stdlib_proof);
        for (const FF& input : inner_public_inputs) {
            outer.set_public_input(input.get_witness_index());
        }
        const size_t recursive_gates = outer.get_num_finalized_gates_inefficient();
        DefaultIO::add_default(outer);

        auto prover_instance = std::make_shared<ProverInstance_<UltraFlavor>>(outer);
        auto verification_key = std::make_shared<UltraFlavor::VerificationKey>(prover_instance->get_precomputed());
        UltraProver_<UltraFlavor> prover(prover_instance, verification_key);
        start = std::chrono::steady_clock::now();
        const HonkProof outer_proof = prover.construct_proof();
        const auto outer_prove_ms = elapsed_ms(start);

        auto vk_and_hash = std::make_shared<UltraFlavor::VKAndHash>(verification_key);
        start = std::chrono::steady_clock::now();
        UltraVerifier_<UltraFlavor, DefaultIO> verifier(vk_and_hash);
        const bool outer_ok = verifier.verify_proof(outer_proof).result;
        const auto outer_verify_ms = elapsed_ms(start);
        if (!outer_ok) {
            throw_or_abort("whir_recursion_bench: the outer proof does not verify");
        }

        std::cout << std::left << std::setw(9) << ("2^" + std::to_string(log_n)) << std::right << std::setw(9)
                  << inner_prove_ms << " ms" << std::setw(9) << inner_verify_ms << " ms" << std::setw(10)
                  << inner_proof.size() * sizeof(fr) << " B" << std::setw(12) << recursive_gates << std::setw(10)
                  << ("2^" + std::to_string(prover_instance->log_dyadic_size())) << std::setw(9) << outer_prove_ms
                  << " ms" << std::setw(9) << outer_verify_ms << " ms" << std::setw(9)
                  << outer_proof.size() * sizeof(fr) << " B\n";
    }
}

/**
 * @brief The dominant term of a *ProveKit* recursive verifier, measured in barretenberg.
 *
 * @details ProveKit's Spartan verifier is not succinct in the circuit it proves. It builds an
 * eq-table over the whole constraint hypercube (`calculate_evaluations_over_boolean_hypercube_for_eq`
 * in `provekit/common/src/utils/sumcheck.rs`) and multiplies the three R1CS matrices by it
 * (`multiply_transposed_by_eq_alpha`), so its work is `2^m + nnz(A) + nnz(B) + nnz(C)` field
 * operations against the *raw* matrices, which it takes as an input.
 *
 * Their own gnark recursive verifier does exactly that in-circuit — `matrix_evaluation.go` loops
 * over every nonzero with `api.Add(ans, api.Mul(value, api.Mul(rowEval[row], colEval[col])))` — so
 * this is not an artifact of one implementation. This probe reproduces that inner loop in bb stdlib
 * and reports the per-nonzero and per-hypercube-entry gate cost, which is what prices the whole
 * approach: the matrix indices are public, so the products are the cost and the loop is linear in
 * the inner circuit.
 */
void provekit_cost()
{
    using FF = stdlib::field_t<UltraCircuitBuilder>;

    std::cout << "\nProveKit's R1CS term, priced in barretenberg\n";
    std::cout << std::string(60, '=') << "\n";

    // The eq-table: 2^m entries, each one multiplication.
    for (const size_t log_m : std::vector<size_t>{ 10, 12 }) {
        UltraCircuitBuilder builder;
        std::vector<FF> r;
        for (size_t i = 0; i < log_m; ++i) {
            r.push_back(FF::from_witness(&builder, fr::random_element()));
        }
        std::vector<FF> table{ FF(1) };
        for (const FF& coordinate : r) {
            std::vector<FF> next(table.size() * 2);
            for (size_t j = 0; j < table.size(); ++j) {
                next[2 * j + 1] = table[j] * coordinate;
                next[2 * j] = table[j] - next[2 * j + 1];
            }
            table = std::move(next);
        }
        const size_t gates = builder.get_num_finalized_gates_inefficient();
        std::cout << "  eq table over 2^" << log_m << ": " << gates << " gates ("
                  << double(gates) / double(size_t(1) << log_m) << " per entry)\n";
    }

    // The matrix loop: one product of two witness evaluations per nonzero, scaled by a public entry.
    for (const size_t nonzeros : std::vector<size_t>{ 10000, 50000 }) {
        UltraCircuitBuilder builder;
        std::vector<FF> row_eval;
        std::vector<FF> col_eval;
        for (size_t i = 0; i < 256; ++i) {
            row_eval.push_back(FF::from_witness(&builder, fr::random_element()));
            col_eval.push_back(FF::from_witness(&builder, fr::random_element()));
        }
        FF accumulator(0);
        for (size_t i = 0; i < nonzeros; ++i) {
            // Matrix structure and values are public, so only the evaluation product is a gate.
            accumulator += row_eval[i % 256] * col_eval[(i * 7) % 256] * bb::fr(uint64_t(i % 5) + 1);
        }
        accumulator.assert_is_not_zero();
        const size_t gates = builder.get_num_finalized_gates_inefficient();
        std::cout << "  " << nonzeros << " nonzeros: " << gates << " gates (" << double(gates) / double(nonzeros)
                  << " per nonzero)\n";
    }
    std::cout << "\n  A ProveKit verification key for the 2^19 passport circuit is 9.1 MB, essentially all\n"
                 "  matrix data, so this term runs to millions of gates and grows linearly with the inner\n"
                 "  circuit. The WHIR opening it also contains costs about 92,000 gates and grows\n"
                 "  logarithmically. Full ProveKit recursion is bounded by Spartan, not by WHIR.\n";
}

} // namespace

int main(int argc, char** argv)
{
    const std::string mode = argc > 1 ? argv[1] : "all";
    if (mode == "provekit-cost") {
        provekit_cost();
        std::cout << std::endl;
        return 0;
    }
    if (mode == "full") {
        srs::init_file_crs_factory(srs::bb_crs_path());
        // Distinct dyadic sizes: the lookup table alone floors the mock circuit at 2^13.
        full({ 100, 20000, 100000 }, /*pow_bits=*/20);
        full({ 100, 20000, 100000 }, /*pow_bits=*/0);
        std::cout << std::endl;
        return 0;
    }
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
