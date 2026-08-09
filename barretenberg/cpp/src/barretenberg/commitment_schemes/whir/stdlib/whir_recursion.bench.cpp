#include "barretenberg/commitment_schemes/whir/stdlib/recursion_harness.hpp"

#include "barretenberg/commitment_schemes/mercury/vela_honk.hpp"
#include "barretenberg/commitment_schemes/whir/stdlib/transparent_honk_recursive_verifier.hpp"
#include "barretenberg/crypto/skyscraper/skyscraper.hpp"
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
 * @brief One WHIR-Honk proof carried all the way to a pairing-based proof, with every stage timed.
 *
 * @details The pipeline a rollup would actually run: an inner statement proved with transparent
 * UltraHonk over WHIR, that whole proof verified inside a circuit, and that circuit proved with a
 * pairing-based scheme whose verifier fits on chain. The inner proof never touches the chain, so
 * its size is a bandwidth question rather than a gas one; what settles is the outer proof.
 *
 * Both outer schemes are measured because they answer different questions. UltraHonk + KZG is the
 * production path (Gemini/Shplonk, ~13 KB). UltraHonk + Vela is the same shell over the Mercury-style
 * univariate PCS, which trades prover work for an opening argument the EVM checks more cheaply.
 */
template <typename Recursive, typename StdlibHasher>
void pipeline_for(const char* label, const std::vector<size_t>& inner_gate_counts, size_t pow_bits)
{
    using InnerHonk = typename Recursive::NativeHonk;
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
    const auto row = [](const char* name, auto value, const char* unit) {
        std::cout << "    " << std::left << std::setw(26) << name << std::right << std::setw(12) << value << " " << unit
                  << "\n";
    };

    for (const size_t num_gates : inner_gate_counts) {
        UltraCircuitBuilder sizing = build_inner(num_gates);
        const size_t log_n = ProverInstance_<UltraFlavor>(sizing).log_dyadic_size();

        std::cout << "\n" << std::string(78, '=') << "\n";
        std::cout << label << " Merkle hash, inner circuit 2^" << log_n << ", lambda=100, rate 2^-4, k=4, k0=1, "
                  << pow_bits << "-bit grind\n";
        std::cout << std::string(78, '=') << "\n";

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

        // ---- Stage 0: the same statement proved directly with UltraHonk + KZG ---------------
        // The baseline the whole pipeline has to justify itself against: if the inner statement can
        // simply be proved with the scheme that settles on chain, the WHIR detour buys nothing.
        {
            UltraCircuitBuilder direct = build_inner(num_gates);
            DefaultIO::add_default(direct);
            auto start_direct = std::chrono::steady_clock::now();
            auto instance = std::make_shared<ProverInstance_<UltraFlavor>>(direct);
            auto vk = std::make_shared<UltraFlavor::VerificationKey>(instance->get_precomputed());
            UltraProver_<UltraFlavor> prover(instance, vk);
            const auto key_ms = elapsed_ms(start_direct);
            start_direct = std::chrono::steady_clock::now();
            const HonkProof proof = prover.construct_proof();
            const auto prove_ms = elapsed_ms(start_direct);
            auto vk_and_hash = std::make_shared<UltraFlavor::VKAndHash>(vk);
            start_direct = std::chrono::steady_clock::now();
            UltraVerifier_<UltraFlavor, DefaultIO> verifier(vk_and_hash);
            const bool ok = verifier.verify_proof(proof).result;
            const auto verify_ms = elapsed_ms(start_direct);
            std::cout << "  [0] baseline: the same statement, UltraHonk + KZG directly (no WHIR)\n";
            row("proving key", key_ms, "ms");
            row("prove", prove_ms, "ms");
            row("verify (native)", verify_ms, "ms");
            row("proof", proof.size() * sizeof(fr), "B");
            row("verified", ok ? "yes" : "NO", "");
            std::cout << std::flush;
        }

        // ---- Stage 1: the inner proof, transparent UltraHonk over WHIR ----------------------
        UltraCircuitBuilder inner_builder = build_inner(num_gates);
        auto start = std::chrono::steady_clock::now();
        auto pk = InnerHonk::create_proving_key(inner_builder, config);
        const auto inner_key_ms = elapsed_ms(start);

        start = std::chrono::steady_clock::now();
        const HonkProof inner_proof = InnerHonk::prove(pk);
        const auto inner_prove_ms = elapsed_ms(start);

        start = std::chrono::steady_clock::now();
        const bool inner_ok = InnerHonk::verify(pk.vk, config, inner_proof);
        const auto inner_verify_ms = elapsed_ms(start);
        if (!inner_ok) {
            throw_or_abort("whir_recursion_bench: the inner proof does not verify");
        }
        std::cout << "  [1] inner: UltraHonk + WHIR (" << label << ")\n";
        row("proving key", inner_key_ms, "ms");
        row("prove", inner_prove_ms, "ms");
        row("verify (native)", inner_verify_ms, "ms");
        row("proof", inner_proof.size() * sizeof(fr), "B");
        std::cout << std::flush;

        // ---- Stage 2: that proof verified inside a circuit -----------------------------------
        UltraCircuitBuilder outer;
        typename StdlibTranscript<UltraCircuitBuilder>::Proof stdlib_proof;
        stdlib_proof.reserve(inner_proof.size());
        for (const fr& element : inner_proof) {
            stdlib_proof.push_back(FF::from_witness(&outer, element));
        }
        start = std::chrono::steady_clock::now();
        const std::vector<FF> inner_public_inputs = Recursive::verify(outer, pk.vk, config, stdlib_proof);
        const auto build_ms = elapsed_ms(start);
        for (const FF& input : inner_public_inputs) {
            outer.set_public_input(input.get_witness_index());
        }
        const size_t recursive_gates = outer.get_num_finalized_gates_inefficient();
        const size_t table_rows = outer.get_tables_size();
        std::cout << "  [2] recursive verifier circuit\n";
        row("build (witness gen)", build_ms, "ms");
        row("gates", recursive_gates, "");
        row("lookup table rows", table_rows, "");
        std::cout << std::flush;

        // ---- Stage 3a: settle with UltraHonk + KZG -------------------------------------------
        {
            UltraCircuitBuilder settle = outer;
            DefaultIO::add_default(settle);
            start = std::chrono::steady_clock::now();
            auto instance = std::make_shared<ProverInstance_<UltraFlavor>>(settle);
            auto vk = std::make_shared<UltraFlavor::VerificationKey>(instance->get_precomputed());
            const auto key_ms = elapsed_ms(start);

            UltraProver_<UltraFlavor> prover(instance, vk);
            start = std::chrono::steady_clock::now();
            const HonkProof proof = prover.construct_proof();
            const auto prove_ms = elapsed_ms(start);

            auto vk_and_hash = std::make_shared<UltraFlavor::VKAndHash>(vk);
            start = std::chrono::steady_clock::now();
            UltraVerifier_<UltraFlavor, DefaultIO> verifier(vk_and_hash);
            const bool ok = verifier.verify_proof(proof).result;
            const auto verify_ms = elapsed_ms(start);

            std::cout << "  [3a] outer: UltraHonk + KZG\n";
            row("circuit",
                instance->dyadic_size(),
                ("rows (2^" + std::to_string(instance->log_dyadic_size()) + ")").c_str());
            row("proving key", key_ms, "ms");
            row("prove", prove_ms, "ms");
            row("verify (native)", verify_ms, "ms");
            row("proof", proof.size() * sizeof(fr), "B");
            row("verified", ok ? "yes" : "NO", "");
            std::cout << std::flush;
        }

        // ---- Stage 3b: settle with UltraHonk + Vela ------------------------------------------
        {
            using VelaHonk = bb::vela::VelaHonk;
            UltraCircuitBuilder settle = outer;
            UltraCircuitBuilder sizing_settle = settle;
            const size_t outer_log_n = ProverInstance_<UltraFlavor>(sizing_settle).log_dyadic_size();
            const auto vela_config = VelaHonk::make_config(outer_log_n);

            start = std::chrono::steady_clock::now();
            auto vela_pk = VelaHonk::create_proving_key(settle, vela_config);
            const auto key_ms = elapsed_ms(start);

            start = std::chrono::steady_clock::now();
            const HonkProof proof = VelaHonk::prove(vela_pk);
            const auto prove_ms = elapsed_ms(start);

            start = std::chrono::steady_clock::now();
            const bool ok = VelaHonk::verify(vela_pk.vk, vela_config, proof);
            const auto verify_ms = elapsed_ms(start);

            std::cout << "  [3b] outer: UltraHonk + Vela\n";
            row("circuit", size_t(1) << outer_log_n, ("rows (2^" + std::to_string(outer_log_n) + ")").c_str());
            row("proving key", key_ms, "ms");
            row("prove", prove_ms, "ms");
            row("verify (native)", verify_ms, "ms");
            row("proof", proof.size() * sizeof(fr), "B");
            row("verified", ok ? "yes" : "NO", "");
            std::cout << std::flush;
        }
    }
}

void pipeline(const std::string& which, const std::vector<size_t>& inner_gate_counts, size_t pow_bits)
{
    using Builder = UltraCircuitBuilder;
    std::cout << "\nWHIR-Honk end to end: inner proof -> recursive verifier -> pairing-based outer proof\n";
    if (which == "all" || which == "poseidon2") {
        pipeline_for<whir::recursion::TransparentHonkRecursiveVerifier<Builder>,
                     whir::recursion::StdlibPoseidon2Hasher<Builder>>("Poseidon2", inner_gate_counts, pow_bits);
    }
    if (which == "all" || which == "skyscraper") {
        pipeline_for<
            whir::recursion::TransparentHonkRecursiveVerifier<Builder,
                                                              whir::recursion::RecursionSkyscraperWhirPcs,
                                                              whir::recursion::StdlibSkyscraperHasher<Builder>>,
            whir::recursion::StdlibSkyscraperHasher<Builder>>("Skyscraper", inner_gate_counts, pow_bits);
    }
    if (which == "all" || which == "blake3") {
        pipeline_for<whir::recursion::TransparentHonkRecursiveVerifier<Builder,
                                                                       whir::recursion::RecursionBlake3sWhirPcs,
                                                                       whir::recursion::StdlibBlake3sHasher<Builder>>,
                     whir::recursion::StdlibBlake3sHasher<Builder>>("Blake3s", inner_gate_counts, pow_bits);
    }
}

/**
 * @brief The inner (prover-side) half of the pipeline for every Merkle hash bb's WHIR can commit
 * with, including the two that have no in-circuit verifier yet.
 * @details Stage [1] alone: how fast the hash is to prove with, what its proof costs on the wire,
 * and how fast it verifies natively. That is the half of the trade-off a recursion-free deployment
 * sees, and it runs the ordering opposite to the in-circuit half.
 */
void inner_hashers(const std::vector<size_t>& inner_gate_counts)
{
    using FF = fr;
    const auto build_inner = [](size_t num_gates) {
        UltraCircuitBuilder builder;
        MockCircuits::add_arithmetic_gates_with_public_inputs(builder, 16);
        MockCircuits::add_arithmetic_gates(builder, num_gates);
        MockCircuits::add_lookup_gates(builder, 2);
        const size_t rom_id = builder.create_ROM_array(4);
        for (size_t i = 0; i < 4; ++i) {
            builder.set_ROM_element(rom_id, i, builder.add_variable(FF(3 * i + 1)));
        }
        builder.read_ROM_array(rom_id, builder.add_variable(FF(2)));
        return builder;
    };
    const auto elapsed_ms = [](auto start) {
        return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count();
    };

    std::cout << "\nInner WHIR-Honk proof by Merkle hash (lambda=100, rate 2^-4, k=4, k0=1, 20-bit grind)\n";
    std::cout << std::string(84, '=') << "\n";
    std::cout << std::left << std::setw(14) << "hash" << std::setw(9) << "inner" << std::right << std::setw(12) << "key"
              << std::setw(12) << "prove" << std::setw(12) << "verify" << std::setw(14) << "proof" << "\n";

    for (const size_t num_gates : inner_gate_counts) {
        UltraCircuitBuilder sizing = build_inner(num_gates);
        const size_t log_n = ProverInstance_<UltraFlavor>(sizing).log_dyadic_size();

        const auto run = [&](const char* name, auto pcs_tag) {
            using Pcs = decltype(pcs_tag);
            using Honk = honk_transparent::TransparentHonk<Pcs, UltraProveKitFlavor>;
            WhirConfig config = WhirConfig::create(log_n,
                                                   /*security_bits=*/100,
                                                   /*log_inv_rate=*/4,
                                                   /*folding_factor_bits=*/4,
                                                   /*final_poly_bits=*/4,
                                                   WhirSoundness::REPAIRED_LIST,
                                                   /*zk=*/false,
                                                   /*max_stack_bits=*/0,
                                                   /*initial_folding_factor_bits=*/1,
                                                   /*pow_bits=*/20);
            config.enable_recursion_profile();

            UltraCircuitBuilder builder = build_inner(num_gates);
            auto start = std::chrono::steady_clock::now();
            auto pk = Honk::create_proving_key(builder, config);
            const auto key_ms = elapsed_ms(start);
            start = std::chrono::steady_clock::now();
            const HonkProof proof = Honk::prove(pk);
            const auto prove_ms = elapsed_ms(start);
            start = std::chrono::steady_clock::now();
            const bool ok = Honk::verify(pk.vk, config, proof);
            const auto verify_ms = elapsed_ms(start);
            if (!ok) {
                throw_or_abort("whir_recursion_bench: inner proof does not verify");
            }
            std::cout << std::left << std::setw(14) << name << std::setw(9) << ("2^" + std::to_string(log_n))
                      << std::right << std::setw(9) << key_ms << " ms" << std::setw(9) << prove_ms << " ms"
                      << std::setw(9) << verify_ms << " ms" << std::setw(11) << proof.size() * sizeof(fr) << " B\n"
                      << std::flush;
        };

        run("Poseidon2", whir::recursion::RecursionWhirPcsFor<Poseidon2CompressionHasher>{});
        run("Skyscraper", whir::recursion::RecursionWhirPcsFor<SkyscraperMerkleHasher>{});
        run("Blake3s", whir::recursion::RecursionWhirPcsFor<Blake3sMerkleHasher>{});
    }
}

/**
 * @brief What a Skyscraper "bar" would cost in circuit, priced before writing one.
 * @details The bar rotates the *canonical* little-endian byte string of an Fr by 16 bytes, S-boxes
 * each byte and reduces. Everything but the canonicity is cheap and lookup-shaped; the canonical
 * decomposition is the term that decides whether the hash is viable in a recursive verifier, so it
 * is worth measuring on its own rather than assuming.
 */
void bar_cost()
{
    using FF = stdlib::field_t<UltraCircuitBuilder>;
    std::cout << "\nSkyscraper bar, priced in barretenberg\n";
    std::cout << std::string(60, '=') << "\n";

    for (const size_t count : std::vector<size_t>{ 1, 2, 10 }) {
        UltraCircuitBuilder builder;
        std::vector<FF> inputs;
        for (size_t i = 0; i < count; ++i) {
            inputs.push_back(FF::from_witness(&builder, fr::random_element()));
        }
        const size_t before = builder.get_num_finalized_gates_inefficient();
        for (const FF& x : inputs) {
            stdlib::byte_array<UltraCircuitBuilder> bytes(x, 32);
            // Keep the decomposition alive so it is not optimised away.
            builder.set_public_input(FF(bytes[0]).normalize().get_witness_index());
        }
        const size_t after = builder.get_num_finalized_gates_inefficient();
        std::cout << "  " << std::setw(3) << count << " canonical 32-byte decompositions: " << std::setw(8)
                  << (after - before) << " gates\n";
    }
    std::cout << "  (the difference between successive rows is the marginal cost; the first row\n"
                 "   also pays for the shared range lists)\n";
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

/**
 * @brief The same statement verified with each Merkle hasher, so the choice is measured.
 * @details Poseidon2 and Blake3s are opposites: Blake3s is the fastest of bb's WHIR hashers to
 * *prove* with and the most expensive to *verify* in a circuit, because in a BN254 circuit it is
 * 32-bit word arithmetic and every field element it absorbs has to be decomposed into 32 bytes
 * first. The configurations here are deliberately small — a realistic Blake3s verifier does not fit
 * on this machine, which is itself the finding.
 */
void hashers()
{
    using Poseidon2Harness =
        whir::recursion::RecursionHarness_<whir::recursion::StdlibPoseidon2Hasher<UltraCircuitBuilder>>;
    using Blake3sHarness =
        whir::recursion::RecursionHarness_<whir::recursion::StdlibBlake3sHasher<UltraCircuitBuilder>>;

    std::cout << "\nMerkle hash choice, same statement both times (1 column, rate 2^-4, k=4, k0=1, no grind)\n";
    std::cout << std::string(92, '=') << "\n";
    std::cout << std::left << std::setw(22) << "configuration" << std::right << std::setw(14) << "poseidon2"
              << std::setw(14) << "blake3s" << std::setw(12) << "blake3s/p2" << std::setw(14) << "p2 proof"
              << std::setw(14) << "b3 proof" << "\n";

    struct Point {
        size_t num_variables;
        size_t security_bits;
    };
    for (const Point& point : std::vector<Point>{ { 8, 32 }, { 8, 64 }, { 10, 32 } }) {
        const auto make = [&](auto harness_tag) {
            using Harness = decltype(harness_tag);
            const WhirConfig config = Harness::config(point.num_variables,
                                                      point.security_bits,
                                                      /*log_inv_rate=*/4,
                                                      /*k=*/4,
                                                      /*k0=*/1,
                                                      /*pow_bits=*/0);
            const auto instance = Harness::make_instance(config, { 1 });
            if (!Harness::verify_natively(instance)) {
                throw_or_abort("whir_recursion_bench: the proof under measurement does not verify");
            }
            auto builder = Harness::build_circuit(instance);
            return std::pair<size_t, size_t>{ builder.get_num_finalized_gates_inefficient(),
                                              instance.proof.size() * sizeof(fr) };
        };
        const auto [p2_gates, p2_proof] = make(Poseidon2Harness{});
        const auto [b3_gates, b3_proof] = make(Blake3sHarness{});
        std::cout << std::left << std::setw(22)
                  << ("m=" + std::to_string(point.num_variables) + ", lambda=" + std::to_string(point.security_bits))
                  << std::right << std::setw(14) << p2_gates << std::setw(14) << b3_gates << std::setw(11) << std::fixed
                  << std::setprecision(1) << double(b3_gates) / double(p2_gates) << "x" << std::setw(12) << p2_proof
                  << " B" << std::setw(12) << b3_proof << " B\n";
    }
    std::cout << "\n  Blake3s digests are two field elements where Poseidon2's is one, so the Blake3s proof also\n"
                 "  carries roughly twice the Merkle-path bytes.\n";
}

} // namespace

int main(int argc, char** argv)
{
    const std::string mode = argc > 1 ? argv[1] : "all";
    if (mode == "hashers") {
        hashers();
        std::cout << std::endl;
        return 0;
    }
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
    if (mode == "pipeline") {
        srs::init_file_crs_factory(srs::bb_crs_path());
        const std::string which = argc > 2 ? argv[2] : "all";
        const size_t inner_gates = argc > 3 ? static_cast<size_t>(std::stoul(argv[3])) : 100;
        pipeline(which, { inner_gates }, /*pow_bits=*/20);
        std::cout << std::endl;
        return 0;
    }
    if (mode == "bar-cost") {
        bar_cost();
        std::cout << std::endl;
        return 0;
    }
    if (mode == "inner-hashers") {
        inner_hashers({ 100, 20000 });
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
