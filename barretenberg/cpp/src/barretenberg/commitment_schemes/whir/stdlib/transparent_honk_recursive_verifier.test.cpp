#include "barretenberg/commitment_schemes/whir/stdlib/transparent_honk_recursive_verifier.hpp"

#include "barretenberg/circuit_checker/circuit_checker.hpp"
#include "barretenberg/stdlib_circuit_builders/mock_circuits.hpp"
#include "barretenberg/ultra_honk/prover_instance.hpp"

#include <gtest/gtest.h>

namespace bb::whir::recursion {

namespace {

using Builder = UltraCircuitBuilder;
using FF = stdlib::field_t<Builder>;
using Verifier = TransparentHonkRecursiveVerifier<Builder>;
using Honk = Verifier::NativeHonk;

/**
 * @brief An inner circuit exercising every relation `UltraProveKitFlavor` keeps.
 * @details Arithmetic, public inputs, a lookup and a ROM read, so the recursive verifier really
 * evaluates all of the flavor's subrelations rather than a degenerate subset. No elliptic,
 * non-native-field or Poseidon2 gates — the flavor refuses those.
 */
Builder build_inner_circuit(size_t num_gates)
{
    Builder builder;
    MockCircuits::add_arithmetic_gates_with_public_inputs(builder, 16);
    MockCircuits::add_arithmetic_gates(builder, num_gates);
    MockCircuits::add_lookup_gates(builder, 2);
    const size_t rom_id = builder.create_ROM_array(4);
    for (size_t i = 0; i < 4; ++i) {
        builder.set_ROM_element(rom_id, i, builder.add_variable(fr(3 * i + 1)));
    }
    const uint32_t read_index = builder.add_variable(fr(2));
    builder.read_ROM_array(rom_id, read_index);
    return builder;
}

/** @brief A proven inner instance: the verification key, the WHIR config and the proof. */
struct Inner {
    WhirConfig config;
    typename Honk::VerificationKey vk;
    HonkProof proof;
};

Inner prove_inner(size_t num_gates = 16, size_t security_bits = 64, size_t pow_bits = 20)
{
    // ProverInstance finalizes (mutates) the builder it consumes, so size with a throwaway copy.
    Builder sizing_builder = build_inner_circuit(num_gates);
    const size_t log_n = ProverInstance_<UltraFlavor>(sizing_builder).log_dyadic_size();

    Builder builder = build_inner_circuit(num_gates);
    WhirConfig config = Honk::make_config(log_n, security_bits, /*log_inv_rate=*/4);
    config = WhirConfig::create(log_n,
                                security_bits,
                                /*log_inv_rate=*/4,
                                /*folding_factor_bits=*/4,
                                /*final_poly_bits=*/4,
                                WhirSoundness::REPAIRED_LIST,
                                /*zk=*/false,
                                /*max_stack_bits=*/0,
                                /*initial_folding_factor_bits=*/1,
                                pow_bits);
    config.enable_recursion_profile();

    auto pk = Honk::create_proving_key(builder, config);
    return { config, pk.vk, Honk::prove(pk) };
}

Builder build_verifier_circuit(const Inner& inner, GateReport* report = nullptr)
{
    Builder builder;
    typename StdlibTranscript<Builder>::Proof proof;
    proof.reserve(inner.proof.size());
    for (const fr& element : inner.proof) {
        proof.push_back(FF::from_witness(&builder, element));
    }
    const std::vector<FF> public_inputs = Verifier::verify(builder, inner.vk, inner.config, proof, report);
    // An aggregator would forward these; here it is enough that they exist and are constrained.
    for (const FF& input : public_inputs) {
        builder.set_public_input(input.get_witness_index());
    }
    return builder;
}

} // namespace

TEST(TransparentHonkRecursiveVerifierTest, VerifiesAWholeWhirHonkProof)
{
    const Inner inner = prove_inner();
    ASSERT_TRUE(Honk::verify(inner.vk, inner.config, inner.proof));

    GateReport report;
    Builder builder = build_verifier_circuit(inner, &report);
    EXPECT_TRUE(CircuitChecker::check(builder));

    info("full WHIR-Honk recursive verifier: ",
         builder.get_num_finalized_gates_inefficient(),
         " gates for a 2^",
         inner.vk.log_dyadic_size,
         " inner circuit (proof ",
         inner.proof.size() * sizeof(fr),
         " bytes)");
    for (const std::string& phase : report.phases()) {
        info("  ", phase, ": ", report.gates(phase));
    }
}

TEST(TransparentHonkRecursiveVerifierTest, VerifiesWithoutGrinding)
{
    const Inner inner = prove_inner(/*num_gates=*/16, /*security_bits=*/64, /*pow_bits=*/0);
    ASSERT_TRUE(Honk::verify(inner.vk, inner.config, inner.proof));
    EXPECT_TRUE(CircuitChecker::check(build_verifier_circuit(inner)));
}

TEST(TransparentHonkRecursiveVerifierTest, VerifiesALargerInnerCircuit)
{
    const Inner inner = prove_inner(/*num_gates=*/20000);
    ASSERT_GT(inner.vk.log_dyadic_size, 13U);
    ASSERT_TRUE(Honk::verify(inner.vk, inner.config, inner.proof));
    Builder builder = build_verifier_circuit(inner);
    EXPECT_TRUE(CircuitChecker::check(builder));
    info("full WHIR-Honk recursive verifier: ",
         builder.get_num_finalized_gates_inefficient(),
         " gates for a 2^",
         inner.vk.log_dyadic_size,
         " inner circuit");
}

// The whole point of a verifier is that a bad proof cannot satisfy it. Corrupting anywhere in the
// proof — the sumcheck univariates, the claimed evaluations, the commitment roots, the WHIR
// openings — must break the circuit.
TEST(TransparentHonkRecursiveVerifierTest, TamperedProofsMakeTheCircuitUnsatisfiable)
{
    const Inner inner = prove_inner();
    ASSERT_TRUE(CircuitChecker::check(build_verifier_circuit(inner)));

    const size_t stride = std::max<size_t>(1, inner.proof.size() / 24);
    size_t checked = 0;
    for (size_t i = 0; i < inner.proof.size(); i += stride) {
        Inner tampered = inner;
        tampered.proof[i] += fr(1);
        EXPECT_FALSE(CircuitChecker::check(build_verifier_circuit(tampered)))
            << "proof element " << i << " is unconstrained";
        ++checked;
    }
    EXPECT_GT(checked, 20U);
}

// The precomputed-column root is a circuit constant, so it is what ties the verifier to one inner
// circuit: a proof of a different circuit cannot satisfy it.
TEST(TransparentHonkRecursiveVerifierTest, ThePrecomputedRootPinsTheInnerCircuit)
{
    Inner inner = prove_inner();
    ASSERT_TRUE(CircuitChecker::check(build_verifier_circuit(inner)));
    inner.vk.precomputed_commitment += fr(1);
    EXPECT_FALSE(CircuitChecker::check(build_verifier_circuit(inner)));
}

TEST(TransparentHonkRecursiveVerifierTest, AggregatesSeveralProofs)
{
    std::vector<Inner> inners;
    for (size_t i = 0; i < 2; ++i) {
        inners.push_back(prove_inner());
    }

    Builder builder;
    for (const Inner& inner : inners) {
        typename StdlibTranscript<Builder>::Proof proof;
        for (const fr& element : inner.proof) {
            proof.push_back(FF::from_witness(&builder, element));
        }
        [[maybe_unused]] const auto public_inputs = Verifier::verify(builder, inner.vk, inner.config, proof);
    }
    EXPECT_TRUE(CircuitChecker::check(builder));
    info("two whole WHIR-Honk proofs: ", builder.get_num_finalized_gates_inefficient(), " gates");
}

} // namespace bb::whir::recursion
