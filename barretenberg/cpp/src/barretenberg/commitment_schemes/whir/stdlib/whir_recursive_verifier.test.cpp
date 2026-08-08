#include "barretenberg/commitment_schemes/whir/stdlib/recursion_harness.hpp"

#include "barretenberg/circuit_checker/circuit_checker.hpp"

#include <gtest/gtest.h>

namespace bb::whir::recursion {

namespace {

using Harness = RecursionHarness;
using Builder = Harness::Builder;
using Instance = Harness::Instance;

/** @brief The layout a transparent UltraHonk proof opens with: four groups, 29 columns, 5 shifted. */
const std::vector<size_t>& honk_groups()
{
    static const std::vector<size_t> groups = { 21, 3, 3, 2 };
    return groups;
}
const std::vector<WhirColumnRef>& honk_shifted()
{
    static const std::vector<WhirColumnRef> shifted = { { 1, 0 }, { 1, 1 }, { 1, 2 }, { 2, 2 }, { 3, 1 } };
    return shifted;
}

WhirConfig test_config(size_t num_variables, size_t pow_bits = 8)
{
    return Harness::config(num_variables, /*security_bits=*/64, /*log_inv_rate=*/2, /*k=*/4, /*k0=*/1, pow_bits);
}

} // namespace

TEST(WhirRecursiveVerifierTest, VerifiesARealProof)
{
    const Instance instance = Harness::make_instance(test_config(10), { 3, 2 }, { { 1, 0 }, { 1, 1 } });
    ASSERT_TRUE(Harness::verify_natively(instance));

    GateReport report;
    Builder builder = Harness::build_circuit(instance, &report);
    EXPECT_TRUE(CircuitChecker::check(builder));

    info("WHIR recursive verifier: ", builder.get_num_finalized_gates_inefficient(), " gates");
    for (const std::string& phase : report.phases()) {
        info("  ", phase, ": ", report.gates(phase));
    }
}

TEST(WhirRecursiveVerifierTest, VerifiesATransparentHonkColumnLayout)
{
    const Instance instance = Harness::make_instance(test_config(10), honk_groups(), honk_shifted());
    ASSERT_TRUE(Harness::verify_natively(instance));
    Builder builder = Harness::build_circuit(instance);
    EXPECT_TRUE(CircuitChecker::check(builder));
    info("WHIR recursive verifier (29 columns, 4 groups): ", builder.get_num_finalized_gates_inefficient(), " gates");
}

TEST(WhirRecursiveVerifierTest, VerifiesASingleColumnProveKitShape)
{
    // ProveKit's WHIR commits one vector; the round-0 leaf is then the fold coset alone.
    const Instance instance = Harness::make_instance(test_config(12), { 1 });
    ASSERT_TRUE(Harness::verify_natively(instance));
    Builder builder = Harness::build_circuit(instance);
    EXPECT_TRUE(CircuitChecker::check(builder));
    info("WHIR recursive verifier (1 column, m=12): ", builder.get_num_finalized_gates_inefficient(), " gates");
}

TEST(WhirRecursiveVerifierTest, WithoutGrindingAndWithoutACap)
{
    WhirConfig config = test_config(10, /*pow_bits=*/0);
    config.merkle_cap_levels = 0;
    const Instance instance = Harness::make_instance(config, { 2, 1 }, { { 1, 0 } });
    ASSERT_TRUE(Harness::verify_natively(instance));
    EXPECT_TRUE(CircuitChecker::check(Harness::build_circuit(instance)));
}

TEST(WhirRecursiveVerifierTest, ZeroIterationSchedule)
{
    const Instance instance = Harness::make_instance(test_config(4), { 2, 1 }, { { 1, 0 } });
    ASSERT_EQ(instance.config.num_iterations(), 0U);
    ASSERT_TRUE(Harness::verify_natively(instance));
    EXPECT_TRUE(CircuitChecker::check(Harness::build_circuit(instance)));
}

TEST(WhirRecursiveVerifierTest, LargeFoldAndHighRate)
{
    const WhirConfig config =
        Harness::config(12, /*security_bits=*/80, /*log_inv_rate=*/4, /*k=*/3, /*k0=*/2, /*pow_bits=*/12);
    const Instance instance = Harness::make_instance(config, { 4, 2 }, { { 1, 0 }, { 1, 1 } });
    ASSERT_TRUE(Harness::verify_natively(instance));
    EXPECT_TRUE(CircuitChecker::check(Harness::build_circuit(instance)));
}

// Query indices, and so the whole shape of the Merkle work, are fresh randomness on every proof:
// colliding queries, an index of zero, an index at the last leaf. Completeness has to hold for all
// of them, not for the one instance a single-shot test happens to draw.
TEST(WhirRecursiveVerifierTest, CompletenessOverManyRandomInstances)
{
    for (size_t i = 0; i < 8; ++i) {
        const Instance instance = Harness::make_instance(test_config(10), { 3, 2 }, { { 1, 0 }, { 1, 1 } });
        ASSERT_TRUE(Harness::verify_natively(instance)) << "instance " << i;
        EXPECT_TRUE(CircuitChecker::check(Harness::build_circuit(instance))) << "instance " << i;
    }
}

// A circuit that stayed satisfiable for a bad proof would verify nothing. Each of these corrupts one
// part of the proof and asserts the circuit becomes unsatisfiable.
TEST(WhirRecursiveVerifierTest, TamperedProofsMakeTheCircuitUnsatisfiable)
{
    const Instance instance = Harness::make_instance(test_config(10), { 2, 1 }, { { 1, 0 } });
    ASSERT_TRUE(Harness::verify_natively(instance));
    ASSERT_TRUE(CircuitChecker::check(Harness::build_circuit(instance)));

    const size_t stride = std::max<size_t>(1, instance.proof.size() / 24);
    size_t checked = 0;
    for (size_t i = 0; i < instance.proof.size(); i += stride) {
        Instance tampered = instance;
        tampered.proof[i] += fr(1);
        EXPECT_FALSE(CircuitChecker::check(Harness::build_circuit(tampered)))
            << "proof element " << i << " is unconstrained";
        ++checked;
    }
    EXPECT_GT(checked, 20U);
}

TEST(WhirRecursiveVerifierTest, TamperedClaimMakesTheCircuitUnsatisfiable)
{
    Instance instance = Harness::make_instance(test_config(10), { 2, 1 }, { { 1, 0 } });
    ASSERT_TRUE(CircuitChecker::check(Harness::build_circuit(instance)));
    instance.unshifted_evaluations[0] += fr(1);
    EXPECT_FALSE(CircuitChecker::check(Harness::build_circuit(instance)));
    instance.unshifted_evaluations[0] -= fr(1);
    instance.shifted_evaluations[0] += fr(1);
    EXPECT_FALSE(CircuitChecker::check(Harness::build_circuit(instance)));
}

TEST(WhirRecursiveVerifierTest, TamperedOpeningPointMakesTheCircuitUnsatisfiable)
{
    Instance instance = Harness::make_instance(test_config(10), { 2, 1 }, { { 1, 0 } });
    ASSERT_TRUE(CircuitChecker::check(Harness::build_circuit(instance)));
    instance.u[0] += fr(1);
    EXPECT_FALSE(CircuitChecker::check(Harness::build_circuit(instance)));
}

// The aggregator setup: the inner commitment roots are constants of the verification circuit rather
// than proof data, so a wrong root is a wrong constant and the circuit cannot be satisfied.
TEST(WhirRecursiveVerifierTest, ConstantRootsAreBoundToTheOpenings)
{
    const Instance instance =
        Harness::make_instance(test_config(10), { 2, 1 }, { { 1, 0 } }, /*roots_on_proof_stream=*/false);
    ASSERT_TRUE(Harness::verify_natively(instance));
    ASSERT_TRUE(CircuitChecker::check(Harness::build_circuit(instance)));

    for (size_t g = 0; g < instance.roots.size(); ++g) {
        Instance tampered = instance;
        tampered.roots[g] += fr(1);
        EXPECT_FALSE(CircuitChecker::check(Harness::build_circuit(tampered))) << "root " << g << " is unconstrained";
    }
}

// Several proofs in one circuit: the aggregation primitive. Each is verified independently, so the
// cost is additive and a bad proof anywhere breaks the whole circuit.
TEST(WhirRecursiveVerifierTest, AggregatesSeveralProofs)
{
    std::vector<Instance> instances;
    for (size_t i = 0; i < 3; ++i) {
        instances.push_back(Harness::make_instance(test_config(10), { 2, 1 }, { { 1, 0 } }));
    }

    Builder builder;
    for (const Instance& instance : instances) {
        Harness::verify_in_circuit(builder, instance);
    }
    EXPECT_TRUE(CircuitChecker::check(builder));

    Builder with_a_bad_proof;
    for (size_t i = 0; i < instances.size(); ++i) {
        Instance instance = instances[i];
        if (i == 1) {
            instance.unshifted_evaluations[0] += fr(1);
        }
        Harness::verify_in_circuit(with_a_bad_proof, instance);
    }
    EXPECT_FALSE(CircuitChecker::check(with_a_bad_proof));
}

} // namespace bb::whir::recursion
