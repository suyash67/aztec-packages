#include "barretenberg/commitment_schemes/bolt/bolt_honk.hpp"

#include "barretenberg/commitment_schemes/brakedown/brakedown_honk.hpp"
#include "barretenberg/stdlib_circuit_builders/mock_circuits.hpp"
#include "barretenberg/transcript/transcript.hpp"

#include <gtest/gtest.h>

namespace bb::bolt {

namespace {

constexpr double LN_2P32 = 32.0 * 0.6931471805599453;
constexpr double LN_BN254 = 175.8; // ln of the BN254 scalar field order

std::vector<fr> random_array(size_t size)
{
    std::vector<fr> array(size);
    for (fr& value : array) {
        value = fr::random_element();
    }
    return array;
}

fr mle(std::span<const fr> array, std::span<const fr> u)
{
    return Polynomial<fr>(array).evaluate_mle(u);
}

fr shifted_mle(std::span<const fr> array, std::span<const fr> u)
{
    std::vector<fr> shifted(array.begin() + 1, array.end());
    shifted.push_back(fr::zero());
    return mle(shifted, u);
}

using Hasher = whir::Blake3sMerkleHasher;
using CK = ligero::LigeroCommitmentKey<Hasher, BoltCodePolicy>;
using Prover = ligero::LigeroProver<Hasher>;
using Verifier = ligero::LigeroVerifier<Hasher>;

// Small lambda keeps the query count (and so the proof) test-sized.
constexpr size_t TEST_SECURITY_BITS = 8;

ligero::LigeroConfig test_config(size_t num_variables)
{
    return make_bolt_config(num_variables, TEST_SECURITY_BITS, /*est_polynomials=*/4);
}

struct Instance {
    std::vector<std::vector<fr>> unshifted_arrays;
    std::vector<std::vector<fr>> to_be_shifted_arrays;
    std::vector<fr> u;
    std::vector<ligero::LigeroGroupData<Hasher>> groups;
    Prover::Claims prover_claims;
    Verifier::Claims verifier_claims;
};

Instance make_instance(const CK& ck, size_t num_unshifted, size_t num_shifted)
{
    const size_t n = size_t(1) << ck.config.num_variables;
    Instance instance;
    instance.u = random_array(ck.config.num_variables);
    for (size_t p = 0; p < num_unshifted; ++p) {
        instance.unshifted_arrays.push_back(random_array(n));
    }
    for (size_t p = 0; p < num_shifted; ++p) {
        std::vector<fr> array = random_array(n);
        array[0] = fr::zero();
        instance.to_be_shifted_arrays.push_back(std::move(array));
    }
    size_t shifted_group = 0;
    if (num_unshifted > 0) {
        instance.groups.push_back(ck.commit_group(instance.unshifted_arrays));
        shifted_group = 1;
    }
    if (num_shifted > 0) {
        instance.groups.push_back(ck.commit_group(instance.to_be_shifted_arrays));
    }
    for (auto& group : instance.groups) {
        instance.prover_claims.groups.push_back(&group);
        instance.verifier_claims.group_num_columns.push_back(group.num_columns());
    }
    for (size_t c = 0; c < num_unshifted; ++c) {
        instance.prover_claims.unshifted.push_back({ 0, c });
        instance.prover_claims.unshifted_evaluations.push_back(mle(instance.unshifted_arrays[c], instance.u));
    }
    for (size_t c = 0; c < num_shifted; ++c) {
        instance.prover_claims.to_be_shifted.push_back({ shifted_group, c });
        instance.prover_claims.shifted_evaluations.push_back(shifted_mle(instance.to_be_shifted_arrays[c], instance.u));
    }
    instance.verifier_claims.unshifted = instance.prover_claims.unshifted;
    instance.verifier_claims.unshifted_evaluations = instance.prover_claims.unshifted_evaluations;
    instance.verifier_claims.to_be_shifted = instance.prover_claims.to_be_shifted;
    instance.verifier_claims.shifted_evaluations = instance.prover_claims.shifted_evaluations;
    return instance;
}

HonkProof prove_instance(const CK& ck, const Instance& instance)
{
    auto transcript = NativeTranscript::test_prover_init_empty();
    Prover::prove(ck, instance.prover_claims, instance.u, transcript);
    return transcript->export_proof();
}

bool verify_proof(const ligero::LigeroConfig& config,
                  const Verifier::Claims& claims,
                  std::span<const fr> u,
                  const HonkProof& proof)
{
    try {
        auto transcript = std::make_shared<NativeTranscript>(proof);
        [[maybe_unused]] auto init = transcript->template receive_from_prover<fr>("Init");
        const BoltCodePolicy code(config);
        return Verifier::verify(config, claims, u, transcript, code);
    } catch (const std::exception&) {
        return false;
    }
}

} // namespace

// The paper states that at q = 2^32, j = 16, k = 128 the root of omega is x0 ~ 0.094114 and so the
// ensemble's distance is at least 0.0941. Reproducing that from Theorem 6.5's formulas is the
// strongest available check on the transcription.
TEST(BoltParamsTest, DistanceBoundMatchesPaperInstantiation)
{
    const double gamma = ldpc_distance_bound(LN_2P32, 16, 128);
    EXPECT_NEAR(gamma, 0.094114, 1e-4);
    EXPECT_GE(gamma, BoltParams::reference_2p32().gamma);
}

// The finding that sets the BN254 parameters apart from the paper's: omega's leading term is
// x*ln(q-1), so a fixed column degree certifies less distance as the field grows. The paper's own
// (j=16, k=128) is worthless at BN254; the degree has to scale with ln q.
TEST(BoltParamsTest, ColumnDegreeMustScaleWithTheFieldSize)
{
    // Same parameters, two field sizes: a 1000x collapse.
    const double at_2p32 = ldpc_distance_bound(LN_2P32, 16, 128);
    const double at_bn254 = ldpc_distance_bound(LN_BN254, 16, 128);
    EXPECT_GT(at_2p32, 0.09);
    EXPECT_LT(at_bn254, 0.001);

    // Raising the degree recovers it, and the shipped preset is what the formula certifies.
    const BoltParams params = BoltParams::bn254();
    const double recovered = ldpc_distance_bound(LN_BN254, params.column_degree, params.check_degree);
    EXPECT_NEAR(recovered, params.gamma, 2e-3);
    // Beating Brakedown on distance is the whole point of paying for the larger degree.
    EXPECT_GT(params.gamma, brakedown::BrakedownParams::distance_7pct().distance());
}

TEST(BoltCodeTest, IsSystematicWithASketchTail)
{
    const size_t n = 1024;
    const BoltCode code(n, BoltParams::bn254());
    const std::vector<fr> message = random_array(n);
    const std::vector<fr> codeword = code.encode(message);

    ASSERT_EQ(codeword.size(), code.codeword_length());
    for (size_t i = 0; i < n; ++i) {
        EXPECT_EQ(codeword[i], message[i]) << "systematic prefix broken at " << i;
    }
    EXPECT_EQ(code.sketch_begin(), n);
    EXPECT_EQ(code.codeword_length(), n + code.sketch_codeword_length());
}

TEST(BoltCodeTest, EncodingIsLinear)
{
    const size_t n = 512;
    const BoltCode code(n, BoltParams::bn254());
    const std::vector<fr> a = random_array(n);
    const std::vector<fr> b = random_array(n);
    const fr scalar = fr::random_element();

    std::vector<fr> combination(n);
    for (size_t i = 0; i < n; ++i) {
        combination[i] = a[i] + (scalar * b[i]);
    }
    const std::vector<fr> encoded_a = code.encode(a);
    const std::vector<fr> encoded_b = code.encode(b);
    const std::vector<fr> encoded_combination = code.encode(combination);
    for (size_t i = 0; i < encoded_a.size(); ++i) {
        EXPECT_EQ(encoded_combination[i], encoded_a[i] + (scalar * encoded_b[i])) << "nonlinear at " << i;
    }
}

TEST(BoltCodeTest, SameSeedGivesTheSameCode)
{
    const std::vector<fr> message = random_array(256);
    const BoltCode first(256, BoltParams::bn254(), 2, /*seed=*/11);
    const BoltCode second(256, BoltParams::bn254(), 2, /*seed=*/11);
    const BoltCode other(256, BoltParams::bn254(), 2, /*seed=*/12);
    EXPECT_EQ(first.encode(message), second.encode(message));
    EXPECT_NE(first.encode(message), other.encode(message));
}

// H must be a parity check matrix: no low-weight message may sketch to zero, or the systematic piece
// of the distance guarantee collapses. Unit and near-unit vectors are the adversarial cases.
TEST(BoltCodeTest, LowWeightMessagesHaveNonzeroSketches)
{
    const size_t n = 1024;
    const BoltCode code(n, BoltParams::bn254());
    for (const size_t position : { size_t(0), size_t(1), n / 2, n - 1 }) {
        std::vector<fr> unit(n, fr::zero());
        unit[position] = fr::one();
        const std::vector<fr> codeword = code.encode(unit);
        bool sketch_nonzero = false;
        for (size_t i = code.sketch_begin(); i < codeword.size(); ++i) {
            sketch_nonzero = sketch_nonzero || !codeword[i].is_zero();
        }
        EXPECT_TRUE(sketch_nonzero) << "unit vector at " << position << " sketched to zero";
    }
}

TEST(BoltHonkTest, SingleUnshiftedCompleteness)
{
    const auto config = test_config(10);
    const CK ck(config);
    const auto instance = make_instance(ck, 1, 0);
    const auto proof = prove_instance(ck, instance);
    EXPECT_TRUE(verify_proof(config, instance.verifier_claims, instance.u, proof));
}

TEST(BoltHonkTest, BatchedWithShiftedCompleteness)
{
    const auto config = test_config(9);
    const CK ck(config);
    const auto instance = make_instance(ck, 3, 2);
    const auto proof = prove_instance(ck, instance);
    EXPECT_TRUE(verify_proof(config, instance.verifier_claims, instance.u, proof));
}

TEST(BoltHonkTest, WrongEvaluationRejected)
{
    const auto config = test_config(10);
    const CK ck(config);
    const auto instance = make_instance(ck, 2, 1);
    const auto proof = prove_instance(ck, instance);

    auto bad_claims = instance.verifier_claims;
    bad_claims.unshifted_evaluations[1] += fr(1);
    EXPECT_FALSE(verify_proof(config, bad_claims, instance.u, proof));

    bad_claims = instance.verifier_claims;
    bad_claims.shifted_evaluations[0] += fr(1);
    EXPECT_FALSE(verify_proof(config, bad_claims, instance.u, proof));
}

TEST(BoltHonkTest, TamperedProofRejected)
{
    const auto config = test_config(10);
    const CK ck(config);
    const auto instance = make_instance(ck, 2, 1);
    const auto proof = prove_instance(ck, instance);

    for (const size_t position : { size_t(1), proof.size() / 4, proof.size() / 2, proof.size() - 2 }) {
        HonkProof tampered = proof;
        tampered[position] += fr(1);
        EXPECT_FALSE(verify_proof(config, instance.verifier_claims, instance.u, tampered)) << "position " << position;
    }
}

// Both stretches must be sampled, each inside its own range: querying only one of them would leave
// the other piece of the piecewise guarantee untested.
TEST(BoltHonkTest, QueryPlanCoversBothStretchesSeparately)
{
    const auto config = test_config(10);
    const BoltCodePolicy code(config);
    const auto plan = code.query_plan(config);
    ASSERT_EQ(plan.size(), 2U);
    EXPECT_EQ(plan[0].begin, 0U);
    EXPECT_EQ(plan[0].length, code.message_length());
    EXPECT_EQ(plan[1].begin, code.message_length());
    EXPECT_EQ(plan[0].length + plan[1].length, code.codeword_length());
    EXPECT_GT(plan[0].num_queries, 0U);
    EXPECT_GT(plan[1].num_queries, 0U);
}

// The piecewise test is what buys Bolt its proof size: testing each stretch against its own distance
// costs far fewer queries than sampling uniformly against the diluted overall distance.
TEST(BoltHonkTest, PiecewiseQueriesBeatTheDilutedAlternative)
{
    const ligero::LigeroConfig config{ 20, 14, 2, 100, 0 };
    const size_t piecewise = BoltCodePolicy::num_queries(config);
    const size_t diluted = BoltCodePolicy::diluted_queries(config, BoltParams::bn254());
    EXPECT_LT(piecewise, diluted);
    // And it must also beat Brakedown, whose single distance is smaller still.
    EXPECT_LT(piecewise, brakedown::BrakedownCodePolicy::provable_num_queries(100, 0.07));
}

class BoltHonkCircuitTest : public ::testing::Test {
  public:
    using Honk = BoltHonk<whir::Blake3sMerkleHasher>;

    static UltraCircuitBuilder build_test_circuit()
    {
        UltraCircuitBuilder builder;
        MockCircuits::add_arithmetic_gates_with_public_inputs(builder, 16);
        MockCircuits::add_arithmetic_gates(builder, 16);
        MockCircuits::add_lookup_gates(builder, 2);
        const size_t rom_id = builder.create_ROM_array(4);
        for (size_t i = 0; i < 4; ++i) {
            builder.set_ROM_element(rom_id, i, builder.add_variable(fr(3 * i + 1)));
        }
        builder.read_ROM_array(rom_id, builder.add_variable(fr(2)));
        return builder;
    }
};

TEST_F(BoltHonkCircuitTest, ProveAndVerify)
{
    UltraCircuitBuilder sizing_builder = build_test_circuit();
    ProverInstance_<UltraFlavor> sizing_instance(sizing_builder);
    const size_t log_n = sizing_instance.log_dyadic_size();

    UltraCircuitBuilder proving_builder = build_test_circuit();
    const auto config = Honk::make_config(log_n, TEST_SECURITY_BITS, 2);
    auto pk = Honk::create_proving_key(proving_builder, config);
    const auto vk = pk.vk;
    const HonkProof proof = Honk::prove(pk);
    EXPECT_TRUE(Honk::verify(vk, config, proof));

    HonkProof tampered = proof;
    tampered[tampered.size() / 2] += fr(1);
    bool accepted = false;
    try {
        accepted = Honk::verify(vk, config, tampered);
    } catch (const std::exception&) {
        accepted = false;
    }
    EXPECT_FALSE(accepted);
}

} // namespace bb::bolt
