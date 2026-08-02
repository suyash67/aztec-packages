#include "barretenberg/commitment_schemes/whir/whir.hpp"

#include "barretenberg/transcript/transcript.hpp"

#include <gtest/gtest.h>

namespace bb::whir {

namespace {

std::vector<fr> random_array(size_t size)
{
    std::vector<fr> array(size);
    for (fr& value : array) {
        value = fr::random_element();
    }
    return array;
}

std::vector<fr> random_point(size_t size)
{
    return random_array(size);
}

fr mle(std::span<const fr> array, std::span<const fr> u)
{
    return Polynomial<fr>(array).evaluate_mle(u);
}

/** @brief MLE of the shifted array (a₁, ..., a_{n-1}, 0) at u — Honk's shifted-claim semantics. */
fr shifted_mle(std::span<const fr> array, std::span<const fr> u)
{
    std::vector<fr> shifted(array.begin() + 1, array.end());
    shifted.push_back(fr::zero());
    return mle(shifted, u);
}

} // namespace

template <typename Hasher> class WhirTest : public ::testing::Test {
  public:
    using CK = WhirCommitmentKey<Hasher>;
    using Prover = WhirProver<Hasher>;
    using Verifier = WhirVerifier<Hasher>;
    using GroupData = WhirGroupData<Hasher>;

    static WhirConfig test_config(size_t num_variables, bool zk = false)
    {
        return WhirConfig::create(num_variables,
                                  /*security_bits=*/64,
                                  /*log_inv_rate=*/2,
                                  /*folding_factor_bits=*/4,
                                  /*final_poly_bits=*/4,
                                  WhirSoundness::CONJECTURED_LIST,
                                  zk);
    }

    struct Instance {
        std::vector<std::vector<fr>> unshifted_arrays;
        std::vector<std::vector<fr>> to_be_shifted_arrays; // constant term zero
        std::vector<fr> u;
        std::vector<WhirGroupData<Hasher>> groups; // one group of unshifted columns, one of shifted
        typename Prover::Claims prover_claims;
        typename Verifier::Claims verifier_claims;
    };

    static Instance make_instance(const CK& ck, size_t num_unshifted, size_t num_shifted)
    {
        const size_t n = size_t(1) << ck.config.num_payload_variables;
        Instance instance;
        instance.u = random_point(ck.config.num_payload_variables);
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
            instance.groups.push_back(
                ck.commit_group(instance.to_be_shifted_arrays, std::vector<bool>(num_shifted, true)));
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
            instance.prover_claims.shifted_evaluations.push_back(
                shifted_mle(instance.to_be_shifted_arrays[c], instance.u));
        }
        instance.verifier_claims.unshifted = instance.prover_claims.unshifted;
        instance.verifier_claims.unshifted_evaluations = instance.prover_claims.unshifted_evaluations;
        instance.verifier_claims.to_be_shifted = instance.prover_claims.to_be_shifted;
        instance.verifier_claims.shifted_evaluations = instance.prover_claims.shifted_evaluations;
        return instance;
    }

    static HonkProof prove_instance(const CK& ck, const Instance& instance)
    {
        auto prover_transcript = NativeTranscript::test_prover_init_empty();
        Prover::prove(ck, instance.prover_claims, instance.u, prover_transcript);
        return prover_transcript->export_proof();
    }

    static bool verify_proof(const WhirConfig& config,
                             const typename Verifier::Claims& claims,
                             std::span<const fr> u,
                             const HonkProof& proof)
    {
        auto verifier_transcript = std::make_shared<NativeTranscript>(proof);
        [[maybe_unused]] auto init = verifier_transcript->template receive_from_prover<fr>("Init");
        return Verifier::verify(config, claims, u, verifier_transcript);
    }
};

using HasherTypes = ::testing::Types<Poseidon2MerkleHasher, Blake3sMerkleHasher>;
TYPED_TEST_SUITE(WhirTest, HasherTypes);

TYPED_TEST(WhirTest, SingleUnshiftedCompleteness)
{
    const WhirConfig config = TestFixture::test_config(10);
    typename TestFixture::CK ck(config);
    const auto instance = TestFixture::make_instance(ck, 1, 0);
    const auto proof = TestFixture::prove_instance(ck, instance);
    EXPECT_TRUE(TestFixture::verify_proof(config, instance.verifier_claims, instance.u, proof));
}

TYPED_TEST(WhirTest, BatchedWithShiftedCompleteness)
{
    const WhirConfig config = TestFixture::test_config(10);
    typename TestFixture::CK ck(config);
    const auto instance = TestFixture::make_instance(ck, 3, 2);
    const auto proof = TestFixture::prove_instance(ck, instance);
    EXPECT_TRUE(TestFixture::verify_proof(config, instance.verifier_claims, instance.u, proof));
}

// A configuration small enough that there are no fold-and-commit iterations: the final phase checks
// the batched virtual oracle (with shift scaling) directly against the clear polynomial.
TYPED_TEST(WhirTest, RepairedSoundnessCompleteness)
{
    const WhirConfig config =
        WhirConfig::create(10, /*security_bits=*/64, /*log_inv_rate=*/2, 4, 4, WhirSoundness::REPAIRED_LIST);
    typename TestFixture::CK ck(config);
    const auto instance = TestFixture::make_instance(ck, 2, 1);
    const auto proof = TestFixture::prove_instance(ck, instance);
    EXPECT_TRUE(TestFixture::verify_proof(config, instance.verifier_claims, instance.u, proof));
}

TYPED_TEST(WhirTest, ZeroIterationEdgeCase)
{
    const WhirConfig config = TestFixture::test_config(4);
    ASSERT_EQ(config.num_iterations(), 0U);
    typename TestFixture::CK ck(config);
    const auto instance = TestFixture::make_instance(ck, 2, 1);
    const auto proof = TestFixture::prove_instance(ck, instance);
    EXPECT_TRUE(TestFixture::verify_proof(config, instance.verifier_claims, instance.u, proof));
}

TYPED_TEST(WhirTest, WrongEvaluationRejected)
{
    const WhirConfig config = TestFixture::test_config(10);
    typename TestFixture::CK ck(config);
    const auto instance = TestFixture::make_instance(ck, 2, 1);
    const auto proof = TestFixture::prove_instance(ck, instance);

    auto bad_claims = instance.verifier_claims;
    bad_claims.unshifted_evaluations[1] += fr(1);
    EXPECT_FALSE(TestFixture::verify_proof(config, bad_claims, instance.u, proof));

    bad_claims = instance.verifier_claims;
    bad_claims.shifted_evaluations[0] += fr(1);
    EXPECT_FALSE(TestFixture::verify_proof(config, bad_claims, instance.u, proof));
}

TYPED_TEST(WhirTest, WrongOpeningPointRejected)
{
    const WhirConfig config = TestFixture::test_config(10);
    typename TestFixture::CK ck(config);
    const auto instance = TestFixture::make_instance(ck, 1, 0);
    const auto proof = TestFixture::prove_instance(ck, instance);

    auto wrong_u = instance.u;
    wrong_u[3] += fr(1);
    EXPECT_FALSE(TestFixture::verify_proof(config, instance.verifier_claims, wrong_u, proof));
}

TYPED_TEST(WhirTest, TamperedProofRejected)
{
    const WhirConfig config = TestFixture::test_config(10);
    typename TestFixture::CK ck(config);
    const auto instance = TestFixture::make_instance(ck, 1, 0);
    const auto proof = TestFixture::prove_instance(ck, instance);

    // Tamper with elements spread across the proof: roots, sumcheck messages, opened values, the
    // final polynomial.
    for (const size_t position : { size_t(1), proof.size() / 4, proof.size() / 2, proof.size() - 2 }) {
        HonkProof tampered = proof;
        tampered[position] += fr(1);
        EXPECT_FALSE(TestFixture::verify_proof(config, instance.verifier_claims, instance.u, tampered))
            << "position " << position;
    }
}

TYPED_TEST(WhirTest, ZkBatchedCompleteness)
{
    const WhirConfig config = TestFixture::test_config(8, /*zk=*/true);
    ASSERT_EQ(config.num_variables, 9U);
    typename TestFixture::CK ck(config);
    const auto instance = TestFixture::make_instance(ck, 2, 1);
    const auto proof = TestFixture::prove_instance(ck, instance);
    EXPECT_TRUE(TestFixture::verify_proof(config, instance.verifier_claims, instance.u, proof));
}

TYPED_TEST(WhirTest, ZkZeroIterationCompleteness)
{
    const WhirConfig config = TestFixture::test_config(6, /*zk=*/true);
    ASSERT_EQ(config.num_iterations(), 0U);
    typename TestFixture::CK ck(config);
    const auto instance = TestFixture::make_instance(ck, 2, 1);
    const auto proof = TestFixture::prove_instance(ck, instance);
    EXPECT_TRUE(TestFixture::verify_proof(config, instance.verifier_claims, instance.u, proof));
}

TYPED_TEST(WhirTest, ZkWrongEvaluationRejected)
{
    const WhirConfig config = TestFixture::test_config(8, /*zk=*/true);
    typename TestFixture::CK ck(config);
    const auto instance = TestFixture::make_instance(ck, 2, 1);
    const auto proof = TestFixture::prove_instance(ck, instance);

    auto bad_claims = instance.verifier_claims;
    bad_claims.shifted_evaluations[0] += fr(1);
    EXPECT_FALSE(TestFixture::verify_proof(config, bad_claims, instance.u, proof));
}

// Salted, blinded commitments to the same payload are distinct: the commitment is hiding.
TYPED_TEST(WhirTest, ZkCommitmentIsHiding)
{
    const WhirConfig config = TestFixture::test_config(8, /*zk=*/true);
    typename TestFixture::CK ck(config);
    const std::vector<fr> payload = random_array(size_t(1) << config.num_payload_variables);
    const auto data_a = ck.commit(std::span<const fr>(payload));
    const auto data_b = ck.commit(std::span<const fr>(payload));
    EXPECT_NE(data_a.tree.root(), data_b.tree.root());
}

TEST(WhirConfigTest, ZkSchedule)
{
    const WhirConfig config = WhirConfig::create(8, 64, 2, 4, 4, WhirSoundness::CONJECTURED_LIST, /*zk=*/true);
    EXPECT_TRUE(config.zk);
    EXPECT_EQ(config.num_variables, 9U);
    EXPECT_EQ(config.num_payload_variables, 8U);
    ASSERT_EQ(config.num_iterations(), 1U);
    // Blinding must cover every round-0 query with slack.
    EXPECT_EQ(config.num_blinding_coefficients, config.rounds[0].num_queries + 8);
}

TEST(WhirConfigTest, QueryCountFormulas)
{
    EXPECT_EQ(WhirConfig::compute_num_queries(100, 2, WhirSoundness::CONJECTURED_LIST), 50U);
    EXPECT_EQ(WhirConfig::compute_num_queries(100, 3, WhirSoundness::CONJECTURED_LIST), 34U);
    EXPECT_EQ(WhirConfig::compute_num_queries(100, 2, WhirSoundness::PROVABLE_LIST), 100U);
    // Unique decoding at rate 1/2 tests distance 1/4: -log2(3/4) ≈ 0.415 bits per query.
    const size_t ud_queries = WhirConfig::compute_num_queries(64, 1, WhirSoundness::UNIQUE_DECODING);
    EXPECT_GE(ud_queries, 154U);
    EXPECT_LE(ud_queries, 156U);
}

// The repaired-conjecture regime (Crites-Stewart, eprint 2025/2046) tests distance δ* with
// H_q(δ*) = 1-ρ instead of the disproved capacity δ = 1-ρ. At 254-bit q the entropy penalty
// h₂(δ)/log₂q costs about one extra query at aggressive rates and vanishes at high rates.
TEST(WhirConfigTest, RepairedQueryCountsExceedDisprovenCapacityCounts)
{
    // r=2: per-query error 1/4 + h₂(0.747)/253 ≈ 2^-1.98 -> 51 queries against 50.
    EXPECT_EQ(WhirConfig::compute_num_queries(100, 2, WhirSoundness::REPAIRED_LIST), 51U);
    // r=5: 2^-4.96 bits per query -> 21 against 20.
    EXPECT_EQ(WhirConfig::compute_num_queries(100, 5, WhirSoundness::REPAIRED_LIST), 21U);
    // r=8 and beyond: the penalty is below the ceiling granularity.
    EXPECT_EQ(WhirConfig::compute_num_queries(100, 8, WhirSoundness::REPAIRED_LIST), 13U);
    EXPECT_EQ(WhirConfig::compute_num_queries(100, 11, WhirSoundness::REPAIRED_LIST), 10U);
    // The repaired count is never below the capacity count and never above Johnson.
    for (size_t rate = 1; rate <= 20; ++rate) {
        const size_t repaired = WhirConfig::compute_num_queries(100, rate, WhirSoundness::REPAIRED_LIST);
        EXPECT_GE(repaired, WhirConfig::compute_num_queries(100, rate, WhirSoundness::CONJECTURED_LIST));
        EXPECT_LE(repaired, WhirConfig::compute_num_queries(100, rate, WhirSoundness::PROVABLE_LIST));
    }
}

// The README.md §6 worked example: m = 20, r₀ = 2, k = 4, λ = 100.
TEST(WhirConfigTest, ScheduleWorkedExample)
{
    const WhirConfig config = WhirConfig::create(20, 100, 2, 4, 4);
    ASSERT_EQ(config.num_iterations(), 4U);
    const std::vector<size_t> expected_rates = { 2, 5, 8, 11 };
    const std::vector<size_t> expected_queries = { 50, 20, 13, 10 };
    for (size_t i = 0; i < 4; ++i) {
        EXPECT_EQ(config.rounds[i].num_variables, 20 - 4 * i);
        EXPECT_EQ(config.rounds[i].log_domain_size, 22 - i);
        EXPECT_EQ(config.rounds[i].log_inv_rate, expected_rates[i]);
        EXPECT_EQ(config.rounds[i].num_queries, expected_queries[i]);
    }
    EXPECT_EQ(config.final_round.num_variables, 4U);
    EXPECT_EQ(config.final_round.log_domain_size, 18U);
    EXPECT_EQ(config.final_round.log_inv_rate, 14U);
    EXPECT_EQ(config.final_round.num_queries, 8U);
}

} // namespace bb::whir
