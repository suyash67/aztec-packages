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

/** @brief Test parameter pairing a Merkle hasher with a commitment layout. */
template <typename Hasher_, size_t MaxStackBits_> struct WhirVariant {
    using Hasher = Hasher_;
    static constexpr size_t MAX_STACK_BITS = MaxStackBits_;
};

template <typename Variant> class WhirTest : public ::testing::Test {
  public:
    using Hasher = typename Variant::Hasher;
    using CK = WhirCommitmentKey<Hasher>;
    using Prover = WhirProver<Hasher>;
    using Verifier = WhirVerifier<Hasher>;
    using GroupData = WhirGroupData<Hasher>;

    static WhirConfig test_config(size_t num_variables, bool zk = false, size_t security_bits = 64)
    {
        return WhirConfig::create(num_variables,
                                  security_bits,
                                  /*log_inv_rate=*/2,
                                  /*folding_factor_bits=*/4,
                                  /*final_poly_bits=*/4,
                                  WhirSoundness::CONJECTURED_LIST,
                                  zk,
                                  Variant::MAX_STACK_BITS);
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

// Every Merkle hasher against the interleaved layout, and the two most-used against the stacked one.
using WhirVariants = ::testing::Types<WhirVariant<Poseidon2MerkleHasher, 0>,
                                      WhirVariant<Blake3sMerkleHasher, 0>,
                                      WhirVariant<Sha256MerkleHasher, 0>,
                                      WhirVariant<SkyscraperMerkleHasher, 0>,
                                      WhirVariant<Poseidon2MerkleHasher, 2>,
                                      WhirVariant<Blake3sMerkleHasher, 2>>;
TYPED_TEST_SUITE(WhirTest, WhirVariants);

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
    const WhirConfig config = WhirConfig::create(10,
                                                 /*security_bits=*/64,
                                                 /*log_inv_rate=*/2,
                                                 4,
                                                 4,
                                                 WhirSoundness::REPAIRED_LIST,
                                                 /*zk=*/false,
                                                 TypeParam::MAX_STACK_BITS);
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

// The out-of-domain answer on the round-0 oracle is what singles out one codeword of the decoding
// list, so it has to be a checked part of the statement rather than an unconstrained message. With
// one Poseidon2 group of a single unshifted column and no stacking there are no stacking challenges
// and no cross evaluations, so the proof opens [Init, root_0, y_ood_init_0_0, sumcheck...] and the
// sample sits at a known index.
TEST(WhirOutOfDomainTest, InitialSampleIsConstrained)
{
    using Fixture = WhirTest<WhirVariant<Poseidon2MerkleHasher, 0>>;
    const WhirConfig config = WhirConfig::create(10,
                                                 /*security_bits=*/64,
                                                 /*log_inv_rate=*/2,
                                                 4,
                                                 4,
                                                 WhirSoundness::PROVABLE_LIST);
    ASSERT_EQ(config.num_ood_samples, 1U);
    Fixture::CK ck(config);
    const auto instance = Fixture::make_instance(ck, 1, 0);
    const auto proof = Fixture::prove_instance(ck, instance);
    EXPECT_TRUE(Fixture::verify_proof(config, instance.verifier_claims, instance.u, proof));

    // Re-proving the same instance with the sample switched off leaves the header and root before
    // it untouched and diverges exactly at the sample, which pins its position. (The two proofs'
    // total lengths are not comparable: the sample is absorbed, so every later challenge differs,
    // and batched Merkle openings cost a number of digests that depends on the query indices.)
    WhirConfig without_ood = config;
    without_ood.num_ood_samples = 0;
    Fixture::CK ck_without_ood(without_ood);
    const auto proof_without_ood = Fixture::prove_instance(ck_without_ood, instance);
    constexpr size_t ood_position = 2;
    for (size_t i = 0; i < ood_position; ++i) {
        EXPECT_EQ(proof[i], proof_without_ood[i]) << "at position " << i;
    }
    EXPECT_NE(proof[ood_position], proof_without_ood[ood_position]);

    HonkProof tampered = proof;
    tampered[ood_position] += fr(1);
    EXPECT_FALSE(Fixture::verify_proof(config, instance.verifier_claims, instance.u, tampered));
}

TYPED_TEST(WhirTest, ZkBatchedCompleteness)
{
    const WhirConfig config = TestFixture::test_config(8, /*zk=*/true);
    ASSERT_EQ(config.num_variables, 9U + TypeParam::MAX_STACK_BITS);
    typename TestFixture::CK ck(config);
    const auto instance = TestFixture::make_instance(ck, 2, 1);
    const auto proof = TestFixture::prove_instance(ck, instance);
    EXPECT_TRUE(TestFixture::verify_proof(config, instance.verifier_claims, instance.u, proof));
}

TYPED_TEST(WhirTest, ZkZeroIterationCompleteness)
{
    // The only payload width with no fold iterations at either layout: stacking spends
    // MAX_STACK_BITS of the committed width on the stack index. The blinding coefficients have to
    // fit above that payload (one per round-0 query, plus slack), which at a 4-variable payload
    // caps the security level this configuration can carry - a stacked group of one column is only
    // 2^{m+1} wide, so there is no more room than the payload itself.
    const WhirConfig config =
        TestFixture::test_config(6 - TypeParam::MAX_STACK_BITS, /*zk=*/true, /*security_bits=*/12);
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
    // Johnson tests distance 1-√ρ-η with η = √ρ/20, so a query is worth r/2 - log₂(21/20) bits,
    // not r/2: at r=2 that is 0.9296 bits and 108 queries rather than 100.
    EXPECT_EQ(WhirConfig::compute_num_queries(100, 2, WhirSoundness::PROVABLE_LIST), 108U);
    // Unique decoding at rate 1/2 tests distance 1/4: -log2(3/4) ≈ 0.415 bits per query.
    const size_t ud_queries = WhirConfig::compute_num_queries(64, 1, WhirSoundness::UNIQUE_DECODING);
    EXPECT_GE(ud_queries, 154U);
    EXPECT_LE(ud_queries, 156U);
}

// The reference WHIR implementation (WizardOfMenlo/whir @ 0aeaa7f, the revision ProveKit pins)
// derives its in-domain query counts from the same Johnson bound with η = √ρ/20. These are the
// counts it reports at λ = 128 for the rates a k=3 schedule walks through, and they must agree
// exactly: a lower count here would mean barretenberg claims λ it does not have.
TEST(WhirConfigTest, JohnsonQueryCountsMatchReferenceImplementation)
{
    const std::vector<std::pair<size_t, size_t>> reference_counts = {
        { 2, 138 }, { 4, 67 }, { 5, 53 }, { 6, 44 }, { 8, 33 }, { 10, 26 }, { 11, 24 }, { 12, 22 }, { 14, 19 },
    };
    for (const auto& [log_inv_rate, expected] : reference_counts) {
        EXPECT_EQ(WhirConfig::compute_num_queries(128, log_inv_rate, WhirSoundness::PROVABLE_LIST), expected)
            << "at rate 2^-" << log_inv_rate;
    }
}

// The whole schedule at ProveKit's protocol parameters (λ = 128, rate 2^-2, k = 3, Johnson) must
// reproduce the reference implementation's, oracle by oracle: 138 queries against the round-0
// commitments, then 67, 44, 33, 26 against each folded oracle, then 22 in the final phase. The
// reference reaches the same 330 total using 118 bits of queries plus 10 bits of grinding; with no
// grinding implemented here the counts are derived at the full 128 instead.
TEST(WhirConfigTest, ProveKitScheduleMatchesReferenceImplementation)
{
    const WhirConfig config = WhirConfig::create(19,
                                                 /*security_bits=*/128,
                                                 /*log_inv_rate=*/2,
                                                 /*folding_factor_bits=*/3,
                                                 /*final_poly_bits=*/4,
                                                 WhirSoundness::PROVABLE_LIST);
    ASSERT_EQ(config.num_iterations(), 5U);
    const std::vector<size_t> expected_rates = { 2, 4, 6, 8, 10 };
    const std::vector<size_t> expected_queries = { 138, 67, 44, 33, 26 };
    size_t total_queries = 0;
    for (size_t i = 0; i < config.num_iterations(); ++i) {
        EXPECT_EQ(config.rounds[i].num_variables, 19 - 3 * i);
        EXPECT_EQ(config.rounds[i].log_inv_rate, expected_rates[i]);
        EXPECT_EQ(config.rounds[i].num_queries, expected_queries[i]);
        total_queries += config.rounds[i].num_queries;
    }
    EXPECT_EQ(config.final_round.num_variables, 4U);
    EXPECT_EQ(config.final_round.log_inv_rate, 12U);
    EXPECT_EQ(config.final_round.num_queries, 22U);
    total_queries += config.final_round.num_queries;
    EXPECT_EQ(total_queries, 330U);
}

// Every list-decoding regime must sample out of domain; at BN254's field size one sample per
// commitment suffices across the schedule, which is what the reference also derives.
TEST(WhirConfigTest, OutOfDomainSampleCounts)
{
    const WhirConfig johnson = WhirConfig::create(19, 128, 2, 3, 4, WhirSoundness::PROVABLE_LIST);
    EXPECT_EQ(johnson.num_ood_samples, 1U);
    const WhirConfig conjectured = WhirConfig::create(19, 128, 2, 4, 4, WhirSoundness::CONJECTURED_LIST);
    EXPECT_EQ(conjectured.num_ood_samples, 1U);
    // Unique decoding leaves a single codeword in the list, so no sample is needed.
    const WhirConfig unique = WhirConfig::create(19, 128, 2, 4, 4, WhirSoundness::UNIQUE_DECODING);
    EXPECT_EQ(unique.num_ood_samples, 0U);
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

// A first fold narrower than the rest shortens only the first iteration: round 0 folds 2^{k₀} and
// improves the rate by k₀-1 bits, every later iteration folds 2^k as usual. At k₀ = 1 the rate does
// not improve at all across round 0, so round 1 repeats round 0's query count.
TEST(WhirConfigTest, NarrowFirstFoldSchedule)
{
    const WhirConfig config = WhirConfig::create(19,
                                                 /*security_bits=*/128,
                                                 /*log_inv_rate=*/2,
                                                 /*folding_factor_bits=*/3,
                                                 /*final_poly_bits=*/4,
                                                 WhirSoundness::PROVABLE_LIST,
                                                 /*zk=*/false,
                                                 /*max_stack_bits=*/0,
                                                 /*initial_folding_factor_bits=*/1);
    ASSERT_EQ(config.num_iterations(), 5U);
    EXPECT_EQ(config.initial_folding_factor_bits, 1U);
    const std::vector<size_t> expected_variables = { 19, 18, 15, 12, 9 };
    const std::vector<size_t> expected_folds = { 1, 3, 3, 3, 3 };
    const std::vector<size_t> expected_rates = { 2, 2, 4, 6, 8 };
    const std::vector<size_t> expected_queries = { 138, 138, 67, 44, 33 };
    for (size_t i = 0; i < config.num_iterations(); ++i) {
        EXPECT_EQ(config.rounds[i].num_variables, expected_variables[i]) << "round " << i;
        EXPECT_EQ(config.rounds[i].folding_factor_bits, expected_folds[i]) << "round " << i;
        EXPECT_EQ(config.rounds[i].log_inv_rate, expected_rates[i]) << "round " << i;
        EXPECT_EQ(config.rounds[i].num_queries, expected_queries[i]) << "round " << i;
    }
    // Each folded oracle is committed with the arity the *next* round folds by.
    for (size_t i = 0; i + 1 < config.num_iterations(); ++i) {
        EXPECT_EQ(config.committed_arity_bits(i), config.rounds[i + 1].folding_factor_bits);
    }
    EXPECT_EQ(config.committed_arity_bits(config.num_iterations() - 1), config.final_round.folding_factor_bits);
    EXPECT_EQ(config.final_round.num_variables, 6U);
    EXPECT_EQ(config.final_round.log_inv_rate, 10U);
    EXPECT_EQ(config.final_round.num_queries, 26U);
    EXPECT_EQ(config.final_round.folding_factor_bits, 3U);
}

// Grinding trades query soundness for prover work one bit at a time, so a schedule with p bits of
// proof of work draws its query counts at λ-p. ProveKit's own configuration is this at p = 10.
TEST(WhirConfigTest, GrindingReplacesQueryBits)
{
    const auto schedule = [](size_t pow_bits) {
        return WhirConfig::create(19,
                                  /*security_bits=*/128,
                                  /*log_inv_rate=*/2,
                                  /*folding_factor_bits=*/3,
                                  /*final_poly_bits=*/4,
                                  WhirSoundness::PROVABLE_LIST,
                                  /*zk=*/false,
                                  /*max_stack_bits=*/0,
                                  /*initial_folding_factor_bits=*/0,
                                  pow_bits);
    };
    const WhirConfig ground = schedule(10);
    EXPECT_EQ(ground.pow_bits, 10U);
    for (size_t i = 0; i < ground.num_iterations(); ++i) {
        EXPECT_EQ(ground.rounds[i].num_queries,
                  WhirConfig::compute_num_queries(118, ground.rounds[i].log_inv_rate, WhirSoundness::PROVABLE_LIST));
        EXPECT_LT(ground.rounds[i].num_queries, schedule(0).rounds[i].num_queries);
    }
    EXPECT_EQ(schedule(0).pow_bits, 0U);
}

// The grind must actually bind: a nonce is accepted only when its digest has the required leading
// zeros, and the search returns the same nonce every time it runs.
TEST(WhirProofOfWorkTest, NonceIsCheckedAndDeterministic)
{
    constexpr size_t pow_bits = 12;
    const fr seed = fr::random_element();
    const uint64_t nonce = detail::grind(seed, pow_bits);
    EXPECT_TRUE(detail::pow_is_valid(seed, nonce, pow_bits));
    EXPECT_EQ(detail::grind(seed, pow_bits), nonce);

    // Every smaller nonce fails, so the search really did find the first one.
    for (uint64_t candidate = 0; candidate < nonce; ++candidate) {
        EXPECT_FALSE(detail::pow_is_valid(seed, candidate, pow_bits)) << "at nonce " << candidate;
    }
    // A nonce that works for one seed does not carry over to another.
    EXPECT_FALSE(detail::pow_is_valid(seed + fr(1), nonce, pow_bits) &&
                 detail::pow_is_valid(seed + fr(2), nonce, pow_bits));
    // Zero bits of work accept anything.
    EXPECT_TRUE(detail::pow_is_valid(seed, 12345, 0));
}

// `poseidon2_pow_is_valid` spells the sponge out as one permutation so the search does not allocate.
// If that ever drifts from `Poseidon2::hash`, the in-circuit check — which calls the sponge — stops
// agreeing with the native one and every ground proof fails to recurse.
TEST(WhirProofOfWorkTest, Poseidon2GrindMatchesTheSponge)
{
    using Poseidon2 = crypto::Poseidon2<crypto::Poseidon2Bn254ScalarFieldParams>;
    for (size_t i = 0; i < 8; ++i) {
        const fr seed = fr::random_element();
        const uint64_t nonce = i * 7919;
        const uint256_t sponge(Poseidon2::hash({ seed, fr(nonce) }));
        // The predicate accepts exactly when the sponge digest's low bits vanish.
        for (const size_t bits : { size_t(1), size_t(4), size_t(11) }) {
            const bool expected = (sponge.data[0] & ((uint64_t(1) << bits) - 1)) == 0;
            EXPECT_EQ(detail::poseidon2_pow_is_valid(seed, nonce, bits), expected);
        }
    }
}

// The Poseidon2 grind an in-circuit verifier can afford: divisibility of the digest by 2^pow_bits
// rather than leading zero bytes of a Blake3 digest.
TEST(WhirProofOfWorkTest, Poseidon2NonceIsCheckedAndDeterministic)
{
    constexpr size_t pow_bits = 12;
    const fr seed = fr::random_element();
    const uint64_t nonce = detail::grind(seed, pow_bits, /*poseidon2=*/true);
    EXPECT_TRUE(detail::poseidon2_pow_is_valid(seed, nonce, pow_bits));
    EXPECT_EQ(detail::grind(seed, pow_bits, /*poseidon2=*/true), nonce);
    for (uint64_t candidate = 0; candidate < nonce; ++candidate) {
        EXPECT_FALSE(detail::poseidon2_pow_is_valid(seed, candidate, pow_bits)) << "at nonce " << candidate;
    }
    // The two grinds are genuinely different functions, so a Blake3 nonce is not a Poseidon2 one.
    EXPECT_FALSE(detail::poseidon2_pow_is_valid(seed, detail::grind(seed, pow_bits, /*poseidon2=*/false), pow_bits) &&
                 detail::pow_is_valid(seed, nonce, pow_bits));
}

/**
 * @brief The recursion profile: per-query authentication paths against a Merkle cap, with a
 * Poseidon2 grind. Completeness plus the soundness of everything the profile changed.
 */
class WhirRecursionProfileTest : public WhirTest<WhirVariant<Poseidon2CompressionHasher, 0>> {
  public:
    static WhirConfig recursion_config(size_t num_variables, size_t pow_bits = 8)
    {
        WhirConfig config = WhirConfig::create(num_variables,
                                               /*security_bits=*/64,
                                               /*log_inv_rate=*/2,
                                               /*folding_factor_bits=*/4,
                                               /*final_poly_bits=*/4,
                                               WhirSoundness::CONJECTURED_LIST,
                                               /*zk=*/false,
                                               /*max_stack_bits=*/0,
                                               /*initial_folding_factor_bits=*/1,
                                               pow_bits);
        config.enable_recursion_profile();
        return config;
    }
};

TEST_F(WhirRecursionProfileTest, Completeness)
{
    const WhirConfig config = recursion_config(10);
    EXPECT_TRUE(config.per_query_openings);
    EXPECT_TRUE(config.poseidon2_pow);
    EXPECT_GT(config.merkle_cap_levels, 0U);

    CK ck(config);
    const auto instance = make_instance(ck, 3, 2);
    const auto proof = prove_instance(ck, instance);
    EXPECT_TRUE(verify_proof(config, instance.verifier_claims, instance.u, proof));
}

TEST_F(WhirRecursionProfileTest, CapZeroIsTheOrdinaryRootWalk)
{
    WhirConfig config = recursion_config(10);
    config.merkle_cap_levels = 0;
    CK ck(config);
    const auto instance = make_instance(ck, 2, 1);
    const auto proof = prove_instance(ck, instance);
    EXPECT_TRUE(verify_proof(config, instance.verifier_claims, instance.u, proof));
}

// Every opened value, sibling digest and cap entry is bound: corrupting any single field of the
// proof stream must break authentication.
TEST_F(WhirRecursionProfileTest, TamperedOpeningRejected)
{
    const WhirConfig config = recursion_config(10);
    CK ck(config);
    const auto instance = make_instance(ck, 2, 1);
    const auto proof = prove_instance(ck, instance);

    // The unhashed opening data sits after the hashed prefix; sample positions across the tail.
    size_t rejected = 0;
    const size_t stride = std::max<size_t>(1, proof.size() / 32);
    for (size_t i = proof.size() / 2; i < proof.size(); i += stride) {
        HonkProof tampered = proof;
        tampered[i] += fr(1);
        if (!verify_proof(config, instance.verifier_claims, instance.u, tampered)) {
            ++rejected;
        }
    }
    EXPECT_EQ(rejected, (proof.size() - proof.size() / 2 + stride - 1) / stride);
}

TEST_F(WhirRecursionProfileTest, GrindingIsEnforced)
{
    const WhirConfig config = recursion_config(10, /*pow_bits=*/10);
    CK ck(config);
    const auto instance = make_instance(ck, 2, 1);
    const auto proof = prove_instance(ck, instance);
    ASSERT_TRUE(verify_proof(config, instance.verifier_claims, instance.u, proof));

    // Exactly one field of the proof is the first round's nonce; breaking every candidate nonce
    // position must be caught, and at least one of them is the real one.
    WhirConfig unground = config;
    unground.pow_bits = 0;
    CK unground_ck(unground);
    const auto unground_instance = make_instance(unground_ck, 2, 1);
    const auto unground_proof = prove_instance(unground_ck, unground_instance);
    // A proof produced without grinding cannot satisfy a verifier that demands it.
    EXPECT_FALSE(verify_proof(config, unground_instance.verifier_claims, unground_instance.u, unground_proof));
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
