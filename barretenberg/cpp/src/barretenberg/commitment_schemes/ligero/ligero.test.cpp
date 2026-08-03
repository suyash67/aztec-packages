#include "barretenberg/commitment_schemes/ligero/ligero.hpp"

#include "barretenberg/transcript/transcript.hpp"

#include <gtest/gtest.h>

namespace bb::ligero {

namespace {

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

/** @brief MLE of the shifted array (a₁, ..., a_{n-1}, 0) at u — Honk's shifted-claim semantics. */
fr shifted_mle(std::span<const fr> array, std::span<const fr> u)
{
    std::vector<fr> shifted(array.begin() + 1, array.end());
    shifted.push_back(fr::zero());
    return mle(shifted, u);
}

} // namespace

template <typename Hasher> class LigeroTest : public ::testing::Test {
  public:
    using CK = LigeroCommitmentKey<Hasher>;
    using Prover = LigeroProver<Hasher>;
    using Verifier = LigeroVerifier<Hasher>;

    static LigeroConfig test_config(size_t num_variables)
    {
        return LigeroConfig::create(num_variables,
                                    /*security_bits=*/64,
                                    /*log_inv_rate=*/2,
                                    /*est_total_polynomials=*/4);
    }

    struct Instance {
        std::vector<std::vector<fr>> unshifted_arrays;
        std::vector<std::vector<fr>> to_be_shifted_arrays;
        std::vector<fr> u;
        std::vector<LigeroGroupData<Hasher>> groups;
        typename Prover::Claims prover_claims;
        typename Verifier::Claims verifier_claims;
    };

    static Instance make_instance(const CK& ck, size_t num_unshifted, size_t num_shifted)
    {
        const size_t n = size_t(1) << ck.config.num_variables;
        Instance instance;
        instance.u = random_array(ck.config.num_variables);
        for (size_t p = 0; p < num_unshifted; ++p) {
            instance.unshifted_arrays.push_back(random_array(n));
        }
        for (size_t p = 0; p < num_shifted; ++p) {
            instance.to_be_shifted_arrays.push_back(random_array(n));
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
        auto transcript = NativeTranscript::test_prover_init_empty();
        Prover::prove(ck, instance.prover_claims, instance.u, transcript);
        return transcript->export_proof();
    }

    static bool verify_proof(const LigeroConfig& config,
                             const typename Verifier::Claims& claims,
                             std::span<const fr> u,
                             const HonkProof& proof)
    {
        auto transcript = std::make_shared<NativeTranscript>(proof);
        [[maybe_unused]] auto init = transcript->template receive_from_prover<fr>("Init");
        const RSCodePolicy code(config);
        return Verifier::verify(config, claims, u, transcript, code);
    }
};

using HasherTypes = ::testing::Types<whir::Poseidon2MerkleHasher, whir::Blake3sMerkleHasher>;
TYPED_TEST_SUITE(LigeroTest, HasherTypes);

TYPED_TEST(LigeroTest, SingleUnshiftedCompleteness)
{
    const LigeroConfig config = TestFixture::test_config(10);
    typename TestFixture::CK ck(config);
    const auto instance = TestFixture::make_instance(ck, 1, 0);
    const auto proof = TestFixture::prove_instance(ck, instance);
    EXPECT_TRUE(TestFixture::verify_proof(config, instance.verifier_claims, instance.u, proof));
}

TYPED_TEST(LigeroTest, BatchedWithShiftedCompleteness)
{
    const LigeroConfig config = TestFixture::test_config(10);
    typename TestFixture::CK ck(config);
    const auto instance = TestFixture::make_instance(ck, 3, 2);
    const auto proof = TestFixture::prove_instance(ck, instance);
    EXPECT_TRUE(TestFixture::verify_proof(config, instance.verifier_claims, instance.u, proof));
}

TYPED_TEST(LigeroTest, WrongEvaluationRejected)
{
    const LigeroConfig config = TestFixture::test_config(10);
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

TYPED_TEST(LigeroTest, TamperedProofRejected)
{
    const LigeroConfig config = TestFixture::test_config(10);
    typename TestFixture::CK ck(config);
    const auto instance = TestFixture::make_instance(ck, 2, 1);
    const auto proof = TestFixture::prove_instance(ck, instance);

    for (const size_t position : { size_t(1), proof.size() / 4, proof.size() / 2, proof.size() - 2 }) {
        HonkProof tampered = proof;
        tampered[position] += fr(1);
        EXPECT_FALSE(TestFixture::verify_proof(config, instance.verifier_claims, instance.u, tampered))
            << "position " << position;
    }
}

TEST(LigeroConfigTest, DefaultSplit)
{
    const LigeroConfig config = LigeroConfig::create(20, 100, 2, 36);
    EXPECT_EQ(config.num_queries, 50U);
    // C* ~ sqrt(50 * 36 * 2^19) ~ 2^14.8 -> 2^15
    EXPECT_EQ(config.log_num_cols, 15U);
    EXPECT_EQ(config.num_rows(), size_t(1) << 5);
}

} // namespace bb::ligero
