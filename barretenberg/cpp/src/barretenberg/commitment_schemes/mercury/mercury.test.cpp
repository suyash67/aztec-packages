#include "barretenberg/commitment_schemes/mercury/mercury.hpp"

#include "barretenberg/srs/global_crs.hpp"
#include "barretenberg/transcript/transcript.hpp"

#include <gtest/gtest.h>

namespace bb::mercury {

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

fr shifted_mle(std::span<const fr> array, std::span<const fr> u)
{
    std::vector<fr> shifted(array.begin() + 1, array.end());
    shifted.push_back(fr::zero());
    return mle(shifted, u);
}

} // namespace

class MercuryTest : public ::testing::Test {
  public:
    struct Instance {
        std::vector<std::vector<fr>> unshifted_arrays;
        std::vector<std::vector<fr>> to_be_shifted_arrays;
        std::vector<fr> u;
        std::vector<MercuryGroupData> groups;
        MercuryProver::Claims prover_claims;
        MercuryVerifier::Claims verifier_claims;
    };

    static Instance make_instance(const MercuryCommitmentKey& ck, size_t num_unshifted, size_t num_shifted)
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
            instance.prover_claims.shifted_evaluations.push_back(
                shifted_mle(instance.to_be_shifted_arrays[c], instance.u));
        }
        instance.verifier_claims.unshifted = instance.prover_claims.unshifted;
        instance.verifier_claims.unshifted_evaluations = instance.prover_claims.unshifted_evaluations;
        instance.verifier_claims.to_be_shifted = instance.prover_claims.to_be_shifted;
        instance.verifier_claims.shifted_evaluations = instance.prover_claims.shifted_evaluations;
        return instance;
    }

    static HonkProof prove_instance(const MercuryCommitmentKey& ck, const Instance& instance)
    {
        auto transcript = NativeTranscript::test_prover_init_empty();
        MercuryProver::prove(ck, instance.prover_claims, instance.u, transcript);
        return transcript->export_proof();
    }

    static bool verify_proof(const MercuryConfig& config,
                             const MercuryVerifier::Claims& claims,
                             std::span<const fr> u,
                             const HonkProof& proof)
    {
        try {
            auto transcript = std::make_shared<NativeTranscript>(proof);
            [[maybe_unused]] auto init = transcript->template receive_from_prover<fr>("Init");
            return MercuryVerifier::verify(config, claims, u, transcript);
        } catch (const std::exception&) {
            // Malformed proof data (e.g. a tampered commitment that is not on the curve) is a
            // rejection, not a crash.
            return false;
        }
    }

  protected:
    static void SetUpTestSuite() { srs::init_file_crs_factory(srs::bb_crs_path()); }
};

TEST_F(MercuryTest, SingleUnshiftedCompleteness)
{
    const MercuryConfig config = MercuryConfig::create(10);
    MercuryCommitmentKey ck(config);
    const auto instance = make_instance(ck, 1, 0);
    const auto proof = prove_instance(ck, instance);
    EXPECT_TRUE(verify_proof(config, instance.verifier_claims, instance.u, proof));
}

TEST_F(MercuryTest, OddVariablesCompleteness)
{
    const MercuryConfig config = MercuryConfig::create(9);
    MercuryCommitmentKey ck(config);
    const auto instance = make_instance(ck, 2, 0);
    const auto proof = prove_instance(ck, instance);
    EXPECT_TRUE(verify_proof(config, instance.verifier_claims, instance.u, proof));
}

TEST_F(MercuryTest, BatchedWithShiftedCompleteness)
{
    const MercuryConfig config = MercuryConfig::create(10);
    MercuryCommitmentKey ck(config);
    const auto instance = make_instance(ck, 3, 2);
    const auto proof = prove_instance(ck, instance);
    EXPECT_TRUE(verify_proof(config, instance.verifier_claims, instance.u, proof));
}

TEST_F(MercuryTest, WrongEvaluationRejected)
{
    const MercuryConfig config = MercuryConfig::create(10);
    MercuryCommitmentKey ck(config);
    const auto instance = make_instance(ck, 2, 1);
    const auto proof = prove_instance(ck, instance);

    auto bad_claims = instance.verifier_claims;
    bad_claims.unshifted_evaluations[1] += fr(1);
    EXPECT_FALSE(verify_proof(config, bad_claims, instance.u, proof));

    bad_claims = instance.verifier_claims;
    bad_claims.shifted_evaluations[0] += fr(1);
    EXPECT_FALSE(verify_proof(config, bad_claims, instance.u, proof));
}

TEST_F(MercuryTest, TamperedProofRejected)
{
    const MercuryConfig config = MercuryConfig::create(10);
    MercuryCommitmentKey ck(config);
    const auto instance = make_instance(ck, 2, 1);
    const auto proof = prove_instance(ck, instance);

    for (const size_t position : { size_t(1), proof.size() / 4, proof.size() / 2, proof.size() - 2 }) {
        HonkProof tampered = proof;
        tampered[position] += fr(1);
        EXPECT_FALSE(verify_proof(config, instance.verifier_claims, instance.u, tampered)) << "position " << position;
    }
}

} // namespace bb::mercury
