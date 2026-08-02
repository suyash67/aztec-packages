#include "barretenberg/commitment_schemes/mercury/chopin_honk.hpp"

#include "barretenberg/srs/global_crs.hpp"
#include "barretenberg/stdlib_circuit_builders/mock_circuits.hpp"
#include "barretenberg/transcript/transcript.hpp"

#include <gtest/gtest.h>

namespace bb::chopin {

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

class ChopinTest : public ::testing::Test {
  public:
    static void SetUpTestSuite() { srs::init_file_crs_factory(srs::bb_crs_path()); }

    struct Instance {
        std::vector<std::vector<fr>> unshifted_arrays;
        std::vector<std::vector<fr>> to_be_shifted_arrays;
        std::vector<fr> u;
        std::vector<ChopinGroupData> groups;
        ChopinProver::Claims prover_claims;
        ChopinVerifier::Claims verifier_claims;
    };

    static Instance make_instance(const ChopinCommitmentKey& ck, size_t num_unshifted, size_t num_shifted)
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

    static HonkProof prove_instance(const ChopinCommitmentKey& ck, const Instance& instance)
    {
        auto transcript = NativeTranscript::test_prover_init_empty();
        ChopinProver::prove(ck, instance.prover_claims, instance.u, transcript);
        return transcript->export_proof();
    }

    static bool verify_proof(const ChopinConfig& config,
                             const ChopinVerifier::Claims& claims,
                             std::span<const fr> u,
                             const HonkProof& proof)
    {
        try {
            auto transcript = std::make_shared<NativeTranscript>(proof);
            [[maybe_unused]] auto init = transcript->template receive_from_prover<fr>("Init");
            return ChopinVerifier::verify(config, claims, u, transcript);
        } catch (const std::exception&) {
            return false;
        }
    }
};

TEST_F(ChopinTest, SingleUnshiftedCompleteness)
{
    const ChopinConfig config = ChopinConfig::create(10);
    ChopinCommitmentKey ck(config);
    const auto instance = make_instance(ck, 1, 0);
    const auto proof = prove_instance(ck, instance);
    EXPECT_TRUE(verify_proof(config, instance.verifier_claims, instance.u, proof));
}

TEST_F(ChopinTest, OddVariablesBatchedWithShiftedCompleteness)
{
    const ChopinConfig config = ChopinConfig::create(9);
    ChopinCommitmentKey ck(config);
    const auto instance = make_instance(ck, 3, 2);
    const auto proof = prove_instance(ck, instance);
    EXPECT_TRUE(verify_proof(config, instance.verifier_claims, instance.u, proof));
}

TEST_F(ChopinTest, WrongEvaluationRejected)
{
    const ChopinConfig config = ChopinConfig::create(10);
    ChopinCommitmentKey ck(config);
    const auto instance = make_instance(ck, 2, 1);
    const auto proof = prove_instance(ck, instance);

    auto bad_claims = instance.verifier_claims;
    bad_claims.unshifted_evaluations[1] += fr(1);
    EXPECT_FALSE(verify_proof(config, bad_claims, instance.u, proof));

    bad_claims = instance.verifier_claims;
    bad_claims.shifted_evaluations[0] += fr(1);
    EXPECT_FALSE(verify_proof(config, bad_claims, instance.u, proof));
}

TEST_F(ChopinTest, TamperedProofRejected)
{
    const ChopinConfig config = ChopinConfig::create(10);
    ChopinCommitmentKey ck(config);
    const auto instance = make_instance(ck, 2, 1);
    const auto proof = prove_instance(ck, instance);

    for (const size_t position : { size_t(1), proof.size() / 4, proof.size() / 2, proof.size() - 2 }) {
        HonkProof tampered = proof;
        tampered[position] += fr(1);
        EXPECT_FALSE(verify_proof(config, instance.verifier_claims, instance.u, tampered)) << "position " << position;
    }
}

// The opening payload is constant-size in N: 10 commitments and 15 field elements for the
// two-chain claim set (Mercury's is larger and grows with the BDFG opening-set structure).
TEST_F(ChopinTest, ProofSizeIsConstant)
{
    constexpr size_t frs_per_commitment = NativeTranscript::Codec::calc_num_fields<Commitment>();
    constexpr size_t frs_init = 1;
    size_t sizes[2] = { 0, 0 };
    size_t index = 0;
    for (const size_t num_variables : { size_t(8), size_t(10) }) {
        const ChopinConfig config = ChopinConfig::create(num_variables);
        ChopinCommitmentKey ck(config);
        const auto instance = make_instance(ck, 2, 1);
        const auto proof = prove_instance(ck, instance);
        ASSERT_TRUE(verify_proof(config, instance.verifier_claims, instance.u, proof));
        sizes[index++] = proof.size();
    }
    EXPECT_EQ(sizes[0], sizes[1]);
    // Roots: 3 columns. Opening: C0_0, C1_0, C0_1, C0s_1, C1_1, S, W, pi_s, pi_q, Q1, Q2 = 11
    // commitments; fields: a_0, a_1, a2_1, w20_1, 4x2 beta evaluations, s_beta, s_binv, y = 15.
    EXPECT_EQ(sizes[0], frs_init + (3 + 11) * frs_per_commitment + 15);
}

class ChopinHonkTest : public ::testing::Test {
  public:
    using Honk = ChopinHonk;

    static void SetUpTestSuite() { srs::init_file_crs_factory(srs::bb_crs_path()); }

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

TEST_F(ChopinHonkTest, ProveAndVerify)
{
    UltraCircuitBuilder sizing_builder = build_test_circuit();
    ProverInstance_<UltraFlavor> sizing_instance(sizing_builder);
    const size_t log_n = sizing_instance.log_dyadic_size();

    UltraCircuitBuilder proving_builder = build_test_circuit();
    const ChopinConfig config = Honk::make_config(log_n);
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

} // namespace bb::chopin
