#include "barretenberg/commitment_schemes/switchfold/switchfold_honk.hpp"

#include "barretenberg/commitment_schemes/ligero/ligero_honk.hpp"
#include "barretenberg/stdlib_circuit_builders/mock_circuits.hpp"
#include "barretenberg/transcript/transcript.hpp"

#include <gtest/gtest.h>

namespace bb::switchfold {

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

template <typename Hasher> class SwitchFoldTest : public ::testing::Test {
  public:
    using CK = SwitchFoldCommitmentKey<Hasher>;
    using Prover = SwitchFoldProver<Hasher>;
    using Verifier = SwitchFoldVerifier<Hasher>;

    static SwitchFoldConfig test_config(size_t num_variables)
    {
        return SwitchFoldConfig::create(num_variables, /*security_bits=*/64, /*log_inv_rate=*/2);
    }

    struct Instance {
        std::vector<std::vector<fr>> unshifted_arrays;
        std::vector<std::vector<fr>> to_be_shifted_arrays;
        std::vector<fr> u;
        std::vector<SwitchFoldGroupData<Hasher>> groups;
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

    static HonkProof prove_instance(const CK& ck, const Instance& instance)
    {
        auto transcript = NativeTranscript::test_prover_init_empty();
        Prover::prove(ck, instance.prover_claims, instance.u, transcript);
        return transcript->export_proof();
    }

    static bool verify_proof(const SwitchFoldConfig& config,
                             const typename Verifier::Claims& claims,
                             std::span<const fr> u,
                             const HonkProof& proof)
    {
        try {
            auto transcript = std::make_shared<NativeTranscript>(proof);
            [[maybe_unused]] auto init = transcript->template receive_from_prover<fr>("Init");
            return Verifier::verify(config, claims, u, transcript);
        } catch (const std::exception&) {
            return false;
        }
    }
};

using HasherTypes = ::testing::Types<whir::Poseidon2MerkleHasher, whir::Blake3sMerkleHasher>;
TYPED_TEST_SUITE(SwitchFoldTest, HasherTypes);

TYPED_TEST(SwitchFoldTest, SingleUnshiftedCompleteness)
{
    const SwitchFoldConfig config = TestFixture::test_config(10);
    typename TestFixture::CK ck(config);
    const auto instance = TestFixture::make_instance(ck, 1, 0);
    const auto proof = TestFixture::prove_instance(ck, instance);
    EXPECT_TRUE(TestFixture::verify_proof(config, instance.verifier_claims, instance.u, proof));
}

TYPED_TEST(SwitchFoldTest, BatchedWithShiftedCompleteness)
{
    const SwitchFoldConfig config = TestFixture::test_config(9);
    typename TestFixture::CK ck(config);
    const auto instance = TestFixture::make_instance(ck, 3, 2);
    const auto proof = TestFixture::prove_instance(ck, instance);
    EXPECT_TRUE(TestFixture::verify_proof(config, instance.verifier_claims, instance.u, proof));
}

TYPED_TEST(SwitchFoldTest, WrongEvaluationRejected)
{
    const SwitchFoldConfig config = TestFixture::test_config(10);
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

TYPED_TEST(SwitchFoldTest, TamperedProofRejected)
{
    const SwitchFoldConfig config = TestFixture::test_config(10);
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

// The descent replaces Ligero's three clear-text combined rows (3 * num_cols field elements) with
// a logarithmic number of levels, so the opening must grow far slower than Ligero's in num_cols.
TYPED_TEST(SwitchFoldTest, OpeningIsSublinearInRowLength)
{
    std::vector<size_t> sizes;
    for (const size_t num_variables : { size_t(9), size_t(11) }) {
        const SwitchFoldConfig config = TestFixture::test_config(num_variables);
        typename TestFixture::CK ck(config);
        const auto instance = TestFixture::make_instance(ck, 2, 1);
        const auto proof = TestFixture::prove_instance(ck, instance);
        ASSERT_TRUE(TestFixture::verify_proof(config, instance.verifier_claims, instance.u, proof));
        sizes.push_back(proof.size());
    }
    // Four times the row length, but the descent only adds two levels plus wider Merkle paths.
    EXPECT_LT(sizes[1], 2 * sizes[0]);
}

class SwitchFoldHonkTest : public ::testing::Test {
  public:
    using Honk = SwitchFoldHonk<whir::Blake3sMerkleHasher>;

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

TEST_F(SwitchFoldHonkTest, ProveAndVerify)
{
    UltraCircuitBuilder sizing_builder = build_test_circuit();
    ProverInstance_<UltraFlavor> sizing_instance(sizing_builder);
    const size_t log_n = sizing_instance.log_dyadic_size();

    UltraCircuitBuilder proving_builder = build_test_circuit();
    const SwitchFoldConfig config = Honk::make_config(log_n, 64, 2);
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

// The whole point of the descent: same commitment as Ligero, materially smaller opening.
TEST_F(SwitchFoldHonkTest, ProofIsSmallerThanLigero)
{
    UltraCircuitBuilder sizing_builder = build_test_circuit();
    const size_t log_n = ProverInstance_<UltraFlavor>(sizing_builder).log_dyadic_size();

    UltraCircuitBuilder sf_builder = build_test_circuit();
    const SwitchFoldConfig sf_config = Honk::make_config(log_n, 64, 2);
    auto sf_pk = Honk::create_proving_key(sf_builder, sf_config);
    const HonkProof sf_proof = Honk::prove(sf_pk);
    ASSERT_TRUE(Honk::verify(sf_pk.vk, sf_config, sf_proof));

    using LigeroHonk = ligero::LigeroHonk<whir::Blake3sMerkleHasher>;
    UltraCircuitBuilder ligero_builder = build_test_circuit();
    const ligero::LigeroConfig ligero_config = LigeroHonk::make_config(log_n, 64, 2);
    auto ligero_pk = LigeroHonk::create_proving_key(ligero_builder, ligero_config);
    const HonkProof ligero_proof = LigeroHonk::prove(ligero_pk);
    ASSERT_TRUE(LigeroHonk::verify(ligero_pk.vk, ligero_config, ligero_proof));

    EXPECT_LT(sf_proof.size(), ligero_proof.size());
}

} // namespace bb::switchfold
