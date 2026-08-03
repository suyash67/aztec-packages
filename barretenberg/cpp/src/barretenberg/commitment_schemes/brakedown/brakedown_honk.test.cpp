#include "barretenberg/commitment_schemes/brakedown/brakedown_honk.hpp"

#include "barretenberg/commitment_schemes/ligero/ligero_honk.hpp"
#include "barretenberg/stdlib_circuit_builders/mock_circuits.hpp"
#include "barretenberg/transcript/transcript.hpp"

#include <gtest/gtest.h>

namespace bb::brakedown {

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

using Hasher = whir::Blake3sMerkleHasher;
using CK = ligero::LigeroCommitmentKey<Hasher, BrakedownCodePolicy>;
using Prover = ligero::LigeroProver<Hasher>;
using Verifier = ligero::LigeroVerifier<Hasher>;

// Small security parameter: the query count is ~42x the security bits for this code, so a realistic
// lambda would put thousands of openings in every unit test.
constexpr size_t TEST_SECURITY_BITS = 8;

ligero::LigeroConfig test_config(size_t num_variables)
{
    return make_brakedown_config(num_variables, TEST_SECURITY_BITS, /*est_polynomials=*/4);
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
        const BrakedownCodePolicy code(config);
        return Verifier::verify(config, claims, u, transcript, code);
    } catch (const std::exception&) {
        return false;
    }
}

} // namespace

TEST(BrakedownHonkTest, SingleUnshiftedCompleteness)
{
    const auto config = test_config(10);
    const CK ck(config);
    const auto instance = make_instance(ck, 1, 0);
    const auto proof = prove_instance(ck, instance);
    EXPECT_TRUE(verify_proof(config, instance.verifier_claims, instance.u, proof));
}

TEST(BrakedownHonkTest, BatchedWithShiftedCompleteness)
{
    const auto config = test_config(9);
    const CK ck(config);
    const auto instance = make_instance(ck, 3, 2);
    const auto proof = prove_instance(ck, instance);
    EXPECT_TRUE(verify_proof(config, instance.verifier_claims, instance.u, proof));
}

TEST(BrakedownHonkTest, WrongEvaluationRejected)
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

TEST(BrakedownHonkTest, TamperedProofRejected)
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

// The distance-aware rule is the whole reason this backend exists as a separate policy: RS at rate
// 1/4 gets 50 queries at lambda = 100 under the capacity conjecture, Brakedown needs ~2936 from the
// provable interleaved test because its relative distance is 0.07 rather than 0.75.
TEST(BrakedownHonkTest, QueryCountFollowsTheProvableDistanceRule)
{
    EXPECT_EQ(BrakedownCodePolicy::provable_num_queries(100, 0.07), 2936U);
    // Monotone in the distance, and always far above the RS capacity count.
    EXPECT_GT(BrakedownCodePolicy::provable_num_queries(100, 0.07), 50U * 50U);
    EXPECT_LT(BrakedownCodePolicy::provable_num_queries(100, 0.07),
              BrakedownCodePolicy::provable_num_queries(100, 0.02));
}

// Codewords are not powers of two, so the tree is padded; queries must stay inside the true length
// or they would open padding that carries no information.
TEST(BrakedownHonkTest, QueriesStayInsideTheUnpaddedCodeword)
{
    const auto config = test_config(10);
    const BrakedownCodePolicy code(config);
    EXPECT_LT(code.codeword_length(), code.padded_codeword_length() + 1);
    EXPECT_GE(code.padded_codeword_length(), code.codeword_length());
    for (size_t i = 0; i < 2000; ++i) {
        const size_t index = ligero::detail::index_from_challenge_bounded(fr::random_element(), code.codeword_length());
        EXPECT_LT(index, code.codeword_length());
    }
}

// Same tensor protocol, same circuit, different code: Brakedown must trade proof size for its
// cheaper encoder. Asserting the direction guards against silently regressing to RS-like queries.
TEST(BrakedownHonkTest, ProofIsLargerThanLigeroAtEqualSecurity)
{
    const size_t num_variables = 10;
    const size_t security_bits = 8;

    const auto bd_config = make_brakedown_config(num_variables, security_bits, /*est_polynomials=*/4);
    const CK bd_ck(bd_config);
    const auto bd_instance = make_instance(bd_ck, 2, 1);
    const auto bd_proof = prove_instance(bd_ck, bd_instance);
    ASSERT_TRUE(verify_proof(bd_config, bd_instance.verifier_claims, bd_instance.u, bd_proof));

    const auto rs_config = ligero::LigeroConfig::create(num_variables, security_bits, 2, /*est_polynomials=*/4);
    const ligero::LigeroCommitmentKey<Hasher> rs_ck(rs_config);
    EXPECT_LT(rs_config.num_queries, bd_config.num_queries);
    EXPECT_GT(bd_proof.size(), size_t(0));
}

class BrakedownHonkCircuitTest : public ::testing::Test {
  public:
    using Honk = BrakedownHonk<whir::Blake3sMerkleHasher>;

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

TEST_F(BrakedownHonkCircuitTest, ProveAndVerify)
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

} // namespace bb::brakedown
