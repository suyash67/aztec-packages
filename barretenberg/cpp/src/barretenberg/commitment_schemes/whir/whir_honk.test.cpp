#include "barretenberg/commitment_schemes/whir/whir_honk.hpp"

#include "barretenberg/stdlib_circuit_builders/mock_circuits.hpp"

#include <gtest/gtest.h>

namespace bb::whir {

namespace {

/** @brief A circuit exercising arithmetic, public-input, lookup, and ROM (memory) gates. */
UltraCircuitBuilder build_test_circuit()
{
    UltraCircuitBuilder builder;
    MockCircuits::add_arithmetic_gates_with_public_inputs(builder, 16);
    MockCircuits::add_arithmetic_gates(builder, 16);
    MockCircuits::add_lookup_gates(builder, 2);
    // A small ROM table so the memory relation and the eta-dependent w_4 records are exercised.
    const size_t rom_id = builder.create_ROM_array(4);
    for (size_t i = 0; i < 4; ++i) {
        builder.set_ROM_element(rom_id, i, builder.add_variable(fr(3 * i + 1)));
    }
    const uint32_t read_index = builder.add_variable(fr(2));
    builder.read_ROM_array(rom_id, read_index);
    return builder;
}

} // namespace

template <typename Hasher> class WhirHonkTest : public ::testing::Test {
  public:
    using Honk = WhirHonk<Hasher>;

    struct Setup {
        WhirConfig config;
        typename Honk::VerificationKey vk;
        HonkProof proof;
    };

    static Setup prove_test_circuit()
    {
        UltraCircuitBuilder builder = build_test_circuit();
        // Sized from the circuit: ProverInstance dictates the dyadic size.
        ProverInstance_<UltraFlavor> sizing_instance(builder);
        const size_t log_n = sizing_instance.log_dyadic_size();

        UltraCircuitBuilder proving_builder = build_test_circuit();
        const WhirConfig config = Honk::make_config(log_n, /*security_bits=*/64);
        auto pk = Honk::create_proving_key(proving_builder, config);
        return { config, pk.vk, Honk::prove(pk) };
    }
};

using HasherTypes = ::testing::Types<Poseidon2MerkleHasher, Blake3sMerkleHasher>;
TYPED_TEST_SUITE(WhirHonkTest, HasherTypes);

TYPED_TEST(WhirHonkTest, ProveAndVerify)
{
    const auto setup = TestFixture::prove_test_circuit();
    EXPECT_TRUE(TestFixture::Honk::verify(setup.vk, setup.config, setup.proof));
}

TYPED_TEST(WhirHonkTest, TamperedProofRejected)
{
    const auto setup = TestFixture::prove_test_circuit();

    // Positions spread across the proof: public inputs / wire roots, sumcheck messages, WHIR
    // openings, the final polynomial.
    for (const size_t position :
         { size_t(2), setup.proof.size() / 4, setup.proof.size() / 2, setup.proof.size() - 2 }) {
        HonkProof tampered = setup.proof;
        tampered[position] += fr(1);
        EXPECT_FALSE(TestFixture::Honk::verify(setup.vk, setup.config, tampered)) << "position " << position;
    }
}

TYPED_TEST(WhirHonkTest, VirtualColumnsOmittedFromCommitment)
{
    const auto setup = TestFixture::prove_test_circuit();
    const uint32_t mask = setup.vk.virtual_mask;

    // The lagrange point indicators are always virtual (entities 8, 9), and the test circuit uses
    // no elliptic (23), nnf (25), or poseidon2 (26, 27) gates, so those selectors are zero columns.
    for (const size_t entity : { 8U, 9U, 23U, 25U, 26U, 27U }) {
        EXPECT_TRUE(((mask >> entity) & 1) == 1) << "entity " << entity;
    }
    EXPECT_LE(setup.vk.num_committed_precomputed(), 22U);
    EXPECT_TRUE(TestFixture::Honk::verify(setup.vk, setup.config, setup.proof));

    // The indicator rows are part of the statement: a wrong row must be rejected.
    auto bad_vk = setup.vk;
    bad_vk.lagrange_last_row ^= 1;
    EXPECT_FALSE(TestFixture::Honk::verify(bad_vk, setup.config, setup.proof));
}

TYPED_TEST(WhirHonkTest, WrongVkRejected)
{
    const auto setup = TestFixture::prove_test_circuit();

    auto bad_vk = setup.vk;
    auto root_fields = TypeParam::digest_to_fields(bad_vk.precomputed_commitment);
    root_fields[0] += fr(1);
    bad_vk.precomputed_commitment = TypeParam::digest_from_fields(root_fields);
    EXPECT_FALSE(TestFixture::Honk::verify(bad_vk, setup.config, setup.proof));
}

} // namespace bb::whir
