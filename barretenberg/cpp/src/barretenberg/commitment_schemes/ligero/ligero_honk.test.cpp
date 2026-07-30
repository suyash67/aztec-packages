#include "barretenberg/commitment_schemes/ligero/ligero_honk.hpp"

#include "barretenberg/stdlib_circuit_builders/mock_circuits.hpp"

#include <gtest/gtest.h>

namespace bb::ligero {

namespace {

/** @brief A circuit exercising arithmetic, public-input, lookup, and ROM (memory) gates. */
UltraCircuitBuilder build_test_circuit()
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

} // namespace

template <typename Hasher> class LigeroHonkTest : public ::testing::Test {
  public:
    using Honk = LigeroHonk<Hasher>;

    struct Setup {
        LigeroConfig config;
        typename Honk::VerificationKey vk;
        HonkProof proof;
    };

    static Setup prove_test_circuit()
    {
        UltraCircuitBuilder sizing_builder = build_test_circuit();
        ProverInstance_<UltraFlavor> sizing_instance(sizing_builder);
        const size_t log_n = sizing_instance.log_dyadic_size();

        UltraCircuitBuilder proving_builder = build_test_circuit();
        const LigeroConfig config = Honk::make_config(log_n, /*security_bits=*/64);
        auto pk = Honk::create_proving_key(proving_builder, config);
        return { config, pk.vk, Honk::prove(pk) };
    }
};

using HasherTypes = ::testing::Types<whir::Poseidon2MerkleHasher, whir::Blake3sMerkleHasher>;
TYPED_TEST_SUITE(LigeroHonkTest, HasherTypes);

TYPED_TEST(LigeroHonkTest, ProveAndVerify)
{
    const auto setup = TestFixture::prove_test_circuit();
    EXPECT_TRUE(TestFixture::Honk::verify(setup.vk, setup.config, setup.proof));
}

TYPED_TEST(LigeroHonkTest, TamperedProofRejected)
{
    const auto setup = TestFixture::prove_test_circuit();
    for (const size_t position :
         { size_t(2), setup.proof.size() / 4, setup.proof.size() / 2, setup.proof.size() - 2 }) {
        HonkProof tampered = setup.proof;
        tampered[position] += fr(1);
        EXPECT_FALSE(TestFixture::Honk::verify(setup.vk, setup.config, tampered)) << "position " << position;
    }
}

TYPED_TEST(LigeroHonkTest, WrongVkRejected)
{
    const auto setup = TestFixture::prove_test_circuit();
    auto bad_vk = setup.vk;
    auto root_fields = TypeParam::digest_to_fields(bad_vk.precomputed_root);
    root_fields[0] += fr(1);
    bad_vk.precomputed_root = TypeParam::digest_from_fields(root_fields);
    EXPECT_FALSE(TestFixture::Honk::verify(bad_vk, setup.config, setup.proof));
}

} // namespace bb::ligero
