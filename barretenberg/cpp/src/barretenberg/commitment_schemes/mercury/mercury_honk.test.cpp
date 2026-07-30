#include "barretenberg/commitment_schemes/mercury/mercury_honk.hpp"

#include "barretenberg/srs/global_crs.hpp"
#include "barretenberg/stdlib_circuit_builders/mock_circuits.hpp"

#include <gtest/gtest.h>

namespace bb::mercury {

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

class MercuryHonkTest : public ::testing::Test {
  public:
    using Honk = MercuryHonk;

    struct Setup {
        MercuryConfig config;
        Honk::VerificationKey vk;
        HonkProof proof;
    };

    static Setup prove_test_circuit()
    {
        UltraCircuitBuilder sizing_builder = build_test_circuit();
        ProverInstance_<UltraFlavor> sizing_instance(sizing_builder);
        const size_t log_n = sizing_instance.log_dyadic_size();

        UltraCircuitBuilder proving_builder = build_test_circuit();
        const MercuryConfig config = Honk::make_config(log_n);
        auto pk = Honk::create_proving_key(proving_builder, config);
        return { config, pk.vk, Honk::prove(pk) };
    }

  protected:
    static void SetUpTestSuite() { srs::init_file_crs_factory(srs::bb_crs_path()); }
};

TEST_F(MercuryHonkTest, ProveAndVerify)
{
    const auto setup = prove_test_circuit();
    EXPECT_TRUE(Honk::verify(setup.vk, setup.config, setup.proof));
}

TEST_F(MercuryHonkTest, TamperedProofRejected)
{
    const auto setup = prove_test_circuit();
    for (const size_t position :
         { size_t(2), setup.proof.size() / 4, setup.proof.size() / 2, setup.proof.size() - 2 }) {
        HonkProof tampered = setup.proof;
        tampered[position] += fr(1);
        bool accepted = false;
        try {
            accepted = Honk::verify(setup.vk, setup.config, tampered);
        } catch (const std::exception&) {
            accepted = false; // malformed proof data (e.g. off-curve point) is a rejection
        }
        EXPECT_FALSE(accepted) << "position " << position;
    }
}

TEST_F(MercuryHonkTest, WrongVkRejected)
{
    const auto setup = prove_test_circuit();
    auto bad_vk = setup.vk;
    std::swap(bad_vk.precomputed_commitment[3], bad_vk.precomputed_commitment[4]);
    EXPECT_FALSE(Honk::verify(bad_vk, setup.config, setup.proof));
}

} // namespace bb::mercury
