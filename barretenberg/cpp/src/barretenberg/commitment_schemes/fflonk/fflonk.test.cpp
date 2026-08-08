#include "barretenberg/commitment_schemes/fflonk/fflonk_honk.hpp"

#include "barretenberg/commitment_schemes/pcs_test_utils.hpp"
#include "barretenberg/stdlib_circuit_builders/mock_circuits.hpp"
#include "barretenberg/srs/global_crs.hpp"
#include "barretenberg/transcript/transcript.hpp"

#include <gtest/gtest.h>

namespace bb::fflonk {

namespace {

/** @brief A circuit exercising arithmetic, public-input, lookup and ROM gates. */
UltraCircuitBuilder build_test_circuit()
{
    UltraCircuitBuilder builder;
    MockCircuits::add_arithmetic_gates_with_public_inputs(builder, 8);
    const fr a_value(1);
    const auto a_idx = builder.add_variable(a_value);
    const auto accumulators =
        plookup::get_lookup_accumulators(plookup::MultiTableId::HONK_DUMMY_MULTI, a_value, fr(0), true);
    builder.create_gates_from_plookup_accumulators(
        plookup::MultiTableId::HONK_DUMMY_MULTI, accumulators, a_idx, std::nullopt);
    const size_t rom_id = builder.create_ROM_array(4);
    for (size_t i = 0; i < 4; ++i) {
        builder.set_ROM_element(rom_id, i, builder.add_variable(fr(3 * i + 1)));
    }
    builder.read_ROM_array(rom_id, builder.add_variable(fr(2)));
    return builder;
}

} // namespace

TEST(FflonkTest, PackedResidueIsTheEvaluationVector)
{
    // g mod (X^t - z) has the columns' evaluations at z as its coefficients: the identity the whole
    // construction rests on.
    constexpr size_t t = 4;
    constexpr size_t n = 8;
    std::vector<std::vector<fr>> columns(t, std::vector<fr>(n));
    for (auto& column : columns) {
        for (fr& value : column) {
            value = fr::random_element();
        }
    }
    std::vector<fr> packed(n * t, fr::zero());
    for (size_t i = 0; i < t; ++i) {
        for (size_t j = 0; j < n; ++j) {
            packed[j * t + i] = columns[i][j];
        }
    }
    const fr z = fr::random_element();
    std::vector<fr> quotient;
    std::vector<fr> residue;
    detail::divide_by_power_minus(packed, t, z, quotient, residue);
    for (size_t i = 0; i < t; ++i) {
        EXPECT_EQ(residue[i], polynomial_arithmetic::evaluate(columns[i].data(), z, n));
    }
}

TEST(FflonkTest, InterpolantMatchesTheResidue)
{
    // The closed form for R(y) agrees with the residue modulo (X^t - z)(X^t - 1/z).
    constexpr size_t t = 4;
    constexpr size_t n = 8;
    std::vector<fr> packed(n * t);
    for (fr& value : packed) {
        value = fr::random_element();
    }
    const fr z = fr::random_element();
    const fr z_inv = z.invert();
    const fr y = fr::random_element();

    std::vector<fr> quotient;
    std::vector<fr> p_evals;
    std::vector<fr> q_evals;
    detail::divide_by_power_minus(packed, t, z, quotient, p_evals);
    detail::divide_by_power_minus(packed, t, z_inv, quotient, q_evals);

    std::vector<fr> quotient_one;
    std::vector<fr> residue_one;
    detail::divide_by_power_minus(packed, t, z, quotient_one, residue_one);
    std::vector<fr> quotient_two;
    std::vector<fr> residue_two;
    detail::divide_by_power_minus(quotient_one, t, z_inv, quotient_two, residue_two);
    // R(X) = residue_one(X) + residue_two(X) * (X^t - z)
    const fr y_pow = y.pow(uint256_t(uint64_t(t)));
    const fr expected =
        detail::evaluate(residue_one, y) + detail::evaluate(residue_two, y) * (y_pow - z);
    EXPECT_EQ(detail::interpolant_at(p_evals, q_evals, z, z_inv, y, t), expected);
}

class FflonkHonkTest : public ::testing::Test {
  public:
    static void SetUpTestSuite() { srs::init_file_crs_factory(srs::bb_crs_path()); }
};

TEST_F(FflonkHonkTest, ProveAndVerify)
{
    UltraCircuitBuilder builder = build_test_circuit();
    UltraCircuitBuilder sizing = build_test_circuit();
    const size_t log_n = ProverInstance_<UltraFlavor>(sizing).log_dyadic_size();
    const auto config = FflonkHonk::make_config(log_n);
    auto pk = FflonkHonk::create_proving_key(builder, config);
    const auto vk = pk.vk;
    const HonkProof proof = FflonkHonk::prove(pk);
    EXPECT_TRUE(FflonkHonk::verify(vk, config, proof));
}

TEST_F(FflonkHonkTest, TamperedProofRejected)
{
    UltraCircuitBuilder builder = build_test_circuit();
    UltraCircuitBuilder sizing = build_test_circuit();
    const size_t log_n = ProverInstance_<UltraFlavor>(sizing).log_dyadic_size();
    const auto config = FflonkHonk::make_config(log_n);
    auto pk = FflonkHonk::create_proving_key(builder, config);
    const auto vk = pk.vk;
    HonkProof proof = FflonkHonk::prove(pk);
    // A tampered proof must not be accepted. Positions that fall inside a group element make
    // deserialization throw rather than return false; either outcome is a rejection.
    for (const size_t position : { size_t(1), proof.size() / 3, proof.size() / 2, proof.size() - 1 }) {
        HonkProof tampered = proof;
        tampered[position] += fr::one();
        bool accepted = false;
        try {
            accepted = FflonkHonk::verify(vk, config, tampered);
        } catch (...) {
        }
        EXPECT_FALSE(accepted) << "position " << position;
    }
}

} // namespace bb::fflonk
