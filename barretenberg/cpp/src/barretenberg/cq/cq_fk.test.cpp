#include "cq_fk.hpp"

#include "barretenberg/cq/cq_domain.hpp"
#include "barretenberg/cq/cq_trusted_setup.hpp"
#include "barretenberg/ecc/scalar_multiplication/scalar_multiplication.hpp"
#include "barretenberg/polynomials/polynomial_arithmetic.hpp"

#include <gtest/gtest.h>

namespace bb::cq {

class CqFkTest : public ::testing::Test {
  public:
    static constexpr size_t N = 32;
    static TestSrs srs;

    static void SetUpTestSuite() { srs = TestSrs::create(N, N + 1); }

    static g1::affine_element commit(std::span<const fr> coeffs)
    {
        return g1::affine_element(scalar_multiplication::pippenger<curve::BN254>(
            PolynomialSpan<const fr>{ 0, coeffs }, srs.g1_powers, /*handle_edge_cases=*/true));
    }
};

TestSrs CqFkTest::srs;

TEST_F(CqFkTest, ToeplitzFftMatchesNaive)
{
    std::vector<fr> coeffs(N);
    for (auto& c : coeffs) {
        c = fr::random_element();
    }
    const std::vector<g1::element> expected = toeplitz_srs_products_naive(coeffs, srs.g1_powers);
    const std::vector<g1::element> actual = toeplitz_srs_products(coeffs, srs.g1_powers);
    ASSERT_EQ(actual.size(), expected.size());
    for (size_t r = 0; r < N; ++r) {
        EXPECT_EQ(actual[r], expected[r]) << "row " << r;
    }
}

TEST_F(CqFkTest, FkMatchesNaiveKzgQuotients)
{
    const SubgroupDomain domain(N);
    std::vector<fr> t_coeffs(N);
    for (auto& c : t_coeffs) {
        c = fr::random_element();
    }
    const std::vector<g1::affine_element> proofs = fk_all_kzg_opening_proofs(t_coeffs, srs.g1_powers);
    ASSERT_EQ(proofs.size(), N);
    const std::vector<fr> roots = domain.root_powers();
    for (size_t i = 0; i < N; ++i) {
        // U_i(X) = (T(X) - T(omega^i)) / (X - omega^i) computed directly.
        std::vector<fr> quotient = t_coeffs;
        quotient[0] -= polynomial_arithmetic::evaluate<fr>(t_coeffs, roots[i]);
        polynomial_arithmetic::factor_roots<fr>(quotient, roots[i]);
        EXPECT_EQ(proofs[i], commit(quotient)) << "opening " << i;
    }
}

TEST_F(CqFkTest, LagrangeCommitmentsMatchNaive)
{
    const SubgroupDomain domain(N);
    const std::vector<g1::affine_element> commitments = lagrange_basis_commitments(srs.g1_powers, N);
    ASSERT_EQ(commitments.size(), N);
    for (size_t i = 0; i < N; ++i) {
        std::vector<fr> lagrange_coeffs(N, fr::zero());
        lagrange_coeffs[i] = fr::one();
        domain.ifft<fr>(lagrange_coeffs);
        EXPECT_EQ(commitments[i], commit(lagrange_coeffs)) << "Lagrange " << i;
    }
}

TEST_F(CqFkTest, LagrangeZeroQuotientsMatchNaive)
{
    const SubgroupDomain domain(N);
    const std::vector<g1::affine_element> commitments = lagrange_zero_quotient_commitments(srs.g1_powers, N);
    ASSERT_EQ(commitments.size(), N);
    for (size_t i = 0; i < N; ++i) {
        std::vector<fr> lagrange_coeffs(N, fr::zero());
        lagrange_coeffs[i] = fr::one();
        domain.ifft<fr>(lagrange_coeffs);
        // (L_i(X) - L_i(0)) / X drops the constant term and shifts down.
        const std::vector<fr> shifted(lagrange_coeffs.begin() + 1, lagrange_coeffs.end());
        EXPECT_EQ(commitments[i], commit(shifted)) << "Lagrange zero quotient " << i;
    }
}

} // namespace bb::cq
