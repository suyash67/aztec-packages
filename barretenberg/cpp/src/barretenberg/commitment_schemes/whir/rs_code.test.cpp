#include "barretenberg/commitment_schemes/whir/rs_code.hpp"

#include "barretenberg/polynomials/polynomial.hpp"

#include <gtest/gtest.h>

namespace bb::whir {

namespace {
std::vector<fr> random_array(size_t size)
{
    std::vector<fr> array(size);
    for (fr& value : array) {
        value = fr::random_element();
    }
    return array;
}
} // namespace

// codeword[i] = A(ω^i): the FFT output is in natural order.
TEST(WhirRSCode, EncodeIsNaturalOrderEvaluation)
{
    constexpr size_t n = 8;
    constexpr size_t domain_size = 32;
    const std::vector<fr> coeffs = random_array(n);
    RSDomains domains;
    const auto& domain = domains.get(domain_size);
    const std::vector<fr> codeword = rs_encode(coeffs, domain);
    ASSERT_EQ(codeword.size(), domain_size);

    fr point = fr::one();
    for (size_t i = 0; i < domain_size; ++i) {
        EXPECT_EQ(codeword[i], polynomial_arithmetic::evaluate(coeffs.data(), point, n)) << "index " << i;
        point *= domain.root;
    }
}

// Folding all variables at u equals the MLE evaluation of the array at u.
TEST(WhirRSCode, FoldArrayIsMlePartialEvaluation)
{
    constexpr size_t m = 4;
    std::vector<fr> array = random_array(size_t(1) << m);
    const Polynomial<fr> reference{ std::span<const fr>(array) };

    std::vector<fr> u;
    for (size_t j = 0; j < m; ++j) {
        u.push_back(fr::random_element());
    }
    for (const fr& u_j : u) {
        fold_array_in_place(array, u_j);
    }
    ASSERT_EQ(array.size(), 1U);
    EXPECT_EQ(array[0], reference.evaluate_mle(u));
}

// fold_coset over codeword values equals evaluating the folded coefficient array at x_base^{2^k}:
// README.md eq. (3.1) iterated k times.
TEST(WhirRSCode, FoldCosetMatchesFoldedEncode)
{
    constexpr size_t m = 6;
    constexpr size_t k = 2;
    constexpr size_t domain_size = size_t(1) << (m + 2);
    const std::vector<fr> coeffs = random_array(size_t(1) << m);
    RSDomains domains;
    const auto& domain = domains.get(domain_size);
    const std::vector<fr> codeword = rs_encode(coeffs, domain);

    const std::vector<fr> alphas = { fr::random_element(), fr::random_element() };
    std::vector<fr> folded = coeffs;
    for (const fr& alpha : alphas) {
        fold_array_in_place(folded, alpha);
    }

    const size_t num_cosets = domain_size >> k;
    const fr eta = domain.root.pow(uint256_t(num_cosets));
    const fr eta_inv = eta.invert();
    for (const size_t j : { size_t(0), size_t(1), size_t(7), num_cosets - 1 }) {
        std::vector<fr> coset_values(size_t(1) << k);
        for (size_t t = 0; t < coset_values.size(); ++t) {
            coset_values[t] = codeword[j + t * num_cosets];
        }
        const fr x_base = domain.root.pow(uint256_t(j));
        const fr folded_value = fold_coset(coset_values, x_base.invert(), eta_inv, alphas);
        const fr expected =
            polynomial_arithmetic::evaluate(folded.data(), x_base.pow(uint256_t(1) << k), folded.size());
        EXPECT_EQ(folded_value, expected) << "coset " << j;
    }
}

// For a to-be-shifted array (zero constant term), the shifted array's codeword is A(x)/x pointwise.
TEST(WhirRSCode, ShiftedCodewordIsDivisionByX)
{
    constexpr size_t n = 16;
    constexpr size_t domain_size = 64;
    std::vector<fr> coeffs = random_array(n);
    coeffs[0] = fr::zero();
    std::vector<fr> shifted(coeffs.begin() + 1, coeffs.end());

    RSDomains domains;
    const auto& domain = domains.get(domain_size);
    const std::vector<fr> codeword = rs_encode(coeffs, domain);
    const std::vector<fr> shifted_codeword = rs_encode(shifted, domain);

    fr point_inv = fr::one();
    for (size_t i = 0; i < domain_size; ++i) {
        EXPECT_EQ(shifted_codeword[i], codeword[i] * point_inv) << "index " << i;
        point_inv *= domain.root_inverse;
    }
}

} // namespace bb::whir
