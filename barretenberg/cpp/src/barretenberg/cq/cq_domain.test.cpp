#include "cq_domain.hpp"

#include "barretenberg/ecc/curves/bn254/g1.hpp"

#include <gtest/gtest.h>

namespace bb::cq {

namespace {

// Reference DFT: evaluations[j] = sum_k coeffs[k] * root^{jk}.
std::vector<fr> naive_dft(std::span<const fr> coeffs, const fr& root)
{
    const size_t n = coeffs.size();
    std::vector<fr> result(n, fr::zero());
    for (size_t j = 0; j < n; ++j) {
        fr point = root.pow(j);
        fr power = fr::one();
        for (size_t k = 0; k < n; ++k) {
            result[j] += coeffs[k] * power;
            power *= point;
        }
    }
    return result;
}

} // namespace

TEST(CqDomain, FftMatchesNaiveDft)
{
    for (const size_t n : { 2UL, 4UL, 8UL, 32UL, 128UL }) {
        SubgroupDomain domain(n);
        std::vector<fr> coeffs(n);
        for (auto& c : coeffs) {
            c = fr::random_element();
        }
        const std::vector<fr> expected = naive_dft(coeffs, domain.root());
        std::vector<fr> actual = coeffs;
        domain.fft<fr>(actual);
        EXPECT_EQ(actual, expected) << "size " << n;
    }
}

TEST(CqDomain, IfftInvertsFft)
{
    const size_t n = 64;
    SubgroupDomain domain(n);
    std::vector<fr> coeffs(n);
    for (auto& c : coeffs) {
        c = fr::random_element();
    }
    std::vector<fr> transformed = coeffs;
    domain.fft<fr>(transformed);
    domain.ifft<fr>(transformed);
    EXPECT_EQ(transformed, coeffs);
}

TEST(CqDomain, IfftInterpolates)
{
    // ifft of evaluations must produce coefficients whose evaluation at each omega^j gives back the value.
    const size_t n = 16;
    SubgroupDomain domain(n);
    std::vector<fr> evaluations(n);
    for (auto& e : evaluations) {
        e = fr::random_element();
    }
    std::vector<fr> coeffs = evaluations;
    domain.ifft<fr>(coeffs);
    for (size_t j = 0; j < n; ++j) {
        const fr point = domain.root().pow(j);
        fr result = fr::zero();
        fr power = fr::one();
        for (size_t k = 0; k < n; ++k) {
            result += coeffs[k] * power;
            power *= point;
        }
        EXPECT_EQ(result, evaluations[j]) << "evaluation " << j;
    }
}

TEST(CqDomain, GroupFftMatchesFieldFft)
{
    // FFT commutes with the homomorphism x -> x*G, so the group FFT of (c_i * G) equals the field FFT mapped to G.
    const size_t n = 32;
    SubgroupDomain domain(n);
    std::vector<fr> coeffs(n);
    std::vector<g1::element> points(n);
    for (size_t i = 0; i < n; ++i) {
        coeffs[i] = fr::random_element();
        points[i] = g1::one * coeffs[i];
    }
    domain.fft<fr>(coeffs);
    domain.fft<g1::element>(points);
    for (size_t i = 0; i < n; ++i) {
        EXPECT_EQ(points[i], g1::one * coeffs[i]) << "index " << i;
    }
}

TEST(CqDomain, GroupIfftMatchesFieldIfft)
{
    const size_t n = 16;
    SubgroupDomain domain(n);
    std::vector<fr> coeffs(n);
    std::vector<g1::element> points(n);
    for (size_t i = 0; i < n; ++i) {
        coeffs[i] = fr::random_element();
        points[i] = g1::one * coeffs[i];
    }
    domain.ifft<fr>(coeffs);
    domain.ifft<g1::element>(points);
    for (size_t i = 0; i < n; ++i) {
        EXPECT_EQ(points[i], g1::one * coeffs[i]) << "index " << i;
    }
}

TEST(CqDomain, GroupFftHandlesInfinity)
{
    // Vectors with points at infinity (used for zero padding in the FK convolution) must transform consistently.
    const size_t n = 8;
    SubgroupDomain domain(n);
    std::vector<fr> coeffs(n, fr::zero());
    std::vector<g1::element> points(n);
    for (auto& p : points) {
        p.self_set_infinity();
    }
    for (size_t i = 0; i < n / 2; ++i) {
        coeffs[i] = fr::random_element();
        points[i] = g1::one * coeffs[i];
    }
    domain.fft<fr>(coeffs);
    domain.fft<g1::element>(points);
    for (size_t i = 0; i < n; ++i) {
        EXPECT_EQ(points[i], g1::one * coeffs[i]) << "index " << i;
    }
}

} // namespace bb::cq
