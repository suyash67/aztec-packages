#include "cq_fk.hpp"

#include "barretenberg/common/assert.hpp"
#include "barretenberg/common/thread.hpp"
#include "barretenberg/cq/cq_domain.hpp"

namespace bb::cq {

namespace {

g1::element infinity_element()
{
    g1::element result;
    result.self_set_infinity();
    return result;
}

std::vector<g1::affine_element> to_affine_batch(std::vector<g1::element>&& points)
{
    g1::element::batch_normalize(points.data(), points.size());
    std::vector<g1::affine_element> result(points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        result[i] = points[i].is_point_at_infinity() ? g1::affine_element::infinity()
                                                     : g1::affine_element(points[i].x, points[i].y);
    }
    return result;
}

} // namespace

std::vector<g1::element> toeplitz_srs_products_naive(std::span<const fr> coeffs,
                                                     std::span<const g1::affine_element> srs)
{
    const size_t n = coeffs.size();
    BB_ASSERT_GTE(srs.size(), n - 1, "toeplitz_srs_products_naive: SRS too small");
    std::vector<g1::element> h(n, infinity_element());
    parallel_for_range(n, [&](size_t start, size_t end) {
        for (size_t r = start; r < end; ++r) {
            g1::element acc = infinity_element();
            for (size_t d = 0; r + 1 + d < n; ++d) {
                acc += g1::element(srs[d]) * coeffs[r + 1 + d];
            }
            h[r] = acc;
        }
    });
    return h;
}

std::vector<g1::element> toeplitz_srs_products(std::span<const fr> coeffs, std::span<const g1::affine_element> srs)
{
    const size_t n = coeffs.size();
    BB_ASSERT(n >= 2 && (n & (n - 1)) == 0, "toeplitz_srs_products: size must be a power of 2");
    BB_ASSERT_GTE(srs.size(), n - 1, "toeplitz_srs_products: SRS too small");
    const size_t m = 2 * n;
    const SubgroupDomain domain(m);

    // Cyclic convolution of the SRS prefix with the reversed coefficient tail:
    //   h_r = sum_d coeffs[r+1+d] * srs[d] = (s_hat (*) c_rev)[2n - 1 - r],
    // where s_hat = (srs[0], ..., srs[n-2], inf, ..., inf) and c_rev places coeffs[n-1], ..., coeffs[1] in
    // positions n+1, ..., 2n-1.
    std::vector<g1::element> s_hat(m, infinity_element());
    for (size_t d = 0; d + 1 < n; ++d) {
        s_hat[d] = g1::element(srs[d]);
    }
    std::vector<fr> c_rev(m, fr::zero());
    for (size_t a = n + 1; a < m; ++a) {
        c_rev[a] = coeffs[m - a];
    }

    domain.fft<g1::element>(s_hat);
    domain.fft<fr>(c_rev);
    parallel_for_range(
        m,
        [&](size_t start, size_t end) {
            for (size_t k = start; k < end; ++k) {
                s_hat[k] = s_hat[k] * c_rev[k];
            }
        },
        /*no_multhreading_if_less_or_equal=*/64);
    domain.ifft<g1::element>(s_hat);

    std::vector<g1::element> h(n);
    for (size_t r = 0; r < n; ++r) {
        h[r] = s_hat[m - 1 - r];
    }
    return h;
}

std::vector<g1::affine_element> fk_all_kzg_opening_proofs(std::span<const fr> t_coeffs,
                                                          std::span<const g1::affine_element> srs)
{
    const size_t n = t_coeffs.size();
    std::vector<g1::element> h = toeplitz_srs_products(t_coeffs, srs);
    const SubgroupDomain domain(n);
    domain.fft<g1::element>(h);
    return to_affine_batch(std::move(h));
}

std::vector<g1::affine_element> lagrange_basis_commitments(std::span<const g1::affine_element> srs, size_t domain_size)
{
    BB_ASSERT_GTE(srs.size(), domain_size, "lagrange_basis_commitments: SRS too small");
    std::vector<g1::element> points(domain_size);
    for (size_t i = 0; i < domain_size; ++i) {
        points[i] = g1::element(srs[i]);
    }
    const SubgroupDomain domain(domain_size);
    domain.ifft<g1::element>(points);
    return to_affine_batch(std::move(points));
}

std::vector<g1::affine_element> lagrange_zero_quotient_commitments(std::span<const g1::affine_element> srs,
                                                                   size_t domain_size)
{
    BB_ASSERT_GTE(srs.size(), domain_size - 1, "lagrange_zero_quotient_commitments: SRS too small");
    std::vector<g1::element> points(domain_size, infinity_element());
    for (size_t i = 0; i + 1 < domain_size; ++i) {
        points[i] = g1::element(srs[i]);
    }
    const SubgroupDomain domain(domain_size);
    domain.ifft<g1::element>(points);
    // (L_i(X) - L_i(0))/X = omega^{-i} * (1/N) * sum_{k=0}^{N-2} omega^{-ik} X^k, so scale entry i by omega^{-i}.
    std::vector<fr> root_inverse_powers(domain_size);
    root_inverse_powers[0] = fr::one();
    for (size_t i = 1; i < domain_size; ++i) {
        root_inverse_powers[i] = root_inverse_powers[i - 1] * domain.root_inverse();
    }
    parallel_for_range(
        domain_size,
        [&](size_t start, size_t end) {
            for (size_t i = start; i < end; ++i) {
                points[i] = points[i] * root_inverse_powers[i];
            }
        },
        /*no_multhreading_if_less_or_equal=*/64);
    return to_affine_batch(std::move(points));
}

} // namespace bb::cq
