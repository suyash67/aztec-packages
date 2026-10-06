#pragma once

#include "barretenberg/common/thread.hpp"
#include "barretenberg/numeric/bitop/get_msb.hpp"
#include "barretenberg/numeric/uint256/uint256.hpp"
#include "barretenberg/polynomials/polynomial.hpp"

#include <algorithm>
#include <span>
#include <vector>

/**
 * @file pippenger_generic.hpp
 * @brief Pippenger's bucket method for curves without a GLV endomorphism split (e.g. Pallas / Vesta).
 *
 * @details Full-width scalars are cut into signed c-bit digits (d in [-2^{c-1}, 2^{c-1}]); every window accumulates
 * its points into 2^{c-1} buckets with mixed additions, reduces the buckets with a running sum, and the windows are
 * combined by Horner's rule. Windows are processed in parallel. Every operation is a complete projective addition, so
 * no input restrictions apply (repeated points and the point at infinity are fine).
 */
namespace bb::scalar_multiplication::generic {

inline size_t window_bits(size_t num_points)
{
    if (num_points < 32) {
        return 3;
    }
    // ~ log2(n) - 2, the usual optimum for the bucket method with a cheap reduction.
    return std::clamp<size_t>(numeric::get_msb(num_points) - 2, 4, 16);
}

template <typename Curve>
typename Curve::Element pippenger(PolynomialSpan<const typename Curve::ScalarField> scalars,
                                  std::span<const typename Curve::AffineElement> points)
{
    using Fr = typename Curve::ScalarField;
    using Element = typename Curve::Element;
    using AffineElement = typename Curve::AffineElement;

    const size_t n = scalars.size();
    BB_ASSERT_LTE(scalars.start_index + n, points.size());
    if (n == 0) {
        return Element::infinity();
    }
    const auto point = [&](size_t i) -> const AffineElement& { return points[scalars.start_index + i]; };

    // Integer representatives of the scalars.
    std::vector<uint256_t> k(n);
    parallel_for_range(n, [&](size_t start, size_t end) {
        for (size_t i = start; i < end; ++i) {
            k[i] = uint256_t(scalars.span[i]);
        }
    });

    const size_t c = window_bits(n);
    const size_t num_bits = Fr::modulus.get_msb() + 1;
    // One extra window absorbs the final carry of the signed recoding.
    const size_t num_windows = ((num_bits + c - 1) / c) + 1;
    const size_t num_buckets = size_t{ 1 } << (c - 1);

    // Signed digits, window-major: digit(w, i) in [-2^{c-1}, 2^{c-1}].
    std::vector<int32_t> digits(num_windows * n);
    parallel_for_range(n, [&](size_t start, size_t end) {
        const uint64_t mask = (uint64_t{ 1 } << c) - 1;
        const int64_t half = int64_t{ 1 } << (c - 1);
        for (size_t i = start; i < end; ++i) {
            int64_t carry = 0;
            for (size_t w = 0; w < num_windows; ++w) {
                const size_t lo = w * c;
                const uint64_t bits = (lo < 256) ? (k[i] >> lo).data[0] & mask : 0;
                int64_t d = static_cast<int64_t>(bits) + carry;
                carry = 0;
                if (d > half) {
                    d -= int64_t{ 1 } << c;
                    carry = 1;
                }
                digits[(w * n) + i] = static_cast<int32_t>(d);
            }
        }
    });

    std::vector<Element> window_sums(num_windows, Element::infinity());
    parallel_for(num_windows, [&](size_t w) {
        std::vector<Element> buckets(num_buckets, Element::infinity());
        const int32_t* d = &digits[w * n];
        for (size_t i = 0; i < n; ++i) {
            const int32_t digit = d[i];
            if (digit > 0) {
                buckets[static_cast<size_t>(digit - 1)] += point(i);
            } else if (digit < 0) {
                buckets[static_cast<size_t>(-digit - 1)] -= point(i);
            }
        }
        Element running = Element::infinity();
        Element sum = Element::infinity();
        for (size_t b = num_buckets; b-- > 0;) {
            running += buckets[b];
            sum += running;
        }
        window_sums[w] = sum;
    });

    Element result = window_sums[num_windows - 1];
    for (size_t w = num_windows - 1; w-- > 0;) {
        for (size_t j = 0; j < c; ++j) {
            result.self_dbl();
        }
        result += window_sums[w];
    }
    return result;
}

} // namespace bb::scalar_multiplication::generic
