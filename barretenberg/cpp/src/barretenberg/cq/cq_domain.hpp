#pragma once

#include "barretenberg/common/assert.hpp"
#include "barretenberg/common/thread.hpp"
#include "barretenberg/ecc/curves/bn254/fr.hpp"
#include "barretenberg/numeric/bitop/get_msb.hpp"

#include <span>
#include <vector>

namespace bb::cq {

/**
 * @brief Multiplicative subgroup H = <omega> of order 2^k in the BN254 scalar field, with radix-2 FFT/IFFT.
 *
 * @details The transforms are templated over the coefficient type so the same code serves two purposes:
 *  - T = fr: ordinary polynomial FFTs (interpolation of tables and witness columns),
 *  - T = g1::element: FFTs over vectors of group elements, required by the Feist-Khovratovich algorithm that
 *    computes all N KZG opening proofs of the table polynomial in O(N log N) group operations.
 * The only operations required of T are addition, subtraction and right-multiplication by fr, which both fr and
 * Jacobian group elements provide.
 *
 * Convention: fft() maps coefficients (a_0, ..., a_{n-1}) to evaluations (A(omega^0), ..., A(omega^{n-1})) in
 * natural order; ifft() is its inverse.
 */
class SubgroupDomain {
  public:
    explicit SubgroupDomain(size_t size)
        : size_(size)
        , log_size_(numeric::get_msb(size))
    {
        BB_ASSERT(size >= 2 && (size & (size - 1)) == 0, "SubgroupDomain size must be a power of 2 and >= 2");
        root_ = fr::get_root_of_unity(log_size_);
        root_inverse_ = root_.invert();
        size_inverse_ = fr(size_).invert();
        forward_twiddles_ = compute_twiddles(root_);
        inverse_twiddles_ = compute_twiddles(root_inverse_);
    }

    size_t size() const { return size_; }
    size_t log_size() const { return log_size_; }
    fr root() const { return root_; }
    fr root_inverse() const { return root_inverse_; }
    fr size_inverse() const { return size_inverse_; }

    /**
     * @brief omega^i for i in [0, size); computed on demand (O(size) scalar muls).
     */
    std::vector<fr> root_powers() const
    {
        std::vector<fr> powers(size_);
        powers[0] = fr::one();
        for (size_t i = 1; i < size_; ++i) {
            powers[i] = powers[i - 1] * root_;
        }
        return powers;
    }

    template <typename T> void fft(std::span<T> data) const { transform(data, forward_twiddles_); }

    template <typename T> void ifft(std::span<T> data) const
    {
        transform(data, inverse_twiddles_);
        parallel_for_range(
            size_,
            [&](size_t start, size_t end) {
                for (size_t i = start; i < end; ++i) {
                    data[i] = data[i] * size_inverse_;
                }
            },
            /*no_multhreading_if_less_or_equal=*/512);
    }

    /**
     * @brief Evaluations of the vanishing polynomial Z_H(X) = X^n - 1 at a point.
     */
    fr evaluate_vanishing(const fr& x) const { return x.pow(size_) - fr::one(); }

  private:
    // twiddles[j] = root^j for j in [0, size/2); precomputed so each FFT stage can index twiddles directly and
    // butterflies parallelize without a sequential twiddle recurrence.
    std::vector<fr> compute_twiddles(const fr& root) const
    {
        std::vector<fr> twiddles(size_ / 2);
        twiddles[0] = fr::one();
        for (size_t i = 1; i < size_ / 2; ++i) {
            twiddles[i] = twiddles[i - 1] * root;
        }
        return twiddles;
    }

    template <typename T> void transform(std::span<T> data, const std::vector<fr>& twiddles) const
    {
        BB_ASSERT_EQ(data.size(), size_, "SubgroupDomain: data size must match domain size");
        bit_reverse_permute(data);
        // Iterative Cooley-Tukey, decimation in time. Stage with half-block size m combines pairs at distance m;
        // the butterfly with in-block index j uses twiddle root^(j * size / (2m)).
        for (size_t m = 1; m < size_; m <<= 1) {
            const size_t twiddle_stride = size_ / (2 * m);
            const size_t num_butterflies = size_ / 2;
            parallel_for_range(
                num_butterflies,
                [&](size_t start, size_t end) {
                    for (size_t b = start; b < end; ++b) {
                        const size_t block = b / m;
                        const size_t j = b % m;
                        const size_t lo = block * 2 * m + j;
                        const size_t hi = lo + m;
                        T odd_term = (j == 0) ? data[hi] : data[hi] * twiddles[j * twiddle_stride];
                        T even_term = data[lo];
                        data[lo] = even_term + odd_term;
                        data[hi] = even_term - odd_term;
                    }
                },
                /*no_multhreading_if_less_or_equal=*/256);
        }
    }

    template <typename T> void bit_reverse_permute(std::span<T> data) const
    {
        for (size_t i = 0; i < size_; ++i) {
            const size_t j = reverse_bits(i);
            if (i < j) {
                std::swap(data[i], data[j]);
            }
        }
    }

    size_t reverse_bits(size_t x) const
    {
        size_t result = 0;
        for (size_t bit = 0; bit < log_size_; ++bit) {
            result = (result << 1) | ((x >> bit) & 1);
        }
        return result;
    }

    size_t size_;
    size_t log_size_;
    fr root_;
    fr root_inverse_;
    fr size_inverse_;
    std::vector<fr> forward_twiddles_;
    std::vector<fr> inverse_twiddles_;
};

} // namespace bb::cq
