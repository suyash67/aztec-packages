#pragma once

#include "barretenberg/common/assert.hpp"
#include "barretenberg/common/thread.hpp"
#include "barretenberg/ecc/curves/bn254/fr.hpp"
#include "barretenberg/polynomials/evaluation_domain.hpp"
#include "barretenberg/polynomials/polynomial_arithmetic.hpp"

#include <map>
#include <memory>
#include <span>
#include <vector>

namespace bb::whir {

/**
 * @brief Lazily-built cache of FFT evaluation domains, one per codeword size used by the schedule.
 */
class RSDomains {
  public:
    const EvaluationDomain<fr>& get(size_t size)
    {
        auto it = domains_.find(size);
        if (it == domains_.end()) {
            auto domain = std::make_unique<EvaluationDomain<fr>>(size);
            domain->compute_lookup_table();
            it = domains_.emplace(size, std::move(domain)).first;
        }
        return *it->second;
    }

  private:
    std::map<size_t, std::unique_ptr<EvaluationDomain<fr>>> domains_;
};

/**
 * @brief Reed-Solomon encode: evaluate the univariate with coefficients `coeffs` over the domain's
 * subgroup, natural order (`codeword[i] = A(ωⁱ)`).
 */
inline std::vector<fr> rs_encode(std::span<const fr> coeffs, const EvaluationDomain<fr>& domain)
{
    BB_ASSERT_LTE(coeffs.size(), domain.size, "coefficients exceed codeword length");
    std::vector<fr> scratch(domain.size, fr::zero());
    std::copy(coeffs.begin(), coeffs.end(), scratch.begin());
    std::vector<fr> codeword(domain.size);
    polynomial_arithmetic::fft_inner_parallel<fr>(
        scratch.data(), codeword.data(), domain, domain.root, domain.get_round_roots());
    return codeword;
}

/**
 * @brief Stride-2 Lagrange fold: a[j] <- (1-α)·a[2j] + α·a[2j+1]; halves the array.
 * @details Simultaneously partial evaluation of the array-as-MLE at X₁ = α and the coefficient fold
 * of the array-as-univariate (README.md §3). Writes go to a fresh buffer: an in-place parallel fold
 * would race, since one thread's output range overlaps another thread's input range.
 */
inline void fold_array_in_place(std::vector<fr>& array, const fr& alpha)
{
    const size_t half = array.size() / 2;
    BB_ASSERT_EQ(array.size(), half * 2, "fold requires even-length array");
    std::vector<fr> folded(half);
    parallel_for_range(half, [&](size_t start, size_t end) {
        for (size_t j = start; j < end; ++j) {
            folded[j] = array[2 * j] + alpha * (array[2 * j + 1] - array[2 * j]);
        }
    });
    array = std::move(folded);
}

/**
 * @brief Fold a 2^k coset of codeword values down to the folded polynomial's evaluation, README.md
 * eq. (3.1) applied k times.
 *
 * @param values coset values in stride order: values[t] = A(x_base·ηᵗ), η a primitive 2^k-th root
 * @param x_base_inv inverse of the first coset point
 * @param eta_inv inverse of η = ω^{N/2^k}, the primitive 2^k-th root of unity
 * @param alphas the k fold challenges, applied low variable first
 * @return Fold(A, α)(x_base^{2^k}) — inversion-free given x_base_inv
 */
inline fr fold_coset(std::span<const fr> values, const fr& x_base_inv, const fr& eta_inv, std::span<const fr> alphas)
{
    static const fr two_inv = fr(2).invert();
    size_t size = values.size();
    BB_ASSERT_EQ(size, size_t(1) << alphas.size(), "coset size must be 2^k");
    std::vector<fr> current(values.begin(), values.end());
    fr base_inv = x_base_inv;
    fr level_eta_inv = eta_inv;
    for (const fr& alpha : alphas) {
        const size_t half = size / 2;
        fr point_inv = base_inv; // (x_base·ηᵗ)⁻¹, updated multiplicatively over t
        for (size_t t = 0; t < half; ++t) {
            const fr& v_pos = current[t];        // A(+x_t)
            const fr& v_neg = current[t + half]; // A(-x_t), since η^{half} = -1 at every level
            const fr even = (v_pos + v_neg) * two_inv;
            const fr odd = (v_pos - v_neg) * two_inv * point_inv;
            current[t] = even + alpha * (odd - even);
            point_inv *= level_eta_inv;
        }
        base_inv = base_inv.sqr();
        level_eta_inv = level_eta_inv.sqr();
        size = half;
    }
    return current[0];
}

} // namespace bb::whir
