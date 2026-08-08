#pragma once

#include "barretenberg/common/assert.hpp"
#include "barretenberg/ecc/curves/bn254/bn254.hpp"
#include "barretenberg/polynomials/evaluation_domain.hpp"
#include "barretenberg/polynomials/polynomial_arithmetic.hpp"

#include <cstddef>
#include <span>
#include <vector>

namespace bb::fflonk_plonk {

using Curve = curve::BN254;
using FF = bb::fr;
using Commitment = Curve::AffineElement;
using GroupElement = Curve::Element;

/**
 * @brief The interleaved packing of a group of columns: g(X) = sum_{i<t} f_i(X^t) X^i.
 *
 * @details The reason this is the right encoding is the residue identity
 * `g mod (X^t - z) = sum_i f_i(z) X^i`: one opening of `g` certifies every column's evaluation at
 * `z`, and the residue's coefficients are exactly those evaluations. Columns may have different
 * lengths; the packing is sized by the longest.
 */
std::vector<FF> pack_columns(std::span<const std::vector<FF>> columns);

/**
 * @brief Quotient and residue of `poly` by `X^t - c`.
 *
 * @details From `poly[k] = quotient[k-t] - c*quotient[k]` the quotient falls out top-down as
 * `quotient[j] = poly[j+t] + c*quotient[j+t]`, and then `residue[i] = poly[i] + c*quotient[i]`.
 * The residue has exactly `t` coefficients.
 */
void divide_by_power_minus(
    std::span<const FF> poly, size_t t, const FF& c, std::vector<FF>& quotient, std::vector<FF>& residue);

/**
 * @brief Synthetic division of `poly` by `X - r`, returning the quotient and the remainder.
 * @details The remainder is `poly(r)`; callers that require exact division must check it.
 */
std::vector<FF> divide_by_linear(std::span<const FF> poly, const FF& r, FF& remainder);

/**
 * @brief Divide `poly` by the vanishing polynomial `X^n - 1`, requiring exact division.
 * @return false if `poly` is not divisible, in which case the result is meaningless.
 */
bool divide_by_vanishing(std::span<const FF> poly, size_t n, std::vector<FF>& quotient);

/** @brief `poly(x)` by Horner, over however many coefficients the span holds. */
FF evaluate(std::span<const FF> poly, const FF& x);

/**
 * @brief Drop structurally-zero leading coefficients, keeping at least one.
 * @details Division routines here size their output by the dividend rather than by the true degree,
 * so callers that need the degree - to bound a commitment, or to check a split fits - must trim.
 */
void trim(std::vector<FF>& poly);

/** @brief Coefficients -> evaluations over `domain`. `coefficients` may be shorter than the domain. */
std::vector<FF> coefficients_to_evaluations(std::span<const FF> coefficients, const EvaluationDomain<FF>& domain);

/** @brief Evaluations over `domain` -> coefficients. `evaluations` must span the whole domain. */
std::vector<FF> evaluations_to_coefficients(std::span<const FF> evaluations, const EvaluationDomain<FF>& domain);

/** @brief `{1, root, root^2, ...}`, `count` of them. */
std::vector<FF> compute_power_table(const FF& root, size_t count);

} // namespace bb::fflonk_plonk
