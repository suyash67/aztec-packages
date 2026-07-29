#pragma once

#include "barretenberg/ecc/curves/bn254/fr.hpp"
#include "barretenberg/ecc/curves/bn254/g1.hpp"

#include <span>
#include <vector>

namespace bb::cq {

/**
 * @brief Feist-Khovratovich batch computation of KZG opening proofs, plus Lagrange-basis SRS transforms.
 *
 * @details Given a polynomial T of degree < N and a monomial SRS { [tau^d]_1 }, the commitments to all N KZG
 * quotients U_i(X) = (T(X) - T(omega^i)) / (X - omega^i) over a subgroup of order N are the size-N group DFT of the
 * "Toeplitz products" h_r = sum_d c_{r+1+d} [tau^d]_1, which are themselves computable with a size-2N group FFT
 * (https://eprint.iacr.org/2023/033, originally Feist-Khovratovich). Total cost O(N log N) group operations
 * instead of the naive N MSMs of size N.
 */

/**
 * @brief h_r = sum_{d >= 0, r+1+d < N} coeffs[r+1+d] * srs[d] for r in [0, N), where N = coeffs.size().
 * Naive O(N^2) reference implementation used for testing.
 */
std::vector<g1::element> toeplitz_srs_products_naive(std::span<const fr> coeffs,
                                                     std::span<const g1::affine_element> srs);

/**
 * @brief Same as toeplitz_srs_products_naive but via a size-2N cyclic convolution: O(N log N) group operations.
 */
std::vector<g1::element> toeplitz_srs_products(std::span<const fr> coeffs, std::span<const g1::affine_element> srs);

/**
 * @brief Commitments [U_i]_1 to the KZG quotients U_i(X) = (T(X) - T(omega^i)) / (X - omega^i) for all i in [0, N),
 * where T is given by its N monomial coefficients and omega generates the subgroup of order N.
 */
std::vector<g1::affine_element> fk_all_kzg_opening_proofs(std::span<const fr> t_coeffs,
                                                          std::span<const g1::affine_element> srs);

/**
 * @brief Commitments [L_i]_1 to the Lagrange basis polynomials of the subgroup of order domain_size.
 * @details L_i(X) = (1/N) sum_k omega^{-ik} X^k, so the vector ([L_0], ..., [L_{N-1}]) is the group IFFT of the
 * first N SRS points.
 */
std::vector<g1::affine_element> lagrange_basis_commitments(std::span<const g1::affine_element> srs, size_t domain_size);

/**
 * @brief Commitments [(L_i(X) - L_i(0)) / X]_1 for all i, used for KZG openings at zero of polynomials given in
 * Lagrange basis. (L_i(X) - 1/N)/X = (omega^{-i}/N) sum_{k=0}^{N-2} omega^{-ik} X^k.
 */
std::vector<g1::affine_element> lagrange_zero_quotient_commitments(std::span<const g1::affine_element> srs,
                                                                   size_t domain_size);

} // namespace bb::cq
