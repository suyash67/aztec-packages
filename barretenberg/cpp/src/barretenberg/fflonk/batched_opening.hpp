#pragma once

#include "barretenberg/commitment_schemes/commitment_key.hpp"
#include "barretenberg/commitment_schemes/pairing_points.hpp"
#include "barretenberg/fflonk/polynomial_utils.hpp"

#include <span>
#include <vector>

namespace bb::fflonk_plonk {

/**
 * @brief How one packed group is opened: `pack` interleaved columns, at one point or at two.
 *
 * @details A group opened at `{xi}` certifies each of its columns at `xi`; a group opened at
 * `{xi, xi*omega}` certifies both, and therefore reveals *every* column it contains at both points.
 * Group membership is a zero-knowledge decision as much as a cost one - a column revealed at a point
 * where no verifier identity pins it down needs a blinder of its own. See PROTOCOL.md section 7.
 */
struct GroupShape {
    size_t pack = 1;
    bool two_point = false;

    /** @brief How many field elements this group contributes to the proof. */
    [[nodiscard]] constexpr size_t num_evaluations() const { return two_point ? 2 * pack : pack; }

    /** @brief `deg Z_g`, where `Z_g` vanishes exactly on the group's opening set. */
    [[nodiscard]] constexpr size_t vanishing_degree() const { return num_evaluations(); }
};

/** @brief One group's claimed evaluations. `at_xi_omega` is empty unless the group is two-point. */
struct GroupEvaluations {
    std::vector<FF> at_xi;
    std::vector<FF> at_xi_omega;
};

/**
 * @brief `Z_g(y)` for every group, where `Z_g = X^t - xi` or `(X^t - xi)(X^t - xi*omega)`.
 */
std::vector<FF> group_vanishing_at(std::span<const GroupShape> shapes, const FF& y, const FF& xi, const FF& xi_omega);

/**
 * @brief `R_g(y)`: the residue of the group's packed polynomial modulo `Z_g`, evaluated at `y`.
 *
 * @details For a one-point group the residue's coefficients *are* the claimed evaluations, so this
 * is a Horner evaluation over them. For a two-point group, with
 * `A(Y) = sum_i f_i(xi) Y^i` and `B(Y) = sum_i f_i(xi*omega) Y^i`, the Chinese-remainder interpolant
 * of the two residues is
 *
 * ```
 * R = [ B (Y^t - xi) - A (Y^t - xi*omega) ] / (xi*omega - xi)
 * ```
 *
 * which at `t = 1` is the straight line through `(xi, f(xi))` and `(xi*omega, f(xi*omega))`.
 *
 * @param inverse_xi_omega_minus_xi `1/(xi*omega - xi)`; unused by one-point groups.
 */
FF residue_at(const GroupShape& shape,
              const GroupEvaluations& evaluations,
              const FF& y,
              const FF& xi,
              const FF& xi_omega,
              const FF& inverse_xi_omega_minus_xi);

/**
 * @brief Rounds 3 to 5 of the protocol: open every group at its own point set, then batch every
 * opening into the two group elements `W` and `W'`.
 *
 * @details The caller drives the transcript, because the order in which the evaluations are absorbed
 * is part of whichever protocol is using this - so the flow is: construct (which divides out the
 * residues), absorb the evaluations and squeeze `nu`, `commit_w`, absorb `W` and squeeze `y`,
 * `commit_w_prime`.
 *
 * Soundness with a product `Z_T = prod_g Z_g` rather than the vanishing polynomial of the union of
 * the point sets - which double-counts whenever two groups share a point - is argued in PROTOCOL.md
 * section 8. It rests on `nu` being drawn after the claimed evaluations are fixed.
 */
class BatchedOpeningProver {
  public:
    /**
     * @param packed_groups one packed polynomial per group, in the order the `nu` powers apply.
     * @param shapes parallel to `packed_groups`.
     */
    BatchedOpeningProver(std::span<const std::vector<FF>* const> packed_groups,
                         std::span<const GroupShape> shapes,
                         const FF& xi,
                         const FF& xi_omega);

    /** @brief The claimed evaluations, in group order. Read them out to build the proof. */
    [[nodiscard]] const std::vector<GroupEvaluations>& evaluations() const { return evaluations_; }

    /** @brief `W = sum_g nu^g (g_g - R_g) / Z_g`. Call once `nu` has been squeezed. */
    [[nodiscard]] std::vector<FF> compute_w(const FF& nu);

    /**
     * @brief `W' = L / (X - y)`, with `L = sum_g nu^g Gamma_g (g_g - R_g(y)) - Z_T(y) W`.
     * @details Aborts if `y` collides with an opening point, which would make `Gamma_g` a division
     * by zero, or if `L` does not vanish at `y`, which would mean this file is wrong.
     */
    [[nodiscard]] std::vector<FF> compute_w_prime(const FF& nu, const FF& y, const std::vector<FF>& w_polynomial);

  private:
    std::span<const std::vector<FF>* const> packed_groups_;
    std::span<const GroupShape> shapes_;
    FF xi_;
    FF xi_omega_;
    std::vector<GroupEvaluations> evaluations_;
    std::vector<std::vector<FF>> quotients_;
};

/** @brief Everything the verifier's fold needs besides the group vanishing values. */
struct OpeningClaim {
    std::span<const Commitment> commitments;
    std::span<const GroupShape> shapes;
    std::span<const GroupEvaluations> evaluations;
    Commitment w = Commitment::infinity();
    Commitment w_prime = Commitment::infinity();
    FF xi = FF::zero();
    FF xi_omega = FF::zero();
    FF nu = FF::zero();
    FF y = FF::zero();
    FF inverse_xi_omega_minus_xi = FF::zero();
};

/**
 * @brief Fold the group commitments into the two points of the pairing check.
 *
 * @details `F = sum_g s_g C_g - (sum_g s_g R_g(y)) [1] - Z_T(y) W + y W'` with
 * `s_g = nu^g Z_T(y)/Z_g(y)`, and the check is `e(F, [1]_2) = e(W', [x]_2)`.
 *
 * The vanishing values are passed in rather than computed here so the caller can fold them into
 * whatever batched inversion it already runs; `inverse_group_vanishing[g]` must be `1/Z_g(y)`, and
 * the caller is responsible for having rejected a zero.
 */
[[nodiscard]] PairingPoints<Curve> fold_opening(const OpeningClaim& claim,
                                                std::span<const FF> group_vanishing,
                                                std::span<const FF> inverse_group_vanishing);

} // namespace bb::fflonk_plonk
