#include "barretenberg/fflonk/batched_opening.hpp"

#include "barretenberg/common/throw_or_abort.hpp"

#include <algorithm>

namespace bb::fflonk_plonk {

namespace {

/** @brief `y^t - c`, the vanishing polynomial of the `t`-th roots of `c`, evaluated at `y`. */
FF power_minus(const FF& y, const size_t t, const FF& c)
{
    return y.pow(static_cast<uint64_t>(t)) - c;
}

/** @brief The top `t` coefficients of a `divide_by_power_minus` quotient are structurally zero. */
void drop_structural_zeros(std::vector<FF>& quotient, const size_t t)
{
    quotient.resize(quotient.size() > t ? quotient.size() - t : static_cast<size_t>(0));
}

} // namespace

std::vector<FF> group_vanishing_at(std::span<const GroupShape> shapes, const FF& y, const FF& xi, const FF& xi_omega)
{
    std::vector<FF> vanishing(shapes.size(), FF::zero());
    for (size_t g = 0; g < shapes.size(); ++g) {
        const FF at_xi = power_minus(y, shapes[g].pack, xi);
        vanishing[g] = shapes[g].two_point ? at_xi * power_minus(y, shapes[g].pack, xi_omega) : at_xi;
    }
    return vanishing;
}

FF residue_at(const GroupShape& shape,
              const GroupEvaluations& evaluations,
              const FF& y,
              const FF& xi,
              const FF& xi_omega,
              const FF& inverse_xi_omega_minus_xi)
{
    const FF a = evaluate(evaluations.at_xi, y);
    if (!shape.two_point) {
        return a;
    }
    const FF b = evaluate(evaluations.at_xi_omega, y);
    const FF y_pow_t = y.pow(static_cast<uint64_t>(shape.pack));
    return (b * (y_pow_t - xi) - a * (y_pow_t - xi_omega)) * inverse_xi_omega_minus_xi;
}

BatchedOpeningProver::BatchedOpeningProver(std::span<const std::vector<FF>* const> packed_groups,
                                           std::span<const GroupShape> shapes,
                                           const FF& xi,
                                           const FF& xi_omega)
    : packed_groups_(packed_groups)
    , shapes_(shapes)
    , xi_(xi)
    , xi_omega_(xi_omega)
{
    BB_ASSERT_EQ(packed_groups.size(), shapes.size(), "one shape per packed group");
    if (xi == xi_omega) {
        throw_or_abort("fflonk: the two opening points coincide");
    }
    const FF inverse_gap = (xi_omega - xi).invert();

    evaluations_.resize(shapes.size());
    quotients_.resize(shapes.size());

    for (size_t g = 0; g < shapes.size(); ++g) {
        const size_t t = shapes_[g].pack;
        const std::vector<FF>& packed = *packed_groups_[g];

        std::vector<FF> quotient;
        divide_by_power_minus(packed, t, xi_, quotient, evaluations_[g].at_xi);
        drop_structural_zeros(quotient, t);

        if (!shapes_[g].two_point) {
            quotients_[g] = std::move(quotient);
            continue;
        }

        // g - R is divisible by (X^t - xi) with cofactor `quotient + (A - B)/(xi*omega - xi)`, and
        // that cofactor is in turn divisible by (X^t - xi*omega): reducing `quotient` modulo the
        // second factor gives exactly (B - A)/(xi*omega - xi), which the correction cancels.
        std::vector<FF> discarded;
        divide_by_power_minus(packed, t, xi_omega_, discarded, evaluations_[g].at_xi_omega);
        discarded = {};

        quotient.resize(std::max(quotient.size(), t), FF::zero());
        for (size_t i = 0; i < t; ++i) {
            quotient[i] += (evaluations_[g].at_xi[i] - evaluations_[g].at_xi_omega[i]) * inverse_gap;
        }

        std::vector<FF> residue;
        divide_by_power_minus(quotient, t, xi_omega_, quotients_[g], residue);
        drop_structural_zeros(quotients_[g], t);
        for (const FF& coefficient : residue) {
            if (!coefficient.is_zero()) {
                throw_or_abort("fflonk: a two-point group's opening does not divide exactly");
            }
        }
    }
}

std::vector<FF> BatchedOpeningProver::compute_w(const FF& nu)
{
    size_t widest = 0;
    for (const std::vector<FF>& quotient : quotients_) {
        widest = std::max(widest, quotient.size());
    }

    std::vector<FF> w_polynomial(widest, FF::zero());
    FF nu_power = FF::one();
    for (const std::vector<FF>& quotient : quotients_) {
        for (size_t i = 0; i < quotient.size(); ++i) {
            w_polynomial[i] += nu_power * quotient[i];
        }
        nu_power *= nu;
    }
    return w_polynomial;
}

std::vector<FF> BatchedOpeningProver::compute_w_prime(const FF& nu, const FF& y, const std::vector<FF>& w_polynomial)
{
    std::vector<FF> vanishing_at_y = group_vanishing_at(shapes_, y, xi_, xi_omega_);
    FF total_vanishing = FF::one();
    for (const FF& value : vanishing_at_y) {
        if (value.is_zero()) {
            throw_or_abort("fflonk: the batching challenge collides with an opening point");
        }
        total_vanishing *= value;
    }
    std::vector<FF> inverse_vanishing = vanishing_at_y;
    FF::batch_invert(inverse_vanishing);
    const FF inverse_gap = (xi_omega_ - xi_).invert();

    size_t widest = w_polynomial.size();
    for (const std::vector<FF>* packed : packed_groups_) {
        widest = std::max(widest, packed->size());
    }

    std::vector<FF> linearization(widest, FF::zero());
    FF constant_term = FF::zero();
    FF nu_power = FF::one();
    for (size_t g = 0; g < shapes_.size(); ++g) {
        const FF scalar = nu_power * total_vanishing * inverse_vanishing[g];
        const std::vector<FF>& packed = *packed_groups_[g];
        for (size_t i = 0; i < packed.size(); ++i) {
            linearization[i] += scalar * packed[i];
        }
        constant_term += scalar * residue_at(shapes_[g], evaluations_[g], y, xi_, xi_omega_, inverse_gap);
        nu_power *= nu;
    }
    linearization[0] -= constant_term;
    for (size_t i = 0; i < w_polynomial.size(); ++i) {
        linearization[i] -= total_vanishing * w_polynomial[i];
    }

    FF remainder = FF::zero();
    std::vector<FF> w_prime = divide_by_linear(linearization, y, remainder);
    if (!remainder.is_zero()) {
        throw_or_abort("fflonk: the batched opening polynomial does not vanish at y");
    }
    return w_prime;
}

PairingPoints<Curve> fold_opening(const OpeningClaim& claim,
                                  std::span<const FF> group_vanishing,
                                  std::span<const FF> inverse_group_vanishing)
{
    const size_t num_groups = claim.shapes.size();
    BB_ASSERT_EQ(claim.commitments.size(), num_groups, "one commitment per group");
    BB_ASSERT_EQ(claim.evaluations.size(), num_groups, "one evaluation set per group");
    BB_ASSERT_EQ(group_vanishing.size(), num_groups, "one vanishing value per group");
    BB_ASSERT_EQ(inverse_group_vanishing.size(), num_groups, "one inverse vanishing value per group");

    FF total_vanishing = FF::one();
    for (const FF& value : group_vanishing) {
        total_vanishing *= value;
    }

    GroupElement accumulator = GroupElement::infinity();
    FF constant_term = FF::zero();
    FF nu_power = FF::one();
    for (size_t g = 0; g < num_groups; ++g) {
        const FF scalar = nu_power * total_vanishing * inverse_group_vanishing[g];
        accumulator += GroupElement(claim.commitments[g]) * scalar;
        constant_term += scalar * residue_at(claim.shapes[g],
                                             claim.evaluations[g],
                                             claim.y,
                                             claim.xi,
                                             claim.xi_omega,
                                             claim.inverse_xi_omega_minus_xi);
        nu_power *= claim.nu;
    }
    accumulator -= GroupElement(Commitment::one()) * constant_term;
    accumulator -= GroupElement(claim.w) * total_vanishing;
    accumulator += GroupElement(claim.w_prime) * claim.y;

    return PairingPoints<Curve>{ Commitment(accumulator), Commitment(-GroupElement(claim.w_prime)) };
}

} // namespace bb::fflonk_plonk
