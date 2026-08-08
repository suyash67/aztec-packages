#include "barretenberg/fflonk/verifier.hpp"

#include "barretenberg/commitment_schemes/pairing_points.hpp"
#include "barretenberg/fflonk/transcript.hpp"

#include <array>

namespace bb::fflonk_plonk {

namespace {

/**
 * @brief On the curve, and not the point at infinity.
 * @details BN254's G1 has prime order, so being on the curve implies being in the subgroup and no
 * cofactor check is needed. Infinity is rejected rather than accepted: an honest prover never
 * produces one - it would mean a committed polynomial was identically zero - and rejecting keeps
 * this verifier byte-for-byte equivalent to the Solidity one, whose `y^2 = x^3 + 3` test has no
 * natural encoding for it.
 */
bool is_valid_point(const Commitment& point)
{
    return !point.is_point_at_infinity() && point.on_curve();
}

} // namespace

VerificationReport verify_detailed(const VerificationKey& key, const Proof& proof, std::span<const FF> public_inputs)
{
    VerificationReport report;

    const size_t n = key.circuit_size;
    if (n == 0 || (n & (n - 1)) != 0) {
        return report;
    }
    if (public_inputs.size() != key.num_public_inputs || public_inputs.size() >= n) {
        return report;
    }
    for (const Commitment& point : { key.c0, proof.c1, proof.c2, proof.c3, proof.w, proof.w_prime }) {
        if (!is_valid_point(point)) {
            return report;
        }
    }

    const auto& e = proof.evaluations;

    // -------------------------------------------------------------------------------------------
    // Fiat-Shamir, replayed in the order the prover wrote it.
    // -------------------------------------------------------------------------------------------
    Transcript transcript;
    transcript.absorb(key.hash());
    for (const FF& public_input : public_inputs) {
        transcript.absorb(public_input);
    }
    transcript.absorb(proof.c1);
    const FF beta = transcript.squeeze();
    const FF gamma = transcript.squeeze();
    transcript.absorb(proof.c2);
    transcript.absorb(proof.c3);
    const FF xi = transcript.squeeze();
    for (const FF& evaluation : e) {
        transcript.absorb(evaluation);
    }
    const FF nu = transcript.squeeze();
    transcript.absorb(proof.w);
    const FF y = transcript.squeeze();

    const FF xi_pow_n = xi.pow(static_cast<uint64_t>(n));
    // xi must miss the domain, or the vanishing polynomial is zero and every quotient claim is free.
    if (xi.is_zero() || xi_pow_n == FF::one()) {
        return report;
    }
    const FF vanishing = xi_pow_n - FF::one();
    const FF xi_omega = xi * key.omega;

    // -------------------------------------------------------------------------------------------
    // One batched inversion for the whole protocol: the Lagrange denominators, the interpolant's
    // slope, and the four Shplonk denominators.
    // -------------------------------------------------------------------------------------------
    const size_t num_public = public_inputs.size();
    const std::vector<FF> group_vanishing = group_vanishing_at(GROUP_SHAPES, y, xi, xi_omega);

    // Layout: [ xi - w^i for i < l ] [ xi - 1 ] [ xi*w - xi ] [ Z_g(y) for g < 4 ]
    const size_t offset_lagrange_first = num_public;
    const size_t offset_interpolant = num_public + 1;
    const size_t offset_groups = num_public + 2;
    std::vector<FF> denominators(num_public + 2 + NUM_GROUPS, FF::zero());
    {
        FF root_power = FF::one();
        for (size_t i = 0; i < num_public; ++i) {
            denominators[i] = xi - root_power;
            root_power *= key.omega;
        }
    }
    denominators[offset_lagrange_first] = xi - FF::one();
    denominators[offset_interpolant] = xi_omega - xi;
    for (size_t g = 0; g < NUM_GROUPS; ++g) {
        denominators[offset_groups + g] = group_vanishing[g];
    }
    for (const FF& denominator : denominators) {
        // xi is off the domain and xi != xi*omega, so only y colliding with an opening point can
        // land here; that is a malformed proof, not a soundness break, and it is rejected.
        if (denominator.is_zero()) {
            return report;
        }
    }
    FF::batch_invert(denominators);

    report.well_formed = true;

    const FF n_inverse = FF(static_cast<uint64_t>(n)).invert();
    const FF lagrange_first = vanishing * n_inverse * denominators[offset_lagrange_first];

    // PI(xi) = - sum_i x_i L_i(xi), with L_i(X) = w^i (X^n - 1) / (n (X - w^i)).
    FF public_input_evaluation = FF::zero();
    {
        FF root_power = FF::one();
        for (size_t i = 0; i < num_public; ++i) {
            public_input_evaluation -= public_inputs[i] * root_power * vanishing * n_inverse * denominators[i];
            root_power *= key.omega;
        }
    }

    // -------------------------------------------------------------------------------------------
    // The three constraint identities, entirely in the field over the claimed evaluations. This is
    // what replaces PlonK's linearisation, and why the verifier's group work does not grow with the
    // number of selectors.
    // -------------------------------------------------------------------------------------------
    {
        const FF gate = e[EVAL_Q_L] * e[EVAL_A] + e[EVAL_Q_R] * e[EVAL_B] + e[EVAL_Q_O] * e[EVAL_C] +
                        e[EVAL_Q_M] * e[EVAL_A] * e[EVAL_B] + e[EVAL_Q_C] + public_input_evaluation;
        report.gate_identity = gate == (e[EVAL_T0_LO] + xi_pow_n * e[EVAL_T0_HI]) * vanishing;
    }
    report.grand_product_start = (e[EVAL_Z] - FF::one()) * lagrange_first == e[EVAL_T1] * vanishing;
    {
        const FF identity_side = (e[EVAL_A] + beta * xi + gamma) * (e[EVAL_B] + beta * key.k1 * xi + gamma) *
                                 (e[EVAL_C] + beta * key.k2 * xi + gamma) * e[EVAL_Z];
        const FF sigma_side = (e[EVAL_A] + beta * e[EVAL_S_1] + gamma) * (e[EVAL_B] + beta * e[EVAL_S_2] + gamma) *
                              (e[EVAL_C] + beta * e[EVAL_S_3] + gamma) * e[EVAL_Z_OMEGA];
        const FF quotient = e[EVAL_T2_LO] + xi_pow_n * e[EVAL_T2_MID] + xi_pow_n * xi_pow_n * e[EVAL_T2_HI];
        report.permutation_identity = identity_side - sigma_side == quotient * vanishing;
    }

    // -------------------------------------------------------------------------------------------
    // The batched opening. R_g is the residue of the group's packed polynomial modulo Z_g, so its
    // coefficients are the claimed evaluations themselves; only the grand product, opened at two
    // points, needs an interpolant, and that one is linear.
    // -------------------------------------------------------------------------------------------
    const std::array<GroupEvaluations, NUM_GROUPS> group_evaluations = {
        GroupEvaluations{ std::vector<FF>(e.begin() + EVAL_Q_L, e.begin() + EVAL_Q_L + PACK_PREPROCESSED), {} },
        GroupEvaluations{ std::vector<FF>(e.begin() + EVAL_A, e.begin() + EVAL_A + PACK_WIRES), {} },
        GroupEvaluations{ { e[EVAL_Z] }, { e[EVAL_Z_OMEGA] } },
        GroupEvaluations{ std::vector<FF>(e.begin() + EVAL_T1, e.begin() + EVAL_T1 + PACK_QUOTIENTS), {} },
    };
    const std::array<Commitment, NUM_GROUPS> commitments = { key.c0, proof.c1, proof.c2, proof.c3 };

    const OpeningClaim claim{ .commitments = commitments,
                              .shapes = GROUP_SHAPES,
                              .evaluations = group_evaluations,
                              .w = proof.w,
                              .w_prime = proof.w_prime,
                              .xi = xi,
                              .xi_omega = xi_omega,
                              .nu = nu,
                              .y = y,
                              .inverse_xi_omega_minus_xi = denominators[offset_interpolant] };

    report.opening =
        fold_opening(claim, group_vanishing, std::span<const FF>(denominators).subspan(offset_groups, NUM_GROUPS))
            .check();

    return report;
}

bool verify(const VerificationKey& key, const Proof& proof, std::span<const FF> public_inputs)
{
    return verify_detailed(key, proof, public_inputs).accepted();
}

bool verify(const VerificationKey& key, std::span<const uint8_t> proof_bytes, std::span<const FF> public_inputs)
{
    Proof proof;
    if (!Proof::from_buffer(proof_bytes, proof)) {
        return false;
    }
    return verify(key, proof, public_inputs);
}

} // namespace bb::fflonk_plonk
