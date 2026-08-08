#include "barretenberg/ultra_fflonk/verifier.hpp"

#include "barretenberg/fflonk/encoding.hpp"
#include "barretenberg/fflonk/transcript.hpp"
#include "barretenberg/honk/library/grand_product_delta.hpp"
#include "barretenberg/ultra_fflonk/relation_batch.hpp"

#include <array>

namespace bb::ultra_fflonk {

VerificationReport verify_detailed(const VerificationKey& key, const Proof& proof, std::span<const FF> public_inputs)
{
    VerificationReport report;

    const size_t n = key.circuit_size;
    if (n < MIN_CIRCUIT_SIZE || (n & (n - 1)) != 0) {
        return report;
    }
    if (public_inputs.size() != key.num_public_inputs || public_inputs.size() >= n) {
        return report;
    }
    for (const Commitment& point : key.preprocessed) {
        if (!fflonk_plonk::is_valid_point(point)) {
            return report;
        }
    }
    for (const Commitment& point :
         { proof.wires, proof.memory, proof.grand_product, proof.quotient, proof.w, proof.w_prime }) {
        if (!fflonk_plonk::is_valid_point(point)) {
            return report;
        }
    }

    // ---------------------------------------------------------------------------------------------
    // Fiat-Shamir, replayed in the order the prover wrote it: one group per Oink round.
    // ---------------------------------------------------------------------------------------------
    fflonk_plonk::Transcript transcript;
    transcript.absorb(key.hash());
    for (const FF& public_input : public_inputs) {
        transcript.absorb(public_input);
    }

    RelationParameters<FF> parameters;
    transcript.absorb(proof.wires);
    parameters.compute_eta_powers(transcript.squeeze());
    parameters.rom_logup_gamma = transcript.squeeze();

    transcript.absorb(proof.memory);
    parameters.compute_beta_powers(transcript.squeeze());
    parameters.gamma = transcript.squeeze();

    transcript.absorb(proof.grand_product);
    const FF alpha = transcript.squeeze();

    transcript.absorb(proof.quotient);
    const FF xi = transcript.squeeze();

    for (const FF& evaluation : proof.evaluations) {
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
    const FF xi_omega = xi * key.omega;

    // The public inputs enter the permutation argument through this correction term, which is the
    // only place the verifier binds them to the witness.
    parameters.public_input_delta = compute_public_input_delta<Flavor>(
        public_inputs, parameters.beta, parameters.gamma, FF(static_cast<uint64_t>(key.pub_inputs_offset)));

    // One batched inversion for the whole protocol: the eight Shplonk denominators and the gap
    // between the two opening points.
    const std::vector<FF> group_vanishing = fflonk_plonk::group_vanishing_at(GROUP_SHAPES, y, xi, xi_omega);
    std::vector<FF> denominators = group_vanishing;
    denominators.push_back(xi_omega - xi);
    for (const FF& denominator : denominators) {
        // xi is off the domain and xi != xi*omega, so only y colliding with an opening point can
        // land here; that is a malformed proof, not a soundness break, and it is rejected.
        if (denominator.is_zero()) {
            return report;
        }
    }
    FF::batch_invert(denominators);

    report.well_formed = true;

    // ---------------------------------------------------------------------------------------------
    // The quotient identity, entirely in the field over the claimed evaluations. The relations are
    // the flavor's own, called through the same accumulators the prover divided out, so there is one
    // arithmetization rather than a prover's and a verifier's copy of it.
    // ---------------------------------------------------------------------------------------------
    const std::vector<GroupEvaluations> evaluations = unflatten_evaluations(proof.evaluations);
    {
        const RunningSumEvaluations running_sums = to_running_sum_evaluations(evaluations);
        const FF numerator = quotient_numerator(
            to_all_values(evaluations), parameters, alpha, running_sums.at_xi, running_sums.at_xi_omega);

        FF quotient_at_xi = FF::zero();
        for (size_t chunk = NUM_QUOTIENT_CHUNKS; chunk-- > 0;) {
            quotient_at_xi = (quotient_at_xi * xi_pow_n) + evaluations[GROUP_QUOTIENT].at_xi[chunk];
        }
        report.quotient = numerator == quotient_at_xi * (xi_pow_n - FF::one());
    }

    // ---------------------------------------------------------------------------------------------
    // The batched opening.
    // ---------------------------------------------------------------------------------------------
    std::array<Commitment, NUM_GROUPS> commitments{};
    for (size_t g = 0; g < NUM_PREPROCESSED_GROUPS; ++g) {
        commitments[g] = key.preprocessed[g];
    }
    const std::array<Commitment, 4> committed = proof.committed_groups();
    for (size_t g = 0; g < committed.size(); ++g) {
        commitments[NUM_PREPROCESSED_GROUPS + g] = committed[g];
    }

    const fflonk_plonk::OpeningClaim claim{ .commitments = commitments,
                                            .shapes = GROUP_SHAPES,
                                            .evaluations = evaluations,
                                            .w = proof.w,
                                            .w_prime = proof.w_prime,
                                            .xi = xi,
                                            .xi_omega = xi_omega,
                                            .nu = nu,
                                            .y = y,
                                            .inverse_xi_omega_minus_xi = denominators[NUM_GROUPS] };

    report.opening =
        fflonk_plonk::fold_opening(claim, group_vanishing, std::span<const FF>(denominators).subspan(0, NUM_GROUPS))
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

} // namespace bb::ultra_fflonk
