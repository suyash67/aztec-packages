#include "blake3vm_verifier.hpp"

#include "barretenberg/commitment_schemes/shplonk/shplemini.hpp"
#include "barretenberg/commitment_schemes/small_subgroup_ipa/small_subgroup_ipa_utils.hpp"
#include "barretenberg/sumcheck/sumcheck.hpp"

namespace bb {

Blake3VMVerifier::Blake3VMVerifier(const std::shared_ptr<VerificationKey>& verification_key,
                                   const std::shared_ptr<Transcript>& transcript)
    : verification_key(verification_key)
    , transcript(transcript)
{}

bool Blake3VMVerifier::verify_proof(const HonkProof& proof)
{
    using Curve = Flavor::Curve;
    using PCS = Flavor::PCS;
    using Commitment = Flavor::Commitment;
    using VerifierCommitments = Flavor::VerifierCommitments;
    using Shplemini = ShpleminiVerifier_<Curve, Flavor::HasZK>;
    using ClaimBatcher = ClaimBatcher_<Curve>;
    using ClaimBatch = ClaimBatcher::Batch;

    const size_t log_n = verification_key->log_circuit_size;

    transcript->load_proof(proof);

    // Bind the VK exactly as the prover does.
    transcript->add_to_hash_buffer("vk_circuit_size", FF(verification_key->circuit_size));
    for (auto [label, commitment] :
         zip_view(Flavor::PrecomputedEntities<FF>::get_labels(), verification_key->get_all())) {
        transcript->add_to_hash_buffer("vk_" + label, commitment);
    }

    VerifierCommitments commitments{ verification_key };

    // Wire commitments.
    {
        const auto labels = Flavor::WitnessEntities<FF>::get_wire_labels();
        for (auto [commitment, label] : zip_view(commitments.get_wires(), labels)) {
            commitment = transcript->template receive_from_prover<Commitment>(label);
        }
    }

    // Logup challenges and inverse commitments.
    const FF beta = transcript->template get_challenge<FF>("beta");
    const FF gamma = transcript->template get_challenge<FF>("gamma");
    relation_parameters.beta = beta;
    relation_parameters.beta_sqr = beta * beta;
    relation_parameters.gamma = gamma;
    for (auto [commitment, label] :
         zip_view(commitments.get_inverses(), Flavor::DerivedWitnessEntities<FF>::get_labels())) {
        commitment = transcript->template receive_from_prover<Commitment>(label);
    }

    // Sumcheck.
    const FF alpha = transcript->template get_challenge<FF>("Sumcheck:alpha");
    SumcheckVerifier<Flavor> sumcheck(transcript, alpha, log_n);
    const std::vector<FF> gate_challenges =
        transcript->template get_dyadic_powers_of_challenge<FF>("Sumcheck:gate_challenge", log_n);
    SumcheckOutput<Flavor> sumcheck_output = sumcheck.verify(relation_parameters, gate_challenges);

    // Shplemini + KZG pairing check.
    ClaimBatcher claim_batcher{
        .unshifted = ClaimBatch{ commitments.get_unshifted(), sumcheck_output.claimed_evaluations.get_unshifted() },
        .shifted = ClaimBatch{ commitments.get_to_be_shifted(), sumcheck_output.claimed_evaluations.get_shifted() }
    };
    std::array<Commitment, NUM_SMALL_IPA_COMMITMENTS> libra_commitments = {};
    auto shplemini_output = Shplemini::compute_batch_opening_claim(claim_batcher,
                                                                   sumcheck_output.challenge,
                                                                   Commitment::one(),
                                                                   transcript,
                                                                   RepeatedCommitmentsData{},
                                                                   libra_commitments,
                                                                   sumcheck_output.claimed_libra_evaluation);
    const auto pairing_points =
        PCS::reduce_verify_batch_opening_claim(std::move(shplemini_output.batch_opening_claim), transcript);

    if (!sumcheck_output.verified) {
        info("Blake3VM Verifier: sumcheck failed");
        return false;
    }
    if (!pairing_points.check()) {
        info("Blake3VM Verifier: pairing check failed");
        return false;
    }
    return true;
}

} // namespace bb
