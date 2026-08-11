#include "blake3vm_prover.hpp"

#include "barretenberg/commitment_schemes/gemini/gemini.hpp"
#include "barretenberg/commitment_schemes/shplonk/shplemini.hpp"
#include "barretenberg/honk/proof_system/logderivative_library.hpp"
#include "barretenberg/sumcheck/sumcheck.hpp"

namespace bb {

Blake3VMProver::Blake3VMProver(const CircuitBuilder& builder, const std::shared_ptr<Transcript>& transcript)
    : key(std::make_shared<ProvingKey>(builder))
    , verification_key(std::make_shared<VerificationKey>(*key))
    , transcript(transcript)
{}

HonkProof Blake3VMProver::construct_proof()
{
    using Curve = Flavor::Curve;
    using PCS = Flavor::PCS;
    using PolynomialBatcher = GeminiProver_<Curve>::PolynomialBatcher;

    auto& polynomials = key->polynomials;
    const size_t circuit_size = key->circuit_size;
    const size_t log_n = key->log_circuit_size;

    // Bind the VK: circuit size and precomputed commitments.
    transcript->add_to_hash_buffer("vk_circuit_size", FF(circuit_size));
    for (auto [label, commitment] :
         zip_view(Flavor::PrecomputedEntities<FF>::get_labels(), verification_key->get_all())) {
        transcript->add_to_hash_buffer("vk_" + label, commitment);
    }

    // Wire commitments.
    {
        auto batch = key->commitment_key.start_batch();
        const auto labels = Flavor::WitnessEntities<FF>::get_wire_labels();
        for (auto [wire, label] : zip_view(polynomials.get_wires(), labels)) {
            batch.add_to_batch(wire, label);
        }
        batch.commit_and_send_to_verifier(transcript);
    }

    // Logup challenges and inverse columns.
    const FF beta = transcript->template get_challenge<FF>("beta");
    const FF gamma = transcript->template get_challenge<FF>("gamma");
    relation_parameters.beta = beta;
    relation_parameters.beta_sqr = beta * beta;
    relation_parameters.gamma = gamma;

    bb::constexpr_for<0, Flavor::Layout::NUM_LOOKUP_SETS, 1>([&]<size_t SET>() {
        compute_logderivative_inverse<FF, Blake3VMLookupRelation<FF, SET>, Flavor::ProverPolynomials, true>(
            polynomials, relation_parameters, 0);
    });
    compute_logderivative_inverse<FF, Blake3VMLinkRelation<FF>, Flavor::ProverPolynomials, true>(
        polynomials, relation_parameters, 0);
    {
        auto batch = key->commitment_key.start_batch();
        for (auto [inverse, label] :
             zip_view(polynomials.get_inverses(), Flavor::DerivedWitnessEntities<FF>::get_labels())) {
            batch.add_to_batch(inverse, label);
        }
        batch.commit_and_send_to_verifier(transcript);
    }

    // Sumcheck.
    const FF alpha = transcript->template get_challenge<FF>("Sumcheck:alpha");
    const std::vector<FF> gate_challenges =
        transcript->template get_dyadic_powers_of_challenge<FF>("Sumcheck:gate_challenge", log_n);
    SumcheckProver<Flavor> sumcheck(
        circuit_size, polynomials, transcript, alpha, gate_challenges, relation_parameters, log_n);
    SumcheckOutput<Flavor> sumcheck_output = sumcheck.prove();

    // Shplemini + KZG opening.
    PolynomialBatcher polynomial_batcher(circuit_size);
    polynomial_batcher.set_unshifted(polynomials.get_unshifted());
    polynomial_batcher.set_to_be_shifted_by_one(polynomials.get_to_be_shifted());
    const auto opening_claim = ShpleminiProver_<Curve>::prove(
        circuit_size, polynomial_batcher, sumcheck_output.challenge, key->commitment_key, transcript);
    PCS::compute_opening_proof(key->commitment_key, opening_claim, transcript);

    return transcript->export_proof();
}

} // namespace bb
