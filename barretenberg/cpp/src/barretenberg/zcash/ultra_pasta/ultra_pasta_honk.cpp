#include "barretenberg/zcash/ultra_pasta/ultra_pasta_honk.hpp"
#include "barretenberg/commitment_schemes/claim_batcher.hpp"
#include "barretenberg/commitment_schemes/shplonk/shplemini.hpp"
#include "barretenberg/commitment_schemes/small_subgroup_ipa/small_subgroup_ipa_impl.hpp"
#include "barretenberg/common/bb_bench.hpp"
#include "barretenberg/flavor/verifier_commitments.hpp"
#include "barretenberg/honk/composer/composer_lib.hpp"
#include "barretenberg/sumcheck/sumcheck.hpp"
#include "barretenberg/trace_to_polynomials/trace_to_polynomials_impl.hpp"
#include "barretenberg/ultra_honk/oink_prover_impl.hpp"
#include "barretenberg/ultra_honk/oink_verifier_impl.hpp"
#include "barretenberg/ultra_honk/prover_instance_impl.hpp"
#include "barretenberg/ultra_honk/verifier_instance.hpp"
#include "barretenberg/zcash/honk/halo2_ipa.hpp"

namespace bb {
/**
 * @brief The table polynomials of an UltraPastaZKFlavor circuit, from the F_p rows of its tables (the BasicTables only
 * index the rows).
 */
template <>
void construct_lookup_table_polynomials<zcash::UltraPastaZKFlavor>(
    const RefArray<zcash::UltraPastaZKFlavor::Polynomial, 4>& table_polynomials,
    const zcash::UltraPastaCircuitBuilder& circuit)
{
    using FF = zcash::UltraPastaZKFlavor::FF;
    size_t offset = circuit.blocks.lookup.trace_offset();
    for (const auto& table : circuit.get_lookup_tables()) {
        const auto& rows = circuit.pasta_table_rows(table.table_index);
        BB_ASSERT_EQ(rows.size(), table.size());
        for (const auto& row : rows) {
            table_polynomials[0].at(offset) = row[0];
            table_polynomials[1].at(offset) = row[1];
            table_polynomials[2].at(offset) = row[2];
            table_polynomials[3].at(offset) = FF(table.table_index);
            offset++;
        }
    }
    BB_ASSERT(offset <= table_polynomials[0].end_index(),
              "construct_lookup_table_polynomials: total lookup table entries exceed polynomial size");
}

template class TraceToPolynomials<zcash::UltraPastaZKFlavor>;
template class ProverInstance_<zcash::UltraPastaZKFlavor>;
template class OinkProver<zcash::UltraPastaZKFlavor>;
template class OinkVerifier<zcash::UltraPastaZKFlavor>;
template class SmallSubgroupIPAProver<zcash::UltraPastaZKFlavor>;
} // namespace bb

namespace bb::zcash {

using Flavor = UltraPastaZKFlavor;
using FF = Flavor::FF;
using Curve = Flavor::Curve;
using Commitment = Flavor::Commitment;

namespace {
// SmallSubgroupIPA commits to polynomials of size up to SUBGROUP_SIZE + 3.
size_t commitment_key_size(size_t n)
{
    return std::max(n, size_t{ 2 } * Curve::SUBGROUP_SIZE);
}
} // namespace

Flavor::Proof ultra_pasta_prove(const std::shared_ptr<UltraPastaProverInstance>& instance,
                                const std::shared_ptr<Flavor::VerificationKey>& vk)
{
    BB_BENCH_NAME("ultra_pasta_prove");
    const size_t n = instance->dyadic_size();
    const size_t log_n = instance->log_dyadic_size();
    auto transcript = std::make_shared<Flavor::Transcript>();
    Flavor::CommitmentKey ck(commitment_key_size(n));

    OinkProver<Flavor> oink(instance, vk, transcript);
    oink.prove();

    const auto gate_challenges = transcript->get_dyadic_powers_of_challenge<FF>("Sumcheck:gate_challenge", log_n);
    SumcheckProver<Flavor> sumcheck(
        n, instance->polynomials, transcript, instance->alpha, gate_challenges, instance->relation_parameters, log_n);
    ZKSumcheckData<Flavor> zk_sumcheck_data(log_n, transcript, ck);
    auto sumcheck_output = sumcheck.prove(zk_sumcheck_data);

    SmallSubgroupIPAProver<Flavor> small_subgroup_ipa(
        zk_sumcheck_data, sumcheck_output.challenge, sumcheck_output.claimed_libra_evaluation, transcript, ck);
    small_subgroup_ipa.prove();

    using PolynomialBatcher = GeminiProver_<Curve>::PolynomialBatcher;
    PolynomialBatcher batcher(n, instance->polynomials.max_end_index());
    batcher.set_unshifted(instance->polynomials.get_unshifted());
    batcher.set_to_be_shifted_by_one(instance->polynomials.get_to_be_shifted());
    auto opening_claim = ShpleminiProver_<Curve>::prove(
        n, batcher, sumcheck_output.challenge, ck, transcript, small_subgroup_ipa.get_witness_polynomials());
    {
        BB_BENCH_NAME("ultra_pasta_prove/halo2_ipa");
        Halo2IPA<Curve>::prove(halo2_vesta_ipa_generators(ck, n), opening_claim, FF(0), transcript);
    }
    return transcript->export_proof();
}

bool ultra_pasta_verify(const std::shared_ptr<Flavor::VerificationKey>& vk,
                        const Flavor::Proof& proof,
                        std::vector<pasta::fp>* public_inputs)
{
    BB_BENCH_NAME("ultra_pasta_verify");
    const size_t log_n = vk->log_circuit_size;
    const size_t n = size_t{ 1 } << log_n;
    auto transcript = std::make_shared<Flavor::Transcript>();
    transcript->load_proof(proof);
    auto instance = std::make_shared<VerifierInstance_<Flavor>>(std::make_shared<Flavor::VKAndHash>(vk));

    OinkVerifier<Flavor> oink(instance, transcript, vk->num_public_inputs);
    oink.verify();
    instance->gate_challenges = transcript->get_dyadic_powers_of_challenge<FF>("Sumcheck:gate_challenge", log_n);
    auto commitments = VerifierCommitmentsConstructor<Flavor>::construct(
        instance->get_vk(), instance->witness_commitments, instance->gemini_masking_commitment);

    SumcheckVerifier<Flavor> sumcheck(transcript, instance->alpha, log_n);
    std::array<Commitment, NUM_SMALL_IPA_COMMITMENTS> libra_commitments = {};
    libra_commitments[0] = transcript->receive_from_prover<Commitment>("Libra:concatenation_commitment");
    auto sumcheck_output = sumcheck.verify(instance->relation_parameters, instance->gate_challenges);
    libra_commitments[1] = transcript->receive_from_prover<Commitment>("Libra:grand_sum_commitment");
    libra_commitments[2] = transcript->receive_from_prover<Commitment>("Libra:quotient_commitment");

    using ClaimBatcher = ClaimBatcher_<Curve>;
    ClaimBatcher claim_batcher{
        .unshifted =
            ClaimBatcher::Batch{ commitments.get_unshifted(), sumcheck_output.claimed_evaluations.get_unshifted() },
        .shifted =
            ClaimBatcher::Batch{ commitments.get_to_be_shifted(), sumcheck_output.claimed_evaluations.get_shifted() },
    };
    Flavor::CommitmentKey ck(commitment_key_size(n));
    const Commitment g1_identity = ck.get_monomial_points()[0];
    auto shplemini_output =
        ShpleminiVerifier_<Curve, true, true>::compute_batch_opening_claim(claim_batcher,
                                                                           sumcheck_output.challenge,
                                                                           g1_identity,
                                                                           transcript,
                                                                           Flavor::REPEATED_COMMITMENTS,
                                                                           libra_commitments,
                                                                           sumcheck_output.claimed_libra_evaluation);
    const auto& batch_claim = shplemini_output.batch_opening_claim;
    std::vector<FF> scalars;
    std::vector<Commitment> points;
    for (size_t i = 0; i < batch_claim.commitments.size(); ++i) {
        if (!batch_claim.commitments[i].is_point_at_infinity()) {
            scalars.push_back(batch_claim.scalars[i]);
            points.push_back(batch_claim.commitments[i]);
        }
    }
    const OpeningClaim<Curve> opening_claim{ { batch_claim.evaluation_point, FF(0) },
                                             Commitment(Halo2IPA<Curve>::msm(scalars, points)) };
    const bool pcs_verified = Halo2IPA<Curve>::verify(halo2_vesta_ipa_generators(ck, n), opening_claim, transcript);
    vinfo("ultra pasta verifier: sumcheck ",
          sumcheck_output.verified,
          ", libra consistency ",
          shplemini_output.consistency_checked,
          ", pcs ",
          pcs_verified);
    if (public_inputs != nullptr) {
        *public_inputs = instance->public_inputs;
    }
    return sumcheck_output.verified && shplemini_output.consistency_checked && pcs_verified;
}

} // namespace bb::zcash
