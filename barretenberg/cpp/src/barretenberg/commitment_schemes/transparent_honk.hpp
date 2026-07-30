#pragma once

#include "barretenberg/flavor/ultra_flavor.hpp"
#include "barretenberg/honk/library/grand_product_delta.hpp"
#include "barretenberg/sumcheck/sumcheck.hpp"
#include "barretenberg/transcript/transcript.hpp"
#include "barretenberg/ultra_honk/oink_prover.hpp"
#include "barretenberg/ultra_honk/prover_instance.hpp"

#include <memory>
#include <vector>

namespace bb::honk_transparent {

/**
 * @brief UltraHonk with a transparent (hash-based) polynomial commitment scheme.
 *
 * @details The arithmetization, trace layout, relations, and sumcheck are exactly `UltraFlavor`'s;
 * every elliptic-curve commitment is replaced by a Merkle root of the `Pcs` backend (WHIR or
 * Ligero), and the Gemini+Shplonk+KZG opening phase by one batched opening of that backend. The
 * Fiat-Shamir schedule mirrors Oink round by round — commitments are absorbed before the challenges
 * that depend on them, with the identical challenge labels — so the relation-parameter derivation
 * and the derived-polynomial computations (`OinkProver`'s static helpers) are shared with the
 * production prover. Sumcheck runs unpadded (`virtual_log_n = log_n`).
 *
 * Columns committed in the same round share one Merkle tree. The five groups:
 *
 * | group | columns (entity order within) | to-be-shifted |
 * |---|---|---|
 * | 0 (verification key) | the 28 precomputed polynomials | - |
 * | 1 `HONK:wires` | w_l, w_r, w_o | all |
 * | 2 `HONK:counts_w4` | lookup_read_counts, lookup_read_tags, w_4 | w_4 |
 * | 3 `HONK:lookup_inverses` | lookup_inverses | - |
 * | 4 `HONK:z_perm` | z_perm | z_perm |
 *
 * The `Pcs` backend supplies: `Config`, `CommitmentKey` (constructible from `Config`, with
 * `commit_group(span of Polynomial*, to_be_shifted flags)`), `GroupData`, `GroupCommitment` (a
 * Merkle digest for hash backends; a vector of curve points for KZG-based backends) with the
 * send/receive/absorb hooks, `ProverClaims`/`VerifierClaims` (the shared claim shape:
 * groups/commitments, `{group, column}` refs, evaluations), `make_config`, `payload_variables`,
 * `prove_opening`, and `verify_opening`.
 *
 * `prove` consumes the proving key's instance (memory records are appended to w_4, derived
 * polynomials are computed in place): create a fresh proving key per proof.
 */
template <typename Pcs> class TransparentHonk {
  public:
    using Flavor = UltraFlavor;
    using FF = Flavor::FF;
    using ProverInstance = ProverInstance_<Flavor>;
    using Oink = OinkProver<Flavor>;
    using Builder = UltraCircuitBuilder;
    using Config = typename Pcs::Config;
    using CommitmentKey = typename Pcs::CommitmentKey;
    using GroupData = typename Pcs::GroupData;
    using GroupCommitment = typename Pcs::GroupCommitment;
    using Transcript = NativeTranscript;

    static constexpr size_t NUM_WITNESS = Flavor::NUM_WITNESS_ENTITIES;
    static_assert(Flavor::NUM_PRECOMPUTED_ENTITIES == 28 && NUM_WITNESS == 8 && Flavor::NUM_SHIFTED_ENTITIES == 5,
                  "Ultra entity layout changed; revisit the commit schedule and claim assembly");

    // Column counts of the five commitment groups, in group order.
    static constexpr std::array<size_t, 5> GROUP_COLUMNS = { 28, 3, 3, 1, 1 };

    /** @brief Transparent verification key: circuit metadata plus the precomputed-group root. */
    struct VerificationKey {
        size_t log_dyadic_size = 0;
        size_t num_public_inputs = 0;
        size_t pub_inputs_offset = 0;
        GroupCommitment precomputed_commitment = {};
    };

    struct ProvingKey {
        std::shared_ptr<ProverInstance> instance;
        std::shared_ptr<CommitmentKey> ck;
        std::shared_ptr<GroupData> precomputed; // group 0
        VerificationKey vk;
    };

    static Config make_config(size_t log_dyadic_size, size_t security_bits = 100, size_t log_inv_rate = 2)
    {
        return Pcs::make_config(log_dyadic_size, security_bits, log_inv_rate);
    }

    static ProvingKey create_proving_key(Builder& circuit, const Config& config)
    {
        ProvingKey pk;
        pk.instance = std::make_shared<ProverInstance>(circuit);
        BB_ASSERT_EQ(pk.instance->log_dyadic_size(), Pcs::payload_variables(config), "config/circuit size mismatch");
        pk.ck = std::make_shared<CommitmentKey>(config);
        std::vector<const Polynomial<fr>*> precomputed_columns;
        for (const auto& polynomial : pk.instance->polynomials.get_precomputed()) {
            precomputed_columns.push_back(&polynomial);
        }
        pk.precomputed = std::make_shared<GroupData>(pk.ck->commit_group(precomputed_columns));
        pk.vk = VerificationKey{ .log_dyadic_size = pk.instance->log_dyadic_size(),
                                 .num_public_inputs = pk.instance->num_public_inputs(),
                                 .pub_inputs_offset = pk.instance->pub_inputs_offset(),
                                 .precomputed_commitment = Pcs::group_commitment(*pk.precomputed) };
        return pk;
    }

    static HonkProof prove(ProvingKey& pk)
    {
        auto transcript = Transcript::test_prover_init_empty();
        ProverInstance& instance = *pk.instance;
        const CommitmentKey& ck = *pk.ck;
        const size_t log_n = instance.log_dyadic_size();

        absorb_vk(transcript, pk.vk);
        for (size_t i = 0; i < instance.num_public_inputs(); ++i) {
            transcript->send_to_verifier("public_input_" + std::to_string(i), instance.public_inputs[i]);
        }

        // Oink round 1: wire commitments (w_4 deferred until memory records are added).
        auto wires =
            commit_and_send(ck,
                            transcript,
                            { &instance.polynomials.w_l(), &instance.polynomials.w_r(), &instance.polynomials.w_o() },
                            { true, true, true },
                            "wires");

        // Oink round 2: memory-record challenges, lookup counts, and the completed fourth wire.
        auto [eta, rom_logup_gamma] =
            transcript->template get_challenges<FF>(std::array<std::string, 2>{ "eta", "rom_logup_gamma" });
        instance.relation_parameters.eta = eta;
        instance.relation_parameters.eta_two = eta * eta;
        instance.relation_parameters.eta_three = instance.relation_parameters.eta_two * eta;
        instance.relation_parameters.rom_logup_gamma = rom_logup_gamma;
        Oink::add_ram_rom_memory_records_to_wire_4(instance);
        Oink::add_rom_logup_inverses_to_wire_4(instance);
        auto counts_w4 = commit_and_send(ck,
                                         transcript,
                                         { &instance.polynomials.lookup_read_counts(),
                                           &instance.polynomials.lookup_read_tags(),
                                           &instance.polynomials.w_4() },
                                         { false, false, true },
                                         "counts_w4");

        // Oink round 3: log-derivative lookup inverses.
        auto [beta, gamma] = transcript->template get_challenges<FF>(std::array<std::string, 2>{ "beta", "gamma" });
        instance.relation_parameters.beta = beta;
        instance.relation_parameters.beta_sqr = beta * beta;
        instance.relation_parameters.beta_cube = instance.relation_parameters.beta_sqr * beta;
        instance.relation_parameters.gamma = gamma;
        Oink::compute_logderivative_inverses(instance);
        auto lookup_inverses =
            commit_and_send(ck, transcript, { &instance.polynomials.lookup_inverses() }, { false }, "lookup_inverses");

        // Oink round 4: permutation grand product (also sets public_input_delta).
        uint32_t z_perm_dup_count = 0;
        Oink::compute_grand_product_polynomial(instance, z_perm_dup_count);
        auto z_perm = commit_and_send(ck, transcript, { &instance.polynomials.z_perm() }, { true }, "z_perm");

        instance.alpha = transcript->template get_challenge<FF>("alpha");
        instance.gate_challenges =
            transcript->template get_dyadic_powers_of_challenge<FF>("Sumcheck:gate_challenge", log_n);

        SumcheckProver<Flavor> sumcheck(instance.dyadic_size(),
                                        instance.polynomials,
                                        transcript,
                                        instance.alpha,
                                        instance.gate_challenges,
                                        instance.relation_parameters,
                                        /*virtual_log_n=*/log_n);
        SumcheckOutput<Flavor> sumcheck_output = sumcheck.prove();

        // Batched opening of every unshifted and to-be-shifted claim at the sumcheck challenge.
        typename Pcs::ProverClaims claims;
        claims.send_roots = false; // all roots are already bound to the transcript above
        claims.groups = { &*pk.precomputed, &wires, &counts_w4, &lookup_inverses, &z_perm };
        append_unshifted_refs(claims.unshifted);
        const auto unshifted_evaluations = sumcheck_output.claimed_evaluations.get_unshifted();
        claims.unshifted_evaluations.assign(unshifted_evaluations.begin(), unshifted_evaluations.end());
        append_shifted_refs(claims.to_be_shifted);
        const auto shifted_evaluations = sumcheck_output.claimed_evaluations.get_shifted();
        claims.shifted_evaluations.assign(shifted_evaluations.begin(), shifted_evaluations.end());

        Pcs::prove_opening(ck, claims, sumcheck_output.challenge, transcript);
        return transcript->export_proof();
    }

    static bool verify(const VerificationKey& vk, const Config& config, const HonkProof& proof)
    {
        auto transcript = std::make_shared<Transcript>(proof);
        [[maybe_unused]] auto init = transcript->template receive_from_prover<FF>("Init");
        const size_t log_n = vk.log_dyadic_size;

        absorb_vk(transcript, vk);
        std::vector<FF> public_inputs(vk.num_public_inputs);
        for (size_t i = 0; i < vk.num_public_inputs; ++i) {
            public_inputs[i] = transcript->template receive_from_prover<FF>("public_input_" + std::to_string(i));
        }

        const GroupCommitment wires = Pcs::receive_group_commitment(transcript, "HONK:wires", GROUP_COLUMNS[1]);

        auto [eta, rom_logup_gamma] =
            transcript->template get_challenges<FF>(std::array<std::string, 2>{ "eta", "rom_logup_gamma" });
        RelationParameters<FF> relation_parameters;
        relation_parameters.eta = eta;
        relation_parameters.eta_two = eta * eta;
        relation_parameters.eta_three = relation_parameters.eta_two * eta;
        relation_parameters.rom_logup_gamma = rom_logup_gamma;
        const GroupCommitment counts_w4 =
            Pcs::receive_group_commitment(transcript, "HONK:counts_w4", GROUP_COLUMNS[2]);

        auto [beta, gamma] = transcript->template get_challenges<FF>(std::array<std::string, 2>{ "beta", "gamma" });
        relation_parameters.beta = beta;
        relation_parameters.beta_sqr = beta * beta;
        relation_parameters.beta_cube = relation_parameters.beta_sqr * beta;
        relation_parameters.gamma = gamma;
        relation_parameters.public_input_delta =
            compute_public_input_delta<Flavor>(public_inputs, beta, gamma, FF(vk.pub_inputs_offset));
        const GroupCommitment lookup_inverses =
            Pcs::receive_group_commitment(transcript, "HONK:lookup_inverses", GROUP_COLUMNS[3]);
        const GroupCommitment z_perm = Pcs::receive_group_commitment(transcript, "HONK:z_perm", GROUP_COLUMNS[4]);

        const FF alpha = transcript->template get_challenge<FF>("alpha");
        const std::vector<FF> gate_challenges =
            transcript->template get_dyadic_powers_of_challenge<FF>("Sumcheck:gate_challenge", log_n);

        SumcheckVerifier<Flavor> sumcheck(transcript, alpha, log_n);
        SumcheckOutput<Flavor> sumcheck_output = sumcheck.verify(relation_parameters, gate_challenges);
        if (!sumcheck_output.verified) {
            return false;
        }

        typename Pcs::VerifierClaims claims;
        claims.group_num_columns.assign(GROUP_COLUMNS.begin(), GROUP_COLUMNS.end());
        Pcs::set_group_commitments(claims,
                                   { vk.precomputed_commitment, wires, counts_w4, lookup_inverses, z_perm });
        append_unshifted_refs(claims.unshifted);
        const auto unshifted_evaluations = sumcheck_output.claimed_evaluations.get_unshifted();
        claims.unshifted_evaluations.assign(unshifted_evaluations.begin(), unshifted_evaluations.end());
        append_shifted_refs(claims.to_be_shifted);
        const auto shifted_evaluations = sumcheck_output.claimed_evaluations.get_shifted();
        claims.shifted_evaluations.assign(shifted_evaluations.begin(), shifted_evaluations.end());

        return Pcs::verify_opening(config, claims, sumcheck_output.challenge, transcript);
    }

  private:
    /**
     * @brief (group, column) of every unshifted entity, in `get_unshifted()` entity order:
     * the 28 precomputed columns, then w_l, w_r, w_o, w_4, z_perm, lookup_inverses,
     * lookup_read_counts, lookup_read_tags.
     */
    template <typename RefVector> static void append_unshifted_refs(RefVector& refs)
    {
        for (size_t c = 0; c < GROUP_COLUMNS[0]; ++c) {
            refs.push_back({ 0, c });
        }
        refs.push_back({ 1, 0 }); // w_l
        refs.push_back({ 1, 1 }); // w_r
        refs.push_back({ 1, 2 }); // w_o
        refs.push_back({ 2, 2 }); // w_4
        refs.push_back({ 4, 0 }); // z_perm
        refs.push_back({ 3, 0 }); // lookup_inverses
        refs.push_back({ 2, 0 }); // lookup_read_counts
        refs.push_back({ 2, 1 }); // lookup_read_tags
    }

    /** @brief (group, column) of the to-be-shifted entities in `get_shifted()` entity order. */
    template <typename RefVector> static void append_shifted_refs(RefVector& refs)
    {
        refs.push_back({ 1, 0 }); // w_l
        refs.push_back({ 1, 1 }); // w_r
        refs.push_back({ 1, 2 }); // w_o
        refs.push_back({ 2, 2 }); // w_4
        refs.push_back({ 4, 0 }); // z_perm
    }

    static void absorb_vk(const std::shared_ptr<Transcript>& transcript, const VerificationKey& vk)
    {
        transcript->add_to_hash_buffer("vk_log_dyadic_size", FF(vk.log_dyadic_size));
        transcript->add_to_hash_buffer("vk_num_public_inputs", FF(vk.num_public_inputs));
        transcript->add_to_hash_buffer("vk_pub_inputs_offset", FF(vk.pub_inputs_offset));
        Pcs::absorb_group_commitment(transcript, "vk_root", vk.precomputed_commitment);
    }

    static GroupData commit_and_send(const CommitmentKey& ck,
                                     const std::shared_ptr<Transcript>& transcript,
                                     const std::vector<const Polynomial<fr>*>& columns,
                                     const std::vector<bool>& to_be_shifted,
                                     const std::string& label)
    {
        GroupData data = ck.commit_group(columns, to_be_shifted);
        Pcs::send_group_commitment(transcript, "HONK:" + label, data);
        return data;
    }
};

} // namespace bb::honk_transparent
