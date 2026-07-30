#pragma once

#include "barretenberg/commitment_schemes/whir/whir.hpp"
#include "barretenberg/flavor/ultra_flavor.hpp"
#include "barretenberg/honk/library/grand_product_delta.hpp"
#include "barretenberg/sumcheck/sumcheck.hpp"
#include "barretenberg/transcript/transcript.hpp"
#include "barretenberg/ultra_honk/oink_prover.hpp"
#include "barretenberg/ultra_honk/prover_instance.hpp"

#include <memory>
#include <vector>

namespace bb::whir {

/**
 * @brief UltraHonk with WHIR as the polynomial commitment scheme.
 *
 * @details A transparent (no SRS, no pairings) UltraHonk: the arithmetization, trace layout,
 * relations, and sumcheck are exactly `UltraFlavor`'s; every elliptic-curve commitment is replaced
 * by a WHIR Merkle root, and the Gemini+Shplonk+KZG opening phase by one batched WHIR opening
 * (README.md §10). The Fiat-Shamir schedule mirrors Oink round by round — commitments are absorbed
 * before the challenges that depend on them, with the identical challenge labels — so the
 * relation-parameter derivation and the derived-polynomial computations (`OinkProver`'s static
 * helpers) are shared with the production prover. Sumcheck runs unpadded (`virtual_log_n = log_n`).
 *
 * `prove` consumes the proving key's instance (memory records are appended to w_4, derived
 * polynomials are computed in place): create a fresh proving key per proof.
 */
template <typename Hasher = Poseidon2MerkleHasher> class WhirHonk {
  public:
    using Flavor = UltraFlavor;
    using FF = Flavor::FF;
    using ProverInstance = ProverInstance_<Flavor>;
    using Oink = OinkProver<Flavor>;
    using Builder = UltraCircuitBuilder;
    using Digest = typename Hasher::Digest;
    using Transcript = NativeTranscript;

    // get_witness() entity order; the per-round commit schedule below covers exactly these.
    static constexpr size_t NUM_WITNESS = Flavor::NUM_WITNESS_ENTITIES;
    static_assert(Flavor::NUM_PRECOMPUTED_ENTITIES == 28 && NUM_WITNESS == 8 && Flavor::NUM_SHIFTED_ENTITIES == 5,
                  "Ultra entity layout changed; revisit the WhirHonk commit schedule and claim assembly");

    /** @brief Transparent verification key: circuit metadata plus the precomputed-polynomial roots. */
    struct VerificationKey {
        size_t log_dyadic_size = 0;
        size_t num_public_inputs = 0;
        size_t pub_inputs_offset = 0;
        std::vector<Digest> precomputed_roots; // EntityId order (= get_precomputed())
    };

    struct ProvingKey {
        std::shared_ptr<ProverInstance> instance;
        std::shared_ptr<WhirCommitmentKey<Hasher>> ck;
        std::vector<WhirProverData<Hasher>> precomputed; // EntityId order
        VerificationKey vk;
    };

    static WhirConfig make_config(size_t log_dyadic_size,
                                  size_t security_bits = 100,
                                  size_t log_inv_rate = 2,
                                  size_t folding_factor_bits = 4,
                                  size_t final_poly_bits = 4)
    {
        return WhirConfig::create(log_dyadic_size, security_bits, log_inv_rate, folding_factor_bits, final_poly_bits);
    }

    static ProvingKey create_proving_key(Builder& circuit, const WhirConfig& config)
    {
        ProvingKey pk;
        pk.instance = std::make_shared<ProverInstance>(circuit);
        BB_ASSERT_EQ(pk.instance->log_dyadic_size(), config.num_payload_variables, "config/circuit size mismatch");
        pk.ck = std::make_shared<WhirCommitmentKey<Hasher>>(config);
        for (const auto& polynomial : pk.instance->polynomials.get_precomputed()) {
            pk.precomputed.push_back(pk.ck->commit(polynomial));
        }
        pk.vk = VerificationKey{ .log_dyadic_size = pk.instance->log_dyadic_size(),
                                 .num_public_inputs = pk.instance->num_public_inputs(),
                                 .pub_inputs_offset = pk.instance->pub_inputs_offset(),
                                 .precomputed_roots = {} };
        for (const auto& data : pk.precomputed) {
            pk.vk.precomputed_roots.push_back(data.tree.root());
        }
        return pk;
    }

    static HonkProof prove(ProvingKey& pk)
    {
        auto transcript = Transcript::test_prover_init_empty();
        ProverInstance& instance = *pk.instance;
        const WhirCommitmentKey<Hasher>& ck = *pk.ck;
        const size_t log_n = instance.log_dyadic_size();

        absorb_vk(transcript, pk.vk);
        for (size_t i = 0; i < instance.num_public_inputs(); ++i) {
            transcript->send_to_verifier("public_input_" + std::to_string(i), instance.public_inputs[i]);
        }

        // Oink round 1: wire commitments (w_4 deferred until memory records are added).
        auto w_l = commit_and_send(ck, transcript, instance.polynomials.w_l(), "w_l", true);
        auto w_r = commit_and_send(ck, transcript, instance.polynomials.w_r(), "w_r", true);
        auto w_o = commit_and_send(ck, transcript, instance.polynomials.w_o(), "w_o", true);

        // Oink round 2: memory-record challenges, lookup counts, and the completed fourth wire.
        auto [eta, rom_logup_gamma] =
            transcript->template get_challenges<FF>(std::array<std::string, 2>{ "eta", "rom_logup_gamma" });
        instance.relation_parameters.eta = eta;
        instance.relation_parameters.eta_two = eta * eta;
        instance.relation_parameters.eta_three = instance.relation_parameters.eta_two * eta;
        instance.relation_parameters.rom_logup_gamma = rom_logup_gamma;
        Oink::add_ram_rom_memory_records_to_wire_4(instance);
        Oink::add_rom_logup_inverses_to_wire_4(instance);
        auto lookup_read_counts =
            commit_and_send(ck, transcript, instance.polynomials.lookup_read_counts(), "lookup_read_counts", false);
        auto lookup_read_tags =
            commit_and_send(ck, transcript, instance.polynomials.lookup_read_tags(), "lookup_read_tags", false);
        auto w_4 = commit_and_send(ck, transcript, instance.polynomials.w_4(), "w_4", true);

        // Oink round 3: log-derivative lookup inverses.
        auto [beta, gamma] = transcript->template get_challenges<FF>(std::array<std::string, 2>{ "beta", "gamma" });
        instance.relation_parameters.beta = beta;
        instance.relation_parameters.beta_sqr = beta * beta;
        instance.relation_parameters.beta_cube = instance.relation_parameters.beta_sqr * beta;
        instance.relation_parameters.gamma = gamma;
        Oink::compute_logderivative_inverses(instance);
        auto lookup_inverses =
            commit_and_send(ck, transcript, instance.polynomials.lookup_inverses(), "lookup_inverses", false);

        // Oink round 4: permutation grand product (also sets public_input_delta).
        uint32_t z_perm_dup_count = 0;
        Oink::compute_grand_product_polynomial(instance, z_perm_dup_count);
        auto z_perm = commit_and_send(ck, transcript, instance.polynomials.z_perm(), "z_perm", true);

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

        // WHIR opening of every unshifted and to-be-shifted claim at the sumcheck challenge point.
        typename WhirProver<Hasher>::Claims claims;
        claims.send_roots = false; // all roots are already bound to the transcript above
        for (const auto& data : pk.precomputed) {
            claims.unshifted.push_back(&data);
        }
        for (const auto* data :
             { &w_l, &w_r, &w_o, &w_4, &z_perm, &lookup_inverses, &lookup_read_counts, &lookup_read_tags }) {
            claims.unshifted.push_back(data);
        }
        const auto unshifted_evaluations = sumcheck_output.claimed_evaluations.get_unshifted();
        claims.unshifted_evaluations.assign(unshifted_evaluations.begin(), unshifted_evaluations.end());
        for (const auto* data : { &w_l, &w_r, &w_o, &w_4, &z_perm }) {
            claims.to_be_shifted.push_back(data);
        }
        const auto shifted_evaluations = sumcheck_output.claimed_evaluations.get_shifted();
        claims.shifted_evaluations.assign(shifted_evaluations.begin(), shifted_evaluations.end());

        WhirProver<Hasher>::prove(ck, claims, sumcheck_output.challenge, transcript);
        return transcript->export_proof();
    }

    static bool verify(const VerificationKey& vk, const WhirConfig& config, const HonkProof& proof)
    {
        auto transcript = std::make_shared<Transcript>(proof);
        [[maybe_unused]] auto init = transcript->template receive_from_prover<FF>("Init");
        const size_t log_n = vk.log_dyadic_size;

        absorb_vk(transcript, vk);
        std::vector<FF> public_inputs(vk.num_public_inputs);
        for (size_t i = 0; i < vk.num_public_inputs; ++i) {
            public_inputs[i] = transcript->template receive_from_prover<FF>("public_input_" + std::to_string(i));
        }

        const Digest w_l = receive_root(transcript, "w_l");
        const Digest w_r = receive_root(transcript, "w_r");
        const Digest w_o = receive_root(transcript, "w_o");

        auto [eta, rom_logup_gamma] =
            transcript->template get_challenges<FF>(std::array<std::string, 2>{ "eta", "rom_logup_gamma" });
        RelationParameters<FF> relation_parameters;
        relation_parameters.eta = eta;
        relation_parameters.eta_two = eta * eta;
        relation_parameters.eta_three = relation_parameters.eta_two * eta;
        relation_parameters.rom_logup_gamma = rom_logup_gamma;
        const Digest lookup_read_counts = receive_root(transcript, "lookup_read_counts");
        const Digest lookup_read_tags = receive_root(transcript, "lookup_read_tags");
        const Digest w_4 = receive_root(transcript, "w_4");

        auto [beta, gamma] = transcript->template get_challenges<FF>(std::array<std::string, 2>{ "beta", "gamma" });
        relation_parameters.beta = beta;
        relation_parameters.beta_sqr = beta * beta;
        relation_parameters.beta_cube = relation_parameters.beta_sqr * beta;
        relation_parameters.gamma = gamma;
        relation_parameters.public_input_delta =
            compute_public_input_delta<Flavor>(public_inputs, beta, gamma, FF(vk.pub_inputs_offset));
        const Digest lookup_inverses = receive_root(transcript, "lookup_inverses");
        const Digest z_perm = receive_root(transcript, "z_perm");

        const FF alpha = transcript->template get_challenge<FF>("alpha");
        const std::vector<FF> gate_challenges =
            transcript->template get_dyadic_powers_of_challenge<FF>("Sumcheck:gate_challenge", log_n);

        SumcheckVerifier<Flavor> sumcheck(transcript, alpha, log_n);
        SumcheckOutput<Flavor> sumcheck_output = sumcheck.verify(relation_parameters, gate_challenges);
        if (!sumcheck_output.verified) {
            return false;
        }

        typename WhirVerifier<Hasher>::Claims claims;
        const auto unshifted_evaluations = sumcheck_output.claimed_evaluations.get_unshifted();
        claims.unshifted_evaluations.assign(unshifted_evaluations.begin(), unshifted_evaluations.end());
        const auto shifted_evaluations = sumcheck_output.claimed_evaluations.get_shifted();
        claims.shifted_evaluations.assign(shifted_evaluations.begin(), shifted_evaluations.end());
        claims.unshifted_roots = vk.precomputed_roots;
        for (const Digest& root :
             { w_l, w_r, w_o, w_4, z_perm, lookup_inverses, lookup_read_counts, lookup_read_tags }) {
            claims.unshifted_roots.push_back(root);
        }
        claims.shifted_roots = { w_l, w_r, w_o, w_4, z_perm };

        return WhirVerifier<Hasher>::verify(config, claims, sumcheck_output.challenge, transcript);
    }

  private:
    static void absorb_vk(const std::shared_ptr<Transcript>& transcript, const VerificationKey& vk)
    {
        transcript->add_to_hash_buffer("vk_log_dyadic_size", FF(vk.log_dyadic_size));
        transcript->add_to_hash_buffer("vk_num_public_inputs", FF(vk.num_public_inputs));
        transcript->add_to_hash_buffer("vk_pub_inputs_offset", FF(vk.pub_inputs_offset));
        for (const Digest& root : vk.precomputed_roots) {
            transcript->add_to_hash_buffer("vk_root", Hasher::digest_to_fields(root));
        }
    }

    static WhirProverData<Hasher> commit_and_send(const WhirCommitmentKey<Hasher>& ck,
                                                  const std::shared_ptr<Transcript>& transcript,
                                                  const Polynomial<fr>& polynomial,
                                                  const std::string& label,
                                                  bool to_be_shifted)
    {
        WhirProverData<Hasher> data = ck.commit(polynomial, to_be_shifted);
        transcript->send_to_verifier("WHIRHONK:" + label, Hasher::digest_to_fields(data.tree.root()));
        return data;
    }

    template <typename TranscriptPtr>
    static Digest receive_root(const TranscriptPtr& transcript, const std::string& label)
    {
        const auto fields =
            transcript->template receive_from_prover<std::array<fr, Hasher::DIGEST_NUM_FIELDS>>("WHIRHONK:" + label);
        return Hasher::digest_from_fields(fields);
    }
};

} // namespace bb::whir
