#pragma once

#include "barretenberg/commitment_schemes/transparent_honk.hpp"
#include "barretenberg/commitment_schemes/whir/stdlib/whir_recursive_verifier.hpp"
#include "barretenberg/commitment_schemes/whir/whir_honk.hpp"
#include "barretenberg/flavor/ultra_provekit_recursive_flavor.hpp"
#include "barretenberg/honk/library/grand_product_delta.hpp"
#include "barretenberg/sumcheck/sumcheck.hpp"

#include <span>
#include <vector>

namespace bb::whir::recursion {

/**
 * @brief The WHIR instantiation a recursive verifier wants, as `whir/stdlib/README.md` §4 measures.
 * @details `InitialFoldingFactorBits = 1` is the one that matters: round 0 is the only round that
 * touches the wide commitment, and its leaf holds `2^k0` values of *every* Honk column, so the
 * default `k0 = k = 4` makes each query hash sixteen values per column instead of two. Johnson
 * soundness (`PROVABLE_LIST`) also roughly doubles the query count against the repaired-list bound
 * for the same lambda. Grinding supplies 20 bits so the queries do not have to.
 */
using RecursionWhirPcs = WhirPcs<Poseidon2CompressionHasher,
                                 /*MaxStackBits=*/0,
                                 WhirSoundness::REPAIRED_LIST,
                                 /*FoldingFactorBits=*/4,
                                 /*InitialFoldingFactorBits=*/1,
                                 /*PowBits=*/20>;

/**
 * @brief Recursive verifier for a whole transparent-Honk-with-WHIR proof.
 *
 * @details The in-circuit mirror of `TransparentHonk<WhirPcs, Flavor>::verify`, end to end: the
 * verification-key absorption, the public inputs, the three witness commitment rounds and their
 * challenges, the full sumcheck over the flavor's relation set, the virtual-column checks, and the
 * WHIR opening of every claimed evaluation. Verifying its circuit is verifying the inner proof.
 *
 * Three things make this cheap enough to be worth doing.
 *
 * *No curve arithmetic anywhere.* A hash-based commitment is a Merkle root — a field element — so
 * the commitments, the transcript and the opening are all native BN254 `Fr`. A KZG-based Honk
 * verifier has to emulate `G1` operations or defer them to a Goblin/ECCVM tail; this one has
 * nothing to defer.
 *
 * *The relation set is bb's own.* `Flavor::Relations` is re-instantiated over `stdlib::field_t`, so
 * `SumcheckVerifier` evaluates exactly the subrelations the native verifier does, from the same
 * source. There is no hand-written copy of the constraint system to drift.
 *
 * *The reduced flavor pays twice.* `UltraProveKitFlavor` drops the elliptic, non-native-field and
 * Poseidon2 relations, which removes four committed columns *and* shortens every sumcheck round
 * univariate from 7 evaluations to 6 — both of which the recursive verifier feels directly.
 */
template <typename Builder, typename Pcs_ = RecursionWhirPcs> class TransparentHonkRecursiveVerifier {
  public:
    using Pcs = Pcs_;
    using NativeFlavor = UltraProveKitFlavor;
    using Flavor = UltraProveKitRecursiveFlavor_<Builder>;
    using FF = stdlib::field_t<Builder>;
    using Transcript = StdlibTranscript<Builder>;
    using NativeHonk = honk_transparent::TransparentHonk<Pcs, NativeFlavor>;
    using NativeVerificationKey = typename NativeHonk::VerificationKey;
    using WhirVerifier = WhirRecursiveVerifier<Builder>;

    static constexpr size_t NUM_PRECOMPUTED = NativeFlavor::NUM_PRECOMPUTED_ENTITIES;

    /**
     * @brief Verify `proof` against `vk` inside `builder`; returns the inner proof's public inputs.
     * @details `vk` is the native verification key: everything in it is public data fixed when the
     * aggregator is built, so its fields enter the circuit as constants. That is also what binds
     * the statement — the precomputed-column root is a constant, so a proof of a different circuit
     * cannot satisfy this one.
     *
     * **The returned public inputs are not bound to anything yet, and binding them is the caller's
     * job.** Left alone, this circuit proves "some proof of this circuit exists", with the inner
     * public inputs free for the prover to choose. An aggregator that means "a proof of *this
     * statement* exists" must forward them — as public inputs of the outer circuit, or into whatever
     * it hashes — which is why they are returned rather than dropped.
     */
    [[nodiscard]] static std::vector<FF> verify(Builder& builder,
                                                const NativeVerificationKey& vk,
                                                const WhirConfig& config,
                                                const typename Transcript::Proof& proof,
                                                GateReport* report = nullptr)
    {
        auto transcript = std::make_shared<Transcript>(proof);
        [[maybe_unused]] auto init = transcript->template receive_from_prover<FF>("Init");
        const size_t log_n = vk.log_dyadic_size;

        Phase phase(builder, report);

        absorb_vk(builder, transcript, vk);
        std::vector<FF> public_inputs;
        public_inputs.reserve(vk.num_public_inputs);
        for (size_t i = 0; i < vk.num_public_inputs; ++i) {
            public_inputs.push_back(transcript->template receive_from_prover<FF>("public_input_" + std::to_string(i)));
        }

        // Three witness commitment rounds, interleaved with the challenges that depend on them.
        const FF wires_root = receive_root(transcript, "HONK:wires");
        auto [eta, rom_logup_gamma] =
            transcript->template get_challenges<FF>(std::array<std::string, 2>{ "eta", "rom_logup_gamma" });
        RelationParameters<FF> relation_parameters;
        relation_parameters.eta = eta;
        relation_parameters.eta_two = eta * eta;
        relation_parameters.eta_three = relation_parameters.eta_two * eta;
        relation_parameters.rom_logup_gamma = rom_logup_gamma;

        const FF counts_w4_root = receive_root(transcript, "HONK:counts_w4");
        auto [beta, gamma] = transcript->template get_challenges<FF>(std::array<std::string, 2>{ "beta", "gamma" });
        relation_parameters.beta = beta;
        relation_parameters.beta_sqr = beta * beta;
        relation_parameters.beta_cube = relation_parameters.beta_sqr * beta;
        relation_parameters.gamma = gamma;
        relation_parameters.public_input_delta =
            compute_public_input_delta<Flavor>(public_inputs, beta, gamma, FF(uint64_t(vk.pub_inputs_offset)));

        const FF inverses_z_perm_root = receive_root(transcript, "HONK:inverses_z_perm");
        phase.mark("honk: transcript + public inputs");

        const FF alpha = transcript->template get_challenge<FF>("alpha");
        const std::vector<FF> gate_challenges =
            transcript->template get_dyadic_powers_of_challenge<FF>("Sumcheck:gate_challenge", log_n);

        // bb's own sumcheck verifier: for a recursive flavor its round checks are `assert_equal`,
        // so a satisfiable circuit is a passing sumcheck.
        SumcheckVerifier<Flavor> sumcheck(transcript, alpha, log_n);
        SumcheckOutput<Flavor> sumcheck_output = sumcheck.verify(relation_parameters, gate_challenges);
        phase.mark("honk: sumcheck + relations");

        // A virtual precomputed column is not committed, so its claimed evaluation has to equal the
        // value the verifier computes for itself.
        const auto unshifted_evaluations = sumcheck_output.claimed_evaluations.get_unshifted();
        for (size_t e = 0; e < NUM_PRECOMPUTED; ++e) {
            if (((vk.virtual_mask >> e) & 1) == 1) {
                unshifted_evaluations[e].assert_equal(virtual_column_evaluation(vk, e, sumcheck_output.challenge),
                                                      "transparent honk virtual column");
            }
        }
        phase.mark("honk: virtual columns");

        typename WhirVerifier::Claims claims;
        claims.group_num_columns = { vk.num_committed_precomputed(),
                                     NativeHonk::WITNESS_GROUP_COLUMNS[0],
                                     NativeHonk::WITNESS_GROUP_COLUMNS[1],
                                     NativeHonk::WITNESS_GROUP_COLUMNS[2] };
        // The precomputed root is a constant of the aggregator; the witness roots are proof data.
        claims.group_roots = { FF(vk.precomputed_commitment), wires_root, counts_w4_root, inverses_z_perm_root };
        NativeHonk::append_unshifted_refs(claims.unshifted, vk.virtual_mask);
        NativeHonk::append_unshifted_evaluations(claims.unshifted_evaluations, unshifted_evaluations, vk.virtual_mask);
        NativeHonk::append_shifted_refs(claims.to_be_shifted);
        const auto shifted_evaluations = sumcheck_output.claimed_evaluations.get_shifted();
        claims.shifted_evaluations.assign(shifted_evaluations.begin(), shifted_evaluations.end());

        WhirVerifier::verify(builder, config, claims, sumcheck_output.challenge, transcript, report);
        return public_inputs;
    }

  private:
    using Phase = typename WhirVerifier::Phase;

    static FF receive_root(const std::shared_ptr<Transcript>& transcript, const std::string& label)
    {
        return transcript->template receive_from_prover<FF>(label);
    }

    /**
     * @brief `TransparentHonk::absorb_vk`: public metadata, fixed when the verifier circuit is built.
     * @details Each field is a circuit constant promoted to a *fixed* witness, because the stdlib
     * sponge refuses constants — a Poseidon2 custom gate reads witness indices. `fix_witness` is what
     * keeps it public data: the value is pinned by a gate, so it is a constant in everything but
     * representation, and a prover cannot move it.
     */
    static void absorb_vk(Builder& builder,
                          const std::shared_ptr<Transcript>& transcript,
                          const NativeVerificationKey& vk)
    {
        const auto fixed = [&](const bb::fr& value) {
            FF element(value);
            element.convert_constant_to_fixed_witness(&builder);
            return element;
        };
        transcript->add_to_hash_buffer("vk_log_dyadic_size", fixed(bb::fr(uint64_t(vk.log_dyadic_size))));
        transcript->add_to_hash_buffer("vk_num_public_inputs", fixed(bb::fr(uint64_t(vk.num_public_inputs))));
        transcript->add_to_hash_buffer("vk_pub_inputs_offset", fixed(bb::fr(uint64_t(vk.pub_inputs_offset))));
        transcript->add_to_hash_buffer("vk_virtual_mask", fixed(bb::fr(uint64_t(vk.virtual_mask))));
        transcript->add_to_hash_buffer("vk_lagrange_first_row", fixed(bb::fr(uint64_t(vk.lagrange_first_row))));
        transcript->add_to_hash_buffer("vk_lagrange_last_row", fixed(bb::fr(uint64_t(vk.lagrange_last_row))));
        transcript->add_to_hash_buffer("vk_root", fixed(vk.precomputed_commitment));
    }

    /**
     * @brief What sumcheck must have claimed for a virtual precomputed entity, in circuit form.
     * @details Identically-zero columns claim zero; the two Lagrange point indicators claim
     * `eq(bits(row), u)`, whose row is fixed by the verification key, so every factor is a choice
     * between `u_i` and `1 - u_i` made at circuit-build time.
     */
    static FF virtual_column_evaluation(const NativeVerificationKey& vk, size_t entity, std::span<const FF> u)
    {
        if (entity != NativeHonk::LAGRANGE_FIRST_IDX && entity != NativeHonk::LAGRANGE_LAST_IDX) {
            return FF(0);
        }
        const size_t row = entity == NativeHonk::LAGRANGE_FIRST_IDX ? vk.lagrange_first_row : vk.lagrange_last_row;
        FF eq(1);
        for (size_t i = 0; i < u.size(); ++i) {
            eq = eq * (((row >> i) & 1) == 1 ? u[i] : FF(1) - u[i]);
        }
        return eq;
    }
};

} // namespace bb::whir::recursion
