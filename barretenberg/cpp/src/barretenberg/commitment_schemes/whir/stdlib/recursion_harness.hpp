#pragma once

#include "barretenberg/commitment_schemes/whir/stdlib/whir_recursive_verifier.hpp"
#include "barretenberg/polynomials/polynomial.hpp"
#include "barretenberg/stdlib_circuit_builders/ultra_circuit_builder.hpp"
#include "barretenberg/transcript/transcript.hpp"

#include <span>
#include <vector>

namespace bb::whir::recursion {

/**
 * @brief Proving WHIR openings of an arbitrary column layout and verifying them in a circuit.
 * @details Shared by the recursive verifier's tests and its gate-count sweep so both drive the
 * identical statement; the sweep would otherwise be measuring a different protocol from the one
 * under test.
 */
template <typename StdlibHasher = StdlibPoseidon2Hasher<UltraCircuitBuilder>> class RecursionHarness_ {
  public:
    using Builder = UltraCircuitBuilder;
    using FF = stdlib::field_t<Builder>;
    using Verifier = WhirRecursiveVerifier<Builder, StdlibHasher>;
    using CircuitTranscript = StdlibTranscript<Builder>;
    // The in-circuit hasher names the native one it mirrors, so the two can never be paired wrongly.
    using Hasher = typename StdlibHasher::NativeHasher;
    using Digest = typename Hasher::Digest;
    using CK = WhirCommitmentKey<Hasher>;
    using NativeProver = WhirProver<Hasher>;
    using NativeVerifier = WhirVerifier<Hasher>;

    /** @brief A proven WHIR instance: the statement, the commitment roots, and the proof. */
    struct Instance {
        WhirConfig config;
        std::vector<fr> u;
        std::vector<size_t> group_num_columns;
        std::vector<Digest> roots;
        std::vector<WhirColumnRef> unshifted;
        std::vector<fr> unshifted_evaluations;
        std::vector<WhirColumnRef> to_be_shifted;
        std::vector<fr> shifted_evaluations;
        HonkProof proof;
        // Whether the roots travel on the proof stream (the standalone protocol) or are part of the
        // statement (the Honk integration, and an aggregator that hard-codes an inner VK root).
        bool roots_on_proof_stream = true;
    };

    /**
     * @brief Commit `group_columns[g]` columns in group `g`, claim every column at one point, and
     * additionally claim the columns named by `shifted` in their shifted form.
     * @details The layout knob is what lets the sweep run both shapes that matter: a transparent
     * UltraHonk opening (four groups, ~30 columns) and ProveKit's (one group, one column).
     */
    static Instance make_instance(const WhirConfig& config,
                                  const std::vector<size_t>& group_columns,
                                  const std::vector<WhirColumnRef>& shifted = {},
                                  bool roots_on_proof_stream = true)
    {
        CK ck(config);
        const size_t n = size_t(1) << config.num_payload_variables;

        Instance instance;
        instance.config = config;
        instance.group_num_columns = group_columns;
        instance.roots_on_proof_stream = roots_on_proof_stream;
        instance.u = random_array(config.num_payload_variables);

        // A column claimed shifted must have a zero constant term: the shifted virtual oracle is the
        // codeword scaled by x^-1, which represents the coefficient-shifted array only then.
        std::vector<std::vector<std::vector<fr>>> arrays(group_columns.size());
        for (size_t g = 0; g < group_columns.size(); ++g) {
            for (size_t c = 0; c < group_columns[g]; ++c) {
                arrays[g].push_back(random_array(n));
            }
        }
        std::vector<std::vector<bool>> is_shifted(group_columns.size());
        for (size_t g = 0; g < group_columns.size(); ++g) {
            is_shifted[g].assign(group_columns[g], false);
        }
        for (const WhirColumnRef& reference : shifted) {
            is_shifted[reference.group][reference.column] = true;
            arrays[reference.group][reference.column][0] = fr::zero();
        }

        std::vector<WhirGroupData<Hasher>> groups;
        for (size_t g = 0; g < group_columns.size(); ++g) {
            groups.push_back(ck.commit_group(arrays[g], is_shifted[g]));
        }

        typename NativeProver::Claims prover_claims;
        prover_claims.send_roots = roots_on_proof_stream;
        for (auto& group : groups) {
            prover_claims.groups.push_back(&group);
            instance.roots.push_back(group.tree.root());
        }
        for (size_t g = 0; g < group_columns.size(); ++g) {
            for (size_t c = 0; c < group_columns[g]; ++c) {
                prover_claims.unshifted.push_back({ g, c });
                prover_claims.unshifted_evaluations.push_back(mle(arrays[g][c], instance.u));
            }
        }
        for (const WhirColumnRef& reference : shifted) {
            prover_claims.to_be_shifted.push_back(reference);
            prover_claims.shifted_evaluations.push_back(
                shifted_mle(arrays[reference.group][reference.column], instance.u));
        }
        instance.unshifted = prover_claims.unshifted;
        instance.unshifted_evaluations = prover_claims.unshifted_evaluations;
        instance.to_be_shifted = prover_claims.to_be_shifted;
        instance.shifted_evaluations = prover_claims.shifted_evaluations;

        auto prover_transcript = NativeTranscript::test_prover_init_empty();
        NativeProver::prove(ck, prover_claims, instance.u, prover_transcript);
        instance.proof = prover_transcript->export_proof();
        return instance;
    }

    static bool verify_natively(const Instance& instance)
    {
        typename NativeVerifier::Claims claims{ .group_num_columns = instance.group_num_columns,
                                                .unshifted = instance.unshifted,
                                                .unshifted_evaluations = instance.unshifted_evaluations,
                                                .to_be_shifted = instance.to_be_shifted,
                                                .shifted_evaluations = instance.shifted_evaluations,
                                                .group_roots = instance.roots_on_proof_stream ? std::vector<Digest>{}
                                                                                              : instance.roots };
        auto transcript = std::make_shared<NativeTranscript>(instance.proof);
        [[maybe_unused]] auto init = transcript->template receive_from_prover<fr>("Init");
        return NativeVerifier::verify(instance.config, claims, instance.u, transcript);
    }

    /** @brief Constrain `instance` inside `builder`; several calls aggregate several proofs. */
    static void verify_in_circuit(Builder& builder, const Instance& instance, GateReport* report = nullptr)
    {
        typename CircuitTranscript::Proof proof;
        proof.reserve(instance.proof.size());
        for (const fr& element : instance.proof) {
            proof.push_back(FF::from_witness(&builder, element));
        }
        auto transcript = std::make_shared<CircuitTranscript>(proof);
        [[maybe_unused]] auto init = transcript->template receive_from_prover<FF>("Init");

        typename Verifier::Claims claims;
        claims.group_num_columns = instance.group_num_columns;
        if (!instance.roots_on_proof_stream) {
            // A root the aggregator knows at build time is a circuit constant, which is what binds
            // it: `assert_equal` against a constant fixes the witness, where two witnesses would
            // only be copy-constrained to each other.
            for (const Digest& root : instance.roots) {
                const auto limbs = Hasher::digest_to_fields(root);
                std::vector<FF> constants(limbs.begin(), limbs.end());
                claims.group_roots.push_back(StdlibHasher::from_fields(builder, constants));
            }
        }
        claims.unshifted = instance.unshifted;
        for (const fr& value : instance.unshifted_evaluations) {
            claims.unshifted_evaluations.push_back(statement_witness(builder, value));
        }
        claims.to_be_shifted = instance.to_be_shifted;
        for (const fr& value : instance.shifted_evaluations) {
            claims.shifted_evaluations.push_back(statement_witness(builder, value));
        }
        std::vector<FF> u;
        u.reserve(instance.u.size());
        for (const fr& coordinate : instance.u) {
            u.push_back(statement_witness(builder, coordinate));
        }

        Verifier::verify(builder, instance.config, claims, u, transcript, report);
    }

    static Builder build_circuit(const Instance& instance, GateReport* report = nullptr)
    {
        Builder builder;
        verify_in_circuit(builder, instance, report);
        return builder;
    }

    /** @brief The recursion profile at the given schedule parameters. */
    static WhirConfig config(size_t num_variables,
                             size_t security_bits,
                             size_t log_inv_rate,
                             size_t k,
                             size_t k0,
                             size_t pow_bits,
                             WhirSoundness soundness = WhirSoundness::CONJECTURED_LIST)
    {
        WhirConfig result = WhirConfig::create(num_variables,
                                               security_bits,
                                               log_inv_rate,
                                               k,
                                               /*final_poly_bits=*/4,
                                               soundness,
                                               /*zk=*/false,
                                               /*max_stack_bits=*/0,
                                               k0,
                                               pow_bits);
        result.enable_recursion_profile();
        return result;
    }

    /**
     * @brief A statement witness: an input with no Fiat-Shamir provenance of its own.
     * @details Claimed evaluations and the opening point reach a recursive verifier from the outer
     * circuit — as public inputs, or bound by earlier rounds in the Honk integration — so they are
     * not free witnesses and are allowed to meet transcript-derived values.
     */
    static FF statement_witness(Builder& builder, const fr& value)
    {
        FF element = FF::from_witness(&builder, value);
        element.set_origin_tag(OriginTag::constant());
        return element;
    }

  private:
    static std::vector<fr> random_array(size_t size)
    {
        std::vector<fr> array(size);
        for (fr& value : array) {
            value = fr::random_element();
        }
        return array;
    }

    static fr mle(std::span<const fr> array, std::span<const fr> u) { return Polynomial<fr>(array).evaluate_mle(u); }

    /** @brief MLE of the shifted array (a₁, ..., a_{n-1}, 0) at u — Honk's shifted-claim semantics. */
    static fr shifted_mle(std::span<const fr> array, std::span<const fr> u)
    {
        std::vector<fr> shifted(array.begin() + 1, array.end());
        shifted.push_back(fr::zero());
        return mle(shifted, u);
    }
};

using RecursionHarness = RecursionHarness_<>;

} // namespace bb::whir::recursion
