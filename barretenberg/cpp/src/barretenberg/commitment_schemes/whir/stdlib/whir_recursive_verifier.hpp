#pragma once

#include "barretenberg/commitment_schemes/whir/whir.hpp"
#include "barretenberg/stdlib/hash/poseidon2/poseidon2_permutation.hpp"
#include "barretenberg/stdlib/primitives/bool/bool.hpp"
#include "barretenberg/stdlib/primitives/field/field.hpp"
#include "barretenberg/stdlib/primitives/memory/rom_table.hpp"
#include "barretenberg/stdlib/primitives/witness/witness.hpp"
#include "barretenberg/transcript/transcript.hpp"

#include <functional>
#include <map>
#include <span>
#include <string>
#include <vector>

namespace bb::whir::recursion {

/**
 * @brief Gate counts attributed to the phases of one in-circuit WHIR verification.
 * @details The whole point of the exercise is knowing where the constraints go, so the verifier
 * takes one of these and stamps the builder's gate count at each phase boundary. Ordering is
 * insertion order; a phase entered twice accumulates.
 */
class GateReport {
  public:
    void begin(size_t gates) { last_ = gates; }
    void mark(const std::string& phase, size_t gates)
    {
        if (!totals_.contains(phase)) {
            order_.push_back(phase);
        }
        totals_[phase] += gates - last_;
        last_ = gates;
    }
    size_t total() const
    {
        size_t sum = 0;
        for (const auto& [_, gates] : totals_) {
            sum += gates;
        }
        return sum;
    }
    const std::vector<std::string>& phases() const { return order_; }
    size_t gates(const std::string& phase) const { return totals_.at(phase); }

  private:
    size_t last_ = 0;
    std::vector<std::string> order_;
    std::map<std::string, size_t> totals_;
};

/**
 * @brief In-circuit Merkle hashing, bit-identical to `bb::whir::Poseidon2CompressionHasher`.
 * @details Both halves of the hasher are a single `Poseidon2Permutation` call per three absorbed
 * values (leaves) or per node, which is what makes it the cheapest of bb's WHIR hashers to verify
 * recursively: Ultra gives Poseidon2 its own custom gates, where a Blake3 or SHA-256 compression
 * costs thousands of constraints.
 */
template <typename Builder> class StdlibPoseidon2Hasher {
  public:
    using FF = stdlib::field_t<Builder>;
    using Permutation = stdlib::Poseidon2Permutation<Builder>;
    static constexpr size_t RATE = 3;

    static FF hash_leaf(Builder& builder, std::span<const FF> values)
    {
        std::array<FF, 4> state{ FF(0), FF(0), FF(0), iv(builder, Hasher::leaf_iv(values.size(), false)) };
        OriginTag absorbed = OriginTag::constant();
        size_t slot = 0;
        auto permute = [&] {
            state = Permutation::permutation(&builder, state);
            retag(state, absorbed);
        };
        for (const FF& value : values) {
            state[slot] += value;
            absorbed = OriginTag(absorbed, value.get_origin_tag());
            if (++slot == RATE) {
                permute();
                slot = 0;
            }
        }
        if (slot != 0 || values.empty()) {
            permute();
        }
        return state[0];
    }

    static FF hash_node(Builder& builder, const FF& left, const FF& right)
    {
        std::array<FF, 4> state{ left, right, FF(0), iv(builder, Hasher::node_iv()) };
        state = Permutation::permutation(&builder, state);
        retag(state, OriginTag(left.get_origin_tag(), right.get_origin_tag()));
        return state[0];
    }

    /** @brief Permutations one leaf of `width` values costs; the cost model's unit. */
    static size_t leaf_permutations(size_t width) { return width == 0 ? 1 : (width + RATE - 1) / RATE; }

  private:
    using Hasher = Poseidon2CompressionHasher;

    /**
     * @brief Give a permuted state the Fiat-Shamir provenance of everything absorbed into it.
     * @details `Poseidon2Permutation` emits fresh witnesses and leaves them untagged — the
     * transcript, its only other caller, tags the challenges it derives itself. A sponge state is a
     * function of what it has absorbed, so it inherits the merged tag; without this the state looks
     * like a free witness the moment the next value is absorbed into it.
     */
    static void retag(std::array<FF, 4>& state, const OriginTag& tag)
    {
        for (FF& element : state) {
            element.set_origin_tag(tag);
        }
    }

    // Poseidon2's custom gates read witness indices, so the capacity separator has to be a witness
    // rather than a constant, exactly as `stdlib::FieldSponge` does with its length IV.
    static FF iv(Builder& builder, const bb::fr& value)
    {
        FF constant(value);
        constant.convert_constant_to_fixed_witness(&builder);
        return constant;
    }
};

/**
 * @brief Recursive verifier for a WHIR batched opening, the in-circuit mirror of
 * `bb::whir::WhirVerifier::verify`.
 *
 * @details Restricted to the shape a transparent Honk proof opens with and the one ProveKit's WHIR
 * produces: interleaved columns (no stacking), no zero knowledge, every claim at the single
 * sumcheck point. `WhirConfig::enable_recursion_profile` must be on, because three of its settings
 * are what make the circuit cheap:
 *
 *  - *per-query authentication paths*. The batched layout walks the distinct leaf set the queries
 *    induce, whose shape depends on the query values; in-circuit that means a witness-indexed read
 *    at every step of every level. Per-query paths make every index a compile-time constant, and
 *    the extra digests are unhashed proof data, which costs a recursive verifier nothing.
 *  - *a Merkle cap*. Each path stops `merkle_cap_levels` below the root and lands in a table the
 *    verifier folds to the root once, trading `2^c - 1` hashes for `c` per query.
 *  - *a Poseidon2 grind*. Blake3 proof of work costs more in-circuit than the queries it removes.
 *
 * Everything is native BN254 Fr arithmetic — no non-native field emulation and no elliptic-curve
 * work at all, which is the structural reason a hash-based scheme recurses cheaply into UltraHonk.
 */
template <typename Builder> class WhirRecursiveVerifier {
  public:
    using FF = stdlib::field_t<Builder>;
    using Bool = stdlib::bool_t<Builder>;
    using Witness = stdlib::witness_t<Builder>;
    using Transcript = StdlibTranscript<Builder>;
    using Hasher = StdlibPoseidon2Hasher<Builder>;
    using RomTable = stdlib::rom_table<Builder>;

    /** @brief The verifier's statement: the same shape as `WhirVerifier::Claims`, in circuit types. */
    struct Claims {
        std::vector<size_t> group_num_columns;
        std::vector<FF> group_roots; // empty reads them off the proof stream instead
        std::vector<WhirColumnRef> unshifted;
        std::vector<FF> unshifted_evaluations;
        std::vector<WhirColumnRef> to_be_shifted;
        std::vector<FF> shifted_evaluations;
    };

    /**
     * @brief Constrain that the proof behind `transcript` opens `claims` at `u`.
     * @details Every check the native verifier answers with `return false` is an in-circuit
     * assertion, so a satisfiable circuit is a verifying proof.
     */
    static void verify(Builder& builder,
                       const WhirConfig& config,
                       const Claims& claims,
                       std::span<const FF> u,
                       const std::shared_ptr<Transcript>& transcript,
                       GateReport* report = nullptr)
    {
        BB_ASSERT(config.per_query_openings, "the recursive verifier needs per-query openings");
        BB_ASSERT(!config.stack_columns, "column stacking is not supported by the recursive verifier");
        BB_ASSERT(!config.zk, "zk WHIR is not supported by the recursive verifier");
        BB_ASSERT_EQ(u.size(), config.num_payload_variables, "opening point size mismatch");

        Phase phase(builder, report);

        // ---- Commitments, claims and the batching challenges --------------------------------
        std::vector<FF> roots = claims.group_roots;
        for (size_t g = roots.size(); g < claims.group_num_columns.size(); ++g) {
            roots.push_back(transcript->template receive_from_prover<FF>(detail::whir_label("root", g)));
        }

        // The non-stacked plan of `build_stacked_plan`: one constituent per claimed column, the
        // unshifted ones first, all claimed at the single point u.
        struct Constituent {
            size_t tree;
            size_t leaf_column;
            bool shifted;
            FF value;
        };
        std::vector<Constituent> constituents;
        for (size_t j = 0; j < claims.unshifted.size(); ++j) {
            constituents.push_back(
                { claims.unshifted[j].group, claims.unshifted[j].column, false, claims.unshifted_evaluations[j] });
        }
        for (size_t j = 0; j < claims.to_be_shifted.size(); ++j) {
            constituents.push_back(
                { claims.to_be_shifted[j].group, claims.to_be_shifted[j].column, true, claims.shifted_evaluations[j] });
        }

        std::vector<FF> z_ood_initial(config.num_ood_samples);
        std::vector<std::vector<FF>> ood_values(config.num_ood_samples);
        for (size_t s = 0; s < config.num_ood_samples; ++s) {
            z_ood_initial[s] = transcript->template get_challenge<FF>(detail::whir_label("z_ood_init", s));
            for (size_t c = 0; c < constituents.size(); ++c) {
                ood_values[s].push_back(
                    transcript->template receive_from_prover<FF>(detail::whir_label("y_ood_init", s, c)));
            }
        }
        const FF rho = transcript->template get_challenge<FF>("WHIR:rho");
        const FF gamma_claims = transcript->template get_challenge<FF>("WHIR:gamma_claims");

        std::vector<FF> rho_powers;
        rho_powers.reserve(constituents.size());
        FF rho_power(1);
        for (size_t c = 0; c < constituents.size(); ++c) {
            rho_powers.push_back(rho_power);
            rho_power = rho_power * rho;
        }

        // σ₀ = ∑_t γᵗ ∑_c ρᶜ R_c(p_t), one weight term per point: the sumcheck claim entering round 0.
        std::vector<Term> terms;
        FF sigma(0);
        FF gamma_power(1);
        {
            FF oracle_value(0);
            for (size_t c = 0; c < constituents.size(); ++c) {
                oracle_value += rho_powers[c] * constituents[c].value;
            }
            sigma += gamma_power * oracle_value;
            terms.push_back(Term::eq_term(gamma_power));
            gamma_power = gamma_power * gamma_claims;
        }
        for (size_t s = 0; s < config.num_ood_samples; ++s) {
            FF oracle_value(0);
            for (size_t c = 0; c < constituents.size(); ++c) {
                oracle_value += rho_powers[c] * ood_values[s][c];
            }
            sigma += gamma_power * oracle_value;
            terms.push_back(Term::pow_term(gamma_power, z_ood_initial[s]));
            gamma_power = gamma_power * gamma_claims;
        }
        phase.mark("claims + batching");

        // ---- Fold-and-commit iterations -----------------------------------------------------
        std::optional<FF> folded_root;
        size_t bound = 0; // sumcheck variables consumed, i.e. the eq term's next coordinate

        for (size_t i = 0; i < config.rounds.size(); ++i) {
            const WhirRound& round = config.rounds[i];
            const size_t k = round.folding_factor_bits;

            std::vector<FF> alphas;
            for (size_t j = 0; j < k; ++j) {
                const auto h =
                    transcript->template receive_from_prover<std::array<FF, 3>>(detail::whir_label("sc", i, j));
                (h[0] + h[1]).assert_equal(sigma, "whir sumcheck round sum");
                const FF alpha = transcript->template get_challenge<FF>(detail::whir_label("alpha", i, j));
                sigma = evaluate_quadratic(h, alpha);
                for (Term& term : terms) {
                    term.bind(alpha, u, bound);
                }
                alphas.push_back(alpha);
                ++bound;
            }
            phase.mark("whir sumcheck");

            const FF next_root = transcript->template receive_from_prover<FF>(detail::whir_label("root_g", i + 1));
            const FF z_ood = transcript->template get_challenge<FF>(detail::whir_label("z_ood", i));
            const FF y_ood = transcript->template receive_from_prover<FF>(detail::whir_label("y_ood", i));

            const size_t index_bits = round.log_domain_size - k;
            check_grinding(builder, transcript, config, detail::whir_label("pow", i), detail::whir_label("nonce", i));
            phase.mark("grinding");

            const std::vector<QueryIndex> indices = draw_indices(
                builder,
                transcript,
                [&](size_t chunk) { return detail::whir_label("query", i, chunk); },
                round.num_queries,
                index_bits);
            phase.mark("query index extraction");

            const bb::fr omega = bb::fr::get_root_of_unity(round.log_domain_size);
            const bb::fr eta_inv = omega.pow(uint256_t(uint64_t(1)) << index_bits).invert();

            std::vector<std::vector<FF>> cosets;
            if (folded_root) {
                cosets = read_single_column_openings(
                    builder, transcript, config, *folded_root, indices, index_bits, k, phase);
            } else {
                cosets = read_round0_openings(builder,
                                              transcript,
                                              config,
                                              roots,
                                              claims,
                                              constituents,
                                              rho_powers,
                                              indices,
                                              omega,
                                              eta_inv,
                                              phase);
            }

            // 5. Fold each queried coset and γ-batch the results into the next round's claim.
            const FF gamma = transcript->template get_challenge<FF>(detail::whir_label("gamma", i));
            FF gamma_power_round = gamma;
            sigma += gamma_power_round * y_ood;
            terms.push_back(Term::pow_term(gamma_power_round, z_ood));
            for (size_t s = 0; s < round.num_queries; ++s) {
                const FF folded = fold_coset(cosets[s], indices[s].power_inverse(omega), eta_inv, alphas);
                gamma_power_round = gamma_power_round * gamma;
                sigma += gamma_power_round * folded;
                // The folded oracle's domain point ω^{idx·2^k} = (ω^{2^k})^{idx}.
                terms.push_back(Term::pow_term(gamma_power_round, indices[s].power(omega.pow(uint64_t(1) << k))));
            }
            phase.mark("coset folding");
            folded_root = next_root;
        }

        // ---- Final phase: the clear polynomial ------------------------------------------------
        const size_t final_size = size_t(1) << config.final_round.num_variables;
        std::vector<FF> final_poly;
        for (size_t j = 0; j < final_size; ++j) {
            final_poly.push_back(transcript->template receive_from_prover<FF>(detail::whir_label("final", j)));
        }
        const size_t final_k = config.final_round.folding_factor_bits;
        const size_t index_bits = config.final_round.log_domain_size - final_k;
        const bb::fr omega = bb::fr::get_root_of_unity(config.final_round.log_domain_size);
        const bb::fr eta = omega.pow(uint256_t(uint64_t(1)) << index_bits);
        const bb::fr eta_inv = eta.invert();

        check_grinding(builder, transcript, config, "WHIR:fpow", "WHIR:fnonce");
        const std::vector<QueryIndex> final_indices = draw_indices(
            builder,
            transcript,
            [](size_t chunk) { return detail::whir_label("fquery", chunk); },
            config.final_round.num_queries,
            index_bits);
        phase.mark("query index extraction");

        std::vector<std::vector<FF>> final_cosets;
        if (folded_root) {
            final_cosets = read_single_column_openings(
                builder, transcript, config, *folded_root, final_indices, index_bits, final_k, phase);
        } else {
            final_cosets = read_round0_openings(builder,
                                                transcript,
                                                config,
                                                roots,
                                                claims,
                                                constituents,
                                                rho_powers,
                                                final_indices,
                                                omega,
                                                eta_inv,
                                                phase);
        }
        for (size_t s = 0; s < final_indices.size(); ++s) {
            const std::vector<FF> expected =
                evaluate_on_coset(final_poly, final_indices[s].power(omega), eta, final_cosets[s].size());
            for (size_t t = 0; t < final_cosets[s].size(); ++t) {
                final_cosets[s][t].assert_equal(expected[t], "whir final oracle consistency");
            }
        }
        phase.mark("final oracle consistency");

        // Every accumulated weight term against the final polynomial.
        FF total(0);
        std::span<const FF> u_tail = u.subspan(bound);
        total += terms[0].coeff * multilinear_at(final_poly, u_tail);
        for (size_t t = 1; t < terms.size(); ++t) {
            total += terms[t].coeff * evaluate_coefficients(final_poly, terms[t].power);
        }
        total.assert_equal(sigma, "whir final claim");
        phase.mark("final claim");
    }

  private:
    /** @brief RAII-free helper that stamps the builder's gate count at each phase boundary. */
    class Phase {
      public:
        Phase(Builder& builder, GateReport* report)
            : builder_(builder)
            , report_(report)
        {
            if (report_ != nullptr) {
                report_->begin(builder_.num_gates());
            }
        }
        void mark(const std::string& name)
        {
            if (report_ != nullptr) {
                report_->mark(name, builder_.num_gates());
            }
        }

      private:
        Builder& builder_;
        GateReport* report_;
    };

    /**
     * @brief A query index: its value, its bits, and the powers of a domain generator it induces.
     * @details The bits are needed anyway — one per Merkle level to place the sibling, and one per
     * level of the square-and-multiply that turns the index into its domain point — so the index is
     * bit-decomposed once and everything else reads those bits.
     */
    struct QueryIndex {
        FF value;
        std::vector<Bool> bits;

        /** @brief `base^index`, from the bits: each factor is a free linear term in one bit. */
        FF power(const bb::fr& base) const
        {
            FF result(1);
            bb::fr square = base;
            for (size_t i = 0; i < bits.size(); ++i) {
                // (square - 1)·bit + 1 is a linear combination, so only the accumulation is a gate.
                result = result * (FF(bits[i]) * (square - bb::fr(1)) + FF(1));
                square = square.sqr();
            }
            return result;
        }
        FF power_inverse(const bb::fr& base) const { return power(base.invert()); }

        /** @brief The index shifted right by `low`, as a field element (the Merkle cap address). */
        FF high_part(size_t low) const
        {
            FF result(0);
            bb::fr scale(1);
            for (size_t i = low; i < bits.size(); ++i) {
                result += FF(bits[i]) * scale;
                scale += scale;
            }
            return result;
        }
    };

    /**
     * @brief A weight term of the running claim, in the two forms the schedule produces.
     * @details Every term ends the protocol with the final round's variables unbound, so the two
     * forms collapse to a scalar and one point: the `eq(u, X)` term contributes the final
     * polynomial's multilinear extension at the unbound tail of u, and a `pow_y(X)` term
     * contributes the final polynomial — read as coefficients — at y^{2^bound}. Binding a variable
     * multiplies the scalar and squares the running power.
     */
    struct Term {
        FF coeff;
        FF power;
        bool is_eq;

        static Term eq_term(const FF& coeff) { return { coeff, FF(0), true }; }
        static Term pow_term(const FF& coeff, const FF& y) { return { coeff, y, false }; }

        void bind(const FF& alpha, std::span<const FF> u, size_t index)
        {
            if (is_eq) {
                // (1-u) + (2u-1)·α
                coeff = coeff * (FF(1) - u[index] + (u[index] + u[index] - FF(1)) * alpha);
            } else {
                coeff = coeff * (FF(1) + (power - FF(1)) * alpha);
                power = power * power;
            }
        }
    };

    /**
     * @brief `count` query indices, the in-circuit mirror of `detail::draw_query_indices`.
     * @details Each challenge is split into `detail::indices_per_challenge` base-2^b digits plus a
     * high remainder: `challenge = Σ_j index_j·2^{jb} + rest·2^{nb}`. The digits are bit-decomposed
     * because the Merkle walk and the domain-point exponentiation both consume them a bit at a
     * time, so the decomposition is work the verifier needs anyway; what the packing removes is one
     * sponge permutation and one pair of range constraints *per query*.
     *
     * `rest` is pinned to `[0, ⌊r/2^{nb}⌋)` by two range constraints. Without the second one the
     * decomposition is not unique — `challenge + r` is below 2^254 for about a third of all
     * challenges — and a prover free to choose between two index sets per challenge would claw back
     * a bit of soundness per query. Pinning `rest` strictly below ⌊r/2^{nb}⌋ rather than at it
     * rejects an honest challenge with probability about `2^{nb}/r`; `detail::QUERY_DIGIT_BITS`
     * bounds `nb` at 128 to keep that at 2^-128, because Fiat-Shamir gives an honest prover no
     * second draw and a rejected challenge means an unprovable proof.
     */
    static std::vector<QueryIndex> draw_indices(Builder& builder,
                                                const std::shared_ptr<Transcript>& transcript,
                                                const std::function<std::string(size_t)>& label,
                                                size_t count,
                                                size_t index_bits)
    {
        BB_ASSERT_LTE(index_bits, size_t(64), "query index wider than the low limb");
        constexpr size_t FIELD_BITS = 254;
        const size_t per_challenge = detail::indices_per_challenge(index_bits);
        const uint64_t mask = (uint64_t(1) << index_bits) - 1;

        std::vector<QueryIndex> indices;
        indices.reserve(count);
        for (size_t chunk = 0; indices.size() < count; ++chunk) {
            const FF challenge = transcript->template get_challenge<FF>(label(chunk));
            const uint256_t challenge_value(challenge.get_value());
            // Hint witnesses derived from a challenge inherit its Fiat-Shamir provenance; they are
            // functions of it, not independent inputs.
            const OriginTag tag = challenge.get_origin_tag();

            FF sum(0);
            bb::fr scale(1);
            const size_t drawn = std::min(per_challenge, count - indices.size());
            for (size_t j = 0; j < drawn; ++j) {
                const uint64_t index_value = ((challenge_value >> (j * index_bits)).data[0]) & mask;
                QueryIndex index;
                FF digit(0);
                bb::fr digit_scale(1); // the digit's own place value, not its place in the challenge
                for (size_t b = 0; b < index_bits; ++b) {
                    Bool bit(Witness(&builder, ((index_value >> b) & 1) == 1));
                    bit.set_origin_tag(tag);
                    digit += FF(bit) * digit_scale;
                    sum += FF(bit) * scale;
                    digit_scale += digit_scale;
                    scale += scale;
                    index.bits.push_back(bit);
                }
                index.value = digit.normalize();
                indices.push_back(std::move(index));
            }
            // Digits beyond the ones this chunk supplies still belong to the challenge's low part:
            // the remainder starts above the whole `per_challenge`-digit window, so a short final
            // chunk leaves the unused digits inside `rest`.
            const size_t digit_bits = drawn * index_bits;
            const uint256_t rest_value = challenge_value >> digit_bits;
            FF rest(Witness(&builder, bb::fr(rest_value)));
            rest.set_origin_tag(tag);
            (sum + rest * bb::fr(uint256_t(1) << digit_bits)).assert_equal(challenge, "query index decomposition");

            const size_t rest_bits = FIELD_BITS - digit_bits;
            rest.create_range_constraint(rest_bits, "query index high part");
            const uint256_t quotient_bound = (uint256_t(bb::fr::modulus) >> digit_bits);
            (FF(bb::fr(quotient_bound - 1)) - rest).create_range_constraint(rest_bits, "query index canonicity");
        }
        return indices;
    }

    /**
     * @brief Constrain the round's Poseidon2 grind: `Poseidon2(seed, nonce) = q·2^pow_bits`.
     * @details One permutation and one range constraint, against tens of thousands of constraints
     * for the Blake3 form. `pow_bits` bits of the query soundness are supplied by this instead of
     * by queries, and a query costs a whole leaf hash plus a path.
     */
    static void check_grinding(Builder& builder,
                               const std::shared_ptr<Transcript>& transcript,
                               const WhirConfig& config,
                               const std::string& seed_label,
                               const std::string& nonce_label)
    {
        const FF seed = transcript->template get_challenge<FF>(seed_label);
        const FF nonce = transcript->template receive_from_prover<FF>(nonce_label);
        if (config.pow_bits == 0) {
            return;
        }
        BB_ASSERT(config.poseidon2_pow, "the recursive verifier only checks the Poseidon2 grind");
        const FF digest = stdlib::poseidon2<Builder>::hash({ seed, nonce });
        const uint256_t digest_value(digest.get_value());
        FF quotient(Witness(&builder, bb::fr(digest_value >> config.pow_bits)));
        quotient.set_origin_tag(digest.get_origin_tag());
        (quotient * bb::fr(uint256_t(1) << config.pow_bits)).assert_equal(digest, "whir proof of work");

        // `quotient * 2^b = digest` has a second solution whenever `digest + r` is divisible by 2^b
        // and the resulting quotient still fits, which would accept one extra digest residue and
        // hand back a bit of the grind. Pinning `quotient` to [0, ⌊r/2^b⌋] keeps `quotient * 2^b`
        // below r, so the solution is unique; an honest quotient is `digest/2^b < r/2^b` and always
        // inside that range, so nothing is lost on the completeness side.
        const size_t quotient_bits = 254 - config.pow_bits;
        quotient.create_range_constraint(quotient_bits, "whir proof of work quotient");
        const uint256_t quotient_bound = uint256_t(bb::fr::modulus) >> config.pow_bits;
        (FF(bb::fr(quotient_bound)) - quotient).create_range_constraint(quotient_bits, "whir proof of work canonicity");
    }

    /**
     * @brief Authenticate one tree's per-query openings and return each query's leaf values.
     * @details The cap is read first and folded to `root` once; each query then hashes its leaf and
     * climbs `depth - cap_levels` levels, placing the transmitted sibling on the side its own index
     * bit dictates, and the result is looked up in the cap. Every array index here is a
     * compile-time constant except the cap address, which is one ROM read.
     */
    static std::vector<std::vector<FF>> authenticate(Builder& builder,
                                                     const std::shared_ptr<Transcript>& transcript,
                                                     const WhirConfig& config,
                                                     const FF& root,
                                                     std::span<const QueryIndex> indices,
                                                     size_t leaf_width,
                                                     size_t depth,
                                                     Phase& phase)
    {
        const size_t cap_levels = config.cap_levels_for(indices.size(), depth);
        const size_t cap_size = size_t(1) << cap_levels;
        const size_t path_length = depth - cap_levels;
        // The cap is addressed by the index bits the path did not consume, so a short bit vector
        // would silently read a constant-zero slot instead of the node the path landed on.
        for (const QueryIndex& index : indices) {
            BB_ASSERT_EQ(index.bits.size(), depth, "query index width must match the tree depth");
        }

        const std::vector<FF> flat =
            transcript->receive_unhashed_from_prover(cap_size + indices.size() * (leaf_width + path_length));
        size_t cursor = 0;
        std::vector<FF> cap(flat.begin(), flat.begin() + static_cast<std::ptrdiff_t>(cap_size));
        cursor += cap_size;

        // Fold the cap to the committed root: 2^c - 1 hashes once, against c per query.
        {
            std::vector<FF> level = cap;
            while (level.size() > 1) {
                std::vector<FF> above(level.size() / 2);
                for (size_t j = 0; j < above.size(); ++j) {
                    above[j] = Hasher::hash_node(builder, level[2 * j], level[2 * j + 1]);
                }
                level = std::move(above);
            }
            level[0].assert_equal(root, "whir merkle cap");
        }
        phase.mark("merkle: cap fold");
        const RomTable cap_table(cap);

        std::vector<std::vector<FF>> values(indices.size());
        for (size_t s = 0; s < indices.size(); ++s) {
            values[s].assign(flat.begin() + static_cast<std::ptrdiff_t>(cursor),
                             flat.begin() + static_cast<std::ptrdiff_t>(cursor + leaf_width));
            cursor += leaf_width;
            FF digest = Hasher::hash_leaf(builder, values[s]);
            phase.mark("merkle: leaf hashing");
            for (size_t level = 0; level < path_length; ++level) {
                const FF& sibling = flat[cursor++];
                const Bool& bit = indices[s].bits[level];
                // Two normalized selects rather than one select and a linear complement: the
                // permutation normalizes whatever it is handed, so deriving `right` as
                // `digest + sibling - left` moves the gate rather than removing it, and measures
                // slightly worse.
                const FF left = FF::conditional_assign(bit, sibling, digest);
                const FF right = FF::conditional_assign(bit, digest, sibling);
                digest = Hasher::hash_node(builder, left, right);
            }
            phase.mark("merkle: path walk");
            digest.assert_equal(cap_table[indices[s].high_part(path_length)], "whir merkle path");
            phase.mark("merkle: cap lookup");
        }
        return values;
    }

    /** @brief The folded oracle's cosets: one column, so a leaf is the coset itself. */
    static std::vector<std::vector<FF>> read_single_column_openings(Builder& builder,
                                                                    const std::shared_ptr<Transcript>& transcript,
                                                                    const WhirConfig& config,
                                                                    const FF& root,
                                                                    std::span<const QueryIndex> indices,
                                                                    size_t index_bits,
                                                                    size_t k,
                                                                    Phase& phase)
    {
        return authenticate(builder, transcript, config, root, indices, size_t(1) << k, index_bits, phase);
    }

    /**
     * @brief Round 0: authenticate every commitment group and assemble the batched oracle's cosets.
     * @details The batched oracle's value at a coset position is the ρ-combination of every
     * constituent's opened value there, with a shifted constituent scaled by the coset point's
     * inverse — the same `x^{-1}` contract the native `round0_coset` applies.
     */
    static std::vector<std::vector<FF>> read_round0_openings(Builder& builder,
                                                             const std::shared_ptr<Transcript>& transcript,
                                                             const WhirConfig& config,
                                                             const std::vector<FF>& roots,
                                                             const Claims& claims,
                                                             const auto& constituents,
                                                             const std::vector<FF>& rho_powers,
                                                             std::span<const QueryIndex> indices,
                                                             const bb::fr& omega,
                                                             const bb::fr& eta_inv,
                                                             Phase& phase)
    {
        const size_t k = config.initial_folding_factor_bits;
        const size_t arity = size_t(1) << k;
        const size_t tree_index_bits = config.num_variables + config.log_inv_rate - k;

        std::vector<std::vector<std::vector<FF>>> group_values;
        for (size_t g = 0; g < roots.size(); ++g) {
            group_values.push_back(authenticate(builder,
                                                transcript,
                                                config,
                                                roots[g],
                                                indices,
                                                claims.group_num_columns[g] * arity,
                                                tree_index_bits,
                                                phase));
        }

        const bool any_shifted = !claims.to_be_shifted.empty();
        std::vector<std::vector<FF>> cosets(indices.size(), std::vector<FF>(arity, FF(0)));
        for (size_t s = 0; s < indices.size(); ++s) {
            std::vector<FF> x_inverses;
            if (any_shifted) {
                FF x_inv = indices[s].power_inverse(omega);
                for (size_t t = 0; t < arity; ++t) {
                    x_inverses.push_back(x_inv);
                    x_inv = x_inv * eta_inv;
                }
            }
            for (size_t c = 0; c < constituents.size(); ++c) {
                const auto& constituent = constituents[c];
                const std::vector<FF>& leaf = group_values[constituent.tree][s];
                for (size_t t = 0; t < arity; ++t) {
                    FF value = rho_powers[c] * leaf[constituent.leaf_column * arity + t];
                    if (constituent.shifted) {
                        value = value * x_inverses[t];
                    }
                    cosets[s][t] += value;
                }
            }
        }
        phase.mark("batched oracle assembly");
        return cosets;
    }

    /** @brief `Fold(A, α)(x_base^{2^k})` from the coset values, `rs_code.hpp`'s `fold_coset`. */
    static FF fold_coset(std::span<const FF> values,
                         const FF& x_base_inv,
                         const bb::fr& eta_inv,
                         std::span<const FF> alphas)
    {
        static const bb::fr two_inv = bb::fr(2).invert();
        std::vector<FF> current(values.begin(), values.end());
        size_t size = current.size();
        FF base_inv = x_base_inv;
        bb::fr level_eta_inv = eta_inv;
        for (const FF& alpha : alphas) {
            const size_t half = size / 2;
            FF point_inv = base_inv;
            for (size_t t = 0; t < half; ++t) {
                const FF even = (current[t] + current[t + half]) * two_inv;
                const FF odd = (current[t] - current[t + half]) * two_inv * point_inv;
                current[t] = even + alpha * (odd - even);
                point_inv = point_inv * level_eta_inv;
            }
            base_inv = base_inv * base_inv;
            level_eta_inv = level_eta_inv.sqr();
            size = half;
        }
        return current[0];
    }

    /** @brief `detail::evaluate_quadratic`: h given at 0, 1, 2, evaluated at alpha. */
    static FF evaluate_quadratic(const std::array<FF, 3>& h, const FF& alpha)
    {
        static const bb::fr two_inv = bb::fr(2).invert();
        const FF a2 = (h[2] - h[1] - h[1] + h[0]) * two_inv;
        const FF a1 = h[1] - h[0] - a2;
        return h[0] + alpha * (a1 + alpha * a2);
    }

    /**
     * @brief `poly` evaluated at every point of the coset {x, xη, ..., xη^{arity-1}}.
     * @details Horner at each point separately would cost `arity * |poly|` multiplications. Writing
     * `P(xη^t) = Σ_j (poly[j]·x^j)·η^{tj}` shares the `poly[j]·x^j` products across the whole coset
     * and leaves each output a linear combination of them with *constant* coefficients, which an
     * addition gate absorbs several at a time. Used for the final oracle-consistency check, where
     * every query compares its whole opened coset against the clear polynomial.
     */
    static std::vector<FF> evaluate_on_coset(std::span<const FF> poly, const FF& x, const bb::fr& eta, size_t arity)
    {
        std::vector<FF> scaled;
        scaled.reserve(poly.size());
        FF power(1);
        for (const FF& coefficient : poly) {
            scaled.push_back(coefficient * power);
            power = power * x;
        }
        std::vector<FF> values;
        values.reserve(arity);
        for (size_t t = 0; t < arity; ++t) {
            std::vector<FF> terms;
            terms.reserve(scaled.size());
            bb::fr eta_power(1); // η^{t·j}
            const bb::fr eta_t = eta.pow(uint64_t(t));
            for (const FF& term : scaled) {
                terms.push_back(term * eta_power);
                eta_power *= eta_t;
            }
            values.push_back(FF::accumulate(terms));
        }
        return values;
    }

    /** @brief The univariate with coefficients `poly` at `point`, by Horner. */
    static FF evaluate_coefficients(std::span<const FF> poly, const FF& point)
    {
        FF accumulator(0);
        for (size_t i = poly.size(); i > 0; --i) {
            accumulator = accumulator * point + poly[i - 1];
        }
        return accumulator;
    }

    /** @brief The multilinear extension of `poly` (evaluation form, low variable first) at `point`. */
    static FF multilinear_at(std::span<const FF> poly, std::span<const FF> point)
    {
        std::vector<FF> current(poly.begin(), poly.end());
        for (const FF& coordinate : point) {
            std::vector<FF> next(current.size() / 2);
            for (size_t j = 0; j < next.size(); ++j) {
                next[j] = current[2 * j] + coordinate * (current[2 * j + 1] - current[2 * j]);
            }
            current = std::move(next);
        }
        return current[0];
    }
};

} // namespace bb::whir::recursion
