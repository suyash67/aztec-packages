#pragma once

#include "barretenberg/common/assert.hpp"
#include "barretenberg/ecc/curves/bn254/fr.hpp"
#include "barretenberg/numeric/uint256/uint256.hpp"

#include <algorithm>
#include <cstddef>
#include <vector>

namespace bb::whir {

/**
 * @brief Proximity regime assumed for the underlying Reed-Solomon codes; determines the per-round
 * query counts. See README.md §6.
 */
enum class WhirSoundness {
    UNIQUE_DECODING,  // distance (1-ρ)/2, no conjecture, most queries
    PROVABLE_LIST,    // Johnson bound 1-√ρ-η with slack η = √ρ/20 (unconditional)
    REPAIRED_LIST,    // list-decoding-capacity bound H_q(δ) = 1-ρ (Crites-Stewart, eprint 2025/2046)
    CONJECTURED_LIST, // capacity 1-ρ (DISPROVED as a conjecture basis by eprint 2025/2046; kept for
                      // comparison with deployed FRI/STIR systems that still assume it)
};

/**
 * @brief Johnson slack η as the rational η = √ρ/`JOHNSON_SLACK_INV_NUMERATOR`, matching the
 * reference WHIR implementation's choice of η = √ρ/20.
 * @details The Johnson bound gives (1-√ρ-η, 1/(2η√ρ))-list-decodability only for η > 0: the list
 * size diverges as η → 0, so η = 0 is not an instantiation of the bound at all. Fixing η pins both
 * the per-query error √ρ+η and the list size, and hence the number of out-of-domain samples needed
 * to single out one list element.
 */
static constexpr uint64_t JOHNSON_SLACK_INV_NUMERATOR = 20;

/**
 * @brief Per-iteration schedule entry. `num_variables`/`log_domain_size` describe the oracle the
 * iteration starts from; `num_queries` is the number of in-domain queries made against it.
 * @details `folding_factor_bits` is the fold applied in this iteration, and hence also the leaf
 * arity of the oracle the iteration queries: a query must open the whole fold coset in one leaf.
 */
struct WhirRound {
    size_t num_variables;
    size_t log_domain_size;
    size_t log_inv_rate;
    size_t num_queries;
    size_t folding_factor_bits;
};

/**
 * @brief WHIR protocol parameters and the derived per-iteration schedule.
 * @details Immutable after `create`; both prover and verifier derive the identical schedule from the
 * same inputs, so none of it is sent through the transcript. See README.md §2 and §6.
 */
struct WhirConfig {
    size_t num_variables;       // committed arrays have 2^num_variables coefficients
    size_t security_bits;       // λ
    size_t log_inv_rate;        // r₀: initial codeword length is 2^{num_variables + r₀}
    size_t folding_factor_bits; // k: every iteration after the first folds 2^k
    size_t final_poly_bits;     // lower bound on the clear final polynomial's log-size
    WhirSoundness soundness;

    // k₀, the first iteration's fold, and hence the leaf arity of the committed oracles: a round-0
    // query reveals 2^{k₀} values *per committed column*, so on a wide commitment (transparent Honk
    // opens ~30 columns at once) this is the dominant term of the whole proof. Splitting it from k
    // buys a proportional cut there while leaving the later, single-column rounds free to fold
    // faster. The reference implementation calls it `initial_folding_factor`.
    size_t initial_folding_factor_bits;

    // Bits of proof of work the prover grinds out before each round's query indices are drawn.
    // Forging a round means re-grinding, so the work replaces query soundness one bit for one bit
    // and the query counts are derived at `security_bits - pow_bits`. The reference implementation
    // and ProveKit both take 10 bits this way; zero here means every bit comes from queries.
    size_t pow_bits = 0;

    // ---- Recursion profile ----------------------------------------------------------------
    // Three settings that leave the protocol's soundness argument untouched but change what an
    // in-circuit verifier has to do. None of them affect the schedule, so they can be flipped on a
    // config `create` already returned.

    // Grind with Poseidon2 rather than Blake3. Blake3 is the right choice natively — the prover runs
    // it 2^pow_bits times — but it costs tens of thousands of constraints per round in-circuit,
    // which wipes out everything grinding buys. Poseidon2 costs one permutation.
    bool poseidon2_pow = false;

    // Send each query its own authentication path instead of one batch per round. The batched form
    // is smaller on the wire, but its walk is over the *distinct* leaf set the queries induce, whose
    // shape is data-dependent; in-circuit that forces witness-indexed reads at every step. Per-query
    // paths make every index a compile-time constant. Openings are unhashed, so the extra digests
    // cost a recursive verifier nothing but proof bytes.
    bool per_query_openings = false;

    // Stop every authentication path `merkle_cap_levels` below the root and commit to that whole
    // level (a "Merkle cap") instead. See `MerkleTree::cap`. Only meaningful with per-query openings.
    size_t merkle_cap_levels = 0;

    /**
     * @brief The cap height this tree actually uses: right-sized for its own query count.
     * @details `merkle_cap_levels` is set from round 0, which makes the most queries against the
     * deepest tree. A later round queries a smaller oracle fewer times, where the same cap would
     * cost more to fold than the path levels it removes, so each tree takes the height optimal for
     * its own shape. Both parties derive it from the schedule, so nothing is transmitted.
     */
    size_t cap_levels_for(size_t num_queries, size_t depth) const
    {
        if (merkle_cap_levels == 0) {
            return 0;
        }
        return std::min(std::min(merkle_cap_levels, depth), best_cap_levels(num_queries, depth));
    }

    /** @brief Cap height minimizing `2^c - 1 + c*t` hashes for `t` queries against a depth-`d` tree. */
    static size_t best_cap_levels(size_t num_queries, size_t depth)
    {
        size_t best = 0;
        size_t best_cost = num_queries * depth;
        for (size_t c = 1; c <= depth; ++c) {
            const size_t cost = (size_t(1) << c) - 1 + num_queries * (depth - c);
            if (cost < best_cost) {
                best_cost = cost;
                best = c;
            }
        }
        return best;
    }

    /** @brief Turn on every recursion-profile setting, sizing the cap from the round-0 query count. */
    void enable_recursion_profile()
    {
        poseidon2_pow = true;
        per_query_openings = true;
        const WhirRound& first = rounds.empty() ? final_round : rounds[0];
        const size_t depth = first.log_domain_size - first.folding_factor_bits;
        merkle_cap_levels = best_cap_levels(first.num_queries, depth);
    }
    // Out-of-domain samples taken against every committed oracle, the initial one included
    // (README.md §4.2, §6). Zero in unique decoding, where the list holds a single codeword.
    size_t num_ood_samples = 1;

    // Column stacking (README.md §4.1): concatenate a commitment group's columns into one taller
    // array with narrow Merkle leaves, so a query costs 2^k values per group rather than 2^k per
    // column. Shrinks the proof several-fold but widens the batched oracle from one column to the
    // whole group, which costs proportionally more prover time. Off by default.
    bool stack_columns = false;

    // Zero-knowledge mode (README.md §8): committed arrays gain one variable over the payload
    // (blinding coefficients live in the high half), leaves are salted, and the opening batches in a
    // fresh random mask polynomial.
    bool zk = false;
    size_t num_payload_variables = 0;     // m: payload polynomials have 2^m coefficients
    size_t num_blinding_coefficients = 0; // q: random coefficients per commitment in zk mode

    std::vector<WhirRound> rounds; // the M fold-and-commit iterations
    WhirRound final_round;         // query schedule for the final clear-polynomial phase

    size_t num_iterations() const { return rounds.size(); }

    /** @brief Leaf arity of the oracle committed at the end of iteration `i` (the one iteration
     * i+1, or the final round, queries). */
    size_t committed_arity_bits(size_t i) const
    {
        return i + 1 < rounds.size() ? rounds[i + 1].folding_factor_bits : final_round.folding_factor_bits;
    }

    /**
     * @brief -log₂(x) of a Q192 fixed-point value x ∈ (0, 1), returned in Q64 fixed point.
     * @details Normalizes x into [1/2, 1) (accumulating the integer part), then extracts 64
     * fractional bits by repeated squaring: each squaring doubles the log, and a fall below 1/2
     * emits a set bit. All arithmetic is integer (96-bit-split truncating products, relative error
     * < 2⁻⁸⁹ over the whole extraction), so every platform derives the same value.
     */
    static uint256_t neg_log2_q192_to_q64(uint256_t x)
    {
        const uint256_t one_q192 = uint256_t(1) << 192;
        const uint256_t mask96 = (uint256_t(1) << 96) - 1;
        BB_ASSERT_GT(x, uint256_t(0), "-log2 of zero");
        uint256_t integer_part = 0;
        while (x < (one_q192 >> 1)) {
            x = x << 1;
            integer_part += 1;
        }
        uint256_t frac = 0;
        for (size_t bit = 0; bit < 64; ++bit) {
            const uint256_t hi = x >> 96;
            const uint256_t lo = x & mask96;
            x = hi * hi + ((hi * lo) >> 95); // (x/2^192)² in Q192; the lo·lo term truncates to 0
            frac = frac << 1;
            if (x < (one_q192 >> 1)) {
                frac = frac | 1;
                x = x << 1;
            }
        }
        return (integer_part << 64) + frac;
    }

    /**
     * @brief Number of queries for a code of log-inverse-rate `log_inv_rate` at `security_bits`.
     * @details CONJECTURED_LIST: t = ⌈λ/r⌉; PROVABLE_LIST: t = ⌈2λ/r⌉ (per-query soundness r/2 bits).
     * REPAIRED_LIST tests the largest distance δ* allowed by the repaired up-to-capacity conjectures
     * of Crites-Stewart (eprint 2025/2046), H_q(δ*) = 1-ρ: per-query error 1-δ* ≈ ρ + h₂(δ*)/log₂q,
     * solved by fixed-point iteration. UNIQUE_DECODING tests distance (1-ρ)/2, i.e. per-query error
     * (1+2⁻ʳ)/2: the count is the least t with ((1+2⁻ʳ)/2)ᵗ ≤ 2^{-λ}. Both fixed-point paths use
     * Q192/Q64 integer arithmetic so that both parties (and all platforms) derive the same integer.
     */
    static size_t compute_num_queries(size_t security_bits, size_t log_inv_rate, WhirSoundness soundness)
    {
        BB_ASSERT_GT(log_inv_rate, size_t(0), "WHIR requires rate < 1");
        switch (soundness) {
        case WhirSoundness::CONJECTURED_LIST:
            return (security_bits + log_inv_rate - 1) / log_inv_rate;
        case WhirSoundness::PROVABLE_LIST: {
            // Johnson: δ = 1-√ρ-η with η = √ρ/J, so the per-query error is √ρ·(1+1/J) and the bits
            // per query are r/2 - log₂(1+1/J). Taking η = 0 (bits per query exactly r/2) would not
            // be an instantiation of the Johnson bound, and understates the count by ~7% at J = 20.
            const uint256_t one_q192 = uint256_t(1) << 192;
            // -log₂(J/(J+1)) = log₂(1+1/J); the floored argument rounds the penalty up, i.e. toward
            // more queries.
            const uint256_t slack_penalty_q64 =
                neg_log2_q192_to_q64((one_q192 * JOHNSON_SLACK_INV_NUMERATOR) / (JOHNSON_SLACK_INV_NUMERATOR + 1));
            const uint256_t half_rate_q64 = uint256_t(log_inv_rate) << 63;
            BB_ASSERT_GT(half_rate_q64, slack_penalty_q64, "Johnson slack exhausts the per-query soundness");
            const uint256_t bits_per_query_q64 = half_rate_q64 - slack_penalty_q64;
            const uint256_t queries = ((uint256_t(security_bits) << 64) + bits_per_query_q64 - 1) / bits_per_query_q64;
            return static_cast<size_t>(queries.data[0]);
        }
        case WhirSoundness::REPAIRED_LIST: {
            // Solve δ* with H_q(δ*) = 1-ρ. For prime q, H_q(δ) = δ log_q(q-1) - δ log_q δ -
            // (1-δ) log_q(1-δ) ≤ δ + h₂(δ)/log₂q, so iterate δ ← 1-ρ - h₂(δ)/L from δ₀ = 1-ρ.
            // L = 253 rounds log₂|Fr| (≈ 253.5) down, overstating the entropy penalty; the h/L
            // division rounds up and a 2⁻³² guard is subtracted from the per-query bits, so every
            // rounding is toward more queries. The iteration is a contraction with factor
            // |h₂'(δ)|/L < 0.02, monotone from above; 25 rounds land within 2⁻¹⁹² of δ*.
            const uint256_t one_q192 = uint256_t(1) << 192;
            const uint256_t rho_q192 = one_q192 >> log_inv_rate;
            const uint256_t field_bits = 253;
            uint256_t delta = one_q192 - rho_q192;
            for (size_t iteration = 0; iteration < 25; ++iteration) {
                const uint256_t eps = one_q192 - delta;
                // h₂(δ) in Q192: δ·(-log₂δ) + (1-δ)·(-log₂(1-δ)), Q96 x Q96 products.
                const uint256_t entropy_q192 = (delta >> 96) * (neg_log2_q192_to_q64(delta) << 32) +
                                               (eps >> 96) * (neg_log2_q192_to_q64(eps) << 32);
                const uint256_t penalty_q192 = (entropy_q192 + field_bits - 1) / field_bits;
                BB_ASSERT_GT(one_q192 - rho_q192, penalty_q192, "repaired-regime distance vanished");
                delta = one_q192 - rho_q192 - penalty_q192;
            }
            const uint256_t bits_per_query_q64 = neg_log2_q192_to_q64(one_q192 - delta) - (uint256_t(1) << 32);
            const uint256_t queries = ((uint256_t(security_bits) << 64) + bits_per_query_q64 - 1) / bits_per_query_q64;
            return static_cast<size_t>(queries.data[0]);
        }
        case WhirSoundness::UNIQUE_DECODING: {
            // The accumulator represents the value (acc / 2^192)·2^{-deficit}; err = (1 + 2⁻ʳ)/2 in
            // Q192. Each step multiplies acc by err (three-term 96-bit-split product, truncation error
            // < 2⁻⁹⁰ relative) and renormalizes acc back into [2^191, 2^192], so `deficit` counts the
            // accumulated -log₂ of the running error probability.
            const uint256_t one_q192 = uint256_t(1) << 192;
            const uint256_t err_q192 = (one_q192 >> 1) + (one_q192 >> (1 + log_inv_rate));
            const uint256_t err_hi = err_q192 >> 96;
            const uint256_t err_lo = err_q192 & ((uint256_t(1) << 96) - 1);
            uint256_t acc = one_q192;
            size_t t = 0;
            size_t deficit = 0;
            while (deficit < security_bits) {
                acc = (acc >> 96) * err_hi + (((acc >> 96) * err_lo) >> 96) +
                      (((acc & ((uint256_t(1) << 96) - 1)) * err_hi) >> 96);
                while (acc < (one_q192 >> 1)) {
                    acc = acc << 1;
                    deficit += 1;
                }
                ++t;
                BB_ASSERT_LT(t, size_t(100000), "UNIQUE_DECODING query count diverged");
            }
            return t;
        }
        }
        return 0; // unreachable
    }

    /**
     * @brief Out-of-domain samples needed against an oracle of `num_variables` variables at rate
     * 2^{-log_inv_rate}, so that the sampled point singles out one codeword of the decoding list.
     * @details STIR lemma 4.5: with a list of size L, the probability that two list elements agree
     * on all s random out-of-domain points is at most (L choose 2)·((n-1)/|F|)^s, so
     * s ≥ (λ + log₂(L choose 2)) / (log₂|F| - log₂(n-1)). Every quantity is replaced by an integer
     * bound in the direction of more samples: log₂|F| ≥ 253 for BN254 Fr, log₂(n-1) ≤ num_variables,
     * and log₂(L choose 2) ≤ 2·log₂L. In unique decoding the list is a single codeword and no sample
     * is required. At BN254's field size this returns 1 across the whole schedule.
     */
    static size_t compute_num_ood_samples(size_t security_bits,
                                          size_t log_inv_rate,
                                          size_t num_variables,
                                          WhirSoundness soundness)
    {
        if (soundness == WhirSoundness::UNIQUE_DECODING) {
            return 0;
        }
        // log₂ of an upper bound on the list size. Johnson gives exactly 1/(2η√ρ) = J·2^{r-1};
        // the capacity regimes have no proven bound, so charge the generous n/ρ.
        size_t log_slack_inv = 0;
        while ((uint64_t(1) << log_slack_inv) < JOHNSON_SLACK_INV_NUMERATOR) {
            ++log_slack_inv;
        }
        const size_t log_list_bound = soundness == WhirSoundness::PROVABLE_LIST ? (log_inv_rate - 1) + log_slack_inv
                                                                                : num_variables + log_inv_rate;
        const size_t field_bits = 253;
        BB_ASSERT_GT(field_bits, num_variables, "polynomial too large to isolate a list element");
        const size_t bits_per_sample = field_bits - num_variables;
        return (security_bits + 2 * log_list_bound + bits_per_sample - 1) / bits_per_sample;
    }

    /**
     * @brief Build the full parameter schedule. Iterations continue while at least `final_poly_bits`
     * variables would remain after folding; the rate improves by k-1 bits per iteration.
     * @param num_payload_variables log-size of the polynomials being opened; in zk mode the
     * committed arrays have one variable more (the blinded high half).
     * @param max_stack_bits column-stacking headroom: when non-zero, groups are stacked and a group
     * of up to 2^max_stack_bits columns fits, so the protocol arrays carry this many variables over
     * the payload. Zero (the default) keeps the one-codeword-per-column layout.
     * @param initial_folding_factor_bits k₀, the first iteration's fold; zero means "same as k".
     */
    static WhirConfig create(size_t num_payload_variables,
                             size_t security_bits = 100,
                             size_t log_inv_rate = 2,
                             size_t folding_factor_bits = 4,
                             size_t final_poly_bits = 4,
                             WhirSoundness soundness = WhirSoundness::CONJECTURED_LIST,
                             bool zk = false,
                             size_t max_stack_bits = 0,
                             size_t initial_folding_factor_bits = 0,
                             size_t pow_bits = 0)
    {
        const size_t num_variables = num_payload_variables + max_stack_bits + (zk ? 1 : 0);
        const size_t k0 = initial_folding_factor_bits == 0 ? folding_factor_bits : initial_folding_factor_bits;
        // Grinding supplies `pow_bits` of each round's soundness; the queries supply the rest.
        BB_ASSERT_LT(pow_bits, security_bits, "proof of work cannot cover the whole security level");
        const size_t query_bits = security_bits - pow_bits;
        BB_ASSERT_GT(folding_factor_bits, size_t(0));
        BB_ASSERT_GTE(num_variables, std::max(folding_factor_bits, k0), "polynomial too small for one fold");
        // BN254 Fr has 2-adicity 28; the initial codeword domain must be a power-of-two subgroup.
        BB_ASSERT_LTE(num_variables + log_inv_rate, size_t(28), "initial domain exceeds field 2-adicity");

        WhirConfig config{ .num_variables = num_variables,
                           .security_bits = security_bits,
                           .log_inv_rate = log_inv_rate,
                           .folding_factor_bits = folding_factor_bits,
                           .final_poly_bits = final_poly_bits,
                           .soundness = soundness,
                           .initial_folding_factor_bits = k0,
                           .pow_bits = pow_bits,
                           .stack_columns = max_stack_bits > 0,
                           .zk = zk,
                           .num_payload_variables = num_payload_variables,
                           .num_blinding_coefficients = 0,
                           .rounds = {},
                           .final_round = {} };

        size_t vars = num_variables;
        size_t log_domain = num_variables + log_inv_rate;
        for (size_t i = 0;; ++i) {
            const size_t k = i == 0 ? k0 : folding_factor_bits;
            if (vars < k || (vars - k) < final_poly_bits) {
                break;
            }
            const size_t rate_bits = log_domain - vars;
            config.rounds.push_back({ .num_variables = vars,
                                      .log_domain_size = log_domain,
                                      .log_inv_rate = rate_bits,
                                      .num_queries = compute_num_queries(query_bits, rate_bits, soundness),
                                      .folding_factor_bits = k });
            vars -= k;
            log_domain -= 1;
            // Leaves of the next oracle group 2^k values; its domain must be large enough.
            BB_ASSERT_GTE(log_domain, folding_factor_bits, "domain too small for leaf grouping");
        }
        const size_t final_rate_bits = log_domain - vars;
        // The final round queries the last committed oracle, whose leaves already group 2^k values -
        // except when no iteration ran at all, where the round-0 commitment's k₀ leaves apply.
        config.final_round = { .num_variables = vars,
                               .log_domain_size = log_domain,
                               .log_inv_rate = final_rate_bits,
                               .num_queries = compute_num_queries(query_bits, final_rate_bits, soundness),
                               .folding_factor_bits = config.rounds.empty() ? k0 : folding_factor_bits };

        // One sample count covers every commitment: the initial oracle and each iteration's folded
        // oracle, taking the worst case over the schedule.
        config.num_ood_samples = compute_num_ood_samples(security_bits, log_inv_rate, num_variables, soundness);
        for (const WhirRound& round : config.rounds) {
            const size_t folded_variables = round.num_variables - round.folding_factor_bits;
            const size_t folded_rate_bits = round.log_inv_rate + round.folding_factor_bits - 1;
            config.num_ood_samples =
                std::max(config.num_ood_samples,
                         compute_num_ood_samples(security_bits, folded_rate_bits, folded_variables, soundness));
        }

        if (zk) {
            // Round-0 queries are the only openings of per-polynomial leaves; q of them must remain
            // information-theoretically blinded (README.md §8), with a small slack margin.
            const size_t round0_queries =
                config.rounds.empty() ? config.final_round.num_queries : config.rounds[0].num_queries;
            config.num_blinding_coefficients = round0_queries + 8;
            // Blinding lives above the payload and the shift contract's zero slot, so it never
            // collides with either. The room there is 2^num_payload_variables in both layouts: an
            // interleaved column is 2^num_variables wide over a 2^m payload, and the narrowest
            // stacked group (one column) is 2^{m+1} wide over the same payload. Sizing this against
            // the committed width instead would overrun exactly those narrow stacked groups.
            BB_ASSERT_LT(config.num_blinding_coefficients + 1,
                         size_t(1) << num_payload_variables,
                         "blinding coefficients do not fit above the payload");
        }
        return config;
    }
};

} // namespace bb::whir
