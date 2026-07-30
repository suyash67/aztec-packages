#pragma once

#include "barretenberg/common/assert.hpp"
#include "barretenberg/ecc/curves/bn254/fr.hpp"
#include "barretenberg/numeric/uint256/uint256.hpp"

#include <cstddef>
#include <vector>

namespace bb::whir {

/**
 * @brief Proximity regime assumed for the underlying Reed-Solomon codes; determines the per-round
 * query counts. See README.md §6.
 */
enum class WhirSoundness {
    UNIQUE_DECODING,  // distance (1-ρ)/2, no conjecture, most queries
    PROVABLE_LIST,    // Johnson bound 1-√ρ
    CONJECTURED_LIST, // capacity 1-ρ (conjectured; standard in deployed FRI/STIR systems)
};

/**
 * @brief Per-iteration schedule entry. `num_variables`/`log_domain_size` describe the oracle the
 * iteration starts from; `num_queries` is the number of in-domain queries made against it.
 */
struct WhirRound {
    size_t num_variables;
    size_t log_domain_size;
    size_t log_inv_rate;
    size_t num_queries;
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
    size_t folding_factor_bits; // k: each iteration folds 2^k
    size_t final_poly_bits;     // lower bound on the clear final polynomial's log-size
    WhirSoundness soundness;

    // Zero-knowledge mode (README.md §8): committed arrays gain one variable over the payload
    // (blinding coefficients live in the high half), leaves are salted, and the opening batches in a
    // fresh random mask polynomial.
    bool zk = false;
    size_t num_payload_variables = 0;     // m: payload polynomials have 2^m coefficients
    size_t num_blinding_coefficients = 0; // q: random coefficients per commitment in zk mode

    std::vector<WhirRound> rounds; // the M fold-and-commit iterations
    WhirRound final_round;         // query schedule for the final clear-polynomial phase

    size_t folding_factor() const { return size_t(1) << folding_factor_bits; }
    size_t num_iterations() const { return rounds.size(); }

    /**
     * @brief Number of queries for a code of log-inverse-rate `log_inv_rate` at `security_bits`.
     * @details CONJECTURED_LIST: t = ⌈λ/r⌉; PROVABLE_LIST: t = ⌈2λ/r⌉ (per-query soundness r/2 bits).
     * UNIQUE_DECODING tests distance (1-ρ)/2, i.e. per-query error (1+2⁻ʳ)/2: the count is the least t
     * with ((1+2⁻ʳ)/2)ᵗ ≤ 2^{-λ}, computed in Q192 fixed point so that both parties (and all
     * platforms) derive the same integer.
     */
    static size_t compute_num_queries(size_t security_bits, size_t log_inv_rate, WhirSoundness soundness)
    {
        BB_ASSERT_GT(log_inv_rate, size_t(0), "WHIR requires rate < 1");
        switch (soundness) {
        case WhirSoundness::CONJECTURED_LIST:
            return (security_bits + log_inv_rate - 1) / log_inv_rate;
        case WhirSoundness::PROVABLE_LIST:
            return (2 * security_bits + log_inv_rate - 1) / log_inv_rate;
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
     * @brief Build the full parameter schedule. Iterations continue while at least `final_poly_bits`
     * variables would remain after folding; the rate improves by k-1 bits per iteration.
     * @param num_payload_variables log-size of the polynomials being opened; in zk mode the
     * committed arrays have one variable more (the blinded high half).
     */
    static WhirConfig create(size_t num_payload_variables,
                             size_t security_bits = 100,
                             size_t log_inv_rate = 2,
                             size_t folding_factor_bits = 4,
                             size_t final_poly_bits = 4,
                             WhirSoundness soundness = WhirSoundness::CONJECTURED_LIST,
                             bool zk = false)
    {
        const size_t num_variables = num_payload_variables + (zk ? 1 : 0);
        BB_ASSERT_GT(folding_factor_bits, size_t(0));
        BB_ASSERT_GTE(num_variables, folding_factor_bits, "polynomial too small for one fold");
        // BN254 Fr has 2-adicity 28; the initial codeword domain must be a power-of-two subgroup.
        BB_ASSERT_LTE(num_variables + log_inv_rate, size_t(28), "initial domain exceeds field 2-adicity");

        WhirConfig config{ .num_variables = num_variables,
                           .security_bits = security_bits,
                           .log_inv_rate = log_inv_rate,
                           .folding_factor_bits = folding_factor_bits,
                           .final_poly_bits = final_poly_bits,
                           .soundness = soundness,
                           .zk = zk,
                           .num_payload_variables = num_payload_variables,
                           .num_blinding_coefficients = 0,
                           .rounds = {},
                           .final_round = {} };

        size_t vars = num_variables;
        size_t log_domain = num_variables + log_inv_rate;
        while (vars >= folding_factor_bits && (vars - folding_factor_bits) >= final_poly_bits) {
            const size_t rate_bits = log_domain - vars;
            config.rounds.push_back({ .num_variables = vars,
                                      .log_domain_size = log_domain,
                                      .log_inv_rate = rate_bits,
                                      .num_queries = compute_num_queries(security_bits, rate_bits, soundness) });
            vars -= folding_factor_bits;
            log_domain -= 1;
            // Leaves of the next oracle group 2^k values; its domain must be large enough.
            BB_ASSERT_GTE(log_domain, folding_factor_bits, "domain too small for leaf grouping");
        }
        const size_t final_rate_bits = log_domain - vars;
        config.final_round = { .num_variables = vars,
                               .log_domain_size = log_domain,
                               .log_inv_rate = final_rate_bits,
                               .num_queries = compute_num_queries(security_bits, final_rate_bits, soundness) };

        if (zk) {
            // Round-0 queries are the only openings of per-polynomial leaves; q of them must remain
            // information-theoretically blinded (README.md §8), with a small slack margin.
            const size_t round0_queries =
                config.rounds.empty() ? config.final_round.num_queries : config.rounds[0].num_queries;
            config.num_blinding_coefficients = round0_queries + 8;
            // Blinding lives in [2^m + 1, 2^{m+1}) so it never collides with the payload or the
            // shift contract's zero slot at 2^m.
            BB_ASSERT_LT(config.num_blinding_coefficients + 1,
                         size_t(1) << num_payload_variables,
                         "blinding coefficients do not fit above the payload");
        }
        return config;
    }
};

} // namespace bb::whir
