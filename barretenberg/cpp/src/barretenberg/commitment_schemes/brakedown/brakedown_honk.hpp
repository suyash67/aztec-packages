#pragma once

#include "barretenberg/commitment_schemes/brakedown/brakedown_code.hpp"
#include "barretenberg/commitment_schemes/ligero/ligero.hpp"
#include "barretenberg/commitment_schemes/transparent_honk.hpp"

#include <cmath>
#include <memory>

namespace bb::brakedown {

using ligero::LigeroColumnRef;
using ligero::LigeroConfig;

/**
 * @brief Ligero's tensor PCS over the Brakedown code instead of Reed-Solomon.
 *
 * @details The tensor protocol is code-agnostic: it needs an encoder and a query count, both of
 * which this policy supplies (see `ligero/ligero.hpp`). Two things differ from the RS policy and
 * both come straight from the code's properties.
 *
 * **Query count.** RS at rate 1/4 is used under the up-to-capacity conjecture, giving
 * `t = ceil(λ/r)` — 50 queries at λ = 100. Brakedown has no such conjecture available and a far
 * smaller distance, so the query count follows the *provable* interleaved proximity test: per-query
 * soundness `1 - δ/3`, hence `t = ceil(λ / -log2(1 - δ/3))`. At δ = 0.07 that is about 2936 queries.
 * This ~59x is the whole story of this backend, and it is why Brakedown trades proof size for
 * encoder speed rather than winning outright.
 *
 * **Codeword length.** Brakedown codewords are not powers of two (the recursion ceilings at every
 * level), so the Merkle tree is built over the padded length while queries are drawn below the true
 * length — padding is never opened and never enters soundness.
 */
struct BrakedownCodePolicy {
    explicit BrakedownCodePolicy(const LigeroConfig& config)
        : params(BrakedownParams::distance_7pct())
        , code(std::make_shared<BrakedownCode>(config.num_cols(), params))
    {}

    size_t message_length() const { return code->message_length(); }
    size_t codeword_length() const { return code->codeword_length(); }
    size_t padded_codeword_length() const { return ligero::detail::next_power_of_two(code->codeword_length()); }

    std::vector<fr> encode(std::span<const fr> row) const { return code->encode(row); }

    /** @brief The provable interleaved test: t = ceil(λ / -log2(1 - δ/3)). */
    static size_t num_queries(const LigeroConfig& config)
    {
        return provable_num_queries(config.security_bits, BrakedownParams::distance_7pct().distance());
    }

    static size_t provable_num_queries(size_t security_bits, double distance)
    {
        const double per_query_bits = -std::log2(1.0 - (distance / 3.0));
        return static_cast<size_t>(std::ceil(static_cast<double>(security_bits) / per_query_bits));
    }

    BrakedownParams params;
    // Shared so the policy stays copyable: the sampled matrices are large and immutable.
    std::shared_ptr<BrakedownCode> code;
};

/**
 * @brief Brakedown parameters: Ligero's shape, but sized for a query count two orders of magnitude
 * larger.
 * @details Ligero balances the leaf width against the transmitted rows via `default_log_num_cols`,
 * which assumes ~50 queries. With ~2936 the balance moves hard toward wider rows (fewer, narrower
 * leaves), so the column count is recomputed with the real query count rather than inherited.
 */
inline LigeroConfig make_brakedown_config(size_t num_variables, size_t security_bits, size_t est_polynomials = 36)
{
    const size_t queries =
        BrakedownCodePolicy::provable_num_queries(security_bits, BrakedownParams::distance_7pct().distance());
    size_t log_num_cols = LigeroConfig::default_log_num_cols(num_variables, est_polynomials, queries);
    log_num_cols = std::min(log_num_cols, num_variables - 1);
    // log_inv_rate is unused by this code (the rate is 1/r, fixed by the parameters) but the shared
    // config carries it; keep it at Ligero's value so the struct stays comparable.
    return { num_variables, log_num_cols, 2, security_bits, queries };
}

/** @brief PCS backend adapter binding Brakedown-over-Ligero into the UltraHonk shell. */
template <typename Hasher_> struct BrakedownPcs {
    using Hasher = Hasher_;
    using Config = LigeroConfig;
    using CommitmentKey = ligero::LigeroCommitmentKey<Hasher, BrakedownCodePolicy>;
    using GroupData = ligero::LigeroGroupData<Hasher>;
    using Digest = typename Hasher::Digest;
    using ProverClaims = typename ligero::LigeroProver<Hasher>::Claims;
    using VerifierClaims = typename ligero::LigeroVerifier<Hasher>::Claims;
    using GroupCommitment = typename Hasher::Digest;

    static Config make_config(size_t log_dyadic_size, size_t security_bits, size_t /*log_inv_rate*/)
    {
        return make_brakedown_config(log_dyadic_size, security_bits);
    }
    static size_t payload_variables(const Config& config) { return config.num_variables; }

    static GroupCommitment group_commitment(const GroupData& data) { return data.tree.root(); }
    template <typename Transcript>
    static void send_group_commitment(const std::shared_ptr<Transcript>& transcript,
                                      const std::string& label,
                                      const GroupData& data)
    {
        transcript->send_to_verifier(label, Hasher::digest_to_fields(data.tree.root()));
    }
    template <typename Transcript>
    static GroupCommitment receive_group_commitment(const std::shared_ptr<Transcript>& transcript,
                                                    const std::string& label,
                                                    size_t /*num_columns*/,
                                                    const Config& /*config*/)
    {
        return Hasher::digest_from_fields(
            transcript->template receive_from_prover<std::array<fr, Hasher::DIGEST_NUM_FIELDS>>(label));
    }
    template <typename Transcript>
    static void absorb_group_commitment(const std::shared_ptr<Transcript>& transcript,
                                        const std::string& label,
                                        const GroupCommitment& commitment)
    {
        transcript->add_to_hash_buffer(label, Hasher::digest_to_fields(commitment));
    }
    static void set_group_commitments(VerifierClaims& claims, std::vector<GroupCommitment> commitments)
    {
        claims.group_roots = std::move(commitments);
    }

    template <typename Transcript>
    static void prove_opening(const CommitmentKey& ck,
                              const ProverClaims& claims,
                              std::span<const fr> u,
                              const std::shared_ptr<Transcript>& transcript)
    {
        ligero::LigeroProver<Hasher>::prove(ck, claims, u, transcript);
    }
    template <typename Transcript>
    static bool verify_opening(const Config& config,
                               const VerifierClaims& claims,
                               std::span<const fr> u,
                               const std::shared_ptr<Transcript>& transcript)
    {
        const BrakedownCodePolicy code(config);
        return ligero::LigeroVerifier<Hasher>::verify(config, claims, u, transcript, code);
    }
};

/** @brief UltraHonk with Ligero's tensor PCS over the Brakedown linear-time code. */
template <typename Hasher = whir::Blake3sMerkleHasher>
using BrakedownHonk = honk_transparent::TransparentHonk<BrakedownPcs<Hasher>>;

} // namespace bb::brakedown
