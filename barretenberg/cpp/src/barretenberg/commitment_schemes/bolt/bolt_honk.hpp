#pragma once

#include "barretenberg/commitment_schemes/bolt/bolt_code.hpp"
#include "barretenberg/commitment_schemes/ligero/ligero.hpp"
#include "barretenberg/commitment_schemes/transparent_honk.hpp"

#include <cmath>
#include <memory>

namespace bb::bolt {

using ligero::LigeroConfig;
using ligero::QuerySegment;

/**
 * @brief Ligero's tensor PCS over Bolt's sketched code, tested piecewise.
 *
 * @details The sketched code `(x, C(Hx))` carries two *separate* distance guarantees rather than one
 * diluted average (paper Claim 3.1): two distinct codewords differ in more than `gamma` of the
 * systematic stretch, **or** more than `delta` of the sketch stretch. So the two stretches are
 * queried independently, each with the count its own distance warrants:
 *
 *     systematic:  t1 = ceil(lambda / -log2(1 - gamma))
 *     sketch:      t2 = ceil(lambda / log_inv_rate)          (RS, capacity conjecture)
 *
 * A cheating prover has to survive both, so the error is the max of the two, and each is driven by a
 * distance far better than the code's overall one. Sampling uniformly against the diluted distance
 * `min(gamma, delta*alpha/rho) / (1 + alpha/rho)` instead would cost several times more queries —
 * `QueryPlanIsPiecewiseAndCheaperThanDiluted` measures the gap. This piecewise test is the reason
 * Bolt beats Brakedown on proof size despite a more expensive encoder.
 */
struct BoltCodePolicy {
    explicit BoltCodePolicy(const LigeroConfig& config)
        : params(BoltParams::bn254())
        , code(std::make_shared<BoltCode>(config.num_cols(), params, config.log_inv_rate))
    {}

    size_t message_length() const { return code->message_length(); }
    size_t codeword_length() const { return code->codeword_length(); }
    size_t padded_codeword_length() const { return ligero::detail::next_power_of_two(code->codeword_length()); }

    std::vector<fr> encode(std::span<const fr> row) const { return code->encode(row); }

    static size_t systematic_queries(size_t security_bits, double gamma)
    {
        return static_cast<size_t>(std::ceil(static_cast<double>(security_bits) / -std::log2(1.0 - gamma)));
    }

    static size_t sketch_queries(const LigeroConfig& config)
    {
        return (config.security_bits + config.log_inv_rate - 1) / config.log_inv_rate;
    }

    /** @brief What a single-distance test would cost, for comparison. */
    static size_t diluted_queries(const LigeroConfig& config, const BoltParams& params)
    {
        const double rho = std::pow(2.0, -static_cast<double>(config.log_inv_rate));
        const double sketch_fraction = params.alpha() / rho;
        const double diluted = std::min(params.gamma, (1.0 - rho) * sketch_fraction) / (1.0 + sketch_fraction);
        return static_cast<size_t>(
            std::ceil(static_cast<double>(config.security_bits) / -std::log2(1.0 - (diluted / 3.0))));
    }

    static size_t num_queries(const LigeroConfig& config)
    {
        return systematic_queries(config.security_bits, BoltParams::bn254().gamma) + sketch_queries(config);
    }

    /** @brief Two segments: the systematic prefix and the encoded sketch. */
    std::vector<QuerySegment> query_plan(const LigeroConfig& config) const
    {
        return { { 0, code->message_length(), systematic_queries(config.security_bits, params.gamma) },
                 { code->sketch_begin(), code->sketch_codeword_length(), sketch_queries(config) } };
    }

    BoltParams params;
    // Shared so the policy stays copyable: the sampled sketch matrix is large and immutable.
    std::shared_ptr<BoltCode> code;
};

/** @brief Bolt parameters: Ligero's shape re-optimized for the piecewise query count. */
inline LigeroConfig make_bolt_config(size_t num_variables, size_t security_bits, size_t est_polynomials = 36)
{
    const LigeroConfig probe{ num_variables, 1, 2, security_bits, 0 };
    const size_t queries = BoltCodePolicy::num_queries(probe);
    size_t log_num_cols = LigeroConfig::default_log_num_cols(num_variables, est_polynomials, queries);
    log_num_cols = std::min(log_num_cols, num_variables - 1);
    return { num_variables, log_num_cols, 2, security_bits, queries };
}

/** @brief PCS backend adapter binding Bolt-over-Ligero into the UltraHonk shell. */
template <typename Hasher_> struct BoltPcs {
    using Hasher = Hasher_;
    using Config = LigeroConfig;
    using CommitmentKey = ligero::LigeroCommitmentKey<Hasher, BoltCodePolicy>;
    using GroupData = ligero::LigeroGroupData<Hasher>;
    using Digest = typename Hasher::Digest;
    using ProverClaims = typename ligero::LigeroProver<Hasher>::Claims;
    using VerifierClaims = typename ligero::LigeroVerifier<Hasher>::Claims;
    using GroupCommitment = typename Hasher::Digest;

    static Config make_config(size_t log_dyadic_size, size_t security_bits, size_t /*log_inv_rate*/)
    {
        return make_bolt_config(log_dyadic_size, security_bits);
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
        const BoltCodePolicy code(config);
        return ligero::LigeroVerifier<Hasher>::verify(config, claims, u, transcript, code);
    }
};

/** @brief UltraHonk with Ligero's tensor PCS over Bolt's sketched code. */
template <typename Hasher = whir::Blake3sMerkleHasher>
using BoltHonk = honk_transparent::TransparentHonk<BoltPcs<Hasher>>;

} // namespace bb::bolt
