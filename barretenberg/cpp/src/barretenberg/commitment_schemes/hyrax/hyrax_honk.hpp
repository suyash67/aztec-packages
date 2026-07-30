#pragma once

#include "barretenberg/commitment_schemes/hyrax/hyrax.hpp"
#include "barretenberg/commitment_schemes/transparent_honk.hpp"

namespace bb::hyrax {

/** @brief PCS backend adapter binding Hyrax into the UltraHonk shell. */
struct HyraxPcs {
    using Config = HyraxConfig;
    using CommitmentKey = HyraxCommitmentKey;
    using GroupData = HyraxGroupData;
    using ProverClaims = HyraxProver::Claims;
    using VerifierClaims = HyraxVerifier::Claims;
    // Per column, the vector of Pedersen row commitments.
    using GroupCommitment = std::vector<std::vector<Commitment>>;

    static Config make_config(size_t log_dyadic_size, size_t security_bits, size_t log_inv_rate)
    {
        return HyraxConfig::create(log_dyadic_size, security_bits, log_inv_rate);
    }
    static size_t payload_variables(const Config& config) { return config.num_variables; }

    static GroupCommitment group_commitment(const GroupData& data) { return data.commitments; }
    template <typename Transcript>
    static void send_group_commitment(const std::shared_ptr<Transcript>& transcript,
                                      const std::string& label,
                                      const GroupData& data)
    {
        for (size_t c = 0; c < data.commitments.size(); ++c) {
            for (size_t r = 0; r < data.commitments[c].size(); ++r) {
                transcript->send_to_verifier(label + "_" + std::to_string(c) + "_" + std::to_string(r),
                                             data.commitments[c][r]);
            }
        }
    }
    template <typename Transcript>
    static GroupCommitment receive_group_commitment(const std::shared_ptr<Transcript>& transcript,
                                                    const std::string& label,
                                                    size_t num_columns,
                                                    const Config& config)
    {
        GroupCommitment commitments(num_columns);
        for (size_t c = 0; c < num_columns; ++c) {
            for (size_t r = 0; r < config.num_rows(); ++r) {
                commitments[c].push_back(transcript->template receive_from_prover<Commitment>(
                    label + "_" + std::to_string(c) + "_" + std::to_string(r)));
            }
        }
        return commitments;
    }
    template <typename Transcript>
    static void absorb_group_commitment(const std::shared_ptr<Transcript>& transcript,
                                        const std::string& label,
                                        const GroupCommitment& commitments)
    {
        for (const auto& column : commitments) {
            for (const Commitment& commitment : column) {
                transcript->add_to_hash_buffer(label, commitment);
            }
        }
    }
    static void set_group_commitments(VerifierClaims& claims, std::vector<GroupCommitment> commitments)
    {
        claims.group_commitments = std::move(commitments);
    }

    template <typename Transcript>
    static void prove_opening(const CommitmentKey& ck,
                              const ProverClaims& claims,
                              std::span<const fr> u,
                              const std::shared_ptr<Transcript>& transcript)
    {
        HyraxProver::prove(ck, claims, u, transcript);
    }
    template <typename Transcript>
    static bool verify_opening(const Config& config,
                               const VerifierClaims& claims,
                               std::span<const fr> u,
                               const std::shared_ptr<Transcript>& transcript)
    {
        const HyraxCommitmentKey ck(config);
        return HyraxVerifier::verify(config, claims, u, transcript, ck);
    }
};

/** @brief UltraHonk with Hyrax as the polynomial commitment scheme. */
using HyraxHonk = honk_transparent::TransparentHonk<HyraxPcs>;

} // namespace bb::hyrax
