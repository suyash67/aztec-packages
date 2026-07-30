#pragma once

#include "barretenberg/commitment_schemes/dory/dory.hpp"
#include "barretenberg/commitment_schemes/transparent_honk.hpp"

namespace bb::dory {

/** @brief PCS backend adapter binding the Dory-style two-tier PCS into the UltraHonk shell. */
struct DoryPcs {
    using Config = DoryConfig;
    using CommitmentKey = DoryCommitmentKey;
    using GroupData = DoryGroupData;
    using ProverClaims = DoryProver::Claims;
    using VerifierClaims = DoryVerifier::Claims;
    // One GT element per column.
    using GroupCommitment = std::vector<fq12>;

    static Config make_config(size_t log_dyadic_size, size_t security_bits, size_t log_inv_rate)
    {
        return DoryConfig::create(log_dyadic_size, security_bits, log_inv_rate);
    }
    static size_t payload_variables(const Config& config) { return config.num_variables; }

    static GroupCommitment group_commitment(const GroupData& data) { return data.commitments; }
    template <typename Transcript>
    static void send_group_commitment(const std::shared_ptr<Transcript>& transcript,
                                      const std::string& label,
                                      const GroupData& data)
    {
        for (size_t c = 0; c < data.commitments.size(); ++c) {
            transcript->send_to_verifier(label + "_" + std::to_string(c), detail::fq12_to_fields(data.commitments[c]));
        }
    }
    template <typename Transcript>
    static GroupCommitment receive_group_commitment(const std::shared_ptr<Transcript>& transcript,
                                                    const std::string& label,
                                                    size_t num_columns,
                                                    const Config& /*config*/)
    {
        GroupCommitment commitments;
        for (size_t c = 0; c < num_columns; ++c) {
            commitments.push_back(detail::fq12_from_fields(
                transcript->template receive_from_prover<std::array<fr, 24>>(label + "_" + std::to_string(c))));
        }
        return commitments;
    }
    template <typename Transcript>
    static void absorb_group_commitment(const std::shared_ptr<Transcript>& transcript,
                                        const std::string& label,
                                        const GroupCommitment& commitments)
    {
        for (const fq12& commitment : commitments) {
            transcript->add_to_hash_buffer(label, detail::fq12_to_fields(commitment));
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
        DoryProver::prove(ck, claims, u, transcript);
    }
    template <typename Transcript>
    static bool verify_opening(const Config& config,
                               const VerifierClaims& claims,
                               std::span<const fr> u,
                               const std::shared_ptr<Transcript>& transcript)
    {
        const DoryCommitmentKey ck(config);
        return DoryVerifier::verify(config, claims, u, transcript, ck);
    }
};

/** @brief UltraHonk with the Dory-style PCS as the polynomial commitment scheme. */
using DoryHonk = honk_transparent::TransparentHonk<DoryPcs>;

} // namespace bb::dory
