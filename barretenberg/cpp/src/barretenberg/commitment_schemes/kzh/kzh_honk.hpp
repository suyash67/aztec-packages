#pragma once

#include "barretenberg/commitment_schemes/kzh/kzh.hpp"
#include "barretenberg/commitment_schemes/transparent_honk.hpp"

namespace bb::kzh {

/** @brief PCS backend adapter binding KZH2 into the UltraHonk shell. */
struct KzhPcs {
    using Config = KzhConfig;
    using CommitmentKey = KzhCommitmentKey;
    using GroupData = KzhGroupData;
    using ProverClaims = KzhProver::Claims;
    using VerifierClaims = KzhVerifier::Claims;
    using GroupCommitment = std::vector<Commitment>;

    static Config make_config(size_t log_dyadic_size, size_t security_bits, size_t log_inv_rate)
    {
        return KzhConfig::create(log_dyadic_size, security_bits, log_inv_rate);
    }
    static size_t payload_variables(const Config& config) { return config.num_variables; }

    static GroupCommitment group_commitment(const GroupData& data) { return data.commitments; }
    template <typename Transcript>
    static void send_group_commitment(const std::shared_ptr<Transcript>& transcript,
                                      const std::string& label,
                                      const GroupData& data)
    {
        for (size_t c = 0; c < data.commitments.size(); ++c) {
            transcript->send_to_verifier(label + "_" + std::to_string(c), data.commitments[c]);
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
            commitments.push_back(
                transcript->template receive_from_prover<Commitment>(label + "_" + std::to_string(c)));
        }
        return commitments;
    }
    template <typename Transcript>
    static void absorb_group_commitment(const std::shared_ptr<Transcript>& transcript,
                                        const std::string& label,
                                        const GroupCommitment& commitments)
    {
        for (const Commitment& commitment : commitments) {
            transcript->add_to_hash_buffer(label, commitment);
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
        KzhProver::prove(ck, claims, u, transcript);
    }
    template <typename Transcript>
    static bool verify_opening(const Config& config,
                               const VerifierClaims& claims,
                               std::span<const fr> u,
                               const std::shared_ptr<Transcript>& transcript)
    {
        // Deterministic test-only setup: the verifier derives the same key (see KzhCommitmentKey).
        const KzhCommitmentKey ck(config);
        return KzhVerifier::verify(config, claims, u, transcript, ck);
    }
};

/** @brief UltraHonk with KZH2 as the polynomial commitment scheme. */
using KzhHonk = honk_transparent::TransparentHonk<KzhPcs>;

} // namespace bb::kzh
