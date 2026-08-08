#pragma once

#include "barretenberg/commitment_schemes/fflonk/fflonk.hpp"
#include "barretenberg/commitment_schemes/transparent_honk.hpp"

namespace bb::fflonk {

/** @brief PCS backend adapter binding the fflonk opening into the UltraHonk shell. */
struct FflonkPcs {
    using Config = FflonkConfig;
    using CommitmentKey = FflonkCommitmentKey;
    using GroupData = FflonkGroupData;
    using ProverClaims = FflonkProver::Claims;
    using VerifierClaims = FflonkVerifier::Claims;
    // One packed commitment per commitment round, whatever the round's column count.
    using GroupCommitment = Commitment;

    static Config make_config(size_t log_dyadic_size, size_t security_bits, size_t log_inv_rate)
    {
        return FflonkConfig::create(log_dyadic_size, security_bits, log_inv_rate);
    }
    static size_t payload_variables(const Config& config) { return config.num_variables; }

    static GroupCommitment group_commitment(const GroupData& data) { return data.commitment; }

    template <typename Transcript>
    static void send_group_commitment(const std::shared_ptr<Transcript>& transcript,
                                      const std::string& label,
                                      const GroupData& data)
    {
        transcript->send_to_verifier(label, data.commitment);
    }
    template <typename Transcript>
    static GroupCommitment receive_group_commitment(const std::shared_ptr<Transcript>& transcript,
                                                    const std::string& label,
                                                    size_t /*num_columns*/,
                                                    const Config& /*config*/)
    {
        return transcript->template receive_from_prover<Commitment>(label);
    }
    template <typename Transcript>
    static void absorb_group_commitment(const std::shared_ptr<Transcript>& transcript,
                                        const std::string& label,
                                        const GroupCommitment& commitment)
    {
        transcript->add_to_hash_buffer(label, commitment);
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
        FflonkProver::prove(ck, claims, u, transcript);
    }
    template <typename Transcript>
    static bool verify_opening(const Config& config,
                               const VerifierClaims& claims,
                               std::span<const fr> u,
                               const std::shared_ptr<Transcript>& transcript)
    {
        return FflonkVerifier::verify(config, claims, u, transcript);
    }
};

/** @brief UltraHonk with the fflonk opening as the polynomial commitment scheme. */
using FflonkHonk = honk_transparent::TransparentHonk<FflonkPcs>;

} // namespace bb::fflonk
