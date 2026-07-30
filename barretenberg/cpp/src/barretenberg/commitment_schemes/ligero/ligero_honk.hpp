#pragma once

#include "barretenberg/commitment_schemes/ligero/ligero.hpp"
#include "barretenberg/commitment_schemes/transparent_honk.hpp"

namespace bb::ligero {

/** @brief PCS backend adapter binding Ligero into the transparent-UltraHonk shell. */
template <typename Hasher_> struct LigeroPcs {
    using Hasher = Hasher_;
    using Config = LigeroConfig;
    using CommitmentKey = LigeroCommitmentKey<Hasher>;
    using GroupData = LigeroGroupData<Hasher>;
    using Digest = typename Hasher::Digest;
    using ProverClaims = typename LigeroProver<Hasher>::Claims;
    using VerifierClaims = typename LigeroVerifier<Hasher>::Claims;

    static Config make_config(size_t log_dyadic_size, size_t security_bits, size_t log_inv_rate)
    {
        return LigeroConfig::create(log_dyadic_size, security_bits, log_inv_rate, /*est_total_polynomials=*/36);
    }
    static size_t payload_variables(const Config& config) { return config.num_variables; }

    using GroupCommitment = typename Hasher::Digest;

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
        LigeroProver<Hasher>::prove(ck, claims, u, transcript);
    }
    template <typename Transcript>
    static bool verify_opening(const Config& config,
                               const VerifierClaims& claims,
                               std::span<const fr> u,
                               const std::shared_ptr<Transcript>& transcript)
    {
        RSDomains domains;
        return LigeroVerifier<Hasher>::verify(config, claims, u, transcript, domains);
    }
};

/** @brief UltraHonk with Ligero as the polynomial commitment scheme (README.md). */
template <typename Hasher = whir::Blake3sMerkleHasher>
using LigeroHonk = honk_transparent::TransparentHonk<LigeroPcs<Hasher>>;

} // namespace bb::ligero
