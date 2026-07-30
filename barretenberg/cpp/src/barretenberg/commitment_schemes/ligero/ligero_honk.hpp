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
