#pragma once

#include "barretenberg/commitment_schemes/transparent_honk.hpp"
#include "barretenberg/commitment_schemes/whir/whir.hpp"

namespace bb::whir {

/** @brief PCS backend adapter binding WHIR into the transparent-UltraHonk shell. */
template <typename Hasher_> struct WhirPcs {
    using Hasher = Hasher_;
    using Config = WhirConfig;
    using CommitmentKey = WhirCommitmentKey<Hasher>;
    using GroupData = WhirGroupData<Hasher>;
    using Digest = typename Hasher::Digest;
    using ProverClaims = typename WhirProver<Hasher>::Claims;
    using VerifierClaims = typename WhirVerifier<Hasher>::Claims;

    static Config make_config(size_t log_dyadic_size, size_t security_bits, size_t log_inv_rate)
    {
        return WhirConfig::create(log_dyadic_size, security_bits, log_inv_rate);
    }
    static size_t payload_variables(const Config& config) { return config.num_payload_variables; }

    template <typename Transcript>
    static void prove_opening(const CommitmentKey& ck,
                              const ProverClaims& claims,
                              std::span<const fr> u,
                              const std::shared_ptr<Transcript>& transcript)
    {
        WhirProver<Hasher>::prove(ck, claims, u, transcript);
    }
    template <typename Transcript>
    static bool verify_opening(const Config& config,
                               const VerifierClaims& claims,
                               std::span<const fr> u,
                               const std::shared_ptr<Transcript>& transcript)
    {
        return WhirVerifier<Hasher>::verify(config, claims, u, transcript);
    }
};

/** @brief UltraHonk with WHIR as the polynomial commitment scheme (README.md §10). */
template <typename Hasher = Poseidon2MerkleHasher> using WhirHonk = honk_transparent::TransparentHonk<WhirPcs<Hasher>>;

} // namespace bb::whir
