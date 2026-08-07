#pragma once

#include "barretenberg/commitment_schemes/transparent_honk.hpp"
#include "barretenberg/commitment_schemes/whir/whir.hpp"

namespace bb::whir {

/**
 * @brief PCS backend adapter binding WHIR into the transparent-UltraHonk shell.
 * @tparam MaxStackBits 0 keeps one codeword per column (fast prover, larger proof); non-zero stacks
 * each commitment group into one narrow-leaf oracle sized for up to 2^MaxStackBits columns, which
 * shrinks the proof several-fold at a proportional prover cost. See `WhirGroupData`.
 * @tparam Soundness the proximity regime the query schedule is derived from. The default is the
 * unconditional Johnson bound, which rests on no conjecture.
 * @tparam FoldingFactorBits k, the number of variables folded per iteration.
 * @tparam InitialFoldingFactorBits k₀, the first iteration's fold; 0 means "same as k". A round-0
 * query reveals 2^{k₀} values of *every* committed column, so on the ~30-column Honk commitment this
 * is the single largest term in the proof; the later rounds run on one column and are free to fold
 * faster.
 */
template <typename Hasher_,
          size_t MaxStackBits = 0,
          WhirSoundness Soundness = WhirSoundness::PROVABLE_LIST,
          size_t FoldingFactorBits = 4,
          size_t InitialFoldingFactorBits = 0,
          size_t PowBits = 0>
struct WhirPcs {
    using Hasher = Hasher_;
    using Config = WhirConfig;
    using CommitmentKey = WhirCommitmentKey<Hasher>;
    using GroupData = WhirGroupData<Hasher>;
    using Digest = typename Hasher::Digest;
    using ProverClaims = typename WhirProver<Hasher>::Claims;
    using VerifierClaims = typename WhirVerifier<Hasher>::Claims;

    using GroupCommitment = typename Hasher::Digest;

    static Config make_config(size_t log_dyadic_size, size_t security_bits, size_t log_inv_rate)
    {
        return WhirConfig::create(log_dyadic_size,
                                  security_bits,
                                  log_inv_rate,
                                  FoldingFactorBits,
                                  /*final_poly_bits=*/4,
                                  Soundness,
                                  /*zk=*/false,
                                  MaxStackBits,
                                  InitialFoldingFactorBits,
                                  PowBits);
    }
    static size_t payload_variables(const Config& config) { return config.num_payload_variables; }

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

/**
 * @brief WhirHonk with per-group column stacking: roughly a third of the proof size, at several
 * times the prover cost. Choose it only when proof size dominates.
 * @details The headroom covers the widest commitment group, the (at most 28) precomputed columns.
 */
template <typename Hasher = Poseidon2MerkleHasher>
using WhirStackedHonk = honk_transparent::TransparentHonk<WhirPcs<Hasher, 5>>;

/**
 * @brief WhirHonk on ProveKit's own WHIR parameters: the Johnson bound with k = 3.
 * @details ProveKit runs the reference WHIR at λ = 128, rate 2^-2, k = 3 and 10 bits of per-round
 * grinding. Grinding is not implemented here, so at the same λ this schedule draws every bit of
 * soundness from queries and is strictly the more conservative of the two.
 */
template <typename Hasher = SkyscraperMerkleHasher>
using WhirProveKitHonk =
    honk_transparent::TransparentHonk<WhirPcs<Hasher, 0, WhirSoundness::PROVABLE_LIST, /*FoldingFactorBits=*/3>>;

/**
 * @brief WhirProveKitHonk with the first fold narrowed to 2, the rest left at ProveKit's 8.
 * @details Round-0 queries are the only ones that touch the wide commitment, and each reveals
 * 2^{k₀} values of all ~30 opened columns; the later rounds run against a single-column oracle where
 * a larger fold is what shortens the schedule. Trading k₀ = 3 for k₀ = 1 costs one extra iteration
 * and two extra Merkle levels on the round-0 trees, and removes three quarters of the proof's
 * largest term.
 */
template <typename Hasher = SkyscraperMerkleHasher>
using WhirNarrowHonk = honk_transparent::TransparentHonk<
    WhirPcs<Hasher, 0, WhirSoundness::PROVABLE_LIST, /*FoldingFactorBits=*/3, /*InitialFoldingFactorBits=*/1>>;

/**
 * @brief WhirNarrowHonk with 20 bits of per-round grinding: the smallest proof at λ = 128.
 * @details Grinding is the cheapest soundness bit there is here — one nonce in the proof and a
 * Blake3 search the honest prover runs once per round, against a query that costs the wide
 * commitment's whole leaf width. ProveKit takes 10 bits this way; 20 is still milliseconds.
 */
template <typename Hasher = SkyscraperMerkleHasher>
using WhirCompactHonk = honk_transparent::TransparentHonk<WhirPcs<Hasher,
                                                                  0,
                                                                  WhirSoundness::PROVABLE_LIST,
                                                                  /*FoldingFactorBits=*/3,
                                                                  /*InitialFoldingFactorBits=*/1,
                                                                  /*PowBits=*/20>>;

/**
 * @brief WhirHonk under the up-to-capacity conjecture, which eprint 2025/2046 disproves.
 * @details Retained only to reproduce the published FRI/STIR-style figures that assume it; it is not
 * a sound configuration to deploy. See README.md §6.
 */
template <typename Hasher = Poseidon2MerkleHasher>
using WhirConjecturedHonk = honk_transparent::TransparentHonk<WhirPcs<Hasher, 0, WhirSoundness::CONJECTURED_LIST>>;

} // namespace bb::whir
