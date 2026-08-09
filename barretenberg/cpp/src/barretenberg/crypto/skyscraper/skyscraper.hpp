#pragma once

#include "barretenberg/ecc/curves/bn254/fr.hpp"

#include <span>
#include <utility>

namespace bb::crypto::skyscraper {

/**
 * @brief Skyscraper-v1 permutation and compression over BN254 Fr (ePrint 2025/058).
 *
 * @details Matches ProveKit's reference implementation (worldfnd/ProveKit skyscraper/core,
 * `Sky_BN254_1`): a two-element Feistel state with nine double-rounds — squarings scaled by
 * sigma^-1 in rounds {0,1,2,4,6,7,8} and byte-level "bar" rounds (rotate the 32 canonical
 * little-endian bytes by 16, apply the Chi-like S-box to each byte, reduce mod r) in rounds
 * {3,5} — with 18 round constants. Compression is Davies-Meyer on the left branch:
 * compress(l, r) = permute(l, r).first + l.
 */

/** @brief The Chi-like byte S-box, Table 3 of the paper. */
uint8_t sbox(uint8_t v);

/** @brief sigma^-1, the scaling applied to every squaring round. */
const fr& sigma_inv();

/** @brief Round constant `index` of the 18 the permutation uses. */
const fr& round_constant(size_t index);

/** @brief The full 9-double-round permutation, Figure 2.a. */
std::pair<fr, fr> permute(const fr& left, const fr& right);

/** @brief Two-to-one compression: left feed-forward over the permutation. */
fr compress(const fr& left, const fr& right);

/**
 * @brief Left-fold of `compress` over a sequence, ProveKit's leaf-hash convention:
 * h = values[0]; h = compress(h, values[i]) for i >= 1. A single value hashes to itself.
 */
fr fold_compress(std::span<const fr> values);

} // namespace bb::crypto::skyscraper
