#pragma once

#include "barretenberg/stdlib/primitives/field/field.hpp"

#include <span>
#include <utility>

namespace bb::stdlib::skyscraper {

/**
 * @brief In-circuit Skyscraper-v1 over BN254 Fr, bit-identical to `crypto::skyscraper`.
 *
 * @details The permutation is nine double-rounds over a two-element Feistel state: fourteen
 * squaring half-rounds and four "bar" half-rounds. The split in cost is extreme and it is the whole
 * story of the hash in a circuit.
 *
 * A squaring half-round is `r + l^2 * sigma_inv + c`, which is exactly what one Ultra arithmetic
 * gate expresses (`q_m*l*l + q_1*r + q_c`), so all fourteen cost one gate each.
 *
 * A bar half-round rotates the *canonical* little-endian byte string of the state by 16 bytes,
 * applies a byte S-box, and reads the result back as a field element. Two things make that
 * affordable: the S-box is a 256-row lookup, so each of the 32 bytes costs one gate, and the
 * canonical decomposition is `byte_array(field_t, 32)`, whose marginal cost is small once a circuit
 * has paid for the shared range lists. Canonicity is not optional - a field element has several
 * 32-byte representations below 2^256, and each would give a different bar output.
 */
template <typename Builder> class Skyscraper {
  public:
    using field_ct = field_t<Builder>;

    /** @brief Two-to-one compression: Davies-Meyer on the left branch of the permutation. */
    static field_ct compress(const field_ct& left, const field_ct& right);

    /** @brief The nine-double-round permutation. */
    static std::pair<field_ct, field_ct> permute(const field_ct& left, const field_ct& right);

    /** @brief One bar application: canonical bytes -> rotate 16 -> S-box -> field element. */
    static field_ct bar(const field_ct& x);

    /**
     * @brief Left-fold of `compress` over a sequence, matching `crypto::skyscraper::fold_compress`
     * and hence ProveKit's leaf-hash convention. A single value hashes to itself.
     */
    static field_ct fold_compress(std::span<const field_ct> values);

  private:
    static constexpr size_t NUM_BYTES = 32;

    /** @brief Pin a bar read's byte decomposition to the canonical representative below `r`. */
    static void enforce_canonical(const field_ct& value, const field_ct& high);
};

} // namespace bb::stdlib::skyscraper
