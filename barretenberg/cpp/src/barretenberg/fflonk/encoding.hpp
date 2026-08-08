#pragma once

#include "barretenberg/fflonk/polynomial_utils.hpp"

#include <cstdint>
#include <span>
#include <vector>

namespace bb::fflonk_plonk {

/**
 * @brief The wire format shared by every fflonk proof in this repository: 32-byte big-endian words,
 * a group element as `(x, y)` with the point at infinity as `(0, 0)`.
 *
 * @details Parsing rejects rather than repairs. A non-canonical coordinate or scalar - one at or
 * above the modulus - would be silently reduced by the field constructor, which would give a single
 * proof many encodings and break the binding between a proof's bytes and its transcript.
 */
inline void append_word(std::vector<uint8_t>& buffer, const uint256_t& word)
{
    for (size_t limb = 4; limb-- > 0;) {
        for (size_t byte = 8; byte-- > 0;) {
            buffer.push_back(static_cast<uint8_t>((word.data[limb] >> (8 * byte)) & 0xff));
        }
    }
}

inline uint256_t read_word(std::span<const uint8_t> buffer, size_t offset)
{
    uint256_t word{ 0 };
    for (size_t limb = 4; limb-- > 0;) {
        uint64_t value = 0;
        for (size_t byte = 0; byte < 8; ++byte) {
            value = (value << 8) | buffer[offset++];
        }
        word.data[limb] = value;
    }
    return word;
}

inline void append_point(std::vector<uint8_t>& buffer, const Commitment& point)
{
    if (point.is_point_at_infinity()) {
        append_word(buffer, uint256_t(0));
        append_word(buffer, uint256_t(0));
    } else {
        append_word(buffer, static_cast<uint256_t>(point.x));
        append_word(buffer, static_cast<uint256_t>(point.y));
    }
}

/**
 * @brief Read a group element, rejecting anything a verifier must not accept.
 * @details BN254's G1 has prime order, so being on the curve implies being in the subgroup and no
 * separate cofactor check is needed. Infinity is rejected: an honest prover never produces one - it
 * would mean a committed polynomial was identically zero - and rejecting keeps the native verifier
 * equivalent to a Solidity one, whose `y^2 = x^3 + 3` test has no natural encoding for it.
 */
inline bool read_point(std::span<const uint8_t> buffer, size_t& offset, Commitment& point)
{
    const uint256_t x = read_word(buffer, offset);
    const uint256_t y = read_word(buffer, offset + 32);
    offset += 64;

    if (x >= Curve::BaseField::modulus || y >= Curve::BaseField::modulus) {
        return false;
    }

    point = Commitment(Curve::BaseField(x), Curve::BaseField(y));
    return point.on_curve() && !point.is_point_at_infinity();
}

/**
 * @brief On the curve, and not the point at infinity - the same test `read_point` applies, for a
 * commitment that arrived as a struct rather than as bytes.
 */
inline bool is_valid_point(const Commitment& point)
{
    return !point.is_point_at_infinity() && point.on_curve();
}

/** @brief Read a scalar, rejecting a non-canonical residue. */
inline bool read_scalar(std::span<const uint8_t> buffer, size_t& offset, FF& scalar)
{
    const uint256_t word = read_word(buffer, offset);
    offset += 32;
    if (word >= FF::modulus) {
        return false;
    }
    scalar = FF(word);
    return true;
}

} // namespace bb::fflonk_plonk
