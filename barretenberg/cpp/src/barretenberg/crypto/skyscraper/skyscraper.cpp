#include "barretenberg/crypto/skyscraper/skyscraper.hpp"

#include "barretenberg/common/assert.hpp"
#include "barretenberg/numeric/uint256/uint256.hpp"

#include <array>
#include <cstring>

namespace bb::crypto::skyscraper {

namespace {

// sigma^-1 for Sky_BN254_1 (skyscraper.sage: Sky_BN254_1.sigma_inv).
const fr SIGMA_INV = fr(uint256_t("0x15ebf95182c5551cc8260de4aeb85d5d090ef5a9e111ec87dc5ba0056db1194e"));

// Round constants for BN254-Fr, t = 1, canonical little-endian limbs (ProveKit
// skyscraper/core/src/constants.rs, generated from the reference sage implementation).
const std::array<fr, 18> ROUND_CONSTANTS = {
    fr(0),
    fr(uint256_t(0x903c4324270bd744ULL, 0x873125f708a7d269ULL, 0x081dd27906c83855ULL, 0x276b1823ea6d7667ULL)),
    fr(uint256_t(0x7ac8edbb4b378d71ULL, 0xe29d79f3d99e2cb7ULL, 0x751417914c1a5a18ULL, 0x0cf02bd758a484a6ULL)),
    fr(uint256_t(0xfa7adc6769e5bc36ULL, 0x1c3f8e297cca387dULL, 0x0eb7730d63481db0ULL, 0x25b0e03f18ede544ULL)),
    fr(uint256_t(0x57847e652f03cfb7ULL, 0x33440b9668873404ULL, 0x955a32e849af80bcULL, 0x002882fcbe14ae70ULL)),
    fr(uint256_t(0x979231396257d4d7ULL, 0x29989c3e1b37d3c1ULL, 0x12ef02b47f1277baULL, 0x039ad8571e2b7a9cULL)),
    fr(uint256_t(0xb5b48465abbb7887ULL, 0xa72a6bc5e6ba2d2bULL, 0x4cd48043712f7b29ULL, 0x1142d5410fc1fc1aULL)),
    fr(uint256_t(0x7ab2c156059075d3ULL, 0x17cb3594047999b2ULL, 0x44f2c93598f289f7ULL, 0x1d78439f69bc0becULL)),
    fr(uint256_t(0x05d7a965138b8edbULL, 0x36ef35a3d55c48b1ULL, 0x8ddfb8a1ac6f1628ULL, 0x258588a508f4ff82ULL)),
    fr(uint256_t(0x1596fb9afccb49e9ULL, 0x9a7367d69a09a95bULL, 0x9bc43f6984e4c157ULL, 0x13087879d2f514feULL)),
    fr(uint256_t(0x295ccd233b4109faULL, 0xe1d72f89ed868012ULL, 0x2e9e1eea4bc88a8eULL, 0x17dadee898c45232ULL)),
    fr(uint256_t(0x9a8590b4aa1f486fULL, 0xb75834b430e9130eULL, 0xb8e90b1034d5de31ULL, 0x295c6d1546e7f4a6ULL)),
    fr(uint256_t(0x850adcb74c6eb892ULL, 0x07699ef305b92fc3ULL, 0x4ef96a2ba1720f2dULL, 0x1288ca0e1d3ed446ULL)),
    fr(uint256_t(0x01960f9349d1b5eeULL, 0x8ccad30769371c69ULL, 0xe5c81e8991c98662ULL, 0x17563b4d1ae023f3ULL)),
    fr(uint256_t(0x6ba01e9476b32917ULL, 0xa1cb0a3add977bc9ULL, 0x86815a945815f030ULL, 0x2869043be91a1eeaULL)),
    fr(uint256_t(0x81776c885511d976ULL, 0x7475d34f47f414e7ULL, 0x5d090056095d96cfULL, 0x14941f0aff59e79aULL)),
    fr(uint256_t(0xbc40b4fd8fc8c034ULL, 0xbb7142c3cce4fd48ULL, 0x318356758a39005aULL, 0x1ce337a190f4379fULL)),
    fr(0),
};

constexpr uint8_t rotl8(uint8_t v, unsigned n)
{
    return static_cast<uint8_t>(static_cast<uint8_t>(v << n) | static_cast<uint8_t>(v >> (8U - n)));
}

constexpr uint8_t sbox_of(uint8_t v)
{
    const auto chi = static_cast<uint8_t>(rotl8(static_cast<uint8_t>(~v), 1) & rotl8(v, 2) & rotl8(v, 3));
    return rotl8(static_cast<uint8_t>(v ^ chi), 1);
}

// The S-box is a bijection on a byte, so it is a 256-entry table. Evaluating it costs four rotates
// and three bitwise ops; a compression applies it 128 times, which is enough for the table to be
// worth its cache line.
constexpr std::array<uint8_t, 256> SBOX_TABLE = [] {
    std::array<uint8_t, 256> table{};
    for (size_t i = 0; i < 256; ++i) {
        table[i] = sbox_of(static_cast<uint8_t>(i));
    }
    return table;
}();

/** @brief The S-box applied to each of a limb's eight bytes. */
constexpr uint64_t sbox_limb(uint64_t limb)
{
    uint64_t out = 0;
    for (size_t byte = 0; byte < 8; ++byte) {
        out |= static_cast<uint64_t>(SBOX_TABLE[(limb >> (8 * byte)) & 0xFFULL]) << (8 * byte);
    }
    return out;
}

/**
 * @brief Rotate the canonical little-endian byte string by 16, S-box each byte, reduce mod r.
 * @details Rotating by 16 bytes is exactly swapping the two 128-bit halves, so it is a relabelling
 * of the limbs rather than anything the bytes have to be materialised for. Working on the limbs also
 * makes the routine endianness-independent, which the byte-copy version was not.
 */
fr bar(const fr& x)
{
    const uint256_t canonical(x);
    return fr(uint256_t(sbox_limb(canonical.data[2]),
                        sbox_limb(canonical.data[3]),
                        sbox_limb(canonical.data[0]),
                        sbox_limb(canonical.data[1])));
}

void square_round(size_t round, fr& l, fr& r)
{
    r += l.sqr() * SIGMA_INV + ROUND_CONSTANTS[round];
    std::swap(l, r);
    r += l.sqr() * SIGMA_INV + ROUND_CONSTANTS[round + 1];
    std::swap(l, r);
}

void bar_round(size_t round, fr& l, fr& r)
{
    r += bar(l) + ROUND_CONSTANTS[round];
    std::swap(l, r);
    r += bar(l) + ROUND_CONSTANTS[round + 1];
    std::swap(l, r);
}

} // namespace

const fr& sigma_inv()
{
    return SIGMA_INV;
}

const fr& round_constant(size_t index)
{
    BB_ASSERT_LT(index, ROUND_CONSTANTS.size(), "skyscraper round constant out of range");
    return ROUND_CONSTANTS[index];
}

uint8_t sbox(uint8_t v)
{
    return SBOX_TABLE[v];
}

std::pair<fr, fr> permute(const fr& left, const fr& right)
{
    fr l = left;
    fr r = right;
    square_round(0, l, r);
    square_round(2, l, r);
    square_round(4, l, r);
    bar_round(6, l, r);
    square_round(8, l, r);
    bar_round(10, l, r);
    square_round(12, l, r);
    square_round(14, l, r);
    square_round(16, l, r);
    return { l, r };
}

fr compress(const fr& left, const fr& right)
{
    return permute(left, right).first + left;
}

fr fold_compress(std::span<const fr> values)
{
    BB_ASSERT_GT(values.size(), static_cast<size_t>(0), "skyscraper fold of an empty sequence");
    fr acc = values[0];
    for (size_t i = 1; i < values.size(); ++i) {
        acc = compress(acc, values[i]);
    }
    return acc;
}

} // namespace bb::crypto::skyscraper
