#include "skyscraper.hpp"

#include "barretenberg/crypto/skyscraper/skyscraper.hpp"
#include "barretenberg/stdlib/primitives/byte_array/byte_array.hpp"
#include "barretenberg/stdlib/primitives/circuit_builders/circuit_builders.hpp"
#include "barretenberg/stdlib/primitives/plookup/plookup.hpp"

namespace bb::stdlib::skyscraper {

namespace {
constexpr size_t NUM_BYTES = 32;
constexpr size_t ROTATION_BYTES = 16;
} // namespace

/**
 * @brief Rotate the canonical little-endian byte string by 16, S-box each byte, read it back.
 *
 * @details `byte_array(field_t, 32)` is what makes this sound: a field element has several 32-byte
 * representations below 2^256 and each would give a different result, so the decomposition has to be
 * the canonical one. Its constructor imposes that. The 32 S-box reads then cost one gate each and
 * range-constrain their own input byte, and the rotation is free - it is just where each substituted
 * byte is placed in the recomposition.
 */
template <typename Builder> field_t<Builder> Skyscraper<Builder>::bar(const field_ct& x)
{
    using byte_array_ct = byte_array<Builder>;

    // `byte_array` is big-endian; Skyscraper's byte string is the canonical little-endian one.
    byte_array_ct bytes(x, NUM_BYTES);

    std::vector<field_ct> terms;
    terms.reserve(NUM_BYTES);
    for (size_t i = 0; i < NUM_BYTES; ++i) {
        const field_ct byte = bytes[NUM_BYTES - 1 - i];
        const field_ct substituted =
            plookup_read<Builder>::read_from_1_to_2_table(plookup::MultiTableId::SKYSCRAPER_SBOX_MULTI, byte);
        // rotated[j] = le[(j + 16) mod 32], so byte i of the input lands at position (i - 16) mod 32.
        const size_t position = (i + NUM_BYTES - ROTATION_BYTES) % NUM_BYTES;
        terms.push_back(substituted * bb::fr(uint256_t(1) << (8 * position)));
    }
    return field_ct::accumulate(terms);
}

template <typename Builder>
std::pair<field_t<Builder>, field_t<Builder>> Skyscraper<Builder>::permute(const field_ct& left, const field_ct& right)
{
    const bb::fr sigma_inv = crypto::skyscraper::sigma_inv();

    field_ct l = left;
    field_ct r = right;

    // `r + l^2 * sigma_inv + c`. Scaling by sigma_inv and adding the constant are selector changes,
    // so `madd` puts the whole half-round in one gate.
    const auto square_half_round = [&](size_t round) {
        r = l.madd(l * sigma_inv, r + crypto::skyscraper::round_constant(round));
        std::swap(l, r);
    };
    const auto bar_half_round = [&](size_t round) {
        r = r + bar(l) + crypto::skyscraper::round_constant(round);
        std::swap(l, r);
    };
    const auto square_round = [&](size_t round) {
        square_half_round(round);
        square_half_round(round + 1);
    };
    const auto bar_round = [&](size_t round) {
        bar_half_round(round);
        bar_half_round(round + 1);
    };

    square_round(0);
    square_round(2);
    square_round(4);
    bar_round(6);
    square_round(8);
    bar_round(10);
    square_round(12);
    square_round(14);
    square_round(16);
    return { l, r };
}

template <typename Builder> field_t<Builder> Skyscraper<Builder>::compress(const field_ct& left, const field_ct& right)
{
    return permute(left, right).first + left;
}

template <typename Builder> field_t<Builder> Skyscraper<Builder>::fold_compress(std::span<const field_ct> values)
{
    BB_ASSERT_GT(values.size(), static_cast<size_t>(0), "skyscraper fold of an empty sequence");
    field_ct acc = values[0];
    for (size_t i = 1; i < values.size(); ++i) {
        acc = compress(acc, values[i]);
    }
    return acc;
}

template class Skyscraper<UltraCircuitBuilder>;
template class Skyscraper<MegaCircuitBuilder>;

} // namespace bb::stdlib::skyscraper
