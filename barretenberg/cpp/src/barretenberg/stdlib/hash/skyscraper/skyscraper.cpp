#include "skyscraper.hpp"

#include "barretenberg/crypto/skyscraper/skyscraper.hpp"
#include "barretenberg/stdlib/primitives/circuit_builders/circuit_builders.hpp"
#include "barretenberg/stdlib/primitives/plookup/plookup.hpp"
#include "barretenberg/stdlib/primitives/witness/witness.hpp"

namespace bb::stdlib::skyscraper {

/**
 * @brief Rotate the canonical little-endian byte string by 16, S-box each byte, read it back.
 *
 * @details One `SKYSCRAPER_BAR` read does all three: its key column accumulates the byte
 * decomposition of `x` (so the decomposition is bound without a separate `byte_array`), its value
 * column accumulates the substituted bytes already rotated into place, and the top accumulator of
 * each is the value we want. That is 32 gates for what otherwise costs a decomposition, 32 separate
 * reads and a 32-term recomposition.
 *
 * What the read does *not* give is canonicity. Its key accumulator binds
 * `sum b_i 2^(8i) == x` in the field, and a field element has four such 32-byte representations
 * below 2^256 - each of which the bar would map somewhere different, so the prover could choose.
 * `enforce_canonical` pins it to the representation below the modulus, reusing the two 128-bit
 * halves the read has already accumulated rather than decomposing anything a second time.
 */
template <typename Builder> field_t<Builder> Skyscraper<Builder>::bar(const field_ct& x)
{
    const auto lookup = plookup_read<Builder>::get_lookup_accumulators(plookup::MultiTableId::SKYSCRAPER_BAR, x);
    enforce_canonical(x, lookup[plookup::ColumnIdx::C1][NUM_BYTES / 2]);
    return lookup[plookup::ColumnIdx::C2][0];
}

/**
 * @brief Prove that the byte decomposition behind a `SKYSCRAPER_BAR` read is the canonical one.
 *
 * @details `high` is the read's own accumulator at slice 16, so it is the top 128 bits of the
 * decomposition and `value - 2^128 * high` is the bottom 128 bits; both are already pinned below
 * 2^128 by the slices themselves. What is left is to show the pair is at most `r - 1`, which is the
 * comparison `byte_array(field_t, 32)` makes: shift the low-half difference up by 2^128 so its
 * 129th bit is the borrow, then require the high-half difference minus the borrow's complement to
 * be a non-negative 128-bit number.
 */
template <typename Builder> void Skyscraper<Builder>::enforce_canonical(const field_ct& value, const field_ct& high)
{
    constexpr uint256_t modulus_minus_one = bb::fr::modulus - 1;
    constexpr uint256_t s_lo = modulus_minus_one.slice(0, 128);
    constexpr uint256_t s_hi = modulus_minus_one.slice(128, 256);
    constexpr uint256_t shift = uint256_t(1) << 128;

    if (value.is_constant()) {
        // A constant is already the canonical representative; there is no prover to constrain.
        return;
    }
    Builder* ctx = value.get_context();

    const field_ct low = value - (high * bb::fr(shift));

    const field_ct diff_lo = -low + bb::fr(s_lo) + bb::fr(shift);
    const uint256_t diff_lo_value(diff_lo.get_value());
    field_ct borrow = witness_t<Builder>(ctx, bb::fr(diff_lo_value >> 128));
    field_ct diff_lo_lo = witness_t<Builder>(ctx, bb::fr(diff_lo_value & (shift - 1)));
    // Both are functions of `value`, so they inherit its Fiat-Shamir provenance. Left untagged they
    // would look like free witnesses the moment they are compared against `diff_lo`, which carries
    // the transcript's origin.
    borrow.set_origin_tag(value.get_origin_tag());
    diff_lo_lo.set_origin_tag(value.get_origin_tag());
    borrow.create_range_constraint(1, "skyscraper: canonicity borrow is not a bit");
    diff_lo_lo.create_range_constraint(128, "skyscraper: canonicity low limb exceeds 128 bits");
    diff_lo.assert_equal(diff_lo_lo + borrow * bb::fr(shift), "skyscraper: canonicity low split");

    // borrow == 1 exactly when the low half does not already force the comparison, so the high half
    // must then be strictly smaller.
    const field_ct diff_hi = -high + bb::fr(s_hi) - (field_ct(1) - borrow);
    diff_hi.create_range_constraint(128, "skyscraper: value is not a canonical field element");
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
