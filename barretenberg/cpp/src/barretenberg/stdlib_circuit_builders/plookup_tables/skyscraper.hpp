#pragma once

#include "barretenberg/crypto/skyscraper/skyscraper.hpp"

#include "types.hpp"

namespace bb::plookup::skyscraper_tables {

/**
 * Skyscraper's "bar" rotates the canonical little-endian byte string of a field element by 16 bytes,
 * applies the Chi-like S-box of ePrint 2025/058 Table 3 to each byte, and reads the result back.
 * All three steps are one MultiTable read: the key column accumulates the byte decomposition of the
 * input, the value column accumulates the substituted bytes, and the rotation is nothing more than
 * where each value-column coefficient places its byte.
 */
static constexpr uint64_t SBOX_DOMAIN = 256;
static constexpr size_t NUM_SLICES = 32;
static constexpr size_t ROTATION_SLICES = 16;

/**
 * @brief `sbox(key) << shift`.
 * @details Slice 0's contribution lands at byte 16 of the output, and a MultiTable's value column is
 * only defined up to division by its first coefficient - so the only way to get an unrotated,
 * unscaled bar out of the accumulator is for slice 0's table to emit its byte already in place.
 */
template <size_t shift> inline std::array<bb::fr, 2> get_sbox_values_from_key(const std::array<uint64_t, 2> key)
{
    const uint256_t substituted(crypto::skyscraper::sbox(static_cast<uint8_t>(key[0])));
    return { bb::fr(substituted << shift), bb::fr(0) };
}

/** @brief The S-box as a 1-to-2 map: byte in column 1, substituted byte in column 2. */
template <size_t shift> inline BasicTable generate_sbox_table(BasicTableId id, const size_t table_index)
{
    BasicTable table;
    table.id = id;
    table.table_index = table_index;
    table.use_twin_keys = false;

    for (uint64_t i = 0; i < SBOX_DOMAIN; ++i) {
        table.column_1.emplace_back(i);
        table.column_2.emplace_back(uint256_t(crypto::skyscraper::sbox(static_cast<uint8_t>(i))) << shift);
        table.column_3.emplace_back(0);
    }

    table.get_values_from_key = &get_sbox_values_from_key<shift>;
    table.column_1_step_size = SBOX_DOMAIN;
    table.column_2_step_size = SBOX_DOMAIN;
    table.column_3_step_size = 0;
    return table;
}

/**
 * @brief The whole bar as one 32-slice read.
 *
 * @details Byte i of the input becomes byte `(i + 16) mod 32` of the output, so the value column's
 * coefficient for slice i is `2^(8 * ((i + 16) mod 32))`, normalised by slice 0's `2^128` - which
 * slice 0's own table already carries. The key column's coefficients are the plain `2^(8i)`, so its
 * top accumulator is the input itself and the read binds the decomposition for free.
 *
 * @note The read constrains the decomposition only as a field element, and a field element has
 * several 32-byte representations below 2^256. The caller must additionally prove the decomposition
 * is the canonical one - see `stdlib::skyscraper::Skyscraper::bar`.
 */
inline MultiTable get_bar_table(const MultiTableId id = SKYSCRAPER_BAR)
{
    std::vector<bb::fr> key_coefficients;
    std::vector<bb::fr> value_coefficients;
    std::vector<bb::fr> unused_coefficients;
    std::vector<BasicTableId> basic_table_ids;
    std::vector<MultiTable::table_out (*)(MultiTable::table_in)> table_values;
    std::vector<uint64_t> slice_sizes;

    for (size_t i = 0; i < NUM_SLICES; ++i) {
        key_coefficients.emplace_back(uint256_t(1) << (8 * i));
        // Slice 0's own table already emits its byte at position 16, so its coefficient is 1 - the
        // accumulator divides through by it, and scaling here as well would double-count.
        value_coefficients.emplace_back(i == 0 ? uint256_t(1)
                                               : uint256_t(1) << (8 * ((i + ROTATION_SLICES) % NUM_SLICES)));
        unused_coefficients.emplace_back(1);
        slice_sizes.emplace_back(SBOX_DOMAIN);
        basic_table_ids.emplace_back(i == 0 ? SKYSCRAPER_SBOX_SHIFT128 : SKYSCRAPER_SBOX);
        table_values.emplace_back(i == 0 ? &get_sbox_values_from_key<8 * ROTATION_SLICES>
                                         : &get_sbox_values_from_key<0>);
    }

    MultiTable table(key_coefficients, value_coefficients, unused_coefficients);
    table.id = id;
    table.slice_sizes = std::move(slice_sizes);
    table.basic_table_ids = std::move(basic_table_ids);
    table.get_table_values = std::move(table_values);
    return table;
}

} // namespace bb::plookup::skyscraper_tables
