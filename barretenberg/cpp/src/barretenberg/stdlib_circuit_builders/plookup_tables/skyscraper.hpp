#pragma once

#include "barretenberg/crypto/skyscraper/skyscraper.hpp"

#include "types.hpp"

namespace bb::plookup::skyscraper_tables {

/**
 * Skyscraper's "bar" applies the Chi-like S-box of ePrint 2025/058 Table 3 to each of the 32
 * canonical little-endian bytes of a field element. A 256-row table turns that into one lookup gate
 * per byte, which is what makes the hash affordable in a circuit at all.
 */
static constexpr uint64_t SBOX_DOMAIN = 256;

inline std::array<bb::fr, 2> get_sbox_values_from_key(const std::array<uint64_t, 2> key)
{
    return { bb::fr(uint256_t(crypto::skyscraper::sbox(static_cast<uint8_t>(key[0])))), bb::fr(0) };
}

/** @brief The S-box as a 1-to-2 map: byte in column 1, S-boxed byte in column 2. */
inline BasicTable generate_sbox_table(BasicTableId id, const size_t table_index)
{
    BasicTable table;
    table.id = id;
    table.table_index = table_index;
    table.use_twin_keys = false;

    for (uint64_t i = 0; i < SBOX_DOMAIN; ++i) {
        table.column_1.emplace_back(i);
        table.column_2.emplace_back(uint256_t(crypto::skyscraper::sbox(static_cast<uint8_t>(i))));
        table.column_3.emplace_back(0);
    }

    table.get_values_from_key = &get_sbox_values_from_key;
    table.column_1_step_size = SBOX_DOMAIN;
    table.column_2_step_size = SBOX_DOMAIN;
    table.column_3_step_size = 0;
    return table;
}

/** @brief A single-byte read: one lookup gate, with the byte range-constrained by the table itself. */
inline MultiTable get_sbox_multitable(const MultiTableId id = SKYSCRAPER_SBOX_MULTI)
{
    MultiTable table(SBOX_DOMAIN, 0, 0, 1);
    table.id = id;
    table.slice_sizes = { SBOX_DOMAIN };
    table.basic_table_ids = { SKYSCRAPER_SBOX };
    table.get_table_values = { &get_sbox_values_from_key };
    return table;
}

} // namespace bb::plookup::skyscraper_tables
