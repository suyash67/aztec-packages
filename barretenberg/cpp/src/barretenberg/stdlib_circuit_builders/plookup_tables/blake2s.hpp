// === AUDIT STATUS ===
// internal:    { status: Complete, auditors: [Nishat], commit: 66052c96cc754339ac3f2761f341f150130555b3}
// external_1:  { status: not started, auditors: [], commit: }
// external_2:  { status: not started, auditors: [], commit: }
// =====================

#pragma once

#include "barretenberg/numeric/bitop/rotate.hpp"

#include "sparse.hpp"
#include "types.hpp"

namespace bb::plookup::blake2s_tables {

/**
 * Blake2s/Blake3s work on 32-bit words, which we split into four 8-bit slices. Each slice pair is looked up in a
 * 256 x 256 (2^16-row) XOR table, so a 32-bit XOR costs 4 lookup gates.
 */
static constexpr uint64_t BITS_PER_SLICE = 8UL;
static constexpr uint64_t SLICE_SIZE = (1UL << BITS_PER_SLICE);
static constexpr size_t NUM_SLICES = 32 / BITS_PER_SLICE;

/**
 * This function performs the operation ROTR^{k}(a ^ b) on `bits_per_slice`-wide inputs, where the rotation is taken
 * over the full 32-bit word. Note that for a slice value x < 2^8 and k in {8, 12, 16}, `rotate32(x, k)` is just
 * `x << (32 - k)`, which is how a slice gets placed at its final position in the rotated word.
 */
template <uint64_t bits_per_slice, uint64_t num_rotated_output_bits>
inline std::array<bb::fr, 2> get_xor_rotate_values_from_key(const std::array<uint64_t, 2> key)
{
    return { uint256_t(numeric::rotate32(uint32_t(key[0]) ^ uint32_t(key[1]), uint32_t(num_rotated_output_bits))),
             0ULL };
}

/**
 * Generates a basic 32-bit (XOR + ROTR) lookup table.
 */
template <uint64_t bits_per_slice, uint64_t num_rotated_output_bits>
inline BasicTable generate_xor_rotate_table(BasicTableId id, const size_t table_index)
{
    const uint64_t base = 1UL << bits_per_slice;
    BasicTable table;
    table.id = id;
    table.table_index = table_index;
    table.use_twin_keys = true;

    for (uint64_t i = 0; i < base; ++i) {
        for (uint64_t j = 0; j < base; ++j) {
            table.column_1.emplace_back(i);
            table.column_2.emplace_back(j);
            table.column_3.emplace_back(
                uint256_t(numeric::rotate32(uint32_t(i) ^ uint32_t(j), uint32_t(num_rotated_output_bits))));
        }
    }

    table.get_values_from_key = &get_xor_rotate_values_from_key<bits_per_slice, num_rotated_output_bits>;

    table.column_1_step_size = base;
    table.column_2_step_size = base;
    table.column_3_step_size = base;

    return table;
}

/**
 * @brief Assemble a 4-slice MultiTable computing ROTR^{k}(a ^ b) on 32-bit inputs.
 *
 * @details Write the 32-bit XOR result as u = s0 + 2^8.s1 + 2^16.s2 + 2^24.s3. The MultiTable's third column
 * accumulates to `sum_i c_i * raw_i` where `raw_i` is what the i-th basic table returns for slice i, and c_0 is
 * pinned to 1 by the accumulator construction (see plookup_tables.cpp). Since the accumulator is only ever defined
 * up to division by c_0, the *only* way to get an unscaled rotated word out of the lookup is for slice 0's basic
 * table to already emit its contribution at its final bit position - which is exactly what a ROTR^{k} basic table
 * does for a single slice. Every rotation below therefore uses a dedicated slice-0 table and plain XOR elsewhere,
 * with one exception: ROTR^{12} straddles slice 1 at a nibble boundary, so slice 1 needs its own ROTR^{4} table.
 *
 * Getting an unscaled result matters beyond aesthetics: a scaled `field_t` has to be normalized (one gate) before it
 * can be used as a lookup key, and in Blake every rotation output feeds a later XOR.
 *
 * -------------------------------------------------------------------------
 * | slice | ROTR_16      | ROTR_12      | ROTR_8       | ROTR_7           |
 * |-------|--------------|--------------|--------------|------------------|
 * | s0    | ROTR16, c=1  | ROTR12, c=1  | ROTR8,  c=1  | ROTR7, c=1       |
 * | s1    | XOR,  c=2^24 | ROTR4,  c=1  | XOR,    c=1  | XOR,   c=2       |
 * | s2    | XOR,  c=1    | XOR,  c=2^4  | XOR,  c=2^8  | XOR,   c=2^9     |
 * | s3    | XOR,  c=2^8  | XOR,  c=2^12 | XOR,  c=2^16 | XOR,   c=2^17    |
 * -------------------------------------------------------------------------
 */
using slice_fn = MultiTable::table_out (*)(MultiTable::table_in);

inline MultiTable get_blake2s_xor_rotate_table(const MultiTableId id,
                                               const BasicTableId slice_0_table,
                                               const slice_fn slice_0_values,
                                               const BasicTableId slice_1_table,
                                               const slice_fn slice_1_values,
                                               const std::array<uint64_t, NUM_SLICES>& column_3_shifts)
{
    std::vector<bb::fr> key_coefficients{ bb::fr(1), bb::fr(1 << 8), bb::fr(1 << 16), bb::fr(uint256_t(1) << 24) };
    std::vector<bb::fr> column_3_coefficients;
    for (const auto& shift : column_3_shifts) {
        column_3_coefficients.emplace_back(bb::fr(uint256_t(1) << shift));
    }

    MultiTable table(key_coefficients, key_coefficients, column_3_coefficients);
    table.id = id;
    table.slice_sizes = { SLICE_SIZE, SLICE_SIZE, SLICE_SIZE, SLICE_SIZE };
    table.basic_table_ids = { slice_0_table, slice_1_table, BLAKE_XOR8, BLAKE_XOR8 };
    table.get_table_values = { slice_0_values,
                               slice_1_values,
                               &get_xor_rotate_values_from_key<BITS_PER_SLICE, 0>,
                               &get_xor_rotate_values_from_key<BITS_PER_SLICE, 0> };
    return table;
}

/**
 * Generates a 4-slice MultiTable for the 32-bit operation (a ^ b).
 */
inline MultiTable get_blake2s_xor_table(const MultiTableId id = BLAKE_XOR)
{
    return get_blake2s_xor_rotate_table(id,
                                        BLAKE_XOR8,
                                        &get_xor_rotate_values_from_key<BITS_PER_SLICE, 0>,
                                        BLAKE_XOR8,
                                        &get_xor_rotate_values_from_key<BITS_PER_SLICE, 0>,
                                        { 0, 8, 16, 24 });
}

/**
 * Generates a 4-slice MultiTable for the 32-bit operation ROTR^{16}(a ^ b).
 */
inline MultiTable get_blake2s_xor_rotate_16_table(const MultiTableId id = BLAKE_XOR_ROTATE_16)
{
    return get_blake2s_xor_rotate_table(id,
                                        BLAKE_XOR8_ROTATE16,
                                        &get_xor_rotate_values_from_key<BITS_PER_SLICE, 16>,
                                        BLAKE_XOR8,
                                        &get_xor_rotate_values_from_key<BITS_PER_SLICE, 0>,
                                        { 0, 24, 0, 8 });
}

/**
 * Generates a 4-slice MultiTable for the 32-bit operation ROTR^{12}(a ^ b).
 */
inline MultiTable get_blake2s_xor_rotate_12_table(const MultiTableId id = BLAKE_XOR_ROTATE_12)
{
    return get_blake2s_xor_rotate_table(id,
                                        BLAKE_XOR8_ROTATE12,
                                        &get_xor_rotate_values_from_key<BITS_PER_SLICE, 12>,
                                        BLAKE_XOR8_ROTATE4,
                                        &get_xor_rotate_values_from_key<BITS_PER_SLICE, 4>,
                                        { 0, 0, 4, 12 });
}

/**
 * Generates a 4-slice MultiTable for the 32-bit operation ROTR^{8}(a ^ b).
 */
inline MultiTable get_blake2s_xor_rotate_8_table(const MultiTableId id = BLAKE_XOR_ROTATE_8)
{
    return get_blake2s_xor_rotate_table(id,
                                        BLAKE_XOR8_ROTATE8,
                                        &get_xor_rotate_values_from_key<BITS_PER_SLICE, 8>,
                                        BLAKE_XOR8,
                                        &get_xor_rotate_values_from_key<BITS_PER_SLICE, 0>,
                                        { 0, 0, 8, 16 });
}

/**
 * Generates a 4-slice MultiTable for the 32-bit operation ROTR^{7}(a ^ b).
 */
inline MultiTable get_blake2s_xor_rotate_7_table(const MultiTableId id = BLAKE_XOR_ROTATE_7)
{
    return get_blake2s_xor_rotate_table(id,
                                        BLAKE_XOR8_ROTATE7,
                                        &get_xor_rotate_values_from_key<BITS_PER_SLICE, 7>,
                                        BLAKE_XOR8,
                                        &get_xor_rotate_values_from_key<BITS_PER_SLICE, 0>,
                                        { 0, 1, 9, 17 });
}

} // namespace bb::plookup::blake2s_tables
