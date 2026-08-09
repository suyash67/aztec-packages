// === AUDIT STATUS ===
// internal:    { status: Complete, auditors: [Nishat], commit: 66052c96cc754339ac3f2761f341f150130555b3}
// external_1:  { status: not started, auditors: [], commit: }
// external_2:  { status: not started, auditors: [], commit: }
// =====================

#pragma once
#include "barretenberg/stdlib/hash/hash_utils.hpp"
#include "barretenberg/stdlib/primitives/plookup/plookup.hpp"
#include "barretenberg/stdlib_circuit_builders/plookup_tables/plookup_tables.hpp"

namespace bb::stdlib::blake_util {

using namespace bb::plookup;

// constants
enum blake_constant { BLAKE_STATE_SIZE = 16 };

constexpr uint8_t MSG_SCHEDULE_BLAKE3[7][16] = {
    { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 }, { 2, 6, 3, 10, 7, 0, 4, 13, 1, 11, 12, 5, 9, 14, 15, 8 },
    { 3, 4, 10, 12, 13, 2, 7, 14, 6, 5, 9, 0, 11, 15, 8, 1 }, { 10, 7, 12, 9, 14, 3, 13, 15, 4, 0, 11, 2, 5, 8, 1, 6 },
    { 12, 13, 9, 11, 15, 10, 14, 8, 7, 2, 5, 3, 0, 1, 6, 4 }, { 9, 14, 11, 5, 8, 12, 15, 1, 13, 3, 0, 10, 2, 6, 4, 7 },
    { 11, 15, 5, 0, 1, 9, 8, 6, 14, 10, 2, 12, 3, 4, 7, 13 },
};

constexpr uint8_t MSG_SCHEDULE_BLAKE2[10][16] = {
    { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 }, { 14, 10, 4, 8, 9, 15, 13, 6, 1, 12, 0, 2, 11, 7, 5, 3 },
    { 11, 8, 12, 0, 5, 2, 15, 13, 10, 14, 3, 6, 7, 1, 9, 4 }, { 7, 9, 3, 1, 13, 12, 11, 14, 2, 6, 5, 10, 4, 0, 15, 8 },
    { 9, 0, 5, 7, 2, 4, 10, 15, 14, 1, 11, 12, 6, 8, 3, 13 }, { 2, 12, 6, 10, 0, 11, 8, 3, 4, 13, 7, 5, 15, 14, 1, 9 },
    { 12, 5, 1, 15, 14, 13, 4, 10, 0, 7, 6, 3, 9, 2, 8, 11 }, { 13, 11, 7, 14, 12, 1, 3, 9, 5, 0, 15, 4, 8, 6, 2, 10 },
    { 6, 15, 14, 9, 11, 3, 0, 8, 12, 2, 13, 7, 1, 4, 10, 5 }, { 10, 2, 8, 4, 7, 6, 1, 5, 15, 11, 9, 14, 3, 12, 13, 0 },
};

/**
 * @brief The `G` mixing function of Blake2s and Blake3s: four modular additions interleaved with four
 * XOR-and-rotate-right steps.
 *
 * @details Cost per call (UltraCircuitBuilder):
 *
 *   4 XOR-rotates x 4 lookup gates (32-bit words in 8-bit slices)  = 16
 *   2 three-term modular additions x 2 gates                       =  4
 *   2 two-term modular additions x 1 gate                          =  2
 *                                                                    --
 *                                                                    22
 * plus four range constraints on the addition overflows, which amortize to well under a gate each since all of Blake
 * shares two range lists (1-bit and 2-bit).
 *
 * +-----------+--------------+-----------------------+---------------------------+
 * |           |  calls to G  | gate count for rounds | lookups outside the rounds|
 * |-----------|--------------|-----------------------|---------------------------|
 * |  Blake2s  |      80      |        80 * 22        |          20 * 4           |
 * |  Blake3s  |      56      |        56 * 22        |           8 * 4           |
 * +-----------+--------------+-----------------------+---------------------------+
 *
 * Every value handed to a lookup here must be a reduced 32-bit word: the 8-bit slicing leaves no headroom above
 * bit 31, so an out-of-range key has no valid decomposition and the lookup simply fails. That is why each addition
 * is reduced by `add_normalize_unsafe` rather than left to overflow. The reduction is sound despite its name: the
 * overflow witness is range-constrained, and the subsequent lookup pins the result below 2^32, which leaves the
 * prover exactly one admissible choice of overflow.
 *
 * Inputs: - A pointer to a 16-word `state`,
 *         - indices a, b, c, d,
 *         - addition messages x and y
 */
template <typename Builder>
void g(field_t<Builder> state[BLAKE_STATE_SIZE],
       size_t a,
       size_t b,
       size_t c,
       size_t d,
       field_t<Builder> x,
       field_t<Builder> y)
{
    using plookup_read_pt = plookup_read<Builder>;

    // For simplicity, state[a] is written as `a' in comments.
    // a = a + b + x. Three 32-bit summands overflow by at most 2 bits.
    state[a] = hash_utils::add_normalize_unsafe(state[a], state[b] + x, /*overflow_bits=*/2);

    // d = (d ^ a).ror(16)
    state[d] =
        plookup_read_pt::get_lookup_accumulators(BLAKE_XOR_ROTATE_16, state[d], state[a], true)[ColumnIdx::C3][0];

    // c = c + d
    state[c] = hash_utils::add_normalize_unsafe(state[c], state[d], /*overflow_bits=*/1);

    // b = (b ^ c).ror(12)
    state[b] =
        plookup_read_pt::get_lookup_accumulators(BLAKE_XOR_ROTATE_12, state[b], state[c], true)[ColumnIdx::C3][0];

    // a = a + b + y
    state[a] = hash_utils::add_normalize_unsafe(state[a], state[b] + y, /*overflow_bits=*/2);

    // d = (d ^ a).ror(8)
    state[d] = plookup_read_pt::get_lookup_accumulators(BLAKE_XOR_ROTATE_8, state[d], state[a], true)[ColumnIdx::C3][0];

    // c = c + d
    state[c] = hash_utils::add_normalize_unsafe(state[c], state[d], /*overflow_bits=*/1);

    // b = (b ^ c).ror(7)
    state[b] = plookup_read_pt::get_lookup_accumulators(BLAKE_XOR_ROTATE_7, state[b], state[c], true)[ColumnIdx::C3][0];
}

/*
 * This is the round function used in Blake2s and Blake3s for Ultra.
 * Inputs: - 16-word state
 *         - 16-word msg
 *         - round number
 *         - which_blake to choose Blake2 or Blake3 (false -> Blake2)
 */
template <typename Builder>
void round_fn(field_t<Builder> state[BLAKE_STATE_SIZE],
              field_t<Builder> msg[BLAKE_STATE_SIZE],
              size_t round,
              const bool which_blake = false)
{
    // Select the message schedule based on the round.
    const uint8_t* schedule = which_blake ? MSG_SCHEDULE_BLAKE3[round] : MSG_SCHEDULE_BLAKE2[round];

    // Mix the columns.
    g<Builder>(state, 0, 4, 8, 12, msg[schedule[0]], msg[schedule[1]]);
    g<Builder>(state, 1, 5, 9, 13, msg[schedule[2]], msg[schedule[3]]);
    g<Builder>(state, 2, 6, 10, 14, msg[schedule[4]], msg[schedule[5]]);
    g<Builder>(state, 3, 7, 11, 15, msg[schedule[6]], msg[schedule[7]]);

    // Mix the rows.
    g<Builder>(state, 0, 5, 10, 15, msg[schedule[8]], msg[schedule[9]]);
    g<Builder>(state, 1, 6, 11, 12, msg[schedule[10]], msg[schedule[11]]);
    g<Builder>(state, 2, 7, 8, 13, msg[schedule[12]], msg[schedule[13]]);
    g<Builder>(state, 3, 4, 9, 14, msg[schedule[14]], msg[schedule[15]]);
}

} // namespace bb::stdlib::blake_util
