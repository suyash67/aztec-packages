#pragma once
/**
 * @file blake3vm_circuit_builder.hpp
 * @brief Trace generation for the Blake3VM, a standalone Honk VM that proves batches of Blake3s hashes.
 *
 * The VM exists to make recursive verification of hash-heavy proofs (e.g. WHIR-Honk with Blake3s Merkle
 * trees) affordable: instead of paying ~2,600 UltraHonk gates per Blake3s compression inside a recursive
 * verifier, the hashing is delegated to this VM, whose own proof is verified with sumcheck + Shplemini/KZG
 * cost independent of the number of hashes. See README.md for the full design; the trace layout summary:
 *
 * - Row 0 is the genesis boundary row (all-zero witness; carries compression 0's block parameters).
 * - Compression c occupies rows [1 + 60c, 60(c+1)]: 56 G-function rows (round r ∈ [0,7), position
 *   p ∈ [0,8)), then 4 output rows computing outᵢ = vᵢ ⊕ vᵢ₊₈ (two words per row).
 * - The 4th output row doubles as the boundary row: it constrains the next row's state to the next
 *   compression's initial state (IV or chained CV) and carries that compression's block_len/flags/chain.
 * - One ghost row after the last boundary absorbs its state-initialisation constraint.
 * - Rows [0, 2^16) of three precomputed columns hold the 8-bit XOR table (x, y, x ⊕ y); every byte-valued
 *   witness is constrained through log-derivative lookups into it.
 *
 * All 32-bit words are decomposed little-endian into 4 byte limbs. Rotations by 16/8 are byte
 * permutations; rotations by 12/7 split one byte at the rotation boundary (nibbles / 7-bit + msb).
 */

#include <array>
#include <cstdint>
#include <span>
#include <vector>

#include "barretenberg/common/assert.hpp"
#include "barretenberg/crypto/blake3s/blake3s.hpp"

namespace bb {

/** @brief Fixed layout constants shared by the trace builder, flavor, and relations. */
struct Blake3VMTraceLayout {
    static constexpr size_t NUM_ROUNDS = 7;
    static constexpr size_t NUM_G_PER_ROUND = 8;
    static constexpr size_t NUM_G_ROWS = NUM_ROUNDS * NUM_G_PER_ROUND; // 56
    static constexpr size_t NUM_OUT_ROWS = 4;
    static constexpr size_t ROWS_PER_COMPRESSION = NUM_G_ROWS + NUM_OUT_ROWS; // 60
    static constexpr size_t TABLE_ROWS = 1UL << 16;                           // 8-bit XOR table
    static constexpr size_t NUM_LOOKUP_SETS = 6;

    // State indices (a, b, c, d) used by the G function at each of the 8 positions within a round.
    static constexpr std::array<size_t, 8> POS_A = { 0, 1, 2, 3, 0, 1, 2, 3 };
    static constexpr std::array<size_t, 8> POS_B = { 4, 5, 6, 7, 5, 6, 7, 4 };
    static constexpr std::array<size_t, 8> POS_C = { 8, 9, 10, 11, 10, 11, 8, 9 };
    static constexpr std::array<size_t, 8> POS_D = { 12, 13, 14, 15, 15, 12, 13, 14 };

    // In-place message permutation applied between rounds: m'[j] = m[MSG_PERM[j]].
    // Equals blake3::MSG_SCHEDULE[1]; with it, round r position p always reads mx = m[2p], my = m[2p+1].
    static constexpr std::array<size_t, 16> MSG_PERM = { 2, 6, 3, 10, 7, 0, 4, 13, 1, 11, 12, 5, 9, 14, 15, 8 };

    static constexpr size_t row_of_g(size_t compression, size_t round, size_t pos)
    {
        return 1 + ROWS_PER_COMPRESSION * compression + NUM_G_PER_ROUND * round + pos;
    }
    static constexpr size_t row_of_out(size_t compression, size_t t)
    {
        return 1 + ROWS_PER_COMPRESSION * compression + NUM_G_ROWS + t;
    }
    // Boundary row carrying compression c's block parameters (genesis row for c = 0, else the
    // previous compression's 4th output row).
    static constexpr size_t boundary_row_of(size_t compression) { return ROWS_PER_COMPRESSION * compression; }
};

/** @brief Witness values for one G-function row. */
struct Blake3VMGRow {
    std::array<uint32_t, 16> v; // state before this G function
    uint32_t a1, c1, b1, a2, c2, b2;
    uint8_t ca1, cc1, ca2, cc2; // addition carries: ca* ∈ {0,1,2}, cc* ∈ {0,1}
    std::array<uint8_t, 4> din_b, bin_b, a1_b, z1_b, c1_b, z2_b, b1_b, a2_b, z3_b, c2_b, z4_b, mx_b, my_b;
    uint8_t z2_nib_lo, z2_nib_hi; // z2_b[1] = z2_nib_lo + 16·z2_nib_hi (rotation by 12)
    uint8_t z4_lo7, z4_msb;       // z4_b[0] = z4_lo7 + 128·z4_msb (rotation by 7)
    std::array<uint32_t, 16> m;   // message words in the current (permuted) order
};

/** @brief Data for one compression's output/boundary block. */
struct Blake3VMOutBlock {
    std::array<uint32_t, 16> v; // final state after round 7 (constant across the 4 output rows)
    std::array<uint32_t, 8> out;
};

/** @brief Block parameters carried on the boundary row preceding each compression. */
struct Blake3VMBoundary {
    uint32_t block_len = 0;
    uint32_t flags = 0;
    bool chain = false; // if set, initial CV = previous compression's output words
};

/**
 * @brief Builds the Blake3VM witness trace for a batch of Blake3s hashes.
 * @details Accepts the same inputs as bb's native `blake3::blake3s` (≤ 1024 bytes, single chunk,
 * 32-byte output) and reproduces its block splitting and flag schedule exactly. Each 64-byte block
 * becomes one compression; multi-block messages chain their CV through the boundary rows.
 */
class Blake3VMCircuitBuilder {
  public:
    using Layout = Blake3VMTraceLayout;

    std::vector<Blake3VMGRow> g_rows;         // 56 per compression
    std::vector<Blake3VMOutBlock> out_blocks; // 1 per compression
    std::vector<Blake3VMBoundary> boundaries; // boundaries[c] = params of compression c
    // Read multiplicities into the XOR table, one array per lookup set.
    std::array<std::vector<uint32_t>, Layout::NUM_LOOKUP_SETS> read_counts;

    Blake3VMCircuitBuilder()
    {
        for (auto& counts : read_counts) {
            counts.resize(Layout::TABLE_ROWS, 0);
        }
        // Genesis row reads: lookup sets 0 and 1 are active on every boundary/output-type row,
        // including row 0, where all byte columns are zero — 4 reads each of table row (0,0).
        read_counts[0][0] += 4;
        read_counts[1][0] += 4;
    }

    size_t num_compressions() const { return out_blocks.size(); }

    // Trace rows used, including genesis and ghost rows (excluding the XOR table floor).
    size_t num_active_rows() const { return 1 + Layout::ROWS_PER_COMPRESSION * num_compressions() + 1; }

    size_t min_circuit_size() const { return std::max(num_active_rows(), Layout::TABLE_ROWS); }

    /**
     * @brief Append the compressions for one Blake3s hash of `input` (≤ 1024 bytes).
     * @return The 32-byte digest, asserted equal to `blake3::blake3s(input)`.
     */
    std::array<uint8_t, 32> add_hash(std::span<const uint8_t> input)
    {
        BB_ASSERT_LTE(input.size(), 1024U, "Blake3VM supports single-chunk inputs of at most 1024 bytes");

        // Mirror blake3_hasher_update/finalize: full blocks are compressed while *more* input remains;
        // the final (possibly partial, possibly full) block is compressed with CHUNK_END | ROOT.
        const size_t num_blocks = input.size() <= blake3::BLAKE3_BLOCK_LEN
                                      ? 1
                                      : (input.size() + blake3::BLAKE3_BLOCK_LEN - 1) / blake3::BLAKE3_BLOCK_LEN;

        std::array<uint32_t, 8> cv;
        std::copy(blake3::IV.begin(), blake3::IV.end(), cv.begin());

        for (size_t b = 0; b < num_blocks; ++b) {
            const size_t offset = b * blake3::BLAKE3_BLOCK_LEN;
            const size_t len = std::min(input.size() - offset, static_cast<size_t>(blake3::BLAKE3_BLOCK_LEN));
            std::array<uint8_t, blake3::BLAKE3_BLOCK_LEN> block{};
            std::copy(input.begin() + static_cast<std::ptrdiff_t>(offset),
                      input.begin() + static_cast<std::ptrdiff_t>(offset + len),
                      block.begin());

            uint32_t flags = 0;
            if (b == 0) {
                flags |= blake3::CHUNK_START;
            }
            if (b + 1 == num_blocks) {
                flags |= blake3::CHUNK_END | blake3::ROOT;
            }

            Blake3VMBoundary boundary{ .block_len = static_cast<uint32_t>(len), .flags = flags, .chain = (b != 0) };
            const auto out = add_compression(cv, block, boundary);
            std::copy(out.begin(), out.begin() + 8, cv.begin());
        }

        std::array<uint8_t, 32> digest;
        for (size_t i = 0; i < 8; ++i) {
            blake3::store32(&digest[4 * i], cv[i]);
        }

        std::array<uint8_t, 32> expected;
        blake3::blake3s(input, expected);
        BB_ASSERT(digest == expected, "Blake3VM trace output disagrees with the blake3s reference");
        return digest;
    }

  private:
    static std::array<uint8_t, 4> to_bytes(uint32_t w)
    {
        return { static_cast<uint8_t>(w),
                 static_cast<uint8_t>(w >> 8),
                 static_cast<uint8_t>(w >> 16),
                 static_cast<uint8_t>(w >> 24) };
    }

    void count_read(size_t set, uint8_t x, uint8_t y) { read_counts[set][static_cast<size_t>(x) * 256 + y] += 1; }
    // Scaled range check: a k-bit value v is read as the table row (v·2^{8-k}, 0).
    void count_scaled_read(size_t set, uint8_t scaled) { count_read(set, scaled, 0); }

    /** @brief Run one compression, emitting 56 G rows + 1 output block, and return the 16 output words. */
    std::array<uint32_t, 16> add_compression(const std::array<uint32_t, 8>& cv,
                                             const std::array<uint8_t, blake3::BLAKE3_BLOCK_LEN>& block,
                                             const Blake3VMBoundary& boundary)
    {
        boundaries.push_back(boundary);

        std::array<uint32_t, 16> state = { cv[0],
                                           cv[1],
                                           cv[2],
                                           cv[3],
                                           cv[4],
                                           cv[5],
                                           cv[6],
                                           cv[7],
                                           blake3::IV[0],
                                           blake3::IV[1],
                                           blake3::IV[2],
                                           blake3::IV[3],
                                           0,
                                           0,
                                           boundary.block_len,
                                           boundary.flags };

        std::array<uint32_t, 16> m;
        for (size_t i = 0; i < 16; ++i) {
            m[i] = blake3::load32(&block[4 * i]);
        }

        // Reference state for the paranoid cross-check below.
        blake3::state_array ref_state;
        blake3::compress_pre(ref_state,
                             { cv[0], cv[1], cv[2], cv[3], cv[4], cv[5], cv[6], cv[7] },
                             block.data(),
                             static_cast<uint8_t>(boundary.block_len),
                             static_cast<uint8_t>(boundary.flags));

        for (size_t round = 0; round < Layout::NUM_ROUNDS; ++round) {
            for (size_t pos = 0; pos < Layout::NUM_G_PER_ROUND; ++pos) {
                add_g_row(state, m, pos);
            }
            if (round + 1 < Layout::NUM_ROUNDS) {
                std::array<uint32_t, 16> permuted;
                for (size_t j = 0; j < 16; ++j) {
                    permuted[j] = m[Layout::MSG_PERM[j]];
                }
                m = permuted;
            }
        }
        BB_ASSERT(state == ref_state, "Blake3VM G-row trace disagrees with the blake3s round function");

        Blake3VMOutBlock out_block;
        out_block.v = state;
        for (size_t i = 0; i < 8; ++i) {
            out_block.out[i] = state[i] ^ state[i + 8];
        }
        out_blocks.push_back(out_block);

        // Output-row lookups: row t computes out_{2t} (set 0) and out_{2t+1} (set 1), byte by byte.
        for (size_t t = 0; t < Layout::NUM_OUT_ROWS; ++t) {
            const auto x_even = to_bytes(state[2 * t]);
            const auto y_even = to_bytes(state[2 * t + 8]);
            const auto x_odd = to_bytes(state[2 * t + 1]);
            const auto y_odd = to_bytes(state[2 * t + 9]);
            for (size_t i = 0; i < 4; ++i) {
                count_read(0, x_even[i], y_even[i]);
                count_read(1, x_odd[i], y_odd[i]);
            }
        }

        std::array<uint32_t, 16> full_out;
        for (size_t i = 0; i < 8; ++i) {
            full_out[i] = state[i] ^ state[i + 8];
            full_out[i + 8] = state[i + 8] ^ cv[i];
        }
        return full_out;
    }

    void add_g_row(std::array<uint32_t, 16>& state, const std::array<uint32_t, 16>& m, size_t pos)
    {
        const uint32_t a = state[Layout::POS_A[pos]];
        const uint32_t b = state[Layout::POS_B[pos]];
        const uint32_t c = state[Layout::POS_C[pos]];
        const uint32_t d = state[Layout::POS_D[pos]];
        const uint32_t mx = m[2 * pos];
        const uint32_t my = m[2 * pos + 1];

        Blake3VMGRow row;
        row.v = state;
        row.m = m;
        row.din_b = to_bytes(d);
        row.bin_b = to_bytes(b);
        row.mx_b = to_bytes(mx);
        row.my_b = to_bytes(my);

        const uint64_t t1 = static_cast<uint64_t>(a) + b + mx;
        row.ca1 = static_cast<uint8_t>(t1 >> 32);
        row.a1 = static_cast<uint32_t>(t1);
        row.a1_b = to_bytes(row.a1);

        const uint32_t z1 = d ^ row.a1;
        row.z1_b = to_bytes(z1);
        const uint32_t d1 = blake3::rotr32(z1, 16);

        const uint64_t t2 = static_cast<uint64_t>(c) + d1;
        row.cc1 = static_cast<uint8_t>(t2 >> 32);
        row.c1 = static_cast<uint32_t>(t2);
        row.c1_b = to_bytes(row.c1);

        const uint32_t z2 = b ^ row.c1;
        row.z2_b = to_bytes(z2);
        row.z2_nib_lo = row.z2_b[1] & 0x0F;
        row.z2_nib_hi = row.z2_b[1] >> 4;
        row.b1 = blake3::rotr32(z2, 12);
        row.b1_b = to_bytes(row.b1);

        const uint64_t t3 = static_cast<uint64_t>(row.a1) + row.b1 + my;
        row.ca2 = static_cast<uint8_t>(t3 >> 32);
        row.a2 = static_cast<uint32_t>(t3);
        row.a2_b = to_bytes(row.a2);

        const uint32_t z3 = d1 ^ row.a2;
        row.z3_b = to_bytes(z3);
        const uint32_t d2 = blake3::rotr32(z3, 8);

        const uint64_t t4 = static_cast<uint64_t>(row.c1) + d2;
        row.cc2 = static_cast<uint8_t>(t4 >> 32);
        row.c2 = static_cast<uint32_t>(t4);
        row.c2_b = to_bytes(row.c2);

        const uint32_t z4 = row.b1 ^ row.c2;
        row.z4_b = to_bytes(z4);
        row.z4_lo7 = row.z4_b[0] & 0x7F;
        row.z4_msb = row.z4_b[0] >> 7;
        row.b2 = blake3::rotr32(z4, 7);

        // The relations recombine rotated words from byte limbs; check both algebraic identities here.
        BB_ASSERT_EQ(row.b1,
                     static_cast<uint32_t>(row.z2_nib_hi) + (static_cast<uint32_t>(row.z2_b[2]) << 4) +
                         (static_cast<uint32_t>(row.z2_b[3]) << 12) + (static_cast<uint32_t>(row.z2_b[0]) << 20) +
                         (static_cast<uint32_t>(row.z2_nib_lo) << 28));
        BB_ASSERT_EQ(row.b2,
                     static_cast<uint32_t>(row.z4_msb) + (static_cast<uint32_t>(row.z4_b[1]) << 1) +
                         (static_cast<uint32_t>(row.z4_b[2]) << 9) + (static_cast<uint32_t>(row.z4_b[3]) << 17) +
                         (static_cast<uint32_t>(row.z4_lo7) << 25));

        // Lookup reads, exactly mirroring the relation read terms.
        for (size_t i = 0; i < 4; ++i) {
            count_read(0, row.din_b[i], row.a1_b[i]);
            count_read(1, row.bin_b[i], row.c1_b[i]);
            count_read(2, row.z1_b[(i + 2) % 4], row.a2_b[i]);
            count_read(3, row.b1_b[i], row.c2_b[i]);
        }
        count_scaled_read(4, static_cast<uint8_t>(row.z2_nib_lo << 4));
        count_scaled_read(4, static_cast<uint8_t>(row.z2_nib_hi << 4));
        count_scaled_read(4, static_cast<uint8_t>(row.z4_lo7 << 1));
        count_read(4, row.mx_b[0], row.mx_b[1]);
        count_read(5, row.mx_b[2], row.mx_b[3]);
        count_read(5, row.my_b[0], row.my_b[1]);
        count_read(5, row.my_b[2], row.my_b[3]);

        g_rows.push_back(row);

        state[Layout::POS_A[pos]] = row.a2;
        state[Layout::POS_B[pos]] = row.b2;
        state[Layout::POS_C[pos]] = row.c2;
        state[Layout::POS_D[pos]] = static_cast<uint32_t>(d2);
    }
};

} // namespace bb
