#pragma once

#include "barretenberg/numeric/bitop/get_msb.hpp"
#include "barretenberg/zcash/halo2/gates.hpp"
#include "barretenberg/zcash/halo2/plonkish_builder.hpp"
#include "barretenberg/zcash/primitives/sinsemilla.hpp"

#include <string>
#include <vector>

/**
 * @file trace.hpp
 * @brief The Honk-anchored execution trace of a halo2-style table, and a native checker of all its constraints.
 *
 * @details halo2 row r of the table lands at trace row `r + row_offset`. The first `row_offset` rows are reserved for
 * zero rows and ZK masking (see OrchardFlavor::TRACE_OFFSET). The selector of a gate with ANCHOR_SHIFT = 1 (a gate that
 * queries Rotation::prev()) is moved one row up so that all of its advice queries become offsets 0, 1, 2 from the
 * selector row. Lookup table columns hold the Sinsemilla S table (idx, x, y) in rows [row_offset, row_offset + 1024)
 * marked by `q_table`; the 10-bit range table is table_idx on the same rows.
 */
namespace bb::zcash::halo2 {

inline constexpr size_t SINSEMILLA_TABLE_SIZE = size_t{ 1 } << 10;

template <typename Cycle> struct AnchoredTrace {
    using FF = typename Cycle::FF;
    using Table = typename PlonkishBuilder<FF>::Table;

    size_t num_rows = 0;   // power of two
    size_t row_offset = 0; // trace row of halo2 row 0
    std::array<std::vector<FF>, NUM_ADVICE> advice;
    std::array<std::vector<FF>, NUM_FIXED> fixed;
    std::array<std::vector<FF>, NUM_SELECTORS> selectors;
    std::array<std::vector<FF>, 3> table; // idx, x, y
    std::vector<FF> q_table;
    // Copy constraints and public-input cells in trace coordinates. Fixed column f0 is the only fixed column taking
    // part in copies; columns are numbered a0..a9 -> 0..9, f0 -> 10, f5..f7 -> 11..13 (see permutation_column).
    std::vector<std::pair<std::pair<size_t, size_t>, std::pair<size_t, size_t>>> copies;
    std::vector<std::pair<size_t, size_t>> public_input_cells;
    std::vector<FF> public_inputs;

    // Index of a halo2 column in the 14 permutation columns.
    static size_t permutation_column(const Column& c)
    {
        if (c.kind == ColumnKind::Advice) {
            return c.index;
        }
        BB_ASSERT(c.kind == ColumnKind::Fixed);
        if (c.index == F0) {
            return NUM_ADVICE;
        }
        BB_ASSERT(c.index >= F5 && c.index <= F7, "fixed column is not equality-enabled");
        return NUM_ADVICE + 1 + (c.index - F5);
    }

    // The value of permutation column `col` at `row`.
    const FF& permutation_value(size_t col, size_t row) const
    {
        if (col < NUM_ADVICE) {
            return advice[col][row];
        }
        if (col == NUM_ADVICE) {
            return fixed[F0][row];
        }
        return fixed[F5 + (col - NUM_ADVICE - 1)][row];
    }

    /**
     * @param min_tail_rows rows that must stay free after the last used row (the gates read up to two rows ahead).
     */
    static AnchoredTrace build(const Table& t, size_t row_offset, size_t min_tail_rows = 2)
    {
        AnchoredTrace out;
        out.row_offset = row_offset;
        const size_t used = row_offset + std::max(t.num_rows, SINSEMILLA_TABLE_SIZE) + min_tail_rows;
        out.num_rows = size_t{ 1 } << numeric::get_msb(used - 1) << 1;
        if ((out.num_rows >> 1) >= used) {
            out.num_rows >>= 1;
        }
        const size_t n = out.num_rows;
        for (auto& col : out.advice) {
            col.assign(n, FF(0));
        }
        for (auto& col : out.fixed) {
            col.assign(n, FF(0));
        }
        for (auto& col : out.selectors) {
            col.assign(n, FF(0));
        }
        for (auto& col : out.table) {
            col.assign(n, FF(0));
        }
        out.q_table.assign(n, FF(0));

        for (size_t c = 0; c < NUM_ADVICE; ++c) {
            for (size_t r = 0; r < t.num_rows; ++r) {
                out.advice[c][r + row_offset] = t.advice[c][r];
            }
        }
        for (size_t c = 0; c < NUM_FIXED; ++c) {
            for (size_t r = 0; r < t.num_rows; ++r) {
                out.fixed[c][r + row_offset] = t.fixed[c][r];
            }
        }
        for (size_t s = 0; s < NUM_SELECTORS; ++s) {
            for (size_t r = 0; r < t.num_rows; ++r) {
                if (t.selectors[s][r] != 0) {
                    BB_ASSERT_GTE(r + row_offset, ANCHOR_SHIFT[s]);
                    out.selectors[s][r + row_offset - ANCHOR_SHIFT[s]] = FF(1);
                }
            }
        }
        const auto& S = Sinsemilla<Cycle>::S_table();
        for (size_t i = 0; i < SINSEMILLA_TABLE_SIZE; ++i) {
            out.table[0][row_offset + i] = FF(i);
            out.table[1][row_offset + i] = S[i].x;
            out.table[2][row_offset + i] = S[i].y;
            out.q_table[row_offset + i] = FF(1);
        }
        for (const auto& [a, b] : t.copies) {
            out.copies.emplace_back(std::make_pair(permutation_column(a.first), a.second + row_offset),
                                    std::make_pair(permutation_column(b.first), b.second + row_offset));
        }
        for (const auto& [col, row] : t.public_input_cells) {
            out.public_input_cells.emplace_back(permutation_column(col), row + row_offset);
        }
        out.public_inputs = t.public_inputs;
        return out;
    }

    // Native accessor at one anchor row.
    struct RowView {
        const AnchoredTrace& t;
        size_t row;
        FF adv(size_t col, size_t off) const { return t.advice[col][row + off]; }
        FF fix(size_t col) const { return t.fixed[col][row]; }
        FF sel(size_t s) const { return t.selectors[s][row]; }
        bool active(size_t s) const { return !t.selectors[s][row].is_zero(); }
    };

    /**
     * @brief Checks every gate, lookup and copy constraint.
     * @return a description of the first few failures (empty if the trace is valid).
     */
    std::vector<std::string> check() const
    {
        std::vector<std::string> failures;
        auto fail = [&](std::string msg) {
            if (failures.size() < 20) {
                failures.push_back(std::move(msg));
            }
        };
        for (size_t row = 0; row + 2 < num_rows; ++row) {
            RowView view{ *this, row };
            OrchardGates<Cycle>::template evaluate<FF>(view, [&](size_t idx, const FF& value) {
                if (!value.is_zero()) {
                    fail("gate constraint " + std::to_string(idx) + " (selector " +
                         std::to_string(CONSTRAINT_SELECTOR[idx]) + ") fails at row " + std::to_string(row));
                }
            });
            if (!selectors[Q_LOOKUP][row].is_zero()) {
                const uint256_t v(OrchardGates<Cycle>::template range_lookup_value<FF>(view));
                if (v >= SINSEMILLA_TABLE_SIZE) {
                    fail("range lookup fails at row " + std::to_string(row));
                }
            }
            const auto& S = Sinsemilla<Cycle>::S_table();
            for (size_t chip = 0; chip < 2; ++chip) {
                const size_t q = (chip == 0) ? Q_SINSEMILLA1_1 : Q_SINSEMILLA1_2;
                if (!selectors[q][row].is_zero()) {
                    const auto [m, x, y] = OrchardGates<Cycle>::template sinsemilla_lookup_value<FF>(view, chip);
                    const uint256_t mi(m);
                    if (mi >= SINSEMILLA_TABLE_SIZE || S[mi.data[0]].x != x || S[mi.data[0]].y != y) {
                        fail("sinsemilla lookup " + std::to_string(chip) + " fails at row " + std::to_string(row));
                    }
                }
            }
        }
        for (const auto& [a, b] : copies) {
            if (permutation_value(a.first, a.second) != permutation_value(b.first, b.second)) {
                fail("copy constraint fails between column " + std::to_string(a.first) + " row " +
                     std::to_string(a.second) + " and column " + std::to_string(b.first) + " row " +
                     std::to_string(b.second));
            }
        }
        for (size_t i = 0; i < public_input_cells.size(); ++i) {
            const auto& [col, row] = public_input_cells[i];
            if (permutation_value(col, row) != public_inputs[i]) {
                fail("public input " + std::to_string(i) + " mismatch");
            }
        }
        return failures;
    }
};

} // namespace bb::zcash::halo2
