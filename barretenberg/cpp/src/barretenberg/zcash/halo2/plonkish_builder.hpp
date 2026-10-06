#pragma once

#include "barretenberg/common/assert.hpp"
#include "barretenberg/common/throw_or_abort.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <map>
#include <numeric>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace bb::zcash::halo2 {

/**
 * @brief A column of a PLONKish (halo2-style) table.
 * @details Advice columns hold witness values, fixed columns hold circuit constants, selector columns are boolean fixed
 * columns that gate custom constraints. Table columns of a lookup argument are not represented here (they are static).
 */
enum class ColumnKind : uint8_t { Advice, Fixed, Selector };

struct Column {
    ColumnKind kind = ColumnKind::Advice;
    uint16_t index = 0;

    bool operator==(const Column& other) const = default;
    auto operator<=>(const Column& other) const = default;
};

inline Column advice(size_t i)
{
    return { ColumnKind::Advice, static_cast<uint16_t>(i) };
}
inline Column fixed(size_t i)
{
    return { ColumnKind::Fixed, static_cast<uint16_t>(i) };
}
inline Column selector(size_t i)
{
    return { ColumnKind::Selector, static_cast<uint16_t>(i) };
}

/**
 * @brief A cell of a region, addressed relative to the region's (not yet known) start row.
 */
struct CellRef {
    uint32_t region = 0;
    Column column;
    uint32_t row = 0;
};

/**
 * @brief A cell together with the value assigned to it (halo2's AssignedCell).
 */
template <typename FF> struct AssignedCell {
    CellRef cell;
    FF value;
};

/**
 * @brief Builder for a halo2-style circuit: regions laid out by a first-fit floor planner, with copy constraints,
 * constants and public inputs.
 *
 * @details Gadgets assign cells inside regions at offsets relative to the region start, exactly like halo2's `Region`.
 * Synthesis is single-pass; after all regions are recorded `finalize()` runs the floor planner (halo2's V1 strategy:
 * regions sorted by advice area, largest first, each placed at the lowest row at which all of its columns are free),
 * places constants in free rows of the constants column, and places public inputs in free rows of the public-input
 * column. The resulting table is consumed by the Honk flavor (`OrchardFlavor`) and by the native constraint checker.
 */
template <typename FF_> class PlonkishBuilder {
  public:
    using FF = FF_;
    using Cell = AssignedCell<FF>;

    struct Config {
        size_t num_advice = 0;
        size_t num_fixed = 0;
        size_t num_selectors = 0;
        // Fixed column holding constants for `assign_advice_from_constant` / `constrain_constant`.
        size_t constants_column = 0;
        // Advice column whose free rows hold public inputs.
        size_t public_input_column = 0;
    };

    struct RegionRecord {
        std::string name;
        std::vector<std::tuple<Column, uint32_t, FF>> assignments; // advice and fixed values
        std::vector<std::pair<uint16_t, uint32_t>> selectors;      // (selector, row)
        std::vector<Column> columns;                               // distinct columns used
        uint32_t height = 0;
        uint32_t start = 0; // set by the floor planner
    };

    class Region {
      public:
        Region(PlonkishBuilder& builder, uint32_t id)
            : builder_(builder)
            , id_(id)
        {}

        Cell assign_advice(Column col, size_t offset, const FF& value)
        {
            BB_ASSERT(col.kind == ColumnKind::Advice);
            return assign(col, offset, value);
        }

        Cell assign_fixed(Column col, size_t offset, const FF& value)
        {
            BB_ASSERT(col.kind == ColumnKind::Fixed);
            return assign(col, offset, value);
        }

        void enable_selector(Column sel, size_t offset)
        {
            BB_ASSERT(sel.kind == ColumnKind::Selector);
            auto& r = record();
            r.selectors.emplace_back(sel.index, static_cast<uint32_t>(offset));
            touch(sel, offset);
        }

        // Assign `src`'s value to (col, offset) and constrain the two cells to be equal.
        Cell copy_advice(const Cell& src, Column col, size_t offset)
        {
            Cell dst = assign_advice(col, offset, src.value);
            builder_.constrain_equal(src, dst);
            return dst;
        }

        Cell assign_advice_from_constant(Column col, size_t offset, const FF& constant)
        {
            Cell dst = assign_advice(col, offset, constant);
            builder_.constrain_constant(dst, constant);
            return dst;
        }

        Cell assign_advice_from_instance(size_t instance_row, Column col, size_t offset)
        {
            BB_ASSERT_LT(instance_row, builder_.public_inputs_.size());
            Cell dst = assign_advice(col, offset, builder_.public_inputs_[instance_row]);
            builder_.constrain_instance(dst, instance_row);
            return dst;
        }

        uint32_t id() const { return id_; }

      private:
        PlonkishBuilder& builder_;
        uint32_t id_;

        RegionRecord& record() { return builder_.regions_[id_]; }

        void touch(Column col, size_t offset)
        {
            auto& r = record();
            r.height = std::max(r.height, static_cast<uint32_t>(offset + 1));
            if (std::find(r.columns.begin(), r.columns.end(), col) == r.columns.end()) {
                r.columns.push_back(col);
            }
        }

        Cell assign(Column col, size_t offset, const FF& value)
        {
            auto& r = record();
            r.assignments.emplace_back(col, static_cast<uint32_t>(offset), value);
            touch(col, offset);
            return Cell{ CellRef{ id_, col, static_cast<uint32_t>(offset) }, value };
        }
    };

    explicit PlonkishBuilder(Config config)
        : config_(config)
    {}

    const Config& config() const { return config_; }

    // Public inputs must be declared before synthesis so that regions can read them.
    void set_public_inputs(std::vector<FF> public_inputs) { public_inputs_ = std::move(public_inputs); }
    const std::vector<FF>& public_inputs() const { return public_inputs_; }

    template <typename Fn> auto assign_region(const std::string& name, Fn&& fn)
    {
        const auto id = static_cast<uint32_t>(regions_.size());
        regions_.push_back(RegionRecord{ .name = name });
        Region region(*this, id);
        return fn(region);
    }

    void constrain_equal(const Cell& a, const Cell& b)
    {
        BB_ASSERT(a.value == b.value, "copy constraint between unequal cells");
        copies_.emplace_back(a.cell, b.cell);
    }

    void constrain_constant(const Cell& cell, const FF& constant)
    {
        BB_ASSERT(cell.value == constant, "constant constraint on a cell with another value");
        constants_.emplace_back(cell.cell, constant);
    }

    void constrain_instance(const Cell& cell, size_t instance_row)
    {
        BB_ASSERT_LT(instance_row, public_inputs_.size());
        BB_ASSERT(cell.value == public_inputs_[instance_row], "instance constraint on a cell with another value");
        instance_constraints_.emplace_back(cell.cell, instance_row);
    }

    /**
     * @brief The materialized table.
     * @details `num_rows` is the number of used rows (before any padding to a power of two). Cell (column, row) of the
     * table is at absolute row `row`; consumers add their own row offset.
     */
    struct Table {
        size_t num_rows = 0;
        std::vector<std::vector<FF>> advice;
        std::vector<std::vector<FF>> fixed;
        std::vector<std::vector<uint8_t>> selectors;
        // Copy constraints between absolute cells (only advice and fixed columns).
        std::vector<std::pair<std::pair<Column, uint32_t>, std::pair<Column, uint32_t>>> copies;
        // Absolute (column, row) of each public input.
        std::vector<std::pair<Column, uint32_t>> public_input_cells;
        std::vector<FF> public_inputs;
    };

    Table finalize()
    {
        place_regions();
        Table t;
        size_t num_rows = 0;
        for (const auto& r : regions_) {
            num_rows = std::max(num_rows, static_cast<size_t>(r.start + r.height));
        }
        // Constants: one fresh cell per constant in a free row of the constants column (halo2 V1 behaviour).
        const Column const_col = fixed(config_.constants_column);
        std::vector<std::pair<uint32_t, FF>> constant_cells;
        for (const auto& [cell, value] : constants_) {
            const uint32_t row = allocate_free_row(const_col);
            constant_cells.emplace_back(row, value);
            num_rows = std::max(num_rows, static_cast<size_t>(row + 1));
        }
        // Public inputs: one cell per public input in a free row of the public-input column.
        const Column pi_col = advice(config_.public_input_column);
        std::vector<uint32_t> pi_rows;
        for (size_t i = 0; i < public_inputs_.size(); ++i) {
            const uint32_t row = allocate_free_row(pi_col);
            pi_rows.push_back(row);
            num_rows = std::max(num_rows, static_cast<size_t>(row + 1));
        }
        t.num_rows = num_rows;
        t.advice.assign(config_.num_advice, std::vector<FF>(num_rows, FF(0)));
        t.fixed.assign(config_.num_fixed, std::vector<FF>(num_rows, FF(0)));
        t.selectors.assign(config_.num_selectors, std::vector<uint8_t>(num_rows, 0));

        for (const auto& r : regions_) {
            for (const auto& [col, row, value] : r.assignments) {
                auto& dst = (col.kind == ColumnKind::Advice) ? t.advice[col.index] : t.fixed[col.index];
                dst[r.start + row] = value;
            }
            for (const auto& [sel, row] : r.selectors) {
                t.selectors[sel][r.start + row] = 1;
            }
        }
        auto absolute = [&](const CellRef& c) { return std::make_pair(c.column, regions_[c.region].start + c.row); };
        for (const auto& [a, b] : copies_) {
            t.copies.emplace_back(absolute(a), absolute(b));
        }
        for (size_t i = 0; i < constants_.size(); ++i) {
            const auto& [row, value] = constant_cells[i];
            t.fixed[config_.constants_column][row] = value;
            t.copies.emplace_back(absolute(constants_[i].first), std::make_pair(const_col, row));
        }
        for (size_t i = 0; i < public_inputs_.size(); ++i) {
            t.advice[config_.public_input_column][pi_rows[i]] = public_inputs_[i];
            t.public_input_cells.emplace_back(pi_col, pi_rows[i]);
        }
        for (const auto& [cell, idx] : instance_constraints_) {
            t.copies.emplace_back(absolute(cell), std::make_pair(pi_col, pi_rows[idx]));
        }
        t.public_inputs = public_inputs_;
        return t;
    }

    const std::vector<RegionRecord>& regions() const { return regions_; }

  private:
    Config config_;
    std::vector<RegionRecord> regions_;
    std::vector<std::pair<CellRef, CellRef>> copies_;
    std::vector<std::pair<CellRef, FF>> constants_;
    std::vector<std::pair<CellRef, size_t>> instance_constraints_;
    std::vector<FF> public_inputs_;

    // Per-column row occupancy bitmaps (grown on demand).
    std::map<Column, std::vector<uint8_t>> occupied_;

    // Returns true if rows [start, end) of `col` are free; otherwise sets `next_candidate` past the first conflict.
    bool is_free(const Column& col, uint32_t start, uint32_t end, uint32_t& next_candidate) const
    {
        auto it = occupied_.find(col);
        if (it == occupied_.end()) {
            return true;
        }
        const auto& bits = it->second;
        const uint32_t limit = std::min<uint32_t>(end, static_cast<uint32_t>(bits.size()));
        for (uint32_t row = start; row < limit; ++row) {
            if (bits[row] != 0) {
                uint32_t skip = row + 1;
                while (skip < bits.size() && bits[skip] != 0) {
                    ++skip;
                }
                next_candidate = std::max(next_candidate, skip);
                return false;
            }
        }
        return true;
    }

    void occupy(const Column& col, uint32_t start, uint32_t end)
    {
        auto& bits = occupied_[col];
        if (bits.size() < end) {
            bits.resize(std::max<size_t>(end, bits.size() * 2), 0);
        }
        std::fill(bits.begin() + start, bits.begin() + end, uint8_t{ 1 });
    }

    uint32_t allocate_free_row(const Column& col)
    {
        uint32_t row = 0;
        while (true) {
            uint32_t next = row;
            if (is_free(col, row, row + 1, next)) {
                occupy(col, row, row + 1);
                return row;
            }
            row = next;
        }
    }

    void place_regions()
    {
        std::vector<size_t> order(regions_.size());
        std::iota(order.begin(), order.end(), 0);
        auto advice_area = [&](size_t i) {
            const auto& r = regions_[i];
            const auto n = std::count_if(
                r.columns.begin(), r.columns.end(), [](const Column& c) { return c.kind == ColumnKind::Advice; });
            return static_cast<size_t>(n) * r.height;
        };
        std::stable_sort(
            order.begin(), order.end(), [&](size_t a, size_t b) { return advice_area(a) > advice_area(b); });
        for (size_t idx : order) {
            auto& r = regions_[idx];
            if (r.height == 0) {
                continue;
            }
            uint32_t start = 0;
            while (true) {
                uint32_t next = start;
                bool ok = true;
                for (const auto& col : r.columns) {
                    if (!is_free(col, start, start + r.height, next)) {
                        ok = false;
                    }
                }
                if (ok) {
                    break;
                }
                start = next;
            }
            r.start = start;
            for (const auto& col : r.columns) {
                occupy(col, start, start + r.height);
            }
        }
    }
};

} // namespace bb::zcash::halo2
