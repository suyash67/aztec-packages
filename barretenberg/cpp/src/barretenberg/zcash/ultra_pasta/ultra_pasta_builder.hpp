#pragma once

#include "barretenberg/ecc/curves/pasta/pasta.hpp"
#include "barretenberg/honk/execution_trace/execution_trace_block.hpp"
#include "barretenberg/honk/execution_trace/ultra_execution_trace.hpp"
#include "barretenberg/stdlib_circuit_builders/ultra_circuit_builder.hpp"

#include <map>

/**
 * @file ultra_pasta_builder.hpp
 * @brief The Ultra circuit builder over the Pallas base field (the Vesta scalar field), so that Ultra circuits can be
 * committed to with Vesta.
 */
namespace bb::zcash {

/**
 * @brief The Ultra execution trace blocks with Pallas base field entries. Same blocks, in the same order, as
 * `UltraTraceBlockData` / `UltraExecutionTraceBlocks`.
 */
class UltraPastaExecutionTraceBlocks {
  public:
    using FF = pasta::fp;
    using BlockBase = ExecutionTraceBlock<FF, 4>;
    static constexpr size_t NUM_WIRES = BlockBase::NUM_WIRES;
    static constexpr size_t TRACE_OFFSET = UltraExecutionTraceBlocks::TRACE_OFFSET;
    static constexpr size_t NUM_BLOCKS = 9;

    BlockBase pub_inputs{};
    BlockBase lookup{ GateKind::Lookup };
    BlockBase arithmetic{ GateKind::Arith };
    BlockBase delta_range{ GateKind::DeltaRange };
    BlockBase elliptic{ GateKind::Elliptic };
    BlockBase memory{ GateKind::Memory };
    BlockBase nnf{ GateKind::Nnf };
    BlockBase poseidon2_external{ GateKind::Poseidon2Ext };
    BlockBase poseidon2_internal{ GateKind::Poseidon2Int };

    std::vector<std::string_view> get_labels() const
    {
        return {
            "pub_inputs", "lookup", "arithmetic",         "delta_range",        "elliptic",
            "memory",     "nnf",    "poseidon2_external", "poseidon2_internal",
        };
    }

    auto get()
    {
        return RefArray(std::array<BlockBase*, NUM_BLOCKS>{ &pub_inputs,
                                                            &lookup,
                                                            &arithmetic,
                                                            &delta_range,
                                                            &elliptic,
                                                            &memory,
                                                            &nnf,
                                                            &poseidon2_external,
                                                            &poseidon2_internal });
    }

    auto get() const
    {
        return RefArray(std::array<const BlockBase*, NUM_BLOCKS>{ &pub_inputs,
                                                                  &lookup,
                                                                  &arithmetic,
                                                                  &delta_range,
                                                                  &elliptic,
                                                                  &memory,
                                                                  &nnf,
                                                                  &poseidon2_external,
                                                                  &poseidon2_internal });
    }

    void summarize() const
    {
        info("Gate blocks summary:");
        const auto labels = get_labels();
        size_t i = 0;
        for (const auto& block : get()) {
            info(labels[i++], ": ", block.size());
        }
    }

    void compute_offsets(size_t trace_offset = TRACE_OFFSET)
    {
        auto offset = static_cast<uint32_t>(trace_offset + NUM_ZERO_ROWS);
        for (auto& block : get()) {
            block.trace_offset_ = offset;
            offset += static_cast<uint32_t>(block.size());
        }
    }

    size_t get_total_content_size()
    {
        size_t total_size = 0;
        for (const auto& block : get()) {
            total_size += block.size();
        }
        return total_size;
    }

    bool operator==(const UltraPastaExecutionTraceBlocks& other) const = default;
};

/**
 * @brief The Ultra circuit builder over F_p, with lookup tables whose entries are elements of F_p.
 * @details barretenberg's plookup::BasicTable stores BN254 scalar field entries. A table registered here keeps its F_p
 * rows in the builder, and its BasicTable holds only the row index (in the first column, zeros elsewhere), which is
 * all that the read-count bookkeeping of the prover uses. The table polynomials of UltraPastaZKFlavor are written from
 * the F_p rows (see construct_lookup_table_polynomials in ultra_pasta_honk.cpp).
 */
class UltraPastaCircuitBuilder : public UltraCircuitBuilder_<UltraPastaExecutionTraceBlocks> {
  public:
    using Base = UltraCircuitBuilder_<UltraPastaExecutionTraceBlocks>;
    using FF = Base::FF;
    using TableRow = std::array<FF, 3>;
    using Base::Base;

    /**
     * @brief Registers a lookup table with the given rows.
     */
    plookup::BasicTable& create_pasta_table(std::vector<TableRow> rows)
    {
        plookup::BasicTable table{};
        table.id = plookup::BasicTableId::HONK_DUMMY_BASIC1;
        table.use_twin_keys = false;
        table.get_values_from_key = nullptr;
        for (size_t i = 0; i < rows.size(); ++i) {
            table.column_1.emplace_back(i);
            table.column_2.emplace_back(0);
            table.column_3.emplace_back(0);
        }
        plookup::BasicTable* registered = this->register_basic_lookup_table(std::move(table));
        pasta_tables_.emplace(registered->table_index, std::move(rows));
        return *registered;
    }

    /**
     * @brief Constrains (key, value_1, value_2) to be the row `row_index` of `table`.
     */
    void create_pasta_lookup(
        plookup::BasicTable& table, size_t row_index, uint32_t key, uint32_t value_1, uint32_t value_2)
    {
        const auto& row = pasta_table_rows(table.table_index).at(row_index);
        if ((this->get_variable(key) != row[0] || this->get_variable(value_1) != row[1] ||
             this->get_variable(value_2) != row[2]) &&
            !this->failed()) {
            this->failure("pasta lookup: values are not the table row");
        }
        plookup::BasicTable::LookupEntry entry;
        entry.key = { uint256_t(row_index), 0 };
        this->create_lookup_gate(key, value_1, value_2, table, entry);
    }

    const std::vector<TableRow>& pasta_table_rows(size_t table_index) const { return pasta_tables_.at(table_index); }

  private:
    std::map<size_t, std::vector<TableRow>> pasta_tables_;
};

} // namespace bb::zcash
