#include "barretenberg/zcash/ultra_pasta/ultra_pasta_builder.hpp"
#include "barretenberg/stdlib_circuit_builders/circuit_builder_base_impl.hpp"
#include "barretenberg/stdlib_circuit_builders/rom_ram_logic_impl.hpp"
#include "barretenberg/stdlib_circuit_builders/ultra_circuit_builder_impl.hpp"

namespace bb {
template class CircuitBuilderBase<pasta::fp>;
template class RomRamLogic_<zcash::UltraPastaExecutionTraceBlocks>;
template class UltraCircuitBuilder_<zcash::UltraPastaExecutionTraceBlocks>;
} // namespace bb
