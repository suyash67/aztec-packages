// === AUDIT STATUS ===
// internal:    { status: Complete, auditors: [Raju], commit: 05a381f8b31ae4648e480f1369e911b148216e8b}
// external_1:  { status: Complete, auditors: [Sherlock], commit: e6694849223 }
// external_2:  { status: not started, auditors: [], commit: }
// =====================

#include "barretenberg/stdlib_circuit_builders/rom_ram_logic_impl.hpp"

namespace bb {

template class RomRamLogic_<UltraExecutionTraceBlocks>;
template class RomRamLogic_<MegaExecutionTraceBlocks>;

} // namespace bb
