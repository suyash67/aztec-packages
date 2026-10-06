// === AUDIT STATUS ===
// internal:    { status: Complete, auditors: [Luke, Raju], commit: }
// external_1:  { status: not started, auditors: [], commit: }
// external_2:  { status: not started, auditors: [], commit: }
// =====================

#include "barretenberg/stdlib_circuit_builders/ultra_circuit_builder_impl.hpp"

namespace bb {

template class UltraCircuitBuilder_<UltraExecutionTraceBlocks>;
template class UltraCircuitBuilder_<MegaExecutionTraceBlocks>;

} // namespace bb
