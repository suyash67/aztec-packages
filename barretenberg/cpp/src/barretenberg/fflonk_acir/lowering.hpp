#pragma once

#include "barretenberg/fflonk/circuit_builder.hpp"
#include "barretenberg/stdlib_circuit_builders/ultra_circuit_builder.hpp"

namespace bb::fflonk_acir {

/**
 * @brief Rewrite an Ultra circuit into the three-wire fflonk arithmetization.
 *
 * @details The three-wire system in `fflonk/` has no ACIR frontend of its own, and writing a second
 * one would mean a second place for the frontend to be wrong. Instead the ACIR frontend that already
 * exists produces an `UltraCircuitBuilder`, and this pass rewrites its trace.
 *
 * Each Ultra arithmetic row
 *
 * ```
 * q_m w_l w_r + q_l w_l + q_r w_r + q_o w_o + q_4 w_4 + q_c [+ w_4(next row)] = 0
 * ```
 *
 * spans five wire slots, which three wires cannot hold, so it becomes two or three three-wire rows
 * threaded by fresh intermediate variables. Copy constraints carry over by construction: two Ultra
 * slots holding the same variable map to the same three-wire variable, because the map is keyed on
 * `real_variable_index` - the representative after Ultra's own equality merging.
 *
 * **This covers arithmetic gates only.** Lookups, range constraints, elliptic additions, RAM/ROM and
 * Poseidon2 have no three-wire encoding that is anything but astronomically expensive - a lookup
 * "expanded into explicit constraints" costs the size of the table - so a circuit using them is
 * rejected here rather than silently lowered into something enormous. `--scheme ultra_fflonk` is the
 * backend for those circuits; this one is for arithmetic-only programs that want the smallest
 * possible on-chain verifier and the smallest possible trusted code base.
 *
 * @throws if the circuit uses a gate this pass does not cover.
 */
fflonk_plonk::CircuitBuilder lower(UltraCircuitBuilder& circuit);

} // namespace bb::fflonk_acir
