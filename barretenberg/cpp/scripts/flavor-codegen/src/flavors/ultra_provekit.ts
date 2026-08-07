import { flavor } from "../flavor.js";
import * as R from "../relations/index.js";

// UltraFlavor restricted to the gate kinds a ProveKit-style circuit actually uses. Elliptic,
// non-native-field and both Poseidon2 relations are dropped: circuits compiled from Noir that do
// no in-circuit curve arithmetic, bignum arithmetic or Poseidon2 hashing leave those trace blocks
// empty, so their selectors are identically zero and their subrelations vanish on every row.
//
// Dropping them shrinks the layout by four precomputed selectors and, because the Poseidon2
// subrelations are the only degree-6 ones in Ultra, drops MAX_PARTIAL_RELATION_LENGTH from 7 to 6
// — one fewer evaluation per sumcheck round polynomial and a shorter edge extension in the hot
// loop. The relation order of the columns it keeps matches ultra.ts, so the trace layout is
// UltraCircuitBuilder's unchanged (`emitsTrace: false` — this flavor reuses the Ultra trace).
export const UltraProveKit = flavor({
    name: "UltraProveKitFlavor",
    family: "ultra_provekit",
    generatedClassName: "UltraProveKitFlavor_Generated",
    relations: [
        // Same ordering rationale as ultra.ts: UltraPermutation first (no gate block) keeps the
        // to-be-shifted witnesses contiguous, LogDerivLookup before Arithmetic keeps the lookup
        // block early in the trace.
        R.UltraPermutationRelation,
        R.LogDerivLookupRelation,
        R.ArithmeticRelation,
        R.DeltaRangeConstraintRelation,
        R.MemoryRelation,
    ],
    composites: {
        selectors: ["non_gate_selectors", "gate_selectors"],
    },
    traceExtraBlocks: ["pub_inputs"],
    emitsTrace: false,
});
