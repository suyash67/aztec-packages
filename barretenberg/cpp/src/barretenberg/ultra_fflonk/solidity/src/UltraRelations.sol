// SPDX-License-Identifier: Apache-2.0
pragma solidity >=0.8.21;

import "./UltraFflonkTypes.sol";
import "./relations/PermutationRelation.sol";
import "./relations/LogDerivLookupRelation.sol";
import "./relations/ArithmeticRelation.sol";
import "./relations/DeltaRangeRelation.sol";
import "./relations/EllipticRelation.sol";
import "./relations/MemoryRelation.sol";
import "./relations/NonNativeFieldRelation.sol";
import "./relations/Poseidon2ExternalRelation.sol";
import "./relations/Poseidon2InternalRelation.sol";

/**
 * @title The Ultra relation set, evaluated at one point.
 *
 * @dev The order below is `UltraFlavor::Relations_`, which is generated from
 * `scripts/flavor-codegen/src/flavors/ultra.ts`. It decides which power of `alpha` batches which
 * subrelation, so it is protocol: reordering it here silently invalidates every proof. The mapping is
 * pinned by `RelationVectors.t.sol`, which checks all thirty-one values against barretenberg's own.
 *
 * Two of the thirty-one are "linearly dependent" — Sumcheck only requires them to sum to zero across
 * the trace, not to vanish at each row. They are still evaluated here; it is the *verifier* that
 * routes them into the running-sum identities instead of the per-row batch.
 */
library UltraRelations {
    uint256 internal constant NUM_SUBRELATIONS = 31;

    /// Global indices of the subrelations Sumcheck enforces as a trace sum, in relation order.
    uint256 internal constant SUMMED_LOOKUP = 4; // LogDerivLookup's lookup identity
    uint256 internal constant SUMMED_ROM_LOGUP = 21; // Memory's ROM-LogUp sum

    /// @dev Writes every subrelation's value at the evaluation point into `out`.
    function evaluate(uint256[41] memory e, UltraFflonkTypes.Params memory p, uint256[NUM_SUBRELATIONS] memory out)
        internal
        pure
    {
        {
            uint256[3] memory o;
            PermutationRelation.accumulate(e, p, o);
            out[0] = o[0];
            out[1] = o[1];
            out[2] = o[2];
        }
        {
            uint256[3] memory o;
            LogDerivLookupRelation.accumulate(e, p, o);
            out[3] = o[0];
            out[4] = o[1];
            out[5] = o[2];
        }
        {
            uint256[2] memory o;
            ArithmeticRelation.accumulate(e, p, o);
            out[6] = o[0];
            out[7] = o[1];
        }
        {
            uint256[4] memory o;
            DeltaRangeRelation.accumulate(e, p, o);
            out[8] = o[0];
            out[9] = o[1];
            out[10] = o[2];
            out[11] = o[3];
        }
        {
            uint256[2] memory o;
            EllipticRelation.accumulate(e, p, o);
            out[12] = o[0];
            out[13] = o[1];
        }
        {
            uint256[8] memory o;
            MemoryRelation.accumulate(e, p, o);
            for (uint256 i = 0; i < 8; ++i) {
                out[14 + i] = o[i];
            }
        }
        {
            uint256[1] memory o;
            NonNativeFieldRelation.accumulate(e, p, o);
            out[22] = o[0];
        }
        {
            uint256[4] memory o;
            Poseidon2ExternalRelation.accumulate(e, p, o);
            for (uint256 i = 0; i < 4; ++i) {
                out[23 + i] = o[i];
            }
        }
        {
            uint256[4] memory o;
            Poseidon2InternalRelation.accumulate(e, p, o);
            for (uint256 i = 0; i < 4; ++i) {
                out[27 + i] = o[i];
            }
        }
    }
}
