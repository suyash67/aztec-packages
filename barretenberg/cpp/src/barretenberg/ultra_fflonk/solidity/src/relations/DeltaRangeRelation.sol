// SPDX-License-Identifier: Apache-2.0
pragma solidity >=0.8.21;

import "../UltraFflonkTypes.sol";

/// @dev Port of relations/delta_range_constraint_relation.hpp — each of the four consecutive wire
/// differences across the row (and into the next row) lies in {0, 1, 2, 3}, which is what makes the
/// sorted-witness range check work.
library DeltaRangeRelation {
    uint256 internal constant R = UltraFflonkTypes.R;

    /// @dev Writes each subrelation's value at xi into `out`, in SUBRELATION_PARTIAL_LENGTHS order.
    function accumulate(uint256[41] memory e, UltraFflonkTypes.Params memory p, uint256[4] memory out) internal pure {
        uint256 qDeltaRange = e[UltraFflonkTypes.E_Q_DELTA_RANGE];

        // Subrelation 0: D_0 = w_2 - w_1 is in {0,1,2,3}.
        {
            uint256 delta = addmod(e[UltraFflonkTypes.E_W_R], R - e[UltraFflonkTypes.E_W_L], R);
            out[0] = gatedRangeCheck(delta, qDeltaRange);
        }

        // Subrelation 1: D_1 = w_3 - w_2 is in {0,1,2,3}.
        {
            uint256 delta = addmod(e[UltraFflonkTypes.E_W_O], R - e[UltraFflonkTypes.E_W_R], R);
            out[1] = gatedRangeCheck(delta, qDeltaRange);
        }

        // Subrelation 2: D_2 = w_4 - w_3 is in {0,1,2,3}.
        {
            uint256 delta = addmod(e[UltraFflonkTypes.E_W_4], R - e[UltraFflonkTypes.E_W_O], R);
            out[2] = gatedRangeCheck(delta, qDeltaRange);
        }

        // Subrelation 3: D_3 = w_1_shift - w_4 links this row's last wire to the next row's first.
        {
            uint256 delta = addmod(e[UltraFflonkTypes.E_W_L_SHIFT], R - e[UltraFflonkTypes.E_W_4], R);
            out[3] = gatedRangeCheck(delta, qDeltaRange);
        }

        p; // the relation reads no relation parameters
    }

    /// @dev q * D(D-1)(D-2)(D-3), evaluated via the identity T*(T+2) == D(D-1)(D-2)(D-3) for T = (D-3)*D.
    function gatedRangeCheck(uint256 delta, uint256 qDeltaRange) private pure returns (uint256) {
        uint256 t = mulmod(addmod(delta, R - 3, R), delta, R);
        return mulmod(mulmod(t, addmod(t, 2, R), R), qDeltaRange, R);
    }
}
