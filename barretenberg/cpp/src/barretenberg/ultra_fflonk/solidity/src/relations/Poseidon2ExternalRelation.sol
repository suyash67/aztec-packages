// SPDX-License-Identifier: Apache-2.0
pragma solidity >=0.8.21;

import "../UltraFflonkTypes.sol";

/// @dev Port of relations/poseidon2_external_relation.hpp — one Poseidon2 external round: add the
/// round constants, apply the x^5 S-box to all four state words, multiply by the external MDS
/// matrix M_E, and constrain the result to equal the next row's witness state.
library Poseidon2ExternalRelation {
    uint256 internal constant R = UltraFflonkTypes.R;

    /// @dev The Poseidon2 S-box, x -> x^5.
    function sbox(uint256 x) private pure returns (uint256) {
        uint256 x2 = mulmod(x, x, R);
        uint256 x4 = mulmod(x2, x2, R);
        return mulmod(x4, x, R);
    }

    /// @dev Writes each subrelation's value at xi into `out`, in SUBRELATION_PARTIAL_LENGTHS order.
    function accumulate(uint256[41] memory e, UltraFflonkTypes.Params memory p, uint256[4] memory out) internal pure {
        p; // the external round reads no relation parameters

        // S-box inputs are the current state plus this round's constants, which are carried in the
        // arithmetic selectors q_l, q_r, q_o, q_4 on Poseidon2 external rows.
        uint256[4] memory u;
        u[0] = sbox(addmod(e[UltraFflonkTypes.E_W_L], e[UltraFflonkTypes.E_Q_L], R));
        u[1] = sbox(addmod(e[UltraFflonkTypes.E_W_R], e[UltraFflonkTypes.E_Q_R], R));
        u[2] = sbox(addmod(e[UltraFflonkTypes.E_W_O], e[UltraFflonkTypes.E_Q_O], R));
        u[3] = sbox(addmod(e[UltraFflonkTypes.E_W_4], e[UltraFflonkTypes.E_Q_4], R));

        // v = M_E * u for the external matrix
        //   [5 7 1 3; 4 6 1 1; 1 3 5 7; 1 1 4 6],
        // built from shared partial sums exactly as the C++ does.
        uint256[4] memory v;
        {
            uint256 t0 = addmod(u[0], u[1], R); // u1 + u2
            uint256 t1 = addmod(u[2], u[3], R); // u3 + u4
            uint256 t2 = addmod(addmod(u[1], u[1], R), t1, R); // 2u2 + u3 + u4
            uint256 t3 = addmod(addmod(u[3], u[3], R), t0, R); // u1 + u2 + 2u4

            uint256 v2 = addmod(t0, t0, R);
            v2 = addmod(v2, v2, R);
            v2 = addmod(v2, t2, R); // 4u1 + 6u2 + u3 + u4

            uint256 v4 = addmod(t1, t1, R);
            v4 = addmod(v4, v4, R);
            v4 = addmod(v4, t3, R); // u1 + u2 + 4u3 + 6u4

            v[0] = addmod(t3, v2, R); // 5u1 + 7u2 + u3 + 3u4
            v[1] = v2;
            v[2] = addmod(t2, v4, R); // u1 + 3u2 + 5u3 + 7u4
            v[3] = v4;
        }

        // Each subrelation gates "round output equals the next row's state word" on the
        // Poseidon2 external selector.
        uint256 q = e[UltraFflonkTypes.E_Q_POSEIDON2_EXTERNAL];
        out[0] = mulmod(q, addmod(v[0], R - e[UltraFflonkTypes.E_W_L_SHIFT], R), R);
        out[1] = mulmod(q, addmod(v[1], R - e[UltraFflonkTypes.E_W_R_SHIFT], R), R);
        out[2] = mulmod(q, addmod(v[2], R - e[UltraFflonkTypes.E_W_O_SHIFT], R), R);
        out[3] = mulmod(q, addmod(v[3], R - e[UltraFflonkTypes.E_W_4_SHIFT], R), R);
    }
}
