// SPDX-License-Identifier: Apache-2.0
pragma solidity >=0.8.21;

import "../UltraFflonkTypes.sol";

/// @dev Port of relations/ultra_arithmetic_relation.hpp — the width-4 Ultra arithmetic gate, whose
/// two subrelations encode four different identities selected by q_arith in {0, 1, 2, 3}.
library ArithmeticRelation {
    uint256 internal constant R = UltraFflonkTypes.R;

    /// -1/2 in the BN254 scalar field: the C++ `FF(-2).invert()`, i.e. (R - 1) / 2.
    /// It rescales the multiplicative term so that (q_arith - 3) * q_m * (-1/2) equals q_m at
    /// q_arith == 1, twice q_m at q_arith == 2, and vanishes at q_arith == 3.
    uint256 internal constant NEG_HALF = 0x183227397098d014dc2822db40c0ac2e9419f4243cdcb848a1f0fac9f8000000;

    /// @dev Writes each subrelation's value at xi into `out`, in SUBRELATION_PARTIAL_LENGTHS order.
    function accumulate(uint256[41] memory e, UltraFflonkTypes.Params memory p, uint256[2] memory out) internal pure {
        p; // this relation reads no relation parameters

        uint256 qArith = e[UltraFflonkTypes.E_Q_ARITH];
        // (q_arith - 1) gates the extra w_4_shift term in the first subrelation and, together with
        // q_arith itself, disables the second subrelation for q_arith in {0, 1}.
        uint256 qArithSub1 = addmod(qArith, R - 1, R);

        // Subrelation 0: the gate identity itself. The whole expression is scaled by q_arith, so it
        // is vacuous on rows where the arithmetic gate is switched off.
        {
            uint256 acc;
            // Multiplicative term q_m * w_l * w_r, weighted by (-1/2) * (q_arith - 3).
            {
                uint256 mulTerm = mulmod(mulmod(e[UltraFflonkTypes.E_W_R], e[UltraFflonkTypes.E_W_L], R), NEG_HALF, R);
                uint256 qmWeight = mulmod(addmod(qArith, R - 3, R), e[UltraFflonkTypes.E_Q_M], R);
                acc = mulmod(mulTerm, qmWeight, R);
            }
            // Linear terms sum_{i in {l, r, o, 4}} q_i * w_i plus the constant selector.
            {
                uint256 lin = mulmod(e[UltraFflonkTypes.E_Q_L], e[UltraFflonkTypes.E_W_L], R);
                lin = addmod(lin, mulmod(e[UltraFflonkTypes.E_Q_R], e[UltraFflonkTypes.E_W_R], R), R);
                lin = addmod(lin, mulmod(e[UltraFflonkTypes.E_Q_O], e[UltraFflonkTypes.E_W_O], R), R);
                lin = addmod(lin, mulmod(e[UltraFflonkTypes.E_Q_4], e[UltraFflonkTypes.E_W_4], R), R);
                lin = addmod(lin, e[UltraFflonkTypes.E_Q_C], R);
                // Carry from the next row, enabled only for q_arith >= 2.
                lin = addmod(lin, mulmod(qArithSub1, e[UltraFflonkTypes.E_W_4_SHIFT], R), R);
                acc = addmod(acc, lin, R);
            }
            out[0] = mulmod(acc, qArith, R);
        }

        // Subrelation 1: active only at q_arith == 3, where q_m is repurposed as an additive term
        // and w_l_shift is constrained to w_l + w_4 + q_m.
        {
            uint256 sum = addmod(e[UltraFflonkTypes.E_W_L], e[UltraFflonkTypes.E_W_4], R);
            sum = addmod(sum, R - e[UltraFflonkTypes.E_W_L_SHIFT], R);
            sum = addmod(sum, e[UltraFflonkTypes.E_Q_M], R);

            uint256 gate = mulmod(sum, addmod(qArith, R - 2, R), R);
            out[1] = mulmod(gate, mulmod(qArithSub1, qArith, R), R);
        }
    }
}
