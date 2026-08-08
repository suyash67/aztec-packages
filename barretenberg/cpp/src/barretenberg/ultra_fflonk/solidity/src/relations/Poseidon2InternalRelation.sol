// SPDX-License-Identifier: Apache-2.0
pragma solidity >=0.8.21;

import "../UltraFflonkTypes.sol";

/// @dev Port of relations/poseidon2_internal_relation.hpp — one Poseidon2 internal (partial) round:
/// the S-box is applied to the first state element only, the internal MDS matrix M_I is applied to
/// the resulting state, and the four outputs must equal the next row's wires.
library Poseidon2InternalRelation {
    uint256 internal constant R = UltraFflonkTypes.R;

    /// Diagonal entries of the internal matrix M_I, from
    /// crypto/poseidon2/poseidon2_params.hpp::Poseidon2Bn254ScalarFieldParams. That header stores
    /// `internal_matrix_diagonal_minus_one[i]` = D_i - 1; the values below are copied verbatim and
    /// D1 is the +1 of the first one (equal, as expected, to `internal_matrix[0][0]`).
    uint256 internal constant D1 = 0x10dc6e9c006ea38b04b1e03b4bd9490c0d03f98929ca1d7fb56821fd19d3b6e8;
    uint256 internal constant D2_MINUS_1 = 0x0c28145b6a44df3e0149b3d0a30b3bb599df9756d4dd9b84a86b38cfb45a740b;
    uint256 internal constant D3_MINUS_1 = 0x00544b8338791518b2c7645a50392798b21f75bb60e3596170067d00141cac15;
    uint256 internal constant D4_MINUS_1 = 0x222c01175718386f2e2e82eb122789e352e105a3b8fa852613bc534433ee428b;

    /// @dev Writes each subrelation's value at xi into `out`, in SUBRELATION_PARTIAL_LENGTHS order.
    function accumulate(uint256[41] memory e, UltraFflonkTypes.Params memory p, uint256[4] memory out) internal pure {
        // The whole relation is gated on the internal-round selector.
        uint256 qPoseidon = e[UltraFflonkTypes.E_Q_POSEIDON2_INTERNAL];

        // u_1 = (w_1 + c_0)^5, the S-box output for the only element the internal round touches.
        // The round constant c_0 for this round is carried by the q_l selector. Every subrelation
        // uses it already multiplied by the selector.
        uint256 scaledU1;
        {
            uint256 s1 = addmod(e[UltraFflonkTypes.E_W_L], e[UltraFflonkTypes.E_Q_L], R);
            uint256 s1Pow4 = mulmod(s1, s1, R);
            s1Pow4 = mulmod(s1Pow4, s1Pow4, R);
            scaledU1 = mulmod(mulmod(s1Pow4, s1, R), qPoseidon, R);
        }

        // u_2 + u_3 + u_4: the off-diagonal part of M_I is all-ones, so this sum appears in every row.
        // The untouched state elements pass through the S-box layer unchanged, so u_k = w_k for k > 1.
        uint256 partialSum;
        {
            uint256 s = addmod(e[UltraFflonkTypes.E_W_R], e[UltraFflonkTypes.E_W_O], R);
            partialSum = addmod(s, e[UltraFflonkTypes.E_W_4], R);
        }

        // Row 1: D_1*u_1 + u_2 + u_3 + u_4 == w_1_shift.
        {
            uint256 rest = addmod(partialSum, R - e[UltraFflonkTypes.E_W_L_SHIFT], R);
            out[0] = addmod(mulmod(scaledU1, D1, R), mulmod(rest, qPoseidon, R), R);
        }

        // Row 2: u_1 + D_2*u_2 + u_3 + u_4 == w_2_shift, written as (D_2 - 1)*u_2 + partialSum so the
        // shared sum can be reused.
        {
            uint256 rest = addmod(mulmod(e[UltraFflonkTypes.E_W_R], D2_MINUS_1, R), partialSum, R);
            rest = addmod(rest, R - e[UltraFflonkTypes.E_W_R_SHIFT], R);
            out[1] = addmod(mulmod(rest, qPoseidon, R), scaledU1, R);
        }

        // Row 3: u_1 + u_2 + D_3*u_3 + u_4 == w_3_shift.
        {
            uint256 rest = addmod(mulmod(e[UltraFflonkTypes.E_W_O], D3_MINUS_1, R), partialSum, R);
            rest = addmod(rest, R - e[UltraFflonkTypes.E_W_O_SHIFT], R);
            out[2] = addmod(mulmod(rest, qPoseidon, R), scaledU1, R);
        }

        // Row 4: u_1 + u_2 + u_3 + D_4*u_4 == w_4_shift.
        {
            uint256 rest = addmod(mulmod(e[UltraFflonkTypes.E_W_4], D4_MINUS_1, R), partialSum, R);
            rest = addmod(rest, R - e[UltraFflonkTypes.E_W_4_SHIFT], R);
            out[3] = addmod(mulmod(rest, qPoseidon, R), scaledU1, R);
        }

        p; // the relation reads no relation parameters
    }
}
