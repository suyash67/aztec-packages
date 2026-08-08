// SPDX-License-Identifier: Apache-2.0
pragma solidity >=0.8.21;

import "../UltraFflonkTypes.sol";

/// @dev Port of relations/permutation_relation.hpp — the copy-constraint grand product over the four wires.
library PermutationRelation {
    uint256 internal constant R = UltraFflonkTypes.R;

    /// @dev Writes each subrelation's value at xi into `out`, in SUBRELATION_PARTIAL_LENGTHS order.
    function accumulate(uint256[41] memory e, UltraFflonkTypes.Params memory p, uint256[3] memory out) internal pure {
        // (0) Grand product construction: the running product advances by the ratio of the id-labelled
        // wire tuple to the sigma-permuted one. `lagrange_first` seeds the product at the first row
        // (where z_perm is constrained to 0 by subrelation (2), making the factor 1), and the
        // public-input delta closes it at the last row.
        {
            uint256 left = mulmod(
                addmod(e[UltraFflonkTypes.E_Z_PERM], e[UltraFflonkTypes.E_LAGRANGE_FIRST], R), numerator(e, p), R
            );

            uint256 publicInputTerm = addmod(
                mulmod(e[UltraFflonkTypes.E_LAGRANGE_LAST], p.publicInputDelta, R),
                e[UltraFflonkTypes.E_Z_PERM_SHIFT],
                R
            );

            out[0] = addmod(left, R - mulmod(publicInputTerm, denominator(e, p), R), R);
        }

        // (1) Left-shiftability: z_perm must vanish immediately after the last active row.
        out[1] = mulmod(e[UltraFflonkTypes.E_LAGRANGE_LAST], e[UltraFflonkTypes.E_Z_PERM_SHIFT], R);

        // (2) Initialization: z_perm is 0 at the first row, so that (z_perm + L_first) starts the
        // product at 1. Without it a prover could pick an arbitrary starting value.
        out[2] = mulmod(e[UltraFflonkTypes.E_LAGRANGE_FIRST], e[UltraFflonkTypes.E_Z_PERM], R);
    }

    /// @dev Grand product numerator: the wires tagged by their own row identity, prod_i (w_i + id_i*beta + gamma).
    function numerator(uint256[41] memory e, UltraFflonkTypes.Params memory p) private pure returns (uint256) {
        uint256 acc;
        {
            uint256 t1 =
                addmod(mulmod(e[UltraFflonkTypes.E_ID_1], p.beta, R), addmod(e[UltraFflonkTypes.E_W_L], p.gamma, R), R);
            uint256 t2 =
                addmod(mulmod(e[UltraFflonkTypes.E_ID_2], p.beta, R), addmod(e[UltraFflonkTypes.E_W_R], p.gamma, R), R);
            acc = mulmod(t1, t2, R);
        }
        {
            uint256 t3 =
                addmod(mulmod(e[UltraFflonkTypes.E_ID_3], p.beta, R), addmod(e[UltraFflonkTypes.E_W_O], p.gamma, R), R);
            uint256 t4 =
                addmod(mulmod(e[UltraFflonkTypes.E_ID_4], p.beta, R), addmod(e[UltraFflonkTypes.E_W_4], p.gamma, R), R);
            acc = mulmod(acc, mulmod(t3, t4, R), R);
        }
        return acc;
    }

    /// @dev Grand product denominator: the same wires tagged by their permutation image,
    /// prod_i (w_i + sigma_i*beta + gamma).
    function denominator(uint256[41] memory e, UltraFflonkTypes.Params memory p) private pure returns (uint256) {
        uint256 acc;
        {
            uint256 t5 = addmod(
                mulmod(e[UltraFflonkTypes.E_SIGMA_1], p.beta, R), addmod(e[UltraFflonkTypes.E_W_L], p.gamma, R), R
            );
            uint256 t6 = addmod(
                mulmod(e[UltraFflonkTypes.E_SIGMA_2], p.beta, R), addmod(e[UltraFflonkTypes.E_W_R], p.gamma, R), R
            );
            acc = mulmod(t5, t6, R);
        }
        {
            uint256 t7 = addmod(
                mulmod(e[UltraFflonkTypes.E_SIGMA_3], p.beta, R), addmod(e[UltraFflonkTypes.E_W_O], p.gamma, R), R
            );
            uint256 t8 = addmod(
                mulmod(e[UltraFflonkTypes.E_SIGMA_4], p.beta, R), addmod(e[UltraFflonkTypes.E_W_4], p.gamma, R), R
            );
            acc = mulmod(acc, mulmod(t7, t8, R), R);
        }
        return acc;
    }
}
