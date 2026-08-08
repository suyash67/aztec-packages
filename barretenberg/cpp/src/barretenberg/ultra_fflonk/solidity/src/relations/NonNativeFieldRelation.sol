// SPDX-License-Identifier: Apache-2.0
pragma solidity >=0.8.21;

import "../UltraFflonkTypes.sol";

/// @dev Port of relations/non_native_field_relation.hpp — bigfield (non-native) arithmetic over
/// 68-bit limbs: the three limb-product gates and the two 14-bit sublimb accumulation gates,
/// all gated by q_nnf.
library NonNativeFieldRelation {
    uint256 internal constant R = UltraFflonkTypes.R;

    /// 2^68, the bigfield limb radix (`FF LIMB_SIZE(uint256_t(1) << 68)`).
    uint256 private constant LIMB_SIZE = 0x100000000000000000;

    /// 2^14, the sublimb radix used by the limb accumulation gates
    /// (`FF SUBLIMB_SHIFT(uint256_t(1) << 14)`).
    uint256 private constant SUBLIMB_SHIFT = 0x4000;

    /// @dev Writes each subrelation's value at xi into `out`, in SUBRELATION_PARTIAL_LENGTHS order.
    function accumulate(uint256[41] memory e, UltraFflonkTypes.Params memory p, uint256[1] memory out) internal pure {
        p; // the non-native field relation reads no relation parameters

        // Subrelation 0: the single combined non-native field identity. The two families of gates
        // (limb products, selected by q_2; limb accumulation, selected by q_3) are summed and the
        // whole thing is switched on by q_nnf.
        uint256 identity = addmod(productIdentity(e), limbAccumulatorIdentity(e), R);
        out[0] = mulmod(identity, e[UltraFflonkTypes.E_Q_NNF], R);
    }

    /// @dev The three bigfield product gates, each selected by its own secondary selector and the
    /// family as a whole gated by q_2. They check that the cross-terms of a 4-limb x 4-limb
    /// multiplication accumulate correctly with their 2^68 weights.
    function productIdentity(uint256[41] memory e) private pure returns (uint256) {
        // The cross-term shared by all three gates: w_1 * w_2' + w_1' * w_2.
        uint256 limbSubproduct = addmod(
            mulmod(e[UltraFflonkTypes.E_W_L], e[UltraFflonkTypes.E_W_R_SHIFT], R),
            mulmod(e[UltraFflonkTypes.E_W_L_SHIFT], e[UltraFflonkTypes.E_W_R], R),
            R
        );

        // Gate 2 consumes the un-scaled subproduct, so it must be evaluated before the rescale below.
        uint256 acc = mulmod(productGate2(e, limbSubproduct), e[UltraFflonkTypes.E_Q_4], R);

        // Gates 1 and 3 share the scaled subproduct (w_1 * w_2' + w_1' * w_2) * 2^68 + w_1' * w_2'.
        limbSubproduct = mulmod(limbSubproduct, LIMB_SIZE, R);
        limbSubproduct =
            addmod(limbSubproduct, mulmod(e[UltraFflonkTypes.E_W_L_SHIFT], e[UltraFflonkTypes.E_W_R_SHIFT], R), R);

        // Gate 1 (selector q_3): scaled subproduct minus the low accumulators w_3 + w_4.
        {
            uint256 t = addmod(limbSubproduct, R - addmod(e[UltraFflonkTypes.E_W_O], e[UltraFflonkTypes.E_W_4], R), R);
            acc = addmod(acc, mulmod(t, e[UltraFflonkTypes.E_Q_O], R), R);
        }

        // Gate 3 (selector q_m): scaled subproduct + w_4 - w_3' - w_4'.
        {
            uint256 t = addmod(limbSubproduct, e[UltraFflonkTypes.E_W_4], R);
            t = addmod(t, R - addmod(e[UltraFflonkTypes.E_W_O_SHIFT], e[UltraFflonkTypes.E_W_4_SHIFT], R), R);
            acc = addmod(acc, mulmod(t, e[UltraFflonkTypes.E_Q_M], R), R);
        }

        return mulmod(acc, e[UltraFflonkTypes.E_Q_R], R);
    }

    /// @dev Bigfield product gate 2 (selector q_4), before its selector is applied:
    ///   limb_subproduct + (w_1 * w_4 + w_2 * w_3 - w_3') * 2^68 - w_4'
    function productGate2(uint256[41] memory e, uint256 limbSubproduct) private pure returns (uint256) {
        uint256 t = addmod(
            mulmod(e[UltraFflonkTypes.E_W_L], e[UltraFflonkTypes.E_W_4], R),
            mulmod(e[UltraFflonkTypes.E_W_R], e[UltraFflonkTypes.E_W_O], R),
            R
        );
        t = addmod(t, R - e[UltraFflonkTypes.E_W_O_SHIFT], R);
        t = mulmod(t, LIMB_SIZE, R);
        t = addmod(t, R - e[UltraFflonkTypes.E_W_4_SHIFT], R);
        return addmod(t, limbSubproduct, R);
    }

    /// @dev The two limb accumulation gates, gated as a family by q_3. Each reconstructs a 68-bit
    /// limb from five 14-bit sublimbs laid out across the current and next row, in Horner form.
    function limbAccumulatorIdentity(uint256[41] memory e) private pure returns (uint256) {
        uint256 acc;

        // Gate 1 (selector q_4): ((((w_2' * 2^14 + w_1') * 2^14 + w_3) * 2^14 + w_2) * 2^14 + w_1) - w_4
        {
            uint256 t = mulmod(e[UltraFflonkTypes.E_W_R_SHIFT], SUBLIMB_SHIFT, R);
            t = addmod(t, e[UltraFflonkTypes.E_W_L_SHIFT], R);
            t = mulmod(t, SUBLIMB_SHIFT, R);
            t = addmod(t, e[UltraFflonkTypes.E_W_O], R);
            t = mulmod(t, SUBLIMB_SHIFT, R);
            t = addmod(t, e[UltraFflonkTypes.E_W_R], R);
            t = mulmod(t, SUBLIMB_SHIFT, R);
            t = addmod(t, e[UltraFflonkTypes.E_W_L], R);
            t = addmod(t, R - e[UltraFflonkTypes.E_W_4], R);
            acc = mulmod(t, e[UltraFflonkTypes.E_Q_4], R);
        }

        // Gate 2 (selector q_m): ((((w_3' * 2^14 + w_2') * 2^14 + w_1') * 2^14 + w_4) * 2^14 + w_3) - w_4'
        {
            uint256 t = mulmod(e[UltraFflonkTypes.E_W_O_SHIFT], SUBLIMB_SHIFT, R);
            t = addmod(t, e[UltraFflonkTypes.E_W_R_SHIFT], R);
            t = mulmod(t, SUBLIMB_SHIFT, R);
            t = addmod(t, e[UltraFflonkTypes.E_W_L_SHIFT], R);
            t = mulmod(t, SUBLIMB_SHIFT, R);
            t = addmod(t, e[UltraFflonkTypes.E_W_4], R);
            t = mulmod(t, SUBLIMB_SHIFT, R);
            t = addmod(t, e[UltraFflonkTypes.E_W_O], R);
            t = addmod(t, R - e[UltraFflonkTypes.E_W_4_SHIFT], R);
            acc = addmod(acc, mulmod(t, e[UltraFflonkTypes.E_Q_M], R), R);
        }

        return mulmod(acc, e[UltraFflonkTypes.E_Q_O], R);
    }
}
