// SPDX-License-Identifier: Apache-2.0
pragma solidity >=0.8.21;

import "../UltraFflonkTypes.sol";

/// @dev Port of relations/memory_relation.hpp — RAM/ROM memory gates: record fingerprints, sorted-trace
/// consistency for RAM and ROM, RAM timestamp jumps, and the single-value ROM LogUp argument.
library MemoryRelation {
    uint256 internal constant R = UltraFflonkTypes.R;

    /// @dev Writes each subrelation's value at xi into `out`, in SUBRELATION_PARTIAL_LENGTHS order.
    function accumulate(uint256[41] memory e, UltraFflonkTypes.Params memory p, uint256[8] memory out) internal pure {
        uint256 qMemory = e[UltraFflonkTypes.E_Q_MEMORY];

        // Record fingerprint residual `q_c + eta*w_l + eta_two*w_r + eta_three*w_o - w_4`. On an honest
        // RAM/ROM access gate this is zero; on a sorted RAM entry the same expression is the negated
        // access type, so it doubles as the value the boolean access check is applied to.
        uint256 memoryRecordCheck = recordCheckResidual(e, p);

        // Sorted traces are ordered by index, so the index delta between adjacent rows must be 0 or 1.
        // `negIndexDelta = w_l - w_l_shift` is 0 or -1, hence `indexDeltaIsZero` is the indicator of a
        // repeated index and `negIndexDelta^2 + negIndexDelta` vanishes exactly on {0, -1}.
        uint256 indexDeltaIsZero;
        uint256 indexIncreasesByZeroOrOne;
        {
            uint256 negIndexDelta = addmod(e[UltraFflonkTypes.E_W_L], R - e[UltraFflonkTypes.E_W_L_SHIFT], R);
            indexDeltaIsZero = addmod(negIndexDelta, 1, R);
            indexIncreasesByZeroOrOne = addmod(mulmod(negIndexDelta, negIndexDelta, R), negIndexDelta, R);
        }

        // ROM sorted-trace subrelations, selected by the ROM_CONSISTENCY_CHECK bitpattern q_1 = q_2 = 1.
        {
            uint256 romGate = mulmod(mulmod(e[UltraFflonkTypes.E_Q_L], e[UltraFflonkTypes.E_Q_R], R), qMemory, R);
            // If two adjacent ROM entries share an index, their full records (and hence both stored
            // values) must agree.
            uint256 recordDelta = addmod(e[UltraFflonkTypes.E_W_4_SHIFT], R - e[UltraFflonkTypes.E_W_4], R);
            out[1] = mulmod(mulmod(indexDeltaIsZero, recordDelta, R), romGate, R);
            out[2] = mulmod(indexIncreasesByZeroOrOne, romGate, R);
        }

        // Negated access type of the *next* sorted row, read off its own record fingerprint. The final
        // sorted RAM entry is left unselected, so checking the next row's type uniformly covers it.
        uint256 negNextGateAccessType;
        {
            uint256 acc = mulmod(e[UltraFflonkTypes.E_W_O_SHIFT], p.etaThree, R);
            acc = addmod(acc, mulmod(e[UltraFflonkTypes.E_W_R_SHIFT], p.etaTwo, R), R);
            acc = addmod(acc, mulmod(e[UltraFflonkTypes.E_W_L_SHIFT], p.eta, R), R);
            negNextGateAccessType = addmod(acc, R - e[UltraFflonkTypes.E_W_4_SHIFT], R);
        }

        // RAM sorted-trace subrelations, selected by the RAM_CONSISTENCY_CHECK bitpattern q_3 = 1.
        uint256 q3ByMemory = mulmod(e[UltraFflonkTypes.E_Q_O], qMemory, R);
        {
            // (indices match) && (next access is a read) => values match.
            uint256 valueDelta = addmod(e[UltraFflonkTypes.E_W_O_SHIFT], R - e[UltraFflonkTypes.E_W_O], R);
            uint256 gated = mulmod(indexDeltaIsZero, valueDelta, R);
            gated = mulmod(gated, addmod(negNextGateAccessType, 1, R), R);
            out[3] = mulmod(gated, q3ByMemory, R);
        }
        out[4] = mulmod(indexIncreasesByZeroOrOne, q3ByMemory, R);
        {
            // The next gate's access type is a read or a write, and nothing else.
            uint256 nextAccessTypeIsBoolean =
                addmod(mulmod(negNextGateAccessType, negNextGateAccessType, R), negNextGateAccessType, R);
            out[5] = mulmod(nextAccessTypeIsBoolean, q3ByMemory, R);
        }

        // Subrelation 0 bundles the three q_memory-gated identities that share the record fingerprint,
        // plus the (already q_3-gated) RAM access-type boolean check.
        {
            // ROM read: the record witness is consistent with (index, value1, value2). q_1 = q_2 = 1.
            uint256 identity =
                mulmod(memoryRecordCheck, mulmod(e[UltraFflonkTypes.E_Q_L], e[UltraFflonkTypes.E_Q_R], R), R);
            {
                // RAM timestamp check (q_1 = q_4 = 1): w_o holds the timestamp delta when the index
                // repeats, and zero when the index advances.
                uint256 timestampDelta = addmod(e[UltraFflonkTypes.E_W_R_SHIFT], R - e[UltraFflonkTypes.E_W_R], R);
                uint256 timestampCheck =
                    addmod(mulmod(indexDeltaIsZero, timestampDelta, R), R - e[UltraFflonkTypes.E_W_O], R);
                uint256 timestampGate = mulmod(e[UltraFflonkTypes.E_Q_4], e[UltraFflonkTypes.E_Q_L], R);
                identity = addmod(identity, mulmod(timestampCheck, timestampGate, R), R);
            }
            // RAM/ROM access gate (q_1 = q_m = 1): the record fingerprint must match w_4.
            {
                uint256 accessGate = mulmod(e[UltraFflonkTypes.E_Q_M], e[UltraFflonkTypes.E_Q_L], R);
                identity = addmod(identity, mulmod(memoryRecordCheck, accessGate, R), R);
            }
            identity = mulmod(identity, qMemory, R);

            // On sorted RAM rows the record residual is the negated access type; force it boolean.
            uint256 accessCheck = addmod(mulmod(memoryRecordCheck, memoryRecordCheck, R), memoryRecordCheck, R);
            out[0] = addmod(identity, mulmod(accessCheck, q3ByMemory, R), R);
        }

        // Single-value ROM tables use a LogUp argument rather than the sorted-trace permutation.
        out[6] = romLogupInverseCheck(e, p);
        out[7] = romLogupSum(e, p);
    }

    /// @dev `q_c + eta*w_l + eta_two*w_r + eta_three*w_o - w_4`: the record fingerprint residual.
    function recordCheckResidual(uint256[41] memory e, UltraFflonkTypes.Params memory p)
        private
        pure
        returns (uint256)
    {
        uint256 acc = mulmod(e[UltraFflonkTypes.E_W_O], p.etaThree, R);
        acc = addmod(acc, mulmod(e[UltraFflonkTypes.E_W_R], p.etaTwo, R), R);
        acc = addmod(acc, mulmod(e[UltraFflonkTypes.E_W_L], p.eta, R), R);
        acc = addmod(acc, e[UltraFflonkTypes.E_Q_C], R);
        return addmod(acc, R - e[UltraFflonkTypes.E_W_4], R);
    }

    /// @dev LogUp fingerprint denominator `rom_logup_gamma + w_l + eta*w_r + eta_two*q_c`, batching
    /// (index, value, ROM array id) with powers of eta under an independent additive challenge.
    function romLogupDenominator(uint256[41] memory e, UltraFflonkTypes.Params memory p)
        private
        pure
        returns (uint256)
    {
        uint256 acc = addmod(e[UltraFflonkTypes.E_W_L], mulmod(e[UltraFflonkTypes.E_W_R], p.eta, R), R);
        acc = addmod(acc, mulmod(e[UltraFflonkTypes.E_Q_C], p.etaTwo, R), R);
        return addmod(acc, p.romLogupGamma, R);
    }

    /// @dev Subrelation 6: on every ROM-LogUp row (table entry q_2*(1-q_1) or read access q_4*(1-q_1)),
    /// w_4 must be the inverse of the fingerprint denominator.
    function romLogupInverseCheck(uint256[41] memory e, UltraFflonkTypes.Params memory p)
        private
        pure
        returns (uint256)
    {
        uint256 oneMinusQ1 = addmod(1, R - e[UltraFflonkTypes.E_Q_L], R);
        uint256 qLogupAny = addmod(e[UltraFflonkTypes.E_Q_R], e[UltraFflonkTypes.E_Q_4], R);
        qLogupAny = mulmod(qLogupAny, oneMinusQ1, R);

        uint256 inverseIsCorrect = addmod(mulmod(e[UltraFflonkTypes.E_W_4], romLogupDenominator(e, p), R), R - 1, R);
        return mulmod(mulmod(e[UltraFflonkTypes.E_Q_MEMORY], qLogupAny, R), inverseIsCorrect, R);
    }

    /// @dev Subrelation 7 (linearly dependent, summed across the trace): read rows contribute +1/denom,
    /// table rows contribute -m_i/denom, where the multiplicity m_i lives in w_o.
    function romLogupSum(uint256[41] memory e, UltraFflonkTypes.Params memory) private pure returns (uint256) {
        uint256 oneMinusQ1 = addmod(1, R - e[UltraFflonkTypes.E_Q_L], R);
        uint256 qLogupTable = mulmod(e[UltraFflonkTypes.E_Q_R], oneMinusQ1, R);
        uint256 qLogupRead = mulmod(e[UltraFflonkTypes.E_Q_4], oneMinusQ1, R);

        uint256 contribution = addmod(qLogupRead, R - mulmod(qLogupTable, e[UltraFflonkTypes.E_W_O], R), R);
        contribution = mulmod(contribution, e[UltraFflonkTypes.E_W_4], R);
        return mulmod(contribution, e[UltraFflonkTypes.E_Q_MEMORY], R);
    }
}
