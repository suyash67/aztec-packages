// SPDX-License-Identifier: Apache-2.0
pragma solidity >=0.8.21;

import "../UltraFflonkTypes.sol";

/// @dev Port of relations/logderiv_lookup_relation.hpp — the log-derivative lookup argument: reads
/// from a table of up to three columns, expressed as a sum of inverses of the read and table terms.
library LogDerivLookupRelation {
    uint256 internal constant R = UltraFflonkTypes.R;

    /// @dev Writes each subrelation's value at xi into `out`, in SUBRELATION_PARTIAL_LENGTHS order.
    function accumulate(uint256[41] memory e, UltraFflonkTypes.Params memory p, uint256[3] memory out) internal pure {
        uint256 tableTerm = computeTableTerm(e, p);
        uint256 lookupTerm = computeLookupTerm(e, p);
        uint256 inverses = e[UltraFflonkTypes.E_LOOKUP_INVERSES];

        // Subrelation 0: the polynomial of inverses I is correctly formed, i.e.
        // I * lookup_term * table_term == inverse_exists. The right-hand side rather than 1 because I
        // is set to 0 on rows that neither perform a read nor hold read table data.
        {
            uint256 readTags = e[UltraFflonkTypes.E_LOOKUP_READ_TAGS];
            uint256 readSelector = e[UltraFflonkTypes.E_Q_LOOKUP];
            // inverse_exists = 1 - (1 - read_tag) * (1 - is_read_gate), expanded assuming both are boolean
            uint256 inverseExists = addmod(addmod(readTags, readSelector, R), R - mulmod(readTags, readSelector, R), R);
            uint256 product = mulmod(mulmod(lookupTerm, tableTerm, R), inverses, R);
            out[0] = addmod(product, R - inverseExists, R);
        }

        // Subrelation 1: the lookup identity itself. Summed over the trace this is
        // sum(q_lookup / lookup_term - read_count / table_term), rearranged over the common
        // denominator I = 1 / (lookup_term * table_term). Linearly dependent: it only vanishes when
        // aggregated across all rows, not row by row.
        {
            uint256 readSelector = e[UltraFflonkTypes.E_Q_LOOKUP];
            uint256 readCounts = e[UltraFflonkTypes.E_LOOKUP_READ_COUNTS];
            uint256 tmp = addmod(mulmod(readSelector, tableTerm, R), R - mulmod(readCounts, lookupTerm, R), R);
            out[1] = mulmod(tmp, inverses, R);
        }

        // Subrelation 2: read_tag is boolean. Required because inverse_exists above is only the
        // intended OR when read_tag is constrained to {0, 1}; otherwise a prover could pick any value
        // for it on rows where the read selector is off.
        {
            uint256 readTag = e[UltraFflonkTypes.E_LOOKUP_READ_TAGS];
            out[2] = addmod(mulmod(readTag, readTag, R), R - readTag, R);
        }
    }

    /// @dev table_1 + gamma + table_2 * beta + table_3 * beta^2 + table_4 * beta^3, where table_1..3
    /// are the lookup table's columns and table_4 is the table's unique identifier.
    function computeTableTerm(uint256[41] memory e, UltraFflonkTypes.Params memory p)
        private
        pure
        returns (uint256 result)
    {
        result = mulmod(e[UltraFflonkTypes.E_TABLE_2], p.beta, R);
        result = addmod(result, mulmod(e[UltraFflonkTypes.E_TABLE_3], p.betaSqr, R), R);
        result = addmod(result, mulmod(e[UltraFflonkTypes.E_TABLE_4], p.betaCube, R), R);
        result = addmod(result, e[UltraFflonkTypes.E_TABLE_1], R);
        result = addmod(result, p.gamma, R);
    }

    /// @dev The read term, built from the same challenge powers as the table term. The looked-up
    /// values are not the wires themselves: wires hold successive accumulators, and column i's entry
    /// is recovered as w_i + (negated step size) * w_i_shift. The step sizes live in q_r/q_m/q_c and
    /// the table identifier in q_o.
    function computeLookupTerm(uint256[41] memory e, UltraFflonkTypes.Params memory p)
        private
        pure
        returns (uint256 result)
    {
        uint256 derivedEntry2 = addmod(
            mulmod(e[UltraFflonkTypes.E_Q_M], e[UltraFflonkTypes.E_W_R_SHIFT], R), e[UltraFflonkTypes.E_W_R], R
        );
        uint256 derivedEntry3 =
            addmod(mulmod(e[UltraFflonkTypes.E_Q_C], e[UltraFflonkTypes.E_W_O_SHIFT], R), e[UltraFflonkTypes.E_W_O], R);

        result = addmod(mulmod(derivedEntry2, p.beta, R), mulmod(derivedEntry3, p.betaSqr, R), R);

        {
            uint256 derivedEntry1 = addmod(
                mulmod(e[UltraFflonkTypes.E_Q_R], e[UltraFflonkTypes.E_W_L_SHIFT], R),
                addmod(e[UltraFflonkTypes.E_W_L], p.gamma, R),
                R
            );
            uint256 tableIndexEntry = mulmod(e[UltraFflonkTypes.E_Q_O], p.betaCube, R);
            result = addmod(result, addmod(derivedEntry1, tableIndexEntry, R), R);
        }
    }
}
