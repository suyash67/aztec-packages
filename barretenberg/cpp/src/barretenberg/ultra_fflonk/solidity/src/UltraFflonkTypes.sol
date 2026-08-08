// SPDX-License-Identifier: Apache-2.0
pragma solidity >=0.8.21;

/**
 * @title The entity layout and relation parameters the Ultra relation set reads.
 *
 * @dev The indices are `flavor/generated/ultra_flavor_generated.hpp`'s `EntityId`, verbatim. They are
 * protocol: the prover packs its groups in this order and the verifier reassembles the row from the
 * proof's evaluations into an array with this layout, so a reordering here is a reordering of the
 * proof.
 */
library UltraFflonkTypes {
    /// BN254 scalar field.
    uint256 internal constant R = 21888242871839275222246405745257275088548364400416034343698204186575808495617;

    uint256 internal constant NUM_ENTITIES = 41;

    // Precomputed (0 .. 27).
    uint256 internal constant E_SIGMA_1 = 0;
    uint256 internal constant E_SIGMA_2 = 1;
    uint256 internal constant E_SIGMA_3 = 2;
    uint256 internal constant E_SIGMA_4 = 3;
    uint256 internal constant E_ID_1 = 4;
    uint256 internal constant E_ID_2 = 5;
    uint256 internal constant E_ID_3 = 6;
    uint256 internal constant E_ID_4 = 7;
    uint256 internal constant E_LAGRANGE_FIRST = 8;
    uint256 internal constant E_LAGRANGE_LAST = 9;
    uint256 internal constant E_Q_LOOKUP = 10;
    uint256 internal constant E_TABLE_1 = 11;
    uint256 internal constant E_TABLE_2 = 12;
    uint256 internal constant E_TABLE_3 = 13;
    uint256 internal constant E_TABLE_4 = 14;
    uint256 internal constant E_Q_M = 15;
    uint256 internal constant E_Q_R = 16;
    uint256 internal constant E_Q_O = 17;
    uint256 internal constant E_Q_C = 18;
    uint256 internal constant E_Q_L = 19;
    uint256 internal constant E_Q_4 = 20;
    uint256 internal constant E_Q_ARITH = 21;
    uint256 internal constant E_Q_DELTA_RANGE = 22;
    uint256 internal constant E_Q_ELLIPTIC = 23;
    uint256 internal constant E_Q_MEMORY = 24;
    uint256 internal constant E_Q_NNF = 25;
    uint256 internal constant E_Q_POSEIDON2_EXTERNAL = 26;
    uint256 internal constant E_Q_POSEIDON2_INTERNAL = 27;

    // Witness (28 .. 35).
    uint256 internal constant E_W_L = 28;
    uint256 internal constant E_W_R = 29;
    uint256 internal constant E_W_O = 30;
    uint256 internal constant E_W_4 = 31;
    uint256 internal constant E_Z_PERM = 32;
    uint256 internal constant E_LOOKUP_INVERSES = 33;
    uint256 internal constant E_LOOKUP_READ_COUNTS = 34;
    uint256 internal constant E_LOOKUP_READ_TAGS = 35;

    // Shifted (36 .. 40): the same column read at the next row.
    uint256 internal constant E_W_L_SHIFT = 36;
    uint256 internal constant E_W_R_SHIFT = 37;
    uint256 internal constant E_W_O_SHIFT = 38;
    uint256 internal constant E_W_4_SHIFT = 39;
    uint256 internal constant E_Z_PERM_SHIFT = 40;

    /**
     * @dev `RelationParameters<FF>`, restricted to the fields the Ultra set reads. `beta_quartic` and
     * the Translator/ECCVM fields are not used by any relation here and are therefore absent.
     */
    struct Params {
        uint256 eta;
        uint256 etaTwo;
        uint256 etaThree;
        uint256 romLogupGamma;
        uint256 beta;
        uint256 gamma;
        uint256 publicInputDelta;
        uint256 betaSqr;
        uint256 betaCube;
    }
}
