// SPDX-License-Identifier: Apache-2.0
pragma solidity >=0.8.21;

import "../UltraFflonkTypes.sol";

/// @dev Port of relations/elliptic_relation.hpp — short-Weierstrass point addition/subtraction and
/// doubling: given P1 = (x1, y1) and P2 = (x2, y2), constrains P3 = P1 ± P2, or P3 = 2*P1.
library EllipticRelation {
    uint256 internal constant R = UltraFflonkTypes.R;

    /// Curve constant `b` of the embedded curve, i.e. `EllipticRelationImpl::get_curve_b()`.
    /// The Ultra flavor's FF is the BN254 scalar field, which equals `grumpkin::fq`, so the
    /// second branch is taken and this is `grumpkin::g1::curve_b`. That value is stored in
    /// grumpkin.hpp as Montgomery-form limbs {0xdd7056026000005a, 0x223fa97acb319311,
    /// 0xcc388229877910c0, 0x34394632b724eaa}; converting out of Montgomery form gives R - 17,
    /// matching the Grumpkin curve equation y^2 = x^3 - 17.
    uint256 internal constant CURVE_B = 0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593effffff0;

    /// @dev Writes each subrelation's value at xi into `out`, in SUBRELATION_PARTIAL_LENGTHS order.
    function accumulate(uint256[41] memory e, UltraFflonkTypes.Params memory, uint256[2] memory out) internal pure {
        out[0] = xCoordinate(e);
        out[1] = yCoordinate(e);
    }

    /// @dev Subrelation 0: the x-coordinate of P3. The addition branch enforces
    /// (x3 + x2 + x1)(x2 - x1)^2 - (y2^2 + y1^2 - 2*q_sign*y2*y1) = 0, the cleared-denominator form of
    /// x3 = lambda^2 - x1 - x2. The doubling branch enforces (x3 + 2*x1)*4*y1^2 - 9*x1*(y1^2 - b) = 0,
    /// which uses the on-curve identity x1^3 = y1^2 - b to keep the degree down. The two branches are
    /// selected by q_elliptic*(1 - q_is_double) and q_elliptic*q_is_double respectively.
    function xCoordinate(uint256[41] memory e) private pure returns (uint256) {
        uint256 x1 = e[UltraFflonkTypes.E_W_R];
        uint256 y1 = e[UltraFflonkTypes.E_W_O];
        uint256 y1Sqr = mulmod(y1, y1, R);
        uint256 x1Mul3 = addmod(x1, addmod(x1, x1, R), R);

        // (x3 - x1), (x3 + 2*x1) and (x3 + x2 + x1), built up from the shifted wires.
        uint256 x2SubX1 = addmod(e[UltraFflonkTypes.E_W_L_SHIFT], R - x1, R);
        uint256 x3PlusTwoX1 = addmod(addmod(e[UltraFflonkTypes.E_W_R_SHIFT], R - x1, R), x1Mul3, R);
        uint256 x3PlusX2PlusX1 = addmod(x3PlusTwoX1, x2SubX1, R);

        uint256 xAddIdentity;
        {
            uint256 y2 = e[UltraFflonkTypes.E_W_4_SHIFT];
            uint256 y2MulQSign = mulmod(y2, e[UltraFflonkTypes.E_Q_L], R);
            uint256 acc = mulmod(x3PlusX2PlusX1, mulmod(x2SubX1, x2SubX1, R), R);
            acc = addmod(acc, R - addmod(mulmod(y2, y2, R), y1Sqr, R), R);
            xAddIdentity = addmod(acc, mulmod(addmod(y2MulQSign, y2MulQSign, R), y1, R), R);
        }

        uint256 xDoubleIdentity;
        {
            // 3*x1*(y1^2 - b) stands in for 3*x1^4; tripling it gives the 9*x1^4 term.
            uint256 xPow4Mul3 = mulmod(addmod(y1Sqr, R - CURVE_B, R), x1Mul3, R);
            uint256 x1Pow4Mul9 = addmod(xPow4Mul3, addmod(xPow4Mul3, xPow4Mul3, R), R);
            uint256 y1SqrMul4 = addmod(y1Sqr, y1Sqr, R);
            y1SqrMul4 = addmod(y1SqrMul4, y1SqrMul4, R);
            xDoubleIdentity = addmod(mulmod(x3PlusTwoX1, y1SqrMul4, R), R - x1Pow4Mul9, R);
        }

        uint256 qElliptic = e[UltraFflonkTypes.E_Q_ELLIPTIC];
        uint256 doubleScaling = mulmod(qElliptic, e[UltraFflonkTypes.E_Q_M], R);
        // q_elliptic * (q_is_double - 1), i.e. the negation of the addition-branch selector.
        uint256 negNotDoubleScaling = addmod(doubleScaling, R - qElliptic, R);

        return addmod(mulmod(xDoubleIdentity, doubleScaling, R), R - mulmod(xAddIdentity, negNotDoubleScaling, R), R);
    }

    /// @dev Subrelation 1: the y-coordinate of P3, i.e. y3 = lambda*(x1 - x3) - y1 with the denominator
    /// cleared. The addition branch enforces (y3 + y1)(x2 - x1) + (x3 - x1)(q_sign*y2 - y1) = 0; the
    /// doubling branch enforces (y3 + y1)(2*y1) - (3*x1^2)(x1 - x3) = 0, accumulated in negated form.
    function yCoordinate(uint256[41] memory e) private pure returns (uint256) {
        uint256 x1 = e[UltraFflonkTypes.E_W_R];
        uint256 y1 = e[UltraFflonkTypes.E_W_O];
        uint256 y1PlusY3 = addmod(y1, e[UltraFflonkTypes.E_W_O_SHIFT], R);
        uint256 x2SubX1 = addmod(e[UltraFflonkTypes.E_W_L_SHIFT], R - x1, R);
        uint256 x3SubX1 = addmod(e[UltraFflonkTypes.E_W_R_SHIFT], R - x1, R);

        uint256 yAddIdentity;
        {
            // q_sign * y2 - y1, the numerator of lambda for the addition case.
            uint256 yDiff = addmod(mulmod(e[UltraFflonkTypes.E_W_4_SHIFT], e[UltraFflonkTypes.E_Q_L], R), R - y1, R);
            yAddIdentity = addmod(mulmod(y1PlusY3, x2SubX1, R), mulmod(x3SubX1, yDiff, R), R);
        }

        uint256 negYDoubleIdentity;
        {
            uint256 x1SqrMul3 = mulmod(addmod(x1, addmod(x1, x1, R), R), x1, R);
            negYDoubleIdentity = addmod(mulmod(x1SqrMul3, x3SubX1, R), mulmod(addmod(y1, y1, R), y1PlusY3, R), R);
        }

        uint256 qElliptic = e[UltraFflonkTypes.E_Q_ELLIPTIC];
        uint256 doubleScaling = mulmod(qElliptic, e[UltraFflonkTypes.E_Q_M], R);
        uint256 negNotDoubleScaling = addmod(doubleScaling, R - qElliptic, R);

        // Both branches enter negated, so sum them and negate once.
        uint256 sum =
            addmod(mulmod(yAddIdentity, negNotDoubleScaling, R), mulmod(negYDoubleIdentity, doubleScaling, R), R);
        return addmod(0, R - sum, R);
    }
}
