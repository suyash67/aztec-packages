// SPDX-License-Identifier: Apache-2.0
pragma solidity >=0.8.21;

/**
 * @title fflonk opening verifier (ePrint 2021/1167) for transparent UltraHonk.
 *
 * @notice What fflonk buys on chain, and what it does not.
 *
 * Every KZG-based multilinear scheme here ends by batching the ~36 committed columns into one
 * commitment with a challenge drawn after the commitments are fixed, which forces a multi-scalar
 * multiplication over all of them - 50 points for Shplemini, 35 for Vela. fflonk removes that MSM
 * entirely. A commitment round of `t` columns is committed as one interleaved polynomial
 * `g(X) = sum_i f_i(X^t) X^i`, and `g mod (X^t - z) = sum_i f_i(z) X^i`, so one opening certifies
 * every column's evaluation at `z` and the residue's coefficients *are* those evaluations. Batching
 * moves out of the group and into the field.
 *
 * What is left is a BDFG21 batch over the rounds' point sets: eight scalar multiplications - four
 * round commitments, the Laurent auxiliary `h`, the quotient `W`, the linearization `W'` and the
 * generator - and one pairing check, whatever the column count.
 *
 * @dev The interpolant has a closed form. With `Z_r(X) = (X^t - z)(X^t - 1/z)` the residue is
 * `A + B X^t` where `B = (P - Q)/(z - 1/z)` and `A = P - zB`, for `P` and `Q` the evaluation vectors
 * at `z` and `1/z`. So `R_r(y)` costs two Horner passes of length `t` and no interpolation.
 */
contract FflonkOpeningVerifier {
    uint256 internal constant P = 21888242871839275222246405745257275088548364400416034343698204186575808495617;
    uint256 internal constant Q = 21888242871839275222246405745257275088696311157297823662689037894645226208583;
    uint256 internal constant SCRATCH_SIZE = 0x180;

    struct G1Point {
        uint256 x;
        uint256 y;
    }

    struct Opening {
        uint256 logN;
        uint256 numUnshifted;
        uint256 numShifted;
        uint256 previousChallenge;
        uint256[] packs; // per round, then 1 for h
        uint256[] evaluations; // the claimed sumcheck evaluations
        uint256[] claimGroup;
        uint256[] claimColumn;
        uint256[] claimEntity;
        uint256[] u; // the sumcheck challenge
        G1Point[] commitments; // one per round, then h
        uint256[] pairs; // (p, q) per packed slot, rounds in order then h
        G1Point cW;
        G1Point cWp;
    }

    struct Ctx {
        uint256 rho;
        uint256 lambda;
        uint256 z;
        uint256 zInv;
        uint256 diffInv; // 1/(z - 1/z)
        uint256 nu;
        uint256 y;
        uint256 scratch;
    }

    function verify(Opening memory o) public view returns (bool) {
        (bool ok,) = verifyWithGas(o);
        return ok;
    }

    /// `verify`, also reporting the gas spent on elliptic-curve work alone.
    function verifyWithGas(Opening memory o) public view returns (bool, uint256) {
        Ctx memory c;
        assembly {
            let free := mload(0x40)
            mstore(0x40, add(free, SCRATCH_SIZE))
            mstore(add(c, 0xe0), free)
        }
        _challenges(o, c);
        require(_laurentHolds(o, c), "fflonk: multilinear claim");

        // Shplonk scalars: nu^i * Z_T(y) / Z_i(y), with all the inversions batched into one modexp.
        uint256 rounds = o.packs.length;
        uint256[] memory scalars = new uint256[](rounds);
        uint256 zTotal = _shplonkScalars(o, c, scalars);

        uint256 ecStart = gasleft();
        G1Point memory acc;
        for (uint256 i = 0; i < rounds; i++) {
            acc = _mulAdd(acc, o.commitments[i], scalars[i], c.scratch);
        }
        acc = _mulAdd(acc, G1Point(1, 2), P - _constantTerm(o, c, scalars), c.scratch);
        acc = _mulAdd(acc, o.cW, P - zTotal, c.scratch);
        acc = _mulAdd(acc, o.cWp, c.y, c.scratch);
        bool ok = _pairing(acc, G1Point(o.cWp.x, Q - (o.cWp.y % Q)), c.scratch);
        return (ok, ecStart - gasleft());
    }

    function _challenges(Opening memory o, Ctx memory c) internal view {
        c.rho = _hashArray(o.previousChallenge, o.evaluations);
        c.lambda = _hash(c.rho);
        c.z = _hash3(c.lambda, o.commitments[o.packs.length - 1].x, o.commitments[o.packs.length - 1].y, c.scratch);
        require(c.z != 0 && c.z != 1 && c.z != P - 1, "degenerate z");
        c.nu = _hashArray(c.z, o.pairs);
        c.y = _hash3(c.nu, o.cW.x, o.cW.y, c.scratch);

        // 1/z and 1/(z - 1/z) from a single inversion: t = 1/(z(z^2-1)) gives 1/z = t(z^2-1) and
        // 1/(z - 1/z) = t z^2.
        uint256 zSquaredMinusOne = addmod(mulmod(c.z, c.z, P), P - 1, P);
        uint256 t = _invert(mulmod(c.z, zSquaredMinusOne, P), c.scratch);
        c.zInv = mulmod(t, zSquaredMinusOne, P);
        c.diffInv = mulmod(t, mulmod(c.z, c.z, P), P);
    }

    /// The Laurent identity that turns the certified univariate evaluations into the multilinear
    /// claim - Vela's reduction, but every term computed in the field from certified values.
    function _laurentHolds(Opening memory o, Ctx memory c) internal pure returns (bool) {
        uint256[4] memory v; // v0a, v1a, v0b, v1b
        uint256 yCombined;
        uint256 power = 1;
        for (uint256 i = 0; i < o.claimGroup.length; i++) {
            uint256 slot = _slot(o, o.claimGroup[i], o.claimColumn[i]);
            uint256 offset = i < o.numUnshifted ? 0 : 2;
            v[offset] = addmod(v[offset], mulmod(power, o.pairs[slot], P), P);
            v[offset + 1] = addmod(v[offset + 1], mulmod(power, o.pairs[slot + 1], P), P);
            uint256 weight = i < o.numUnshifted ? power : mulmod(c.lambda, power, P);
            yCombined = addmod(yCombined, mulmod(weight, o.evaluations[o.claimEntity[i]], P), P);
            power = mulmod(power, c.rho, P);
        }
        uint256 gZ = addmod(v[0], mulmod(mulmod(c.lambda, c.zInv, P), v[2], P), P);
        uint256 gZInv = addmod(v[1], mulmod(mulmod(c.lambda, c.z, P), v[3], P), P);
        uint256 left = addmod(mulmod(gZ, _tensor(o.u, c.zInv), P), mulmod(gZInv, _tensor(o.u, c.z), P), P);
        left = addmod(left, P - addmod(yCombined, yCombined, P), P);
        uint256 hSlot = o.pairs.length - 2;
        uint256 right = addmod(mulmod(c.z, o.pairs[hSlot], P), mulmod(c.zInv, o.pairs[hSlot + 1], P), P);
        return left == right;
    }

    /// `nu^i Z_T(y) / Z_i(y)` for every round, and `Z_T(y)`. One modexp covers all the inversions.
    function _shplonkScalars(Opening memory o, Ctx memory c, uint256[] memory scalars)
        internal
        view
        returns (uint256 zTotal)
    {
        uint256 rounds = o.packs.length;
        uint256[] memory zAtY = new uint256[](rounds);
        zTotal = 1;
        for (uint256 i = 0; i < rounds; i++) {
            uint256 yPow = _pow(c.y, o.packs[i]);
            zAtY[i] = mulmod(addmod(yPow, P - c.z, P), addmod(yPow, P - c.zInv, P), P);
            zTotal = mulmod(zTotal, zAtY[i], P);
        }
        // Montgomery's trick: with Z_T = prod Z_i, the scalar nu^i Z_T / Z_i is nu^i times the
        // product of the other Z_j, which a forward and backward pass produce without any inversion.
        uint256[] memory prefix = new uint256[](rounds + 1);
        prefix[0] = 1;
        for (uint256 i = 0; i < rounds; i++) {
            prefix[i + 1] = mulmod(prefix[i], zAtY[i], P);
        }
        uint256 suffix = 1;
        uint256 nuPower = 1;
        for (uint256 i = rounds; i > 0; i--) {
            scalars[i - 1] = mulmod(prefix[i - 1], suffix, P);
            suffix = mulmod(suffix, zAtY[i - 1], P);
        }
        for (uint256 i = 0; i < rounds; i++) {
            scalars[i] = mulmod(scalars[i], nuPower, P);
            nuPower = mulmod(nuPower, c.nu, P);
        }
    }

    /// `sum_i scalar_i R_i(y)`, the generator's coefficient in the Shplonk accumulation.
    function _constantTerm(Opening memory o, Ctx memory c, uint256[] memory scalars)
        internal
        pure
        returns (uint256 total)
    {
        uint256 base = 0;
        for (uint256 i = 0; i < o.packs.length; i++) {
            uint256 t = o.packs[i];
            uint256 a = 0;
            uint256 b = 0;
            // Horner over A and B, whose coefficients follow from the two evaluation vectors.
            for (uint256 j = t; j > 0; j--) {
                uint256 slot = base + (j - 1) * 2;
                uint256 bj = mulmod(addmod(o.pairs[slot], P - o.pairs[slot + 1], P), c.diffInv, P);
                uint256 aj = addmod(o.pairs[slot], P - mulmod(c.z, bj, P), P);
                a = addmod(mulmod(a, c.y, P), aj, P);
                b = addmod(mulmod(b, c.y, P), bj, P);
            }
            total = addmod(total, mulmod(scalars[i], addmod(a, mulmod(_pow(c.y, t), b, P), P), P), P);
            base += t * 2;
        }
    }

    /// Index of `(p, q)` for a column, in the flat evaluation-pair array.
    function _slot(Opening memory o, uint256 group, uint256 column) internal pure returns (uint256 slot) {
        for (uint256 g = 0; g < group; g++) {
            slot += o.packs[g] * 2;
        }
        slot += column * 2;
    }

    function _tensor(uint256[] memory r, uint256 x) internal pure returns (uint256 result) {
        result = 1;
        uint256 power = x;
        for (uint256 i = 0; i < r.length; i++) {
            result = mulmod(result, addmod(addmod(mulmod(r[i], power, P), P - r[i], P), 1, P), P);
            power = mulmod(power, power, P);
        }
    }

    function _pow(uint256 base, uint256 exponent) internal pure returns (uint256 result) {
        result = 1;
        uint256 square = base;
        uint256 e = exponent;
        while (e > 0) {
            if (e & 1 == 1) {
                result = mulmod(result, square, P);
            }
            square = mulmod(square, square, P);
            e >>= 1;
        }
    }

    function _mulAdd(G1Point memory acc, G1Point memory point, uint256 scalar, uint256 scratch)
        internal
        view
        returns (G1Point memory out)
    {
        assembly {
            mstore(scratch, mload(acc))
            mstore(add(scratch, 0x20), mload(add(acc, 0x20)))
            mstore(add(scratch, 0x40), mload(point))
            mstore(add(scratch, 0x60), mload(add(point, 0x20)))
            mstore(add(scratch, 0x80), scalar)
            let ok := staticcall(gas(), 7, add(scratch, 0x40), 0x60, add(scratch, 0x40), 0x40)
            ok := and(ok, staticcall(gas(), 6, scratch, 0x80, scratch, 0x40))
            if iszero(ok) { revert(0, 0) }
            out := mload(0x40)
            mstore(0x40, add(out, 0x40))
            mstore(out, mload(scratch))
            mstore(add(out, 0x20), mload(add(scratch, 0x20)))
        }
    }

    function _hashArray(uint256 previous, uint256[] memory data) internal pure returns (uint256 result) {
        assembly {
            let saved := mload(data)
            mstore(data, previous)
            result := mod(keccak256(data, add(0x20, mul(saved, 0x20))), P)
            mstore(data, saved)
        }
    }

    function _hash(uint256 a) internal pure returns (uint256 result) {
        assembly {
            mstore(0x00, a)
            result := mod(keccak256(0x00, 0x20), P)
        }
    }

    function _hash3(uint256 a, uint256 b, uint256 d, uint256 scratch) internal pure returns (uint256 result) {
        assembly {
            mstore(scratch, a)
            mstore(add(scratch, 0x20), b)
            mstore(add(scratch, 0x40), d)
            result := mod(keccak256(scratch, 0x60), P)
        }
    }

    function _invert(uint256 a, uint256 scratch) internal view returns (uint256 result) {
        uint256 p = P;
        assembly {
            mstore(scratch, 0x20)
            mstore(add(scratch, 0x20), 0x20)
            mstore(add(scratch, 0x40), 0x20)
            mstore(add(scratch, 0x60), a)
            mstore(add(scratch, 0x80), sub(p, 2))
            mstore(add(scratch, 0xa0), p)
            if iszero(staticcall(gas(), 0x05, scratch, 0xc0, scratch, 0x20)) { revert(0, 0) }
            result := mload(scratch)
        }
    }

    function _pairing(G1Point memory lhs, G1Point memory rhs, uint256 scratch) internal view returns (bool ok) {
        assembly {
            mstore(scratch, mload(lhs))
            mstore(add(scratch, 0x20), mload(add(lhs, 0x20)))
            mstore(add(scratch, 0x40), 0x198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c2)
            mstore(add(scratch, 0x60), 0x1800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed)
            mstore(add(scratch, 0x80), 0x090689d0585ff075ec9e99ad690c3395bc4b313370b38ef355acdadcd122975b)
            mstore(add(scratch, 0xa0), 0x12c85ea5db8c6deb4aab71808dcb408fe3d1e7690c43d37b4ce6cc0166fa7daa)
            mstore(add(scratch, 0xc0), mload(rhs))
            mstore(add(scratch, 0xe0), mload(add(rhs, 0x20)))
            mstore(add(scratch, 0x100), 0x260e01b251f6f1c7e7ff4e580791dee8ea51d87a358e038b4efe30fac09383c1)
            mstore(add(scratch, 0x120), 0x0118c4d5b837bcc2bc89b5b398b5974e9f5944073b32078b7e231fec938883b0)
            mstore(add(scratch, 0x140), 0x04fc6369f7110fe3d25156c1bb9a72859cf2a04641f99ba4ee413c80da6a5fe4)
            mstore(add(scratch, 0x160), 0x22febda3c0c0632a56475b4214e5615e11e6dd3f96e6cea2854a87d4dacc5e55)
            ok := staticcall(gas(), 8, scratch, 0x180, scratch, 0x20)
            ok := and(ok, mload(scratch))
        }
    }
}
