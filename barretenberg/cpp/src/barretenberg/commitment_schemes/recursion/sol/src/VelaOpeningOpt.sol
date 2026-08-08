// SPDX-License-Identifier: Apache-2.0
pragma solidity >=0.8.21;

import {VelaOpeningVerifier} from "./VelaOpening.sol";

/**
 * @title Vela opening verifier, with the optimizations barretenberg's `--optimized` Honk verifier uses.
 *
 * @notice Same protocol and same inputs as `VelaOpeningVerifier`; only the implementation changes, so
 * the pair measures what assembly buys rather than what a different algorithm would.
 *
 * Four techniques carry essentially all of it, taken from
 * `barretenberg/sol/src/honk/optimised/honk-optimized.sol.template`:
 *
 *  1. **Precompiles called straight into a fixed scratch region.** The accumulator, the staged point
 *     and the scalar sit at fixed offsets chosen so that `ecMul`'s output region is `ecAdd`'s second
 *     operand. No `abi.encodePacked`, no `abi.decode`, no allocation per curve operation - the whole
 *     MSM allocates once.
 *  2. **One keccak over contiguous memory.** A challenge's preimage is `previous || round buffer`,
 *     and the round buffer is already a memory array; writing the previous challenge into the word
 *     ahead of its data makes the preimage contiguous, so the hash is one `keccak256(ptr, len)`
 *     instead of a quadratic chain of `abi.encodePacked` copies.
 *  3. **One modular inversion for the whole opening.** Vela needs `1/z` and `1/(z - 1/z)`; both come
 *     from inverting `z(z^2-1)` once, since `t = 1/(z(z^2-1))` gives `1/z = t(z^2-1)` and
 *     `1/(z - 1/z) = t z^2`. That is Montgomery's trick specialized to this pair, and it halves the
 *     `modexp` calls.
 *  4. **A preassembled pairing input.** The 384-byte argument is written once into scratch rather
 *     than built by `abi.encodePacked`.
 *
 * The optimizations that dominate the Honk verifier but do *not* apply here are the ones aimed at its
 * field work: unrolling `log n` sumcheck rounds, inlining the verification key, and batching ~40
 * barycentric and Shplemini denominators into a single inversion. Vela's opening has almost no field
 * work to fold - its cost is the MSM and the pairing.
 */
contract VelaOpeningVerifierOpt {
    uint256 internal constant P = 21888242871839275222246405745257275088548364400416034343698204186575808495617;
    uint256 internal constant Q = 21888242871839275222246405745257275088696311157297823662689037894645226208583;

    // Scratch layout, mirroring the template's: the accumulator occupies [0x00, 0x40) of the region,
    // the staged point and scalar [0x40, 0xa0). `ecMul` reads [0x40, 0xa0) and writes [0x40, 0x80);
    // `ecAdd` then reads [0x00, 0x80) and writes [0x00, 0x40). One buffer, no moves.
    uint256 internal constant SCRATCH_SIZE = 0x180;

    /// The opening's challenges and the two inverses derived from them. Held in memory so the
    /// verifier keeps one stack slot rather than a dozen.
    struct Ctx {
        uint256 rho;
        uint256 lambda;
        uint256 z;
        uint256 zInv;
        uint256 diffInv;
        uint256 alpha;
        uint256 alphaSquared;
        uint256 zeta;
        uint256 scratch;
    }

    function verify(VelaOpeningVerifier.Opening memory o) public view returns (bool) {
        (bool ok,) = verifyWithGas(o);
        return ok;
    }

    /// `verify`, additionally reporting the gas spent on elliptic-curve work alone, so it can be set
    /// against the same split of the generated Shplemini verifier.
    function verifyWithGas(VelaOpeningVerifier.Opening memory o) public view returns (bool, uint256) {
        Ctx memory c;
        assembly {
            let free := mload(0x40)
            mstore(0x40, add(free, SCRATCH_SIZE))
            mstore(add(c, 0x100), free)
        }
        _challenges(o, c);

        uint256[] memory rhoPowers = new uint256[](o.numUnshifted + o.numShifted);
        uint256 w1 = _recoverW1(o, c, _batchScalars(o, c, rhoPowers));

        // C_batched = C_A + alpha C_B + alpha^2 C_h, with C_A and C_B the rho-combinations. Folding
        // alpha into C_B's scalars keeps it to a single accumulator and saves a scalar multiplication.
        VelaOpeningVerifier.G1Point memory acc;
        uint256 ecStart = gasleft();
        acc = _msm(acc, o.commitments, rhoPowers, 0, o.numUnshifted, 1, c.scratch);
        acc = _msm(acc, o.commitments, rhoPowers, o.numUnshifted, o.numShifted, c.alpha, c.scratch);
        acc = _mulAdd(acc, o.cH, c.alphaSquared, c.scratch);

        // P0 = C_batched - Z(zeta) C_q - R(zeta) [1] + zeta pi_L. Negating a scalar is free where
        // negating a point is not, so the subtractions ride on the same accumulate step.
        (uint256 rAtZeta, uint256 zAtZeta) = _line(o, c, w1);
        acc = _mulAdd(acc, o.cQ, P - zAtZeta, c.scratch);
        acc = _mulAdd(acc, VelaOpeningVerifier.G1Point(1, 2), P - rAtZeta, c.scratch);
        acc = _mulAdd(acc, o.piL, c.zeta, c.scratch);
        bool ok = _pairing(acc, VelaOpeningVerifier.G1Point(o.piL.x, Q - (o.piL.y % Q)), c.scratch);
        return (ok, ecStart - gasleft());
    }

    /// The Fiat-Shamir challenges, plus 1/z and 1/(z - 1/z) from a single inversion.
    function _challenges(VelaOpeningVerifier.Opening memory o, Ctx memory c) internal view {
        c.rho = _hashArray(o.previousChallenge, o.evaluations);
        c.lambda = _hash(c.rho);
        c.z = _hash3(c.lambda, o.cH.x, o.cH.y, c.scratch);
        require(c.z != 0 && c.z != 1 && c.z != P - 1, "degenerate z");
        c.alpha = _hash6(c.z, o.v0a, o.v1a, o.v0b, o.v1b, o.w0);
        c.alphaSquared = mulmod(c.alpha, c.alpha, P);
        c.zeta = _hash3(c.alpha, o.cQ.x, o.cQ.y, c.scratch);

        uint256 zSquaredMinusOne = addmod(mulmod(c.z, c.z, P), P - 1, P);
        uint256 t = _invert(mulmod(c.z, zSquaredMinusOne, P), c.scratch);
        c.zInv = mulmod(t, zSquaredMinusOne, P);
        c.diffInv = mulmod(t, mulmod(c.z, c.z, P), P);
    }

    /// R(zeta) for the line through the batched polynomial's two evaluations, and Z(zeta).
    function _line(VelaOpeningVerifier.Opening memory o, Ctx memory c, uint256 w1)
        internal
        pure
        returns (uint256 rAtZeta, uint256 zAtZeta)
    {
        uint256 batched0 =
            addmod(addmod(o.v0a, mulmod(c.alpha, o.v0b, P), P), mulmod(c.alphaSquared, o.w0, P), P);
        uint256 batched1 =
            addmod(addmod(o.v1a, mulmod(c.alpha, o.v1b, P), P), mulmod(c.alphaSquared, w1, P), P);
        uint256 slope = mulmod(addmod(batched0, P - batched1, P), c.diffInv, P);
        rAtZeta = addmod(addmod(batched0, P - mulmod(slope, c.z, P), P), mulmod(slope, c.zeta, P), P);
        zAtZeta = mulmod(addmod(c.zeta, P - c.z, P), addmod(c.zeta, P - c.zInv, P), P);
    }

    /// The rho powers, and the combined claim value they weight (the shifted half by lambda too).
    function _batchScalars(VelaOpeningVerifier.Opening memory o, Ctx memory c, uint256[] memory rhoPowers)
        internal
        pure
        returns (uint256 yCombined)
    {
        uint256 power = 1;
        for (uint256 i = 0; i < rhoPowers.length; i++) {
            rhoPowers[i] = power;
            uint256 weight = i < o.numUnshifted ? power : mulmod(c.lambda, power, P);
            yCombined = addmod(yCombined, mulmod(weight, o.evaluations[o.claimEntity[i]], P), P);
            power = mulmod(power, c.rho, P);
        }
    }

    function _recoverW1(VelaOpeningVerifier.Opening memory o, Ctx memory c, uint256 yCombined)
        internal
        pure
        returns (uint256)
    {
        uint256 gZ = addmod(o.v0a, mulmod(mulmod(c.lambda, c.zInv, P), o.v0b, P), P);
        uint256 gZInv = addmod(o.v1a, mulmod(mulmod(c.lambda, c.z, P), o.v1b, P), P);
        uint256 inner =
            addmod(mulmod(gZ, _tensor(o.u, c.zInv), P), mulmod(gZInv, _tensor(o.u, c.z), P), P);
        inner = addmod(inner, P - addmod(yCombined, yCombined, P), P);
        inner = addmod(inner, P - mulmod(c.z, o.w0, P), P);
        return mulmod(c.z, inner, P);
    }

    function _tensor(uint256[] memory r, uint256 x) internal pure returns (uint256 result) {
        result = 1;
        uint256 power = x;
        for (uint256 i = 0; i < r.length; i++) {
            result = mulmod(result, addmod(addmod(mulmod(r[i], power, P), P - r[i], P), 1, P), P);
            power = mulmod(power, power, P);
        }
    }

    /// `acc + sum_{i in [start, start+count)} (weight * scalars[i]) * points[i]`.
    function _msm(
        VelaOpeningVerifier.G1Point memory acc,
        VelaOpeningVerifier.G1Point[] memory points,
        uint256[] memory scalars,
        uint256 start,
        uint256 count,
        uint256 weight,
        uint256 scratch
    ) internal view returns (VelaOpeningVerifier.G1Point memory out) {
        uint256 p = P;
        assembly {
            mstore(scratch, mload(acc))
            mstore(add(scratch, 0x20), mload(add(acc, 0x20)))
            let ok := 1
            let pointBase := add(add(points, 0x20), mul(start, 0x20))
            let scalarBase := add(add(scalars, 0x20), mul(start, 0x20))
            for { let i := 0 } lt(i, count) { i := add(i, 1) } {
                let point := mload(add(pointBase, mul(i, 0x20)))
                mstore(add(scratch, 0x40), mload(point))
                mstore(add(scratch, 0x60), mload(add(point, 0x20)))
                mstore(add(scratch, 0x80), mulmod(mload(add(scalarBase, mul(i, 0x20))), weight, p))
                ok := and(ok, staticcall(gas(), 7, add(scratch, 0x40), 0x60, add(scratch, 0x40), 0x40))
                ok := and(ok, staticcall(gas(), 6, scratch, 0x80, scratch, 0x40))
            }
            if iszero(ok) { revert(0, 0) }
            out := mload(0x40)
            mstore(0x40, add(out, 0x40))
            mstore(out, mload(scratch))
            mstore(add(out, 0x20), mload(add(scratch, 0x20)))
        }
    }

    /// `acc + scalar * point`.
    function _mulAdd(
        VelaOpeningVerifier.G1Point memory acc,
        VelaOpeningVerifier.G1Point memory point,
        uint256 scalar,
        uint256 scratch
    ) internal view returns (VelaOpeningVerifier.G1Point memory out) {
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

    /// `keccak(previous || data) mod P`, hashing the array in place.
    function _hashArray(uint256 previous, uint256[] memory data) internal pure returns (uint256 result) {
        assembly {
            let saved := mload(data) // the length slot, which becomes the previous-challenge word
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

    /// Hashes into the caller's scratch region: 0x40 holds the free-memory pointer and must not be
    /// used as a hash buffer.
    function _hash3(uint256 a, uint256 b, uint256 c, uint256 scratch) internal pure returns (uint256 result) {
        assembly {
            mstore(scratch, a)
            mstore(add(scratch, 0x20), b)
            mstore(add(scratch, 0x40), c)
            result := mod(keccak256(scratch, 0x60), P)
        }
    }

    function _hash6(uint256 a, uint256 b, uint256 c, uint256 d, uint256 e, uint256 f)
        internal
        pure
        returns (uint256)
    {
        return uint256(keccak256(abi.encodePacked(a, b, c, d, e, f))) % P;
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

    function _pairing(
        VelaOpeningVerifier.G1Point memory lhs,
        VelaOpeningVerifier.G1Point memory rhs,
        uint256 scratch
    ) internal view returns (bool ok) {
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
