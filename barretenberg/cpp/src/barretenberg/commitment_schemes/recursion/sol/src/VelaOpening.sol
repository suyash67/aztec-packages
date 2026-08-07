// SPDX-License-Identifier: Apache-2.0
pragma solidity >=0.8.21;

/**
 * @title Vela opening verifier (ePrint 2026/1438 §2.2) for transparent UltraHonk.
 *
 * @notice The step that differs between UltraHonk with Shplemini+KZG and UltraHonk with Vela.
 * Everything before it - the transcript, sumcheck and relation evaluation - is identical for the two
 * because they share a flavor, so this contract plus the generated Shplemini verifier are what a gas
 * comparison has to weigh against each other.
 *
 * @dev The multilinear claim f(r) = y is the constant coefficient of the Laurent product
 * H(X) = f_v(X)·T_r(1/X). Its inversion-symmetric residual yields one ordinary polynomial h whose
 * commitment plus evaluations at z and 1/z certify the claim; the verifier recovers the fourth
 * evaluation h(1/z) from the Laurent identity and closes with a single degree-2 vanishing check,
 * linearized at a further challenge because barretenberg's SRS publishes only [1]_2 and [tau]_2.
 *
 * The whole claim batch (every unshifted column, then every to-be-shifted one) collapses into one
 * commitment: the rho-combination of the unshifted claims and lambda·X^-1 times that of the shifted
 * ones. The elliptic-curve work is therefore one MSM over the claim commitments plus three further
 * scalar multiplications and one pairing.
 */
contract VelaOpeningVerifier {
    uint256 internal constant P = 21888242871839275222246405745257275088548364400416034343698204186575808495617;
    uint256 internal constant Q = 21888242871839275222246405745257275088696311157297823662689037894645226208583;

    struct G1Point {
        uint256 x;
        uint256 y;
    }

    /// Everything the opening step consumes, already in memory - the shape `verifyShplemini` is
    /// handed in the generated verifier, so the two measurements start from the same place.
    struct Opening {
        uint256 logN;
        uint256 numUnshifted;
        uint256 numShifted;
        uint256 previousChallenge; // transcript state entering the opening phase
        uint256[] evaluations; // every claimed sumcheck evaluation (the round buffer rho absorbs)
        uint256[] claimEntity; // claim index -> entity index into `evaluations`
        uint256[] u; // the sumcheck challenge
        G1Point[] commitments; // claim commitments, unshifted then shifted
        G1Point cH;
        uint256 v0a;
        uint256 v1a;
        uint256 v0b;
        uint256 v1b;
        uint256 w0;
        G1Point cQ;
        G1Point piL;
    }

    /// The opening phase's Fiat-Shamir challenges, drawn in transcript order.
    struct Challenges {
        uint256 rho;
        uint256 lambda;
        uint256 z;
        uint256 zInv;
        uint256 alpha;
        uint256 zeta;
    }

    function verify(Opening memory o) public view returns (bool) {
        Challenges memory ch = _challenges(o);
        (G1Point memory cA, G1Point memory cB, uint256 yCombined) = _batch(o, ch);
        uint256 w1 = _recoverW1(o, ch, yCombined);
        return _close(o, ch, cA, cB, w1);
    }

    /// challenge = keccak(previous || round buffer) mod P, the convention bb's `KeccakTranscript`
    /// and the generated Honk verifier both use.
    function _challenges(Opening memory o) internal view returns (Challenges memory ch) {
        ch.rho = _hashWith(o.previousChallenge, o.evaluations);
        ch.lambda = _hash1(ch.rho);
        ch.z = _hash3(ch.lambda, o.cH.x, o.cH.y);
        ch.alpha = _hash6(ch.z, o.v0a, o.v1a, o.v0b, o.v1b, o.w0);
        ch.zeta = _hash3(ch.alpha, o.cQ.x, o.cQ.y);
        require(ch.z != 0 && ch.z != 1 && ch.z != P - 1, "degenerate z");
        ch.zInv = _invert(ch.z);
    }

    /// The rho-combination of the unshifted and to-be-shifted claim commitments, and the combined
    /// claim value (the shifted half weighted by lambda).
    function _batch(Opening memory o, Challenges memory ch)
        internal
        view
        returns (G1Point memory cA, G1Point memory cB, uint256 yCombined)
    {
        uint256 rhoPower = 1;
        for (uint256 i = 0; i < o.numUnshifted; i++) {
            cA = _add(cA, _mul(o.commitments[i], rhoPower));
            yCombined = addmod(yCombined, mulmod(rhoPower, o.evaluations[o.claimEntity[i]], P), P);
            rhoPower = mulmod(rhoPower, ch.rho, P);
        }
        for (uint256 l = 0; l < o.numShifted; l++) {
            uint256 idx = o.numUnshifted + l;
            cB = _add(cB, _mul(o.commitments[idx], rhoPower));
            uint256 weighted = mulmod(mulmod(ch.lambda, rhoPower, P), o.evaluations[o.claimEntity[idx]], P);
            yCombined = addmod(yCombined, weighted, P);
            rhoPower = mulmod(rhoPower, ch.rho, P);
        }
    }

    /// The alpha-batched polynomial, the line through its two evaluations, and the pairing check
    /// e(P0, [1]_2) e(P1, [tau]_2) == 1 that certifies the degree-2 vanishing.
    function _close(
        Opening memory o,
        Challenges memory ch,
        G1Point memory cA,
        G1Point memory cB,
        uint256 w1
    ) internal view returns (bool) {
        uint256 alphaSquared = mulmod(ch.alpha, ch.alpha, P);
        G1Point memory cBatched = _add(_add(cA, _mul(cB, ch.alpha)), _mul(o.cH, alphaSquared));
        (uint256 rAtZeta, uint256 zAtZeta) = _line(o, ch, w1, alphaSquared);
        return _finalPairing(o, ch, cBatched, rAtZeta, zAtZeta);
    }

    /// R(zeta) for the line through (z, batched(z)) and (1/z, batched(1/z)), and Z(zeta).
    function _line(Opening memory o, Challenges memory ch, uint256 w1, uint256 alphaSquared)
        internal
        view
        returns (uint256 rAtZeta, uint256 zAtZeta)
    {
        uint256 batched0 =
            addmod(addmod(o.v0a, mulmod(ch.alpha, o.v0b, P), P), mulmod(alphaSquared, o.w0, P), P);
        uint256 batched1 =
            addmod(addmod(o.v1a, mulmod(ch.alpha, o.v1b, P), P), mulmod(alphaSquared, w1, P), P);
        uint256 slope =
            mulmod(addmod(batched0, P - batched1, P), _invert(addmod(ch.z, P - ch.zInv, P)), P);
        rAtZeta = addmod(addmod(batched0, P - mulmod(slope, ch.z, P), P), mulmod(slope, ch.zeta, P), P);
        zAtZeta = mulmod(addmod(ch.zeta, P - ch.z, P), addmod(ch.zeta, P - ch.zInv, P), P);
    }

    function _finalPairing(
        Opening memory o,
        Challenges memory ch,
        G1Point memory cBatched,
        uint256 rAtZeta,
        uint256 zAtZeta
    ) internal view returns (bool) {
        G1Point memory p0 = _add(cBatched, _neg(_mul(o.cQ, zAtZeta)));
        p0 = _add(p0, _neg(_mul(G1Point(1, 2), rAtZeta)));
        p0 = _add(p0, _mul(o.piL, ch.zeta));
        return _pairing(p0, _neg(o.piL));
    }

    /// w1 = z (g(z) T_r(1/z) + g(1/z) T_r(z) - 2Y - z w0): the paper's recovery of the fourth
    /// evaluation, which is what keeps it out of the proof.
    function _recoverW1(Opening memory o, Challenges memory ch, uint256 yCombined)
        internal
        pure
        returns (uint256)
    {
        uint256 gZ = addmod(o.v0a, mulmod(mulmod(ch.lambda, ch.zInv, P), o.v0b, P), P);
        uint256 gZInv = addmod(o.v1a, mulmod(mulmod(ch.lambda, ch.z, P), o.v1b, P), P);
        uint256 inner =
            addmod(mulmod(gZ, _tensor(o.u, ch.zInv), P), mulmod(gZInv, _tensor(o.u, ch.z), P), P);
        inner = addmod(inner, P - addmod(yCombined, yCombined, P), P);
        inner = addmod(inner, P - mulmod(ch.z, o.w0, P), P);
        return mulmod(ch.z, inner, P);
    }

    /// T_r(x) = prod_k ((1 - r_k) + r_k x^(2^k)).
    function _tensor(uint256[] memory r, uint256 x) internal pure returns (uint256 result) {
        result = 1;
        uint256 power = x;
        for (uint256 i = 0; i < r.length; i++) {
            uint256 term = addmod(mulmod(r[i], power, P), P - r[i], P);
            result = mulmod(result, addmod(term, 1, P), P);
            power = mulmod(power, power, P);
        }
    }

    function _hashWith(uint256 previous, uint256[] memory data) internal pure returns (uint256) {
        bytes memory buffer = abi.encodePacked(previous);
        for (uint256 i = 0; i < data.length; i++) {
            buffer = abi.encodePacked(buffer, data[i]);
        }
        return uint256(keccak256(buffer)) % P;
    }

    function _hash1(uint256 a) internal pure returns (uint256) {
        return uint256(keccak256(abi.encodePacked(a))) % P;
    }

    function _hash3(uint256 a, uint256 b, uint256 c) internal pure returns (uint256) {
        return uint256(keccak256(abi.encodePacked(a, b, c))) % P;
    }

    function _hash6(uint256 a, uint256 b, uint256 c, uint256 d, uint256 e, uint256 f)
        internal
        pure
        returns (uint256)
    {
        return uint256(keccak256(abi.encodePacked(a, b, c, d, e, f))) % P;
    }

    function _invert(uint256 a) internal view returns (uint256 result) {
        bytes memory input = abi.encodePacked(uint256(32), uint256(32), uint256(32), a, P - 2, P);
        (bool ok, bytes memory output) = address(0x05).staticcall(input);
        require(ok, "modexp");
        result = abi.decode(output, (uint256));
    }

    /// Point addition, treating (0, 0) as the point at infinity (the precompile's convention).
    function _add(G1Point memory a, G1Point memory b) internal view returns (G1Point memory r) {
        (bool ok, bytes memory output) = address(0x06).staticcall(abi.encodePacked(a.x, a.y, b.x, b.y));
        require(ok, "ecAdd");
        (r.x, r.y) = abi.decode(output, (uint256, uint256));
    }

    function _mul(G1Point memory a, uint256 s) internal view returns (G1Point memory r) {
        (bool ok, bytes memory output) = address(0x07).staticcall(abi.encodePacked(a.x, a.y, s));
        require(ok, "ecMul");
        (r.x, r.y) = abi.decode(output, (uint256, uint256));
    }

    function _neg(G1Point memory a) internal pure returns (G1Point memory) {
        if (a.x == 0 && a.y == 0) {
            return a;
        }
        return G1Point(a.x, Q - (a.y % Q));
    }

    function _pairing(G1Point memory lhs, G1Point memory rhs) internal view returns (bool) {
        bytes memory input = abi.encodePacked(
            lhs.x,
            lhs.y,
            // [1]_2
            uint256(0x198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c2),
            uint256(0x1800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed),
            uint256(0x090689d0585ff075ec9e99ad690c3395bc4b313370b38ef355acdadcd122975b),
            uint256(0x12c85ea5db8c6deb4aab71808dcb408fe3d1e7690c43d37b4ce6cc0166fa7daa),
            rhs.x,
            rhs.y,
            // [tau]_2
            uint256(0x260e01b251f6f1c7e7ff4e580791dee8ea51d87a358e038b4efe30fac09383c1),
            uint256(0x0118c4d5b837bcc2bc89b5b398b5974e9f5944073b32078b7e231fec938883b0),
            uint256(0x04fc6369f7110fe3d25156c1bb9a72859cf2a04641f99ba4ee413c80da6a5fe4),
            uint256(0x22febda3c0c0632a56475b4214e5615e11e6dd3f96e6cea2854a87d4dacc5e55)
        );
        (bool ok, bytes memory output) = address(0x08).staticcall(input);
        return ok && abi.decode(output, (bool));
    }
}
