// SPDX-License-Identifier: Apache-2.0
pragma solidity >=0.8.21;

/**
 * @title The EVM cost of an inner-product-argument verifier's size-n MSM.
 *
 * @notice A Bulletproofs-style IPA needs no trusted setup, but its verifier is O(n) in the circuit
 * size: after folding `log n` rounds it must certify the prover's claimed final generator
 * `G_0 = <s, G>` with one multi-scalar multiplication over all `n` generators
 * (`IpaVerifier::verify_impl` in `commitment_schemes/pedersen_ipa/`). KZG and Vela instead close with
 * a pairing against a constant-size key, so their verifiers touch only the ~40 claim commitments.
 *
 * @dev This measures the MSM the way a Solidity verifier would have to run it — the `ecMul` and
 * `ecAdd` precompiles in a loop — at sizes small enough to execute, so the cost per generator is
 * measured rather than assumed and the figure for a real circuit follows by extrapolation. The
 * per-generator cost is flat by construction (one `ecMul` plus one `ecAdd` each), so the linear model
 * is exact up to loop overhead, which the intercept absorbs.
 *
 * It deliberately does *not* price obtaining the generators. A transparent scheme cannot ship them in
 * a trusted key: they have to be hash-derived on chain, held in contract code (over the 24 KB EIP-170
 * limit past n = 384), or sent as calldata (64 bytes each). All three are themselves O(n) and none is
 * cheaper than the MSM.
 */
contract IpaMsmCost {
    uint256 internal constant P = 21888242871839275222246405745257275088548364400416034343698204186575808495617;

    struct G1Point {
        uint256 x;
        uint256 y;
    }

    /// `sum_i s_i G_i` for `n` generators, with the scalars and points derived the way a verifier's
    /// would be: the scalars from the round challenges, the generators from a seed.
    function msm(uint256 n, uint256 seed) external view returns (uint256 gasUsed, G1Point memory acc) {
        G1Point[] memory generators = new G1Point[](n);
        uint256[] memory scalars = new uint256[](n);
        for (uint256 i = 0; i < n; i++) {
            // Any fixed-cost stand-in for hash-to-curve: the measurement is of the MSM, and the real
            // derivation is strictly more expensive.
            generators[i] = _mul(G1Point(1, 2), uint256(keccak256(abi.encodePacked(seed, i))) % P);
            scalars[i] = uint256(keccak256(abi.encodePacked(i, seed))) % P;
        }

        uint256 start = gasleft();
        for (uint256 i = 0; i < n; i++) {
            acc = _add(acc, _mul(generators[i], scalars[i]));
        }
        gasUsed = start - gasleft();
    }

    /// The same MSM written the way the generated Honk verifier writes its own: one scratch buffer,
    /// `staticcall` straight into it, no per-call memory allocation. This is the floor a production
    /// IPA verifier could reach, so it brackets the naive figure from below.
    function msmAssembly(uint256 n, uint256 seed) external view returns (uint256 gasUsed) {
        G1Point[] memory generators = new G1Point[](n);
        uint256[] memory scalars = new uint256[](n);
        for (uint256 i = 0; i < n; i++) {
            generators[i] = _mul(G1Point(1, 2), uint256(keccak256(abi.encodePacked(seed, i))) % P);
            scalars[i] = uint256(keccak256(abi.encodePacked(i, seed))) % P;
        }

        uint256 start = gasleft();
        bool success = true;
        assembly {
            let free := mload(0x40)
            // [free, free+0x40) accumulates; [free+0x40, free+0xa0) stages each term.
            mstore(free, 0)
            mstore(add(free, 0x20), 0)
            let points := add(generators, 0x20)
            let values := add(scalars, 0x20)
            for { let i := 0 } lt(i, n) { i := add(i, 1) } {
                let point := mload(add(points, mul(i, 0x20)))
                mstore(add(free, 0x40), mload(point))
                mstore(add(free, 0x60), mload(add(point, 0x20)))
                mstore(add(free, 0x80), mload(add(values, mul(i, 0x20))))
                success := and(success, staticcall(gas(), 7, add(free, 0x40), 0x60, add(free, 0x40), 0x40))
                success := and(success, staticcall(gas(), 6, free, 0x80, free, 0x40))
            }
        }
        gasUsed = start - gasleft();
        require(success, "msm");
    }

    /// The `log n` folding rounds: two scalar multiplications and two additions each.
    function folding(uint256 logN, uint256 seed) external view returns (uint256 gasUsed) {
        G1Point memory acc = G1Point(1, 2);
        uint256 start = gasleft();
        for (uint256 i = 0; i < logN; i++) {
            uint256 u = uint256(keccak256(abi.encodePacked(seed, i))) % P;
            uint256 uInv = uint256(keccak256(abi.encodePacked(i, seed))) % P;
            acc = _add(acc, _mul(G1Point(1, 2), u));
            acc = _add(acc, _mul(G1Point(1, 2), uInv));
        }
        gasUsed = start - gasleft();
    }

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
}
