// SPDX-License-Identifier: Apache-2.0
pragma solidity >=0.8.21;

import "./UltraFflonkTypes.sol";
import "./UltraRelations.sol";

/**
 * @title fflonk over the arithmetization Noir compiles to: an L1 verifier for UltraHonk circuits.
 *
 * @notice The circuit is UltraHonk's — four wires, plookup, the custom gates, the whole ACIR
 * frontend. Only the proof system on top of it is different: a univariate quotient and a batched
 * fflonk opening rather than sumcheck and Shplemini. That trades a much heavier prover for a verifier
 * whose group work is a fixed eleven scalar multiplications and one pairing, and a proof that is a
 * fixed 2,112 bytes, whatever the circuit contains.
 *
 * @dev The mirror of `barretenberg/ultra_fflonk`; `PROTOCOL.md` there is the specification both
 * implement. The relation set below is a port of `barretenberg/relations/*.hpp`, and every one of its
 * thirty-one subrelations is checked against barretenberg's own evaluation by `RelationVectors.t.sol`
 * — the batched identity folds them into a single field element, inside which one wrong term would
 * otherwise be invisible.
 */
contract UltraFflonkVerifier {
    /// BN254 base field.
    uint256 internal constant Q = 21888242871839275222246405745257275088696311157297823662689037894645226208583;
    /// BN254 scalar field. Spelled out rather than aliased: inline assembly can only reference a
    /// direct number constant, and the constructor pins it to `UltraFflonkTypes.R`.
    uint256 internal constant R = 21888242871839275222246405745257275088548364400416034343698204186575808495617;

    uint256 internal constant PROOF_LENGTH = 2112;
    uint256 internal constant NUM_EVALUATIONS = 54;
    uint256 internal constant NUM_QUOTIENT_CHUNKS = 6;
    uint256 internal constant NUM_PREPROCESSED_GROUPS = 4;
    uint256 internal constant PACK_PREPROCESSED = 7;

    /// @dev `PERMUTATION_ARGUMENT_VALUE_SEPARATOR` from `barretenberg/constants.hpp`.
    uint256 internal constant PERMUTATION_SEPARATOR = 1 << 28;

    // Byte offsets into the proof blob: four committed groups, then W and W'.
    uint256 internal constant OFF_C_WIRES = 0;
    uint256 internal constant OFF_C_MEMORY = 64;
    uint256 internal constant OFF_C_GRAND_PRODUCT = 128;
    uint256 internal constant OFF_C_QUOTIENT = 192;
    uint256 internal constant OFF_W = 256;
    uint256 internal constant OFF_WP = 320;
    uint256 internal constant OFF_EVALUATIONS = 384;

    // Offsets into the flat evaluation vector. Each group contributes its `xi` openings, then its
    // `xi*omega` ones; the order is protocol.
    uint256 internal constant FLAT_WIRES_XI = 28;
    uint256 internal constant FLAT_WIRES_XI_OMEGA = 31;
    uint256 internal constant FLAT_MEMORY_XI = 34;
    uint256 internal constant FLAT_MEMORY_XI_OMEGA = 37;
    uint256 internal constant FLAT_GRAND_PRODUCT_XI = 40;
    uint256 internal constant FLAT_GRAND_PRODUCT_XI_OMEGA = 44;
    uint256 internal constant FLAT_QUOTIENT = 48;

    // The BN254 G2 generator and the SRS's [x]_2, in the (c1, c0) ordering the precompile expects.
    uint256 internal constant G2_GEN_X_C1 = 0x198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c2;
    uint256 internal constant G2_GEN_X_C0 = 0x1800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed;
    uint256 internal constant G2_GEN_Y_C1 = 0x090689d0585ff075ec9e99ad690c3395bc4b313370b38ef355acdadcd122975b;
    uint256 internal constant G2_GEN_Y_C0 = 0x12c85ea5db8c6deb4aab71808dcb408fe3d1e7690c43d37b4ce6cc0166fa7daa;
    uint256 internal constant G2_TAU_X_C1 = 0x260e01b251f6f1c7e7ff4e580791dee8ea51d87a358e038b4efe30fac09383c1;
    uint256 internal constant G2_TAU_X_C0 = 0x0118c4d5b837bcc2bc89b5b398b5974e9f5944073b32078b7e231fec938883b0;
    uint256 internal constant G2_TAU_Y_C1 = 0x04fc6369f7110fe3d25156c1bb9a72859cf2a04641f99ba4ee413c80da6a5fe4;
    uint256 internal constant G2_TAU_Y_C0 = 0x22febda3c0c0632a56475b4214e5615e11e6dd3f96e6cea2854a87d4dacc5e55;

    // The verification key, inlined into the bytecode at deployment.
    uint256 public immutable CIRCUIT_SIZE;
    uint256 public immutable LOG_CIRCUIT_SIZE;
    uint256 public immutable NUM_PUBLIC_INPUTS;
    uint256 public immutable PUB_INPUTS_OFFSET;
    uint256 public immutable OMEGA;
    uint256 public immutable VK_HASH;
    uint256 internal immutable P0X;
    uint256 internal immutable P0Y;
    uint256 internal immutable P1X;
    uint256 internal immutable P1Y;
    uint256 internal immutable P2X;
    uint256 internal immutable P2Y;
    uint256 internal immutable P3X;
    uint256 internal immutable P3Y;

    /// @dev Challenges and the values derived from them, kept off the stack.
    struct Context {
        uint256 alpha;
        uint256 xi;
        uint256 nu;
        uint256 y;
        uint256 xiPowN;
        uint256 vanishing;
        uint256 xiOmega;
        uint256 totalVanishing;
    }

    error InvalidVerificationKey();

    /**
     * @param circuitSize the number of rows, a power of two
     * @param numPublicInputs how many public inputs the circuit exposes
     * @param pubInputsOffset the row the public inputs start at, which the permutation delta reads
     * @param omega a primitive `circuitSize`-th root of unity
     * @param preprocessed the four preprocessed group commitments, as (x, y) pairs
     */
    constructor(
        uint256 circuitSize,
        uint256 numPublicInputs,
        uint256 pubInputsOffset,
        uint256 omega,
        uint256[8] memory preprocessed
    ) {
        // The relations are written against the library's modulus and this contract's assembly against
        // its own literal; they have to be the same field.
        assert(R == UltraFflonkTypes.R);

        if (circuitSize < 8 || (circuitSize & (circuitSize - 1)) != 0) revert InvalidVerificationKey();
        if (numPublicInputs >= circuitSize) revert InvalidVerificationKey();
        if (pubInputsOffset >= circuitSize) revert InvalidVerificationKey();
        if (omega == 0 || omega >= R) revert InvalidVerificationKey();
        for (uint256 i = 0; i < NUM_PREPROCESSED_GROUPS; ++i) {
            if (!_onCurve(preprocessed[2 * i], preprocessed[2 * i + 1])) revert InvalidVerificationKey();
        }

        uint256 logN = 0;
        while ((uint256(1) << logN) != circuitSize) {
            logN++;
        }

        // omega must have order exactly circuitSize, or the domain this verifier's vanishing
        // polynomial and shift are taken over is not the one the prover committed to.
        {
            uint256 half = omega;
            for (uint256 i = 0; i + 1 < logN; ++i) {
                half = mulmod(half, half, R);
            }
            if (half == 1 || mulmod(half, half, R) != 1) revert InvalidVerificationKey();
        }

        CIRCUIT_SIZE = circuitSize;
        LOG_CIRCUIT_SIZE = logN;
        NUM_PUBLIC_INPUTS = numPublicInputs;
        PUB_INPUTS_OFFSET = pubInputsOffset;
        OMEGA = omega;
        P0X = preprocessed[0];
        P0Y = preprocessed[1];
        P1X = preprocessed[2];
        P1Y = preprocessed[3];
        P2X = preprocessed[4];
        P2Y = preprocessed[5];
        P3X = preprocessed[6];
        P3Y = preprocessed[7];
        VK_HASH = uint256(
            keccak256(
                abi.encodePacked(
                    circuitSize,
                    numPublicInputs,
                    pubInputsOffset,
                    preprocessed[0],
                    preprocessed[1],
                    preprocessed[2],
                    preprocessed[3],
                    preprocessed[4],
                    preprocessed[5],
                    preprocessed[6],
                    preprocessed[7]
                )
            )
        ) % R;
    }

    /**
     * @notice Verify a proof of this circuit for these public inputs.
     * @dev Never reverts: a malformed proof, a bad point or a failed check all return false.
     */
    function verify(bytes calldata proof, uint256[] calldata publicInputs) external view returns (bool) {
        if (proof.length != PROOF_LENGTH) return false;
        if (publicInputs.length != NUM_PUBLIC_INPUTS) return false;

        uint256[NUM_EVALUATIONS] memory ev;
        uint256[12] memory points; // C_wires, C_memory, C_grandProduct, C_quotient, W, W' as (x, y)
        if (!_parse(proof, publicInputs, ev, points)) return false;

        Context memory ctx;
        UltraFflonkTypes.Params memory params;
        if (!_challenges(proof, publicInputs, points, ctx, params)) return false;

        // Everything the protocol inverts, in one vector: the five distinct Shplonk denominators, the
        // gap between the two opening points, and the permutation delta's denominator.
        uint256[7] memory inverses;
        if (!_denominators(publicInputs, ctx, params, inverses)) return false;

        if (!_checkQuotient(ev, ctx, params)) return false;
        return _checkOpening(ev, points, inverses, ctx);
    }

    /**
     * @dev Fiat-Shamir: `state <- keccak256(state || absorbed words) mod R`, from `state = 0`, one
     * group per round of `ultra_honk/oink_prover.cpp` so the challenge schedule is Honk's.
     */
    function _challenges(
        bytes calldata proof,
        uint256[] calldata publicInputs,
        uint256[12] memory points,
        Context memory ctx,
        UltraFflonkTypes.Params memory params
    ) internal view returns (bool) {
        params.eta = _challengeWires(proof, publicInputs);
        params.etaTwo = mulmod(params.eta, params.eta, R);
        params.etaThree = mulmod(params.etaTwo, params.eta, R);
        params.romLogupGamma = uint256(keccak256(abi.encodePacked(params.eta))) % R;

        params.beta = uint256(keccak256(abi.encodePacked(params.romLogupGamma, points[2], points[3]))) % R;
        params.betaSqr = mulmod(params.beta, params.beta, R);
        params.betaCube = mulmod(params.betaSqr, params.beta, R);
        params.gamma = uint256(keccak256(abi.encodePacked(params.beta))) % R;

        ctx.alpha = uint256(keccak256(abi.encodePacked(params.gamma, points[4], points[5]))) % R;
        ctx.xi = uint256(keccak256(abi.encodePacked(ctx.alpha, points[6], points[7]))) % R;
        ctx.nu = _challengeEvaluations(proof, ctx.xi);
        ctx.y = uint256(keccak256(abi.encodePacked(ctx.nu, points[8], points[9]))) % R;

        uint256 xiPowN = ctx.xi;
        for (uint256 i = 0; i < LOG_CIRCUIT_SIZE; ++i) {
            xiPowN = mulmod(xiPowN, xiPowN, R);
        }
        // xi must miss the domain, or the vanishing polynomial is zero and every quotient claim is free.
        if (ctx.xi == 0 || xiPowN == 1) return false;

        ctx.xiPowN = xiPowN;
        ctx.vanishing = addmod(xiPowN, R - 1, R);
        ctx.xiOmega = mulmod(ctx.xi, OMEGA, R);
        return true;
    }

    /**
     * @dev The Shplonk denominators, the opening-point gap, and the permutation delta's denominator,
     * inverted together by Montgomery's trick. Writes the delta into `params` as a side effect,
     * because it is the only consumer of the last inverse.
     *
     * The four preprocessed groups have the same shape, so they share one denominator; only five
     * distinct `Z_g(y)` exist.
     */
    function _denominators(
        uint256[] calldata publicInputs,
        Context memory ctx,
        UltraFflonkTypes.Params memory params,
        uint256[7] memory inverses
    ) internal view returns (bool) {
        {
            uint256 y2 = mulmod(ctx.y, ctx.y, R);
            uint256 y3 = mulmod(y2, ctx.y, R);
            uint256 y4 = mulmod(y2, y2, R);
            uint256 y6 = mulmod(y3, y3, R);
            uint256 y7 = mulmod(y6, ctx.y, R);

            inverses[0] = addmod(y7, R - ctx.xi, R); // preprocessed: y^7 - xi
            inverses[1] = mulmod(addmod(y3, R - ctx.xi, R), addmod(y3, R - ctx.xiOmega, R), R); // wires
            inverses[2] = inverses[1]; // memory: same pack and point set
            inverses[3] = mulmod(addmod(y4, R - ctx.xi, R), addmod(y4, R - ctx.xiOmega, R), R); // grand product
            inverses[4] = addmod(y6, R - ctx.xi, R); // quotient: y^6 - xi
        }
        inverses[5] = addmod(ctx.xiOmega, R - ctx.xi, R);

        // The public inputs reach the relations only through this correction term, which is the
        // verifier's whole binding of them to the witness.
        uint256 numerator = 1;
        {
            uint256 offset = PUB_INPUTS_OFFSET;
            uint256 numeratorAcc =
                addmod(params.gamma, mulmod(params.beta, addmod(PERMUTATION_SEPARATOR, offset, R), R), R);
            uint256 denominatorAcc = addmod(params.gamma, R - mulmod(params.beta, addmod(offset, 1, R), R), R);
            uint256 denominator = 1;
            uint256 length = publicInputs.length;
            for (uint256 i = 0; i < length; ++i) {
                numerator = mulmod(numerator, addmod(numeratorAcc, publicInputs[i], R), R);
                denominator = mulmod(denominator, addmod(denominatorAcc, publicInputs[i], R), R);
                if (i + 1 < length) {
                    numeratorAcc = addmod(numeratorAcc, params.beta, R);
                    denominatorAcc = addmod(denominatorAcc, R - params.beta, R);
                }
            }
            inverses[6] = denominator;
        }

        if (!_batchInvert(inverses)) return false;
        params.publicInputDelta = mulmod(numerator, inverses[6], R);
        return true;
    }

    /**
     * @dev The batched quotient identity, entirely in the field over the claimed evaluations.
     *
     * `sum_i alpha^i R_i + sum_j alpha^{31+j} (S_j(xi omega) - S_j(xi) - f_j) = T(xi) Z_H(xi)`.
     *
     * The `alpha` power advances across *every* subrelation, including the two Sumcheck enforces as a
     * trace sum rather than per row, so a subrelation's power depends only on the flavor's relation
     * order. Those two are held back and enter through their running sums instead.
     */
    function _checkQuotient(
        uint256[NUM_EVALUATIONS] memory ev,
        Context memory ctx,
        UltraFflonkTypes.Params memory params
    ) internal pure returns (bool) {
        uint256 numerator;
        {
            uint256[31] memory sub;
            {
                uint256[41] memory ent;
                _entities(ev, ent);
                UltraRelations.evaluate(ent, params, sub);
            }

            uint256 alphaPower = 1;
            for (uint256 i = 0; i < UltraRelations.NUM_SUBRELATIONS; ++i) {
                if (i != UltraRelations.SUMMED_LOOKUP && i != UltraRelations.SUMMED_ROM_LOGUP) {
                    numerator = addmod(numerator, mulmod(alphaPower, sub[i], R), R);
                }
                alphaPower = mulmod(alphaPower, ctx.alpha, R);
            }

            // The two running sums, in relation order: the lookup identity then the ROM-LogUp sum.
            uint256[2] memory summed = [sub[UltraRelations.SUMMED_LOOKUP], sub[UltraRelations.SUMMED_ROM_LOGUP]];
            for (uint256 j = 0; j < 2; ++j) {
                uint256 shifted = ev[FLAT_GRAND_PRODUCT_XI_OMEGA + 2 + j];
                uint256 current = ev[FLAT_GRAND_PRODUCT_XI + 2 + j];
                uint256 identity = addmod(addmod(shifted, R - current, R), R - summed[j], R);
                numerator = addmod(numerator, mulmod(alphaPower, identity, R), R);
                alphaPower = mulmod(alphaPower, ctx.alpha, R);
            }
        }

        uint256 quotientAtXi;
        for (uint256 c = NUM_QUOTIENT_CHUNKS; c > 0; --c) {
            quotientAtXi = addmod(mulmod(quotientAtXi, ctx.xiPowN, R), ev[FLAT_QUOTIENT + c - 1], R);
        }
        return numerator == mulmod(quotientAtXi, ctx.vanishing, R);
    }

    /**
     * @dev Reassemble the flavor's per-row container from the proof's evaluations.
     *
     * The five shifted entities come from the `xi*omega` half of the two-point groups, which is the
     * whole reason those groups are opened twice. Three evaluations are not read here — the lookup
     * read counts, read tags and inverses at `xi*omega` — because no relation reads them shifted;
     * they are revealed only because they share a group with a column that is.
     */
    function _entities(uint256[NUM_EVALUATIONS] memory ev, uint256[41] memory ent) internal pure {
        for (uint256 i = 0; i < 28; ++i) {
            ent[i] = ev[i];
        }
        ent[UltraFflonkTypes.E_W_L] = ev[FLAT_WIRES_XI];
        ent[UltraFflonkTypes.E_W_R] = ev[FLAT_WIRES_XI + 1];
        ent[UltraFflonkTypes.E_W_O] = ev[FLAT_WIRES_XI + 2];
        ent[UltraFflonkTypes.E_W_L_SHIFT] = ev[FLAT_WIRES_XI_OMEGA];
        ent[UltraFflonkTypes.E_W_R_SHIFT] = ev[FLAT_WIRES_XI_OMEGA + 1];
        ent[UltraFflonkTypes.E_W_O_SHIFT] = ev[FLAT_WIRES_XI_OMEGA + 2];

        ent[UltraFflonkTypes.E_LOOKUP_READ_COUNTS] = ev[FLAT_MEMORY_XI];
        ent[UltraFflonkTypes.E_LOOKUP_READ_TAGS] = ev[FLAT_MEMORY_XI + 1];
        ent[UltraFflonkTypes.E_W_4] = ev[FLAT_MEMORY_XI + 2];
        ent[UltraFflonkTypes.E_W_4_SHIFT] = ev[FLAT_MEMORY_XI_OMEGA + 2];

        ent[UltraFflonkTypes.E_LOOKUP_INVERSES] = ev[FLAT_GRAND_PRODUCT_XI];
        ent[UltraFflonkTypes.E_Z_PERM] = ev[FLAT_GRAND_PRODUCT_XI + 1];
        ent[UltraFflonkTypes.E_Z_PERM_SHIFT] = ev[FLAT_GRAND_PRODUCT_XI_OMEGA + 1];
    }

    /**
     * @dev `F = sum_g s_g C_g - (sum_g s_g R_g(y)) [1] - Z_T(y) W + y W'`, then
     * `e(F, [1]_2) = e(W', [x]_2)`, with `s_g = nu^g Z_T(y)/Z_g(y)`.
     *
     * Eleven scalar multiplications and one pairing, whatever the circuit contains.
     */
    function _checkOpening(
        uint256[NUM_EVALUATIONS] memory ev,
        uint256[12] memory points,
        uint256[7] memory inverses,
        Context memory ctx
    ) internal view returns (bool) {
        // Z_T(y) = prod_g Z_g(y); the preprocessed denominator appears four times.
        {
            uint256 inverseTotal = mulmod(inverses[0], inverses[0], R);
            inverseTotal = mulmod(inverseTotal, inverseTotal, R);
            inverseTotal = mulmod(inverseTotal, mulmod(inverses[1], inverses[2], R), R);
            inverseTotal = mulmod(inverseTotal, mulmod(inverses[3], inverses[4], R), R);
            ctx.totalVanishing = _inverse(inverseTotal);
        }

        uint256[3] memory accumulator;
        uint256 constantTerm;
        {
            uint256 nuPower = 1;
            for (uint256 g = 0; g < 8; ++g) {
                uint256 denominatorIndex = g < NUM_PREPROCESSED_GROUPS ? 0 : g - NUM_PREPROCESSED_GROUPS + 1;
                uint256 scalar = mulmod(mulmod(nuPower, ctx.totalVanishing, R), inverses[denominatorIndex], R);
                constantTerm = addmod(constantTerm, mulmod(scalar, _residueAt(ev, g, inverses[5], ctx), R), R);
                if (!_accumulateGroup(accumulator, points, g, scalar)) return false;
                nuPower = mulmod(nuPower, ctx.nu, R);
            }
        }

        if (!_accumulate(accumulator, 1, 2, R - constantTerm)) return false;
        if (!_accumulate(accumulator, points[8], points[9], R - ctx.totalVanishing)) return false;
        if (!_accumulate(accumulator, points[10], points[11], ctx.y)) return false;

        // -W'. Negating zero has to stay zero, or Q would not be a canonical coordinate and the
        // precompile would reject a proof the native verifier accepts.
        uint256 negatedY = points[11] == 0 ? 0 : Q - points[11];
        return _pairing(accumulator[0], accumulator[1], points[10], negatedY);
    }

    /**
     * @dev `R_g(y)`, the residue of the group's packed polynomial modulo `Z_g`, evaluated at `y`.
     *
     * For a one-point group the residue's coefficients *are* the claimed evaluations, so this is a
     * Horner evaluation over them. For a two-point group it is the Chinese-remainder interpolant of
     * the two residues, `[B(y)(y^t - xi) - A(y)(y^t - xi omega)] / (xi omega - xi)`.
     */
    function _residueAt(uint256[NUM_EVALUATIONS] memory ev, uint256 group, uint256 inverseGap, Context memory ctx)
        internal
        pure
        returns (uint256)
    {
        if (group < NUM_PREPROCESSED_GROUPS) {
            return _horner(ev, group * PACK_PREPROCESSED, PACK_PREPROCESSED, ctx.y);
        }
        if (group == 7) {
            return _horner(ev, FLAT_QUOTIENT, NUM_QUOTIENT_CHUNKS, ctx.y);
        }

        (uint256 firstXi, uint256 firstXiOmega, uint256 pack) = group == 4
            ? (FLAT_WIRES_XI, FLAT_WIRES_XI_OMEGA, uint256(3))
            : (group == 5
                    ? (FLAT_MEMORY_XI, FLAT_MEMORY_XI_OMEGA, uint256(3))
                    : (FLAT_GRAND_PRODUCT_XI, FLAT_GRAND_PRODUCT_XI_OMEGA, uint256(4)));

        uint256 a = _horner(ev, firstXi, pack, ctx.y);
        uint256 b = _horner(ev, firstXiOmega, pack, ctx.y);
        uint256 yPowT = _power(ctx.y, pack);
        uint256 left = mulmod(b, addmod(yPowT, R - ctx.xi, R), R);
        uint256 right = mulmod(a, addmod(yPowT, R - ctx.xiOmega, R), R);
        return mulmod(addmod(left, R - right, R), inverseGap, R);
    }

    function _horner(uint256[NUM_EVALUATIONS] memory ev, uint256 first, uint256 count, uint256 y)
        internal
        pure
        returns (uint256 result)
    {
        for (uint256 i = count; i > 0; --i) {
            result = addmod(mulmod(result, y, R), ev[first + i - 1], R);
        }
    }

    function _power(uint256 base, uint256 exponent) internal pure returns (uint256 result) {
        result = 1;
        for (uint256 i = 0; i < exponent; ++i) {
            result = mulmod(result, base, R);
        }
    }

    /// @dev The preprocessed groups' commitments live in the verification key, not the proof.
    function _accumulateGroup(uint256[3] memory accumulator, uint256[12] memory points, uint256 group, uint256 scalar)
        internal
        view
        returns (bool)
    {
        if (group == 0) return _accumulate(accumulator, P0X, P0Y, scalar);
        if (group == 1) return _accumulate(accumulator, P1X, P1Y, scalar);
        if (group == 2) return _accumulate(accumulator, P2X, P2Y, scalar);
        if (group == 3) return _accumulate(accumulator, P3X, P3Y, scalar);
        uint256 base = 2 * (group - NUM_PREPROCESSED_GROUPS);
        return _accumulate(accumulator, points[base], points[base + 1], scalar);
    }

    /**
     * @dev Read the proof, rejecting non-canonical scalars and points that are not on the curve.
     *
     * Read with `calldataload` rather than by slicing: at sixty-six words the per-slice bounds
     * checking would be the largest non-cryptographic cost in the verifier. The length was already
     * checked by the caller, which is what makes the unchecked reads safe.
     */
    function _parse(
        bytes calldata proof,
        uint256[] calldata publicInputs,
        uint256[NUM_EVALUATIONS] memory ev,
        uint256[12] memory points
    ) internal pure returns (bool) {
        for (uint256 i = 0; i < publicInputs.length; ++i) {
            if (publicInputs[i] >= R) return false;
        }

        uint256 canonical;
        assembly {
            let cursor := proof.offset
            for { let i := 0 } lt(i, 12) { i := add(i, 1) } {
                mstore(add(points, mul(i, 0x20)), calldataload(add(cursor, mul(i, 0x20))))
            }
            cursor := add(cursor, OFF_EVALUATIONS)
            canonical := 1
            for { let i := 0 } lt(i, NUM_EVALUATIONS) { i := add(i, 1) } {
                let value := calldataload(add(cursor, mul(i, 0x20)))
                // A scalar at or above the modulus must be rejected rather than reduced, or a proof
                // would have many encodings and the transcript would not bind it.
                if iszero(lt(value, R)) { canonical := 0 }
                mstore(add(ev, mul(i, 0x20)), value)
            }
        }
        if (canonical != 1) return false;

        for (uint256 i = 0; i < 6; ++i) {
            // The point at infinity is not accepted. An honest prover never produces one - it would
            // mean a committed polynomial was identically zero - and rejecting keeps this verifier and
            // the native one identical on every input.
            if (!_onCurve(points[2 * i], points[2 * i + 1])) return false;
        }
        return true;
    }

    /// @dev keccak256(0 || vkHash || publicInputs || C_wires), the first challenge.
    function _challengeWires(bytes calldata proof, uint256[] calldata publicInputs)
        internal
        view
        returns (uint256 eta)
    {
        uint256 vkHash = VK_HASH;
        assembly {
            let length := mul(publicInputs.length, 0x20)
            let ptr := mload(0x40)
            mstore(ptr, 0)
            mstore(add(ptr, 0x20), vkHash)
            calldatacopy(add(ptr, 0x40), publicInputs.offset, length)
            calldatacopy(add(add(ptr, 0x40), length), add(proof.offset, OFF_C_WIRES), 0x40)
            mstore(0x40, add(add(ptr, 0x80), length))
            eta := mod(keccak256(ptr, add(0x80, length)), R)
        }
    }

    /// @dev keccak256(xi || the fifty-four evaluations), copied straight out of calldata.
    function _challengeEvaluations(bytes calldata proof, uint256 xi) internal pure returns (uint256 nu) {
        assembly {
            let ptr := mload(0x40)
            mstore(ptr, xi)
            calldatacopy(add(ptr, 0x20), add(proof.offset, OFF_EVALUATIONS), mul(NUM_EVALUATIONS, 0x20))
            mstore(0x40, add(ptr, add(0x20, mul(NUM_EVALUATIONS, 0x20))))
            nu := mod(keccak256(ptr, add(0x20, mul(NUM_EVALUATIONS, 0x20))), R)
        }
    }

    /**
     * @dev `accumulator += (px, py) * scalar`, laid out so that ecMul writes its result exactly where
     * ecAdd expects its second operand. Two staticcalls, no copying.
     */
    function _accumulate(uint256[3] memory accumulator, uint256 px, uint256 py, uint256 scalar)
        internal
        view
        returns (bool ok)
    {
        assembly {
            let ptr := mload(0x40)
            mstore(ptr, mload(accumulator))
            mstore(add(ptr, 0x20), mload(add(accumulator, 0x20)))
            mstore(add(ptr, 0x40), px)
            mstore(add(ptr, 0x60), py)
            mstore(add(ptr, 0x80), scalar)

            ok := staticcall(gas(), 0x07, add(ptr, 0x40), 0x60, add(ptr, 0x40), 0x40)
            ok := and(ok, staticcall(gas(), 0x06, ptr, 0x80, ptr, 0x40))

            mstore(accumulator, mload(ptr))
            mstore(add(accumulator, 0x20), mload(add(ptr, 0x20)))
        }
    }

    /// @dev e(a, [1]_2) * e(b, [x]_2) == 1: two pairs, 384 bytes, written straight into scratch.
    function _pairing(uint256 ax, uint256 ay, uint256 bx, uint256 by) internal view returns (bool ok) {
        assembly {
            let ptr := mload(0x40)
            mstore(ptr, ax)
            mstore(add(ptr, 0x20), ay)
            mstore(add(ptr, 0x40), G2_GEN_X_C1)
            mstore(add(ptr, 0x60), G2_GEN_X_C0)
            mstore(add(ptr, 0x80), G2_GEN_Y_C1)
            mstore(add(ptr, 0xa0), G2_GEN_Y_C0)
            mstore(add(ptr, 0xc0), bx)
            mstore(add(ptr, 0xe0), by)
            mstore(add(ptr, 0x100), G2_TAU_X_C1)
            mstore(add(ptr, 0x120), G2_TAU_X_C0)
            mstore(add(ptr, 0x140), G2_TAU_Y_C1)
            mstore(add(ptr, 0x160), G2_TAU_Y_C0)
            mstore(0x40, add(ptr, 0x180))
            ok := staticcall(gas(), 0x08, ptr, 0x180, ptr, 0x20)
            ok := and(ok, eq(mload(ptr), 1))
        }
    }

    /// @dev Montgomery's trick: one modular exponentiation for the whole vector.
    function _batchInvert(uint256[7] memory values) internal view returns (bool) {
        uint256[7] memory prefix;
        uint256 running = 1;
        for (uint256 i = 0; i < 7; ++i) {
            if (values[i] == 0) return false;
            running = mulmod(running, values[i], R);
            prefix[i] = running;
        }

        uint256 inverse = _inverse(running);
        for (uint256 i = 7; i > 1; --i) {
            uint256 value = values[i - 1];
            values[i - 1] = mulmod(inverse, prefix[i - 2], R);
            inverse = mulmod(inverse, value, R);
        }
        values[0] = inverse;
        return true;
    }

    /// @dev `value^(R-2) mod R` via the modexp precompile.
    function _inverse(uint256 value) internal view returns (uint256 result) {
        assembly {
            let ptr := mload(0x40)
            mstore(ptr, 0x20)
            mstore(add(ptr, 0x20), 0x20)
            mstore(add(ptr, 0x40), 0x20)
            mstore(add(ptr, 0x60), value)
            mstore(add(ptr, 0x80), sub(R, 2))
            mstore(add(ptr, 0xa0), R)
            if iszero(staticcall(gas(), 0x05, ptr, 0xc0, ptr, 0x20)) { revert(0, 0) }
            result := mload(ptr)
        }
    }

    /// @dev `y^2 == x^3 + 3` over the base field, with canonical coordinates. Rejects infinity.
    function _onCurve(uint256 x, uint256 y) internal pure returns (bool) {
        if (x >= Q || y >= Q) return false;
        if (x == 0 && y == 0) return false;
        return mulmod(y, y, Q) == addmod(mulmod(mulmod(x, x, Q), x, Q), 3, Q);
    }
}
