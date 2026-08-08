// SPDX-License-Identifier: Apache-2.0
pragma solidity >=0.8.21;

import "forge-std/Test.sol";
import "../src/UltraFflonkVerifier.sol";
import "./Fixture.sol";

/// @dev A same-signature no-op, so the gas measurement can subtract the external-call overhead.
contract Probe {
    function noop(bytes calldata, uint256[] calldata) external pure returns (bool) {
        return true;
    }
}

/**
 * @notice The on-chain verifier against proofs barretenberg produced and its own verifier accepted.
 *
 * @dev Every fixture here was verified natively before it was written out, so a passing acceptance
 * test is two independent implementations agreeing on one proof, and each rejection test is both of
 * them rejecting the same tampering. What this suite cannot show on its own is that the relation set
 * is right — a wrong relation would reject honest proofs, which looks the same as a wrong proof.
 * `RelationVectors.t.sol` is what covers that.
 */
contract UltraFflonkVerifierTest is Test {
    Probe internal probe;

    function setUp() public {
        probe = new Probe();
    }

    function _deploy(Fixture.Case memory c) internal returns (UltraFflonkVerifier) {
        return new UltraFflonkVerifier(c.circuitSize, c.numPublicInputs, c.pubInputsOffset, c.omega, c.preprocessed);
    }

    // -----------------------------------------------------------------------------------------
    // Acceptance
    // -----------------------------------------------------------------------------------------

    function testAcceptsEveryFixture() public {
        for (uint256 i = 0; i < Fixture.COUNT; ++i) {
            Fixture.Case memory c = Fixture.get(i);
            assertTrue(_deploy(c).verify(c.proof, c.publicInputs), c.label);
        }
    }

    function testProofSizeIsFixed() public pure {
        for (uint256 i = 0; i < Fixture.COUNT; ++i) {
            assertEq(Fixture.get(i).proof.length, 2112);
        }
    }

    // -----------------------------------------------------------------------------------------
    // Rejection
    // -----------------------------------------------------------------------------------------

    function testRejectsATamperedEvaluation() public {
        Fixture.Case memory c = Fixture.get(1); // the every-gate circuit
        UltraFflonkVerifier verifier = _deploy(c);

        for (uint256 k = 0; k < 54; ++k) {
            bytes memory tampered = c.proof;
            uint256 offset = 384 + (k * 32) + 31;
            tampered[offset] = bytes1(uint8(tampered[offset]) ^ 0x01);
            assertFalse(verifier.verify(tampered, c.publicInputs), "evaluation");
        }
    }

    function testRejectsATamperedCommitment() public {
        Fixture.Case memory c = Fixture.get(1);
        UltraFflonkVerifier verifier = _deploy(c);

        for (uint256 point = 0; point < 6; ++point) {
            bytes memory tampered = c.proof;
            uint256 offset = (point * 64) + 31;
            tampered[offset] = bytes1(uint8(tampered[offset]) ^ 0x01);
            assertFalse(verifier.verify(tampered, c.publicInputs), "commitment");
        }
    }

    function testRejectsSwappedCommitments() public {
        Fixture.Case memory c = Fixture.get(1);
        UltraFflonkVerifier verifier = _deploy(c);

        bytes memory swapped = c.proof;
        for (uint256 i = 0; i < 64; ++i) {
            bytes1 held = swapped[i];
            swapped[i] = swapped[64 + i];
            swapped[64 + i] = held;
        }
        assertFalse(verifier.verify(swapped, c.publicInputs));
    }

    function testRejectsWrongPublicInputs() public {
        Fixture.Case memory c = Fixture.get(0);
        UltraFflonkVerifier verifier = _deploy(c);

        uint256[] memory wrong = c.publicInputs;
        wrong[0] = addmod(wrong[0], 1, UltraFflonkTypes.R);
        assertFalse(verifier.verify(c.proof, wrong));

        uint256[] memory tooFew = new uint256[](c.publicInputs.length - 1);
        assertFalse(verifier.verify(c.proof, tooFew));
    }

    function testRejectsAProofForAnotherCircuit() public {
        Fixture.Case memory first = Fixture.get(0);
        Fixture.Case memory second = Fixture.get(2);

        // Same public-input count and offset, different circuit: only the key can tell them apart.
        assertFalse(_deploy(first).verify(second.proof, second.publicInputs));
        assertFalse(_deploy(second).verify(first.proof, first.publicInputs));
    }

    function testRejectsAMalformedProof() public {
        Fixture.Case memory c = Fixture.get(0);
        UltraFflonkVerifier verifier = _deploy(c);

        assertFalse(verifier.verify(new bytes(0), c.publicInputs));
        assertFalse(verifier.verify(new bytes(2112), c.publicInputs));

        bytes memory truncated = new bytes(2111);
        for (uint256 i = 0; i < 2111; ++i) {
            truncated[i] = c.proof[i];
        }
        assertFalse(verifier.verify(truncated, c.publicInputs));
    }

    function testRejectsNonCanonicalEncodings() public {
        Fixture.Case memory c = Fixture.get(0);
        UltraFflonkVerifier verifier = _deploy(c);

        // A scalar at or above the modulus must be rejected rather than silently reduced, or a proof
        // would have many encodings and the transcript would not bind it.
        bytes memory nonCanonical = c.proof;
        for (uint256 i = 0; i < 32; ++i) {
            nonCanonical[384 + i] = bytes1(uint8(0xff));
        }
        assertFalse(verifier.verify(nonCanonical, c.publicInputs));

        // A public input likewise.
        uint256[] memory publics = c.publicInputs;
        publics[0] = UltraFflonkTypes.R;
        assertFalse(verifier.verify(c.proof, publics));
    }

    function testRejectsAnInvalidVerificationKey() public {
        Fixture.Case memory c = Fixture.get(0);

        vm.expectRevert(UltraFflonkVerifier.InvalidVerificationKey.selector);
        new UltraFflonkVerifier(c.circuitSize + 1, c.numPublicInputs, c.pubInputsOffset, c.omega, c.preprocessed);

        // An omega whose order is not the circuit size would put this verifier's vanishing polynomial
        // and shift on a different domain than the prover committed to.
        vm.expectRevert(UltraFflonkVerifier.InvalidVerificationKey.selector);
        new UltraFflonkVerifier(c.circuitSize, c.numPublicInputs, c.pubInputsOffset, 1, c.preprocessed);

        uint256[8] memory offCurve = c.preprocessed;
        offCurve[0] = addmod(offCurve[0], 1, UltraFflonkTypes.R);
        vm.expectRevert(UltraFflonkVerifier.InvalidVerificationKey.selector);
        new UltraFflonkVerifier(c.circuitSize, c.numPublicInputs, c.pubInputsOffset, c.omega, offCurve);
    }

    function testFuzzRejectsAnyByteChange(uint256 index, uint8 mask) public {
        Fixture.Case memory c = Fixture.get(0);
        UltraFflonkVerifier verifier = _deploy(c);

        vm.assume(mask != 0);
        index = index % c.proof.length;

        bytes memory tampered = c.proof;
        tampered[index] = bytes1(uint8(tampered[index]) ^ mask);
        assertFalse(verifier.verify(tampered, c.publicInputs));
    }

    // -----------------------------------------------------------------------------------------
    // Cost
    // -----------------------------------------------------------------------------------------

    function testGas() public {
        console.log("circuit | rows | public inputs | execution | calldata | transaction");
        for (uint256 i = 0; i < Fixture.COUNT; ++i) {
            Fixture.Case memory c = Fixture.get(i);
            UltraFflonkVerifier verifier = _deploy(c);

            // The proof is in memory, not storage: reading 2,112 bytes back out of storage would add
            // ~166k of cold SLOADs that a real transaction, whose proof arrives as calldata, never pays.
            uint256 mark = gasleft();
            probe.noop(c.proof, c.publicInputs);
            uint256 overhead = mark - gasleft();

            uint256 start = gasleft();
            bool ok = verifier.verify(c.proof, c.publicInputs);
            uint256 execution = start - gasleft() - overhead;
            assertTrue(ok, c.label);

            uint256 calldataGas = 0;
            for (uint256 b = 0; b < c.proof.length; ++b) {
                calldataGas += c.proof[b] == 0 ? 4 : 16;
            }
            for (uint256 k = 0; k < c.publicInputs.length; ++k) {
                bytes32 word = bytes32(c.publicInputs[k]);
                for (uint256 b = 0; b < 32; ++b) {
                    calldataGas += word[b] == 0 ? 4 : 16;
                }
            }

            console.log(c.label);
            console.log("  rows          ", c.circuitSize);
            console.log("  public inputs ", c.publicInputs.length);
            console.log("  execution gas ", execution);
            console.log("  calldata gas  ", calldataGas);
            console.log("  transaction   ", execution + calldataGas + 21000);
        }
    }
}
