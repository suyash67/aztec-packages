// SPDX-License-Identifier: Apache-2.0
pragma solidity >=0.8.21;

import {Test, console} from "forge-std/Test.sol";
import {HonkVerifierInstrumented} from "./VerifierOptInstrumented.sol";

/// @notice Splits the `--optimized` UltraHonk verifier into the field work every commitment scheme
/// shares and the elliptic-curve work Shplemini+KZG spends on its opening, so the Vela opening can be
/// compared against a like-for-like figure rather than against the unoptimized verifier.
contract OptSplitTest is Test {
    HonkVerifierInstrumented public verifier;

    function setUp() public {
        verifier = new HonkVerifierInstrumented();
    }

    function testSplit() public view {
        string memory dir = vm.envString("PROOF_DIR");
        bytes memory proof = vm.readFileBinary(string.concat(dir, "/proof"));
        bytes memory rawInputs = vm.readFileBinary(string.concat(dir, "/public_inputs"));
        bytes32[] memory publicInputs = new bytes32[](rawInputs.length / 32);
        for (uint256 i = 0; i < publicInputs.length; i++) {
            bytes32 word;
            uint256 offset = 32 + i * 32;
            assembly {
                word := mload(add(rawInputs, offset))
            }
            publicInputs[i] = word;
        }
        uint256 before = gasleft();
        (bool ok, uint256 ecGas) = verifier.verifyWithGas(proof, publicInputs);
        uint256 total = before - gasleft();
        require(ok, "verification failed");
        console.log("opt_total_gas", total);
        console.log("opt_opening_ec_gas", ecGas);
        console.log("opt_shared_gas", total - ecGas);
    }
}
