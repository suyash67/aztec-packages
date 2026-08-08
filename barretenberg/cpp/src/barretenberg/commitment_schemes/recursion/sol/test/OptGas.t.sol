// SPDX-License-Identifier: Apache-2.0
pragma solidity >=0.8.21;

import {Test, console} from "forge-std/Test.sol";
import {HonkVerifier as OptVerifier} from "../src/VerifierOpt.sol";

/// @notice Gas of the `--optimized` UltraHonk verifier - the one Aztec generates for its L1 rollup
/// root proof - on the same 2^19 proof the unoptimized verifier was measured against.
contract OptGasTest is Test {
    OptVerifier public verifier;

    function setUp() public {
        verifier = new OptVerifier();
    }

    function testOptimizedVerifyGas() public view {
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
        bool ok = verifier.verify(proof, publicInputs);
        uint256 executionGas = before - gasleft();
        require(ok, "optimized verification failed");
        console.log("opt_proof_bytes", proof.length);
        console.log("opt_execution_gas", executionGas);
    }
}
