// SPDX-License-Identifier: Apache-2.0
pragma solidity >=0.8.21;

import {Test, console} from "forge-std/Test.sol";
import {HonkVerifier} from "../src/Verifier.sol";

/// @notice Gas of one on-chain UltraHonk (Shplemini + KZG) verification.
/// @dev The proof and public inputs are the ones `bb prove -t evm-no-zk` wrote for the recursive
/// verifier circuit named by PROOF_DIR; the harness reports both the execution gas the call burns
/// and the full transaction cost a sender would pay, which on a proof this size is calldata-heavy.
contract HonkGasTest is Test {
    HonkVerifier public verifier;

    function setUp() public {
        verifier = new HonkVerifier();
    }

    function testVerifyGas() public view {
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
        require(ok, "verification failed");

        uint256 calldataGas = _calldataGas(proof) + _calldataGasWords(publicInputs);
        console.log("proof_bytes", proof.length);
        console.log("public_inputs", publicInputs.length);
        console.log("execution_gas", executionGas);
        console.log("calldata_gas", calldataGas);
        console.log("total_tx_gas", executionGas + calldataGas + 21000);
    }

    /// EIP-2028 calldata pricing: 4 gas per zero byte, 16 per non-zero byte.
    function _calldataGas(bytes memory data) internal pure returns (uint256 total) {
        for (uint256 i = 0; i < data.length; i++) {
            total += data[i] == 0 ? 4 : 16;
        }
    }

    function _calldataGasWords(bytes32[] memory words) internal pure returns (uint256 total) {
        for (uint256 i = 0; i < words.length; i++) {
            for (uint256 b = 0; b < 32; b++) {
                total += words[i][b] == 0 ? 4 : 16;
            }
        }
    }
}
