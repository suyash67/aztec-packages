// SPDX-License-Identifier: Apache-2.0
pragma solidity >=0.8.21;

import {Test, console} from "forge-std/Test.sol";
import {HonkVerifier, Honk, Transcript, TranscriptLib, LOG_N, VK_HASH, NUMBER_OF_PUBLIC_INPUTS} from
    "../src/Verifier.sol";
import {VelaOpeningVerifier} from "../src/VelaOpening.sol";
import {VelaOpeningVerifierOpt} from "../src/VelaOpeningOpt.sol";

/// @notice Exposes the generated verifier's opening step so its gas can be separated from the
/// transcript/sumcheck work the two commitment schemes share.
contract ShpleminiHarness is HonkVerifier {
    function shpleminiGas(bytes calldata proof, bytes32[] calldata publicInputs)
        external
        view
        returns (uint256 openingGas, uint256 sharedGas)
    {
        uint256 start = gasleft();
        Honk.VerificationKey memory vk = loadVerificationKey();
        Honk.Proof memory p = TranscriptLib.loadProof(proof, LOG_N);
        Transcript memory t =
            TranscriptLib.generateTranscript(p, publicInputs, VK_HASH, NUMBER_OF_PUBLIC_INPUTS, LOG_N);
        t.relationParameters.publicInputsDelta = computePublicInputDelta(
            publicInputs, p.pairingPointObject, t.relationParameters.beta, t.relationParameters.gamma, 5
        );
        require(verifySumcheck(p, t), "sumcheck");
        sharedGas = start - gasleft();

        start = gasleft();
        bool ok = verifyShplemini(p, vk, t);
        openingGas = start - gasleft();
        require(ok, "shplemini");
    }
}

/// @notice On-chain cost of one UltraHonk verification, split into the shared Honk work and the
/// polynomial-commitment opening, for Shplemini+KZG and for Vela.
contract PcsGasTest is Test {
    ShpleminiHarness public shplemini;
    VelaOpeningVerifier public vela;
    VelaOpeningVerifierOpt public velaOpt;

    function setUp() public {
        shplemini = new ShpleminiHarness();
        vela = new VelaOpeningVerifier();
        velaOpt = new VelaOpeningVerifierOpt();
    }

    function testShpleminiGas() public view {
        string memory dir = vm.envString("PROOF_DIR");
        bytes memory proof = vm.readFileBinary(string.concat(dir, "/proof"));
        bytes memory rawInputs = vm.readFileBinary(string.concat(dir, "/public_inputs"));
        bytes32[] memory publicInputs = new bytes32[](rawInputs.length / 32);
        for (uint256 i = 0; i < publicInputs.length; i++) {
            publicInputs[i] = _word(rawInputs, i);
        }
        (uint256 openingGas, uint256 sharedGas) = shplemini.shpleminiGas(proof, publicInputs);
        console.log("kzg_proof_bytes", proof.length);
        console.log("kzg_shared_gas", sharedGas);
        console.log("kzg_opening_gas", openingGas);
        console.log("kzg_calldata_gas", _calldataGas(proof));
    }

    function testVelaGas() public view {
        string memory dir = vm.envString("PROOF_DIR");
        bytes memory raw = vm.readFileBinary(string.concat(dir, "/vela_opening.bin"));
        VelaOpeningVerifier.Opening memory o;
        o.logN = uint256(_word(raw, 0));
        o.numUnshifted = uint256(_word(raw, 1));
        o.numShifted = uint256(_word(raw, 2));
        uint256 numAll = uint256(_word(raw, 3));
        o.previousChallenge = uint256(_word(raw, 4));

        uint256 cursor = 5;
        o.evaluations = new uint256[](numAll);
        for (uint256 i = 0; i < numAll; i++) {
            o.evaluations[i] = uint256(_word(raw, cursor++));
        }
        uint256 numClaims = o.numUnshifted + o.numShifted;
        o.claimEntity = new uint256[](numClaims);
        for (uint256 i = 0; i < numClaims; i++) {
            o.claimEntity[i] = uint256(_word(raw, cursor++));
        }
        o.u = new uint256[](o.logN);
        for (uint256 i = 0; i < o.logN; i++) {
            o.u[i] = uint256(_word(raw, cursor++));
        }
        o.commitments = new VelaOpeningVerifier.G1Point[](numClaims);
        for (uint256 i = 0; i < numClaims; i++) {
            o.commitments[i].x = uint256(_word(raw, cursor++));
            o.commitments[i].y = uint256(_word(raw, cursor++));
        }
        o.cH = VelaOpeningVerifier.G1Point(uint256(_word(raw, cursor)), uint256(_word(raw, cursor + 1)));
        cursor += 2;
        o.v0a = uint256(_word(raw, cursor++));
        o.v1a = uint256(_word(raw, cursor++));
        o.v0b = uint256(_word(raw, cursor++));
        o.v1b = uint256(_word(raw, cursor++));
        o.w0 = uint256(_word(raw, cursor++));
        o.cQ = VelaOpeningVerifier.G1Point(uint256(_word(raw, cursor)), uint256(_word(raw, cursor + 1)));
        cursor += 2;
        o.piL = VelaOpeningVerifier.G1Point(uint256(_word(raw, cursor)), uint256(_word(raw, cursor + 1)));
        cursor += 2;
        require(cursor * 32 == raw.length, "vela export length");

        uint256 start = gasleft();
        bool ok = vela.verify(o);
        uint256 openingGas = start - gasleft();
        require(ok, "vela opening failed");

        start = gasleft();
        (bool okOpt, uint256 ecGas) = velaOpt.verifyWithGas(o);
        uint256 openingGasOpt = start - gasleft();
        require(okOpt, "optimized vela opening failed");
        console.log("vela_opening_gas_optimized", openingGasOpt);
        console.log("vela_opening_ec_gas_optimized", ecGas);
        bytes memory velaProof = vm.readFileBinary(string.concat(dir, "/vela_opening.bin.proof"));
        console.log("vela_claims", numClaims);
        console.log("vela_proof_bytes", velaProof.length);
        console.log("vela_opening_gas", openingGas);
        console.log("vela_calldata_gas", _calldataGas(velaProof));
    }

    function _word(bytes memory data, uint256 index) internal pure returns (bytes32 word) {
        uint256 offset = 32 + index * 32;
        assembly {
            word := mload(add(data, offset))
        }
    }

    /// EIP-2028 calldata pricing: 4 gas per zero byte, 16 per non-zero byte.
    function _calldataGas(bytes memory data) internal pure returns (uint256 total) {
        for (uint256 i = 0; i < data.length; i++) {
            total += data[i] == 0 ? 4 : 16;
        }
    }
}
