// SPDX-License-Identifier: Apache-2.0
pragma solidity >=0.8.21;

import {Test, console} from "forge-std/Test.sol";
import {FflonkOpeningVerifier} from "../src/FflonkOpening.sol";

/// @notice On-chain cost of the fflonk opening, against the Shplemini and Vela figures for the same
/// claim set. The packing collapses the batching MSM to one scalar multiplication per commitment
/// round, so the group work no longer scales with the column count.
contract FflonkGasTest is Test {
    FflonkOpeningVerifier public fflonk;

    function setUp() public {
        fflonk = new FflonkOpeningVerifier();
    }

    function testFflonkGas() public view {
        string memory dir = vm.envString("PROOF_DIR");
        bytes memory raw = vm.readFileBinary(string.concat(dir, "/fflonk_opening.bin"));
        FflonkOpeningVerifier.Opening memory o;

        uint256 cursor = 0;
        o.logN = _w(raw, cursor++);
        uint256 numGroups = _w(raw, cursor++);
        uint256 numAll = _w(raw, cursor++);
        o.numUnshifted = _w(raw, cursor++);
        o.numShifted = _w(raw, cursor++);
        o.previousChallenge = _w(raw, cursor++);

        o.packs = new uint256[](numGroups + 1);
        for (uint256 i = 0; i < numGroups + 1; i++) {
            o.packs[i] = _w(raw, cursor++);
        }
        o.evaluations = new uint256[](numAll);
        for (uint256 i = 0; i < numAll; i++) {
            o.evaluations[i] = _w(raw, cursor++);
        }
        uint256 numClaims = o.numUnshifted + o.numShifted;
        o.claimGroup = new uint256[](numClaims);
        for (uint256 i = 0; i < numClaims; i++) {
            o.claimGroup[i] = _w(raw, cursor++);
        }
        o.claimColumn = new uint256[](numClaims);
        for (uint256 i = 0; i < numClaims; i++) {
            o.claimColumn[i] = _w(raw, cursor++);
        }
        o.claimEntity = new uint256[](numClaims);
        for (uint256 i = 0; i < numClaims; i++) {
            o.claimEntity[i] = _w(raw, cursor++);
        }
        o.u = new uint256[](o.logN);
        for (uint256 i = 0; i < o.logN; i++) {
            o.u[i] = _w(raw, cursor++);
        }
        // The round commitments, then h's: `packs` and `commitments` run in step.
        o.commitments = new FflonkOpeningVerifier.G1Point[](numGroups + 1);
        for (uint256 i = 0; i < numGroups + 1; i++) {
            o.commitments[i].x = _w(raw, cursor++);
            o.commitments[i].y = _w(raw, cursor++);
        }
        uint256 totalPack = 0;
        for (uint256 i = 0; i < o.packs.length; i++) {
            totalPack += o.packs[i];
        }
        o.pairs = new uint256[](totalPack * 2);
        for (uint256 i = 0; i < totalPack * 2; i++) {
            o.pairs[i] = _w(raw, cursor++);
        }
        o.cW = FflonkOpeningVerifier.G1Point(_w(raw, cursor), _w(raw, cursor + 1));
        cursor += 2;
        o.cWp = FflonkOpeningVerifier.G1Point(_w(raw, cursor), _w(raw, cursor + 1));
        cursor += 2;
        require(cursor * 32 == raw.length, "fflonk export length");

        uint256 start = gasleft();
        (bool ok, uint256 ecGas) = fflonk.verifyWithGas(o);
        uint256 openingGas = start - gasleft();
        require(ok, "fflonk opening failed");

        bytes memory proof = vm.readFileBinary(string.concat(dir, "/fflonk_opening.bin.proof"));
        console.log("fflonk_packs_total", totalPack);
        console.log("fflonk_proof_bytes", proof.length);
        console.log("fflonk_opening_gas", openingGas);
        console.log("fflonk_opening_ec_gas", ecGas);
        console.log("fflonk_calldata_gas", _calldataGas(proof));
    }

    function _w(bytes memory data, uint256 index) internal pure returns (uint256 word) {
        uint256 offset = 32 + index * 32;
        assembly {
            word := mload(add(data, offset))
        }
    }

    function _calldataGas(bytes memory data) internal pure returns (uint256 total) {
        for (uint256 i = 0; i < data.length; i++) {
            total += data[i] == 0 ? 4 : 16;
        }
    }
}
