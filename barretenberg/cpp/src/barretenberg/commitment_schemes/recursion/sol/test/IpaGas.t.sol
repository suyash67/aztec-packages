// SPDX-License-Identifier: Apache-2.0
pragma solidity >=0.8.21;

import {Test, console} from "forge-std/Test.sol";
import {IpaMsmCost} from "../src/IpaMsmCost.sol";

/// @notice What an inner-product argument would cost to verify on Ethereum, measured per generator.
contract IpaGasTest is Test {
    IpaMsmCost public ipa;

    function setUp() public {
        ipa = new IpaMsmCost();
    }

    function testMsmCostPerGenerator() public view {
        uint256[4] memory sizes = [uint256(256), 512, 1024, 2048];
        uint256[4] memory used;
        for (uint256 i = 0; i < sizes.length; i++) {
            (uint256 gasUsed,) = ipa.msm(sizes[i], 7);
            used[i] = gasUsed;
            console.log("msm_n", sizes[i]);
            console.log("msm_gas", gasUsed);
        }
        // Slope between the two largest points: the marginal cost of one more generator.
        uint256 slope = (used[3] - used[2]) / (sizes[3] - sizes[2]);
        console.log("gas_per_generator_naive", slope);
        console.log("projected_2_pow_19_naive", slope * (1 << 19));

        uint256 asm1024 = ipa.msmAssembly(1024, 7);
        uint256 asm2048 = ipa.msmAssembly(2048, 7);
        console.log("msm_asm_gas_1024", asm1024);
        console.log("msm_asm_gas_2048", asm2048);
        uint256 asmSlope = (asm2048 - asm1024) / 1024;
        console.log("gas_per_generator_asm", asmSlope);
        console.log("projected_2_pow_19_asm", asmSlope * (1 << 19));
    }

    function testFoldingCost() public view {
        console.log("folding_gas_log_n_19", ipa.folding(19, 7));
    }
}
