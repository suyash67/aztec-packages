// SPDX-License-Identifier: Apache-2.0
pragma solidity >=0.8.21;

import "forge-std/Test.sol";
import "../src/UltraRelations.sol";
import "./RelationVectors.sol";

/**
 * @notice The Solidity relation set against barretenberg's own, subrelation by subrelation.
 *
 * @dev The verifier folds all thirty-one subrelations into one field element with powers of `alpha`,
 * and a proof either verifies or does not — inside which a single mistranslated term is invisible,
 * and indistinguishable from a hundred of them. This test is the only thing in the suite that can say
 * *which* relation is wrong, so it is the one that makes the port reviewable at all.
 *
 * The rows are uniformly random rather than drawn from a real trace. That is deliberate: this is a
 * differential test between two evaluations of the same polynomial, not a satisfiability test, so
 * random inputs exercise every term and no term can be masked by a selector that happens to be zero.
 */
contract RelationVectorsTest is Test {
    function _params(uint256[9] memory raw) internal pure returns (UltraFflonkTypes.Params memory p) {
        p.eta = raw[0];
        p.etaTwo = raw[1];
        p.etaThree = raw[2];
        p.romLogupGamma = raw[3];
        p.beta = raw[4];
        p.gamma = raw[5];
        p.publicInputDelta = raw[6];
        p.betaSqr = raw[7];
        p.betaCube = raw[8];
    }

    function testEverySubrelationMatchesBarretenberg() public {
        for (uint256 i = 0; i < RelationVectors.COUNT; ++i) {
            RelationVectors.Case memory c = RelationVectors.get(i);

            uint256[31] memory got;
            UltraRelations.evaluate(c.e, _params(c.params), got);

            for (uint256 k = 0; k < RelationVectors.NUM_SUBRELATIONS; ++k) {
                assertEq(
                    got[k], c.expected[k], string.concat("vector ", vm.toString(i), ", subrelation ", vm.toString(k))
                );
            }
        }
    }

    /// @notice A vector whose expected values are all zero would make the test above vacuous.
    function testTheVectorsAreNonTrivial() public pure {
        uint256 nonZero = 0;
        for (uint256 i = 0; i < RelationVectors.COUNT; ++i) {
            RelationVectors.Case memory c = RelationVectors.get(i);
            for (uint256 k = 0; k < RelationVectors.NUM_SUBRELATIONS; ++k) {
                if (c.expected[k] != 0) nonZero++;
            }
        }
        assertEq(nonZero, RelationVectors.COUNT * RelationVectors.NUM_SUBRELATIONS);
    }
}
