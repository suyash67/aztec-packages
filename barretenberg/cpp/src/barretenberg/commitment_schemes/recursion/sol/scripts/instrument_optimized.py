#!/usr/bin/env python3
"""Test-only copy of the generated `--optimized` verifier that reports its elliptic-curve gas.

The optimized verifier is a single assembly block that returns from inside the block, so its phases
cannot be timed from Solidity. This inserts a `gas()` reading at the boundary between the field work
(transcript, sumcheck, relations, batching scalars) and the batch multiplication that opens the
commitments, and returns the difference alongside the verification result. 0x9000 is the first word
above the generated memory map, whose top is 0x8e60.
"""
import pathlib
import sys

source = pathlib.Path(sys.argv[1]).read_text()
source = source.replace("contract HonkVerifier is IVerifier {", "contract HonkVerifierInstrumented {")
source = source.replace("import {IVerifier} from", "// import {IVerifier} from")
source = source.replace("    function verify(\n", "    function verifyWithGas(\n", 1)
source = source.replace(
    "    )\n        public\n        view\n        override\n        returns (bool)\n    {",
    "    )\n        public\n        view\n        returns (bool, uint256)\n    {",
)
source = source.replace(
    """            {
                // The initial accumulator = 1 * shplonk_q
                mcopy(ACCUMULATOR, SHPLONK_Q_X_LOC, 0x40)
            }""",
    """            {
                mstore(0x9000, gas())
                // The initial accumulator = 1 * shplonk_q
                mcopy(ACCUMULATOR, SHPLONK_Q_X_LOC, 0x40)
            }""",
)
source = source.replace(
    """                {
                    mstore(0x00, 0x01)
                    return(0x00, 0x20) // Proof succeeded!
                }""",
    """                {
                    let ecGas := sub(mload(0x9000), gas())
                    mstore(0x00, 0x01)
                    mstore(0x20, ecGas)
                    return(0x00, 0x40) // Proof succeeded!
                }""",
)
if "verifyWithGas" not in source or "0x9000" not in source:
    raise SystemExit("instrumentation did not apply; the generated verifier's shape changed")
pathlib.Path(sys.argv[2]).write_text(source)
