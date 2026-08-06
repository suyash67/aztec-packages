#!/usr/bin/env python3
"""Run one fresh-process ProveKit prove+verify rep and emit a raw-runs JSON line.

Usage: measure_provekit.py <provekit-cli> <package-dir> <pkp> <pkv>

The package dir must contain Prover.toml. Timings come from provekit-cli's own tracing spans (the
top-level `prove:`/`verify:` durations), so key-file loading is excluded — matching the UltraHonk
rows, whose prove_ms excludes proving-key setup. ProveKit's prove span still includes its witness
solve (~3% here); the UltraHonk rows consume a pre-solved witness. Peak RSS is the prove child's
process-wide maximum, which does include the key load.
"""

import json
import os
import re
import resource
import subprocess
import sys
import tempfile

SPAN_UNITS = {"μs": 1e-3, "ms": 1.0, "s": 1e3}


def top_span_ms(stderr: str, label: str) -> float:
    """Largest `<label>: <value> <unit> duration` span in the tracing tree (the outermost one)."""
    # provekit-cli separates value and unit with U+202F (narrow no-break space).
    spans = re.findall(rf"\b{label}: ([\d.]+)[\s ]+(μs|ms|s) duration", re.sub(r"\x1b\[[0-9;]*m", "", stderr))
    if not spans:
        sys.exit(f"no '{label}:' span in provekit-cli output; tracing format changed?")
    return max(float(v) * SPAN_UNITS[u] for v, u in spans)


def main() -> None:
    cli, pkg, pkp, pkv = sys.argv[1:5]
    proof = os.path.join(tempfile.mkdtemp(), "proof.np")

    run = subprocess.run([cli, "prove", "--prover", pkp, "--input", os.path.join(pkg, "Prover.toml"), "--out", proof],
                         cwd=pkg, check=True, capture_output=True, text=True)
    prove_ms = top_span_ms(run.stdout + run.stderr, "prove")
    peak_rss = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    if sys.platform != "darwin":
        peak_rss *= 1024

    run = subprocess.run([cli, "verify", "--verifier", pkv, "--proof", proof],
                         cwd=pkg, check=True, capture_output=True, text=True)
    verify_ms = top_span_ms(run.stdout + run.stderr, "verify")

    print(json.dumps({
        "pcs": "provekit",
        "log_n": 0,
        "num_gates": 0,
        "circuit_ms": 0,
        "pk_ms": 0,
        "prove_ms": round(prove_ms, 2),
        "verify_ms": round(verify_ms, 2),
        "proof_bytes": os.path.getsize(proof),
        "peak_rss_bytes": peak_rss,
        "verified": True,
    }))


if __name__ == "__main__":
    main()
