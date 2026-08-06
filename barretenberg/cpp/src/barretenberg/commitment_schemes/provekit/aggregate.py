#!/usr/bin/env python3
"""Reduce raw pcs_acir_bench runs to per-backend medians and re-embed them into the HTML report.

Usage: aggregate.py raw_runs.jsonl results.json provekit_pcs_report.html

Every metric is the median across fresh-process repetitions; ratios are medians divided by the
same-file KZG medians, so all normalized numbers come from the same sweep (absolute wall clock on
this machine is not reproducible across sessions, ratios are).
"""

import json
import platform
import statistics
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

def _machine_description() -> str:
    # platform.machine() and uname both report the calling process's architecture, which lies
    # under Rosetta; hw.optional.arm64 reports the hardware.
    def sysctl(key: str) -> str:
        return subprocess.run(["sysctl", "-n", key], capture_output=True, text=True, check=False).stdout.strip()

    cpu = sysctl("machdep.cpu.brand_string")
    arch = "arm64" if sysctl("hw.optional.arm64") == "1" else platform.machine()
    base = (cpu + " " if cpu else "") + arch
    return base + " (DISABLE_ASM build; compare ratios, not absolutes)"


MARKER_START = "/*RESULTS_JSON_START*/"
MARKER_END = "/*RESULTS_JSON_END*/"

METRICS = ["prove_ms", "verify_ms", "proof_bytes", "peak_rss_bytes"]


def main() -> None:
    raw_path, results_path, html_path = sys.argv[1], sys.argv[2], sys.argv[3]
    runs = [json.loads(line) for line in Path(raw_path).read_text().splitlines() if line.strip()]
    if not runs:
        sys.exit("no runs in " + raw_path)
    bad = [r for r in runs if not r.get("verified")]
    if bad:
        sys.exit(f"{len(bad)} runs failed verification: {sorted({r['pcs'] for r in bad})}")

    by_pcs: dict[str, list[dict]] = {}
    for run in runs:
        by_pcs.setdefault(run["pcs"], []).append(run)

    backends = {}
    for pcs, pcs_runs in by_pcs.items():
        entry = {"reps": len(pcs_runs), "log_n": pcs_runs[0]["log_n"], "num_gates": pcs_runs[0]["num_gates"]}
        for metric in METRICS + ["pk_ms", "circuit_ms"]:
            entry[metric] = statistics.median(r[metric] for r in pcs_runs)
        backends[pcs] = entry

    if "kzg" not in backends:
        sys.exit("KZG baseline missing; ratios need a same-run KZG measurement")
    kzg = backends["kzg"]
    for entry in backends.values():
        for metric in METRICS:
            entry["ratio_" + metric] = entry[metric] / kzg[metric] if kzg[metric] else None

    try:
        commit = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True, check=True
        ).stdout.strip()
    except Exception:
        commit = "unknown"

    results = {
        "circuit": {
            "name": "ProveKit complete_age_check (World ID RSA ePassport age check)",
            "source": "github.com/worldfnd/ProveKit noir-examples/noir-passport-monolithic",
            "log_n": kzg["log_n"],
            "num_gates": kzg["num_gates"],
        },
        "machine": _machine_description(),
        "commit": commit,
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "reps": min(e["reps"] for e in backends.values()),
        "backends": backends,
    }
    Path(results_path).write_text(json.dumps(results, indent=1) + "\n")

    html_file = Path(html_path)
    if html_file.exists():
        html = html_file.read_text()
        start = html.index(MARKER_START) + len(MARKER_START)
        end = html.index(MARKER_END)
        html_file.write_text(html[:start] + "\n" + json.dumps(results) + "\n" + html[end:])
        print(f"embedded results into {html_path}")
    else:
        print(f"warning: {html_path} not found; wrote {results_path} only")


if __name__ == "__main__":
    main()
