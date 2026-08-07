#!/usr/bin/env python3
"""Aggregate hash_driver.sh's raw runs into a per-hash comparison table.

Usage: hash_report.py <raw.jsonl>

Reports medians across reps. Absolute wall clock on a laptop is not reproducible across sessions,
so the ratio columns (each side against its own fastest hash) are the numbers to trust; see
PCS_REVIEW.md on benchmarking method.
"""

import json
import statistics
import sys

HASHES = ["skyscraper", "sha256", "blake3", "poseidon2"]
SIDES = ["barretenberg", "provekit"]


def median(rows, key):
    values = [r[key] for r in rows if key in r]
    return statistics.median(values) if values else float("nan")


def main() -> None:
    rows = [json.loads(line) for line in open(sys.argv[1]) if line.strip()]
    grouped = {
        (side, hash_): [r for r in rows if r["side"] == side and r["hash"] == hash_]
        for side in SIDES
        for hash_ in HASHES
    }

    for side in SIDES:
        present = [h for h in HASHES if grouped[(side, h)]]
        if not present:
            continue
        proves = {h: median(grouped[(side, h)], "prove_ms") for h in present}
        best = min(proves.values())
        print(f"\n=== {side} (reps: {len(grouped[(side, present[0])])}) ===")
        print(f"{'hash':<12}{'prove_s':>10}{'xfastest':>10}{'verify_ms':>11}{'proof_KiB':>11}{'RSS_GiB':>9}")
        for h in present:
            g = grouped[(side, h)]
            print(
                f"{h:<12}{proves[h] / 1000:>10.2f}{proves[h] / best:>10.2f}"
                f"{median(g, 'verify_ms'):>11.1f}{median(g, 'proof_bytes') / 1024:>11.1f}"
                f"{median(g, 'peak_rss_bytes') / 2**30:>9.2f}"
            )
        if not all(r.get("verified", True) for h in present for r in grouped[(side, h)]):
            print("  WARNING: at least one run failed verification")


if __name__ == "__main__":
    main()
