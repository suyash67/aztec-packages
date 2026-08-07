#!/usr/bin/env bash
# Compares barretenberg's transparent UltraHonk+WHIR against ProveKit's own prover across Merkle
# hash functions, at matched WHIR parameters.
#
# Both sides run 128-bit security under the Johnson bound (`unique_decoding: false`), rate 2^-2 and
# folding factor 3 — no conjectural assumption on either side. They differ in where the bits come
# from: ProveKit takes 118 from queries plus 10 bits of per-round grinding, barretenberg implements
# no grinding and takes all 128 from queries, which is the more conservative of the two. Everything
# else about the schedule is identical, and `WhirConfigTest.ProveKitScheduleMatchesReferenceImplementation`
# pins that.
#
# The hash selects the Merkle tree compression on both sides (and, on ProveKit's side, the
# Fiat-Shamir sponge and public-input binding too — its HashConfig ties all three together).
#
# Env: REPS (default 3), PROVEKIT_CLI (required), BENCH_BIN, WORK_DIR.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

REPS="${REPS:-3}"
WORK_DIR="${WORK_DIR:-$PWD/workdir}"
BENCH_BIN="${BENCH_BIN:-$PWD/../../../../build-arm64/bin/pcs_acir_bench}"
PROVEKIT_CLI="${PROVEKIT_CLI:-$WORK_DIR/ProveKit/target/release/provekit-cli}"
SECURITY_BITS=128
RAW="$PWD/hash_runs.jsonl"

CIRCUIT_DIR="$WORK_DIR/complete_age_check/target"
[[ -f "$CIRCUIT_DIR/complete_age_check.json" ]] || { echo "run driver.sh first to fetch the circuit" >&2; exit 1; }
[[ -x "$PROVEKIT_CLI" ]] || { echo "no provekit-cli at $PROVEKIT_CLI" >&2; exit 1; }

# ProveKit's hash name -> barretenberg's backend suffix.
declare -A BB_SUFFIX=([skyscraper]=sky [sha256]=sha256 [blake3]=b3 [poseidon2]=p2)
HASHES=(skyscraper sha256 blake3 poseidon2)

# provekit-cli `prepare` recompiles the package with ProveKit's own Noir frontend and overwrites
# target/, clobbering the nargo artifact barretenberg consumes — keep it in a package copy.
PKG_DIR="$WORK_DIR/provekit-pkg"
if [[ ! -d "$PKG_DIR" ]]; then
    cp -R "$(dirname "$CIRCUIT_DIR")" "$PKG_DIR"
    rm -rf "$PKG_DIR/target"
fi

# Phase A: one scheme per hash. Keys are cached; this is the slow part (minutes per hash).
for hash in "${HASHES[@]}"; do
    keys="$WORK_DIR/provekit-keys-$hash"
    mkdir -p "$keys"
    if [[ ! -f "$keys/scheme.pkp" ]]; then
        echo "--- provekit prepare hash=$hash"
        (cd "$PKG_DIR" && "$PROVEKIT_CLI" prepare --hash "$hash" --pkp "$keys/scheme.pkp" --pkv "$keys/scheme.pkv")
    fi
done

# Phase B: measured reps, strictly sequential — a concurrent run skews these badly.
: > "$RAW"
for rep in $(seq 1 "$REPS"); do
    for hash in "${HASHES[@]}"; do
        echo "--- rep $rep/$REPS  barretenberg  hash=$hash"
        PCS_ACIR_BENCH_SECURITY_BITS="$SECURITY_BITS" "$BENCH_BIN" \
            -b "$CIRCUIT_DIR/complete_age_check.json" -w "$CIRCUIT_DIR/witness.gz" \
            --pcs "whir-provekit-${BB_SUFFIX[$hash]}" --flavor provekit |
            python3 -c "import sys,json;d=json.load(sys.stdin);d['side']='barretenberg';d['hash']='$hash';print(json.dumps(d))" >> "$RAW"
        tail -1 "$RAW"

        echo "--- rep $rep/$REPS  provekit      hash=$hash"
        python3 measure_provekit.py "$PROVEKIT_CLI" "$PKG_DIR" \
            "$WORK_DIR/provekit-keys-$hash/scheme.pkp" "$WORK_DIR/provekit-keys-$hash/scheme.pkv" |
            python3 -c "import sys,json;d=json.load(sys.stdin);d['side']='provekit';d['hash']='$hash';print(json.dumps(d))" >> "$RAW"
        tail -1 "$RAW"
    done
done

python3 hash_report.py "$RAW"
