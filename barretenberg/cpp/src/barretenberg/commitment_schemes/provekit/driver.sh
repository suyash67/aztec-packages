#!/usr/bin/env bash
# Benchmarks the ProveKit World ID passport circuit (complete_age_check) across every PCS backend
# reachable through TransparentHonk, plus the production KZG UltraHonk baseline.
#
# Each (backend, rep) pair runs in a fresh pcs_acir_bench process so peak_rss_bytes is a valid
# per-run maximum. Raw per-run JSON lines land in raw_runs.jsonl; aggregate.py reduces them to
# medians and same-run ratios against KZG in results.json, and re-embeds that JSON into
# provekit_pcs_report.html.
#
# Usage: ./driver.sh [reps]        (default 5)
# Env:   BUILD_DIR   cmake build dir (default: ../../../../build-arm64 relative to this script)
#        CIRCUIT_DIR directory containing complete_age_check.json + witness.gz (default: auto-fetch)

set -euo pipefail
cd "$(dirname "$0")"

REPS="${1:-5}"
BACKENDS=(kzg mercury whir whir-p2 whir-sky whir-sky-stacked ligero hyrax kzh2 ipa dory)
BUILD_DIR="${BUILD_DIR:-$(cd ../../../.. && pwd)/build-arm64}"
BENCH_BIN="$BUILD_DIR/bin/pcs_acir_bench"
WORK_DIR="$PWD/workdir"

if [[ ! -x "$BENCH_BIN" ]]; then
    echo "building pcs_acir_bench..."
    cmake --build "$BUILD_DIR" --target pcs_acir_bench
fi
# Snapshot the binary: a concurrent ninja run can delete/replace it mid-sweep.
mkdir -p "$WORK_DIR"
cp "$BENCH_BIN" "$WORK_DIR/pcs_acir_bench.snapshot"
BENCH_BIN="$WORK_DIR/pcs_acir_bench.snapshot"

# --- circuit: ProveKit's noir-passport-monolithic/complete_age_check ---
if [[ -z "${CIRCUIT_DIR:-}" ]]; then
    CIRCUIT_DIR="$WORK_DIR/complete_age_check/target"
    if [[ ! -f "$CIRCUIT_DIR/complete_age_check.json" || ! -f "$CIRCUIT_DIR/witness.gz" ]]; then
        echo "fetching + compiling ProveKit circuit (needs git, nargo >= 1.0.0-beta.26 on PATH)..."
        rm -rf "$WORK_DIR/ProveKit"
        git clone --depth 1 --filter=blob:none --sparse https://github.com/worldfnd/ProveKit.git "$WORK_DIR/ProveKit"
        git -C "$WORK_DIR/ProveKit" sparse-checkout set noir-examples/noir-passport-monolithic
        SRC="$WORK_DIR/ProveKit/noir-examples/noir-passport-monolithic/complete_age_check"
        (cd "$SRC" && nargo compile && nargo execute witness)
        rm -rf "$WORK_DIR/complete_age_check"
        cp -R "$SRC" "$WORK_DIR/complete_age_check"
    fi
fi
BYTECODE="$CIRCUIT_DIR/complete_age_check.json"
WITNESS="$CIRCUIT_DIR/witness.gz"
[[ -f "$BYTECODE" && -f "$WITNESS" ]] || { echo "missing $BYTECODE or $WITNESS" >&2; exit 1; }

RAW="$WORK_DIR/raw_runs.jsonl"
: > "$RAW"
echo "circuit: $BYTECODE"
echo "reps: $REPS  backends: ${BACKENDS[*]}"

for rep in $(seq 1 "$REPS"); do
    for pcs in "${BACKENDS[@]}"; do
        echo "--- rep $rep/$REPS  pcs=$pcs"
        if ! "$BENCH_BIN" -b "$BYTECODE" -w "$WITNESS" --pcs "$pcs" >> "$RAW"; then
            echo "FAILED: pcs=$pcs rep=$rep" >&2
            exit 1
        fi
        tail -1 "$RAW"
    done
done

# Optional comparison row: ProveKit's own prover (Spartan+WHIR over R1CS, Skyscraper hash) on the
# same Noir package. Point PROVEKIT_CLI at a release build of provekit-cli from
# github.com/worldfnd/ProveKit (cargo build --release --bin provekit-cli); skipped when unset.
PROVEKIT_CLI="${PROVEKIT_CLI:-}"
if [[ -x "$PROVEKIT_CLI" ]]; then
    # provekit-cli `prepare` recompiles the package with ProveKit's own Noir frontend and
    # overwrites target/, clobbering the nargo artifact bb consumes — run it in a package copy.
    SRC_PKG="$(dirname "$CIRCUIT_DIR")" # the Noir package holding Nargo.toml + Prover.toml
    PKG_DIR="$WORK_DIR/provekit-pkg"
    if [[ ! -d "$PKG_DIR" ]]; then
        cp -R "$SRC_PKG" "$PKG_DIR"
        rm -rf "$PKG_DIR/target"
    fi
    PK_KEYS="$WORK_DIR/provekit-keys"
    mkdir -p "$PK_KEYS"
    if [[ ! -f "$PK_KEYS/scheme.pkp" ]]; then
        echo "--- provekit prepare (one-time key generation)"
        (cd "$PKG_DIR" && "$PROVEKIT_CLI" prepare --pkp "$PK_KEYS/scheme.pkp" --pkv "$PK_KEYS/scheme.pkv")
    fi
    for rep in $(seq 1 "$REPS"); do
        echo "--- rep $rep/$REPS  pcs=provekit (native)"
        python3 measure_provekit.py "$PROVEKIT_CLI" "$PKG_DIR" "$PK_KEYS/scheme.pkp" "$PK_KEYS/scheme.pkv" >> "$RAW"
        tail -1 "$RAW"
    done
else
    echo "PROVEKIT_CLI not set/executable; skipping the ProveKit native comparison row"
fi

python3 aggregate.py "$RAW" results.json provekit_pcs_report.html
echo "wrote results.json and provekit_pcs_report.html"
