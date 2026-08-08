#!/usr/bin/env bash
# Recursive verification of WHIR-Honk and Ligero-Honk in Noir: build, measure, benchmark.
#
#   ./driver.sh sizes      sweep PCS parameters and report each verifier circuit's gate count
#   ./driver.sh outer      prove/verify the two verifier circuits with UltraHonk+KZG and +Vela
#   ./driver.sh gas        Ethereum gas for both openings and the shared Honk work
#   ./driver.sh ipa-gas    what a transparent (IPA) outer PCS would cost natively and on chain
#   ./driver.sh all        all of the above
#
# Env: BB_BUILD (default ../../../../build-arm64), NARGO (default ~/.nargo/bin/nargo),
#      WORK (default ./workdir), REPS (default 5).
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BB_BUILD="${BB_BUILD:-$here/../../../../build-arm64}"
NARGO="${NARGO:-$HOME/.nargo/bin/nargo}"
FORGE="${FORGE:-$HOME/.foundry/bin/forge}"
WORK="${WORK:-$here/workdir}"
REPS="${REPS:-5}"
export PATH="$(dirname "$FORGE"):$PATH"

exporter="$BB_BUILD/bin/pcs_noir_export_bench"
acir_bench="$BB_BUILD/bin/pcs_acir_bench"
vela_export="$BB_BUILD/bin/vela_evm_export_bench"
bb="$BB_BUILD/bin/bb"
noir="$here/noir"

mkdir -p "$WORK"

# Regenerate one verifier crate's params/inputs, compile it, and print its gate count.
# usage: measure_size <label> <crate> <exporter args...>
measure_size() {
  local label="$1" crate="$2"
  shift 2
  local out="$WORK/$label"
  mkdir -p "$out"
  "$exporter" "$@" --out "$out" 2>"$out/export.log"
  cp "$out/params.nr" "$noir/$crate/src/params.nr"
  cp "$out/Prover.toml" "$noir/$crate/Prover.toml"
  (cd "$noir/$crate" && "$NARGO" execute --silence-warnings >/dev/null)
  (cd "$noir/$crate" && "$NARGO" compile --silence-warnings)
  local gates
  gates=$("$bb" gates -b "$noir/target/${crate}.json" 2>/dev/null | python3 -c \
    'import json,sys; d=json.load(sys.stdin)["functions"][0]; print(f"{d["acir_opcodes"]} {d["circuit_size"]}")')
  echo "$label $(tail -1 "$out/export.log") | acir_opcodes+circuit_size $gates"
}

cmd_sizes() {
  echo "== recursive verifier circuit sizes =="
  measure_size whir-r4-repaired-k4-k0_2 whir_verifier \
    --pcs whir --log-n 10 --rate 4 --soundness repaired --k 4 --k0 2
  measure_size whir-r4-provable-k4-k0_2 whir_verifier \
    --pcs whir --log-n 10 --rate 4 --soundness provable --k 4 --k0 2
  measure_size whir-r2-repaired-k4-k0_2 whir_verifier \
    --pcs whir --log-n 10 --rate 2 --soundness repaired --k 4 --k0 2
  measure_size whir-r6-repaired-k4-k0_2 whir_verifier \
    --pcs whir --log-n 10 --rate 6 --soundness repaired --k 4 --k0 2
  measure_size whir-r4-repaired-k4-k0_1 whir_verifier \
    --pcs whir --log-n 10 --rate 4 --soundness repaired --k 4 --k0 1
  measure_size whir-m12-r4-repaired whir_verifier \
    --pcs whir --log-n 12 --rate 4 --soundness repaired --k 4 --k0 2
  measure_size whir-m14-r4-repaired whir_verifier \
    --pcs whir --log-n 14 --rate 4 --soundness repaired --k 4 --k0 2
  measure_size ligero-r4 ligero_verifier --pcs ligero --log-n 10 --rate 4
  measure_size ligero-r2 ligero_verifier --pcs ligero --log-n 10 --rate 2
  measure_size ligero-m12-r4 ligero_verifier --pcs ligero --log-n 12 --rate 4
}

cmd_outer() {
  echo "== outer proving: UltraHonk over each recursive verifier =="
  : >"$WORK/outer_runs.jsonl"
  for _ in $(seq "$REPS"); do
    for crate in ligero_verifier whir_verifier; do
      for pcs in kzg vela; do
        "$acir_bench" -b "$noir/target/$crate.json" -w "$noir/target/$crate.gz" --pcs "$pcs" 2>/dev/null |
          python3 -c "import json,sys; d=json.load(sys.stdin); d['circuit']='$crate'; print(json.dumps(d))" \
            >>"$WORK/outer_runs.jsonl"
      done
    done
  done
  python3 - "$WORK/outer_runs.jsonl" <<'PY'
import collections, json, statistics, sys
rows = [json.loads(line) for line in open(sys.argv[1])]
groups = collections.defaultdict(list)
for row in rows:
    groups[(row['circuit'], row['pcs'])].append(row)
print(f"{'circuit':17} {'pcs':6} {'prove_ms':>9} {'verify_ms':>10} {'proof_B':>8} {'peak_MiB':>9}")
for key in sorted(groups):
    runs = groups[key]
    print(f"{key[0]:17} {key[1]:6} {statistics.median(r['prove_ms'] for r in runs):9.1f} "
          f"{statistics.median(r['verify_ms'] for r in runs):10.2f} {runs[0]['proof_bytes']:8} "
          f"{statistics.median(r['peak_rss_bytes'] for r in runs) / 1048576:9.0f}")
PY
}

cmd_gas() {
  echo "== Ethereum gas =="
  local crate=ligero_verifier
  local out="$WORK/gas"
  mkdir -p "$out"
  "$bb" write_vk -b "$noir/target/$crate.json" -o "$out" -t evm-no-zk >/dev/null 2>&1
  "$bb" prove -b "$noir/target/$crate.json" -w "$noir/target/$crate.gz" -k "$out/vk" -o "$out" \
    -t evm-no-zk >/dev/null 2>&1
  "$bb" write_solidity_verifier -k "$out/vk" -o "$here/sol/src/Verifier.sol" -t evm-no-zk >/dev/null 2>&1
  # The `--optimized` backend is the one Aztec deploys for its L1 rollup verifier.
  "$bb" write_solidity_verifier -k "$out/vk" -o "$here/sol/src/VerifierOpt.sol" -t evm-no-zk --optimized \
    >/dev/null 2>&1
  python3 "$here/sol/scripts/instrument_optimized.py" "$here/sol/src/VerifierOpt.sol" \
    "$here/sol/test/VerifierOptInstrumented.sol"
  "$vela_export" -b "$noir/target/$crate.json" -w "$noir/target/$crate.gz" -o "$out/vela_opening.bin" 2>&1 |
    tail -1
  "$BB_BUILD/bin/fflonk_evm_export_bench" -b "$noir/target/$crate.json" -w "$noir/target/$crate.gz" \
    -o "$out/fflonk_opening.bin" 2>&1 | tail -1
  (cd "$here/sol" && PROOF_DIR="$out" "$FORGE" test -vv 2>&1 |
    grep -E "kzg_|vela_|fflonk_|opt_|execution_gas|calldata_gas|total_tx_gas|proof_bytes|public_inputs")
}

cmd_ipa_gas() {
  echo "== IPA: the on-chain cost of an O(n) opening =="
  for crate in ligero_verifier whir_verifier; do
    "$acir_bench" -b "$noir/target/$crate.json" -w "$noir/target/$crate.gz" --pcs ipa 2>/dev/null
  done
  (cd "$here/sol" && "$FORGE" test --match-contract IpaGasTest -vv --gas-limit 100000000000 2>&1 |
    grep -E "msm_|gas_per_generator|projected_|folding_gas")
}

case "${1:-all}" in
sizes) cmd_sizes ;;
ipa-gas) cmd_ipa_gas ;;
outer) cmd_outer ;;
gas) cmd_gas ;;
all)
  cmd_sizes
  cmd_outer
  cmd_gas
  cmd_ipa_gas
  ;;
*)
  echo "usage: $0 <sizes|outer|gas|ipa-gas|all>" >&2
  exit 1
  ;;
esac
