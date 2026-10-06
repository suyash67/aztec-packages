# Zcash Orchard reference prover

Baseline numbers and test vectors for barretenberg's Orchard Action circuits (`barretenberg/cpp/src/barretenberg/zcash`),
taken from the production Zcash prover: the `orchard` crate (release 0.16.0) on `halo2_proofs` 0.4 and
`pasta_curves` 0.6.

## Setup

The crate depends on a local checkout of `orchard` at `./orchard`:

```bash
cd barretenberg/cpp/scripts/zcash_orchard_reference
git clone https://github.com/zcash/orchard.git && git -C orchard checkout 0.16.0
# Only needed for the test-vector export (src/bin/vectors.rs):
git -C orchard apply ../orchard-export.patch
```

## Binaries

- `cargo run --release -- 1,2,4,8,16 > results/zcash_mt.json` builds Orchard bundles with exactly `n` Actions (one
  output each, dummy spends, no padding) and times `create_proof` and single-proof `verify` with the post-NU6.2
  circuit (`OrchardCircuitVersion::FixedPostNu6_2`), as `orchard/benches/circuit.rs` does. It records the serialized
  proof length and the proving + verifying key generation time. `RAYON_NUM_THREADS=1` gives the single-threaded
  numbers (`results/zcash_st.json`).
- `cargo run --release --bin vectors > ../../src/barretenberg/zcash/test_vectors/orchard_vectors.hpp` regenerates the
  C++ test vectors (Pasta hash-to-curve outputs of every Orchard generator, and complete Action witnesses with their
  public inputs), using the export patch.
- `cargo run --release --bin params` prints halo2's Vesta IPA parameters (`G_0`, `G_1`, `G_15`, `W`, `U` for k = 4),
  against which `zcash/honk/pasta_crs.test.cpp` checks barretenberg's generators.

## barretenberg numbers

`results/bb_*.jsonl` are produced by barretenberg's benchmark binary (one JSON object per line):

```bash
cd barretenberg/cpp && cmake --preset default && cmake --build build --target zcash_action_bench
for system in halo2-honk-pasta halo2-honk-bn254 ultra-zk-pasta ultra-zk-bn254; do
  HARDWARE_CONCURRENCY=4 ./build/bin/zcash_action_bench $system 1,2,4,8,16 5 > scripts/zcash_orchard_reference/results/bb_${system}_mt.jsonl
  HARDWARE_CONCURRENCY=1 ./build/bin/zcash_action_bench $system 1,2,4 3 > scripts/zcash_orchard_reference/results/bb_${system}_st.jsonl
done
```

`ultra-zk-bn254` needs the BN254 SRS in `~/.bb-crs`. The Pasta systems use halo2's transparent generators, derived on
the fly.

## Report

`report/make_report.py results report/orchard-benchmarks.html` renders the HTML report (charts and tables) from the
result files, using `report/report_template.html`.
