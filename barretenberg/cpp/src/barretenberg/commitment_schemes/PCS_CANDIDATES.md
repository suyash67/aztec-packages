# UltraHonk PCS backends: implemented set and evaluated candidates

Four backends are implemented and benchmarked against each other on identical circuits through
the shared `TransparentHonk` shell (see `whir/`, `ligero/`, `mercury/` and `whir_honk.bench.cpp`):

| Backend | Family | Trust | Strength (measured) |
|---|---|---|---|
| Gemini+Shplonk+KZG | pairing | SRS | balanced incumbent; 13.3 KiB proofs |
| Mercury | pairing, constant-size | SRS | smallest proofs (8.4–9.9 KiB), flat ~4 ms verify |
| WHIR | hash, RS folding | transparent | flat ~12–16 ms verify, ~1.2 MiB proofs, recursion-friendly |
| Ligero | hash, tensor √N | transparent | fastest prover through 2^18 |

## Evaluated and rejected as PCS backends: Binius, Flock

Binius / FRI-Binius ([2024/504](https://eprint.iacr.org/2024/504)) and Flock
([2026/1329](https://eprint.iacr.org/2026/1329)) are proof systems over **binary tower fields**
whose performance comes from the arithmetization, not from a detachable commitment scheme. They
cannot fill UltraHonk's PCS slot:

- **Field mismatch is structural.** UltraHonk's polynomials have BN254 Fr coefficients and its
  sumcheck emits Fr evaluation claims. A binary-field PCS commits to multilinears over
  GF(2^k) towers; no homomorphism connects the two. The only bridge — bit-decomposing each Fr
  coefficient and proving the Fr-weighted recombination (with modular reduction) inside the
  binary system — turns the opening argument into large-prime-field emulation over GF(2), the
  workload binary systems are worst at. Raw committed bits stay comparable, but the evaluation
  argument inflates by orders of magnitude and the protocol would be novel research, strictly
  dominated by the native-Fr hash backends already implemented.
- **Their advantage is Boolean-native witnesses.** Flock reports 42k SHA-256 compressions/s
  on one M4 Max core (9x Binius64, 500x EC SNARKs) because hash circuits are bitwise and
  binary witnesses cost 1 bit each. UltraHonk traces are Fr-native (Poseidon, Pedersen, EC,
  field arithmetic); expressed over GF(2), each Fr multiplication costs thousands of Boolean
  gates, inverting the comparison.

**What is worth taking:** (i) Flock's lincheck/zerocheck optimizations are field-agnostic
sumcheck engineering, relevant to Honk's sumcheck rather than its PCS; (ii) the realistic
integration of either system into an Aztec-like stack is a *Boolean coprocessor* — proving
Keccak/SHA/Blake batches natively and connecting through a recursion bridge — whose open cost is
verifying a binary-tower proof inside an Fr circuit (tower arithmetic emulation), a proof-system
integration study rather than a `commitment_schemes/` addition.

## Next candidates, prioritized

The implemented set spans {pairing, hash} x {logarithmic, constant, folding, tensor} at BN254.
Three cells remain genuinely open, plus one unmeasured axis:

1. **IPA/Bulletproofs over BN254 G1** — the {transparent x small-proof} corner (the Halo2 trade:
   ~2-3 KiB proofs, no setup, O(N) verifier). Cheapest addition: bb's IPA template and Shplemini's
   IPA support exist (Grumpkin/ECCVM); the work is an independent-generator setup (hash-to-curve —
   the powers-of-tau SRS cannot serve as IPA generators) plus a backend adapter.
2. **Dory** ([2020/1274](https://eprint.iacr.org/2020/1274)) — the strongest missing combination:
   transparent AND O(log n) verifier AND O(log n) proof (~10-20 KiB), linear prover paid in
   pairing-group operations (expect ~5-10x the KZG prover). Needs inner-pairing-product argument
   machinery and GT multi-exponentiation on top of bb's Fq12 arithmetic. Highest research value of
   the remaining points.
3. **KZH / KZH-Fold** ([2025/144](https://eprint.iacr.org/2025/144)) — sublinear opening for both
   parties with *native sublinear accumulation* (reported 50x decider improvement over Nova).
   Accumulation-friendliness is the axis this comparison has not measured and the most
   Aztec-relevant one (Chonk/IVC): KZG-family backends accumulate homomorphically, hash backends
   do not fold natively. Adding KZH would turn the suite into a folding-oriented comparison.
4. **Samaritan** ([2025/419](https://eprint.iacr.org/2025/419)) — intra-corner A/B against
   Mercury; quick on the existing infrastructure, low marginal insight.
5. **Brakedown with expander codes** ([2021/1043](https://eprint.iacr.org/2021/1043)) — FFT-free
   linear-time encoding under our Ligero structure; extends the prover-optimal frontier
   (predictably wins the prover column at all sizes, worse proof constants).
6. **zip proof compression** ([2025/1446](https://eprint.iacr.org/2025/1446)) — not a PCS: a
   compression layer for hash-based proofs that directly attacks WHIR/Ligero's MiB proof sizes.

Tracked but not planned: class-group constant-size transparent PCS
([2025/1233](https://eprint.iacr.org/2025/1233), 10 group elements, shortest transparent PCS —
no GUO arithmetic in bb, seconds-scale provers); Zeromorph/PST/HyperKZG (the Shplemini corner,
already represented); BaseFold/STIR (dominated by WHIR; useful only as ablation baselines);
Hyrax (dominated by Dory/KZH in its corner); Greyhound/lattice (wrong ring);
FRIttata ([2025/1285](https://eprint.iacr.org/2025/1285), distributed proving — an orthogonal
axis).
