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

**Other tracked candidates:** Samaritan ([2025/419](https://eprint.iacr.org/2025/419), Mercury's
concurrent sibling — implemented design point covered by Mercury); Greyhound/lattice PCS
(post-quantum with structured verification, wrong-field and immature tooling for BN254 Honk);
FRIttata ([2025/1285](https://eprint.iacr.org/2025/1285), distributed-prover FRI — orthogonal
axis: prover distribution rather than a new cost profile).
