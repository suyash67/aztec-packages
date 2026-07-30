# Small fields in barretenberg: M31 and the cost of leaving BN254

This note answers "can bb easily switch UltraHonk from BN254 Fr to a small field like M31?" with
a working M31/QM31 implementation (`m31.hpp`), measurements on this machine (`m31_bench`), and an
inventory of what in the stack actually depends on Fr. Summary: the *field* is a day of work (done
here); the *switch* is not a configuration change but a second proving stack, because every layer
above the field — codes, transcript, sumcheck soundness, the gate set, and the Noir frontend — is
load-dependent on Fr's properties.

## 1. What is implemented

`bb::small_field::m31` (p = 2^31 − 1, 4-byte values, shift-add Mersenne reduction, no Montgomery
form), `cm31` = F_p[i], and `qm31` = CM31[u]/(u² − (2+i)) — the ~124-bit "secure field" that
Fiat-Shamir challenges must be drawn from when witnesses are 31-bit (stwo/Plonky3 convention).
Field-axiom, edge-case, and extension-structure tests.

## 2. Measurements (M4 Pro, single-threaded, arm64 build without Fr assembly)

| Workload (2^20 elements) | M31 | mixed M31xQM31 | QM31 | BN254 Fr | ratio to Fr |
|---|---|---|---|---|---|
| multiplication stream | 555 Mmul/s | — | 61 Mmul/s | 70 Mmul/s | **7.9x** (M31); **0.87x** (QM31) |
| sumcheck round kernel | 2.77 ms | 7.45 ms | — | 27.3 ms | **9.9x** / **3.7x** |
| Blake3 over the values | 7.5 ms (4 MiB) | — | — | 60.2 ms (32 MiB) | **8.0x** |

Three facts to hold together:

1. **Base-field work is ~8–10x cheaper** — the small-field promise is real where computation
   stays in M31 (first sumcheck round, commitment hashing, encoding).
2. **The challenge field gives it back.** A QM31 multiplication costs about one Fr
   multiplication. Every sumcheck round after the first folds with QM31 challenges, so the
   realistic sumcheck gain is the *mixed* row: **~3.7x**, not 10x.
3. **The Fr baseline here is handicapped** (this build disables the x86 field assembly); on
   deployment hardware the ratios shrink by roughly 2x. Conversely an M31 implementation with
   NEON/AVX batching would gain a similar factor. Treat the ratios as the right order, not
   final values.

And one accounting caveat that dominates everything: these are per-*element* comparisons. An Fr
element carries 254 bits of payload to M31's 31. Circuits with genuinely small-valued semantics
(hashing, u32 arithmetic — the Flock/stwo regime, see `../PCS_CANDIDATES.md`) really do get the
element-count for free; Fr-semantic circuits (Poseidon over Fr, elliptic-curve gates, field
arithmetic — the Aztec workload) would expand ~8x in element count when re-arithmetized, eating
the per-element advantage before it starts.

## 3. Port-surface inventory: what in bb assumes Fr

| Layer | Fr dependence | M31 status |
|---|---|---|
| field type | `field<Params>` is 4x64 Montgomery-specific | done: standalone class (easy) |
| FFT / RS codes | `SupportsFFT` needs high 2-adicity; Fr has 2^28 | **fundamental**: M31's multiplicative group has 2-adicity 1 — no radix-2 domains exist. All RS-code backends (WHIR, Ligero, FRI) require the circle FFT over the order-2^31 circle group (Circle STARKs, [2024/278](https://eprint.iacr.org/2024/278)): new domains, new twiddle structure, circle-variant folding; circle-WHIR is research, not a port |
| EC backends | KZG/Mercury need a pairing curve with the field as scalar field | **impossible at useful security** — no such curve for a 31-bit field; both pairing backends are simply gone |
| transcript | `FrCodec` and Poseidon2 parameters are BN254-specific | moderate: M31 codec + an M31 sponge (Poseidon2-M31 exists in the literature) or a byte-hash transcript |
| sumcheck soundness | one field everywhere; 254-bit challenges | **rearchitecture**: 31-bit challenges are unsound; every protocol needs the Base/Challenge split (M31 witnesses, QM31 challenges) threaded through `Sumcheck`, all relations, `Univariate`, barycentric tables — the type signatures of the entire proving core |
| gate set / flavor | Poseidon2 gates hardcode Fr round constants; the elliptic gate is Grumpkin arithmetic in Fr; NNF gates emulate Fq; range tables assume Fr decomposition | new arithmetization: an "UltraM31Flavor" shares nothing semantic with Ultra — this is building a stwo/Plonky3-class system |
| frontend | Noir/ACIR compile to Fr circuits | ecosystem break: no Noir programs run on an M31 backend |

## 4. Verdict

"Easily switch": **no**. The four implemented backends cover the design space *available at
BN254*: switching fields removes the two pairing backends outright and forces the two hash
backends onto circle-group codes that do not exist in bb. Above the PCS, the honest cost is the
Base/Challenge-field rearchitecture of sumcheck and a new gate set, i.e. a parallel proving stack
with a different frontend — at which point one is re-deriving stwo rather than porting Honk.

What the measurements support instead: (i) the mixed-field sumcheck row (~3.7x) and the hashing
row (8x) quantify the ceiling a small-field stack offers on hash-heavy Boolean workloads,
consistent with Flock/Binius64's positioning; (ii) for bb the leverage path remains BN254 Honk
with the transparent backends for PCS properties, plus small-field systems as recursion-bridged
coprocessors for bulk Boolean work (`../PCS_CANDIDATES.md`).
