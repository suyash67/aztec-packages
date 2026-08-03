# Brakedown — a linear-time, field-agnostic linear code

Implementation of the Spielman-style recursive code of Brakedown
([eprint 2021/1043](https://eprint.iacr.org/2021/1043), Algorithm 1). This is a *code*, not a PCS:
it is the encoder that `ligero/`, `switchfold/` and the Lightning/Bolt constructions are all
parameterized over, and the only one in the tree that needs no field structure beyond size.

## 1. The construction

```
Enc_n(x):
    y = x · A            A is n × ⌈αn⌉,      c_n nonzeros per row
    z = Enc(y)           recursive, on the contracted message
    v = z · B            B is |z| × (rn−n−|z|),  d_n nonzeros per row
    return (x, z, v)
```

Systematic, linear, rate `1/r`, relative distance `β/r`. Each level costs `n·c_n + |z|·d_n`
multiplications and the message contracts by `α`, so the total is a geometric series: `O(n)`, with
no FFT and no smooth subgroup. `c_n` and `d_n` come from the paper's Equations (7) and (8) and are
computed per level, since both depend on the level's length.

The recursion bottoms out at `base_length` with a systematic Reed–Solomon-style parity block. At
that size a dense block is free asymptotically, and being MDS its relative distance `(r−1)/r`
comfortably exceeds the `β/r` the recursion needs from it.

Matrices are sampled deterministically from `(seed, params, length)` via a Blake3 stream, so prover
and verifier derive bit-identical codes without transmitting anything.

## 2. Parameters

`BrakedownParams` exposes the six rows of the paper's Figure 2 as named presets, each chosen so the
code misses its distance with probability at most `2^-100` for messages up to `2^30` over a field of
at least `2^127` elements — BN254 Fr is larger, so the bound carries.

| preset | α | β | r | distance | rate | c_n | d_n |
|---|---|---|---|---|---|---|---|
| `distance_2pct` | 0.1195 | 0.0284 | 1.42 | 0.02 | 0.704 | 6 | 33 |
| `distance_3pct` | 0.138 | 0.0444 | 1.47 | 0.03 | 0.68 | 7 | 26 |
| `distance_4pct` | 0.178 | 0.061 | 1.521 | 0.04 | 0.65 | 7 | 22 |
| `distance_5pct` | 0.2 | 0.082 | 1.64 | 0.05 | 0.60 | 8 | 19 |
| `distance_6pct` | 0.211 | 0.097 | 1.616 | 0.06 | 0.61 | 9 | 21 |
| `distance_7pct` | 0.238 | 0.1205 | 1.72 | 0.07 | 0.58 | 10 | 23 |

`BrakedownParamsTest.SparsitiesMatchPaperFigure2` re-derives these from Equations (7) and (8).
Every `c_n` reproduces exactly. For `d_n` the test asserts the real-valued bound rather than its
ceiling: Figure 2 quotes α, β and r to three or four significant figures, and at those rounded
values two rows evaluate to 26.0033 and 22.0014 against published 26 and 22 — the true optima sit a
hair below integers the ceiling rounds up. The shipped integer is always the paper's value or one
above, never below, and one extra nonzero per row only helps distance.

## 3. Why the matrices are stored column-major

The construction specifies `A` and `B` row-wise, but the product the encoder needs is
`vector · matrix`. Row-major storage makes that a scatter (`out[index] += …`), which races — so it
cannot be parallelized — and whose random writes miss cache. Transposing to compressed-column form
at build time turns each output entry into an independent gather.

This is not a micro-optimization; it decides whether the code is competitive at all. Measured
encode time at 2^16 went from **29.1 ms to 6.3 ms**, a 4.6x improvement, and flipped the comparison
against Reed–Solomon from a loss into a win.

## 4. Against Reed–Solomon

Encoding a message of length `n`, `distance_7pct` against `whir/rs_code.hpp` at rate 1/4 (what
Ligero uses today), on an Apple M4 Pro arm64 build with `DISABLE_ASM=1`:

| n | Brakedown | RS (rate 1/4) | speedup | Brakedown symbols | RS symbols |
|---|---|---|---|---|---|
| 2^12 | 0.98 ms | 1.11 ms | 1.13x | 7 046 | 16 384 |
| 2^14 | 1.97 ms | 2.27 ms | 1.15x | 28 181 | 65 536 |
| 2^16 | 6.32 ms | 7.69 ms | 1.22x | 112 722 | 262 144 |

The measured cost is 26.9 multiplications per message symbol at 2^16, against the 25.5n the paper
predicts for this row — close enough to treat the implementation as faithful. The speedup widens
with `n`, which is the `Θ(n)` versus `Θ(n log n)` crossover doing its work, and Brakedown emits 2.3x
fewer symbols, so its Merkle trees are correspondingly smaller.

**But the distance is the catch, and it dominates.** Brakedown's relative distance is 0.07 against
Reed–Solomon's 0.75 at rate 1/4 under the capacity conjecture. Under the provable interleaved
proximity test, per-query soundness is `1 − δ/3`, so at λ = 100:

- RS at rate 1/4: `⌈100/2⌉` = **50 queries**
- Brakedown at δ = 0.07: `100 / −log₂(1 − 0.07/3)` ≈ **2936 queries**

Roughly 59x more openings, which is why Brakedown's own paper reports opening 6593 columns and why
its proofs are measured in tens of megabytes.

## 5. End to end: `BrakedownHonk`

`brakedown_honk.hpp` runs Ligero's tensor PCS over this code (the protocol is code-agnostic; see
`ligero/RSCodePolicy`), so `LigeroHonk` and `BrakedownHonk` are a controlled A/B — same trace, same
tensor protocol, same Merkle hasher, only the row code and the query rule differ. Measured with the
Blake3s hasher:

| | Ligero (RS) | Brakedown | ratio |
|---|---|---|---|
| prove 2^14 | 40.3 ms | 78.0 ms | **1.9x slower** |
| verify 2^14 | 8.07 ms | 174 ms | **21.6x slower** |
| proof 2^14 | 833 KiB | 20.2 MiB | **24.8x larger** |

**The faster encoder does not even win the prover.** Opening 2936 Merkle paths costs more than the
~20% saved on encoding, so Brakedown loses on all three axes at once. The 59x query ratio shows up
as ~25x proof size (the combined rows and sumcheck do not scale with queries, which dilutes it) and
~22x verification, since the verifier authenticates every one of those openings.

This is the honest baseline the linear-time-code line has to beat, and it is why the module's value
is as a reusable encoder rather than as a competitive backend. It also locates the real target
precisely: not encoding speed, but distance.

## 6. What it unblocks

Lightning ([2026/258](https://eprint.iacr.org/2026/258)) and Bolt
([2026/310](https://eprint.iacr.org/2026/310)) both need a constant-relative-distance linear-time
base code; Lightning applies its sparsity-preserving compression `C_L(m) = m ‖ C_D(mA)` on top of
one, and Bolt's sketched random LDPCs are a sibling family. This module is the prerequisite both
were blocked on (see `../PCS_CANDIDATES.md` §7). Their contribution is precisely to keep this
encoding speed while recovering distance — the axis the table above shows is the binding one.
