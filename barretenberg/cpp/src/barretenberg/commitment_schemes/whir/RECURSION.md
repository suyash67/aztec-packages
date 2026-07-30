# Recursive WHIR in barretenberg

This note analyzes an in-circuit WHIR verifier: its cost model, how it compares with the current
recursive Shplemini (Gemini+Shplonk+KZG) verifier, the parameter and layout choices that recursion
forces, and a design sketch for a `stdlib` implementation. Protocol notation is README.md §2; the
verifier being modeled is `WhirVerifier` as instantiated by `WhirHonk` (Ultra claim shape: 36
committed columns, 5 of them also opened shifted).

Headline conclusions, established in §2–§4:

- **(R1)** An in-circuit WHIR verifier is hashing-dominated and entirely native-field: no biggroup,
  no bigfield, no pairing. With the Poseidon2 hasher every protocol step is field arithmetic plus
  sponge permutations (73 Ultra gates each; ~26–28 with Mega's compressed Poseidon2 encoding).
- **(R2)** The naive per-polynomial-tree layout is impractical for recursion (~5M gates at
  $m = 20$); the cost is round-0 openings. Three layout changes — one shared tree per commitment
  round, deduplicated shifted openings, and excluding Merkle-authenticated data from Fiat-Shamir
  re-absorption — bring it to ~0.5–1.2M Ultra gates depending on rate and grinding.
- **(R3)** That is competitive with a non-goblinized recursive Shplemini verifier (whose ~60-point
  BN254 MSM must be emulated in non-native arithmetic), but not with Goblin-deferred EC recursion,
  which keeps only op-queue writes in-circuit. The case for recursive WHIR is therefore not raw
  gate count today: it is transparency (no SRS), a plausibly post-quantum hash-only stack, and
  removal of the ECCVM/Translator tail from the recursion chain.

## 1. Cost anchors

| Quantity | Value | Source |
|---|---|---|
| Poseidon2 permutation, Ultra encoding | 73 gates | `stdlib/hash/poseidon2/poseidon2.test.cpp` `gate_count` |
| Poseidon2 permutation, Mega compressed encoding | $28 P - 2$ gates for $P$ permutations | same |
| Sponge rate | 3 field elements per permutation | Poseidon2 t=4 sponge |
| Hash of an $N$-element leaf | $\lceil N/3 \rceil$ permutations | sponge absorption |
| Merkle path step | 1 permutation (2-element input) | `MerkleTree::verify` |
| Native field mul-add | ~1 gate | Ultra arithmetic gate |

## 2. Verifier hashing model

For the schedule of README.md §6 the verifier's hashing is

$$H \;=\; \underbrace{t_0 \cdot \Big(\sum_{\text{trees } T} \big(\lceil \text{leaf}_T/3 \rceil + \text{depth}_T\big)\Big)}_{\text{round-0 openings}} \;+\; \underbrace{\sum_{i \ge 1} t_i\,(\lceil 2^k/3 \rceil + \text{depth}_i)}_{\text{inner and final openings}} \;+\; H_{\text{FS}},$$

where $H_{\text{FS}}$ is transcript absorption. Worked at $m = 20$, $r_0 = 2$, $k = 4$,
$\lambda = 100$ (schedule $t = \{50, 20, 13, 10\}$, final $8$; depths $\{18, 17, 16, 15\}$,
final $14$):

| Layout | Round-0 perms | Inner+final perms | $H_{\text{FS}}$ perms | Total perms | Ultra gates |
|---|---|---|---|---|---|
| Per-poly trees, shifted dups, full absorption (current `WhirHonk`) | $50 \cdot 41 \cdot (6{+}18) = 49{,}200$ | $1{,}116$ | $\approx 24{,}000$ | $\approx 74{,}300$ | $\approx 5.4$M |
| Shared trees + dedup + minimal absorption | $50 \cdot 284 = 14{,}200$ | $1{,}116$ | $\approx 300$ | $\approx 15{,}600$ | $\approx 1.14$M |
| same, $r_0 = 3$ ($t_0 = 34$) | $9{,}656$ | $\approx 1{,}050$ | $\approx 300$ | $\approx 11{,}000$ | $\approx 0.80$M |
| same, $r_0 = 4$ + 20-bit grinding ($t_0 = 20$) | $5{,}680$ | $\approx 700$ | $\approx 300$ | $\approx 6{,}700$ | $\approx 0.49$M |

The three layout changes behind rows 2–4:

1. **One tree per commitment round.** Five trees (precomputed 28 columns; wires 3; counts+w_4 3;
   lookup_inverses; z_perm) instead of 36: a round-0 query then costs
   $\lceil 448/3 \rceil + 2\lceil 48/3 \rceil + 6 + 6 = 194$ leaf permutations plus $5 \cdot 18 = 90$
   path permutations $= 284$. Value hashing is irreducible — every opened value must enter a leaf
   hash — so column count is the fundamental driver; the change eliminates the 36-fold path
   duplication.
2. **Deduplicate shifted openings.** A to-be-shifted column's shifted value is derived from the
   same opened leaf ($A(x)/x$); opening its tree twice per query is pure waste.
3. **Absorb only unbound prover messages into Fiat-Shamir.** Opened values and Merkle paths are
   already bound to pre-challenge roots by path verification; re-absorbing them (as the current
   uniform transcript flow does) costs $\approx$ proof-size/96 B permutations for no soundness
   benefit. Only roots, sumcheck univariates, OOD answers, and the final polynomial need absorption.

Field-operation costs are negligible beside hashing: per query one inversion, $k \cdot 2^{k-1}$
fold mul-adds, and an $\omega^{\text{idx}}$ exponentiation by bit-decomposition
($\le \text{depth}$ conditional muls); plus weight-term bookkeeping
($O(\text{terms} \cdot (k + 2^{m_{\text{fin}}}))$) — in total under ~50k gates at these parameters,
all native.

The zk mode (README.md §8) adds one salt element per round-0 opened leaf
($t_0 \cdot 6$ absorptions $\approx 100$ permutations) and one extra column (the mask) — a few
percent.

## 3. Comparison with recursive Shplemini

The recursive Shplemini+KZG verifier's dominant cost is one batched MSM of
$\approx \text{columns} + \log n + O(1) \approx 60$ BN254 points (`Flavor::FINAL_PCS_MSM_SIZE`),
plus the pairing left to the accumulator. Two regimes:

- **Non-goblinized (biggroup):** each point contributes thousands of gates of non-native Fq
  emulation; the MSM alone is order $10^5$–$10^6$ gates. Recursive WHIR at ~0.5–1.2M gates is in
  the same band, with no SRS and no pairing anywhere.
- **Goblinized (current Aztec stack):** in-circuit cost collapses to op-queue writes, with the real
  EC work deferred to ECCVM+Translator and one final native MSM/pairing. In-circuit, Goblin wins
  decisively; end-to-end it retains the ECCVM/Translator provers and the trusted-setup pairing
  check. A WHIR-based recursion chain removes all of that in exchange for the hashing above and
  larger proofs between layers (§9 benchmarks).

A Mega-arithmetized verifier circuit re-prices every permutation from 73 to ~26–28 gates, putting
the optimized configurations at roughly 180–420k gates — below the non-goblinized EC verifier.

## 4. stdlib design sketch

`WhirVerifier` ports mechanically because every operation is already native circuit arithmetic:

- **Transcript:** the existing `StdlibTranscript` (Poseidon2 Fiat-Shamir) replaces
  `NativeTranscript`; digests are `field_t` directly (`Poseidon2MerkleHasher` is the only sensible
  in-circuit hasher — Blake3s costs thousands of gates per compression).
- **Merkle verification:** `stdlib::poseidon2` sponge per leaf/node; path direction bits come from
  the range-constrained bit decomposition of the query index challenge, selecting left/right with
  conditional swaps (2 muls per level).
- **Coset folds and weights:** `field_t` arithmetic as in `rs_code.hpp` / `weights.hpp`;
  $\omega^{\text{idx}}$ via square-and-multiply over the decomposed index bits; the single
  inversion per query is one witness-and-constrain gate.
- **Config:** the schedule is a compile-time/native-computed constant; no in-circuit derivation.
- **Witness volume:** the proof stream (hundreds of KiB at $m = 20$ with the round-layout
  optimizations pending) enters as ~tens of thousands of witnesses — well within builder limits.

Prerequisite before a stdlib port is worthwhile: implement the §2 layout changes natively
(shared per-round trees, dedup, minimal absorption), since they change the proof format the
in-circuit verifier consumes.

## 5. Recommendation

For Aztec's current Goblin-based stack, recursive WHIR does not beat in-circuit EC deferral on
gates. It becomes the right choice when any of these matter: removing the trusted setup end-to-end,
a post-quantum-plausible recursion story, eliminating the ECCVM/Translator provers, or a uniform
hash-only verifier (e.g. final verification layers on hash-friendly L1s). The concrete path:
land the round-layout optimizations, re-benchmark proof size, then port `WhirVerifier` to stdlib
with the Poseidon2 hasher and measure a real recursive-verifier circuit against
`UltraRecursiveFlavor`'s.
