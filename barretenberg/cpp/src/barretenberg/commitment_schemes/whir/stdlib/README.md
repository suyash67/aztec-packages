# Recursive verification of WHIR in barretenberg

A WHIR opening verifier written directly in barretenberg's circuit stdlib, together with the
protocol-level settings that make it cheap and the measurements that justify each one. A WHIR proof
is hundreds of kilobytes and will never go on Ethereum; the way to use WHIR in a rollup is to verify
it inside an UltraHonk circuit and settle the 13 KB pairing-based proof instead.

| | |
|---|---|
| `whir_recursive_verifier.hpp` | the WHIR opening verifier, the mirror of `bb::whir::WhirVerifier::verify` |
| `transparent_honk_recursive_verifier.hpp` | the *whole* proof: Honk shell + sumcheck + relations + the opening |
| `flavor/ultra_provekit_recursive_flavor.hpp` | `UltraProveKitFlavor` over `stdlib::field_t`, so bb's own sumcheck verifier runs in-circuit |
| `recursion_harness.hpp` | proves WHIR openings of a chosen column layout and verifies them in a circuit |
| `whir_recursive_verifier.test.cpp` | completeness, soundness (tampering every part of the proof), aggregation |
| `whir_recursion.bench.cpp` | the gate sweeps below, and the end-to-end outer proof |

## 1. Why WHIR recurses cheaply into UltraHonk

Everything a WHIR verifier does is BN254 `Fr` arithmetic and Poseidon2 hashing. Both are native to an
UltraHonk circuit over the same field: there is no non-native field emulation, no elliptic-curve
arithmetic, and nothing deferred to a Goblin/ECCVM tail. That is the structural reason to reach for a
hash-based scheme when the goal is aggregation, and it is why the circuit below is a few hundred
thousand gates rather than a few million.

The cost is concentrated in one place. **Poseidon2 permutations for Merkle authentication are 85 % of
the circuit**, so the entire design problem is "how many field elements does the verifier have to
hash", which decomposes into

```
permutations ≈ Σ_rounds  queries × [ Σ_groups ⌈(columns_g · 2^k₀)/3⌉   +   groups × (depth − cap) ]
                                    └── leaf hashing ──┘                └── path walk ──┘
```

Every lever below attacks one factor of that expression.

## 2. The recursion profile

Three settings on `WhirConfig`, turned on together by `enable_recursion_profile()`. None of them
touches the soundness argument — the schedule, the query counts and the proximity regime are
unchanged — they only change what the verifier has to do.

**Per-query authentication paths** (`per_query_openings`). The batched layout `WhirProver` normally
uses walks the *distinct ascending* leaf set the query indices induce. That set's shape depends on
the query values, so in-circuit it forces a witness-indexed read at every step of every level, and
the walk has to be unrolled to its worst case anyway — which is exactly the per-query cost, without
the sharing. Sending each query its own path makes every array index a compile-time constant. The
extra digests are *unhashed* proof data, so they cost a recursive verifier nothing but bytes.

**A Merkle cap** (`merkle_cap_levels`). Each path stops `c` levels below the root; the prover sends
that whole level and the verifier folds it to the root once. That trades `2^c − 1` hashes, paid once
per tree, for `c` hashes on every query, which pays whenever `2^c − 1 < c·t`. `cap_levels_for` picks
the optimum for each tree's own query count and depth — a later round queries a smaller oracle fewer
times, and reusing round 0's cap there costs more than it saves.

**A Poseidon2 grind** (`poseidon2_pow`). Grinding is the cheapest way to buy soundness bits, and
soundness bits are queries. But bb's grind is Blake3 — the right choice natively, since the prover
runs it 2^b times — and in-circuit a Blake3s compression measures 4,330 gates with another 2,865 per
field element absorbed just to decompose it (§5), so each grind check runs to roughly ten thousand
gates and the three in one of these schedules would cost about what the removed queries were worth. The Poseidon2 form accepts when `Poseidon2(seed, nonce)` is divisible by `2^b`, which is
one permutation plus two range constraints on the quotient: 1.3 % of the circuit for a fifth off the
query count. It is the one lever with a real prover-side price — see §4.

The second range constraint on the quotient is not decorative. `quotient · 2^b = digest` has a
second solution whenever `digest + r` is divisible by `2^b` and the quotient still fits, which would
accept one extra digest residue and hand a bit of the grind back to the prover; pinning the quotient
to `[0, ⌊r/2^b⌋]` keeps the product below `r`, and an honest quotient is always inside that range.

A fourth change is in the protocol proper: **query indices are cut several to a challenge**
(`detail::draw_query_indices`). A challenge carries ~254 bits and an index needs 12–20 of them, so one
challenge per query wasted a sponge permutation *and* a canonicity argument per query. Packing them
took that phase from 8,050 gates to 2,944.

## 3. Circuit-side choices

**A capacity-separated hasher.** `Poseidon2CompressionHasher` puts the leaf/node domain separator in
the sponge's capacity instead of spending a rate slot on a tag element. A leaf of `n` values costs
`⌈n/3⌉` permutations instead of `⌈(n+1)/3⌉`, and a node is one bare permutation. Separation is at
least as strong: the capacity is never touched by absorbed data, where a tag element shares the rate
with it.

**Sound index extraction.** `challenge = Σ_j index_j·2^{jb} + rest·2^{nb}`, with the digits
bit-decomposed and `rest` pinned to `[0, ⌊r/2^{nb}⌋)` by *two* range constraints. The second one is
not optional: `challenge + r` is below 2^254 for about a third of all challenges, so without it a
prover chooses between two index sets per challenge and claws back a soundness bit per query.

Pinning `rest` strictly below `⌊r/2^{nb}⌋` rather than at it costs completeness, not soundness: it
makes an honest proof unprovable with probability about `2^{nb}/r`, and Fiat-Shamir gives the prover
no second draw. That is why `detail::QUERY_DIGIT_BITS` caps `nb` at 128 rather than filling the
challenge. At 240 bits the figure is 2^-15.3 — roughly one honest proof in five thousand would be
unprovable, which is not a rounding error. At 128 bits it is 2^-128, and the extra challenges cost
545 gates (0.4 %).

The digits are bit-decomposed rather than range-constrained because the bits are needed anyway: one
per Merkle level to place the sibling, and one per level of the square-and-multiply that turns an
index into its domain point `ω^idx`. Each factor of that product, `(ω^{2^i} − 1)·bit + 1`, is a linear
term in one bit, so only the accumulation is a gate.

**Coset evaluation instead of repeated Horner.** The final oracle-consistency check compares a whole
opened coset `{x, xη, …}` against the clear polynomial. Horner at each point costs `arity × |poly|`
multiplications; writing `P(xη^t) = Σ_j (poly[j]·x^j)·η^{tj}` shares the products across the coset and
leaves each output a constant-coefficient linear combination, which an addition gate absorbs several
at a time. That phase went from 12,156 gates to 2,916.

**Measured, not assumed.** One tempting micro-optimization — deriving the right Merkle child as
`digest + sibling − left` so only one multiplication places the pair — measures *worse* (32,448 →
32,864), because the permutation normalizes whatever it is handed and the gate moves rather than
disappears. `field_t::conditional_assign` is already at the floor.

## 4. Results

All figures are gate counts of a real UltraHonk circuit verifying a real WHIR proof, on an M4 Pro
(`DISABLE_ASM=1` arm64 build). Regenerate with `whir_recursion_bench [levers|schedule|columns|outer N]`.

### What each lever is worth

Transparent-Honk column layout — four commitment groups holding 21 committed precomputed columns and
8 witness columns, 5 of them also claimed shifted — at the parameters the Noir verifier's `whir-r4-repaired-k4-k0_2` family was generated at:
m = 10, λ = 100, rate 2^-4, k = 4, repaired-list soundness.

| configuration | gates | vs base | round-0 queries | proof (Fr) |
|---|---:|---:|---:|---:|
| per-query paths only (the base circuit) | 235,887 | 100 % | 26 | 5,187 |
| + Merkle cap | 196,539 | 83 % | 26 | 4,680 |
| + Poseidon2 grinding, 20 bits | 193,433 | 82 % | 21 | 4,219 |
| + k₀ = 1 (narrow first fold) | 215,618 | 91 % | 26 | 3,969 |
| cap + grinding | 163,941 | 69 % | 21 | 3,840 |
| **cap + grinding + k₀ = 1** | **144,455** | **61 %** | 21 | 2,819 |
| cap + grinding (24 bits) + k₀ = 1 | 138,546 | 59 % | 20 | 2,689 |

The bold row is the recommendation. The last row is the gate-minimal point but not the one to ship:
24 bits of Poseidon2 grinding costs the prover about fifteen times the search for the last 4 % of
gates — see "Rate and grinding" below.

The levers compose almost independently, because each divides a different factor of the cost
expression in §1. `k₀ = 1` on its own even looks like a regression (91 %) — halving the round-0 leaf
width also adds a tree level, and the extra path work only turns into a win once a cap is absorbing
it.

### Where the constraints go

At the recommended configuration above (144,455 gates):

| phase | gates | share |
|---|---:|---:|
| merkle: path walk | 62,712 | 43.4 % |
| merkle: leaf hashing | 46,626 | 32.3 % |
| merkle: cap fold | 12,580 | 8.7 % |
| query index extraction | 3,263 | 2.3 % |
| coset folding | 3,237 | 2.2 % |
| batched oracle assembly | 3,234 | 2.2 % |
| final claim | 3,019 | 2.1 % |
| final oracle consistency | 2,940 | 2.0 % |
| grinding | 1,924 | 1.3 % |
| claims + batching | 1,299 | 0.9 % |
| whir sumcheck | 1,083 | 0.7 % |
| merkle: cap lookup | 780 | 0.5 % |

Everything that is not hashing is under 15 %. A verifier optimized past this point is optimizing
Poseidon2 permutations, not arithmetic.

(The phases are stamped from the builder's in-progress gate count and sum to about 98.7 % of the
total; the remainder is what finalization adds — the sorted lists backing the range constraints and
the lookup tables.)

### Rate and grinding

Lengthening the codeword buys query count directly. m = 10, λ = 100, k = 4, k₀ = 1.

| | no grinding | 20-bit grind |
|---|---:|---:|
| rate 2^-2 | 259,924 (51 q) | 216,664 (41 q) |
| rate 2^-3 | 201,926 (34 q) | 166,638 (27 q) |
| rate 2^-4 | 173,729 (26 q) | 144,455 (21 q) |
| rate 2^-6 | 138,021 (17 q) | **116,183 (14 q)** |

Rate is a prover-side cost — 2^-6 quadruples the round-0 codeword against 2^-4 — and cheap for a
rollup that proves once and aggregates many. **Grinding is not as cheap as it looks**, and it is the
one place where the recursion profile makes the prover materially worse.

Measured trial rates on this machine (14 threads, `DISABLE_ASM=1`, so a build with field assembly
would be several times faster):

| | trials/s | 16-bit grind | 20-bit | 24-bit |
|---|---:|---:|---:|---:|
| Blake3 (bb's default) | 31.8 M | 0.002 s | 0.03 s | 0.53 s |
| Poseidon2 (recursion profile) | 0.70 M | 0.09 s | 1.5 s | 24 s |

Poseidon2 is **45x slower to grind** than Blake3 — the same algebraic structure that makes it cheap to
*verify* in a circuit makes it expensive to *search* natively. These schedules grind three times per
proof (two fold rounds and the final round), so 24 bits costs the prover about 72 s against a WHIR
prover that otherwise runs in seconds. **20 bits is the practical setting**: about 4.5 s of grinding
for 17 % off the circuit, where the last 4 % that 24 bits buys costs fifteen times the search. The
tables here quote 24 bits because they rank circuits, not provers.

(`poseidon2_pow_is_valid` spells the two-element sponge out as a single permutation rather than
calling `Poseidon2::hash`, which allocates a vector per trial; that alone was worth 4.9x on the
search.)

### Column count is the other big dial

m = 12, λ = 100, rate 2^-4, k₀ = 1, 20-bit grind:

| shape | gates | proof (Fr) |
|---|---:|---:|
| ProveKit: 1 column, 1 group | **92,094** | 1,211 |
| 4 columns, 1 group | 95,607 | 1,340 |
| 8 columns, 2 groups | 120,075 | 1,756 |
| transparent Honk: 29 columns, 4 groups | 183,034 | 3,149 |

This is the strongest argument for targeting ProveKit-shaped proofs: Spartan reduces the statement to
**one** committed vector, where transparent Honk opens ~30 columns in four commitment phases. The
group count matters as much as the column count — four groups means four independent authentication
paths per query, and paths are the largest single phase.

### Growth in the inner circuit is mild

λ = 100, rate 2^-4, k = 4, k₀ = 1, 20-bit grind, transparent-Honk layout:

| inner circuit | 2^8 | 2^10 | 2^12 | 2^14 |
|---|---:|---:|---:|---:|
| gates | 135,632 | 144,455 | 183,034 | 195,848 |

Doubling the inner circuit adds one Merkle level per query and one sumcheck round — logarithmic, as
the protocol promises. A verifier for a 2^20 inner circuit is not far above these numbers.

### The outer proof

Verifications compose by construction — each proof is checked independently against its own
transcript, so a bad proof anywhere makes the whole circuit unsatisfiable. Cost is additive with a
little amortization (the range-constraint tables and the lookup machinery are shared):

| proofs in one circuit | gates | gates/proof | smallest outer circuit |
|---:|---:|---:|---|
| 1 | 92,094 | 92,094 | 2^17 |
| 2 | 180,817 | 90,408 | 2^18 |
| 4 | 358,268 | 89,567 | 2^19 |

At ~90k gates per proof a 2^21 outer circuit holds about **23 ProveKit-shape WHIR proofs**, or about
11 transparent-Honk ones.

`whir_recursion_bench outer N` then proves that circuit with UltraHonk + KZG. The output is a
13,120-byte pairing-based proof whose size does not depend on N.

| proofs | outer circuit | key generation | prove | verify | proof |
|---:|---|---:|---:|---:|---:|
| 1 | 2^17 rows | 102 ms | 266 ms | 3 ms | 13,120 B |
| 4 | 2^19 rows | 303 ms | 923 ms | 5 ms | 13,120 B |
| 8 | 2^20 rows | 746 ms | 2,496 ms | 6 ms | 13,120 B |

About 310 ms of outer proving per aggregated WHIR proof, and neither the proof size nor the
verification cost moves with N.

On Ethereum that proof is checked by the `--optimized` Solidity verifier for **630,067 execution gas
/ 781,543 as a transaction** — measured on a 2^19 proof in `commitment_schemes/recursion/README.md`
§6, and nearly flat in circuit size (each extra `log n` adds one Gemini fold commitment to the MSM,
about 6.3k gas). So the whole cost of settling N WHIR proofs on L1 is one ~780k-gas transaction plus
the outer prover's time, against a WHIR proof that could not be posted at any price.

## 5. What is left

Three levers are priced but not implemented.

**Ternary Merkle trees.** A binary node hashes two children into a rate-3 permutation and wastes a
third of the rate; a ternary node uses all of it, cutting the depth by `log₂3 = 1.585`. The index
arithmetic is the price: base-3 digits need a booleanity-plus-sum argument and a three-way mux per
level, about 10 extra gates against 2 for the binary swap. Modelling that against the measured
per-level cost gives roughly **18 % off the path-walk and cap phases, so 9–10 % overall**. The
domain-point exponentiation gets *cheaper*, since `ω^idx = Π_j (ω^{3^j})^{d_j}` and a three-way
constant select is a linear term.

**A hard cap.** The cap is currently folded to the root in-circuit at `2^c − 1` permutations per
tree. Absorbing the cap into the transcript as the commitment instead costs `⌈2^c/3⌉` — for `c = 5`,
11 permutations against 31. That removes most of the 9.1 % the cap fold costs, at the price of making
the commitment a vector of digests rather than one root, which ripples into `WhirGroupData` and the
transparent-Honk verification key.

**Per-group column stacking.** Round 0's leaf width is `columns × 2^k₀`, and 21 of the 29 columns are
the *precomputed* ones — the verification key, fixed before any proof exists. Stacking just that
group into one tall codeword makes its leaf two values instead of forty-two, and the encoding is paid
once at key-generation time rather than per proof. The tree gets deeper (the stacked codeword is 21×
longer), so the win is partial. Group 0 is about a third of the circuit — roughly 70 % of the leaf
hashing and a quarter of the path walk — and the model puts the saving at about **36 % of that, so
~12 % overall**. `WhirConfig::stack_columns` is currently all-or-nothing and costs ~7× prover time
when applied to witness columns, so this needs per-group stacking to be worth it.

**Skyscraper** is the wrong in-circuit hash despite being the best native one (it is the fastest of
bb's four WHIR hashers to prove with). Its `bar` needs a byte decomposition and S-box lookups per
field element — roughly 2–4× Poseidon2's ~23 gates per absorbed element — and Ultra gives Poseidon2
dedicated custom gates that nothing else gets. The native ranking and the in-circuit ranking are
opposites, and only the in-circuit one matters here. **Blake3s** is the same story, measured: see
below.

### Blake3s, measured rather than asserted

`StdlibBlake3sHasher` is a second in-circuit hasher, bit-identical to `bb::whir::Blake3sMerkleHasher`
(a test pins both leaf and node against it), so the same verifier runs on Blake3s-committed proofs
and the choice can be priced. The unit costs on this machine:

| | gates |
|---|---:|
| Poseidon2 permutation — a Merkle node, or three absorbed values | 75 |
| Blake3s of 65 bytes — a Merkle node | 4,330 |
| Blake3s of 769 bytes — one 24-value leaf chunk | 28,575 |
| `field_t` → its canonical 32 bytes | 2,865 |

Two separate penalties. A Blake3s node compression is **58× a Poseidon2 permutation**, because in a
BN254 circuit it is 32-bit word arithmetic rather than field arithmetic. And every field element it
absorbs must first be decomposed into 32 bytes with a canonicity check, which alone costs 38× a whole
Poseidon2 permutation — so a leaf of `n` values pays `2,865n` before any hashing happens.

End to end, the same statement verified both ways (`whir_recursion_bench hashers`):

| configuration | Poseidon2 | Blake3s | ratio | Poseidon2 proof | Blake3s proof |
|---|---:|---:|---:|---:|---:|
| m = 8, λ = 32 | 31,482 | 719,321 | 22.8× | 12,608 B | 16,256 B |
| m = 8, λ = 64 | 51,696 | 1,291,389 | 25.0× | 19,776 B | 25,984 B |
| m = 10, λ = 32 | 30,032 | 1,104,605 | 36.8× | 14,656 B | 20,416 B |

The ratio *grows* with the tree depth, because the deeper the tree the more of the cost is node
hashing, where Blake3s is worst (58× against the leaf side's ~30×). A Blake3s digest is also two
field elements to Poseidon2's one, so its proof carries about twice the Merkle-path bytes.

The table stops at deliberately small parameters because that is where it stops being measurable:
extrapolating the m = 10 ratio, the transparent-Honk verifier of §4 would be roughly 4–5 M gates on
Blake3s against 144,455 on Poseidon2, and the whole-proof verifier of §6 would be larger still — past
what one outer proof holds on this machine. `TransparentHonkRecursiveVerifier` is therefore fixed to
Poseidon2 rather than templated on the hasher: the opening comparison above already isolates the
choice, and the Honk shell around it hashes with the Poseidon2 transcript either way, so a
whole-proof comparison would measure the same ratio diluted.

So the hash choice inverts between the two sides of the pipeline: bb's WHIR defaults to Blake3s
because it is the fastest to *prove* with, and a recursive verifier must use Poseidon2 because it is
20–40× cheaper to *verify*. If proofs are going to be aggregated, that decision belongs to the
verifier.

### On turning the verifier into lookups

The natural thought, given that 85 % of the circuit is hashing, is to replace the Merkle work with
table queries — cq being the obvious candidate, since its prover is independent of table size. Three
things get in the way, and they are worth stating because they rule out the whole family rather than
one scheme.

*A Merkle path is not a membership claim.* A lookup argument proves "every element of this vector
appears in that table". An authentication path is a chain of hashes whose structure depends on the
query index, and the intermediate digests are not members of any fixed set — they are computed. There
is no table to look into.

*The codeword is committed, not public.* The one statement here that does look like a lookup is
"value `v` sits at index `i` of the codeword". If the codeword were a public table, `t` queries would
indeed be `t` lookups. But the codeword is exactly the object the commitment hides, and committing it
so that a lookup argument can address it means a polynomial commitment — which is where we came in.

*A pairing-based verifier is the expensive part in-circuit, not the cheap part.* cq's verifier costs
pairings. Inside a BN254 circuit those are non-native field arithmetic — thousands of gates per
operation, well past the Merkle hashing they would replace. cq's table-size independence is a
*prover*-side property; an in-circuit verifier inherits none of it. The same objection applies to any
scheme whose verifier lives on the curve rather than in the field.

What lookups *do* buy here is smaller and real: bb's plookup range constraints are what make the
index-canonicity argument and the grinding quotient affordable — together 3.6 % of the circuit for
checks that a naive bit decomposition would price at ten times that. The place a bigger table could
still pay is the domain-point exponentiation `ω^idx`: a 2^k-entry table of `ω` powers would replace k
multiplications with one lookup. At 2.2 % of the circuit for coset folding and exponentiation
combined, it is not where the remaining gates are.

## 6. Verifying a whole proof, not just the opening

`TransparentHonkRecursiveVerifier` is the complete thing: verification-key absorption, public
inputs, the three witness commitment rounds and their challenges, the full sumcheck over the
flavor's relation set, the virtual-column checks, and then the WHIR opening. Satisfying its circuit
is verifying the inner proof.

It reuses bb's own `SumcheckVerifier` rather than restating it. `UltraProveKitRecursiveFlavor_`
re-instantiates `Flavor::Relations` over `stdlib::field_t`, and because the flavor is registered
with the `IsRecursiveFlavor` concept, the sumcheck round checks become `assert_equal` instead of
native comparisons. The relation set the circuit evaluates is therefore the same source the native
verifier uses, with no hand-written copy to drift.

**The Honk shell is nearly free.** For a 2^13 inner circuit at λ = 64 the whole verifier is 117,314
gates, of which the shell is 4,825:

| phase | gates |
|---|---:|
| honk: transcript + public inputs | 1,007 |
| honk: sumcheck + relations (all 20 subrelations, 13 rounds) | 3,794 |
| honk: virtual columns | 24 |
| *everything else — the WHIR opening* | *112,489* |

That is the useful headline: **verifying a whole WHIR-Honk proof costs about 4 % more than verifying
its WHIR opening alone.** Sumcheck is `log n` rounds of a degree-6 univariate plus one evaluation of
the relation set — a few thousand gates — against the hundreds of thousands the opening's Merkle
work costs. Which also means the levers in §4 are the levers for the whole thing.

One thing the circuit does **not** do for you: `verify` returns the inner proof's public inputs, and
binding them is the caller's job. Left unbound, the circuit proves "a proof of this inner circuit
exists" with the public inputs free for the prover to choose. An aggregator that means "a proof of
*this statement* exists" has to forward them to its own public inputs or hash them in — which is why
they are returned rather than dropped.

### The whole pipeline, measured

Inner proofs are transparent UltraHonk with WHIR on `UltraProveKitFlavor` (λ = 100, rate 2^-4, k = 4,
k₀ = 1, Poseidon2 Merkle); the outer proof is ordinary UltraHonk + KZG. Same machine as §4.
`whir_recursion_bench full` regenerates it.

With a 20-bit grind:

| inner circuit | inner prove | inner verify | inner proof | verifier gates | outer circuit | outer prove | outer verify | outer proof |
|---|---:|---:|---:|---:|---|---:|---:|---:|
| 2^13 | 4,921 ms | 23 ms | 112,832 B | 187,655 | 2^18 | 497 ms | 3 ms | 13,120 B |
| 2^15 | 4,530 ms | 27 ms | 122,880 B | 218,860 | 2^18 | 579 ms | 4 ms | 13,120 B |
| 2^17 | 17,014 ms | 28 ms | 136,416 B | 242,490 | 2^18 | 618 ms | 3 ms | 13,120 B |

Without grinding:

| inner circuit | inner prove | inner verify | inner proof | verifier gates | outer circuit | outer prove | outer verify | outer proof |
|---|---:|---:|---:|---:|---|---:|---:|---:|
| 2^13 | 603 ms | 28 ms | 136,448 B | 225,889 | 2^18 | 596 ms | 3 ms | 13,120 B |
| 2^15 | 2,385 ms | 33 ms | 148,416 B | 263,518 | 2^19 | 732 ms | 4 ms | 13,120 B |
| 2^17 | 9,521 ms | 33 ms | 164,704 B | 291,595 | 2^19 | 831 ms | 4 ms | 13,120 B |

Three things to read off it. The verifier grows **logarithmically** in the inner circuit — 2^13 to
2^17 is sixteen times the circuit for 1.29x the gates. A 112 KB inner proof becomes a 13,120-byte
outer one whatever the inner size, and the outer verification is 3 ms rather than 23–28 ms. And
grinding is worth about 17 % of the verifier circuit but is the dominant term in inner proving at
these sizes — 4.3 s of the 4.9 s at 2^13 — which is the §4 trade seen from the prover's side.

## 7. What full ProveKit recursion would cost

The obvious next target is a whole ProveKit proof rather than a whole WHIR-Honk one. It does not
work, and the reason is Spartan rather than WHIR.

**ProveKit's verifier is not succinct in the circuit it proves.** `WhirR1CSVerifier::verify` takes
the R1CS as an argument and calls `multiply_transposed_by_eq_alpha`, which builds an eq-table over
the whole constraint hypercube and multiplies all three matrices by it — `2^m + nnz(A) + nnz(B) +
nnz(C)` field operations against the raw matrices. That is why their verifying key for the passport
circuit is **9.1 MB** and their native verify is ~371 ms where bb's WHIR-Honk verify is 23–33 ms.

This is not one implementation's choice. ProveKit's *own* gnark recursive verifier does the same
thing in-circuit: `recursive-verifier/app/circuit/matrix_evaluation.go` loops over every nonzero with
`api.Add(ans, api.Mul(cell.value, api.Mul(rowEval[cell.row], colEval[cell.column])))`, and its CLI
takes "the R1CS JSON file describing the constraint system of the inner circuit" as a required input.

Priced in barretenberg (`whir_recursion_bench provekit-cost`), that inner loop costs **2.0 gates per
nonzero** and an eq-table **2.0 gates per hypercube entry** — the matrix structure is public, so only
the product of the two evaluation tables is a gate.

The row table is indexed by constraint, so it is `2^m` entries and there is no folding variant for
it: at 2^19 constraints that is **~1.05M gates before a single matrix entry is touched**. (The
column table can be small — `evaluateFoldedR1CSMatrixExtension` wraps column indices modulo a
blinding-sized mask — so the column side is not the problem.) On top of that the matrices themselves
cost 2 gates per nonzero, and ProveKit's verifying key for the passport circuit is 9.1 MB of them.

So a full ProveKit recursive verifier is Θ(inner circuit size): larger than the circuit it verifies,
and tied to one specific circuit. The WHIR opening inside it — the part this directory implements —
is 92,094 gates at m = 12 (§4) and grows logarithmically, so it is one to two orders of magnitude
below the Spartan term at any realistic size. **WHIR is not what makes ProveKit hard to recurse.**
The fix is on their side: a sparse commitment to the matrices so the verifier stops reading them.
ProveKit has the machinery (`provekit/spark`, a Spartan-style SPARK compiler), but it is wired into
the prover only — `grep -rl spark provekit/verifier/src` is empty. Until it reaches the verifier,
aggregating ProveKit proofs means aggregating the WHIR opening and re-proving the R1CS, which is not
aggregation.

## 8. Scope, and running it

The verifier covers the shape transparent Honk and ProveKit both use: interleaved columns (no
stacking), no zero knowledge, every claim at the single sumcheck point. `zk`, `stack_columns` and
multi-point claim plans are rejected with an assertion rather than silently mis-verified.

`WhirRecursiveVerifier` verifies the opening alone; `TransparentHonkRecursiveVerifier` (§6) verifies
a whole proof. The §4 tables measure the former, so they are not directly comparable to the Noir
verifier's 501,884 in `commitment_schemes/recursion/`, which bundles the Honk shell — §6 gives the
comparable figure.

```bash
cd barretenberg/cpp
cmake --build build --target commitment_schemes_tests whir_recursion_bench
./build/bin/commitment_schemes_tests --gtest_filter='WhirRecursiveVerifierTest.*:TransparentHonkRecursiveVerifierTest.*:WhirRecursionProfileTest.*'
./build/bin/whir_recursion_bench                 # every table in §4
./build/bin/whir_recursion_bench full            # the whole pipeline, §6
./build/bin/whir_recursion_bench provekit-cost   # the ProveKit term, §7
./build/bin/whir_recursion_bench outer 8         # aggregate 8 openings and prove the result
```

The tests cover completeness over fresh randomness, aggregation of several proofs in one circuit, and
soundness: corrupting any sampled element of the proof stream, any claimed evaluation, any coordinate
of the opening point, or any commitment root makes the circuit unsatisfiable. For the whole-proof
verifier the same tamper sweep covers the sumcheck univariates and the claimed evaluations, and the
precomputed-column root — a circuit constant — is what pins the verifier to one inner circuit.
