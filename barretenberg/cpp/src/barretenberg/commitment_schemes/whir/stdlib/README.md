# Recursive verification of WHIR in barretenberg

A WHIR opening verifier written directly in barretenberg's circuit stdlib, together with the
protocol-level settings that make it cheap and the measurements that justify each one. A WHIR proof
is hundreds of kilobytes and will never go on Ethereum; the way to use WHIR in a rollup is to verify
it inside an UltraHonk circuit and settle the 13 KB pairing-based proof instead.

| | |
|---|---|
| `whir_recursive_verifier.hpp` | the in-circuit verifier, the mirror of `bb::whir::WhirVerifier::verify` |
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
runs it 2^b times — and a Blake3 compression is tens of thousands of constraints in-circuit. There is
one grind per round, so seven of them would cost more than the ~30,000 gates the removed queries were
worth. The Poseidon2 form accepts when `Poseidon2(seed, nonce)` is divisible by `2^b`, which is one
permutation plus two range constraints on the quotient: 1.3 % of the circuit for a 15 % cut in
queries.

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
| per-query paths only (the base circuit) | 235,639 | 100 % | 26 | 5,187 |
| + Merkle cap | 196,291 | 83 % | 26 | 4,680 |
| + Poseidon2 grinding, 20 bits | 193,151 | 82 % | 21 | 4,219 |
| + k₀ = 1 (narrow first fold) | 215,073 | 91 % | 26 | 3,969 |
| cap + grinding | 163,659 | 69 % | 21 | 3,840 |
| cap + grinding + k₀ = 1 | 143,892 | 61 % | 21 | 2,819 |
| **cap + grinding (24 bits) + k₀ = 1** | **137,967** | **59 %** | 20 | 2,689 |

The levers compose almost independently, because each divides a different factor of the cost
expression in §1.

### Where the constraints go

At the best configuration above (137,967 gates):

| phase | gates | share |
|---|---:|---:|
| merkle: path walk | 59,592 | 43.2 % |
| merkle: leaf hashing | 44,212 | 32.0 % |
| merkle: cap fold | 12,580 | 9.1 % |
| query index extraction | 3,117 | 2.3 % |
| coset folding | 3,090 | 2.2 % |
| batched oracle assembly | 3,080 | 2.2 % |
| final claim | 2,889 | 2.1 % |
| final oracle consistency | 2,695 | 2.0 % |
| grinding | 1,830 | 1.3 % |
| claims + batching | 1,299 | 0.9 % |
| whir sumcheck | 1,071 | 0.8 % |
| merkle: cap lookup | 751 | 0.5 % |

Everything that is not hashing is under 15 %. A verifier optimized past this point is optimizing
Poseidon2 permutations, not arithmetic.

### Rate and grinding

Lengthening the codeword buys query count directly. m = 10, λ = 100, k = 4, k₀ = 1.

| | no grinding | 24-bit grind |
|---|---:|---:|
| rate 2^-2 | 259,610 (51 q) | 206,483 (39 q) |
| rate 2^-3 | 202,110 (34 q) | 159,756 (26 q) |
| rate 2^-4 | 173,184 (26 q) | 137,967 (20 q) |
| rate 2^-6 | 138,131 (17 q) | **110,440 (13 q)** |

Rate is a prover-side cost (2^-6 quadruples the round-0 codeword against 2^-4) and grinding is a
one-time 2^24 hash search, so both are cheap for a rollup that proves once and aggregates many.

### Column count is the other big dial

m = 12, λ = 100, rate 2^-4, k₀ = 1, 24-bit grind:

| shape | gates | proof (Fr) |
|---|---:|---:|
| ProveKit: 1 column, 1 group | **87,530** | 1,155 |
| 4 columns, 1 group | 90,880 | 1,278 |
| 8 columns, 2 groups | 114,302 | 1,676 |
| transparent Honk: 29 columns, 4 groups | 174,523 | 3,007 |

This is the strongest argument for targeting ProveKit-shaped proofs: Spartan reduces the statement to
**one** committed vector, where transparent Honk opens ~30 columns in four commitment phases. The
group count matters as much as the column count — four groups means four independent authentication
paths per query, and paths are the largest single phase.

### Growth in the inner circuit is mild

λ = 100, rate 2^-4, k = 4, k₀ = 1, 24-bit grind, transparent-Honk layout:

| inner circuit | 2^8 | 2^10 | 2^12 | 2^14 |
|---|---:|---:|---:|---:|
| gates | 129,679 | 137,967 | 174,523 | 186,222 |

Doubling the inner circuit adds one Merkle level per query and one sumcheck round — logarithmic, as
the protocol promises. A verifier for a 2^20 inner circuit is not far above these numbers.

### The outer proof

Verifications compose by construction — each proof is checked independently against its own
transcript, so a bad proof anywhere makes the whole circuit unsatisfiable. Cost is additive with a
little amortization (the range-constraint tables and the lookup machinery are shared):

| proofs in one circuit | gates | gates/proof | smallest outer circuit |
|---:|---:|---:|---|
| 1 | 87,530 | 87,530 | 2^17 |
| 2 | 171,583 | 85,791 | 2^18 |
| 4 | 339,689 | 84,922 | 2^19 |

At ~85k gates per proof a 2^21 outer circuit holds about **24 ProveKit-shape WHIR proofs**, or about
12 transparent-Honk ones.

`whir_recursion_bench outer N` then proves that circuit with UltraHonk + KZG. The output is a
13,120-byte pairing-based proof whose size does not depend on N.

| proofs | outer circuit | key generation | prove | verify | proof |
|---:|---|---:|---:|---:|---:|
| 1 | 2^17 rows | 183 ms | 556 ms | 5 ms | 13,120 B |
| 4 | 2^19 rows | 1,000 ms | 2,864 ms | 7 ms | 13,120 B |
| 8 | 2^20 rows | 1,799 ms | 5,166 ms | 5 ms | 13,120 B |

Roughly 640 ms of outer proving per aggregated WHIR proof, and the proof and verification cost do not
move.

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
opposites, and only the in-circuit one matters here.

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

## 6. Scope

The verifier covers the shape transparent Honk and ProveKit both use: interleaved columns (no
stacking), no zero knowledge, every claim at the single sumcheck point. `zk`, `stack_columns` and
multi-point claim plans are rejected with an assertion rather than silently mis-verified. The Honk
shell around the opening — sumcheck and the flavor's subrelations — is not included here; the Noir
verifier in `commitment_schemes/recursion/` has it, and its cost is a constant independent of the PCS.
