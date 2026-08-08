# Recursive verification of WHIR-Honk and Ligero-Honk in Noir

Two Noir circuits that verify a transparent UltraHonk proof — one for the WHIR polynomial commitment
scheme, one for Ligero — together with the tooling to prove those circuits with UltraHonk (KZG and
Vela) and to price the result on Ethereum.

Everything here verifies real proofs: the Noir circuits are driven by proofs `WhirHonk` and
`LigeroHonk` actually produced, and every number below was measured, not modelled.

## 1. What is here

| Path | What it is |
|---|---|
| `noir/honk_common/` | Noir library: Poseidon2, the Fiat-Shamir transcript, the Honk shell (sumcheck + all twenty subrelations), Merkle authentication |
| `noir/whir_verifier/` | the WHIR-Honk recursive verifier circuit |
| `noir/ligero_verifier/` | the Ligero-Honk recursive verifier circuit |
| `pcs_noir_export.bench.cpp` | proves the inner circuit and writes each verifier's `params.nr` and `Prover.toml` |
| `vela_evm_export.bench.cpp` | proves an ACIR circuit with UltraHonk + Vela under a Keccak transcript and exports its opening argument |
| `fflonk/` (in `commitment_schemes/`) | the fflonk opening as a PCS backend, with `fflonk_evm_export.bench.cpp` for its Keccak-transcript export |
| `sol/` | Foundry harness: gas for the generated UltraHonk verifier (default and `--optimized`), the Vela opening in both plain and assembly form, and an IPA verifier's size-n MSM |
| `driver.sh` | `sizes`, `outer`, `gas`, `ipa-gas` — regenerates every measurement below |

The inner circuit is a small UltraHonk circuit that exercises every relation the flavor keeps —
arithmetic gates, public inputs, a lookup, and a ROM table — so the recursive verifier really
evaluates all twenty subrelations rather than a degenerate subset. It runs on `UltraProveKitFlavor`
(24 precomputed + 8 witness + 5 shifted entities, `BATCHED_RELATION_PARTIAL_LENGTH = 6 + 1`).

## 2. How the Noir verifiers are structured

Noir needs every array length and loop bound at compile time, and a WHIR or Ligero schedule fixes
those only once its parameters are chosen. So `params.nr` is **generated**: the exporter writes each
proof family's structural constants (round query counts, tree depths, leaf widths, claim layout,
domain roots) and the circuit is compiled against them. One compiled circuit verifies any proof of
that family; a parameter sweep is a regeneration, which is exactly what `driver.sh sizes` does.

Three things make the port faithful rather than approximate.

**The transcript is bit-identical to `bb::NativeTranscript`.** A challenge is
`Poseidon2(previous_challenge ‖ round_buffer)`, and barretenberg's Poseidon2 sponge (rate 3,
capacity 1, IV = `len << 64`) is reproduced in `poseidon.nr` as a chunked absorb so every index stays
a compile-time constant. It reproduces bb's own `HashConsistencyCheck` test vector exactly.

**Values the prover sends "unhashed" are witnesses, not a parsed byte stream.** WHIR and Ligero send
Merkle openings without absorbing them — they are bound to the transcript by a root that was absorbed
earlier. Feeding them as separately-sized witness arrays instead of one flat proof array keeps every
array index a compile-time constant, which is what avoids paying a RAM lookup for each of the tens of
thousands of opened values.

**The batched Merkle walk uses a checked hint.** `MerkleTree::verify_batch` climbs the tree over the
*distinct ascending* leaf set the query indices induce, which is data-dependent — in-circuit that
would mean sorting witnesses. Instead the prover supplies the leaf set and each query's slot in it,
and the circuit checks the set is strictly ascending and that every query index sits at its claimed
slot. That is all soundness needs: each query's coset is still authenticated under the root at its
own leaf index, and a longer leaf set would only add authentication work.

Every array is sized for the worst case (all queries distinct) with the real count carried as a
witness, so a batch with repeated indices simply leaves the tail of each level inactive.

## 3. Circuit sizes

Inner circuit 2^10 rows, λ = 100, Poseidon2 Merkle hashing (the only sensible in-circuit choice — a
Blake3 or SHA-256 leaf costs thousands of gates per compression where Poseidon2 costs one
permutation).

**WHIR-Honk.** `queries` is the round-0 count, the term that dominates: each one reveals
`2^k0` values of every committed column.

| configuration | round-0 queries | arity 2^k0 | tree depth | inner proof | ACIR opcodes | UltraHonk gates |
|---|---:|---:|---:|---:|---:|---:|
| `whir-r4-repaired-k4-k0_2` | 26 | 4 | 12 | 137 KiB | 202,196 | 501,884 |
| `whir-r4-provable-k4-k0_2` | 52 | 4 | 12 | 264 KiB | 405,454 | 993,352 |
| `whir-r2-repaired-k4-k0_2` | 51 | 4 | 10 | 236 KiB | 330,361 | 833,498 |
| `whir-r6-repaired-k4-k0_2` | 17 | 4 | 14 | 105 KiB | 153,080 | 376,713 |
| `whir-r4-repaired-k4-k0_1` | 26 | 2 | 13 | 102 KiB | 213,876 | 487,688 |
| `whir-m12-r4-repaired` | 26 | 4 | 14 | 150 KiB | 242,181 | 576,910 |
| `whir-m14-r4-repaired` | 26 | 4 | 16 | 167 KiB | 289,964 | 694,326 |

**Ligero-Honk.** No rounds: `t` queries, each opening one column of every commitment group.

| configuration | queries | matrix rows x cols | index bits | inner proof | ACIR opcodes | UltraHonk gates |
|---|---:|---:|---:|---:|---:|---:|
| `ligero-r4` | 25 | 2 x 512 | 13 | 137 KiB | 66,517 | 345,768 |
| `ligero-r2` | 50 | 2 x 512 | 11 | 211 KiB | 125,617 | 578,941 |
| `ligero-m12-r4` | 25 | 4 x 1024 | 14 | 234 KiB | 110,966 | 479,109 |

Reading the sweep: WHIR's cost is dominated by the round-0 openings, which reveal `2^k0` values of
every one of the ~30 committed columns per query. The two levers that matter are therefore the query
count (rate and soundness regime) and `k0` — not `k`, which only governs the later single-column
rounds. Ligero has no rounds at all: its cost is `t` queries × one column of every group, so it grows
as √N in the inner circuit and shrinks directly with the query count.

## 4. Proving the recursive verifiers

Both verifier circuits land at 2^19 UltraHonk rows. Medians of 5 fresh-process runs on an M4 Pro
(`DISABLE_ASM=1` arm64 build); absolute times on this machine are not comparable across sessions, but
the KZG/Vela ratio within one run is.

| verifier circuit | outer PCS | prove | verify | proof | peak RSS |
|---|---|---:|---:|---:|---:|
| Ligero-Honk | UltraHonk + KZG | 1016 ms | 5.89 ms | 13,408 B | 832 MiB |
| Ligero-Honk | UltraHonk + Vela | 1348 ms | 4.33 ms | 8,320 B | 1,742 MiB |
| WHIR-Honk | UltraHonk + KZG | 1209 ms | 5.41 ms | 13,408 B | 1,230 MiB |
| WHIR-Honk | UltraHonk + Vela | 1385 ms | 5.59 ms | 8,320 B | 2,366 MiB |

Vela's proof is 8,320 bytes against KZG's 13,408 — 38 % smaller — because it replaces Gemini's
`log n` fold commitments and evaluations with three group elements and five field elements. It pays
for that with a slower prover (the Laurent-product expansion is `μ` shift-and-add passes over the
whole coefficient table) and roughly twice the peak memory.

## 5. Ethereum gas

The generated Solidity verifier (`bb write_solidity_verifier -t evm-no-zk`) gives a real end-to-end
number for UltraHonk + Shplemini + KZG. For Vela, `VelaOpening.sol` implements the opening step and is
driven by a real Keccak-transcript Vela proof of the same circuit, so both openings are measured, not
modelled. The transcript, sumcheck and relation evaluation before the opening are identical for the
two schemes — same flavor, same `log n` — so the totals compose.

Measured on the Ligero-Honk verifier circuit's proof (2^19 rows, 9 public inputs). The WHIR-Honk
circuit is the same size with the same public-input count, and the verifier's work does not depend on
the verification key's contents, so its figures are the same to within noise.

| | Shplemini + KZG | Vela |
|---|---:|---:|
| shared Honk work (transcript, sumcheck, relations) | 895,736 | 895,736 |
| opening argument | 729,504 | **499,685** |
| execution gas, total | 1,625,240 | **1,395,421** |
| proof calldata | 8,832 B / 130,476 gas | 7,616 B / **110,720 gas** |
| transaction gas (21,000 + calldata + execution) | 1,776,716 | **1,527,141** |

The full generated verifier, called end to end rather than decomposed, costs 1,635,967 execution gas —
0.7 % above the sum of the two components, which is the external-call and calldata-copy overhead the
decomposition leaves out.

Vela therefore saves **~230k gas on the opening and ~20k on calldata**, about 14 % of the whole
verification.

Two caveats worth stating with these numbers. The generated Shplemini verifier does its MSM in
hand-written assembly; `VelaOpening.sol` calls the `ecAdd`/`ecMul` precompiles through ordinary
Solidity, allocating memory per call. The Vela figure is therefore an upper bound, and the gap is if
anything understated. And the shared component is measured on the production UltraHonk transcript,
which differs from `TransparentHonk`'s by a handful of extra keccaks over the verification-key
metadata — noise at this scale.

The structural reason Vela wins: its opening is one MSM over the claim commitments (31 points here)
plus three further scalar multiplications and one pairing, where Shplemini needs
`NUM_ENTITIES + log n + 2 ≈ 62` points for the same claim set.

## 6. The optimized verifier, and what Aztec actually pays on L1

`bb write_solidity_verifier` has two backends. Section 5 measured the default one. Aztec's L1 rollup
verifier is the other: `noir-projects/fnd/noir-protocol-circuits/bootstrap.sh` generates it with
`--optimized` from the `rollup_root` circuit's key, and `l1-contracts/bootstrap.sh` copies the result
to `generated/HonkVerifier.sol`. On the same 2^19 proof:

| verifier | execution gas | + calldata | transaction gas |
|---|---:|---:|---:|
| default (`Verifier.sol`) | 1,635,967 | 130,476 | 1,787,443 |
| `--optimized` | **630,067** | 130,476 | **781,543** |

That is the ~800k figure: 2.6x better than the default verifier, and it is the number that matters
because it is the one deployed.

**What the published Aztec benchmark does and does not include.** `l1-contracts/gas_benchmark.md`
reports `submitEpochRootProof` at 697,655 average (896,101 with validators), but
`test/builder/RollupBuilder.sol` deploys `MockVerifier` — whose `verify` returns `true` immediately.
Those figures are rollup bookkeeping only. A production epoch proof pays them **plus** the optimized
Honk verifier, so the real cost is on the order of 1.3M gas without validators and 1.5M with them.
The verifier's own cost is nearly flat in circuit size — each extra `log n` adds one Gemini fold
commitment to the MSM (~6.3k gas) and one sumcheck round of field work — so 2^19 is representative.

### Which optimizations these are

From `barretenberg/sol/src/honk/optimised/honk-optimized.sol.template`:

1. **The whole verifier is one `assembly` block over a fixed memory map.** Named offset constants
   replace structs, so nothing is allocated, ABI-encoded or copied between steps.
2. **Precompiles are called straight into scratch.** `ACCUMULATOR` at 0x00, the staged point at
   `G1_LOCATION`, the scalar at `SCALAR_LOCATION` — laid out so `ecMul`'s output region *is* `ecAdd`'s
   second operand. An MSM step is two `staticcall`s and an `mcopy`, with no memory churn.
3. **One batched inversion for the entire protocol.** Every barycentric denominator, the
   public-input-delta denominator and all the Shplemini denominators accumulate into a running
   product, one `modexp` inverts it, and a backward pass extracts each inverse — one precompile call
   instead of about forty.
4. **The proof is `calldatacopy`'d once** into contiguous memory and hashed in place with
   `keccak256(ptr, len)`; the transcript never rebuilds a buffer.
5. **Everything is unrolled** — `log n` sumcheck rounds, the powers of the evaluation challenge, the
   inverse collection — with per-index memory constants generated into the template.
6. **Constants are literals**: barycentric Lagrange denominators, `NEG_HALF_MODULO_P`, and the
   verification key itself, injected at generation time.
7. **Scratch space is aliased across phases** that provably do not overlap in time.

### The same optimizations applied to Vela

`VelaOpeningOpt.sol` is the same protocol as `VelaOpening.sol` with techniques 2, 3 and 4 applied —
they are the ones that bite here. Techniques 1, 5 and 6 target the sumcheck and relation evaluation,
which the opening step does not contain. Vela's inversion batching is a special case worth naming: it
needs `1/z` and `1/(z - 1/z)`, and both fall out of inverting `z(z^2-1)` once, since
`t = 1/(z(z^2-1))` gives `1/z = t(z^2-1)` and `1/(z - 1/z) = t z^2`.

Measured with the same split — elliptic-curve work separated from the field work that precedes it —
on the same proof:

| | Shplemini + KZG | Vela |
|---|---:|---:|
| shared Honk work (transcript, sumcheck, relations, batching scalars) | 155,803 | 155,803 |
| opening, field work | *(counted in shared)* | 70,973 |
| opening, elliptic-curve work | 474,433 | **346,330** |
| full verification | 630,236 | **573,106** |
| proof calldata | 8,832 B / 130,476 | 7,616 B / **110,720** |
| transaction gas | 781,712 | **704,826** |

Vela's opening does 27 % less curve work: 35 scalar multiplications (31 claim commitments, `C_h`,
`C_q`, `[1]`, `pi_L`) against Shplemini's ~50 once the Gemini fold commitments and the shifted-column
reuse are counted. Naive Solidity cost it 499,720; the assembly rewrite brings it to 417,303, so the
techniques transfer at about the same 1.2x the Honk verifier gets on its own curve work.

The composed Vela total is **conservative**. Shplemini's opening also does substantial *field* work —
the Gemini fold scalars, the ~40-way batch inversion, 41 rho powers — and the checkpoint leaves all of
it in the shared bucket, while Vela's 70,973 of equivalent work is charged against Vela. The true gap
is wider than the 57k the table shows.

**So: yes, the optimizations transfer, and they were worth about 83k gas on Vela.** But the headline
is that at equal optimization the schemes are closer than the unoptimized comparison in section 5
suggests — 630k against 573k execution, ~10 % of the transaction — because the generated verifier's
hand-written MSM already extracts most of what assembly can give. Vela's durable advantages are the
smaller point count and the 1,216 fewer calldata bytes, not implementation slack.

## 7. fflonk as the final wrapper

fflonk's reputation is a PlonK number: verifiers around 200k gas. Most of that comes from PlonK's
*linearisation*, which leaves the verifier an MSM over a handful of selector commitments. Honk has no
linearisation - its verifier runs `log n` sumcheck rounds and thirty-one subrelations (155,803 gas
here, measured, and untouchable by any commitment scheme) and holds ~36 column commitments. So fflonk
cannot turn an UltraHonk wrapper into a 200k verifier. What it *can* do is delete the one part that
scales with the column count, and that is worth measuring.

`commitment_schemes/fflonk/` implements it as a PCS backend. A commitment round of `t` columns is
committed as the single interleaved polynomial `g(X) = sum_i f_i(X^t) X^i`, and

    g mod (X^t - z) = sum_i f_i(z) X^i

so one opening of `g` certifies *every* column's evaluation at `z`, and the residue's coefficients are
exactly those evaluations - no root of unity, no inverse DFT. Opening at `z` and `1/z` together lets
Vela's Laurent reduction carry the multilinear claim, with `Z_r(X) = (X^t - z)(X^t - 1/z)` and the
interpolant in closed form (`B = (P - Q)/(z - 1/z)`, `A = P - zB`). A BDFG21 batch over the rounds'
point sets closes it. Batching moves out of the group and into the field: what was a 35- or 50-point
MSM becomes **eight scalar multiplications** - four round commitments, the Laurent auxiliary, the
quotient, the linearization and the generator - whatever the column count.

Measured on the same 2^19 proof, at the same optimization level, with the same split:

| | Shplemini + KZG | Vela | fflonk |
|---|---:|---:|---:|
| opening, elliptic-curve work | 474,433 | 346,330 | **225,187** |
| opening, total | *(in shared)* | 417,303 | 363,361 |
| full verification | 630,236 | 573,106 | **519,164** |
| proof calldata | 8,832 B / 130,476 | 7,616 B / 110,720 | 9,888 B / 134,592 |
| **transaction gas** | 781,712 | 704,826 | **674,756** |

And the price, natively, on the same circuit:

| | prove | verify | proof | peak RSS |
|---|---:|---:|---:|---:|
| UltraHonk + KZG | 1,016 ms | 5.89 ms | 8,832 B | 832 MiB |
| UltraHonk + Vela | 1,348 ms | 4.33 ms | 7,616 B | 1,742 MiB |
| UltraHonk + fflonk | **7,609 ms** | **2.48 ms** | 10,272 B | **5,311 MiB** |

That is the whole trade in one line: fflonk gives the cheapest L1 verification of the three and the
fastest native verifier, and pays for it with a 7.5x prover and 6x the memory. Packing a round of `t`
columns multiplies its committed degree by `t`, so the precomputed round - 26 columns, packed to 32 -
commits at 2^24 and every Shplonk polynomial pass runs at that degree. For a final wrapper, proven
once and verified forever, that is usually the right side of the trade; for anything proven often it
is not.

The packing factor is a dial rather than a constant. Splitting the precomputed round into four
sub-rounds of eight would cut the prover's degree fourfold and cost three more scalar multiplications
(~19k gas), which is a better point on the curve for most uses. The two calldata regressions are
structural, though: fflonk sends `sum_r t_r` evaluations per opening point (43 here, against Vela's
five), which is why its proof is the largest of the three even as its verifier is the cheapest.

## 6. What a transparent outer PCS would cost: IPA

A Bulletproofs-style inner-product argument needs no trusted setup, which is the one thing KZG and
Vela cannot offer. The price is that its verifier is **O(n) in the circuit size**: after folding
`log n` rounds it must certify the prover's claimed final generator `G_0 = <s, G>` with a single
multi-scalar multiplication over all `n` generators (`IpaVerifier::verify_impl` in
`commitment_schemes/pedersen_ipa/`). KZG and Vela instead close against a constant-size key, so their
verifiers touch only the ~40 claim commitments however large the circuit is.

Natively, over the same 2^19 recursive-verifier circuits (`pcs_acir_bench --pcs ipa`):

| | prove | verify | proof | peak RSS |
|---|---:|---:|---:|---:|
| Ligero-Honk verifier, UltraHonk + IPA | 2,453 ms | 174 ms | 12,832 B | 1,786 MiB |
| WHIR-Honk verifier, UltraHonk + IPA | 2,688 ms | 169 ms | 12,832 B | 2,319 MiB |

The verifier is already ~30x slower than KZG's 5.9 ms off chain. On chain that linear term is fatal.
`IpaMsmCost.sol` measures the MSM the way a Solidity verifier would have to run it, at sizes small
enough to execute, and reads the marginal cost off the slope (`driver.sh ipa-gas`):

| MSM style | gas per generator | at n = 2^19 |
|---|---:|---:|
| plain Solidity precompile calls | 14,878 | 7.80 x 10^9 |
| scratch-buffer assembly, as the generated Honk verifier writes its own | 6,558 | **3.44 x 10^9** |

The assembly figure is essentially the precompile floor — `ecMul` 6,000 plus `ecAdd` 150 — so
3.44 billion gas is a lower bound, not an estimate that better engineering can move. Against it the
`log n` folding rounds (310,770 gas measured at `log n = 19`) and the shared Honk work (895,736)
round to zero.

That is **115x the 30M block gas limit** (76x at 45M), and about 2,100x the cost of the same
verification under Vela. Two further obstacles compound it rather than substitute for it:

- the verifier needs all `n` generators. At 2^19 that is 33.5 MB — far past EIP-170's 24 KB contract
  limit, so they cannot be baked into code; as calldata they would cost 536M gas by themselves, and
  hash-to-curve derivation on chain is more expensive still;
- the cost scales with the *circuit*, not the claim count, so it grows exactly where recursion is
  meant to help.

The Halo-style escape is to defer the MSM rather than perform it — accumulate `G_0` claims and
discharge them in a later proof — which is what barretenberg does with IPA over Grumpkin in the
rollup. That keeps IPA transparent and cheap in-circuit, but the deferred claim has to be settled
somewhere, and on Ethereum that somewhere is a pairing-based scheme again. As a *direct* L1 verifier
IPA is not a candidate at these sizes; as the inner half of a cycle it is.

## 6. Reproducing

```bash
cd barretenberg/cpp
cmake --build build-arm64 --target pcs_noir_export_bench --target pcs_acir_bench --target vela_evm_export_bench
src/barretenberg/commitment_schemes/recursion/driver.sh all
```

Needs `nargo` (tested on 1.0.0-beta.26) and a bb build. The gas section additionally needs Foundry
and, once, `git clone --depth 1 https://github.com/foundry-rs/forge-std sol/lib/forge-std`;
`sol/src/Verifier.sol` is generated by `driver.sh gas`, not committed.

Each `params.nr` in the repository is the one its section-3 headline row was measured with, so
`nargo execute` in either verifier crate reproduces a passing run straight away.
