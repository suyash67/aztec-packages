# fflonk: a plonkish proof system for the last hop to Ethereum

A complete SNARK — circuit builder, preprocessing, prover, native verifier, Solidity verifier —
built for one purpose: to be the proof that actually gets verified on L1.

**237,795 gas** for a whole transaction, flat in circuit size, from a universal setup.

`PROTOCOL.md` is the specification. The C++ prover and verifier and the Solidity verifier are three
implementations of it; where they disagree, that document is what needs fixing first.

## 1. The problem this solves

A project proving with a hash-based commitment scheme — WHIR, Ligero, Brakedown, anything
FRI-shaped — cannot put that proof on Ethereum. The verifier is hundreds of thousands of gas of
hashing at best and unbounded at worst. So it wraps: a second circuit proves "I know a valid inner
proof", and *that* proof goes on chain.

Today the wrapper is almost always Groth16. Groth16 is cheap, and it needs a **circuit-specific
trusted setup** — every change to the wrapper circuit means a new ceremony, with everything that
implies for coordination, auditing and trust.

fflonk ([ePrint 2021/1167](https://eprint.iacr.org/2021/1167)) is the same order of gas from a
**universal** setup: one powers-of-tau serves every circuit, and changing the wrapper is a
recompile. That is the whole trade, and this module implements it end to end.

## 2. What it costs

Measured with `forge test`, one deployment per circuit shape, the proof arriving in memory rather
than storage (a real transaction's proof is calldata, and reading 928 bytes back out of storage
would add ~73k of cold SLOADs that never happen on chain). The external-call overhead of the
measurement harness is subtracted; the calldata figure is 16 gas per non-zero byte, 4 per zero.

| circuit | public inputs | execution gas | calldata gas | **transaction gas** |
|---|---:|---:|---:|---:|
| 2^8 rows | 1 | 189,142 | 16,424 | **226,566** |
| 2^12 rows | 4 | 198,392 | 17,948 | **237,340** |
| 2^16 rows | 16 | 216,611 | 24,020 | **261,631** |
| 2^19 rows | 4 | 198,847 | 17,948 | **237,795** |

**The cost does not grow with the circuit.** 2^19 rows costs the same as 2^12 — the only size-
dependent work is `log n` squarings to reach `ξ^n`. What does move the number is the public-input
count, at roughly 1.1k gas each for the Lagrange evaluation.

Against the other verifiers measured in this repository, on the same 2^19 circuit
(`commitment_schemes/recursion/README.md` §6–7):

| verifier | execution gas | transaction gas |
|---|---:|---:|
| UltraHonk + Shplemini/KZG, `--optimized` (what Aztec deploys) | 630,236 | 781,712 |
| UltraHonk + Vela | 573,106 | 704,826 |
| UltraHonk + fflonk *opening* | 519,164 | 674,756 |
| **fflonk over plonkish (this module)** | **198,847** | **237,795** |

**3.2× cheaper than the deployed optimized Honk verifier**, and 2.6× cheaper than the best that
swapping Honk's commitment scheme can achieve.

### Why the opening backend could not get here, and this can

`commitment_schemes/fflonk/` already implements fflonk as a PCS for UltraHonk, and its own
measurement says 519k. The floor there is Honk itself: 155,803 gas of sumcheck and relation
evaluation that no commitment scheme touches, plus ~36 column commitments to batch.

fflonk's reputation is a *PlonK* number, and it comes from PlonK's linearisation. The specific trick
here goes one step further than linearisation: because every preprocessed selector's evaluation is
*opened* rather than reconstructed in the group, all three constraint identities are checked **in
the field**. There is no linearisation MSM at all. The verifier's entire group budget is seven
scalar multiplications and one pairing — 42,900 and 113,000 gas respectively, and neither grows with
the number of selectors, the number of gates, or `n`.

### Where the 198k goes

Measured with `gasleft()` checkpoints on the 2^12 case (whose clean total is 198,392; the
checkpoints' own logging accounts for the ~1k the parts overshoot by):

| | gas |
|---|---:|
| the pairing, two pairs | 113,604 |
| the rest of the batched opening: 7 scalar multiplications, 7 additions, the Shplonk scalars | 59,953 |
| one batched inversion — Montgomery's trick, one `modexp` | 12,511 |
| parsing and validating the proof: 29 words, 5 on-curve checks | 6,956 |
| five keccak challenges | 3,052 |
| the three constraint identities | 3,339 |

The pairing alone is 57% of the verifier, and at 113,604 against a precompile price of
`45,000 + 2 × 34,000` it is within 0.5% of the EVM's floor. Of the opening's 59,953, the `ecMul` and
`ecAdd` precompiles are 43,050 — so 79% of the whole verifier is precompile cost that no amount of
further engineering can remove.

## 3. What it costs to produce

Native, on an M4 Pro (`DISABLE_ASM=1` arm64 build), one export run on an otherwise-idle machine.
Absolute wall clock on this machine is not reproducible across sessions - a concurrent build tripled
these figures - so treat the ratios between rows as the durable part:

| circuit | preprocess | prove | verify | proof |
|---|---:|---:|---:|---:|
| 2^8 rows | 6 ms | 19 ms | 0.85 ms | 928 B |
| 2^12 rows | 24 ms | 75 ms | 0.86 ms | 928 B |
| 2^16 rows | 209 ms | 923 ms | 1.17 ms | 928 B |
| 2^19 rows | 1,605 ms | 7,578 ms | 1.18 ms | 928 B |

Peak RSS for the whole four-case export was 4.1 GiB, dominated by the 2^19 case.

**The proof is 928 bytes whatever the circuit** — five group elements and nineteen scalars — plus 32
bytes per public input.

The prover is the price. fflonk packs a group of `t` polynomials into one commitment of `t` times
the degree, and the quotient needs an FFT domain of `8n`; at 2^19 rows that is a 2^22 FFT and an SRS
of 2^22 points. Against UltraHonk+KZG's 1,016 ms and 832 MiB on a circuit of the same row count,
this prover is ~7.5× slower and ~5× the memory. For something proven once and verified forever, that
is the right side of the trade; for anything proven often it is not.

The packing factor is a dial rather than a constant — see `PROTOCOL.md` §9.

## 4. How the two verifiers are kept honest

The C++ and Solidity verifiers are independent implementations, and the test suites are built so
that "it was rejected" is never accepted as evidence on its own.

**Every proof element except `W'` is absorbed into the transcript**, so tampering with any of them
breaks the *pairing*. That means a tamper test cannot, by itself, show that the three constraint
identities are checked at all. So `verify_detailed` reports each check separately and none of them
short-circuits, and the suite pins down each one:

- `ArbitraryPolynomialsOpenCorrectlyAndAreRejectedByTheConstraints` commits and opens **random**
  polynomials, honestly. The batched opening and the pairing accept them — asserted, not assumed.
  Only the three identities reject, which is the statement that they are what makes a proof mean
  anything.
- `LyingAboutAPublicInputIsCaughtByTheGateIdentity` builds a proof that is internally perfect —
  every opening honest, the pairing satisfied — for a prover claiming a public input its witness
  does not contain. Nothing but the gate identity can reject it.
- `EachIdentityConstrainsExactlyTheEvaluationsItShould` perturbs each of the nineteen evaluations in
  turn and asserts exactly which identities notice, so an identity that had silently become vacuous
  would show up.

Those tests were validated by mutation: making each of the four checks accept unconditionally
(`gate_identity`, `grand_product_start`, `permutation_identity`, `opening`) makes tests fail, and
each mutation fails a different set.

On the Solidity side, every fixture is produced by `fflonk_solidity_export_bench`, which verifies
the proof with barretenberg's own verifier before writing it out. A passing acceptance test is
therefore two implementations agreeing on the same proof, and each rejection test is both of them
rejecting the same tampering. A 256-run fuzz test asserts that flipping *any* single byte of the
proof is rejected.

The prover refuses to emit an unsatisfiable proof: every quotient division is checked for an exact
remainder, the grand product is checked to close, and the three identities are re-checked against
the prover's own claimed evaluations before the proof is returned.

## 5. Using it

```cpp
#include "barretenberg/fflonk/prover.hpp"
#include "barretenberg/fflonk/verifier.hpp"

using namespace bb::fflonk_plonk;

CircuitBuilder builder;
const uint32_t x = builder.add_public_input(FF(3));
const uint32_t y = builder.add_variable(FF(4));
const uint32_t z = builder.create_mul(x, y);
builder.fix_variable(z, FF(12));

const ProvingKey key = preprocess(builder);
const Proof proof = prove(key);
const bool ok = verify(key.verification_key, proof.to_buffer(), key.trace.public_inputs);
```

`verify` never throws and trusts nothing: it enforces the encoding (exact length, canonical scalars,
on-curve points) as well as the protocol.

Deploying the on-chain verifier takes the verification key as constructor arguments, so one contract
source serves every circuit:

```solidity
new FflonkVerifier(circuitSize, numPublicInputs, omega, k1, k2, c0x, c0y);
```

## 6. Reproducing

```bash
cd barretenberg/cpp
cmake --build build --target fflonk_tests fflonk_solidity_export_bench
./build/bin/fflonk_tests

./build/bin/fflonk_solidity_export_bench            # regenerates solidity/test/Fixture.sol
cd src/barretenberg/fflonk/solidity && forge test -vv
```

The export takes `--shapes "8:1,12:4,16:16,19:4"` (`log_n:public_inputs`) and `-o <path>`. The gas
harness needs Foundry and borrows `forge-std` from
`commitment_schemes/recursion/sol/lib`; clone it there once if it is missing.

## 7. What is not here

**The arithmetization is vanilla plonkish** — three wires, `q_M q_L q_R q_O q_C`, copy constraints,
public inputs. No lookups, no custom gates, no ACIR frontend. A wrapper circuit for a hash-based PCS
verifier has to be written against `CircuitBuilder` today; compiling one from Noir would mean an
ACIR → fflonk frontend, or widening the arithmetization to four wires and plookup, which changes the
group layout but not the shape of the verifier. That is the main thing standing between this module
and being drop-in for the use case in §1.

**The prover has not been optimised.** The `8n` quotient domain is the simple, obviously-correct
choice (`PROTOCOL.md` §9 derives why `4n` does not fit with full blinding); it is also where the
memory goes. Splitting the preprocessed group in two would drop every degree to `5n` and shrink the
SRS with it, for one more scalar multiplication (~6k gas) in the verifier.

**Zero knowledge is argued, not proven mechanically.** `PROTOCOL.md` §7 states exactly what each
opening reveals, what pins it down, and which blinder covers the rest; the layout was chosen so that
argument closes, which is why the grand product sits alone in its own group rather than in the
classic 3-group fflonk layout.
