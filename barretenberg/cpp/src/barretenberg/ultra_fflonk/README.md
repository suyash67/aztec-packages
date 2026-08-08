# ultra_fflonk: fflonk over the arithmetization Noir already compiles to

A proof system with **UltraHonk's circuit and fflonk's verifier shape**. The circuit builder, the
ACIR frontend, the lookup tables, the custom gates and the constraint set are UltraHonk's, taken
unchanged; the proof on top of them is a univariate quotient and a batched fflonk opening instead of
sumcheck and Shplemini.

`PROTOCOL.md` is the specification. `prover.cpp` and `verifier.cpp` are two implementations of it.

```bash
bb prove  --scheme ultra_fflonk --write_vk -b target/program.json -w target/witness.gz -o out
bb verify --scheme ultra_fflonk -k out/vk -p out/proof -i out/public_inputs
```

## 1. Why not just use `fflonk/`

The three-wire system in `../fflonk/` verifies a whole Ethereum transaction for 237,795 gas, but its
arithmetization is vanilla plonkish: three wires, five selectors, copy constraints. No lookups, no
custom gates, no frontend. Getting a Noir program onto it means either 3–5× the rows or a second
ACIR frontend to write and audit — and lookups have no three-wire encoding that is not astronomically
expensive, because "expand the lookup into explicit constraints" costs the size of the table.

`UltraCircuitBuilder` is already the four-wire, plookup, custom-gate arithmetization every Noir
program compiles to. Replacing only the proof system on top of it means the frontend, every
constraint type and the whole gate-count profile carry over for free, and the only new cryptographic
surface is `PROTOCOL.md`.

Both are exposed. `../fflonk_acir/` lowers arithmetic and Poseidon2 gates onto the three-wire system
for `--scheme fflonk`, which is the right choice when a program is built from those alone and the
cheapest possible verifier is what matters. Anything using integer types, comparisons, bitwise
operations, the SHA/Blake/Keccak hashes or dynamically indexed arrays belongs here — see
`../fflonk_acir/README.md` for exactly why each of those cannot cross.

## 2. What it does, measured

A real Noir program (`complete_age_check`, 130,233 ACIR opcodes → **427,716 gates**, 2^19 rows, 11
public inputs), proven and verified through `bb` on an M4 Pro (`DISABLE_ASM=1` arm64 build):

| | |
|---|---|
| proof size | **2,112 bytes**, flat in circuit size |
| verification key | 384 bytes — four group elements, whatever the circuit contains |
| peak prover RSS | 4.14 GiB |
| native verification | ~25 ms |

The proof is six group elements and fifty-four scalars regardless of the circuit; the verifier's
group work is eleven scalar multiplications and one pairing, and its field work is one evaluation of
the nine Ultra relations.

### On chain

`solidity/` is the verifier, and `bb write_solidity_verifier --scheme ultra_fflonk` emits it with the
key baked in. Measured with `forge test`, the proof arriving in memory rather than storage, with the
harness's own external-call overhead subtracted:

| circuit | rows | public inputs | execution gas | calldata gas | **transaction gas** |
|---|---:|---:|---:|---:|---:|
| arithmetic | 2^5 | 4 | 272,841 | 28,928 | **322,769** |
| every gate type | 2^14 | 1 | 276,206 | 34,244 | **331,450** |
| arithmetic | 2^13 | 4 | 277,861 | 28,844 | **327,705** |
| `complete_age_check`, a real Noir program | 2^19 | 11 | 282,381 | 33,604 | **336,985** |

**The cost does not grow with the circuit** — 2^19 rows costs 3% more than 2^5, and that is the
public inputs, not the size. Against the other verifiers measured in this repository:

| verifier | transaction gas |
|---|---:|
| UltraHonk + Shplemini/KZG, `--optimized` (what Aztec deploys) | 781,712 |
| UltraHonk + fflonk *opening* | 674,756 |
| **fflonk over the Ultra arithmetization (this module)** | **336,985** |
| fflonk over three-wire plonkish (`../fflonk/`, no lookups) | 248,109 |

**2.3× cheaper than the deployed optimized Honk verifier.** The gap to the three-wire system —
about 75k of execution and 19k of calldata — is what the lookup-bearing arithmetization costs at the
verifier: thirty-one subrelations instead of three, eight groups instead of four, 2,112 proof bytes
instead of 928. In exchange the circuit needs 3–5× fewer rows and every Noir construct works.

## 3. How the verifier is kept honest

Every proof element except `W'` is absorbed into the transcript, so tampering with any of them breaks
the *pairing*. That means "the proof was rejected" is not, on its own, evidence that the constraint
system is checked at all. So `verify_detailed` reports each check separately, none of them
short-circuits, and the suite pins each one down:

- **`ArbitraryPolynomialsOpenCorrectlyAndAreRejectedByTheQuotientIdentity`** commits and opens
  **random** polynomials, honestly. The batched opening and the pairing accept them — asserted, not
  assumed. Only the quotient identity rejects, which is the statement that the relations are what
  make a proof mean anything.
- **`EveryEvaluationTheRelationsReadIsConstrained`** perturbs each of the fifty-four claimed
  evaluations in turn and asserts exactly which ones the quotient identity notices. It notices all of
  them except three — the lookup read counts, read tags and inverses at `xi*omega`, which no relation
  reads and which are only revealed because they share a two-point group. Naming that set exactly is
  what makes the test sharp.
- **`TheQuotientIdentityReadsThePublicInputsThroughTheirDelta`** builds a proof that is internally
  perfect — every opening honest, the pairing satisfied, the same transcript on both sides — but where
  the verifier derives a different `public_input_delta`. Only the quotient identity can reject it.
- **`ForgedOpeningsCannotSubstituteThePreprocessedGroups`** confirms the verification key really binds
  the preprocessed columns: substituting one makes the *opening* fail, not merely the identity.

Those tests were validated by mutation. Making each check accept unconditionally makes tests fail,
and each mutation fails a different set:

| mutation | tests that fail |
|---|---|
| `report.quotient = true` | `ArbitraryPolynomials…`, `EveryEvaluation…`, `TheQuotientIdentityReadsThePublicInputs…` |
| `report.opening = true` | `RejectsATamperedCommitment`, `ForgedOpeningsCannotSubstitute…`, `EveryEvaluation…` |
| encoding checks always accept | `RejectsMalformedEncodings` (in this suite and in `fflonk/`'s) |

The prover refuses to emit an unsatisfiable proof: the quotient division is checked for an exact
remainder, each running sum is checked to close, the realised quotient degree is checked against the
group layout, and the batched identity is re-checked against the prover's own claimed evaluations
before the proof is returned.

## 4. The one thing to know about the design

Sumcheck imposes two different obligations on a subrelation, and only one of them is a quotient
claim. Twenty-nine of the Ultra set's thirty-one subrelations must vanish at every row — those divide
by `Z_H`. Two must only *sum to zero across the trace*, which dividing by `Z_H` would prove the wrong
thing about.

Each of those two gets a committed running sum `S` constrained by `S(omega X) − S(X) − f(X) = 0`,
whose sum over the domain telescopes around the cycle and leaves exactly `sum_rows f = 0`. The count
is derived from the flavor rather than hard-coded, so adding a third linearly dependent subrelation
breaks the build instead of silently going unenforced.

Everything else — degrees, the group layout, the SRS bound, the zero-knowledge argument — is in
`PROTOCOL.md` §3–§8, derived rather than assumed.

## 5. How the Solidity verifier is kept honest

It is a hand port of nine C++ relation templates, and the batched identity folds all thirty-one
subrelations into a single field element — inside which one mistranslated term is invisible, and
indistinguishable from a hundred of them. A proof either verifies or does not, so the fixture tests
alone cannot tell a correct port from a uniformly broken one.

So `ultra_fflonk_solidity_export_bench` exports random rows together with barretenberg's own value
for **each subrelation separately**, and `RelationVectors.t.sol` checks all thirty-one. That is the
only test in the suite that can say *which* relation is wrong, and it is what makes the port
reviewable. On top of it, every fixture is a proof barretenberg's own verifier accepted before it was
written out, so acceptance is two implementations agreeing; and a 256-run fuzz test asserts that
flipping any single byte of the proof is rejected.

## 6. What is not here

- **The prover is the price.** 4.14 GiB and minutes at 2^19 rows: the quotient needs an `8n` FFT
  domain, thirty-eight columns have to reach it, and the nine relations are evaluated at every one of
  its `8n` points. Working one coset at a time keeps memory at `38n` field elements rather than
  `38·8n`, which is the difference between "large" and "impossible", but this is a proof system for
  something proven once and verified forever.
- **Zero knowledge is argued, not proven mechanically** — `PROTOCOL.md` §8 states what each opening
  reveals and which blinder covers it; there is no simulator.
- **No bb.js backend.** The scheme is reachable from `bb` but not yet from `ts/bb.js`.

## 7. Reproducing

```bash
cd barretenberg/cpp
cmake --build build --target ultra_fflonk_tests ultra_fflonk_solidity_export_bench
./build/bin/ultra_fflonk_tests

./build/bin/ultra_fflonk_solidity_export_bench   # regenerates solidity/test/{RelationVectors,Fixture}.sol
cd src/barretenberg/ultra_fflonk/solidity && forge test -vv
```

The gas harness needs Foundry and borrows `forge-std` from
`commitment_schemes/recursion/sol/lib`; clone it there once if it is missing. `generate_contract_header.sh`
reflattens the contract into `dsl/acir_proofs/ultra_fflonk_contract.hpp` after any change to `solidity/src/`,
which is what `bb write_solidity_verifier` emits.
