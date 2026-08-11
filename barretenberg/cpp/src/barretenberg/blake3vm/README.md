# Blake3VM

A standalone Honk VM, in the style of the ECCVM and Translator, that proves batches of Blake3s
hashes. One proof attests to every compression in the batch; its verification cost is independent
of the batch size.

## Why

A recursive verifier for a hash-based proof system pays for every hash it re-executes. For
WHIR-Honk with Blake3s Merkle trees, a single 64-byte Merkle node costs about 2,620 UltraHonk gates
even with 8-bit XOR lookup tables, and a full recursive verifier lands at 3–4M gates — almost all
of it Merkle-path hashing. Delegating the hashing to a machine specialised for the Blake3
compression function removes that cost from the recursion: the intended flow is

1. a client produces a WHIR-Honk proof (transparent, hash-based, Blake3s Merkle trees);
2. a sequencer verifies it natively (milliseconds) for a soft confirmation;
3. the sequencer proves the verification recursively, delegating every Blake3s compression to a
   Blake3VM proof, so the recursion circuit keeps only the (cheap) non-hash verifier logic plus
   the committed chunk column that links it to the VM;
4. many such proofs are aggregated into one.

This module implements the standalone machine (step 3's hash side): trace builder, flavor,
relations, prover, and verifier, together with the complete linking argument that binds a consuming
circuit's hash claims to the trace — the VM's dense claim column, the consuming circuit's committed
chunk column, and the commitment equality between them.

Measured on an M4 Pro (native arm64 build with `DISABLE_ASM=1`): 1,000 64-byte compressions fill a
2^16-row trace and prove in **1.14 s**, verify natively in **6 ms**, with a **29,184-byte** proof
(`Blake3VMTests.DISABLED_LargeBatchProveAndVerify`). The same hashing in UltraHonk is a ~2^22
circuit.

## The pipeline, measured

`whir_settlement.test.cpp` runs every stage of that flow on one statement and prices the recursion
step both ways. Nothing in it is modelled: the client proof is produced and verified, both
recursion circuits are built, the VM proves the hash workload the delegated circuit consumed, the
settled proofs are produced and verified, and the link's commitment equality is checked between
them. It takes about 25 s and peaks near 6 GB, almost all of it the in-circuit-Blake3s baseline.

| Stage | Measurement |
|---|---|
| client proves (2^13 statement, zk WHIR) | 3.5 s, 189,440-byte proof |
| sequencer verifies natively → soft confirmation | 8 ms |
| recursion, Blake3s in circuit | **5,604,393 gates** (5,288,073 = 94% Merkle hashing) |
| recursion, hashing delegated and linked | **166,021 gates** — 33x smaller |
| ↳ of which the link's chunk side | 124,236 gates committing 37,204 chunks |
| Blake3VM over the same workload | 1,541 calls / 3,880 compressions, 2^18 trace, 6.0 s prove, 11 ms verify, 31,296-byte proof |
| settle the delegated circuit (MegaHonk) | 2^18 outer, 696 ms prove, 7 ms verify, 17,696-byte proof |
| link check `C_calldata == C_link` | accepted, both commitments read from the two proofs |
| aggregate 2 client proofs | 2^19 outer, 164,511 gates per client proof |

The delegated circuit is measured with `DelegatingBlake3sHasher`, which keeps the verifier's
interface and its digests but takes each digest from an oracle, witnesses the hashed bytes eight at
a time, and commits them in the circuit's databus calldata column — the circuit's half of the
linking argument. Two effects compound: the compressions leave the circuit, and a digest becomes
the two field elements it travels the transcript as instead of 32 bytes, which drops the Merkle
path's conditional swap from thirty-two selects per level to two. The settlement stage verifies
both proofs and checks that the circuit's calldata commitment equals the VM proof's link column
commitment. The baseline is not proved in this test — a separate disabled test settles it, priced
against the delegated path below.

Two limits the run makes concrete. The client proof is zk, but the recursion stage re-proves the
same statement without it, because the in-circuit verifier rejects salted leaves
(`whir_recursive_verifier.hpp` asserts `!config.zk`). And the VM's 6.0 s dominates the sequencer's
delegated path, so the VM prover, not the recursion circuit, is what to optimise next.

### The recursion step, proved both ways

The pipeline test leaves the 5.6M-gate baseline unproved to stay fast;
`DISABLED_BaselineInCircuitRecursionProof` settles it for real (run with
`--gtest_also_run_disabled_tests`, ~40 s). Same statement, same machine, both paths measured with
real proofs end to end:

| | Blake3s in circuit | delegated + linked |
|---|---|---|
| recursion circuit | 5,604,393 gates → 2^23 rows | 166,021 gates → 2^18 rows |
| circuit proof | 20.2 s (UltraHonk) | 0.70 s (MegaHonk) |
| Blake3VM proof | — | 6.0 s (2^18 trace) |
| **sequencer proves, total** | **20.2 s** | **6.7 s** |
| verification | 6 ms | 18 ms (both proofs + the link's group-element equality) |
| settled proof bytes | 13,632 | 48,992 (17,696 + 31,296) |
| outer-circuit rows per client proof | 5,604,393 | 164,511 — **34x more proofs per aggregator** |

The delegated path pays with a second proof to verify and ~35 KB of extra proof bytes; it wins
everywhere a sequencer is actually constrained:

- **3x prover time today, with the headroom on one side.** 90% of the delegated path is the VM
  prover, which is unoptimised — nearly every VM column is byte-valued, so small-scalar commitment
  hints and per-column sizing attack exactly where the 6.0 s goes. The circuit proof itself is
  already 29x faster than the baseline's, and the two proofs are independent, so a two-core
  sequencer sees 6.0 s wall clock.
- **The baseline's 20.2 s buys one client proof.** A 2^23 proving key is spent entirely on a single
  verification; batching N clients means N such proofs (2^25 for four clients). The delegated
  aggregator fits 34x more client proofs per outer circuit, and the VM amortises the other way:
  one VM proof covers the whole batch's hashing, with verification cost independent of batch size.
- **Memory follows circuit size.** The baseline settlement peaks near 6 GB for one client proof;
  the delegated outer circuit is a 2^18 MegaHonk key.

## Parameters

| Symbol | Meaning | Value |
|---|---|---|
| field / PCS | proof system base | BN254, KZG, non-ZK sumcheck |
| rows per compression | 56 G rows + 4 output rows | 60 (`Blake3VMTraceLayout`) |
| XOR table | rows of (x, y, x ⊕ y), x, y ∈ [0, 256) | 2^16 (the minimum trace size) |
| precomputed columns | selectors + table + Lagrange + link index | 21 |
| witness columns | 41 to-be-shifted + 78 working + 7 logup inverses | 126 |
| entities in sumcheck | precomputed + witness + 41 shifts | 188 |
| lookup reads per G row | across six logup sets | 23 |
| max relation length | the 4-read lookup subrelations | 8 (batched: 9) |

## Trace layout

Row 0 is the genesis boundary row (all-zero witness). Compression `c` occupies rows
`[1 + 60c, 60(c+1)]`: one row per G function (round r ∈ [0,7), position p ∈ [0,8)), then four
output rows. A ghost row after the last compression absorbs the final boundary constraint. The
XOR table occupies rows [0, 2^16) of its own precomputed columns and coexists with witness rows.

Each G row holds the full 16-word state `v_0..v_15` and message `m_0..m_15` (words), the byte
decompositions of every XOR operand/result, the addition results and carries, and the rotation
split pieces. State flows row to row through shifted columns:

- On a G row at position p, the four words `(A[p], B[p], C[p], D[p])` take the G outputs
  `(a2, b2, c2, d2)`; the other twelve copy. Each word's role is fixed by its index
  (j < 4 ↦ a, j ∈ [4,8) ↦ b, j ∈ [8,12) ↦ c, j ∈ [12,16) ↦ d), so each `v_j` has exactly two
  update positions per round (`Blake3VMWiringRelationImpl::UPDATE_POS`).
- Message columns copy on positions 0–6 and permute by `MSG_PERM` (= `blake3::MSG_SCHEDULE[1]`)
  on position-7 rows of rounds 0–5, so round r position p always reads `mx = m_{2p}`,
  `my = m_{2p+1}`.
- Output row t computes `out_{2t} = v_{2t} ⊕ v_{2t+8}` and `out_{2t+1} = v_{2t+1} ⊕ v_{2t+9}`,
  reusing the set-0/1 lookup byte columns. Output words copy across the four rows.
- The fourth output row is the boundary row: it pins the next row's state to the next
  compression's initial state — `IV` or, when `chain = 1`, the current output words (multi-block
  messages chain their CV this way) — and to that compression's `blk_len`/`blk_flags` in
  `v_14`/`v_15`. The counter words `v_12, v_13` are pinned to 0 (single-chunk inputs ≤ 1024 bytes,
  matching `blake3::blake3s`).

### Rotations and range checks

All words are little-endian byte limbs. Every byte is bound to [0, 256) by an XOR-table read; the
table read (x, y, z) also proves z = x ⊕ y. Rotations by 16 and 8 are byte permutations.
Rotation by 12 splits z2's byte 1 into nibbles; rotation by 7 splits z4's byte 0 into a 7-bit part
and its top bit. A k-bit piece v is range-checked by reading the table row `(v·2^{8−k}, 0)`.
Additions mod 2^32 carry explicitly (`carry ∈ {0,1}` two-term, `{0,1,2}` three-term) with the
result bounded by its byte decomposition. Message words are byte-decomposed and range-checked at
every use (`mx_b`/`my_b` with the `junk_*` XOR outputs discarded).

The relation inventory lives in `relations/blake3vm/blake3vm_relations.hpp` (G algebra, wiring,
zero row), `relations/blake3vm/blake3vm_lookup_relation.hpp` (the six logup sets), and
`relations/blake3vm/blake3vm_link_relation.hpp` (the claim exposure and its gather); each file's
header comment is the constraint-level reference.

## Transcript schedule

| Step | Labels |
|---|---|
| VK binding (hash buffer) | `vk_circuit_size`, `vk_<precomputed entity name>` |
| 119 wire commitments | the witness entity names (`v_0` … `counts_5`) |
| logup challenges | `beta`, `gamma` (β² is derived, not drawn) |
| 7 inverse commitments | `inv_0` … `inv_5`, `inv_link` |
| sumcheck | `Sumcheck:alpha`, `Sumcheck:gate_challenge`, then the standard non-ZK rounds |
| opening | Shplemini + KZG standard labels |

`Blake3VMProver::construct_proof` and `Blake3VMVerifier::verify_proof` implement the two sides.

## Soundness notes

- Logup soundness does not require the read-count columns to be small: a read denominator absent
  from the table side leaves an unmatched pole regardless of the count values.
- The shifted polynomial views require a zero row 0; `Blake3VMZeroRowRelation` enforces it for all
  40 to-be-shifted columns, and `lagrange_first · chain = 0` stops the genesis compression from
  chaining an arbitrary CV.
- `blk_len`/`blk_flags` are range-checked at first use (they enter round 0 as XOR operands whose
  bytes are table reads). Their *semantics* (that `blk_len` is the true byte length, flags follow
  the chunk schedule) are the batch author's claim, exposed per compression. The linking argument
  binds them through the digest claims: the parameters enter the compression the VM proves, so a
  batch with a different schedule produces a different digest for the same message bytes, and a
  consuming circuit whose digests must chain to a committed Merkle root cannot absorb that
  difference without a Blake3s collision.
- The VK is per-circuit-size: the precomputed selectors depend only on the trace length. A
  deployment would pin a fixed size the way the ECCVM does.

## The linking argument

Delegation on its own proves nothing: the delegated circuit reads Merkle siblings off the proof
stream and only ever hands them to the hash oracle, so with the compressions gone, tampering with a
sibling makes the inner proof fail native verification while the circuit stays satisfiable when the
oracle replays the honest digests. In-circuit hashing was the only thing constraining that sibling.

The link must establish that the sequence of `(input, digest)` pairs the circuit consumed is the
sequence the VM proved. Three of the four ways to do that are unsound or unaffordable, and ruling
them out pins the design.

**A fingerprint under a challenge fixed before the circuit is built is forgeable.** With claims
`c_i` and a Horner fingerprint `F = Σ cᵢ γⁱ`, a prover that learns `γ` before choosing the `c_i`
solves one linear equation in one free claim. Every digest the oracle answers is such a free claim:
the prover picks fake digests to satisfy the Merkle root checks, then tunes one of them to land `F`
on the target. This holds however `γ` is derived — from the inner proof, from the VM's
commitments, or from a hash of both — as long as the circuit's own claims are not committed first.

**Deriving the challenge in-circuit costs more than it saves.** `γ = Poseidon2(claims)` is
unforgeable because `γ` depends on the claims, but it is 3,530 permutations for 10,590 elements,
about 265,000 gates, and the VM cannot recompute a Poseidon2 hash to match it.

**Publishing the claims defeats the purpose.** 10,590 field elements is a 338 KB settled proof, and
an aggregator would have to absorb every one of them.

What remains is the standard construction, and it is the one Goblin already uses to bind an Ultra
circuit to the ECCVM: **the claims live in a committed column of each proof, and the two are
compared**. Since both proofs are KZG over BN254 against the same SRS, and a commitment is
`Σ vᵢ·Gᵢ`, two columns carrying the same values at the same indices have the *same group element*.
The check is `C_circuit == C_vm` — no shared challenge, no in-circuit hashing, and binding follows
from KZG binding alone.

The claim encoding is the hashed bytes themselves, eight at a time: entry `j` is
`w_lo + 2³²·w_hi` for the word pair the trace exposes on one row. Chunking at eight bytes rather
than at field elements is what keeps the VM side row-local, since the one-byte domain tag offsets
every field element against the VM's 32-bit words.

Both halves are built:

- **The VM's half.** `claims` carries the sequence, `Blake3VMClaimRelation` pins each row's
  contribution to the trace data that row already holds, and `Blake3VMLinkRelation` gathers them
  into the dense `link_value` column. `Blake3VMTests.ClaimVectorIsTheHashedData` rebuilds the
  hashed byte stream from the claims, and `TamperedLinkColumnFails` shows a link column that does
  not carry the trace's claims cannot prove.
- **The circuit's half.** `DelegatingBlake3sHasher` witnesses each call's eight-byte chunks, proves
  they recombine to the field elements it holds, and appends every chunk to the Mega databus
  calldata column. A tagged buffer offsets every piece (16-byte digest halves, 32-byte leaf values)
  one byte against the chunk grid, so each piece carries a `(7 bytes, 8-byte middles, 1 byte)`
  split; both split parts are range-constrained. The ranges are what make the split sound: the
  chunks are pinned to the VM's byte-built claims by the commitment equality, and range-checked
  splits make the byte decomposition unique, so the recombination pins the field element itself —
  without them a prover can shift one piece and a neighbour's split in tandem and claim different
  field elements recombine to the same honest chunks. Chained buffers open with the previous digest
  instead of the tag, stay 8-byte aligned, and need neither splits nor ranges.
- **The comparison.** Both proofs are KZG over the same SRS and both columns start at
  `NUM_DISABLED_ROWS_IN_SUMCHECK` — the offset bb's databus polynomials use for cross-flavor
  commitment compatibility, which is the same property this link relies on — so
  `C_calldata == C_link` as group elements is the entire check, performed by whoever verifies both
  proofs. Nothing about the claims crosses the transcript.

`WhirSettlementTests.TamperedSiblingDivergesTheLinkCommitments` is the argument end to end: the
honest circuit's calldata commitment equals the VM's link commitment, while a tampered sibling —
still satisfiable in-circuit, since the oracle's digests are free witnesses — diverges from the
link column of every VM proof a prover could produce: the honest workload differs in the tampered
message chunks, and the circuit's actual calls hash to digests that differ from the replayed honest
ones. `FullPipelineWithAndWithoutTheVM` performs the check on commitments read from the two real
proofs.

Measured, the circuit's half costs 124,236 gates for 37,204 chunks (busread gates, splits, and
range checks), putting the linked recursion circuit at 166,021 gates — against 5,604,393 with the
hashing in circuit.

## Roadmap

- **Cheaper chunk side**: the link's 124,236 gates are three quarters of the delegated circuit —
  one busread gate per chunk plus the splits and range checks. A dedicated append-only bus gate
  (no read machinery), batched range checks, and skipping the zero-padding chunks' busreads would
  each take a measurable slice.
- **zk on the recursion path**: salted leaves are one extra absorbed element per leaf, which the
  delegated hasher would carry for free; the blinding and the extra variable in the WHIR config are
  the real work.
- **Recursive Blake3VM verifier**: sumcheck over 184 entities plus a ~131-term Shplemini MSM,
  independent of the batch size; needed for step 4 of the flow.
- **ZK**: the machine is sequencer-side, so the non-ZK sumcheck suffices for now; Libra masking can
  be added the way the Translator does if the trace ever carries private data.
- **Prover speed**: most columns are byte-valued, so commitment hints for small scalars, per-column
  physical sizing, and relation `skip()` methods are the next wins.
