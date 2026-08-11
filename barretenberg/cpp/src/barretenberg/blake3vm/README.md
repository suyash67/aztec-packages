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
   Blake3VM proof, so the UltraHonk circuit keeps only the (cheap) non-hash verifier logic;
4. many such proofs are aggregated into one.

This module implements the standalone machine (step 3's hash side): trace builder, flavor,
relations, prover, and verifier, together with the VM's half of the linking argument that binds a
consuming circuit's hash claims to the trace.

Measured on an M4 Pro (native arm64 build with `DISABLE_ASM=1`): 1,000 64-byte compressions fill a
2^16-row trace and prove in **1.14 s**, verify natively in **6 ms**, with a **29,184-byte** proof
(`Blake3VMTests.DISABLED_LargeBatchProveAndVerify`). The same hashing in UltraHonk is a ~2^22
circuit.

## The pipeline, measured

`whir_settlement.test.cpp` runs every stage of that flow on one statement and prices the recursion
step both ways. Nothing in it is modelled: the client proof is produced and verified, both
recursion circuits are built, the VM proves the hash workload the delegated circuit consumed, and
the settled proofs are produced and verified. It takes about 20 s and peaks near 6 GB, almost all
of it the in-circuit-Blake3s baseline.

| Stage | Measurement |
|---|---|
| client proves (2^13 statement, zk WHIR) | 4.8 s, 189,440-byte proof |
| sequencer verifies natively → soft confirmation | 10 ms |
| recursion, Blake3s in circuit | **5,604,393 gates** (5,288,073 = 94% Merkle hashing) |
| recursion, hashing delegated | **49,553 gates** — 113x smaller |
| ↳ of which the link fingerprint | 10,590 gates binding 10,590 field elements |
| Blake3VM over the same workload | 1,541 calls / 3,880 compressions, 2^18 trace, 6.4 s prove, 10 ms verify, 30,080-byte proof |
| settle the delegated circuit | 2^16 outer, 244 ms prove, 5 ms verify, 13,664-byte proof |
| aggregate 2 client proofs | 2^17 outer, 48,027 gates per client proof |

The delegated circuit is measured with `DelegatingBlake3sHasher`, which keeps the verifier's
interface and its digests but takes each digest from an oracle and binds the call's inputs and
output into a Horner fingerprint. Two effects compound: the compressions leave the circuit, and a
digest becomes the two field elements it travels the transcript as instead of 32 bytes, which
drops the Merkle path's conditional swap from thirty-two selects per level to two.

**That circuit is not sound on its own.** Nothing yet forces the oracle's digests to be the hashes
the VM proved; the fingerprint is the hook a linking argument would use, and it is what the 10,590
gates buy. The baseline is not proved in the test either — 5.6M gates is a 2^23 outer circuit,
against 2^16 for the delegated one.

Two limits the run makes concrete. The client proof is zk, but the recursion stage re-proves the
same statement without it, because the in-circuit verifier rejects salted leaves
(`whir_recursive_verifier.hpp` asserts `!config.zk`). And the VM's 6.4 s dominates the sequencer's
delegated path, so the VM prover, not the recursion circuit, is what to optimise next.

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
  the chunk schedule) are the batch author's claim, exposed per compression; a consumer of the VM
  must bind them, which is the linking argument's job.
- The VK is per-circuit-size: the precomputed selectors depend only on the trace length. A
  deployment would pin a fixed size the way the ECCVM does.

## The linking argument

`whir_settlement.test.cpp` demonstrates what is missing: tampering with a Merkle sibling in the
query openings makes the inner proof fail native verification, while the delegated circuit stays
satisfiable when the oracle replays the honest digests. In-circuit hashing was the only thing
constraining that sibling.

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

The VM's half is built: `claims` carries the sequence, `Blake3VMClaimRelation` pins each row's
contribution to the trace data that row already holds, and `Blake3VMLinkRelation` gathers them into
the dense `link_value` column. `Blake3VMTests.ClaimVectorIsTheHashedData` rebuilds the hashed byte
stream from the claims, and `TamperedLinkColumnFails` shows a link column that does not carry the
trace's claims cannot prove.

The claim encoding is the hashed bytes themselves, eight at a time: entry `j` is
`w_lo + 2³²·w_hi` for the word pair the trace exposes on one row. Chunking at eight bytes rather
than at field elements is what keeps the VM side row-local, since the one-byte domain tag offsets
every field element against the VM's 32-bit words.

What remains is the consuming circuit's half:

- **Emit the same vector.** The delegating hasher builds each call's byte buffer already, so it can
  witness the eight-byte chunks directly. It must then prove they recombine to the field elements
  it holds — the buffer is a run of pieces of mixed width (a 1-byte tag, 16-byte digest halves,
  32-byte leaf values), each offset one byte by the tag, so each piece boundary needs a
  `(1 byte, 7 bytes)` split with the byte range-constrained. About 2-3 gates per piece.
- **Commit it.** `MegaCircuitBuilder`'s databus gives a committed column of circuit witnesses. Its
  polynomial starts at `NUM_DISABLED_ROWS_IN_SUMCHECK` rather than 0 — an offset bb maintains for
  cross-flavor commitment compatibility, which is the same property this link relies on — so the
  VM's link section must start there too.
- **Compare.** `C_calldata == C_link` as group elements, checked by whoever verifies both proofs.
  `DelegatedCircuitIsForgeableWithoutTheLink` then becomes the green test: the tampered sibling
  changes the circuit's claims, so the two commitments diverge.

Costed at ~42,500 chunks and ~15,000 pieces, the circuit's half is roughly 87,000 gates, putting
the linked recursion circuit near 140,000 — against 5,604,393 with the hashing in circuit.

## Roadmap

- **Linking argument**: expose the per-compression IO (message words at the first G row;
  `blk_len`/`blk_flags`/`chain` on the preceding boundary row; output words on the output rows) to
  a consuming circuit via a databus-style log-derivative argument over a shared transcript, so a
  Noir/UltraHonk recursive WHIR verifier can consume hash results without hashing. The VM side is
  the missing half: it must recompute the fingerprint `whir_settlement.test.cpp` builds, over the
  32-byte message the circuit absorbs rather than the 64-byte blocks the trace holds.
- **Cheaper binding**: the fingerprint costs one gate per field element, and the leaf values it
  absorbs are already transcript data the circuit holds. Binding a leaf by a value the verifier
  computes anyway would remove most of the 10,590 gates.
- **zk on the recursion path**: salted leaves are one extra absorbed element per leaf, which the
  delegated hasher would carry for free; the blinding and the extra variable in the WHIR config are
  the real work.
- **Recursive Blake3VM verifier**: sumcheck over 184 entities plus a ~131-term Shplemini MSM,
  independent of the batch size; needed for step 4 of the flow.
- **ZK**: the machine is sequencer-side, so the non-ZK sumcheck suffices for now; Libra masking can
  be added the way the Translator does if the trace ever carries private data.
- **Prover speed**: most columns are byte-valued, so commitment hints for small scalars, per-column
  physical sizing, and relation `skip()` methods are the next wins.
