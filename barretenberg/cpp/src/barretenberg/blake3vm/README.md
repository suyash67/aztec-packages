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
relations, prover, and verifier. The linking argument that connects a consuming circuit's hash
claims to the VM trace is future work (see Roadmap).

Measured on an M4 Pro (native arm64 build with `DISABLE_ASM=1`): 1,000 64-byte compressions fill a
2^16-row trace and prove in **1.14 s**, verify natively in **6 ms**, with a **29,184-byte** proof
(`Blake3VMTests.DISABLED_LargeBatchProveAndVerify`). The same hashing in UltraHonk is a ~2^22
circuit.

## Parameters

| Symbol | Meaning | Value |
|---|---|---|
| field / PCS | proof system base | BN254, KZG, non-ZK sumcheck |
| rows per compression | 56 G rows + 4 output rows | 60 (`Blake3VMTraceLayout`) |
| XOR table | rows of (x, y, x ⊕ y), x, y ∈ [0, 256) | 2^16 (the minimum trace size) |
| precomputed columns | selectors + table + Lagrange | 19 |
| witness columns | 40 to-be-shifted + 79 working + 6 logup inverses | 125 |
| entities in sumcheck | precomputed + witness + 40 shifts | 184 |
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
zero row) and `relations/blake3vm/blake3vm_lookup_relation.hpp` (the six logup sets); each file's
header comment is the constraint-level reference.

## Transcript schedule

| Step | Labels |
|---|---|
| VK binding (hash buffer) | `vk_circuit_size`, `vk_<precomputed entity name>` |
| 119 wire commitments | the witness entity names (`v_0` … `counts_5`) |
| logup challenges | `beta`, `gamma` (β² is derived, not drawn) |
| 6 inverse commitments | `inv_0` … `inv_5` |
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

## Roadmap

- **Linking argument**: expose the per-compression IO (message words at the first G row;
  `blk_len`/`blk_flags`/`chain` on the preceding boundary row; output words on the output rows) to
  a consuming circuit via a databus-style log-derivative argument over a shared transcript, so a
  Noir/UltraHonk recursive WHIR verifier can consume hash results without hashing.
- **Recursive Blake3VM verifier**: sumcheck over 184 entities plus a ~131-term Shplemini MSM,
  independent of the batch size; needed for step 4 of the flow.
- **ZK**: the machine is sequencer-side, so the non-ZK sumcheck suffices for now; Libra masking can
  be added the way the Translator does if the trace ever carries private data.
- **Prover speed**: most columns are byte-valued, so commitment hints for small scalars, per-column
  physical sizing, and relation `skip()` methods are the next wins.
