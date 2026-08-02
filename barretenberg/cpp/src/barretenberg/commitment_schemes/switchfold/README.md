# SwitchFold — recursive code switching over Reed–Solomon

Implementation of SwitchFold ([eprint 2026/1489](https://eprint.iacr.org/2026/1489)) instantiated
over the RS code sequence already present in `whir/rs_code.hpp`, sharing Ligero's commitment.

## 1. What it changes

Ligero's proof has two parts: the query openings, and the **combined rows sent in the clear**. The
second is `3 · num_cols` field elements — `O(√N)` — and at 2^18 it is about 60% of the proof.
SwitchFold sends nothing of row length. The combined rows become a fresh message that is
*recursively code-switched* down a sequence of RS codes with geometrically shrinking block lengths,
until a 16-element base message remains and is transmitted directly.

The commitment is untouched, so `LigeroHonk` and `SwitchFoldHonk` are a controlled A/B: identical
trace, identical Merkle tree, identical payload queries, different opening argument.

## 2. Why Reed–Solomon collapses two of the paper's three modules

The paper's framework has three modules: recursive code switching, accumulation of
generator-matrix MLE claims, and a final recursive opening of those accumulated claims. Modules two
and three exist because Brakedown's generator matrix is sparse-but-large, so the code-switching
identity

```
Ξ · Enc[C_j](m_j) = Ξ G_j m_j
```

leaves an unresolved claim about `G_j` that must be committed, accumulated across openings, and
eventually opened.

For an RS code the generator row at domain point `ω^s` is `(1, ω^s, ω^{2s}, …)` — the **pow tensor**.
So `(G_j m_j)[s] = m̂_j(ω^s) = ⟨m_j, pow(ω^s)⟩` is already an ordinary tensor claim on the next
message, of exactly the kind the recursion carries. No generator-matrix commitment, no accumulator,
no auxiliary foldable code. This is the paper's own point that SwitchFold "decouples PCS design from
code design", read in the direction it does not emphasise: choosing a structured code buys back
protocol simplicity, at the cost of the linear-time encoding that motivated Brakedown.

## 3. The descent

State entering a level: a message `m` of length `k`, and a list of pending claims about it, each a
sum of `WeightTerm`s (`whir/weights.hpp`) with a claimed value.

1. **Commit.** Reshape `m` into `2^log_switch_factor` columns of length `k / 2^log_switch_factor`,
   RS-encode each, and Merkle-commit the interleaved codeword (one leaf per domain position, holding
   every column's value). Send the root.
2. **Query.** Draw `num_queries` positions and open them.
3. **Batch and collapse.** Draw `μ`, batch the pending claims into a single inner product
   `⟨m, W⟩ = V`, and run one sumcheck (`log k` rounds, two field elements each) reducing it to a
   single tensor claim `⟨m, ts(x)⟩ = v` — the paper's Fact 1.
4. **Split.** `ts(x)` factors as `ts(x_low) ⊗ ts(x_high)` because the low `log(k/t)` variables index
   within a column. So `⟨m', ts(x_low)⟩ = v` where `m' = Σ_c ts(x_high)[c] · column_c` — and `m'` is
   the next level's message.
5. **Carry the queries.** For each queried position `s`, `Enc(m')[s]` is the same
   `ts(x_high)`-combination of the opened leaf, which the verifier computes itself. That yields
   `⟨m', pow(ω^s)⟩ = Enc(m')[s]`, one code-switching claim per query, pending for the next level.

At the base the message is sent in the clear and every pending claim is checked directly.

## 4. The four segments

The Honk claim set arrives as Ligero's rank-2 shift decomposition: three combined rows `w_u`
(unshifted), `w_s` (shifted) and `w_s2` (the shifted carry). They are laid out as the four segments
of one descent message of length `4 · num_cols`, the segment index occupying the top two variables.
This keeps a **single** descent, and every claim stays a tensor:

| claim | weight |
|---|---|
| the evaluation claim, unshifted part | `eq(u_lo)` on segment 0 |
| its shifted part | `shr(eq(u_lo))` on segment 1 — `log num_cols` eq terms, see below |
| its carry term | `eq(0)` scaled by `a[C-1]` on segment 2 |
| payload query `s`, segment `σ` | `pow(ω^s)` on segment `σ` |

Pinning a segment bit is an eq factor at the constant 0 or 1, which is why `WeightTerm` needs only
one extra primitive (`append_eq_variable`) to express all of these.

`shr(eq(u))[c] = eq(u, c-1)` is not a tensor, but it is a sum of `log num_cols` of them. Adding one
to `c` flips a maximal run of low 1-bits to 0 and the next 0-bit to 1, so the vector splits by carry
position `k`: bits below `k` pinned to 0, bit `k` pinned to 1, bits above `k` carrying the ordinary
eq factor.

## 5. Parameters

Because nothing of row length is transmitted, the proof-size optimum moves. Ligero picks `num_cols`
to balance the leaf width (`num_rows × num_polynomials` values per query) against the transmitted
rows; with the second term gone, `num_cols` goes as large as BN254's 2-adicity allows and the leaves
shrink accordingly. `SwitchFoldConfig::create` therefore takes `log_num_cols = num_variables - 1`,
backing off only if the descent codeword — which carries two extra segment variables — would exceed
2^28.

`log_switch_factor` (4 columns per level) and `log_base_length` (a 16-element base) trade level
count against per-level query cost.

## 6. Costs and what the measurements say

Against `LigeroHonk` at 2^18 with the Blake3s hasher and an identical commitment: proof size falls
**3.5x** (matching the paper's headline against Brakedown) and verification roughly halves. The
prover is slower, for two reasons that should not be conflated:

- The paper's linear prover comes from Brakedown's linear-time encoder. Over RS the descent pays
  FFT encoding at every level.
- The batched claim carries one term per query per segment, and expanding those into the sumcheck
  weight table is multiplication-heavy. This arm64 build has `DISABLE_ASM=1`, where a BN254 Fr
  multiplication costs ~4.2 additions, so that phase is penalised relative to Ligero's
  FFT-and-hash profile. The table build is parallelised over terms with private accumulators; it
  was 2.6x slower before that.

See `../PCS_CANDIDATES.md` §5 and §7 for the numbers in context and for what a linear-time code
module would change.
