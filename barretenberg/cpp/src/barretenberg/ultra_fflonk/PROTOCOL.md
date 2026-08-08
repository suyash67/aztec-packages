# fflonk over the Ultra arithmetization

A SNARK with UltraHonk's circuit and fflonk's verifier: the same `UltraCircuitBuilder` trace, the
same `relations/*.hpp` constraint set, the same ACIR frontend — proven with a univariate quotient and
a batched fflonk opening instead of sumcheck and Shplemini.

This is the specification. `prover.cpp` and `verifier.cpp` are two implementations of it; where they
disagree, this document is what needs fixing first. The three-wire system in `../fflonk/` is a
separate, self-contained proof system with its own `PROTOCOL.md`; what is shared between them is the
packing, the transcript and the batched opening, in `../fflonk/{polynomial_utils,transcript,
batched_opening}.hpp`.

## 1. Why this exists

`../fflonk/` gets a whole Ethereum transaction down to 237,795 gas, but its arithmetization is
vanilla plonkish — three wires, five selectors, copy constraints — with no lookups, no custom gates
and no ACIR frontend. A Noir program compiled to that arithmetization needs 3–5× the rows, and
somebody has to write and audit the frontend that gets it there.

`UltraCircuitBuilder` is already the four-wire, plookup, custom-gate arithmetization that every Noir
program compiles to. Replacing only the proof system on top of it means the frontend, every
constraint type and the whole gate-count profile carry over unchanged, and the *only* new
cryptographic surface is the argument in this document.

The price, paid once per proof rather than once per verification, is in §7.

## 2. Arithmetization

Fix `n = 2^k`, `ω` a primitive `n`-th root of unity, `H = {ω^i}`, `Z_H(X) = X^n − 1`. Row `i` of the
Honk trace is identified with `ω^i`.

The 41 entities are exactly `flavor/generated/ultra_flavor_generated.hpp`:

| kind | count | columns |
|---|---:|---|
| precomputed | 28 | `σ_1..σ_4`, `id_1..id_4`, `L_first`, `L_last`, `q_lookup`, `table_1..table_4`, `q_m`, `q_r`, `q_o`, `q_c`, `q_l`, `q_4`, `q_arith`, `q_delta_range`, `q_elliptic`, `q_memory`, `q_nnf`, `q_poseidon2_external`, `q_poseidon2_internal` |
| witness | 8 | `w_l`, `w_r`, `w_o`, `w_4`, `z_perm`, `lookup_inverses`, `lookup_read_counts`, `lookup_read_tags` |
| shifted | 5 | `w_l`, `w_r`, `w_o`, `w_4`, `z_perm`, each read at the next row |

**Shifts are cyclic here**: `p_shift` is the polynomial `p(ωX)`, so at row `n−1` it wraps to `p[0]`.
Honk's shift is `p_shift[i] = p[i+1]` with `p_shift[n−1] = 0`. The two agree because barretenberg
allocates every to-be-shifted polynomial with `Polynomial::shiftable`, which fixes `p[0] = 0`. That
identity is what lets this proof system consume Honk's polynomials untouched, and it is the one place
where a change to the Honk allocation would silently break this one.

**Relations.** The nine relations of `UltraFlavor::Relations` contribute 31 subrelations, called
through the same `accumulate` entry point Sumcheck uses. `relation_batch.hpp` is the only place the
arithmetization is read, and both the prover (at every point of the quotient domain) and the verifier
(once, at `ξ`) go through it — so "the verifier checks the same relations the prover divided out" is
a property of the code, not of a review.

Sumcheck imposes two different obligations on a subrelation, and only one of them is a quotient
claim:

- 29 are **linearly independent**: they must vanish at every row. These divide by `Z_H`.
- 2 are **linearly dependent**: they must only *sum to zero across the trace*. These do not.

The two are `LogDerivLookupRelation`'s lookup identity and `MemoryRelation`'s ROM-LogUp sum. §4
handles them.

**Public inputs** are not a separate identity. They enter the permutation argument through
`public_input_delta`, computed by the verifier from the claimed values, `β`, `γ` and the verification
key's `pub_inputs_offset` — the same `compute_public_input_delta` the Honk verifier uses. That is the
only place the public inputs reach the constraint system, and it is a *multiplicative* correction the
prover's own quotient has to satisfy, so unlike the three-wire system a prover cannot produce an
internally consistent proof for a public input its witness does not carry: its quotient stops
dividing.

## 3. Polynomial groups

fflonk commits a group of `t` polynomials as one interleaved polynomial

```
g(X) = Σ_{i<t} f_i(X^t) · X^i        i.e.  g[j·t + i] = f_i[j]
```

with the residue identity `g mod (X^t − ζ) = Σ_i f_i(ζ) X^i`, so one opening of `g` certifies every
`f_i(ζ)` at once and the residue's coefficients *are* those evaluations.

Eight groups, in this exact order (the ordering is protocol; it is the order the `ν` powers apply):

| group | `t` | columns | opened at | committed |
|---|---:|---|---|---|
| `P0` | 7 | precomputed `0..6` | `{ξ}` | preprocessing (in the VK) |
| `P1` | 7 | precomputed `7..13` | `{ξ}` | preprocessing (in the VK) |
| `P2` | 7 | precomputed `14..20` | `{ξ}` | preprocessing (in the VK) |
| `P3` | 7 | precomputed `21..27` | `{ξ}` | preprocessing (in the VK) |
| `W` | 3 | `w_l, w_r, w_o` | `{ξ, ξω}` | round 1 |
| `M` | 3 | `lookup_read_counts, lookup_read_tags, w_4` | `{ξ, ξω}` | round 2 |
| `Z` | 4 | `lookup_inverses, z_perm, S_0, S_1` | `{ξ, ξω}` | round 3 |
| `T` | 6 | `t_0 .. t_5` | `{ξ}` | round 4 |

with vanishing polynomials `Z_g = X^t − ξ` for the one-point groups and
`Z_g = (X^t − ξ)(X^t − ξω)` for the two-point ones.

The precomputed columns are tiled seven at a time in the flavor's own order; nothing in this protocol
names an individual selector, so adding one to the flavor changes only `NUM_PRECOMPUTED` and the
tiling.

**Why these groups.** Round boundaries are forced: a group can only hold columns that are committed
before the same challenge. Within a round the choice is between fewer, wider groups (cheaper
verifier, higher prover degree and a larger SRS) and more, narrower ones. Splitting the 28
preprocessed columns four ways rather than committing them as one group of 28 drops the widest packed
degree from `28n` to `7n` — and the SRS with it — for three extra scalar multiplications.

**Two-point groups reveal every member at both points.** `lookup_read_counts`, `lookup_read_tags` and
`lookup_inverses` are not read shifted by any relation; they are revealed at `ξω` only because they
share a group with a column that is. That costs three field elements of proof and one blinder each,
and it is cheaper than the extra group that separating them would need. It is also *checkable*:
`EveryEvaluationTheRelationsReadIsConstrained` asserts that those three evaluations, and no others,
leave the quotient identity unmoved.

### Residues of a two-point group

With `A(Y) = Σ_i f_i(ξ) Y^i` and `B(Y) = Σ_i f_i(ξω) Y^i`, the residue of `g` modulo
`(X^t − ξ)(X^t − ξω)` is the Chinese-remainder interpolant

```
R(Y) = [ B(Y)(Y^t − ξ) − A(Y)(Y^t − ξω) ] / (ξω − ξ)
```

At `t = 1` this is the straight line through `(ξ, f(ξ))` and `(ξω, f(ξω))`, which is what the
three-wire system uses for its grand product. The prover divides in two steps:
`g = q₁(X^t − ξ) + A`, then `q₁ + (A − B)/(ξω − ξ)` is divisible by `X^t − ξω`, exactly — because
`q₁ mod (X^t − ξω) = (B − A)/(ξω − ξ)`.

## 4. The trace-sum subrelations

A quotient argument can say "this vanishes at every row". It cannot say "this sums to zero over the
rows" — that is a different and strictly weaker claim, and dividing by `Z_H` would prove the wrong
thing.

So each linearly dependent subrelation `f_j` gets a committed running sum `S_j`, defined by

```
S_j[0] = 0,   S_j[i+1] = S_j[i] + f_j(ω^i)
```

and constrained by the single identity

```
S_j(ωX) − S_j(X) − f_j(X) = 0     for all X ∈ H
```

Summing that over `H` telescopes around the cycle — this is where the shift being *cyclic* is
essential — and leaves `Σ_i f_j(ω^i) = 0`. So the mere existence of an `S_j` satisfying the identity
is the statement Sumcheck was making, and nothing further is needed: in particular `S_j` is **not**
pinned to start at zero, because soundness comes from the sum and not from the starting value. The
freedom in the additive constant is harmless and the prover fixes it at zero.

Both `S_j` are committed in round 3, alongside `z_perm`: their summands depend on `η`,
`rom_logup_gamma`, `β` and `γ`, all of which are known by then, and on nothing later. In particular
they do not depend on `α`, which is what keeps them out of the quotient group.

## 5. The protocol

`vk_hash = keccak256(n ‖ ℓ ‖ pub_inputs_offset ‖ P0.x ‖ P0.y ‖ … ‖ P3.x ‖ P3.y)`, every field a
32-byte big-endian word.

The transcript is `state ← keccak256(state ‖ absorbed words) mod r`, starting from `state = 0`;
group elements absorb as `(x, y)`. This is `../fflonk/transcript.hpp` unchanged.

The rounds mirror `ultra_honk/oink_prover.cpp` exactly — one fflonk group per Oink round, so the
challenge schedule is Honk's:

```
absorb vk_hash, then x_0 … x_{ℓ−1}
round 1   prover sends C_W  = commit(W)
          η ← squeeze,  rom_logup_gamma ← squeeze
round 2   memory records and ROM-LogUp inverses land in w_4
          prover sends C_M  = commit(M)
          β ← squeeze,  γ ← squeeze
round 3   lookup inverses, grand product and the two running sums are built
          prover sends C_Z  = commit(Z)
          α ← squeeze
round 4   prover sends C_T  = commit(T)
          ξ ← squeeze
round 5   prover sends the 54 evaluations, group by group, each group's ξ openings then its ξω ones
          ν ← squeeze
round 6   prover sends W = Σ_g ν^g (g_g − R_g) / Z_g
          y ← squeeze
round 7   prover sends W' = L / (X − y),  where
              L(X) = Σ_g ν^g Γ_g (g_g(X) − R_g(y)) − Z_T(y)·W(X)
              Z_T = Π_g Z_g,  Γ_g = Z_T(y)/Z_g(y)
```

**Proof:** `C_W, C_M, C_Z, C_T, W, W'` and 54 field elements — 6 G1 points and 54 scalars,
**2,112 bytes**, plus 32 bytes per public input.

## 6. Verification

Reject unless every scalar is `< r`, every point is on the curve, `ξ ∉ H`, `ξ ≠ 0`, and no `Z_g(y)`
is zero. Then, writing `e` for the claimed evaluations assembled back into the flavor's per-row
container — unshifted entities from the `ξ` openings, the five shifted entities from the `ξω`
openings of `W`, `M` and `Z` — check the single quotient identity

```
  Σ_{i ∈ per-row} α^i · R_i(e)
+ Σ_{j<2}        α^{31+j} · ( S_j(ξω) − S_j(ξ) − f_j(e) )
=  ( Σ_{c<6} ξ^{cn} · t_c )  ·  Z_H(ξ)
```

where the `α` power advances across *every* subrelation in the flavor's relation order, including the
two that are held back for a running sum, so a subrelation's power depends only on that order.

Individual relations do **not** vanish at `ξ` — only the batch does. That is the substantive
difference from the three-wire system's three separate identities, and it is why the per-identity
sensitivity test there becomes a per-evaluation one here.

Then the batched opening, with `s_g = ν^g · Z_T(y)/Z_g(y)`:

```
F = Σ_g s_g·C_g  −  (Σ_g s_g·R_g(y))·[1]  −  Z_T(y)·W  +  y·W'
e(F, [1]_2) = e(W', [x]_2)
```

Group work: eight group commitments, `[1]`, `W` and `W'` — **eleven scalar multiplications and one
pairing**, whatever the circuit contains. Field work is one evaluation of the nine relations.

### Soundness of the product `Z_T`

BDFG21 defines `Z_T` as the vanishing polynomial of the *union* of the point sets; this uses the
product, which double-counts when two groups share a point. The argument in `../fflonk/PROTOCOL.md`
§8 carries over verbatim: `ν` is drawn after the claimed evaluations are fixed, the pairing forces
`Q(y) = 0` for random `y` hence `Q ≡ 0`, and the residue of `Σ_g ν^g (f_g − R_g)/Z_g` at any pole `s`
is a nonzero polynomial in `ν` unless `f_g(s) = R_g(s)` for every group containing `s`. Overlap
changes which terms appear in the residue; it does not let them cancel.

## 7. Degrees, the SRS and what it costs

Derived here rather than assumed, and re-checked in code: `prover.cpp` trims the realised quotient
and aborts if it does not fit the six chunks.

**Column degrees.** Precomputed columns interpolate `n` values, so `deg ≤ n−1`. Every witness column
carries two blinders — `p + (b_0X + b_1)Z_H`, which changes nothing on `H` — so `deg ≤ n+1`.

**Numerator.** `UltraFlavor::MAX_PARTIAL_RELATION_LENGTH == 7`, i.e. the maximum subrelation degree
in the entities is 6. A degree-6 monomial in witness columns reaches `6(n+1) = 6n+6`; the
running-sum identities are lower. So

```
deg N ≤ 6n + 6
```

**Quotient domain.** The inverse FFT that recovers `N`'s coefficients is exact only if
`deg N < D·n`. `6n+6 < 8n` for `n > 3`, and `6n+6 ≥ 4n` always, so `D = 8` is the smallest power of
two that works — the same factor the three-wire system needs, for a different reason.

**Quotient split.** `deg T ≤ 5n+6`, which is `5n+7` coefficients, so six chunks of `n` suffice for
`n ≥ 7`. `NUM_QUOTIENT_CHUNKS = 6`.

**Minimum circuit size** is therefore 8 — the next power of two above the `n > 3` and `n ≥ 7` bounds,
and enough to hold the four disabled rows at the top of an Ultra trace.

**Packed degrees and the SRS.** A group of `t` columns of degree `d` packs to degree `t·d + t−1`:

| group | packed degree |
|---|---|
| `P0..P3` | `7(n−1) + 6 = 7n − 1` |
| `W`, `M` | `3(n+1) + 2 = 3n + 5` |
| `Z` | `4(n+1) + 3 = 4n + 7` |
| `T` | `6n + 5` |

`W` and `W'` inherit the largest, so **the SRS must hold `7n` points**. The implementation asks for
`8n`, matching the quotient domain, which is one allocation rather than two bounds to keep in step.

**Prover cost.** Per proof: `36` inverse FFTs of size `n` to reach coefficient form, `38 × 8`
forward FFTs of size `n` to reach the quotient domain one coset at a time, `8n` evaluations of the
nine relations, one inverse FFT of size `8n`, and six MSMs totalling `3n + 3n + 4n + 6n + 7n + 7n =
30n` points. Peak live memory is dominated by `38n` field elements of coset evaluations plus `38n` of
coefficients plus `8n` for the numerator.

Working one coset at a time is what keeps that at `38n` rather than `38·8n`: the quotient domain is
eight cosets of `H`, and multiplying a coset point by `ω` moves it one position along the *same*
coset — so `p(ωX)` on a coset is that coset's evaluations rotated by one, and the shifted entities
cost no FFTs at all.

## 8. Zero knowledge

Twenty-one blinding scalars: two for each of the eleven committed witness columns, and five in the
kernel of the quotient split.

Every witness group is opened at exactly two points, so every witness column is revealed at exactly
two — and two blinders make a pair of evaluations off `H` uniformly distributed, by the same 2×2
determinant argument as `../fflonk/PROTOCOL.md` §7: the matrix `[[ξ Z_H(ξ), Z_H(ξ)], [ξω Z_H(ξω),
Z_H(ξω)]]` has determinant `Z_H(ξ)Z_H(ξω)(ξ − ξω) ≠ 0`.

The six quotient chunks are revealed at one point each. Five are free and one is pinned by the
quotient identity given everything else, so five blinders suffice; they are placed inside the kernel
of the split map — `t_c += b_c X^n`, `t_{c+1} −= b_c` — which leaves `Σ_c X^{cn} t_c` unchanged as a
polynomial.

Like the three-wire system, **this is argued rather than proven mechanically**. The argument is that
each opening is either uniformly random given its blinders or determined by the verifier's own
identity from openings that are; what is not established here is a simulator.

## 9. What is not here

- **No Solidity verifier.** `../fflonk/solidity/` does not carry over: this scheme's on-chain
  verifier has to evaluate the nine Ultra relations at `ξ`, which is a different contract and the
  bulk of the remaining work. Until it exists, the gas figure for this scheme is an estimate and is
  not reported as a measurement.
- **The prover has not been optimised.** The coset-at-a-time quotient is the memory-shaped choice,
  not the fast one; the relation evaluation at `8n` rows is the dominant term and is a straight port
  of the row-wise accumulators rather than anything vectorised.
- **`lookup_read_counts`, `lookup_read_tags` and `lookup_inverses` are revealed at `ξω` for
  nothing.** Splitting them into their own one-point group would trade one scalar multiplication for
  three field elements and three blinders; the current layout takes the other side of that trade.
