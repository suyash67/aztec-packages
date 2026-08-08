# fflonk over a plonkish arithmetization

A complete SNARK — circuit builder, preprocessing, prover, native verifier, Solidity verifier —
whose only purpose is to be verified on Ethereum as cheaply as possible.

This is the specification. The C++ prover (`prover.cpp`), the C++ verifier (`verifier.cpp`) and the
Solidity verifier (`solidity/`) are three independent implementations of what is written here; where
they disagree, this document is what needs fixing first.

## 1. Why this exists

A project proving with a hash-based commitment scheme cannot put that proof on L1 — the verifier is
hundreds of thousands of gas of hashing at best, and unbounded at worst. The standard escape is to
wrap: prove "I know a valid inner proof" in a second circuit whose own verifier is cheap. Today that
wrapper is almost always Groth16, which is cheap (~200k gas) but needs a **circuit-specific trusted
setup** — every change to the wrapper circuit means a new ceremony.

fflonk (Gabizon–Williamson, [ePrint 2021/1167](https://eprint.iacr.org/2021/1167)) is the same order
of gas with a **universal** setup: the same powers-of-tau works for every circuit. That is the whole
trade this module implements.

### Why not the fflonk *opening* under Honk

`commitment_schemes/fflonk/` already implements fflonk as a PCS backend for UltraHonk, and
`commitment_schemes/recursion/README.md` §7 measures it: 519,164 gas of execution, against 630,236
for Shplemini+KZG. The floor there is Honk itself — 155,803 gas of sumcheck and relation evaluation
that no commitment scheme can touch, plus ~36 column commitments to batch.

fflonk's ~200k reputation is a *PlonK* number, and it comes from PlonK's linearisation: a PlonK
verifier's only group work is a small multi-scalar multiplication and one pairing. That is what this
module builds, and it is why it is a separate proof system rather than another PCS backend.

The specific trick that makes it cheap here: because every preprocessed polynomial's evaluation at
the challenge is *opened* rather than reconstructed in the group, the verifier checks all three
constraint identities **entirely in the field**. There is no linearisation MSM at all — the group
work is one batched opening, five scalar multiplications and one pairing, independent of the number
of selectors.

## 2. Arithmetization

Fix `n = 2^k`, `ω` a primitive `n`-th root of unity, `H = {ω^i}`, `Z_H(X) = X^n − 1`.

A circuit is `n` rows of

```
q_M(x)·a(x)·b(x) + q_L(x)·a(x) + q_R(x)·b(x) + q_O(x)·c(x) + q_C(x) + PI(x) = 0   for all x ∈ H
```

plus copy constraints: an arbitrary partition of the `3n` wire slots into classes that must hold
equal values.

**Public inputs.** The `ℓ` public inputs occupy rows `0 … ℓ−1` on wire `a`, each row carrying
`q_L = 1` and every other selector zero, so the row reads `x_i + PI(ω^i) = 0`. Accordingly

```
PI(X) = − Σ_{i<ℓ} x_i · L_i(X),        L_i(X) = ω^i (X^n − 1) / (n (X − ω^i))
```

**Copy constraints.** Slot `(j, i)` — wire `j ∈ {0,1,2}`, row `i` — is identified with the field
element `k_j · ω^i`, where `k_0 = 1` and `k_1`, `k_2` are chosen so that `H`, `k_1 H`, `k_2 H` are
pairwise disjoint (asserted at preprocessing time, which is what makes the identification
injective). `σ` is the permutation of the `3n` slots that cycles each equality class, and
`S_σj(ω^i) = id(σ(j, i))`.

## 3. Polynomial groups

fflonk commits a *group* of `t` polynomials as one interleaved polynomial

```
g(X) = Σ_{i<t} f_i(X^t) · X^i        i.e.  g[j·t + i] = f_i[j]
```

The identity this module leans on throughout is

```
g mod (X^t − ζ) = Σ_{i<t} f_i(ζ) · X^i
```

so a single opening of `g` on the `t`-th roots of `ζ` certifies **every** `f_i(ζ)` at once, and the
residue's coefficients *are* those evaluations. Nothing ever has to construct a `t`-th root of `ζ`,
take a root of unity of order `t`, or run an inverse DFT — which is what lets `t` be any positive
integer, not just one for which `ζ^{1/t}` is guaranteed to exist.

Four groups, in this exact column order (the ordering is part of the protocol; the Solidity verifier
hard-codes it):

| group | `t` | columns | opened at | committed |
|---|---:|---|---|---|
| `G0` | 8 | `q_L, q_R, q_O, q_M, q_C, S_σ1, S_σ2, S_σ3` | `{ξ}` | preprocessing (in the VK) |
| `G1` | 5 | `a, b, c, T0_lo, T0_hi` | `{ξ}` | round 1 |
| `G2` | 1 | `z` | `{ξ, ξω}` | round 2 |
| `G3` | 4 | `T1, T2_lo, T2_mid, T2_hi` | `{ξ}` | round 2 |

with vanishing polynomials

```
Z_0(X) = X^8 − ξ      Z_1(X) = X^5 − ξ      Z_2(X) = (X − ξ)(X − ξω)      Z_3(X) = X^4 − ξ
```

`z` sits alone in `G2` because it is the only polynomial that has to be opened at two points, and
grouping anything with it would reveal that thing at `ξω` too — where no verifier identity pins it
down, so it could not be simulated. See §7.

## 4. Quotients

Three separate quotients, one per identity, rather than one `α`-combined quotient. This costs
nothing (they are opened anyway, inside groups that already exist) and it removes `α` from the
protocol entirely — which is what lets `G3` be committed in the same round as `G2`.

```
T0 = [ q_M·a·b + q_L·a + q_R·b + q_O·c + q_C + PI ] / Z_H                            deg ≤ 2n+1
T1 = [ (z − 1)·L_0 ] / Z_H  =  (z − 1) / (n·(X − 1))                                 deg ≤ n+1
T2 = [ f(X)·z(X) − g(X)·z(ωX) ] / Z_H                                                deg ≤ 3n+5
    f = (a + βX + γ)(b + βk_1X + γ)(c + βk_2X + γ)
    g = (a + βS_σ1 + γ)(b + βS_σ2 + γ)(c + βS_σ3 + γ)
```

`T1`'s closed form follows from `L_0(X) = (X^n − 1)/(n(X − 1))`; it needs no FFT, just one synthetic
division of `z − 1` by `X − 1` (exact because `z(1) = 1`).

`T0` is split as `T0 = T0_lo + X^n·T0_hi`, and `T2` as `T2 = T2_lo + X^n·T2_mid + X^{2n}·T2_hi`.

## 5. The protocol

`vk_hash = keccak256(n ‖ ℓ ‖ k_1 ‖ k_2 ‖ C0.x ‖ C0.y)`, every field a 32-byte big-endian word.

The transcript is `state ← keccak256(state ‖ absorbed words) mod r`, starting from `state = 0`.
Group elements absorb as `(x, y)`, with the point at infinity absorbing as `(0, 0)`.

```
absorb vk_hash, then x_0 … x_{ℓ−1}
round 1   prover sends C1 = commit(G1)
          β ← squeeze,  γ ← squeeze
round 2   prover sends C2 = commit(G2), C3 = commit(G3)
          ξ ← squeeze
round 3   prover sends the 19 evaluations, in this order:
              q_L, q_R, q_O, q_M, q_C, S_σ1, S_σ2, S_σ3      (G0 at ξ)
              a, b, c, T0_lo, T0_hi                          (G1 at ξ)
              z, zω                                          (G2 at ξ and ξω)
              T1, T2_lo, T2_mid, T2_hi                       (G3 at ξ)
          ν ← squeeze
round 4   prover sends W = Σ_i ν^i (g_i − R_i) / Z_i
          y ← squeeze
round 5   prover sends W' = L / (X − y),  where
              L(X) = Σ_i ν^i Γ_i (g_i(X) − R_i(y)) − Z_T(y)·W(X)
              Z_T = Π_i Z_i,  Γ_i = Z_T(y)/Z_i(y)
```

`R_i` is the residue of `g_i` modulo `Z_i`, i.e. the interpolant of the claimed evaluations:

```
R_0(Y) = q_L + q_R Y + q_O Y² + q_M Y³ + q_C Y⁴ + S_σ1 Y⁵ + S_σ2 Y⁶ + S_σ3 Y⁷
R_1(Y) = a + b Y + c Y² + T0_lo Y³ + T0_hi Y⁴
R_2(Y) = A + B Y,     B = (zω − z)/(ξω − ξ),   A = z − ξB
R_3(Y) = T1 + T2_lo Y + T2_mid Y² + T2_hi Y³
```

**Proof:** `C1, C2, C3, W, W'` and 19 field elements — 5 G1 points and 19 scalars, 928 bytes, plus
32 bytes per public input.

## 6. Verification

Reject unless every scalar is `< r` and every point is on the curve. Then, with
`Z_H(ξ) = ξ^n − 1` and `L_0(ξ) = Z_H(ξ)/(n(ξ − 1))`, reject unless `ξ ∉ H` and `ξ ≠ 0`, and check

```
(gate)    q_L a + q_R b + q_O c + q_M a b + q_C + PI(ξ)  =  (T0_lo + ξ^n T0_hi) · Z_H(ξ)
(start)   (z − 1) · L_0(ξ)                               =  T1 · Z_H(ξ)
(perm)    (a + βξ + γ)(b + βk_1ξ + γ)(c + βk_2ξ + γ)·z
            − (a + βS_σ1 + γ)(b + βS_σ2 + γ)(c + βS_σ3 + γ)·zω
                                                         =  (T2_lo + ξ^n T2_mid + ξ^{2n} T2_hi) · Z_H(ξ)
```

all three in the field, over the *claimed* evaluations. Then the batched opening, with
`s_i = ν^i · Z_T(y)/Z_i(y)`:

```
F = Σ_i s_i·C_i  −  (Σ_i s_i·R_i(y))·[1]  −  Z_T(y)·W
e(F + y·W', [1]_2) = e(W', [x]_2)
```

Group work: five scalar multiplications (`C0`, `C1`, `C2`, `C3`, `[1]` — `W` and `W'` fold in as
one more each, so seven in all) and one pairing, whatever the circuit contains.

## 7. Zero knowledge

Twelve blinding scalars. `a, b, c` take two each, `z` takes three, and the quotient splits are
randomised inside the kernel of the split map — `T0_lo += b_10 X^n`, `T0_hi −= b_10`, and
`T2_lo += b_11 X^n`, `T2_mid −= b_11 + b_12 X^n`, `T2_hi −= b_12` — which leaves `T0` and `T2`
unchanged as polynomials.

What a verifier learns, and why each is simulatable:

| revealed | count | pinned by | free | blinders |
|---|---:|---|---:|---:|
| `a, b, c` at `ξ` | 3 | — | 3 | 6 |
| `T0_lo, T0_hi` at `ξ` | 2 | the gate identity | 1 | 1 |
| `z` at `ξ, ξω` | 2 | — | 2 | 3 |
| `T1` at `ξ` | 1 | the start identity, given `z(ξ)` | 0 | — |
| `T2_lo, T2_mid, T2_hi` at `ξ` | 3 | the permutation identity | 2 | 2 |

Every group is opened at one point except `G2`, which holds only `z` — and `z`'s three blinders
cover its two openings because the 2×2 matrix `[[ξ Z_H(ξ), Z_H(ξ)], [ξω Z_H(ξω), Z_H(ξω)]]` has
determinant `Z_H(ξ)Z_H(ξω)(ξ − ξω) ≠ 0`.

This is why `z` is not grouped with the quotients, as it is in the classic 3-group fflonk layout.
There, `T2` is revealed at `ξω` as well, where no verifier identity determines it and the split
kernel's evaluations at two points span only a 4-dimensional space against the 5 that need hiding.
Separating `z` costs one scalar multiplication and 64 bytes; it buys a ZK argument that closes.

## 8. Why the batched opening is sound with overlapping point sets

BDFG21 defines `Z_T` as the vanishing polynomial of the *union* of the point sets. This module uses
the product `Z_T = Π_i Z_i`, which double-counts if two groups share a point — and they can:
`S_0 ∩ S_1 ≠ ∅` exactly when `ξ^9 = 1`, and `ξ ∈ S_0` exactly when `ξ^7 = 1`.

The product is still sound, because `ν` is drawn after the claimed evaluations are fixed. The
verifier's pairing check forces `Q(y) = 0` for random `y`, where

```
Q(X) = Σ_i ν^i (Z_T/Z_i)(X)·(f_i(X) − R_i(X)) − Z_T(X)·W(X)
```

so `Q ≡ 0`, so `Σ_i ν^i (f_i − R_i)/Z_i = W` as rational functions — the left side has no poles. At a
point `s`, the residue of the left side is `Σ_{i: s ∈ S_i} ν^i (f_i(s) − R_i(s))/Z_i'(s)`, a nonzero
polynomial in `ν` unless `f_i(s) = R_i(s)` for every group containing `s`. Since `ν` is uniform and
independent of those values, that happens with probability at most `deg/r`. Overlap changes which
terms appear in the residue; it does not let them cancel.

The verifier still has to reject `Z_i(y) = 0` — checked explicitly, since it would otherwise be a
division by zero rather than a soundness break.

## 9. Degrees and the SRS

With blinding: `deg a, b, c ≤ n+1`, `deg z ≤ n+2`, so the `T2` numerator reaches `4n+5` and the
quotient FFT runs on a domain of size `8n`. Packed degrees are `8n−1`, `5n+9`, `n+2`, `4n+23`, and
`W`, `W'` inherit the largest, so **the SRS must hold `8n` points** — the price fflonk charges for
folding a group of `t` columns into one commitment.

Per proof the prover runs five MSMs totalling `5n + n + 4n + 8n + 8n = 26n` points, which is the
same order as UltraHonk's ~17`n`; `G0` is preprocessed and does not recur.

The packing factor is a dial. Splitting `G0` into two groups of four drops every degree to `5n` and
the SRS with it, for one more scalar multiplication (~6k gas) in the verifier. This module takes the
cheap-verifier end of that trade, because the wrapper is proven once and verified forever.
