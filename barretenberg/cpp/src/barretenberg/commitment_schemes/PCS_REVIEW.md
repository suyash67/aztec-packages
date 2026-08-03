# Polynomial commitment schemes in barretenberg: correctness, efficiency, and the pluggable-backend question

A review of every polynomial commitment scheme in `commitment_schemes/`, written to answer one
question: **can a Noir-like frontend expose the proof system's commitment scheme as a user-selectable
backend?** The answer is developed in §7; §§1–6 are the evidence.

Scope: 15 end-to-end UltraHonk backends, the shared shell they plug into, the two support protocols
(`small_subgroup_ipa`, `triple_ipa`), and the small-field study. Every scheme was read against its
source paper. Measurements are from this tree on an Apple M4 Pro (arm64, `DISABLE_ASM=1`); see §4.4
for the methodology, which uses same-run ratios rather than absolute times.

Status of the tree at review time: `commitment_schemes_tests` is **246 passing, 2 skipped** (the two
skips are pre-existing ECCVM translation-masking tests, disabled before this review).

**The tree moved under this review.** `BrakedownHonk`, the `QuerySegment` code-policy hook, and a full
`Bolt` backend all landed while it was being written. Where that changes a finding, the report says so
rather than describing a tree that no longer exists — §3.2 in particular.

---

## 0. Summary of findings

**Correctness.** No soundness break was found in the protocol logic of any backend. Every scheme
matches its paper's algebraic structure, the Fiat–Shamir orderings are correct (every challenge is
drawn strictly after the messages it must bind), and the three different Honk-shift constructions are
each sound. Seven findings are recorded in §3, of which the material ones are:

- **Four backends are not deployable as written, by construction, not by bug** (§3.1). KZH2, KZH3,
  CHOPIN and Dory ship *test-only* setups whose trapdoors are either sampled in-process (KZH2/KZH3),
  fixed to a public constant (CHOPIN's `τ_Y = 2`), or known discrete logs of the G2 generator
  (Dory). Each is clearly labelled in the source; none is a mistake. But it means their measured
  numbers are performance models, not security claims.
- **The Reed–Solomon query count rests on a conjecture that has since been disproved** (§3.2), the same
  one WHIR carries a repaired preset for. The shared query-plan hook this needed landed mid-review; the
  repaired RS rule behind it has not, and SwitchFold is not wired to the hook.
- **Ligero-family batching doubles as the proximity test, which is only sound because `u` is a
  Fiat–Shamir challenge** (§3.3). Correct inside `TransparentHonk`; unsound if the standalone API is
  called at an attacker-chosen point, which nothing in the API documents or enforces.
- **Brakedown's base-case code is not MDS** (§3.4), contradicting the stated property its distance
  argument relies on. A Cauchy matrix fixes it in one function.

**Efficiency.** A single shared defect dominated the whole experimental suite: **the ρ-batched random
linear combination of the ~41 Honk claims was written as a serial nested loop in every non-production
backend**, while production Gemini does it with the parallel `add_scaled_batch`. Fixing that, plus
six related serial loops and one algorithmically wrong IPA inner loop, gives the measured improvements
in §4.4. Headline, all with the proof bytes and the verifier untouched: **the IPA backend is 13.3x
faster**, **CHOPIN now proves faster than the KZG baseline** (0.81x, previously 1.05x) while producing
a smaller proof, **Vela's prover penalty drops from 1.70x to 1.18x**, and **SwitchFold improves 16%**.

**Zero-knowledge is the real gap.** Only the production Shplemini path and WHIR have zk. The other
thirteen backends are non-hiding: they leak the witness outright (§5). This is the single biggest
obstacle to the product, and it is not uniform work — each family needs a different masking design.

**Verdict (§7).** The pluggable-backend product is *architecturally* real and already demonstrated:
`TransparentHonk<Pcs>` proves that everything above the commitment layer — arithmetization,
relations, sumcheck, Fiat–Shamir schedule — is genuinely PCS-independent, and 15 backends ride it
unmodified. What is not yet real is the *security surface*: zk on 13 backends, ceremonies or
hash-to-curve for 4, the repaired query rule, and recursive verifiers. Concrete roadmap in §7.3.

---

## 1. The PCS slot: what a backend actually has to supply

### 1.1 The interface

`transparent_honk.hpp` defines the contract. A backend supplies eight things:

| Concept | Meaning |
|---|---|
| `Config` | Size/security parameters, derived identically by both parties, never transmitted |
| `CommitmentKey` | Constructed from `Config`; supplies `commit_group(columns, to_be_shifted)` |
| `GroupData` | Prover-side commitment state for one round's columns |
| `GroupCommitment` | What the verifier holds: a Merkle digest, a `vector<G1>`, or a `vector<GT>` |
| `ProverClaims` / `VerifierClaims` | `{group, column}` references plus evaluations |
| `prove_opening` / `verify_opening` | The opening argument at one multilinear point |
| send/receive/absorb hooks | Transcript plumbing for the commitment type |
| `payload_variables` | The polynomial log-size the config expects |

Everything else — the 28 precomputed columns, the Ultra relation set, the memory-record and
log-derivative derivations, the grand product, sumcheck — is shared verbatim with production
`OinkProver`/`SumcheckProver`. **This is the load-bearing evidence for the product question**: the
seam between "proof system" and "commitment scheme" is real and narrow.

### 1.2 The claim set a backend must open

After sumcheck at challenge `u ∈ F^n`, UltraHonk hands the PCS **41 claims over 36 committed
columns**: 36 unshifted evaluations `P_i(u) = v_i`, and 5 shifted evaluations `shift(P_j)(u) = w_j`
where `shift(P)[i] = P[i+1]`. Columns committed in the same Oink round share one group (5 groups:
28 precomputed, 3 wires, 3 counts+w4, 1 lookup_inverses, 1 z_perm).

### 1.3 bb's array-is-both-bases convention

A polynomial is one dense array of length `N = 2^n`, read simultaneously as

- the **evaluation table** of the multilinear extension over `{0,1}^n`, and
- the **coefficient vector** of a univariate of degree `< N`.

The stride-2 Lagrange fold `a'[j] = (1-α)a[2j] + α·a[2j+1]` is *simultaneously* partial multilinear
evaluation at the lowest variable and the univariate even/odd fold. No basis conversion ever happens.
Gemini, WHIR's folding, Mercury's twin, Vela's Laurent twin and SwitchFold's descent all exploit this;
it is why a 14-backend suite can share one data layout and one `eq_tensor`.

### 1.4 Three ways to handle the Honk shift

Every backend has to open `shift(P)` given a commitment to `P`, and the suite contains three genuinely
different solutions. This is the part of the design space that a pluggable-backend product must get
right, because it is where a new backend most easily goes wrong.

**(a) Univariate division contract** — Gemini/Shplemini, WHIR, Mercury, Vela.
The univariate twin of `shift(P)` is exactly `P(X)/X`, an *exact* division because the Honk shift
contract guarantees `P[0] = 0`. Mercury multiplies its fold identity through by `X`
(`B(X) = Q(X)(X^b − α) + X·g_B(X)` with `Q := X·q_B`) so the verifier can use the *unshifted*
homomorphic commitment. Vela folds the shifted chain into the same Laurent twin as
`g(X) = f_A(X) + λ·X^{-1}f_B(X)`. Precondition: `P[0] = 0`, asserted at prove time.

**(b) Rank-2 tensor decomposition** — Ligero, SwitchFold, Hyrax, KZH2, Dory (rank 3 for KZH3).
The shifted-eq vector is not a tensor, but it is a rank-2 sum of tensors. With `j = rC + c`:

```
Q'[r][c] = b[r]·a[c-1]        (c ≥ 1)
Q'[r][0] = b[r-1]·a[C-1]      (r ≥ 1)
Q'[0][0] = 0
i.e.  Q' = b ⊗ shr(a) + shr(b) ⊗ (a_{C-1}·e_0)
```

Costs one extra combined row (`w_s2`), of which only entry 0 is used. **No zero-constant-term
precondition** — the `(0,0)` term is dropped structurally. I verified the rank-3 extension in
`kzh3.hpp` term by term: the shift of the flattened index `k = (i₁d₂+i₂)d₃+i₃` crosses both the `d₃`
and `d₂` boundaries, contributing `T3` under shifted-`a₃` weights, `a₃[d₃−1]·T3'[0]`, and
`a₂[d₂−1]a₃[d₃−1]·T3''[0]`. Correct.

**(c) Successor-kernel closed form** — Pedersen IPA, Dory's scalar side.
The IPA fold reduces `b` to a single scalar `Σ_i b_i·Π_{j∈bits(i)} m_j`. For `b = eq(u)` that is
`Π_j((1−u_j) + u_j m_j)`; for the shifted-eq vector it telescopes to

```
Σ_k [Π_{j<k} u_j] (1−u_k) m_k [Π_{j>k} ((1−u_j) + u_j m_j)]
```

(derived by splitting on the lowest set bit of the index). The verifier never touches a shifted
commitment. I re-derived both closed forms independently; both are correct.

---

## 2. The backends

Fifteen, grouped by trust model and the resource each optimizes. Prove times are given as a **ratio to
the KZG baseline** (see §4.4 for why); verify and proof size are absolute, at 2^18, and unaffected by
this review's changes.

| Backend | Family | Trust | Prove (×KZG) | Verify | Proof | zk |
|---|---|---|---|---|---|---|
| Gemini+Shplonk+KZG | pairing | SRS ceremony | 1.00x | 4.7 ms | 13.3 KiB | **yes** |
| Mercury | pairing, constant-size | SRS ceremony | 1.14x | 4.2 ms | 9.4 KiB | no |
| Vela | pairing, reciprocal | SRS ceremony | 1.18x | **3.5 ms** | **7.8 KiB** | no |
| CHOPIN | pairing, bivariate KZG | **test-only SRS** | **0.81x** | 4.7 ms | 9.2 KiB | no |
| KZH2 | pairing, √N | **test-only SRS** | ~14x | 216 ms | 183.3 KiB | no |
| KZH3 | pairing, ∛N | **test-only SRS** | ~46x | 79.7 ms | 55.3 KiB | no |
| Dory | pairing, two-tier | **test-only G2** | ~12x | 362 ms | 94.3 KiB | no |
| IPA (Pedersen) | DL | transparent | ~5.6x | 3596 ms | 16.4 KiB | no |
| Hyrax | DL, Pedersen rows | transparent | ~11x | 243 ms | 566.3 KiB | no |
| WHIR | hash, RS folding | transparent | 1.08x | 13.5 ms | 1219.7 KiB | **yes** |
| SwitchFold | hash, code switching | transparent | 2.50x | 13.6 ms | 772.9 KiB | no |
| Ligero | hash, tensor √N | transparent | **0.72x** | 20.6 ms | 2692.6 KiB | no |
| BrakedownHonk | hash, linear-time code | transparent | 1.9x Ligero @2^14 | 21.6x Ligero | 20.2 MiB @2^14 | no |
| Bolt | hash, sketched LDPC | transparent | see `bolt/README.md` | — | — | no |

Support modules, not backends: `small_subgroup_ipa` (the zk-sumcheck and ECCVM-translation inner
product), `triple_ipa` (ECCVM's three-claim reduction on Grumpkin), `small_field` (M31/QM31 study).

### 2.1 Gemini + Shplonk + KZG (Shplemini) — the production baseline

**Math.** Gemini reduces `k` multilinear claims at `u` to `d+1` *univariate* claims. Batch
`A₀ = F + G/X` with `F = Σρ^i f_i`, `G = Σρ^{k+j} g_j`. Fold `A_{l+1}(X) = (1−u_l)even(A_l)(X) +
u_l·odd(A_l)(X)`; commit `A_1…A_{d−1}`. Open `A₀₊` at `r`, `A₀₋` at `−r`, and `A_l` at `−r^{2^l}`.
The verifier reconstructs the positive evaluations from the Gemini relation, so only the negative ones
are sent. Shplonk collapses the `d+1` single-point claims to one point `z` via the quotient
`Q(X) = Σ_j ν^{j−1}(f_j(X) − v_j)/(X − x_j)`; KZG closes with one pairing. The verifier's whole job is
**one MSM** over `[Q], [f_1..f_n], [1]` plus a pairing.

**Honk fit.** Native: this *is* the Honk PCS. Shift via the univariate contract (a); the `1/X` becomes
the `r^{-1}` factor in the shifted Shplonk scalar.

**zk.** The only fully worked-out story in the tree: a sparse `gemini_masking_poly` with `2d` random
coefficients on the tail-halving support `{E−1,E−2, N/2,N/2−1, N/4,N/4−1, …}`, plus `small_subgroup_ipa`
for the Libra univariates. `SHPLEMINI_ZK_MASKING.md` gives the reduction and points at a Lean
development for the linear algebra. There is a real edge case handled explicitly: if the Gemini
challenge `r` lands in the small subgroup used by `small_subgroup_ipa`, prover evaluations leak, and
`gemini_impl.hpp` aborts rather than proving.

**Audit.** Production code, internally audited (`ipa.hpp` carries an audit banner; `gemini_impl.hpp`
one planned). The prover is well optimized: SIMD `fold_stride2`, extent tracking that skips the
provably-zero tail of each fold, `parallel_for_heuristic` throughout.

**Future work.** None specific to this review.

### 2.2 Mercury (ePrint 2025/385)

**Math.** With `N = 2^n`, `t = ⌈n/2⌉`, `b = 2^t`, the array reshapes to `rows × b`. The claim
`f̃(u) = v` factors through the column-combined polynomial `h(X) = Σ_j eq_j(u_lo)·f_j(X)` as four
sub-claims: (1) the binomial fold `f(X) = q(X)(X^b − α) + g(X)`, checked at random `ζ` through the
quotient `H = (f − (ζ^b−α)q − g(ζ))/(X−ζ)`; (2) the degree bound `deg g < b` via `D(X) = X^{b−1}g(1/X)`;
(3) `⟨P_{u_lo}, g⟩ = h(α)`; (4) `⟨P_{u_hi}, h⟩ = v`. The two inner products batch by `γ` into one
Laurent identity

```
Σ_k γ^k (p_k(X)r_k(1/X) + p_k(1/X)r_k(X)) = 2Σγ^k c_k + X·S(X) + X^{-1}·S(1/X)
```

whose constant term extracts the inner products. All small-polynomial openings batch into two group
elements with BDFG. Proof: 13 G1 + 14 F, **constant in `N` and in the claim count**.

**Audit.** I verified the Laurent accumulator coefficient-by-coefficient: `[X^k](p(X)r(1/X) +
p(1/X)r(X)) = Σ_{i≥k}(p_i r_{i−k} + p_{i−k} r_i)`, which is what `accumulate_ip_terms` computes, and
the verifier's `lhs == 2·constant + ζS(ζ) + ζ^{-1}S(1/ζ)` matches. **No degree bound on `S` is
needed** and none is present — correctly, because `X·S(X) + X^{-1}S(1/X)` has zero constant term for
*any* `S`, so the identity still forces the claim. The degree bound on `g` is sound: `ζ` is drawn
after `C_D`, and `D(X) = X^{b−1}g(1/X)` as Laurent polynomials forces `deg g ≤ b−1` because `D` has no
negative exponents. Shift handling is contract (a) with the `X`-multiplied fold identity; verified.

**zk.** Not implemented. The paper's §5 gives the design (bounded-degree blinders on `h`, `g`, `S`).

**Future work.** zk; recursive verifier.

### 2.3 Vela (ePrint 2026/1438)

**Math.** The multilinear claim is the constant coefficient of the Laurent product
`H(X) = f_v(X)·T_r(1/X)` where `T_r(X) = Π_k((1−r_k) + r_k X^{2^k})`. The inversion-symmetric residual
`D(X) = H(X) + H(1/X) − 2y` has zero constant coefficient **exactly when the claim holds**, and its
symmetry (`D_k = D_{−k}`) yields a single ordinary polynomial `h` with `D(X) = X·h(X) + X^{-1}h(1/X)`.
One commitment to `h` plus evaluations at `z` and `1/z` certify the claim — and the fourth evaluation
`h(1/z)` is **never transmitted**: the verifier recovers it as

```
w₁ = z·( g(z)T_r(1/z) + g(1/z)T_r(z) − 2Y − z·w₀ )
```

which is precisely the identity being proved, so a lying prover's recovered `w₁` fails the subsequent
two-point opening.

**Honk fit.** The whole batched claim set collapses into *one* identity: the shifted twin is exactly
`X^{-1}f_B(X)`, so the combined twin `g = f_A + λX^{-1}f_B` carries both chains under one `h`.

**Audit.** I checked the Laurent table's exponent range: `g` spans `[−1, N−1]`, `T_r(1/X)` spans
`[−(N−1), 0]`, so the product spans `[−N, N−1]` and the table covers exactly that — **no truncation**.
The in-place forward pass reads `table[e+2^k]` before writing it, so it was correct serially. The
recovery formula and the degree-2 linearization (forced by bb publishing only `[1]₂, [τ]₂`) are both
sound: the check `E(ζ) = 0` at random `ζ` drawn after `C_q` forces `F − R = Z·q` identically. The
prover's rejection of `z ∈ {0, 1, −1}` is necessary (the two-point interpolation divides by `z − 1/z`)
and mirrored on the verifier.

**Deviation from the paper** (already documented): the paper's two-point check needs `[τ²]₂`; bb has
only `[1]₂, [τ]₂`, so the degree-2 vanishing is linearized at one further challenge, costing one extra
G1 (3 G1 + 5 F here vs the paper's 2 G1 + 3 F unbatched).

**zk.** Not implemented.

**Future work.** zk; the `O(N log N)` Laurent expansion is inherent to never sending `h(1/z)` — see
§4 for its parallelization.

### 2.4 CHOPIN (ePrint 2026/480)

**Math.** Reads the claim as a bilinear form over the coefficient matrix of a *bivariate* twin
`f(X,Y)`, and opens it with a column restriction `w = F·Ψ_zR`, one α-fold, a γ-batched Lagrange IPA
(Mercury's Laurent accumulator, reused directly), a δ/z multi-polynomial multi-point univariate batch
proof (paper Fig. 7), and one bivariate KZG opening at `(α,β)`. The X-quotient of that last opening is
the **single** size-N MSM, against Mercury's fold quotient plus per-chain BDFG quotients.

**Honk fit.** The shift is a verifier-side query decomposition: chain B opens one extra column
restriction `w2 = F·shr(Ψ_zR)` plus its value at zero, and the batch proof additionally opens `w2` at
`X = 0`. The verifier's shifted weights use `shr_lo(β) = β·Ψ_lo(β) − eq_lo[M₁−1]·β^{M₁}`, the standard
shift-in-the-evaluation-domain identity; verified.

**Audit.** Protocol logic is faithful. The `t_points` set grows to 4 (`{α, β, 1/β, 0}`) when shifted
claims are present, and the prover and verifier compute the complement products consistently. The
three pairing identities (batch-proof `pi_s`, quotient `pi_q`, bivariate opening) are ξ-batched into
one 3-term multipairing against `[1]₂, [τ]₂, [2]₂`.

**Setup caveat (§3.1).** `τ_Y` is fixed to the public constant 2 so the grid derives from the ceremony
X-row by doublings and `[τ_Y]₂ = 2[1]₂`. That is a **complete break of the Y-direction binding** for an
adversary — the module says so. Honest-prover cost is unaffected, which is the point of the shortcut.

**Future work.** A real second trapdoor (i.e. a bivariate ceremony); zk; recursion.

### 2.5 KZH2 / KZH3 (ePrint 2025/144, Fig. 2 and App. C)

**Math.** KZH2 splits the array into `R × C`; the commitment is `C = Σ_{i,j} f_{i,j}(τ_i A_j)`. The
opening sends the row commitments `D_i = Σ_j f_{i,j}A_j` and the combined row `w = Σ_i b_i·f_{i,·}`.
The verifier checks (1) `e(C, V) = Π_i e(D_i, V_i)` with `V_i = τ_i V`; (2) `Commit_A(w) = Σ_i b_i D_i`;
(3) `⟨w, a⟩ = v`. Proof `O(√N)`, verify one `R+1` multipairing plus two MSMs.

KZH3 adds a dimension: slice commitments `D1` (dim 1), contracted row commitments `D2` (dim 2), and the
final `d₃`-vector, bound by a two-level pairing chain. Proof drops to `O(∛N)` and verification with it.

**Honk fit.** Rank-2 (KZH2) / rank-3 (KZH3) query decomposition, shift term-by-term verified above.

**Audit.** Both are faithful to the paper's structure. KZH3's extra chain-B data is properly bound:
`T3'` via `⟨T3', A⟩ = ⟨shr(a₂), D2⟩` against the *same* `D2` layer, and `T3''` via a second pairing
chain on `C1'' = ⟨shr(a₁), D1⟩` followed by `Commit_A(T3'') = D2''₀`.

**Setup caveat (§3.1).** Both sample their own trapdoors in-process. Deployment needs a ceremony.

**Cost.** KZH3's proof gain widens with size (3.3x at 2^18, 4.1x at 2^20 vs KZH2) but so does its
prover (1.6x at 2^14 → 4.8x at 2^20), because the slice-commitment layer's count grows with `d1`. The
crossover where the smaller proof pays for the prover sits well below 2^18. KZH-k exists to shrink an
*accumulator and decider*, which is a different objective from raw proving throughput.

### 2.6 Dory (ePrint 2020/1274, adapted)

**Math.** The commitment is AFGHO: `T = Π_i e(P_i, Γ2_i)` with `P_i` the Pedersen row commitments — one
GT element per polynomial, *transparent* (no trapdoor). The opening runs `log R` inner-pairing fold
rounds: each round sends the GT cross terms `T_eo`, `T_oe` of `T`'s even/odd split and the G1 cross
terms of every Q-claim, and folds `rows' = rows_even + x·rows_odd`, `Γ2' = Γ2_even + x^{-1}Γ2_odd`. The
final check is one pairing `e(P_final, Γ2_final) = T_acc` plus the scalar-side closed forms.

**Audit.** The fold algebra is consistent between prover and verifier
(`T' = T·T_eo^{x^{-1}}·T_oe^{x}`, `Q' = Q + xQ₊ + x^{-1}Q₋`), and the b-side uses the same
successor-kernel closed form as the IPA backend (contract (c)), which I verified.

**Setup caveat (§3.1).** `Γ2_i = s_i·G2` with `s_i` public deterministic scalars. bb has no
cofactor-cleared hash-to-G2, so the module substitutes known discrete logs — **which destroys AFGHO
binding**: knowing the `s_i` lets an adversary find GT collisions. Documented as test-only. Real
deployment needs hash-to-G2, at identical runtime cost.

**Cost.** The commitment does `R` Miller loops *per polynomial*: at 2^18 with 36 columns that is ~18k
pairings, and it is the whole prover cost. This is inherent to AFGHO, not an implementation artifact.

**Future work.** Hash-to-G2; GT-Pippenger + cyclotomic squaring for the verifier's ρ-combination
(§4.4); zk.

### 2.7 IPA over Pedersen generators (Bulletproofs)

**Math.** The classical inner-product argument over hash-derived generators, run directly on the
multilinear array. Per chain, `log n` rounds: `L = ⟨a_hi, G_lo⟩ + ⟨a_hi, b_lo⟩U`, symmetrically `R`;
challenge `x` folds `a' = a_lo + x·a_hi`, `b' = b_lo + x^{-1}b_hi`, `G' = G_lo + x^{-1}G_hi`. Two chains
share the transcript: `b = eq(u)` unshifted, `b = shifted-eq(u)` shifted, with the verifier using the
closed forms of contract (c). Verification is `O(n)` (the folded-generator MSM).

**Audit.** The fold algebra is exactly Bulletproofs; I re-derived
`⟨a',G'⟩ = ⟨a,G⟩ + x^{-1}⟨a_lo,G_hi⟩ + x⟨a_hi,G_lo⟩` and the matching `⟨a',b'⟩` relation. `U` is drawn
after the commitments (challenge `IPA:x_u`), as Bulletproofs requires.

**Non-hiding.** No blinders. This is the plain, non-zk IPA, not the zero-knowledge variant.

**Efficiency (fixed in this review — §4.2).** The round loop was computing `L` and `R` with
`generators[2t] * a[2t+1]` — **2N individual elliptic-curve scalar multiplications**, not an MSM. bb's
production `ipa.hpp` uses `pippenger_unsafe` on contiguous lo/hi halves for precisely this reason. I
converted the backend from even/odd to lo/hi splits (which makes both MSM operands contiguous), routed
`L`/`R` through Pippenger, and batched the generator fold. See §4.2 for the measurement.

**Future work.** Adopt production `ipa.hpp`'s two-round-fused `batch_two_round_fold` and short 127-bit
challenge schedule (a further large constant factor); zk blinders.

### 2.8 Hyrax (ePrint 2017/1132)

**Math.** Pedersen row commitments over an `R × C` split; the prover sends the combined rows in the
clear and the verifier checks `Commit(w_x)` against the homomorphically folded row commitments, then
evaluates the rank-2 claim equation. Deterministic — no queries, no extra challenges.

**Deviation from the paper.** This is Hyrax's *uncompressed* variant. The paper composes a Bulletproofs
inner-product argument on top of the `√N` combined row to reach an `O(log N)` proof; that step is not
implemented, which is why the measured proof is 566 KiB rather than a few KiB. It is also non-hiding
(the paper's Pedersen commitments carry blinders). Both are honest simplifications, but they mean the
row labelled "Hyrax" in the comparison table is not the scheme's headline design point.

**Future work.** Compose with the (now fast) IPA backend to recover the log-size proof; hiding blinders.

### 2.9 WHIR (ePrint 2024/1586)

**Math.** Per iteration: `k` sumcheck rounds against a weight polynomial (array and weight table fold
at each challenge); commit the folded polynomial on the *halved* domain; one out-of-domain sample
`z_ood` with `y_ood = ĝ(z_ood)`; `t` in-domain queries into the previous oracle, each opening a
`2^k`-coset whose fold value becomes a claim about the new oracle at `x^{2^k}`; γ-batch the OOD and
query claims into the next iteration's weight. The final polynomial is sent in the clear and checked
against the queried cosets plus the accumulated weighted sum.

The rate improves by `k−1` bits per iteration (`log_domain` drops 1, `num_variables` drops `k`), which
is what makes WHIR's query counts fall off round by round.

**Honk fit.** Batched by ρ into one virtual round-0 oracle; shifted columns enter with a `x^{-1}` factor
per coset position (contract (a) in the evaluation domain). Columns of one Oink round share one Merkle
tree with **column-major interleaved leaves**, so one authentication path opens every column at a query
index — this is why WHIR's proof is 1.2 MiB rather than 3.4.

**Audit.** `fold_coset` is correct: `even = (v₊+v₋)/2`, `odd = (v₊−v₋)/(2x)`, result `(1−α)even + α·odd`
— matching bb's Lagrange fold rather than the classical `A_e + αA_o`, consistently with §1.3. The
Merkle leaf grouping (`leaf j ← {j + t·num_leaves}`) is exactly the coset `{ω^j·η^t}` that `fold_coset`
assumes. The schedule invariant `rounds[i].log_domain = rounds[i−1].log_domain − 1` holds, so querying
`folded_trees[i−1]` at round `i` targets the right oracle.

**Soundness regime.** `whir_config.hpp` implements four regimes. The default `CONJECTURED_LIST`
(`t = ⌈λ/r⌉`) rests on the up-to-capacity mutual-correlated-agreement conjecture that Crites–Stewart
(ePrint 2025/2046) **disprove**. The module ships a `REPAIRED_LIST` preset solving `H_q(δ*) = 1−ρ` by
deterministic integer fixed-point iteration with every rounding pushed toward more queries. At BN254 the
repair costs exactly one extra query per aggressive-rate round and **+2.0% proof size** at every size
from 2^12 to 2^20. Good work; the same repair is missing from Ligero/SwitchFold (§3.2).

**zk.** Implemented, and the only backend besides production Shplemini that has it: committed arrays
gain one variable (blinding coefficients in the high half at offset `2^m` or `2^m+1` for to-be-shifted
columns, so the shift contract's zero slot survives), leaves are salted, and the opening batches in a
fresh uniformly random mask polynomial whose claimed evaluation is revealed. `num_blinding_coefficients
= round0_queries + 8` — round-0 queries are the only openings of per-polynomial leaves, so that many
must remain information-theoretically blinded, with margin.

**Future work.** Grinding (proof-of-work) to cut query counts; stdlib recursive verifier
(`RECURSION.md` estimates ~0.5–1.2M gates with the shared-tree layout); UltraZK+WHIR wiring.

### 2.10 Ligero (AHIV17 / ePrint 2022/1608, Brakedown's tensor query)

**Math.** Reshape to `R × C`, RS-encode every row to `C·2^r`, and interleave into one Merkle tree with
leaf `j` = column `j` of all rows. The evaluation is the tensor sandwich `P(u) = bᵀMa`. The opening
sends three combined rows (`w_u`, `w_s`, `w_s2` — the rank-2 shift decomposition), then `t` sampled
columns, each checked as `Enc(w_x)[j] = B_xᵀ·col_j`.

**Honk fit.** Contract (b). One commit-and-query round, no folding at all — the prover-optimal corner.

**Audit.** The claim equation and the three encodings match. The rows travel unhashed with a digest
absorbed for Fiat–Shamir, so the column challenges are properly bound.

**Two findings**, §3.2 (the RS query count still uses the disproved-conjecture rule) and §3.3 (batching
doubles as the proximity test — sound only because `u` is a challenge).

**zk.** Not implemented; the README sketches the design (one random masking row per group, salted
leaves, random column padding).

### 2.11 SwitchFold (ePrint 2026/1489)

**Math.** Keeps Ligero's commitment and replaces the clear-text combined rows — the `O(√N)` term that
dominates Ligero's proof — with a recursive *code-switching* descent. The three combined rows become
segments of one message `M` of length `4·num_cols` (segment index in the top two variables, so every
claim about `M` is a tensor). Each level commits `M` reshaped under the next shorter RS code, opens it at
that level's queries, collapses all pending claims to a single tensor claim with one inner-product
sumcheck, and splits that claim into the column combination that becomes the next message. This level's
queries become the next level's code-switching claims. A 16-element base message is sent in the clear.

**Instantiating over Reed–Solomon collapses two of the paper's three modules**: the RS generator-matrix
row at a domain point *is* the pow tensor, so the code-switching claim needs neither a commitment to the
generator matrix nor the paper's accumulation scheme. Because nothing of row length is transmitted, the
proof-optimal shape moves: `num_cols` goes as large as BN254's 2-adicity allows.

**Audit.** I verified the tensor bookkeeping end to end. `append_shifted_eq_terms` decomposes
`shr(eq(u))` into `log C` eq terms by carry position — term `k` pins bits `<k` to 0 and bit `k` to 1, so
it fires exactly when `k` is the lowest set bit of the index, giving `Π_{j<k}u_j·(1−u_k)·Π_{j>k}eq`.
Correct. The claim split (low `log_sub` challenges bind the within-column index, high ones give the
column combination) matches the `message[c·sub_length + i]` layout under LSB-first ordering. The derived
code-switching values `Σ_c comb[c]·leaf[c] = Enc(next)[index]` by linearity. Correct.

**Result.** Reproduces the paper's headline 3.5x proof-size reduction over Ligero at 2^18 and does
better at 2^20 (5.9x), because Ligero's rows grow as `√N` while the descent only adds levels. It is
also 1.6x smaller than WHIR at comparable verify time — the best hash-based proof size in the suite.

**Future work.** §3.2 (query rule); §4.5 (a 4x reduction in the dominant prover cost, quantified but
not implemented).

### 2.12 BrakedownHonk (ePrint 2021/1043)

`brakedown_code.hpp` implements the Spielman-style recursive linear-time code; `brakedown_honk.hpp` runs
Ligero's tensor PCS over it, making `LigeroHonk` vs `BrakedownHonk` a controlled A/B — same trace, same
tensor protocol, same hasher, only the row code and query rule differ.

**The result is the honest baseline the linear-time-code line has to beat, and it loses on all three
axes**: at 2^14, 1.9x slower prover, 21.6x slower verifier, 24.8x larger proof. The encoder *is* faster
(1.13–1.22x over the RS FFT, emitting 2.3x fewer symbols), but Brakedown's relative distance is 0.07
against RS's 0.75, so the provable interleaved test needs ~2936 queries against 50 — and opening 2936
Merkle paths costs more than the ~20% saved on encoding.

**This correctly relocates the target for Lightning/Bolt: not encoding speed, but distance.**

**Finding (§3.4).** The base-case code is claimed MDS; it is not.

### 2.12a Bolt (ePrint 2026/310) — landed during this review

`bolt/` implements Bolt's sketched code `C_H(x) = (x, C(Hx))` — a sparse random-LDPC parity-check `H`
applied to compress the message, with the expensive base code (RS here) applied only to the short
sketch — and `BoltHonk` runs Ligero's tensor PCS over it.

The interesting contribution is not the encoder but the **piecewise distance guarantee** (paper Claim
3.1): two distinct codewords differ either in more than `γ` of the systematic stretch *or* in more than
`δ` of the sketch stretch, which is strictly more than the diluted overall distance
`min(γ, δα/ρ)/(1+α/ρ)`. So the two stretches are queried *independently*, each at the count its own
distance warrants (428 and 50 at λ = 100), and a cheating prover must survive both. This is what the
new `QuerySegment` hook exists for.

The module's own headline finding is one this review endorses and would have predicted: **the LDPC
distance bound collapses at BN254.** `ω`'s leading term is `x·ln(q−1)`, so a fixed column degree
certifies less distance as the field grows — `j = 16, k = 128` gives `γ = 0.0941` at `q = 2^32` and
`γ = 0.00006` at BN254; recovering `γ ≈ 0.15` needs `j = 64, k = 256`. That is the same structural
observation as §2.13's small-field conclusion, arriving from the coding-theory side: **the schemes that
look best in the literature are tuned to small fields, and the tuning does not transfer.** Any product
claiming "pick your backend" over BN254 has to re-derive parameters per backend rather than adopt the
paper's.

### 2.13 Support: small_subgroup_ipa, triple_ipa, small_field

`small_subgroup_ipa` proves `⟨F, G⟩ = s` for small vectors over a multiplicative subgroup, with a
grand-sum polynomial `A`, a quotient `Q`, and a linear-time verifier. It carries the zk-sumcheck (Libra)
and the ECCVM translation-masking consistency. Its README contains an unusually good piece of analysis:
**a fifth opening `A(1) = 0` is required**, because at `X = 1` the boundary condition and the `j = 0`
recurrence collapse into a single equation with a one-dimensional kernel, and the explicit homogeneous
perturbation `δ_A(1) = δ, δ_A(g^j) = −δ/(g−1), δ_s = −δ/(g−1)` forges `s' = s + δ_s`. The README also
states plainly that the local `check_consistency` is *not* sufficient on its own — soundness needs the
Shplemini batched opening. That is exactly the kind of caveat a pluggable-backend product must carry
into its API contract, and it is the model the other modules should follow.

`triple_ipa` reduces ECCVM's three structured claims (eq, shifted-eq, pow) to one ordinary IPA claim on
Grumpkin. Note the documented non-property: it does **not** assert that `F'` is the shift of `F`; it only
proves the inner product against the shifted-eq tensor. That is fine for ECCVM's use but is a trap for
reuse.

`small_field` answers "can bb move UltraHonk to M31?" with a working M31/CM31/QM31 implementation and
measurements. The honest conclusion is worth quoting into the product question: base-field work is
8–10x cheaper, but **the challenge field gives most of it back** (a QM31 multiplication costs about one
Fr multiplication, and every sumcheck round after the first folds with QM31 challenges), so the realistic
sumcheck gain is ~3.7x, not 10x — and the switch is a second proving stack, not a configuration flag.

---

## 3. Correctness findings

### 3.1 Four backends ship deliberately insecure setups

| Backend | What is substituted | Consequence |
|---|---|---|
| KZH2 | `τ_i` sampled in-process from a Blake3 seed | Prover knows the trapdoors; binding gone |
| KZH3 | `μ1_i`, `μ2_i` likewise | Same |
| CHOPIN | `τ_Y = 2`, `[τ_Y]₂ = 2·[1]₂` published | Y-direction binding gone |
| Dory | `Γ2_i = s_i·G2` with public `s_i` | AFGHO binding gone (GT collisions computable) |

All four are labelled `@warning Test-only` in the source and the honest-prover cost profile is
unaffected, which is the point. **This is not a defect — it is a scope boundary that the review needs
to state loudly**, because these four rows in the comparison table are performance models. Dory's fix
is mechanical (cofactor-cleared hash-to-G2, same runtime); KZH2/KZH3/CHOPIN need real ceremonies, and
CHOPIN's is a *bivariate* ceremony that does not exist today.

### 3.2 The RS query rule still rests on a disproved conjecture

`RSCodePolicy::num_queries` returns `t = ⌈λ/r⌉` — the RS up-to-capacity rule. WHIR ships a
`REPAIRED_LIST` preset for exactly this bound because Crites–Stewart (ePrint 2025/2046) disproved the
mutual-correlated-agreement conjecture it rests on, and the same argument applies verbatim to the
interleaved-RS proximity test Ligero and SwitchFold use. At BN254 the repair costs about one extra
query (50 → 51) and ~2% proof size — cheap to fix, expensive to leave undocumented.

**Half of this finding was closed while the review was being written.** The shared,
distance-parameterized query hook that `PCS_CANDIDATES.md` §"Remaining candidates" item 3 asked for now
exists: `ligero/QuerySegment` plus `Code::query_plan(config)`, with `BrakedownCodePolicy` using the
*provable* rule `t = ⌈λ / −log₂(1 − δ/3)⌉` and `BoltCodePolicy` returning two segments so each stretch
of its piecewise-distance code is queried at the count its own distance warrants. That is a better
design than the one this review would have recommended.

What remains is narrower and still worth doing:

1. **`RSCodePolicy` has no repaired-capacity variant.** The hook is there; the corrected rule is not.
   WHIR's `compute_num_queries(REPAIRED_LIST)` already implements it in deterministic integer
   arithmetic and could be called directly.
2. **SwitchFold does not use the hook at all** — it still reads `config.num_queries` in four places
   (`switchfold.hpp:364, 484, 604, 665`), so it cannot be run over a non-RS code and will not pick up
   the repaired rule when it lands.

### 3.3 The Ligero-family proximity test rides on `u` being a Fiat–Shamir challenge

Ligero, SwitchFold and Hyrax use a *single* combination for both jobs: the row weights are
`ρ^i·eq_r(u_hi)`, and the same combined rows serve as the proximity test and the evaluation extraction.
Classical Ligero/Brakedown run a *separate* proximity test with a uniformly random vector, precisely
because the eq tensor is not uniform.

This is sound here — tensor-structured randomness supports proximity gaps (Diamond–Posen, ePrint
2023/630) — **but only when `u` is uniformly random and drawn after the commitments.** Inside
`TransparentHonk` it is: `u` is the sumcheck challenge. But `LigeroProver::prove(ck, claims, u, transcript)`
accepts an arbitrary caller-supplied point, and nothing in the signature, the doc-comment, or the README
says the point must be a challenge. A user of a pluggable-PCS product who opens at a *public, fixed*
point — a perfectly reasonable thing to want from a PCS — gets no proximity guarantee at all.

Two actions: document the precondition in the API, and (for honesty in the query count) account for the
Diamond–Posen soundness loss, which is a `log N`-ish factor, in `num_queries`.

### 3.4 Brakedown's base-case code is not MDS

`build_base` constructs the parity block as `V_{ji} = p_j^i` with `p_j = j+2` — a Vandermonde in the
*coefficient* basis. The README and the code comment both claim this is MDS with relative distance
`(r−1)/r`, which the recursion's distance argument relies on. It is not:

- `[I | V]` is MDS iff every square submatrix of `V` is nonsingular. Square submatrices of a Vandermonde
  are *generalized* Vandermonde matrices, which are nonsingular over ℝ but can be singular over `F_p`.
- The code's actual codewords are `(coeffs(f), f(p_0), …, f(p_{parity−1}))` — the systematic part holds
  the *coefficients*, not evaluations, so this is not a Reed–Solomon code at all. A nonzero `f` of
  degree `< k` has at most `k−1` roots, so the provable minimum distance is only
  `1 + max(0, parity − (k−1))`. With `base_length = 32` and `parity = 24` at `r = 1.72`, that bound is
  **1**, against the `d = 25` a true MDS `[56, 32]` code would give.

To be precise about severity: the code's *actual* distance is probably fine — a low-weight codeword
needs a sparse coefficient vector whose polynomial vanishes on almost all of the 24 tiny evaluation
points `{2, …, 25}`, which is very unlikely to exist. But Brakedown's entire value proposition is a
distance that holds with *provable* failure probability `2^-100`, and the recursion's distance argument
consumes the base case's distance as a hypothesis. An unproven base case makes the whole chain unproven.

**Fix:** use a Cauchy matrix, `C_{ij} = 1/(y_i − z_j)` with disjoint `{y}`, `{z}`. Cauchy matrices are
superregular (every square submatrix is nonsingular), so `[I | C]` is genuinely MDS and the systematic
generalized-Reed–Solomon distance `(r−1)/r` holds. One function, no cost change.

Also minor: `build_base` constructs a `SeededPrng` and then discards it (`static_cast<void>(prng)`) —
dead code left from an earlier random-points design.

### 3.5 Query-index collisions are not deduplicated

Every hash-based backend derives query indices as `challenge.data[0] & mask` independently per query.
Duplicates are possible, and duplicated queries do not add soundness. This is standard practice in
deployed FRI implementations and the i.i.d.-sampling bound `ρ^t` is the correct one for sampling *with*
replacement, so this is not a bug — but it is worth a comment, because a reader checking `t` distinct
queries against the analysis will not find them.

### 3.6 `index_from_challenge_bounded` has documented modulo bias

`ligero::detail::index_from_challenge_bounded` uses `% bound` for non-power-of-two codeword lengths
(Brakedown). The comment states the bias is below 2^-40 at these sizes, which is right, but the
constant is worth recomputing if the module is ever used with a much smaller security parameter.

### 3.7 Mercury/Vela/CHOPIN assert rather than reject on the shift precondition

`BB_ASSERT_EQ(column[0], fr::zero(), "to-be-shifted polynomial must have zero constant term")` is a
prover-side assertion, which is stripped in release builds. That is correct for these backends (an
honest prover always satisfies it, and a dishonest one only breaks its own proof), but if the PCS layer
is ever exposed to untrusted *prover-side* input in a product, it should be a real check.

---

## 4. Efficiency: what was wrong, what was fixed, what remains

### 4.1 The shared defect: the ρ-batched linear combination was serial everywhere

Every backend starts by forming the ρ-combination of the ~41 claims:

```cpp
for (each claim) {
    for (size_t k = 0; k < n; ++k) { array[k] += rho_power * column[k]; }
    rho_power *= rho;
}
```

That is a **serial `num_claims × N` pass** — 10.7M multiply-adds at 2^18. Production Gemini does the
same job with `bb::add_scaled_batch`, which splits the *output* range across threads and loops the
sources inside each chunk. The experimental backends could not reuse it directly because it requires
`Polynomial` operands and they hold `std::vector<fr>`, so the loop was open-coded — serially — in
Mercury, Vela, CHOPIN, KZH2, KZH3, Dory and the IPA backend.

The size of this is not incidental: at 2^18 it is a **large fraction of the prover gap between the
constant-size pairing backends and the KZG baseline**, and KZG only avoids it because production
Gemini already parallelizes the identical work.

**Fix.** `utils/batch_accumulate.hpp` adds the `std::vector`-backed counterpart:

```cpp
struct ScaledTerm { const fr* source; fr scalar; bool shifted; };
void accumulate_scaled(std::span<fr> out, std::span<const ScaledTerm> terms);
```

one parallel pass over the output range, terms looped inside each chunk (disjoint writes, no locking,
each output cache line touched once per chunk rather than once per term). Applied to Mercury, Vela
(three sites, including the α-batch), CHOPIN, KZH2, KZH3, Dory and IPA.

### 4.2 The IPA backend was not using MSMs at all

`pedersen_ipa.hpp`'s round loop computed the cross terms as

```cpp
for (size_t t = 0; t < half; ++t) {
    left  += generators[2 * t]     * a[2 * t + 1];
    right += generators[2 * t + 1] * a[2 * t];
}
```

— **2N individual elliptic-curve scalar multiplications** across the whole reduction, where an MSM is
called for. bb's production `ipa.hpp` uses `pippenger_unsafe` on contiguous lo/hi halves for exactly
this. The even/odd split was what blocked it: Pippenger needs contiguous point and scalar spans.

**Fix.** Converted the backend from even/odd to **lo/hi** splits, which is also the classical
Bulletproofs presentation. Both cross-term operands become contiguous, so each is one
`pippenger_unsafe` call; the generator fold `G'[t] = G_lo[t] + x^{-1}G_hi[t]` is now a parallel pass
followed by one `batch_normalize` back to affine.

Because round `j` now binds the *highest* remaining variable rather than the lowest, the verifier's
per-variable fold multipliers are the round-challenge inverses **in reverse order**; this is the new
`detail::fold_multipliers`, threaded into both the folded-generator scalars and the eq / shifted-eq
closed forms. All four `PedersenIpaTest` cases (including the two negative tests) and `IpaHonkTest`
pass unchanged.

### 4.3 Six more serial loops

- **SwitchFold's combined-row build** was the serial version of a loop Ligero already parallelizes —
  `3 × 36 × N` multiply-adds, ~28M at 2^18. Now split over the *column* range so each thread owns a
  disjoint slice of all three segments (no private accumulators, no merge).
- **SwitchFold's `inner_product_round`** was a serial sumcheck round; now `parallel_for_range` with
  private accumulators, matching WHIR's `sumcheck_round_univariate`.
- **Hyrax's combined-row build** — same loop, same fix.
- **Vela's `laurent_product`** is the `O(N log N)` shift-and-add expansion, `log N` serial passes over a
  `2N` table. Parallelized with an explicit double buffer: the pass reads `table[e + 2^k]`, so an
  in-place parallel version would race (the same trap already fixed once in WHIR's array fold). This is
  the largest single contributor to Vela's 1.70x → 1.18x.
- **Mercury's column combination** `h[row] = Σ_j eq_lo[j]·f[row·b + j]` — an `O(N)` serial pass; now
  parallel over rows.
- **Mercury's Laurent accumulator** `accumulate_ip_terms` is `O(b²) + O(rows²) = O(N)`. Each shift writes
  only `S[shift−1]`, so the shifts are **independent** — parallelized directly with no accumulator merge.
  (Work per shift falls as `len − shift`, so index-ordered chunks are unbalanced by under 2x, which is
  not worth a custom schedule.)

### 4.4 Measured effect

**Methodology, and its limits.** The review machine was shared with another build/edit session for part
of the measurement window, so absolute wall-clock is not comparable between the committed baseline
(`pcs_bench_results.json`, single run) and the post-fix runs. Every figure below is therefore a **ratio
to the `ultra_honk_kzg_prove` baseline measured in the same run**, which cancels machine-condition
drift. **`LigeroHonk` is the control**: its code is untouched by this review, so the movement in its
ratio is the noise floor. The post-fix figures come from a 5-repetition median taken while the machine
was quiet; a later 3-repetition run under contention moved the control by 26% and was discarded rather
than reported. Proof bytes and verifier work are unchanged by every fix, so these are pure prover-side
wins — the proof-size column of the comparison table is unaffected.

All at 2^18, prove, ratio to same-run KZG:

| Backend | before | after | change |
|---|---|---|---|
| **Ligero** (control, code untouched) | 0.70x | 0.72x | +3% ← **noise floor** |
| CHOPIN | 1.05x | **0.81x** | −23% |
| Mercury | 1.27x | **1.14x** | −10% |
| Vela | 1.70x | **1.18x** | −30% |
| SwitchFold | 2.97x | **2.50x** | −16% |
| IPA (Pedersen) | 82.0x | **≈5.6x** | **13.3x faster** |
| Hyrax, KZH2, KZH3, Dory | — | — | unchanged within noise |

The IPA figure is 44 726 ms → 3351 ms in comparable single-run conditions; the change is large enough
that it survives the contention (CPU time in the noisy run was 3391 ms). The four unchanged backends
are exactly the ones whose provers are dominated by `N` elliptic-curve scalar multiplications per
column (and `R` Miller loops per column for Dory), where the ρ-batch was never more than a few percent
— that they did *not* move is a consistency check on the attribution, not a disappointment.

**Two results worth calling out.** CHOPIN now proves **faster than the KZG baseline** while producing a
9.2 KiB proof against KZG's 13.3, and Vela's prover penalty for the smallest-proof-in-the-suite drops
from 1.70x to 1.18x. Both were previously read as "pay prover time for proof size"; that trade is now
much weaker, which changes their standing in §7.3.

### 4.5 Quantified but not implemented

**SwitchFold's weight tables do 4x more work than necessary.** `WeightTerm::accumulate_table` builds the
full `2^{log C + 2}` tensor for every term, but `select_segment` appends two eq factors at the *constants*
0 and 1, which zero three quarters of the tensor — and those zeros are still computed (as multiplications
by zero) and added. With ~150 of the ~169 level-0 terms being per-query, segment-restricted `pow_weight`
terms, restricting the build to the live quarter is a straight 4x on the dominant cost, worth ~195 ms of
the remaining prover at 2^18. The general form is: partition each term's variables into *free* and
*pinned* (exactly one of `a`, `a+b` zero), build the tensor over the free variables only, and scatter into
the strided support. That also collapses the `log C` shifted-eq terms from `log C` full tables to ~1.

**Dory's verifier spends its time in GT exponentiation.** The ρ-combination `Π_p T_p^{ρ^p}` is 41
square-and-multiply exponentiations of 254 bits in `fq12` — ~15.6k `fq12` operations. Cyclotomic squaring
(a standard ~2x for GT) plus a bucket-method multi-exponentiation over GT (another ~4x at these counts)
would take the 362 ms verifier down substantially. Neither primitive exists in bb's `fq12` today.

**Ligero's verifier re-encodes three full rows.** Three size-`C·2^r` FFTs per verification. Direct
evaluation at only the `t` queried positions is *not* a win at the current shape (≈4.9M vs ≈3.3M
multiplications), so this is correctly left alone — recorded so it is not re-investigated.

**Mercury's binomial-fold recurrence is parallelizable and was left alone.** `q[k−b] = f[k] + α·q[k]`
looks sequential but decomposes into `b` independent chains, one per residue class mod `b` — so it
parallelizes over `b = √N ≈ 512` chains at 2^18. Together with the two remaining `O(N)` passes (the
`numerator` construction in round 4, and the `divide_by_linear` calls, which are genuinely sequential
recurrences of length `N`) this is the last few percent of Mercury's prover. Recorded for completeness;
the measured payoff is small relative to the risk of touching the fold.

**IPA can go substantially further.** The current fix uses one MSM per cross term per round. Production
`ipa.hpp` additionally fuses round pairs (`batch_two_round_fold`, deferring the SRS fold so the second
round's MSMs run against the pre-fold SRS) and uses a rescaled fold with 127-bit short challenges so
each `batch_mul` multiplies by a half-width scalar. Adopting both would be another large constant
factor on top of the 13.3x.

**KZH2/KZH3/Hyrax/Dory commit costs are inherent.** Their provers are dominated by `N` elliptic-curve
scalar multiplications per committed column (and `R` Miller loops per column for Dory). No
implementation change moves these; they are the schemes' cost.

---

## 5. Zero-knowledge: the decisive gap

| Backend | zk status |
|---|---|
| Gemini+Shplonk+KZG | **Implemented.** Sparse `2d`-coefficient masking poly + `small_subgroup_ipa` for Libra; proof sketch in `SHPLEMINI_ZK_MASKING.md` with a Lean development for the linear algebra |
| WHIR | **Implemented.** Blinded high half, salted Merkle leaves, batched random mask polynomial, `q = round0_queries + 8` blinding coefficients |
| Mercury, Vela, CHOPIN, KZH2, KZH3, Dory, Hyrax, IPA, Ligero, SwitchFold, Brakedown, Bolt | **None.** Non-hiding commitments, unmasked openings |

Thirteen of fifteen backends leak the witness. For a product this is not a footnote — a user who selects
"Mercury" and gets a non-zk proof has a security failure, not a performance trade.

The work is also **not uniform**, which matters for scoping:

- *Pairing/KZG family* (Mercury, Vela, CHOPIN): bounded-degree blinders on the small polynomials plus
  degree slack in the SRS. Each paper gives the design; each needs its own argument because each opens a
  different set of auxiliary polynomials.
- *`√N`/`∛N` families* (KZH2, KZH3, Hyrax, Dory): the combined rows are sent **in the clear**. Hiding
  them requires either hiding Pedersen commitments plus a zk inner-product argument (Hyrax's actual
  design), or masking rows — a structural change, not a wrapper.
- *Hash family* (Ligero, SwitchFold, Brakedown, Bolt): WHIR's recipe (mask polynomial + salted leaves +
  blinded region) transfers in shape, but Ligero's combined rows and SwitchFold's descent messages are
  transmitted in the clear and each needs its own masking-row design and query-count accounting.

**A pluggable-backend product cannot ship a backend without zk, and it cannot make zk a shared layer.**
This is the single largest item in §7.3.

One structural note that makes this less bleak than it looks: WHIR's zk design is *entirely inside the
PCS*, requiring no change to `TransparentHonk` or to anything above it — extra blinding coefficients in
a region the claim does not see, salted leaves, and one batched mask polynomial. That is evidence that
"zk-capable PCS" can be a property of the backend rather than a property of the proof system, which is
exactly what a pluggable design needs. It just has to actually be implemented thirteen more times.

---

## 6. Recursion

Only the production KZG/IPA path has stdlib (in-circuit) verifiers today. `whir/RECURSION.md` models
WHIR's: ~73 gates per Poseidon2 permutation in Ultra, ~5.4M gates for a naive layout, ~0.5–1.2M with the
shared-per-round trees, deduplicated openings and minimal Fiat–Shamir absorption that the native
implementation already has.

For the other backends the picture divides cleanly:

- **Pairing backends** (Mercury, Vela, CHOPIN, KZH2, KZH3) recurse like KZG: the verifier is an MSM plus
  a pairing, and bb already has the `biggroup`/`bigfield` machinery. Vela and Mercury are the cheapest
  (constant-size proof, `O(log N)` field work). KZH2/KZH3 need `√N`/`∛N` in-circuit MSM scalars, which is
  the point of the KZH line for accumulation.
- **Dory** needs in-circuit `fq12` arithmetic and `log R` GT operations — expensive and not currently
  available.
- **Hash backends** need in-circuit Poseidon2 Merkle verification (the `Poseidon2MerkleHasher` variant
  exists precisely for this) plus the sumcheck and weight bookkeeping. WHIR and SwitchFold are the
  plausible ones; Ligero's `t·P·R` opened values make the in-circuit cost scale with the proof, which is
  MiB-scale.
- **IPA/Hyrax** verify in `O(N)` group operations, which does not recurse at reasonable cost without the
  accumulation/deferral machinery `ipa.hpp` already implements for ECCVM (`NativeAccumulator`,
  `batch_verify_accumulators`).

---

## 7. Can a Noir frontend expose a selectable PCS backend?

### 7.1 What the evidence says yes to

The **architecture question is answered, and answered positively**. `TransparentHonk<Pcs>` is a real
demonstration that:

- the boundary between the proof system and the commitment scheme is narrow (§1.1) and does not leak;
- fifteen schemes spanning four trust models, three commitment types (Merkle digest, `vector<G1>`,
  `vector<GT>`) and three different shift constructions all fit behind it, unmodified above the PCS;
- the Fiat–Shamir schedule, relation parameters, derived polynomials and sumcheck are genuinely shared
  with the production prover, so backends do not fork the proof system;
- the trade-offs a user would actually select on are real and large — at 2^18 the suite spans 7.8 KiB to
  2.7 MiB in proof size, 3.5 ms to 3.6 s in verification, and 0.72x to ~46x the KZG prover;
- the abstraction keeps absorbing new structure rather than fighting it. Three things landed *during*
  this review — a linear-time code, a code-policy parameterization, and a piecewise-distance query plan
  for a sketched-LDPC backend — and none required a change above the PCS layer. An abstraction that
  survives contact with genuinely different cryptography is the thing you cannot fake.

That is a strong result and it is the hard part of the product. A `#[backend(vela)]` attribute on a Noir
program is not a fantasy.

### 7.2 What blocks it today

Ranked by how much they cost to fix, most binding first.

1. **Zero-knowledge on 13 of 15 backends** (§5). Not a shared layer; ~3 distinct designs plus per-scheme
   arguments. This is the gating item.
2. **Four backends have no real setup** (§3.1). Dory needs hash-to-G2 (mechanical). KZH2/KZH3 need a
   ceremony. CHOPIN needs a *bivariate* ceremony that does not exist.
3. **Recursive verifiers exist for one backend family** (§6). A product whose proofs cannot be verified
   inside another Aztec circuit is a much narrower product.
4. **The RS query rule still rests on a disproved conjecture** (§3.2). Now small work — the hook exists —
   but it is a correctness-of-parameters issue, which is exactly the class of bug a user selecting a
   backend cannot audit for themselves.
5. **API preconditions are undocumented** (§3.3, §3.7). A pluggable PCS is used by people who did not
   read the module; "the opening point must be a Fiat–Shamir challenge" has to be in the type system or
   the doc-comment, not only in the caller.
6. **Parameters do not transfer from the papers.** Bolt's LDPC distance collapses from `γ = 0.094` to
   `γ = 0.00006` moving from `2^32` to BN254 (§2.12a); WHIR's repaired-conjecture penalty scales as
   `1/log₂ q`; the M31 study finds the challenge field eats most of the small-field win (§2.13). A
   backend menu has to own per-backend parameter derivation at *its* field, not cite the paper's table.
7. **Proof-size-driven backends need proof compression.** The hash family sits at 0.77–2.7 MiB. ePrint
   2025/1446 (`zip`) is tracked in `PCS_CANDIDATES.md` and attacks exactly this.

### 7.3 A roadmap that would make it shippable

**Phase 1 — make the honest set small and real.** Pick three backends to productize rather than
fifteen: **Gemini+Shplonk+KZG** (has zk, has recursion, is the default), **Vela** (smallest proof,
fastest verifier, and — importantly — a real ceremony already available, since it reuses bb's existing
KZG commitment key), and **WHIR** (has zk, transparent, recursion designed). Ship those three behind the
selector and gate the remaining twelve as research backends. This converts the zk gap from "13 backends"
to "1 backend" (Vela), which is a scoped piece of work rather than a programme.

Note that Vela's case strengthened during this review: its prover penalty was 1.70x KZG and is now 1.18x
(§4.4), so "smallest proof and fastest verifier" no longer costs much on the third axis.

**Phase 2 — factor the cross-cutting concerns out of the backends.**
- Finish the query rule: `RSCodePolicy` should return the repaired-capacity count (WHIR already computes
  it), and SwitchFold should consume `query_plan` instead of `config.num_queries`. The hook is done.
- A masking interface in the `Pcs` concept, so `TransparentHonk` can *require* zk rather than each
  backend optionally providing it. WHIR's design (§5) shows this can live entirely inside the backend.
- A precondition contract on the opening point, enforced by construction — e.g. the point can only be
  obtained from the transcript, so §3.3's unsound usage is unrepresentable.
- A per-backend parameter-derivation story at BN254 rather than the paper's field (§7.2 item 6).

**Phase 3 — widen.** Add backends only when they clear the Phase-2 bar. On current evidence the best
next candidates are **CHOPIN** (now the fastest prover in the pairing corner, *if* a bivariate ceremony
is solvable), **SwitchFold** (best hash-based proof size, verifier already competitive) and **Mercury**
(constant proof, well-understood zk design), followed by **Titan** from `PCS_CANDIDATES.md`'s list.

**Phase 4 — the frontend.** Only at this point does the Noir-side work (backend attribute, per-backend
proof types, verifier codegen) become the bottleneck rather than the cryptography.

### 7.4 The one-sentence answer

**Yes, the product is possible — the pluggable architecture is already built and validated across
fifteen schemes — but what exists today is a research bench, not a backend menu: the gap is not
architecture, it is that thirteen of the fifteen backends have no zero-knowledge, four have no real
setup, and fourteen have no recursive verifier.**

---

## Appendix A: files changed by this review

| File | Change |
|---|---|
| `utils/batch_accumulate.hpp` | New: parallel `std::vector`-backed RLC accumulator |
| `mercury/mercury.hpp` | Parallel ρ-batch (both chains); parallel column-combination `h`; parallel Laurent accumulator (shifts are independent) |
| `mercury/vela.hpp` | Parallel ρ-batch, parallel α-batch, parallel double-buffered `laurent_product` |
| `mercury/chopin.hpp` | Parallel ρ-batch (both chains) |
| `kzh/kzh.hpp`, `kzh/kzh3.hpp` | Parallel ρ-batch |
| `dory/dory.hpp` | Parallel ρ-batch (array side; the row-commitment side is EC work) |
| `hyrax/hyrax.hpp` | Column-parallel combined-row build |
| `switchfold/switchfold.hpp` | Column-parallel combined-row build; parallel sumcheck round |
| `pedersen_ipa/pedersen_ipa.hpp` | lo/hi split; Pippenger cross terms; batched generator fold; reversed fold multipliers |

## Appendix B: papers audited against

Mercury [2025/385], Vela/Carina [2026/1438], CHOPIN [2026/480], KZH/KZH-k [2025/144],
Dory [2020/1274], Bulletproofs/IPA, Hyrax [2017/1132], WHIR [2024/1586],
Crites–Stewart proximity-gap disproof [2025/2046], Ligero [2022/1608], Brakedown [2021/1043],
SwitchFold [2026/1489], Bolt [2026/310], BDFG multi-point batching [2020/081],
Diamond–Posen tensor proximity [2023/630].
