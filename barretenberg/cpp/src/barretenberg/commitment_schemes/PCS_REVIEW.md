# Polynomial commitment schemes in barretenberg: correctness, efficiency, and the pluggable-backend question

A review of every polynomial commitment scheme in `commitment_schemes/`, written to answer one
question: **can a Noir-like frontend expose the proof system's commitment scheme as a user-selectable
backend?** The answer is developed in §7; §§1–6 are the evidence.

> A colour-coded companion to this document, with the central equations annotated term by term, is
> `PCS_REVIEW.html` — self-contained, open it in a browser.

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
  fixed to a public constant (CHOPIN's $\tau_Y = 2$), or known discrete logs of the $\mathbb{G}_2$
  generator (Dory). Each is clearly labelled in the source; none is a mistake. But it means their
  measured numbers are performance models, not security claims.
- **The Reed–Solomon query count rests on a conjecture that has since been disproved** (§3.2), the same
  one WHIR carries a repaired preset for. The shared query-plan hook this needed landed mid-review; the
  repaired RS rule behind it has not, and SwitchFold is not wired to the hook.
- **Ligero-family batching doubles as the proximity test, which is only sound because $\vec{u}$ is a
  Fiat–Shamir challenge** (§3.3). Correct inside `TransparentHonk`; unsound if the standalone API is
  called at an attacker-chosen point, which nothing in the API documents or enforces.
- **Brakedown's base-case code is not MDS** (§3.4), contradicting the stated property its distance
  argument relies on. A Cauchy matrix fixes it in one function.

**Efficiency.** A single shared defect dominated the whole experimental suite: **the $\rho$-batched
random linear combination of the $\approx 41$ Honk claims was written as a serial nested loop in every
non-production backend**, while production Gemini does it with the parallel `add_scaled_batch`. Fixing
that, plus six related serial loops and one algorithmically wrong IPA inner loop, gives the measured
improvements in §4.4. Headline: **the IPA backend's prover went from $82\times$ KZG to
$\approx 2.2\times$** (rebuilt on the production `ipa.hpp` core — §2.7, §4.2 — which also cut its
verify from 3.6 s to ~0.1 s and its proof below KZG's), **CHOPIN now proves faster than the KZG
baseline** ($0.81\times$, previously $1.05\times$) while producing a smaller proof, **Vela's prover
penalty drops from $1.70\times$ to $1.18\times$**, and **SwitchFold improves 16%**. Apart from the
IPA rebuild, every fix leaves proof bytes and verifier work untouched.

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
boundary between "proof system" and "commitment scheme" is real and narrow.

### 1.2 The claim set a backend must open

After sumcheck at challenge $\vec{u} \in \mathbb{F}^n$, UltraHonk hands the PCS **41 claims over 36
committed columns**:

$$
\underbrace{P_i(\vec{u}) = v_i}_{i \,=\, 1,\ldots,36}
\qquad\text{and}\qquad
\underbrace{\mathsf{shift}(P_j)(\vec{u}) = w_j}_{j \,=\, 1,\ldots,5},
\qquad
\mathsf{shift}(P)[i] := P[i+1].
$$

Columns committed in the same Oink round share one group — 5 groups: 28 precomputed, 3 wires,
3 counts $+\, w_4$, 1 `lookup_inverses`, 1 `z_perm`.

### 1.3 bb's array-is-both-bases convention

A polynomial is one dense array $a$ of length $N = 2^n$, read simultaneously as

- the **evaluation table** of the multilinear extension $\tilde{a}$ over $\{0,1\}^n$, and
- the **coefficient vector** of a univariate $A(X) = \sum_i a_i X^i$ of degree $< N$.

The stride-2 Lagrange fold

$$
a'[j] \;=\; (1-\alpha)\,a[2j] \;+\; \alpha\,a[2j+1]
$$

is *simultaneously* partial multilinear evaluation at the lowest variable and the univariate even/odd
fold. Splitting $A(X) = A_e(X^2) + X A_o(X^2)$, the folded array's univariate is exactly

$$
A'(Y) \;=\; (1-\alpha) A_e(Y) \;+\; \alpha A_o(Y),
$$

so one array and one loop serve both readings, and no basis conversion ever happens. Gemini, WHIR's
folding, Mercury's twin, Vela's Laurent twin and SwitchFold's descent all exploit this; it is why a
15-backend suite can share one data layout and one `eq_tensor`.

### 1.4 Three ways to handle the Honk shift

Every backend has to open $\mathsf{shift}(P)$ given a commitment to $P$, and the suite contains three
genuinely different solutions. This is the part of the design space that a pluggable-backend product
must get right, because it is where a new backend most easily goes wrong.

#### (a) Univariate division contract — Gemini/Shplemini, WHIR, Mercury, Vela

The univariate twin of $\mathsf{shift}(P)$ is exactly $P(X)/X$, an *exact* division because the Honk
shift contract guarantees $P[0] = 0$. Mercury multiplies its fold identity through by $X$,

$$
B(X) \;=\; Q(X)\,(X^b - \alpha) \;+\; X\,g_B(X),
\qquad
Q(X) \,:=\, X\,q_B(X),
$$

so the verifier can use the *unshifted* homomorphic commitment $C_B$. Vela folds the shifted chain into
the same Laurent twin,

$$
g(X) \;=\; f_A(X) \;+\; \lambda\,X^{-1} f_B(X).
$$

Precondition: $P[0] = 0$, asserted at prove time.

#### (b) Rank-2 tensor decomposition — Ligero, SwitchFold, Hyrax, KZH2, Dory (rank 3 for KZH3)

The shifted-eq vector is not a tensor, but it is a rank-2 sum of tensors. Reshaping the flat index as
$j = rC + c$ over an $R \times C$ matrix, with

$$
a = \bigotimes_{i < c_{\log}} (1-u_i,\, u_i) \in \mathbb{F}^C,
\qquad
b = \bigotimes_{i \ge c_{\log}} (1-u_i,\, u_i) \in \mathbb{F}^R,
$$

the shifted weight matrix is

$$
Q'_{r,c} \;=\;
\begin{cases}
b_r\, a_{c-1}, & c \ge 1,\\[4pt]
b_{r-1}\, a_{C-1}, & c = 0,\ r \ge 1,\\[4pt]
0, & c = r = 0,
\end{cases}
\qquad\text{i.e.}\qquad
Q' \;=\; \underbrace{b \otimes \mathsf{shr}(a)}_{\text{main term}}
\;+\; \underbrace{\mathsf{shr}(b) \otimes \big(a_{C-1}\, e_0\big)}_{\text{carry term}} .
$$

The carry term costs one extra combined row ($w_{s2}$), of which only entry $0$ is ever used. **No
zero-constant-term precondition** is needed — the $(0,0)$ entry is dropped structurally. I verified the
rank-3 extension in `kzh3.hpp` term by term: the shift of the flattened index
$k = (i_1 d_2 + i_2) d_3 + i_3$ crosses both the $d_3$ and $d_2$ boundaries, giving

$$
\big\langle \mathsf{shift}(f),\, \mathsf{eq}(\vec{u}) \big\rangle
\;=\;
\underbrace{\sum_{i_3 \ge 1} T_3[i_3]\, a_3[i_3 - 1]}_{\text{within the } d_3 \text{ row}}
\;+\;
\underbrace{a_3[d_3-1]\, T_3'[0]}_{d_3 \text{ boundary}}
\;+\;
\underbrace{a_2[d_2-1]\, a_3[d_3-1]\, T_3''[0]}_{d_2 \text{ boundary}} .
$$

Correct.

#### (c) Successor-kernel closed form — Pedersen IPA, Dory's scalar side

The IPA fold reduces the vector $b$ to a single scalar
$\sum_i b_i \prod_{j \in \mathsf{bits}(i)} m_j$, with $m_j$ the multiplier of variable $j$. For
$b = \mathsf{eq}(\vec{u})$ that collapses to a product:

$$
\Big\langle \mathsf{eq}(\vec{u}),\ \bigotimes_j (1,\, m_j) \Big\rangle
\;=\; \prod_j \big( (1-u_j) + u_j m_j \big).
$$

For the shifted-eq vector ($b_0 = 0$, $b_i = \mathsf{eq}_{i-1}(\vec{u})$) it telescopes. Split on $k$,
the position of the **lowest set bit** of the index — so that $i-1$ clears bit $k$ and sets every bit
below it:

$$
\sum_k
\underbrace{\Big[\prod_{j<k} u_j\Big]}_{\text{bits below } k \text{ set}}
\underbrace{(1-u_k)}_{\text{bit } k \text{ cleared}}
\underbrace{m_k}_{\text{fold factor}}
\underbrace{\Big[\prod_{j>k} \big((1-u_j) + u_j m_j\big)\Big]}_{\text{bits above } k \text{ free}} .
$$

The verifier never touches a shifted commitment. I re-derived both closed forms independently; both are
correct.

---

## 2. The backends

Fifteen, grouped by trust model and the resource each optimizes. Prove times are given as a **ratio to
the KZG baseline** (see §4.4 for why); verify and proof size are absolute, at $N = 2^{18}$. IPA's
verify and proof reflect its rebuild on the production core (§2.7, §4.2); every other backend's are
unaffected by this review's changes.

| Backend | Family | Trust | Prove ($\times$KZG) | Verify | Proof | zk |
|---|---|---|---|---|---|---|
| Gemini+Shplonk+KZG | pairing | SRS ceremony | $1.00\times$ | 4.7 ms | 13.3 KiB | **yes** |
| Mercury | pairing, constant-size | SRS ceremony | $1.14\times$ | 4.2 ms | 9.4 KiB | no |
| Vela | pairing, reciprocal | SRS ceremony | $1.18\times$ | **3.5 ms** | **7.8 KiB** | no |
| CHOPIN | pairing, bivariate KZG | **test-only SRS** | $\mathbf{0.81\times}$ | 4.7 ms | 9.2 KiB | no |
| KZH2 | pairing, $\sqrt{N}$ | **test-only SRS** | $\sim\!14\times$ | 216 ms | 183.3 KiB | no |
| KZH3 | pairing, $\sqrt[3]{N}$ | **test-only SRS** | $\sim\!46\times$ | 79.7 ms | 55.3 KiB | no |
| Dory | pairing, two-tier | **test-only $\mathbb{G}_2$** | $\sim\!12\times$ | 362 ms | 94.3 KiB | no |
| IPA (Pedersen) | DL | transparent | $\sim\!2.2\times$ | 107 ms | 12.0 KiB | no |
| Hyrax | DL, Pedersen rows | transparent | $\sim\!11\times$ | 243 ms | 566.3 KiB | no |
| WHIR | hash, RS folding | transparent | $1.08\times$ | 13.5 ms | 1219.7 KiB | **yes** |
| SwitchFold | hash, code switching | transparent | $2.50\times$ | 13.6 ms | 772.9 KiB | no |
| Ligero | hash, tensor $\sqrt{N}$ | transparent | $\mathbf{0.72\times}$ | 20.6 ms | 2692.6 KiB | no |
| BrakedownHonk | hash, linear-time code | transparent | $1.9\times$ Ligero @ $2^{14}$ | $21.6\times$ Ligero | 20.2 MiB @ $2^{14}$ | no |
| Bolt | hash, sketched LDPC | transparent | see `bolt/README.md` | — | — | no |

Support modules, not backends: `small_subgroup_ipa` (the zk-sumcheck and ECCVM-translation inner
product), `triple_ipa` (ECCVM's three-claim reduction on Grumpkin), `small_field` (M31/QM31 study).

### 2.1 Gemini + Shplonk + KZG (Shplemini) — the production baseline

**Math.** Gemini reduces $k$ multilinear claims at $\vec{u}$ to $d+1$ *univariate* claims. Batch

$$
A_0 \;=\; F + G/X,
\qquad
F = \sum_i \rho^i f_i,
\qquad
G = \sum_j \rho^{\,k+j} g_j ,
$$

then fold and commit $A_1, \ldots, A_{d-1}$:

$$
A_{l+1}(X) \;=\; (1 - u_l)\,\mathsf{even}(A_l)(X) \;+\; u_l\,\mathsf{odd}(A_l)(X).
$$

Open $A_0^{+}$ at $r$, $A_0^{-}$ at $-r$, and $A_l$ at $-r^{2^l}$. The verifier reconstructs the
positive evaluations from the Gemini relation, so only the negative ones are sent. Shplonk collapses the
$d+1$ single-point claims to one point $z$ via

$$
Q(X) \;=\; \sum_j \nu^{\,j-1}\,\frac{f_j(X) - v_j}{X - x_j},
$$

and KZG closes with one pairing. The verifier's whole job is **one MSM** over
$[Q], [f_1], \ldots, [f_n], [1]$ plus a pairing.

**Honk fit.** Native — this *is* the Honk PCS. Shift via contract (a); the $1/X$ becomes the $r^{-1}$
factor in the shifted Shplonk scalar:

$$
s_{\text{unshifted}} = \frac{1}{z-r} + \frac{\nu}{z+r},
\qquad
s_{\text{shifted}} = \frac{1}{r}\left(\frac{1}{z-r} - \frac{\nu}{z+r}\right).
$$

**zk.** The only fully worked-out story in the tree: a sparse masking polynomial
$M(X) = \sum_{s \in S} c_s X^s$ with $|S| = 2d$ i.i.d. uniform coefficients on the tail-halving support

$$
S = \big\{ E-1,\, E-2,\ \ N/2,\, N/2 - 1,\ \ N/4,\, N/4 - 1,\ \ \ldots,\ 2,\, 1 \big\},
$$

plus `small_subgroup_ipa` for the Libra univariates. `SHPLEMINI_ZK_MASKING.md` gives the reduction and
points at a Lean development for the linear algebra. There is a real edge case handled explicitly: if
the Gemini challenge $r$ lands in the small subgroup used by `small_subgroup_ipa`, prover evaluations
leak, and `gemini_impl.hpp` aborts rather than proving.

**Audit.** Production code, internally audited (`ipa.hpp` carries an audit banner; `gemini_impl.hpp`
one planned). The prover is well optimized: SIMD `fold_stride2`, extent tracking that skips the
provably-zero tail of each fold, `parallel_for_heuristic` throughout.

**Future work.** None specific to this review.

### 2.2 Mercury (ePrint 2025/385)

**Math.** With $N = 2^n$, $t = \lceil n/2 \rceil$, $b = 2^t$, the array reshapes to
$\mathsf{rows} \times b$. The claim $\tilde{f}(\vec{u}) = v$ factors through the column-combined
polynomial

$$
h(X) \;=\; \sum_j \mathsf{eq}_j(\vec{u}_{\mathrm{lo}})\, f_j(X)
$$

into four sub-claims:

1. **Fold correctness.** $f(X) = q(X)(X^b - \alpha) + g(X)$ for a verifier challenge $\alpha$, checked
   at a random $\zeta$ through the quotient
   $H(X) = \big(f(X) - (\zeta^b - \alpha) q(X) - g(\zeta)\big) / (X - \zeta)$ and the pairing

$$
e\big(C_f - (\zeta^b - \alpha) C_q - g(\zeta)[1]_1,\ [1]_2\big) \;=\; e\big(C_H,\ [\tau - \zeta]_2\big).
$$

2. **Degree bound** $\deg g < b$, via $D(X) = X^{b-1} g(1/X)$ — committable only if the bound holds.
3. $\langle P_{\vec{u}_{\mathrm{lo}}},\, g \rangle = h(\alpha)$, and
4. $\langle P_{\vec{u}_{\mathrm{hi}}},\, h \rangle = v$.

The two inner products batch by $\gamma$ into one Laurent identity:

$$
\underbrace{\sum_k \gamma^k \big( p_k(X)\,r_k(1/X) + p_k(1/X)\,r_k(X) \big)}_{\text{symmetric Laurent product}}
\;=\;
\underbrace{2 \sum_k \gamma^k c_k}_{\text{the inner products}}
\;+\;
\underbrace{X\,S(X) + X^{-1} S(1/X)}_{\text{has no constant term, for any } S} ,
$$

whose constant term therefore extracts exactly the inner products. All small-polynomial openings batch
into two group elements with BDFG. Proof: $13\,\mathbb{G}_1 + 14\,\mathbb{F}$, **constant in $N$ and in
the claim count**.

**Audit.** I verified the Laurent accumulator coefficient-by-coefficient:

$$
[X^k]\big( p(X) r(1/X) + p(1/X) r(X) \big) \;=\; \sum_{i \ge k} \big( p_i r_{i-k} + p_{i-k}\, r_i \big),
$$

which is what `accumulate_ip_terms` computes, and the verifier's check

$$
\mathrm{lhs} \;\overset{?}{=}\; 2\,\mathrm{const} \;+\; \zeta\, S(\zeta) \;+\; \zeta^{-1} S(1/\zeta)
$$

matches. **No degree bound on $S$ is needed** and none is present — correctly, because
$X S(X) + X^{-1} S(1/X)$ has zero constant term for *any* $S$, so the identity still forces the claim.
The degree bound on $g$ is sound: $\zeta$ is drawn after $C_D$, and $D(X) = X^{b-1} g(1/X)$ as Laurent
polynomials forces $\deg g \le b-1$, because $D$ has no negative exponents. Shift handling is contract
(a) with the $X$-multiplied fold identity; verified.

**zk.** Not implemented. The paper's §5 gives the design (bounded-degree blinders on $h$, $g$, $S$).

**Future work.** zk; recursive verifier.

### 2.3 Vela (ePrint 2026/1438)

**Math.** The multilinear claim is the **constant coefficient** of a Laurent product. With $f_v$ the
univariate twin (coefficients $=$ the evaluation table) and
$T_{\vec{r}}(X) = \prod_k \big( (1-r_k) + r_k X^{2^k} \big)$:

$$
H(X) \;=\; f_v(X) \cdot T_{\vec{r}}(1/X),
\qquad
[X^0]\,H \;=\; \sum_i a_i\, \mathsf{eq}_i(\vec{r}) \;=\; y .
$$

The inversion-symmetric residual has zero constant coefficient **exactly when the claim holds**:

$$
D(X) \;=\; H(X) + H(1/X) - 2y,
\qquad
D_0 = 2\big( [X^0] H - y \big),
\qquad
D_k = D_{-k} .
$$

That symmetry then yields a *single* ordinary polynomial $h$ of degree $\le N - 2$ with

$$
D(X) \;=\; X\,h(X) \;+\; X^{-1} h(1/X),
\qquad
h(X) = \sum_{k \ge 1} D_k X^{k-1} .
$$

One commitment to $h$ plus evaluations at $z$ and $1/z$ certify the claim — and the fourth evaluation
$h(1/z)$ is **never transmitted**. The verifier recovers it by rearranging
$D(z) = z\,w_0 + z^{-1} w_1$:

$$
w_1 \;=\; z\Big(
\underbrace{g(z)\,T_{\vec{r}}(1/z) + g(1/z)\,T_{\vec{r}}(z)}_{\text{the verifier computes this}}
\;-\; 2Y \;-\; z\,w_0 \Big),
$$

which is *precisely the identity being proved* — so a lying prover's **recovered** $w_1$ fails the
subsequent two-point opening. That is the whole trick: the constraint is enforced by making the verifier
derive, rather than receive, the value it would otherwise have to check.

**Honk fit.** The whole batched claim set collapses into *one* identity: the shifted twin is exactly
$X^{-1} f_B(X)$, so the combined twin $g = f_A + \lambda X^{-1} f_B$ carries both chains under one $h$.

**Audit.** I checked the Laurent table's exponent range: $g$ spans $[-1,\, N-1]$ and $T_{\vec{r}}(1/X)$
spans $[-(N-1),\, 0]$, so the product spans $[-N,\, N-1]$ and the table covers exactly that —
**no truncation**. The in-place forward pass reads `table[e + 2^k]` before writing it, so it was correct
serially. The recovery formula and the degree-2 linearization (forced by bb publishing only $[1]_2$ and
$[\tau]_2$) are both sound: with

$$
E(X) \;=\; F(X) - Z(X)\,q(X) - R(X),
\qquad
Z(X) = (X-z)(X - 1/z),
$$

the check $E(\zeta) = 0$ at a random $\zeta$ drawn after $C_q$ forces $F - R = Z \cdot q$ identically.
The prover's rejection of $z \in \{0, 1, -1\}$ is necessary (the two-point interpolation divides by
$z - 1/z$) and is mirrored on the verifier.

**Deviation from the paper** (already documented): the paper's two-point check needs $[\tau^2]_2$; bb
has only $[1]_2, [\tau]_2$, so the degree-2 vanishing is linearized at one further challenge, costing
one extra $\mathbb{G}_1$ element — $3\,\mathbb{G}_1 + 5\,\mathbb{F}$ here against the paper's
$2\,\mathbb{G}_1 + 3\,\mathbb{F}$ unbatched.

**zk.** Not implemented.

**Future work.** zk; the $O(N \log N)$ Laurent expansion is inherent to never sending $h(1/z)$ — see §4
for its parallelization.

### 2.4 CHOPIN (ePrint 2026/480)

**Math.** Reads the claim as a bilinear form over the coefficient matrix of a *bivariate* twin
$f(X, Y)$, and opens it with a column restriction $w = F\,\Psi_{zR}$, one $\alpha$-fold, a
$\gamma$-batched Lagrange IPA (Mercury's Laurent accumulator, reused directly), a $\delta/z$
multi-polynomial multi-point univariate batch proof (paper Fig. 7), and one bivariate KZG opening at
$(\alpha, \beta)$. The $X$-quotient of that last opening is the **single** size-$N$ MSM, against
Mercury's fold quotient plus per-chain BDFG quotients.

**Honk fit.** The shift is a verifier-side query decomposition: chain B opens one extra column
restriction $w_2 = F\,\mathsf{shr}(\Psi_{zR})$ plus its value at zero, and the batch proof additionally
opens $w_2$ at $X = 0$. The verifier's shifted weights use the standard shift-in-the-evaluation-domain
identity

$$
\mathsf{shr}\Psi_{\mathrm{lo}}(\beta) \;=\; \beta\,\Psi_{\mathrm{lo}}(\beta) \;-\; \mathsf{eq}_{\mathrm{lo}}[M_1 - 1]\,\beta^{M_1};
$$

verified.

**Audit.** Protocol logic is faithful. The opening-point set grows to
$T = \{\alpha,\, \beta,\, 1/\beta,\, 0\}$ when shifted claims are present, and prover and verifier
compute the complement products $Z_{T \setminus S_t}(z)$ consistently. The three pairing identities
(batch proof $\pi_s$, quotient $\pi_q$, bivariate opening) are $\xi$-batched into one 3-term
multipairing against $[1]_2, [\tau]_2, [2]_2$.

**Setup caveat (§3.1).** $\tau_Y$ is fixed to the public constant $2$, so the grid derives from the
ceremony $X$-row by pointwise doublings and $[\tau_Y]_2 = 2[1]_2$. That is a **complete break of the
$Y$-direction binding** for an adversary — the module says so. Honest-prover cost is unaffected, which
is the point of the shortcut.

**Future work.** A real second trapdoor (i.e. a bivariate ceremony); zk; recursion.

### 2.5 KZH2 / KZH3 (ePrint 2025/144, Fig. 2 and App. C)

**Math.** KZH2 splits the array into $R \times C$; the commitment and its opening data are

$$
C \;=\; \sum_{i,j} f_{i,j}\,(\tau_i A_j),
\qquad
D_i \;=\; \sum_j f_{i,j} A_j,
\qquad
w \;=\; \sum_i b_i\, f_{i,\cdot} .
$$

The verifier checks three things:

$$
\underbrace{e(C, V) = \prod_i e(D_i,\, V_i)}_{\text{binding},\ V_i = \tau_i V}
\qquad
\underbrace{\mathsf{Commit}_A(w) = \sum_i b_i D_i}_{\text{row combination}}
\qquad
\underbrace{\langle w,\, a \rangle = v}_{\text{the claim}}
$$

Proof $O(\sqrt{N})$; verify one $(R+1)$-term multipairing plus two MSMs.

KZH3 adds a dimension: slice commitments $D_1$ (dim 1), contracted row commitments $D_2$ (dim 2), and
the final $d_3$-vector, bound by a two-level pairing chain. Proof drops to $O(\sqrt[3]{N})$ and
verification with it.

**Honk fit.** Rank-2 (KZH2) / rank-3 (KZH3) query decomposition; the rank-3 shift is verified term by
term in §1.4(b).

**Audit.** Both are faithful to the paper's structure. KZH3's extra chain-B data is properly bound:
$T_3'$ via $\langle T_3',\, A \rangle = \langle \mathsf{shr}(a_2),\, D_2 \rangle$ against the *same*
$D_2$ layer, and $T_3''$ via a second pairing chain on
$C_1'' = \langle \mathsf{shr}(a_1),\, D_1 \rangle$ followed by $\mathsf{Commit}_A(T_3'') = D_2''[0]$.

**Setup caveat (§3.1).** Both sample their own trapdoors in-process. Deployment needs a ceremony.

**Cost.** KZH3's proof gain widens with size ($3.3\times$ at $2^{18}$, $4.1\times$ at $2^{20}$ against
KZH2) but so does its prover ($1.6\times$ at $2^{14}$ rising to $4.8\times$ at $2^{20}$), because the
slice-commitment layer's count grows with $d_1$. The crossover where the smaller proof pays for the
prover sits well below $2^{18}$. KZH-$k$ exists to shrink an *accumulator and decider*, which is a
different objective from raw proving throughput.

### 2.6 Dory (ePrint 2020/1274, adapted)

**Math.** The commitment is AFGHO — one $\mathbb{G}_T$ element per polynomial, *transparent* (no
trapdoor):

$$
T \;=\; \prod_i e\big( P_i,\, \Gamma_{2,i} \big),
\qquad
P_i = \text{Pedersen commitment to row } i .
$$

The opening runs $\log R$ inner-pairing fold rounds. Each round sends the $\mathbb{G}_T$ cross terms of
$T$'s even/odd split and the $\mathbb{G}_1$ cross terms of every $Q$-claim, then folds

$$
T' = T \cdot T_{eo}^{\,x^{-1}} \cdot T_{oe}^{\,x},
\qquad
Q' = Q + x\,Q_{+} + x^{-1} Q_{-},
\qquad
\Gamma_2' = \Gamma_{2,\mathrm{even}} + x^{-1}\Gamma_{2,\mathrm{odd}} .
$$

The final check is one pairing $e(P_{\mathrm{final}},\, \Gamma_{2,\mathrm{final}}) = T_{\mathrm{acc}}$
plus the scalar-side closed forms.

**Audit.** The fold algebra is consistent between prover and verifier, and the $b$-side uses the same
successor-kernel closed form as the IPA backend (contract (c)), which I verified.

**Setup caveat (§3.1).** $\Gamma_{2,i} = s_i\,G_2$ with $s_i$ public deterministic scalars. bb has no
cofactor-cleared hash-to-$\mathbb{G}_2$, so the module substitutes known discrete logs — **which
destroys AFGHO binding**: knowing the $s_i$ lets an adversary find $\mathbb{G}_T$ collisions. Documented
as test-only. Real deployment needs hash-to-$\mathbb{G}_2$, at identical runtime cost.

**Cost.** The commitment does $R$ Miller loops *per polynomial*: at $2^{18}$ with 36 columns that is
$\approx 18\,000$ pairings, and it is the whole prover cost. This is inherent to AFGHO, not an
implementation artifact.

**Future work.** Hash-to-$\mathbb{G}_2$; $\mathbb{G}_T$-Pippenger $+$ cyclotomic squaring for the
verifier's $\rho$-combination (§4.5); zk.

### 2.7 IPA over Pedersen generators (Bulletproofs)

**Math.** The classical inner-product argument over hash-derived BN254 generators, run directly on
the multilinear array through the shared production core (`IPA<Curve, log n>` in `ipa.hpp`, the same
code path as the ECCVM's Grumpkin IPA). The two claim families are first merged into a **single** IPA
claim — `triple_ipa`'s reduction with the pow tensor omitted. With $F$ the $\rho$-batched unshifted
combination ($v_F$ its batched evaluation) and $F'$ the $\rho$-batched to-be-shifted combination
($v_{sh}$), the prover sends the cross-sum

$$
c \;=\; \langle F,\, b_{sh} \rangle + \langle F',\, \mathsf{eq}(\vec{u}) \rangle
$$

(label `IPA:cross_F_shift`), and the challenges $\zeta_F, \zeta_{sh}$ merge witness and tensor:

$$
A = \zeta_F F + \zeta_{sh} F',
\qquad
b = \zeta_F\, \mathsf{eq}(\vec{u}) + \zeta_{sh}\, \mathsf{shifted\text{-}eq}(\vec{u}),
\qquad
\langle A, b \rangle = \zeta_F^2 v_F + \zeta_{sh}^2 v_{sh} + \zeta_F \zeta_{sh}\, c .
$$

One core run (`compute_inner_product_proof_internal`, over the derived generators in place of the
SRS) then opens $\langle A, b \rangle$: $\log n$ rounds of

$$
L = \langle a_{\mathrm{lo}},\, G_{\mathrm{hi}} \rangle + \langle a_{\mathrm{lo}},\, b_{\mathrm{hi}} \rangle\, U,
\qquad
R = \langle a_{\mathrm{hi}},\, G_{\mathrm{lo}} \rangle + \langle a_{\mathrm{hi}},\, b_{\mathrm{lo}} \rangle\, U,
$$

with the challenge $u$ folding $a' = u^{-1} a_{\mathrm{lo}} + a_{\mathrm{hi}}$,
$b' = u\, b_{\mathrm{lo}} + b_{\mathrm{hi}}$, $G' = u\, G_{\mathrm{lo}} + G_{\mathrm{hi}}$ (the
rescaled form of the classical fold; the running scale $\prod u$ is removed from the final $G_0$ and
$a_0$, so the transcript equals the classical one). The core contributes its three prover
optimizations: 127-bit short round challenges (`get_short_challenge`), the rescaled fold above —
each generator is multiplied by the short raw challenge, never a full-width inverse — and the fused
two-round generator fold (`batch_two_round_fold`, one batch-affine pass per two rounds; see
`ipa/ELEMENT_IMPL_FOLD.md`).

**Verification** is $O(n)$: the transcript reduction (`read_inner_product_transcript_data`) checks
the group relation against the prover-claimed $G_0$, with the $b_0$ contraction in $O(\log n)$ via
the closed forms of contract (c) (`ShiftedEqPolynomial::evaluate_eq_folded` / `evaluate_folded`,
$\zeta$-weighted); the linear work is the one size-$n$ MSM certifying
$G_0 = \langle \vec{s}, \vec{G} \rangle$. That check is the Halo amortization point — the same claim
shape the production core defers to an accumulator and discharges across proofs in one batched MSM
(`verify_accumulator` / `batch_verify_accumulators`); this backend checks it inline, one MSM per
proof.

**Audit.** The fold algebra is exactly Bulletproofs; I re-derived

$$
\langle a',\, G' \rangle
\;=\; \langle a,\, G \rangle
\;+\; x^{-1} \langle a_{\mathrm{lo}},\, G_{\mathrm{hi}} \rangle
\;+\; x\, \langle a_{\mathrm{hi}},\, G_{\mathrm{lo}} \rangle
$$

and the matching $\langle a',\, b' \rangle$ relation, so that $P' = P + xL + x^{-1}R$ closes (the
core's rescaling is a change of variables on top of this). The $\zeta$-merge is the same
diagonal-plus-cross-sum identity verified for `triple_ipa` (§2.13), restricted to two tensors. $U$
is drawn after the commitments (challenge `IPA:generator_challenge`), as Bulletproofs requires; the
opening point, combined commitment and combined evaluation are absorbed into the hash buffer before
it.

**Non-hiding.** No blinders. This is the plain, non-zk IPA, not the zero-knowledge variant.

**Future work.** zk blinders; a deferred-accumulator verify entry point so many proofs share one
$G_0$ MSM (the core already ships the machinery).

### 2.8 Hyrax (ePrint 2017/1132)

**Math.** Pedersen row commitments over an $R \times C$ split; the prover sends the combined rows in
the clear and the verifier checks $\mathsf{Commit}(w_x)$ against the homomorphically folded row
commitments, then evaluates the rank-2 claim equation

$$
\sum_c w_u[c]\, a[c]
\;+\; \sum_{c \ge 1} w_s[c]\, a[c-1]
\;+\; a[C-1]\, w_{s2}[0]
\;\overset{?}{=}\; \text{expected}.
$$

Deterministic — no queries, no extra challenges.

**Deviation from the paper.** This is Hyrax's *uncompressed* variant. The paper composes a Bulletproofs
inner-product argument on top of the $\sqrt{N}$ combined row to reach an $O(\log N)$ proof; that step is
not implemented, which is why the measured proof is 566 KiB rather than a few KiB. It is also
non-hiding (the paper's Pedersen commitments carry blinders). Both are honest simplifications, but they
mean the row labelled "Hyrax" in the comparison table is not the scheme's headline design point.

**Future work.** Compose with the (now fast) IPA backend to recover the log-size proof; hiding blinders.

### 2.9 WHIR (ePrint 2024/1586)

**Math.** Per iteration: $k$ sumcheck rounds against a weight polynomial (array and weight table fold at
each challenge); commit the folded polynomial on the *halved* domain; one out-of-domain sample
$z_{\mathrm{ood}}$ with $y_{\mathrm{ood}} = \hat{g}(z_{\mathrm{ood}})$; $t$ in-domain queries into the
previous oracle, each opening a $2^k$-coset whose fold value becomes a claim about the new oracle at
$x^{2^k}$; then $\gamma$-batch the OOD and query claims into the next iteration's weight. The final
polynomial is sent in the clear and checked against the queried cosets plus the accumulated weighted sum.

The coset fold is the classical FRI identity, re-expressed in bb's Lagrange convention:

$$
A_e(x^2) = \frac{A(x) + A(-x)}{2},
\qquad
A_o(x^2) = \frac{A(x) - A(-x)}{2x},
\qquad
\mathsf{Fold}(A, \alpha)(x^2) = (1-\alpha)\,A_e(x^2) + \alpha\, A_o(x^2).
$$

The rate improves by $k-1$ bits per iteration ($\log |L|$ drops by 1 while the variable count drops by
$k$), which is what makes WHIR's query counts fall off round by round.

**Honk fit.** Batched by $\rho$ into one virtual round-0 oracle; shifted columns enter with an $x^{-1}$
factor per coset position (contract (a), in the evaluation domain). Columns of one Oink round share one
Merkle tree with **column-major interleaved leaves**, so one authentication path opens every column at a
query index — this is why WHIR's proof is 1.2 MiB rather than 3.4.

**Audit.** `fold_coset` is correct and matches bb's Lagrange fold rather than the classical
$A_e + \alpha A_o$, consistently with §1.3. The Merkle leaf grouping
$\text{leaf } j \leftarrow \{\, j + t\cdot\mathsf{num\_leaves} \,\}$ is exactly the coset
$\{\omega^j \eta^t\}$ that `fold_coset` assumes. The schedule invariant
$\log|L_i| = \log|L_{i-1}| - 1$ holds, so querying `folded_trees[i-1]` at round $i$ targets the right
oracle.

**Soundness regime.** `whir_config.hpp` implements four regimes. The default `CONJECTURED_LIST`,
$t = \lceil \lambda / r \rceil$, rests on the up-to-capacity mutual-correlated-agreement conjecture that
Crites–Stewart (ePrint 2025/2046) **disprove**. The module ships a `REPAIRED_LIST` preset solving

$$
H_q(\delta^{*}) \;=\; 1 - \rho
$$

by deterministic integer fixed-point iteration with every rounding pushed toward more queries. At BN254
the repair costs exactly one extra query per aggressive-rate round and **$+2.0\%$ proof size** at every
size from $2^{12}$ to $2^{20}$. Good work; the same repair is missing from Ligero/SwitchFold (§3.2).

**zk.** Implemented, and the only backend besides production Shplemini that has it: committed arrays
gain one variable (blinding coefficients in the high half at offset $2^m$, or $2^m + 1$ for
to-be-shifted columns, so the shift contract's zero slot survives), leaves are salted, and the opening
batches in a fresh uniformly random mask polynomial whose claimed evaluation is revealed. The blinding
budget is $q = t_{\text{round-0}} + 8$ — round-0 queries are the only openings of per-polynomial leaves,
so that many must remain information-theoretically blinded, with margin.

**Future work.** Grinding (proof-of-work) to cut query counts; stdlib recursive verifier
(`RECURSION.md` estimates $\approx 0.5$–$1.2$M gates with the shared-tree layout); UltraZK+WHIR wiring.

### 2.10 Ligero (AHIV17 / ePrint 2022/1608, Brakedown's tensor query)

**Math.** Reshape to $R \times C$, RS-encode every row to length $C \cdot 2^r$, and interleave into one
Merkle tree with leaf $j$ $=$ column $j$ of all rows. The evaluation is the tensor sandwich

$$
P(\vec{u}) \;=\; b^{\top} M\, a,
\qquad
a = \bigotimes_{i < c_{\log}} (1-u_i,\, u_i),
\qquad
b = \bigotimes_{i \ge c_{\log}} (1-u_i,\, u_i).
$$

The opening sends three combined rows ($w_u$, $w_s$, $w_{s2}$ — the rank-2 shift decomposition), then
$t$ sampled columns, each checked as

$$
\mathsf{Enc}(w_x)[j] \;=\; B_x^{\top}\,\mathsf{col}_j .
$$

**Honk fit.** Contract (b). One commit-and-query round, no folding at all — the prover-optimal corner.

**Audit.** The claim equation and the three encodings match. The rows travel unhashed with a digest
absorbed for Fiat–Shamir, so the column challenges are properly bound.

**Two findings**: §3.2 (the RS query count still uses the disproved-conjecture rule) and §3.3 (batching
doubles as the proximity test — sound only because $\vec{u}$ is a challenge).

**zk.** Not implemented; the README sketches the design (one random masking row per group, salted
leaves, random column padding).

### 2.11 SwitchFold (ePrint 2026/1489)

**Math.** Keeps Ligero's commitment and replaces the clear-text combined rows — the $O(\sqrt{N})$ term
that dominates Ligero's proof — with a recursive *code-switching* descent. The three combined rows
become the segments of one message $M$ of length $4C$, the segment index sitting in the top two
variables so that **every** claim about $M$ is a tensor. Each level commits $M$ reshaped under the next
shorter RS code, opens it at that level's queries, collapses all pending claims to a single tensor claim
with one inner-product sumcheck, and splits that claim into the column combination that becomes the next
message. This level's queries become the next level's code-switching claims. A 16-element base message is
sent in the clear.

**Instantiating over Reed–Solomon collapses two of the paper's three modules**: the RS generator-matrix
row at a domain point *is* the pow tensor, so the code-switching claim

$$
\Xi \cdot \mathsf{Enc}[C](m) \;=\; \Xi \cdot G \cdot m
$$

needs neither a commitment to the generator matrix nor the paper's accumulation scheme. Because nothing
of row length is transmitted, the proof-optimal shape moves: $C$ goes as large as BN254's 2-adicity
allows.

**Audit.** I verified the tensor bookkeeping end to end. `append_shifted_eq_terms` decomposes
$\mathsf{shr}(\mathsf{eq}(\vec{u}))$ into $\log C$ eq terms by carry position — term $k$ pins bits $< k$
to $0$ and bit $k$ to $1$, so it fires exactly when $k$ is the lowest set bit of the index, giving

$$
\Big[\prod_{j<k} u_j\Big](1 - u_k)\prod_{j>k} \mathsf{eq}_{c_j}(u_j)
\;=\; \mathsf{eq}_{c-1}(\vec{u}).
$$

Correct. The claim split (low $\log(\mathsf{sub})$ challenges bind the within-column index, high ones
give the column combination) matches the `message[c*sub_length + i]` layout under LSB-first ordering.
The derived code-switching values satisfy

$$
\sum_c \mathsf{comb}[c]\cdot \mathsf{leaf}[c] \;=\; \mathsf{Enc}(\mathsf{next})[\mathsf{index}]
$$

by linearity. Correct.

**Result.** Reproduces the paper's headline $3.5\times$ proof-size reduction over Ligero at $2^{18}$ and
does better at $2^{20}$ ($5.9\times$), because Ligero's rows grow as $\sqrt{N}$ while the descent only
adds levels. It is also $1.6\times$ smaller than WHIR at comparable verify time — the best hash-based
proof size in the suite.

**Future work.** §3.2 (query rule); §4.5 (a $4\times$ reduction in the dominant prover cost, quantified
but not implemented).

### 2.12 BrakedownHonk (ePrint 2021/1043)

`brakedown_code.hpp` implements the Spielman-style recursive linear-time code,

$$
\mathsf{Enc}_n(x):\qquad
y = xA, \quad
z = \mathsf{Enc}(y), \quad
v = zB, \quad
\mathsf{Enc}(x) = (x,\, z,\, v),
$$

with $A$ of size $n \times \lceil \alpha n \rceil$ and $c_n$ nonzeros per row, and $B$ of size
$|z| \times (rn - n - |z|)$ with $d_n$ nonzeros per row. Each level costs $n c_n + |z| d_n$
multiplications while the message contracts by $\alpha$, so the total is a geometric series: $O(n)$, no
FFT, no smooth subgroup. `brakedown_honk.hpp` runs Ligero's tensor PCS over it, making `LigeroHonk`
against `BrakedownHonk` a controlled A/B — same trace, same tensor protocol, same hasher, only the row
code and query rule differ.

**The result is the honest baseline the linear-time-code line has to beat, and it loses on all three
axes**: at $2^{14}$, $1.9\times$ slower prover, $21.6\times$ slower verifier, $24.8\times$ larger proof.
The encoder *is* faster ($1.13$–$1.22\times$ over the RS FFT, emitting $2.3\times$ fewer symbols), but
Brakedown's relative distance is $0.07$ against RS's $0.75$, so the provable interleaved test needs

$$
t \;=\; \left\lceil \frac{\lambda}{-\log_2\!\big(1 - \delta/3\big)} \right\rceil \approx 2936
\qquad\text{against}\qquad
\left\lceil \frac{\lambda}{r} \right\rceil = 50,
$$

and opening 2936 Merkle paths costs more than the $\approx 20\%$ saved on encoding.

**This correctly relocates the target for Lightning/Bolt: not encoding speed, but distance.**

**Finding (§3.4).** The base-case code is claimed MDS; it is not.

### 2.12a Bolt (ePrint 2026/310) — landed during this review

`bolt/` implements Bolt's sketched code

$$
C_H(x) \;=\; \big(\, x,\ C(Hx) \,\big),
$$

a sparse random-LDPC parity-check $H$ (of size $\alpha n \times n$, $\alpha = j/k$) compressing the
message, with the expensive base code $C$ (Reed–Solomon here) applied only to the short sketch.
`BoltHonk` runs Ligero's tensor PCS over it.

The interesting contribution is not the encoder but the **piecewise distance guarantee** (paper Claim
3.1): two distinct codewords differ either in more than $\gamma$ of the systematic stretch *or* in more
than $\delta$ of the sketch stretch. That is strictly more information than the diluted overall distance

$$
\frac{\min\!\big(\gamma,\ \delta\alpha/\rho\big)}{1 + \alpha/\rho},
$$

which averages the two and is dominated by the worse one. So the two stretches are queried
*independently*, each at the count its own distance warrants ($428$ and $50$ at $\lambda = 100$), and a
cheating prover must survive both. This is what the new `QuerySegment` hook exists for.

The module's own headline finding is one this review endorses and would have predicted: **the LDPC
distance bound collapses at BN254.** The leading term of $\omega_{q,j,k}$ is $x\ln(q-1)$, so a fixed
column degree certifies less and less distance as the field grows:

| $q$ | $j$ | $k$ | $\alpha$ | $\gamma$ |
|---|---|---|---|---|
| $2^{32}$ | 16 | 128 | 0.125 | **0.0941** |
| BN254 | 16 | 128 | 0.125 | **0.00006** |
| BN254 | 64 | 256 | 0.25 | **0.1497** |

That is the same structural observation as §2.13's small-field conclusion, arriving from the
coding-theory side: **the schemes that look best in the literature are tuned to small fields, and the
tuning does not transfer.** Any product claiming "pick your backend" over BN254 has to re-derive
parameters per backend rather than adopt the paper's.

### 2.13 Support: small_subgroup_ipa, triple_ipa, small_field

`small_subgroup_ipa` proves $\langle F,\, G \rangle = s$ for small vectors over a multiplicative
subgroup $H$ of size $m+1$, with a grand-sum polynomial $A$, a quotient $Q$, and a linear-time verifier.
The identity checked at a random $r$ is

$$
Z_H(r)\, Q(r)
\;=\; L_1(r)\, A(r)
\;+\; (r - g^{-1})\big( A(gr) - A(r) - F(r) G(r) \big)
\;+\; L_{|H|}(r)\big( A(r) - s \big).
$$

Its README contains an unusually good piece of analysis: **a fifth opening $A(1) = 0$ is required**,
because at $X = 1$ the boundary condition and the $j = 0$ recurrence collapse into the single equation

$$
A(1) + (1 - g^{-1})\big( A(g) - A(1) - F(1) G(1) \big) = 0
$$

with a one-dimensional kernel. The explicit homogeneous perturbation

$$
\delta_A(1) = \delta,
\qquad
\delta_A(g^j) = \frac{-\delta}{g-1}\ (j \ge 1),
\qquad
\delta_s = \frac{-\delta}{g-1}
$$

satisfies the identity and forges $s' = s + \delta_s$. The README also states plainly that the local
`check_consistency` is *not* sufficient on its own — soundness needs the Shplemini batched opening. That
is exactly the kind of caveat a pluggable-backend product must carry into its API contract, and it is
the model the other modules should follow.

`triple_ipa` reduces ECCVM's three structured claims (eq, shifted-eq, pow) to one ordinary IPA claim on
Grumpkin. Note the documented non-property: it does **not** assert that $F'$ is the shift of $F$; it
only proves the inner product against the shifted-eq tensor. That is fine for ECCVM's use but is a trap
for reuse.

`small_field` answers "can bb move UltraHonk to M31?" with a working M31/CM31/QM31 implementation and
measurements. The honest conclusion is worth quoting into the product question: base-field work is
$8$–$10\times$ cheaper, but **the challenge field gives most of it back** — a QM31 multiplication costs
about one $\mathbb{F}_r$ multiplication, and every sumcheck round after the first folds with QM31
challenges — so the realistic sumcheck gain is $\approx 3.7\times$, not $10\times$, and the switch is a
second proving stack, not a configuration flag.

---

## 3. Correctness findings

### 3.1 Four backends ship deliberately insecure setups

| Backend | What is substituted | Consequence |
|---|---|---|
| KZH2 | $\tau_i$ sampled in-process from a Blake3 seed | Prover knows the trapdoors; binding gone |
| KZH3 | $\mu^{(1)}_i$, $\mu^{(2)}_i$ likewise | Same |
| CHOPIN | $\tau_Y = 2$, with $[\tau_Y]_2 = 2[1]_2$ published | $Y$-direction binding gone |
| Dory | $\Gamma_{2,i} = s_i G_2$ with public $s_i$ | AFGHO binding gone ($\mathbb{G}_T$ collisions computable) |

All four are labelled `@warning Test-only` in the source and the honest-prover cost profile is
unaffected, which is the point. **This is not a defect — it is a scope boundary that the review needs to
state loudly**, because these four rows in the comparison table are performance models. Dory's fix is
mechanical (cofactor-cleared hash-to-$\mathbb{G}_2$, same runtime); KZH2/KZH3/CHOPIN need real
ceremonies, and CHOPIN's is a *bivariate* ceremony that does not exist today.

### 3.2 The RS query rule still rests on a disproved conjecture

`RSCodePolicy::num_queries` returns $t = \lceil \lambda / r \rceil$ — the RS up-to-capacity rule. WHIR
ships a `REPAIRED_LIST` preset for exactly this bound because Crites–Stewart (ePrint 2025/2046)
disproved the mutual-correlated-agreement conjecture it rests on. Their counterexamples live between the
list-decoding capacity bound and the code's capacity: they construct $u^{(0)}, u^{(1)}$ with $u^{(1)}$
farther than $\delta$ from the code such that $u^{(0)} + \lambda u^{(1)}$ is $\delta$-close for *every*
$\lambda$, whenever $\delta > 1 - H_q(\rho)$. The minimally repaired conjectures replace

$$
\delta \le 1 - \rho - \eta
\qquad\text{by}\qquad
H_q(\delta) \le 1 - \rho - \eta .
$$

The same argument applies verbatim to the interleaved-RS proximity test Ligero and SwitchFold use. At
BN254 the repair costs about one extra query ($50 \to 51$) and $\approx 2\%$ proof size — cheap to fix,
expensive to leave undocumented.

**Half of this finding was closed while the review was being written.** The shared,
distance-parameterized query hook that `PCS_CANDIDATES.md` §"Remaining candidates" item 3 asked for now
exists: `ligero/QuerySegment` plus `Code::query_plan(config)`, with `BrakedownCodePolicy` using the
*provable* rule $t = \lceil \lambda / -\log_2(1 - \delta/3) \rceil$ and `BoltCodePolicy` returning two
segments so each stretch of its piecewise-distance code is queried at the count its own distance
warrants. That is a better design than the one this review would have recommended.

What remains is narrower and still worth doing:

1. **`RSCodePolicy` has no repaired-capacity variant.** The hook is there; the corrected rule is not.
   WHIR's `compute_num_queries(REPAIRED_LIST)` already implements it in deterministic integer
   arithmetic and could be called directly.
2. **SwitchFold does not use the hook at all** — it still reads `config.num_queries` in four places
   (`switchfold.hpp:364, 484, 604, 665`), so it cannot be run over a non-RS code and will not pick up
   the repaired rule when it lands.

### 3.3 The Ligero-family proximity test rides on $\vec{u}$ being a Fiat–Shamir challenge

Ligero, SwitchFold and Hyrax use a *single* combination for both jobs: the row weights are
$\rho^i\,\mathsf{eq}_r(\vec{u}_{\mathrm{hi}})$, and the same combined rows serve as the proximity test
and the evaluation extraction. Classical Ligero/Brakedown run a *separate* proximity test with a
uniformly random vector, precisely because the eq tensor is not uniform.

This is sound here — tensor-structured randomness supports proximity gaps (Diamond–Posen, ePrint
2023/630) — **but only when $\vec{u}$ is uniformly random and drawn after the commitments.** Inside
`TransparentHonk` it is: $\vec{u}$ is the sumcheck challenge. But
`LigeroProver::prove(ck, claims, u, transcript)` accepts an arbitrary caller-supplied point, and nothing
in the signature, the doc-comment, or the README says the point must be a challenge. A user of a
pluggable-PCS product who opens at a *public, fixed* point — a perfectly reasonable thing to want from a
PCS — gets no proximity guarantee at all.

Two actions: document the precondition in the API, and (for honesty in the query count) account for the
Diamond–Posen soundness loss, which is a $\log N$-ish factor, in `num_queries`.

### 3.4 Brakedown's base-case code is not MDS

`build_base` constructs the parity block as $V_{ji} = p_j^{\,i}$ with $p_j = j + 2$ — a Vandermonde in
the *coefficient* basis. The README and the code comment both claim this is MDS with relative distance
$(r-1)/r$, which the recursion's distance argument relies on. It is not:

- $[\,I \mid V\,]$ is MDS iff every square submatrix of $V$ is nonsingular. Square submatrices of a
  Vandermonde are *generalized* Vandermonde matrices, which are nonsingular over $\mathbb{R}$ but can be
  singular over $\mathbb{F}_p$.
- The code's actual codewords are
  $\big( \mathsf{coeffs}(f),\ f(p_0),\ \ldots,\ f(p_{\mathsf{parity}-1}) \big)$ — the systematic part
  holds the *coefficients*, not evaluations, so this is not a Reed–Solomon code at all. A nonzero $f$ of
  degree $< k$ has at most $k-1$ roots, so the provable minimum distance is only

$$
d \;\ge\; 1 + \max\big(0,\ \mathsf{parity} - (k-1)\big).
$$

With $k = 32$ and $\mathsf{parity} = 24$ at $r = 1.72$, that bound is **1**, against the $d = 25$ a true
MDS $[56, 32]$ code would give.

To be precise about severity: the code's *actual* distance is probably fine — a low-weight codeword
needs a sparse coefficient vector whose polynomial vanishes on almost all of the 24 tiny evaluation
points $\{2, \ldots, 25\}$, which is very unlikely to exist. But Brakedown's entire value proposition is
a distance that holds with *provable* failure probability $2^{-100}$, and the recursion's distance
argument consumes the base case's distance as a hypothesis. An unproven base case makes the whole chain
unproven.

**Fix:** use a Cauchy matrix, $C_{ij} = 1/(y_i - z_j)$ with disjoint $\{y\}$, $\{z\}$. Cauchy matrices
are superregular (every square submatrix is nonsingular), so $[\,I \mid C\,]$ is genuinely MDS and the
systematic generalized-Reed–Solomon distance $(r-1)/r$ holds. One function, no cost change.

Also minor: `build_base` constructs a `SeededPrng` and then discards it (`static_cast<void>(prng)`) —
dead code left from an earlier random-points design.

### 3.5 Query-index collisions are not deduplicated

Every hash-based backend derives query indices as `challenge.data[0] & mask` independently per query.
Duplicates are possible, and duplicated queries do not add soundness. This is standard practice in
deployed FRI implementations and the i.i.d.-sampling bound $\rho^t$ is the correct one for sampling
*with* replacement, so this is not a bug — but it is worth a comment, because a reader checking $t$
distinct queries against the analysis will not find them.

### 3.6 `index_from_challenge_bounded` has documented modulo bias

`ligero::detail::index_from_challenge_bounded` uses `% bound` for non-power-of-two codeword lengths
(Brakedown). The comment states the bias is below $2^{-40}$ at these sizes, which is right, but the
constant is worth recomputing if the module is ever used with a much smaller security parameter.

### 3.7 Mercury/Vela/CHOPIN assert rather than reject on the shift precondition

`BB_ASSERT_EQ(column[0], fr::zero(), ...)` is a prover-side assertion, which is stripped in release
builds. That is correct for these backends (an honest prover always satisfies $P[0] = 0$, and a
dishonest one only breaks its own proof), but if the PCS layer is ever exposed to untrusted
*prover-side* input in a product, it should be a real check.

---

## 4. Efficiency: what was wrong, what was fixed, what remains

### 4.1 The shared defect: the $\rho$-batched linear combination was serial everywhere

Every backend starts by forming the $\rho$-combination of the $\approx 41$ claims,

$$
F \;=\; \sum_i \rho^{\,i} f_i ,
$$

which was written as:

```cpp
for (each claim) {
    for (size_t k = 0; k < n; ++k) { array[k] += rho_power * column[k]; }
    rho_power *= rho;
}
```

That is a **serial $\mathsf{num\_claims} \times N$ pass** — $10.7$M multiply-adds at $2^{18}$.
Production Gemini does the same job with `bb::add_scaled_batch`, which splits the *output* range across
threads and loops the sources inside each chunk. The experimental backends could not reuse it directly
because it requires `Polynomial` operands and they hold `std::vector<fr>`, so the loop was open-coded —
serially — in Mercury, Vela, CHOPIN, KZH2, KZH3, Dory and the IPA backend.

The size of this is not incidental: at $2^{18}$ it is a **large fraction of the prover gap between the
constant-size pairing backends and the KZG baseline**, and KZG only avoids it because production Gemini
already parallelizes the identical work.

**Fix.** `utils/batch_accumulate.hpp` adds the `std::vector`-backed counterpart:

```cpp
struct ScaledTerm { const fr* source; fr scalar; bool shifted; };
void accumulate_scaled(std::span<fr> out, std::span<const ScaledTerm> terms);
```

one parallel pass over the output range, terms looped inside each chunk (disjoint writes, no locking,
each output cache line touched once per chunk rather than once per term). Applied to Mercury, Vela
(three sites, including the $\alpha$-batch), CHOPIN, KZH2, KZH3, Dory and IPA.

### 4.2 The IPA backend was not using MSMs at all

`pedersen_ipa.hpp`'s round loop computed the cross terms as

```cpp
for (size_t t = 0; t < half; ++t) {
    left  += generators[2 * t]     * a[2 * t + 1];
    right += generators[2 * t + 1] * a[2 * t];
}
```

— **$2N$ individual elliptic-curve scalar multiplications** across the whole reduction, where an MSM is
called for. bb's production `ipa.hpp` uses `pippenger_unsafe` on contiguous lo/hi halves for exactly
this. The even/odd split was what blocked it: Pippenger needs contiguous point and scalar spans.

**Fix, first pass.** Converted the backend from even/odd to **lo/hi** splits, which is also the
classical Bulletproofs presentation. Both cross-term operands become contiguous, so each is one
`pippenger_unsafe` call; the generator fold became a parallel pass followed by one `batch_normalize`
back to affine.

**Fix, second pass.** Even with Pippenger cross terms, a hand-rolled round loop pays a full-width
scalar multiplication per generator in the fold ($x^{-1}$ is a 254-bit scalar) — across the rounds
that is $\approx n$ full scalar muls, several times the cost of all the round MSMs combined. The
backend now routes its (merged, see §2.7) claim through the production core
`IPA::compute_inner_product_proof_internal` via a generator-span overload, so it inherits the core's
short-challenge schedule, rescaled fold, and fused two-round `batch_two_round_fold` instead of
re-implementing the reduction. The verifier correspondingly consumes
`read_inner_product_transcript_data` with the $\zeta$-weighted closed forms of §1.4(c)
(`ShiftedEqPolynomial`), and the hash-derived generator vector is memoized per size so re-creating a
verification key does not re-run the O(n) hash-to-curve derivation. All four `PedersenIpaTest` cases
(including the two negative tests) and `IpaHonkTest` pass.

### 4.3 Six more serial loops

- **SwitchFold's combined-row build** was the serial version of a loop Ligero already parallelizes —
  $3 \times 36 \times N$ multiply-adds, $\approx 28$M at $2^{18}$. Now split over the *column* range so
  each thread owns a disjoint slice of all three segments (no private accumulators, no merge).
- **SwitchFold's `inner_product_round`** was a serial sumcheck round; now `parallel_for_range` with
  private accumulators, matching WHIR's `sumcheck_round_univariate`.
- **Hyrax's combined-row build** — same loop, same fix.
- **Vela's `laurent_product`** is the $O(N \log N)$ shift-and-add expansion, $\log N$ serial passes over
  a $2N$ table. Parallelized with an explicit double buffer: the pass computes

$$
\mathsf{table}'[e] \;=\; (1 - r_k)\,\mathsf{table}[e] \;+\; r_k\,\mathsf{table}[e + 2^k],
$$

  reading ahead of the write cursor, so an in-place parallel version would race (the same trap already
  fixed once in WHIR's array fold). This is the largest single contributor to Vela's
  $1.70\times \to 1.18\times$.
- **Mercury's column combination**
  $h[\mathsf{row}] = \sum_j \mathsf{eq}_{\mathrm{lo}}[j]\, f[\mathsf{row}\cdot b + j]$ — an $O(N)$
  serial pass; now parallel over rows.
- **Mercury's Laurent accumulator** `accumulate_ip_terms` is $O(b^2) + O(\mathsf{rows}^2) = O(N)$. Each
  shift writes only $S[\mathsf{shift}-1]$, so the shifts are **independent** — parallelized directly
  with no accumulator merge. (Work per shift falls as $\mathsf{len} - \mathsf{shift}$, so index-ordered
  chunks are unbalanced by under $2\times$, which is not worth a custom schedule.)

### 4.4 Measured effect

**Methodology, and its limits.** The review machine was shared with another build/edit session for part
of the measurement window, so absolute wall-clock is not comparable between the committed baseline
(`pcs_bench_results.json`, single run) and the post-fix runs. Every figure below is therefore a **ratio
to the `ultra_honk_kzg_prove` baseline measured in the same run**, which cancels machine-condition
drift. **`LigeroHonk` is the control**: its code is untouched by this review, so the movement in its
ratio is the noise floor. The post-fix figures come from a 5-repetition median taken while the machine
was quiet; a later 3-repetition run under contention moved the control by 26% and was discarded rather
than reported. Proof bytes and verifier work are unchanged by every fix, so these are pure prover-side
wins — the proof-size column of §2's table is unaffected.

All at $2^{18}$, prove, ratio to same-run KZG:

| Backend | before | after | change |
|---|---|---|---|
| **Ligero** (control, code untouched) | $0.70\times$ | $0.72\times$ | $+3\%$ ← **noise floor** |
| CHOPIN | $1.05\times$ | $\mathbf{0.81\times}$ | $-23\%$ |
| Mercury | $1.27\times$ | $\mathbf{1.14\times}$ | $-10\%$ |
| Vela | $1.70\times$ | $\mathbf{1.18\times}$ | $-30\%$ |
| SwitchFold | $2.97\times$ | $\mathbf{2.50\times}$ | $-16\%$ |
| IPA (Pedersen) | $82.0\times$ | $\mathbf{\approx 2.2\times}$ | **~37× faster** |
| Hyrax, KZH2, KZH3, Dory | — | — | unchanged within noise |

The IPA row has two stages. The MSM fix of §4.2 took it from $82.0\times$ to $\approx 5.6\times$
($44\,726\ \mathrm{ms} \to 3351\ \mathrm{ms}$ in comparable single-run conditions — large enough to
survive the contention; CPU time in the noisy run was 3391 ms). The rebuild on the production core
(§2.7, §4.2) then took it to $\approx 2.2\times$: same-run 3-repetition medians of 2026-08-06 give
IPA/KZG prove ratios $1.8$–$2.5\times$ across $2^{12}$–$2^{20}$ ($1837/819\ \mathrm{ms}$ at
$2^{18}$), verify $3596 \to 121\ \mathrm{ms}$ at $2^{18}$ (single MSM instead of two, and no
per-verify generator derivation), proof $16.4 \to 12.0$ KiB (one reduction chain instead of two).
The IPA entries in `pcs_bench_results.json` / `pcs_report.html` are these medians transported into
the committed baseline's frame via the same-size same-run KZG ratio. The four
unchanged backends are exactly the ones whose provers are dominated by $N$ elliptic-curve scalar
multiplications per column (and $R$ Miller loops per column for Dory), where the $\rho$-batch was never
more than a few percent — that they did *not* move is a consistency check on the attribution, not a
disappointment.

**Two results worth calling out.** CHOPIN now proves **faster than the KZG baseline** while producing a
9.2 KiB proof against KZG's 13.3, and Vela's prover penalty for the smallest-proof-in-the-suite drops
from $1.70\times$ to $1.18\times$. Both were previously read as "pay prover time for proof size"; that
trade is now much weaker, which changes their standing in §7.3.

### 4.5 Quantified but not implemented

**SwitchFold's weight tables do $4\times$ more work than necessary.** `WeightTerm::accumulate_table`
builds the full $2^{\log C + 2}$ tensor for every term, but `select_segment` appends two eq factors at
the *constants* $0$ and $1$. An eq factor at a constant is a point mass — $(a, b) = (1, -1)$ zeroes the
high half, $(a, b) = (0, 1)$ the low half — so three quarters of the tensor are identically zero, and
those zeros are still computed (as multiplications by zero) and added. With $\approx 150$ of the
$\approx 169$ level-0 terms being per-query, segment-restricted `pow_weight` terms, restricting the
build to the live quarter is a straight $4\times$ on the dominant cost, worth $\approx 195$ ms of the
remaining prover at $2^{18}$. The general form: partition each term's variables into *free* and *pinned*
(exactly one of $a$, $a+b$ zero), build the tensor over the free variables only, and scatter into the
strided support. That also collapses the $\log C$ shifted-eq terms from $\log C$ full tables to
$\approx 1$.

**Dory's verifier spends its time in $\mathbb{G}_T$ exponentiation.** The $\rho$-combination
$\prod_p T_p^{\rho^p}$ is 41 square-and-multiply exponentiations of 254 bits in $\mathbb{F}_{q^{12}}$ —
$\approx 15\,600$ `fq12` operations. Cyclotomic squaring (a standard $\approx 2\times$ for
$\mathbb{G}_T$) plus a bucket-method multi-exponentiation over $\mathbb{G}_T$ (another $\approx 4\times$
at these counts) would take the 362 ms verifier down substantially. Neither primitive exists in bb's
`fq12` today.

**Ligero's verifier re-encodes three full rows.** Three size-$C 2^r$ FFTs per verification. Direct
evaluation at only the $t$ queried positions is *not* a win at the current shape ($\approx 4.9$M against
$\approx 3.3$M multiplications), so this is correctly left alone — recorded so it is not
re-investigated.

**Mercury's binomial-fold recurrence is parallelizable and was left alone.** The recurrence

$$
q[k - b] \;=\; f[k] + \alpha\, q[k]
$$

looks sequential but decomposes into $b$ independent chains, one per residue class mod $b$ — so it
parallelizes over $b = \sqrt{N} \approx 512$ chains at $2^{18}$. Together with the two remaining $O(N)$
passes (the `numerator` construction in round 4, and the `divide_by_linear` calls, which are genuinely
sequential recurrences of length $N$) this is the last few percent of Mercury's prover. Recorded for
completeness; the measured payoff is small relative to the risk of touching the fold.

**IPA can go substantially further.** The current fix uses one MSM per cross term per round. Production
`ipa.hpp` additionally fuses round pairs (`batch_two_round_fold`, deferring the SRS fold so the second
round's MSMs run against the pre-fold SRS) and uses a rescaled fold
$G' = u\,G_{\mathrm{lo}} + G_{\mathrm{hi}}$ with 127-bit short challenges, so each `batch_mul`
multiplies by a half-width scalar. Adopting both would be another large constant factor on top of the
$13.3\times$.

**KZH2/KZH3/Hyrax/Dory commit costs are inherent.** Their provers are dominated by $N$ elliptic-curve
scalar multiplications per committed column (and $R$ Miller loops per column for Dory). No
implementation change moves these; they are the schemes' cost.

---

## 5. Zero-knowledge: the decisive gap

| Backend | zk status |
|---|---|
| Gemini+Shplonk+KZG | **Implemented.** Sparse $2d$-coefficient masking poly $+$ `small_subgroup_ipa` for Libra; proof sketch in `SHPLEMINI_ZK_MASKING.md` with a Lean development for the linear algebra |
| WHIR | **Implemented.** Blinded high half, salted Merkle leaves, batched random mask polynomial, $q = t_{\text{round-0}} + 8$ blinding coefficients |
| Mercury, Vela, CHOPIN, KZH2, KZH3, Dory, Hyrax, IPA, Ligero, SwitchFold, Brakedown, Bolt | **None.** Non-hiding commitments, unmasked openings |

Thirteen of fifteen backends leak the witness. For a product this is not a footnote — a user who selects
"Mercury" and gets a non-zk proof has a security failure, not a performance trade.

The work is also **not uniform**, which matters for scoping:

- *Pairing/KZG family* (Mercury, Vela, CHOPIN): bounded-degree blinders on the small polynomials plus
  degree slack in the SRS. Each paper gives the design; each needs its own argument because each opens a
  different set of auxiliary polynomials.
- *$\sqrt{N}$ / $\sqrt[3]{N}$ families* (KZH2, KZH3, Hyrax, Dory): the combined rows are sent **in the
  clear**. Hiding them requires either hiding Pedersen commitments plus a zk inner-product argument
  (Hyrax's actual design), or masking rows — a structural change, not a wrapper.
- *Hash family* (Ligero, SwitchFold, Brakedown, Bolt): WHIR's recipe (mask polynomial $+$ salted leaves
  $+$ blinded region) transfers in shape, but Ligero's combined rows and SwitchFold's descent messages
  are transmitted in the clear and each needs its own masking-row design and query-count accounting.

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
WHIR's: $\approx 73$ gates per Poseidon2 permutation in Ultra, $\approx 5.4$M gates for a naive layout,
$\approx 0.5$–$1.2$M with the shared-per-round trees, deduplicated openings and minimal Fiat–Shamir
absorption that the native implementation already has.

For the other backends the picture divides cleanly:

- **Pairing backends** (Mercury, Vela, CHOPIN, KZH2, KZH3) recurse like KZG: the verifier is an MSM plus
  a pairing, and bb already has the `biggroup`/`bigfield` machinery. Vela and Mercury are the cheapest
  (constant-size proof, $O(\log N)$ field work). KZH2/KZH3 need $\sqrt{N}$ / $\sqrt[3]{N}$ in-circuit
  MSM scalars, which is the point of the KZH line for accumulation.
- **Dory** needs in-circuit $\mathbb{F}_{q^{12}}$ arithmetic and $\log R$ $\mathbb{G}_T$ operations —
  expensive and not currently available.
- **Hash backends** need in-circuit Poseidon2 Merkle verification (the `Poseidon2MerkleHasher` variant
  exists precisely for this) plus the sumcheck and weight bookkeeping. WHIR and SwitchFold are the
  plausible ones; Ligero's $t \cdot P \cdot R$ opened values make the in-circuit cost scale with the
  proof, which is MiB-scale.
- **IPA/Hyrax** verify in $O(N)$ group operations, which does not recurse at reasonable cost without the
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
- the trade-offs a user would actually select on are real and large — at $2^{18}$ the suite spans
  7.8 KiB to 2.7 MiB in proof size, 3.5 ms to 3.6 s in verification, and $0.72\times$ to
  $\approx 46\times$ the KZG prover;
- the abstraction keeps absorbing new structure rather than fighting it. Three things landed *during*
  this review — a linear-time code, a code-policy parameterization, and a piecewise-distance query plan
  for a sketched-LDPC backend — and none required a change above the PCS layer. An abstraction that
  survives contact with genuinely different cryptography is the thing you cannot fake.

That is a strong result and it is the hard part of the product. A `#[backend(vela)]` attribute on a Noir
program is not a fantasy.

### 7.2 What blocks it today

Ranked by how much they cost to fix, most binding first.

1. **Zero-knowledge on 13 of 15 backends** (§5). Not a shared layer; $\approx 3$ distinct designs plus
   per-scheme arguments. This is the gating item.
2. **Four backends have no real setup** (§3.1). Dory needs hash-to-$\mathbb{G}_2$ (mechanical).
   KZH2/KZH3 need a ceremony. CHOPIN needs a *bivariate* ceremony that does not exist.
3. **Recursive verifiers exist for one backend family** (§6). A product whose proofs cannot be verified
   inside another Aztec circuit is a much narrower product.
4. **The RS query rule still rests on a disproved conjecture** (§3.2). Now small work — the hook exists —
   but it is a correctness-of-parameters issue, which is exactly the class of bug a user selecting a
   backend cannot audit for themselves.
5. **API preconditions are undocumented** (§3.3, §3.7). A pluggable PCS is used by people who did not
   read the module; "the opening point must be a Fiat–Shamir challenge" has to be in the type system or
   the doc-comment, not only in the caller.
6. **Parameters do not transfer from the papers.** Bolt's LDPC distance collapses from
   $\gamma = 0.094$ to $\gamma = 0.00006$ moving from $2^{32}$ to BN254 (§2.12a); WHIR's
   repaired-conjecture penalty scales as $1/\log_2 q$; the M31 study finds the challenge field eats most
   of the small-field win (§2.13). A backend menu has to own per-backend parameter derivation at *its*
   field, not cite the paper's table.
7. **Proof-size-driven backends need proof compression.** The hash family sits at 0.77–2.7 MiB. ePrint
   2025/1446 (`zip`) is tracked in `PCS_CANDIDATES.md` and attacks exactly this.

### 7.3 A roadmap that would make it shippable

**Phase 1 — make the honest set small and real.** Pick three backends to productize rather than
fifteen: **Gemini+Shplonk+KZG** (has zk, has recursion, is the default), **Vela** (smallest proof,
fastest verifier, and — importantly — a real ceremony already available, since it reuses bb's existing
KZG commitment key), and **WHIR** (has zk, transparent, recursion designed). Ship those three behind the
selector and gate the remaining twelve as research backends. This converts the zk gap from "13 backends"
to "1 backend" (Vela), which is a scoped piece of work rather than a programme.

Note that Vela's case strengthened during this review: its prover penalty was $1.70\times$ KZG and is
now $1.18\times$ (§4.4), so "smallest proof and fastest verifier" no longer costs much on the third axis.

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
| `mercury/mercury.hpp` | Parallel $\rho$-batch (both chains); parallel column-combination $h$; parallel Laurent accumulator (shifts are independent) |
| `mercury/vela.hpp` | Parallel $\rho$-batch, parallel $\alpha$-batch, parallel double-buffered `laurent_product` |
| `mercury/chopin.hpp` | Parallel $\rho$-batch (both chains) |
| `kzh/kzh.hpp`, `kzh/kzh3.hpp` | Parallel $\rho$-batch |
| `dory/dory.hpp` | Parallel $\rho$-batch (array side; the row-commitment side is EC work) |
| `hyrax/hyrax.hpp` | Column-parallel combined-row build |
| `switchfold/switchfold.hpp` | Column-parallel combined-row build; parallel sumcheck round |
| `pedersen_ipa/pedersen_ipa.hpp` | lo/hi split; Pippenger cross terms; batched generator fold; reversed fold multipliers |
| `PCS_REVIEW.md`, `PCS_REVIEW.html` | This report |

## Appendix B: papers audited against

Mercury [2025/385], Vela/Carina [2026/1438], CHOPIN [2026/480], KZH/KZH-k [2025/144],
Dory [2020/1274], Bulletproofs/IPA, Hyrax [2017/1132], WHIR [2024/1586],
Crites–Stewart proximity-gap disproof [2025/2046], Ligero [2022/1608], Brakedown [2021/1043],
SwitchFold [2026/1489], Bolt [2026/310], BDFG multi-point batching [2020/081],
Diamond–Posen tensor proximity [2023/630].
