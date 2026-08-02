# WHIR: a hash-based multilinear PCS for Honk

This module implements WHIR [1] as a multilinear polynomial commitment scheme over the BN254
scalar field, with the interfaces and conventions of the Honk proving stack. It provides:

- a transparent commitment (a Merkle root over a Reed–Solomon codeword) with no structured
  reference string and no elliptic-curve operations;
- an opening argument for batched multilinear evaluation claims at a common point $\vec u$,
  including Honk's shifted-by-one claims, matching the claim shape produced by sumcheck
  (the same shape `Shplemini` consumes);
- an optional zero-knowledge mode (hiding commitments, leakage-free openings);
- plug-in Merkle hashers: Blake3s (fast native proving) and Poseidon2 (recursion-friendly).

Security properties, established by the cited analyses and the parameter schedule of §6:

- **(C1) Binding.** A commitment binds the prover to a unique polynomial of the claimed size,
  under collision resistance of the Merkle hash and the proximity soundness of the WHIR IOP
  in the chosen soundness regime (§6).
- **(C2) Evaluation soundness.** A verifier accepting an opening at point $\vec u$ with values
  $(v_j)$ implies, with error at most $2^{-\lambda}$ plus negligible field terms, that the
  committed polynomials evaluate to $(v_j)$ at $\vec u$.
- **(C3) Zero knowledge** (zk mode, §8). Commitments are statistically hiding and the opening
  transcript is simulatable given only the claimed evaluations.

## 1. Why WHIR — the hash-based PCS landscape

Hash-based PCS constructions descend from two lineages that converged in 2023–2024:

- **RS proximity testing:** FRI [4] folds an RS codeword by 2 per round on a fixed-rate
  schedule; the verifier does $O(\lambda \cdot \log^2 d)$ hashing. STIR [5] restores the rate
  after each fold by re-committing the folded polynomial on a domain only half the size,
  so the code's rate improves by $2^{k-1}$ per round and later rounds need fewer queries.
- **Tensor/linear codes:** Ligero and Brakedown commit to a matrix of coefficients and open
  via random row combinations; proofs are $O(\sqrt d)$ but the multilinear structure comes
  for free. BaseFold [6] showed that FRI-style folding *is* multilinear partial evaluation if
  the fold challenges are shared with a sumcheck, turning FRI into a multilinear PCS.

WHIR [1] combines both: it is STIR's rate-improving iteration applied to *constrained*
Reed–Solomon codes, where a sumcheck over a weighted hypercube sum rides along with the
folding. One protocol simultaneously proves proximity and a multilinear evaluation claim.
At 100-bit security, the WHIR paper reports 63 KiB proofs and 360 µs verification for
$d = 2^{22}$, the smallest proofs and fastest verification in the hash-based family [1, §7].
The 2026 SoK [3] surveys the full design space and confirms WHIR as the current frontier of
the RS line. That is the scheme worth benchmarking against the Gemini+Shplonk+KZG stack.

Trade-offs to keep in mind when reading benchmarks (§9): hash-based commitment is
$O(N \log N)$ field ops (FFT) plus $O(N)$ hashing versus one $O(N)$ MSM for KZG; proofs are
tens of KiB versus KZG's constant few hundred bytes; verification is hashing-bound rather
than pairing-bound; there is no trusted setup and (with a suitable hash) a plausibly
post-quantum security story.

## 2. Notation and parameters

| Symbol | Meaning | Typical value |
|---|---|---|
| $\mathbb F$ | BN254 scalar field (2-adicity 28, $\lvert\mathbb F\rvert \approx 2^{254}$) | — |
| $m$ | number of multilinear variables; polynomial size $n = 2^m$ | 12–20 |
| $r_0$ | initial log-inverse-rate; codeword length $N_0 = 2^{m + r_0}$ | 1–3 |
| $k$ | folding factor per iteration (`folding_factor_bits`) | 4 |
| $\lambda$ | target security level (`security_level_bits`) | 100 |
| $M$ | number of fold-and-commit iterations | $\lfloor (m - m_{\text{fin}})/k \rfloor$ |
| $m_i$ | variables remaining at iteration $i$: $m_i = m - ik$ | — |
| $N_i$ | oracle domain size at iteration $i$: $N_i = N_0 / 2^i$ | — |
| $r_i$ | log-inverse-rate at iteration $i$: $r_i = r_0 + i(k-1)$ | — |
| $t_i$ | number of in-domain queries at iteration $i$ (§6) | 8–100 |
| $m_{\text{fin}}$ | variables of the final clear polynomial (`final_poly_bits`) | 4–6 |
| $L_i$ | evaluation domain: the order-$N_i$ subgroup $\langle \omega_i \rangle$, $\omega_{i+1} = \omega_i^2$ | — |
| $\gamma, \rho$ | claim-combination / polynomial-batching challenges | — |

Parameters live in `WhirConfig`; the derived per-iteration schedule (rates, query counts,
domain sizes) is computed by `WhirConfig::create` and asserted against the field's 2-adicity.

## 3. The basis convention

Honk stores a polynomial as an array $a = (a_0, \dots, a_{n-1})$ and uses it in two roles at
once: sumcheck reads $a$ as the evaluations of a multilinear $P$ on the hypercube
($P(\operatorname{bits}(i)) = a_i$, low variable first), while the commitment scheme reads
$a$ as univariate coefficients $A(X) = \sum_i a_i X^i$. Gemini ties the two roles together
through the stride-2 fold

$$\operatorname{fold}(a, \alpha)_j \;=\; (1-\alpha)\,a_{2j} + \alpha\, a_{2j+1},$$

which is simultaneously (i) partial evaluation of $P$ at $X_1 = \alpha$ and (ii) the
coefficient array of $(1-\alpha) A_e + \alpha A_o$, where $A(X) = A_e(X^2) + X A_o(X^2)$.
This module keeps exactly that convention (the kernel is the shared `fold_stride2`), so WHIR
plugs into Honk without any basis conversion. Two consequences:

**Pointwise foldability.** Given codeword values at $\pm x$, the folded codeword value at
$x^2$ is

$$\operatorname{fold}(A,\alpha)(x^2) \;=\; (1-\alpha)\,\frac{A(x) + A(-x)}{2} \;+\; \alpha\,\frac{A(x) - A(-x)}{2x}. \tag{3.1}$$

**Weighted-sum form of every claim.** Both claim types the protocol manipulates are sums
over the hypercube of $a$ against a weight that is a product of per-variable affine factors:

$$\text{MLE claim } P(\vec u) = v: \qquad \sum_{b \in \{0,1\}^m} \operatorname{eq}(\vec u, b)\, a_b = v, \qquad \operatorname{eq}(\vec u, X) = \prod_j \big((1-u_j)(1-X_j) + u_j X_j\big),$$

$$\text{univariate claim } A(y) = \beta: \qquad \sum_{b \in \{0,1\}^m} \operatorname{pow}_y(b)\, a_b = \beta, \qquad \operatorname{pow}_y(X) = \prod_j \big((1-X_j) + X_j\, y^{2^{j-1}}\big). \tag{3.2}$$

The second identity is $A(y) = \sum_i a_i y^i$ with $y^i = \prod_j (y^{2^{j-1}})^{b_j}$ for
$i = \operatorname{bits}^{-1}(b)$. Univariate evaluations of an oracle (out-of-domain samples
and in-domain queries) therefore become sumcheck-compatible weights on the *same* array the
prover already holds. Weight terms are represented uniformly by `WeightTerm` (a scalar
coefficient and one affine factor per variable); tables and partial evaluations are shared
between the two types.

**Shift contract.** Honk's shifted-by-one claims use the same contract as `Shplemini`: a
to-be-shifted polynomial has $a_0 = 0$, its shift has array $a'_i = a_{i+1}$, and on the
codeword $A'(x) = A(x)/x$ pointwise (all $x \in L_0$ are nonzero). The verifier derives
shifted-oracle values from unshifted leaf openings with one multiplication by $x^{-1}$; no
separate commitment exists for a shift.

The protocol below is WHIR [1] with the fold matrix $(1-\alpha, \alpha)$ in place of
$(1, \alpha)$ and weights adjusted per (3.2). Both folds are affine lines through the value
pair, so the mutual-correlated-agreement analysis [1, §4] applies unchanged (proximity gaps
hold for affine lines [7]).

## 4. Protocol

### 4.1 Commitment

`WhirCommitmentKey::commit_group`: polynomials committed together (a "group" — e.g. all
columns of one Honk round) are padded to $n = 2^m$, FFT'd onto $L_0$ ($N_0 = 2^{m+r_0}$
points, natural order), and bound by one Merkle tree whose leaf $j$ hashes, column-major,
every column's values on fold-coset $j$,

$$\operatorname{coset}(j) = \{\, j + t \cdot N_0/2^k \;:\; t \in [2^k] \,\}, \qquad j \in [N_0 / 2^k],$$

so one opening authenticates everything needed to fold every column of the group at index
$j$. The commitment is the root. The prover retains the coefficient arrays and tree
(`WhirGroupData`); a single polynomial is a group of one.

### 4.2 Opening statement

The prover opens claims $\{P_j(\vec u) = v_j\}$ for committed $P_j$ and shifted claims
$\{P^{\text{shift}}_l(\vec u) = v^{\text{shift}}_l\}$, all at one point $\vec u$ (the
sumcheck challenge in Honk). After all roots are in the transcript the verifier sends
batching challenge $\rho$; both sides form the virtual initial oracle and claim

$$F(x) = \sum_j \rho^{\,j} A_j(x) \;+\; \sum_l \rho^{\,n_u + l}\, \frac{A_l(x)}{x}, \qquad
\sigma_0 = \sum_j \rho^{\,j} v_j + \sum_l \rho^{\,n_u+l} v^{\text{shift}}_l,$$

with initial weight list $W_0 = \{\operatorname{eq}(\vec u, \cdot)\}$. Soundness of this
batching is by correlated agreement [7], as in batch-FRI. The single-polynomial case is the
same with an empty $\rho$-combination.

### 4.3 Iteration $i = 0, \dots, M-1$ (claim on $\hat g_i$, oracle on $L_i$; $\hat g_0 = F$)

Implemented by `WhirProver::prove_round` / mirrored in `WhirVerifier::verify_round`.

1. **Sumcheck ($k$ rounds).** For the claim
   $\sum_{b} W_i(b)\, \hat g_i(b) = \sigma$, round $j$ sends the degree-2 univariate
   $h_j(X) = \sum_{b} \hat g_i(\alpha_{<j}, X, b)\, W_i(\alpha_{<j}, X, b)$ as three
   evaluations $h_j(0), h_j(1), h_j(2)$; the verifier checks $h_j(0) + h_j(1) = \sigma_{j-1}$
   and sets $\sigma_j = h_j(\alpha_j)$ for challenge $\alpha_j$. Prover-side, array and
   weight table both fold by `fold_stride2` at $\alpha_j$.
2. **Fold-and-commit.** The prover commits $g_{i+1} = \operatorname{fold}(g_i, \vec\alpha)$
   — coefficient array of length $2^{m_{i+1}}$ — as a fresh codeword on
   $L_{i+1} = L_i^2$ of size $N_{i+1} = N_i/2$ (coset-grouped leaves as in §4.1).
3. **Out-of-domain sample.** Verifier sends $z \leftarrow \mathbb F$; prover replies
   $y^{\text{ood}} = A_{g_{i+1}}(z)$. This pins the codeword to a unique nearby polynomial
   (list-decoding disambiguation [1, §5]).
4. **In-domain queries.** The verifier derives $t_i$ indices
   $q_1, \dots, q_{t_i} \in [N_i / 2^k]$ from transcript challenges. For each, the prover
   opens leaf $q_s$ of *every* constituent tree of the round-$i$ oracle (all $P_j$ trees at
   $i=0$; the single $g_i$ tree for $i>0$). The verifier checks the Merkle paths,
   assembles the virtual value (RLC and $x^{-1}$ scaling at $i = 0$), and computes
   $y_s = \operatorname{fold}(g_i, \vec\alpha)(x_s)$ at $x_s = \omega_i^{\,2^k q_s}$ by $k$
   applications of (3.1) over the coset.
5. **Combine.** Verifier sends $\gamma$. The next claim on $\hat g_{i+1}$ is

$$W_{i+1} = W_i(\vec\alpha, \cdot) \;\cup\; \{\gamma \cdot \operatorname{pow}_z\} \;\cup\; \{\gamma^{1+s} \cdot \operatorname{pow}_{x_s}\}_{s=1}^{t_i}, \qquad
\sigma_{i+1} = \sigma_k + \gamma\, y^{\text{ood}} + \sum_{s} \gamma^{1+s} y_s,$$

where $\sigma_k$ is the running sumcheck value after step 1 and $W_i(\vec\alpha,\cdot)$ is
the partial evaluation of every term of $W_i$ at the sumcheck challenges. Prover-side the
combined weight *table* over $\{0,1\}^{m_{i+1}}$ is updated in $O(2^{m_{i+1}})$ per new term
(tensor accumulation, `WeightTerm::accumulate_table`).

### 4.4 Final phase

After iteration $M-1$ the claim concerns $\hat g_M$ with $m_{\text{fin}} = m - Mk$
variables, oracle on $L_M$. The prover sends the array of $\hat g_M$ in the clear
(`WHIR:final_poly`, $2^{m_{\text{fin}}}$ field elements). The verifier:

1. derives $t_M$ final query indices, opens leaves of the $g_M$ tree, and checks every
   opened codeword value against direct evaluation of the clear polynomial,
   $A_{g_M}(\omega_M^{\,\text{idx}})$;
2. checks the outstanding claim directly:
   $\sum_{b \in \{0,1\}^{m_{\text{fin}}}} W_M(b)\, \hat g_M(b) = \sigma_M$, evaluating each
   weight term of $W_M$ (partially evaluated at all past sumcheck challenges) over the small
   cube.

### 4.5 Transcript schedule

All messages go through the standard Honk transcript (Poseidon2 Fiat–Shamir for
`NativeTranscript`); challenges are full-width field elements. Query indices are reduced
from field challenges modulo the (power-of-two) index space, which introduces no usable
bias at 254 bits. Opened leaves and Merkle paths travel in the proof stream without
Fiat-Shamir absorption (`send_unhashed_to_verifier`): they are bound by the roots absorbed
before the query indices were drawn, so re-absorbing them adds no soundness and would
dominate verifier hashing.

| Order | Label | Direction | Content |
|---|---|---|---|
| 1 | `WHIR:root_<poly>` | P→V | commitment roots (or from earlier Honk rounds) |
| 2 | `WHIR:rho` | V→P | batching challenge (batched case) |
| per iter $i$, round $j$ | `WHIR:sumcheck_<i>_<j>` | P→V | $h_j(0), h_j(1), h_j(2)$ |
| ″ | `WHIR:alpha_<i>_<j>` | V→P | fold/sumcheck challenge |
| ″ | `WHIR:root_g<i+1>` | P→V | folded-oracle root |
| ″ | `WHIR:z_ood_<i>` / `WHIR:y_ood_<i>` | V→P / P→V | OOD point / answer |
| ″ | `WHIR:query_<i>_<s>` | V→P | index challenges |
| ″ | `WHIR:answers_<i>` | P→V | opened leaves + Merkle paths |
| ″ | `WHIR:gamma_<i>` | V→P | combination challenge |
| end | `WHIR:final_poly` | P→V | $2^{m_{\text{fin}}}$ coefficients |
| ″ | `WHIR:final_query_<s>` | V→P | final index challenges |
| ″ | `WHIR:final_answers` | P→V | opened leaves + paths |

## 5. Correctness

Completeness is arithmetic: (3.1) makes the verifier's coset folds equal the prover's
`fold_stride2` outputs; (3.2) makes OOD/query answers equal weighted hypercube sums; the
sumcheck telescopes $\sigma_0 \to \sigma_M$. The unit tests exercise each identity in
isolation (`whir.test.cpp`, `rs_code.test.cpp`, `weights.test.cpp`) and end-to-end
completeness plus tamper rejection for every message type.

## 6. Soundness regimes and the parameter schedule

Per-query soundness depends on the proximity regime assumed for RS codes of rate
$\rho_i = 2^{-r_i}$ [1, §5; 3]:

| Regime (`WhirSoundness`) | Distance tested | Bits per query at rate $2^{-r}$ | Query count $t_i$ |
|---|---|---|---|
| `UNIQUE_DECODING` | $(1-\rho)/2$ | $-\log_2\!\big(\tfrac{1+2^{-r}}{2}\big)$ | $\lceil \lambda / \text{bits} \rceil$ |
| `PROVABLE_LIST` (Johnson) | $1 - \sqrt\rho$ | $r/2$ | $\lceil 2\lambda / r \rceil$ |
| `REPAIRED_LIST` | $\delta^*$: $H_q(\delta^*) = 1-\rho$ | $-\log_2(1-\delta^*)$ | $\lceil \lambda / \text{bits} \rceil$ |
| `CONJECTURED_LIST` (capacity) | $1 - \rho$ | $r$ | $\lceil \lambda / r \rceil$ |

The rate improves every iteration, $r_{i+1} = r_i + (k-1)$, because degree divides by $2^k$
while the domain only halves — this is what makes later rounds cheap and total queries much
lower than FRI at the same $\lambda$. OOD samples (1 per iteration; 0 in unique decoding)
and the $\gamma$/$\rho$/sumcheck field terms add $O(\text{poly}(2^m)/|\mathbb F|)$ error,
negligible at 254 bits for $\lambda \le 128$. Grinding (query-phase proof-of-work) is a
standard further reduction of $t_i$; it is not implemented, and all benchmark numbers are
grinding-free.

### The up-to-capacity conjecture is false; the repaired regime

Crites-Stewart ([eprint 2025/2046](https://eprint.iacr.org/2025/2046)) disprove the
correlated-agreement, mutual-correlated-agreement, and list-decodability up-to-capacity
conjectures that `CONJECTURED_LIST` rests on (for WHIR specifically, Conjecture 4.12 — the
mutual correlated agreement adaptation of [BCI+23] Conjecture 8.4). The counterexamples live
between the list-decoding capacity bound and the code's capacity: they construct words
$u^{(0)}, u^{(1)}$ with $u^{(1)}$ farther than $\delta$ from the code such that
$u^{(0)} + \lambda u^{(1)}$ is $\delta$-close for *every* $\lambda$, whenever
$\delta > 1 - H_q(\rho)$ (with $H_q$ the $q$-ary entropy). Their minimally repaired
conjectures replace $\delta \le 1-\rho-\eta$ by $H_q(\delta) \le 1-\rho-\eta$.

`REPAIRED_LIST` implements that repair: it tests the largest distance $\delta^*$ with
$H_q(\delta^*) = 1-\rho$, i.e. per-query error
$1-\delta^* \approx \rho + h_2(\delta^*)/\log_2 q$. At $\log_2 q \approx 254$ (BN254 Fr) the
entropy penalty is a fraction of a bit per query: at $r=2$ it costs 51 queries instead of 50,
at $r=5$ 21 instead of 20, and from $r=8$ the ceiling absorbs it. The schedule solver runs in
the same deterministic Q192/Q64 fixed-point style as `UNIQUE_DECODING`, with every rounding
(the $\log_2 q = 253$ floor, the penalty division, a $2^{-32}$ bits-per-query guard) pushed
toward more queries. The default preset remains `CONJECTURED_LIST` for comparability with
deployed FRI/STIR systems; the benchmark reports both (`whir_honk_repaired_*`), and the
measured deltas are ~2% proof size and noise-level verify time.

Worked example (`CONJECTURED_LIST`, $\lambda = 100$, $m = 20$, $r_0 = 2$, $k = 4$,
$m_{\text{fin}} = 4$, so $M = 4$):

| $i$ | $m_i$ | $N_i$ | $r_i$ | $t_i$ |
|---|---|---|---|---|
| 0 | 20 | $2^{22}$ | 2 | 50 |
| 1 | 16 | $2^{21}$ | 5 | 20 |
| 2 | 12 | $2^{20}$ | 8 | 13 |
| 3 | 8 | $2^{19}$ | 11 | 10 |
| final | 4 | $2^{18}$ | 14 | 8 |

## 7. Costs

Counts for a single polynomial of size $n = 2^m$ (batched opening shares all per-iteration
costs across the batch; only round-0 leaf openings scale with the number of trees):

| Cost | Amount |
|---|---|
| Prover commit | one size-$N_0$ FFT + $N_0/2^k$ leaf hashes + $N_0/2^k$ tree hashes |
| Prover open | per iteration: $O(2^{m_i})$ sumcheck/weight field ops, size-$N_{i+1}$ FFT, tree build |
| Proof size | round 0: $t_0 \cdot 32\,\text{B} \cdot \big(C \cdot 2^k + \sum_{\text{groups}} d \cdot \log_2(N_0/2^k)\big)$ ($C$ total columns, $d$ digest fields); rounds $i \ge 1$: $t_i \cdot 32\,\text{B} \cdot \big(2^k + d \log_2(N_i/2^k)\big)$; + $3 \cdot 32\,\text{B} \cdot k M$ (sumcheck) + $2^{m_{\text{fin}}} \cdot 32\,\text{B}$ |
| Verifier | $\sum_i t_i \cdot (\log_2(N_i/2^k) + 1)$ hashes + $O\big(\sum_i t_i (2^k + m_i)\big)$ field ops |

For the §6 example the proof is $\approx 90$ KiB and the verifier computes $\approx 1{,}700$
hashes; the corresponding Shplemini+KZG proof is $\approx 4$ KiB with a pairing-bound
verifier. Measured numbers live in `whir.bench.cpp` results (see the benchmark section of
the PR).

## 8. Zero-knowledge mode

zk mode makes (C3) hold with three mechanisms, following the standard RS-IOP recipe
(Aurora [8]; STARK zk note [2]; binary-field analogue [9]):

1. **Salted leaves.** Leaf $j$ hashes $(\text{salt}_j \,\|\, \text{values})$ with a fresh
   random 256-bit salt, revealed only for opened leaves: unopened leaves are hiding, so the
   Merkle tree is a statistically hiding vector commitment.
2. **Codeword blinding.** Each committed array is extended by one variable
   ($m \to m+1$): the low half is the payload, and $q = t_0 + 8$ uniformly random
   coefficients are placed in the high half (offset $2^m$ for unshifted polynomials;
   offset $2^m + 1$ for to-be-shifted polynomials, keeping slot $2^m$ zero so the shift
   contract and shifted-claim values are unchanged). Every claim at $\vec u$ lifts to
   $(\vec u, 0)$, whose eq-weight vanishes on the blinded half, so claimed values are
   unchanged. Any $t_0$ codeword positions of a blinded polynomial are jointly uniform
   (Vandermonde), and round-0 queries are the only openings of per-polynomial leaves.
3. **Mask polynomial.** The opening prover commits one extra fully-random blinded
   polynomial $R$, reveals its (independent, uniform) claimed value at $(\vec u, 0)$, and
   includes it in the $\rho$-combination. The virtual oracle $F$ is then uniform, so every
   message derived from $F$ — all sumcheck univariates, folded oracles, OOD answers, the
   final clear polynomial — is independent of the witness.

Simulation argument (honest-verifier, statistical): sample $F$ uniform conditioned on the
public claim $\sigma_0$; run the honest protocol on $F$ for all post-$\rho$ messages;
answer round-0 leaf openings with values jointly uniform conditioned on consistency with
$F$ (valid because opened positions number at most $t_0 < q$); salts cover all unopened
leaves. Fiat–Shamir preserves zk. Costs: one extra committed polynomial, and each first-round
codeword doubles in length ($N_0 \to 2 N_0$).

The Honk sumcheck phase has its own masking (Libra); zk-WHIR covers the PCS phase only.
§10 discusses the combination.

## 9. Benchmark methodology

`whir.bench.cpp` measures commit, open, verify and proof size for
$m \in \{12, \dots, 20\}$, both hashers, and $r_0 \in \{1, 2, 3\}$, against
Shplemini+KZG on identical claim sets (same polynomial count and sizes as an Ultra trace).
Native numbers on this machine carry a caveat: the arm64 build disables the x86 field
assembly, which slows field ops (and thus KZG MSMs) more than it slows Blake3 hashing;
relative conclusions should be re-checked on the remote benchmark machine.

## 10. UltraHonk integration

`whir_honk.hpp` (`WhirHonk`) is UltraHonk with WHIR as the PCS: `UltraFlavor`'s
arithmetization, relations, and sumcheck are unchanged, every commitment is a WHIR Merkle
root, and the opening phase is one batched WHIR run.

- **Commitments and VK.** The transparent `WhirHonk::VerificationKey` holds circuit metadata
  and the precomputed polynomials' roots. Witness commitments follow Oink's round schedule
  with the identical challenge labels, so the derived-polynomial computations are the shared
  `OinkProver` static helpers and the relation parameters bind to hash commitments.
- **Opening phase.** `WhirHonk::prove` assembles the §4.2 statement from
  `polynomials.get_unshifted()` / `get_to_be_shifted()` and
  `sumcheck_output.claimed_evaluations` — the identical inputs Shplemini consumes — with the
  roots already transcript-bound (`Claims::send_roots = false`). The verifier replaces the
  Shplemini batch-mul + pairing with `WhirVerifier` (no `PairingPoints`, no SRS anywhere).
- **Per-round tree groups.** Polynomials committed in the same Honk round share one tree
  (five groups: the 28 precomputed columns, wires, counts+w_4, lookup_inverses, z_perm), so
  a round-0 query opens one path per commitment round rather than per polynomial, and a
  shifted claim reuses its column's opened values with $x^{-1}$ scaling instead of a second
  opening. RECURSION.md §2 quantifies the effect.
- **ZK flavors.** Witness masking rows and Libra sumcheck masking carry over unchanged;
  the PCS phase uses §8. The `SmallSubgroupIPA` sub-protocol reduces to standard opening
  claims, which fold into the same batched WHIR statement.

## 11. Recursion

An in-circuit WHIR verifier is hashing-dominated: with Poseidon2 leaves/nodes it needs
$\approx \sum_i t_i \log_2 N_i$ in-circuit permutations (order $10^3$ for §6 parameters)
plus $O(10^4)$ field ops — no non-native group arithmetic at all, unlike the recursive
Shplemini verifier whose cost is a large biggroup MSM. The quantitative analysis and a
stdlib design sketch live in `RECURSION.md`.

## References

1. G. Arnon, A. Chiesa, G. Fenzi, E. Yogev, *WHIR: Reed–Solomon Proximity Testing with
   Super-Fast Verification*, [ePrint 2024/1586](https://eprint.iacr.org/2024/1586).
2. U. Haböck, Al Kindi, *A note on adding zero-knowledge to STARKs*,
   [ePrint 2024/1037](https://eprint.iacr.org/2024/1037).
3. *SoK: Hash-Based Polynomial Commitments and Low-Degree Tests: From FRI to Basefold,
   STIR, and WHIR*, [ePrint 2026/1367](https://eprint.iacr.org/2026/1367).
4. E. Ben-Sasson, I. Bentov, Y. Horesh, M. Riabzev, *Fast Reed–Solomon Interactive Oracle
   Proofs of Proximity* (FRI), ICALP 2018.
5. G. Arnon, A. Chiesa, G. Fenzi, E. Yogev, *STIR: Reed–Solomon Proximity Testing with
   Fewer Queries*, [ePrint 2024/390](https://eprint.iacr.org/2024/390).
6. H. Zeilberger, B. Chen, B. Fisch, *BaseFold: Efficient Field-Agnostic Polynomial
   Commitment Schemes*, [ePrint 2023/1705](https://eprint.iacr.org/2023/1705).
7. E. Ben-Sasson, D. Carmon, Y. Ishai, S. Kopparty, S. Saraf, *Proximity Gaps for
   Reed–Solomon Codes*, [ePrint 2020/654](https://eprint.iacr.org/2020/654).
8. E. Ben-Sasson, A. Chiesa, M. Riabzev, N. Spooner, M. Virza, N. Ward, *Aurora:
   Transparent Succinct Arguments for R1CS*, [ePrint 2018/828](https://eprint.iacr.org/2018/828).
9. B. E. Diamond, *Zero-Knowledge Polynomial Commitment in Binary Fields*,
   [ePrint 2025/1015](https://eprint.iacr.org/2025/1015).

## Correspondence map

| Spec item | Code |
|---|---|
| parameters, query schedule (§2, §6) | `whir_config.hpp`: `WhirConfig`, `WhirRoundSchedule` |
| coset-grouped Merkle tree (§4.1) | `merkle_tree.hpp`: `MerkleTree<Hasher>` |
| hashers | `merkle_tree.hpp`: `Poseidon2MerkleHasher`, `Blake3sMerkleHasher` |
| RS encode, domains, coset fold (3.1) | `rs_code.hpp` |
| weight terms, tables, partial evals (3.2) | `weights.hpp`: `WeightTerm` |
| commit (§4.1) | `whir.hpp`: `WhirCommitmentKey::commit` |
| prover iterations (§4.3–4.4) | `whir.hpp`: `WhirProver` |
| verifier (§4.3–4.4) | `whir.hpp`: `WhirVerifier` |
| batched statement (§4.2) | `whir.hpp`: `WhirProver::Claims`, `WhirVerifier::Claims` |
| zk mode (§8) | `WhirConfig::zk` schedule fields, blinding in `WhirCommitmentKey`, salting in `merkle_tree.hpp` |
| UltraHonk integration (§10) | `whir_honk.hpp`: `WhirHonk` |
| benchmarks (§9) | `whir.bench.cpp` |
| recursion analysis (§11) | `RECURSION.md` |
