# Ligero: a tensor-code multilinear PCS for Honk

This module implements a Ligero-style [1, 2] multilinear polynomial commitment scheme over the
BN254 scalar field: the prover-optimal corner of the hash-based design space, complementing WHIR
(`../whir/`) at the other end (WHIR: small proofs, iterated folding; Ligero: one commit-and-query
round, larger proofs, no proximity folding at all). It shares the WHIR module's Merkle tree,
RS-encoding, transcript, and Honk integration infrastructure, and plugs into UltraHonk through the
same `TransparentHonk` shell.

Claims, established by the interleaved-RS proximity analysis of [1] and the tensor-query argument
of [2]:

- **(C1) Binding**: the Merkle root binds the prover to an interleaved codeword within unique
  decoding distance, hence to unique row polynomials, under hash collision resistance.
- **(C2) Evaluation soundness**: an accepted opening implies the claimed multilinear evaluations
  with error 2^{-λ} + negl (rate-dependent per-query soundness, §4).

## 1. Scheme

A polynomial of $2^m$ evaluations is reshaped row-major into an $R \times C$ matrix $M$
($R = 2^{m - c_{\log}}$ rows, $C = 2^{c_{\log}}$ columns, the low $c_{\log}$ variables indexing
columns). Its multilinear evaluation at $\vec u$ is the tensor sandwich

$$P(\vec u) \;=\; b^\top M\, a, \qquad a = \bigotimes_{i < c_{\log}} (1-u_i, u_i) \in \mathbb F^C, \quad b = \bigotimes_{i \ge c_{\log}} (1-u_i, u_i) \in \mathbb F^R.$$

**Commit** (`LigeroCommitmentKey::commit_group`): every row of every polynomial in the group is
RS-encoded to length $C' = C \cdot 2^{r}$ (rate $2^{-r}$), and one Merkle tree commits the
interleaved codeword with leaf $j$ = column $j$ of all rows (`MerkleTree` with arity 1). The
commitment is the root.

**Open** (batched, at one point $\vec u$, claims ordered as in the WHIR module): after the
batching challenge $\rho$, the prover sends three combined rows

$$w_u = B_u^\top M_{\text{stack}}, \qquad w_s = B_s^\top M_{\text{stack}}, \qquad w_{s2} = B_{s2}^\top M_{\text{stack}},$$

where $M_{\text{stack}}$ stacks all committed rows and the per-row batching vectors are assembled
from $b$: for the $i$-th unshifted claim on a polynomial, $\rho^i b$ on its rows into $B_u$; for
the $l$-th shifted claim, $\rho^{n_u+l} b$ into $B_s$ and $\rho^{n_u+l} \operatorname{shr}(b)$
into $B_{s2}$ ($\operatorname{shr}$ prepends a zero). The rows travel unhashed; a digest of
$(w_u \| w_s \| w_{s2})$ is absorbed for Fiat-Shamir. The verifier then samples $t$ column indices;
the prover opens each group's tree at those leaves (unhashed, root-bound).

**Verify**: (i) each opened column satisfies, for $x \in \{u, s, s2\}$,
$\operatorname{Enc}(w_x)_j = B_x^\top \operatorname{col}_j$ — one full RS encoding of each $w_x$
plus $t \cdot 3$ inner products; (ii) Merkle paths; (iii) the claim equation

$$\sum_i \rho^i v_i + \sum_l \rho^{n_u+l} v^{\text{shift}}_l \;=\; \langle w_u, a\rangle + \langle w_s, \operatorname{shr}(a)\rangle + a_{C-1} \cdot w_{s2}[0].$$

## 2. Shifted claims: the rank-2 decomposition

Honk's shifted claim is $\operatorname{shift}(P)(\vec u) = \sum_{j \ge 1} P_j\, \operatorname{eq}_{j-1}(\vec u)$,
an inner product with the shifted eq vector $Q'$, which is not a tensor. Reshaped, with
$j = rC + c$:

$$Q'_{r,c} = \begin{cases} b_r\, a_{c-1} & c \ge 1 \\ b_{r-1}\, a_{C-1} & c = 0,\ r \ge 1 \\ 0 & c = r = 0, \end{cases}
\qquad\text{i.e.}\qquad Q' = b \otimes \operatorname{shr}(a) \;+\; \operatorname{shr}(b) \otimes (a_{C-1}\, e_0).$$

A rank-2 sum of tensors: the first term reuses the $b$-direction (pairs with
$\operatorname{shr}(a)$ through $w_s$), the second needs the one extra combined row $w_{s2}$ and
contributes only its entry at column 0. Unlike the WHIR/Shplemini univariate contract, no
zero-constant-term condition on the polynomial is needed.

## 3. Parameters

`LigeroConfig`: $m$, column log-size $c_{\log}$, log-inverse-rate $r$, security λ, query count
$t = \lceil \lambda / r \rceil$ (up-to-capacity conjecture for interleaved RS proximity, matching
the WHIR module's default regime for comparability; the proven Ligero bound needs a constant
factor more queries [1]). `default_log_num_cols` picks $C^* \approx \sqrt{t \cdot P \cdot 2^{m-1}}$
($P$ = total committed polynomials), the minimizer of proof size
$3C + t \cdot P \cdot R$ field elements plus paths.

## 4. Costs (vs WHIR, same Ultra batch)

Prover: one row-encoding pass ($P \cdot R$ FFTs of size $C'$), one tree of $P \cdot 2^{m+r}$
values, and a single $3$-scalar pass over the data for the combined rows — no folding rounds, no
iterated FFTs, no inner sumcheck. Verifier: $3$ size-$C'$ encodings, $t(3 P R)$ multiplications,
$t$ column hashes + paths. Proof: $3C + t \cdot P \cdot R$ field elements + $t$ paths — a few MiB
at the Ultra column count, the price of the square-root regime. zk (not implemented; design): one
random masking row per group added to each combined row, salted leaves, and random column padding,
as in [1, §4].

## References

1. S. Ames, C. Hazay, Y. Ishai, M. Venkitasubramaniam, *Ligero: Lightweight Sublinear Arguments
   Without a Trusted Setup*, CCS 2017 / [ePrint 2022/1608](https://eprint.iacr.org/2022/1608).
2. A. Golovnev, J. Lee, S. Setty, J. Thaler, R. Wahby, *Brakedown: Linear-time and field-agnostic
   SNARKs for R1CS*, [ePrint 2021/1043](https://eprint.iacr.org/2021/1043).
3. Comparison context: [SoK: Hash-Based Polynomial Commitments](https://eprint.iacr.org/2026/1367),
   [Mercury: constant proof size, no prover FFTs](https://eprint.iacr.org/2025/385) — the
   pairing-based candidate documented as follow-up work.
