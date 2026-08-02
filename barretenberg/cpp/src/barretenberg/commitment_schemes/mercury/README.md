# Mercury: a constant-proof-size multilinear PCS over KZG

This module implements Mercury [1] over BN254 KZG: the pairing-based corner of the backend set,
competing directly with Gemini+Shplonk+KZG (`Shplemini`) at the same trust model (powers-of-τ SRS,
pairing check) while producing **constant-size opening proofs** with a **linear, FFT-free prover**.
It plugs into UltraHonk through the same `TransparentHonk` shell as WHIR and Ligero.

Claims, established by [1] and the multi-point batching of [2]:

- **(C1)** Commitments are standard univariate KZG commitments of the evaluation array (bb's
  array-is-both-bases convention is Mercury's native convention: hypercube values are univariate
  coefficients).
- **(C2)** An accepted opening implies the claimed multilinear evaluations under the same
  assumptions as KZG (AGM/q-DLOG), with proof size independent of the polynomial size.

## 1. Scheme (single claim)

Following [1] and the exposition of [3]: with $N = 2^n$, $t = \lceil n/2 \rceil$, $b = 2^t$, the
array reshapes to $\text{rows} = N/b$ by $b$ columns ($j = \text{row} \cdot b + \text{col}$, columns
indexed by the low $t$ variables). The claim $\tilde f(\vec u) = v$ factors through
$h(X) = \sum_j \operatorname{eq}_j(\vec u_{\text{lo}}) f_j(X)$ (the column-combined polynomial,
$f_j$ the column polynomials) as four sub-claims:

1. **Fold correctness**: $f(X) = q(X)(X^b - \alpha) + g(X)$ for verifier challenge $\alpha$, where
   $g$'s coefficients are $f_j(\alpha)$ — one binomial division, checked at a random $\zeta$
   through the quotient $H(X) = (f(X) - (\zeta^b-\alpha)q(X) - g(\zeta))/(X-\zeta)$ and the pairing
   $e(C_f - (\zeta^b-\alpha)C_q - g(\zeta)[1]_1,\ [1]_2) = e(C_H, [\tau-\zeta]_2)$.
2. **Degree bound** $\deg g < b$: via $D(X) = X^{b-1}g(1/X)$ (committable only if the bound holds).
3. **$\langle P_{u_{\text{lo}}}, g \rangle = h(\alpha)$** and 4. **$\langle P_{u_{\text{hi}}}, h \rangle = v$**:
   two inner products, batched by $\gamma$ into one polynomial identity
   $\sum \gamma^k (p_k(X) r_k(1/X) + p_k(1/X) r_k(X)) = 2\sum\gamma^k c_k + X S(X) + X^{-1} S(1/X)$,
   whose constant term extracts the inner products.

All small-polynomial openings ($g, h, S, D$ at subsets of $\{\zeta, 1/\zeta, \alpha\}$) are batched
into two group elements with the BDFG multi-point argument [2]. The verifier evaluates the tensor
polynomials in $O(\log N)$: $P_{u}(X) = \prod_i (u_i X^{2^i} + 1 - u_i)$.

## 2. Honk batching: two chains

**Unshifted claims** batch homomorphically: the verifier computes
$C_A = \sum_i \rho^i C_i$ itself and one chain proves the combined claim
$\langle F_A, \operatorname{eq}(\vec u)\rangle = \sum_i \rho^i v_i$.

**Shifted claims** run a second chain on the virtual polynomial $f' = B(X)/X$ (the shifted
array, where $B = \sum_l \rho^{n_u+l} B_l$ over the to-be-shifted polynomials, each with zero
constant term). $f'$ is never committed: its fold identity multiplies through by $X$,

$$B(X) = Q(X)(X^b - \alpha) + X g_B(X), \qquad Q(X) := X q_B(X),$$

so the prover commits $Q$ and the pairing check uses the *unshifted* homomorphic commitment
$C_B$ with a $\zeta$-scaled remainder term:
$e(C_B - (\zeta^b-\alpha)C_Q - \zeta g_B(\zeta)[1]_1,\ [1]_2) = e(C_{H_B}, [\tau-\zeta]_2)$.
Both chains then use the standard $\operatorname{eq}$ tensors — no shifted-query machinery.

The chains share $\alpha, \gamma, \zeta, \beta, z$; their four inner products share one $S$; all
seven small polynomials share one BDFG batch; and the three pairing identities (two fold checks +
BDFG) combine with one more challenge into a single pairing check. Proof: 13 G1 + 14 F
(~1.2 KiB in bb field-element serialization), constant in $N$ and in the claim count.

Costs: prover $\approx 4N + O(\sqrt N)$ MSM scalars (the $q, H$ pairs of both chains) plus $O(N)$
field operations with no FFT; verifier: one claim-count-sized commitment RLC, $O(\log N)$ field
work, one pairing check. zk (not implemented; design): standard KZG-style masking of $h, g, S$
with bounded-degree blinders, as in [1, §5].

## Siblings in this directory

Two later schemes reuse this module's commitment key, tensor evaluation, Laurent accumulator and
exact division, and are benchmarked against it in `../pcs_report.html`:

- **`vela.hpp`** — Vela [4]. Same constant-term encoding, but the inversion-symmetric residual
  $D(X) = H(X) + H(1/X) - 2y$ collapses the opening to a single auxiliary polynomial $h$ opened at
  $z$ and $1/z$; the fourth evaluation $h(1/z)$ is recovered by the verifier from the Laurent
  identity rather than transmitted. Smallest proof and fastest verifier in the suite, at the cost
  of $O(N\log N)$ structured field work.
- **`chopin.hpp`** — CHOPIN [5]. Reads the claim as a bilinear form over a *bivariate* twin, so the
  fold consistency check is a cross-evaluation rather than a univariate division. That removes the
  high-degree quotient commitment and leaves a single size-$N$ MSM, at one extra pairing.

Both are described in `../PCS_CANDIDATES.md` §1–2, including the deviations bb's two-element G2
SRS forces on them.

## References

1. L. Eagen, A. Gabizon, *MERCURY: A multilinear Polynomial Commitment Scheme with constant proof
   size and no prover FFTs*, [ePrint 2025/385](https://eprint.iacr.org/2025/385).
2. D. Boneh, J. Drake, B. Fisch, A. Gabizon, *Efficient polynomial commitment schemes for multiple
   points and polynomials*, [ePrint 2020/081](https://eprint.iacr.org/2020/081).
3. J. Xie, Y. Guo, *Mercury Notes* (sec-bit/mle-pcs), the protocol exposition this implementation
   follows round by round.
4. Y. Zhang, *Vela and Carina: Fast Pairing-Based Multilinear Polynomial Commitments from Reciprocal
   Polynomials*, [ePrint 2026/1438](https://eprint.iacr.org/2026/1438).
5. J. Belohorec, P. Hubáček, A. Kalsta, K. Mašková, *CHOPIN: Optimal Pairing-Based Multilinear
   Polynomial Commitments from Bivariate KZG*, [ePrint 2026/480](https://eprint.iacr.org/2026/480).
