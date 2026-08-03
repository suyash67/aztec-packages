# UltraHonk PCS backends: implemented set and evaluated candidates

Fifteen PCS backends are implemented and benchmarked against each other on identical circuits through
the shared `TransparentHonk` shell (`transparent_honk.hpp`); the sweep lives in
`whir/whir_honk.bench.cpp` and its output in `pcs_bench_results.json` / `pcs_report.html`.

Bold entries were added in the 2026 state-of-the-art pass described below. All figures are from one
sweep on an Apple M4 Pro (arm64, `DISABLE_ASM=1`); single-run wall clock, so treat differences under
a few percent as noise.

| Backend | Family | Trust | Prove 2^18 | Verify 2^18 | Proof 2^18 |
|---|---|---|---|---|---|
| **Vela** | pairing, reciprocal | SRS | 929 ms | **3.5 ms** | **7.8 KiB** |
| **CHOPIN** | pairing, bivariate KZG | SRS | 573 ms | 4.7 ms | 9.2 KiB |
| Mercury | pairing, constant-size | SRS | 691 ms | 4.2 ms | 9.4 KiB |
| Gemini+Shplonk+KZG | pairing | SRS | **545 ms** | 4.7 ms | 13.3 KiB |
| IPA (Pedersen) | DL | transparent | 44 726 ms | 3596 ms | 16.4 KiB |
| **KZH3** | pairing, ∛N opening | SRS | 25 167 ms | 79.7 ms | 55.3 KiB |
| Dory | pairing, two-tier | transparent | 6655 ms | 362 ms | 94.3 KiB |
| KZH2 | pairing, √N opening | SRS | 7837 ms | 216 ms | 183.3 KiB |
| Hyrax | DL, Pedersen rows | transparent | 6087 ms | 243 ms | 566.3 KiB |
| **SwitchFold** | hash, code switching | transparent | 1619 ms | 13.6 ms | 772.9 KiB |
| WHIR | hash, RS folding | transparent | 591 ms | 13.5 ms | 1219.7 KiB |
| Ligero | hash, tensor √N | transparent | **382 ms** | 20.6 ms | 2692.6 KiB |
| Bolt | hash, sketched code | transparent | see §7 | see §7 | see §7 |
| Brakedown | hash, Spielman code | transparent | see §7 | see §7 | see §7 |

The two linear-time-code backends are compared against Ligero at 2^14 in §7 rather than here: their
query counts are one to two orders of magnitude larger, so a 2^18 row would be dominated by that
single parameter and would invite the wrong comparison.

Reading the table: **Vela is the smallest proof and the fastest verifier in the suite**, KZG is the
fastest pairing prover and Ligero the fastest prover overall, and SwitchFold now occupies the
hash-based corner Ligero used to own on proof size. The full 2^12–2^20 sweep, including the
Poseidon2 variants of every hash backend and the repaired-conjecture WHIR preset, is in
`pcs_bench_results.json`; `pcs_report.html` renders it (regenerate both with `pcs_report.py`).

## 1. Vela (eprint 2026/1438)

Vela reads the multilinear claim `f(r) = y` as the constant coefficient of the Laurent product
`H(X) = f_v(X)·T_r(1/X)`, where `f_v` is the univariate twin whose coefficients are the evaluation
table. The inversion-symmetric residual `D = H(X) + H(1/X) - 2y` has zero constant coefficient
exactly when the claim holds, and its symmetry yields a single ordinary polynomial `h` with
`D = X·h(X) + X⁻¹·h(1/X)`. One commitment to `h` plus evaluations at `z` and `1/z` certify the
claim, and the fourth evaluation `h(1/z)` is never transmitted — the verifier recovers it from the
Laurent identity.

The batched Honk claim set collapses into a single identity. The shifted twin of a to-be-shifted
polynomial is exactly `X⁻¹·f_B(X)` (the division is exact because the shift contract forces a zero
constant term), so the combined twin `f_A + λ·X⁻¹·f_B` carries both chains under one `h`. `h` is
built by `μ` shift-and-add passes over the Laurent coefficient table — `O(N log N)` field
operations, no FFT — reusing Mercury's commitment key, tensor evaluation and exact division.

**Deviation from the paper.** The paper's two-point opening at `{z, 1/z}` checks
`e(C, [1]₂) = e(π, [Z(τ)]₂)` with `Z` of degree two, which needs `[τ²]₂`. bb's verifier SRS
publishes only `[1]₂` and `[τ]₂`, so the degree-2 vanishing check is linearized at one further
challenge. That costs one extra G1 element: 3 G1 + 5 F here against the paper's 2 G1 + 3 F for a
single unbatched claim.

Vela confirms the paper's trade on all three axes. At 2^18 it produces the smallest proof in the
suite (7.8 KiB against Mercury's 9.4 and KZG's 13.3) and the fastest verification (3.5 ms against
4.2 and 4.7), paid for with prover-side structured field work: 929 ms against Mercury's 691. The
`O(N log N)` Laurent expansion is the price of never sending the fourth evaluation.

## 2. CHOPIN (eprint 2026/480)

CHOPIN expresses the multilinear evaluation as a bilinear form over the coefficient matrix of a
bivariate twin `f(X, Y)` and opens it with a column restriction, one α-fold, a γ-batched Lagrange
IPA (Mercury's Laurent accumulator, reused directly from `mercury/`), a δ/z multi-polynomial
multi-point univariate KZG batch proof (paper Fig. 7), and one bivariate KZG opening at `(α, β)`.
The X-quotient of that last opening is the **single** size-N MSM — against Mercury's fold quotient
plus per-chain BDFG quotients, which in the batched Honk setting are four size-N MSMs. The Honk
shift is a verifier-side query decomposition: chain B opens one extra column restriction under
shifted weights plus its value at zero.

The test-only bivariate SRS keeps `τ_X` as bb's real ceremony trapdoor and fixes `τ_Y = 2`, so the
grid derives from the ceremony X-row by pointwise doublings and `[τ_Y]₂ = 2·[1]₂` is public. The
honest prover's cost profile is unchanged by this: commitments and the opening quotient are flat
Pippenger MSMs over the materialized grid, exactly matching the univariate-KZG commit path Mercury
is benchmarked against.

The paper's claimed two-fold prover speedup does not fully materialize here. CHOPIN proves in 573 ms
at 2^18 against Mercury's 691 (0.83x) and 2342 ms against 2446 at 2^20 (0.96x), with the expected
one-extra-pairing verifier cost (4.7 ms against 4.2) and a marginally smaller proof (9.2 against
9.4 KiB). The reason is that the paper counts opening MSMs, whereas in the batched Honk setting a
large and growing share of the prover is the shared commitment and sumcheck phases, which CHOPIN
does not touch. The MSM reduction is real; it is simply diluted.

## 3. KZH3 (eprint 2025/144, Appendix C)

`kzh/` implemented KZH2 (the two-dimensional scheme of Figure 2). The paper's Appendix C generalizes
it to KZH-k; `kzh/kzh3.hpp` adds the k = 3 case: a three-dimensional tensor split with slice
commitments `D1`, contracted row commitments `D2` and a final `d₃`-length vector, bound by a
two-level pairing chain. The Honk shift extends the KZH2 backend's verifier-side query decomposition
to rank 3 — the shift of the flattened index crosses both the `d₃` and `d₂` boundaries, so chain B
carries two extra final-layer vectors and one extra slice-commitment layer.

Proof size drops from `O(√N)` to `O(∛N)` and verification with it, and both gains widen with size:
at 2^18 the proof is 55.3 KiB against KZH2's 183.3 (3.3x) and verification 79.7 ms against 216
(2.7x); at 2^20, 87.8 KiB against 359.8 (4.1x) and 124 ms against 407 (3.3x).

The prover cost also widens, and by more than the small-size measurement suggests: 1.6x at 2^14
(3.1 s against 1.9 s) but **3.2x at 2^18** (25.2 s against 7.8 s) and 4.8x at 2^20. The extra
dimension adds a whole layer of slice commitments whose count grows with `d1`, so the crossover
where KZH3's smaller proof is worth its prover sits well below 2^18. KZH3 is the right choice when
proof size or decider cost dominates — its stated purpose in the paper, where the k-dimensional
generalization exists to shrink the accumulator and decider — and the wrong one when raw proving
throughput does.

## 4. WHIR soundness audit against Crites–Stewart (eprint 2025/2046)

Crites and Stewart **disprove** the up-to-capacity correlated-agreement, mutual-correlated-agreement
and list-decodability conjectures. WHIR specifically relies on the mutual correlated agreement
adaptation (their Conjecture 4.12, from [BCI+23] Conjecture 8.4), which is among the disproved set.
The counterexamples live between the list-decoding capacity bound and the code's capacity: they
construct `u⁽⁰⁾, u⁽¹⁾` with `u⁽¹⁾` farther than `δ` from the code such that `u⁽⁰⁾ + λu⁽¹⁾` is
`δ`-close for *every* `λ`, whenever `δ > 1 - H_q(ρ)`. Their minimally repaired conjectures replace
`δ ≤ 1-ρ-η` by `H_q(δ) ≤ 1-ρ-η`.

`whir/whir_config.hpp`'s default `CONJECTURED_LIST` regime (`t = ⌈λ/r⌉`) rests squarely on the
disproved statement. The audit added a `REPAIRED_LIST` preset implementing the repaired bound: the
largest testable distance `δ*` with `H_q(δ*) = 1-ρ`, solved by a deterministic Q192/Q64 fixed-point
iteration (integer-only `-log₂` by repeated squaring) with every rounding pushed toward more
queries.

**The delta is small.** At BN254's 254-bit field the entropy penalty `h₂(δ)/log₂q` costs exactly one
extra query per aggressive-rate round — 51 and 21 queries at rates `r = 2, 5` against 50 and 20, and
nothing at all from `r = 8` upward. End to end that is **+2.0% proof size at every size from 2^12 to
2^20** (1097.0 → 1119.2 KiB at 2^12, 1266.8 → 1292.0 at 2^20) with verify time unchanged inside
run-to-run noise. The default preset stays `CONJECTURED_LIST` for comparability with deployed
FRI/STIR systems; `whir_honk_repaired_*` benchmarks both, and both are documented in
`whir/README.md` §6.

The practical conclusion: WHIR's parameters were built on a conjecture that is now false, but at a
254-bit prime field the repair is nearly free. Small-field WHIR deployments would pay much more,
since the penalty scales as `1/log₂q`.

## 5. SwitchFold (eprint 2026/1489)

SwitchFold keeps Ligero's interleaved commitment and replaces the clear-text combined rows — the
`O(√N)` term that dominates Ligero's proof — with a recursive code-switching descent. The three
combined rows of the rank-2 shift decomposition become the four segments of one descent message, so
every claim about it is a tensor. Each level commits the message reshaped under the next shorter RS
code, opens it at that level's queries, collapses the pending claims to a single tensor claim with
one inner-product sumcheck (paper Fact 1), and splits that claim into the column combination that
becomes the next message. This level's queries become the next level's code-switching claims, until
a 16-element base message is sent in the clear.

**Instantiating over Reed–Solomon collapses two of the paper's three modules.** The RS
generator-matrix row at a domain point *is* the pow tensor, so the code-switching claim
`Ξ·Enc[C](m) = Ξ·G·m` needs neither a commitment to the generator matrix nor the accumulation
scheme that Brakedown's sparse generator forces. What the paper spends its second and third modules
on is free here. Since nothing of row length is transmitted, the proof-optimal shape also moves:
`num_cols` goes as large as BN254's 2-adicity allows, which shrinks the payload leaves.

Measured against LigeroHonk with an identical commitment, SwitchFold reproduces the paper's headline
**3.5x proof-size reduction** at 2^18 (772.9 KiB against 2692.6) and does better as the trace grows
— **5.9x at 2^20** (878.0 against 5144.8), because Ligero's transmitted rows grow as `√N` while the
descent only adds levels. Verification improves 1.5x (13.6 ms against 20.6), and SwitchFold's proof
is also 1.6x smaller than WHIR's at comparable verify time, making it the best hash-based proof size
in the suite.

The prover costs 1619 ms against Ligero's 382, for two reasons worth separating. The paper's linear
prover comes from Brakedown's linear-time encoder, which bb does not have (see §7) — over RS the
descent pays FFT encoding at every level. And the batched claim carries one term per query per
segment, so building the sumcheck weight table is multiplication-heavy, which this no-assembly
arm64 build penalizes relative to Ligero's FFT-and-hash profile; parallelizing that build over
terms with private accumulators already took the prover down 2.6x from where it started.

## 6. Titan (eprint 2026/908) — designed, not implemented

Titan is a transparent DL-only PCS: Pedersen row commitments feed a group polynomial whose
commitment is an IOPP over groups (WHIR/BaseFold adapted to prime-order groups) rather than Dory's
pairing-based two-tier commitment, with a compressed sigma protocol and delegated generator folding
for the inner layer. It is the strongest remaining unimplemented candidate and the design work is
done; two findings are worth recording.

**The Honk shift cannot be derived from the base oracle.** The natural shortcut is to derive the
shifted group polynomial's codeword as `(c_B[s] - B₀)·ω⁻ˢ` from the committed one, so no second
commitment is needed. This is **unsound**: a corrupted `B₀` perturbs the derived word by
`Δ·x⁻¹`, and `x⁻¹` restricted to the domain is a codeword of degree `|L|-1`, so the corrupted word
sits at Hamming distance as low as 1 from a valid codeword and the query test cannot see it. The
correct construction commits the shifted group polynomial as an extra codeword column (the
`to_be_shifted` flag that `commit_group` already carries) and enforces
`ĉ_B(X) = B₀ + X·ĉ_D(X)` pointwise at the query positions; that residual has degree ≤ p, so it fails
on a `1-ρ` fraction of the domain and `t` queries catch it with probability `1-ρᵗ`.

**The cost is the real obstacle.** Titan's commitment is Pedersen rows (`N` scalar multiplications
per committed column) *plus* the outer RS encoding of the group polynomial (`N/ρ` more). Hyrax pays
only the first term and already measures ~6.0 s at 2^18 in this suite; KZH2 measures 7.8 s and Dory
6.6 s. Titan therefore lands in the 20–30 s range at 2^18 on this build — an order of magnitude
past the transparent hash backends, and past the group backends it is meant to dominate. The
encoding can be restructured so the group work is `|L₁|` MSMs of size `q` (encode the `q` matrix
*columns* with cheap field FFTs first, then one MSM per domain point), which keeps it at the
paper's `O(n/ρ)` rather than worse, but does not change the order of magnitude.

The paper's own comparison is against Dory and Hyrax on Pasta curves with a Rust implementation; the
claim that Titan dominates them is plausible and orthogonal to the constant-factor problem above.
Implementing it remains worthwhile, but as a scoped piece of work rather than a variation on an
existing module.

## 7. Lightning (2026/258) and Bolt (2026/310) — viability at BN254 Fr

The hypothesis under test was that these prover-optimal code-based schemes are tuned to small or
binary fields and that their advantage — in particular Bolt's replacement of multiplications by
additions — would not survive at a 254-bit prime field. **The measurements do not support that
hypothesis for either scheme.**

**Field-operation cost ratio.** A microbenchmark over 2^22 random BN254 Fr elements on this build
(Apple M4 Pro, arm64, `DISABLE_ASM=1`) gives **4.09 ns per addition and 17.13 ns per multiplication,
a ratio of 4.19x**. This is the number that decides Bolt's trade. Over GF(2^32), where Bolt is
benchmarked, an addition is an XOR and the ratio is one to two orders of magnitude — so the trade
is far less dramatic at BN254, exactly as suspected. But it does not invert: Bolt's stated
commitment cost of `(3+ε)·N` additions is ~0.7N multiplication-equivalents, while RS encoding of
the same data costs `Θ(N log N)` multiplications. At 2^18 with Ligero's shape that is roughly 8.4M
multiplications, or ~35M addition-equivalents, against Bolt's ~0.8M additions — a 40x reduction in
field work even at BN254.

**Is encoding actually the bottleneck?** Yes, on this build. Comparing LigeroHonk under the two
hashers at 2^18 (Blake3 389 ms, Poseidon2 2938 ms) and attributing the 2549 ms difference to
hashing, the Blake3 configuration spends roughly 90% of its prover outside hashing — i.e. in
encoding and row combination. A 40x cheaper encoder would therefore translate into a large real
speedup rather than being absorbed by Merkle hashing.

**So what blocks them?** Not the field, and not the PCS framework — bb now has the code-switching
framework, via the SwitchFold module of §5. What is missing is the **code family**: Lightning needs
a base code with constant relative distance (Spielman/expander style, as in Brakedown) to which it
applies its sparsity-preserving compression `C_L(m) = m ‖ C_D(mA)`; Bolt needs its sketched random
LDPC codes, whose large-field variant the paper treats in its §6.2. Neither exists in bb, and
implementing one means bringing in a new code family together with certified distance parameters —
work that is substantial, independent of the PCS layer, and shared between the two schemes.

**The prerequisite now exists.** `brakedown/` implements the Spielman-style code of
[2021/1043](https://eprint.iacr.org/2021/1043) — see that module's README for the construction, the
Figure 2 parameter validation, and the measurements. Two results from it change the picture above:

- *The encoder is genuinely faster, but only modestly.* Against `rs_encode` at rate 1/4 it wins
  1.13x at 2^12 widening to 1.22x at 2^16, at 26.9 multiplications per symbol (the paper predicts
  25.5n), emitting 2.3x fewer symbols. The `Θ(n)` versus `Θ(n log n)` gap is real and widening, but
  bb's RS encoder is a well-optimized parallel FFT, so the constant factors nearly cancel at the
  sizes this suite covers. The 40x field-work reduction estimated above is an *operation count*,
  not a wall-clock speedup.
- *Distance, not encoding speed, is the binding constraint.* Brakedown's relative distance is 0.07
  against RS's 0.75 at rate 1/4, so the provable interleaved proximity test needs about 2936
  queries at λ = 100 against RS's 50 — a 59x increase in openings. This is why Brakedown's own
  paper opens 6593 columns and reports proofs in the tens of megabytes, and it is why the module
  ships as a reusable code rather than as another `*Honk` backend: dropping it into Ligero or
  SwitchFold would trade a ~20% encode win for a ~59x query-side proof blowup.

**Bolt is now implemented too** (`bolt/`), and it settles the question the original hypothesis posed.

Bolt's code is `C_H(x) = (x, C(Hx))` — a sparse LDPC sketch fed to a base code — and its real
contribution is a *piecewise* distance guarantee: distinct codewords differ either in more than `γ`
of the systematic stretch or more than `δ` of the sketch stretch. `BoltHonk` tests the two stretches
independently rather than sampling uniformly against the diluted average, which is what
`ligero/QuerySegment` exists for.

**The additions-vs-multiplications trade does not survive at BN254 — but not for the reason
expected.** `γ` is the root of the paper's `ω_{q,j,k}`; the implementation reproduces the paper's own
reference point (`q = 2^32, j = 16, k = 128 → 0.094114`). The identical parameters at BN254 give
`γ = 0.00006`. The cause is structural: `ω`'s leading term is `x·ln(q−1)`, so a fixed column degree
certifies less distance as the field grows and the degree must scale with `ln q ≈ 176`. Recovering
`γ = 0.15` needs `j = 64`, i.e. **64 multiplications per message symbol**. Bolt's advertised
`(3+ε)N` field *additions* is the Boolean-`H` instantiation of its §6.1, valid over bits; it does not
carry to a 254-bit prime. So the earlier conclusion — that the operation mix left Bolt viable — was
right about the arithmetic and wrong about the parameters.

**What survives is the distance, and that is the half worth having.** At 2^14, the same tensor
protocol over three codes:

| | code | distance | queries | prove | verify | proof |
|---|---|---|---|---|---|---|
| Ligero | RS, rate 1/4 | 0.75 (conj.) | 50 | **41.2 ms** | **9.4 ms** | **833 KiB** |
| Bolt | sketched, piecewise | 0.15 / 0.75 | 428 + 50 | 86.3 ms | 96.3 ms | 3.94 MiB |
| Brakedown | Spielman | 0.07 | 2936 | 78.9 ms | 156 ms | 20.2 MiB |

Bolt's proof is **5.1x smaller than Brakedown's** with a 1.6x faster verifier — the piecewise test
doing exactly what it should — but still 4.7x Ligero's. The ordering on the binding axis is
`Ligero < Bolt < Brakedown`.

**Conclusion for the whole linear-time-code line at BN254.** It is not yet competitive with a tuned
Reed–Solomon encoder, and the obstacle is the certified distance these codes can offer at a 254-bit
field — not their encoding speed, which was the original hypothesis and which the measurements
retire. Lightning would not change this: it explicitly trades 2.4x proof size for its 2.7x prover
win, the wrong direction when proof size binds. The levers that would change the answer are a
tighter distance analysis (RS gets its 0.75 from a *conjecture*; these codes are held to
union-bound-style proofs), a smaller sketch ratio at fixed `γ`, or composing a sketched code with
the `switchfold/` descent so the transmitted-row and query terms are attacked together. See
`bolt/README.md` §5.

## 8. Carina (eprint 2026/1438) — not implemented

Carina is Vela's sibling and shares §1's paper. It applies the constant-term reduction once along
each coordinate of a `√N × √N` tensor embedding and proves the committed bivariate together with
two auxiliary polynomials with a single grid opening, reaching `N + 2√N - 6` opening MSM scalars —
close to the theoretical floor — with a 4 G1 + 8 F proof.

It was not implemented for a concrete reason: its final check is
`e(C_E, [1]₂) = e(Π_X, [Z_A(τ)]₂)·e(Π_Y, [Z_B(σ)]₂)` with both `Z_A` and `Z_B` of degree two,
requiring five fixed G2 elements — `[1]₂, [τ]₂, [τ²]₂, [σ]₂, [σ²]₂`. bb publishes two. The
`σ = 2` trick used for CHOPIN's second trapdoor supplies `[σ]₂` and `[σ²]₂` for free, but `[τ²]₂`
is unavailable for the real ceremony trapdoor and would need the same linearization Vela uses,
applied to the X coordinate. That is tractable but it is a second deviation stacked on a
substantially more complex protocol, and Vela already occupies the corner of the design space
Carina's sibling targets. CHOPIN (§2) covers the prover-oriented corner Carina competes in.

## Previously evaluated and rejected: Binius, Flock

Binius / FRI-Binius ([2024/504](https://eprint.iacr.org/2024/504)) and Flock
([2026/1329](https://eprint.iacr.org/2026/1329)) are proof systems over **binary tower fields**
whose performance comes from the arithmetization, not from a detachable commitment scheme. They
cannot fill UltraHonk's PCS slot:

- **Field mismatch is structural.** UltraHonk's polynomials have BN254 Fr coefficients and its
  sumcheck emits Fr evaluation claims. A binary-field PCS commits to multilinears over
  GF(2^k) towers; no homomorphism connects the two. The only bridge — bit-decomposing each Fr
  coefficient and proving the Fr-weighted recombination (with modular reduction) inside the
  binary system — turns the opening argument into large-prime-field emulation over GF(2), the
  workload binary systems are worst at.
- **Their advantage is Boolean-native witnesses.** UltraHonk traces are Fr-native (Poseidon,
  Pedersen, EC, field arithmetic); expressed over GF(2), each Fr multiplication costs thousands of
  Boolean gates, inverting the comparison.

## Remaining candidates, prioritized

1. **A sketched code under the code-switching descent.** `switchfold/` removes Ligero's
   transmitted-row term but does nothing for the query term, which is exactly what dominates Bolt and
   Brakedown (§7). Running the descent over `bolt/`'s code attacks both at once, and all three pieces
   now exist, so this is the cheapest remaining experiment with a real chance of moving the
   transparent frontier.
2. **Titan** ([2026/908](https://eprint.iacr.org/2026/908)) — designed (§6); needs the group-IOPP
   and compressed-sigma machinery plus a resolution of the constant-factor problem.
3. **zip proof compression** ([2025/1446](https://eprint.iacr.org/2025/1446)) — not a PCS: a
   compression layer for hash-based proofs. It attacks the MiB-scale proofs of WHIR, Ligero,
   SwitchFold, Bolt and Brakedown simultaneously, which given §7's conclusion is arguably better
   leverage than any further code.
4. **Orion** ([2022/1010](https://eprint.iacr.org/2022/1010)) — the one construction with both a
   linear prover and polylog proof and verify. Its expander testing targets precisely the
   small-distance problem §7 identifies, and its proof composition the query term. Needs a careful
   read of the current version first: the original expander-testing soundness argument had a gap that
   was later patched.
5. **Carina** ([2026/1438](https://eprint.iacr.org/2026/1438)) — §8; needs the X-coordinate
   linearization.
6. **Samaritan** ([2025/419](https://eprint.iacr.org/2025/419)) — intra-corner A/B against Mercury
   and Vela; low marginal insight now that Vela and CHOPIN occupy that corner.

Tracked but not planned: class-group constant-size transparent PCS
([2025/1233](https://eprint.iacr.org/2025/1233), no GUO arithmetic in bb, seconds-scale provers);
Zeromorph/PST/HyperKZG (the Shplemini corner, already represented); BaseFold/STIR (dominated by
WHIR; useful only as ablation baselines); Greyhound/lattice (wrong ring);
FRIttata ([2025/1285](https://eprint.iacr.org/2025/1285), distributed proving — an orthogonal axis).
