# Bolt — sketched codes, tested piecewise

Implementation of Bolt's sketched code ([eprint 2026/310](https://eprint.iacr.org/2026/310) §3, §6.2)
and `BoltHonk`, Ligero's tensor PCS over it. The interesting part is not the encoder — it is the
*piecewise* distance guarantee, which is what actually buys proof size.

## 1. The code

```
C_H(x) = (x, C(Hx))
```

`H` is a sparse parity-check matrix from the random `(j, k, q)`-LDPC ensemble — `n` variable nodes of
degree `j`, `nj/k` check nodes of degree `k`, connected by a uniformly random socket bijection with
an independent random nonzero field weight per check socket. So `H` is `αn × n` with `α = j/k`, and
the expensive base code `C` (Reed–Solomon here) is applied only to the short sketch `Hx`. The result
is systematic: positions `[0, n)` are the message, `[n, n + αn/ρ)` the encoded sketch.

## 2. Piecewise distance is the whole point

Claim 3.1 of the paper: for two distinct codewords, **either** they differ in more than `γ` of the
systematic stretch **or** in more than `δ` of the sketch stretch — where `γ` is the LDPC's certified
distance and `δ` the base code's. That is strictly more information than the code's overall distance

```
min(γ, δ·α/ρ) / (1 + α/ρ)
```

which averages the two and is dominated by the worse one. So `BoltHonk` queries the two stretches
*independently*, each with the count its own distance warrants:

| stretch | distance | queries at λ=100 |
|---|---|---|
| systematic `[0, n)` | `γ = 0.15` | `⌈100 / −log₂(1−γ)⌉` = 428 |
| sketch `[n, …)` | RS at rate 1/4 | `⌈100/2⌉` = 50 |

A cheating prover must survive both, so the error is the max of the two. Sampling uniformly against
the diluted distance instead costs several times more; `PiecewiseQueriesBeatTheDilutedAlternative`
asserts the gap. `ligero/QuerySegment` is the hook that makes this expressible — Reed–Solomon and
Brakedown return one segment, Bolt returns two.

## 3. The field size changes the parameters, and it is the main finding

`γ` is the root of the paper's `ω_{q,j,k}` (Theorem 6.5, from Yang et al.), implemented in
`ldpc_distance_bound`. `BoltParamsTest.DistanceBoundMatchesPaperInstantiation` reproduces the
paper's own reference point — `q = 2^32, j = 16, k = 128 → x₀ ≈ 0.094114` — which is the check that
the transcription is right.

Run the same formula at BN254 and it collapses:

| q | j | k | α | γ |
|---|---|---|---|---|
| 2^32 | 16 | 128 | 0.125 | **0.0941** |
| BN254 | 16 | 128 | 0.125 | **0.00006** |
| BN254 | 24 | 128 | 0.1875 | 0.00375 |
| BN254 | **64** | **256** | 0.25 | **0.1497** |
| BN254 | 128 | 512 | 0.25 | 0.2419 |

The cause is structural: `ω`'s leading term is `x·ln(q−1)` from `H_q`, so a fixed column degree
certifies less and less distance as the field grows — the degree has to scale with `ln q ≈ 176`.

**This is why Bolt's headline cost does not carry to BN254.** The advertised `(3+ε)N` field
*additions* is the Boolean-`H` instantiation of §6.1, where `j = 3` and the 0/1 entries make the
sketch addition-only. Over a 254-bit prime the same theorem needs random field weights at `j = 64`,
i.e. **64 multiplications per message symbol** — against Brakedown's 27 for a whole encoding, and
against 3 additions. The encoder advantage is gone.

What survives is the distance: `γ = 0.15` against Brakedown's `0.07`. Since distance is the binding
constraint on proof size (`../brakedown/README.md` §5), that is the half worth having.

## 4. Measured

`LigeroHonk`, `BrakedownHonk` and `BoltHonk` are the same tensor protocol over three different codes,
so they are a clean three-way A/B. At 2^14 with the Blake3s hasher:

| | code | distance | queries | prove | verify | proof |
|---|---|---|---|---|---|---|
| Ligero | RS, rate 1/4 | 0.75 (conj.) | 50 | **41.2 ms** | **9.4 ms** | **833 KiB** |
| Bolt | sketched, piecewise | 0.15 / 0.75 | 428 + 50 | 86.3 ms | 96.3 ms | 3.94 MiB |
| Brakedown | Spielman | 0.07 | 2936 | 78.9 ms | 156 ms | 20.2 MiB |

**Bolt's proof is 5.1x smaller than Brakedown's** and its verifier 1.6x faster, which is the
piecewise test doing exactly what it is supposed to. It does not overtake Reed–Solomon: still 4.7x
Ligero's proof and 10x its verification, because 428 systematic queries is a lot next to 50, and
because a 64-multiplication sketch is not cheap.

So the ordering on the axis that matters is `Ligero < Bolt < Brakedown`. Bolt closes most of the gap
Brakedown opened without closing all of it. The honest summary: at BN254 the linear-time-code line
is not yet competitive with a well-tuned RS encoder, and the reason is the certified distance these
codes can offer at a 254-bit field — not their encoding speed, which was the original hypothesis.

## 5. What would change the answer

- **A better distance bound.** The `ω` bound is a union-bound-style argument; a tighter analysis (or
  a conjecture, as RS enjoys via capacity) would cut the 428 systematic queries directly.
- **A base code with a shorter sketch.** `α = 0.25` at `j = 64` means the sketch is a quarter of the
  message. Driving `α` down while holding `γ` shrinks the encoded tail.
- **Composing with the code-switching descent.** `switchfold/` removes the *transmitted-row* term
  from Ligero's proof; it does nothing for the query term, which is what dominates here. Running the
  descent over a sketched code would attack both at once, and is the natural next experiment.
