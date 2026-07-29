# cq: lookups with prover cost independent of table size

UltraHonk's lookup argument materializes every table in the execution trace: the `table_1..table_4` polynomials
concatenate all used basic tables, the trace must be at least as long as the total number of table rows, and the
log-derivative relation does work proportional to the full (dyadic) circuit size. A table of $2^{20}$ entries
therefore costs the prover $2^{20}$ trace rows on every proof, whether the circuit performs two lookups or two
thousand. This module implements the cq lookup argument [1], in which a table of size $N$ is preprocessed once for
$O(N \log N)$ group operations, after which each proof costs the prover $O(n \log n)$ field operations and $O(n)$
group operations for $n$ lookups — with no dependence on $N$ anywhere in proving or verification.

cq is the state of the art in the sequence of preprocessing-based lookup arguments Caulk [2] → Caulk+ [3] →
Flookup [4] → Baloo [5] → cq: constant-size proofs, pairing-based verification compatible with the BN254 KZG stack
used by UltraHonk, and a prover that only performs sparse multi-scalar multiplications over cached commitments.
The alternative family (Lasso [8]) replaces preprocessing with sumcheck over decomposable tables; it wins for
structured tables too large to materialize at all (e.g. $2^{32}$ two-operand tables) but requires different
machinery and does not produce constant-size pairing-checkable proofs. A survey of the design space is [9].

## Claims

For a table $t \in \mathbb{F}^N$ and witness $f \in \mathbb{F}^n$ (both sizes powers of two, $n \le N$), the
protocol proves $\{f_j\} \subseteq \{t_i\}$ as multisets of values. Following the log-derivative formulation [6],
this is equivalent (with overwhelming probability over a challenge $\beta$) to the existence of multiplicities
$m_i$ with

$$\sum_{i=0}^{N-1} \frac{m_i}{\beta + t_i} \;=\; \sum_{j=0}^{n-1} \frac{1}{\beta + f_j}. \tag{S}$$

The two sides are established by three labelled identities over the subgroups $H = \langle\omega\rangle$ of order
$N$ and $V = \langle\nu\rangle$ of order $n$, with $Z_H(X) = X^N - 1$ and $Z_V(X) = X^n - 1$:

$$A(X)\,(T(X) + \beta) - m(X) \;=\; Q_A(X)\, Z_H(X) \tag{C1}$$

$$P_B(X) \;=\; B_0(X)\cdot X^{N-n+1} \tag{C2}$$

$$B(X)\,(f(X) + \beta) - 1 \;=\; Q_B(X)\, Z_V(X) \tag{C3}$$

where $A$ interpolates $m_i/(\beta + t_i)$ over $H$, $B$ interpolates $1/(\beta + f_j)$ over $V$, and
$B_0(X) = (B(X) - B(0))/X$. Because the sum of a polynomial's values over a subgroup is the subgroup order times
its constant coefficient (given degree less than the order), (S) reduces to the linkage $N \cdot A(0) = n \cdot
B(0)$: the verifier learns $A(0)$ from a KZG opening at zero and *derives* $B(0) := N \cdot A(0) / n$, so the
linkage holds by construction and (C3) at a random point pins $B$ to its intended values.

## Notation

| Symbol | Meaning | Code |
|---|---|---|
| $N$ | table size (power of 2) | `CqProvingKey::table_size` |
| $n$ | number of lookups (power of 2, $2 \le n \le N$) | `num_lookups` |
| $K$ | number of table columns | `CqProvingKey::num_columns()` |
| $H, V$ | subgroups of order $N$, $n$ | `SubgroupDomain` |
| $t^k, f^k$ | table / witness column $k$ | `columns`, `lookup_columns` |
| $T^k$ | polynomial interpolating $t^k$ over $H$ | `column_monomials[k]` |
| $\theta$ | column-combination challenge | `"CQ:theta"` |
| $\beta$ | log-derivative challenge | `"CQ:beta"` |
| $m$ | multiplicity polynomial over $H$ (sparse) | `m_values` on `support_indices` |
| $A, Q_A$ | table-side rational values and quotient | `a_values`, `q_a_commitment` |
| $B, B_0, Q_B$ | witness-side values, shifted poly, quotient | `b_monomial`, `b_shifted`, `q_b` |
| $P_B$ | degree-shift $B_0 \cdot X^{N-n+1}$ | `degree_shift_commitment` |
| $\gamma, \eta, \rho$ | evaluation / opening-batch / check-batch challenges | `"CQ:gamma"`, `"CQ:eta"`, `"CQ:rho"` |
| $\tau$ | SRS trapdoor | `TestSrs` (tests only) |

## Preprocessing (`CqProvingKey::create`)

Once per table, for each column $k$, and shared across all future proofs:

1. $T^k \leftarrow \mathrm{IFFT}_H(t^k)$, with commitments $[T^k]_1$ and $[T^k]_2$ (the G2 commitment is what lets
   the verifier pair against a committed table).
2. All $N$ cached quotients
   $$[Q^k_i]_1 = \left[\frac{L_i(X)(T^k(X) - t^k_i)}{Z_H(X)}\right]_1 = \frac{\omega^i}{N}\left[\frac{T^k(X) - T^k(\omega^i)}{X - \omega^i}\right]_1,$$
   i.e. scaled KZG opening proofs of $T^k$ at every point of $H$, computed in $O(N \log N)$ group operations with
   the Feist–Khovratovich algorithm [7] (`fk_all_kzg_opening_proofs`): a Toeplitz matrix–vector product with the
   SRS realized as a size-$2N$ cyclic convolution (`toeplitz_srs_products`), followed by a size-$N$ group FFT.
   The group transforms are the field FFT templated over Jacobian points (`SubgroupDomain::fft<g1::element>`).
3. Column-independent Lagrange caches: $[L_i]_1$ for committing to sparse Lagrange-basis polynomials
   (`lagrange_basis_commitments` — the group IFFT of the SRS), and $[(L_i(X) - L_i(0))/X]_1$ for KZG openings at
   zero of such polynomials (`lagrange_zero_quotient_commitments`).
4. A row-index map for multiplicity counting (`find_row`).

The verification key (`CqVerificationKey`) stores $[T^k]_1, [T^k]_2$, $[Z_H]_2 = [\tau^N - 1]_2$ and the G2 powers
needed for the degree check.

## Prover (`CqProver::construct_proof`)

Transcript schedule (`NativeTranscript`), with literal labels; every challenge is drawn after all messages it must
bind:

| Round | Sends | Then derives |
|---|---|---|
| 0 | `CQ:table_size`, `CQ:num_lookups`, `CQ:num_columns`, `CQ:T_k` (each $[T^k]_1$) | — |
| 1 | `CQ:m` $= [m]_1$, `CQ:f_k` $= [f^k]_1$ | `CQ:theta`, then `CQ:beta` |
| 2 | `CQ:A`, `CQ:Q_A`, `CQ:B_0`, `CQ:Q_B`, `CQ:P_B` | `CQ:gamma` |
| 3 | `CQ:b0_gamma`, `CQ:f_gamma`, `CQ:a0` | `CQ:eta` |
| 3 | `CQ:pi_gamma`, `CQ:pi_zero` | (`CQ:rho`, verifier side) |

Round 1. Multiplicities are counted on witness *rows*; $[m]_1 = \sum_{i \in S} m_i [L_i]_1$ is a sparse MSM over
the cached Lagrange commitments, where $S$ (the support) has at most $n$ elements. Each witness column is
interpolated over $V$ and committed. After $\theta$, rows are reduced to scalars: $t_i(\theta) = \sum_k \theta^k
t^k_i$ (computed only on $S$) and $f_j(\theta) = \sum_k \theta^k f^k_j$.

Round 2, table side. $A_i = m_i / (\beta + t_i(\theta))$ on $S$ via one batched inversion;
$[A]_1 = \sum_{i \in S} A_i [L_i]_1$ and, using the homomorphic combination of per-column caches,
$$[Q_A]_1 = \sum_k \theta^k \sum_{i \in S} A_i\, [Q^k_i]_1,$$
$K$ sparse MSMs of size $|S|$. This is the cached-quotients step that removes the table size from the prover: an
uncached prover would need an $O(N \log N)$ polynomial division per proof to produce $Q_A$.

Round 2, witness side. $B_j = 1/(\beta + f_j(\theta))$ via batched inversion, $B \leftarrow \mathrm{IFFT}_V$,
$Q_B = (B \cdot (f + \beta) - 1)/Z_V$ via a size-$2n$ product (`polynomial_product`) and exact division by
$X^n - 1$ (`divide_by_vanishing`). $[P_B]_1$ is committed with the SRS span shifted by $N - n + 1$
(`PolynomialSpan{ big_n - n + 1, b_shifted }`), which is possible precisely when $\deg B_0 \le n - 2$.

Round 3. Evaluations $B_0(\gamma)$, $f(\gamma)$ and $a_0 = A(0)$ are sent. The batched KZG witness for
$(B_0, f, Q_B)$ at $\gamma$ is one quotient (`factor_roots`), and the opening of $A$ at zero is the sparse MSM
$\pi_0 = \sum_{i \in S} A_i [(L_i - L_i(0))/X]_1$ over the cached zero-quotients.

## Verifier (`CqVerifier::verify_proof`)

After recomputing challenges and homomorphically combining $[f]_1 = \sum_k \theta^k [f^k]_1$ and
$[T(\theta)]_2 = \sum_k \theta^k [T^k]_2$, three pairing products (7 pairings total):

1. (C1) $\;e([A], [T(\theta)]_2)\; e(\beta [A] - [m], [1]_2)\; e(-[Q_A], [Z_H]_2) = 1$.
2. (C2) $\;e([B_0], [\tau^{N-n+1}]_2)\; e(-[P_B], [1]_2) = 1$; with the SRS truncated at degree $N-1$ this bounds
   $\deg B_0 \le n-2$, hence $\deg B \le n-1$ as required by the sum linkage.
3. (C3 + opening of $A$ at 0, batched by $\rho$): with $b_0 := N a_0 / n$, $B(\gamma) := \gamma B_0(\gamma) + b_0$
   and $Q_B(\gamma)$ *derived* from (C3), verify the two KZG openings as
   $$e\big(C - v[1]_1 + \gamma \pi_\gamma + \rho([A] - a_0 [1]_1),\, [1]_2\big)\; e\big(-(\pi_\gamma + \rho\, \pi_0),\, [\tau]_2\big) = 1,$$
   where $C = [B_0] + \eta [f] + \eta^2 [Q_B]$ and $v = B_0(\gamma) + \eta f(\gamma) + \eta^2 Q_B(\gamma)$.

## Costs

| | field ops | group ops | pairings | data |
|---|---|---|---|---|
| Preprocess (per column) | $O(N \log N)$ | $O(N \log N)$ (3 group FFTs, 2 MSMs) | — | stores $(2{+}K) N$ G1 + $N{+}1$ G2 |
| Prove | $O(K n \log n)$ | $K{+}3$ sparse MSMs $\le n$, $K{+}4$ dense MSMs $\le n$ | — | proof: $2K{+}8$ G1, 3 $\mathbb{F}$ |
| Verify | $O(\log n)$ | $O(K)$ | 7 | — |

Nothing on the prove or verify rows depends on $N$. Measured with `cq.bench.cpp` on an M4 Pro (arm64,
`DISABLE_ASM` build, so absolute times are conservative), single column:

| $N$ | preprocess (once) | prove, $n = 256$ | prove, $n = 1024$ |
|---|---|---|---|
| $2^{10}$ | 0.18 s | 14.3 ms | 20.8 ms |
| $2^{12}$ | 0.83 s | 14.7 ms | 24.0 ms |
| $2^{14}$ | 3.96 s | 15.9 ms | 20.7 ms |
| $2^{16}$ | 18.2 s | 14.2 ms | 21.2 ms |

Verification: 2.55 ms. For scale: an UltraHonk circuit forced to $2^{16}$ rows just to hold such a table takes
seconds to prove on the same machine, on every proof.

## Soundness notes

- $[m]$, $[f^k]$ and the (preprocessed) table commitments enter the transcript before $\beta$ and $\theta$, as
  required for the Schwartz–Zippel arguments over (S) and over the row combination.
- $m$ needs no degree bound or integrality proof: (C1) only pins the values of $A$ against *some* committed $m$,
  and (S) for a random $\beta$ already forces every pole $1/(\beta + f_j)$ to be matched by a table pole.
- $\deg A < N$ is enforced by the SRS itself: the G1 SRS for a size-$N$ table must contain exactly the powers
  $[\tau^0]_1..[\tau^{N-1}]_1$ (cq specifies this truncation; `TestSrs` reproduces it). $\deg B_0$ gets the
  explicit check (C2) because its bound $n-2$ is stricter than the SRS bound.
- $B(0)$ is derived, never sent, so the sum linkage cannot be misstated; the verifier rejects if $\gamma \in V$.
- Extraction is in the algebraic group model, as for cq [1] and KZG generally.
- This is the non-zero-knowledge variant; cq's zk variant (blinding $m$, $A$, $B$ and the quotients) is future
  work, as is batching the three pairing products into one.

## Trusted setup

cq requires G2 powers of $\tau$ up to degree $N$. Public powers-of-tau ceremonies publish these; the SRS files
shipped with barretenberg contain only $[\tau]_2$, so `TestSrs` synthesizes a fresh $\tau$ in-process — test-only
by construction, since the trapdoor is known. A deployment would load G1 and G2 powers from a ceremony transcript
and install the G1 prefix as the global CRS.

## Relation to UltraHonk tables

`cq_plookup.test.cpp` proves lookups against unmodified UltraHonk `plookup::BasicTable` data (the 6-bit XOR table
and the SHA-256 majority normalization table), using the multi-column reduction for their 3-column rows. The
remaining engineering to use cq inside an UltraHonk proof is the witness linkage: UltraHonk commits wire values in
evaluation basis, while cq's $[f^k]$ commit the interpolating polynomial over $V$. Two options:

1. Commit $f^k$ as here and add one consistency opening to the existing Shplemini batch: the multilinear
   evaluation of the wire at the sumcheck point and $f^k(\gamma)$ both open commitments to the same value vector.
2. Run cq's witness side in Lagrange basis over $V$ with a preprocessed Lagrange SRS for the wire domain, making
   $[f^k]$ *equal* to the wire commitment.

Either way the verifier gains 7 pairings (or, batched, 2 extra pairing-point contributions that fold into the
existing `PairingPoints` aggregation), and the recursive verifier needs the G2 points as circuit constants.

## Applications

The tables UltraHonk uses today are small precisely because they must fit in the trace. With cost moved to
preprocessing, table size is bounded only by preprocessing time and storage ($\approx (2{+}K) N$ G1 points, i.e.
$\approx 5$ GB at $N = 2^{24}$ with $K = 3$), so $N = 2^{20}..2^{24}$ becomes realistic. The advantage over
current gate counts comes from three distinct mechanisms:

1. **Trace decoupling** — table rows leave the circuit. Today a circuit pays the summed size of every table it
   touches, then rounds up to a power of two.
2. **Width scaling** — the number of lookups per operation shrinks roughly linearly in the chunk width while the
   table size grows exponentially; cq moves the feasible frontier from $2^{12}$ to $\approx 2^{24}$.
3. **Semantic fusion** — a table row encodes an arbitrary function of its key, so a whole gadget (a hash-round
   subfunction, a VM instruction, an activation function, a Merkle membership check) collapses to one lookup.

Baseline gate counts for this branch, measured by comparing finalized `UltraCircuitBuilder` sizes of a circuit
containing one instance of a primitive against an empty circuit (arm64 build):

| primitive | gates today |
|---|---|
| SHA-256 one compression | 6,703 + **33,944 table rows** in the trace (tables shared by all uses in a circuit) |
| AES-128 CBC, per block | $\approx$ 1,057 |
| Poseidon2 2-to-1 hash | 74 |
| `uint32` XOR (marginal) | 6 (plus the shared 4,096-row table) |
| single 20-bit range constraint | 2,751 unamortized (one-time sorted-list machinery; marginal cost 2–3) |

### Static-set membership: registries without Merkle paths

A membership check against a preprocessed set is one table row; the same check via a depth-32 Poseidon2 Merkle
path costs $32 \times 74 \approx 2{,}400$ gates — a factor of several hundred per check. Candidate sets are
anything public, large, and slow-changing: a contract-class registry, the set of approved protocol-circuit VK
hashes (checked by Merkle proof in kernels today), token or precompile registries, credential allowlists. With an
index column the same argument proves *positional* reads $(i, v_i)$ of a committed $2^{24}$-cell array — a batch
of $k$ random reads costs $k$ rows and one $\approx 10$-G1 argument, against $k$ Merkle paths of 30+ hashes each
(relevant to data-availability-sampling-style statements). Non-membership in a *sorted* static set also reduces
to one lookup: preprocess the table of adjacent pairs $(t_i, t_{i+1})$ and prove $t_i < x < t_{i+1}$ with one row
plus two comparisons (denylist screening).

The update model: changing any entry invalidates all cached quotients, and re-preprocessing is $O(N \log N)$
(18 s at $N = 2^{16}$ on the machine of the Costs table, extrapolating to minutes at $2^{20}$) — acceptable at
registry-update cadence. The preprocessing is publicly recomputable by anyone, so publishing a table's caches has
the same trust profile as publishing the table. Snapshots of *historic* chain state (a finalized epoch's archive-tree leaves, an L1
block-hash history) are static in exactly this sense: an epoch table turns archival membership proofs — light
clients, bridges — into single rows, with one certifying SNARK per epoch attesting that $[T]_2$ commits to the
same leaves as the epoch's Merkle root. For incremental maintenance, the cached quotients are exactly the
all-positions KZG opening proofs of the aSVC vector commitment [12] (FK is how aSVC computes them), so aSVC's
update keys apply: after a change at position $j$, any single cached proof updates in $O(1)$, letting a *user*
maintain their own membership witness under a batch of $B$ table changes in $O(B)$ — while a maintainer of all
$N$ caches still prefers one FK re-run per epoch.

#### Case study: World ID / Semaphore anonymity sets

World ID keeps its $\approx 2^{24}$ identity commitments in a depth-30 Poseidon Merkle tree on Ethereum; a
membership proof is a Groth16/Semaphore circuit containing the 30-hash Merkle path plus the nullifier derivation,
and users must fetch a fresh path from the operator's inclusion-proof service as the tree changes. As a
single-column cq table, the path disappears from the circuit: the user-side SNARK shrinks to the nullifier and
commitment-opening hashes ($\approx 3$ Poseidons, roughly a tenth of the constraints), membership itself becomes
one blinded lookup, and the per-user witness is three cached group elements ($\approx 200$ bytes, self-updatable
via update keys) instead of a 960-byte re-fetched path. Semacaulk [11] validated exactly this architecture on
Caulk-generation machinery — precomputed membership quotients, blinded per use, cleanly separated from the
nullifier circuit — and cq is the same design with the table-size-independent prover. Two caveats are structural:
the membership side must run the zk variant of cq (the base protocol reveals $f(\gamma)$ and the support), and
the set is dynamic, which the epoch model above absorbs at the operator's existing batch-insertion cadence. Where
the swap is strongest is not the single on-chain verify (Groth16's 3 pairings are already cheap) but composition
and batch settings: $k$ users' membership claims are one cq argument with $n = k$ rows, and a recursive verifier
inside another proof system replaces 30 in-circuit hashes with one lookup row plus pairing accumulation. The same
shape covers Privacy-Pools-style association sets and any Semaphore group.

### A universal ALU table for VMs

One table with rows $(\mathrm{opcode} \,\|\, a \,\|\, b) \to \mathrm{out}$ over 8-bit operands is $16 \cdot
2^{16} = 2^{20}$ rows and serves every bitwise, arithmetic, comparison and shift instruction at once; u32
operations decompose into 4 chunk lookups plus carry glue. Output columns are nearly free (they combine
homomorphically under $\theta$), so the same row can carry $(\mathrm{result}, \mathrm{carry}, \mathrm{zero},
\mathrm{overflow})$ — flags that each cost their own relation in a hand-built VM arithmetization. This collapses
a large fraction of an AVM-style instruction set into one preprocessed table plus a handful of rows per
instruction, the construction Jolt [10] builds from Lasso, here obtained natively in BN254/KZG with no new proof
machinery. For a *fixed* program, bytecode itself is a static table $(\mathrm{pc}, \mathrm{insn})$, making
instruction fetch one lookup per cycle at the cost of per-program preprocessing (seconds at $2^{16}$).

### zkTLS: AES-GCM

Software AES's T-table formulation becomes provable directly: four $2^8$-row tables mapping a byte to its fused
SubBytes+MixColumns 32-bit column contribution, folded with wide-chunk XOR lookups — an estimated 150–300 rows
per block against the measured 1,057. For GHASH, the multiplier $H$ is session-derived so no table can depend on
it, but $\mathrm{GF}(2^{128})$ multiplication decomposes into carryless $8 \times 8$-bit chunk products (a
$2^{16}$ table) combined by Karatsuba — hundreds of rows against thousands of gates. TLS transcripts are
megabyte-scale AEAD, so a 3–5$\times$ block cost reduction moves end-to-end proving time by the same factor.

### Hash functions: decoupling now, wider slices next

The measured SHA-256 numbers make the trace-decoupling mechanism concrete: the tables (33,944 rows) are five
times the compute (6,703 gates), so a prove-one-preimage circuit shrinks from dyadic $2^{16}$ to $2^{13}$ before
any table is widened. Widening then attacks the compute: the sparse normalization tables process 2–3 bits per
lookup (choose: $28^2 = 784$ rows, majority: $16^3 = 4096$ rows), so a 32-bit normalization costs 11–16 lookups;
doubling slice widths (majority $16^6 \approx 2^{24}$, choose $28^4 \approx 2^{19.2}$) halves the lookups of
every normalization step, an estimated 6,703 → 3,500–4,000 gates per compression. The same applies to Blake
(11-bit XOR slices: the measured 6-gate u32 XOR becomes 3) and Keccak's base-11 normalization tables.

### Byte-format proofs: base64, JWTs, regex

Base64 fits cq exactly: 4 characters decode to 3 bytes and $64^4 = 2^{24}$, so JWT/email-style payload decoding
costs one lookup per 3 output bytes instead of $\approx 10$ gates per byte. Regex and other DFA claims step
$(\mathrm{state}, \mathrm{byte}_1, \mathrm{byte}_2) \to \mathrm{state}$ through a $2^8$-state automaton at 2
input bytes per row ($2^{24}$ table), replacing several gates per byte; canonicalization maps (lowercasing,
UTF-8 validation, hex/RLP decoding) are the same shape.

### Non-native arithmetic: the range-check economy

RSA and ECDSA circuits are dominated by limb-decomposition range checks attached to every non-native
multiplication. Each such check currently rides the 14-bit sorted-list machinery (2,751 gates before
amortization, 2–3 after); a $2^{20}$ range table prices any 20-bit check at one row with *no* fixed machinery,
which both cuts the amortized cost and removes the fixed cost entirely from small circuits. A $2^{16}$
reciprocal-seed table (one Newton step from a 16-bit seed) similarly shortcuts integer division and modular
reduction digits. Estimated effect on an RSA-2048 verification (a passport-credential workload): 1.5–2.5$\times$
from ranges alone, on top of the hash-side savings above. Field-to-binary decomposition — everywhere in scalar
arithmetic — reads 20 bits per row from a decomposition table instead of one gate per bit.

### zkML: nonlinearities and quantization

Fixed-point activation functions (GELU, sigmoid, exp, rsqrt) at 16–20-bit precision are unary tables: one row
each, against hundreds of gates of polynomial approximation plus range checks today. Quantization of 24-bit
accumulators to int8 (scale, clamp, round) is one $2^{24}$ lookup. Narrow floating point (fp8 directly at
$2^{16}$; bf16 with mantissa/exponent handling split between tables and gates) brings table-priced arithmetic to
inference circuits whose current cost is dominated by exactly these steps.

### Fixed-base elliptic-curve windows

The secp256r1 fixed-base path chains 32 per-window tables of $2^8$ points; $2^{16}$–$2^{20}$-point windows halve
or better the window additions per scalar multiplication. This matters where Goblin/ECCVM is not available to
absorb the curve arithmetic (plain UltraHonk deployments); inside the Aztec stack Goblin already amortizes these.

### Erasure-coding and DA encoding proofs

Reed–Solomon encoding checks over $\mathrm{GF}(2^8)$/$\mathrm{GF}(2^{16})$ reduce to chunked carryless
multiplications against fixed generator constants — the same $2^{16}$ carryless-product table as GHASH, with the
per-symbol constants foldable into per-position tables. Speculative relative to the entries above, but the
operation profile (many small-field multiplications by fixed constants) is the best possible fit for
preprocessed lookups.

### One universe table, and lookups as a deferred side-protocol

Because $N$ never enters proving cost, every table in the system can be concatenated — with a table-id column —
into one universe table, so a proof carries one cq argument ($\approx 10$ G1) regardless of how many distinct
tables the circuit touches, and no circuit pays trace rows for any of them. Pushing the same observation across
proofs: cq claims can be *deferred* through an IVC chain the way Goblin defers EC operations — accumulate
$(\mathrm{table{:}id}, \mathrm{row})$ claims across Chonk steps and discharge a single argument for the whole
chain at the end, removing lookup machinery from every kernel step. Deferral also shifts work out of witness
generation: a looked-up row is a precomputed result, not a recomputed one.

### Boundaries

- Tables must be static and instance-independent: no live nullifier/note trees, no RAM, no tables keyed by a
  session or public key (the GHASH $H$-powers and pubkey-window ideas fail this test; historic snapshots pass).
- Practical ceiling $N \approx 2^{24}..2^{26}$ from cache storage and SRS length (public powers-of-tau reach
  $2^{28}$). Two-operand 16-bit tables ($2^{32}$) stay out of reach; that regime belongs to decomposable-table
  arguments (Lasso [8]), and a hybrid — virtual wide tables whose chunks are cq tables — is the natural bridge.
- A lookup still costs about one trace row of witness wiring; cq removes the table's cost and the per-table
  argument overhead, while fusion is what reduces lookup counts.
- This implementation is non-zk, and in-circuit (recursive) verification needs the pairing checks folded into the
  existing `PairingPoints` aggregation with $n$ fixed per circuit (the degree-check G2 point depends on $n$).

## References

1. cq: Cached quotients for fast lookups — L. Eagen, D. Fiore, A. Gabizon. https://eprint.iacr.org/2022/1763
2. Caulk — A. Zapico, V. Buterin, D. Khovratovich, M. Maller, A. Nitulescu, M. Simkin. https://eprint.iacr.org/2022/621
3. Caulk+ — J. Posen, A. Kattis. https://eprint.iacr.org/2022/957
4. Flookup — A. Gabizon, D. Khovratovich. https://eprint.iacr.org/2022/1447
5. Baloo — A. Zapico, A. Gabizon, D. Khovratovich, M. Maller, C. Ràfols. https://eprint.iacr.org/2022/1565
6. Multivariate lookups based on logarithmic derivatives — U. Haböck. https://eprint.iacr.org/2022/1530
7. Fast amortized KZG proofs — D. Feist, D. Khovratovich. https://eprint.iacr.org/2023/033
8. Unlocking the lookup singularity with Lasso — S. Setty, J. Thaler, R. Wahby. https://eprint.iacr.org/2023/1216
9. SoK: Lookup Table Arguments. https://eprint.iacr.org/2025/1876
10. Jolt: SNARKs for Virtual Machines via Lookups — A. Arun, S. Setty, J. Thaler. https://eprint.iacr.org/2023/1217
11. Semacaulk — Geometry research (K. W. J. Koh et al.). https://github.com/geometryxyz/semacaulk
12. Aggregatable Subvector Commitments for Stateless Cryptocurrencies — A. Tomescu, I. Abraham, V. Buterin,
    J. Drake, D. Feist, D. Khovratovich. https://eprint.iacr.org/2020/527
