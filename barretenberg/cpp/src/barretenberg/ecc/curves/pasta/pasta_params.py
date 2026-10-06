#!/usr/bin/env python3
"""
Derives every constant barretenberg's `field<Params>` and `group<...>` templates need for the Pasta cycle
(Pallas / Vesta) and prints them as C++ (`pasta_params.py > /tmp/out.hpp`). The printed blocks are pasted into
pasta.hpp. Running the script also checks:
  * the cube roots of unity are paired so that [lambda]P = (beta * x, y) on both curves, where the Fp cube root is
    beta for Pallas and lambda for Vesta (and symmetrically for Fq), as barretenberg stores one cube root per field;
  * the short GLV lattice basis is too long for barretenberg's 128-bit split (see USE_ENDOMORPHISM in pasta.hpp).

Pasta parameters (https://github.com/zcash/pasta):
  p = 0x40000000000000000000000000000000224698fc094cf91b992d30ed00000001   (Pallas base, Vesta scalar)
  q = 0x40000000000000000000000000000000224698fc0994a8dd8c46eb2100000001   (Pallas scalar, Vesta base)
  Pallas: y^2 = x^3 + 5 over F_p, order q.   Vesta: y^2 = x^3 + 5 over F_q, order p.
  Generators (as in pasta_curves): (-1, 2) on both curves.
"""
import sys

P = 0x40000000000000000000000000000000224698FC094CF91B992D30ED00000001
Q = 0x40000000000000000000000000000000224698FC0994A8DD8C46EB2100000001
B = 5
MULTIPLICATIVE_GENERATOR = 5  # pasta_curves GENERATOR for both fields


def limbs(x, bits, n):
    return [(x >> (bits * i)) & ((1 << bits) - 1) for i in range(n)]


def two_adicity(m):
    s = 0
    t = m - 1
    while t % 2 == 0:
        t //= 2
        s += 1
    return s, t


def nontrivial_cube_roots(m):
    # m - 1 is divisible by 3 for both Pasta primes
    assert (m - 1) % 3 == 0
    for g in range(2, 100):
        w = pow(g, (m - 1) // 3, m)
        if w != 1:
            return [w, w * w % m]
    raise ValueError("no cube root found")


# ---------------------------------------------------------------- curve arithmetic (affine, Python ints)
def ec_add(P1, P2, m):
    if P1 is None:
        return P2
    if P2 is None:
        return P1
    (x1, y1), (x2, y2) = P1, P2
    if x1 == x2:
        if (y1 + y2) % m == 0:
            return None
        lam = 3 * x1 * x1 * pow(2 * y1, -1, m) % m
    else:
        lam = (y2 - y1) * pow(x2 - x1, -1, m) % m
    x3 = (lam * lam - x1 - x2) % m
    return (x3, (lam * (x1 - x3) - y1) % m)


def ec_mul(k, Pt, m):
    acc = None
    while k:
        if k & 1:
            acc = ec_add(acc, Pt, m)
        Pt = ec_add(Pt, Pt, m)
        k >>= 1
    return acc


# ---------------------------------------------------------------- GLV lattice (same routine as ../../fields/endomorphism_scalars.py)
def bb_short_basis(lam, n):
    from math import isqrt
    approx_sqrt = isqrt(n)
    remainder, prev_remainder = lam, n
    coeff, prev_coeff = 1, 0
    while remainder >= approx_sqrt:
        quot = prev_remainder // remainder
        prev_remainder, remainder = remainder, prev_remainder - quot * remainder
        prev_coeff, coeff = coeff, prev_coeff - quot * coeff
    vec_before = (-prev_remainder, prev_coeff)
    quot = prev_remainder // remainder
    r_after = prev_remainder - quot * remainder
    s_after = prev_coeff - quot * coeff
    a1, b1 = (-remainder, coeff)
    a2, b2 = (-r_after, s_after) if (r_after**2 + s_after**2 < prev_remainder**2 + prev_coeff**2) else vec_before
    assert (a1 + lam * b1) % n == 0 and (a2 + lam * b2) % n == 0
    return a1, b1, a2, b2


def mont(x, m, rbits=256):
    return x * (1 << rbits) % m


def emit_field(name, m, cube_root, schema):
    s, t = two_adicity(m)
    root = pow(MULTIPLICATIVE_GENERATOR, t, m)  # primitive 2^s-th root of unity
    assert pow(root, 1 << (s - 1), m) == m - 1
    coset = MULTIPLICATIVE_GENERATOR
    out = []
    w = out.append
    w(f"    // modulus = {hex(m)}")
    for i, v in enumerate(limbs(m, 64, 4)):
        w(f"    static constexpr uint64_t modulus_{i} = {hex(v)}ULL;")
    for i, v in enumerate(limbs(pow(2, 512, m), 64, 4)):
        w(f"    static constexpr uint64_t r_squared_{i} = {hex(v)}ULL;")
    w(f"    static constexpr uint64_t r_inv = {hex((-pow(m, -1, 1 << 64)) % (1 << 64))}ULL;")
    for i, v in enumerate(limbs(pow(2, -64, m), 64, 4)):
        w(f"    static constexpr uint64_t r_inv_{i} = {hex(v)}ULL;")
    for i, v in enumerate(limbs(mont(cube_root, m), 64, 4)):
        w(f"    static constexpr uint64_t cube_root_{i} = {hex(v)}ULL;")
    for i, v in enumerate(limbs(mont(root, m), 64, 4)):
        w(f"    static constexpr uint64_t primitive_root_{i} = {hex(v)}ULL;")
    for i, v in enumerate(limbs(mont(coset, m), 64, 4)):
        w(f"    static constexpr uint64_t coset_generator_{i} = {hex(v)}ULL;")
    for i, v in enumerate(limbs(m, 29, 9)):
        w(f"    static constexpr uint64_t modulus_wasm_{i} = {hex(v)};")
    for i, v in enumerate(limbs(pow(2, 522, m), 64, 4)):
        w(f"    static constexpr uint64_t r_squared_wasm_{i} = {hex(v)}ULL;")
    for i, v in enumerate(limbs(pow(2, -29, m), 29, 9)):
        w(f"    static constexpr uint64_t r_inv_wasm_{i} = {hex(v)};")
    for i, v in enumerate(limbs(mont(cube_root, m, 261), 64, 4)):
        w(f"    static constexpr uint64_t cube_root_wasm_{i} = {hex(v)}ULL;")
    for i, v in enumerate(limbs(mont(root, m, 261), 64, 4)):
        w(f"    static constexpr uint64_t primitive_root_wasm_{i} = {hex(v)}ULL;")
    for i, v in enumerate(limbs(mont(coset, m, 261), 64, 4)):
        w(f"    static constexpr uint64_t coset_generator_wasm_{i} = {hex(v)}ULL;")
    # barretenberg's MSM GLV path needs both split halves to fit in 128 bits. For Pasta the short lattice basis has a
    # 128-bit coordinate, and the floor-rounded split used by field::split_into_endomorphism_scalars then overflows
    # 128 bits for a constant fraction of scalars, so the curves run with USE_ENDOMORPHISM = false and the splitting
    # constants are left zero.
    a1, b1, a2, b2 = bb_short_basis(cube_root, m)
    assert max(abs(a1), abs(b1), abs(a2), abs(b2)).bit_length() == 128
    for nm in ("endo_g1_lo", "endo_g1_mid", "endo_g1_hi", "endo_g2_lo", "endo_g2_mid", "endo_minus_b1_lo",
               "endo_minus_b1_mid", "endo_b2_lo", "endo_b2_mid"):
        w(f"    static constexpr uint64_t {nm} = 0;")
    w(f"    // 2-adicity = {s}")
    return f"// ---- {name} ({schema})\n" + "\n".join(out)


def main():
    # Pick cube roots so that one constant per field serves both as beta (base-field role) and lambda (scalar role).
    gen_p = (P - 1, 2)  # Pallas generator over F_p
    gen_q = (Q - 1, 2)  # Vesta generator over F_q
    assert (gen_p[1] ** 2 - gen_p[0] ** 3 - B) % P == 0
    assert (gen_q[1] ** 2 - gen_q[0] ** 3 - B) % Q == 0
    chosen = None
    for zp in nontrivial_cube_roots(P):
        for zq in nontrivial_cube_roots(Q):
            pallas_ok = ec_mul(zq, gen_p, P) == (zp * gen_p[0] % P, gen_p[1])
            vesta_ok = ec_mul(zp, gen_q, Q) == (zq * gen_q[0] % Q, gen_q[1])
            if pallas_ok and vesta_ok:
                chosen = (zp, zq)
    assert chosen is not None, "no consistent pair of cube roots"
    zp, zq = chosen
    print(f"// Fp cube root (pasta_curves Fp::ZETA = 0x12ccca83...): {hex(zp)}", file=sys.stderr)
    print(f"// Fq cube root: {hex(zq)}", file=sys.stderr)
    print(emit_field("Fp: Pallas base field / Vesta scalar field", P, zp, "pallas_fq"))
    print()
    print(emit_field("Fq: Pallas scalar field / Vesta base field", Q, zq, "pallas_fr"))


if __name__ == "__main__":
    main()
