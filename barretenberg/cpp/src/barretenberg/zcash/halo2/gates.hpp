#pragma once

#include "barretenberg/zcash/halo2/layout.hpp"
#include "barretenberg/zcash/primitives/cycle.hpp"
#include "barretenberg/zcash/primitives/poseidon.hpp"

#include <array>
#include <cstddef>

/**
 * @file gates.hpp
 * @brief The custom gates of the halo2 Orchard Action circuit (orchard 0.16 / halo2_gadgets 0.6), written once and
 * evaluated both natively (constraint checker) and on sumcheck univariates (Honk relation).
 *
 * @details Every halo2 gate `q_s * [p_1, ..., p_m]` contributes m polynomial identities `q_s * p_j = 0`. halo2 queries
 * advice cells at rotations -1, 0, +1 relative to the selector row. Honk only queries "current" and "shifted" values,
 * so each gate is re-anchored at its first queried row: a gate that queries Rotation::prev() is anchored one row
 * earlier
 * (`ANCHOR_SHIFT[s] = 1`) and its selector polynomial is stored one row earlier as well. Advice cells are then read at
 * offsets 0, 1, 2 from the anchor (current, shift, shift-by-two). Fixed columns are only ever queried at
 * Rotation::cur() by gates with ANCHOR_SHIFT = 0.
 *
 * The evaluator takes an input accessor `in` with
 *   T   in.adv(size_t column, size_t offset)  advice column at anchor + offset (offset in {0, 1, 2})
 *   T   in.fix(size_t column)                 fixed column at the anchor
 *   T   in.sel(size_t selector)               selector at the anchor
 *   bool in.active(size_t selector)           false only if the selector is identically zero on the inputs
 * and an emitter `emit(size_t index, const T& value)` that receives every constraint already multiplied by its
 * selector. Constraints of inactive gates are not emitted; indices are stable regardless.
 */
namespace bb::zcash::halo2 {

using namespace layout;

/// Number of rows by which each selector (and its gate) is anchored before halo2's selector row.
inline constexpr std::array<size_t, NUM_SELECTORS> ANCHOR_SHIFT = [] {
    std::array<size_t, NUM_SELECTORS> shift{};
    for (size_t s : { Q_MUL_HI_2,
                      Q_MUL_HI_3,
                      Q_MUL_LO_2,
                      Q_MUL_LO_3,
                      Q_MUL_DECOMPOSE_VAR,
                      Q_MUL_OVERFLOW,
                      Q_MUL_FIXED_BASE_FIELD,
                      Q_POSEIDON_PAD_AND_ADD,
                      Q_BITSHIFT }) {
        shift[s] = 1;
    }
    return shift;
}();

/**
 * @brief Number of polynomial identities of each gate. The Pasta circuit is orchard's; the BN254/Grumpkin port replaces
 * the gadgets that depend on the Pasta moduli (canonicity checks, variable-base multiplication overflow check), which
 * changes the identities of those gates (see OrchardGates::evaluate).
 */
constexpr std::array<size_t, NUM_SELECTORS> gate_sizes(bool pasta)
{
    std::array<size_t, NUM_SELECTORS> n{};
    n[Q_ORCHARD] = 4;
    n[Q_ADD_FIELD] = 1;
    n[Q_POINT] = 2;
    n[Q_POINT_NON_ID] = 1;
    n[Q_ADD_INCOMPLETE] = 2;
    n[Q_ADD] = 12;
    n[Q_MUL_HI_1] = n[Q_MUL_LO_1] = 1;
    n[Q_MUL_HI_2] = n[Q_MUL_LO_2] = 6;
    n[Q_MUL_HI_3] = n[Q_MUL_LO_3] = 4;
    n[Q_MUL_DECOMPOSE_VAR] = 2;
    n[Q_MUL_LSB] = 3;
    n[Q_MUL_OVERFLOW] = pasta ? 5 : 7;
    n[Q_MUL_FIXED_RUNNING_SUM] = 4;
    n[Q_MUL_FIXED_FULL] = 4;
    n[Q_MUL_FIXED_SHORT] = 4;
    n[Q_MUL_FIXED_BASE_FIELD] = pasta ? 8 : 5;
    n[Q_POSEIDON_FULL] = 3;
    n[Q_POSEIDON_PARTIAL] = 4;
    n[Q_POSEIDON_PAD_AND_ADD] = 3;
    n[Q_SINSEMILLA1_1] = n[Q_SINSEMILLA1_2] = 2;
    n[Q_SINSEMILLA4_1] = n[Q_SINSEMILLA4_2] = 1;
    n[Q_MERKLE_DECOMPOSE_1] = n[Q_MERKLE_DECOMPOSE_2] = 4;
    n[Q_SWAP_1] = n[Q_SWAP_2] = 3;
    n[Q_LOOKUP] = 0;
    n[Q_RUNNING] = 0;
    n[Q_BITSHIFT] = 1;
    n[Q_COMMIT_IVK] = pasta ? 14 : 17;
    n[Q_NOTECOMMIT_B] = 3;
    n[Q_NOTECOMMIT_D] = 3;
    n[Q_NOTECOMMIT_E] = 1;
    n[Q_NOTECOMMIT_G] = 2;
    n[Q_NOTECOMMIT_H] = 2;
    n[Q_NOTECOMMIT_G_D] = pasta ? 5 : 7;
    n[Q_NOTECOMMIT_PK_D] = pasta ? 4 : 8;
    n[Q_NOTECOMMIT_VALUE] = 1;
    n[Q_NOTECOMMIT_RHO] = pasta ? 4 : 8;
    n[Q_NOTECOMMIT_PSI] = pasta ? 5 : 9;
    n[Q_Y_CANON] = pasta ? 7 : 8;
    return n;
}

/// Gate sizes, offsets and the selector of every constraint, for the Pasta circuit or the BN254 port.
template <bool PASTA> struct GateMetadata {
    static constexpr std::array<size_t, NUM_SELECTORS> GATE_SIZE = gate_sizes(PASTA);

    static constexpr size_t NUM_GATE_CONSTRAINTS = [] {
        size_t n = 0;
        for (size_t s : GATE_SIZE) {
            n += s;
        }
        return n;
    }();

    /// Index of the first constraint of each gate.
    static constexpr std::array<size_t, NUM_SELECTORS> GATE_OFFSET = [] {
        std::array<size_t, NUM_SELECTORS> out{};
        size_t n = 0;
        for (size_t s = 0; s < NUM_SELECTORS; ++s) {
            out[s] = n;
            n += GATE_SIZE[s];
        }
        return out;
    }();

    /// The selector gating each constraint.
    static constexpr std::array<size_t, NUM_GATE_CONSTRAINTS> CONSTRAINT_SELECTOR = [] {
        std::array<size_t, NUM_GATE_CONSTRAINTS> out{};
        size_t i = 0;
        for (size_t s = 0; s < NUM_SELECTORS; ++s) {
            for (size_t j = 0; j < GATE_SIZE[s]; ++j) {
                out[i++] = s;
            }
        }
        return out;
    }();
};

/// Highest total degree (selector included) over all gates: the fixed-base Lagrange interpolation and the 3-bit window
/// range checks have degree 8 in the advice cells.
inline constexpr size_t MAX_GATE_DEGREE = 9;

// Pieces of the sinsemilla configurations needed by the gates and the lookup.
struct SinsemillaColumns {
    size_t x_a, x_p, bits, lambda_1, lambda_2;
    size_t fixed_y_q, q_s2;
};
inline constexpr std::array<SinsemillaColumns, 2> SINSEMILLA_COLUMNS = {
    SinsemillaColumns{ 0, 1, 2, 3, 4, F0, Q_SINSEMILLA2_1 },
    SinsemillaColumns{ 5, 6, 7, 8, 9, F1, Q_SINSEMILLA2_2 },
};
inline constexpr size_t RANGE_CHECK_COLUMN = 9;

template <typename Cycle> struct OrchardGates {
    using FF = typename Cycle::FF;
    static constexpr bool IS_PASTA = std::is_same_v<Cycle, PastaCycle>;
    using Metadata = GateMetadata<IS_PASTA>;
    static constexpr auto GATE_SIZE = Metadata::GATE_SIZE;
    static constexpr auto GATE_OFFSET = Metadata::GATE_OFFSET;
    static constexpr size_t NUM_GATE_CONSTRAINTS = Metadata::NUM_GATE_CONSTRAINTS;
    static constexpr auto CONSTRAINT_SELECTOR = Metadata::CONSTRAINT_SELECTOR;
    using Poseidon = PoseidonP128Pow5T3<FF>;

    static FF two_pow(size_t k) { return FF(uint256_t(1) << k); }
    static FF t_p() { return FF(uint256_t(FF::modulus) - (uint256_t(1) << 254)); }
    static FF t_q() { return FF(uint256_t(Cycle::EmbeddedScalar::modulus) - (uint256_t(1) << 254)); }

    // ---- BN254 / Grumpkin port constants
    // A field element x < r encoded on 255 bits is canonical iff bit 254 is 0 and, when bits 252..253 are both set,
    // bits 250..251 are 0 and the low 250 bits are below t_r = r - 3 * 2^252 (~2^246.6).
    static uint256_t t_r_uint() { return uint256_t(FF::modulus) - (uint256_t(3) << 252); }
    static FF t_r() { return FF(t_r_uint()); }
    // Number of 10-bit words of the canonicity range check (L + 2^250 - t_r < 2^250).
    static constexpr size_t CANONICITY_WORDS = IS_PASTA ? 13 : 25;
    // Variable-base mul computes [k + 2^254] T for the 255-bit k it decomposes; on Grumpkin k = alpha + t' with
    // t' = 2 q - 2^254, and the overflow check proves t' <= k < t' + r over the integers.
    static uint256_t var_mul_offset()
    {
        const uint256_t q(Cycle::EmbeddedScalar::modulus);
        return IS_PASTA ? q - (uint256_t(1) << 254) : (q + q) - (uint256_t(1) << 254);
    }
    static uint256_t var_mul_upper() { return var_mul_offset() + uint256_t(FF::modulus) - 1; }
    static FF lo130(const uint256_t& x) { return FF(x.slice(0, 130)); }
    static FF hi130(const uint256_t& x) { return FF(x.slice(130, 256)); }
    // Indicator of t == 3 for t in [0, 4).
    template <typename T> static T is_three(const T& t)
    {
        static const FF inv6 = FF(6).invert();
        return t * (t - FF(1)) * (t - FF(2)) * inv6;
    }

    template <typename T> static T range_check(const T& word, size_t range)
    {
        T acc = word;
        for (size_t i = 1; i < range; ++i) {
            acc *= (FF(i) - word);
        }
        return acc;
    }
    template <typename T> static T bool_check(const T& x) { return range_check(x, 2); }
    // a * b + (1 - a) * c
    template <typename T> static T ternary(const T& a, const T& b, const T& c) { return a * (b - c) + c; }

    /**
     * @brief The log-derivative lookup read terms of a row (anchored at the row).
     * @details Range check (halo2 LookupRangeCheckConfig): active iff q_lookup; reads
     *   q_running * (z_cur - 2^K z_next) + (1 - q_running) * z_cur.
     * Sinsemilla (one per chip): active iff q_sinsemilla1; reads (m, x_p, y_p) with
     *   m = z_cur - q_run * 2^K * z_next,  q_run = q_s2 - q_s2 (q_s2 - 1),
     *   y_p = Y_A / 2 - lambda_1 (x_a - x_p).
     */
    template <typename T, typename In> static T range_lookup_value(const In& in)
    {
        const T q_running = in.sel(Q_RUNNING);
        const T z_cur = in.adv(RANGE_CHECK_COLUMN, 0);
        const T z_next = in.adv(RANGE_CHECK_COLUMN, 1);
        return z_cur - q_running * z_next * two_pow(10);
    }

    template <typename T, typename In> static std::array<T, 3> sinsemilla_lookup_value(const In& in, size_t chip)
    {
        const auto& c = SINSEMILLA_COLUMNS[chip];
        const T q_s2 = in.fix(c.q_s2);
        const T q_s3 = q_s2 * (q_s2 - FF(1));
        const T q_run = q_s2 - q_s3;
        const T word = in.adv(c.bits, 0) - q_run * in.adv(c.bits, 1) * two_pow(10);
        const T x_a = in.adv(c.x_a, 0);
        const T x_p = in.adv(c.x_p, 0);
        const T l1 = in.adv(c.lambda_1, 0);
        const T l2 = in.adv(c.lambda_2, 0);
        const T x_r = l1.sqr() - x_a - x_p;
        const T Y_A = (l1 + l2) * (x_a - x_r);
        const T y_p = Y_A * FF(2).invert() - l1 * (x_a - x_p);
        return { word, x_p, y_p };
    }

    template <typename T, typename In, typename Emit> static void evaluate(const In& in, Emit&& emit)
    {
        const FF one(1);
        const FF two(2);
        const FF two_inv = FF(2).invert();
        const FF curve_b = Cycle::Group::curve_b;

        size_t emitted = 0;
        // Runs `body(constrain)` if selector `s` is active; `constrain(expr)` emits `q_s * expr`.
        auto gate = [&](size_t s, auto&& body) {
            if (in.active(s)) {
                const T q = in.sel(s);
                size_t idx = GATE_OFFSET[s];
                body([&](const T& expr) { emit(idx++, q * expr); });
                BB_ASSERT_EQ(idx - GATE_OFFSET[s], GATE_SIZE[s]);
            }
            emitted += GATE_SIZE[s];
        };
        // Advice cell at halo2 rotation `rot` for a gate with anchor shift `base`.
        auto A = [&](size_t base, size_t col, int rot) -> T {
            return in.adv(col, static_cast<size_t>(static_cast<int>(base) + rot));
        };

        // ---- Orchard circuit checks
        gate(Q_ORCHARD, [&](auto c) {
            const T v_old = A(0, 0, 0);
            const T v_new = A(0, 1, 0);
            const T magnitude = A(0, 2, 0);
            const T sign = A(0, 3, 0);
            const T root = A(0, 4, 0);
            const T anchor = A(0, 5, 0);
            const T enable_spend = A(0, 6, 0);
            const T enable_output = A(0, 7, 0);
            c(v_old - v_new - magnitude * sign);
            c(v_old * (root - anchor));
            c(v_old * (one - enable_spend));
            c(v_new * (one - enable_output));
        });
        // ---- Field element addition c = a + b (a = a7, b = a8, c = a6)
        gate(Q_ADD_FIELD, [&](auto c) { c(A(0, 7, 0) + A(0, 8, 0) - A(0, 6, 0)); });

        // ---- ECC: witness point (x = a0, y = a1)
        auto curve_eqn = [&](const T& x, const T& y) { return y.sqr() - x.sqr() * x - curve_b; };
        gate(Q_POINT, [&](auto c) {
            const T x = A(0, 0, 0);
            const T y = A(0, 1, 0);
            const T eq = curve_eqn(x, y);
            c(x * eq);
            c(y * eq);
        });
        gate(Q_POINT_NON_ID, [&](auto c) { c(curve_eqn(A(0, 0, 0), A(0, 1, 0))); });

        // ---- ECC: incomplete addition (x_p = a0, y_p = a1, x_qr = a2, y_qr = a3)
        gate(Q_ADD_INCOMPLETE, [&](auto c) {
            const T x_p = A(0, 0, 0);
            const T y_p = A(0, 1, 0);
            const T x_q = A(0, 2, 0);
            const T y_q = A(0, 3, 0);
            const T x_r = A(0, 2, 1);
            const T y_r = A(0, 3, 1);
            const T dx = x_p - x_q;
            c((x_r + x_q + x_p) * dx * dx - (y_p - y_q).sqr());
            c((y_r + y_q) * dx - (y_p - y_q) * (x_q - x_r));
        });

        // ---- ECC: complete addition (lambda = a4, alpha = a5, beta = a6, gamma = a7, delta = a8)
        gate(Q_ADD, [&](auto c) {
            const T x_p = A(0, 0, 0);
            const T y_p = A(0, 1, 0);
            const T x_q = A(0, 2, 0);
            const T y_q = A(0, 3, 0);
            const T x_r = A(0, 2, 1);
            const T y_r = A(0, 3, 1);
            const T lambda = A(0, 4, 0);
            const T alpha = A(0, 5, 0);
            const T beta = A(0, 6, 0);
            const T gamma = A(0, 7, 0);
            const T delta = A(0, 8, 0);

            const T x_q_minus_x_p = x_q - x_p;
            const T x_p_minus_x_r = x_p - x_r;
            const T y_q_plus_y_p = y_q + y_p;
            const T if_alpha = x_q_minus_x_p * alpha;
            const T if_beta = x_p * beta;
            const T if_gamma = x_q * gamma;
            const T if_delta = y_q_plus_y_p * delta;

            c(x_q_minus_x_p * (x_q_minus_x_p * lambda - (y_q - y_p)));
            c((one - if_alpha) * (y_p * two * lambda - x_p.sqr() * FF(3)));
            const T nonexceptional_x_r = lambda.sqr() - x_p - x_q - x_r;
            const T nonexceptional_y_r = lambda * x_p_minus_x_r - y_p - y_r;
            const T x_p_x_q = x_p * x_q;
            c(x_p_x_q * x_q_minus_x_p * nonexceptional_x_r);
            c(x_p_x_q * x_q_minus_x_p * nonexceptional_y_r);
            c(x_p_x_q * y_q_plus_y_p * nonexceptional_x_r);
            c(x_p_x_q * y_q_plus_y_p * nonexceptional_y_r);
            c((one - if_beta) * (x_r - x_q));
            c((one - if_beta) * (y_r - y_q));
            c((one - if_gamma) * (x_r - x_p));
            c((one - if_gamma) * (y_r - y_p));
            const T not_alpha_delta = one - if_alpha - if_delta;
            c(not_alpha_delta * x_r);
            c(not_alpha_delta * y_r);
        });

        // ---- ECC: variable-base mul, incomplete double-and-add (hi: z = a9, x_a = a3, lambda = a4, a5;
        //      lo: z = a6, x_a = a7, lambda = a8, a2; both: x_p = a0, y_p = a1)
        struct DoubleAndAdd {
            size_t z, x_a, lambda_1, lambda_2;
            size_t q_1, q_2, q_3;
        };
        for (const DoubleAndAdd& d : { DoubleAndAdd{ 9, 3, 4, 5, Q_MUL_HI_1, Q_MUL_HI_2, Q_MUL_HI_3 },
                                       DoubleAndAdd{ 6, 7, 8, 2, Q_MUL_LO_1, Q_MUL_LO_2, Q_MUL_LO_3 } }) {
            auto x_r = [&](size_t base, int rot) {
                return A(base, d.lambda_1, rot).sqr() - A(base, d.x_a, rot) - A(base, 0, rot);
            };
            auto y_a = [&](size_t base, int rot) {
                return (A(base, d.lambda_1, rot) + A(base, d.lambda_2, rot)) * (A(base, d.x_a, rot) - x_r(base, rot)) *
                       two_inv;
            };
            auto for_loop = [&](auto& c, size_t base, const T& y_a_next) {
                const T z_cur = A(base, d.z, 0);
                const T z_prev = A(base, d.z, -1);
                const T x_a_cur = A(base, d.x_a, 0);
                const T x_a_next = A(base, d.x_a, 1);
                const T x_p_cur = A(base, 0, 0);
                const T y_p_cur = A(base, 1, 0);
                const T lambda1_cur = A(base, d.lambda_1, 0);
                const T lambda2_cur = A(base, d.lambda_2, 0);
                const T y_a_cur = y_a(base, 0);
                const T k = z_cur - z_prev * two;
                c(bool_check(k));
                c(lambda1_cur * (x_a_cur - x_p_cur) - y_a_cur + (k * two - one) * y_p_cur);
                c(lambda2_cur.sqr() - x_a_next - x_r(base, 0) - x_a_cur);
                c(lambda2_cur * (x_a_cur - x_a_next) - y_a_cur - y_a_next);
            };
            gate(d.q_1, [&](auto c) { c(A(0, d.lambda_1, 0) - y_a(0, 1)); });
            gate(d.q_2, [&](auto c) {
                c(A(1, 0, 0) - A(1, 0, 1));
                c(A(1, 1, 0) - A(1, 1, 1));
                for_loop(c, 1, y_a(1, 1));
            });
            gate(d.q_3, [&](auto c) { for_loop(c, 1, A(1, d.lambda_1, 1)); });
        }

        // ---- ECC: variable-base mul, complete bits (z_complete = a9, y_p = a1)
        gate(Q_MUL_DECOMPOSE_VAR, [&](auto c) {
            const T k = A(1, 9, 1) - A(1, 9, -1) * two;
            const T base_y = A(1, 9, 0);
            const T y_p = A(1, 1, -1);
            c(bool_check(k));
            c(ternary(k, base_y - y_p, base_y + y_p));
        });
        // ---- ECC: variable-base mul, LSB
        gate(Q_MUL_LSB, [&](auto c) {
            const T z_1 = A(0, 9, 0);
            const T z_0 = A(0, 9, 1);
            const T x_p = A(0, 0, 0);
            const T y_p = A(0, 1, 0);
            const T base_x = A(0, 0, 1);
            const T base_y = A(0, 1, 1);
            const T lsb = z_0 - z_1 * two;
            c(bool_check(lsb));
            c(ternary(lsb, x_p, x_p - base_x));
            c(ternary(lsb, y_p, y_p + base_y));
        });
        // ---- ECC: variable-base mul, overflow check (a6, a7, a8)
        gate(Q_MUL_OVERFLOW, [&](auto c) {
            if constexpr (IS_PASTA) {
                const T z_0 = A(1, 6, -1);
                const T z_130 = A(1, 6, 0);
                const T eta = A(1, 6, 1);
                const T k_254 = A(1, 7, -1);
                const T alpha = A(1, 7, 0);
                const T s_minus_lo_130 = A(1, 7, 1);
                const T s = A(1, 8, 0);
                c(s - (alpha + k_254 * two_pow(130)));
                c(z_0 - alpha - t_q());
                c(k_254 * (z_130 - two_pow(124)));
                c(k_254 * s_minus_lo_130);
                c((one - k_254) * (one - z_130 * eta) * s_minus_lo_130);
            } else {
                // rows (anchor = halo2 row - 1): [z_0, b_1, b_2], [z_130, alpha, d1_lo], [d1_hi, d2_lo, d2_hi]
                const T z_0 = A(1, 6, -1);
                const T b_1 = A(1, 7, -1);
                const T b_2 = A(1, 8, -1);
                const T z_130 = A(1, 6, 0);
                const T alpha = A(1, 7, 0);
                const T d1_lo = A(1, 8, 0);
                const T d1_hi = A(1, 6, 1);
                const T d2_lo = A(1, 7, 1);
                const T d2_hi = A(1, 8, 1);
                const uint256_t t = var_mul_offset();
                const uint256_t u = var_mul_upper();
                const T k_lo = z_0 - z_130 * two_pow(130);
                c(z_0 - alpha - FF(t));
                c(k_lo - lo130(t) + b_1 * two_pow(130) - d1_lo);
                c(z_130 - hi130(t) - b_1 - d1_hi);
                c(lo130(u) - k_lo + b_2 * two_pow(130) - d2_lo);
                c(hi130(u) - z_130 - b_2 - d2_hi);
                c(bool_check(b_1));
                c(bool_check(b_2));
            }
        });

        // ---- ECC: fixed-base mul (window = a4, u = a5, x_p = a0, y_p = a1, Lagrange coefficients f0..f7, z)
        auto coords_check = [&](auto& c, const T& window) {
            const T x_p = A(0, 0, 0);
            const T y_p = A(0, 1, 0);
            const T u = A(0, 5, 0);
            T interpolated_x = in.fix(F0);
            T window_pow = window;
            for (size_t k = 1; k < 8; ++k) {
                interpolated_x += window_pow * in.fix(F0 + k);
                if (k + 1 < 8) {
                    window_pow *= window;
                }
            }
            c(interpolated_x - x_p);
            c(u.sqr() - y_p - in.fix(FIXED_Z));
            c(curve_eqn(x_p, y_p));
        };
        gate(Q_MUL_FIXED_RUNNING_SUM, [&](auto c) {
            const T word = A(0, 4, 0) - A(0, 4, 1) * FF(8);
            c(range_check(word, 8));
            coords_check(c, word);
        });
        gate(Q_MUL_FIXED_FULL, [&](auto c) {
            const T window = A(0, 4, 0);
            coords_check(c, window);
            c(range_check(window, 8));
        });
        gate(Q_MUL_FIXED_SHORT, [&](auto c) {
            const T y_p = A(0, 1, 0);
            const T y_a = A(0, 3, 0);
            const T last_window = A(0, 5, 0);
            const T sign = A(0, 4, 0);
            c(bool_check(last_window));
            c(sign.sqr() - one);
            c((y_p - y_a) * (y_p + y_a));
            c(sign * y_p - y_a);
        });
        gate(Q_MUL_FIXED_BASE_FIELD, [&](auto c) {
            if constexpr (IS_PASTA) {
                const T alpha = A(1, 6, -1);
                const T z_84_alpha = A(1, 8, -1);
                const T alpha_0 = alpha - z_84_alpha * two_pow(252);
                const T alpha_1 = A(1, 7, 0);
                const T alpha_2 = A(1, 8, 0);
                const T alpha_0_prime = A(1, 6, 0);
                const T z_13_alpha_0_prime = A(1, 6, 1);
                const T z_44_alpha = A(1, 7, 1);
                const T z_43_alpha = A(1, 8, 1);
                c(alpha_2 * alpha_1);
                c(alpha_2 * (z_44_alpha - z_84_alpha * two_pow(120)));
                c(alpha_2 * bool_check(z_43_alpha - z_44_alpha * FF(8)));
                c(alpha_2 * z_13_alpha_0_prime);
                c(range_check(alpha_1, 4));
                c(bool_check(alpha_2));
                c(z_84_alpha - (alpha_1 + alpha_2 * FF(4)));
                c(alpha_0_prime - (alpha_0 + two_pow(130) - t_p()));
            } else {
                const T alpha = A(1, 6, -1);
                const T z_84_alpha = A(1, 8, -1);
                const T alpha_0 = alpha - z_84_alpha * two_pow(252);
                const T alpha_1 = A(1, 7, 0);
                const T alpha_2 = A(1, 8, 0);
                const T alpha_0_prime = A(1, 6, 0);
                const T z_25_alpha_0_prime = A(1, 6, 1);
                c(range_check(alpha_1, 4));
                c(alpha_2);
                c(z_84_alpha - (alpha_1 + alpha_2 * FF(4)));
                c(alpha_0_prime - (alpha_0 + two_pow(250) - t_r()));
                c(is_three(alpha_1) * z_25_alpha_0_prime);
            }
        });

        // ---- Poseidon (state a6..a8, partial sbox a5, rc_a f2..f4, rc_b f5..f7)
        {
            const auto& p = Poseidon::params();
            constexpr std::array<size_t, 3> state = { 6, 7, 8 };
            constexpr std::array<size_t, 3> rc_a = { F2, F3, F4 };
            constexpr std::array<size_t, 3> rc_b = { F5, F6, F7 };
            auto pow5 = [](const T& x) {
                const T x2 = x.sqr();
                return x2.sqr() * x;
            };
            gate(Q_POSEIDON_FULL, [&](auto c) {
                std::array<T, 3> sboxed{ pow5(A(0, state[0], 0) + in.fix(rc_a[0])),
                                         pow5(A(0, state[1], 0) + in.fix(rc_a[1])),
                                         pow5(A(0, state[2], 0) + in.fix(rc_a[2])) };
                for (size_t next_idx = 0; next_idx < 3; ++next_idx) {
                    c(sboxed[0] * p.mds[next_idx][0] + sboxed[1] * p.mds[next_idx][1] + sboxed[2] * p.mds[next_idx][2] -
                      A(0, state[next_idx], 1));
                }
            });
            gate(Q_POSEIDON_PARTIAL, [&](auto c) {
                const T cur_0 = A(0, state[0], 0);
                const T mid_0 = A(0, 5, 0);
                const T s1 = A(0, state[1], 0) + in.fix(rc_a[1]);
                const T s2 = A(0, state[2], 0) + in.fix(rc_a[2]);
                auto mid = [&](size_t i) { return mid_0 * p.mds[i][0] + s1 * p.mds[i][1] + s2 * p.mds[i][2]; };
                std::array<T, 3> next_state{ A(0, state[0], 1), A(0, state[1], 1), A(0, state[2], 1) };
                auto next = [&](size_t i) {
                    return next_state[0] * p.mds_inv[i][0] + next_state[1] * p.mds_inv[i][1] +
                           next_state[2] * p.mds_inv[i][2];
                };
                c(pow5(cur_0 + in.fix(rc_a[0])) - mid_0);
                c(pow5(mid(0) + in.fix(rc_b[0])) - next(0));
                for (size_t i = 1; i < 3; ++i) {
                    c(mid(i) + in.fix(rc_b[i]) - next(i));
                }
            });
            gate(Q_POSEIDON_PAD_AND_ADD, [&](auto c) {
                for (size_t i = 0; i < 2; ++i) {
                    c(A(1, state[i], -1) + A(1, state[i], 0) - A(1, state[i], 1));
                }
                c(A(1, state[2], -1) - A(1, state[2], 1));
            });
        }

        // ---- Sinsemilla (one configuration per chip)
        constexpr std::array<std::array<size_t, 2>, 2> sinsemilla_selectors = { {
            { Q_SINSEMILLA1_1, Q_SINSEMILLA4_1 },
            { Q_SINSEMILLA1_2, Q_SINSEMILLA4_2 },
        } };
        for (size_t chip = 0; chip < 2; ++chip) {
            const auto& cfg = SINSEMILLA_COLUMNS[chip];
            auto x_r = [&](int rot) { return A(0, cfg.lambda_1, rot).sqr() - A(0, cfg.x_a, rot) - A(0, cfg.x_p, rot); };
            auto Y_A = [&](int rot) {
                return (A(0, cfg.lambda_1, rot) + A(0, cfg.lambda_2, rot)) * (A(0, cfg.x_a, rot) - x_r(rot));
            };
            gate(sinsemilla_selectors[chip][0], [&](auto c) {
                const T q_s2 = in.fix(cfg.q_s2);
                const T q_s3 = q_s2 * (q_s2 - one);
                const T lambda_1_next = A(0, cfg.lambda_1, 1);
                const T lambda_2_cur = A(0, cfg.lambda_2, 0);
                const T x_a_cur = A(0, cfg.x_a, 0);
                const T x_a_next = A(0, cfg.x_a, 1);
                c(lambda_2_cur.sqr() - (x_a_next + x_r(0) + x_a_cur));
                const T lhs = lambda_2_cur * FF(4) * (x_a_cur - x_a_next);
                const T rhs = Y_A(0) * two + (two - q_s3) * Y_A(1) + q_s3 * two * lambda_1_next;
                c(lhs - rhs);
            });
            gate(sinsemilla_selectors[chip][1], [&](auto c) { c(in.fix(cfg.fixed_y_q) * two - Y_A(0)); });
        }

        // ---- Merkle chips: decomposition check (advices of the respective Sinsemilla chip)
        constexpr std::array<size_t, 2> merkle_decompose = { Q_MERKLE_DECOMPOSE_1, Q_MERKLE_DECOMPOSE_2 };
        constexpr std::array<size_t, 2> merkle_swap = { Q_SWAP_1, Q_SWAP_2 };
        for (size_t chip = 0; chip < 2; ++chip) {
            const size_t a0 = (chip == 0) ? 0 : 5;
            gate(merkle_decompose[chip], [&](auto c) {
                const T l_whole = A(0, a0 + 4, 1);
                const T a_whole = A(0, a0 + 0, 0);
                const T b_whole = A(0, a0 + 1, 0);
                const T c_whole = A(0, a0 + 2, 0);
                const T left_node = A(0, a0 + 3, 0);
                const T right_node = A(0, a0 + 4, 0);
                const T a_1 = A(0, a0 + 0, 1);
                const T a_0 = a_whole - a_1 * two_pow(10);
                const T z1_b = A(0, a0 + 1, 1);
                const T b_1 = A(0, a0 + 2, 1);
                const T b_2 = A(0, a0 + 3, 1);
                const T b_0 = b_whole - z1_b * two_pow(10);
                c(a_0 - l_whole);
                c(a_1 + (b_0 + b_1 * two_pow(10)) * two_pow(240) - left_node);
                c(b_2 + c_whole * two_pow(5) - right_node);
                c(z1_b - (b_1 + b_2 * two_pow(5)));
            });
            gate(merkle_swap[chip], [&](auto c) {
                const T a = A(0, a0 + 0, 0);
                const T b = A(0, a0 + 1, 0);
                const T a_swapped = A(0, a0 + 2, 0);
                const T b_swapped = A(0, a0 + 3, 0);
                const T swap = A(0, a0 + 4, 0);
                c(a_swapped - ternary(swap, b, a));
                c(b_swapped - ternary(swap, a, b));
                c(bool_check(swap));
            });
        }

        // ---- Lookup range check: short lookup bitshift (running sum column a9)
        gate(Q_BITSHIFT, [&](auto c) {
            const T word = A(1, RANGE_CHECK_COLUMN, -1);
            const T shifted_word = A(1, RANGE_CHECK_COLUMN, 0);
            const T inv_two_pow_s = A(1, RANGE_CHECK_COLUMN, 1);
            c(word * two_pow(10) * inv_two_pow_s - shifted_word);
        });

        // ---- CommitIvk canonicity (a0..a8)
        gate(Q_COMMIT_IVK, [&](auto c) {
            if constexpr (IS_PASTA) {
                const T ak = A(0, 0, 0);
                const T nk = A(0, 0, 1);
                const T a = A(0, 1, 0);
                const T b_whole = A(0, 2, 0);
                const T cc = A(0, 1, 1);
                const T d_whole = A(0, 2, 1);
                const T b_0 = A(0, 3, 0);
                const T b_1 = A(0, 4, 0);
                const T b_2 = A(0, 5, 0);
                const T d_0 = A(0, 3, 1);
                const T d_1 = A(0, 4, 1);
                c(bool_check(b_1));
                c(bool_check(d_1));
                c(b_whole - (b_0 + b_1 * two_pow(4) + b_2 * two_pow(5)));
                c(d_whole - (d_0 + d_1 * two_pow(9)));
                c(a + b_0 * two_pow(250) + b_1 * two_pow(254) - ak);
                c(b_2 + cc * two_pow(5) + d_0 * two_pow(245) + d_1 * two_pow(254) - nk);
                // ak canonicity
                c(b_1 * b_0);
                c(b_1 * A(0, 6, 0));
                c(a + two_pow(130) - t_p() - A(0, 7, 0));
                c(b_1 * A(0, 8, 0));
                // nk canonicity
                c(d_1 * d_0);
                c(d_1 * A(0, 6, 1));
                c(b_2 + cc * two_pow(5) + two_pow(140) - t_p() - A(0, 7, 1));
                c(d_1 * A(0, 8, 1));
            } else {
                const T ak = A(0, 0, 0);
                const T nk = A(0, 0, 1);
                const T a = A(0, 1, 0);
                const T b_whole = A(0, 2, 0);
                const T cc = A(0, 1, 1);
                const T d_whole = A(0, 2, 1);
                const T b_0 = A(0, 3, 0);
                const T b_1 = A(0, 4, 0);
                const T b_2 = A(0, 5, 0);
                const T d_0 = A(0, 3, 1);
                const T d_1 = A(0, 4, 1);
                const T d_lo = A(0, 5, 1);
                // ak: t_hi = bits 252..253, b_0 - 4 t_hi = bits 250..251
                const T t_hi_ak = A(0, 6, 0);
                const T a_prime = A(0, 7, 0);
                const T z25_a_prime = A(0, 8, 0);
                const T t_lo_ak = b_0 - t_hi_ak * FF(4);
                // nk: d_0 = d_lo + 2^5 (t_lo + 4 t_hi)
                const T t_hi_nk = A(0, 6, 1);
                const T l_prime = A(0, 7, 1);
                const T z25_l_prime = A(0, 8, 1);
                const T t_lo_nk = A(0, 9, 1);
                c(b_whole - (b_0 + b_1 * two_pow(4) + b_2 * two_pow(5)));
                c(d_whole - (d_0 + d_1 * two_pow(9)));
                c(a + b_0 * two_pow(250) + b_1 * two_pow(254) - ak);
                c(b_2 + cc * two_pow(5) + d_0 * two_pow(245) + d_1 * two_pow(254) - nk);
                c(b_1);
                c(d_1);
                c(range_check(t_hi_ak, 4));
                c(range_check(t_lo_ak, 4));
                c(is_three(t_hi_ak) * t_lo_ak);
                c(is_three(t_hi_ak) * z25_a_prime);
                c(a + two_pow(250) - t_r() - a_prime);
                c(range_check(t_hi_nk, 4));
                c(range_check(t_lo_nk, 4));
                c(d_0 - (d_lo + t_lo_nk * two_pow(5) + t_hi_nk * two_pow(7)));
                c(is_three(t_hi_nk) * t_lo_nk);
                c(is_three(t_hi_nk) * z25_l_prime);
                c(b_2 + cc * two_pow(5) + d_lo * two_pow(245) + two_pow(250) - t_r() - l_prime);
            }
        });

        // ---- NoteCommit (col_l = a6, col_m = a7, col_r = a8, col_z = a9)
        constexpr size_t L = 6;
        constexpr size_t M = 7;
        constexpr size_t R = 8;
        constexpr size_t Z = 9;
        gate(Q_NOTECOMMIT_B, [&](auto c) {
            const T b = A(0, L, 0);
            const T b_0 = A(0, M, 0);
            const T b_1 = A(0, R, 0);
            const T b_2 = A(0, M, 1);
            const T b_3 = A(0, R, 1);
            c(bool_check(b_1));
            c(bool_check(b_2));
            c(b - (b_0 + b_1 * two_pow(4) + b_2 * two_pow(5) + b_3 * two_pow(6)));
        });
        gate(Q_NOTECOMMIT_D, [&](auto c) {
            const T d = A(0, L, 0);
            const T d_0 = A(0, M, 0);
            const T d_1 = A(0, R, 0);
            const T d_2 = A(0, M, 1);
            const T d_3 = A(0, R, 1);
            c(bool_check(d_0));
            c(bool_check(d_1));
            c(d - (d_0 + d_1 * two + d_2 * two_pow(2) + d_3 * two_pow(10)));
        });
        gate(Q_NOTECOMMIT_E, [&](auto c) { c(A(0, L, 0) - (A(0, M, 0) + A(0, R, 0) * two_pow(6))); });
        gate(Q_NOTECOMMIT_G, [&](auto c) {
            const T g = A(0, L, 0);
            const T g_0 = A(0, M, 0);
            const T g_1 = A(0, L, 1);
            const T g_2 = A(0, M, 1);
            c(bool_check(g_0));
            c(g - (g_0 + g_1 * two + g_2 * two_pow(10)));
        });
        gate(Q_NOTECOMMIT_H, [&](auto c) {
            const T h = A(0, L, 0);
            const T h_0 = A(0, M, 0);
            const T h_1 = A(0, R, 0);
            c(bool_check(h_1));
            c(h - (h_0 + h_1 * two_pow(5)));
        });
        gate(Q_NOTECOMMIT_G_D, [&](auto c) {
            if constexpr (IS_PASTA) {
                const T gd_x = A(0, L, 0);
                const T b_0 = A(0, M, 0);
                const T b_1 = A(0, M, 1);
                const T a = A(0, R, 0);
                const T a_prime = A(0, R, 1);
                const T z13_a = A(0, Z, 0);
                const T z13_a_prime = A(0, Z, 1);
                c(a + b_0 * two_pow(250) + b_1 * two_pow(254) - gd_x);
                c(a + two_pow(130) - t_p() - a_prime);
                c(b_1 * b_0);
                c(b_1 * z13_a);
                c(b_1 * z13_a_prime);
            } else {
                const T gd_x = A(0, L, 0);
                const T b_0 = A(0, M, 0);
                const T b_1 = A(0, M, 1);
                const T a = A(0, R, 0);
                const T a_prime = A(0, R, 1);
                const T t_hi = A(0, Z, 0);
                const T z25_a_prime = A(0, Z, 1);
                const T t_lo = b_0 - t_hi * FF(4);
                c(a + b_0 * two_pow(250) + b_1 * two_pow(254) - gd_x);
                c(a + two_pow(250) - t_r() - a_prime);
                c(b_1);
                c(range_check(t_hi, 4));
                c(range_check(t_lo, 4));
                c(is_three(t_hi) * t_lo);
                c(is_three(t_hi) * z25_a_prime);
            }
        });
        gate(Q_NOTECOMMIT_PK_D, [&](auto c) {
            if constexpr (IS_PASTA) {
                const T pkd_x = A(0, L, 0);
                const T b_3 = A(0, M, 0);
                const T d_0 = A(0, M, 1);
                const T cc = A(0, R, 0);
                const T b3_c_prime = A(0, R, 1);
                const T z13_c = A(0, Z, 0);
                const T z14_b3_c_prime = A(0, Z, 1);
                c(b_3 + cc * two_pow(4) + d_0 * two_pow(254) - pkd_x);
                c(b_3 + cc * two_pow(4) + two_pow(140) - t_p() - b3_c_prime);
                c(d_0 * z13_c);
                c(d_0 * z14_b3_c_prime);
            } else {
                // rows: [x, low, mid, z24_mid], [t_lo, top, l_prime, z25_l_prime], [t_hi, w_lo]
                const T x = A(0, L, 0);
                const T low = A(0, M, 0);
                const T mid = A(0, R, 0);
                const T z24_mid = A(0, Z, 0);
                const T t_lo = A(0, L, 1);
                const T top = A(0, M, 1);
                const T l_prime = A(0, R, 1);
                const T z25_l_prime = A(0, Z, 1);
                const T t_hi = A(0, L, 2);
                const T w_lo = A(0, M, 2);
                c(low + mid * two_pow(4) + top * two_pow(254) - x);
                c(top);
                c(z24_mid - (w_lo + t_lo * two_pow(6) + t_hi * two_pow(8)));
                c(range_check(t_lo, 4));
                c(range_check(t_hi, 4));
                c(is_three(t_hi) * t_lo);
                c(is_three(t_hi) * z25_l_prime);
                c(low + mid * two_pow(4) - (t_lo + t_hi * FF(4)) * two_pow(250) + two_pow(250) - t_r() - l_prime);
            }
        });
        gate(Q_NOTECOMMIT_VALUE, [&](auto c) {
            const T value = A(0, L, 0);
            const T d_2 = A(0, M, 0);
            const T d_3 = A(0, R, 0);
            const T e_0 = A(0, Z, 0);
            c(d_2 + d_3 * two_pow(8) + e_0 * two_pow(58) - value);
        });
        gate(Q_NOTECOMMIT_RHO, [&](auto c) {
            if constexpr (IS_PASTA) {
                const T rho = A(0, L, 0);
                const T e_1 = A(0, M, 0);
                const T g_0 = A(0, M, 1);
                const T f = A(0, R, 0);
                const T e1_f_prime = A(0, R, 1);
                const T z13_f = A(0, Z, 0);
                const T z14_e1_f_prime = A(0, Z, 1);
                c(e_1 + f * two_pow(4) + g_0 * two_pow(254) - rho);
                c(e_1 + f * two_pow(4) + two_pow(140) - t_p() - e1_f_prime);
                c(g_0 * z13_f);
                c(g_0 * z14_e1_f_prime);
            } else {
                // rows: [x, low, mid, z24_mid], [t_lo, top, l_prime, z25_l_prime], [t_hi, w_lo]
                const T x = A(0, L, 0);
                const T low = A(0, M, 0);
                const T mid = A(0, R, 0);
                const T z24_mid = A(0, Z, 0);
                const T t_lo = A(0, L, 1);
                const T top = A(0, M, 1);
                const T l_prime = A(0, R, 1);
                const T z25_l_prime = A(0, Z, 1);
                const T t_hi = A(0, L, 2);
                const T w_lo = A(0, M, 2);
                c(low + mid * two_pow(4) + top * two_pow(254) - x);
                c(top);
                c(z24_mid - (w_lo + t_lo * two_pow(6) + t_hi * two_pow(8)));
                c(range_check(t_lo, 4));
                c(range_check(t_hi, 4));
                c(is_three(t_hi) * t_lo);
                c(is_three(t_hi) * z25_l_prime);
                c(low + mid * two_pow(4) - (t_lo + t_hi * FF(4)) * two_pow(250) + two_pow(250) - t_r() - l_prime);
            }
        });
        gate(Q_NOTECOMMIT_PSI, [&](auto c) {
            if constexpr (IS_PASTA) {
                const T psi = A(0, L, 0);
                const T h_0 = A(0, L, 1);
                const T g_1 = A(0, M, 0);
                const T h_1 = A(0, M, 1);
                const T g_2 = A(0, R, 0);
                const T g1_g2_prime = A(0, R, 1);
                const T z13_g = A(0, Z, 0);
                const T z13_g1_g2_prime = A(0, Z, 1);
                c(g_1 + g_2 * two_pow(9) + h_0 * two_pow(249) + h_1 * two_pow(254) - psi);
                c(g_1 + g_2 * two_pow(9) + two_pow(130) - t_p() - g1_g2_prime);
                c(h_1 * h_0);
                c(h_1 * z13_g);
                c(h_1 * z13_g1_g2_prime);
            } else {
                // rows: [psi, g_1, g_2, t_hi], [h_0, h_1, l_prime, z25_l_prime], [h_lsb, t_lo]
                const T psi = A(0, L, 0);
                const T g_1 = A(0, M, 0);
                const T g_2 = A(0, R, 0);
                const T t_hi = A(0, Z, 0);
                const T h_0 = A(0, L, 1);
                const T h_1 = A(0, M, 1);
                const T l_prime = A(0, R, 1);
                const T z25_l_prime = A(0, Z, 1);
                const T h_lsb = A(0, L, 2);
                const T t_lo = A(0, M, 2);
                c(g_1 + g_2 * two_pow(9) + h_0 * two_pow(249) + h_1 * two_pow(254) - psi);
                c(h_1);
                c(bool_check(h_lsb));
                c(range_check(t_lo, 4));
                c(range_check(t_hi, 4));
                c(h_0 - (h_lsb + t_lo * FF(2) + t_hi * FF(8)));
                c(is_three(t_hi) * t_lo);
                c(is_three(t_hi) * z25_l_prime);
                c(g_1 + g_2 * two_pow(9) + h_lsb * two_pow(249) + two_pow(250) - t_r() - l_prime);
            }
        });
        gate(Q_Y_CANON, [&](auto c) {
            if constexpr (IS_PASTA) {
                const T y = A(0, 5, 0);
                const T lsb = A(0, 6, 0);
                const T k_0 = A(0, 7, 0);
                const T k_2 = A(0, 8, 0);
                const T k_3 = A(0, 9, 0);
                const T j = A(0, 5, 1);
                const T z1_j = A(0, 6, 1);
                const T z13_j = A(0, 7, 1);
                const T j_prime = A(0, 8, 1);
                const T z13_j_prime = A(0, 9, 1);
                c(bool_check(k_3));
                c(j - (lsb + k_0 * two + z1_j * two_pow(10)));
                c(y - (j + k_2 * two_pow(250) + k_3 * two_pow(254)));
                c(j + two_pow(130) - t_p() - j_prime);
                c(k_3 * k_2);
                c(k_3 * z13_j);
                c(k_3 * z13_j_prime);
            } else {
                const T y = A(0, 5, 0);
                const T lsb = A(0, 6, 0);
                const T k_0 = A(0, 7, 0);
                const T k_2 = A(0, 8, 0);
                const T k_3 = A(0, 9, 0);
                const T j = A(0, 5, 1);
                const T z1_j = A(0, 6, 1);
                const T t_hi = A(0, 7, 1);
                const T j_prime = A(0, 8, 1);
                const T z25_j_prime = A(0, 9, 1);
                const T t_lo = k_2 - t_hi * FF(4);
                c(k_3);
                c(j - (lsb + k_0 * two + z1_j * two_pow(10)));
                c(y - (j + k_2 * two_pow(250) + k_3 * two_pow(254)));
                c(j + two_pow(250) - t_r() - j_prime);
                c(range_check(t_hi, 4));
                c(range_check(t_lo, 4));
                c(is_three(t_hi) * t_lo);
                c(is_three(t_hi) * z25_j_prime);
            }
        });

        BB_ASSERT_EQ(emitted, NUM_GATE_CONSTRAINTS);
    }
};

} // namespace bb::zcash::halo2
