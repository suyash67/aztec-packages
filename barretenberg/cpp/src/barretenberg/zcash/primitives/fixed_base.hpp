#pragma once

#include "barretenberg/common/assert.hpp"
#include "barretenberg/common/thread.hpp"
#include "cycle.hpp"

#include <array>
#include <span>
#include <vector>

namespace bb::zcash {

/**
 * @brief Precomputed data for halo2's windowed fixed-base scalar multiplication (halo2_gadgets ecc/chip/constants.rs).
 * @details The scalar is split into 3-bit windows k_w. Window w < W-1 uses the points
 *   M_w(k) = [(k + 2) * 8^w] B,   k in [0, 8),
 * and the last window uses M_{W-1}(k) = [k * 8^{W-1} - Σ_{j<W-1} 2^{3j+1}] B, so that the offsets cancel. In-circuit,
 * x(M_w(k)) is the evaluation at k of the degree-7 interpolant with coefficients `lagrange_coeffs[w]`, and
 * y(M_w(k)) is pinned by u_w(k)^2 = y + z_w, where z_w is chosen so that -y + z_w is a non-square for every k.
 */
template <typename Cycle> struct FixedBase {
    using FF = typename Cycle::FF;
    using AffineElement = typename Cycle::AffineElement;
    using Element = typename Cycle::Element;
    using Scalar = typename Cycle::EmbeddedScalar;

    static constexpr size_t WINDOW_BITS = 3;
    static constexpr size_t H = size_t{ 1 } << WINDOW_BITS;

    AffineElement generator;
    size_t num_windows = 0;
    std::vector<std::array<AffineElement, H>> window_table;
    std::vector<std::array<FF, H>> lagrange_coeffs;
    std::vector<uint64_t> z;
    std::vector<std::array<FF, H>> u;

    static Scalar last_window_offset(size_t num_windows)
    {
        Scalar sum = 0;
        for (size_t j = 0; j + 1 < num_windows; ++j) {
            sum += Scalar(uint256_t(1) << (WINDOW_BITS * j + 1));
        }
        return sum;
    }

    static std::vector<std::array<AffineElement, H>> compute_window_table(const AffineElement& base, size_t num_windows)
    {
        std::vector<std::array<AffineElement, H>> table(num_windows);
        const Scalar eight = Scalar(H);
        const Scalar offset = last_window_offset(num_windows);
        parallel_for(num_windows, [&](size_t w) {
            const Scalar eight_pow_w = eight.pow(w);
            for (size_t k = 0; k < H; ++k) {
                Scalar s = (w + 1 < num_windows) ? Scalar(k + 2) * eight_pow_w : Scalar(k) * eight_pow_w - offset;
                table[w][k] = AffineElement(Element(base) * s);
            }
        });
        return table;
    }

    // Coefficients c_0..c_7 (monomial basis) of the polynomial through (k, xs[k]), k = 0..7.
    static std::array<FF, H> interpolate(const std::array<FF, H>& xs)
    {
        std::array<FF, H> coeffs{};
        for (size_t j = 0; j < H; ++j) {
            // basis polynomial L_j(X) = Π_{m≠j} (X - m) / (j - m), accumulated in monomial form
            std::array<FF, H> basis{};
            basis[0] = 1;
            size_t deg = 0;
            FF denom = 1;
            for (size_t m = 0; m < H; ++m) {
                if (m == j) {
                    continue;
                }
                // multiply basis by (X - m)
                for (size_t d = deg + 1; d > 0; --d) {
                    basis[d] = basis[d - 1] - basis[d] * FF(m);
                }
                basis[0] = -basis[0] * FF(m);
                ++deg;
                denom *= FF(j) - FF(m);
            }
            const FF scale = xs[j] * denom.invert();
            for (size_t d = 0; d < H; ++d) {
                coeffs[d] += basis[d] * scale;
            }
        }
        return coeffs;
    }

    static FixedBase make(const AffineElement& generator, size_t num_windows, std::span<const uint64_t> zs)
    {
        BB_ASSERT_EQ(zs.size(), num_windows);
        FixedBase fb;
        fb.generator = generator;
        fb.num_windows = num_windows;
        fb.window_table = compute_window_table(generator, num_windows);
        fb.lagrange_coeffs.resize(num_windows);
        fb.u.resize(num_windows);
        fb.z.assign(zs.begin(), zs.end());
        parallel_for(num_windows, [&](size_t w) {
            std::array<FF, H> xs{};
            for (size_t k = 0; k < H; ++k) {
                xs[k] = fb.window_table[w][k].x;
            }
            fb.lagrange_coeffs[w] = interpolate(xs);
            for (size_t k = 0; k < H; ++k) {
                const FF y = fb.window_table[w][k].y;
                auto [is_square, root] = (y + FF(fb.z[w])).sqrt();
                BB_ASSERT(is_square, "fixed-base z value does not make y + z a square");
                fb.u[w][k] = root;
            }
        });
        return fb;
    }

    /**
     * @brief Smallest z such that, for every window point, y + z is a square and -y + z is not (find_zs_and_us).
     * @details Expected ~2^16 candidates per window; used to derive z-values for curves without published constants.
     */
    static std::vector<uint64_t> find_zs(const AffineElement& generator, size_t num_windows)
    {
        const auto table = compute_window_table(generator, num_windows);
        std::vector<uint64_t> zs(num_windows);
        parallel_for(num_windows, [&](size_t w) {
            for (uint64_t z = 0;; ++z) {
                bool ok = true;
                for (size_t k = 0; k < H && ok; ++k) {
                    const FF y = table[w][k].y;
                    ok = !is_quadratic_residue(-y + FF(z)) && is_quadratic_residue(y + FF(z));
                }
                if (ok) {
                    zs[w] = z;
                    return;
                }
            }
        });
        return zs;
    }

    static bool is_quadratic_residue(const FF& x)
    {
        if (x.is_zero()) {
            return true;
        }
        return x.pow((FF::modulus - 1) >> 1) == FF::one();
    }
};

} // namespace bb::zcash
