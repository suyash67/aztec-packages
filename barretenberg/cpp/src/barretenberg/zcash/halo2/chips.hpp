#pragma once

#include "barretenberg/zcash/halo2/layout.hpp"
#include "barretenberg/zcash/halo2/plonkish_builder.hpp"
#include "barretenberg/zcash/primitives/orchard.hpp"

#include <array>
#include <optional>
#include <vector>

/**
 * @file chips.hpp
 * @brief Witness generation for the halo2_gadgets 0.6 chips used by the Orchard Action circuit.
 *
 * @details Each function is a port of the corresponding halo2 chip's `assign`/`assign_region`: the same regions, the
 * same cell positions (relative to the region), the same selectors and the same copy constraints. Gate polynomials live
 * in `relations.hpp`. Function-level comments name the halo2 source.
 */
namespace bb::zcash::halo2 {

using namespace layout;

template <typename FF> FF inv0(const FF& x)
{
    return x.is_zero() ? FF(0) : x.invert();
}

// halo2_gadgets utilities::bitrange_subset
template <typename FF> FF bitrange_subset(const FF& value, size_t start, size_t end)
{
    const uint256_t v(value);
    return FF(v.slice(start, end));
}

// halo2_gadgets utilities::decompose_word: little-endian windows of `window_bits` bits over the first
// `word_num_bits` bits (zero-padded to a multiple of the window size).
template <typename FF> std::vector<uint64_t> decompose_word(const FF& word, size_t word_num_bits, size_t window_bits)
{
    const uint256_t v(word);
    const size_t num_windows = (word_num_bits + window_bits - 1) / window_bits;
    std::vector<uint64_t> out(num_windows);
    for (size_t i = 0; i < num_windows; ++i) {
        uint64_t w = 0;
        for (size_t b = 0; b < window_bits; ++b) {
            const size_t bit = i * window_bits + b;
            if (bit < word_num_bits && v.get_bit(bit)) {
                w |= uint64_t{ 1 } << b;
            }
        }
        out[i] = w;
    }
    return out;
}

/**
 * @brief All chips of the Action circuit over one builder.
 */
template <typename Cycle> class Chips {
  public:
    using FF = typename Cycle::FF;
    using Scalar = typename Cycle::EmbeddedScalar;
    using AffineElement = typename Cycle::AffineElement;
    using Element = typename Cycle::Element;
    using Builder = PlonkishBuilder<FF>;
    using Region = typename Builder::Region;
    using Cell = typename Builder::Cell;
    using O = Orchard<Cycle>;
    using FB = FixedBase<Cycle>;

    static constexpr size_t K = 10; // Sinsemilla / lookup word size
    static constexpr size_t H = 8;  // fixed-base window table size
    static constexpr size_t FIXED_BASE_WINDOW_SIZE = 3;
    static constexpr size_t L_SCALAR_SHORT = 64;
    static constexpr size_t NUM_WINDOWS = 85;
    static constexpr size_t NUM_WINDOWS_SHORT = 22;
    static constexpr size_t SCALAR_NUM_BITS = 255;
    static constexpr size_t BASE_NUM_BITS = 255;
    static constexpr size_t NUM_COMPLETE_BITS = 3;
    static constexpr size_t INCOMPLETE_LEN = SCALAR_NUM_BITS - 1 - NUM_COMPLETE_BITS;
    static constexpr size_t INCOMPLETE_HI_LEN = INCOMPLETE_LEN / 2;
    static constexpr size_t INCOMPLETE_LO_LEN = INCOMPLETE_LEN - INCOMPLETE_HI_LEN;

    // T_P = p - 2^254 and T_Q = q - 2^254 for the Pasta moduli (circuit field and embedded-curve scalar field).
    static uint256_t t_p() { return uint256_t(FF::modulus) - (uint256_t(1) << 254); }
    static uint256_t t_q() { return uint256_t(Scalar::modulus) - (uint256_t(1) << 254); }
    static FF two_pow(size_t k) { return FF(uint256_t(1) << k); }

    struct Point {
        Cell x;
        Cell y;
        AffineElement value() const
        {
            if (x.value.is_zero() && y.value.is_zero()) {
                return AffineElement::infinity();
            }
            return AffineElement(x.value, y.value);
        }
    };

    // A cell together with the bit length it is known to fit in (halo2's RangeConstrained<AssignedCell>).
    struct RangeConstrainedCell {
        Cell cell;
        size_t num_bits;
    };
    // An unassigned value with a bit length (halo2's RangeConstrained<Value>).
    struct RangeConstrainedValue {
        FF value;
        size_t num_bits;
    };

    struct MessagePiece {
        Cell cell;
        size_t num_words;
    };

    struct SinsemillaConfig {
        size_t x_a, x_p, bits, lambda_1, lambda_2; // advice columns
        size_t witness_pieces;                     // advice column
        size_t fixed_y_q;                          // fixed column
        size_t q_s2;                               // fixed column q_sinsemilla2
        size_t q_s1, q_s4, q_decompose, q_swap;    // selectors (Sinsemilla and the Merkle chip on top of it)
        std::array<size_t, 5> advices() const { return { x_a, x_p, bits, lambda_1, lambda_2 }; }
    };

    static SinsemillaConfig sinsemilla_config_1()
    {
        return {
            0, 1, 2, 3, 4, 6, F0, Q_SINSEMILLA2_1, Q_SINSEMILLA1_1, Q_SINSEMILLA4_1, Q_MERKLE_DECOMPOSE_1, Q_SWAP_1
        };
    }
    static SinsemillaConfig sinsemilla_config_2()
    {
        return {
            5, 6, 7, 8, 9, 7, F1, Q_SINSEMILLA2_2, Q_SINSEMILLA1_2, Q_SINSEMILLA4_2, Q_MERKLE_DECOMPOSE_2, Q_SWAP_2
        };
    }

    explicit Chips(Builder& builder, bool anchored_base = true)
        : b_(builder)
        , anchored_base_(anchored_base)
    {}

    Builder& builder() { return b_; }

    // ---------------------------------------------------------------- utilities

    // orchard circuit/gadget.rs assign_free_advice
    Cell assign_free_advice(size_t column, const FF& value)
    {
        return b_.assign_region("load private", [&](Region& r) { return r.assign_advice(advice(column), 0, value); });
    }

    // ---------------------------------------------------------------- lookup range check (running sum in a9)
    static constexpr size_t RANGE_CHECK_COLUMN = 9;

    std::vector<Cell> range_check(Region& r, const Cell& z0, size_t num_words, bool strict)
    {
        const auto words = decompose_word(z0.value, num_words * K, K);
        const FF inv_two_pow_k = FF(uint64_t{ 1 } << K).invert();
        std::vector<Cell> zs{ z0 };
        Cell z = z0;
        for (size_t idx = 0; idx < num_words; ++idx) {
            r.enable_selector(selector(Q_LOOKUP), idx);
            r.enable_selector(selector(Q_RUNNING), idx);
            const FF z_val = (z.value - FF(words[idx])) * inv_two_pow_k;
            z = r.assign_advice(advice(RANGE_CHECK_COLUMN), idx + 1, z_val);
            zs.push_back(z);
        }
        if (strict) {
            b_.constrain_constant(zs.back(), FF(0));
        }
        return zs;
    }

    std::vector<Cell> copy_check(const Cell& element, size_t num_words, bool strict)
    {
        return b_.assign_region("words range check", [&](Region& r) {
            Cell z0 = r.copy_advice(element, advice(RANGE_CHECK_COLUMN), 0);
            return range_check(r, z0, num_words, strict);
        });
    }

    std::vector<Cell> witness_check(const FF& value, size_t num_words, bool strict)
    {
        return b_.assign_region("Witness element", [&](Region& r) {
            Cell z0 = r.assign_advice(advice(RANGE_CHECK_COLUMN), 0, value);
            return range_check(r, z0, num_words, strict);
        });
    }

    void short_range_check(Region& r, const Cell& element, size_t num_bits)
    {
        r.enable_selector(selector(Q_LOOKUP), 0);
        r.enable_selector(selector(Q_LOOKUP), 1);
        r.enable_selector(selector(Q_BITSHIFT), 1);
        r.assign_advice(advice(RANGE_CHECK_COLUMN), 1, element.value * FF(uint64_t{ 1 } << (K - num_bits)));
        r.assign_advice_from_constant(advice(RANGE_CHECK_COLUMN), 2, FF(uint64_t{ 1 } << num_bits).invert());
    }

    Cell witness_short_check(const FF& element, size_t num_bits)
    {
        BB_ASSERT_LTE(num_bits, K);
        return b_.assign_region("Range check bits", [&](Region& r) {
            Cell e = r.assign_advice(advice(RANGE_CHECK_COLUMN), 0, element);
            short_range_check(r, e, num_bits);
            return e;
        });
    }

    // RangeConstrained::witness_short
    RangeConstrainedCell witness_short(const FF& value, size_t start, size_t end)
    {
        BB_ASSERT_LT(end - start, K);
        return { witness_short_check(bitrange_subset(value, start, end), end - start), end - start };
    }

    static RangeConstrainedValue bitrange_of(const FF& value, size_t start, size_t end)
    {
        return { bitrange_subset(value, start, end), end - start };
    }

    // ---------------------------------------------------------------- Sinsemilla (sinsemilla/chip.rs)

    MessagePiece witness_message_piece(const SinsemillaConfig& cfg, const FF& value, size_t num_words)
    {
        Cell c = b_.assign_region("witness message piece",
                                  [&](Region& r) { return r.assign_advice(advice(cfg.witness_pieces), 0, value); });
        return { c, num_words };
    }

    // MessagePiece::from_subpieces
    MessagePiece from_subpieces(const SinsemillaConfig& cfg, const std::vector<RangeConstrainedValue>& subpieces)
    {
        FF acc = 0;
        size_t bits = 0;
        for (const auto& s : subpieces) {
            acc += s.value * two_pow(bits);
            bits += s.num_bits;
        }
        BB_ASSERT_EQ(bits % K, size_t{ 0 });
        return witness_message_piece(cfg, acc, bits / K);
    }

    struct HashOutput {
        Point point;
        std::vector<std::vector<Cell>> zs; // running sums per piece
    };

    // SinsemillaChip::hash_message with a public Q (allow_init_from_private_point = false)
    HashOutput hash_to_point(const SinsemillaConfig& cfg,
                             const AffineElement& Q,
                             const std::vector<MessagePiece>& pieces)
    {
        return b_.assign_region("hash_to_point", [&](Region& r) {
            size_t offset = 0;
            // public_q_initialization
            r.enable_selector(selector(cfg.q_s4), offset);
            r.assign_fixed(fixed(cfg.fixed_y_q), offset, Q.y);
            Cell x_a = r.assign_advice_from_constant(advice(cfg.x_a), offset, Q.x);
            FF y_a = Q.y;

            const auto& table = Sinsemilla<Cycle>::S_table();
            const FF inv_2_k = FF(uint64_t{ 1 } << K).invert();
            std::vector<std::vector<Cell>> zs_sum;
            for (size_t idx = 0; idx < pieces.size(); ++idx) {
                const auto& piece = pieces[idx];
                const bool final_piece = idx + 1 == pieces.size();
                for (size_t row = 0; row < piece.num_words; ++row) {
                    r.enable_selector(selector(cfg.q_s1), offset + row);
                }
                for (size_t row = 0; row + 1 < piece.num_words; ++row) {
                    r.assign_fixed(fixed(cfg.q_s2), offset + row, FF(1));
                }
                r.assign_fixed(fixed(cfg.q_s2), offset + piece.num_words - 1, FF(final_piece ? 2 : 0));

                const auto words = decompose_word(piece.cell.value, K * piece.num_words, K);
                std::vector<Cell> zs;
                zs.push_back(r.copy_advice(piece.cell, advice(cfg.bits), offset));
                FF z = piece.cell.value;
                for (size_t i = 0; i + 1 < words.size(); ++i) {
                    z = (z - FF(words[i])) * inv_2_k;
                    zs.push_back(r.assign_advice(advice(cfg.bits), offset + i + 1, z));
                }
                for (size_t row = 0; row < words.size(); ++row) {
                    const auto& gen = table[words[row]];
                    const FF x_p = gen.x;
                    const FF y_p = gen.y;
                    r.assign_advice(advice(cfg.x_p), offset + row, x_p);
                    const FF lambda_1 = (y_a - y_p) * (x_a.value - x_p).invert();
                    r.assign_advice(advice(cfg.lambda_1), offset + row, lambda_1);
                    const FF x_r = lambda_1.sqr() - x_a.value - x_p;
                    const FF lambda_2 = y_a * FF(2) * (x_a.value - x_r).invert() - lambda_1;
                    r.assign_advice(advice(cfg.lambda_2), offset + row, lambda_2);
                    const FF x_a_new = lambda_2.sqr() - x_a.value - x_r;
                    Cell x_a_cell = r.assign_advice(advice(cfg.x_a), offset + row + 1, x_a_new);
                    y_a = lambda_2 * (x_a.value - x_a_new) - y_a;
                    x_a = x_a_cell;
                }
                offset += piece.num_words;
                zs_sum.push_back(std::move(zs));
            }
            Cell y_a_cell = r.assign_advice(advice(cfg.lambda_1), offset, y_a);
            r.assign_advice(advice(cfg.lambda_2), offset, FF(0));
            r.assign_advice(advice(cfg.x_p), offset, FF(0));
            BB_ASSERT(!x_a.value.is_zero() && !y_a.is_zero());
            return HashOutput{ Point{ x_a, y_a_cell }, std::move(zs_sum) };
        });
    }

    // ---------------------------------------------------------------- ECC: witness_point.rs

    Point witness_point(const AffineElement& p)
    {
        return b_.assign_region("witness point", [&](Region& r) {
            r.enable_selector(selector(Q_POINT), 0);
            const bool id = p.is_point_at_infinity();
            Cell x = r.assign_advice(advice(0), 0, id ? FF(0) : p.x);
            Cell y = r.assign_advice(advice(1), 0, id ? FF(0) : p.y);
            return Point{ x, y };
        });
    }

    Point witness_point_non_id(const AffineElement& p)
    {
        BB_ASSERT(!p.is_point_at_infinity());
        return b_.assign_region("witness non-identity point", [&](Region& r) {
            r.enable_selector(selector(Q_POINT_NON_ID), 0);
            Cell x = r.assign_advice(advice(0), 0, p.x);
            Cell y = r.assign_advice(advice(1), 0, p.y);
            return Point{ x, y };
        });
    }

    void constrain_equal(const Point& a, const Point& b)
    {
        b_.constrain_equal(a.x, b.x);
        b_.constrain_equal(a.y, b.y);
    }

    // ---------------------------------------------------------------- ECC: add_incomplete.rs (a0, a1, a2, a3)

    Point add_incomplete_region(Region& r, const Point& p, const Point& q, size_t offset)
    {
        r.enable_selector(selector(Q_ADD_INCOMPLETE), offset);
        BB_ASSERT(p.x.value != q.x.value, "incomplete addition exceptional case");
        r.copy_advice(p.x, advice(0), offset);
        r.copy_advice(p.y, advice(1), offset);
        r.copy_advice(q.x, advice(2), offset);
        r.copy_advice(q.y, advice(3), offset);
        const FF lambda = (q.y.value - p.y.value) * (q.x.value - p.x.value).invert();
        const FF x_r = lambda.sqr() - p.x.value - q.x.value;
        const FF y_r = lambda * (p.x.value - x_r) - p.y.value;
        Cell xr = r.assign_advice(advice(2), offset + 1, x_r);
        Cell yr = r.assign_advice(advice(3), offset + 1, y_r);
        return Point{ xr, yr };
    }

    // ---------------------------------------------------------------- ECC: add.rs (complete addition, a0..a8)

    Point add_region(Region& r, const Point& p, const Point& q, size_t offset)
    {
        r.enable_selector(selector(Q_ADD), offset);
        r.copy_advice(p.x, advice(0), offset);
        r.copy_advice(p.y, advice(1), offset);
        r.copy_advice(q.x, advice(2), offset);
        r.copy_advice(q.y, advice(3), offset);
        const FF x_p = p.x.value;
        const FF y_p = p.y.value;
        const FF x_q = q.x.value;
        const FF y_q = q.y.value;
        const FF alpha = inv0(x_q - x_p);
        r.assign_advice(advice(5), offset, alpha);
        r.assign_advice(advice(6), offset, inv0(x_p));
        r.assign_advice(advice(7), offset, inv0(x_q));
        r.assign_advice(advice(8), offset, x_q == x_p ? inv0(y_q + y_p) : FF(0));
        FF lambda = 0;
        if (x_q != x_p) {
            lambda = (y_q - y_p) * alpha;
        } else if (!y_p.is_zero()) {
            lambda = x_p.sqr() * FF(3) * (y_p * FF(2)).invert();
        }
        r.assign_advice(advice(4), offset, lambda);
        FF x_r;
        FF y_r;
        if (x_p.is_zero()) {
            x_r = x_q;
            y_r = y_q;
        } else if (x_q.is_zero()) {
            x_r = x_p;
            y_r = y_p;
        } else if (x_q == x_p && y_q == -y_p) {
            x_r = 0;
            y_r = 0;
        } else {
            x_r = lambda.sqr() - x_p - x_q;
            y_r = lambda * (x_p - x_r) - y_p;
        }
        Cell xr = r.assign_advice(advice(2), offset + 1, x_r);
        Cell yr = r.assign_advice(advice(3), offset + 1, y_r);
        return Point{ xr, yr };
    }

    Point add(const Point& p, const Point& q)
    {
        return b_.assign_region("complete point addition", [&](Region& r) { return add_region(r, p, q, 0); });
    }

    // ---------------------------------------------------------------- ECC: mul_fixed.rs (shared by all variants)

    // The window scalars (as embedded-curve scalars and as integers) of a fixed-base scalar.
    struct FixedWindows {
        std::vector<uint64_t> k; // window values in [0, 8)
    };

    void assign_fixed_constants(Region& r, size_t offset, const FB& base, size_t coords_check_toggle)
    {
        for (size_t w = 0; w < base.num_windows; ++w) {
            r.enable_selector(selector(coords_check_toggle), offset + w);
            for (size_t k = 0; k < H; ++k) {
                r.assign_fixed(fixed(F0 + k), offset + w, base.lagrange_coeffs[w][k]);
            }
            r.assign_fixed(fixed(FIXED_Z), offset + w, FF(base.z[w]));
        }
    }

    Point process_window(Region& r, size_t offset, size_t w, uint64_t k, const FB& base)
    {
        // The window point M_w(k) is precomputed in the window table (it equals [scalar] B for the scalar of
        // process_lower_bits / process_msb).
        const AffineElement mul_b = base.window_table[w][k];
        BB_ASSERT(!mul_b.x.is_zero() && !mul_b.y.is_zero());
        Cell x = r.assign_advice(advice(0), offset + w, mul_b.x);
        Cell y = r.assign_advice(advice(1), offset + w, mul_b.y);
        r.assign_advice(advice(5), offset + w, base.u[w][k]);
        return Point{ x, y };
    }

    // assign_region_inner: returns (acc after all incomplete additions, last window point)
    std::pair<Point, Point> mul_fixed_inner(
        Region& r, size_t offset, const std::vector<uint64_t>& windows, const FB& base, size_t coords_check_toggle)
    {
        BB_ASSERT_EQ(windows.size(), base.num_windows);
        assign_fixed_constants(r, offset, base, coords_check_toggle);
        Point acc = process_window(r, offset, 0, windows[0], base);
        for (size_t w = 1; w + 1 < base.num_windows; ++w) {
            Point mul_b = process_window(r, offset, w, windows[w], base);
            acc = add_incomplete_region(r, mul_b, acc, offset + w);
        }
        Point mul_b = process_window(r, offset, base.num_windows - 1, windows.back(), base);
        return { acc, mul_b };
    }

    // decompose_running_sum.rs: copy_decompose with WINDOW_NUM_BITS = 3 on the window column a4
    std::vector<Cell> running_sum_decompose(
        Region& r, size_t offset, const Cell& alpha, bool strict, size_t word_num_bits, size_t num_windows)
    {
        Cell z0 = r.copy_advice(alpha, advice(4), offset);
        for (size_t idx = 0; idx < num_windows; ++idx) {
            r.enable_selector(selector(Q_MUL_FIXED_RUNNING_SUM), offset + idx);
        }
        const auto words = decompose_word(z0.value, word_num_bits, FIXED_BASE_WINDOW_SIZE);
        BB_ASSERT_EQ(words.size(), num_windows);
        const FF inv8 = FF(8).invert();
        std::vector<Cell> zs{ z0 };
        Cell z = z0;
        for (size_t i = 0; i < num_windows; ++i) {
            z = r.assign_advice(advice(4), offset + i + 1, (z.value - FF(words[i])) * inv8);
            zs.push_back(z);
        }
        if (strict) {
            b_.constrain_constant(zs.back(), FF(0));
        }
        return zs;
    }

    static std::vector<uint64_t> windows_from_running_sum(const std::vector<Cell>& zs)
    {
        std::vector<uint64_t> out;
        for (size_t i = 0; i + 1 < zs.size(); ++i) {
            out.push_back(uint256_t(zs[i].value - zs[i + 1].value * FF(8)).data[0]);
        }
        return out;
    }

    // mul_fixed/full_width.rs
    Point mul_fixed_full(const Scalar& scalar, const FB& base)
    {
        auto [acc, mul_b] = b_.assign_region("Full-width fixed-base mul (incomplete addition)", [&](Region& r) {
            for (size_t idx = 0; idx < NUM_WINDOWS; ++idx) {
                r.enable_selector(selector(Q_MUL_FIXED_FULL), idx);
            }
            const auto windows = decompose_word(scalar, SCALAR_NUM_BITS, FIXED_BASE_WINDOW_SIZE);
            for (size_t idx = 0; idx < windows.size(); ++idx) {
                r.assign_advice(advice(4), idx, FF(windows[idx]));
            }
            return mul_fixed_inner(r, 0, windows, base, Q_MUL_FIXED_FULL);
        });
        return b_.assign_region("Full-width fixed-base mul (last window, complete addition)",
                                [&](Region& r) { return add_region(r, mul_b, acc, 0); });
    }

    // mul_fixed/short.rs (magnitude in [0, 2^64), sign in {-1, 1})
    Point mul_fixed_short(const Cell& magnitude, const Cell& sign, const FB& base)
    {
        std::vector<Cell> running_sum;
        auto [acc, mul_b] = b_.assign_region("Short fixed-base mul (incomplete addition)", [&](Region& r) {
            running_sum = running_sum_decompose(r, 0, magnitude, true, L_SCALAR_SHORT, NUM_WINDOWS_SHORT);
            return mul_fixed_inner(r, 0, windows_from_running_sum(running_sum), base, Q_MUL_FIXED_RUNNING_SUM);
        });
        return b_.assign_region("Short fixed-base mul (most significant word)", [&](Region& r) {
            Point magnitude_mul = add_region(r, mul_b, acc, 0);
            const size_t offset = 1;
            r.copy_advice(sign, advice(4), offset);
            r.copy_advice(running_sum[21], advice(5), offset);
            const FF y_val = (sign.value == -FF(1)) ? -magnitude_mul.y.value : magnitude_mul.y.value;
            r.enable_selector(selector(Q_MUL_FIXED_SHORT), offset);
            Cell y_var = r.assign_advice(advice(1), offset, y_val);
            return Point{ magnitude_mul.x, y_var };
        });
    }

    // mul_fixed/base_field_elem.rs
    Point mul_fixed_base_field(const Cell& scalar, const FB& base)
    {
        std::vector<Cell> running_sum;
        auto [acc, mul_b] = b_.assign_region("Base-field elem fixed-base mul (incomplete addition)", [&](Region& r) {
            running_sum = running_sum_decompose(r, 0, scalar, true, BASE_NUM_BITS, NUM_WINDOWS);
            return mul_fixed_inner(r, 0, windows_from_running_sum(running_sum), base, Q_MUL_FIXED_RUNNING_SUM);
        });
        Point result = b_.assign_region("Base-field elem fixed-base mul (complete addition)",
                                        [&](Region& r) { return add_region(r, mul_b, acc, 0); });

        const Cell& alpha = running_sum[0];
        const Cell z_43_alpha = running_sum[43];
        const Cell z_44_alpha = running_sum[44];
        const Cell z_84_alpha = running_sum[84];
        const FF alpha_0 = alpha.value - z_84_alpha.value * two_pow(252);
        const FF alpha_0_prime = alpha_0 + two_pow(130) - FF(t_p());
        const auto zs = witness_check(alpha_0_prime, 13, false);
        const Cell alpha_0_prime_cell = zs[0];
        const Cell z_13_alpha_0_prime = zs[13];

        b_.assign_region("Canonicity checks", [&](Region& r) {
            r.enable_selector(selector(Q_MUL_FIXED_BASE_FIELD), 1);
            r.copy_advice(alpha, advice(6), 0);
            r.copy_advice(z_84_alpha, advice(8), 0);
            r.copy_advice(alpha_0_prime_cell, advice(6), 1);
            r.assign_advice(advice(7), 1, bitrange_subset(alpha.value, 252, 254));
            r.assign_advice(advice(8), 1, bitrange_subset(alpha.value, 254, 255));
            r.copy_advice(z_13_alpha_0_prime, advice(6), 2);
            r.copy_advice(z_44_alpha, advice(7), 2);
            r.copy_advice(z_43_alpha, advice(8), 2);
            return 0;
        });
        return result;
    }

    // ---------------------------------------------------------------- ECC: variable-base mul (mul.rs)

    struct IncompleteConfig {
        size_t z, x_a, x_p, y_p, lambda_1, lambda_2;
        size_t q_mul_1, q_mul_2, q_mul_3;
    };
    static IncompleteConfig hi_config() { return { 9, 3, 0, 1, 4, 5, Q_MUL_HI_1, Q_MUL_HI_2, Q_MUL_HI_3 }; }
    static IncompleteConfig lo_config() { return { 6, 7, 0, 1, 8, 2, Q_MUL_LO_1, Q_MUL_LO_2, Q_MUL_LO_3 }; }

    struct IncompleteOutput {
        Cell x_a;
        Cell y_a;
        std::vector<Cell> zs;
    };

    // mul/incomplete.rs double_and_add
    IncompleteOutput double_and_add(Region& r,
                                    const IncompleteConfig& cfg,
                                    size_t offset,
                                    const Point& base,
                                    const std::vector<bool>& bits,
                                    const Cell& acc_x,
                                    const Cell& acc_y,
                                    const Cell& acc_z)
    {
        const size_t num_bits = bits.size();
        const FF x_p = base.x.value;
        const FF y_p_base = base.y.value;
        r.enable_selector(selector(cfg.q_mul_1), offset);
        for (size_t idx = 0; idx + 1 < num_bits; ++idx) {
            r.enable_selector(selector(cfg.q_mul_2), offset + 1 + idx);
        }
        r.enable_selector(selector(cfg.q_mul_3), offset + 1 + num_bits - 1);

        Cell z = r.copy_advice(acc_z, advice(cfg.z), offset);
        Cell x_a = r.copy_advice(acc_x, advice(cfg.x_a), offset + 1);
        Cell y_a_cell = r.copy_advice(acc_y, advice(cfg.lambda_1), offset);
        FF y_a = y_a_cell.value;
        offset += 1;

        std::vector<Cell> zs;
        for (size_t row = 0; row < num_bits; ++row) {
            const bool k = bits[row];
            z = r.assign_advice(advice(cfg.z), row + offset, FF(2) * z.value + FF(k ? 1 : 0));
            zs.push_back(z);
            if (anchored_base_ && row == 0) {
                r.copy_advice(base.x, advice(cfg.x_p), row + offset);
                r.copy_advice(base.y, advice(cfg.y_p), row + offset);
            } else {
                r.assign_advice(advice(cfg.x_p), row + offset, x_p);
                r.assign_advice(advice(cfg.y_p), row + offset, y_p_base);
            }
            const FF y_p = k ? y_p_base : -y_p_base;
            const FF lambda1 = (y_a - y_p) * (x_a.value - x_p).invert();
            r.assign_advice(advice(cfg.lambda_1), row + offset, lambda1);
            const FF x_r = lambda1.sqr() - x_a.value - x_p;
            const FF lambda2 = y_a * FF(2) * (x_a.value - x_r).invert() - lambda1;
            r.assign_advice(advice(cfg.lambda_2), row + offset, lambda2);
            const FF x_a_new = lambda2.sqr() - x_a.value - x_r;
            y_a = lambda2 * (x_a.value - x_a_new) - y_a;
            x_a = r.assign_advice(advice(cfg.x_a), row + offset + 1, x_a_new);
        }
        Cell y_a_final = r.assign_advice(advice(cfg.lambda_1), offset + num_bits, y_a);
        return { x_a, y_a_final, zs };
    }

    // mul.rs decompose_for_scalar_mul: bits of k = alpha + t_q, most significant first
    static std::vector<bool> decompose_for_scalar_mul(const FF& alpha)
    {
        const uint256_t k = uint256_t(alpha) + t_q();
        std::vector<bool> bits(SCALAR_NUM_BITS);
        for (size_t i = 0; i < SCALAR_NUM_BITS; ++i) {
            bits[SCALAR_NUM_BITS - 1 - i] = k.get_bit(i);
        }
        return bits;
    }

    // mul.rs Config::assign: [alpha] base for alpha a circuit-field element
    Point mul_var(const Cell& alpha, const Point& base)
    {
        std::vector<Cell> zs_all;
        Point result = b_.assign_region("variable-base scalar mul", [&](Region& r) {
            const auto bits = decompose_for_scalar_mul(alpha.value);
            std::vector<bool> bits_hi(bits.begin(), bits.begin() + INCOMPLETE_HI_LEN);
            std::vector<bool> bits_lo(bits.begin() + INCOMPLETE_HI_LEN, bits.begin() + INCOMPLETE_LEN);
            std::vector<bool> bits_complete(bits.begin() + INCOMPLETE_LEN, bits.begin() + INCOMPLETE_LEN + 3);
            const bool lsb = bits[SCALAR_NUM_BITS - 1];

            size_t offset = 0;
            Point acc = add_region(r, base, base, offset);
            offset += 1;
            Cell z_init = r.assign_advice_from_constant(advice(hi_config().z), offset, FF(0));
            auto hi = double_and_add(r, hi_config(), offset, base, bits_hi, acc.x, acc.y, z_init);
            auto lo = double_and_add(r, lo_config(), offset, base, bits_lo, hi.x_a, hi.y_a, hi.zs.back());
            offset += INCOMPLETE_LO_LEN + 2;

            // complete.rs assign_region (z_complete = a9, add config)
            Point acc_c{ lo.x_a, lo.y_a };
            std::vector<Cell> zs_complete;
            {
                for (size_t row = 0; row < NUM_COMPLETE_BITS; ++row) {
                    r.enable_selector(selector(Q_MUL_DECOMPOSE_VAR), 2 * row + offset + 1);
                }
                Cell z = r.copy_advice(lo.zs.back(), advice(9), offset);
                for (size_t iter = 0; iter < NUM_COMPLETE_BITS; ++iter) {
                    const size_t row = 2 * iter;
                    const bool k = bits_complete[iter];
                    z = r.assign_advice(advice(9), row + offset + 2, z.value * FF(2) + FF(k ? 1 : 0));
                    zs_complete.push_back(z);
                    Cell base_y = r.copy_advice(base.y, advice(9), row + offset + 1);
                    Cell y_p = r.assign_advice(advice(1), row + offset, k ? base_y.value : -base_y.value);
                    Point U{ base.x, y_p };
                    Point tmp_acc = add_region(r, U, acc_c, row + offset);
                    acc_c = add_region(r, acc_c, tmp_acc, row + offset + 1);
                }
            }
            offset += NUM_COMPLETE_BITS * 2;

            // process_lsb
            r.enable_selector(selector(Q_MUL_LSB), offset);
            const Cell& z_1 = zs_complete.back();
            Cell z_0 = r.assign_advice(advice(9), offset + 1, z_1.value * FF(2) + FF(lsb ? 1 : 0));
            r.copy_advice(base.x, advice(0), offset + 1);
            r.copy_advice(base.y, advice(1), offset + 1);
            Cell x = r.assign_advice(advice(0), offset, lsb ? FF(0) : base.x.value);
            Cell y = r.assign_advice(advice(1), offset, lsb ? FF(0) : -base.y.value);
            Point res = add_region(r, Point{ x, y }, acc_c, offset);

            // zs = [z_init, hi..., lo..., complete..., z_0] reversed
            zs_all.push_back(z_init);
            zs_all.insert(zs_all.end(), hi.zs.begin(), hi.zs.end());
            zs_all.insert(zs_all.end(), lo.zs.begin(), lo.zs.end());
            zs_all.insert(zs_all.end(), zs_complete.begin(), zs_complete.end());
            zs_all.push_back(z_0);
            BB_ASSERT_EQ(zs_all.size(), SCALAR_NUM_BITS + 1);
            std::reverse(zs_all.begin(), zs_all.end());
            return res;
        });
        overflow_check(alpha, zs_all);
        return result;
    }

    // mul/overflow.rs (advices a6, a7, a8)
    void overflow_check(const Cell& alpha, const std::vector<Cell>& zs)
    {
        const Cell& k_254 = zs[254];
        const FF s_val = alpha.value + k_254.value * two_pow(130);
        Cell s = b_.assign_region("s = alpha + k_254 * 2^130",
                                  [&](Region& r) { return r.assign_advice(advice(6), 0, s_val); });
        const auto s_zs = copy_check(s, 13, false);
        const Cell s_minus_lo_130 = s_zs.back();
        b_.assign_region("overflow check", [&](Region& r) {
            r.enable_selector(selector(Q_MUL_OVERFLOW), 1);
            r.copy_advice(zs[0], advice(6), 0);
            r.copy_advice(zs[130], advice(6), 1);
            r.assign_advice(advice(6), 2, inv0(zs[130].value));
            r.copy_advice(k_254, advice(7), 0);
            r.copy_advice(alpha, advice(7), 1);
            r.copy_advice(s_minus_lo_130, advice(7), 2);
            r.copy_advice(s, advice(8), 1);
            return 0;
        });
    }

    // ---------------------------------------------------------------- Poseidon Pow5 chip (state a6..a8, sbox a5)

    using Poseidon = PoseidonP128Pow5T3<FF>;
    static constexpr std::array<size_t, 3> POSEIDON_STATE = { 6, 7, 8 };
    static constexpr size_t POSEIDON_PARTIAL_SBOX = 5;
    static constexpr std::array<size_t, 3> POSEIDON_RC_A = { F2, F3, F4 };
    static constexpr std::array<size_t, 3> POSEIDON_RC_B = { F5, F6, F7 };

    // Poseidon::Hash<ConstantLength<2>>::hash(nk, rho)
    Cell poseidon_hash(const Cell& a, const Cell& b)
    {
        const auto& p = Poseidon::params();
        std::array<Cell, 3> state = b_.assign_region("initial state for domain ConstantLength<2>", [&](Region& r) {
            return std::array<Cell, 3>{
                r.assign_advice_from_constant(advice(POSEIDON_STATE[0]), 0, FF(0)),
                r.assign_advice_from_constant(advice(POSEIDON_STATE[1]), 0, FF(0)),
                r.assign_advice_from_constant(advice(POSEIDON_STATE[2]), 0, Poseidon::constant_length_capacity(2)),
            };
        });
        // add_input
        state = b_.assign_region("add input for domain ConstantLength<2>", [&](Region& r) {
            r.enable_selector(selector(Q_POSEIDON_PAD_AND_ADD), 1);
            std::array<Cell, 3> init{};
            for (size_t i = 0; i < 3; ++i) {
                init[i] = r.copy_advice(state[i], advice(POSEIDON_STATE[i]), 0);
            }
            std::array<Cell, 2> input{ r.copy_advice(a, advice(POSEIDON_STATE[0]), 1),
                                       r.copy_advice(b, advice(POSEIDON_STATE[1]), 1) };
            return std::array<Cell, 3>{
                r.assign_advice(advice(POSEIDON_STATE[0]), 2, init[0].value + input[0].value),
                r.assign_advice(advice(POSEIDON_STATE[1]), 2, init[1].value + input[1].value),
                r.assign_advice(advice(POSEIDON_STATE[2]), 2, init[2].value),
            };
        });
        // permute
        state = b_.assign_region("permute state", [&](Region& r) {
            std::array<Cell, 3> s{};
            for (size_t i = 0; i < 3; ++i) {
                s[i] = r.copy_advice(state[i], advice(POSEIDON_STATE[i]), 0);
            }
            auto mds_mul = [&](const std::array<FF, 3>& v) {
                std::array<FF, 3> out{};
                for (size_t i = 0; i < 3; ++i) {
                    out[i] = p.mds[i][0] * v[0] + p.mds[i][1] * v[1] + p.mds[i][2] * v[2];
                }
                return out;
            };
            auto load_rc = [&](const std::array<size_t, 3>& cols, size_t round, size_t offset) {
                for (size_t i = 0; i < 3; ++i) {
                    r.assign_fixed(fixed(cols[i]), offset, p.round_constants[round][i]);
                }
            };
            auto store_state = [&](const std::array<FF, 3>& v, size_t offset) {
                for (size_t i = 0; i < 3; ++i) {
                    s[i] = r.assign_advice(advice(POSEIDON_STATE[i]), offset, v[i]);
                }
            };
            auto full_round = [&](size_t round, size_t offset) {
                r.enable_selector(selector(Q_POSEIDON_FULL), offset);
                load_rc(POSEIDON_RC_A, round, offset);
                std::array<FF, 3> v{};
                for (size_t i = 0; i < 3; ++i) {
                    v[i] = Poseidon::sbox(s[i].value + p.round_constants[round][i]);
                }
                store_state(mds_mul(v), offset + 1);
            };
            auto partial_round = [&](size_t round, size_t offset) {
                r.enable_selector(selector(Q_POSEIDON_PARTIAL), offset);
                load_rc(POSEIDON_RC_A, round, offset);
                std::array<FF, 3> v{ Poseidon::sbox(s[0].value + p.round_constants[round][0]),
                                     s[1].value + p.round_constants[round][1],
                                     s[2].value + p.round_constants[round][2] };
                r.assign_advice(advice(POSEIDON_PARTIAL_SBOX), offset, v[0]);
                auto mid = mds_mul(v);
                load_rc(POSEIDON_RC_B, round + 1, offset);
                std::array<FF, 3> w{ Poseidon::sbox(mid[0] + p.round_constants[round + 1][0]),
                                     mid[1] + p.round_constants[round + 1][1],
                                     mid[2] + p.round_constants[round + 1][2] };
                store_state(mds_mul(w), offset + 1);
            };
            const size_t half_full = Poseidon::FULL_ROUNDS / 2;
            const size_t half_partial = Poseidon::PARTIAL_ROUNDS / 2;
            for (size_t rnd = 0; rnd < half_full; ++rnd) {
                full_round(rnd, rnd);
            }
            for (size_t rnd = 0; rnd < half_partial; ++rnd) {
                partial_round(half_full + 2 * rnd, half_full + rnd);
            }
            for (size_t rnd = 0; rnd < half_full; ++rnd) {
                full_round(half_full + 2 * half_partial + rnd, half_full + half_partial + rnd);
            }
            return s;
        });
        return state[0];
    }

    // ---------------------------------------------------------------- AddChip (a = a7, b = a8, c = a6)

    Cell add_field(const Cell& a, const Cell& b)
    {
        return b_.assign_region("c = a + b", [&](Region& r) {
            r.enable_selector(selector(Q_ADD_FIELD), 0);
            r.copy_advice(a, advice(7), 0);
            r.copy_advice(b, advice(8), 0);
            return r.assign_advice(advice(6), 0, a.value + b.value);
        });
    }

    // ---------------------------------------------------------------- Merkle (sinsemilla/merkle/chip.rs, cond_swap.rs)

    std::pair<Cell, Cell> cond_swap(const SinsemillaConfig& cfg, const Cell& a, const FF& b, bool swap)
    {
        const auto adv = cfg.advices();
        return b_.assign_region("swap", [&](Region& r) {
            r.enable_selector(selector(cfg.q_swap), 0);
            Cell a_cell = r.copy_advice(a, advice(adv[0]), 0);
            Cell b_cell = r.assign_advice(advice(adv[1]), 0, b);
            r.assign_advice(advice(adv[4]), 0, FF(swap ? 1 : 0));
            Cell a_swapped = r.assign_advice(advice(adv[2]), 0, swap ? b_cell.value : a_cell.value);
            Cell b_swapped = r.assign_advice(advice(adv[3]), 0, swap ? a_cell.value : b_cell.value);
            return std::make_pair(a_swapped, b_swapped);
        });
    }

    Cell merkle_hash_layer(
        const SinsemillaConfig& cfg, const AffineElement& Q, size_t l, const Cell& left, const Cell& right)
    {
        const auto adv = cfg.advices();
        MessagePiece a = from_subpieces(cfg, { { FF(l), 10 }, bitrange_of(left.value, 0, 240) });
        const auto b_0 = bitrange_of(left.value, 240, 250);
        const auto b_1 = witness_short(left.value, 250, BASE_NUM_BITS);
        const auto b_2 = witness_short(right.value, 0, 5);
        MessagePiece b =
            from_subpieces(cfg, { b_0, { b_1.cell.value, b_1.num_bits }, { b_2.cell.value, b_2.num_bits } });
        MessagePiece c = from_subpieces(cfg, { bitrange_of(right.value, 5, BASE_NUM_BITS) });
        auto h = hash_to_point(cfg, Q, { a, b, c });
        const Cell z1_a = h.zs[0][1];
        const Cell z1_b = h.zs[1][1];
        b_.assign_region("Check piece decomposition", [&](Region& r) {
            r.enable_selector(selector(cfg.q_decompose), 0);
            r.assign_advice_from_constant(advice(adv[4]), 1, FF(l));
            r.copy_advice(a.cell, advice(adv[0]), 0);
            r.copy_advice(b.cell, advice(adv[1]), 0);
            r.copy_advice(c.cell, advice(adv[2]), 0);
            r.copy_advice(left, advice(adv[3]), 0);
            r.copy_advice(right, advice(adv[4]), 0);
            r.copy_advice(z1_a, advice(adv[0]), 1);
            r.copy_advice(z1_b, advice(adv[1]), 1);
            r.copy_advice(b_1.cell, advice(adv[2]), 1);
            r.copy_advice(b_2.cell, advice(adv[3]), 1);
            return 0;
        });
        return h.point.x;
    }

    // MerklePath::calculate_root with two Merkle chips (PAR = 2): layers 0..15 on chip 1, 16..31 on chip 2.
    Cell merkle_root(const Cell& leaf, const std::array<FF, 32>& path, uint32_t pos)
    {
        const auto Q = O::constants().q_merkle_crh;
        Cell node = leaf;
        for (size_t l = 0; l < 32; ++l) {
            const auto cfg = (l < 16) ? sinsemilla_config_1() : sinsemilla_config_2();
            auto [left, right] = cond_swap(cfg, node, path[l], ((pos >> l) & 1) != 0);
            node = merkle_hash_layer(cfg, Q, l, left, right);
        }
        return node;
    }

    // ---------------------------------------------------------------- Sinsemilla commitments

    // CommitDomain::commit: [r] R (full-width fixed-base mul), then HashToPoint, then complete addition
    std::pair<Point, std::vector<std::vector<Cell>>> commit(const SinsemillaConfig& cfg,
                                                            const AffineElement& Q,
                                                            const FB& r_base,
                                                            const std::vector<MessagePiece>& pieces,
                                                            const Scalar& r)
    {
        Point blind = mul_fixed_full(r, r_base);
        auto h = hash_to_point(cfg, Q, pieces);
        Point commitment = add(h.point, blind);
        return { commitment, std::move(h.zs) };
    }

  private:
    Builder& b_;
    bool anchored_base_;
};

} // namespace bb::zcash::halo2
