#include "barretenberg/zcash/ultra_pasta/action_circuit_ultra_pasta.hpp"

namespace bb::zcash::ultra_pasta {

namespace {
using FF = ActionCircuitUltraPasta::FF;
using Scalar = ActionCircuitUltraPasta::Scalar;
using AffineElement = ActionCircuitUltraPasta::AffineElement;
using Element = ActionCircuitUltraPasta::Element;

constexpr size_t K = 10;
constexpr size_t LO_BITS = 130;  // bits 0..129 of a canonical field element
constexpr size_t TOP_BIT = 254;  // bit 254
constexpr size_t MAX_SLICE = 10; // range constraints of up to 10 bits are direct, longer ones use 10-bit limbs
constexpr size_t WINDOW_BITS = 3;
constexpr size_t WINDOW_SIZE = 8;
constexpr size_t NUM_FULL_WINDOWS = 85;  // windows of a 255-bit scalar
constexpr size_t LAST_WINDOW = 84;       // the top window of a 255-bit scalar
constexpr size_t NUM_SHORT_WINDOWS = 22; // windows of a 64-bit scalar (the last one holds bit 63 only)

FF pow2(size_t k)
{
    return FF(uint256_t(1) << k);
}

// t_p = p - 2^254
const FF& t_p()
{
    static const FF t = FF(uint256_t(FF::modulus) - (uint256_t(1) << TOP_BIT));
    return t;
}

AffineElement add(const AffineElement& p, const AffineElement& q)
{
    return AffineElement(Element(p) + Element(q));
}

} // namespace

ActionCircuitUltraPasta::ActionCircuitUltraPasta(Builder& builder)
    : builder_(builder)
{
    const auto& table = O::Sins::S_table();
    std::vector<Builder::TableRow> rows;
    rows.reserve(table.size());
    for (size_t i = 0; i < table.size(); ++i) {
        rows.push_back({ FF(i), table[i].x, table[i].y });
    }
    s_table_ = &builder_.create_pasta_table(std::move(rows));
}

uint32_t ActionCircuitUltraPasta::linear_combination(const std::vector<Term>& terms, const FF& constant_term)
{
    const uint32_t zero = builder_.zero_idx();
    FF total = constant_term;
    for (const auto& t : terms) {
        total += t.coeff * value(t.var);
    }
    const uint32_t out = witness(total);
    auto term = [&](size_t i) { return i < terms.size() ? terms[i] : Term{ zero, FF(0) }; };
    if (terms.size() <= 3) {
        const auto a = term(0);
        const auto b = term(1);
        const auto c = term(2);
        builder_.create_big_add_gate({ a.var, b.var, c.var, out, a.coeff, b.coeff, c.coeff, FF(-1), constant_term });
        return out;
    }
    // Row layout: the first row holds four terms; every later row holds A_j = -(partial sum so far) in w_4 (written
    // by the previous row through w_4_shift); middle rows add three terms, the last row two terms and the output.
    FF partial = constant_term;
    for (size_t i = 0; i < 4; ++i) {
        partial += terms[i].coeff * value(terms[i].var);
    }
    uint32_t acc = witness(-partial);
    builder_.create_big_add_gate({ terms[0].var,
                                   terms[1].var,
                                   terms[2].var,
                                   terms[3].var,
                                   terms[0].coeff,
                                   terms[1].coeff,
                                   terms[2].coeff,
                                   terms[3].coeff,
                                   constant_term },
                                 /*include_next_gate_w_4=*/true);
    size_t i = 4;
    while (terms.size() - i > 2) {
        const auto a = term(i);
        const auto b = term(i + 1);
        const auto c = term(i + 2);
        partial += a.coeff * value(a.var) + b.coeff * value(b.var) + c.coeff * value(c.var);
        const uint32_t next = witness(-partial);
        builder_.create_big_add_gate({ a.var, b.var, c.var, acc, a.coeff, b.coeff, c.coeff, FF(-1), FF(0) },
                                     /*include_next_gate_w_4=*/true);
        acc = next;
        i += 3;
    }
    const auto a = term(i);
    const auto b = term(i + 1);
    builder_.create_big_add_gate({ a.var, b.var, out, acc, a.coeff, b.coeff, FF(-1), FF(-1), FF(0) });
    return out;
}

uint32_t ActionCircuitUltraPasta::mul(uint32_t a, uint32_t b)
{
    const uint32_t out = witness(value(a) * value(b));
    builder_.create_arithmetic_gate({ a, b, out, FF(1), FF(0), FF(0), FF(-1), FF(0) });
    return out;
}

void ActionCircuitUltraPasta::range_constrain(uint32_t var, size_t num_bits)
{
    if (num_bits <= MAX_SLICE) {
        builder_.create_small_range_constraint(var, (uint64_t{ 1 } << num_bits) - 1, "range");
    } else {
        builder_.create_limbed_range_constraint(var, num_bits, MAX_SLICE, "range");
    }
}

uint32_t ActionCircuitUltraPasta::boolean_witness(bool bit)
{
    const uint32_t var = witness(FF(bit ? 1 : 0));
    range_constrain(var, 1);
    return var;
}

ActionCircuitUltraPasta::Point ActionCircuitUltraPasta::witness_point(const AffineElement& p)
{
    BB_ASSERT(p.on_curve() && !p.is_point_at_infinity());
    const Point out{ witness(p.x), witness(p.y) };
    // y^2 = x^3 + 5
    const uint32_t x2 = mul(out.x, out.x);
    const uint32_t x3 = mul(x2, out.x);
    builder_.create_arithmetic_gate({ out.y, out.y, x3, FF(1), FF(0), FF(0), FF(-1), -pallas::g1::curve_b });
    return out;
}

ActionCircuitUltraPasta::Point ActionCircuitUltraPasta::incomplete_add(const Point& p, const Point& q)
{
    const AffineElement pv = point_value(p);
    const AffineElement qv = point_value(q);
    BB_ASSERT(pv.x != qv.x, "incomplete addition exceptional case");
    const AffineElement r = add(pv, qv);
    const Point out{ witness(r.x), witness(r.y) };
    builder_.create_ecc_add_gate(
        { .x1 = p.x, .y1 = p.y, .x2 = q.x, .y2 = q.y, .x3 = out.x, .y3 = out.y, .is_addition = true });
    return out;
}

ActionCircuitUltraPasta::Point ActionCircuitUltraPasta::dbl(const Point& p)
{
    const AffineElement r(Element(point_value(p)).dbl());
    const Point out{ witness(r.x), witness(r.y) };
    builder_.create_ecc_dbl_gate({ .x1 = p.x, .y1 = p.y, .x3 = out.x, .y3 = out.y });
    return out;
}

ActionCircuitUltraPasta::Point ActionCircuitUltraPasta::complete_add(const Point& p, const Point& q)
{
    const AffineElement pv = point_value(p);
    const AffineElement qv = point_value(q);
    const bool same_x = pv.x == qv.x;
    BB_ASSERT(!same_x || pv.y == qv.y, "complete addition of P and -P");
    const FF lambda = same_x ? (pv.x.sqr() * FF(3)) / (pv.y + pv.y) : (qv.y - pv.y) / (qv.x - pv.x);
    const FF x3 = lambda.sqr() - pv.x - qv.x;
    const FF y3 = lambda * (pv.x - x3) - pv.y;

    // d = x_q - x_p, e = 1 - d * inv(d) (e = 1 iff d = 0), d * e = 0
    const uint32_t d = linear_combination({ { q.x, FF(1) }, { p.x, FF(-1) } });
    const uint32_t inv = witness(same_x ? FF(0) : value(d).invert());
    const uint32_t e = witness(FF(same_x ? 1 : 0));
    builder_.create_arithmetic_gate({ d, inv, e, FF(1), FF(0), FF(0), FF(1), FF(-1) });
    builder_.create_arithmetic_gate({ d, e, builder_.zero_idx(), FF(1), FF(0), FF(0), FF(0), FF(0) });
    // s = y_q - y_p; e * s = 0 excludes P = -Q; lambda * d = s
    const uint32_t s = linear_combination({ { q.y, FF(1) }, { p.y, FF(-1) } });
    builder_.create_arithmetic_gate({ e, s, builder_.zero_idx(), FF(1), FF(0), FF(0), FF(0), FF(0) });
    const uint32_t lam = witness(lambda);
    builder_.create_arithmetic_gate({ lam, d, s, FF(1), FF(0), FF(0), FF(-1), FF(0) });
    // e * (2 y_p lambda - 3 x_p^2) = 0
    const uint32_t xp2 = mul(p.x, p.x);
    const uint32_t w = witness(FF(2) * pv.y * lambda - FF(3) * value(xp2));
    builder_.create_big_mul_add_gate({ p.y, lam, xp2, w, FF(2), FF(0), FF(0), FF(-3), FF(-1), FF(0) });
    builder_.create_arithmetic_gate({ e, w, builder_.zero_idx(), FF(1), FF(0), FF(0), FF(0), FF(0) });
    // x3 = lambda^2 - x_p - x_q, y3 = lambda (x_p - x3) - y_p
    const uint32_t sum_x = linear_combination({ { p.x, FF(1) }, { q.x, FF(1) } });
    const Point out{ witness(x3), witness(y3) };
    builder_.create_big_mul_add_gate({ lam, lam, sum_x, out.x, FF(1), FF(0), FF(0), FF(-1), FF(-1), FF(0) });
    const uint32_t dx = linear_combination({ { p.x, FF(1) }, { out.x, FF(-1) } });
    builder_.create_big_mul_add_gate({ lam, dx, p.y, out.y, FF(1), FF(0), FF(0), FF(-1), FF(-1), FF(0) });
    return out;
}

ActionCircuitUltraPasta::Point ActionCircuitUltraPasta::conditional_negate(const Point& p, uint32_t sign)
{
    return { p.x, mul(sign, p.y) };
}

void ActionCircuitUltraPasta::assert_equal(const Point& p, const Point& q)
{
    builder_.assert_equal(p.x, q.x, "point x");
    builder_.assert_equal(p.y, q.y, "point y");
}

std::vector<uint32_t> ActionCircuitUltraPasta::message_words(const std::vector<Segment>& segments)
{
    size_t total_bits = 0;
    for (const auto& s : segments) {
        total_bits += s.num_bits;
    }
    const size_t num_words = (total_bits + K - 1) / K;
    std::vector<std::vector<Term>> word_terms(num_words);

    struct Slice {
        uint32_t var;
        size_t start; // within the segment
    };
    std::vector<std::vector<Slice>> segment_slices(segments.size());
    size_t offset = 0; // global bit offset of the segment
    for (size_t si = 0; si < segments.size(); ++si) {
        const auto& seg = segments[si];
        const uint256_t v(value(seg.value));
        std::vector<size_t> cuts{ 0 };
        for (size_t b = 1; b < seg.num_bits; ++b) {
            if ((offset + b) % K == 0 || (seg.canonical && (b == LO_BITS || b == TOP_BIT))) {
                cuts.push_back(b);
            }
        }
        cuts.push_back(seg.num_bits);
        for (size_t i = 0; i + 1 < cuts.size(); ++i) {
            const size_t start = cuts[i];
            const size_t width = cuts[i + 1] - start;
            const FF slice_value(v.slice(start, start + width));
            const size_t global = offset + start;
            const bool whole_word = width == K && global % K == 0;
            uint32_t var = 0;
            if (seg.is_constant) {
                var = constant(slice_value);
            } else {
                var = witness(slice_value);
                // A whole word is range-checked by its S-table read.
                if (!whole_word) {
                    range_constrain(var, width);
                }
            }
            segment_slices[si].push_back({ var, start });
            word_terms[global / K].push_back({ var, pow2(global % K) });
        }
        offset += seg.num_bits;
    }

    std::vector<uint32_t> words;
    words.reserve(num_words);
    for (const auto& terms : word_terms) {
        words.push_back(terms.size() == 1 && terms[0].coeff == FF(1) ? terms[0].var : linear_combination(terms));
    }

    for (size_t si = 0; si < segments.size(); ++si) {
        const auto& seg = segments[si];
        if (seg.is_constant) {
            continue;
        }
        const auto& slices = segment_slices[si];
        if (seg.canonical) {
            BB_ASSERT_EQ(seg.num_bits, O::L_BASE);
            std::vector<Term> lo;
            std::vector<Term> mid;
            uint32_t top = builder_.zero_idx();
            for (const auto& s : slices) {
                if (s.start < LO_BITS) {
                    lo.push_back({ s.var, pow2(s.start) });
                } else if (s.start < TOP_BIT) {
                    mid.push_back({ s.var, pow2(s.start - LO_BITS) });
                } else {
                    top = s.var;
                }
            }
            const uint32_t lo_var = linear_combination(lo);
            const uint32_t mid_var = linear_combination(mid);
            builder_.create_big_add_gate(
                { lo_var, mid_var, top, seg.value, FF(1), pow2(LO_BITS), pow2(TOP_BIT), FF(-1), FF(0) });
            assert_canonical(lo_var, mid_var, top, LO_BITS);
        } else {
            std::vector<Term> terms;
            for (const auto& s : slices) {
                terms.push_back({ s.var, pow2(s.start) });
            }
            terms.push_back({ seg.value, FF(-1) });
            builder_.assert_equal(linear_combination(terms), builder_.zero_idx(), "message recomposition");
        }
    }
    return words;
}

ActionCircuitUltraPasta::Point ActionCircuitUltraPasta::hash_to_point(const AffineElement& q,
                                                                      const std::vector<uint32_t>& words)
{
    const auto& table = O::Sins::S_table();
    Point acc = constant_point(q);
    for (const uint32_t word : words) {
        const size_t m = static_cast<size_t>(uint256_t(value(word)).data[0]);
        BB_ASSERT_LT(m, table.size());
        const Point s{ witness(table[m].x), witness(table[m].y) };
        builder_.create_pasta_lookup(*s_table_, m, word, s.x, s.y);
        const Point t = incomplete_add(acc, s);
        acc = incomplete_add(t, acc);
    }
    return acc;
}

void ActionCircuitUltraPasta::assert_canonical(uint32_t lo, uint32_t mid, uint32_t top, size_t lo_bits)
{
    // top = 1 => mid = 0 and lo < t_p, i.e. lo + 2^lo_bits - t_p < 2^lo_bits.
    builder_.create_arithmetic_gate({ top, mid, builder_.zero_idx(), FF(1), FF(0), FF(0), FF(0), FF(0) });
    const FF shift = pow2(lo_bits) - t_p();
    const uint32_t z = witness(value(top) * (value(lo) + shift));
    builder_.create_arithmetic_gate({ top, lo, z, FF(1), shift, FF(0), FF(-1), FF(0) });
    range_constrain(z, lo_bits);
}

uint32_t ActionCircuitUltraPasta::y_lsb(uint32_t y)
{
    const uint256_t v(value(y));
    const uint32_t lsb = boolean_witness(v.get_bit(0));
    const uint32_t rest = witness(FF(v.slice(1, LO_BITS)));
    const uint32_t mid = witness(FF(v.slice(LO_BITS, TOP_BIT)));
    const uint32_t top = boolean_witness(v.get_bit(TOP_BIT));
    range_constrain(rest, LO_BITS - 1);
    range_constrain(mid, TOP_BIT - LO_BITS);
    const uint32_t lo = linear_combination({ { lsb, FF(1) }, { rest, FF(2) } });
    builder_.create_big_add_gate({ lo, mid, top, y, FF(1), pow2(LO_BITS), pow2(TOP_BIT), FF(-1), FF(0) });
    assert_canonical(lo, mid, top, LO_BITS);
    return lsb;
}

std::vector<uint32_t> ActionCircuitUltraPasta::scalar_bits(const uint256_t& k, size_t num_bits)
{
    std::vector<uint32_t> bits;
    bits.reserve(num_bits);
    for (size_t i = 0; i < num_bits; ++i) {
        bits.push_back(boolean_witness(k.get_bit(i)));
    }
    return bits;
}

std::vector<uint32_t> ActionCircuitUltraPasta::field_bits(uint32_t x)
{
    const auto bits = scalar_bits(uint256_t(value(x)), O::L_BASE);
    std::vector<Term> lo;
    std::vector<Term> mid;
    for (size_t i = 0; i < TOP_BIT; ++i) {
        if (i < LO_BITS) {
            lo.push_back({ bits[i], pow2(i) });
        } else {
            mid.push_back({ bits[i], pow2(i - LO_BITS) });
        }
    }
    const uint32_t lo_var = linear_combination(lo);
    const uint32_t mid_var = linear_combination(mid);
    builder_.create_big_add_gate(
        { lo_var, mid_var, bits[TOP_BIT], x, FF(1), pow2(LO_BITS), pow2(TOP_BIT), FF(-1), FF(0) });
    assert_canonical(lo_var, mid_var, bits[TOP_BIT], LO_BITS);
    return bits;
}

std::vector<uint32_t> ActionCircuitUltraPasta::scalar_windows(const uint256_t& k, size_t num_windows)
{
    BB_ASSERT_LTE(k.get_msb() + 1, num_windows * WINDOW_BITS);
    std::vector<uint32_t> windows;
    windows.reserve(num_windows);
    for (size_t w = 0; w < num_windows; ++w) {
        windows.push_back(witness(FF(k.slice(w * WINDOW_BITS, (w + 1) * WINDOW_BITS))));
    }
    return windows;
}

std::vector<uint32_t> ActionCircuitUltraPasta::field_windows(uint32_t x)
{
    // 85 windows; the canonical split is at bit 129 (window 43) and the top window is (bits 252, 253) + 4 * bit 254.
    constexpr size_t SPLIT = 43;
    constexpr size_t LO = SPLIT * WINDOW_BITS;
    const uint256_t v(value(x));
    auto windows = scalar_windows(v, NUM_FULL_WINDOWS);
    const uint32_t top_lo = witness(FF(v.slice(LAST_WINDOW * WINDOW_BITS, TOP_BIT)));
    const uint32_t top_bit = boolean_witness(v.get_bit(TOP_BIT));
    range_constrain(top_lo, TOP_BIT - (LAST_WINDOW * WINDOW_BITS));
    builder_.create_big_add_gate(
        { top_lo, top_bit, windows[LAST_WINDOW], builder_.zero_idx(), FF(1), FF(4), FF(-1), FF(0), FF(0) });
    std::vector<Term> lo;
    std::vector<Term> mid;
    for (size_t w = 0; w < LAST_WINDOW; ++w) {
        if (w < SPLIT) {
            lo.push_back({ windows[w], pow2(w * WINDOW_BITS) });
        } else {
            mid.push_back({ windows[w], pow2((w * WINDOW_BITS) - LO) });
        }
    }
    mid.push_back({ top_lo, pow2((LAST_WINDOW * WINDOW_BITS) - LO) });
    const uint32_t lo_var = linear_combination(lo);
    const uint32_t mid_var = linear_combination(mid);
    builder_.create_big_add_gate({ lo_var, mid_var, top_bit, x, FF(1), pow2(LO), pow2(TOP_BIT), FF(-1), FF(0) });
    assert_canonical(lo_var, mid_var, top_bit, LO);
    return windows;
}

const std::vector<plookup::BasicTable*>& ActionCircuitUltraPasta::fixed_base_tables(const AffineElement& base,
                                                                                    size_t num_windows,
                                                                                    size_t top_size)
{
    auto& tables = fixed_base_tables_[{ uint256_t(base.x), num_windows }];
    if (!tables.empty()) {
        return tables;
    }
    // Window w < W: rows (k, [(k + 2) 8^w] B) for k in [0, 8). Top window: rows (k, [k 8^W - offset] B) for k in
    // [0, top_size), offset = sum_{w < W} 2 8^w.
    Element window_base(base);
    Scalar offset(0);
    Scalar eight_pow_w(1);
    for (size_t w = 0; w < num_windows; ++w) {
        std::vector<Builder::TableRow> rows;
        for (size_t k = 0; k < WINDOW_SIZE; ++k) {
            const AffineElement p(window_base * Scalar(k + 2));
            rows.push_back({ FF(k), p.x, p.y });
        }
        tables.push_back(&builder_.create_pasta_table(std::move(rows)));
        offset += eight_pow_w + eight_pow_w;
        eight_pow_w *= Scalar(WINDOW_SIZE);
        window_base = window_base.dbl().dbl().dbl();
    }
    std::vector<Builder::TableRow> rows;
    for (size_t k = 0; k < top_size; ++k) {
        const AffineElement p(Element(base) * ((Scalar(k) * eight_pow_w) - offset));
        BB_ASSERT(!p.is_point_at_infinity());
        rows.push_back({ FF(k), p.x, p.y });
    }
    tables.push_back(&builder_.create_pasta_table(std::move(rows)));
    return tables;
}

ActionCircuitUltraPasta::FixedBaseParts ActionCircuitUltraPasta::fixed_base_mul_parts(
    const AffineElement& base, const std::vector<uint32_t>& windows, size_t top_size)
{
    // Window w < W adds [(k_w + 2) 8^w] B. Before window w the accumulated multiple lies in [2 (8^w - 1) / 7,
    // 9 (8^w - 1) / 7] and the added one in [2 8^w, 9 8^w], so the incomplete additions are exceptional-free while
    // 11 * 8^w < q, i.e. for W <= 84.
    const size_t num_windows = windows.size() - 1;
    BB_ASSERT_LTE(num_windows, LAST_WINDOW);
    const auto& tables = fixed_base_tables(base, num_windows, top_size);
    auto lookup = [&](size_t w) {
        const size_t k = static_cast<size_t>(uint256_t(value(windows[w])).data[0]);
        const auto& row = builder_.pasta_table_rows(tables[w]->table_index).at(k);
        const Point p{ witness(row[1]), witness(row[2]) };
        builder_.create_pasta_lookup(*tables[w], k, windows[w], p.x, p.y);
        return p;
    };
    Point acc = lookup(0);
    for (size_t w = 1; w < num_windows; ++w) {
        acc = incomplete_add(acc, lookup(w));
    }
    return { acc, lookup(num_windows) };
}

ActionCircuitUltraPasta::Point ActionCircuitUltraPasta::fixed_base_mul(const AffineElement& base,
                                                                       const std::vector<uint32_t>& windows)
{
    const auto parts = fixed_base_mul_parts(base, windows, WINDOW_SIZE);
    return complete_add(parts.acc, parts.top);
}

ActionCircuitUltraPasta::Point ActionCircuitUltraPasta::variable_base_mul(const Point& t,
                                                                          const std::vector<uint32_t>& bits)
{
    // Ladder over the m = 254 bits of s >> 1, most significant first: Acc <- (Acc + S) + Acc with S = (2 b - 1) T,
    // from Acc = [2] T, gives [2^254 + 2 (s >> 1) + 1] T. Before step j the multiple a_j lies in [2^j + 1, 3 2^j - 1],
    // so both additions are exceptional-free while 6 * 2^j < q, i.e. for j < 252; the last two steps use complete
    // additions.
    constexpr size_t M = 254;
    constexpr size_t NUM_INCOMPLETE = 252;
    BB_ASSERT_EQ(bits.size(), M + 1);

    // D = [2^254] T and D1 = D + T, for the final correction.
    Point d = t;
    for (size_t i = 0; i < M; ++i) {
        d = dbl(d);
    }
    const Point d1 = incomplete_add(d, t);

    Point acc = dbl(t);
    for (size_t j = 0; j < M; ++j) {
        const uint32_t b = bits[M - j];
        const FF y_s = (value(b) == FF(1)) ? value(t.y) : -value(t.y);
        const Point s{ t.x, witness(y_s) };
        // y_s = (2 b - 1) y_t
        builder_.create_big_mul_add_gate(
            { b, t.y, s.y, builder_.zero_idx(), FF(2), FF(0), FF(-1), FF(-1), FF(0), FF(0) });
        if (j < NUM_INCOMPLETE) {
            const Point r = incomplete_add(acc, s);
            acc = incomplete_add(r, acc);
        } else {
            const Point r = complete_add(acc, s);
            acc = complete_add(r, acc);
        }
    }

    // [s] T = Acc - [2^254 + 1 - b_0] T, with [2^254 + 1 - b_0] T = b_0 ? D : D1.
    const uint32_t b0 = bits[0];
    const bool lsb = value(b0) == FF(1);
    const Point& c = lsb ? d : d1;
    const uint32_t dx = linear_combination({ { d.x, FF(1) }, { d1.x, FF(-1) } });
    const uint32_t dy = linear_combination({ { d.y, FF(1) }, { d1.y, FF(-1) } });
    const Point neg_c{ witness(value(c.x)), witness(-value(c.y)) };
    // x_c = x_d1 + b0 (x_d - x_d1); -y_c = -(y_d1 + b0 (y_d - y_d1))
    builder_.create_big_mul_add_gate({ b0, dx, neg_c.x, d1.x, FF(1), FF(0), FF(0), FF(-1), FF(1), FF(0) });
    builder_.create_big_mul_add_gate({ b0, dy, neg_c.y, d1.y, FF(1), FF(0), FF(0), FF(1), FF(1), FF(0) });
    return complete_add(acc, neg_c);
}

uint32_t ActionCircuitUltraPasta::poseidon_hash(uint32_t a, uint32_t b)
{
    using P = O::Poseidon;
    const auto& params = P::params();
    std::array<uint32_t, 3> s{ a, b, constant(P::constant_length_capacity(2)) };
    // (x + rc)^5 in three gates
    auto sbox = [&](uint32_t x, const FF& rc) {
        const FF t = value(x) + rc;
        const uint32_t t2 = witness(t.sqr());
        builder_.create_arithmetic_gate({ x, x, t2, FF(1), rc + rc, FF(0), FF(-1), rc.sqr() });
        const uint32_t t4 = mul(t2, t2);
        const uint32_t t5 = witness(value(t4) * t);
        builder_.create_arithmetic_gate({ t4, x, t5, FF(1), rc, FF(0), FF(-1), FF(0) });
        return t5;
    };
    size_t r = 0;
    auto full_round = [&]() {
        std::array<uint32_t, 3> y;
        for (size_t j = 0; j < 3; ++j) {
            y[j] = sbox(s[j], params.round_constants[r][j]);
        }
        for (size_t i = 0; i < 3; ++i) {
            s[i] = linear_combination(
                { { y[0], params.mds[i][0] }, { y[1], params.mds[i][1] }, { y[2], params.mds[i][2] } });
        }
        ++r;
    };
    auto partial_round = [&]() {
        const auto& rc = params.round_constants[r];
        const uint32_t y0 = sbox(s[0], rc[0]);
        std::array<uint32_t, 3> out;
        for (size_t i = 0; i < 3; ++i) {
            out[i] =
                linear_combination({ { y0, params.mds[i][0] }, { s[1], params.mds[i][1] }, { s[2], params.mds[i][2] } },
                                   (params.mds[i][1] * rc[1]) + (params.mds[i][2] * rc[2]));
        }
        s = out;
        ++r;
    };
    for (size_t i = 0; i < P::FULL_ROUNDS / 2; ++i) {
        full_round();
    }
    for (size_t i = 0; i < P::PARTIAL_ROUNDS; ++i) {
        partial_round();
    }
    for (size_t i = 0; i < P::FULL_ROUNDS / 2; ++i) {
        full_round();
    }
    return s[0];
}

ActionCircuitUltraPasta::Point ActionCircuitUltraPasta::note_commit(const Point& g_d,
                                                                    const Point& pk_d,
                                                                    uint32_t v,
                                                                    uint32_t rho,
                                                                    uint32_t psi,
                                                                    const std::vector<uint32_t>& rcm_windows)
{
    const auto& k = O::constants();
    const uint32_t g_d_lsb = y_lsb(g_d.y);
    const uint32_t pk_d_lsb = y_lsb(pk_d.y);
    const auto words = message_words({
        { g_d.x, O::L_BASE, true },
        { g_d_lsb, 1 },
        { pk_d.x, O::L_BASE, true },
        { pk_d_lsb, 1 },
        { v, O::L_VALUE },
        { rho, O::L_BASE, true },
        { psi, O::L_BASE, true },
    });
    const Point hash = hash_to_point(k.q_note_commit, words);
    return complete_add(hash, fixed_base_mul(k.note_commit_r, rcm_windows));
}

uint32_t ActionCircuitUltraPasta::commit_ivk(uint32_t ak_x, uint32_t nk, const std::vector<uint32_t>& rivk_windows)
{
    const auto& k = O::constants();
    const auto words = message_words({ { ak_x, O::L_BASE, true }, { nk, O::L_BASE, true } });
    const Point hash = hash_to_point(k.q_commit_ivk, words);
    return complete_add(hash, fixed_base_mul(k.commit_ivk_r, rivk_windows)).x;
}

uint32_t ActionCircuitUltraPasta::merkle_root(uint32_t leaf, const std::array<FF, O::MERKLE_DEPTH>& path, uint32_t pos)
{
    const auto& q = O::constants().q_merkle_crh;
    uint32_t node = leaf;
    for (size_t l = 0; l < O::MERKLE_DEPTH; ++l) {
        const uint32_t layer = constant(FF(l));
        const uint32_t sibling = witness(path[l]);
        const uint32_t right_child = boolean_witness(((pos >> l) & 1) != 0);
        // left = node + b (sibling - node), right = node + sibling - left
        const uint32_t diff = linear_combination({ { sibling, FF(1) }, { node, FF(-1) } });
        const uint32_t left = witness(value(node) + (value(right_child) * value(diff)));
        builder_.create_big_mul_add_gate({ right_child, diff, left, node, FF(1), FF(0), FF(0), FF(-1), FF(1), FF(0) });
        const uint32_t right = linear_combination({ { node, FF(1) }, { sibling, FF(1) }, { left, FF(-1) } });
        // As in halo2's Merkle chip, the 255-bit encodings of the nodes are not required to be canonical: a
        // non-canonical encoding changes the message, which only helps the prover with a Sinsemilla collision.
        const auto words = message_words({ { layer, K, false, true }, { left, O::L_BASE }, { right, O::L_BASE } });
        node = hash_to_point(q, words).x;
    }
    return node;
}

void ActionCircuitUltraPasta::synthesize_action(const Witness& w, std::span<const FF> public_inputs)
{
    using namespace halo2::layout;
    const auto& k = O::constants();

    const uint32_t psi_old = witness(w.psi_old);
    const uint32_t rho_old = witness(w.rho_old);
    const Point cm_old{ witness(w.cm_old.x), witness(w.cm_old.y) };
    const Point g_d_old = witness_point(w.g_d_old);
    const Point ak = witness_point(w.ak);
    const uint32_t nk = witness(w.nk);
    // v_old and v_new are 64-bit segments of the note commitments' messages, which range-constrains them.
    const uint32_t v_old = witness(FF(w.v_old));
    const uint32_t v_new = witness(FF(w.v_new));

    // Public inputs, in instance order.
    std::array<uint32_t, NUM_PUBLIC_INPUTS_PER_ACTION> pis{};
    for (size_t i = 0; i < NUM_PUBLIC_INPUTS_PER_ACTION; ++i) {
        pis[i] = builder_.add_public_variable(public_inputs[i]);
    }
    auto public_input = [&](size_t idx, uint32_t var) {
        builder_.assert_equal(var, pis[idx], "public input " + std::to_string(idx));
    };

    // Merkle path validity
    const uint32_t root = merkle_root(cm_old.x, w.path, w.pos);

    // Value commitment integrity: cv_net = [sign * magnitude] V + [rcv] R, with v_old - v_new = sign * magnitude.
    // [magnitude] V = acc_V + top_V is split so that magnitude = 0 needs no addition of opposite points:
    // cv_net = ([rcv] R + sign * acc_V) + sign * top_V.
    {
        const auto [magnitude_u64, negative] = O::magnitude_sign(w.v_old, w.v_new);
        const auto magnitude_windows = scalar_windows(uint256_t(magnitude_u64), NUM_SHORT_WINDOWS);
        std::vector<Term> terms;
        for (size_t w = 0; w < NUM_SHORT_WINDOWS; ++w) {
            terms.push_back({ magnitude_windows[w], pow2(w * WINDOW_BITS) });
        }
        const uint32_t magnitude = linear_combination(terms);
        const uint32_t sign = witness(negative ? FF(-1) : FF(1));
        builder_.create_arithmetic_gate({ sign, sign, builder_.zero_idx(), FF(1), FF(0), FF(0), FF(0), FF(-1) });
        builder_.create_big_mul_add_gate({ sign, magnitude, v_old, v_new, FF(1), FF(0), FF(0), FF(-1), FF(1), FF(0) });
        // The top window of the magnitude is bit 63, so its table has two rows.
        const auto v_parts = fixed_base_mul_parts(k.value_commit_v, magnitude_windows, 2);
        const Point blind = fixed_base_mul(k.value_commit_r, scalar_windows(uint256_t(w.rcv), NUM_FULL_WINDOWS));
        const Point partial = complete_add(blind, conditional_negate(v_parts.acc, sign));
        const Point cv_net = complete_add(partial, conditional_negate(v_parts.top, sign));
        public_input(CV_NET_X, cv_net.x);
        public_input(CV_NET_Y, cv_net.y);
    }

    // Nullifier integrity: nf = Extract([Poseidon(nk, rho) + psi] K + cm)
    {
        const uint32_t prf = poseidon_hash(nk, rho_old);
        const uint32_t scalar = linear_combination({ { prf, FF(1) }, { psi_old, FF(1) } });
        const Point nf = complete_add(fixed_base_mul(k.nullifier_k, field_windows(scalar)), cm_old);
        public_input(NF_OLD, nf.x);
    }

    // Spend authority: rk = [alpha] G + ak
    {
        const Point rk =
            complete_add(fixed_base_mul(k.spend_auth_g, scalar_windows(uint256_t(w.alpha), NUM_FULL_WINDOWS)), ak);
        public_input(RK_X, rk.x);
        public_input(RK_Y, rk.y);
    }

    // Diversified address integrity: pk_d_old = [CommitIvk(ak, nk)] g_d_old
    const uint32_t ivk = commit_ivk(ak.x, nk, scalar_windows(uint256_t(w.rivk), NUM_FULL_WINDOWS));
    const Point pk_d_old = variable_base_mul(g_d_old, field_bits(ivk));

    // Old note commitment integrity
    {
        const Point cm = note_commit(
            g_d_old, pk_d_old, v_old, rho_old, psi_old, scalar_windows(uint256_t(w.rcm_old), NUM_FULL_WINDOWS));
        assert_equal(cm, cm_old);
    }

    // New note commitment integrity (rho_new = nf_old)
    {
        const Point g_d_new = witness_point(w.g_d_new);
        const Point pk_d_new = witness_point(w.pk_d_new);
        const uint32_t psi_new = witness(w.psi_new);
        const Point cm = note_commit(
            g_d_new, pk_d_new, v_new, pis[NF_OLD], psi_new, scalar_windows(uint256_t(w.rcm_new), NUM_FULL_WINDOWS));
        public_input(CMX, cm.x);
    }

    // Orchard circuit checks: v_old = 0 or root = anchor; v_old = 0 or enable_spend; v_new = 0 or enable_output.
    {
        const uint32_t root_minus_anchor = linear_combination({ { root, FF(1) }, { pis[ANCHOR], FF(-1) } });
        builder_.create_arithmetic_gate(
            { v_old, root_minus_anchor, builder_.zero_idx(), FF(1), FF(0), FF(0), FF(0), FF(0) });
        builder_.create_arithmetic_gate(
            { v_old, pis[ENABLE_SPEND], builder_.zero_idx(), FF(-1), FF(1), FF(0), FF(0), FF(0) });
        builder_.create_arithmetic_gate(
            { v_new, pis[ENABLE_OUTPUT], builder_.zero_idx(), FF(-1), FF(1), FF(0), FF(0), FF(0) });
    }
}

void ActionCircuitUltraPasta::build(Builder& builder,
                                    const std::vector<Witness>& witnesses,
                                    const std::vector<FF>& public_inputs)
{
    using halo2::layout::NUM_PUBLIC_INPUTS_PER_ACTION;
    BB_ASSERT_EQ(public_inputs.size(), witnesses.size() * NUM_PUBLIC_INPUTS_PER_ACTION);
    ActionCircuitUltraPasta circuit(builder);
    for (size_t i = 0; i < witnesses.size(); ++i) {
        circuit.synthesize_action(
            witnesses[i],
            std::span<const FF>(public_inputs).subspan(i * NUM_PUBLIC_INPUTS_PER_ACTION, NUM_PUBLIC_INPUTS_PER_ACTION));
    }
}

} // namespace bb::zcash::ultra_pasta
