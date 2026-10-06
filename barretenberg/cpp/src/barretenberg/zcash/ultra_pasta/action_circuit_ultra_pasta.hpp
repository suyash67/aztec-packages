#pragma once

#include "barretenberg/zcash/halo2/layout.hpp"
#include "barretenberg/zcash/primitives/orchard.hpp"
#include "barretenberg/zcash/ultra_pasta/ultra_pasta_builder.hpp"

#include <map>
#include <span>
#include <vector>

/**
 * @file action_circuit_ultra_pasta.hpp
 * @brief The Orchard Action statement over the Pallas base field as an UltraHonk circuit written directly with Ultra
 * gates (arithmetic, elliptic, lookup, delta-range), without barretenberg's stdlib.
 *
 * @details The statement is Orchard's (Orchard<PastaCycle> is the native reference); the gadgets are:
 *
 *   - Sinsemilla: the S table is a lookup table of rows (m, x(S[m]), y(S[m])). A hash step is one lookup (which also
 *     range-checks the 10-bit word m) and Acc <- (Acc + S[m]) + Acc with two chained Ultra elliptic addition gates
 *     (one row each). As in halo2's Sinsemilla chip the incomplete-addition exceptional cases are not excluded
 *     (negligible, see the Sinsemilla security argument).
 *   - message decomposition: 10-bit words are linear combinations of slices of the message fields; slices that are
 *     not a whole word are range-constrained. Field elements that must be encoded canonically (NoteCommit, CommitIvk)
 *     are split at bits 130 and 254 and checked against p = 2^254 + t_p as in halo2: if bit 254 is set, bits 130..253
 *     are zero and the low 130 bits are below t_p.
 *   - fixed-base scalar multiplication: halo2's 3-bit windows. Window w < W looks up [(k_w + 2) 8^w] B in a table of
 *     the base (which range-checks k_w) and adds it with one elliptic gate; the offsets keep the accumulator away from
 *     the exceptional cases of incomplete addition. The last window looks up [k_W 8^W - offsets] B, added with a
 *     complete addition. Full-width scalars use 85 windows, the 64-bit value magnitude 22 (the last one being bit 63).
 *   - variable-base scalar multiplication [s] T: the bits of s >> 1 drive the ladder Acc <- (Acc + (+-T)) + Acc from
 *     Acc = [2] T (two chained elliptic gates per bit; the last two bits use complete additions), which yields
 *     [2^254 + 2 (s >> 1) + 1] T; [2^254 + 1 - lsb(s)] T is then subtracted with a complete addition.
 *   - complete addition: an arithmetic-gate gadget handling P = Q (13 gates); P = -Q is unsatisfiable.
 *   - Poseidon P128Pow5T3: x^5 S-boxes in three multiplication gates, the MDS mix in one gate per state element.
 * Public inputs (10 per Action, in instance order) are the circuit's public inputs.
 */
namespace bb::zcash::ultra_pasta {

using Builder = UltraPastaCircuitBuilder;

class ActionCircuitUltraPasta {
  public:
    using Cycle = PastaCycle;
    using O = Orchard<Cycle>;
    using FF = Cycle::FF;
    using Scalar = Cycle::EmbeddedScalar;
    using AffineElement = Cycle::AffineElement;
    using Element = Cycle::Element;
    using Witness = O::ActionWitness;

    // A Pallas point in the circuit (never the identity).
    struct Point {
        uint32_t x;
        uint32_t y;
    };
    struct Term {
        uint32_t var;
        FF coeff;
    };
    // A field of a Sinsemilla message: the `num_bits` low bits of `value`. `canonical` requires value < p as an integer
    // (for a full field element encoded on 255 bits); `is_constant` marks a constant value (its slices are constants).
    struct Segment {
        uint32_t value;
        size_t num_bits;
        bool canonical = false;
        bool is_constant = false;
    };

    /**
     * @brief Adds every Action to the builder; the public inputs are added in order (10 per Action).
     */
    static void build(Builder& builder, const std::vector<Witness>& witnesses, const std::vector<FF>& public_inputs);

    explicit ActionCircuitUltraPasta(Builder& builder);

    void synthesize_action(const Witness& w, std::span<const FF> public_inputs);

    // Arithmetic
    FF value(uint32_t var) const { return builder_.get_variable(var); }
    uint32_t witness(const FF& value) { return builder_.add_variable(value); }
    uint32_t constant(const FF& value) { return builder_.put_constant_variable(value); }
    // A new variable equal to sum_i coeff_i var_i + constant_term (three terms per gate, chained through w_4).
    uint32_t linear_combination(const std::vector<Term>& terms, const FF& constant_term = FF(0));
    uint32_t mul(uint32_t a, uint32_t b);
    void range_constrain(uint32_t var, size_t num_bits);
    uint32_t boolean_witness(bool bit);

    // Pallas arithmetic
    AffineElement point_value(const Point& p) const { return AffineElement(value(p.x), value(p.y)); }
    Point constant_point(const AffineElement& p) { return { constant(p.x), constant(p.y) }; }
    Point witness_point(const AffineElement& p);
    Point incomplete_add(const Point& p, const Point& q);
    Point complete_add(const Point& p, const Point& q);
    Point dbl(const Point& p);
    // (x, sign * y) for sign in {1, -1}.
    Point conditional_negate(const Point& p, uint32_t sign);
    void assert_equal(const Point& p, const Point& q);

    // Sinsemilla
    std::vector<uint32_t> message_words(const std::vector<Segment>& segments);
    Point hash_to_point(const AffineElement& q, const std::vector<uint32_t>& words);
    // The 255-bit integer lo + 2^lo_bits mid + 2^254 top (lo < 2^lo_bits, mid < 2^(254 - lo_bits), top a bit, and
    // 126 <= lo_bits) is below p.
    void assert_canonical(uint32_t lo, uint32_t mid, uint32_t top, size_t lo_bits);
    // The least-significant bit of the canonical representative of y.
    uint32_t y_lsb(uint32_t y);

    // Scalar multiplication. Bits and 3-bit windows are little-endian.
    std::vector<uint32_t> scalar_bits(const uint256_t& k, size_t num_bits);
    // The 255 bits of the canonical representative of x.
    std::vector<uint32_t> field_bits(uint32_t x);
    // Window witnesses of k (range-constrained by the fixed-base table lookups they are used in).
    std::vector<uint32_t> scalar_windows(const uint256_t& k, size_t num_windows);
    // The 85 windows of the canonical representative of x.
    std::vector<uint32_t> field_windows(uint32_t x);
    // The window tables of a fixed base (created on first use).
    const std::vector<plookup::BasicTable*>& fixed_base_tables(const AffineElement& base,
                                                               size_t num_windows,
                                                               size_t top_size);
    // [k] B = acc + top, with acc the sum over all windows but the last (offsets included) and top the last window's
    // point (its table has top_size rows).
    struct FixedBaseParts {
        Point acc;
        Point top;
    };
    FixedBaseParts fixed_base_mul_parts(const AffineElement& base,
                                        const std::vector<uint32_t>& windows,
                                        size_t top_size);
    Point fixed_base_mul(const AffineElement& base, const std::vector<uint32_t>& windows);
    Point variable_base_mul(const Point& t, const std::vector<uint32_t>& bits);

    uint32_t poseidon_hash(uint32_t a, uint32_t b);
    Point note_commit(const Point& g_d,
                      const Point& pk_d,
                      uint32_t v,
                      uint32_t rho,
                      uint32_t psi,
                      const std::vector<uint32_t>& rcm_windows);
    uint32_t commit_ivk(uint32_t ak_x, uint32_t nk, const std::vector<uint32_t>& rivk_windows);
    uint32_t merkle_root(uint32_t leaf, const std::array<FF, O::MERKLE_DEPTH>& path, uint32_t pos);

  private:
    Builder& builder_;
    plookup::BasicTable* s_table_ = nullptr;
    std::map<std::pair<uint256_t, size_t>, std::vector<plookup::BasicTable*>> fixed_base_tables_;
};

} // namespace bb::zcash::ultra_pasta
