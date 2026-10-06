#include "action_circuit_ultra.hpp"
#include "barretenberg/stdlib/primitives/field/field_utils.hpp"
#include "barretenberg/stdlib/primitives/witness/witness.hpp"

namespace bb::zcash::ultra {

namespace {
using witness_ct = stdlib::witness_t<Builder>;
constexpr size_t K = 10;
// Field elements in Sinsemilla messages are encoded on 255 bits; the canonical split is at bit 250 so that the high
// part (bits 250..253, bit 254 being zero for every field element) fits validate_split_in_field.
constexpr size_t CANONICAL_SPLIT = 250;

field_ct constant(const fr& value)
{
    return field_ct(value);
}

field_ct pow_of_two(size_t k)
{
    return constant(fr(uint256_t(1) << k));
}
} // namespace

ActionCircuitUltra::ActionCircuitUltra(Builder& builder)
    : builder_(builder)
{
    const auto& table = Sinsemilla<Cycle>::S_table();
    std::vector<std::array<field_ct, 2>> entries;
    entries.reserve(table.size());
    for (const auto& p : table) {
        entries.push_back({ constant(p.x), constant(p.y) });
    }
    s_table_ = twin_rom_table_ct(entries);
}

std::vector<field_ct> ActionCircuitUltra::message_words(const std::vector<Segment>& segments)
{
    size_t total_bits = 0;
    for (const auto& s : segments) {
        total_bits += s.num_bits;
    }
    const size_t num_words = (total_bits + K - 1) / K;
    std::vector<field_ct> words(num_words, constant(0));

    size_t offset = 0; // global bit offset of the segment
    for (const auto& seg : segments) {
        const uint256_t value(seg.value.get_value());
        // Slice boundaries inside the segment: word boundaries, and the canonical split.
        std::vector<size_t> cuts{ 0 };
        for (size_t b = 1; b < seg.num_bits; ++b) {
            if ((offset + b) % K == 0 || (seg.canonical && b == CANONICAL_SPLIT)) {
                cuts.push_back(b);
            }
        }
        cuts.push_back(seg.num_bits);

        field_ct recomposed = constant(0);
        field_ct lo = constant(0);
        field_ct hi = constant(0);
        for (size_t i = 0; i + 1 < cuts.size(); ++i) {
            const size_t start = cuts[i];
            const size_t width = cuts[i + 1] - start;
            const fr slice_value(value.slice(start, start + width));
            field_ct slice =
                seg.value.is_constant() ? constant(slice_value) : field_ct::from_witness(&builder_, slice_value);
            if (!slice.is_constant()) {
                slice.create_range_constraint(width, "sinsemilla message slice");
            }
            recomposed += slice * pow_of_two(start);
            if (seg.canonical) {
                if (start < CANONICAL_SPLIT) {
                    lo += slice * pow_of_two(start);
                } else {
                    hi += slice * pow_of_two(start - CANONICAL_SPLIT);
                }
            }
            const size_t global = offset + start;
            words[global / K] += slice * pow_of_two(global % K);
        }
        if (!seg.value.is_constant()) {
            recomposed.assert_equal(seg.value, "sinsemilla message recomposition");
            if (seg.canonical) {
                // hi holds bits 250..254; requiring hi < 2^4 forces bit 254 = 0, then (lo, hi) < r.
                hi = hi.normalize();
                hi.create_range_constraint(fr::modulus.get_msb() + 1 - CANONICAL_SPLIT, "canonical high bits");
                stdlib::validate_split_in_field_unsafe(lo.normalize(), hi, CANONICAL_SPLIT, fr::modulus);
            }
        }
        offset += seg.num_bits;
    }
    for (auto& w : words) {
        w = w.normalize();
    }
    return words;
}

cycle_group_ct ActionCircuitUltra::hash_to_point(const AffineElement& q, const std::vector<field_ct>& words)
{
    using Element = Cycle::Element;
    // The accumulator A and the S points are witnesses chained through Ultra elliptic addition gates; as in halo2's
    // Sinsemilla chip the incomplete-addition exceptional cases are not checked (they occur with negligible
    // probability, see the Sinsemilla security analysis).
    uint32_t acc_x = builder_.put_constant_variable(q.x);
    uint32_t acc_y = builder_.put_constant_variable(q.y);
    AffineElement acc = q;
    auto add =
        [&](uint32_t x1, uint32_t y1, const AffineElement& p1, uint32_t x2, uint32_t y2, const AffineElement& p2) {
            const AffineElement p3(Element(p1) + Element(p2));
            const uint32_t x3 = builder_.add_variable(p3.x);
            const uint32_t y3 = builder_.add_variable(p3.y);
            builder_.create_ecc_add_gate(
                { .x1 = x1, .y1 = y1, .x2 = x2, .y2 = y2, .x3 = x3, .y3 = y3, .is_addition = true });
            return std::make_tuple(x3, y3, p3);
        };
    for (const auto& word : words) {
        const auto [sx, sy] = s_table_[word];
        const AffineElement s(sx.get_value(), sy.get_value());
        // A constant word (e.g. a Merkle layer index) reads a constant point.
        auto index = [&](const field_ct& v) {
            return v.is_constant() ? builder_.put_constant_variable(v.get_value()) : v.normalize().get_witness_index();
        };
        const auto [tx, ty, t] = add(acc_x, acc_y, acc, index(sx), index(sy), s);
        std::tie(acc_x, acc_y, acc) = add(tx, ty, t, acc_x, acc_y, acc);
    }
    return cycle_group_ct(field_ct::from_witness_index(&builder_, acc_x),
                          field_ct::from_witness_index(&builder_, acc_y),
                          /*assert_on_curve=*/false);
}

field_ct ActionCircuitUltra::poseidon_hash(const field_ct& a, const field_ct& b)
{
    using P = PoseidonP128Pow5T3<fr>;
    const auto& params = P::params();
    std::array<field_ct, 3> s{ a, b, constant(P::constant_length_capacity(2)) };
    auto sbox = [](const field_ct& x) {
        const field_ct x2 = x.sqr();
        return x2.sqr() * x;
    };
    auto mds = [&]() {
        std::array<field_ct, 3> out;
        for (size_t i = 0; i < 3; ++i) {
            out[i] = (s[0] * params.mds[i][0] + s[1] * params.mds[i][1] + s[2] * params.mds[i][2]).normalize();
        }
        s = out;
    };
    size_t r = 0;
    auto full_round = [&]() {
        for (size_t j = 0; j < 3; ++j) {
            s[j] = sbox(s[j] + params.round_constants[r][j]);
        }
        mds();
        ++r;
    };
    for (size_t i = 0; i < P::FULL_ROUNDS / 2; ++i) {
        full_round();
    }
    for (size_t i = 0; i < P::PARTIAL_ROUNDS; ++i) {
        s[0] = sbox(s[0] + params.round_constants[r][0]);
        s[1] = s[1] + params.round_constants[r][1];
        s[2] = s[2] + params.round_constants[r][2];
        mds();
        ++r;
    }
    for (size_t i = 0; i < P::FULL_ROUNDS / 2; ++i) {
        full_round();
    }
    return s[0];
}

cycle_group_ct ActionCircuitUltra::fixed_base_mul(const AffineElement& base, const cycle_scalar_ct& scalar)
{
    return cycle_group_ct::fixed_batch_mul({ cycle_group_ct(base) }, std::vector<cycle_scalar_ct>{ scalar });
}

cycle_scalar_ct ActionCircuitUltra::base_field_scalar(const field_ct& x)
{
    auto [lo, hi] = stdlib::split_unique(x, cycle_scalar_ct::LO_BITS);
    return cycle_scalar_ct(lo, hi);
}

field_ct ActionCircuitUltra::y_lsb(const field_ct& y)
{
    auto [lsb, rest] = stdlib::split_unique(y, 1);
    return lsb;
}

cycle_group_ct ActionCircuitUltra::note_commit(const cycle_group_ct& g_d,
                                               const cycle_group_ct& pk_d,
                                               const field_ct& value,
                                               const field_ct& rho,
                                               const field_ct& psi,
                                               const cycle_scalar_ct& rcm)
{
    const auto& k = O::constants();
    const auto words = message_words({
        { g_d.x(), O::L_BASE, true },
        { y_lsb(g_d.y()), 1 },
        { pk_d.x(), O::L_BASE, true },
        { y_lsb(pk_d.y()), 1 },
        { value, O::L_VALUE },
        { rho, O::L_BASE, true },
        { psi, O::L_BASE, true },
    });
    const auto hash = hash_to_point(k.q_note_commit, words);
    return hash + fixed_base_mul(k.note_commit_r, rcm);
}

field_ct ActionCircuitUltra::commit_ivk(const field_ct& ak, const field_ct& nk, const cycle_scalar_ct& rivk)
{
    const auto& k = O::constants();
    const auto words = message_words({ { ak, O::L_BASE, true }, { nk, O::L_BASE, true } });
    const auto commitment = hash_to_point(k.q_commit_ivk, words) + fixed_base_mul(k.commit_ivk_r, rivk);
    return commitment.x();
}

field_ct ActionCircuitUltra::merkle_root(const field_ct& leaf,
                                         const std::array<FF, O::MERKLE_DEPTH>& path,
                                         uint32_t pos)
{
    const auto& q = O::constants().q_merkle_crh;
    field_ct node = leaf;
    for (size_t l = 0; l < O::MERKLE_DEPTH; ++l) {
        const field_ct sibling = field_ct::from_witness(&builder_, path[l]);
        const bool_ct right_child(witness_ct(&builder_, ((pos >> l) & 1) != 0));
        const field_ct left = field_ct::conditional_assign(right_child, sibling, node);
        const field_ct right = field_ct::conditional_assign(right_child, node, sibling);
        // As in halo2's Merkle chip, the 255-bit encodings of the nodes are not required to be canonical: a
        // non-canonical encoding changes the message, which only helps the prover with a Sinsemilla collision.
        const auto words = message_words({ { constant(fr(l)), K }, { left, O::L_BASE }, { right, O::L_BASE } });
        node = hash_to_point(q, words).x();
    }
    return node;
}

void ActionCircuitUltra::synthesize_action(const Witness& w, std::span<const FF> public_inputs)
{
    using namespace halo2::layout;
    const auto& k = O::constants();
    Builder* b = &builder_;

    const field_ct psi_old = field_ct::from_witness(b, w.psi_old);
    const field_ct rho_old = field_ct::from_witness(b, w.rho_old);
    const auto cm_old = cycle_group_ct::from_witness(b, w.cm_old);
    const auto g_d_old = cycle_group_ct::from_witness(b, w.g_d_old);
    const auto ak = cycle_group_ct::from_witness(b, w.ak);
    const field_ct nk = field_ct::from_witness(b, w.nk);
    const field_ct v_old = field_ct::from_witness(b, fr(w.v_old));
    const field_ct v_new = field_ct::from_witness(b, fr(w.v_new));
    v_old.create_range_constraint(O::L_VALUE, "v_old");
    v_new.create_range_constraint(O::L_VALUE, "v_new");

    // Public inputs, in instance order.
    std::array<field_ct, NUM_PUBLIC_INPUTS_PER_ACTION> pis;
    for (size_t i = 0; i < NUM_PUBLIC_INPUTS_PER_ACTION; ++i) {
        pis[i] = field_ct::from_witness(b, public_inputs[i]);
        pis[i].set_public();
    }
    auto public_input = [&](size_t idx, const field_ct& value) {
        value.assert_equal(pis[idx], "public input " + std::to_string(idx));
    };

    // Merkle path validity
    const field_ct root = merkle_root(cm_old.x(), w.path, w.pos);

    // Value commitment integrity: cv_net = [v_old - v_new] V + [rcv] R
    const auto [magnitude_u64, negative] = O::magnitude_sign(w.v_old, w.v_new);
    const field_ct magnitude = field_ct::from_witness(b, fr(magnitude_u64));
    magnitude.create_range_constraint(O::L_VALUE, "magnitude");
    const field_ct sign = field_ct::from_witness(b, negative ? -fr(1) : fr(1));
    (sign * sign).assert_equal(constant(1), "sign is +-1");
    {
        const auto magnitude_mul = fixed_base_mul(k.value_commit_v, cycle_scalar_ct(magnitude, constant(0)));
        const cycle_group_ct commitment(magnitude_mul.x(), magnitude_mul.y() * sign, /*assert_on_curve=*/false);
        const auto blind = fixed_base_mul(k.value_commit_r, cycle_scalar_ct::from_witness(b, w.rcv));
        const auto cv_net = commitment + blind;
        public_input(CV_NET_X, cv_net.x());
        public_input(CV_NET_Y, cv_net.y());
    }

    // Nullifier integrity: nf = Extract([Poseidon(nk, rho) + psi] K + cm)
    {
        const field_ct scalar = poseidon_hash(nk, rho_old) + psi_old;
        const auto nf = fixed_base_mul(k.nullifier_k, base_field_scalar(scalar)) + cm_old;
        public_input(NF_OLD, nf.x());
    }

    // Spend authority: rk = [alpha] G + ak
    {
        const auto rk = fixed_base_mul(k.spend_auth_g, cycle_scalar_ct::from_witness(b, w.alpha)) + ak;
        public_input(RK_X, rk.x());
        public_input(RK_Y, rk.y());
    }

    // Diversified address integrity: pk_d_old = [CommitIvk(ak, nk)] g_d_old
    const auto pk_d_old = cycle_group_ct::from_witness(b, w.pk_d_old);
    {
        const field_ct ivk = commit_ivk(ak.x(), nk, cycle_scalar_ct::from_witness(b, w.rivk));
        auto derived = cycle_group_ct::batch_mul({ g_d_old }, std::vector<cycle_scalar_ct>{ base_field_scalar(ivk) });
        derived.x().assert_equal(pk_d_old.x(), "pk_d_old x");
        derived.y().assert_equal(pk_d_old.y(), "pk_d_old y");
    }

    // Old note commitment integrity
    {
        const auto cm =
            note_commit(g_d_old, pk_d_old, v_old, rho_old, psi_old, cycle_scalar_ct::from_witness(b, w.rcm_old));
        cm.x().assert_equal(cm_old.x(), "cm_old x");
        cm.y().assert_equal(cm_old.y(), "cm_old y");
    }

    // New note commitment integrity (rho_new = nf_old)
    {
        const auto g_d_new = cycle_group_ct::from_witness(b, w.g_d_new);
        const auto pk_d_new = cycle_group_ct::from_witness(b, w.pk_d_new);
        const field_ct psi_new = field_ct::from_witness(b, w.psi_new);
        const field_ct& nf_old = pis[NF_OLD];
        const auto cm =
            note_commit(g_d_new, pk_d_new, v_new, nf_old, psi_new, cycle_scalar_ct::from_witness(b, w.rcm_new));
        public_input(CMX, cm.x());
    }

    // Orchard circuit checks
    {
        const field_ct& anchor = pis[ANCHOR];
        const field_ct& enable_spend = pis[ENABLE_SPEND];
        const field_ct& enable_output = pis[ENABLE_OUTPUT];
        (v_old - v_new).assert_equal(magnitude * sign, "v_old - v_new = magnitude * sign");
        (v_old * (root - anchor)).assert_equal(constant(0), "v_old = 0 or root = anchor");
        (v_old * (constant(1) - enable_spend)).assert_equal(constant(0), "v_old = 0 or enable_spend");
        (v_new * (constant(1) - enable_output)).assert_equal(constant(0), "v_new = 0 or enable_output");
    }
}

void ActionCircuitUltra::build(Builder& builder,
                               const std::vector<Witness>& witnesses,
                               const std::vector<FF>& public_inputs)
{
    BB_ASSERT_EQ(public_inputs.size(), witnesses.size() * halo2::layout::NUM_PUBLIC_INPUTS_PER_ACTION);
    ActionCircuitUltra circuit(builder);
    for (size_t i = 0; i < witnesses.size(); ++i) {
        circuit.synthesize_action(
            witnesses[i],
            std::span<const FF>(public_inputs)
                .subspan(i * halo2::layout::NUM_PUBLIC_INPUTS_PER_ACTION, halo2::layout::NUM_PUBLIC_INPUTS_PER_ACTION));
    }
}

} // namespace bb::zcash::ultra
