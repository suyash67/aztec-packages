#pragma once

#include "barretenberg/zcash/halo2/chips.hpp"

#include <vector>

/**
 * @file action_circuit.hpp
 * @brief The Orchard Action circuit (orchard 0.16.0 `Circuit::synthesize`) on the halo2-style builder.
 *
 * @details `synthesize_action` assigns one Action, in the same order as orchard's `synthesize_base`:
 *   witness shared values -> Merkle path -> value commitment -> nullifier -> spend authority ->
 *   CommitIvk + [ivk] g_d_old -> old NoteCommit -> new NoteCommit -> checks gate (-> NU6.3 cross-address rows).
 * Several Actions are synthesized into one builder (one proof covers the bundle); their public inputs occupy
 * consecutive blocks of 10 instance rows.
 */
namespace bb::zcash::halo2 {

enum class CircuitVersion {
    // orchard FixedPostNu6_2: anchored variable-base mul, no cross-address rows (what orchard's benchmarks build).
    FixedPostNu6_2,
    // orchard PostNu6_3: FixedPostNu6_2 plus the 4 cross-address rows enforcing instance row 9.
    PostNu6_3,
};

template <typename Cycle> class ActionCircuit {
  public:
    using C = Chips<Cycle>;
    using FF = typename C::FF;
    using Scalar = typename C::Scalar;
    using AffineElement = typename C::AffineElement;
    using Cell = typename C::Cell;
    using Point = typename C::Point;
    using Region = typename C::Region;
    using MessagePiece = typename C::MessagePiece;
    using RangeConstrainedCell = typename C::RangeConstrainedCell;
    using RangeConstrainedValue = typename C::RangeConstrainedValue;
    using SinsemillaConfig = typename C::SinsemillaConfig;
    using O = Orchard<Cycle>;
    using Witness = typename O::ActionWitness;
    using Builder = PlonkishBuilder<FF>;

    static typename Builder::Config builder_config()
    {
        return { .num_advice = NUM_ADVICE,
                 .num_fixed = NUM_FIXED,
                 .num_selectors = NUM_SELECTORS,
                 .constants_column = F0,
                 .public_input_column = 0 };
    }

    /**
     * @brief Builds the halo2-style table for `witnesses.size()` Actions.
     * @param public_inputs the expected public inputs (10 per Action), e.g. from Orchard<Cycle>::evaluate.
     */
    static typename Builder::Table build(const std::vector<Witness>& witnesses,
                                         const std::vector<FF>& public_inputs,
                                         CircuitVersion version = CircuitVersion::FixedPostNu6_2)
    {
        BB_ASSERT_EQ(public_inputs.size(), witnesses.size() * NUM_PUBLIC_INPUTS_PER_ACTION);
        Builder builder(builder_config());
        builder.set_public_inputs(public_inputs);
        C chips(builder, /*anchored_base=*/true);
        for (size_t i = 0; i < witnesses.size(); ++i) {
            synthesize_action(chips, witnesses[i], i * NUM_PUBLIC_INPUTS_PER_ACTION, version);
        }
        return builder.finalize();
    }

    static void synthesize_action(C& c, const Witness& w, size_t pi_offset, CircuitVersion version)
    {
        Builder& b = c.builder();
        const auto& k = O::constants();

        // Witness private inputs that are used across multiple checks.
        Cell psi_old = c.assign_free_advice(0, w.psi_old);
        Cell rho_old = c.assign_free_advice(0, w.rho_old);
        Point cm_old = c.witness_point(w.cm_old);
        Point g_d_old = c.witness_point_non_id(w.g_d_old);
        Point ak_P = c.witness_point_non_id(w.ak);
        Cell nk = c.assign_free_advice(0, w.nk);
        Cell v_old = c.assign_free_advice(0, FF(w.v_old));
        Cell v_new = c.assign_free_advice(0, FF(w.v_new));

        // Merkle path validity.
        Cell root = c.merkle_root(cm_old.x, w.path, w.pos);

        // Value commitment integrity.
        const auto [magnitude_u64, negative] = O::magnitude_sign(w.v_old, w.v_new);
        Cell magnitude = c.assign_free_advice(9, FF(magnitude_u64));
        Cell sign = c.assign_free_advice(9, negative ? -FF(1) : FF(1));
        {
            Point commitment = c.mul_fixed_short(magnitude, sign, k.fb_value_commit_v);
            Point blind = c.mul_fixed_full(w.rcv, k.fb_value_commit_r);
            Point cv_net = c.add(commitment, blind);
            b.constrain_instance(cv_net.x, pi_offset + CV_NET_X);
            b.constrain_instance(cv_net.y, pi_offset + CV_NET_Y);
        }

        // Nullifier integrity.
        Cell nf_old = [&]() {
            Cell hash = c.poseidon_hash(nk, rho_old);
            Cell scalar = c.add_field(hash, psi_old);
            Point product = c.mul_fixed_base_field(scalar, k.fb_nullifier_k);
            Point nf = c.add(cm_old, product);
            b.constrain_instance(nf.x, pi_offset + NF_OLD);
            return nf.x;
        }();

        // Spend authority.
        {
            Point alpha_commitment = c.mul_fixed_full(w.alpha, k.fb_spend_auth_g);
            Point rk = c.add(alpha_commitment, ak_P);
            b.constrain_instance(rk.x, pi_offset + RK_X);
            b.constrain_instance(rk.y, pi_offset + RK_Y);
        }

        // Diversified address integrity.
        Point pk_d_old = [&]() {
            Cell ivk = commit_ivk(c, ak_P.x, nk, w.rivk);
            Point derived_pk_d_old = c.mul_var(ivk, g_d_old);
            Point pk_d = c.witness_point_non_id(w.pk_d_old);
            c.constrain_equal(derived_pk_d_old, pk_d);
            return pk_d;
        }();

        // Old note commitment integrity.
        {
            Point derived_cm_old =
                note_commit(c, C::sinsemilla_config_1(), g_d_old, pk_d_old, v_old, rho_old, psi_old, w.rcm_old);
            c.constrain_equal(derived_cm_old, cm_old);
        }

        Point g_d_new = c.witness_point_non_id(w.g_d_new);
        Point pk_d_new = c.witness_point_non_id(w.pk_d_new);

        // New note commitment integrity.
        {
            Cell psi_new = c.assign_free_advice(0, w.psi_new);
            Point cm_new =
                note_commit(c, C::sinsemilla_config_2(), g_d_new, pk_d_new, v_new, nf_old, psi_new, w.rcm_new);
            b.constrain_instance(cm_new.x, pi_offset + CMX);
        }

        // Orchard circuit checks gate.
        b.assign_region("Orchard circuit checks", [&](Region& r) {
            r.copy_advice(v_old, advice(0), 0);
            r.copy_advice(v_new, advice(1), 0);
            r.copy_advice(magnitude, advice(2), 0);
            r.copy_advice(sign, advice(3), 0);
            r.copy_advice(root, advice(4), 0);
            r.assign_advice_from_instance(pi_offset + ANCHOR, advice(5), 0);
            r.assign_advice_from_instance(pi_offset + ENABLE_SPEND, advice(6), 0);
            r.assign_advice_from_instance(pi_offset + ENABLE_OUTPUT, advice(7), 0);
            r.enable_selector(selector(Q_ORCHARD), 0);
            return 0;
        });

        if (version == CircuitVersion::PostNu6_3) {
            b.assign_region("post-NU 6.3 cross-address checks", [&](Region& r) {
                const std::array<std::pair<Cell, Cell>, 4> checks = { std::make_pair(g_d_old.x, g_d_new.x),
                                                                      std::make_pair(g_d_old.y, g_d_new.y),
                                                                      std::make_pair(pk_d_old.x, pk_d_new.x),
                                                                      std::make_pair(pk_d_old.y, pk_d_new.y) };
                for (size_t offset = 0; offset < checks.size(); ++offset) {
                    Cell dca = r.assign_advice_from_instance(pi_offset + DISABLE_CROSS_ADDRESS, advice(0), offset);
                    r.assign_advice_from_constant(advice(1), offset, FF(0));
                    r.copy_advice(dca, advice(2), offset);
                    r.assign_advice_from_constant(advice(3), offset, FF(1));
                    r.copy_advice(checks[offset].first, advice(4), offset);
                    r.copy_advice(checks[offset].second, advice(5), offset);
                    r.assign_advice_from_constant(advice(6), offset, FF(1));
                    r.assign_advice_from_constant(advice(7), offset, FF(1));
                    r.copy_advice(dca, advice(8), offset);
                    r.copy_advice(dca, advice(9), offset);
                    r.enable_selector(selector(Q_ORCHARD), offset);
                }
                return 0;
            });
        }
    }

    // ------------------------------------------------------------------ CommitIvk (commit_ivk.rs)

    static std::pair<Cell, Cell> canonicity_check(C& c, const FF& value, size_t num_words)
    {
        const auto zs = c.witness_check(value, num_words, false);
        return { zs[0], zs[num_words] };
    }

    static FF two_pow(size_t k) { return C::two_pow(k); }
    static FF t_p() { return FF(C::t_p()); }

    static Cell commit_ivk(C& c, const Cell& ak, const Cell& nk, const Scalar& rivk)
    {
        Builder& b = c.builder();
        const auto cfg = C::sinsemilla_config_1();
        MessagePiece a = c.from_subpieces(cfg, { C::bitrange_of(ak.value, 0, 250) });
        RangeConstrainedCell b_0 = c.witness_short(ak.value, 250, 254);
        RangeConstrainedValue b_1 = C::bitrange_of(ak.value, 254, 255);
        RangeConstrainedCell b_2 = c.witness_short(nk.value, 0, 5);
        MessagePiece bp = c.from_subpieces(cfg, { { b_0.cell.value, 4 }, b_1, { b_2.cell.value, 5 } });
        MessagePiece cp = c.from_subpieces(cfg, { C::bitrange_of(nk.value, 5, 245) });
        RangeConstrainedCell d_0 = c.witness_short(nk.value, 245, 254);
        RangeConstrainedValue d_1 = C::bitrange_of(nk.value, 254, 255);
        MessagePiece dp = c.from_subpieces(cfg, { { d_0.cell.value, 9 }, d_1 });

        const auto& k = O::constants();
        auto [commitment, zs] = c.commit(cfg, k.q_commit_ivk, k.fb_commit_ivk_r, { a, bp, cp, dp }, rivk);
        const Cell ivk = commitment.x;
        const Cell z13_a = zs[0][13];
        const Cell z13_c = zs[2][13];

        auto [a_prime, z13_a_prime] = canonicity_check(c, a.cell.value + two_pow(130) - t_p(), 13);
        auto [b2_c_prime, z14_b2_c_prime] =
            canonicity_check(c, b_2.cell.value + cp.cell.value * two_pow(5) + two_pow(140) - t_p(), 14);

        b.assign_region("Assign cells used in canonicity gate", [&](Region& r) {
            r.enable_selector(selector(Q_COMMIT_IVK), 0);
            r.copy_advice(ak, advice(0), 0);
            r.copy_advice(a.cell, advice(1), 0);
            r.copy_advice(bp.cell, advice(2), 0);
            r.copy_advice(b_0.cell, advice(3), 0);
            r.assign_advice(advice(4), 0, b_1.value);
            r.copy_advice(b_2.cell, advice(5), 0);
            r.copy_advice(z13_a, advice(6), 0);
            r.copy_advice(a_prime, advice(7), 0);
            r.copy_advice(z13_a_prime, advice(8), 0);

            r.copy_advice(nk, advice(0), 1);
            r.copy_advice(cp.cell, advice(1), 1);
            r.copy_advice(dp.cell, advice(2), 1);
            r.copy_advice(d_0.cell, advice(3), 1);
            r.assign_advice(advice(4), 1, d_1.value);
            r.copy_advice(z13_c, advice(6), 1);
            r.copy_advice(b2_c_prime, advice(7), 1);
            r.copy_advice(z14_b2_c_prime, advice(8), 1);
            return 0;
        });
        return ivk;
    }

    // ------------------------------------------------------------------ NoteCommit (note_commit.rs)
    // col_l = a6, col_m = a7, col_r = a8, col_z = a9; y canonicity uses a5..a9.
    static constexpr size_t COL_L = 6;
    static constexpr size_t COL_M = 7;
    static constexpr size_t COL_R = 8;
    static constexpr size_t COL_Z = 9;

    static RangeConstrainedCell y_canonicity(C& c, const Cell& y, const RangeConstrainedValue& lsb)
    {
        Builder& b = c.builder();
        RangeConstrainedCell k_0 = c.witness_short(y.value, 1, 10);
        RangeConstrainedValue k_1 = C::bitrange_of(y.value, 10, 250);
        RangeConstrainedCell k_2 = c.witness_short(y.value, 250, 254);
        RangeConstrainedValue k_3 = C::bitrange_of(y.value, 254, 255);
        const FF j_val = lsb.value + FF(2) * k_0.cell.value + two_pow(10) * k_1.value;
        const auto zs = c.witness_check(j_val, 25, true);
        const Cell j = zs[0];
        const Cell z1_j = zs[1];
        const Cell z13_j = zs[13];
        auto [j_prime, z13_j_prime] = canonicity_check(c, j.value + two_pow(130) - t_p(), 13);
        return b.assign_region("y canonicity", [&](Region& r) {
            r.enable_selector(selector(Q_Y_CANON), 0);
            r.copy_advice(y, advice(5), 0);
            Cell lsb_cell = r.assign_advice(advice(6), 0, lsb.value);
            r.copy_advice(k_0.cell, advice(7), 0);
            r.copy_advice(k_2.cell, advice(8), 0);
            r.assign_advice(advice(9), 0, k_3.value);
            r.copy_advice(j, advice(5), 1);
            r.copy_advice(z1_j, advice(6), 1);
            r.copy_advice(z13_j, advice(7), 1);
            r.copy_advice(j_prime, advice(8), 1);
            r.copy_advice(z13_j_prime, advice(9), 1);
            return RangeConstrainedCell{ lsb_cell, lsb.num_bits };
        });
    }

    static Point note_commit(C& c,
                             const SinsemillaConfig& cfg,
                             const Point& g_d,
                             const Point& pk_d,
                             const Cell& value,
                             const Cell& rho,
                             const Cell& psi,
                             const Scalar& rcm)
    {
        Builder& b = c.builder();
        auto rv = [](const RangeConstrainedCell& x) { return RangeConstrainedValue{ x.cell.value, x.num_bits }; };

        MessagePiece a = c.from_subpieces(cfg, { C::bitrange_of(g_d.x.value, 0, 250) });
        // DecomposeB
        RangeConstrainedCell b_0 = c.witness_short(g_d.x.value, 250, 254);
        RangeConstrainedValue b_1 = C::bitrange_of(g_d.x.value, 254, 255);
        RangeConstrainedValue b_2 = C::bitrange_of(g_d.y.value, 0, 1);
        RangeConstrainedCell b_3 = c.witness_short(pk_d.x.value, 0, 4);
        MessagePiece bp = c.from_subpieces(cfg, { rv(b_0), b_1, b_2, rv(b_3) });
        MessagePiece cp = c.from_subpieces(cfg, { C::bitrange_of(pk_d.x.value, 4, 254) });
        // DecomposeD
        RangeConstrainedValue d_0 = C::bitrange_of(pk_d.x.value, 254, 255);
        RangeConstrainedValue d_1 = C::bitrange_of(pk_d.y.value, 0, 1);
        RangeConstrainedCell d_2 = c.witness_short(value.value, 0, 8);
        RangeConstrainedValue d_3 = C::bitrange_of(value.value, 8, 58);
        MessagePiece dp = c.from_subpieces(cfg, { d_0, d_1, rv(d_2), d_3 });
        // DecomposeE
        RangeConstrainedCell e_0 = c.witness_short(value.value, 58, 64);
        RangeConstrainedCell e_1 = c.witness_short(rho.value, 0, 4);
        MessagePiece ep = c.from_subpieces(cfg, { rv(e_0), rv(e_1) });
        MessagePiece fp = c.from_subpieces(cfg, { C::bitrange_of(rho.value, 4, 254) });
        // DecomposeG
        RangeConstrainedValue g_0 = C::bitrange_of(rho.value, 254, 255);
        RangeConstrainedCell g_1 = c.witness_short(psi.value, 0, 9);
        RangeConstrainedValue g_2 = C::bitrange_of(psi.value, 9, 249);
        MessagePiece gp = c.from_subpieces(cfg, { g_0, rv(g_1), g_2 });
        // DecomposeH
        RangeConstrainedCell h_0 = c.witness_short(psi.value, 249, 254);
        RangeConstrainedValue h_1 = C::bitrange_of(psi.value, 254, 255);
        MessagePiece hp = c.from_subpieces(cfg, { rv(h_0), h_1, { FF(0), 4 } });

        RangeConstrainedCell b_2_cell = y_canonicity(c, g_d.y, b_2);
        RangeConstrainedCell d_1_cell = y_canonicity(c, pk_d.y, d_1);

        const auto& k = O::constants();
        auto [cm, zs] = c.commit(cfg, k.q_note_commit, k.fb_note_commit_r, { a, bp, cp, dp, ep, fp, gp, hp }, rcm);
        const Cell z13_a = zs[0][13];
        const Cell z13_c = zs[2][13];
        const Cell z1_d = zs[3][1];
        const Cell z13_f = zs[5][13];
        const Cell z1_g = zs[6][1];
        const Cell g_2_cell = z1_g;
        const Cell z13_g = zs[6][13];

        auto [a_prime, z13_a_prime] = canonicity_check(c, a.cell.value + two_pow(130) - t_p(), 13);
        auto [b3_c_prime, z14_b3_c_prime] =
            canonicity_check(c, b_3.cell.value + two_pow(4) * cp.cell.value + two_pow(140) - t_p(), 14);
        auto [e1_f_prime, z14_e1_f_prime] =
            canonicity_check(c, e_1.cell.value + two_pow(4) * fp.cell.value + two_pow(140) - t_p(), 14);
        auto [g1_g2_prime, z13_g1_g2_prime] =
            canonicity_check(c, g_1.cell.value + two_pow(9) * g_2_cell.value + two_pow(130) - t_p(), 13);

        // DecomposeB::assign
        Cell b_1_cell = b.assign_region("NoteCommit MessagePiece b", [&](Region& r) {
            r.enable_selector(selector(Q_NOTECOMMIT_B), 0);
            r.copy_advice(bp.cell, advice(COL_L), 0);
            r.copy_advice(b_0.cell, advice(COL_M), 0);
            Cell b_1c = r.assign_advice(advice(COL_R), 0, b_1.value);
            r.copy_advice(b_2_cell.cell, advice(COL_M), 1);
            r.copy_advice(b_3.cell, advice(COL_R), 1);
            return b_1c;
        });
        // DecomposeD::assign
        Cell d_0_cell = b.assign_region("NoteCommit MessagePiece d", [&](Region& r) {
            r.enable_selector(selector(Q_NOTECOMMIT_D), 0);
            r.copy_advice(dp.cell, advice(COL_L), 0);
            Cell d_0c = r.assign_advice(advice(COL_M), 0, d_0.value);
            r.copy_advice(d_1_cell.cell, advice(COL_R), 0);
            r.copy_advice(d_2.cell, advice(COL_M), 1);
            r.copy_advice(z1_d, advice(COL_R), 1);
            return d_0c;
        });
        // DecomposeE::assign
        b.assign_region("NoteCommit MessagePiece e", [&](Region& r) {
            r.enable_selector(selector(Q_NOTECOMMIT_E), 0);
            r.copy_advice(ep.cell, advice(COL_L), 0);
            r.copy_advice(e_0.cell, advice(COL_M), 0);
            r.copy_advice(e_1.cell, advice(COL_R), 0);
            return 0;
        });
        // DecomposeG::assign
        Cell g_0_cell = b.assign_region("NoteCommit MessagePiece g", [&](Region& r) {
            r.enable_selector(selector(Q_NOTECOMMIT_G), 0);
            r.copy_advice(gp.cell, advice(COL_L), 0);
            Cell g_0c = r.assign_advice(advice(COL_M), 0, g_0.value);
            r.copy_advice(g_1.cell, advice(COL_L), 1);
            r.copy_advice(z1_g, advice(COL_M), 1);
            return g_0c;
        });
        // DecomposeH::assign
        Cell h_1_cell = b.assign_region("NoteCommit MessagePiece h", [&](Region& r) {
            r.enable_selector(selector(Q_NOTECOMMIT_H), 0);
            r.copy_advice(hp.cell, advice(COL_L), 0);
            r.copy_advice(h_0.cell, advice(COL_M), 0);
            return r.assign_advice(advice(COL_R), 0, h_1.value);
        });
        // GdCanonicity::assign
        b.assign_region("NoteCommit input g_d", [&](Region& r) {
            r.copy_advice(g_d.x, advice(COL_L), 0);
            r.copy_advice(b_0.cell, advice(COL_M), 0);
            r.copy_advice(b_1_cell, advice(COL_M), 1);
            r.copy_advice(a.cell, advice(COL_R), 0);
            r.copy_advice(a_prime, advice(COL_R), 1);
            r.copy_advice(z13_a, advice(COL_Z), 0);
            r.copy_advice(z13_a_prime, advice(COL_Z), 1);
            r.enable_selector(selector(Q_NOTECOMMIT_G_D), 0);
            return 0;
        });
        // PkdCanonicity::assign
        b.assign_region("NoteCommit input pk_d", [&](Region& r) {
            r.copy_advice(pk_d.x, advice(COL_L), 0);
            r.copy_advice(b_3.cell, advice(COL_M), 0);
            r.copy_advice(d_0_cell, advice(COL_M), 1);
            r.copy_advice(cp.cell, advice(COL_R), 0);
            r.copy_advice(b3_c_prime, advice(COL_R), 1);
            r.copy_advice(z13_c, advice(COL_Z), 0);
            r.copy_advice(z14_b3_c_prime, advice(COL_Z), 1);
            r.enable_selector(selector(Q_NOTECOMMIT_PK_D), 0);
            return 0;
        });
        // ValueCanonicity::assign
        b.assign_region("NoteCommit input value", [&](Region& r) {
            r.copy_advice(value, advice(COL_L), 0);
            r.copy_advice(d_2.cell, advice(COL_M), 0);
            r.copy_advice(z1_d, advice(COL_R), 0);
            r.copy_advice(e_0.cell, advice(COL_Z), 0);
            r.enable_selector(selector(Q_NOTECOMMIT_VALUE), 0);
            return 0;
        });
        // RhoCanonicity::assign
        b.assign_region("NoteCommit input rho", [&](Region& r) {
            r.copy_advice(rho, advice(COL_L), 0);
            r.copy_advice(e_1.cell, advice(COL_M), 0);
            r.copy_advice(g_0_cell, advice(COL_M), 1);
            r.copy_advice(fp.cell, advice(COL_R), 0);
            r.copy_advice(e1_f_prime, advice(COL_R), 1);
            r.copy_advice(z13_f, advice(COL_Z), 0);
            r.copy_advice(z14_e1_f_prime, advice(COL_Z), 1);
            r.enable_selector(selector(Q_NOTECOMMIT_RHO), 0);
            return 0;
        });
        // PsiCanonicity::assign
        b.assign_region("NoteCommit input psi", [&](Region& r) {
            r.copy_advice(psi, advice(COL_L), 0);
            r.copy_advice(h_0.cell, advice(COL_L), 1);
            r.copy_advice(g_1.cell, advice(COL_M), 0);
            r.copy_advice(h_1_cell, advice(COL_M), 1);
            r.copy_advice(z1_g, advice(COL_R), 0);
            r.copy_advice(g1_g2_prime, advice(COL_R), 1);
            r.copy_advice(z13_g, advice(COL_Z), 0);
            r.copy_advice(z13_g1_g2_prime, advice(COL_Z), 1);
            r.enable_selector(selector(Q_NOTECOMMIT_PSI), 0);
            return 0;
        });
        return cm;
    }
};

} // namespace bb::zcash::halo2
