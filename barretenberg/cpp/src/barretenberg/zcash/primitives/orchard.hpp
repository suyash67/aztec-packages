#pragma once

#include "barretenberg/common/assert.hpp"
#include "cycle.hpp"
#include "fixed_base.hpp"
#include "orchard_fixed_base_z.hpp"
#include "orchard_fixed_base_z_bn254.hpp"
#include "poseidon.hpp"
#include "sinsemilla.hpp"

#include <array>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace bb::zcash {

/**
 * @brief Orchard Action constants and native (out-of-circuit) evaluation of the Action relation.
 * @details This is the reference against which both circuit implementations are checked. Definitions follow the Zcash
 * protocol specification §4.17 (Action statement), §5.4.8 (commitments) and §5.4.2 (PRF^nf); the comments use the
 * same names as the attached circuit description.
 */
template <typename Cycle> struct Orchard {
    using FF = typename Cycle::FF;
    using AffineElement = typename Cycle::AffineElement;
    using Element = typename Cycle::Element;
    using Scalar = typename Cycle::EmbeddedScalar;
    using Sins = Sinsemilla<Cycle>;
    using Poseidon = PoseidonP128Pow5T3<FF>;

    static constexpr size_t MERKLE_DEPTH = 32;
    // I2LEBSP length used for every field element in Sinsemilla messages (ℓ_P = 255 for Pallas).
    static constexpr size_t L_BASE = 255;
    static constexpr size_t L_VALUE = 64;
    static constexpr size_t NUM_WINDOWS = 85;
    static constexpr size_t NUM_WINDOWS_SHORT = 22;

    static constexpr const char* NOTE_COMMIT_DOMAIN = "z.cash:Orchard-NoteCommit";
    static constexpr const char* COMMIT_IVK_DOMAIN = "z.cash:Orchard-CommitIvk";
    static constexpr const char* MERKLE_CRH_DOMAIN = "z.cash:Orchard-MerkleCRH";

    /**
     * @brief Generators and fixed-base tables. For Pasta these are Orchard's constants (with Orchard's published
     * z-values); for other cycles the same domain strings are hashed with the cycle's GroupHash and z-values are
     * searched for.
     */
    struct Constants {
        AffineElement spend_auth_g;   // G
        AffineElement nullifier_k;    // K
        AffineElement value_commit_v; // V
        AffineElement value_commit_r; // R
        AffineElement note_commit_r;  // R_nc
        AffineElement commit_ivk_r;   // R_ivk
        AffineElement q_merkle_crh;
        AffineElement q_note_commit;
        AffineElement q_commit_ivk;
    };

    /**
     * @brief Window tables of the fixed bases for halo2's fixed-base multiplication (built on first use: for cycles
     * other than Pasta the z-values are searched for, which is slow).
     */
    struct FixedBases {
        FixedBase<Cycle> fb_spend_auth_g;
        FixedBase<Cycle> fb_value_commit_r;
        FixedBase<Cycle> fb_note_commit_r;
        FixedBase<Cycle> fb_commit_ivk_r;
        FixedBase<Cycle> fb_nullifier_k;    // base-field-element scalar
        FixedBase<Cycle> fb_value_commit_v; // short signed scalar
    };

    static const Constants& constants()
    {
        static const Constants c = []() {
            Constants c;
            c.spend_auth_g = group_hash<Cycle>("z.cash:Orchard", "G");
            c.nullifier_k = group_hash<Cycle>("z.cash:Orchard", "K");
            c.value_commit_v = group_hash<Cycle>("z.cash:Orchard-cv", "v");
            c.value_commit_r = group_hash<Cycle>("z.cash:Orchard-cv", "r");
            c.note_commit_r = Sins::R(NOTE_COMMIT_DOMAIN);
            c.commit_ivk_r = Sins::R(COMMIT_IVK_DOMAIN);
            c.q_merkle_crh = Sins::Q(MERKLE_CRH_DOMAIN);
            c.q_note_commit = Sins::Q(std::string(NOTE_COMMIT_DOMAIN) + "-M");
            c.q_commit_ivk = Sins::Q(std::string(COMMIT_IVK_DOMAIN) + "-M");
            return c;
        }();
        return c;
    }

    static const FixedBases& fixed_bases()
    {
        static const FixedBases f = []() {
            const auto& k = constants();
            FixedBases f;
            // Pasta: Orchard's published z-values; BN254: the generated table (fixed_base_z_gen.test.cpp); other cycles
            // search for them.
            auto zs = [&](const AffineElement& g,
                          size_t windows,
                          std::span<const uint64_t> pasta,
                          std::span<const uint64_t> bn254) {
                if constexpr (std::is_same_v<Cycle, PastaCycle>) {
                    return std::vector<uint64_t>(pasta.begin(), pasta.end());
                } else if constexpr (std::is_same_v<Cycle, Bn254Cycle>) {
                    return std::vector<uint64_t>(bn254.begin(), bn254.end());
                } else {
                    return FixedBase<Cycle>::find_zs(g, windows);
                }
            };
            namespace pz = orchard_constants;
            namespace bz = orchard_constants_bn254;
            f.fb_spend_auth_g = FixedBase<Cycle>::make(
                k.spend_auth_g, NUM_WINDOWS, zs(k.spend_auth_g, NUM_WINDOWS, pz::SPEND_AUTH_G_Z, bz::SPEND_AUTH_G_Z));
            f.fb_value_commit_r =
                FixedBase<Cycle>::make(k.value_commit_r,
                                       NUM_WINDOWS,
                                       zs(k.value_commit_r, NUM_WINDOWS, pz::VALUE_COMMIT_R_Z, bz::VALUE_COMMIT_R_Z));
            f.fb_note_commit_r =
                FixedBase<Cycle>::make(k.note_commit_r,
                                       NUM_WINDOWS,
                                       zs(k.note_commit_r, NUM_WINDOWS, pz::NOTE_COMMIT_R_Z, bz::NOTE_COMMIT_R_Z));
            f.fb_commit_ivk_r = FixedBase<Cycle>::make(
                k.commit_ivk_r, NUM_WINDOWS, zs(k.commit_ivk_r, NUM_WINDOWS, pz::COMMIT_IVK_R_Z, bz::COMMIT_IVK_R_Z));
            f.fb_nullifier_k = FixedBase<Cycle>::make(
                k.nullifier_k, NUM_WINDOWS, zs(k.nullifier_k, NUM_WINDOWS, pz::NULLIFIER_K_Z, bz::NULLIFIER_K_Z));
            f.fb_value_commit_v = FixedBase<Cycle>::make(
                k.value_commit_v,
                NUM_WINDOWS_SHORT,
                zs(k.value_commit_v, NUM_WINDOWS_SHORT, pz::VALUE_COMMIT_V_Z, bz::VALUE_COMMIT_V_Z));
            return f;
        }();
        return f;
    }

    // P* = I2LEBSP_255(x(P)) || ỹ(P)
    static void append_point_repr(std::vector<bool>& bits, const AffineElement& p)
    {
        append_bits(bits, p.x, L_BASE);
        bits.push_back(uint256_t(p.y).get_bit(0));
    }

    static std::optional<AffineElement> note_commit(const AffineElement& g_d,
                                                    const AffineElement& pk_d,
                                                    uint64_t v,
                                                    const FF& rho,
                                                    const FF& psi,
                                                    const Scalar& rcm)
    {
        std::vector<bool> bits;
        append_point_repr(bits, g_d);
        append_point_repr(bits, pk_d);
        append_bits(bits, uint256_t(v), L_VALUE);
        append_bits(bits, rho, L_BASE);
        append_bits(bits, psi, L_BASE);
        return Sins::commit(NOTE_COMMIT_DOMAIN, bits, rcm);
    }

    static std::optional<FF> commit_ivk(const FF& ak, const FF& nk, const Scalar& rivk)
    {
        std::vector<bool> bits;
        append_bits(bits, ak, L_BASE);
        append_bits(bits, nk, L_BASE);
        return Sins::short_commit(COMMIT_IVK_DOMAIN, bits, rivk);
    }

    // MerkleCRH(l, left, right) = Extract_P(SinsemillaHashToPoint("z.cash:Orchard-MerkleCRH", l_[10] || left || right))
    static FF merkle_crh(size_t layer, const FF& left, const FF& right)
    {
        std::vector<bool> bits;
        append_bits(bits, uint256_t(layer), 10);
        append_bits(bits, left, L_BASE);
        append_bits(bits, right, L_BASE);
        auto p = Sins::hash_to_point(constants().q_merkle_crh, bits);
        BB_ASSERT(p.has_value());
        return p->x;
    }

    static FF merkle_root(const FF& leaf, const std::array<FF, MERKLE_DEPTH>& path, uint32_t pos)
    {
        FF node = leaf;
        for (size_t l = 0; l < MERKLE_DEPTH; ++l) {
            const bool right_child = ((pos >> l) & 1) != 0;
            node = right_child ? merkle_crh(l, path[l], node) : merkle_crh(l, node, path[l]);
        }
        return node;
    }

    // Embedded-curve scalar with the same integer representative as a circuit-field element (ToScalar).
    static Scalar to_scalar(const FF& x) { return Scalar(uint256_t(x)); }

    // nf = Extract_P([PRF^nf_nk(rho) + psi] K + cm)
    static FF derive_nullifier(const FF& nk, const FF& rho, const FF& psi, const AffineElement& cm)
    {
        const FF scalar = Poseidon::hash(nk, rho) + psi;
        return AffineElement(Element(constants().nullifier_k) * to_scalar(scalar) + Element(cm)).x;
    }

    struct ActionWitness {
        std::array<FF, MERKLE_DEPTH> path{};
        uint32_t pos = 0;
        AffineElement g_d_old;
        AffineElement pk_d_old;
        uint64_t v_old = 0;
        FF rho_old;
        FF psi_old;
        Scalar rcm_old;
        AffineElement cm_old;
        Scalar alpha;
        AffineElement ak;
        FF nk;
        Scalar rivk;
        AffineElement g_d_new;
        AffineElement pk_d_new;
        uint64_t v_new = 0;
        FF psi_new;
        Scalar rcm_new;
        Scalar rcv;
    };

    // Rows 0..9 of the instance column.
    struct PublicInputs {
        FF anchor;
        AffineElement cv_net;
        FF nf_old;
        AffineElement rk;
        FF cmx;
        bool enable_spend = true;
        bool enable_output = true;
        bool disable_cross_address = false;

        std::array<FF, 10> to_field_elements() const
        {
            return { anchor,
                     cv_net.x,
                     cv_net.y,
                     nf_old,
                     rk.x,
                     rk.y,
                     cmx,
                     FF(enable_spend ? 1 : 0),
                     FF(enable_output ? 1 : 0),
                     FF(disable_cross_address ? 1 : 0) };
        }
    };

    // v_old - v_new = magnitude * sign
    static std::pair<uint64_t, bool> magnitude_sign(uint64_t v_old, uint64_t v_new)
    {
        return v_old >= v_new ? std::make_pair(v_old - v_new, false) : std::make_pair(v_new - v_old, true);
    }

    static AffineElement value_commit(uint64_t v_old, uint64_t v_new, const Scalar& rcv)
    {
        const auto [magnitude, negative] = magnitude_sign(v_old, v_new);
        Scalar v_net = Scalar(magnitude);
        if (negative) {
            v_net = -v_net;
        }
        const auto& c = constants();
        return AffineElement(Element(c.value_commit_v) * v_net + Element(c.value_commit_r) * rcv);
    }

    /**
     * @brief Evaluate the Action statement natively.
     * @return the public inputs implied by the witness, or nullopt if the witness violates the statement (the parts
     * that are not determined by the public inputs: old-note commitment integrity and diversified-address integrity).
     */
    static std::optional<PublicInputs> evaluate(const ActionWitness& w,
                                                bool enable_spend = true,
                                                bool enable_output = true)
    {
        const auto& c = constants();
        // 1. Old note commitment integrity
        auto cm_old = note_commit(w.g_d_old, w.pk_d_old, w.v_old, w.rho_old, w.psi_old, w.rcm_old);
        if (!cm_old || *cm_old != w.cm_old) {
            return std::nullopt;
        }
        // 6. Diversified address integrity: pk_d_old = [ivk] g_d_old, ivk = CommitIvk(Extract(ak), nk)
        auto ivk = commit_ivk(w.ak.x, w.nk, w.rivk);
        if (!ivk || AffineElement(Element(w.g_d_old) * to_scalar(*ivk)) != w.pk_d_old) {
            return std::nullopt;
        }
        // 8. Enable flags
        if ((w.v_old != 0 && !enable_spend) || (w.v_new != 0 && !enable_output)) {
            return std::nullopt;
        }
        PublicInputs pi;
        // 2. Merkle path validity (the anchor is only bound when v_old != 0)
        pi.anchor = merkle_root(w.cm_old.x, w.path, w.pos);
        // 3. Value commitment integrity
        pi.cv_net = value_commit(w.v_old, w.v_new, w.rcv);
        // 4. Nullifier integrity
        pi.nf_old = derive_nullifier(w.nk, w.rho_old, w.psi_old, w.cm_old);
        // 5. Spend authority
        pi.rk = AffineElement(Element(c.spend_auth_g) * w.alpha + Element(w.ak));
        // 7. New note commitment integrity, rho_new = nf_old
        auto cm_new = note_commit(w.g_d_new, w.pk_d_new, w.v_new, pi.nf_old, w.psi_new, w.rcm_new);
        if (!cm_new) {
            return std::nullopt;
        }
        pi.cmx = cm_new->x;
        pi.enable_spend = enable_spend;
        pi.enable_output = enable_output;
        return pi;
    }

    /**
     * @brief A random valid Action: a real spend of a fresh note at a random tree position, and a real output.
     * @details Keys follow Orchard's key schedule shape (ak, nk, rivk -> ivk -> pk_d = [ivk] g_d), with g_d a random
     * point rather than DiversifyHash(d), which the circuit does not check.
     */
    template <typename Rng> static ActionWitness random_witness(Rng& engine)
    {
        ActionWitness w;
        auto random_point = [&]() { return AffineElement(Element::random_element(&engine)); };
        w.ak = random_point();
        w.nk = FF::random_element(&engine);
        w.rivk = Scalar::random_element(&engine);
        auto ivk = commit_ivk(w.ak.x, w.nk, w.rivk);
        BB_ASSERT(ivk.has_value());
        w.g_d_old = random_point();
        w.pk_d_old = AffineElement(Element(w.g_d_old) * to_scalar(*ivk));
        w.v_old = engine.get_random_uint64() >> 8;
        w.rho_old = FF::random_element(&engine);
        w.psi_old = FF::random_element(&engine);
        w.rcm_old = Scalar::random_element(&engine);
        auto cm = note_commit(w.g_d_old, w.pk_d_old, w.v_old, w.rho_old, w.psi_old, w.rcm_old);
        BB_ASSERT(cm.has_value());
        w.cm_old = *cm;
        for (auto& node : w.path) {
            node = FF::random_element(&engine);
        }
        w.pos = engine.get_random_uint32();
        w.alpha = Scalar::random_element(&engine);
        w.g_d_new = random_point();
        w.pk_d_new = random_point();
        w.v_new = engine.get_random_uint64() >> 8;
        w.psi_new = FF::random_element(&engine);
        w.rcm_new = Scalar::random_element(&engine);
        w.rcv = Scalar::random_element(&engine);
        return w;
    }
};

} // namespace bb::zcash
