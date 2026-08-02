#pragma once

#include "barretenberg/commitment_schemes/whir/weights.hpp"
#include "barretenberg/common/assert.hpp"
#include "barretenberg/crypto/blake3s/blake3s.hpp"
#include "barretenberg/ecc/curves/bn254/bn254.hpp"
#include "barretenberg/ecc/curves/bn254/pairing.hpp"
#include "barretenberg/polynomials/polynomial.hpp"

#include <cstring>
#include <memory>
#include <span>
#include <string>
#include <vector>

namespace bb::kzh3 {

using Curve = curve::BN254;
using Commitment = Curve::AffineElement;
using GroupElement = Curve::Element;
using G2Affine = Curve::G2AffineElement;

/** @brief Reference to one polynomial of one committed group (duck-type shared with the suite). */
struct Kzh3ColumnRef {
    size_t group;
    size_t column;
};

/**
 * @brief KZH3 parameters: the three-dimensional tensor split (KZH-k with k=3, eprint 2025/144
 * App. C.1). Dimension 3 covers the lowest variables, dimension 1 the highest, so a coefficient
 * index decomposes as k = (i1*d2 + i2)*d3 + i3. Proof size is d1+d2 group elements plus d3 field
 * elements per opening chain; the verifier's pairing work scales with d1+d2, so the split keeps
 * d1 <= d2 <= d3.
 */
struct Kzh3Config {
    size_t num_variables;
    size_t log_d1;
    size_t log_d2;
    size_t log_d3;

    size_t d1() const { return size_t(1) << log_d1; }
    size_t d2() const { return size_t(1) << log_d2; }
    size_t d3() const { return size_t(1) << log_d3; }

    static Kzh3Config create(size_t num_variables, size_t /*security_bits*/ = 100, size_t /*log_inv_rate*/ = 0)
    {
        BB_ASSERT_GT(num_variables, size_t(2));
        return { num_variables, num_variables / 3, (num_variables + 1) / 3, (num_variables + 2) / 3 };
    }
};

/** @brief Prover-side commitment to a group: dense arrays and one KZH3 commitment per column. */
struct Kzh3GroupData {
    std::vector<std::vector<fr>> coefficients;
    std::vector<Commitment> commitments;

    size_t num_columns() const { return coefficients.size(); }
};

namespace detail {

inline std::string kzh3_label(const std::string& name, size_t i)
{
    return "KZH3:" + name + "_" + std::to_string(i);
}

inline fr derive_trapdoor(const char* tag, size_t index)
{
    std::vector<uint8_t> seed(tag, tag + std::strlen(tag));
    for (size_t byte = 0; byte < 8; ++byte) {
        seed.push_back(static_cast<uint8_t>((index >> (8 * byte)) & 0xff));
    }
    const auto digest = blake3::blake3s(seed);
    uint256_t value(0);
    std::memcpy(value.data, digest.data(), 32);
    return fr(value);
}

} // namespace detail

/**
 * @brief KZH3 structured setup: G1 generators A_{i3} for one innermost row, dimension trapdoors
 * mu1_{i1}, mu2_{i2} with V1_{i1} = mu1_{i1} V and V2_{i2} = mu2_{i2} V in G2. The commitment
 * C = sum_{i1,i2,i3} f_{i1,i2,i3} (mu1_{i1} mu2_{i2} A_{i3}) binds through the two-level pairing
 * chain e(C, V) = prod_{i1} e(D1_{i1}, V1_{i1}) with D1_{i1} = sum_{i2,i3} f (mu2_{i2} A_{i3}).
 *
 * @warning Test-only setup: this key SAMPLES the trapdoors mu (a real deployment runs a ceremony
 * and discards them) and uses them to commit via the nested MSM shortcut, which computes the
 * identical group elements a production prover derives from the published H grids at the same
 * asymptotic cost.
 */
class Kzh3CommitmentKey {
  public:
    explicit Kzh3CommitmentKey(const Kzh3Config& config)
        : config(config)
        , generators(
              Curve::Group::derive_generators(std::vector<uint8_t>{ 'b', 'b', '_', 'k', 'z', 'h', '3' }, config.d3()))
    {
        mu1.resize(config.d1());
        v1_powers.resize(config.d1());
        for (size_t i = 0; i < config.d1(); ++i) {
            mu1[i] = detail::derive_trapdoor("kzh3_mu1", i);
            v1_powers[i] = G2Affine(bb::g2::element(bb::g2::affine_one) * mu1[i]);
        }
        mu2.resize(config.d2());
        v2_powers.resize(config.d2());
        for (size_t i = 0; i < config.d2(); ++i) {
            mu2[i] = detail::derive_trapdoor("kzh3_mu2", i);
            v2_powers[i] = G2Affine(bb::g2::element(bb::g2::affine_one) * mu2[i]);
        }
    }

    /** @brief D1_{i1} = sum_{i2} mu2_{i2} * (sum_{i3} array[(i1*d2+i2)*d3+i3] A_{i3}) for one i1. */
    Commitment slice_commitment(std::span<const fr> array, size_t i1) const
    {
        const size_t d2 = config.d2();
        const size_t d3 = config.d3();
        std::vector<Commitment> inner(d2);
        for (size_t i2 = 0; i2 < d2; ++i2) {
            std::vector<fr> row(array.begin() + static_cast<std::ptrdiff_t>((i1 * d2 + i2) * d3),
                                array.begin() + static_cast<std::ptrdiff_t>((i1 * d2 + i2 + 1) * d3));
            inner[i2] = Commitment::batch_mul(std::span<const Commitment>(generators), std::span<fr>(row));
        }
        std::vector<fr> mu2_copy = mu2;
        return Commitment::batch_mul(std::span<const Commitment>(inner), std::span<fr>(mu2_copy));
    }

    Kzh3GroupData commit_group(std::vector<std::vector<fr>> payload_columns,
                               const std::vector<bool>& /*to_be_shifted*/ = {}) const
    {
        const size_t n = size_t(1) << config.num_variables;
        Kzh3GroupData data;
        for (auto& column : payload_columns) {
            BB_ASSERT_LTE(column.size(), n, "polynomial too large for the configured size");
            column.resize(n, fr::zero());
            std::vector<Commitment> slices(config.d1());
            for (size_t i1 = 0; i1 < config.d1(); ++i1) {
                slices[i1] = slice_commitment(column, i1);
            }
            std::vector<fr> mu1_copy = mu1;
            data.commitments.push_back(
                Commitment::batch_mul(std::span<const Commitment>(slices), std::span<fr>(mu1_copy)));
            data.coefficients.push_back(std::move(column));
        }
        return data;
    }

    Kzh3GroupData commit_group(std::span<const Polynomial<fr>* const> polynomials,
                               const std::vector<bool>& to_be_shifted = {}) const
    {
        std::vector<std::vector<fr>> payload_columns;
        for (const Polynomial<fr>* polynomial : polynomials) {
            std::vector<fr> column(polynomial->end_index(), fr::zero());
            for (size_t i = polynomial->start_index(); i < polynomial->end_index(); ++i) {
                column[i] = (*polynomial)[i];
            }
            payload_columns.push_back(std::move(column));
        }
        return commit_group(std::move(payload_columns), to_be_shifted);
    }

    Kzh3Config config;
    std::vector<Commitment> generators; // A_{i3}
    std::vector<fr> mu1;                // dimension-1 trapdoors (test-only)
    std::vector<fr> mu2;                // dimension-2 trapdoors (test-only)
    std::vector<G2Affine> v1_powers;    // V1_{i1} = mu1_{i1} V
    std::vector<G2Affine> v2_powers;    // V2_{i2} = mu2_{i2} V
};

/**
 * @brief KZH3 opening prover for the batched Honk claim set: chain A carries the unshifted
 * rho-combination, chain B the (unshifted!) to-be-shifted combination whose shift lives in the
 * verifier-side rank-3 query decomposition (the KZH2 backend's rank-2 trick extended one level:
 * the shift-by-one of the flattened index crosses both the d3 and d2 boundaries, so chain B sends
 * two extra final-layer vectors and one extra slice-commitment layer).
 *
 * Per chain the prover sends the KZH-k opening data: slice commitments D1 (dim 1), the contracted
 * tensor's row commitments D2 (dim 2), and the final vector T3 (dim 3). Chain B additionally sends
 * T3' (final vector under the shifted-dim-2 weights), D2'' and T3'' (the recursive opening of the
 * shifted-dim-1 contraction at row (0,0)).
 */
class Kzh3Prover {
  public:
    struct Claims {
        std::vector<const Kzh3GroupData*> groups;
        std::vector<Kzh3ColumnRef> unshifted;
        std::vector<fr> unshifted_evaluations;
        std::vector<Kzh3ColumnRef> to_be_shifted;
        std::vector<fr> shifted_evaluations;
        bool send_roots = true;
    };

    template <typename Transcript>
    static void prove(const Kzh3CommitmentKey& ck,
                      const Claims& claims,
                      std::span<const fr> u,
                      const std::shared_ptr<Transcript>& transcript)
    {
        const Kzh3Config& config = ck.config;
        const size_t n = size_t(1) << config.num_variables;
        const size_t d1 = config.d1();
        const size_t d2 = config.d2();
        const size_t d3 = config.d3();
        BB_ASSERT_EQ(u.size(), config.num_variables, "opening point size mismatch");

        if (claims.send_roots) {
            for (size_t g = 0; g < claims.groups.size(); ++g) {
                for (size_t c = 0; c < claims.groups[g]->num_columns(); ++c) {
                    transcript->send_to_verifier(detail::kzh3_label("root", g) + "_" + std::to_string(c),
                                                 claims.groups[g]->commitments[c]);
                }
            }
        }
        const fr rho = transcript->template get_challenge<fr>("KZH3:rho");

        std::vector<std::vector<fr>> chains;
        chains.emplace_back(n, fr::zero());
        fr rho_power = fr::one();
        for (size_t i = 0; i < claims.unshifted.size(); ++i) {
            const auto& column = claims.groups[claims.unshifted[i].group]->coefficients[claims.unshifted[i].column];
            for (size_t k = 0; k < n; ++k) {
                chains[0][k] += rho_power * column[k];
            }
            rho_power *= rho;
        }
        if (!claims.to_be_shifted.empty()) {
            chains.emplace_back(n, fr::zero());
            for (size_t l = 0; l < claims.to_be_shifted.size(); ++l) {
                const auto& column =
                    claims.groups[claims.to_be_shifted[l].group]->coefficients[claims.to_be_shifted[l].column];
                for (size_t k = 0; k < n; ++k) {
                    chains[1][k] += rho_power * column[k];
                }
                rho_power *= rho;
            }
        }

        const std::vector<fr> a1 = whir::eq_tensor(u.subspan(config.log_d3 + config.log_d2));
        const std::vector<fr> a2 = whir::eq_tensor(u.subspan(config.log_d3, config.log_d2));

        for (size_t chain = 0; chain < chains.size(); ++chain) {
            const std::vector<fr>& array = chains[chain];

            for (size_t i1 = 0; i1 < d1; ++i1) {
                transcript->send_to_verifier(detail::kzh3_label("D1", chain) + "_" + std::to_string(i1),
                                             ck.slice_commitment(array, i1));
            }

            // T2 = <T, a1>: contract dimension 1.
            std::vector<fr> t2(d2 * d3, fr::zero());
            for (size_t i1 = 0; i1 < d1; ++i1) {
                const fr* slice = array.data() + i1 * d2 * d3;
                for (size_t idx = 0; idx < d2 * d3; ++idx) {
                    t2[idx] += a1[i1] * slice[idx];
                }
            }
            for (size_t i2 = 0; i2 < d2; ++i2) {
                std::vector<fr> row(t2.begin() + static_cast<std::ptrdiff_t>(i2 * d3),
                                    t2.begin() + static_cast<std::ptrdiff_t>((i2 + 1) * d3));
                const Commitment d2_commitment =
                    Commitment::batch_mul(std::span<const Commitment>(ck.generators), std::span<fr>(row));
                transcript->send_to_verifier(detail::kzh3_label("D2", chain) + "_" + std::to_string(i2), d2_commitment);
            }

            // T3 = <T2, a2>: contract dimension 2.
            std::vector<fr> t3(d3, fr::zero());
            for (size_t i2 = 0; i2 < d2; ++i2) {
                for (size_t i3 = 0; i3 < d3; ++i3) {
                    t3[i3] += a2[i2] * t2[i2 * d3 + i3];
                }
            }
            for (size_t i3 = 0; i3 < d3; ++i3) {
                transcript->send_to_verifier(detail::kzh3_label("T3", chain) + "_" + std::to_string(i3), t3[i3]);
            }

            if (chain == 1) {
                // T3' = <T2, shr(a2)>: the shifted-dim-2 weights against the same D2 layer.
                std::vector<fr> t3_shift2(d3, fr::zero());
                for (size_t i2 = 1; i2 < d2; ++i2) {
                    for (size_t i3 = 0; i3 < d3; ++i3) {
                        t3_shift2[i3] += a2[i2 - 1] * t2[i2 * d3 + i3];
                    }
                }
                for (size_t i3 = 0; i3 < d3; ++i3) {
                    transcript->send_to_verifier(detail::kzh3_label("T3p", chain) + "_" + std::to_string(i3),
                                                 t3_shift2[i3]);
                }

                // T2'' = <T, shr(a1)>, opened recursively at row (0, 0).
                std::vector<fr> t2_shift1(d2 * d3, fr::zero());
                for (size_t i1 = 1; i1 < d1; ++i1) {
                    const fr* slice = array.data() + i1 * d2 * d3;
                    for (size_t idx = 0; idx < d2 * d3; ++idx) {
                        t2_shift1[idx] += a1[i1 - 1] * slice[idx];
                    }
                }
                for (size_t i2 = 0; i2 < d2; ++i2) {
                    std::vector<fr> row(t2_shift1.begin() + static_cast<std::ptrdiff_t>(i2 * d3),
                                        t2_shift1.begin() + static_cast<std::ptrdiff_t>((i2 + 1) * d3));
                    const Commitment d2s_commitment =
                        Commitment::batch_mul(std::span<const Commitment>(ck.generators), std::span<fr>(row));
                    transcript->send_to_verifier(detail::kzh3_label("D2s", chain) + "_" + std::to_string(i2),
                                                 d2s_commitment);
                }
                for (size_t i3 = 0; i3 < d3; ++i3) {
                    transcript->send_to_verifier(detail::kzh3_label("T3s", chain) + "_" + std::to_string(i3),
                                                 t2_shift1[i3]);
                }
            }
        }
    }
};

/** @brief The mirror of `Kzh3Prover`: two multipairings and the layered MSM identities per chain. */
class Kzh3Verifier {
  public:
    struct Claims {
        std::vector<size_t> group_num_columns;
        std::vector<std::vector<Commitment>> group_commitments; // when non-empty, transcript-bound
        std::vector<Kzh3ColumnRef> unshifted;
        std::vector<fr> unshifted_evaluations;
        std::vector<Kzh3ColumnRef> to_be_shifted;
        std::vector<fr> shifted_evaluations;
    };

    template <typename Transcript>
    static bool verify(const Kzh3Config& config,
                       const Claims& claims,
                       std::span<const fr> u,
                       const std::shared_ptr<Transcript>& transcript,
                       const Kzh3CommitmentKey& ck)
    {
        const size_t d1 = config.d1();
        const size_t d2 = config.d2();
        const size_t d3 = config.d3();
        BB_ASSERT_EQ(u.size(), config.num_variables, "opening point size mismatch");

        std::vector<std::vector<Commitment>> group_commitments = claims.group_commitments;
        if (group_commitments.empty()) {
            for (size_t g = 0; g < claims.group_num_columns.size(); ++g) {
                std::vector<Commitment> commitments;
                for (size_t c = 0; c < claims.group_num_columns[g]; ++c) {
                    commitments.push_back(transcript->template receive_from_prover<Commitment>(
                        detail::kzh3_label("root", g) + "_" + std::to_string(c)));
                }
                group_commitments.push_back(std::move(commitments));
            }
        }
        const fr rho = transcript->template get_challenge<fr>("KZH3:rho");

        struct ChainView {
            GroupElement commitment = GroupElement::infinity();
            fr claimed_evaluation = fr::zero();
        };
        std::vector<ChainView> chains(claims.to_be_shifted.empty() ? 1 : 2);
        fr rho_power = fr::one();
        for (size_t i = 0; i < claims.unshifted.size(); ++i) {
            chains[0].commitment +=
                GroupElement(group_commitments[claims.unshifted[i].group][claims.unshifted[i].column]) * rho_power;
            chains[0].claimed_evaluation += rho_power * claims.unshifted_evaluations[i];
            rho_power *= rho;
        }
        if (chains.size() > 1) {
            for (size_t l = 0; l < claims.to_be_shifted.size(); ++l) {
                chains[1].commitment +=
                    GroupElement(group_commitments[claims.to_be_shifted[l].group][claims.to_be_shifted[l].column]) *
                    rho_power;
                chains[1].claimed_evaluation += rho_power * claims.shifted_evaluations[l];
                rho_power *= rho;
            }
        }

        const std::vector<fr> a1 = whir::eq_tensor(u.subspan(config.log_d3 + config.log_d2));
        const std::vector<fr> a2 = whir::eq_tensor(u.subspan(config.log_d3, config.log_d2));
        const std::vector<fr> a3 = whir::eq_tensor(u.subspan(0, config.log_d3));
        std::vector<fr> a1_shift(d1, fr::zero());
        for (size_t i1 = 1; i1 < d1; ++i1) {
            a1_shift[i1] = a1[i1 - 1];
        }
        std::vector<fr> a2_shift(d2, fr::zero());
        for (size_t i2 = 1; i2 < d2; ++i2) {
            a2_shift[i2] = a2[i2 - 1];
        }

        for (size_t chain = 0; chain < chains.size(); ++chain) {
            std::vector<Commitment> d1_commitments(d1);
            for (size_t i1 = 0; i1 < d1; ++i1) {
                d1_commitments[i1] = transcript->template receive_from_prover<Commitment>(
                    detail::kzh3_label("D1", chain) + "_" + std::to_string(i1));
            }
            std::vector<Commitment> d2_commitments(d2);
            for (size_t i2 = 0; i2 < d2; ++i2) {
                d2_commitments[i2] = transcript->template receive_from_prover<Commitment>(
                    detail::kzh3_label("D2", chain) + "_" + std::to_string(i2));
            }
            std::vector<fr> t3(d3);
            for (size_t i3 = 0; i3 < d3; ++i3) {
                t3[i3] = transcript->template receive_from_prover<fr>(detail::kzh3_label("T3", chain) + "_" +
                                                                      std::to_string(i3));
            }

            // Check 1 (binding, dim 1): e(-C, V) * prod_{i1} e(D1_{i1}, V1_{i1}) == 1.
            if (!pairing_chain_check(Commitment(-chains[chain].commitment), d1_commitments, ck.v1_powers)) {
                return false;
            }

            // Check 2 (binding, dim 2): C1 = <a1, D1>; e(-C1, V) * prod_{i2} e(D2_{i2}, V2_{i2}) == 1.
            std::vector<fr> a1_copy = a1;
            const Commitment c1 =
                Commitment::batch_mul(std::span<const Commitment>(d1_commitments), std::span<fr>(a1_copy));
            if (!pairing_chain_check(Commitment(-GroupElement(c1)), d2_commitments, ck.v2_powers)) {
                return false;
            }

            // Check 3 (final layer): <T3, A> == <a2, D2>.
            if (!final_layer_check(ck, t3, d2_commitments, a2)) {
                return false;
            }

            if (chain == 0) {
                // Check 4 (claim): <T3, a3> == claimed evaluation.
                fr total = fr::zero();
                for (size_t i3 = 0; i3 < d3; ++i3) {
                    total += t3[i3] * a3[i3];
                }
                if (total != chains[chain].claimed_evaluation) {
                    return false;
                }
            } else {
                std::vector<fr> t3_shift2(d3);
                for (size_t i3 = 0; i3 < d3; ++i3) {
                    t3_shift2[i3] = transcript->template receive_from_prover<fr>(detail::kzh3_label("T3p", chain) +
                                                                                 "_" + std::to_string(i3));
                }
                std::vector<Commitment> d2s_commitments(d2);
                for (size_t i2 = 0; i2 < d2; ++i2) {
                    d2s_commitments[i2] = transcript->template receive_from_prover<Commitment>(
                        detail::kzh3_label("D2s", chain) + "_" + std::to_string(i2));
                }
                std::vector<fr> t3_shift1(d3);
                for (size_t i3 = 0; i3 < d3; ++i3) {
                    t3_shift1[i3] = transcript->template receive_from_prover<fr>(detail::kzh3_label("T3s", chain) +
                                                                                 "_" + std::to_string(i3));
                }

                // Chain B check 3': <T3', A> == <shr(a2), D2> (same D2 layer, shifted weights).
                if (!final_layer_check(ck, t3_shift2, d2_commitments, a2_shift)) {
                    return false;
                }

                // Chain B checks for T'' = <T, shr(a1)>: pairing-bind D2'' to C1'' = <shr(a1), D1>,
                // then bind the row-(0) vector T3'' to D2''_0.
                std::vector<fr> a1_shift_copy = a1_shift;
                const Commitment c1s =
                    Commitment::batch_mul(std::span<const Commitment>(d1_commitments), std::span<fr>(a1_shift_copy));
                if (!pairing_chain_check(Commitment(-GroupElement(c1s)), d2s_commitments, ck.v2_powers)) {
                    return false;
                }
                std::vector<fr> t3_shift1_copy = t3_shift1;
                const Commitment t3s_commitment =
                    Commitment::batch_mul(std::span<const Commitment>(ck.generators), std::span<fr>(t3_shift1_copy));
                if (t3s_commitment != d2s_commitments[0]) {
                    return false;
                }

                // Chain B claim: the rank-3 shift decomposition. The shift-by-one of the flattened
                // index contributes (i) T3 under the shifted-dim-3 weights, (ii) the dim-2 boundary
                // term via T3'[0], and (iii) the dim-1 boundary term via T3''[0].
                fr total = fr::zero();
                for (size_t i3 = 1; i3 < d3; ++i3) {
                    total += t3[i3] * a3[i3 - 1];
                }
                total += a3[d3 - 1] * t3_shift2[0];
                total += a2[d2 - 1] * a3[d3 - 1] * t3_shift1[0];
                if (total != chains[chain].claimed_evaluation) {
                    return false;
                }
            }
        }
        return true;
    }

  private:
    static bool pairing_chain_check(const Commitment& negated_parent,
                                    const std::vector<Commitment>& children,
                                    const std::vector<G2Affine>& v_powers)
    {
        std::vector<Commitment> pairing_g1;
        std::vector<G2Affine> pairing_g2;
        pairing_g1.push_back(negated_parent);
        pairing_g2.push_back(G2Affine(bb::g2::affine_one));
        for (size_t i = 0; i < children.size(); ++i) {
            pairing_g1.push_back(children[i]);
            pairing_g2.push_back(v_powers[i]);
        }
        return pairing::reduced_ate_pairing_batch(pairing_g1.data(), pairing_g2.data(), pairing_g1.size()) ==
               fq12::one();
    }

    static bool final_layer_check(const Kzh3CommitmentKey& ck,
                                  const std::vector<fr>& t3,
                                  const std::vector<Commitment>& d2_commitments,
                                  const std::vector<fr>& weights)
    {
        std::vector<fr> weights_copy = weights;
        const Commitment expected =
            Commitment::batch_mul(std::span<const Commitment>(d2_commitments), std::span<fr>(weights_copy));
        std::vector<fr> t3_copy = t3;
        const Commitment computed =
            Commitment::batch_mul(std::span<const Commitment>(ck.generators), std::span<fr>(t3_copy));
        return computed == expected;
    }
};

} // namespace bb::kzh3
