#pragma once

#include "barretenberg/commitment_schemes/commitment_key.hpp"
#include "barretenberg/commitment_schemes/mercury/mercury.hpp"
#include "barretenberg/commitment_schemes/mercury/vela.hpp"
#include "barretenberg/commitment_schemes/pairing_points.hpp"
#include "barretenberg/commitment_schemes/utils/batch_accumulate.hpp"
#include "barretenberg/common/assert.hpp"
#include "barretenberg/ecc/curves/bn254/bn254.hpp"
#include "barretenberg/polynomials/polynomial.hpp"
#include "barretenberg/polynomials/polynomial_arithmetic.hpp"

#include <memory>
#include <span>
#include <string>
#include <vector>

namespace bb::fflonk {

using Curve = curve::BN254;
using Commitment = Curve::AffineElement;
using GroupElement = Curve::Element;

/**
 * @brief fflonk-style packed commitments (ePrint 2021/1167) for a batch of multilinear claims.
 *
 * @details The verifier's dominant cost in every KZG-based multilinear scheme here is the
 * multi-scalar multiplication that batches the ~36 committed columns into one commitment, because the
 * batching challenge is drawn after the commitments are fixed. fflonk removes it. A group of `t`
 * columns is committed as the single interleaved polynomial
 *
 *      g(X) = sum_{i<t} f_i(X^t) * X^i,   deg g < t*n,
 *
 * and the identity `g mod (X^t - z) = sum_i f_i(z) X^i` means one opening of `g` on the `t`-th roots
 * of `z` certifies *every* column's evaluation at `z` at once - and the residue's coefficients are
 * exactly those evaluations, so the verifier never has to touch a root of unity or an inverse DFT.
 * Batching then happens in the field, over the certified evaluations, instead of in the group.
 *
 * Opening at `z` and `1/z` together lets the same machinery carry Vela's Laurent reduction
 * (ePrint 2026/1438) from multilinear claims to univariate ones: with
 * `Z_r(X) = (X^{t_r} - z)(X^{t_r} - 1/z)` the residue splits as `A + B X^{t_r}` with
 * `B = (P - Q)/(z - 1/z)` and `A = P - zB`, where `P` and `Q` are the two evaluation vectors. Both
 * are closed forms in `t_r` field operations.
 *
 * A BDFG21 (Shplonk) batch over the groups' point sets closes with one quotient, one linearization
 * and one pairing check, so the verifier's whole group-element budget is
 * `#groups + [h] + [W] + [W'] + [1]` - eight points here, against thirty-five for the same claim set
 * under Vela.
 *
 * The price is on the prover: the committed degree of a group is `t` times the payload, so the SRS
 * and the commitment MSM grow with the packing factor. That is the tradeoff this backend exposes -
 * verifier group operations fall as `1/t`, prover degree rises as `t`.
 */
struct FflonkConfig {
    size_t num_variables;

    static FflonkConfig create(size_t num_variables, size_t /*security_bits*/ = 100, size_t /*log_inv_rate*/ = 0)
    {
        BB_ASSERT_GT(num_variables, size_t(1));
        return { num_variables };
    }
};

/** @brief One commitment round: the columns, their interleaved packing, and its commitment. */
struct FflonkGroupData {
    std::vector<std::vector<fr>> coefficients; // the payload columns, dense
    std::vector<fr> packed;                    // g(X) = sum_i f_i(X^t) X^i
    Commitment commitment;
    size_t pack = 1; // t, the next power of two at or above the column count

    size_t num_columns() const { return coefficients.size(); }
};

namespace detail {

inline std::string fflonk_label(const std::string& name)
{
    return "FFLONK:" + name;
}
inline std::string fflonk_label(const std::string& name, size_t i, size_t j)
{
    return "FFLONK:" + name + "_" + std::to_string(i) + "_" + std::to_string(j);
}

inline size_t next_power_of_two(size_t value)
{
    size_t power = 1;
    while (power < value) {
        power <<= 1;
    }
    return power;
}

/**
 * @brief Quotient and residue of `poly` by `X^t - c`.
 * @details From `f[k] = q[k-t] - c q[k]`, the quotient falls out top-down as
 * `q[j] = f[j+t] + c q[j+t]`, and the residue is `r[i] = f[i] + c q[i]`.
 */
inline void divide_by_power_minus(
    const std::vector<fr>& poly, size_t t, const fr& c, std::vector<fr>& quotient, std::vector<fr>& residue)
{
    const size_t size = poly.size();
    quotient.assign(size, fr::zero());
    for (size_t j = size; j-- > 0;) {
        if (j + t < size) {
            quotient[j] = poly[j + t] + c * quotient[j + t];
        }
    }
    residue.assign(t, fr::zero());
    for (size_t i = 0; i < t && i < size; ++i) {
        residue[i] = poly[i] + c * quotient[i];
    }
}

/** @brief `poly` at `x`, Horner. */
inline fr evaluate(std::span<const fr> poly, const fr& x)
{
    fr accumulator = fr::zero();
    for (size_t i = poly.size(); i-- > 0;) {
        accumulator = accumulator * x + poly[i];
    }
    return accumulator;
}

/**
 * @brief `R_r(y)`, the interpolant of a packed group on `{t-th roots of z} u {t-th roots of 1/z}`.
 * @details `P` and `Q` are the group's column evaluations at `z` and `1/z`, which are exactly the
 * residues of `g` modulo `X^t - z` and `X^t - 1/z`. The interpolant is `A + B X^t` with
 * `B = (P - Q)/(z - 1/z)` and `A = P - zB`.
 */
inline fr interpolant_at(
    std::span<const fr> p_evals, std::span<const fr> q_evals, const fr& z, const fr& z_inv, const fr& y, size_t t)
{
    const fr denominator_inv = (z - z_inv).invert();
    std::vector<fr> a(t);
    std::vector<fr> b(t);
    for (size_t i = 0; i < t; ++i) {
        b[i] = (p_evals[i] - q_evals[i]) * denominator_inv;
        a[i] = p_evals[i] - z * b[i];
    }
    return evaluate(a, y) + y.pow(uint256_t(uint64_t(t))) * evaluate(b, y);
}

} // namespace detail

/** @brief Commitment key: packs a round's columns into one polynomial and commits it. */
class FflonkCommitmentKey {
  public:
    explicit FflonkCommitmentKey(const FflonkConfig& config)
        : config(config)
        // The widest packing this claim shape uses is the precomputed round; size the SRS for it.
        , kzg_ck(std::make_shared<CommitmentKey<Curve>>(size_t(1) << (config.num_variables + 5)))
    {}

    FflonkGroupData commit_group(std::vector<std::vector<fr>> payload_columns,
                                 const std::vector<bool>& /*to_be_shifted*/ = {}) const
    {
        const size_t n = size_t(1) << config.num_variables;
        FflonkGroupData data;
        data.pack = detail::next_power_of_two(payload_columns.size());
        for (auto& column : payload_columns) {
            BB_ASSERT_LTE(column.size(), n, "polynomial too large for the configured size");
            column.resize(n, fr::zero());
            data.coefficients.push_back(std::move(column));
        }
        // g[j*t + i] = f_i[j]: the interleaving that makes `g mod (X^t - z)` the evaluation vector.
        data.packed.assign(n * data.pack, fr::zero());
        for (size_t i = 0; i < data.coefficients.size(); ++i) {
            for (size_t j = 0; j < n; ++j) {
                data.packed[j * data.pack + i] = data.coefficients[i][j];
            }
        }
        data.commitment = kzg_ck->commit(Polynomial<fr>(std::span<const fr>(data.packed)));
        return data;
    }

    FflonkGroupData commit_group(std::span<const Polynomial<fr>* const> polynomials,
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

    FflonkConfig config;
    std::shared_ptr<CommitmentKey<Curve>> kzg_ck;
};

/** @brief Reference to one column of one committed round. */
struct FflonkColumnRef {
    size_t group;
    size_t column;
};

/**
 * @brief fflonk opening prover: certify every column's evaluations at {z, 1/z}, then close the
 * multilinear claim in the field.
 */
class FflonkProver {
  public:
    struct Claims {
        std::vector<const FflonkGroupData*> groups;
        std::vector<FflonkColumnRef> unshifted;
        std::vector<fr> unshifted_evaluations;
        std::vector<FflonkColumnRef> to_be_shifted;
        std::vector<fr> shifted_evaluations;
        bool send_roots = true;
    };

    template <typename Transcript>
    static void prove(const FflonkCommitmentKey& ck,
                      const Claims& claims,
                      std::span<const fr> r,
                      const std::shared_ptr<Transcript>& transcript)
    {
        const size_t n = size_t(1) << r.size();
        const bool has_shifted = !claims.to_be_shifted.empty();

        if (claims.send_roots) {
            for (size_t g = 0; g < claims.groups.size(); ++g) {
                transcript->send_to_verifier(detail::fflonk_label("root", g, 0), claims.groups[g]->commitment);
            }
        }
        const fr rho = transcript->template get_challenge<fr>("FFLONK:rho");
        const fr lambda = transcript->template get_challenge<fr>("FFLONK:lambda");

        // The Laurent residual h of the rho-combined twin, exactly as Vela builds it. Only h needs
        // the combination as a *polynomial*; every other use of it is over certified evaluations.
        std::vector<fr> f_a(n, fr::zero());
        std::vector<fr> f_b(n, fr::zero());
        fr y_combined = fr::zero();
        fr rho_power = fr::one();
        std::vector<pcs_utils::ScaledTerm> terms_a;
        std::vector<pcs_utils::ScaledTerm> terms_b;
        for (size_t i = 0; i < claims.unshifted.size(); ++i) {
            const auto& column = claims.groups[claims.unshifted[i].group]->coefficients[claims.unshifted[i].column];
            terms_a.push_back({ column.data(), rho_power, false });
            y_combined += rho_power * claims.unshifted_evaluations[i];
            rho_power *= rho;
        }
        for (size_t l = 0; l < claims.to_be_shifted.size(); ++l) {
            const auto& column =
                claims.groups[claims.to_be_shifted[l].group]->coefficients[claims.to_be_shifted[l].column];
            terms_b.push_back({ column.data(), rho_power, false });
            y_combined += lambda * rho_power * claims.shifted_evaluations[l];
            rho_power *= rho;
        }
        pcs_utils::accumulate_scaled(f_a, terms_a);
        pcs_utils::accumulate_scaled(f_b, terms_b);

        std::vector<fr> g(n + 1, fr::zero());
        for (size_t k = 0; k < n; ++k) {
            g[k + 1] += f_a[k];
        }
        if (has_shifted) {
            for (size_t k = 0; k < n; ++k) {
                g[k] += lambda * f_b[k];
            }
        }
        const std::vector<fr> table = vela::detail::laurent_product(g, /*shift_low=*/1, r);
        const size_t offset = n;
        std::vector<fr> h(n, fr::zero());
        for (size_t k = 1; k <= n; ++k) {
            fr d_k = fr::zero();
            if (offset + k < table.size()) {
                d_k += table[offset + k];
            }
            if (offset >= k) {
                d_k += table[offset - k];
            }
            h[k - 1] = d_k;
        }
        transcript->send_to_verifier(detail::fflonk_label("C_h"), ck.kzg_ck->commit(Polynomial<fr>(std::span<const fr>(h))));

        const fr z = transcript->template get_challenge<fr>("FFLONK:z");
        BB_ASSERT(z != fr::zero() && z != fr::one() && z != -fr::one(), "degenerate fflonk challenge");
        const fr z_inv = z.invert();

        // Every column at z and 1/z: these are the coefficients of the packed polynomial's residues
        // modulo (X^t - z) and (X^t - 1/z), so sending them *is* sending the openings.
        for (size_t gi = 0; gi < claims.groups.size(); ++gi) {
            const FflonkGroupData& group = *claims.groups[gi];
            for (size_t c = 0; c < group.pack; ++c) {
                const fr at_z = c < group.num_columns()
                                    ? polynomial_arithmetic::evaluate(group.coefficients[c].data(), z, n)
                                    : fr::zero();
                const fr at_z_inv = c < group.num_columns()
                                        ? polynomial_arithmetic::evaluate(group.coefficients[c].data(), z_inv, n)
                                        : fr::zero();
                transcript->send_to_verifier(detail::fflonk_label("p", gi, c), at_z);
                transcript->send_to_verifier(detail::fflonk_label("q", gi, c), at_z_inv);
            }
        }
        transcript->send_to_verifier(detail::fflonk_label("h_z"), polynomial_arithmetic::evaluate(h.data(), z, n));
        transcript->send_to_verifier(detail::fflonk_label("h_z_inv"),
                                     polynomial_arithmetic::evaluate(h.data(), z_inv, n));

        // Shplonk (BDFG21) over the groups' point sets plus h's: one quotient, one linearization.
        std::vector<const std::vector<fr>*> polynomials;
        std::vector<size_t> packs;
        for (const FflonkGroupData* group : claims.groups) {
            polynomials.push_back(&group->packed);
            packs.push_back(group->pack);
        }
        polynomials.push_back(&h);
        packs.push_back(1);

        const fr nu = transcript->template get_challenge<fr>("FFLONK:nu");
        size_t widest = 0;
        for (const std::vector<fr>* poly : polynomials) {
            widest = std::max(widest, poly->size());
        }
        std::vector<fr> w_poly(widest, fr::zero());
        std::vector<std::vector<fr>> residues(polynomials.size());
        fr nu_power = fr::one();
        for (size_t i = 0; i < polynomials.size(); ++i) {
            std::vector<fr> quotient_one;
            std::vector<fr> residue_one;
            detail::divide_by_power_minus(*polynomials[i], packs[i], z, quotient_one, residue_one);
            // (P - R) / ((X^t - z)(X^t - 1/z)): divide by the first factor, then the second. The
            // residue modulo the product is recovered by the verifier's closed form, so only the
            // quotient is needed here.
            std::vector<fr> quotient_two;
            std::vector<fr> residue_two;
            detail::divide_by_power_minus(quotient_one, packs[i], z_inv, quotient_two, residue_two);
            for (size_t k = 0; k < quotient_two.size(); ++k) {
                w_poly[k] += nu_power * quotient_two[k];
            }
            nu_power *= nu;
        }
        transcript->send_to_verifier(detail::fflonk_label("C_w"),
                                     ck.kzg_ck->commit(Polynomial<fr>(std::span<const fr>(w_poly))));

        const fr y = transcript->template get_challenge<fr>("FFLONK:y");
        fr z_total = fr::one();
        std::vector<fr> z_at_y(polynomials.size());
        for (size_t i = 0; i < polynomials.size(); ++i) {
            const fr y_pow = y.pow(uint256_t(uint64_t(packs[i])));
            z_at_y[i] = (y_pow - z) * (y_pow - z_inv);
            z_total *= z_at_y[i];
        }
        std::vector<fr> l_poly(widest, fr::zero());
        nu_power = fr::one();
        fr constant = fr::zero();
        for (size_t i = 0; i < polynomials.size(); ++i) {
            const fr scalar = nu_power * z_total * z_at_y[i].invert();
            for (size_t k = 0; k < polynomials[i]->size(); ++k) {
                l_poly[k] += scalar * (*polynomials[i])[k];
            }
            std::vector<fr> p_evals(packs[i]);
            std::vector<fr> q_evals(packs[i]);
            {
                std::vector<fr> quotient;
                detail::divide_by_power_minus(*polynomials[i], packs[i], z, quotient, p_evals);
                detail::divide_by_power_minus(*polynomials[i], packs[i], z_inv, quotient, q_evals);
            }
            constant += scalar * detail::interpolant_at(p_evals, q_evals, z, z_inv, y, packs[i]);
            nu_power *= nu;
        }
        l_poly[0] -= constant;
        for (size_t k = 0; k < w_poly.size(); ++k) {
            l_poly[k] -= z_total * w_poly[k];
        }
        vela::detail::divide_by_linear(l_poly, y);
        transcript->send_to_verifier(detail::fflonk_label("C_wp"),
                                     ck.kzg_ck->commit(Polynomial<fr>(std::span<const fr>(l_poly))));
    }
};

/** @brief The mirror of `FflonkProver`: eight scalar multiplications and one pairing check. */
class FflonkVerifier {
  public:
    struct Claims {
        std::vector<size_t> group_num_columns;
        std::vector<Commitment> group_commitments;
        std::vector<FflonkColumnRef> unshifted;
        std::vector<fr> unshifted_evaluations;
        std::vector<FflonkColumnRef> to_be_shifted;
        std::vector<fr> shifted_evaluations;
    };

    template <typename Transcript>
    static bool verify(const FflonkConfig& /*config*/,
                       const Claims& claims,
                       std::span<const fr> r,
                       const std::shared_ptr<Transcript>& transcript)
    {
        const size_t num_groups = claims.group_num_columns.size();
        const bool has_shifted = !claims.to_be_shifted.empty();

        std::vector<Commitment> commitments = claims.group_commitments;
        if (commitments.empty()) {
            for (size_t g = 0; g < num_groups; ++g) {
                commitments.push_back(
                    transcript->template receive_from_prover<Commitment>(detail::fflonk_label("root", g, 0)));
            }
        }
        const fr rho = transcript->template get_challenge<fr>("FFLONK:rho");
        const fr lambda = transcript->template get_challenge<fr>("FFLONK:lambda");
        const Commitment c_h = transcript->template receive_from_prover<Commitment>(detail::fflonk_label("C_h"));

        const fr z = transcript->template get_challenge<fr>("FFLONK:z");
        if (z == fr::zero() || z == fr::one() || z == -fr::one()) {
            return false;
        }
        const fr z_inv = z.invert();

        std::vector<size_t> packs;
        std::vector<std::vector<fr>> p_evals(num_groups + 1);
        std::vector<std::vector<fr>> q_evals(num_groups + 1);
        for (size_t g = 0; g < num_groups; ++g) {
            packs.push_back(detail::next_power_of_two(claims.group_num_columns[g]));
            for (size_t c = 0; c < packs[g]; ++c) {
                p_evals[g].push_back(transcript->template receive_from_prover<fr>(detail::fflonk_label("p", g, c)));
                q_evals[g].push_back(transcript->template receive_from_prover<fr>(detail::fflonk_label("q", g, c)));
            }
        }
        packs.push_back(1);
        p_evals[num_groups].push_back(transcript->template receive_from_prover<fr>(detail::fflonk_label("h_z")));
        q_evals[num_groups].push_back(transcript->template receive_from_prover<fr>(detail::fflonk_label("h_z_inv")));

        // The multilinear claim, closed entirely in the field over the certified evaluations.
        fr v0_a = fr::zero();
        fr v1_a = fr::zero();
        fr v0_b = fr::zero();
        fr v1_b = fr::zero();
        fr y_combined = fr::zero();
        fr rho_power = fr::one();
        for (size_t i = 0; i < claims.unshifted.size(); ++i) {
            const auto& ref = claims.unshifted[i];
            v0_a += rho_power * p_evals[ref.group][ref.column];
            v1_a += rho_power * q_evals[ref.group][ref.column];
            y_combined += rho_power * claims.unshifted_evaluations[i];
            rho_power *= rho;
        }
        for (size_t l = 0; l < claims.to_be_shifted.size(); ++l) {
            const auto& ref = claims.to_be_shifted[l];
            v0_b += rho_power * p_evals[ref.group][ref.column];
            v1_b += rho_power * q_evals[ref.group][ref.column];
            y_combined += lambda * rho_power * claims.shifted_evaluations[l];
            rho_power *= rho;
        }
        const fr t_z = mercury::detail::tensor_poly_eval(r, z);
        const fr t_z_inv = mercury::detail::tensor_poly_eval(r, z_inv);
        fr g_z = v0_a;
        fr g_z_inv = v1_a;
        if (has_shifted) {
            g_z += lambda * z_inv * v0_b;
            g_z_inv += lambda * z * v1_b;
        }
        const fr w0 = p_evals[num_groups][0];
        const fr w1 = q_evals[num_groups][0];
        if (g_z * t_z_inv + g_z_inv * t_z - (y_combined + y_combined) != z * w0 + z_inv * w1) {
            return false;
        }

        // Shplonk: one accumulation over the group commitments, h, the quotient and the generator.
        const fr nu = transcript->template get_challenge<fr>("FFLONK:nu");
        const Commitment c_w = transcript->template receive_from_prover<Commitment>(detail::fflonk_label("C_w"));
        const fr y = transcript->template get_challenge<fr>("FFLONK:y");
        const Commitment c_wp = transcript->template receive_from_prover<Commitment>(detail::fflonk_label("C_wp"));

        fr z_total = fr::one();
        std::vector<fr> z_at_y(packs.size());
        for (size_t i = 0; i < packs.size(); ++i) {
            const fr y_pow = y.pow(uint256_t(uint64_t(packs[i])));
            z_at_y[i] = (y_pow - z) * (y_pow - z_inv);
            z_total *= z_at_y[i];
        }

        GroupElement accumulator = GroupElement::infinity();
        fr constant = fr::zero();
        fr nu_power = fr::one();
        for (size_t i = 0; i < packs.size(); ++i) {
            const fr scalar = nu_power * z_total * z_at_y[i].invert();
            const Commitment& commitment = i < num_groups ? commitments[i] : c_h;
            accumulator += GroupElement(commitment) * scalar;
            constant += scalar * detail::interpolant_at(p_evals[i], q_evals[i], z, z_inv, y, packs[i]);
            nu_power *= nu;
        }
        accumulator -= GroupElement(Commitment::one()) * constant;
        accumulator -= GroupElement(c_w) * z_total;
        accumulator += GroupElement(c_wp) * y;

        PairingPoints<Curve> pairing_points{ Commitment(accumulator), Commitment(-GroupElement(c_wp)) };
        return pairing_points.check();
    }
};

} // namespace bb::fflonk
