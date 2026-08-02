#pragma once

#include "barretenberg/commitment_schemes/commitment_key.hpp"
#include "barretenberg/commitment_schemes/mercury/mercury.hpp"
#include "barretenberg/commitment_schemes/whir/weights.hpp"
#include "barretenberg/common/assert.hpp"
#include "barretenberg/ecc/curves/bn254/bn254.hpp"
#include "barretenberg/ecc/curves/bn254/pairing.hpp"
#include "barretenberg/ecc/scalar_multiplication/scalar_multiplication.hpp"
#include "barretenberg/polynomials/polynomial.hpp"
#include "barretenberg/polynomials/polynomial_arithmetic.hpp"
#include "barretenberg/srs/global_crs.hpp"

#include <memory>
#include <span>
#include <string>
#include <vector>

namespace bb::chopin {

using Curve = curve::BN254;
using Commitment = Curve::AffineElement;
using GroupElement = Curve::Element;
using G2Affine = Curve::G2AffineElement;

/** @brief Reference to one polynomial of one committed group (duck-type shared with the suite). */
struct ChopinColumnRef {
    size_t group;
    size_t column;
};

/**
 * @brief CHOPIN parameters (eprint 2026/480): the bivariate split. A committed array of size
 * N = 2^n is the coefficient matrix of a bivariate polynomial f(X, Y) with X-degree < M1 (the low
 * `num_col_vars` variables) and Y-degree < M2 (the high variables); flat index k = i + j*M1.
 */
struct ChopinConfig {
    size_t num_variables;
    size_t num_col_vars; // m1: X dimension covers the low variables, as in Mercury

    size_t num_cols() const { return size_t(1) << num_col_vars; }                   // M1
    size_t num_rows() const { return size_t(1) << (num_variables - num_col_vars); } // M2

    static ChopinConfig create(size_t num_variables, size_t /*security_bits*/ = 100, size_t /*log_inv_rate*/ = 0)
    {
        BB_ASSERT_GT(num_variables, size_t(1));
        return { num_variables, (num_variables + 1) / 2 };
    }
};

/** @brief Prover-side commitment to a group: dense arrays and their bivariate-KZG commitments. */
struct ChopinGroupData {
    std::vector<std::vector<fr>> coefficients;
    std::vector<Commitment> commitments;

    size_t num_columns() const { return coefficients.size(); }
};

namespace detail {

inline std::string chopin_label(const std::string& name, size_t i)
{
    return "CHOPIN:" + name + "_" + std::to_string(i);
}

/** @brief In-place multiply by (X - root); the buffer has slack for the degree increase. */
inline void multiply_by_linear(std::vector<fr>& poly, const fr& root)
{
    for (size_t d = poly.size(); d-- > 1;) {
        poly[d] = poly[d - 1] - root * poly[d];
    }
    poly[0] = -root * poly[0];
}

/** @brief Coefficients of the Lagrange interpolant through up to four points. */
inline std::vector<fr> interpolant_coefficients(std::span<const fr> xs, std::span<const fr> ys)
{
    std::vector<fr> coefficients(xs.size(), fr::zero());
    for (size_t i = 0; i < xs.size(); ++i) {
        fr denominator = fr::one();
        for (size_t j = 0; j < xs.size(); ++j) {
            if (j != i) {
                denominator *= (xs[i] - xs[j]);
            }
        }
        const fr scale = ys[i] * denominator.invert();
        std::vector<fr> basis = { fr::one() };
        for (size_t j = 0; j < xs.size(); ++j) {
            if (j != i) {
                basis.push_back(fr::zero());
                for (size_t d = basis.size(); d-- > 1;) {
                    basis[d] = basis[d - 1] - xs[j] * basis[d];
                }
                basis[0] = -xs[j] * basis[0];
            }
        }
        for (size_t d = 0; d < basis.size(); ++d) {
            coefficients[d] += scale * basis[d];
        }
    }
    return coefficients;
}

} // namespace detail

/**
 * @brief CHOPIN commitment key: bivariate KZG over the grid H[i + j*M1] = tau_Y^j [tau^i]_1.
 *
 * @details tau_X is the ceremony trapdoor of bb's BN254 SRS (unknown, real). tau_Y is fixed to 2
 * so the grid derives from the ceremony X-row by pointwise doublings and [tau_Y]_2 = 2 [1]_2 is
 * public — a TEST-ONLY shortcut standing in for the second independent trapdoor a production
 * ceremony would generate. The honest prover/verifier cost profile is identical to a real
 * bivariate SRS: commitments and the opening quotient are flat size-N Pippenger MSMs over the
 * materialized grid, exactly matching the univariate-KZG commit path Mercury benchmarks against.
 */
class ChopinCommitmentKey {
  public:
    explicit ChopinCommitmentKey(const ChopinConfig& config)
        : config(config)
        , kzg_ck(std::make_shared<CommitmentKey<Curve>>(config.num_cols() + 8))
    {
        const size_t m1 = config.num_cols();
        const size_t m2 = config.num_rows();
        std::span<const Commitment> x_row = kzg_ck->get_monomial_points();
        grid.resize(m1 * m2);
        for (size_t i = 0; i < m1; ++i) {
            grid[i] = x_row[i];
        }
        std::vector<GroupElement> row(m1);
        for (size_t j = 1; j < m2; ++j) {
            for (size_t i = 0; i < m1; ++i) {
                row[i] = GroupElement(grid[(j - 1) * m1 + i]).dbl();
            }
            GroupElement::batch_normalize(row.data(), m1);
            for (size_t i = 0; i < m1; ++i) {
                grid[j * m1 + i] = Commitment(row[i].x, row[i].y);
            }
        }
        y_row.resize(m2 + 1);
        y_row[0] = Commitment::one();
        for (size_t j = 1; j <= m2; ++j) {
            y_row[j] = Commitment(GroupElement(y_row[j - 1]).dbl());
        }
    }

    Commitment commit_grid(std::span<const fr> coefficients) const
    {
        BB_ASSERT_LTE(coefficients.size(), grid.size());
        return Commitment(scalar_multiplication::pippenger_unsafe<Curve>(PolynomialSpan<const fr>(0, coefficients),
                                                                         std::span<const Commitment>(grid)));
    }

    Commitment commit_x(std::span<const fr> coefficients) const
    {
        size_t length = coefficients.size();
        while (length > 1 && coefficients[length - 1] == fr::zero()) {
            --length;
        }
        return kzg_ck->commit(Polynomial<fr>(coefficients.subspan(0, length)));
    }

    Commitment commit_y(std::span<const fr> coefficients) const
    {
        BB_ASSERT_LTE(coefficients.size(), y_row.size());
        return Commitment(scalar_multiplication::pippenger_unsafe<Curve>(PolynomialSpan<const fr>(0, coefficients),
                                                                         std::span<const Commitment>(y_row)));
    }

    ChopinGroupData commit_group(std::vector<std::vector<fr>> payload_columns,
                                 const std::vector<bool>& /*to_be_shifted*/ = {}) const
    {
        const size_t n = size_t(1) << config.num_variables;
        ChopinGroupData data;
        for (auto& column : payload_columns) {
            BB_ASSERT_LTE(column.size(), n, "polynomial too large for the configured size");
            column.resize(n, fr::zero());
            data.commitments.push_back(commit_grid(column));
            data.coefficients.push_back(std::move(column));
        }
        return data;
    }

    ChopinGroupData commit_group(std::span<const Polynomial<fr>* const> polynomials,
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

    ChopinConfig config;
    std::shared_ptr<CommitmentKey<Curve>> kzg_ck; // ceremony X-row and small-poly commitments
    std::vector<Commitment> grid;                 // H[i + j*M1] = 2^j [tau^i]_1
    std::vector<Commitment> y_row;                // [2^j]_1
};

/**
 * @brief One CHOPIN opening chain: chain A carries the unshifted rho-combination, chain B the
 * (unshifted!) to-be-shifted combination. The Honk shift lives in a verifier-side query
 * decomposition: with column restrictions w = F Psi_zR and w2 = F shr(Psi_zR),
 * shifted_eval = <w, shr(Psi_zL)> + eqL[M1-1] * w2[0], so chain B opens one extra column
 * restriction and the batch proof additionally opens it at 0.
 */
struct ChopinChain {
    std::vector<fr> array; // combined coefficient array, length N
    fr claimed_evaluation = fr::zero();
    bool is_shifted = false;

    std::vector<fr> w;       // F * Psi_zR, length M1
    std::vector<fr> w2;      // chain B only: F * shr(Psi_zR), length M1
    std::vector<fr> f_alpha; // per-column evaluations at alpha, length M2
};

/**
 * @brief Batched CHOPIN opening prover (eprint 2026/480 Fig. 5, extended to the two-chain Honk
 * claim set): column restrictions, one alpha-fold per chain, one gamma-batched Lagrange-IPA
 * witness S, one delta/z multi-polynomial multi-point univariate batch proof (Fig. 7), and one
 * mu-combined bivariate KZG opening at (alpha, beta) whose X-quotient is the single size-N MSM.
 */
class ChopinProver {
  public:
    struct Claims {
        std::vector<const ChopinGroupData*> groups;
        std::vector<ChopinColumnRef> unshifted;
        std::vector<fr> unshifted_evaluations;
        std::vector<ChopinColumnRef> to_be_shifted;
        std::vector<fr> shifted_evaluations;
        bool send_roots = true;
    };

    template <typename Transcript>
    static void prove(const ChopinCommitmentKey& ck,
                      const Claims& claims,
                      std::span<const fr> u,
                      const std::shared_ptr<Transcript>& transcript)
    {
        const ChopinConfig& config = ck.config;
        const size_t n = size_t(1) << config.num_variables;
        const size_t m1 = config.num_cols();
        const size_t m2 = config.num_rows();
        BB_ASSERT_EQ(u.size(), config.num_variables, "opening point size mismatch");

        if (claims.send_roots) {
            for (size_t g = 0; g < claims.groups.size(); ++g) {
                for (size_t c = 0; c < claims.groups[g]->num_columns(); ++c) {
                    transcript->send_to_verifier(detail::chopin_label("root", g) + "_" + std::to_string(c),
                                                 claims.groups[g]->commitments[c]);
                }
            }
        }
        const fr rho = transcript->template get_challenge<fr>("CHOPIN:rho");

        std::vector<ChopinChain> chains;
        {
            ChopinChain chain_a;
            chain_a.array.assign(n, fr::zero());
            fr rho_power = fr::one();
            for (size_t i = 0; i < claims.unshifted.size(); ++i) {
                const auto& column = claims.groups[claims.unshifted[i].group]->coefficients[claims.unshifted[i].column];
                for (size_t k = 0; k < n; ++k) {
                    chain_a.array[k] += rho_power * column[k];
                }
                chain_a.claimed_evaluation += rho_power * claims.unshifted_evaluations[i];
                rho_power *= rho;
            }
            chains.push_back(std::move(chain_a));
            if (!claims.to_be_shifted.empty()) {
                ChopinChain chain_b;
                chain_b.is_shifted = true;
                chain_b.array.assign(n, fr::zero());
                for (size_t l = 0; l < claims.to_be_shifted.size(); ++l) {
                    const auto& column =
                        claims.groups[claims.to_be_shifted[l].group]->coefficients[claims.to_be_shifted[l].column];
                    for (size_t k = 0; k < n; ++k) {
                        chain_b.array[k] += rho_power * column[k];
                    }
                    chain_b.claimed_evaluation += rho_power * claims.shifted_evaluations[l];
                    rho_power *= rho;
                }
                chains.push_back(std::move(chain_b));
            }
        }

        const std::vector<fr> eq_lo = whir::eq_tensor(u.subspan(0, config.num_col_vars));
        const std::vector<fr> eq_hi = whir::eq_tensor(u.subspan(config.num_col_vars));

        // Column restrictions w = F Psi_zR (and w2 = F shr(Psi_zR) for chain B).
        for (size_t c = 0; c < chains.size(); ++c) {
            ChopinChain& chain = chains[c];
            chain.w.assign(m1, fr::zero());
            if (chain.is_shifted) {
                chain.w2.assign(m1, fr::zero());
            }
            for (size_t j = 0; j < m2; ++j) {
                const fr* column = chain.array.data() + j * m1;
                for (size_t i = 0; i < m1; ++i) {
                    chain.w[i] += eq_hi[j] * column[i];
                    if (chain.is_shifted && j > 0) {
                        chain.w2[i] += eq_hi[j - 1] * column[i];
                    }
                }
            }
            transcript->send_to_verifier(detail::chopin_label("C0", c), ck.commit_x(chain.w));
            if (chain.is_shifted) {
                transcript->send_to_verifier(detail::chopin_label("C0s", c), ck.commit_x(chain.w2));
            }
        }

        const fr alpha = transcript->template get_challenge<fr>("CHOPIN:alpha");

        // Alpha-fold: f_alpha[j] = column_j(alpha), plus the evaluations tying w (and w2) to it.
        for (size_t c = 0; c < chains.size(); ++c) {
            ChopinChain& chain = chains[c];
            chain.f_alpha.resize(m2);
            for (size_t j = 0; j < m2; ++j) {
                chain.f_alpha[j] = polynomial_arithmetic::evaluate(chain.array.data() + j * m1, alpha, m1);
            }
            transcript->send_to_verifier(detail::chopin_label("C1", c), ck.commit_x(chain.f_alpha));
            transcript->send_to_verifier(detail::chopin_label("a", c),
                                         polynomial_arithmetic::evaluate(chain.w.data(), alpha, m1));
            if (chain.is_shifted) {
                transcript->send_to_verifier(detail::chopin_label("a2", c),
                                             polynomial_arithmetic::evaluate(chain.w2.data(), alpha, m1));
                transcript->send_to_verifier(detail::chopin_label("w20", c), chain.w2[0]);
            }
        }

        const fr gamma = transcript->template get_challenge<fr>("CHOPIN:gamma");

        // Gamma-batched Lagrange-IPA witness S (Lemma 2 of the paper; Mercury's accumulator).
        std::vector<fr> shr_lo(m1, fr::zero());
        for (size_t i = 1; i < m1; ++i) {
            shr_lo[i] = eq_lo[i - 1];
        }
        std::vector<fr> shr_hi(m2, fr::zero());
        for (size_t j = 1; j < m2; ++j) {
            shr_hi[j] = eq_hi[j - 1];
        }
        std::vector<fr> s_poly(std::max(m1, m2), fr::zero());
        fr gamma_power = fr::one();
        for (auto& chain : chains) {
            mercury::detail::accumulate_ip_terms(chain.w, chain.is_shifted ? shr_lo : eq_lo, gamma_power, s_poly);
            gamma_power *= gamma;
            mercury::detail::accumulate_ip_terms(chain.f_alpha, eq_hi, gamma_power, s_poly);
            gamma_power *= gamma;
            if (chain.is_shifted) {
                mercury::detail::accumulate_ip_terms(chain.f_alpha, shr_hi, gamma_power, s_poly);
                gamma_power *= gamma;
            }
        }
        transcript->send_to_verifier(std::string("CHOPIN:S"), ck.commit_x(s_poly));

        const fr beta = transcript->template get_challenge<fr>("CHOPIN:beta");
        const fr beta_inv = beta.invert();

        // Evaluations of the committed univariate polynomials at beta and 1/beta.
        for (size_t c = 0; c < chains.size(); ++c) {
            const ChopinChain& chain = chains[c];
            transcript->send_to_verifier(detail::chopin_label("w_beta", c),
                                         polynomial_arithmetic::evaluate(chain.w.data(), beta, m1));
            transcript->send_to_verifier(detail::chopin_label("w_binv", c),
                                         polynomial_arithmetic::evaluate(chain.w.data(), beta_inv, m1));
            transcript->send_to_verifier(detail::chopin_label("fa_beta", c),
                                         polynomial_arithmetic::evaluate(chain.f_alpha.data(), beta, m2));
            transcript->send_to_verifier(detail::chopin_label("fa_binv", c),
                                         polynomial_arithmetic::evaluate(chain.f_alpha.data(), beta_inv, m2));
        }
        transcript->send_to_verifier(std::string("CHOPIN:s_beta"),
                                     polynomial_arithmetic::evaluate(s_poly.data(), beta, s_poly.size()));
        transcript->send_to_verifier(std::string("CHOPIN:s_binv"),
                                     polynomial_arithmetic::evaluate(s_poly.data(), beta_inv, s_poly.size()));

        // Univariate multi-polynomial multi-point batch proof (paper Fig. 7).
        const bool has_shifted = chains.size() > 1;
        const auto [batch_polys, batch_sets] = batch_layout(chains, s_poly, alpha, beta, beta_inv, has_shifted);
        std::vector<fr> t_points = { alpha, beta, beta_inv };
        if (has_shifted) {
            t_points.push_back(fr::zero());
        }

        const fr delta = transcript->template get_challenge<fr>("CHOPIN:delta");
        const size_t max_len = std::max(m1, m2) + t_points.size();
        std::vector<fr> p_poly(max_len, fr::zero());
        fr delta_power = fr::one();
        for (size_t t = 0; t < batch_polys.size(); ++t) {
            const std::vector<fr>& poly = *batch_polys[t];
            const std::vector<fr>& xs = batch_sets[t];
            std::vector<fr> ys;
            for (const fr& x : xs) {
                ys.push_back(polynomial_arithmetic::evaluate(poly.data(), x, poly.size()));
            }
            std::vector<fr> term(poly.begin(), poly.end());
            const std::vector<fr> ell = detail::interpolant_coefficients(xs, ys);
            for (size_t d = 0; d < ell.size(); ++d) {
                term[d] -= ell[d];
            }
            term.resize(max_len, fr::zero());
            for (const fr& point : t_points) {
                if (std::find(xs.begin(), xs.end(), point) == xs.end()) {
                    detail::multiply_by_linear(term, point);
                }
            }
            for (size_t d = 0; d < term.size(); ++d) {
                p_poly[d] += delta_power * term[d];
            }
            delta_power *= delta;
        }
        for (const fr& point : t_points) {
            mercury::detail::divide_by_linear(p_poly, point);
        }
        transcript->send_to_verifier(std::string("CHOPIN:W"), ck.commit_x(p_poly));

        const fr z = transcript->template get_challenge<fr>("CHOPIN:z");
        std::vector<fr> s_batch(std::max(m1, m2), fr::zero());
        delta_power = fr::one();
        for (size_t t = 0; t < batch_polys.size(); ++t) {
            const std::vector<fr>& poly = *batch_polys[t];
            const std::vector<fr>& xs = batch_sets[t];
            std::vector<fr> ys;
            for (const fr& x : xs) {
                ys.push_back(polynomial_arithmetic::evaluate(poly.data(), x, poly.size()));
            }
            fr z_complement = fr::one();
            for (const fr& point : t_points) {
                if (std::find(xs.begin(), xs.end(), point) == xs.end()) {
                    z_complement *= (z - point);
                }
            }
            const fr ell_z = mercury::detail::interpolate_at(xs, ys, z);
            const fr scale = delta_power * z_complement;
            for (size_t d = 0; d < poly.size(); ++d) {
                s_batch[d] += scale * poly[d];
            }
            s_batch[0] -= scale * ell_z;
            delta_power *= delta;
        }
        const fr y_value = polynomial_arithmetic::evaluate(s_batch.data(), z, s_batch.size());
        transcript->send_to_verifier(std::string("CHOPIN:y"), y_value);

        s_batch[0] -= y_value;
        mercury::detail::divide_by_linear(s_batch, z);
        transcript->send_to_verifier(std::string("CHOPIN:pi_s"), ck.commit_x(s_batch));

        fr z_t = fr::one();
        for (const fr& point : t_points) {
            z_t *= (z - point);
        }
        std::vector<fr> q_sh(p_poly.begin(), p_poly.end());
        for (fr& coefficient : q_sh) {
            coefficient *= z_t;
        }
        q_sh[0] -= y_value;
        mercury::detail::divide_by_linear(q_sh, z);
        transcript->send_to_verifier(std::string("CHOPIN:pi_q"), ck.commit_x(q_sh));

        // Mu-combined bivariate KZG opening of the chains at (alpha, beta): the single size-N MSM.
        const fr mu = transcript->template get_challenge<fr>("CHOPIN:mu");
        std::vector<fr> combined = chains[0].array;
        std::vector<fr> combined_alpha = chains[0].f_alpha;
        if (has_shifted) {
            for (size_t k = 0; k < n; ++k) {
                combined[k] += mu * chains[1].array[k];
            }
            for (size_t j = 0; j < m2; ++j) {
                combined_alpha[j] += mu * chains[1].f_alpha[j];
            }
        }
        // q1(X, Y): per-column exact division of (column - column(alpha)) by (X - alpha).
        for (size_t j = 0; j < m2; ++j) {
            fr* column = combined.data() + j * m1;
            column[0] -= combined_alpha[j];
            fr acc = fr::zero();
            for (size_t i = m1; i-- > 0;) {
                const fr coefficient = column[i];
                column[i] = acc;
                acc = coefficient + alpha * acc;
            }
            BB_ASSERT_EQ(acc, fr::zero(), "bivariate X-division not exact");
        }
        transcript->send_to_verifier(std::string("CHOPIN:Q1"), ck.commit_grid(combined));

        const fr v_combined = polynomial_arithmetic::evaluate(combined_alpha.data(), beta, m2);
        combined_alpha[0] -= v_combined;
        mercury::detail::divide_by_linear(combined_alpha, beta);
        transcript->send_to_verifier(std::string("CHOPIN:Q2"), ck.commit_y(combined_alpha));

        [[maybe_unused]] const fr pairing_batch = transcript->template get_challenge<fr>("CHOPIN:pairing_batch");
    }

  private:
    /** @brief The batch-proof polynomial list and per-polynomial opening sets, in transcript order. */
    static std::pair<std::vector<const std::vector<fr>*>, std::vector<std::vector<fr>>> batch_layout(
        const std::vector<ChopinChain>& chains,
        const std::vector<fr>& s_poly,
        const fr& alpha,
        const fr& beta,
        const fr& beta_inv,
        bool has_shifted)
    {
        std::vector<const std::vector<fr>*> polys;
        std::vector<std::vector<fr>> sets;
        polys.push_back(&chains[0].w);
        sets.push_back({ alpha, beta, beta_inv });
        polys.push_back(&chains[0].f_alpha);
        sets.push_back({ beta, beta_inv });
        if (has_shifted) {
            polys.push_back(&chains[1].w);
            sets.push_back({ alpha, beta, beta_inv });
            polys.push_back(&chains[1].w2);
            sets.push_back({ fr::zero(), alpha });
            polys.push_back(&chains[1].f_alpha);
            sets.push_back({ beta, beta_inv });
        }
        polys.push_back(&s_poly);
        sets.push_back({ beta, beta_inv });
        return { polys, sets };
    }
};

/** @brief Batched CHOPIN opening verifier; the mirror of `ChopinProver`. */
class ChopinVerifier {
  public:
    struct Claims {
        std::vector<size_t> group_num_columns;
        std::vector<std::vector<Commitment>> group_commitments; // when non-empty, transcript-bound
        std::vector<ChopinColumnRef> unshifted;
        std::vector<fr> unshifted_evaluations;
        std::vector<ChopinColumnRef> to_be_shifted;
        std::vector<fr> shifted_evaluations;
    };

    template <typename Transcript>
    static bool verify(const ChopinConfig& config,
                       const Claims& claims,
                       std::span<const fr> u,
                       const std::shared_ptr<Transcript>& transcript)
    {
        const size_t m1 = config.num_cols();
        const size_t m2 = config.num_rows();
        BB_ASSERT_EQ(u.size(), config.num_variables, "opening point size mismatch");

        std::vector<std::vector<Commitment>> group_commitments = claims.group_commitments;
        if (group_commitments.empty()) {
            for (size_t g = 0; g < claims.group_num_columns.size(); ++g) {
                std::vector<Commitment> commitments;
                for (size_t c = 0; c < claims.group_num_columns[g]; ++c) {
                    commitments.push_back(transcript->template receive_from_prover<Commitment>(
                        detail::chopin_label("root", g) + "_" + std::to_string(c)));
                }
                group_commitments.push_back(std::move(commitments));
            }
        }
        const fr rho = transcript->template get_challenge<fr>("CHOPIN:rho");

        struct ChainView {
            GroupElement commitment = GroupElement::infinity();
            fr claimed_evaluation = fr::zero();
            bool is_shifted = false;
            Commitment c0, c0s, c1;
            fr a, a2, w20;
            fr w_beta, w_binv, fa_beta, fa_binv;
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
            chains[1].is_shifted = true;
            for (size_t l = 0; l < claims.to_be_shifted.size(); ++l) {
                chains[1].commitment +=
                    GroupElement(group_commitments[claims.to_be_shifted[l].group][claims.to_be_shifted[l].column]) *
                    rho_power;
                chains[1].claimed_evaluation += rho_power * claims.shifted_evaluations[l];
                rho_power *= rho;
            }
        }
        const bool has_shifted = chains.size() > 1;

        const std::vector<fr> eq_lo = whir::eq_tensor(u.subspan(0, config.num_col_vars));
        const std::vector<fr> eq_hi = whir::eq_tensor(u.subspan(config.num_col_vars));

        for (size_t c = 0; c < chains.size(); ++c) {
            chains[c].c0 = transcript->template receive_from_prover<Commitment>(detail::chopin_label("C0", c));
            if (chains[c].is_shifted) {
                chains[c].c0s = transcript->template receive_from_prover<Commitment>(detail::chopin_label("C0s", c));
            }
        }
        const fr alpha = transcript->template get_challenge<fr>("CHOPIN:alpha");
        for (size_t c = 0; c < chains.size(); ++c) {
            chains[c].c1 = transcript->template receive_from_prover<Commitment>(detail::chopin_label("C1", c));
            chains[c].a = transcript->template receive_from_prover<fr>(detail::chopin_label("a", c));
            if (chains[c].is_shifted) {
                chains[c].a2 = transcript->template receive_from_prover<fr>(detail::chopin_label("a2", c));
                chains[c].w20 = transcript->template receive_from_prover<fr>(detail::chopin_label("w20", c));
            }
        }
        const fr gamma = transcript->template get_challenge<fr>("CHOPIN:gamma");
        const Commitment c_s = transcript->template receive_from_prover<Commitment>(std::string("CHOPIN:S"));
        const fr beta = transcript->template get_challenge<fr>("CHOPIN:beta");
        const fr beta_inv = beta.invert();
        for (size_t c = 0; c < chains.size(); ++c) {
            chains[c].w_beta = transcript->template receive_from_prover<fr>(detail::chopin_label("w_beta", c));
            chains[c].w_binv = transcript->template receive_from_prover<fr>(detail::chopin_label("w_binv", c));
            chains[c].fa_beta = transcript->template receive_from_prover<fr>(detail::chopin_label("fa_beta", c));
            chains[c].fa_binv = transcript->template receive_from_prover<fr>(detail::chopin_label("fa_binv", c));
        }
        const fr s_beta = transcript->template receive_from_prover<fr>(std::string("CHOPIN:s_beta"));
        const fr s_binv = transcript->template receive_from_prover<fr>(std::string("CHOPIN:s_binv"));

        // Weight-polynomial evaluations: psi twins and their shifted variants.
        const fr psi_lo_beta = mercury::detail::tensor_poly_eval(u.subspan(0, config.num_col_vars), beta);
        const fr psi_lo_binv = mercury::detail::tensor_poly_eval(u.subspan(0, config.num_col_vars), beta_inv);
        const fr psi_hi_beta = mercury::detail::tensor_poly_eval(u.subspan(config.num_col_vars), beta);
        const fr psi_hi_binv = mercury::detail::tensor_poly_eval(u.subspan(config.num_col_vars), beta_inv);
        const fr beta_m1 = beta.pow(uint256_t(uint64_t(m1)));
        const fr binv_m1 = beta_inv.pow(uint256_t(uint64_t(m1)));
        const fr beta_m2 = beta.pow(uint256_t(uint64_t(m2)));
        const fr binv_m2 = beta_inv.pow(uint256_t(uint64_t(m2)));
        const fr shr_lo_beta = beta * psi_lo_beta - eq_lo[m1 - 1] * beta_m1;
        const fr shr_lo_binv = beta_inv * psi_lo_binv - eq_lo[m1 - 1] * binv_m1;
        const fr shr_hi_beta = beta * psi_hi_beta - eq_hi[m2 - 1] * beta_m2;
        const fr shr_hi_binv = beta_inv * psi_hi_binv - eq_hi[m2 - 1] * binv_m2;

        // Check 1: the gamma-batched Lagrange-IPA identity at beta (paper Fig. 4 step 5).
        fr lhs = fr::zero();
        fr constant = fr::zero();
        fr gamma_power = fr::one();
        for (const ChainView& chain : chains) {
            const fr w_weight_beta = chain.is_shifted ? shr_lo_beta : psi_lo_beta;
            const fr w_weight_binv = chain.is_shifted ? shr_lo_binv : psi_lo_binv;
            const fr w_value = chain.is_shifted
                                   ? chain.claimed_evaluation - eq_lo[m1 - 1] * chain.w20 // <w, shr(Psi_zL)>
                                   : chain.claimed_evaluation;
            lhs += gamma_power * (chain.w_beta * w_weight_binv + chain.w_binv * w_weight_beta);
            constant += gamma_power * w_value;
            gamma_power *= gamma;
            lhs += gamma_power * (chain.fa_beta * psi_hi_binv + chain.fa_binv * psi_hi_beta);
            constant += gamma_power * chain.a;
            gamma_power *= gamma;
            if (chain.is_shifted) {
                lhs += gamma_power * (chain.fa_beta * shr_hi_binv + chain.fa_binv * shr_hi_beta);
                constant += gamma_power * chain.a2;
                gamma_power *= gamma;
            }
        }
        if (lhs != constant + constant + beta * s_beta + beta_inv * s_binv) {
            return false;
        }

        // Batch-proof bookkeeping mirroring the prover's layout.
        std::vector<Commitment> batch_commitments;
        std::vector<std::vector<fr>> batch_xs;
        std::vector<std::vector<fr>> batch_ys;
        for (const ChainView& chain : chains) {
            batch_commitments.push_back(chain.c0);
            batch_xs.push_back({ alpha, beta, beta_inv });
            batch_ys.push_back({ chain.a, chain.w_beta, chain.w_binv });
            if (chain.is_shifted) {
                batch_commitments.push_back(chain.c0s);
                batch_xs.push_back({ fr::zero(), alpha });
                batch_ys.push_back({ chain.w20, chain.a2 });
            }
            batch_commitments.push_back(chain.c1);
            batch_xs.push_back({ beta, beta_inv });
            batch_ys.push_back({ chain.fa_beta, chain.fa_binv });
        }
        batch_commitments.push_back(c_s);
        batch_xs.push_back({ beta, beta_inv });
        batch_ys.push_back({ s_beta, s_binv });
        // The prover's layout orders chain B as (w, w2, f_alpha); ours pushed (w, w2, f_alpha) too.
        std::vector<fr> t_points = { alpha, beta, beta_inv };
        if (has_shifted) {
            t_points.push_back(fr::zero());
        }

        const fr delta = transcript->template get_challenge<fr>("CHOPIN:delta");
        const Commitment w_commitment = transcript->template receive_from_prover<Commitment>(std::string("CHOPIN:W"));
        const fr z = transcript->template get_challenge<fr>("CHOPIN:z");
        const fr y_value = transcript->template receive_from_prover<fr>(std::string("CHOPIN:y"));
        const Commitment pi_s = transcript->template receive_from_prover<Commitment>(std::string("CHOPIN:pi_s"));
        const Commitment pi_q = transcript->template receive_from_prover<Commitment>(std::string("CHOPIN:pi_q"));
        const fr mu = transcript->template get_challenge<fr>("CHOPIN:mu");
        const Commitment q1 = transcript->template receive_from_prover<Commitment>(std::string("CHOPIN:Q1"));
        const Commitment q2 = transcript->template receive_from_prover<Commitment>(std::string("CHOPIN:Q2"));
        const fr xi = transcript->template get_challenge<fr>("CHOPIN:pairing_batch");

        // C_batch = sum_t delta^{t-1} Z_{T \ S_t}(z) (C_t - [ell_t(z)]_1).
        GroupElement c_batch = GroupElement::infinity();
        fr scalar_accumulator = fr::zero();
        fr delta_power = fr::one();
        fr z_t = fr::one();
        for (const fr& point : t_points) {
            z_t *= (z - point);
        }
        for (size_t t = 0; t < batch_commitments.size(); ++t) {
            fr z_complement = fr::one();
            for (const fr& point : t_points) {
                if (std::find(batch_xs[t].begin(), batch_xs[t].end(), point) == batch_xs[t].end()) {
                    z_complement *= (z - point);
                }
            }
            const fr ell_z = mercury::detail::interpolate_at(batch_xs[t], batch_ys[t], z);
            const fr scale = delta_power * z_complement;
            c_batch += GroupElement(batch_commitments[t]) * scale;
            scalar_accumulator += scale * ell_z;
            delta_power *= delta;
        }
        c_batch -= GroupElement(Commitment::one()) * scalar_accumulator;

        // Pairing terms, xi-batched into e(A1, [1]) e(A_tau, [tau]) e(A_tau_y, [2]_2) == 1.
        GroupElement a_one = GroupElement::infinity();
        GroupElement a_tau = GroupElement::infinity();
        GroupElement a_tau_y = GroupElement::infinity();
        fr xi_power = fr::one();
        // (i) e(C_batch - y [1] + z pi_s, [1]) = e(pi_s, [tau]).
        a_one += (c_batch - GroupElement(Commitment::one()) * y_value + GroupElement(pi_s) * z) * xi_power;
        a_tau -= GroupElement(pi_s) * xi_power;
        xi_power *= xi;
        // (ii) e(Z_T(z) W - y [1] + z pi_q, [1]) = e(pi_q, [tau]).
        a_one +=
            (GroupElement(w_commitment) * z_t - GroupElement(Commitment::one()) * y_value + GroupElement(pi_q) * z) *
            xi_power;
        a_tau -= GroupElement(pi_q) * xi_power;
        xi_power *= xi;
        // (iii) bivariate: e(C_mu - v [1] + alpha Q1 + beta Q2, [1]) = e(Q1, [tau]) e(Q2, [2]_2).
        GroupElement c_mu = chains[0].commitment;
        fr v_combined = chains[0].fa_beta;
        if (has_shifted) {
            c_mu += chains[1].commitment * mu;
            v_combined += mu * chains[1].fa_beta;
        }
        a_one +=
            (c_mu - GroupElement(Commitment::one()) * v_combined + GroupElement(q1) * alpha + GroupElement(q2) * beta) *
            xi_power;
        a_tau -= GroupElement(q1) * xi_power;
        a_tau_y -= GroupElement(q2) * xi_power;

        const G2Affine g2_one(bb::g2::affine_one);
        const G2Affine g2_tau = srs::get_crs_factory<Curve>()->get_verifier_crs()->get_g2x();
        const G2Affine g2_two(bb::g2::element(bb::g2::affine_one) + bb::g2::element(bb::g2::affine_one));
        std::array<Commitment, 3> pairing_g1{ Commitment(a_one), Commitment(a_tau), Commitment(a_tau_y) };
        std::array<G2Affine, 3> pairing_g2{ g2_one, g2_tau, g2_two };
        return pairing::reduced_ate_pairing_batch(pairing_g1.data(), pairing_g2.data(), 3) == fq12::one();
    }
};

} // namespace bb::chopin
