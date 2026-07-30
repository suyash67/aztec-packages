#pragma once

#include "barretenberg/commitment_schemes/commitment_key.hpp"
#include "barretenberg/commitment_schemes/pairing_points.hpp"
#include "barretenberg/commitment_schemes/whir/weights.hpp"
#include "barretenberg/common/assert.hpp"
#include "barretenberg/ecc/curves/bn254/bn254.hpp"
#include "barretenberg/polynomials/polynomial.hpp"
#include "barretenberg/polynomials/polynomial_arithmetic.hpp"

#include <map>
#include <memory>
#include <span>
#include <string>
#include <vector>

namespace bb::mercury {

using Curve = curve::BN254;
using Commitment = Curve::AffineElement;

/** @brief Reference to one polynomial of one committed group (duck-type shared with WHIR/Ligero). */
struct MercuryColumnRef {
    size_t group;
    size_t column;
};

/** @brief Mercury parameters: the matrix split. See README.md §1. */
struct MercuryConfig {
    size_t num_variables; // n: polynomials have 2^n evaluations
    size_t num_col_vars;  // t = ceil(n/2): matrix is 2^{n-t} rows x 2^t cols

    size_t num_cols() const { return size_t(1) << num_col_vars; }
    size_t num_rows() const { return size_t(1) << (num_variables - num_col_vars); }

    static MercuryConfig create(size_t num_variables, size_t /*security_bits*/ = 100, size_t /*log_inv_rate*/ = 0)
    {
        BB_ASSERT_GT(num_variables, size_t(1));
        return { num_variables, (num_variables + 1) / 2 };
    }
};

/** @brief Prover-side commitment to a group: dense arrays and their individual KZG commitments. */
struct MercuryGroupData {
    std::vector<std::vector<fr>> coefficients;
    std::vector<Commitment> commitments;

    size_t num_columns() const { return coefficients.size(); }
};

namespace detail {

inline std::string mercury_label(const std::string& name, size_t i)
{
    return "MERCURY:" + name + "_" + std::to_string(i);
}

/** @brief P_u(x) = prod_i (u_i x^{2^i} + 1 - u_i): the eq-tensor polynomial evaluated at x. */
inline fr tensor_poly_eval(std::span<const fr> u, const fr& x)
{
    fr result = fr::one();
    fr power = x;
    for (const fr& u_i : u) {
        result *= (u_i * power + fr::one() - u_i);
        power = power.sqr();
    }
    return result;
}

/** @brief In-place exact division by (X - root); the caller guarantees divisibility. */
inline void divide_by_linear(std::vector<fr>& poly, const fr& root)
{
    fr acc = fr::zero();
    for (size_t k = poly.size(); k-- > 0;) {
        const fr coefficient = poly[k];
        poly[k] = acc;
        acc = coefficient + root * acc;
    }
    BB_ASSERT_EQ(acc, fr::zero(), "division by (X - root) is not exact");
}

/**
 * @brief Accumulate the batched inner-product polynomial S (README.md §1):
 * S[k-1] += scale * ([X^k] + [X^{-k}]) of p(X)r(1/X) + p(1/X)r(X), for k >= 1.
 */
inline void accumulate_ip_terms(std::span<const fr> p, std::span<const fr> r, const fr& scale, std::vector<fr>& s)
{
    BB_ASSERT_EQ(p.size(), r.size());
    const size_t len = p.size();
    for (size_t shift = 1; shift < len; ++shift) {
        fr sum = fr::zero();
        for (size_t i = shift; i < len; ++i) {
            sum += p[i] * r[i - shift] + p[i - shift] * r[i];
        }
        // p(X)r(1/X) contributes sum at X^{+shift} and (by symmetry of the added mirror term)
        // the same value at X^{-shift}; both land in S per the identity.
        s[shift - 1] += scale * sum;
    }
}

/** @brief Lagrange interpolation value at z through up to three points. */
inline fr interpolate_at(std::span<const fr> xs, std::span<const fr> ys, const fr& z)
{
    fr result = fr::zero();
    for (size_t i = 0; i < xs.size(); ++i) {
        fr term = ys[i];
        for (size_t j = 0; j < xs.size(); ++j) {
            if (j != i) {
                term *= (z - xs[j]) * (xs[i] - xs[j]).invert();
            }
        }
        result += term;
    }
    return result;
}

} // namespace detail

/**
 * @brief One Mercury fold chain (README.md §2): the per-combined-polynomial data flowing through
 * the protocol. Chain A carries the unshifted batch (standard fold identity); chain B carries the
 * shifted batch as the virtual polynomial B(X)/X (fold identity multiplied through by X).
 */
struct MercuryChain {
    std::vector<fr> array;   // the (virtual) combined coefficient array, length N
    std::vector<fr> raw;     // chain B only: the unshifted combined array (X * array)
    bool is_shifted = false; // selects the pairing form
    fr claimed_evaluation = fr::zero();

    std::vector<fr> h; // column-combined polynomial, length rows
    std::vector<fr> g; // fold remainder, length b
    std::vector<fr> q; // fold quotient (chain B: Q = X*q_B), length N
    fr h_alpha = fr::zero();
};

/** @brief Transparent-shell-compatible commitment key wrapping bb's KZG commitment key. */
class MercuryCommitmentKey {
  public:
    explicit MercuryCommitmentKey(const MercuryConfig& config)
        : config(config)
        , kzg_ck(std::make_shared<CommitmentKey<Curve>>(size_t(1) << config.num_variables))
    {}

    MercuryGroupData commit_group(std::vector<std::vector<fr>> payload_columns,
                                  const std::vector<bool>& /*to_be_shifted*/ = {}) const
    {
        const size_t n = size_t(1) << config.num_variables;
        MercuryGroupData data;
        for (auto& column : payload_columns) {
            BB_ASSERT_LTE(column.size(), n, "polynomial too large for the configured size");
            data.commitments.push_back(kzg_ck->commit(Polynomial<fr>(std::span<const fr>(column))));
            column.resize(n, fr::zero());
            data.coefficients.push_back(std::move(column));
        }
        return data;
    }

    MercuryGroupData commit_group(std::span<const Polynomial<fr>* const> polynomials,
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

    MercuryConfig config;
    std::shared_ptr<CommitmentKey<Curve>> kzg_ck;
};

/**
 * @brief Batched Mercury opening prover: two fold chains sharing challenges, one batched
 * inner-product polynomial, one BDFG multi-point opening (README.md §2).
 */
class MercuryProver {
  public:
    struct Claims {
        std::vector<const MercuryGroupData*> groups;
        std::vector<MercuryColumnRef> unshifted;
        std::vector<fr> unshifted_evaluations;
        std::vector<MercuryColumnRef> to_be_shifted;
        std::vector<fr> shifted_evaluations;
        bool send_roots = true;
    };

    template <typename Transcript>
    static void prove(const MercuryCommitmentKey& ck,
                      const Claims& claims,
                      std::span<const fr> u,
                      const std::shared_ptr<Transcript>& transcript)
    {
        const MercuryConfig& config = ck.config;
        const size_t n = size_t(1) << config.num_variables;
        const size_t b = config.num_cols();
        const size_t rows = config.num_rows();
        BB_ASSERT_EQ(u.size(), config.num_variables, "opening point size mismatch");

        if (claims.send_roots) {
            for (size_t g = 0; g < claims.groups.size(); ++g) {
                for (size_t c = 0; c < claims.groups[g]->num_columns(); ++c) {
                    transcript->send_to_verifier(detail::mercury_label("root", g) + "_" + std::to_string(c),
                                                 claims.groups[g]->commitments[c]);
                }
            }
        }
        const fr rho = transcript->template get_challenge<fr>("MERCURY:rho");

        // Combined arrays: chain A directly; chain B as the shifted virtual polynomial raw/X.
        std::vector<MercuryChain> chains;
        {
            MercuryChain chain_a;
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
                MercuryChain chain_b;
                chain_b.is_shifted = true;
                chain_b.raw.assign(n, fr::zero());
                for (size_t l = 0; l < claims.to_be_shifted.size(); ++l) {
                    const auto& column =
                        claims.groups[claims.to_be_shifted[l].group]->coefficients[claims.to_be_shifted[l].column];
                    BB_ASSERT_EQ(column[0], fr::zero(), "to-be-shifted polynomial must have zero constant term");
                    for (size_t k = 0; k < n; ++k) {
                        chain_b.raw[k] += rho_power * column[k];
                    }
                    chain_b.claimed_evaluation += rho_power * claims.shifted_evaluations[l];
                    rho_power *= rho;
                }
                chain_b.array.assign(n, fr::zero());
                for (size_t k = 0; k + 1 < n; ++k) {
                    chain_b.array[k] = chain_b.raw[k + 1];
                }
                chains.push_back(std::move(chain_b));
            }
        }

        // Round 1: column-combined polynomials h = sum_j eq_j(u_lo) f_j.
        const std::vector<fr> eq_lo = whir::eq_tensor(u.subspan(0, config.num_col_vars));
        for (size_t c = 0; c < chains.size(); ++c) {
            MercuryChain& chain = chains[c];
            chain.h.assign(rows, fr::zero());
            for (size_t row = 0; row < rows; ++row) {
                fr acc = fr::zero();
                const fr* slice = chain.array.data() + row * b;
                for (size_t j = 0; j < b; ++j) {
                    acc += eq_lo[j] * slice[j];
                }
                chain.h[row] = acc;
            }
            transcript->send_to_verifier(detail::mercury_label("h", c), commit(ck, chain.h));
        }

        // Round 2: fold at alpha. f = q (X^b - alpha) + g; chain B commits Q = X q_B.
        const fr alpha = transcript->template get_challenge<fr>("MERCURY:alpha");
        for (size_t c = 0; c < chains.size(); ++c) {
            MercuryChain& chain = chains[c];
            std::vector<fr> quotient(n, fr::zero());
            std::vector<fr> remainder(chain.array.begin(), chain.array.begin() + static_cast<std::ptrdiff_t>(b));
            for (size_t k = n; k-- > b;) {
                // f_k = q_{k-b} - alpha q_k  =>  q_{k-b} = f_k + alpha q_k
                quotient[k - b] = chain.array[k] + alpha * quotient[k];
            }
            for (size_t k = 0; k < b; ++k) {
                remainder[k] += alpha * quotient[k];
            }
            chain.g = std::move(remainder);
            if (chain.is_shifted) {
                // Q = X * q_B
                chain.q.assign(n, fr::zero());
                for (size_t k = 0; k + 1 < n; ++k) {
                    chain.q[k + 1] = quotient[k];
                }
            } else {
                chain.q = std::move(quotient);
            }
            chain.h_alpha = polynomial_arithmetic::evaluate(chain.h.data(), alpha, chain.h.size());
            transcript->send_to_verifier(detail::mercury_label("g", c), commit(ck, chain.g));
            transcript->send_to_verifier(detail::mercury_label("q", c), commit(ck, chain.q));
            transcript->send_to_verifier(detail::mercury_label("h_alpha", c), chain.h_alpha);
        }

        // Round 3: batched inner-product polynomial S and degree-bound reversals D.
        const fr gamma = transcript->template get_challenge<fr>("MERCURY:gamma");
        const std::vector<fr> p_hi = whir::eq_tensor(u.subspan(config.num_col_vars));
        std::vector<fr> s_poly(std::max(b, rows), fr::zero());
        fr gamma_power = fr::one();
        for (auto& chain : chains) {
            detail::accumulate_ip_terms(chain.g, eq_lo, gamma_power, s_poly);
            gamma_power *= gamma;
            detail::accumulate_ip_terms(chain.h, p_hi, gamma_power, s_poly);
            gamma_power *= gamma;
        }
        transcript->send_to_verifier(std::string("MERCURY:S"), commit(ck, s_poly));
        for (size_t c = 0; c < chains.size(); ++c) {
            std::vector<fr> reversal(b);
            for (size_t k = 0; k < b; ++k) {
                reversal[k] = chains[c].g[b - 1 - k];
            }
            transcript->send_to_verifier(detail::mercury_label("D", c), commit(ck, reversal));
        }

        // Round 4: evaluations at zeta and 1/zeta; fold-consistency quotients H.
        const fr zeta = transcript->template get_challenge<fr>("MERCURY:zeta");
        const fr zeta_inv = zeta.invert();
        const fr zeta_pow_b = zeta.pow(uint256_t(uint64_t(b)));
        for (size_t c = 0; c < chains.size(); ++c) {
            MercuryChain& chain = chains[c];
            const fr g_zeta = polynomial_arithmetic::evaluate(chain.g.data(), zeta, b);
            const fr g_zeta_inv = polynomial_arithmetic::evaluate(chain.g.data(), zeta_inv, b);
            const fr h_zeta = polynomial_arithmetic::evaluate(chain.h.data(), zeta, rows);
            const fr h_zeta_inv = polynomial_arithmetic::evaluate(chain.h.data(), zeta_inv, rows);
            transcript->send_to_verifier(detail::mercury_label("g_zeta", c), g_zeta);
            transcript->send_to_verifier(detail::mercury_label("g_zeta_inv", c), g_zeta_inv);
            transcript->send_to_verifier(detail::mercury_label("h_zeta", c), h_zeta);
            transcript->send_to_verifier(detail::mercury_label("h_zeta_inv", c), h_zeta_inv);
            transcript->send_to_verifier(detail::mercury_label("D_zeta", c),
                                         zeta.pow(uint256_t(uint64_t(b - 1))) * g_zeta_inv);

            // H = (base - (zeta^b - alpha) q - remainder_term) / (X - zeta), where base/remainder
            // follow the chain's pairing form (README.md §2).
            std::vector<fr> numerator = chain.is_shifted ? chain.raw : chain.array;
            const fr factor = zeta_pow_b - alpha;
            for (size_t k = 0; k < n; ++k) {
                numerator[k] -= factor * chain.q[k];
            }
            numerator[0] -= chain.is_shifted ? zeta * g_zeta : g_zeta;
            detail::divide_by_linear(numerator, zeta);
            transcript->send_to_verifier(detail::mercury_label("H", c), commit(ck, numerator));
        }
        const fr s_zeta = polynomial_arithmetic::evaluate(s_poly.data(), zeta, s_poly.size());
        const fr s_zeta_inv = polynomial_arithmetic::evaluate(s_poly.data(), zeta_inv, s_poly.size());
        transcript->send_to_verifier(std::string("MERCURY:S_zeta"), s_zeta);
        transcript->send_to_verifier(std::string("MERCURY:S_zeta_inv"), s_zeta_inv);

        // Round 5: BDFG batch quotient q_m over T = {zeta, 1/zeta, alpha}.
        const fr beta = transcript->template get_challenge<fr>("MERCURY:beta");
        // Small-poly list with per-poly opening sets: g_c, h_c, D_c per chain, then S.
        std::vector<std::vector<fr>> smalls;
        std::vector<int> set_kind; // 0: {zeta,1/zeta}; 1: {alpha,zeta,1/zeta}; 2: {zeta}
        for (auto& chain : chains) {
            smalls.push_back(chain.g);
            set_kind.push_back(0);
            smalls.push_back(chain.h);
            set_kind.push_back(1);
            std::vector<fr> reversal(b);
            for (size_t k = 0; k < b; ++k) {
                reversal[k] = chain.g[b - 1 - k];
            }
            smalls.push_back(std::move(reversal));
            set_kind.push_back(2);
        }
        smalls.push_back(s_poly);
        set_kind.push_back(0);

        // m(X) = sum_k beta^k Z_{T\S_k}(X) (p_k(X) - p_k*(X)); q_m = m / Z_T via three exact
        // linear divisions.
        const size_t max_len = std::max(b, rows) + 2;
        std::vector<fr> m_poly(max_len + 1, fr::zero());
        fr beta_power = fr::one();
        for (size_t k = 0; k < smalls.size(); ++k) {
            const std::vector<fr>& p = smalls[k];
            std::vector<fr> term(p.size() + 2, fr::zero());
            const auto [xs, ys] = opening_points(set_kind[k], p, zeta, zeta_inv, alpha);
            // p(X) - p*(X): subtract the interpolant through the opening points.
            std::vector<fr> shifted_p(p.begin(), p.end());
            subtract_interpolant(shifted_p, xs, ys);
            multiply_by_complement(term, shifted_p, set_kind[k], zeta_inv, alpha);
            for (size_t i = 0; i < term.size(); ++i) {
                m_poly[i] += beta_power * term[i];
            }
            beta_power *= beta;
        }
        detail::divide_by_linear(m_poly, zeta);
        detail::divide_by_linear(m_poly, zeta_inv);
        detail::divide_by_linear(m_poly, alpha);
        transcript->send_to_verifier(std::string("MERCURY:q_m"), commit(ck, m_poly));

        // Round 6: linearization quotient q_L at z.
        const fr z = transcript->template get_challenge<fr>("MERCURY:z");
        const fr z_t_z = (z - zeta) * (z - zeta_inv) * (z - alpha);
        std::vector<fr> l_poly(max_len + 1, fr::zero());
        beta_power = fr::one();
        for (size_t k = 0; k < smalls.size(); ++k) {
            const std::vector<fr>& p = smalls[k];
            const auto [xs, ys] = opening_points(set_kind[k], p, zeta, zeta_inv, alpha);
            const fr complement_z = complement_at(set_kind[k], z, zeta_inv, alpha);
            const fr interpolant_z = detail::interpolate_at(xs, ys, z);
            const fr scale = beta_power * complement_z;
            for (size_t i = 0; i < p.size(); ++i) {
                l_poly[i] += scale * p[i];
            }
            l_poly[0] -= scale * interpolant_z;
            beta_power *= beta;
        }
        for (size_t i = 0; i < m_poly.size(); ++i) {
            l_poly[i] -= z_t_z * m_poly[i];
        }
        detail::divide_by_linear(l_poly, z);
        transcript->send_to_verifier(std::string("MERCURY:q_L"), commit(ck, l_poly));

        // Pairing-batch challenge: drawn on both sides to keep the transcripts aligned.
        [[maybe_unused]] const fr pairing_batch = transcript->template get_challenge<fr>("MERCURY:pairing_batch");
    }

  private:
    template <typename Container> static Commitment commit(const MercuryCommitmentKey& ck, const Container& poly)
    {
        return ck.kzg_ck->commit(Polynomial<fr>(std::span<const fr>(poly)));
    }

    static std::pair<std::vector<fr>, std::vector<fr>> opening_points(
        int set_kind, const std::vector<fr>& p, const fr& zeta, const fr& zeta_inv, const fr& alpha)
    {
        std::vector<fr> xs;
        if (set_kind == 0) {
            xs = { zeta, zeta_inv };
        } else if (set_kind == 1) {
            xs = { alpha, zeta, zeta_inv };
        } else {
            xs = { zeta };
        }
        std::vector<fr> ys;
        for (const fr& x : xs) {
            ys.push_back(polynomial_arithmetic::evaluate(p.data(), x, p.size()));
        }
        return { xs, ys };
    }

    /** @brief p -= interpolant(xs, ys) (degree < |xs|), in place on the low coefficients. */
    static void subtract_interpolant(std::vector<fr>& p, const std::vector<fr>& xs, const std::vector<fr>& ys)
    {
        // Lagrange to coefficients for up to three points.
        std::vector<fr> coeffs(xs.size(), fr::zero());
        for (size_t i = 0; i < xs.size(); ++i) {
            // scale * prod_{j != i} (X - xs[j])
            fr denom = fr::one();
            for (size_t j = 0; j < xs.size(); ++j) {
                if (j != i) {
                    denom *= (xs[i] - xs[j]);
                }
            }
            const fr scale = ys[i] * denom.invert();
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
                coeffs[d] += scale * basis[d];
            }
        }
        for (size_t d = 0; d < coeffs.size(); ++d) {
            p[d] -= coeffs[d];
        }
    }

    /** @brief term = Z_{T\S}(X) * p(X) for the set kind's complement in T = {zeta, 1/zeta, alpha}. */
    static void multiply_by_complement(
        std::vector<fr>& term, const std::vector<fr>& p, int set_kind, const fr& zeta_inv, const fr& alpha)
    {
        std::copy(p.begin(), p.end(), term.begin());
        if (set_kind == 0) {
            // Z = (X - alpha)
            multiply_by_linear(term, alpha);
        } else if (set_kind == 2) {
            // Z = (X - alpha)(X - 1/zeta)
            multiply_by_linear(term, alpha);
            multiply_by_linear(term, zeta_inv);
        } // set_kind 1: Z = 1
    }

    /** @brief In-place multiply by (X - root); the buffer has slack for the degree increase. */
    static void multiply_by_linear(std::vector<fr>& poly, const fr& root)
    {
        for (size_t d = poly.size(); d-- > 1;) {
            poly[d] = poly[d - 1] - root * poly[d];
        }
        poly[0] = -root * poly[0];
    }

    static fr complement_at(int set_kind, const fr& z, const fr& zeta_inv, const fr& alpha)
    {
        if (set_kind == 0) {
            return z - alpha;
        }
        if (set_kind == 2) {
            return (z - alpha) * (z - zeta_inv);
        }
        return fr::one();
    }

    friend class MercuryVerifier;
};

/** @brief Batched Mercury opening verifier; the mirror of `MercuryProver`. */
class MercuryVerifier {
  public:
    struct Claims {
        std::vector<size_t> group_num_columns;
        std::vector<std::vector<Commitment>> group_commitments; // when non-empty, transcript-bound
        std::vector<MercuryColumnRef> unshifted;
        std::vector<fr> unshifted_evaluations;
        std::vector<MercuryColumnRef> to_be_shifted;
        std::vector<fr> shifted_evaluations;
    };

    template <typename Transcript>
    static bool verify(const MercuryConfig& config,
                       const Claims& claims,
                       std::span<const fr> u,
                       const std::shared_ptr<Transcript>& transcript)
    {
        using GroupElement = Curve::Element;
        const size_t b = config.num_cols();
        BB_ASSERT_EQ(u.size(), config.num_variables, "opening point size mismatch");

        std::vector<std::vector<Commitment>> group_commitments = claims.group_commitments;
        if (group_commitments.empty()) {
            for (size_t g = 0; g < claims.group_num_columns.size(); ++g) {
                std::vector<Commitment> commitments;
                for (size_t c = 0; c < claims.group_num_columns[g]; ++c) {
                    commitments.push_back(transcript->template receive_from_prover<Commitment>(
                        detail::mercury_label("root", g) + "_" + std::to_string(c)));
                }
                group_commitments.push_back(std::move(commitments));
            }
        }
        const fr rho = transcript->template get_challenge<fr>("MERCURY:rho");

        // Homomorphic combined commitments and claimed evaluations per chain.
        struct ChainView {
            GroupElement commitment;
            fr claimed_evaluation = fr::zero();
            bool is_shifted = false;
            fr h_alpha, g_zeta, g_zeta_inv, h_zeta, h_zeta_inv, d_zeta;
            Commitment c_h, c_g, c_q, c_d, c_big_h;
        };
        std::vector<ChainView> chains;
        {
            ChainView chain_a;
            chain_a.commitment = GroupElement::infinity();
            fr rho_power = fr::one();
            for (size_t i = 0; i < claims.unshifted.size(); ++i) {
                chain_a.commitment +=
                    GroupElement(group_commitments[claims.unshifted[i].group][claims.unshifted[i].column]) * rho_power;
                chain_a.claimed_evaluation += rho_power * claims.unshifted_evaluations[i];
                rho_power *= rho;
            }
            chains.push_back(chain_a);
            if (!claims.to_be_shifted.empty()) {
                ChainView chain_b;
                chain_b.is_shifted = true;
                chain_b.commitment = GroupElement::infinity();
                for (size_t l = 0; l < claims.to_be_shifted.size(); ++l) {
                    chain_b.commitment +=
                        GroupElement(group_commitments[claims.to_be_shifted[l].group][claims.to_be_shifted[l].column]) *
                        rho_power;
                    chain_b.claimed_evaluation += rho_power * claims.shifted_evaluations[l];
                    rho_power *= rho;
                }
                chains.push_back(chain_b);
            }
        }

        for (size_t c = 0; c < chains.size(); ++c) {
            chains[c].c_h = transcript->template receive_from_prover<Commitment>(detail::mercury_label("h", c));
        }
        const fr alpha = transcript->template get_challenge<fr>("MERCURY:alpha");
        for (size_t c = 0; c < chains.size(); ++c) {
            chains[c].c_g = transcript->template receive_from_prover<Commitment>(detail::mercury_label("g", c));
            chains[c].c_q = transcript->template receive_from_prover<Commitment>(detail::mercury_label("q", c));
            chains[c].h_alpha = transcript->template receive_from_prover<fr>(detail::mercury_label("h_alpha", c));
        }
        const fr gamma = transcript->template get_challenge<fr>("MERCURY:gamma");
        const Commitment c_s = transcript->template receive_from_prover<Commitment>(std::string("MERCURY:S"));
        for (size_t c = 0; c < chains.size(); ++c) {
            chains[c].c_d = transcript->template receive_from_prover<Commitment>(detail::mercury_label("D", c));
        }
        const fr zeta = transcript->template get_challenge<fr>("MERCURY:zeta");
        const fr zeta_inv = zeta.invert();
        const fr zeta_pow_b = zeta.pow(uint256_t(uint64_t(b)));
        for (size_t c = 0; c < chains.size(); ++c) {
            ChainView& chain = chains[c];
            chain.g_zeta = transcript->template receive_from_prover<fr>(detail::mercury_label("g_zeta", c));
            chain.g_zeta_inv = transcript->template receive_from_prover<fr>(detail::mercury_label("g_zeta_inv", c));
            chain.h_zeta = transcript->template receive_from_prover<fr>(detail::mercury_label("h_zeta", c));
            chain.h_zeta_inv = transcript->template receive_from_prover<fr>(detail::mercury_label("h_zeta_inv", c));
            chain.d_zeta = transcript->template receive_from_prover<fr>(detail::mercury_label("D_zeta", c));
            chain.c_big_h = transcript->template receive_from_prover<Commitment>(detail::mercury_label("H", c));
        }
        const fr s_zeta = transcript->template receive_from_prover<fr>(std::string("MERCURY:S_zeta"));
        const fr s_zeta_inv = transcript->template receive_from_prover<fr>(std::string("MERCURY:S_zeta_inv"));
        const fr beta = transcript->template get_challenge<fr>("MERCURY:beta");
        const Commitment c_qm = transcript->template receive_from_prover<Commitment>(std::string("MERCURY:q_m"));
        const fr z = transcript->template get_challenge<fr>("MERCURY:z");
        const Commitment c_ql = transcript->template receive_from_prover<Commitment>(std::string("MERCURY:q_L"));
        const fr pairing_batch = transcript->template get_challenge<fr>("MERCURY:pairing_batch");

        // Check 1: batched inner-product identity at zeta.
        const fr p_lo_zeta = detail::tensor_poly_eval(u.subspan(0, config.num_col_vars), zeta);
        const fr p_lo_zeta_inv = detail::tensor_poly_eval(u.subspan(0, config.num_col_vars), zeta_inv);
        const fr p_hi_zeta = detail::tensor_poly_eval(u.subspan(config.num_col_vars), zeta);
        const fr p_hi_zeta_inv = detail::tensor_poly_eval(u.subspan(config.num_col_vars), zeta_inv);
        fr lhs = fr::zero();
        fr constant = fr::zero();
        fr gamma_power = fr::one();
        for (const ChainView& chain : chains) {
            lhs += gamma_power * (chain.g_zeta * p_lo_zeta_inv + chain.g_zeta_inv * p_lo_zeta);
            constant += gamma_power * chain.h_alpha;
            gamma_power *= gamma;
            lhs += gamma_power * (chain.h_zeta * p_hi_zeta_inv + chain.h_zeta_inv * p_hi_zeta);
            constant += gamma_power * chain.claimed_evaluation;
            gamma_power *= gamma;
        }
        if (lhs != constant + constant + zeta * s_zeta + zeta_inv * s_zeta_inv) {
            return false;
        }

        // Check 2: degree bounds deg(g) < b.
        const fr zeta_pow_b_minus_1 = zeta.pow(uint256_t(uint64_t(b - 1)));
        for (const ChainView& chain : chains) {
            if (chain.d_zeta != zeta_pow_b_minus_1 * chain.g_zeta_inv) {
                return false;
            }
        }

        // Check 3: fold-consistency pairings and the BDFG opening, batched into one pairing.
        GroupElement p0 = GroupElement::infinity();
        GroupElement p1 = GroupElement::infinity();
        fr batch_power = fr::one();
        const fr fold_factor = zeta_pow_b - alpha;
        for (const ChainView& chain : chains) {
            const fr remainder_term = chain.is_shifted ? zeta * chain.g_zeta : chain.g_zeta;
            p0 += (chain.commitment - GroupElement(chain.c_q) * fold_factor - GroupElement::one() * remainder_term +
                   GroupElement(chain.c_big_h) * zeta) *
                  batch_power;
            p1 -= GroupElement(chain.c_big_h) * batch_power;
            batch_power *= pairing_batch;
        }

        // BDFG: F = sum_k beta^k Z_{T\S_k}(z) (C_k - p_k*(z) [1]) - Z_T(z) C_qm; identity
        // e(F + z q_L, [1]) = e(q_L, [tau]).
        const fr z_t_z = (z - zeta) * (z - zeta_inv) * (z - alpha);
        GroupElement f_point = GroupElement::infinity();
        fr scalar_accumulator = fr::zero();
        fr beta_power = fr::one();
        for (const ChainView& chain : chains) {
            // g: set {zeta, 1/zeta}
            add_bdfg_term(f_point,
                          scalar_accumulator,
                          GroupElement(chain.c_g),
                          beta_power * (z - alpha),
                          detail::interpolate_at(
                              std::vector<fr>{ zeta, zeta_inv }, std::vector<fr>{ chain.g_zeta, chain.g_zeta_inv }, z));
            beta_power *= beta;
            // h: set {alpha, zeta, 1/zeta}
            add_bdfg_term(f_point,
                          scalar_accumulator,
                          GroupElement(chain.c_h),
                          beta_power,
                          detail::interpolate_at(std::vector<fr>{ alpha, zeta, zeta_inv },
                                                 std::vector<fr>{ chain.h_alpha, chain.h_zeta, chain.h_zeta_inv },
                                                 z));
            beta_power *= beta;
            // D: set {zeta}
            add_bdfg_term(f_point,
                          scalar_accumulator,
                          GroupElement(chain.c_d),
                          beta_power * (z - alpha) * (z - zeta_inv),
                          chain.d_zeta);
            beta_power *= beta;
        }
        add_bdfg_term(
            f_point,
            scalar_accumulator,
            GroupElement(c_s),
            beta_power * (z - alpha),
            detail::interpolate_at(std::vector<fr>{ zeta, zeta_inv }, std::vector<fr>{ s_zeta, s_zeta_inv }, z));
        f_point -= GroupElement::one() * scalar_accumulator;
        f_point -= GroupElement(c_qm) * z_t_z;

        p0 += (f_point + GroupElement(c_ql) * z) * batch_power;
        p1 -= GroupElement(c_ql) * batch_power;

        PairingPoints<Curve> pairing_points{ Commitment(p0), Commitment(p1) };
        return pairing_points.check();
    }

  private:
    using GroupElement = Curve::Element;

    /** @brief f_point += complement_z * (C - interpolant_z [1]); the [1] parts accumulate in scalar. */
    static void add_bdfg_term(GroupElement& f_point,
                              fr& scalar_accumulator,
                              const GroupElement& commitment,
                              const fr& complement_z,
                              const fr& interpolant_z)
    {
        f_point += commitment * complement_z;
        scalar_accumulator += complement_z * interpolant_z;
    }
};

} // namespace bb::mercury
