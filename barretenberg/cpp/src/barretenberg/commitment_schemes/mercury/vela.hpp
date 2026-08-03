#pragma once

#include "barretenberg/commitment_schemes/commitment_key.hpp"
#include "barretenberg/commitment_schemes/mercury/mercury.hpp"
#include "barretenberg/commitment_schemes/pairing_points.hpp"
#include "barretenberg/commitment_schemes/utils/batch_accumulate.hpp"
#include "barretenberg/common/assert.hpp"
#include "barretenberg/common/thread.hpp"
#include "barretenberg/ecc/curves/bn254/bn254.hpp"
#include "barretenberg/polynomials/polynomial.hpp"
#include "barretenberg/polynomials/polynomial_arithmetic.hpp"

#include <memory>
#include <span>
#include <string>
#include <vector>

namespace bb::vela {

using Curve = curve::BN254;
using Commitment = Curve::AffineElement;
using GroupElement = Curve::Element;

// Vela reuses Mercury's commitment layer verbatim: the commitment to a multilinear is the KZG
// commitment to the univariate twin f_v(X) = sum_i a_i X^i whose coefficients are the evaluation
// table. Only the opening argument differs.
using VelaColumnRef = mercury::MercuryColumnRef;
using VelaCommitmentKey = mercury::MercuryCommitmentKey;
using VelaGroupData = mercury::MercuryGroupData;

/** @brief Vela parameters. Only the payload size matters; the matrix split is Mercury-only. */
struct VelaConfig {
    size_t num_variables;

    static VelaConfig create(size_t num_variables, size_t /*security_bits*/ = 100, size_t /*log_inv_rate*/ = 0)
    {
        BB_ASSERT_GT(num_variables, size_t(1));
        return { num_variables };
    }

    operator mercury::MercuryConfig() const { return { num_variables, (num_variables + 1) / 2 }; }
};

namespace detail {

inline std::string vela_label(const std::string& name)
{
    return "VELA:" + name;
}

/**
 * @brief Expand the Laurent product H(X) = g(X)·T_r(X^{-1}) into its coefficient table.
 *
 * @details g is the combined twin, supplied as coefficients over exponents `[-shift_low, len-1-
 * shift_low]`; T_r(X) = prod_k ((1-r_k) + r_k X^{2^k}) has 2^mu coefficients but only mu factors,
 * so the product is mu shift-and-add passes rather than a convolution: O(N log N) field operations
 * and no FFT (paper §2.1). The returned table is indexed by `exponent + offset` where
 * `offset = shift_low + N - 1`, covering exponents [-offset, N-1-shift_low].
 *
 * Each pass writes into a fresh buffer rather than updating in place: `table[e]` reads `table[e+2^k]`,
 * so a parallel in-place pass would race (one thread's writes overlap another's reads).
 */
inline std::vector<fr> laurent_product(std::span<const fr> g, size_t shift_low, std::span<const fr> r)
{
    const size_t n = size_t(1) << r.size();
    const size_t offset = shift_low + n - 1;
    // Exponents run from -(offset) to (g.size() - 1 - shift_low).
    const size_t table_size = offset + g.size() - shift_low;
    std::vector<fr> table(table_size, fr::zero());
    for (size_t i = 0; i < g.size(); ++i) {
        table[offset - shift_low + i] = g[i];
    }
    std::vector<fr> next(table_size);
    // Multiply by (1 - r_k) + r_k X^{-2^k}: new[e] = (1-r_k) cur[e] + r_k cur[e + 2^k].
    for (size_t k = 0; k < r.size(); ++k) {
        const size_t step = size_t(1) << k;
        const fr one_minus = fr::one() - r[k];
        parallel_for_range(table_size, [&](size_t start, size_t end) {
            for (size_t e = start; e < end; ++e) {
                const fr high = (e + step < table_size) ? table[e + step] : fr::zero();
                next[e] = one_minus * table[e] + r[k] * high;
            }
        });
        table.swap(next);
    }
    return table;
}

/** @brief In-place exact division by (X - root). */
inline void divide_by_linear(std::vector<fr>& poly, const fr& root)
{
    mercury::detail::divide_by_linear(poly, root);
}

} // namespace detail

/**
 * @brief Vela opening prover (eprint 2026/1438 §2.2).
 *
 * @details The multilinear claim f(r) = y is the constant coefficient of the Laurent product
 * H(X) = f_v(X)·T_r(X^{-1}). Its inversion-symmetric residual D(X) = H(X) + H(X^{-1}) - 2y has
 * zero constant coefficient exactly when the claim holds, and symmetry then yields a single
 * ordinary polynomial h(X) of degree <= N-2 with D(X) = X h(X) + X^{-1} h(X^{-1}). One commitment
 * to h plus evaluations at z and 1/z therefore certifies the claim: three of the four evaluations
 * are transmitted and the fourth, h(1/z), is recovered by the verifier from the Laurent identity.
 *
 * The batched Honk claim set folds into one identity. With f_A the rho-combination of the
 * unshifted claims and f_B that of the to-be-shifted ones, the shifted twin is X^{-1} f_B(X)
 * (the shift contract's zero constant term makes the division exact), so the single combined twin
 * g(X) = f_A(X) + lambda X^{-1} f_B(X) carries both chains and only one h is committed.
 *
 * @note Deviation from the paper, forced by bb's verifier SRS: the paper opens the alpha-batched
 * polynomial at {z, 1/z} with one quotient and checks e(C, [1]_2) = e(pi, [Z(tau)]_2), which needs
 * [tau^2]_2. bb publishes only [1]_2 and [tau]_2, so the degree-2 vanishing check is linearized at
 * a further challenge, adding one G1 element to the proof (3 G1 + 5 F here against the paper's
 * 2 G1 + 3 F for a single unbatched claim).
 */
class VelaProver {
  public:
    struct Claims {
        std::vector<const VelaGroupData*> groups;
        std::vector<VelaColumnRef> unshifted;
        std::vector<fr> unshifted_evaluations;
        std::vector<VelaColumnRef> to_be_shifted;
        std::vector<fr> shifted_evaluations;
        bool send_roots = true;
    };

    template <typename Transcript>
    static void prove(const VelaCommitmentKey& ck,
                      const Claims& claims,
                      std::span<const fr> r,
                      const std::shared_ptr<Transcript>& transcript)
    {
        const size_t num_variables = r.size();
        const size_t n = size_t(1) << num_variables;
        const bool has_shifted = !claims.to_be_shifted.empty();

        if (claims.send_roots) {
            for (size_t g = 0; g < claims.groups.size(); ++g) {
                for (size_t c = 0; c < claims.groups[g]->num_columns(); ++c) {
                    transcript->send_to_verifier(detail::vela_label("root_") + std::to_string(g) + "_" +
                                                     std::to_string(c),
                                                 claims.groups[g]->commitments[c]);
                }
            }
        }
        const fr rho = transcript->template get_challenge<fr>("VELA:rho");
        const fr lambda = transcript->template get_challenge<fr>("VELA:lambda");

        // Combined twins and the combined claimed value.
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
            BB_ASSERT_EQ(column[0], fr::zero(), "to-be-shifted polynomial must have zero constant term");
            terms_b.push_back({ column.data(), rho_power, false });
            y_combined += lambda * rho_power * claims.shifted_evaluations[l];
            rho_power *= rho;
        }
        pcs_utils::accumulate_scaled(f_a, terms_a);
        pcs_utils::accumulate_scaled(f_b, terms_b);

        // g(X) = f_A(X) + lambda X^{-1} f_B(X), stored over exponents [-1, N-1].
        std::vector<fr> g(n + 1, fr::zero());
        for (size_t k = 0; k < n; ++k) {
            g[k + 1] += f_a[k];
        }
        if (has_shifted) {
            for (size_t k = 0; k < n; ++k) {
                g[k] += lambda * f_b[k];
            }
        }
        const std::vector<fr> table = detail::laurent_product(g, /*shift_low=*/1, r);
        const size_t offset = n; // exponent e sits at table[e + n]

        // h(X) = sum_{k>=1} d_k X^{k-1} with d_k = H_k + H_{-k}; the k = 0 coefficient is
        // 2(claim - y) and is not committed.
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
        transcript->send_to_verifier(detail::vela_label("C_h"), commit(ck, h));

        const fr z = transcript->template get_challenge<fr>("VELA:z");
        BB_ASSERT(z != fr::zero() && z != fr::one() && z != -fr::one(), "degenerate Vela challenge");
        const fr z_inv = z.invert();

        const fr v0_a = polynomial_arithmetic::evaluate(f_a.data(), z, n);
        const fr v1_a = polynomial_arithmetic::evaluate(f_a.data(), z_inv, n);
        const fr w0 = polynomial_arithmetic::evaluate(h.data(), z, n);
        transcript->send_to_verifier(detail::vela_label("v0_a"), v0_a);
        transcript->send_to_verifier(detail::vela_label("v1_a"), v1_a);
        fr v0_b = fr::zero();
        fr v1_b = fr::zero();
        if (has_shifted) {
            v0_b = polynomial_arithmetic::evaluate(f_b.data(), z, n);
            v1_b = polynomial_arithmetic::evaluate(f_b.data(), z_inv, n);
            transcript->send_to_verifier(detail::vela_label("v0_b"), v0_b);
            transcript->send_to_verifier(detail::vela_label("v1_b"), v1_b);
        }
        transcript->send_to_verifier(detail::vela_label("w0"), w0);

        // The verifier recovers w1 = h(1/z) from the Laurent identity; the prover recomputes it the
        // same way so that both sides batch the identical value.
        const fr t_z = mercury::detail::tensor_poly_eval(r, z);
        const fr t_z_inv = mercury::detail::tensor_poly_eval(r, z_inv);
        const fr w1 = recover_w1(z, z_inv, v0_a, v1_a, v0_b, v1_b, w0, y_combined, t_z, t_z_inv, lambda, has_shifted);

        // alpha-batch {f_A, f_B, h} and open the combination at {z, 1/z}.
        const fr alpha = transcript->template get_challenge<fr>("VELA:alpha");
        std::vector<fr> batched(n, fr::zero());
        {
            std::vector<pcs_utils::ScaledTerm> batch_terms;
            fr alpha_power = fr::one();
            batch_terms.push_back({ f_a.data(), alpha_power, false });
            alpha_power *= alpha;
            if (has_shifted) {
                batch_terms.push_back({ f_b.data(), alpha_power, false });
            }
            alpha_power *= alpha;
            batch_terms.push_back({ h.data(), alpha_power, false });
            pcs_utils::accumulate_scaled(batched, batch_terms);
        }

        const auto [r0, r1] = interpolant(z,
                                          z_inv,
                                          batched_value(alpha, v0_a, v0_b, w0, has_shifted),
                                          batched_value(alpha, v1_a, v1_b, w1, has_shifted));
        std::vector<fr> numerator = batched;
        numerator[0] -= r0;
        numerator[1] -= r1;
        detail::divide_by_linear(numerator, z);
        detail::divide_by_linear(numerator, z_inv);
        transcript->send_to_verifier(detail::vela_label("C_q"), commit(ck, numerator));

        // Linearization: L(X) = batched(X) - Z(zeta) q(X) satisfies L(zeta) = R(zeta).
        const fr zeta = transcript->template get_challenge<fr>("VELA:zeta");
        const fr z_at_zeta = (zeta - z) * (zeta - z_inv);
        std::vector<fr> linear = batched;
        for (size_t k = 0; k < numerator.size(); ++k) {
            linear[k] -= z_at_zeta * numerator[k];
        }
        linear[0] -= r0 + r1 * zeta;
        detail::divide_by_linear(linear, zeta);
        transcript->send_to_verifier(detail::vela_label("pi_L"), commit(ck, linear));
    }

    /** @brief w1 = z (g(z) T_r(1/z) + g(1/z) T_r(z) - 2Y - z w0), the paper's recovery formula. */
    static fr recover_w1(const fr& z,
                         const fr& z_inv,
                         const fr& v0_a,
                         const fr& v1_a,
                         const fr& v0_b,
                         const fr& v1_b,
                         const fr& w0,
                         const fr& y_combined,
                         const fr& t_z,
                         const fr& t_z_inv,
                         const fr& lambda,
                         bool has_shifted)
    {
        fr g_z = v0_a;
        fr g_z_inv = v1_a;
        if (has_shifted) {
            g_z += lambda * z_inv * v0_b;
            g_z_inv += lambda * z * v1_b;
        }
        return z * (g_z * t_z_inv + g_z_inv * t_z - (y_combined + y_combined) - z * w0);
    }

    /** @brief The alpha-batched value of {f_A, f_B, h} at one point. */
    static fr batched_value(const fr& alpha, const fr& a, const fr& b, const fr& h_value, bool has_shifted)
    {
        fr result = a;
        fr power = alpha;
        if (has_shifted) {
            result += power * b;
        }
        power *= alpha;
        result += power * h_value;
        return result;
    }

    /** @brief Coefficients (r0, r1) of the line through (z, y_z) and (1/z, y_z_inv). */
    static std::pair<fr, fr> interpolant(const fr& z, const fr& z_inv, const fr& y_z, const fr& y_z_inv)
    {
        const fr slope = (y_z - y_z_inv) * (z - z_inv).invert();
        return { y_z - slope * z, slope };
    }

  private:
    template <typename Container> static Commitment commit(const VelaCommitmentKey& ck, const Container& poly)
    {
        return ck.kzg_ck->commit(Polynomial<fr>(std::span<const fr>(poly)));
    }
};

/** @brief Batched Vela opening verifier: one recovery, one linear combination, one pairing check. */
class VelaVerifier {
  public:
    struct Claims {
        std::vector<size_t> group_num_columns;
        std::vector<std::vector<Commitment>> group_commitments; // when non-empty, transcript-bound
        std::vector<VelaColumnRef> unshifted;
        std::vector<fr> unshifted_evaluations;
        std::vector<VelaColumnRef> to_be_shifted;
        std::vector<fr> shifted_evaluations;
    };

    template <typename Transcript>
    static bool verify(const VelaConfig& /*config*/,
                       const Claims& claims,
                       std::span<const fr> r,
                       const std::shared_ptr<Transcript>& transcript)
    {
        const bool has_shifted = !claims.to_be_shifted.empty();

        std::vector<std::vector<Commitment>> group_commitments = claims.group_commitments;
        if (group_commitments.empty()) {
            for (size_t g = 0; g < claims.group_num_columns.size(); ++g) {
                std::vector<Commitment> commitments;
                for (size_t c = 0; c < claims.group_num_columns[g]; ++c) {
                    commitments.push_back(transcript->template receive_from_prover<Commitment>(
                        detail::vela_label("root_") + std::to_string(g) + "_" + std::to_string(c)));
                }
                group_commitments.push_back(std::move(commitments));
            }
        }
        const fr rho = transcript->template get_challenge<fr>("VELA:rho");
        const fr lambda = transcript->template get_challenge<fr>("VELA:lambda");

        GroupElement c_a = GroupElement::infinity();
        GroupElement c_b = GroupElement::infinity();
        fr y_combined = fr::zero();
        fr rho_power = fr::one();
        for (size_t i = 0; i < claims.unshifted.size(); ++i) {
            c_a += GroupElement(group_commitments[claims.unshifted[i].group][claims.unshifted[i].column]) * rho_power;
            y_combined += rho_power * claims.unshifted_evaluations[i];
            rho_power *= rho;
        }
        for (size_t l = 0; l < claims.to_be_shifted.size(); ++l) {
            c_b += GroupElement(group_commitments[claims.to_be_shifted[l].group][claims.to_be_shifted[l].column]) *
                   rho_power;
            y_combined += lambda * rho_power * claims.shifted_evaluations[l];
            rho_power *= rho;
        }

        const Commitment c_h = transcript->template receive_from_prover<Commitment>(detail::vela_label("C_h"));
        const fr z = transcript->template get_challenge<fr>("VELA:z");
        if (z == fr::zero() || z == fr::one() || z == -fr::one()) {
            return false;
        }
        const fr z_inv = z.invert();

        const fr v0_a = transcript->template receive_from_prover<fr>(detail::vela_label("v0_a"));
        const fr v1_a = transcript->template receive_from_prover<fr>(detail::vela_label("v1_a"));
        fr v0_b = fr::zero();
        fr v1_b = fr::zero();
        if (has_shifted) {
            v0_b = transcript->template receive_from_prover<fr>(detail::vela_label("v0_b"));
            v1_b = transcript->template receive_from_prover<fr>(detail::vela_label("v1_b"));
        }
        const fr w0 = transcript->template receive_from_prover<fr>(detail::vela_label("w0"));

        const fr t_z = mercury::detail::tensor_poly_eval(r, z);
        const fr t_z_inv = mercury::detail::tensor_poly_eval(r, z_inv);
        const fr w1 =
            VelaProver::recover_w1(z, z_inv, v0_a, v1_a, v0_b, v1_b, w0, y_combined, t_z, t_z_inv, lambda, has_shifted);

        const fr alpha = transcript->template get_challenge<fr>("VELA:alpha");
        const Commitment c_q = transcript->template receive_from_prover<Commitment>(detail::vela_label("C_q"));
        const fr zeta = transcript->template get_challenge<fr>("VELA:zeta");
        const Commitment pi_l = transcript->template receive_from_prover<Commitment>(detail::vela_label("pi_L"));

        // Batched commitment and interpolant, mirroring the prover's alpha layout.
        GroupElement c_batched = c_a;
        fr alpha_power = alpha;
        if (has_shifted) {
            c_batched += c_b * alpha_power;
        }
        alpha_power *= alpha;
        c_batched += GroupElement(c_h) * alpha_power;

        const auto [r0, r1] = VelaProver::interpolant(z,
                                                      z_inv,
                                                      VelaProver::batched_value(alpha, v0_a, v0_b, w0, has_shifted),
                                                      VelaProver::batched_value(alpha, v1_a, v1_b, w1, has_shifted));

        // e(C_L - R(zeta)[1] + zeta pi_L, [1]_2) = e(pi_L, [tau]_2).
        const fr z_at_zeta = (zeta - z) * (zeta - z_inv);
        const fr r_at_zeta = r0 + r1 * zeta;
        GroupElement p0 = c_batched - GroupElement(c_q) * z_at_zeta - GroupElement(Commitment::one()) * r_at_zeta +
                          GroupElement(pi_l) * zeta;
        GroupElement p1 = -GroupElement(pi_l);

        PairingPoints<Curve> pairing_points{ Commitment(p0), Commitment(p1) };
        return pairing_points.check();
    }
};

} // namespace bb::vela
