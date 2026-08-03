#pragma once

#include "barretenberg/commitment_schemes/utils/batch_accumulate.hpp"
#include "barretenberg/commitment_schemes/whir/weights.hpp"
#include "barretenberg/common/assert.hpp"
#include "barretenberg/common/thread.hpp"
#include "barretenberg/ecc/curves/bn254/bn254.hpp"
#include "barretenberg/ecc/scalar_multiplication/scalar_multiplication.hpp"
#include "barretenberg/polynomials/polynomial.hpp"

#include <memory>
#include <span>
#include <string>
#include <vector>

namespace bb::pedersen_ipa {

using Curve = curve::BN254;
using Commitment = Curve::AffineElement;
using GroupElement = Curve::Element;

/** @brief Reference to one polynomial of one committed group (duck-type shared with the suite). */
struct IpaColumnRef {
    size_t group;
    size_t column;
};

/** @brief Parameters: only the size; the scheme is transparent and unstructured. */
struct IpaConfig {
    size_t num_variables;

    static IpaConfig create(size_t num_variables, size_t /*security_bits*/ = 100, size_t /*log_inv_rate*/ = 0)
    {
        return { num_variables };
    }
};

/** @brief Prover-side commitment to a group: dense arrays and one Pedersen commitment per column. */
struct IpaGroupData {
    std::vector<std::vector<fr>> coefficients;
    std::vector<Commitment> commitments;

    size_t num_columns() const { return coefficients.size(); }
};

/**
 * @brief Transparent commitment key: hash-derived independent generators (the powers-of-tau SRS
 * cannot serve as Pedersen/IPA generators — its points are DL-related).
 */
class IpaCommitmentKey {
  public:
    explicit IpaCommitmentKey(const IpaConfig& config)
        : config(config)
        , generators(Curve::Group::derive_generators(std::vector<uint8_t>{ 'b', 'b', '_', 'i', 'p', 'a' },
                                                     size_t(1) << config.num_variables))
        , auxiliary_generator(
              Curve::Group::derive_generators(std::vector<uint8_t>{ 'b', 'b', '_', 'i', 'p', 'a', '_', 'u' }, 1)[0])
    {}

    Commitment commit(std::span<const fr> coefficients) const
    {
        std::vector<fr> scalars(coefficients.begin(), coefficients.end());
        scalars.resize(size_t(1) << config.num_variables, fr::zero());
        return Commitment::batch_mul(std::span<const Commitment>(generators), std::span<fr>(scalars));
    }

    IpaGroupData commit_group(std::vector<std::vector<fr>> payload_columns,
                              const std::vector<bool>& /*to_be_shifted*/ = {}) const
    {
        const size_t n = size_t(1) << config.num_variables;
        IpaGroupData data;
        for (auto& column : payload_columns) {
            BB_ASSERT_LTE(column.size(), n, "polynomial too large for the configured size");
            data.commitments.push_back(commit(column));
            column.resize(n, fr::zero());
            data.coefficients.push_back(std::move(column));
        }
        return data;
    }

    IpaGroupData commit_group(std::span<const Polynomial<fr>* const> polynomials,
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

    IpaConfig config;
    std::vector<Commitment> generators;
    Commitment auxiliary_generator;
};

namespace detail {

inline std::string ipa_label(const std::string& name, size_t i, size_t j)
{
    return "IPA:" + name + "_" + std::to_string(i) + "_" + std::to_string(j);
}

/**
 * @brief Per-variable fold multipliers from the round challenge inverses.
 * @details Round j splits the current vectors into contiguous lo/hi halves, so it binds the highest
 * remaining variable: variable j of the index is folded by the challenge of round `ell-1-j`.
 */
inline std::vector<fr> fold_multipliers(std::span<const fr> challenge_inverses)
{
    return { challenge_inverses.rbegin(), challenge_inverses.rend() };
}

/**
 * @brief Final b-scalar of the eq(u) vector after the IPA folds:
 * <eq(u), tensor(1, m_j)> = prod_j ((1-u_j) + u_j m_j), m_j the multiplier of variable j.
 */
inline fr eq_fold_scalar(std::span<const fr> u, std::span<const fr> multipliers)
{
    fr result = fr::one();
    for (size_t j = 0; j < u.size(); ++j) {
        result *= (fr::one() - u[j]) + u[j] * multipliers[j];
    }
    return result;
}

/**
 * @brief Final b-scalar of the shifted-eq vector (b_0 = 0, b_i = eq_{i-1}(u)): the successor
 * kernel's telescoping closed form,
 * sum_k [prod_{j<k} u_j] (1-u_k) m_k [prod_{j>k} ((1-u_j) + u_j m_j)].
 */
inline fr shifted_eq_fold_scalar(std::span<const fr> u, std::span<const fr> multipliers)
{
    const size_t ell = u.size();
    // Suffix products of ((1-u_j) + u_j m_j)
    std::vector<fr> suffix(ell + 1, fr::one());
    for (size_t j = ell; j-- > 0;) {
        suffix[j] = suffix[j + 1] * ((fr::one() - u[j]) + u[j] * multipliers[j]);
    }
    fr result = fr::zero();
    fr prefix = fr::one(); // prod_{j<k} u_j
    for (size_t k = 0; k < ell; ++k) {
        result += prefix * (fr::one() - u[k]) * multipliers[k] * suffix[k + 1];
        prefix *= u[k];
    }
    return result;
}

} // namespace detail

/**
 * @brief Multilinear Bulletproofs inner-product argument over hash-derived BN254 generators
 * (non-hiding): the {transparent, small-proof, linear-verifier} corner. Two chains share the
 * transcript: the unshifted rho-combination against b = eq(u), and the to-be-shifted
 * rho-combination against b = shifted-eq(u), whose folded scalar the verifier computes with the
 * successor-kernel closed form — no shifted commitments and no univariate reduction needed.
 *
 * Per chain, log n rounds: L_j = <a_hi, G_lo> + x_U <a_hi, b_lo> U and symmetrically R_j;
 * challenge x_j folds a' = a_lo + x_j a_hi, b' = b_lo + x_j^{-1} b_hi,
 * G' = G_lo + x_j^{-1} G_hi. The verifier recomputes the folded generator with one size-n MSM
 * (the linear-verifier trade) and checks the final one-point identity.
 *
 * The lo/hi (rather than even/odd) split is what makes the round cost an MSM: both operands of each
 * cross term are contiguous, so `pippenger_unsafe` applies directly. Round j binds the *highest*
 * remaining variable, so the verifier's per-variable fold multipliers are the round challenges in
 * reverse order — see `fold_multipliers`.
 */
class IpaProver {
  public:
    struct Claims {
        std::vector<const IpaGroupData*> groups;
        std::vector<IpaColumnRef> unshifted;
        std::vector<fr> unshifted_evaluations;
        std::vector<IpaColumnRef> to_be_shifted;
        std::vector<fr> shifted_evaluations;
        bool send_roots = true;
    };

    template <typename Transcript>
    static void prove(const IpaCommitmentKey& ck,
                      const Claims& claims,
                      std::span<const fr> u,
                      const std::shared_ptr<Transcript>& transcript)
    {
        const size_t n = size_t(1) << ck.config.num_variables;
        BB_ASSERT_EQ(u.size(), ck.config.num_variables, "opening point size mismatch");

        if (claims.send_roots) {
            for (size_t g = 0; g < claims.groups.size(); ++g) {
                for (size_t c = 0; c < claims.groups[g]->num_columns(); ++c) {
                    transcript->send_to_verifier(detail::ipa_label("root", g, c), claims.groups[g]->commitments[c]);
                }
            }
        }
        const fr rho = transcript->template get_challenge<fr>("IPA:rho");

        // Chain arrays: unshifted combination and (unshifted!) to-be-shifted combination — the
        // shift lives entirely in chain 1's b-vector.
        std::vector<std::vector<fr>> chains;
        chains.emplace_back(n, fr::zero());
        std::vector<pcs_utils::ScaledTerm> terms_a;
        std::vector<pcs_utils::ScaledTerm> terms_b;
        fr rho_power = fr::one();
        for (size_t i = 0; i < claims.unshifted.size(); ++i) {
            const auto& column = claims.groups[claims.unshifted[i].group]->coefficients[claims.unshifted[i].column];
            terms_a.push_back({ column.data(), rho_power, false });
            rho_power *= rho;
        }
        if (!claims.to_be_shifted.empty()) {
            chains.emplace_back(n, fr::zero());
            for (size_t l = 0; l < claims.to_be_shifted.size(); ++l) {
                const auto& column =
                    claims.groups[claims.to_be_shifted[l].group]->coefficients[claims.to_be_shifted[l].column];
                terms_b.push_back({ column.data(), rho_power, false });
                rho_power *= rho;
            }
        }
        pcs_utils::accumulate_scaled(chains[0], terms_a);
        if (chains.size() > 1) {
            pcs_utils::accumulate_scaled(chains[1], terms_b);
        }

        const fr x_u = transcript->template get_challenge<fr>("IPA:x_u");
        for (size_t chain = 0; chain < chains.size(); ++chain) {
            std::vector<fr> a = std::move(chains[chain]);
            std::vector<fr> b = (chain == 0) ? whir::eq_tensor(u) : shifted_eq_vector(u);
            std::vector<Commitment> generators(ck.generators.begin(), ck.generators.end());
            const GroupElement aux = GroupElement(ck.auxiliary_generator) * x_u;
            std::vector<GroupElement> folded;

            for (size_t round = 0; a.size() > 1; ++round) {
                const size_t half = a.size() / 2;
                // L = <a_hi, G_lo> + <a_hi, b_lo> aux ; R = <a_lo, G_hi> + <a_lo, b_hi> aux.
                // Both operands of each cross term are contiguous, so each is one Pippenger MSM.
                fr ip_left = fr::zero();
                fr ip_right = fr::zero();
                for (size_t t = 0; t < half; ++t) {
                    ip_left += a[half + t] * b[t];
                    ip_right += a[t] * b[half + t];
                }
                GroupElement left = scalar_multiplication::pippenger_unsafe<Curve>(
                    PolynomialSpan<const fr>(0, { a.data() + half, half }),
                    std::span<const Commitment>(generators.data(), half));
                GroupElement right = scalar_multiplication::pippenger_unsafe<Curve>(
                    PolynomialSpan<const fr>(0, { a.data(), half }),
                    std::span<const Commitment>(generators.data() + half, half));
                left += aux * ip_left;
                right += aux * ip_right;
                transcript->send_to_verifier(detail::ipa_label("L", chain, round), Commitment(left));
                transcript->send_to_verifier(detail::ipa_label("R", chain, round), Commitment(right));
                const fr x = transcript->template get_challenge<fr>(detail::ipa_label("x", chain, round));
                const fr x_inv = x.invert();
                for (size_t t = 0; t < half; ++t) {
                    a[t] = a[t] + x * a[half + t];
                    b[t] = b[t] + x_inv * b[half + t];
                }
                // G'[t] = G_lo[t] + x^{-1} G_hi[t], batch-normalized back to affine so the next
                // round's MSMs can consume it directly.
                folded.resize(half);
                parallel_for_range(half, [&](size_t start, size_t end) {
                    for (size_t t = start; t < end; ++t) {
                        folded[t] = GroupElement(generators[half + t]) * x_inv + generators[t];
                    }
                });
                GroupElement::batch_normalize(folded.data(), half);
                for (size_t t = 0; t < half; ++t) {
                    generators[t] = Commitment(folded[t].x, folded[t].y);
                }
                a.resize(half);
                b.resize(half);
                generators.resize(half);
            }
            transcript->send_to_verifier(detail::ipa_label("a_final", chain, 0), a[0]);
        }
    }

    static std::vector<fr> shifted_eq_vector(std::span<const fr> u)
    {
        const std::vector<fr> eq = whir::eq_tensor(u);
        std::vector<fr> shifted(eq.size(), fr::zero());
        for (size_t i = 1; i < eq.size(); ++i) {
            shifted[i] = eq[i - 1];
        }
        return shifted;
    }
};

/** @brief The mirror of `IpaProver`; the folded-generator MSM makes verification O(n). */
class IpaVerifier {
  public:
    struct Claims {
        std::vector<size_t> group_num_columns;
        std::vector<std::vector<Commitment>> group_commitments; // when non-empty, transcript-bound
        std::vector<IpaColumnRef> unshifted;
        std::vector<fr> unshifted_evaluations;
        std::vector<IpaColumnRef> to_be_shifted;
        std::vector<fr> shifted_evaluations;
    };

    template <typename Transcript>
    static bool verify(const IpaConfig& config,
                       const Claims& claims,
                       std::span<const fr> u,
                       const std::shared_ptr<Transcript>& transcript,
                       const IpaCommitmentKey& ck)
    {
        const size_t ell = config.num_variables;
        BB_ASSERT_EQ(u.size(), ell, "opening point size mismatch");

        std::vector<std::vector<Commitment>> group_commitments = claims.group_commitments;
        if (group_commitments.empty()) {
            for (size_t g = 0; g < claims.group_num_columns.size(); ++g) {
                std::vector<Commitment> commitments;
                for (size_t c = 0; c < claims.group_num_columns[g]; ++c) {
                    commitments.push_back(
                        transcript->template receive_from_prover<Commitment>(detail::ipa_label("root", g, c)));
                }
                group_commitments.push_back(std::move(commitments));
            }
        }
        const fr rho = transcript->template get_challenge<fr>("IPA:rho");

        struct ChainView {
            GroupElement commitment = GroupElement::infinity();
            fr claimed_evaluation = fr::zero();
            bool is_shifted = false;
        };
        std::vector<ChainView> chains(claims.to_be_shifted.empty() ? 1 : 2);
        chains[0].is_shifted = false;
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

        const fr x_u = transcript->template get_challenge<fr>("IPA:x_u");
        const GroupElement aux = GroupElement(ck.auxiliary_generator) * x_u;

        for (size_t chain = 0; chain < chains.size(); ++chain) {
            GroupElement accumulator = chains[chain].commitment + aux * chains[chain].claimed_evaluation;
            std::vector<fr> challenges(ell);
            std::vector<fr> challenge_inverses(ell);
            for (size_t round = 0; round < ell; ++round) {
                const Commitment left =
                    transcript->template receive_from_prover<Commitment>(detail::ipa_label("L", chain, round));
                const Commitment right =
                    transcript->template receive_from_prover<Commitment>(detail::ipa_label("R", chain, round));
                const fr x = transcript->template get_challenge<fr>(detail::ipa_label("x", chain, round));
                challenges[round] = x;
                challenge_inverses[round] = x.invert();
                accumulator += GroupElement(left) * challenges[round] + GroupElement(right) * challenge_inverses[round];
            }
            const fr a_final = transcript->template receive_from_prover<fr>(detail::ipa_label("a_final", chain, 0));

            // Folded generator: G_final = sum_i (prod_{j : bit_j(i)=1} m_j) G_i — one size-n MSM,
            // with m_j the multiplier of variable j under the lo/hi round schedule.
            const std::vector<fr> multipliers = detail::fold_multipliers(challenge_inverses);
            const size_t n = size_t(1) << ell;
            std::vector<fr> scalars(n, fr::one());
            for (size_t j = 0; j < ell; ++j) {
                const size_t bit = size_t(1) << j;
                for (size_t i = 0; i < n; ++i) {
                    if ((i & bit) != 0) {
                        scalars[i] *= multipliers[j];
                    }
                }
            }
            const Commitment folded_generator =
                Commitment::batch_mul(std::span<const Commitment>(ck.generators), std::span<fr>(scalars));

            const fr b_final = chains[chain].is_shifted ? detail::shifted_eq_fold_scalar(u, multipliers)
                                                        : detail::eq_fold_scalar(u, multipliers);
            const GroupElement expected = (GroupElement(folded_generator) + aux * b_final) * a_final;
            if (Commitment(accumulator) != Commitment(expected)) {
                return false;
            }
        }
        return true;
    }
};

} // namespace bb::pedersen_ipa
