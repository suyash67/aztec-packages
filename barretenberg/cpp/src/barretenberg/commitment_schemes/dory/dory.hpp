#pragma once

#include "barretenberg/commitment_schemes/utils/batch_accumulate.hpp"
#include "barretenberg/commitment_schemes/whir/weights.hpp"
#include "barretenberg/common/assert.hpp"
#include "barretenberg/ecc/curves/bn254/bn254.hpp"
#include "barretenberg/ecc/curves/bn254/fq12.hpp"
#include "barretenberg/ecc/curves/bn254/pairing.hpp"
#include "barretenberg/polynomials/polynomial.hpp"

#include <array>
#include <cstring>
#include <map>
#include <memory>
#include <span>
#include <string>
#include <vector>

namespace bb::dory {

using Curve = curve::BN254;
using Commitment = Curve::AffineElement;
using GroupElement = Curve::Element;
using G2Affine = Curve::G2AffineElement;
using G2Element = bb::g2::element;

/** @brief Reference to one polynomial of one committed group (duck-type shared with the suite). */
struct DoryColumnRef {
    size_t group;
    size_t column;
};

/** @brief Dory parameters: the matrix split (columns = low variables). */
struct DoryConfig {
    size_t num_variables;
    size_t log_num_cols;

    size_t num_rows() const { return size_t(1) << (num_variables - log_num_cols); }
    size_t num_cols() const { return size_t(1) << log_num_cols; }
    size_t log_num_rows() const { return num_variables - log_num_cols; }

    static DoryConfig create(size_t num_variables, size_t /*security_bits*/ = 100, size_t /*log_inv_rate*/ = 0)
    {
        BB_ASSERT_GT(num_variables, size_t(1));
        return { num_variables, (num_variables + 1) / 2 };
    }
};

/**
 * @brief Prover-side commitment to a group: dense arrays, per-column AFGHO commitment T (one GT
 * element), and the per-column Pedersen row commitments P_i that T aggregates.
 */
struct DoryGroupData {
    std::vector<std::vector<fr>> coefficients;
    std::vector<fq12> commitments;                        // T_p = prod_i e(P_{p,i}, Gamma2_i)
    std::vector<std::vector<Commitment>> row_commitments; // [column][row]

    size_t num_columns() const { return coefficients.size(); }
};

namespace detail {

inline std::string dory_label(const std::string& name, size_t i, size_t j)
{
    return "DORY:" + name + "_" + std::to_string(i) + "_" + std::to_string(j);
}

/** @brief fq12 <-> 24 fr limbs (128-bit halves of each of the 12 fq coordinates). */
inline std::array<fr, 24> fq12_to_fields(const fq12& value)
{
    std::array<fr, 24> fields;
    const std::array<const fq*, 12> coords = {
        &value.c0.c0.c0, &value.c0.c0.c1, &value.c0.c1.c0, &value.c0.c1.c1, &value.c0.c2.c0, &value.c0.c2.c1,
        &value.c1.c0.c0, &value.c1.c0.c1, &value.c1.c1.c0, &value.c1.c1.c1, &value.c1.c2.c0, &value.c1.c2.c1,
    };
    for (size_t i = 0; i < 12; ++i) {
        const uint256_t canonical(*coords[i]);
        uint256_t lo(0);
        uint256_t hi(0);
        std::memcpy(lo.data, canonical.data, 16);
        std::memcpy(hi.data, reinterpret_cast<const uint8_t*>(canonical.data) + 16, 16);
        fields[2 * i] = fr(lo);
        fields[2 * i + 1] = fr(hi);
    }
    return fields;
}

inline fq12 fq12_from_fields(const std::array<fr, 24>& fields)
{
    fq12 value;
    const std::array<fq*, 12> coords = {
        &value.c0.c0.c0, &value.c0.c0.c1, &value.c0.c1.c0, &value.c0.c1.c1, &value.c0.c2.c0, &value.c0.c2.c1,
        &value.c1.c0.c0, &value.c1.c0.c1, &value.c1.c1.c0, &value.c1.c1.c1, &value.c1.c2.c0, &value.c1.c2.c1,
    };
    for (size_t i = 0; i < 12; ++i) {
        const uint256_t lo(fields[2 * i]);
        const uint256_t hi(fields[2 * i + 1]);
        uint256_t canonical(0);
        std::memcpy(canonical.data, lo.data, 16);
        std::memcpy(reinterpret_cast<uint8_t*>(canonical.data) + 16, hi.data, 16);
        *coords[i] = fq(canonical);
    }
    return value;
}

/** @brief GT exponentiation (square-and-multiply; fq12 has no pow member). */
inline fq12 gt_pow(const fq12& base, const uint256_t& exponent)
{
    fq12 result = fq12::one();
    fq12 running = base;
    uint256_t remaining = exponent;
    while (remaining != uint256_t(0)) {
        if ((remaining.data[0] & 1) != 0) {
            result *= running;
        }
        running = running.sqr();
        remaining = remaining >> 1;
    }
    return result;
}

/** @brief Final b-scalar of eq(u_hi) after stride-2 folds with the given odd-half multipliers. */
inline fr eq_fold_scalar(std::span<const fr> u, std::span<const fr> multipliers)
{
    fr result = fr::one();
    for (size_t j = 0; j < u.size(); ++j) {
        result *= (fr::one() - u[j]) + u[j] * multipliers[j];
    }
    return result;
}

/** @brief Final b-scalar of shifted-eq(u_hi) (successor kernel; see pedersen_ipa). */
inline fr shifted_eq_fold_scalar(std::span<const fr> u, std::span<const fr> multipliers)
{
    const size_t ell = u.size();
    std::vector<fr> suffix(ell + 1, fr::one());
    for (size_t j = ell; j-- > 0;) {
        suffix[j] = suffix[j + 1] * ((fr::one() - u[j]) + u[j] * multipliers[j]);
    }
    fr result = fr::zero();
    fr prefix = fr::one();
    for (size_t k = 0; k < ell; ++k) {
        result += prefix * (fr::one() - u[k]) * multipliers[k] * suffix[k + 1];
        prefix *= u[k];
    }
    return result;
}

} // namespace detail

/**
 * @brief Transparent Dory-style commitment key: G1 generators for one matrix row and G2 generators
 * for the row direction.
 *
 * @warning Test-only G2 derivation: the Gamma2_i are scalar multiples of the G2 generator with
 * public deterministic scalars (bb has no cofactor-cleared hash-to-G2). A production setup derives
 * them by hashing to G2, at identical runtime cost; the performance profile measured with this key
 * is that of the real scheme.
 */
class DoryCommitmentKey {
  public:
    explicit DoryCommitmentKey(const DoryConfig& config)
        : config(config)
        , g1_generators(Curve::Group::derive_generators(std::vector<uint8_t>{ 'b', 'b', '_', 'd', 'o', 'r', 'y' },
                                                        config.num_cols()))
    {
        const size_t num_rows = config.num_rows();
        g2_generators.reserve(num_rows);
        for (size_t i = 0; i < num_rows; ++i) {
            std::vector<uint8_t> seed = { 'd', 'o', 'r', 'y', '_', 'g', '2' };
            for (size_t byte = 0; byte < 8; ++byte) {
                seed.push_back(static_cast<uint8_t>((i >> (8 * byte)) & 0xff));
            }
            const auto digest = blake3::blake3s(seed);
            uint256_t value(0);
            std::memcpy(value.data, digest.data(), 32);
            g2_generators.push_back(G2Affine(G2Element(bb::g2::affine_one) * fr(value)));
        }
    }

    DoryGroupData commit_group(std::vector<std::vector<fr>> payload_columns,
                               const std::vector<bool>& /*to_be_shifted*/ = {}) const
    {
        const size_t n = size_t(1) << config.num_variables;
        const size_t num_rows = config.num_rows();
        const size_t num_cols = config.num_cols();
        DoryGroupData data;
        for (auto& column : payload_columns) {
            BB_ASSERT_LTE(column.size(), n, "polynomial too large for the configured size");
            column.resize(n, fr::zero());
            std::vector<Commitment> rows(num_rows);
            for (size_t r = 0; r < num_rows; ++r) {
                std::vector<fr> row(column.begin() + static_cast<std::ptrdiff_t>(r * num_cols),
                                    column.begin() + static_cast<std::ptrdiff_t>((r + 1) * num_cols));
                rows[r] = Commitment::batch_mul(std::span<const Commitment>(g1_generators), std::span<fr>(row));
            }
            data.commitments.push_back(pairing::reduced_ate_pairing_batch(rows.data(), g2_generators.data(), num_rows));
            data.row_commitments.push_back(std::move(rows));
            data.coefficients.push_back(std::move(column));
        }
        return data;
    }

    DoryGroupData commit_group(std::span<const Polynomial<fr>* const> polynomials,
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

    DoryConfig config;
    std::vector<Commitment> g1_generators; // Gamma1, one matrix row
    std::vector<G2Affine> g2_generators;   // Gamma2, the row direction
};

/**
 * @brief Two-tier Dory-style opening (Lee, ePrint 2020/1274 adapted): per chain, the verifier
 * forms Q_0 = Commit_Gamma1(w) from the opened combined row w, and log(R) inner-pairing fold
 * rounds bind Q_0 to the chain's GT commitment T: each round the prover sends the GT cross terms
 * of T's split and the G1 cross terms of every Q-claim; at the end a single pairing checks the
 * surviving row commitment against the folded Gamma2. The b-side folds scalar-wise with the
 * eq/successor closed forms. Chain A carries the unshifted batch (weights b = eq(u_hi)); chain B
 * the to-be-shifted batch with two Q-claims (b and shr(b), the rank-2 shift decomposition).
 */
class DoryProver {
  public:
    struct Claims {
        std::vector<const DoryGroupData*> groups;
        std::vector<DoryColumnRef> unshifted;
        std::vector<fr> unshifted_evaluations;
        std::vector<DoryColumnRef> to_be_shifted;
        std::vector<fr> shifted_evaluations;
        bool send_roots = true;
    };

    template <typename Transcript>
    static void prove(const DoryCommitmentKey& ck,
                      const Claims& claims,
                      std::span<const fr> u,
                      const std::shared_ptr<Transcript>& transcript)
    {
        const DoryConfig& config = ck.config;
        const size_t n = size_t(1) << config.num_variables;
        const size_t num_rows = config.num_rows();
        const size_t num_cols = config.num_cols();
        BB_ASSERT_EQ(u.size(), config.num_variables, "opening point size mismatch");

        if (claims.send_roots) {
            for (size_t g = 0; g < claims.groups.size(); ++g) {
                for (size_t c = 0; c < claims.groups[g]->num_columns(); ++c) {
                    transcript->send_to_verifier(detail::dory_label("root", g, c),
                                                 detail::fq12_to_fields(claims.groups[g]->commitments[c]));
                }
            }
        }
        const fr rho = transcript->template get_challenge<fr>("DORY:rho");

        // Combined arrays and combined row-commitment vectors per chain.
        struct Chain {
            std::vector<fr> array;
            std::vector<GroupElement> rows;
            size_t num_q_claims = 1;
        };
        std::vector<Chain> chains;
        {
            Chain chain_a;
            chain_a.array.assign(n, fr::zero());
            chain_a.rows.assign(num_rows, GroupElement::infinity());
            std::vector<pcs_utils::ScaledTerm> terms_a;
            fr rho_power = fr::one();
            for (size_t i = 0; i < claims.unshifted.size(); ++i) {
                const auto& group = *claims.groups[claims.unshifted[i].group];
                const size_t column = claims.unshifted[i].column;
                terms_a.push_back({ group.coefficients[column].data(), rho_power, false });
                for (size_t r = 0; r < num_rows; ++r) {
                    chain_a.rows[r] += GroupElement(group.row_commitments[column][r]) * rho_power;
                }
                rho_power *= rho;
            }
            pcs_utils::accumulate_scaled(chain_a.array, terms_a);
            chains.push_back(std::move(chain_a));
            if (!claims.to_be_shifted.empty()) {
                Chain chain_b;
                chain_b.num_q_claims = 2;
                chain_b.array.assign(n, fr::zero());
                chain_b.rows.assign(num_rows, GroupElement::infinity());
                std::vector<pcs_utils::ScaledTerm> terms_b;
                for (size_t l = 0; l < claims.to_be_shifted.size(); ++l) {
                    const auto& group = *claims.groups[claims.to_be_shifted[l].group];
                    const size_t column = claims.to_be_shifted[l].column;
                    terms_b.push_back({ group.coefficients[column].data(), rho_power, false });
                    for (size_t r = 0; r < num_rows; ++r) {
                        chain_b.rows[r] += GroupElement(group.row_commitments[column][r]) * rho_power;
                    }
                    rho_power *= rho;
                }
                pcs_utils::accumulate_scaled(chain_b.array, terms_b);
                chains.push_back(std::move(chain_b));
            }
        }

        const std::vector<fr> b = whir::eq_tensor(u.subspan(config.log_num_cols));
        for (size_t chain_index = 0; chain_index < chains.size(); ++chain_index) {
            Chain& chain = chains[chain_index];

            // The opened combined rows (w for weights b; w2 for shr(b) on chain B).
            std::vector<fr> w(num_cols, fr::zero());
            std::vector<fr> w2(num_cols, fr::zero());
            for (size_t r = 0; r < num_rows; ++r) {
                const fr* row = chain.array.data() + r * num_cols;
                for (size_t c = 0; c < num_cols; ++c) {
                    w[c] += b[r] * row[c];
                    if (chain.num_q_claims == 2 && r > 0) {
                        w2[c] += b[r - 1] * row[c];
                    }
                }
            }
            for (size_t c = 0; c < num_cols; ++c) {
                transcript->send_to_verifier(detail::dory_label("w", chain_index, c), w[c]);
            }
            if (chain.num_q_claims == 2) {
                for (size_t c = 0; c < num_cols; ++c) {
                    transcript->send_to_verifier(detail::dory_label("w2", chain_index, c), w2[c]);
                }
            }

            // Inner-pairing fold rounds on the row commitments.
            std::vector<GroupElement> rows = std::move(chain.rows);
            std::vector<G2Element> gamma2(ck.g2_generators.begin(), ck.g2_generators.end());
            std::vector<fr> b_fold(b.begin(), b.end());
            std::vector<fr> b2_fold(num_rows, fr::zero());
            for (size_t r = 1; r < num_rows; ++r) {
                b2_fold[r] = b[r - 1];
            }
            for (size_t round = 0; rows.size() > 1; ++round) {
                const size_t half = rows.size() / 2;
                // Cross terms of T: T_eo = prod e(P_even, Gamma2_odd); T_oe = prod e(P_odd, Gamma2_even).
                std::vector<Commitment> even_rows(half);
                std::vector<Commitment> odd_rows(half);
                std::vector<G2Affine> even_gamma(half);
                std::vector<G2Affine> odd_gamma(half);
                for (size_t t = 0; t < half; ++t) {
                    even_rows[t] = Commitment(rows[2 * t]);
                    odd_rows[t] = Commitment(rows[2 * t + 1]);
                    even_gamma[t] = G2Affine(gamma2[2 * t]);
                    odd_gamma[t] = G2Affine(gamma2[2 * t + 1]);
                }
                const fq12 cross_eo = pairing::reduced_ate_pairing_batch(even_rows.data(), odd_gamma.data(), half);
                const fq12 cross_oe = pairing::reduced_ate_pairing_batch(odd_rows.data(), even_gamma.data(), half);
                transcript->send_to_verifier(detail::dory_label("T_eo", chain_index, round),
                                             detail::fq12_to_fields(cross_eo));
                transcript->send_to_verifier(detail::dory_label("T_oe", chain_index, round),
                                             detail::fq12_to_fields(cross_oe));
                // Cross terms of the Q-claims.
                GroupElement q_plus = GroupElement::infinity();
                GroupElement q_minus = GroupElement::infinity();
                GroupElement q2_plus = GroupElement::infinity();
                GroupElement q2_minus = GroupElement::infinity();
                for (size_t t = 0; t < half; ++t) {
                    q_plus += rows[2 * t + 1] * b_fold[2 * t];
                    q_minus += rows[2 * t] * b_fold[2 * t + 1];
                    if (chain.num_q_claims == 2) {
                        q2_plus += rows[2 * t + 1] * b2_fold[2 * t];
                        q2_minus += rows[2 * t] * b2_fold[2 * t + 1];
                    }
                }
                transcript->send_to_verifier(detail::dory_label("Q_plus", chain_index, round), Commitment(q_plus));
                transcript->send_to_verifier(detail::dory_label("Q_minus", chain_index, round), Commitment(q_minus));
                if (chain.num_q_claims == 2) {
                    transcript->send_to_verifier(detail::dory_label("Q2_plus", chain_index, round),
                                                 Commitment(q2_plus));
                    transcript->send_to_verifier(detail::dory_label("Q2_minus", chain_index, round),
                                                 Commitment(q2_minus));
                }
                const fr x = transcript->template get_challenge<fr>(detail::dory_label("x", chain_index, round));
                const fr x_inv = x.invert();
                for (size_t t = 0; t < half; ++t) {
                    rows[t] = rows[2 * t] + rows[2 * t + 1] * x;
                    gamma2[t] = gamma2[2 * t] + gamma2[2 * t + 1] * x_inv;
                    b_fold[t] = b_fold[2 * t] + x_inv * b_fold[2 * t + 1];
                    b2_fold[t] = b2_fold[2 * t] + x_inv * b2_fold[2 * t + 1];
                }
                rows.resize(half);
                gamma2.resize(half);
                b_fold.resize(half);
                b2_fold.resize(half);
            }
            transcript->send_to_verifier(detail::dory_label("P_final", chain_index, 0), Commitment(rows[0]));
        }
    }
};

/** @brief The mirror of `DoryProver`: two pairings, O(sqrt N) group work, GT fold bookkeeping. */
class DoryVerifier {
  public:
    struct Claims {
        std::vector<size_t> group_num_columns;
        std::vector<std::vector<fq12>> group_commitments; // when non-empty, transcript-bound
        std::vector<DoryColumnRef> unshifted;
        std::vector<fr> unshifted_evaluations;
        std::vector<DoryColumnRef> to_be_shifted;
        std::vector<fr> shifted_evaluations;
    };

    template <typename Transcript>
    static bool verify(const DoryConfig& config,
                       const Claims& claims,
                       std::span<const fr> u,
                       const std::shared_ptr<Transcript>& transcript,
                       const DoryCommitmentKey& ck)
    {
        const size_t num_rows = config.num_rows();
        const size_t num_cols = config.num_cols();
        BB_ASSERT_EQ(u.size(), config.num_variables, "opening point size mismatch");

        std::vector<std::vector<fq12>> group_commitments = claims.group_commitments;
        if (group_commitments.empty()) {
            for (size_t g = 0; g < claims.group_num_columns.size(); ++g) {
                std::vector<fq12> commitments;
                for (size_t c = 0; c < claims.group_num_columns[g]; ++c) {
                    commitments.push_back(
                        detail::fq12_from_fields(transcript->template receive_from_prover<std::array<fr, 24>>(
                            detail::dory_label("root", g, c))));
                }
                group_commitments.push_back(std::move(commitments));
            }
        }
        const fr rho = transcript->template get_challenge<fr>("DORY:rho");

        struct ChainView {
            fq12 commitment = fq12::one();
            fr claimed_evaluation = fr::zero();
            size_t num_q_claims = 1;
        };
        std::vector<ChainView> chains(claims.to_be_shifted.empty() ? 1 : 2);
        fr rho_power = fr::one();
        for (size_t i = 0; i < claims.unshifted.size(); ++i) {
            chains[0].commitment *= detail::gt_pow(
                group_commitments[claims.unshifted[i].group][claims.unshifted[i].column], uint256_t(rho_power));
            chains[0].claimed_evaluation += rho_power * claims.unshifted_evaluations[i];
            rho_power *= rho;
        }
        if (chains.size() > 1) {
            chains[1].num_q_claims = 2;
            for (size_t l = 0; l < claims.to_be_shifted.size(); ++l) {
                chains[1].commitment *=
                    detail::gt_pow(group_commitments[claims.to_be_shifted[l].group][claims.to_be_shifted[l].column],
                                   uint256_t(rho_power));
                chains[1].claimed_evaluation += rho_power * claims.shifted_evaluations[l];
                rho_power *= rho;
            }
        }

        const std::vector<fr> a = whir::eq_tensor(u.subspan(0, config.log_num_cols));
        for (size_t chain_index = 0; chain_index < chains.size(); ++chain_index) {
            ChainView& chain = chains[chain_index];
            std::vector<fr> w(num_cols);
            for (size_t c = 0; c < num_cols; ++c) {
                w[c] = transcript->template receive_from_prover<fr>(detail::dory_label("w", chain_index, c));
            }
            std::vector<fr> w2(num_cols, fr::zero());
            if (chain.num_q_claims == 2) {
                for (size_t c = 0; c < num_cols; ++c) {
                    w2[c] = transcript->template receive_from_prover<fr>(detail::dory_label("w2", chain_index, c));
                }
            }

            // Q_0 = Commit(w) over the G1 generators; the folds must land on b_final * P_final.
            std::vector<fr> w_copy = w;
            GroupElement q = GroupElement(
                Commitment::batch_mul(std::span<const Commitment>(ck.g1_generators), std::span<fr>(w_copy)));
            std::vector<fr> w2_copy = w2;
            GroupElement q2 = chain.num_q_claims == 2
                                  ? GroupElement(Commitment::batch_mul(std::span<const Commitment>(ck.g1_generators),
                                                                       std::span<fr>(w2_copy)))
                                  : GroupElement::infinity();

            fq12 t_accumulator = chain.commitment;
            std::vector<G2Element> gamma2(ck.g2_generators.begin(), ck.g2_generators.end());
            std::vector<fr> multipliers;
            for (size_t round = 0; (num_rows >> round) > 1; ++round) {
                const fq12 cross_eo =
                    detail::fq12_from_fields(transcript->template receive_from_prover<std::array<fr, 24>>(
                        detail::dory_label("T_eo", chain_index, round)));
                const fq12 cross_oe =
                    detail::fq12_from_fields(transcript->template receive_from_prover<std::array<fr, 24>>(
                        detail::dory_label("T_oe", chain_index, round)));
                const Commitment q_plus = transcript->template receive_from_prover<Commitment>(
                    detail::dory_label("Q_plus", chain_index, round));
                const Commitment q_minus = transcript->template receive_from_prover<Commitment>(
                    detail::dory_label("Q_minus", chain_index, round));
                Commitment q2_plus = Commitment::infinity();
                Commitment q2_minus = Commitment::infinity();
                if (chain.num_q_claims == 2) {
                    q2_plus = transcript->template receive_from_prover<Commitment>(
                        detail::dory_label("Q2_plus", chain_index, round));
                    q2_minus = transcript->template receive_from_prover<Commitment>(
                        detail::dory_label("Q2_minus", chain_index, round));
                }
                const fr x = transcript->template get_challenge<fr>(detail::dory_label("x", chain_index, round));
                const fr x_inv = x.invert();
                multipliers.push_back(x_inv);

                // T' = T * T_eo^{x^{-1}} * T_oe^{x}; Q' = Q + x Q_plus + x^{-1} Q_minus.
                t_accumulator *= detail::gt_pow(cross_eo, uint256_t(x_inv)) * detail::gt_pow(cross_oe, uint256_t(x));
                q += GroupElement(q_plus) * x + GroupElement(q_minus) * x_inv;
                if (chain.num_q_claims == 2) {
                    q2 += GroupElement(q2_plus) * x + GroupElement(q2_minus) * x_inv;
                }
                const size_t half = gamma2.size() / 2;
                for (size_t t = 0; t < half; ++t) {
                    gamma2[t] = gamma2[2 * t] + gamma2[2 * t + 1] * x_inv;
                }
                gamma2.resize(half);
            }
            const Commitment p_final =
                transcript->template receive_from_prover<Commitment>(detail::dory_label("P_final", chain_index, 0));

            // Final checks: pairing binds P_final to the folded T; the Q-claims close scalar-side.
            const G2Affine gamma2_final(gamma2[0]);
            if (pairing::reduced_ate_pairing(p_final, gamma2_final) != t_accumulator) {
                return false;
            }
            const auto u_hi = u.subspan(config.log_num_cols);
            const fr b_final = detail::eq_fold_scalar(u_hi, multipliers);
            if (Commitment(q) != Commitment(GroupElement(p_final) * b_final)) {
                return false;
            }
            if (chain.num_q_claims == 2) {
                const fr b2_final = detail::shifted_eq_fold_scalar(u_hi, multipliers);
                if (Commitment(q2) != Commitment(GroupElement(p_final) * b2_final)) {
                    return false;
                }
            }

            // Claim equation (rank-2 shift decomposition on chain B).
            fr total = fr::zero();
            if (chain.num_q_claims == 1) {
                for (size_t c = 0; c < num_cols; ++c) {
                    total += w[c] * a[c];
                }
            } else {
                for (size_t c = 1; c < num_cols; ++c) {
                    total += w[c] * a[c - 1];
                }
                total += a[num_cols - 1] * w2[0];
            }
            if (total != chain.claimed_evaluation) {
                return false;
            }
        }
        return true;
    }
};

} // namespace bb::dory
