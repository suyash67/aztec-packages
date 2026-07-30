#pragma once

#include "barretenberg/commitment_schemes/whir/merkle_tree.hpp"
#include "barretenberg/commitment_schemes/whir/rs_code.hpp"
#include "barretenberg/commitment_schemes/whir/weights.hpp"
#include "barretenberg/commitment_schemes/whir/whir_config.hpp"
#include "barretenberg/common/assert.hpp"
#include "barretenberg/common/thread.hpp"
#include "barretenberg/polynomials/polynomial.hpp"
#include "barretenberg/polynomials/polynomial_arithmetic.hpp"

#include <deque>
#include <mutex>
#include <span>
#include <string>
#include <vector>

namespace bb::whir {

/**
 * @brief Prover-side result of committing a group of polynomials ("columns") into one shared Merkle
 * tree (README.md §4.1, §10): the dense coefficient arrays and the tree over their codewords with
 * column-major interleaved leaves. The commitment sent to the verifier is the single `tree.root()`.
 */
template <typename Hasher> struct WhirGroupData {
    std::vector<std::vector<fr>> coefficients;
    MerkleTree<Hasher> tree;

    size_t num_columns() const { return coefficients.size(); }
};

/** @brief Reference to one column of one committed group, shared by prover and verifier claims. */
struct WhirColumnRef {
    size_t group;
    size_t column;
};

/**
 * @brief Transparent commitment key: the parameter schedule plus cached FFT domains. Replaces the
 * SRS-based `CommitmentKey` of KZG; there is no toxic waste and no setup beyond FFT root tables.
 */
template <typename Hasher> class WhirCommitmentKey {
  public:
    explicit WhirCommitmentKey(const WhirConfig& config)
        : config(config)
    {}

    /**
     * @brief Commit a group of payload columns (each of length <= 2^m, zero-padded) into one tree.
     * @details In zk mode each column's array doubles: the payload occupies the low half and
     * `config.num_blinding_coefficients` fresh random coefficients occupy the high half (offset 2^m,
     * or 2^m + 1 for a to-be-shifted column so the shift contract's zero slot is preserved), and the
     * leaves are salted. See README.md §8.
     *
     * @param payload_columns payload coefficient arrays; consumed
     * @param to_be_shifted per-column flag; empty means all false
     */
    WhirGroupData<Hasher> commit_group(std::vector<std::vector<fr>> payload_columns,
                                       const std::vector<bool>& to_be_shifted = {}) const
    {
        BB_ASSERT(to_be_shifted.empty() || to_be_shifted.size() == payload_columns.size(),
                  "per-column shift flags must match the column count");
        const size_t payload_size = size_t(1) << config.num_payload_variables;
        std::vector<std::vector<fr>> dense;
        dense.reserve(payload_columns.size());
        for (size_t c = 0; c < payload_columns.size(); ++c) {
            BB_ASSERT_LTE(payload_columns[c].size(), payload_size, "polynomial too large for the configured size");
            std::vector<fr> column = std::move(payload_columns[c]);
            column.resize(size_t(1) << config.num_variables, fr::zero());
            add_blinding(column, !to_be_shifted.empty() && to_be_shifted[c]);
            dense.push_back(std::move(column));
        }
        return commit_dense_group(std::move(dense));
    }

    /** @brief Commit a group of Honk polynomials (virtual zeros outside their spans are honored). */
    WhirGroupData<Hasher> commit_group(std::span<const Polynomial<fr>* const> polynomials,
                                       const std::vector<bool>& to_be_shifted = {}) const
    {
        std::vector<std::vector<fr>> payload_columns;
        payload_columns.reserve(polynomials.size());
        for (const Polynomial<fr>* polynomial : polynomials) {
            BB_ASSERT_LTE(polynomial->end_index(),
                          size_t(1) << config.num_payload_variables,
                          "polynomial too large for the configured size");
            std::vector<fr> column(polynomial->end_index(), fr::zero());
            for (size_t i = polynomial->start_index(); i < polynomial->end_index(); ++i) {
                column[i] = (*polynomial)[i];
            }
            payload_columns.push_back(std::move(column));
        }
        return commit_group(std::move(payload_columns), to_be_shifted);
    }

    /** @brief Single-polynomial convenience: a group of one column. */
    WhirGroupData<Hasher> commit(std::span<const fr> coefficients, bool to_be_shifted = false) const
    {
        std::vector<std::vector<fr>> columns;
        columns.emplace_back(coefficients.begin(), coefficients.end());
        return commit_group(std::move(columns), { to_be_shifted });
    }

    /** @brief Commit a group of one full-width uniformly random column (the zk mask, README.md §8). */
    WhirGroupData<Hasher> commit_random_mask() const
    {
        std::vector<fr> mask(size_t(1) << config.num_variables);
        for (fr& value : mask) {
            value = fr::random_element();
        }
        std::vector<std::vector<fr>> columns;
        columns.push_back(std::move(mask));
        return commit_dense_group(std::move(columns));
    }

    WhirConfig config;
    mutable RSDomains domains;

  private:
    WhirGroupData<Hasher> commit_dense_group(std::vector<std::vector<fr>> dense) const
    {
        const size_t n = size_t(1) << config.num_variables;
        const auto& domain = domains.get(n << config.log_inv_rate);
        std::vector<std::vector<fr>> codewords;
        codewords.reserve(dense.size());
        for (const auto& column : dense) {
            BB_ASSERT_EQ(column.size(), n);
            codewords.push_back(rs_encode(column, domain, domains.round_roots()));
        }
        MerkleTree<Hasher> tree(std::move(codewords), config.folding_factor_bits, /*salted=*/config.zk);
        return { std::move(dense), std::move(tree) };
    }

    void add_blinding(std::vector<fr>& dense, bool to_be_shifted) const
    {
        if (!config.zk) {
            return;
        }
        const size_t offset = (size_t(1) << config.num_payload_variables) + (to_be_shifted ? 1 : 0);
        for (size_t i = 0; i < config.num_blinding_coefficients; ++i) {
            dense[offset + i] = fr::random_element();
        }
    }
};

namespace detail {

inline std::string whir_label(const std::string& name, size_t i)
{
    return "WHIR:" + name + "_" + std::to_string(i);
}
inline std::string whir_label(const std::string& name, size_t i, size_t j)
{
    return "WHIR:" + name + "_" + std::to_string(i) + "_" + std::to_string(j);
}

/** @brief h(alpha) for a degree-2 univariate given by evaluations at 0, 1, 2. */
inline fr evaluate_quadratic(const std::array<fr, 3>& h, const fr& alpha)
{
    static const fr two_inv = fr(2).invert();
    const fr a0 = h[0];
    const fr a2 = (h[2] - h[1] - h[1] + h[0]) * two_inv; // second difference / 2
    const fr a1 = h[1] - h[0] - a2;
    return a0 + alpha * (a1 + alpha * a2);
}

/** @brief The three evaluations {h(0), h(1), h(2)} of the round univariate of ∑_b f(X, b)·W(X, b). */
inline std::array<fr, 3> sumcheck_round_univariate(std::span<const fr> f, std::span<const fr> w)
{
    BB_ASSERT_EQ(f.size(), w.size());
    const size_t half = f.size() / 2;
    std::array<fr, 3> result{ fr::zero(), fr::zero(), fr::zero() };
    std::mutex result_mutex;
    parallel_for_range(half, [&](size_t start, size_t end) {
        fr h0 = fr::zero();
        fr h1 = fr::zero();
        fr h2 = fr::zero();
        for (size_t t = start; t < end; ++t) {
            const fr& f0 = f[2 * t];
            const fr& f1 = f[2 * t + 1];
            const fr& w0 = w[2 * t];
            const fr& w1 = w[2 * t + 1];
            h0 += f0 * w0;
            h1 += f1 * w1;
            h2 += (f1 + f1 - f0) * (w1 + w1 - w0); // evaluations at X = 2 by linear extension
        }
        std::scoped_lock lock(result_mutex);
        result[0] += h0;
        result[1] += h1;
        result[2] += h2;
    });
    return result;
}

/** @brief Query index in [0, 2^index_bits) derived from a transcript challenge. */
inline size_t index_from_challenge(const fr& challenge, size_t index_bits)
{
    return static_cast<size_t>(uint256_t(challenge).data[0] & ((uint64_t(1) << index_bits) - 1));
}

/**
 * @brief Number of proof-stream field elements of one opening: the leaf values, the salt (zk), and
 * the authentication path digests.
 */
template <typename Hasher> size_t opening_num_fields(size_t num_columns, size_t arity, size_t depth, bool salted)
{
    return num_columns * arity + (salted ? 1 : 0) + depth * Hasher::DIGEST_NUM_FIELDS;
}

} // namespace detail

/**
 * @brief WHIR opening prover for a batch of evaluation claims {P_j(u) = v_j} and shifted claims at a
 * common point u, over columns committed in shared-tree groups. Implements the iteration of
 * README.md §4.3 and the final phase of §4.4; the transcript schedule is §4.5. Openings travel in
 * the proof stream without Fiat-Shamir absorption — they are bound by the absorbed roots.
 */
template <typename Hasher> class WhirProver {
  public:
    using Tree = MerkleTree<Hasher>;

    struct Claims {
        std::vector<const WhirGroupData<Hasher>*> groups;
        std::vector<WhirColumnRef> unshifted;
        std::vector<fr> unshifted_evaluations;
        // columns with a zero constant coefficient; evaluations are of their shifts
        std::vector<WhirColumnRef> to_be_shifted;
        std::vector<fr> shifted_evaluations;
        // false when the group roots are already bound to the transcript by earlier protocol rounds
        // (the Honk integration); true for standalone use
        bool send_roots = true;
    };

    template <typename Transcript>
    static void prove(const WhirCommitmentKey<Hasher>& ck,
                      const Claims& claims,
                      std::span<const fr> u,
                      const std::shared_ptr<Transcript>& transcript)
    {
        const WhirConfig& config = ck.config;
        const size_t k = config.folding_factor_bits;
        const size_t n = size_t(1) << config.num_variables;
        BB_ASSERT_EQ(u.size(), config.num_payload_variables, "opening point size mismatch");
        BB_ASSERT_EQ(claims.unshifted.size(), claims.unshifted_evaluations.size());
        BB_ASSERT_EQ(claims.to_be_shifted.size(), claims.shifted_evaluations.size());

        // In zk mode every claim lifts to (u, 0): the appended top variable selects the payload half
        // of the blinded arrays (README.md §8).
        std::vector<fr> u_ext(u.begin(), u.end());
        if (config.zk) {
            u_ext.push_back(fr::zero());
        }

        if (claims.send_roots) {
            for (size_t g = 0; g < claims.groups.size(); ++g) {
                transcript->send_to_verifier(detail::whir_label("root", g),
                                             Hasher::digest_to_fields(claims.groups[g]->tree.root()));
            }
        }

        // zk: commit a uniformly random mask polynomial and reveal its (independent, uniform)
        // claimed evaluation; batched into F below, it makes every post-ρ message witness-independent.
        std::vector<const WhirGroupData<Hasher>*> groups = claims.groups;
        std::optional<WhirGroupData<Hasher>> mask_group;
        if (config.zk) {
            mask_group = ck.commit_random_mask();
            transcript->send_to_verifier(std::string("WHIR:root_mask"),
                                         Hasher::digest_to_fields(mask_group->tree.root()));
            const fr mask_evaluation =
                Polynomial<fr>(std::span<const fr>(mask_group->coefficients[0])).evaluate_mle(u_ext);
            transcript->send_to_verifier(std::string("WHIR:mask_eval"), mask_evaluation);
            groups.push_back(&*mask_group);
        }

        const fr rho = transcript->template get_challenge<fr>("WHIR:rho");

        // Batched array F = ∑_j ρʲ·a_j + ∑_l ρ^{n_u+l}·shift(a_l) (+ the mask with the last power).
        std::vector<fr> batched(n, fr::zero());
        fr rho_power = fr::one();
        auto accumulate_column = [&](const std::vector<fr>& column, const fr& scalar, bool shifted) {
            parallel_for_range(n - (shifted ? 1 : 0), [&](size_t start, size_t end) {
                for (size_t i = start; i < end; ++i) {
                    batched[i] += scalar * column[i + (shifted ? 1 : 0)];
                }
            });
        };
        for (const WhirColumnRef& ref : claims.unshifted) {
            accumulate_column(claims.groups[ref.group]->coefficients[ref.column], rho_power, false);
            rho_power *= rho;
        }
        for (const WhirColumnRef& ref : claims.to_be_shifted) {
            const std::vector<fr>& column = claims.groups[ref.group]->coefficients[ref.column];
            BB_ASSERT_EQ(column[0], fr::zero(), "to-be-shifted polynomial must have zero constant term");
            accumulate_column(column, rho_power, true);
            rho_power *= rho;
        }
        if (config.zk) {
            accumulate_column(mask_group->coefficients[0], rho_power, false);
        }

        // Weight table of the initial claim: W = eq(u, ·) over the committed variables.
        std::vector<fr> weight_table(n, fr::zero());
        WeightTerm::eq_weight(u_ext, fr::one()).accumulate_table(weight_table);

        std::vector<fr> current = std::move(batched);
        std::deque<Tree> folded_trees; // owns the trees committed during the proof

        // Openings of the current oracle at a query index: all round-0 groups, or the last folded tree.
        auto send_query_openings = [&](size_t idx) {
            if (folded_trees.empty()) {
                for (const WhirGroupData<Hasher>* group : groups) {
                    send_opening(transcript, group->tree.open(idx));
                }
            } else {
                send_opening(transcript, folded_trees.back().open(idx));
            }
        };

        for (size_t i = 0; i < config.rounds.size(); ++i) {
            const WhirRound& round = config.rounds[i];

            // 1. k sumcheck rounds; array and weight table fold at each challenge.
            for (size_t j = 0; j < k; ++j) {
                const auto h = detail::sumcheck_round_univariate(current, weight_table);
                transcript->send_to_verifier(detail::whir_label("sc", i, j), h);
                const fr alpha = transcript->template get_challenge<fr>(detail::whir_label("alpha", i, j));
                fold_array_in_place(current, alpha);
                fold_array_in_place(weight_table, alpha);
            }

            // 2. Commit the folded polynomial on the halved domain.
            const auto& next_domain = ck.domains.get(size_t(1) << (round.log_domain_size - 1));
            std::vector<fr> codeword = rs_encode(current, next_domain, ck.domains.round_roots());
            folded_trees.emplace_back(std::move(codeword), k, /*salted=*/false);
            transcript->send_to_verifier(detail::whir_label("root_g", i + 1),
                                         Hasher::digest_to_fields(folded_trees.back().root()));

            // 3. Out-of-domain sample.
            const fr z_ood = transcript->template get_challenge<fr>(detail::whir_label("z_ood", i));
            const fr y_ood = polynomial_arithmetic::evaluate(current.data(), z_ood, current.size());
            transcript->send_to_verifier(detail::whir_label("y_ood", i), y_ood);

            // 4. In-domain queries against the round-i oracle. Note the openings are needed for step
            // 2 of the NEXT iteration's committed tree only when i = 0... they always target the
            // previous oracle, so they are emitted before the trees rotate below.
            const size_t index_bits = round.log_domain_size - k;
            std::vector<size_t> indices(round.num_queries);
            for (size_t s = 0; s < round.num_queries; ++s) {
                const fr challenge = transcript->template get_challenge<fr>(detail::whir_label("query", i, s));
                indices[s] = detail::index_from_challenge(challenge, index_bits);
            }
            const bool query_folded_oracle = (i > 0);
            for (size_t s = 0; s < round.num_queries; ++s) {
                if (query_folded_oracle) {
                    send_opening(transcript, folded_trees[i - 1].open(indices[s]));
                } else {
                    for (const WhirGroupData<Hasher>* group : groups) {
                        send_opening(transcript, group->tree.open(indices[s]));
                    }
                }
            }

            // 5. Combine: γ-batch the OOD and query claims into the weight for the next iteration.
            const fr gamma = transcript->template get_challenge<fr>(detail::whir_label("gamma", i));
            const size_t next_vars = round.num_variables - k;
            const fr omega = fr::get_root_of_unity(round.log_domain_size);
            fr gamma_power = gamma;
            WeightTerm::pow_weight(z_ood, next_vars, gamma_power).accumulate_table(weight_table);
            for (size_t s = 0; s < round.num_queries; ++s) {
                gamma_power *= gamma;
                const fr folded_point = omega.pow(uint256_t(uint64_t(indices[s])) << k);
                WeightTerm::pow_weight(folded_point, next_vars, gamma_power).accumulate_table(weight_table);
            }
        }

        // Final phase: the remaining polynomial in the clear, plus direct consistency queries.
        for (size_t j = 0; j < current.size(); ++j) {
            transcript->send_to_verifier(detail::whir_label("final", j), current[j]);
        }
        const size_t index_bits = config.final_round.log_domain_size - k;
        for (size_t s = 0; s < config.final_round.num_queries; ++s) {
            const fr challenge = transcript->template get_challenge<fr>(detail::whir_label("fquery", s));
            send_query_openings(detail::index_from_challenge(challenge, index_bits));
        }
    }

  private:
    template <typename Transcript>
    static void send_opening(const std::shared_ptr<Transcript>& transcript, const typename Tree::Opening& opening)
    {
        std::vector<fr> flat;
        flat.reserve(opening.values.size() + 1 + opening.path.size() * Hasher::DIGEST_NUM_FIELDS);
        flat.insert(flat.end(), opening.values.begin(), opening.values.end());
        if (opening.salt) {
            flat.push_back(*opening.salt);
        }
        for (const auto& digest : opening.path) {
            const auto fields = Hasher::digest_to_fields(digest);
            flat.insert(flat.end(), fields.begin(), fields.end());
        }
        transcript->send_unhashed_to_verifier(flat);
    }
};

/**
 * @brief WHIR opening verifier; the mirror of `WhirProver`. Returns false on any failed check
 * (sumcheck consistency, Merkle authentication, final-polynomial consistency, final claim).
 */
template <typename Hasher> class WhirVerifier {
  public:
    using Tree = MerkleTree<Hasher>;
    using Digest = typename Hasher::Digest;

    struct Claims {
        std::vector<size_t> group_num_columns; // leaf layout of each round-0 group
        std::vector<WhirColumnRef> unshifted;
        std::vector<fr> unshifted_evaluations;
        std::vector<WhirColumnRef> to_be_shifted;
        std::vector<fr> shifted_evaluations;
        // When non-empty, the group roots are already known (bound to the transcript by earlier
        // protocol rounds, as in the Honk integration) and are not read from the proof stream.
        std::vector<Digest> group_roots;
    };

    template <typename Transcript>
    static bool verify(const WhirConfig& config,
                       const Claims& claims,
                       std::span<const fr> u,
                       const std::shared_ptr<Transcript>& transcript)
    {
        const size_t k = config.folding_factor_bits;
        BB_ASSERT_EQ(u.size(), config.num_payload_variables, "opening point size mismatch");
        BB_ASSERT_EQ(claims.unshifted.size(), claims.unshifted_evaluations.size());
        BB_ASSERT_EQ(claims.to_be_shifted.size(), claims.shifted_evaluations.size());

        std::vector<fr> u_ext(u.begin(), u.end());
        if (config.zk) {
            u_ext.push_back(fr::zero());
        }

        // Group roots (and the zk mask root/evaluation), then the batching challenge and σ₀.
        std::vector<size_t> group_columns = claims.group_num_columns;
        std::vector<Digest> roots = claims.group_roots;
        if (roots.empty()) {
            for (size_t g = 0; g < group_columns.size(); ++g) {
                const auto fields = transcript->template receive_from_prover<std::array<fr, Hasher::DIGEST_NUM_FIELDS>>(
                    detail::whir_label("root", g));
                roots.push_back(Hasher::digest_from_fields(fields));
            }
        }
        BB_ASSERT_EQ(roots.size(), group_columns.size(), "one root per group required");
        fr mask_evaluation = fr::zero();
        if (config.zk) {
            const auto fields = transcript->template receive_from_prover<std::array<fr, Hasher::DIGEST_NUM_FIELDS>>(
                std::string("WHIR:root_mask"));
            roots.push_back(Hasher::digest_from_fields(fields));
            group_columns.push_back(1);
            mask_evaluation = transcript->template receive_from_prover<fr>(std::string("WHIR:mask_eval"));
        }

        const fr rho = transcript->template get_challenge<fr>("WHIR:rho");
        // Per-column RLC contributions to the virtual round-0 oracle, in ρ-power order.
        struct Contribution {
            WhirColumnRef ref;
            fr scalar;
            bool shifted;
        };
        std::vector<Contribution> contributions;
        fr sigma = fr::zero();
        fr rho_power = fr::one();
        for (size_t j = 0; j < claims.unshifted.size(); ++j) {
            contributions.push_back({ claims.unshifted[j], rho_power, false });
            sigma += rho_power * claims.unshifted_evaluations[j];
            rho_power *= rho;
        }
        for (size_t l = 0; l < claims.to_be_shifted.size(); ++l) {
            contributions.push_back({ claims.to_be_shifted[l], rho_power, true });
            sigma += rho_power * claims.shifted_evaluations[l];
            rho_power *= rho;
        }
        if (config.zk) {
            contributions.push_back({ { group_columns.size() - 1, 0 }, rho_power, false });
            sigma += rho_power * mask_evaluation;
        }

        std::vector<WeightTerm> terms = { WeightTerm::eq_weight(u_ext, fr::one()) };
        std::optional<Digest> folded_root; // the g_i oracle once i >= 1

        for (size_t i = 0; i < config.rounds.size(); ++i) {
            const WhirRound& round = config.rounds[i];

            // 1. Sumcheck rounds.
            std::vector<fr> alphas(k);
            for (size_t j = 0; j < k; ++j) {
                const auto h =
                    transcript->template receive_from_prover<std::array<fr, 3>>(detail::whir_label("sc", i, j));
                if (h[0] + h[1] != sigma) {
                    return false;
                }
                alphas[j] = transcript->template get_challenge<fr>(detail::whir_label("alpha", i, j));
                sigma = detail::evaluate_quadratic(h, alphas[j]);
                for (WeightTerm& term : terms) {
                    term.bind_variable(alphas[j]);
                }
            }

            // 2. Folded-oracle root; 3. OOD sample.
            const auto root_fields =
                transcript->template receive_from_prover<std::array<fr, Hasher::DIGEST_NUM_FIELDS>>(
                    detail::whir_label("root_g", i + 1));
            const Digest next_root = Hasher::digest_from_fields(root_fields);
            const fr z_ood = transcript->template get_challenge<fr>(detail::whir_label("z_ood", i));
            const fr y_ood = transcript->template receive_from_prover<fr>(detail::whir_label("y_ood", i));

            // 4. Queries: authenticate openings, assemble virtual values, fold cosets.
            const size_t index_bits = round.log_domain_size - k;
            std::vector<size_t> indices(round.num_queries);
            for (size_t s = 0; s < round.num_queries; ++s) {
                const fr challenge = transcript->template get_challenge<fr>(detail::whir_label("query", i, s));
                indices[s] = detail::index_from_challenge(challenge, index_bits);
            }
            const fr omega = fr::get_root_of_unity(round.log_domain_size);
            const fr eta_inv = omega.pow(uint256_t(uint64_t(1)) << index_bits).invert();
            std::vector<fr> folded_values(round.num_queries);
            for (size_t s = 0; s < round.num_queries; ++s) {
                std::vector<fr> virtual_values;
                const bool ok =
                    folded_root
                        ? read_folded_opening(transcript, config, *folded_root, indices[s], index_bits, virtual_values)
                        : read_round0_openings(transcript,
                                               config,
                                               roots,
                                               group_columns,
                                               contributions,
                                               indices[s],
                                               index_bits,
                                               omega,
                                               eta_inv,
                                               virtual_values);
                if (!ok) {
                    return false;
                }
                const fr x_base_inv = omega.pow(uint256_t(uint64_t(indices[s]))).invert();
                folded_values[s] = fold_coset(virtual_values, x_base_inv, eta_inv, alphas);
            }

            // 5. Combine.
            const fr gamma = transcript->template get_challenge<fr>(detail::whir_label("gamma", i));
            const size_t next_vars = round.num_variables - k;
            fr gamma_power = gamma;
            sigma += gamma_power * y_ood;
            terms.push_back(WeightTerm::pow_weight(z_ood, next_vars, gamma_power));
            for (size_t s = 0; s < round.num_queries; ++s) {
                gamma_power *= gamma;
                sigma += gamma_power * folded_values[s];
                const fr folded_point = omega.pow(uint256_t(uint64_t(indices[s])) << k);
                terms.push_back(WeightTerm::pow_weight(folded_point, next_vars, gamma_power));
            }

            folded_root = next_root;
        }

        // Final phase: clear polynomial, oracle consistency at queried cosets, and the claim itself.
        const size_t final_size = size_t(1) << config.final_round.num_variables;
        std::vector<fr> final_poly(final_size);
        for (size_t j = 0; j < final_size; ++j) {
            final_poly[j] = transcript->template receive_from_prover<fr>(detail::whir_label("final", j));
        }
        const size_t index_bits = config.final_round.log_domain_size - k;
        const fr omega = fr::get_root_of_unity(config.final_round.log_domain_size);
        const fr eta = omega.pow(uint256_t(uint64_t(1)) << index_bits);
        const fr eta_inv = eta.invert();
        for (size_t s = 0; s < config.final_round.num_queries; ++s) {
            const fr challenge = transcript->template get_challenge<fr>(detail::whir_label("fquery", s));
            const size_t idx = detail::index_from_challenge(challenge, index_bits);
            std::vector<fr> virtual_values;
            const bool ok = folded_root
                                ? read_folded_opening(transcript, config, *folded_root, idx, index_bits, virtual_values)
                                : read_round0_openings(transcript,
                                                       config,
                                                       roots,
                                                       group_columns,
                                                       contributions,
                                                       idx,
                                                       index_bits,
                                                       omega,
                                                       eta_inv,
                                                       virtual_values);
            if (!ok) {
                return false;
            }
            // Every value of the opened coset must match the clear polynomial.
            fr point = omega.pow(uint256_t(uint64_t(idx)));
            for (const fr& value : virtual_values) {
                if (value != polynomial_arithmetic::evaluate(final_poly.data(), point, final_size)) {
                    return false;
                }
                point *= eta;
            }
        }

        fr total = fr::zero();
        for (const WeightTerm& term : terms) {
            total += term.weighted_sum(final_poly);
        }
        return total == sigma;
    }

  private:
    template <typename Transcript>
    static typename Tree::Opening read_opening(
        const std::shared_ptr<Transcript>& transcript, size_t num_columns, size_t arity, size_t depth, bool salted)
    {
        const std::vector<fr> flat = transcript->receive_unhashed_from_prover(
            detail::opening_num_fields<Hasher>(num_columns, arity, depth, salted));
        typename Tree::Opening opening;
        size_t cursor = num_columns * arity;
        opening.values.assign(flat.begin(), flat.begin() + static_cast<std::ptrdiff_t>(cursor));
        if (salted) {
            opening.salt = flat[cursor++];
        }
        opening.path.reserve(depth);
        for (size_t l = 0; l < depth; ++l) {
            opening.path.push_back(
                Hasher::digest_from_fields(std::span<const fr>(flat).subspan(cursor, Hasher::DIGEST_NUM_FIELDS)));
            cursor += Hasher::DIGEST_NUM_FIELDS;
        }
        return opening;
    }

    /** @brief Open every round-0 group at `idx` and assemble the RLC'd virtual coset values. */
    template <typename Transcript, typename ContributionList>
    static bool read_round0_openings(const std::shared_ptr<Transcript>& transcript,
                                     const WhirConfig& config,
                                     const std::vector<Digest>& roots,
                                     const std::vector<size_t>& group_columns,
                                     const ContributionList& contributions,
                                     size_t idx,
                                     size_t index_bits,
                                     const fr& omega,
                                     const fr& eta_inv,
                                     std::vector<fr>& out)
    {
        const size_t arity = size_t(1) << config.folding_factor_bits;
        std::vector<typename Tree::Opening> openings;
        openings.reserve(roots.size());
        for (size_t g = 0; g < roots.size(); ++g) {
            openings.push_back(read_opening(transcript, group_columns[g], arity, index_bits, config.zk));
            if (!Tree::verify(roots[g], idx, openings.back())) {
                return false;
            }
        }

        const bool any_shifted = std::any_of(
            contributions.begin(), contributions.end(), [](const auto& contribution) { return contribution.shifted; });
        std::vector<fr> x_inverses;
        if (any_shifted) {
            x_inverses.resize(arity);
            fr x_inv = omega.pow(uint256_t(uint64_t(idx))).invert();
            for (size_t t = 0; t < arity; ++t) {
                x_inverses[t] = x_inv;
                x_inv *= eta_inv;
            }
        }

        out.assign(arity, fr::zero());
        for (const auto& contribution : contributions) {
            const std::vector<fr>& values = openings[contribution.ref.group].values;
            for (size_t t = 0; t < arity; ++t) {
                fr value = contribution.scalar * values[contribution.ref.column * arity + t];
                if (contribution.shifted) {
                    value *= x_inverses[t];
                }
                out[t] += value;
            }
        }
        return true;
    }

    /** @brief Open the single-column folded oracle g_i at `idx`; out = the coset values. */
    template <typename Transcript>
    static bool read_folded_opening(const std::shared_ptr<Transcript>& transcript,
                                    const WhirConfig& config,
                                    const Digest& root,
                                    size_t idx,
                                    size_t index_bits,
                                    std::vector<fr>& out)
    {
        const size_t arity = size_t(1) << config.folding_factor_bits;
        const typename Tree::Opening opening = read_opening(transcript, 1, arity, index_bits, /*salted=*/false);
        if (!Tree::verify(root, idx, opening)) {
            return false;
        }
        out = opening.values;
        return true;
    }
};

} // namespace bb::whir
