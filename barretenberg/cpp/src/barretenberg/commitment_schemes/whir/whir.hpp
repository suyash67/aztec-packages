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
 * @brief Prover-side result of committing one polynomial: the dense coefficient array, and the
 * Merkle tree over its codeword. The commitment sent to the verifier is `tree.root()`.
 */
template <typename Hasher> struct WhirProverData {
    std::vector<fr> coefficients;
    MerkleTree<Hasher> tree;
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

    /** @brief Commit to a dense coefficient span (length <= 2^m; zero-padded). */
    WhirProverData<Hasher> commit(std::span<const fr> coefficients, bool salted = false) const
    {
        const size_t n = size_t(1) << config.num_variables;
        BB_ASSERT_LTE(coefficients.size(), n, "polynomial too large for the configured size");
        std::vector<fr> dense(n, fr::zero());
        std::copy(coefficients.begin(), coefficients.end(), dense.begin());
        std::vector<fr> codeword = rs_encode(dense, domains.get(n << config.log_inv_rate));
        MerkleTree<Hasher> tree(codeword, config.folding_factor_bits, salted);
        return { std::move(dense), std::move(tree) };
    }

    /** @brief Commit to a Honk polynomial (virtual zeros outside its span are honored). */
    WhirProverData<Hasher> commit(const Polynomial<fr>& polynomial, bool salted = false) const
    {
        const size_t n = size_t(1) << config.num_variables;
        BB_ASSERT_LTE(polynomial.end_index(), n, "polynomial too large for the configured size");
        std::vector<fr> dense(n, fr::zero());
        for (size_t i = polynomial.start_index(); i < polynomial.end_index(); ++i) {
            dense[i] = polynomial[i];
        }
        std::vector<fr> codeword = rs_encode(dense, domains.get(n << config.log_inv_rate));
        MerkleTree<Hasher> tree(codeword, config.folding_factor_bits, salted);
        return { std::move(dense), std::move(tree) };
    }

    WhirConfig config;
    mutable RSDomains domains;
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

} // namespace detail

/**
 * @brief One constituent of the (possibly virtual) round-0 oracle: a commitment plus how its opened
 * values enter the batched codeword F (README.md §4.2): scaled by `batching_scalar`, and by x⁻¹ when
 * the claim is for the shift.
 */
template <typename Hasher> struct WhirOracleComponent {
    const MerkleTree<Hasher>* tree; // prover side; nullptr on the verifier
    typename Hasher::Digest root;   // verifier side
    fr batching_scalar;
    bool is_shifted;
};

/**
 * @brief WHIR opening prover for a batch of evaluation claims {P_j(u) = v_j} and shifted claims at a
 * common point u. Implements the iteration of README.md §4.3 and the final phase of §4.4; the
 * transcript schedule is §4.5.
 */
template <typename Hasher> class WhirProver {
  public:
    using Tree = MerkleTree<Hasher>;

    struct Claims {
        std::vector<const WhirProverData<Hasher>*> unshifted;
        std::vector<fr> unshifted_evaluations;
        // to-be-shifted commitments (constant coefficient zero); evaluations are of their shifts
        std::vector<const WhirProverData<Hasher>*> to_be_shifted;
        std::vector<fr> shifted_evaluations;
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
        BB_ASSERT_EQ(u.size(), config.num_variables, "opening point size mismatch");

        // Commitment roots enter the transcript before the batching challenge. In the Honk
        // integration they are already there (sent during earlier rounds); standalone, send them now.
        for (size_t p = 0; p < claims.unshifted.size(); ++p) {
            transcript->send_to_verifier(detail::whir_label("root_u", p),
                                         Hasher::digest_to_fields(claims.unshifted[p]->tree.root()));
        }
        for (size_t p = 0; p < claims.to_be_shifted.size(); ++p) {
            transcript->send_to_verifier(detail::whir_label("root_s", p),
                                         Hasher::digest_to_fields(claims.to_be_shifted[p]->tree.root()));
        }

        const fr rho = transcript->template get_challenge<fr>("WHIR:rho");

        // Round-0 oracle components and the batched array F = ∑_j ρʲ·a_j + ∑_l ρ^{n_u+l}·shift(a_l).
        std::vector<WhirOracleComponent<Hasher>> components;
        std::vector<fr> batched(n, fr::zero());
        fr rho_power = fr::one();
        for (const WhirProverData<Hasher>* data : claims.unshifted) {
            components.push_back({ &data->tree, {}, rho_power, false });
            const fr scalar = rho_power;
            parallel_for_range(n, [&](size_t start, size_t end) {
                for (size_t i = start; i < end; ++i) {
                    batched[i] += scalar * data->coefficients[i];
                }
            });
            rho_power *= rho;
        }
        for (const WhirProverData<Hasher>* data : claims.to_be_shifted) {
            BB_ASSERT_EQ(data->coefficients[0], fr::zero(), "to-be-shifted polynomial must have zero constant term");
            components.push_back({ &data->tree, {}, rho_power, true });
            const fr scalar = rho_power;
            parallel_for_range(n - 1, [&](size_t start, size_t end) {
                for (size_t i = start; i < end; ++i) {
                    batched[i] += scalar * data->coefficients[i + 1];
                }
            });
            rho_power *= rho;
        }

        // Weight table of the initial claim: W = eq(u, ·) over {0,1}^m.
        std::vector<fr> weight_table(n, fr::zero());
        WeightTerm::eq_weight(u, fr::one()).accumulate_table(weight_table);

        std::vector<fr> current = std::move(batched);
        std::deque<Tree> folded_trees; // owns the trees committed during the proof
        std::vector<const Tree*> query_trees;
        for (const auto& component : components) {
            query_trees.push_back(component.tree);
        }

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
            std::vector<fr> codeword = rs_encode(current, next_domain);
            folded_trees.emplace_back(codeword, k, /*salted=*/false);
            transcript->send_to_verifier(detail::whir_label("root_g", i + 1),
                                         Hasher::digest_to_fields(folded_trees.back().root()));

            // 3. Out-of-domain sample.
            const fr z_ood = transcript->template get_challenge<fr>(detail::whir_label("z_ood", i));
            const fr y_ood = polynomial_arithmetic::evaluate(current.data(), z_ood, current.size());
            transcript->send_to_verifier(detail::whir_label("y_ood", i), y_ood);

            // 4. In-domain queries against the round-i oracle.
            const size_t index_bits = round.log_domain_size - k;
            std::vector<size_t> indices(round.num_queries);
            for (size_t s = 0; s < round.num_queries; ++s) {
                const fr challenge = transcript->template get_challenge<fr>(detail::whir_label("query", i, s));
                indices[s] = detail::index_from_challenge(challenge, index_bits);
            }
            for (size_t s = 0; s < round.num_queries; ++s) {
                for (size_t p = 0; p < query_trees.size(); ++p) {
                    send_opening(transcript, query_trees[p]->open(indices[s]), i, s, p);
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

            query_trees = { &folded_trees.back() };
        }

        // Final phase: the remaining polynomial in the clear, plus direct consistency queries.
        for (size_t j = 0; j < current.size(); ++j) {
            transcript->send_to_verifier(detail::whir_label("final", j), current[j]);
        }
        const size_t index_bits = config.final_round.log_domain_size - k;
        for (size_t s = 0; s < config.final_round.num_queries; ++s) {
            const fr challenge = transcript->template get_challenge<fr>(detail::whir_label("fquery", s));
            const size_t idx = detail::index_from_challenge(challenge, index_bits);
            for (size_t p = 0; p < query_trees.size(); ++p) {
                send_opening(transcript, query_trees[p]->open(idx), config.rounds.size(), s, p);
            }
        }
    }

  private:
    template <typename Transcript>
    static void send_opening(const std::shared_ptr<Transcript>& transcript,
                             const typename Tree::Opening& opening,
                             size_t i,
                             size_t s,
                             size_t p)
    {
        const std::string base = "WHIR:open_" + std::to_string(i) + "_" + std::to_string(s) + "_" + std::to_string(p);
        for (size_t t = 0; t < opening.values.size(); ++t) {
            transcript->send_to_verifier(base + ":v" + std::to_string(t), opening.values[t]);
        }
        if (opening.salt) {
            transcript->send_to_verifier(base + ":salt", *opening.salt);
        }
        for (size_t l = 0; l < opening.path.size(); ++l) {
            transcript->send_to_verifier(base + ":p" + std::to_string(l), Hasher::digest_to_fields(opening.path[l]));
        }
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
        std::vector<fr> unshifted_evaluations;
        std::vector<fr> shifted_evaluations;
    };

    template <typename Transcript>
    static bool verify(const WhirConfig& config,
                       const Claims& claims,
                       std::span<const fr> u,
                       const std::shared_ptr<Transcript>& transcript)
    {
        const size_t k = config.folding_factor_bits;
        BB_ASSERT_EQ(u.size(), config.num_variables, "opening point size mismatch");

        // Commitment roots, then the batching challenge and the batched claim value σ₀.
        std::vector<WhirOracleComponent<Hasher>> components;
        for (size_t p = 0; p < claims.unshifted_evaluations.size(); ++p) {
            const auto fields = transcript->template receive_from_prover<std::array<fr, Hasher::DIGEST_NUM_FIELDS>>(
                detail::whir_label("root_u", p));
            components.push_back({ nullptr, Hasher::digest_from_fields(fields), fr::one(), false });
        }
        for (size_t p = 0; p < claims.shifted_evaluations.size(); ++p) {
            const auto fields = transcript->template receive_from_prover<std::array<fr, Hasher::DIGEST_NUM_FIELDS>>(
                detail::whir_label("root_s", p));
            components.push_back({ nullptr, Hasher::digest_from_fields(fields), fr::one(), true });
        }
        const fr rho = transcript->template get_challenge<fr>("WHIR:rho");
        fr sigma = fr::zero();
        fr rho_power = fr::one();
        size_t component_idx = 0;
        for (const fr& evaluation : claims.unshifted_evaluations) {
            components[component_idx++].batching_scalar = rho_power;
            sigma += rho_power * evaluation;
            rho_power *= rho;
        }
        for (const fr& evaluation : claims.shifted_evaluations) {
            components[component_idx++].batching_scalar = rho_power;
            sigma += rho_power * evaluation;
            rho_power *= rho;
        }

        std::vector<WeightTerm> terms = { WeightTerm::eq_weight(u, fr::one()) };

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
                if (!assemble_virtual_values(transcript,
                                             components,
                                             indices[s],
                                             index_bits,
                                             round.log_domain_size,
                                             omega,
                                             eta_inv,
                                             i,
                                             s,
                                             virtual_values)) {
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

            components = { { nullptr, next_root, fr::one(), false } };
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
            if (!assemble_virtual_values(transcript,
                                         components,
                                         idx,
                                         index_bits,
                                         config.final_round.log_domain_size,
                                         omega,
                                         eta_inv,
                                         config.rounds.size(),
                                         s,
                                         virtual_values)) {
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
    /**
     * @brief Receive and authenticate the openings of every oracle component at leaf `idx`, then
     * assemble the virtual codeword values ∑_p scalar_p · v_p[t] (· x_t⁻¹ for shifted components).
     */
    template <typename Transcript>
    static bool assemble_virtual_values(const std::shared_ptr<Transcript>& transcript,
                                        const std::vector<WhirOracleComponent<Hasher>>& components,
                                        size_t idx,
                                        size_t index_bits,
                                        size_t log_domain_size,
                                        const fr& omega,
                                        const fr& eta_inv,
                                        size_t i,
                                        size_t s,
                                        std::vector<fr>& out)
    {
        const size_t arity = size_t(1) << (log_domain_size - index_bits);
        const size_t depth = index_bits;
        out.assign(arity, fr::zero());

        const bool any_shifted =
            std::any_of(components.begin(), components.end(), [](const auto& c) { return c.is_shifted; });
        std::vector<fr> x_inverses;
        if (any_shifted) {
            x_inverses.resize(arity);
            fr x_inv = omega.pow(uint256_t(uint64_t(idx))).invert();
            for (size_t t = 0; t < arity; ++t) {
                x_inverses[t] = x_inv;
                x_inv *= eta_inv;
            }
        }

        for (size_t p = 0; p < components.size(); ++p) {
            typename Tree::Opening opening = receive_opening(transcript, arity, depth, i, s, p);
            if (!Tree::verify(components[p].root, idx, opening)) {
                return false;
            }
            for (size_t t = 0; t < arity; ++t) {
                fr value = components[p].batching_scalar * opening.values[t];
                if (components[p].is_shifted) {
                    value *= x_inverses[t];
                }
                out[t] += value;
            }
        }
        return true;
    }

    template <typename Transcript>
    static typename Tree::Opening receive_opening(
        const std::shared_ptr<Transcript>& transcript, size_t arity, size_t depth, size_t i, size_t s, size_t p)
    {
        const std::string base = "WHIR:open_" + std::to_string(i) + "_" + std::to_string(s) + "_" + std::to_string(p);
        typename Tree::Opening opening;
        opening.values.resize(arity);
        for (size_t t = 0; t < arity; ++t) {
            opening.values[t] = transcript->template receive_from_prover<fr>(base + ":v" + std::to_string(t));
        }
        opening.path.resize(depth);
        for (size_t l = 0; l < depth; ++l) {
            const auto fields = transcript->template receive_from_prover<std::array<fr, Hasher::DIGEST_NUM_FIELDS>>(
                base + ":p" + std::to_string(l));
            opening.path[l] = Hasher::digest_from_fields(fields);
        }
        return opening;
    }
};

} // namespace bb::whir
