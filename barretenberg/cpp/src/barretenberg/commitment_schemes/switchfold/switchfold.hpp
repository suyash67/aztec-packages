#pragma once

#include "barretenberg/commitment_schemes/ligero/ligero.hpp"
#include "barretenberg/commitment_schemes/whir/merkle_tree.hpp"
#include "barretenberg/commitment_schemes/whir/rs_code.hpp"
#include "barretenberg/commitment_schemes/whir/weights.hpp"
#include "barretenberg/common/assert.hpp"
#include "barretenberg/common/thread.hpp"
#include "barretenberg/numeric/bitop/get_msb.hpp"

#include <map>
#include <mutex>
#include <span>
#include <string>
#include <vector>

namespace bb::switchfold {

using whir::eq_tensor;
using whir::MerkleTree;
using whir::rs_encode;
using whir::RSDomains;
using whir::WeightTerm;

/** @brief Reference to one polynomial of one committed group (duck-type shared with the suite). */
using SwitchFoldColumnRef = ligero::LigeroColumnRef;

/**
 * @brief SwitchFold parameters (eprint 2026/1489).
 *
 * @details The payload commitment is Ligero's: a `num_rows x num_cols` matrix per polynomial with
 * RS-encoded rows, interleaved into one Merkle tree. The opening differs. Ligero transmits the
 * combined rows in the clear, the O(sqrt N) term that dominates its proof; SwitchFold treats them
 * as a fresh message and *recursively code-switches* down a sequence of RS codes whose block
 * lengths shrink by `2^log_switch_factor` per level, until a base message small enough to transmit
 * remains. Because nothing of row length is sent, the proof-size-optimal shape moves: `num_cols`
 * goes as large as the field's 2-adicity allows, which shrinks the payload leaves (each of which
 * holds `num_rows x num_polynomials` values).
 */
struct SwitchFoldConfig {
    size_t num_variables;
    size_t log_num_cols;
    size_t log_inv_rate;
    size_t security_bits;
    size_t num_queries;
    size_t log_switch_factor; // message length divides by 2^this per code-switching level
    size_t log_base_length;   // descent stops once the message is at most 2^this

    size_t num_rows() const { return size_t(1) << (num_variables - log_num_cols); }
    size_t num_cols() const { return size_t(1) << log_num_cols; }
    size_t codeword_length() const { return size_t(1) << (log_num_cols + log_inv_rate); }
    size_t switch_factor() const { return size_t(1) << log_switch_factor; }

    static SwitchFoldConfig create(size_t num_variables, size_t security_bits = 100, size_t log_inv_rate = 2)
    {
        BB_ASSERT_GT(num_variables, size_t(4));
        const size_t num_queries = (security_bits + log_inv_rate - 1) / log_inv_rate;
        // The descent input carries two extra segment variables on top of the row length, and its
        // codeword must stay inside BN254's 2-adicity.
        size_t log_num_cols = num_variables - 1;
        while (log_num_cols + 2 + log_inv_rate > 28) {
            --log_num_cols;
        }
        return { num_variables, log_num_cols, log_inv_rate, security_bits, num_queries, 2, 4 };
    }

    operator ligero::LigeroConfig() const
    {
        return { num_variables, log_num_cols, log_inv_rate, security_bits, num_queries };
    }
};

template <typename Hasher> using SwitchFoldGroupData = ligero::LigeroGroupData<Hasher>;

/** @brief Ligero's commitment key, carrying the SwitchFold parameters the descent needs. */
template <typename Hasher> class SwitchFoldCommitmentKey {
  public:
    explicit SwitchFoldCommitmentKey(const SwitchFoldConfig& config)
        : config(config)
        , inner(config)
    {}

    template <typename Columns>
    SwitchFoldGroupData<Hasher> commit_group(Columns&& columns, const std::vector<bool>& to_be_shifted = {}) const
    {
        return inner.commit_group(std::forward<Columns>(columns), to_be_shifted);
    }

    SwitchFoldConfig config;
    ligero::LigeroCommitmentKey<Hasher> inner;
};

namespace detail {

inline std::string sf_label(const std::string& name, size_t level)
{
    return "SWITCHFOLD:" + name + "_" + std::to_string(level);
}

inline std::string sf_label(const std::string& name, size_t level, size_t index)
{
    return sf_label(name, level) + "_" + std::to_string(index);
}

inline size_t index_from_challenge(const fr& challenge, size_t index_bits)
{
    return static_cast<size_t>(uint256_t(challenge).data[0] & ((uint64_t(1) << index_bits) - 1));
}

/** @brief One claim about the current message: a weight given as a sum of terms, and its value. */
struct Claim {
    std::vector<WeightTerm> terms;
    fr value;
};

/**
 * @brief shr(eq(u))[c] = eq(u, c-1) as a sum of `u.size()` eq terms, appended to `out`.
 * @details Adding one to c flips a maximal run of low 1-bits to 0 and the next 0-bit to 1, so the
 * vector splits by carry position k: bits below k are pinned to 0, bit k is pinned to 1, and bits
 * above k keep the ordinary eq factor. Pinning a bit is an eq factor at the constant 0 or 1.
 */
inline void append_shifted_eq_terms(std::span<const fr> u, const fr& scale, std::vector<WeightTerm>& out)
{
    fr prefix = fr::one(); // prod_{j<k} u_j
    for (size_t k = 0; k < u.size(); ++k) {
        std::vector<fr> point(u.size(), fr::zero());
        point[k] = fr::one();
        for (size_t j = k + 1; j < u.size(); ++j) {
            point[j] = u[j];
        }
        out.push_back(WeightTerm::eq_weight(point, scale * prefix * (fr::one() - u[k])));
        prefix *= u[k];
    }
}

/** @brief Append the two constant eq factors selecting `segment` of a four-segment vector. */
inline void select_segment(WeightTerm& term, size_t segment)
{
    term.append_eq_variable(((segment >> 0) & 1) != 0 ? fr::one() : fr::zero());
    term.append_eq_variable(((segment >> 1) & 1) != 0 ? fr::one() : fr::zero());
}

/** @brief One sumcheck round of sum_b m(b) w(b), reported at X = 0 and X = 2. */
inline std::pair<fr, fr> inner_product_round(std::span<const fr> m, std::span<const fr> w)
{
    const size_t half = m.size() / 2;
    fr h0 = fr::zero();
    fr h2 = fr::zero();
    std::mutex accumulator_mutex;
    parallel_for_range(half, [&](size_t start, size_t end) {
        fr local0 = fr::zero();
        fr local2 = fr::zero();
        for (size_t b = start; b < end; ++b) {
            const fr& m_even = m[2 * b];
            const fr& m_odd = m[(2 * b) + 1];
            const fr& w_even = w[2 * b];
            const fr& w_odd = w[(2 * b) + 1];
            local0 += m_even * w_even;
            local2 += (m_odd + m_odd - m_even) * (w_odd + w_odd - w_even);
        }
        const std::scoped_lock lock(accumulator_mutex);
        h0 += local0;
        h2 += local2;
    });
    return { h0, h2 };
}

/** @brief Evaluate the quadratic through (0,h0), (1,h1), (2,h2) at alpha. */
inline fr interpolate_quadratic(const fr& h0, const fr& h1, const fr& h2, const fr& alpha)
{
    static const fr two_inv = fr(2).invert();
    const fr c2 = ((h0 + h2) * two_inv) - h1;
    const fr c1 = h1 - h0 - c2;
    return h0 + (alpha * (c1 + (alpha * c2)));
}

/** @brief Batch a claim list with powers of mu into one term list and one value. */
inline std::pair<std::vector<WeightTerm>, fr> batch_claims(const std::vector<Claim>& claims, const fr& mu)
{
    std::vector<WeightTerm> terms;
    fr value = fr::zero();
    fr power = fr::one();
    for (const Claim& claim : claims) {
        for (WeightTerm term : claim.terms) {
            term.scale(power);
            terms.push_back(term);
        }
        value += power * claim.value;
        power *= mu;
    }
    return { terms, value };
}

/**
 * @brief Sum every term's table. Terms are split across threads with private accumulators: the
 * batched claim carries one term per query per segment, so this is the prover's dominant cost.
 */
inline std::vector<fr> weight_table(const std::vector<WeightTerm>& terms, size_t size)
{
    std::vector<fr> table(size, fr::zero());
    std::mutex merge_mutex;
    parallel_for_range(terms.size(), [&](size_t start, size_t end) {
        std::vector<fr> local(size, fr::zero());
        std::vector<fr> scratch(size);
        for (size_t t = start; t < end; ++t) {
            terms[t].accumulate_table(local, scratch);
        }
        const std::scoped_lock lock(merge_mutex);
        for (size_t i = 0; i < size; ++i) {
            table[i] += local[i];
        }
    });
    return table;
}

/** @brief sum of every term evaluated at the fully-bound point. */
inline fr weight_value(const std::vector<WeightTerm>& terms, std::span<const fr> point)
{
    fr total = fr::zero();
    for (WeightTerm term : terms) {
        for (const fr& value : point) {
            term.bind_variable(value);
        }
        total += term.evaluate(0);
    }
    return total;
}

/** @brief Per-(group, column) rho scalars for the unshifted and shifted chains. */
struct ChainScalars {
    fr scalar_u;
    fr scalar_s;
};

template <typename Refs>
inline std::map<std::pair<size_t, size_t>, ChainScalars> chain_scalars(const Refs& unshifted,
                                                                       const Refs& to_be_shifted,
                                                                       const fr& rho)
{
    std::map<std::pair<size_t, size_t>, ChainScalars> scalars;
    fr rho_power = fr::one();
    for (const auto& ref : unshifted) {
        scalars[{ ref.group, ref.column }].scalar_u += rho_power;
        rho_power *= rho;
    }
    for (const auto& ref : to_be_shifted) {
        scalars[{ ref.group, ref.column }].scalar_s += rho_power;
        rho_power *= rho;
    }
    return scalars;
}

} // namespace detail

/**
 * @brief SwitchFold: Ligero's commitment, then a recursive code-switching descent.
 *
 * @details The payload phase reuses `ligero/`: the rho-combined claim set yields the three combined
 * rows of the rank-2 shift decomposition (unshifted, shifted, shifted-carry). Those are laid out as
 * the four segments of one descent message `M` of length `4 * num_cols`, the segment index sitting
 * in the top two variables, so that every claim about `M` is a tensor:
 *
 *  - the evaluation claim, whose shifted part expands into `log num_cols` eq terms;
 *  - one code-switching claim per payload query and segment,
 *    `<M, e_segment (x) pow(omega^s)> = Enc(row)[s]`, the right-hand side derived by the verifier
 *    from the opened leaf. For a Reed-Solomon code the generator-matrix row at a domain point *is*
 *    the pow tensor, so this instantiation needs neither a generator-matrix commitment nor the
 *    paper's accumulation module - the two mechanisms Brakedown's sparse generator forces.
 *
 * Each descent level then commits the message reshaped into `2^log_switch_factor` columns under the
 * next, shorter RS code, opens it at that level's queries, batches the pending claims with a
 * challenge and collapses them to a single tensor claim with one inner-product sumcheck (paper
 * Fact 1), and splits that claim into a claim on the column combination - the next level's message.
 * This level's queries become the next level's code-switching claims.
 */
template <typename Hasher> class SwitchFoldProver {
  public:
    using Tree = MerkleTree<Hasher>;

    struct Claims {
        std::vector<const SwitchFoldGroupData<Hasher>*> groups;
        std::vector<SwitchFoldColumnRef> unshifted;
        std::vector<fr> unshifted_evaluations;
        std::vector<SwitchFoldColumnRef> to_be_shifted;
        std::vector<fr> shifted_evaluations;
        bool send_roots = true;
    };

    template <typename Transcript>
    static void prove(const SwitchFoldCommitmentKey<Hasher>& ck,
                      const Claims& claims,
                      std::span<const fr> u,
                      const std::shared_ptr<Transcript>& transcript)
    {
        const SwitchFoldConfig& config = ck.config;
        const size_t num_rows = config.num_rows();
        const size_t num_cols = config.num_cols();
        BB_ASSERT_EQ(u.size(), config.num_variables, "opening point size mismatch");

        if (claims.send_roots) {
            for (size_t g = 0; g < claims.groups.size(); ++g) {
                transcript->send_to_verifier(detail::sf_label("proot", g),
                                             Hasher::digest_to_fields(claims.groups[g]->tree.root()));
            }
        }
        const fr rho = transcript->template get_challenge<fr>("SWITCHFOLD:rho");

        const std::vector<fr> a = eq_tensor(u.subspan(0, config.log_num_cols));
        const std::vector<fr> b = eq_tensor(u.subspan(config.log_num_cols));
        const auto scalars = detail::chain_scalars(claims.unshifted, claims.to_be_shifted, rho);

        // Descent input M = [w_u | w_s | w_s2 | 0]; the segment index is the top two variables.
        // Threads split the *column* range, so each owns a disjoint slice of all three segments and
        // no accumulator merge is needed.
        std::vector<fr> message(4 * num_cols, fr::zero());
        {
            struct RowSource {
                const fr* dense;
                fr scalar_u;
                fr scalar_s;
            };
            std::vector<RowSource> sources;
            sources.reserve(scalars.size());
            for (const auto& [key, chain] : scalars) {
                sources.push_back(
                    { claims.groups[key.first]->coefficients[key.second].data(), chain.scalar_u, chain.scalar_s });
            }
            parallel_for_range(num_cols, [&](size_t start, size_t end) {
                for (const RowSource& source : sources) {
                    for (size_t r = 0; r < num_rows; ++r) {
                        const fr* row = source.dense + (r * num_cols);
                        const fr cu = source.scalar_u * b[r];
                        const fr cs = source.scalar_s * b[r];
                        const fr cs2 = (r > 0) ? source.scalar_s * b[r - 1] : fr::zero();
                        for (size_t c = start; c < end; ++c) {
                            message[c] += cu * row[c];
                            message[num_cols + c] += cs * row[c];
                            message[(2 * num_cols) + c] += cs2 * row[c];
                        }
                    }
                }
            });
        }

        fr expected = fr::zero();
        {
            fr rho_power = fr::one();
            for (const fr& value : claims.unshifted_evaluations) {
                expected += rho_power * value;
                rho_power *= rho;
            }
            for (const fr& value : claims.shifted_evaluations) {
                expected += rho_power * value;
                rho_power *= rho;
            }
        }

        std::vector<detail::Claim> pending;
        pending.push_back(initial_claim(config, a, u, expected));

        // Payload queries: the openings travel, the derived codeword values do not.
        const size_t codeword_bits = config.log_num_cols + config.log_inv_rate;
        const fr payload_root = fr::get_root_of_unity(codeword_bits);
        for (size_t q = 0; q < config.num_queries; ++q) {
            const fr challenge = transcript->template get_challenge<fr>(detail::sf_label("pquery", 0, q));
            const size_t index = detail::index_from_challenge(challenge, codeword_bits);
            std::array<fr, 3> derived{ fr::zero(), fr::zero(), fr::zero() };
            for (size_t g = 0; g < claims.groups.size(); ++g) {
                const auto opening = claims.groups[g]->tree.open(index);
                send_opening(transcript, opening);
                accumulate_payload_values(scalars, b, g, num_rows, opening.values, derived);
            }
            append_query_claims(config, payload_root, index, derived, pending);
        }

        run_descent(config, std::move(message), std::move(pending), transcript);
    }

    /** @brief The evaluation claim, as eq on segment 0, shifted-eq on segment 1, e_0 on segment 2. */
    static detail::Claim initial_claim(const SwitchFoldConfig& config,
                                       const std::vector<fr>& a,
                                       std::span<const fr> u,
                                       const fr& expected)
    {
        const std::span<const fr> u_lo = u.subspan(0, config.log_num_cols);
        detail::Claim claim;
        claim.value = expected;
        WeightTerm unshifted = WeightTerm::eq_weight(u_lo, fr::one());
        detail::select_segment(unshifted, 0);
        claim.terms.push_back(unshifted);

        std::vector<WeightTerm> shifted;
        detail::append_shifted_eq_terms(u_lo, fr::one(), shifted);
        for (WeightTerm& term : shifted) {
            detail::select_segment(term, 1);
            claim.terms.push_back(term);
        }

        WeightTerm carry =
            WeightTerm::eq_weight(std::vector<fr>(config.log_num_cols, fr::zero()), a[config.num_cols() - 1]);
        detail::select_segment(carry, 2);
        claim.terms.push_back(carry);
        return claim;
    }

    /** @brief Fold one group's payload leaf into the three combined-row codeword values. */
    static void accumulate_payload_values(const std::map<std::pair<size_t, size_t>, detail::ChainScalars>& scalars,
                                          const std::vector<fr>& b,
                                          size_t group,
                                          size_t num_rows,
                                          const std::vector<fr>& leaf,
                                          std::array<fr, 3>& derived)
    {
        for (const auto& [key, chain] : scalars) {
            if (key.first != group) {
                continue;
            }
            for (size_t r = 0; r < num_rows; ++r) {
                const fr& value = leaf[(key.second * num_rows) + r];
                derived[0] += chain.scalar_u * b[r] * value;
                derived[1] += chain.scalar_s * b[r] * value;
                if (r > 0) {
                    derived[2] += chain.scalar_s * b[r - 1] * value;
                }
            }
        }
    }

    /** @brief The three per-segment code-switching claims of one payload query. */
    static void append_query_claims(const SwitchFoldConfig& config,
                                    const fr& payload_root,
                                    size_t index,
                                    const std::array<fr, 3>& derived,
                                    std::vector<detail::Claim>& pending)
    {
        const fr point = payload_root.pow(static_cast<uint64_t>(index));
        for (size_t segment = 0; segment < 3; ++segment) {
            WeightTerm term = WeightTerm::pow_weight(point, config.log_num_cols, fr::one());
            detail::select_segment(term, segment);
            pending.push_back({ { term }, derived.at(segment) });
        }
    }

    template <typename Transcript>
    static void send_opening(const std::shared_ptr<Transcript>& transcript, const typename Tree::Opening& opening)
    {
        std::vector<fr> flat;
        flat.reserve(opening.values.size() + (opening.path.size() * Hasher::DIGEST_NUM_FIELDS));
        flat.insert(flat.end(), opening.values.begin(), opening.values.end());
        for (const auto& digest : opening.path) {
            const auto fields = Hasher::digest_to_fields(digest);
            flat.insert(flat.end(), fields.begin(), fields.end());
        }
        transcript->send_unhashed_to_verifier(flat);
    }

  private:
    template <typename Transcript>
    static void run_descent(const SwitchFoldConfig& config,
                            std::vector<fr> message,
                            std::vector<detail::Claim> pending,
                            const std::shared_ptr<Transcript>& transcript)
    {
        RSDomains domains;
        const size_t switch_factor = config.switch_factor();
        const size_t base = size_t(1) << config.log_base_length;
        size_t level = 0;

        while (message.size() > base) {
            const size_t sub_length = message.size() / switch_factor;
            const size_t codeword_length = sub_length << config.log_inv_rate;
            const auto& domain = domains.get(codeword_length);
            std::vector<std::vector<fr>> codewords(switch_factor);
            for (size_t c = 0; c < switch_factor; ++c) {
                const std::span<const fr> column(message.data() + (c * sub_length), sub_length);
                codewords[c] = rs_encode(column, domain, domains.round_roots());
            }
            Tree tree(std::move(codewords), /*log_arity=*/0, /*salted=*/false);
            transcript->send_to_verifier(detail::sf_label("root", level), Hasher::digest_to_fields(tree.root()));

            const size_t index_bits = numeric::get_msb(codeword_length);
            std::vector<size_t> indices;
            std::vector<std::vector<fr>> leaves;
            for (size_t q = 0; q < config.num_queries; ++q) {
                const fr challenge = transcript->template get_challenge<fr>(detail::sf_label("query", level, q));
                const size_t index = detail::index_from_challenge(challenge, index_bits);
                const auto opening = tree.open(index);
                send_opening(transcript, opening);
                indices.push_back(index);
                leaves.push_back(opening.values);
            }

            const fr mu = transcript->template get_challenge<fr>(detail::sf_label("mu", level));
            auto [terms, value] = detail::batch_claims(pending, mu);
            std::vector<fr> table = detail::weight_table(terms, message.size());

            std::vector<fr> current = message;
            fr claim = value;
            std::vector<fr> challenges;
            while (current.size() > 1) {
                const auto [h0, h2] = detail::inner_product_round(current, table);
                transcript->send_to_verifier(detail::sf_label("h0", level, challenges.size()), h0);
                transcript->send_to_verifier(detail::sf_label("h2", level, challenges.size()), h2);
                const fr alpha =
                    transcript->template get_challenge<fr>(detail::sf_label("sc", level, challenges.size()));
                claim = detail::interpolate_quadratic(h0, claim - h0, h2, alpha);
                whir::fold_array_in_place(current, alpha);
                whir::fold_array_in_place(table, alpha);
                challenges.push_back(alpha);
            }
            transcript->send_to_verifier(detail::sf_label("m_eval", level), current[0]);

            // The tensor claim splits: low challenges bind the within-column index, high ones give
            // the column combination that becomes the next level's message.
            const size_t log_sub = numeric::get_msb(sub_length);
            const std::vector<fr> combination = eq_tensor(std::span<const fr>(challenges).subspan(log_sub));
            std::vector<fr> next(sub_length, fr::zero());
            for (size_t c = 0; c < switch_factor; ++c) {
                for (size_t i = 0; i < sub_length; ++i) {
                    next[i] += combination[c] * message[(c * sub_length) + i];
                }
            }

            pending.clear();
            pending.push_back(
                { { WeightTerm::eq_weight(std::span<const fr>(challenges).subspan(0, log_sub), fr::one()) },
                  current[0] });
            for (size_t q = 0; q < indices.size(); ++q) {
                fr derived = fr::zero();
                for (size_t c = 0; c < switch_factor; ++c) {
                    derived += combination[c] * leaves[q][c];
                }
                const fr point = domain.root.pow(static_cast<uint64_t>(indices[q]));
                pending.push_back({ { WeightTerm::pow_weight(point, log_sub, fr::one()) }, derived });
            }

            message = std::move(next);
            ++level;
        }

        for (size_t i = 0; i < message.size(); ++i) {
            transcript->send_to_verifier(detail::sf_label("base", level, i), message[i]);
        }
    }
};

/** @brief The mirror of `SwitchFoldProver`: the same claim bookkeeping, driven by the transcript. */
template <typename Hasher> class SwitchFoldVerifier {
  public:
    using Tree = MerkleTree<Hasher>;
    using Digest = typename Hasher::Digest;
    using Prover = SwitchFoldProver<Hasher>;

    struct Claims {
        std::vector<size_t> group_num_columns;
        std::vector<Digest> group_roots; // when non-empty, transcript-bound
        std::vector<SwitchFoldColumnRef> unshifted;
        std::vector<fr> unshifted_evaluations;
        std::vector<SwitchFoldColumnRef> to_be_shifted;
        std::vector<fr> shifted_evaluations;
    };

    template <typename Transcript>
    static bool verify(const SwitchFoldConfig& config,
                       const Claims& claims,
                       std::span<const fr> u,
                       const std::shared_ptr<Transcript>& transcript)
    {
        const size_t num_rows = config.num_rows();
        BB_ASSERT_EQ(u.size(), config.num_variables, "opening point size mismatch");

        std::vector<Digest> roots = claims.group_roots;
        if (roots.empty()) {
            for (size_t g = 0; g < claims.group_num_columns.size(); ++g) {
                roots.push_back(Hasher::digest_from_fields(
                    transcript->template receive_from_prover<std::array<fr, Hasher::DIGEST_NUM_FIELDS>>(
                        detail::sf_label("proot", g))));
            }
        }
        const fr rho = transcript->template get_challenge<fr>("SWITCHFOLD:rho");

        const std::vector<fr> a = eq_tensor(u.subspan(0, config.log_num_cols));
        const std::vector<fr> b = eq_tensor(u.subspan(config.log_num_cols));
        const auto scalars = detail::chain_scalars(claims.unshifted, claims.to_be_shifted, rho);

        fr expected = fr::zero();
        {
            fr rho_power = fr::one();
            for (const fr& value : claims.unshifted_evaluations) {
                expected += rho_power * value;
                rho_power *= rho;
            }
            for (const fr& value : claims.shifted_evaluations) {
                expected += rho_power * value;
                rho_power *= rho;
            }
        }

        std::vector<detail::Claim> pending;
        pending.push_back(Prover::initial_claim(config, a, u, expected));

        const size_t codeword_bits = config.log_num_cols + config.log_inv_rate;
        const fr payload_root = fr::get_root_of_unity(codeword_bits);
        for (size_t q = 0; q < config.num_queries; ++q) {
            const fr challenge = transcript->template get_challenge<fr>(detail::sf_label("pquery", 0, q));
            const size_t index = detail::index_from_challenge(challenge, codeword_bits);
            std::array<fr, 3> derived{ fr::zero(), fr::zero(), fr::zero() };
            for (size_t g = 0; g < claims.group_num_columns.size(); ++g) {
                const size_t leaf_width = claims.group_num_columns[g] * num_rows;
                typename Tree::Opening opening;
                if (!read_opening(transcript, leaf_width, codeword_bits, index, roots[g], opening)) {
                    return false;
                }
                Prover::accumulate_payload_values(scalars, b, g, num_rows, opening.values, derived);
            }
            Prover::append_query_claims(config, payload_root, index, derived, pending);
        }

        return check_descent(config, std::move(pending), transcript);
    }

  private:
    template <typename Transcript>
    static bool read_opening(const std::shared_ptr<Transcript>& transcript,
                             size_t leaf_width,
                             size_t index_bits,
                             size_t index,
                             const Digest& root,
                             typename Tree::Opening& opening)
    {
        const std::vector<fr> flat =
            transcript->receive_unhashed_from_prover(leaf_width + (index_bits * Hasher::DIGEST_NUM_FIELDS));
        opening.values.assign(flat.begin(), flat.begin() + static_cast<std::ptrdiff_t>(leaf_width));
        opening.path.clear();
        for (size_t level = 0; level < index_bits; ++level) {
            const size_t offset = leaf_width + (level * Hasher::DIGEST_NUM_FIELDS);
            opening.path.push_back(
                Hasher::digest_from_fields(std::span<const fr>(flat).subspan(offset, Hasher::DIGEST_NUM_FIELDS)));
        }
        return Tree::verify(root, index, opening);
    }

    template <typename Transcript>
    static bool check_descent(const SwitchFoldConfig& config,
                              std::vector<detail::Claim> pending,
                              const std::shared_ptr<Transcript>& transcript)
    {
        const size_t switch_factor = config.switch_factor();
        const size_t base = size_t(1) << config.log_base_length;
        size_t length = 4 * config.num_cols();
        size_t level = 0;

        while (length > base) {
            const size_t sub_length = length / switch_factor;
            const size_t codeword_length = sub_length << config.log_inv_rate;
            const size_t index_bits = numeric::get_msb(codeword_length);
            const fr domain_root = fr::get_root_of_unity(index_bits);

            const Digest root = Hasher::digest_from_fields(
                transcript->template receive_from_prover<std::array<fr, Hasher::DIGEST_NUM_FIELDS>>(
                    detail::sf_label("root", level)));

            std::vector<size_t> indices;
            std::vector<std::vector<fr>> leaves;
            for (size_t q = 0; q < config.num_queries; ++q) {
                const fr challenge = transcript->template get_challenge<fr>(detail::sf_label("query", level, q));
                const size_t index = detail::index_from_challenge(challenge, index_bits);
                typename Tree::Opening opening;
                if (!read_opening(transcript, switch_factor, index_bits, index, root, opening)) {
                    return false;
                }
                indices.push_back(index);
                leaves.push_back(opening.values);
            }

            const fr mu = transcript->template get_challenge<fr>(detail::sf_label("mu", level));
            auto [terms, value] = detail::batch_claims(pending, mu);

            fr claim = value;
            std::vector<fr> challenges;
            const size_t num_variables = numeric::get_msb(length);
            for (size_t round = 0; round < num_variables; ++round) {
                const fr h0 = transcript->template receive_from_prover<fr>(detail::sf_label("h0", level, round));
                const fr h2 = transcript->template receive_from_prover<fr>(detail::sf_label("h2", level, round));
                const fr alpha = transcript->template get_challenge<fr>(detail::sf_label("sc", level, round));
                claim = detail::interpolate_quadratic(h0, claim - h0, h2, alpha);
                challenges.push_back(alpha);
            }
            const fr m_eval = transcript->template receive_from_prover<fr>(detail::sf_label("m_eval", level));
            if (m_eval * detail::weight_value(terms, challenges) != claim) {
                return false;
            }

            const size_t log_sub = numeric::get_msb(sub_length);
            const std::vector<fr> combination = eq_tensor(std::span<const fr>(challenges).subspan(log_sub));
            pending.clear();
            pending.push_back(
                { { WeightTerm::eq_weight(std::span<const fr>(challenges).subspan(0, log_sub), fr::one()) }, m_eval });
            for (size_t q = 0; q < indices.size(); ++q) {
                fr derived = fr::zero();
                for (size_t c = 0; c < switch_factor; ++c) {
                    derived += combination[c] * leaves[q][c];
                }
                const fr point = domain_root.pow(static_cast<uint64_t>(indices[q]));
                pending.push_back({ { WeightTerm::pow_weight(point, log_sub, fr::one()) }, derived });
            }

            length = sub_length;
            ++level;
        }

        std::vector<fr> message(length);
        for (size_t i = 0; i < length; ++i) {
            message[i] = transcript->template receive_from_prover<fr>(detail::sf_label("base", level, i));
        }
        for (const detail::Claim& claim : pending) {
            const std::vector<fr> table = detail::weight_table(claim.terms, length);
            fr total = fr::zero();
            for (size_t i = 0; i < length; ++i) {
                total += message[i] * table[i];
            }
            if (total != claim.value) {
                return false;
            }
        }
        return true;
    }
};

} // namespace bb::switchfold
