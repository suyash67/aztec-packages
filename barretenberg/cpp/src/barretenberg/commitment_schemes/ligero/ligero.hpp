#pragma once

#include "barretenberg/commitment_schemes/whir/merkle_tree.hpp"
#include "barretenberg/commitment_schemes/whir/rs_code.hpp"
#include "barretenberg/commitment_schemes/whir/weights.hpp"
#include "barretenberg/common/assert.hpp"
#include "barretenberg/common/thread.hpp"
#include "barretenberg/numeric/bitop/get_msb.hpp"
#include "barretenberg/polynomials/polynomial.hpp"

#include <cmath>
#include <map>
#include <mutex>
#include <span>
#include <string>
#include <utility>
#include <vector>

namespace bb::ligero {

using whir::MerkleTree;
using whir::rs_encode;
using whir::RSDomains;

/** @brief Reference to one polynomial ("column" of the claim batch) of one committed group. */
struct LigeroColumnRef {
    size_t group;
    size_t column;
};

/**
 * @brief Ligero parameters (README.md §3). Both parties derive the identical values from the same
 * inputs; nothing is sent through the transcript.
 */
struct LigeroConfig {
    size_t num_variables; // m: polynomials have 2^m evaluations
    size_t log_num_cols;  // matrix is 2^{m - log_num_cols} rows x 2^{log_num_cols} cols
    size_t log_inv_rate;  // rows encode to 2^{log_num_cols + log_inv_rate} entries
    size_t security_bits; // λ
    size_t num_queries;   // t = ceil(λ / log_inv_rate), up-to-capacity conjecture

    size_t num_rows() const { return size_t(1) << (num_variables - log_num_cols); }
    size_t num_cols() const { return size_t(1) << log_num_cols; }
    size_t codeword_length() const { return size_t(1) << (log_num_cols + log_inv_rate); }

    /** @brief Proof-size-minimizing column count: C* ~ sqrt(t * total_polys * 2^{m-1}). */
    static size_t default_log_num_cols(size_t num_variables, size_t est_total_polynomials, size_t num_queries)
    {
        const double c_star = std::sqrt(static_cast<double>(num_queries * est_total_polynomials) *
                                        std::pow(2.0, double(num_variables) - 1));
        size_t log_cols = static_cast<size_t>(std::lround(std::log2(c_star)));
        log_cols = std::min(log_cols, num_variables - 1);
        return std::max(log_cols, size_t(1));
    }

    static LigeroConfig create(size_t num_variables,
                               size_t security_bits = 100,
                               size_t log_inv_rate = 2,
                               size_t est_total_polynomials = 36)
    {
        BB_ASSERT_GT(num_variables, size_t(1));
        const size_t num_queries = (security_bits + log_inv_rate - 1) / log_inv_rate;
        const size_t log_num_cols = default_log_num_cols(num_variables, est_total_polynomials, num_queries);
        BB_ASSERT_LTE(log_num_cols + log_inv_rate, size_t(28), "row codeword exceeds field 2-adicity");
        return { num_variables, log_num_cols, log_inv_rate, security_bits, num_queries };
    }
};

/**
 * @brief Prover-side commitment to a group of polynomials: dense evaluation arrays plus the Merkle
 * tree over the interleaved row codewords (leaf j = column j of every row).
 */
template <typename Hasher> struct LigeroGroupData {
    std::vector<std::vector<fr>> coefficients; // dense 2^m arrays, one per polynomial
    MerkleTree<Hasher> tree;

    size_t num_columns() const { return coefficients.size(); }
};

namespace detail {

using whir::eq_tensor;

inline std::string ligero_label(const std::string& name, size_t i)
{
    return "LIGERO:" + name + "_" + std::to_string(i);
}

inline size_t index_from_challenge(const fr& challenge, size_t index_bits)
{
    return static_cast<size_t>(uint256_t(challenge).data[0] & ((uint64_t(1) << index_bits) - 1));
}

/**
 * @brief A column index in [0, bound), for codes whose length is not a power of two.
 * @details Reed-Solomon codewords live on a power-of-two subgroup so a bitmask suffices, but the
 * linear-time codes have arbitrary lengths; the modulo bias is below 2^-40 at these sizes.
 */
inline size_t index_from_challenge_bounded(const fr& challenge, size_t bound)
{
    return static_cast<size_t>(uint256_t(challenge).data[0] % bound);
}

inline size_t next_power_of_two(size_t value)
{
    size_t power = 1;
    while (power < value) {
        power <<= 1;
    }
    return power;
}

} // namespace detail

/**
 * @brief One stretch of codeword positions and how many of them to open.
 *
 * @details Most codes get a single segment covering the whole codeword. Codes with a *piecewise*
 * distance guarantee — where a corrupted word must differ substantially inside one specific stretch,
 * rather than merely somewhere overall — are tested stretch by stretch instead, each with the query
 * count its own distance warrants. That is strictly stronger than sampling uniformly against the
 * diluted overall distance: see `brakedown/` for the single-segment case and `bolt/` for the
 * two-segment one.
 */
struct QuerySegment {
    size_t begin;
    size_t length;
    size_t num_queries;
};

/**
 * @brief The Reed-Solomon code policy: what Ligero has always used.
 *
 * @details A code policy supplies the encoder and the query count. Only two things about the code
 * enter the protocol - `encode`, and how many columns must be opened - so parameterizing over this
 * is enough to run the same tensor PCS over a linear-time code (see `brakedown/`). The tree is
 * built over `padded_codeword_length()` leaves because `MerkleTree` requires a power of two, while
 * queries are drawn below `codeword_length()`, so padding is never opened and never enters
 * soundness.
 */
struct RSCodePolicy {
    explicit RSCodePolicy(const LigeroConfig& config)
        : config(config)
    {
        // `encode` is called from `commit_group`'s worker threads, and `RSDomains::get` inserts into
        // a map. Build the one domain this code needs up front so that `encode` is read-only and
        // safe to share; without this the concurrent inserts corrupt the map.
        domains.get(config.codeword_length());
    }

    size_t message_length() const { return config.num_cols(); }
    size_t codeword_length() const { return config.codeword_length(); }
    size_t padded_codeword_length() const { return config.codeword_length(); }

    std::vector<fr> encode(std::span<const fr> row) const
    {
        return rs_encode(row, domains.get(config.codeword_length()), domains.round_roots());
    }

    /** @brief Under the capacity conjecture the per-query error is the rate, giving t = ceil(λ/r). */
    static size_t num_queries(const LigeroConfig& config) { return config.num_queries; }

    /** @brief One segment spanning the whole codeword: RS has no piecewise structure to exploit. */
    std::vector<QuerySegment> query_plan(const LigeroConfig& config) const
    {
        return { { 0, codeword_length(), num_queries(config) } };
    }

    LigeroConfig config;
    mutable RSDomains domains;
};

/**
 * @brief Transparent commitment key: parameters plus the row code.
 * @tparam Code the code policy; defaults to Reed-Solomon, which is what `LigeroHonk` uses.
 */
template <typename Hasher, typename Code = RSCodePolicy> class LigeroCommitmentKey {
  public:
    using CodePolicy = Code;

    explicit LigeroCommitmentKey(const LigeroConfig& config)
        : config(config)
        , code(config)
    {}

    template <typename CodeConfig>
    LigeroCommitmentKey(const LigeroConfig& config, const CodeConfig& code_config)
        : config(config)
        , code(code_config)
    {}

    /**
     * @brief Commit a group of polynomials (each <= 2^m evaluations, zero-padded) into one tree.
     * @param to_be_shifted accepted for interface parity with WHIR; Ligero's shifted claims need no
     * commitment-side preparation (README.md §2)
     */
    LigeroGroupData<Hasher> commit_group(std::vector<std::vector<fr>> payload_columns,
                                         const std::vector<bool>& to_be_shifted = {}) const
    {
        static_cast<void>(to_be_shifted);
        const size_t n = size_t(1) << config.num_variables;
        const size_t num_rows = config.num_rows();
        const size_t num_cols = config.num_cols();
        const size_t padded = code.padded_codeword_length();

        std::vector<std::vector<fr>> dense;
        dense.reserve(payload_columns.size());
        for (auto& column : payload_columns) {
            BB_ASSERT_LTE(column.size(), n, "polynomial too large for the configured size");
            column.resize(n, fr::zero());
            dense.push_back(std::move(column));
        }

        std::vector<std::vector<fr>> row_codewords(dense.size() * num_rows);
        parallel_for_range(row_codewords.size(), [&](size_t start, size_t end) {
            for (size_t i = start; i < end; ++i) {
                const std::span<const fr> row(dense[i / num_rows].data() + ((i % num_rows) * num_cols), num_cols);
                row_codewords[i] = code.encode(row);
                // The tree needs a power-of-two leaf count; padded positions are never queried.
                row_codewords[i].resize(padded, fr::zero());
            }
        });
        MerkleTree<Hasher> tree(std::move(row_codewords), /*log_arity=*/0, /*salted=*/false);
        return { std::move(dense), std::move(tree) };
    }

    /** @brief Commit a group of Honk polynomials (virtual zeros outside their spans honored). */
    LigeroGroupData<Hasher> commit_group(std::span<const Polynomial<fr>* const> polynomials,
                                         const std::vector<bool>& to_be_shifted = {}) const
    {
        std::vector<std::vector<fr>> payload_columns;
        payload_columns.reserve(polynomials.size());
        for (const Polynomial<fr>* polynomial : polynomials) {
            std::vector<fr> column(polynomial->end_index(), fr::zero());
            for (size_t i = polynomial->start_index(); i < polynomial->end_index(); ++i) {
                column[i] = (*polynomial)[i];
            }
            payload_columns.push_back(std::move(column));
        }
        return commit_group(std::move(payload_columns), to_be_shifted);
    }

    LigeroConfig config;
    Code code;
};

/**
 * @brief Batched Ligero opening prover (README.md §1-§2). Combined rows travel unhashed with a
 * hashed digest for Fiat-Shamir; column openings are Merkle-bound to pre-challenge roots.
 */
template <typename Hasher> class LigeroProver {
  public:
    struct Claims {
        std::vector<const LigeroGroupData<Hasher>*> groups;
        std::vector<LigeroColumnRef> unshifted;
        std::vector<fr> unshifted_evaluations;
        std::vector<LigeroColumnRef> to_be_shifted;
        std::vector<fr> shifted_evaluations;
        bool send_roots = true;
    };

    /** @tparam CommitmentKey any `LigeroCommitmentKey` instantiation; the code enters only through
     * its codeword length and query count. */
    template <typename CommitmentKey, typename Transcript>
    static void prove(const CommitmentKey& ck,
                      const Claims& claims,
                      std::span<const fr> u,
                      const std::shared_ptr<Transcript>& transcript)
    {
        const LigeroConfig& config = ck.config;
        const size_t num_rows = config.num_rows();
        const size_t num_cols = config.num_cols();
        BB_ASSERT_EQ(u.size(), config.num_variables, "opening point size mismatch");

        if (claims.send_roots) {
            for (size_t g = 0; g < claims.groups.size(); ++g) {
                transcript->send_to_verifier(detail::ligero_label("root", g),
                                             Hasher::digest_to_fields(claims.groups[g]->tree.root()));
            }
        }
        const fr rho = transcript->template get_challenge<fr>("LIGERO:rho");

        const std::vector<fr> b = detail::eq_tensor(u.subspan(config.log_num_cols));

        // Per-(group, column) scalars of the three batching vectors (README.md §1).
        struct RowTask {
            const std::vector<fr>* dense;
            fr scalar_u;
            fr scalar_s;
        };
        std::map<std::pair<size_t, size_t>, RowTask> tasks;
        fr rho_power = fr::one();
        for (const LigeroColumnRef& ref : claims.unshifted) {
            auto& task = tasks[{ ref.group, ref.column }];
            task.dense = &claims.groups[ref.group]->coefficients[ref.column];
            task.scalar_u += rho_power;
            rho_power *= rho;
        }
        for (const LigeroColumnRef& ref : claims.to_be_shifted) {
            auto& task = tasks[{ ref.group, ref.column }];
            task.dense = &claims.groups[ref.group]->coefficients[ref.column];
            task.scalar_s += rho_power;
            rho_power *= rho;
        }

        // w_u = B_u^T M, w_s = B_s^T M, w_s2 = B_s2^T M in one pass over the data.
        std::vector<fr> w_u(num_cols, fr::zero());
        std::vector<fr> w_s(num_cols, fr::zero());
        std::vector<fr> w_s2(num_cols, fr::zero());
        std::vector<RowTask> task_list;
        task_list.reserve(tasks.size());
        for (const auto& [key, task] : tasks) {
            task_list.push_back(task);
        }
        std::mutex w_mutex;
        parallel_for_range(task_list.size() * num_rows, [&](size_t start, size_t end) {
            std::vector<fr> local_u(num_cols, fr::zero());
            std::vector<fr> local_s(num_cols, fr::zero());
            std::vector<fr> local_s2(num_cols, fr::zero());
            for (size_t i = start; i < end; ++i) {
                const RowTask& task = task_list[i / num_rows];
                const size_t r = i % num_rows;
                const fr* row = task.dense->data() + r * num_cols;
                const fr cu = task.scalar_u * b[r];
                const fr cs = task.scalar_s * b[r];
                const fr cs2 = (r > 0) ? task.scalar_s * b[r - 1] : fr::zero();
                for (size_t c = 0; c < num_cols; ++c) {
                    local_u[c] += cu * row[c];
                    local_s[c] += cs * row[c];
                    local_s2[c] += cs2 * row[c];
                }
            }
            std::scoped_lock lock(w_mutex);
            for (size_t c = 0; c < num_cols; ++c) {
                w_u[c] += local_u[c];
                w_s[c] += local_s[c];
                w_s2[c] += local_s2[c];
            }
        });

        // Stream the rows unhashed; absorb one digest of them for Fiat-Shamir.
        std::vector<fr> combined;
        combined.reserve(3 * num_cols);
        combined.insert(combined.end(), w_u.begin(), w_u.end());
        combined.insert(combined.end(), w_s.begin(), w_s.end());
        combined.insert(combined.end(), w_s2.begin(), w_s2.end());
        transcript->send_unhashed_to_verifier(combined);
        transcript->send_to_verifier(std::string("LIGERO:w_digest"),
                                     Hasher::digest_to_fields(Hasher::hash_leaf(combined, std::nullopt)));

        // Column queries: open every group's tree at each sampled leaf, segment by segment.
        size_t query = 0;
        for (const QuerySegment& segment : ck.code.query_plan(config)) {
            for (size_t s = 0; s < segment.num_queries; ++s, ++query) {
                const fr challenge = transcript->template get_challenge<fr>(detail::ligero_label("col", query));
                const size_t idx = segment.begin + detail::index_from_challenge_bounded(challenge, segment.length);
                for (const LigeroGroupData<Hasher>* group : claims.groups) {
                    const auto opening = group->tree.open(idx);
                    std::vector<fr> flat = opening.values;
                    for (const auto& digest : opening.path) {
                        const auto fields = Hasher::digest_to_fields(digest);
                        flat.insert(flat.end(), fields.begin(), fields.end());
                    }
                    transcript->send_unhashed_to_verifier(flat);
                }
            }
        }
    }
};

/** @brief Batched Ligero opening verifier; the mirror of `LigeroProver`. */
template <typename Hasher> class LigeroVerifier {
  public:
    using Digest = typename Hasher::Digest;
    using Tree = MerkleTree<Hasher>;

    struct Claims {
        std::vector<size_t> group_num_columns; // polynomials per group
        std::vector<LigeroColumnRef> unshifted;
        std::vector<fr> unshifted_evaluations;
        std::vector<LigeroColumnRef> to_be_shifted;
        std::vector<fr> shifted_evaluations;
        std::vector<Digest> group_roots; // when non-empty, roots are already transcript-bound
    };

    template <typename Code, typename Transcript>
    static bool verify(const LigeroConfig& config,
                       const Claims& claims,
                       std::span<const fr> u,
                       const std::shared_ptr<Transcript>& transcript,
                       const Code& code)
    {
        const size_t num_rows = config.num_rows();
        const size_t num_cols = config.num_cols();
        BB_ASSERT_EQ(u.size(), config.num_variables, "opening point size mismatch");

        std::vector<Digest> roots = claims.group_roots;
        if (roots.empty()) {
            for (size_t g = 0; g < claims.group_num_columns.size(); ++g) {
                const auto fields = transcript->template receive_from_prover<std::array<fr, Hasher::DIGEST_NUM_FIELDS>>(
                    detail::ligero_label("root", g));
                roots.push_back(Hasher::digest_from_fields(fields));
            }
        }
        const fr rho = transcript->template get_challenge<fr>("LIGERO:rho");

        const std::vector<fr> a = detail::eq_tensor(u.subspan(0, config.log_num_cols));
        const std::vector<fr> b = detail::eq_tensor(u.subspan(config.log_num_cols));

        // Batched claim value and per-(group,column) scalars.
        std::map<std::pair<size_t, size_t>, std::pair<fr, fr>> scalars; // (scalar_u, scalar_s)
        fr expected = fr::zero();
        fr rho_power = fr::one();
        for (size_t i = 0; i < claims.unshifted.size(); ++i) {
            scalars[{ claims.unshifted[i].group, claims.unshifted[i].column }].first += rho_power;
            expected += rho_power * claims.unshifted_evaluations[i];
            rho_power *= rho;
        }
        for (size_t l = 0; l < claims.to_be_shifted.size(); ++l) {
            scalars[{ claims.to_be_shifted[l].group, claims.to_be_shifted[l].column }].second += rho_power;
            expected += rho_power * claims.shifted_evaluations[l];
            rho_power *= rho;
        }

        // Combined rows and their Fiat-Shamir digest.
        const std::vector<fr> combined = transcript->receive_unhashed_from_prover(3 * num_cols);
        const auto digest_fields = transcript->template receive_from_prover<std::array<fr, Hasher::DIGEST_NUM_FIELDS>>(
            std::string("LIGERO:w_digest"));
        if (Hasher::digest_from_fields(digest_fields) != Hasher::hash_leaf(combined, std::nullopt)) {
            return false;
        }
        const std::span<const fr> w_u(combined.data(), num_cols);
        const std::span<const fr> w_s(combined.data() + num_cols, num_cols);
        const std::span<const fr> w_s2(combined.data() + 2 * num_cols, num_cols);

        // Claim equation (README.md §1): rank-2 shift decomposition.
        fr total = fr::zero();
        for (size_t c = 0; c < num_cols; ++c) {
            total += w_u[c] * a[c];
            if (c >= 1) {
                total += w_s[c] * a[c - 1];
            }
        }
        total += a[num_cols - 1] * w_s2[0];
        if (total != expected) {
            return false;
        }

        // Encodings of the combined rows, for the per-column consistency checks.
        const std::vector<fr> enc_u = code.encode(w_u);
        const std::vector<fr> enc_s = code.encode(w_s);
        const std::vector<fr> enc_s2 = code.encode(w_s2);

        // Column checks: authenticate each opened column and test the three linear combinations.
        // The path length follows the padded tree; queries stay inside the true codeword.
        const size_t index_bits = numeric::get_msb(code.padded_codeword_length());
        std::vector<std::pair<size_t, size_t>> plan; // (begin, length) per query, in transcript order
        for (const QuerySegment& segment : code.query_plan(config)) {
            for (size_t s = 0; s < segment.num_queries; ++s) {
                plan.emplace_back(segment.begin, segment.length);
            }
        }
        for (size_t s = 0; s < plan.size(); ++s) {
            const fr challenge = transcript->template get_challenge<fr>(detail::ligero_label("col", s));
            const size_t idx = plan[s].first + detail::index_from_challenge_bounded(challenge, plan[s].second);
            fr acc_u = fr::zero();
            fr acc_s = fr::zero();
            fr acc_s2 = fr::zero();
            for (size_t g = 0; g < claims.group_num_columns.size(); ++g) {
                const size_t group_rows = claims.group_num_columns[g] * num_rows;
                const std::vector<fr> flat =
                    transcript->receive_unhashed_from_prover(group_rows + index_bits * Hasher::DIGEST_NUM_FIELDS);
                typename Tree::Opening opening;
                opening.values.assign(flat.begin(), flat.begin() + static_cast<std::ptrdiff_t>(group_rows));
                size_t cursor = group_rows;
                for (size_t l = 0; l < index_bits; ++l) {
                    opening.path.push_back(Hasher::digest_from_fields(
                        std::span<const fr>(flat).subspan(cursor, Hasher::DIGEST_NUM_FIELDS)));
                    cursor += Hasher::DIGEST_NUM_FIELDS;
                }
                if (!Tree::verify(roots[g], idx, opening)) {
                    return false;
                }
                for (size_t p = 0; p < claims.group_num_columns[g]; ++p) {
                    const auto it = scalars.find({ g, p });
                    if (it == scalars.end()) {
                        continue;
                    }
                    const auto& [scalar_u, scalar_s] = it->second;
                    for (size_t r = 0; r < num_rows; ++r) {
                        const fr& value = opening.values[p * num_rows + r];
                        acc_u += scalar_u * b[r] * value;
                        acc_s += scalar_s * b[r] * value;
                        if (r > 0) {
                            acc_s2 += scalar_s * b[r - 1] * value;
                        }
                    }
                }
            }
            if (acc_u != enc_u[idx] || acc_s != enc_s[idx] || acc_s2 != enc_s2[idx]) {
                return false;
            }
        }
        return true;
    }
};

} // namespace bb::ligero
