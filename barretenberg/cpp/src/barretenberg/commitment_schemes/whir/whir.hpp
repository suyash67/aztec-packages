#pragma once

#include "barretenberg/commitment_schemes/whir/merkle_tree.hpp"
#include "barretenberg/commitment_schemes/whir/rs_code.hpp"
#include "barretenberg/commitment_schemes/whir/weights.hpp"
#include "barretenberg/commitment_schemes/whir/whir_config.hpp"
#include "barretenberg/common/assert.hpp"
#include "barretenberg/common/thread.hpp"
#include "barretenberg/polynomials/polynomial.hpp"
#include "barretenberg/polynomials/polynomial_arithmetic.hpp"

#include <bit>
#include <deque>
#include <mutex>
#include <optional>
#include <span>
#include <string>
#include <vector>

namespace bb::whir {

/**
 * @brief Prover-side result of committing a group of polynomials ("columns") into one shared Merkle
 * tree (README.md §4.1, §10). Two layouts, chosen by `WhirConfig::stack_columns`:
 *
 *  - *interleaved* (default): one codeword per column, leaf j holding every column's coset-j values.
 *    The batched oracle stays the width of a single column, so the prover is fast; a query costs
 *    `columns · 2^k` values per path.
 *  - *stacked*: the columns are concatenated into one taller array (column j at offset j·2^m),
 *    encoded as a single codeword with narrow single-column leaves. A query costs `2^k` values per
 *    path regardless of the column count, but the batched oracle widens to the group's full size.
 *
 * The commitment sent to the verifier is the single `tree.root()` either way.
 */
template <typename Hasher> struct WhirGroupData {
    // One array per column when interleaved; exactly one (the stacked array) when stacked.
    std::vector<std::vector<fr>> arrays;
    size_t columns = 0; // logical payload columns, before any power-of-two stack padding
    bool stacked = false;
    MerkleTree<Hasher> tree;

    size_t num_columns() const { return columns; }
    /** @brief Columns per Merkle leaf: 1 when stacked, else one slot per column. */
    size_t leaf_columns() const { return stacked ? 1 : columns; }
    size_t num_variables() const { return static_cast<size_t>(std::countr_zero(arrays[0].size())); }
};

/** @brief Reference to one column of one committed group, shared by prover and verifier claims. */
struct WhirColumnRef {
    size_t group;
    size_t column;
};

namespace detail {

inline size_t ceil_log2(size_t n)
{
    size_t bits = 0;
    while ((size_t(1) << bits) < n) {
        ++bits;
    }
    return bits;
}

} // namespace detail

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
     * @brief Commit a group of payload columns (each of length <= 2^m) into one tree, in whichever
     * layout `config.stack_columns` selects (see `WhirGroupData`).
     * @details In zk mode the committed array(s) double: the payload occupies the low half and
     * `config.num_blinding_coefficients` fresh random coefficients occupy the high half (just above
     * the payload, offset by one more when the array is to-be-shifted so the shift contract's zero
     * slot is preserved), and the leaves are salted. See README.md §8.
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
        const size_t num_columns = payload_columns.size();
        for (const auto& column : payload_columns) {
            BB_ASSERT_LTE(column.size(), payload_size, "polynomial too large for the configured size");
        }

        if (!config.stack_columns) {
            std::vector<std::vector<fr>> dense;
            dense.reserve(num_columns);
            for (size_t c = 0; c < num_columns; ++c) {
                std::vector<fr> column = std::move(payload_columns[c]);
                column.resize(size_t(1) << config.num_variables, fr::zero());
                add_blinding(column, config.num_payload_variables, !to_be_shifted.empty() && to_be_shifted[c]);
                dense.push_back(std::move(column));
            }
            return commit_arrays(std::move(dense), num_columns, /*stacked=*/false);
        }

        const size_t stack_bits = detail::ceil_log2(num_columns);
        const size_t stacked_payload_bits = config.num_payload_variables + stack_bits;
        BB_ASSERT_LTE(stacked_payload_bits + (config.zk ? 1 : 0),
                      config.num_variables,
                      "group too wide for the configured stacking headroom");

        std::vector<fr> stacked(size_t(1) << (stacked_payload_bits + (config.zk ? 1 : 0)), fr::zero());
        bool any_shifted = false;
        for (size_t c = 0; c < num_columns; ++c) {
            std::copy(payload_columns[c].begin(),
                      payload_columns[c].end(),
                      stacked.begin() + static_cast<std::ptrdiff_t>(c * payload_size));
            any_shifted = any_shifted || (!to_be_shifted.empty() && to_be_shifted[c]);
        }
        if (any_shifted) {
            // The shifted virtual oracle is the stacked codeword scaled by x^{-1}, which represents
            // the coefficient-shifted array only when the constant coefficient vanishes.
            BB_ASSERT_EQ(stacked[0], fr::zero(), "to-be-shifted group must have a zero constant term");
        }
        add_blinding(stacked, stacked_payload_bits, any_shifted);
        std::vector<std::vector<fr>> arrays;
        arrays.push_back(std::move(stacked));
        return commit_arrays(std::move(arrays), num_columns, /*stacked=*/true);
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
        std::vector<std::vector<fr>> arrays;
        arrays.push_back(std::move(mask));
        return commit_arrays(std::move(arrays), 1, config.stack_columns);
    }

    WhirConfig config;
    mutable RSDomains domains;

  private:
    /**
     * @brief Encode the committed arrays on their own domain (rate `config.log_inv_rate`) and build
     * one tree over them. Arrays narrower than the protocol width live on a correspondingly smaller
     * domain, which only arises when stacking.
     */
    WhirGroupData<Hasher> commit_arrays(std::vector<std::vector<fr>> arrays, size_t num_columns, bool stacked) const
    {
        const size_t size = arrays.at(0).size();
        BB_ASSERT_EQ(size & (size - 1), size_t(0), "committed array must be a power of two");
        const auto& domain = domains.get(size << config.log_inv_rate);
        std::vector<std::vector<fr>> codewords;
        codewords.reserve(arrays.size());
        for (const auto& array : arrays) {
            BB_ASSERT_EQ(array.size(), size, "committed arrays must share one size");
            codewords.push_back(rs_encode(array, domain, domains.round_roots()));
        }
        MerkleTree<Hasher> tree(std::move(codewords), config.folding_factor_bits, /*salted=*/config.zk);
        return { std::move(arrays), num_columns, stacked, std::move(tree) };
    }

    /** @brief Fill the zk blinding coefficients just above a payload of `payload_bits` variables. */
    void add_blinding(std::vector<fr>& array, size_t payload_bits, bool to_be_shifted) const
    {
        if (!config.zk) {
            return;
        }
        const size_t offset = (size_t(1) << payload_bits) + (to_be_shifted ? 1 : 0);
        for (size_t i = 0; i < config.num_blinding_coefficients; ++i) {
            array[offset + i] = fr::random_element();
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

/** @brief MLE of an array (evaluation form, LSB variable first) at `point`; O(size). */
inline fr mle_of_span(std::span<const fr> values, std::span<const fr> point)
{
    BB_ASSERT_EQ(values.size(), size_t(1) << point.size(), "array/point size mismatch");
    std::vector<fr> current(values.begin(), values.end());
    for (const fr& x : point) {
        fold_array_in_place(current, x);
    }
    return current[0];
}

/** @brief eq(bits(index), tau) with bit 0 of `index` matched against tau[0]. */
inline fr eq_at_corner(size_t index, std::span<const fr> tau)
{
    fr eq = fr::one();
    for (size_t i = 0; i < tau.size(); ++i) {
        eq *= ((index >> i) & 1) == 1 ? tau[i] : fr::one() - tau[i];
    }
    return eq;
}

/**
 * @brief The deterministic reduction of per-column claims to per-constituent point claims on the
 * stacked arrays; derived identically by prover and verifier (README.md §10).
 *
 * Interleaved (the default): every claimed column is its own constituent, at stride 0 and claimed
 * at the single point u — so there is exactly one point, no cross evaluations, and one weight term,
 * recovering the plain ρ-batched opening.
 *
 * Stacked: each commitment group contributes up to two constituents of the batched round-0 oracle:
 * its stacked array S, and — when the group carries to-be-shifted claims — the coefficient-shifted
 * array T whose codeword is S's scaled by x^{-1}. A constituent of 2^s columns whose every stack
 * corner is claimed or known-zero carries one claim at the τ-reduced point (0^d, u, τ[0..s)); a
 * constituent with unclaimed non-zero corners carries one claim per claimed corner at
 * (0^d, u, bits(j)). The d leading zeros place the array's stride-2^d spread into the protocol's
 * 2^P-size oracle: spread(A)(x) = eq(0, x_{1..d})·A(x_{d+1..P}).
 */
struct StackedPlan {
    struct Claim {
        size_t point; // index into `points`
        fr value;     // the claimed evaluation of this constituent at that point
    };
    struct Constituent {
        size_t tree;        // round-0 tree index: the groups in order, then the zk mask
        size_t leaf_column; // slot within the tree's leaf: the column when interleaved, else 0
        bool shifted;
        size_t stride_bits; // d = P - (committed variables of the tree)
        std::vector<Claim> claims;
    };

    std::vector<std::vector<fr>> points; // full P-coordinate opening points, deduplicated
    std::vector<Constituent> constituents;

    size_t add_point(std::vector<fr> coords)
    {
        for (size_t t = 0; t < points.size(); ++t) {
            if (points[t] == coords) {
                return t;
            }
        }
        points.push_back(std::move(coords));
        return points.size() - 1;
    }

    bool is_diagonal(size_t point, size_t constituent) const
    {
        for (const Claim& claim : constituents[constituent].claims) {
            if (claim.point == point) {
                return true;
            }
        }
        return false;
    }
};

/**
 * @brief Build the claim plan from the shared claim shape. `mask_evaluation` appends the zk mask
 * as a final full-width constituent claimed at (u, τ, 0); pass its value (prover computes it,
 * verifier reads it from the proof stream).
 */
inline StackedPlan build_stacked_plan(const WhirConfig& config,
                                      std::span<const size_t> group_num_columns,
                                      std::span<const WhirColumnRef> unshifted,
                                      std::span<const fr> unshifted_evaluations,
                                      std::span<const WhirColumnRef> to_be_shifted,
                                      std::span<const fr> shifted_evaluations,
                                      std::span<const fr> u,
                                      std::span<const fr> taus)
{
    const size_t num_groups = group_num_columns.size();
    const size_t zk_bit = config.zk ? 1 : 0;
    const size_t total_vars = config.num_variables;
    const size_t payload_vars = config.num_payload_variables;
    BB_ASSERT_EQ(taus.size(), total_vars - payload_vars - zk_bit, "tau challenge count mismatch");

    StackedPlan plan;
    auto make_point = [&](size_t stride_bits, std::span<const fr> column_coords) {
        std::vector<fr> coords;
        coords.reserve(total_vars);
        coords.insert(coords.end(), stride_bits, fr::zero());
        coords.insert(coords.end(), u.begin(), u.end());
        coords.insert(coords.end(), column_coords.begin(), column_coords.end());
        coords.insert(coords.end(), zk_bit, fr::zero());
        BB_ASSERT_EQ(coords.size(), total_vars, "opening point has the wrong arity");
        return coords;
    };

    if (!config.stack_columns) {
        // One constituent per claimed column, all at the single point u (plus the zk 0 coordinate).
        const size_t point = plan.add_point(make_point(0, {}));
        for (const bool shifted : { false, true }) {
            const auto& refs = shifted ? to_be_shifted : unshifted;
            const auto& evaluations = shifted ? shifted_evaluations : unshifted_evaluations;
            for (size_t j = 0; j < refs.size(); ++j) {
                plan.constituents.push_back({ .tree = refs[j].group,
                                              .leaf_column = refs[j].column,
                                              .shifted = shifted,
                                              .stride_bits = 0,
                                              .claims = { { point, evaluations[j] } } });
            }
        }
        if (config.zk) {
            plan.constituents.push_back({ .tree = num_groups,
                                          .leaf_column = 0,
                                          .shifted = false,
                                          .stride_bits = 0,
                                          .claims = { { point, fr::zero() } } });
        }
        return plan;
    }

    for (size_t g = 0; g < num_groups; ++g) {
        const size_t columns = group_num_columns[g];
        const size_t stack_bits = ceil_log2(columns);
        const size_t stride_bits = total_vars - (payload_vars + stack_bits + zk_bit);
        for (const bool shifted : { false, true }) {
            const auto& refs = shifted ? to_be_shifted : unshifted;
            const auto& evaluations = shifted ? shifted_evaluations : unshifted_evaluations;
            std::vector<std::optional<fr>> corner_claims(size_t(1) << stack_bits);
            size_t num_claims = 0;
            for (size_t j = 0; j < refs.size(); ++j) {
                if (refs[j].group != g) {
                    continue;
                }
                BB_ASSERT_LT(refs[j].column, columns, "claim references a padding column");
                BB_ASSERT(!corner_claims[refs[j].column].has_value(), "duplicate claim for one column");
                corner_claims[refs[j].column] = evaluations[j];
                ++num_claims;
            }
            if (num_claims == 0) {
                continue;
            }
            // Padding corners are zero columns (for T also: the shift of a padding column, whose
            // leaked-in constant comes from the next column, again padding).
            bool all_corners_known = true;
            for (size_t j = 0; j < columns; ++j) {
                all_corners_known = all_corners_known && corner_claims[j].has_value();
            }
            StackedPlan::Constituent constituent{
                .tree = g, .leaf_column = 0, .shifted = shifted, .stride_bits = stride_bits, .claims = {}
            };
            if (all_corners_known) {
                const auto tau_slice = taus.subspan(0, stack_bits);
                fr combined = fr::zero();
                for (size_t j = 0; j < columns; ++j) {
                    combined += eq_at_corner(j, tau_slice) * *corner_claims[j];
                }
                constituent.claims.push_back({ plan.add_point(make_point(stride_bits, tau_slice)), combined });
            } else {
                for (size_t j = 0; j < columns; ++j) {
                    if (!corner_claims[j].has_value()) {
                        continue;
                    }
                    std::vector<fr> column_coords(stack_bits);
                    for (size_t i = 0; i < stack_bits; ++i) {
                        column_coords[i] = ((j >> i) & 1) == 1 ? fr::one() : fr::zero();
                    }
                    constituent.claims.push_back(
                        { plan.add_point(make_point(stride_bits, column_coords)), *corner_claims[j] });
                }
            }
            plan.constituents.push_back(std::move(constituent));
        }
    }
    if (config.zk) {
        plan.constituents.push_back({ .tree = num_groups,
                                      .leaf_column = 0,
                                      .shifted = false,
                                      .stride_bits = 0,
                                      .claims = { { plan.add_point(make_point(0, taus)), fr::zero() } } });
    }
    return plan;
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
            groups.push_back(&*mask_group);
        }

        // Column-stacking challenges, then the reduction of the per-column claims to one point
        // claim per constituent of the batched oracle.
        std::vector<fr> taus(config.num_variables - config.num_payload_variables - (config.zk ? 1 : 0));
        for (size_t i = 0; i < taus.size(); ++i) {
            taus[i] = transcript->template get_challenge<fr>(detail::whir_label("tau", i));
        }
        std::vector<size_t> group_num_columns;
        group_num_columns.reserve(claims.groups.size());
        for (const WhirGroupData<Hasher>* group : claims.groups) {
            group_num_columns.push_back(group->num_columns());
        }
        detail::StackedPlan plan = detail::build_stacked_plan(config,
                                                              group_num_columns,
                                                              claims.unshifted,
                                                              claims.unshifted_evaluations,
                                                              claims.to_be_shifted,
                                                              claims.shifted_evaluations,
                                                              u,
                                                              taus);

        // Committed arrays of the constituents: the stacked array, or its coefficient shift.
        std::vector<std::span<const fr>> constituent_arrays;
        std::deque<std::vector<fr>> shifted_storage;
        for (const auto& constituent : plan.constituents) {
            const WhirGroupData<Hasher>& group = *groups[constituent.tree];
            const std::vector<fr>& base = group.arrays[group.stacked ? 0 : constituent.leaf_column];
            if (!constituent.shifted) {
                constituent_arrays.emplace_back(base);
                continue;
            }
            BB_ASSERT_EQ(base[0], fr::zero(), "to-be-shifted array must have zero constant term");
            std::vector<fr> shifted(base.begin() + 1, base.end());
            shifted.push_back(fr::zero());
            shifted_storage.push_back(std::move(shifted));
            constituent_arrays.emplace_back(shifted_storage.back());
        }

        // R_c(p) of the stride-spread constituent: eq(0, p_low)·A_c(p_high).
        auto constituent_evaluation = [&](size_t c, std::span<const fr> point) {
            const size_t stride_bits = plan.constituents[c].stride_bits;
            fr eq_low = fr::one();
            for (size_t i = 0; i < stride_bits; ++i) {
                eq_low *= fr::one() - point[i];
            }
            return eq_low * detail::mle_of_span(constituent_arrays[c], point.subspan(stride_bits));
        };

        // The zk mask's claimed evaluation at its own point.
        if (config.zk) {
            auto& mask_claim = plan.constituents.back().claims[0];
            mask_claim.value = constituent_evaluation(plan.constituents.size() - 1, plan.points[mask_claim.point]);
            transcript->send_to_verifier(std::string("WHIR:mask_eval"), mask_claim.value);
        }

        // Cross evaluations: every constituent evaluated at every other constituent's points; the
        // subsequent ρ/γ batching makes each one a claim proven by the same WHIR run.
        for (size_t t = 0; t < plan.points.size(); ++t) {
            for (size_t c = 0; c < plan.constituents.size(); ++c) {
                if (!plan.is_diagonal(t, c)) {
                    transcript->send_to_verifier(detail::whir_label("cross", t, c),
                                                 constituent_evaluation(c, plan.points[t]));
                }
            }
        }

        const fr rho = transcript->template get_challenge<fr>("WHIR:rho");
        const fr gamma_claims = transcript->template get_challenge<fr>("WHIR:gamma_claims");

        // Batched array F = ∑_c ρᶜ·spread(A_c): constituent c enters at coefficient stride 2^d.
        std::vector<fr> batched(n, fr::zero());
        fr rho_power = fr::one();
        for (size_t c = 0; c < plan.constituents.size(); ++c) {
            const std::span<const fr> array = constituent_arrays[c];
            const size_t stride_bits = plan.constituents[c].stride_bits;
            const fr scalar = rho_power;
            parallel_for_range(array.size(), [&](size_t start, size_t end) {
                for (size_t a = start; a < end; ++a) {
                    batched[a << stride_bits] += scalar * array[a];
                }
            });
            rho_power *= rho;
        }

        // Weight table of the initial claims: W = ∑_t γᵗ·eq(p_t, ·) over the committed variables.
        std::vector<fr> weight_table(n, fr::zero());
        fr gamma_power = fr::one();
        for (const auto& point : plan.points) {
            WeightTerm::eq_weight(point, gamma_power).accumulate_table(weight_table);
            gamma_power *= gamma_claims;
        }

        std::vector<fr> current = std::move(batched);
        std::deque<Tree> folded_trees; // owns the trees committed during the proof

        // Openings of the current oracle at a query index: all round-0 trees (each at the index
        // reduced to its own leaf count), or the last folded tree.
        auto send_query_openings = [&](size_t idx) {
            if (folded_trees.empty()) {
                for (const WhirGroupData<Hasher>* group : groups) {
                    send_opening(transcript, group->tree.open(idx & (group->tree.num_leaves() - 1)));
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
                        send_opening(transcript, group->tree.open(indices[s] & (group->tree.num_leaves() - 1)));
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

    /** @brief One constituent's round-0 opening contribution: leaf slot, scalar ρᶜ, stride, shift. */
    struct Contribution {
        size_t tree;
        size_t leaf_column;
        fr scalar;
        bool shifted;
        size_t stride_bits;
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

        // Group roots and the zk mask root.
        const std::vector<size_t>& group_columns = claims.group_num_columns;
        std::vector<Digest> roots = claims.group_roots;
        if (roots.empty()) {
            for (size_t g = 0; g < group_columns.size(); ++g) {
                const auto fields = transcript->template receive_from_prover<std::array<fr, Hasher::DIGEST_NUM_FIELDS>>(
                    detail::whir_label("root", g));
                roots.push_back(Hasher::digest_from_fields(fields));
            }
        }
        BB_ASSERT_EQ(roots.size(), group_columns.size(), "one root per group required");
        if (config.zk) {
            const auto fields = transcript->template receive_from_prover<std::array<fr, Hasher::DIGEST_NUM_FIELDS>>(
                std::string("WHIR:root_mask"));
            roots.push_back(Hasher::digest_from_fields(fields));
        }

        // Column-stacking challenges and the claim plan (mirrors the prover), the mask evaluation,
        // then the cross evaluations R_c(p_t) of every constituent at every foreign point.
        std::vector<fr> taus(config.num_variables - config.num_payload_variables - (config.zk ? 1 : 0));
        for (size_t i = 0; i < taus.size(); ++i) {
            taus[i] = transcript->template get_challenge<fr>(detail::whir_label("tau", i));
        }
        detail::StackedPlan plan = detail::build_stacked_plan(config,
                                                              group_columns,
                                                              claims.unshifted,
                                                              claims.unshifted_evaluations,
                                                              claims.to_be_shifted,
                                                              claims.shifted_evaluations,
                                                              u,
                                                              taus);
        if (config.zk) {
            plan.constituents.back().claims[0].value =
                transcript->template receive_from_prover<fr>(std::string("WHIR:mask_eval"));
        }
        // values[t][c] = R_c(p_t): claimed on the diagonal, read from the proof stream elsewhere.
        std::vector<std::vector<fr>> values(plan.points.size(), std::vector<fr>(plan.constituents.size(), fr::zero()));
        for (size_t c = 0; c < plan.constituents.size(); ++c) {
            for (const auto& claim : plan.constituents[c].claims) {
                values[claim.point][c] = claim.value;
            }
        }
        for (size_t t = 0; t < plan.points.size(); ++t) {
            for (size_t c = 0; c < plan.constituents.size(); ++c) {
                if (!plan.is_diagonal(t, c)) {
                    values[t][c] = transcript->template receive_from_prover<fr>(detail::whir_label("cross", t, c));
                }
            }
        }

        const fr rho = transcript->template get_challenge<fr>("WHIR:rho");
        const fr gamma_claims = transcript->template get_challenge<fr>("WHIR:gamma_claims");

        // σ₀ = ∑_t γᵗ·F(p_t) with F(p_t) = ∑_c ρᶜ·R_c(p_t); one weight term per point.
        std::vector<Contribution> contributions;
        std::vector<fr> rho_powers(plan.constituents.size());
        fr rho_power = fr::one();
        for (size_t c = 0; c < plan.constituents.size(); ++c) {
            rho_powers[c] = rho_power;
            contributions.push_back({ plan.constituents[c].tree,
                                      plan.constituents[c].leaf_column,
                                      rho_power,
                                      plan.constituents[c].shifted,
                                      plan.constituents[c].stride_bits });
            rho_power *= rho;
        }
        fr sigma = fr::zero();
        std::vector<WeightTerm> terms;
        fr gamma_power = fr::one();
        for (size_t t = 0; t < plan.points.size(); ++t) {
            fr oracle_value = fr::zero();
            for (size_t c = 0; c < plan.constituents.size(); ++c) {
                oracle_value += rho_powers[c] * values[t][c];
            }
            sigma += gamma_power * oracle_value;
            terms.push_back(WeightTerm::eq_weight(plan.points[t], gamma_power));
            gamma_power *= gamma_claims;
        }

        // Per-tree leaf layout: committed variable count (for leaf-index reduction) and how many
        // columns each leaf carries. Interleaved trees are full width with one slot per column.
        std::vector<size_t> tree_variables;
        std::vector<size_t> tree_leaf_columns;
        for (const size_t columns : group_columns) {
            tree_variables.push_back(config.stack_columns ? config.num_payload_variables + detail::ceil_log2(columns) +
                                                                (config.zk ? 1 : 0)
                                                          : config.num_variables);
            tree_leaf_columns.push_back(config.stack_columns ? 1 : columns);
        }
        if (config.zk) {
            tree_variables.push_back(config.num_variables); // the full-width mask
            tree_leaf_columns.push_back(1);
        }

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
                                               tree_variables,

                                               tree_leaf_columns,
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
                                                       tree_variables,

                                                       tree_leaf_columns,
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

    /**
     * @brief Open every round-0 tree at `idx` (reduced to each tree's own leaf count) and assemble
     * the RLC'd virtual coset values of the batched oracle.
     * @details A constituent with stride d is the univariate substitution A(X^{2^d}): its value at
     * the big-domain point x is the small codeword at x^{2^d}. Coset position t of query `idx`
     * therefore reads the tree's leaf `idx mod 2^{ib}` at slot (q + t·2^d) mod 2^k, where
     * ib = tree index bits and q = idx >> ib; shifted constituents scale by (x_t^{2^d})^{-1}.
     */
    template <typename Transcript>
    static bool read_round0_openings(const std::shared_ptr<Transcript>& transcript,
                                     const WhirConfig& config,
                                     const std::vector<Digest>& roots,
                                     const std::vector<size_t>& tree_variables,
                                     const std::vector<size_t>& tree_leaf_columns,
                                     const std::vector<Contribution>& contributions,
                                     size_t idx,
                                     size_t index_bits,
                                     const fr& omega,
                                     const fr& eta_inv,
                                     std::vector<fr>& out)
    {
        const size_t k = config.folding_factor_bits;
        const size_t arity = size_t(1) << k;
        std::vector<typename Tree::Opening> openings;
        std::vector<size_t> tree_index_bits(roots.size());
        openings.reserve(roots.size());
        for (size_t g = 0; g < roots.size(); ++g) {
            tree_index_bits[g] = tree_variables[g] + config.log_inv_rate - k;
            BB_ASSERT_LTE(tree_index_bits[g], index_bits, "tree larger than the round-0 oracle");
            openings.push_back(read_opening(transcript, tree_leaf_columns[g], arity, tree_index_bits[g], config.zk));
            if (!Tree::verify(roots[g], idx & ((size_t(1) << tree_index_bits[g]) - 1), openings.back())) {
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
            const std::vector<fr>& values = openings[contribution.tree].values;
            const size_t stride_bits = contribution.stride_bits;
            const size_t slot_base = idx >> tree_index_bits[contribution.tree];
            const size_t column_offset = contribution.leaf_column * arity;
            for (size_t t = 0; t < arity; ++t) {
                fr value =
                    contribution.scalar * values[column_offset + ((slot_base + (t << stride_bits)) & (arity - 1))];
                if (contribution.shifted) {
                    // (x_t^{2^d})^{-1} by d squarings of the coset point's inverse.
                    fr y_inv = x_inverses[t];
                    for (size_t i = 0; i < stride_bits; ++i) {
                        y_inv = y_inv.sqr();
                    }
                    value *= y_inv;
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
