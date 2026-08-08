#pragma once

#include "barretenberg/commitment_schemes/whir/merkle_tree.hpp"
#include "barretenberg/commitment_schemes/whir/rs_code.hpp"
#include "barretenberg/commitment_schemes/whir/weights.hpp"
#include "barretenberg/commitment_schemes/whir/whir_config.hpp"
#include "barretenberg/common/assert.hpp"
#include "barretenberg/common/log.hpp"
#include "barretenberg/common/thread.hpp"
#include "barretenberg/polynomials/polynomial.hpp"
#include "barretenberg/polynomials/polynomial_arithmetic.hpp"

#include <bit>
#include <cstdlib>
#include <deque>
#include <functional>
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

/**
 * @brief Attributes the field elements the opening writes to the proof stream to named categories.
 * @details Enabled by setting WHIR_PROOF_BREAKDOWN in the environment; a no-op otherwise. Each
 * category is charged the growth of the proof between two samples of the transcript's length, so the
 * tally is the proof itself rather than a model of it.
 */
class WhirProofAccounting {
  public:
    static bool enabled()
    {
        static const bool on = std::getenv("WHIR_PROOF_BREAKDOWN") != nullptr;
        return on;
    }

    template <typename Transcript>
    explicit WhirProofAccounting(const std::shared_ptr<Transcript>& transcript)
        : size_at_mark_(transcript->get_proof_size())
    {}

    /** @brief Charge everything written since the last mark to `category`. */
    template <typename Transcript> void mark(const std::shared_ptr<Transcript>& transcript, const std::string& category)
    {
        const size_t size = transcript->get_proof_size();
        add(category, size - size_at_mark_);
        size_at_mark_ = size;
    }

    /** @brief Charge `count` field elements to `category` directly, for a hand-split write. */
    void add(const std::string& category, size_t count)
    {
        if (!enabled()) {
            return;
        }
        for (auto& [name, tally] : categories_) {
            if (name == category) {
                tally += count;
                return;
            }
        }
        categories_.emplace_back(category, count);
    }

    /** @brief Discard the pending span, for writes already charged through `add`. */
    template <typename Transcript> void sync(const std::shared_ptr<Transcript>& transcript)
    {
        size_at_mark_ = transcript->get_proof_size();
    }

    void report() const
    {
        if (!enabled()) {
            return;
        }
        size_t total = 0;
        for (const auto& [name, count] : categories_) {
            total += count;
        }
        info("WHIR proof breakdown (field elements, 32 bytes each):");
        for (const auto& [name, count] : categories_) {
            info("  ",
                 name,
                 std::string(std::max<size_t>(1, 22 - std::min<size_t>(21, name.size())), ' '),
                 count,
                 "\t",
                 count * 32,
                 " B\t",
                 total == 0 ? 0 : (100 * count) / total,
                 "%");
        }
        info("  WHIR opening total     ", total, "\t", total * 32, " B");
    }

  private:
    size_t size_at_mark_;
    std::vector<std::pair<std::string, size_t>> categories_;
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
        MerkleTree<Hasher> tree(std::move(codewords), config.initial_folding_factor_bits, /*salted=*/config.zk);
        return { std::move(arrays), num_columns, stacked, std::move(tree) };
    }

    /** @brief Fill the zk blinding coefficients just above a payload of `payload_bits` variables. */
    void add_blinding(std::vector<fr>& array, size_t payload_bits, bool to_be_shifted) const
    {
        if (!config.zk) {
            return;
        }
        const size_t offset = (size_t(1) << payload_bits) + (to_be_shifted ? 1 : 0);
        BB_ASSERT_LTE(offset + config.num_blinding_coefficients, array.size(), "blinding overruns the array");
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
 * @brief Bits of a challenge that query indices are cut from.
 * @details The recursive verifier pins the digit decomposition to the canonical one by requiring the
 * high part to be strictly below ⌊r/2^{QUERY_DIGIT_BITS}⌋, which rejects an honest challenge with
 * probability about 2^{QUERY_DIGIT_BITS}/r. At 240 bits that is 2^-15 — one honest proof in tens of
 * thousands would be unprovable, since Fiat-Shamir gives the prover no second draw. 128 bits puts it
 * at 2^-128 and still fits nine or more indices per challenge.
 */
static constexpr size_t QUERY_DIGIT_BITS = 128;

/**
 * @brief How many query indices are cut from a single challenge.
 * @details A challenge carries ~254 bits and a query index needs `index_bits` of them, so drawing
 * one challenge per query wastes both a sponge permutation and, in a recursive verifier, a whole
 * decomposition per query. The digits are taken from the low end of the canonical representation.
 */
inline size_t indices_per_challenge(size_t index_bits)
{
    return std::max<size_t>(1, QUERY_DIGIT_BITS / index_bits);
}

/** @brief `count` query indices, cut `indices_per_challenge` at a time from successive challenges. */
template <typename Transcript>
std::vector<size_t> draw_query_indices(const std::shared_ptr<Transcript>& transcript,
                                       const std::function<std::string(size_t)>& label,
                                       size_t count,
                                       size_t index_bits)
{
    const size_t per_challenge = indices_per_challenge(index_bits);
    const uint64_t mask = (uint64_t(1) << index_bits) - 1;
    std::vector<size_t> indices;
    indices.reserve(count);
    for (size_t chunk = 0; indices.size() < count; ++chunk) {
        const uint256_t challenge(transcript->template get_challenge<fr>(label(chunk)));
        for (size_t j = 0; j < per_challenge && indices.size() < count; ++j) {
            indices.push_back(static_cast<size_t>(((challenge >> (j * index_bits)).data[0]) & mask));
        }
    }
    return indices;
}

/**
 * @brief Does `nonce` satisfy the round's proof of work against `seed`?
 * @details Blake3 of the seed's canonical limbs followed by the little-endian nonce, accepted when
 * the leading `pow_bits` bits of the digest are zero. Blake3 rather than the transcript's Poseidon2
 * because the prover runs this 2^{pow_bits} times: the search has to be cheap for the honest party
 * and is the same work for a dishonest one.
 */
/**
 * @brief Poseidon2 proof of work: accepted when the low `pow_bits` bits of `Poseidon2(seed, nonce)`
 * are zero.
 * @details The recursion-friendly alternative to the Blake3 form. A recursive verifier checks a
 * Blake3 grind at tens of thousands of constraints per round, which costs more than the queries the
 * grinding removed; this costs one permutation plus a range constraint on the quotient. Divisibility
 * rather than leading zeros because "the witness `q` with `digest = q * 2^b`" is the cheap statement
 * in a field circuit. The digest is uniform over [0, r), so the acceptance probability is 2^-b to
 * within the rounding of r/2^b.
 */
inline bool poseidon2_pow_is_valid(const fr& seed, uint64_t nonce, size_t pow_bits)
{
    if (pow_bits == 0) {
        return true;
    }
    BB_ASSERT_LT(pow_bits, size_t(64), "proof of work is capped at 63 bits");
    using Poseidon2 = crypto::Poseidon2<crypto::Poseidon2Bn254ScalarFieldParams>;
    const uint256_t digest(Poseidon2::hash({ seed, fr(nonce) }));
    return (digest.data[0] & ((uint64_t(1) << pow_bits) - 1)) == 0;
}

inline bool pow_is_valid(const fr& seed, uint64_t nonce, size_t pow_bits)
{
    if (pow_bits == 0) {
        return true;
    }
    BB_ASSERT_LT(pow_bits, size_t(64), "proof of work is capped at 63 bits");
    std::array<uint8_t, 40> input{};
    const uint256_t limbs = uint256_t(seed);
    for (size_t l = 0; l < 4; ++l) {
        for (size_t b = 0; b < 8; ++b) {
            input[l * 8 + b] = static_cast<uint8_t>(limbs.data[l] >> (8 * b));
        }
    }
    for (size_t b = 0; b < 8; ++b) {
        input[32 + b] = static_cast<uint8_t>(nonce >> (8 * b));
    }
    const std::array<uint8_t, 32> digest = blake3::blake3s_constexpr(input.data(), input.size());
    uint64_t leading = 0;
    for (size_t b = 0; b < 8; ++b) {
        leading = (leading << 8) | digest[b];
    }
    return (leading >> (64 - pow_bits)) == 0;
}

/**
 * @brief The least nonce satisfying `pow_is_valid`; the prover's side of the grind.
 * @details Threads take disjoint contiguous chunks of one block at a time and the block's smallest
 * winner is returned, so the nonce is the same whatever the thread count — proving twice gives the
 * same proof.
 */
inline bool pow_is_valid(const fr& seed, uint64_t nonce, size_t pow_bits, bool poseidon2)
{
    return poseidon2 ? poseidon2_pow_is_valid(seed, nonce, pow_bits) : pow_is_valid(seed, nonce, pow_bits);
}

inline uint64_t grind(const fr& seed, size_t pow_bits, bool poseidon2 = false)
{
    if (pow_bits == 0) {
        return 0;
    }
    constexpr uint64_t chunk = 1024;
    constexpr uint64_t no_winner = std::numeric_limits<uint64_t>::max();
    const size_t num_threads = std::max<size_t>(1, get_num_cpus());
    std::vector<uint64_t> found(num_threads);
    for (uint64_t base = 0;; base += chunk * num_threads) {
        std::fill(found.begin(), found.end(), no_winner);
        parallel_for(num_threads, [&](size_t t) {
            const uint64_t start = base + chunk * t;
            for (uint64_t nonce = start; nonce < start + chunk; ++nonce) {
                if (pow_is_valid(seed, nonce, pow_bits, poseidon2)) {
                    found[t] = nonce;
                    return;
                }
            }
        });
        const uint64_t nonce = *std::min_element(found.begin(), found.end());
        if (nonce != no_winner) {
            return nonce;
        }
    }
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

/**
 * @brief The univariate evaluation at `z` of a constituent spread to stride 2^d.
 * @details spread(A)(x) = A(x^{2^d}), so the value at z is the array's own univariate evaluation at
 * z^{2^d}. With d = 0 (the default, non-stacked layout) this is just A(z).
 */
inline fr spread_univariate_evaluation(std::span<const fr> array, const fr& z, size_t stride_bits)
{
    fr point = z;
    for (size_t i = 0; i < stride_bits; ++i) {
        point = point.sqr();
    }
    return polynomial_arithmetic::evaluate(array.data(), point, array.size());
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
        const size_t n = size_t(1) << config.num_variables;
        BB_ASSERT_EQ(u.size(), config.num_payload_variables, "opening point size mismatch");
        BB_ASSERT_EQ(claims.unshifted.size(), claims.unshifted_evaluations.size());
        BB_ASSERT_EQ(claims.to_be_shifted.size(), claims.shifted_evaluations.size());

        WhirProofAccounting acct(transcript);

        if (claims.send_roots) {
            for (size_t g = 0; g < claims.groups.size(); ++g) {
                transcript->send_to_verifier(detail::whir_label("root", g),
                                             Hasher::digest_to_fields(claims.groups[g]->tree.root()));
            }
        }
        acct.mark(transcript, "group roots");

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
        acct.mark(transcript, "cross evaluations");

        // Out-of-domain samples against the round-0 commitments (README.md §4.2). In any
        // list-decoding regime a committed word can be close to several codewords; answering at a
        // random point outside the evaluation domain singles out one of them, which is what makes
        // the claims below well defined. Drawn before ρ, so each committed array is pinned
        // individually rather than only their batch.
        std::vector<fr> z_ood_initial(config.num_ood_samples);
        for (size_t s = 0; s < config.num_ood_samples; ++s) {
            z_ood_initial[s] = transcript->template get_challenge<fr>(detail::whir_label("z_ood_init", s));
            for (size_t c = 0; c < plan.constituents.size(); ++c) {
                transcript->send_to_verifier(detail::whir_label("y_ood_init", s, c),
                                             detail::spread_univariate_evaluation(constituent_arrays[c],
                                                                                  z_ood_initial[s],
                                                                                  plan.constituents[c].stride_bits));
            }
        }
        acct.mark(transcript, "initial OOD answers");

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

        // Weight table of the initial claims: W = ∑_t γᵗ·eq(p_t, ·) over the committed variables,
        // followed by one univariate pow-weight per out-of-domain sample.
        std::vector<fr> weight_table(n, fr::zero());
        fr gamma_power = fr::one();
        for (const auto& point : plan.points) {
            WeightTerm::eq_weight(point, gamma_power).accumulate_table(weight_table);
            gamma_power *= gamma_claims;
        }
        for (const fr& z : z_ood_initial) {
            WeightTerm::pow_weight(z, config.num_variables, gamma_power).accumulate_table(weight_table);
            gamma_power *= gamma_claims;
        }

        std::vector<fr> current = std::move(batched);
        std::deque<Tree> folded_trees; // owns the trees committed during the proof

        // Openings of the current oracle at a round's whole query set: all round-0 trees (each at the
        // indices reduced to its own leaf count), or the last folded tree. One batch per tree, so
        // the queries share every path node they have in common.
        auto send_query_openings = [&](std::span<const size_t> indices, const std::string& category) {
            if (folded_trees.empty()) {
                for (const WhirGroupData<Hasher>* group : groups) {
                    send_batch_opening(
                        transcript, group->tree, reduce_indices(indices, group->tree), acct, category, config);
                }
            } else {
                send_batch_opening(transcript, folded_trees.back(), indices, acct, category, config);
            }
        };

        for (size_t i = 0; i < config.rounds.size(); ++i) {
            const WhirRound& round = config.rounds[i];
            const size_t k = round.folding_factor_bits;

            // 1. k sumcheck rounds; array and weight table fold at each challenge.
            for (size_t j = 0; j < k; ++j) {
                const auto h = detail::sumcheck_round_univariate(current, weight_table);
                transcript->send_to_verifier(detail::whir_label("sc", i, j), h);
                const fr alpha = transcript->template get_challenge<fr>(detail::whir_label("alpha", i, j));
                fold_array_in_place(current, alpha);
                fold_array_in_place(weight_table, alpha);
            }
            acct.mark(transcript, "sumcheck univariates");

            // 2. Commit the folded polynomial on the halved domain.
            const auto& next_domain = ck.domains.get(size_t(1) << (round.log_domain_size - 1));
            std::vector<fr> codeword = rs_encode(current, next_domain, ck.domains.round_roots());
            folded_trees.emplace_back(std::move(codeword), config.committed_arity_bits(i), /*salted=*/false);
            transcript->send_to_verifier(detail::whir_label("root_g", i + 1),
                                         Hasher::digest_to_fields(folded_trees.back().root()));
            acct.mark(transcript, "folded-oracle roots");

            // 3. Out-of-domain sample.
            const fr z_ood = transcript->template get_challenge<fr>(detail::whir_label("z_ood", i));
            const fr y_ood = polynomial_arithmetic::evaluate(current.data(), z_ood, current.size());
            transcript->send_to_verifier(detail::whir_label("y_ood", i), y_ood);
            acct.mark(transcript, "per-round OOD answers");

            // 4. In-domain queries against the round-i oracle. Note the openings are needed for step
            // 2 of the NEXT iteration's committed tree only when i = 0... they always target the
            // previous oracle, so they are emitted before the trees rotate below.
            const size_t index_bits = round.log_domain_size - k;
            // Grind before the indices are drawn: forging this round means redoing the work.
            const fr pow_seed = transcript->template get_challenge<fr>(detail::whir_label("pow", i));
            transcript->send_to_verifier(detail::whir_label("nonce", i),
                                         fr(detail::grind(pow_seed, config.pow_bits, config.poseidon2_pow)));
            const std::vector<size_t> indices = detail::draw_query_indices(
                transcript,
                [&](size_t chunk) { return detail::whir_label("query", i, chunk); },
                round.num_queries,
                index_bits);
            acct.mark(transcript, "grinding nonces");
            if (i > 0) {
                send_batch_opening(transcript, folded_trees[i - 1], indices, acct, "folded-round queries", config);
            } else {
                for (const WhirGroupData<Hasher>* group : groups) {
                    send_batch_opening(
                        transcript, group->tree, reduce_indices(indices, group->tree), acct, "round-0 queries", config);
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
        acct.mark(transcript, "final polynomial");
        const size_t index_bits = config.final_round.log_domain_size - config.final_round.folding_factor_bits;
        const fr final_pow_seed = transcript->template get_challenge<fr>(std::string("WHIR:fpow"));
        transcript->send_to_verifier(std::string("WHIR:fnonce"),
                                     fr(detail::grind(final_pow_seed, config.pow_bits, config.poseidon2_pow)));
        const std::vector<size_t> final_indices = detail::draw_query_indices(
            transcript,
            [](size_t chunk) { return detail::whir_label("fquery", chunk); },
            config.final_round.num_queries,
            index_bits);
        acct.mark(transcript, "grinding nonces");
        send_query_openings(final_indices, "final-round queries");
        if (WhirProofAccounting::enabled()) {
            size_t total_columns = 0;
            for (const WhirGroupData<Hasher>* group : groups) {
                total_columns += group->leaf_columns();
            }
            info("WHIR schedule: ",
                 config.num_variables,
                 " vars, rate 2^-",
                 config.log_inv_rate,
                 ", k0=",
                 config.initial_folding_factor_bits,
                 ", k=",
                 config.folding_factor_bits,
                 ", ",
                 groups.size(),
                 " round-0 trees holding ",
                 total_columns,
                 " leaf columns, ",
                 plan.constituents.size(),
                 " constituents");
            for (size_t i = 0; i < config.rounds.size(); ++i) {
                info("  round ",
                     i,
                     ": ",
                     config.rounds[i].num_queries,
                     " queries, domain 2^",
                     config.rounds[i].log_domain_size,
                     ", rate 2^-",
                     config.rounds[i].log_inv_rate,
                     ", fold 2^",
                     config.rounds[i].folding_factor_bits);
            }
            info("  final: ",
                 config.final_round.num_queries,
                 " queries, domain 2^",
                 config.final_round.log_domain_size,
                 ", rate 2^-",
                 config.final_round.log_inv_rate);
        }
        acct.report();
    }

  private:
    /** @brief The query indices reduced to a tree's own leaf count, for a tree narrower than the
     * round-0 oracle (which only arises when stacking). */
    static std::vector<size_t> reduce_indices(std::span<const size_t> indices, const Tree& tree)
    {
        std::vector<size_t> reduced(indices.begin(), indices.end());
        for (size_t& index : reduced) {
            index &= tree.num_leaves() - 1;
        }
        return reduced;
    }

    /**
     * @brief Authenticate a round's whole query set against one tree in a single batch.
     * @details Values come first (the distinct leaves in ascending index order), then the sibling
     * digests the verifier cannot derive. Both sides recover the leaf set from the query indices,
     * so nothing about the addressing is transmitted.
     */
    template <typename Transcript>
    static void send_batch_opening(const std::shared_ptr<Transcript>& transcript,
                                   const Tree& tree,
                                   std::span<const size_t> indices,
                                   WhirProofAccounting& acct,
                                   const std::string& category,
                                   const WhirConfig& config)
    {
        if (config.per_query_openings) {
            send_capped_openings(transcript, tree, indices, acct, category, config);
            return;
        }
        const typename Tree::BatchOpening opening = tree.open_batch(indices);
        std::vector<fr> flat;
        for (size_t i = 0; i < opening.values.size(); ++i) {
            flat.insert(flat.end(), opening.values[i].begin(), opening.values[i].end());
            if (!opening.salts.empty()) {
                flat.push_back(opening.salts[i]);
            }
        }
        acct.add(category + ": leaf values", flat.size());
        for (const auto& digest : opening.siblings) {
            const auto fields = Hasher::digest_to_fields(digest);
            flat.insert(flat.end(), fields.begin(), fields.end());
        }
        acct.add(category + ": merkle paths", opening.siblings.size() * Hasher::DIGEST_NUM_FIELDS);
        transcript->send_unhashed_to_verifier(flat);
        acct.sync(transcript);
    }

    /**
     * @brief One independent authentication path per query, stopping at the tree's Merkle cap.
     * @details The layout an in-circuit verifier wants (`WhirConfig::per_query_openings`): the cap
     * first, then query `s`'s leaf values and its `depth - cap_levels` sibling digests, in query
     * order and repeated verbatim when two queries collide. Nothing about the ordering or the
     * distinct-leaf structure has to be reconstructed, so every array index is a constant.
     */
    template <typename Transcript>
    static void send_capped_openings(const std::shared_ptr<Transcript>& transcript,
                                     const Tree& tree,
                                     std::span<const size_t> indices,
                                     WhirProofAccounting& acct,
                                     const std::string& category,
                                     const WhirConfig& config)
    {
        const size_t cap_levels = config.cap_levels_for(indices.size(), tree.depth());
        std::vector<fr> flat;
        for (const auto& digest : tree.cap(cap_levels)) {
            const auto fields = Hasher::digest_to_fields(digest);
            flat.insert(flat.end(), fields.begin(), fields.end());
        }
        acct.add(category + ": merkle cap", flat.size());

        size_t value_fields = 0;
        size_t path_fields = 0;
        for (const size_t index : indices) {
            const typename Tree::Opening opening = tree.open_capped(index, cap_levels);
            flat.insert(flat.end(), opening.values.begin(), opening.values.end());
            value_fields += opening.values.size();
            if (opening.salt) {
                flat.push_back(*opening.salt);
                value_fields += 1;
            }
            for (const auto& digest : opening.path) {
                const auto fields = Hasher::digest_to_fields(digest);
                flat.insert(flat.end(), fields.begin(), fields.end());
            }
            path_fields += opening.path.size() * Hasher::DIGEST_NUM_FIELDS;
        }
        acct.add(category + ": leaf values", value_fields);
        acct.add(category + ": merkle paths", path_fields);
        transcript->send_unhashed_to_verifier(flat);
        acct.sync(transcript);
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

        // Out-of-domain samples against the round-0 commitments; ood_values[s][c] is constituent c's
        // univariate evaluation at the s-th sampled point.
        std::vector<fr> z_ood_initial(config.num_ood_samples);
        std::vector<std::vector<fr>> ood_values(config.num_ood_samples,
                                                std::vector<fr>(plan.constituents.size(), fr::zero()));
        for (size_t s = 0; s < config.num_ood_samples; ++s) {
            z_ood_initial[s] = transcript->template get_challenge<fr>(detail::whir_label("z_ood_init", s));
            for (size_t c = 0; c < plan.constituents.size(); ++c) {
                ood_values[s][c] = transcript->template receive_from_prover<fr>(detail::whir_label("y_ood_init", s, c));
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
        for (size_t s = 0; s < config.num_ood_samples; ++s) {
            fr oracle_value = fr::zero();
            for (size_t c = 0; c < plan.constituents.size(); ++c) {
                oracle_value += rho_powers[c] * ood_values[s][c];
            }
            sigma += gamma_power * oracle_value;
            terms.push_back(WeightTerm::pow_weight(z_ood_initial[s], config.num_variables, gamma_power));
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
            const size_t k = round.folding_factor_bits;

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
            const fr pow_seed = transcript->template get_challenge<fr>(detail::whir_label("pow", i));
            const fr nonce = transcript->template receive_from_prover<fr>(detail::whir_label("nonce", i));
            if (!detail::pow_is_valid(pow_seed, uint256_t(nonce).data[0], config.pow_bits, config.poseidon2_pow)) {
                return false;
            }
            const std::vector<size_t> indices = detail::draw_query_indices(
                transcript,
                [&](size_t chunk) { return detail::whir_label("query", i, chunk); },
                round.num_queries,
                index_bits);
            const fr omega = fr::get_root_of_unity(round.log_domain_size);
            const fr eta_inv = omega.pow(uint256_t(uint64_t(1)) << index_bits).invert();
            RoundOpenings openings;
            const bool read_ok =
                folded_root
                    ? read_folded_openings(transcript, config, k, *folded_root, indices, index_bits, openings)
                    : read_round0_openings(
                          transcript, config, roots, tree_variables, tree_leaf_columns, indices, index_bits, openings);
            if (!read_ok) {
                return false;
            }
            std::vector<fr> folded_values(round.num_queries);
            for (size_t s = 0; s < round.num_queries; ++s) {
                std::vector<fr> virtual_values;
                if (folded_root) {
                    const std::span<const fr> coset =
                        openings.coset_of_query(0, s, indices[s], config.per_query_openings);
                    virtual_values.assign(coset.begin(), coset.end());
                } else {
                    round0_coset(config, openings, contributions, s, indices[s], omega, eta_inv, virtual_values);
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
        const size_t final_k = config.final_round.folding_factor_bits;
        const size_t index_bits = config.final_round.log_domain_size - final_k;
        const fr omega = fr::get_root_of_unity(config.final_round.log_domain_size);
        const fr eta = omega.pow(uint256_t(uint64_t(1)) << index_bits);
        const fr eta_inv = eta.invert();
        const fr final_pow_seed = transcript->template get_challenge<fr>(std::string("WHIR:fpow"));
        const fr final_nonce = transcript->template receive_from_prover<fr>(std::string("WHIR:fnonce"));
        if (!detail::pow_is_valid(
                final_pow_seed, uint256_t(final_nonce).data[0], config.pow_bits, config.poseidon2_pow)) {
            return false;
        }
        const std::vector<size_t> final_indices = detail::draw_query_indices(
            transcript,
            [](size_t chunk) { return detail::whir_label("fquery", chunk); },
            config.final_round.num_queries,
            index_bits);
        RoundOpenings final_openings;
        const bool final_read_ok =
            folded_root ? read_folded_openings(
                              transcript, config, final_k, *folded_root, final_indices, index_bits, final_openings)
                        : read_round0_openings(transcript,
                                               config,
                                               roots,
                                               tree_variables,
                                               tree_leaf_columns,
                                               final_indices,
                                               index_bits,
                                               final_openings);
        if (!final_read_ok) {
            return false;
        }
        for (size_t s = 0; s < final_indices.size(); ++s) {
            const size_t idx = final_indices[s];
            std::vector<fr> virtual_values;
            if (folded_root) {
                const std::span<const fr> coset = final_openings.coset_of_query(0, s, idx, config.per_query_openings);
                virtual_values.assign(coset.begin(), coset.end());
            } else {
                round0_coset(config, final_openings, contributions, s, idx, omega, eta_inv, virtual_values);
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
     * @brief A round's authenticated openings: one batch per queried tree.
     * @details Round 0 queries every committed group's tree; later rounds query the single folded
     * oracle. `leaves[t]` is the distinct ascending leaf set batch t addresses, derived from the
     * query indices, so `coset` can find any query's values by binary search.
     */
    struct RoundOpenings {
        std::vector<typename Tree::BatchOpening> batches;
        std::vector<std::vector<size_t>> leaves;
        std::vector<size_t> tree_index_bits; // round 0 only: each tree's own index width

        std::span<const fr> coset(size_t tree, size_t leaf) const
        {
            const auto& list = leaves[tree];
            const auto it = std::lower_bound(list.begin(), list.end(), leaf);
            if (it == list.end() || *it != leaf) {
                throw_or_abort("WHIR: queried leaf missing from the batch");
            }
            return batches[tree].values[static_cast<size_t>(it - list.begin())];
        }

        /** @brief Query `s`'s coset. In per-query mode the openings are already in query order; the
         * batched layout addresses them through the distinct ascending leaf set instead. */
        std::span<const fr> coset_of_query(size_t tree, size_t query, size_t leaf, bool per_query) const
        {
            return per_query ? std::span<const fr>(batches[tree].values[query]) : coset(tree, leaf);
        }
    };

    /**
     * @brief Read one tree's per-query openings: the Merkle cap, then each query's leaf and path.
     * @details The mirror of `WhirProver::send_capped_openings`. Every query is authenticated on its
     * own against the cap, and the cap is folded back to `root`, so the check is exactly as strong
     * as the batched one while touching the proof stream at fixed offsets.
     */
    template <typename Transcript>
    static bool read_capped(const std::shared_ptr<Transcript>& transcript,
                            const Digest& root,
                            std::span<const size_t> indices,
                            size_t num_columns,
                            size_t arity,
                            size_t depth,
                            bool salted,
                            const WhirConfig& config,
                            typename Tree::BatchOpening& opening)
    {
        const size_t cap_levels = config.cap_levels_for(indices.size(), depth);
        const size_t cap_size = size_t(1) << cap_levels;
        const size_t path_length = depth - cap_levels;
        const size_t coset_fields = num_columns * arity + (salted ? 1 : 0);
        const std::vector<fr> flat = transcript->receive_unhashed_from_prover(
            cap_size * Hasher::DIGEST_NUM_FIELDS +
            indices.size() * (coset_fields + path_length * Hasher::DIGEST_NUM_FIELDS));

        size_t cursor = 0;
        std::vector<Digest> cap(cap_size);
        for (size_t j = 0; j < cap_size; ++j) {
            cap[j] = Hasher::digest_from_fields(std::span<const fr>(flat).subspan(cursor, Hasher::DIGEST_NUM_FIELDS));
            cursor += Hasher::DIGEST_NUM_FIELDS;
        }
        if (Tree::root_from_cap(cap) != root) {
            return false;
        }

        opening.values.resize(indices.size());
        for (size_t s = 0; s < indices.size(); ++s) {
            typename Tree::Opening path_opening;
            path_opening.values.assign(flat.begin() + static_cast<std::ptrdiff_t>(cursor),
                                       flat.begin() + static_cast<std::ptrdiff_t>(cursor + num_columns * arity));
            cursor += num_columns * arity;
            if (salted) {
                path_opening.salt = flat[cursor++];
                opening.salts.push_back(*path_opening.salt);
            }
            for (size_t l = 0; l < path_length; ++l) {
                path_opening.path.push_back(
                    Hasher::digest_from_fields(std::span<const fr>(flat).subspan(cursor, Hasher::DIGEST_NUM_FIELDS)));
                cursor += Hasher::DIGEST_NUM_FIELDS;
            }
            if (!Tree::verify_capped(cap, indices[s], path_opening)) {
                return false;
            }
            opening.values[s] = std::move(path_opening.values);
        }
        return true;
    }

    /** @brief Read one tree's batch off the proof stream and authenticate it against `root`. */
    template <typename Transcript>
    static bool read_batch(const std::shared_ptr<Transcript>& transcript,
                           const Digest& root,
                           const std::vector<size_t>& leaves,
                           size_t num_columns,
                           size_t arity,
                           size_t depth,
                           bool salted,
                           typename Tree::BatchOpening& opening)
    {
        const size_t coset_fields = num_columns * arity + (salted ? 1 : 0);
        const size_t num_siblings = Tree::batch_num_siblings(leaves, depth);
        const std::vector<fr> flat = transcript->receive_unhashed_from_prover(leaves.size() * coset_fields +
                                                                              num_siblings * Hasher::DIGEST_NUM_FIELDS);
        size_t cursor = 0;
        opening.values.resize(leaves.size());
        for (size_t i = 0; i < leaves.size(); ++i) {
            opening.values[i].assign(flat.begin() + static_cast<std::ptrdiff_t>(cursor),
                                     flat.begin() + static_cast<std::ptrdiff_t>(cursor + num_columns * arity));
            cursor += num_columns * arity;
            if (salted) {
                opening.salts.push_back(flat[cursor++]);
            }
        }
        opening.siblings.reserve(num_siblings);
        for (size_t l = 0; l < num_siblings; ++l) {
            opening.siblings.push_back(
                Hasher::digest_from_fields(std::span<const fr>(flat).subspan(cursor, Hasher::DIGEST_NUM_FIELDS)));
            cursor += Hasher::DIGEST_NUM_FIELDS;
        }
        return Tree::verify_batch(root, leaves, depth, opening);
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
                                     std::span<const size_t> indices,
                                     size_t index_bits,
                                     RoundOpenings& openings)
    {
        const size_t k = config.initial_folding_factor_bits;
        const size_t arity = size_t(1) << k;
        openings.batches.resize(roots.size());
        openings.leaves.resize(roots.size());
        openings.tree_index_bits.resize(roots.size());
        for (size_t g = 0; g < roots.size(); ++g) {
            openings.tree_index_bits[g] = tree_variables[g] + config.log_inv_rate - k;
            BB_ASSERT_LTE(openings.tree_index_bits[g], index_bits, "tree larger than the round-0 oracle");
            const size_t mask = (size_t(1) << openings.tree_index_bits[g]) - 1;
            std::vector<size_t> reduced(indices.begin(), indices.end());
            for (size_t& index : reduced) {
                index &= mask;
            }
            if (config.per_query_openings) {
                openings.leaves[g] = reduced;
                if (!read_capped(transcript,
                                 roots[g],
                                 reduced,
                                 tree_leaf_columns[g],
                                 arity,
                                 openings.tree_index_bits[g],
                                 config.zk,
                                 config,
                                 openings.batches[g])) {
                    return false;
                }
                continue;
            }
            openings.leaves[g] = Tree::batch_leaves(reduced);
            if (!read_batch(transcript,
                            roots[g],
                            openings.leaves[g],
                            tree_leaf_columns[g],
                            arity,
                            openings.tree_index_bits[g],
                            config.zk,
                            openings.batches[g])) {
                return false;
            }
        }
        return true;
    }

    /** @brief The batched oracle's virtual coset values at query `idx`, from the round-0 openings. */
    static void round0_coset(const WhirConfig& config,
                             const RoundOpenings& openings,
                             const std::vector<Contribution>& contributions,
                             size_t query,
                             size_t idx,
                             const fr& omega,
                             const fr& eta_inv,
                             std::vector<fr>& out)
    {
        const size_t arity = size_t(1) << config.initial_folding_factor_bits;
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
            const size_t tree_bits = openings.tree_index_bits[contribution.tree];
            const std::span<const fr> values = openings.coset_of_query(
                contribution.tree, query, idx & ((size_t(1) << tree_bits) - 1), config.per_query_openings);
            const size_t stride_bits = contribution.stride_bits;
            const size_t slot_base = idx >> tree_bits;
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
    }

    /** @brief Open the single-column folded oracle g_i at a round's whole query set. */
    template <typename Transcript>
    static bool read_folded_openings(const std::shared_ptr<Transcript>& transcript,
                                     const WhirConfig& config,
                                     size_t folding_factor_bits,
                                     const Digest& root,
                                     std::span<const size_t> indices,
                                     size_t index_bits,
                                     RoundOpenings& openings)
    {
        openings.batches.resize(1);
        if (config.per_query_openings) {
            openings.leaves = { std::vector<size_t>(indices.begin(), indices.end()) };
            return read_capped(transcript,
                               root,
                               indices,
                               /*num_columns=*/1,
                               size_t(1) << folding_factor_bits,
                               index_bits,
                               /*salted=*/false,
                               config,
                               openings.batches[0]);
        }
        openings.leaves = { Tree::batch_leaves(indices) };
        return read_batch(transcript,
                          root,
                          openings.leaves[0],
                          /*num_columns=*/1,
                          size_t(1) << folding_factor_bits,
                          index_bits,
                          /*salted=*/false,
                          openings.batches[0]);
    }
};

} // namespace bb::whir
