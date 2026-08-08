#pragma once

#include "barretenberg/common/assert.hpp"
#include "barretenberg/common/thread.hpp"
#include "barretenberg/crypto/blake3s/blake3s.hpp"
#include "barretenberg/crypto/poseidon2/poseidon2.hpp"
#include "barretenberg/crypto/sha256/sha256.hpp"
#include "barretenberg/crypto/skyscraper/skyscraper.hpp"
#include "barretenberg/ecc/curves/bn254/fr.hpp"
#include "barretenberg/numeric/uint256/uint256.hpp"

#include <array>
#include <cstdint>
#include <cstring>
#include <optional>
#include <span>
#include <vector>

namespace bb::whir {

/**
 * @brief Merkle hasher over BN254 Fr leaves with an Fr digest; one Poseidon2 sponge call per node.
 * @details The recursion-friendly choice: an in-circuit verifier pays one Poseidon2 permutation per
 * tree level. Leaves and inner nodes are domain-separated by a tag element (0 and 1 respectively).
 */
struct Poseidon2MerkleHasher {
    using Digest = fr;
    using Poseidon2 = crypto::Poseidon2<crypto::Poseidon2Bn254ScalarFieldParams>;
    static constexpr size_t DIGEST_NUM_FIELDS = 1;

    static Digest hash_leaf(std::span<const fr> values, const std::optional<fr>& salt)
    {
        std::vector<fr> input;
        input.reserve(values.size() + 2);
        input.push_back(fr(0));
        if (salt) {
            input.push_back(*salt);
        }
        input.insert(input.end(), values.begin(), values.end());
        return Poseidon2::hash(input);
    }
    static Digest hash_node(const Digest& left, const Digest& right) { return Poseidon2::hash({ fr(1), left, right }); }
    static std::array<fr, DIGEST_NUM_FIELDS> digest_to_fields(const Digest& digest) { return { digest }; }
    static Digest digest_from_fields(std::span<const fr> fields) { return fields[0]; }
};

/**
 * @brief Poseidon2 Merkle hashing tuned for the number of permutations an in-circuit verifier pays.
 * @details `Poseidon2MerkleHasher` prepends a domain tag element to both leaves and nodes. That tag
 * occupies a rate slot, so a leaf of `n` values costs ⌈(n+1)/3⌉ permutations rather than ⌈n/3⌉ and a
 * node costs one permutation only because two children already fit in the rate. Here the domain
 * separator moves into the sponge's capacity instead: a leaf absorbs its values under the initial
 * state `(0,0,0, LEAF_IV + n)` and a node is a single permutation of `(left, right, 0, NODE_IV)`.
 * Separation is at least as strong — the capacity is never touched by absorbed data, whereas a tag
 * element shares the rate with it — and the leaf of a wide commitment costs one permutation less per
 * three values.
 *
 * A recursive verifier pays roughly 70 constraints per permutation on Ultra (Poseidon2 has its own
 * custom gates), which makes hashing the dominant term of a WHIR verification circuit; every
 * permutation removed here is removed once per query per commitment group.
 */
struct Poseidon2CompressionHasher {
    using Digest = fr;
    using Permutation = crypto::Poseidon2Permutation<crypto::Poseidon2Bn254ScalarFieldParams>;
    static constexpr size_t DIGEST_NUM_FIELDS = 1;
    static constexpr size_t RATE = 3;

    /** @brief Capacity separator for leaves; the absorbed length is added so leaves of different
     * widths cannot collide. Disjoint from `NODE_IV` and from bb's `length << 64` sponge IVs. */
    static fr leaf_iv(size_t length, bool salted)
    {
        return fr(uint256_t(1) << 128) + fr(uint256_t(length)) + (salted ? fr(uint256_t(1) << 64) : fr::zero());
    }
    static fr node_iv() { return fr(uint256_t(3) << 128); }

    static Digest hash_leaf(std::span<const fr> values, const std::optional<fr>& salt)
    {
        std::array<fr, 4> state{ fr::zero(), fr::zero(), fr::zero(), leaf_iv(values.size(), salt.has_value()) };
        size_t slot = 0;
        auto absorb = [&](const fr& value) {
            state[slot] += value;
            if (++slot == RATE) {
                state = Permutation::permutation(state);
                slot = 0;
            }
        };
        if (salt) {
            absorb(*salt);
        }
        for (const fr& value : values) {
            absorb(value);
        }
        if (slot != 0 || values.empty()) {
            state = Permutation::permutation(state);
        }
        return state[0];
    }

    static Digest hash_node(const Digest& left, const Digest& right)
    {
        return Permutation::permutation({ left, right, fr::zero(), node_iv() })[0];
    }

    static std::array<fr, DIGEST_NUM_FIELDS> digest_to_fields(const Digest& digest) { return { digest }; }
    static Digest digest_from_fields(std::span<const fr> fields) { return fields[0]; }
};

/**
 * @brief Merkle hasher over BN254 Fr with the Skyscraper-v1 compression (ePrint 2025/058).
 * @details Follows ProveKit's Merkle conventions exactly for apples-to-apples comparison: a leaf is
 * the left-fold of the two-to-one compression over its values (a single value hashes to itself), a
 * node is one compression, and there is no leaf/node domain separation. A salt, when present, is
 * folded in before the values. Like Poseidon2 this is an algebraic hash with an Fr digest, but each
 * compression is 9 double-rounds of squarings/byte S-boxes instead of a full Poseidon2 sponge.
 */
struct SkyscraperMerkleHasher {
    using Digest = fr;
    static constexpr size_t DIGEST_NUM_FIELDS = 1;

    static Digest hash_leaf(std::span<const fr> values, const std::optional<fr>& salt)
    {
        if (!salt) {
            return crypto::skyscraper::fold_compress(values);
        }
        fr acc = *salt;
        for (const fr& value : values) {
            acc = crypto::skyscraper::compress(acc, value);
        }
        return acc;
    }
    static Digest hash_node(const Digest& left, const Digest& right)
    {
        return crypto::skyscraper::compress(left, right);
    }
    static std::array<fr, DIGEST_NUM_FIELDS> digest_to_fields(const Digest& digest) { return { digest }; }
    static Digest digest_from_fields(std::span<const fr> fields) { return fields[0]; }
};

namespace detail {

/** @brief Append an Fr as its canonical 4x64-bit little-endian limbs. */
inline void append_fr_bytes(std::vector<uint8_t>& buffer, const fr& value)
{
    const uint256_t canonical(value);
    const size_t offset = buffer.size();
    buffer.resize(offset + 32);
    std::memcpy(&buffer[offset], canonical.data, 32);
}

/**
 * @brief A 32-byte digest as two 128-bit field elements.
 * @details Each half is < 2^128 < r, so the mapping digest -> fields is injective.
 */
inline std::array<fr, 2> byte_digest_to_fields(const std::array<uint8_t, 32>& digest)
{
    uint256_t lo(0);
    uint256_t hi(0);
    std::memcpy(lo.data, digest.data(), 16);
    std::memcpy(hi.data, digest.data() + 16, 16);
    return { fr(lo), fr(hi) };
}

inline std::array<uint8_t, 32> byte_digest_from_fields(std::span<const fr> fields)
{
    std::array<uint8_t, 32> digest{};
    const uint256_t lo(fields[0]);
    const uint256_t hi(fields[1]);
    std::memcpy(digest.data(), lo.data, 16);
    std::memcpy(digest.data() + 16, hi.data, 16);
    return digest;
}

} // namespace detail

/**
 * @brief Merkle hasher with 32-byte SHA-256 digests; the conservative, universally-available choice.
 * @details Structurally identical to `Blake3sMerkleHasher` — the same leaf/node domain tags and the
 * same two-128-bit-halves transcript encoding — so a benchmark against it isolates the compression
 * function. Unlike bb's blake3s there is no input-length cap, so a leaf is one hash of its whole
 * buffer rather than a chain of chunks.
 */
struct Sha256MerkleHasher {
    using Digest = std::array<uint8_t, 32>;
    static constexpr size_t DIGEST_NUM_FIELDS = 2;

    static Digest hash_leaf(std::span<const fr> values, const std::optional<fr>& salt)
    {
        std::vector<uint8_t> buffer;
        buffer.reserve(33 + 32 * values.size());
        buffer.push_back(uint8_t(0)); // leaf tag
        if (salt) {
            detail::append_fr_bytes(buffer, *salt);
        }
        for (const fr& value : values) {
            detail::append_fr_bytes(buffer, value);
        }
        return crypto::sha256(buffer);
    }
    static Digest hash_node(const Digest& left, const Digest& right)
    {
        std::vector<uint8_t> input(65);
        input[0] = uint8_t(1); // node tag
        std::memcpy(input.data() + 1, left.data(), 32);
        std::memcpy(input.data() + 33, right.data(), 32);
        return crypto::sha256(input);
    }
    static std::array<fr, DIGEST_NUM_FIELDS> digest_to_fields(const Digest& digest)
    {
        return detail::byte_digest_to_fields(digest);
    }
    static Digest digest_from_fields(std::span<const fr> fields) { return detail::byte_digest_from_fields(fields); }
};

/**
 * @brief Merkle hasher with 32-byte Blake3s digests; the fast-native-proving choice.
 * @details Fr values are absorbed as their canonical 4x64-bit little-endian limbs. bb's blake3s is
 * restricted to inputs under 1024 bytes, so wide leaves (many columns) are hashed as a chain of
 * fixed 24-element chunks, each absorbing the previous chunk's digest; the chunk structure is
 * determined by the leaf length, which is fixed per tree. Leaves and nodes are domain-separated by
 * a tag byte (0x00 and 0x01). A digest crosses the transcript as two 128-bit field elements.
 */
struct Blake3sMerkleHasher {
    using Digest = std::array<uint8_t, 32>;
    static constexpr size_t DIGEST_NUM_FIELDS = 2;
    static constexpr size_t LEAF_CHUNK_VALUES = 24;

    static Digest hash_leaf(std::span<const fr> values, const std::optional<fr>& salt)
    {
        std::vector<uint8_t> buffer;
        buffer.reserve(33 + 32 * (LEAF_CHUNK_VALUES + 1));
        buffer.push_back(uint8_t(0)); // leaf tag
        if (salt) {
            detail::append_fr_bytes(buffer, *salt);
        }
        Digest digest{};
        size_t absorbed = 0;
        while (absorbed < values.size()) {
            const size_t chunk = std::min(LEAF_CHUNK_VALUES, values.size() - absorbed);
            for (size_t t = 0; t < chunk; ++t) {
                detail::append_fr_bytes(buffer, values[absorbed + t]);
            }
            absorbed += chunk;
            digest = to_digest(blake3::blake3s(buffer));
            buffer.assign(digest.begin(), digest.end()); // chain into the next chunk
        }
        return digest;
    }
    static Digest hash_node(const Digest& left, const Digest& right)
    {
        std::vector<uint8_t> input(65);
        input[0] = uint8_t(1); // node tag
        std::memcpy(input.data() + 1, left.data(), 32);
        std::memcpy(input.data() + 33, right.data(), 32);
        return to_digest(blake3::blake3s(input));
    }
    static std::array<fr, DIGEST_NUM_FIELDS> digest_to_fields(const Digest& digest)
    {
        return detail::byte_digest_to_fields(digest);
    }
    static Digest digest_from_fields(std::span<const fr> fields) { return detail::byte_digest_from_fields(fields); }

  private:
    static Digest to_digest(const std::vector<uint8_t>& bytes)
    {
        Digest digest;
        std::memcpy(digest.data(), bytes.data(), 32);
        return digest;
    }
};

/**
 * @brief In-memory Merkle tree over a Reed-Solomon codeword with fold-coset leaf grouping.
 * @details Leaf j hashes the codeword values at the strided positions {j + t*num_leaves : t in
 * [arity]} - the fold coset of index j (README.md §4.1) - so a single authentication path opens
 * everything needed to fold at one query index. With `salted` set, each leaf additionally absorbs a
 * fresh random Fr salt (revealed only on opening), making the commitment statistically hiding
 * (README.md §8).
 */
template <typename Hasher> class MerkleTree {
  public:
    using Digest = typename Hasher::Digest;

    struct Opening {
        std::vector<fr> values; // coset values, stride order t = 0..arity-1
        std::optional<fr> salt;
        std::vector<Digest> path; // sibling digests, leaf level first
    };

    /**
     * @brief Authentication of several leaves at once, sharing everything their paths have in common.
     * @details Independent paths for t leaves of a depth-d tree repeat every node the paths meet
     * above their branch points; near the root they all coincide. This form sends a sibling only
     * where the verifier cannot already derive it from leaves it holds, which for t random leaves
     * costs about t·(d - log₂t) digests instead of t·d. `leaves` (distinct, ascending) is derived
     * from the query indices by both parties and is not part of the transmitted data.
     */
    struct BatchOpening {
        std::vector<std::vector<fr>> values; // one coset per entry of `leaves`
        std::vector<fr> salts;               // parallel to `leaves`, empty when unsalted
        std::vector<Digest> siblings;        // in the canonical order `walk_batch` visits them
    };

    MerkleTree(std::vector<fr> codeword, size_t log_arity, bool salted = false)
        : MerkleTree(single_column(std::move(codeword)), log_arity, salted)
    {}

    /**
     * @brief Tree over several same-length codewords ("columns") sharing leaves: leaf j holds every
     * column's coset-j values (column-major: values[c*arity + t]). One authentication path then
     * opens all columns of a commitment round at a query index.
     * @details The tree takes ownership of the codewords and materializes leaf values on demand in
     * `open` — leaves are strided views of the codewords, and storing them reorganized would double
     * the tree's memory (5+ GiB for an Ultra trace at 2^20).
     */
    MerkleTree(std::vector<std::vector<fr>> codewords, size_t log_arity, bool salted = false)
        : log_arity_(log_arity)
        , num_leaves_(codewords.at(0).size() >> log_arity)
        , codewords_(std::move(codewords))
    {
        const size_t arity = size_t(1) << log_arity;
        for (const auto& codeword : codewords_) {
            BB_ASSERT_EQ(codeword.size(), num_leaves_ * arity, "codewords must share one size, a multiple of arity");
        }
        BB_ASSERT_GT(num_leaves_, size_t(0));
        BB_ASSERT_EQ(num_leaves_ & (num_leaves_ - 1), size_t(0), "number of leaves must be a power of two");

        if (salted) {
            salts_.resize(num_leaves_);
            for (auto& salt : salts_) {
                salt = fr::random_element();
            }
        }

        // Leaf digests
        std::vector<Digest> level(num_leaves_);
        parallel_for_range(num_leaves_, [&](size_t start, size_t end) {
            std::vector<fr> values;
            for (size_t j = start; j < end; ++j) {
                gather_leaf_values(j, values);
                level[j] = Hasher::hash_leaf(values, salted ? std::optional<fr>(salts_[j]) : std::nullopt);
            }
        });

        levels_.push_back(std::move(level));
        while (levels_.back().size() > 1) {
            const std::vector<Digest>& below = levels_.back();
            std::vector<Digest> above(below.size() / 2);
            parallel_for_range(above.size(), [&](size_t start, size_t end) {
                for (size_t j = start; j < end; ++j) {
                    above[j] = Hasher::hash_node(below[2 * j], below[2 * j + 1]);
                }
            });
            levels_.push_back(std::move(above));
        }
    }

    const Digest& root() const { return levels_.back()[0]; }
    size_t num_leaves() const { return num_leaves_; }
    size_t depth() const { return levels_.size() - 1; }

    Opening open(size_t leaf_index) const
    {
        BB_ASSERT_LT(leaf_index, num_leaves_);
        Opening opening;
        gather_leaf_values(leaf_index, opening.values);
        if (!salts_.empty()) {
            opening.salt = salts_[leaf_index];
        }
        size_t index = leaf_index;
        for (size_t level = 0; level < depth(); ++level) {
            opening.path.push_back(levels_[level][index ^ 1]);
            index >>= 1;
        }
        return opening;
    }

    /**
     * @brief The distinct leaves of a query-index multiset, ascending: the batch's addressing.
     * @details Both parties derive this from the indices alone, so it is never transmitted.
     */
    static std::vector<size_t> batch_leaves(std::span<const size_t> leaf_indices)
    {
        std::vector<size_t> leaves(leaf_indices.begin(), leaf_indices.end());
        std::sort(leaves.begin(), leaves.end());
        leaves.erase(std::unique(leaves.begin(), leaves.end()), leaves.end());
        return leaves;
    }

    /** @brief How many sibling digests a batch over `leaves` sends; derivable before reading them. */
    static size_t batch_num_siblings(std::span<const size_t> leaves, size_t depth)
    {
        size_t count = 0;
        walk_batch(leaves, depth, [&](size_t, size_t, bool paired) { count += paired ? 0 : 1; });
        return count;
    }

    /** @brief Open every leaf of `leaf_indices` against this tree, sharing their common path nodes. */
    BatchOpening open_batch(std::span<const size_t> leaf_indices) const
    {
        const std::vector<size_t> leaves = batch_leaves(leaf_indices);
        BatchOpening opening;
        opening.values.resize(leaves.size());
        for (size_t i = 0; i < leaves.size(); ++i) {
            BB_ASSERT_LT(leaves[i], num_leaves_);
            gather_leaf_values(leaves[i], opening.values[i]);
            if (!salts_.empty()) {
                opening.salts.push_back(salts_[leaves[i]]);
            }
        }
        walk_batch(leaves, depth(), [&](size_t level, size_t index, bool paired) {
            if (!paired) {
                opening.siblings.push_back(levels_[level][index ^ 1]);
            }
        });
        return opening;
    }

    /**
     * @brief Recompute the root from a batch opening; the mirror of `open_batch`.
     * @details `depth` must be the tree's, and is climbed in full: stopping as soon as one node is
     * left would end early whenever the queried leaves all sit under one subtree.
     */
    static bool verify_batch(const Digest& root,
                             std::span<const size_t> leaves,
                             size_t depth,
                             const BatchOpening& opening)
    {
        if (opening.values.size() != leaves.size() || leaves.empty()) {
            return false;
        }
        if (!opening.salts.empty() && opening.salts.size() != leaves.size()) {
            return false;
        }
        std::vector<size_t> level(leaves.begin(), leaves.end());
        std::vector<Digest> digests(leaves.size());
        for (size_t i = 0; i < leaves.size(); ++i) {
            digests[i] = Hasher::hash_leaf(opening.values[i],
                                           opening.salts.empty() ? std::nullopt : std::optional<fr>(opening.salts[i]));
        }
        size_t cursor = 0;
        for (size_t l = 0; l < depth; ++l) {
            std::vector<size_t> parents;
            std::vector<Digest> parent_digests;
            for (size_t p = 0; p < level.size();) {
                const bool paired = p + 1 < level.size() && level[p + 1] == (level[p] ^ 1);
                parents.push_back(level[p] >> 1);
                if (paired) {
                    parent_digests.push_back(Hasher::hash_node(digests[p], digests[p + 1]));
                    p += 2;
                    continue;
                }
                if (cursor >= opening.siblings.size()) {
                    return false;
                }
                const Digest& sibling = opening.siblings[cursor++];
                parent_digests.push_back((level[p] & 1) ? Hasher::hash_node(sibling, digests[p])
                                                        : Hasher::hash_node(digests[p], sibling));
                p += 1;
            }
            level = std::move(parents);
            digests = std::move(parent_digests);
        }
        return level.size() == 1 && level[0] == 0 && cursor == opening.siblings.size() && digests[0] == root;
    }

    static bool verify(const Digest& root, size_t leaf_index, const Opening& opening)
    {
        Digest digest = Hasher::hash_leaf(opening.values, opening.salt);
        size_t index = leaf_index;
        for (const Digest& sibling : opening.path) {
            digest = (index & 1) ? Hasher::hash_node(sibling, digest) : Hasher::hash_node(digest, sibling);
            index >>= 1;
        }
        return digest == root;
    }

    /**
     * @brief The 2^`cap_levels` node digests `cap_levels` levels below the root: a "Merkle cap".
     * @details Sending the cap once per commitment and stopping every authentication path there
     * trades `2^c - 1` hashes, paid once, for `c` hashes on every query. For the t random queries a
     * WHIR round makes, that is a net saving whenever `2^c - 1 < c*t`, and the optimum sits near
     * `2^c ≈ t`. It matters because an in-circuit verifier's cost is dominated by hash count, and
     * unlike the batched-opening layout it saves that cost with no data-dependent control flow: the
     * path length is a compile-time constant and the cap is addressed by the query index's high
     * bits.
     */
    std::vector<Digest> cap(size_t cap_levels) const
    {
        BB_ASSERT_LTE(cap_levels, depth(), "cap reaches above the root");
        return levels_[depth() - cap_levels];
    }

    /** @brief Authentication path from `leaf_index` up to (not including) the cap level. */
    Opening open_capped(size_t leaf_index, size_t cap_levels) const
    {
        Opening opening = open(leaf_index);
        BB_ASSERT_LTE(cap_levels, opening.path.size(), "cap reaches above the root");
        opening.path.resize(opening.path.size() - cap_levels);
        return opening;
    }

    /** @brief `verify` against a cap: the path lands on `cap[leaf_index >> path.size()]`. */
    static bool verify_capped(std::span<const Digest> cap, size_t leaf_index, const Opening& opening)
    {
        Digest digest = Hasher::hash_leaf(opening.values, opening.salt);
        size_t index = leaf_index;
        for (const Digest& sibling : opening.path) {
            digest = (index & 1) ? Hasher::hash_node(sibling, digest) : Hasher::hash_node(digest, sibling);
            index >>= 1;
        }
        return index < cap.size() && digest == cap[index];
    }

    /** @brief Fold a cap back to the root; how a verifier binds a cap to a single-digest commitment. */
    static Digest root_from_cap(std::span<const Digest> cap)
    {
        std::vector<Digest> level(cap.begin(), cap.end());
        while (level.size() > 1) {
            std::vector<Digest> above(level.size() / 2);
            for (size_t j = 0; j < above.size(); ++j) {
                above[j] = Hasher::hash_node(level[2 * j], level[2 * j + 1]);
            }
            level = std::move(above);
        }
        return level[0];
    }

  private:
    /**
     * @brief Walk the levels a batch of leaves induces, reporting each known node once.
     * @details `visit(level, index, paired)` sees every node the verifier holds at that level, in
     * ascending index order; `paired` says its sibling is also held, so no digest need be sent for
     * it. Prover and verifier run the identical walk, which is what fixes the sibling order.
     */
    template <typename Visitor> static void walk_batch(std::span<const size_t> leaves, size_t depth, Visitor&& visit)
    {
        std::vector<size_t> level(leaves.begin(), leaves.end());
        for (size_t l = 0; l < depth; ++l) {
            std::vector<size_t> parents;
            for (size_t p = 0; p < level.size();) {
                const bool paired = p + 1 < level.size() && level[p + 1] == (level[p] ^ 1);
                visit(l, level[p], paired);
                parents.push_back(level[p] >> 1);
                p += paired ? 2 : 1;
            }
            level = std::move(parents);
        }
    }

    static std::vector<std::vector<fr>> single_column(std::vector<fr>&& codeword)
    {
        std::vector<std::vector<fr>> columns;
        columns.push_back(std::move(codeword));
        return columns;
    }

    void gather_leaf_values(size_t leaf_index, std::vector<fr>& values) const
    {
        const size_t arity = size_t(1) << log_arity_;
        values.resize(codewords_.size() * arity);
        for (size_t c = 0; c < codewords_.size(); ++c) {
            for (size_t t = 0; t < arity; ++t) {
                values[c * arity + t] = codewords_[c][leaf_index + t * num_leaves_];
            }
        }
    }

    size_t log_arity_;
    size_t num_leaves_;
    std::vector<std::vector<fr>> codewords_;
    std::vector<fr> salts_;
    std::vector<std::vector<Digest>> levels_;
};

} // namespace bb::whir
