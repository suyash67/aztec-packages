#pragma once

#include "barretenberg/common/assert.hpp"
#include "barretenberg/common/thread.hpp"
#include "barretenberg/crypto/blake3s/blake3s.hpp"
#include "barretenberg/crypto/poseidon2/poseidon2.hpp"
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

    static void append_fr_bytes(std::vector<uint8_t>& buffer, const fr& value)
    {
        const uint256_t canonical(value);
        const size_t offset = buffer.size();
        buffer.resize(offset + 32);
        std::memcpy(&buffer[offset], canonical.data, 32);
    }

    static Digest hash_leaf(std::span<const fr> values, const std::optional<fr>& salt)
    {
        std::vector<uint8_t> buffer;
        buffer.reserve(33 + 32 * (LEAF_CHUNK_VALUES + 1));
        buffer.push_back(uint8_t(0)); // leaf tag
        if (salt) {
            append_fr_bytes(buffer, *salt);
        }
        Digest digest{};
        size_t absorbed = 0;
        while (absorbed < values.size()) {
            const size_t chunk = std::min(LEAF_CHUNK_VALUES, values.size() - absorbed);
            for (size_t t = 0; t < chunk; ++t) {
                append_fr_bytes(buffer, values[absorbed + t]);
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
        // Two 128-bit halves; each is < 2^128 < r so the mapping digest -> fields is injective.
        uint256_t lo(0);
        uint256_t hi(0);
        std::memcpy(lo.data, digest.data(), 16);
        std::memcpy(hi.data, digest.data() + 16, 16);
        return { fr(lo), fr(hi) };
    }
    static Digest digest_from_fields(std::span<const fr> fields)
    {
        Digest digest;
        const uint256_t lo(fields[0]);
        const uint256_t hi(fields[1]);
        std::memcpy(digest.data(), lo.data, 16);
        std::memcpy(digest.data() + 16, hi.data, 16);
        return digest;
    }

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

    MerkleTree(std::span<const fr> codeword, size_t log_arity, bool salted = false)
        : MerkleTree(std::vector<std::span<const fr>>{ codeword }, log_arity, salted)
    {}

    /**
     * @brief Tree over several same-length codewords ("columns") sharing leaves: leaf j holds every
     * column's coset-j values (column-major: values[c*arity + t]). One authentication path then
     * opens all columns of a commitment round at a query index.
     */
    MerkleTree(const std::vector<std::span<const fr>>& codewords, size_t log_arity, bool salted = false)
        : log_arity_(log_arity)
        , num_leaves_(codewords.at(0).size() >> log_arity)
    {
        const size_t arity = size_t(1) << log_arity;
        const size_t num_columns = codewords.size();
        for (const auto& codeword : codewords) {
            BB_ASSERT_EQ(codeword.size(), num_leaves_ * arity, "codewords must share one size, a multiple of arity");
        }
        BB_ASSERT_GT(num_leaves_, size_t(0));
        BB_ASSERT_EQ(num_leaves_ & (num_leaves_ - 1), size_t(0), "number of leaves must be a power of two");

        leaf_values_.resize(num_leaves_);
        if (salted) {
            salts_.resize(num_leaves_);
            for (auto& salt : salts_) {
                salt = fr::random_element();
            }
        }

        // Leaf digests
        std::vector<Digest> level(num_leaves_);
        parallel_for_range(num_leaves_, [&](size_t start, size_t end) {
            for (size_t j = start; j < end; ++j) {
                std::vector<fr>& values = leaf_values_[j];
                values.resize(num_columns * arity);
                for (size_t c = 0; c < num_columns; ++c) {
                    for (size_t t = 0; t < arity; ++t) {
                        values[c * arity + t] = codewords[c][j + t * num_leaves_];
                    }
                }
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
        opening.values = leaf_values_[leaf_index];
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

  private:
    size_t log_arity_;
    size_t num_leaves_;
    std::vector<std::vector<fr>> leaf_values_;
    std::vector<fr> salts_;
    std::vector<std::vector<Digest>> levels_;
};

} // namespace bb::whir
