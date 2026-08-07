#include "barretenberg/commitment_schemes/whir/merkle_tree.hpp"

#include <gtest/gtest.h>

namespace bb::whir {

template <typename Hasher> class WhirMerkleTreeTest : public ::testing::Test {};
using HasherTypes =
    ::testing::Types<Poseidon2MerkleHasher, Blake3sMerkleHasher, Sha256MerkleHasher, SkyscraperMerkleHasher>;
TYPED_TEST_SUITE(WhirMerkleTreeTest, HasherTypes);

namespace {
std::vector<fr> random_codeword(size_t size)
{
    std::vector<fr> codeword(size);
    for (fr& value : codeword) {
        value = fr::random_element();
    }
    return codeword;
}
} // namespace

TYPED_TEST(WhirMerkleTreeTest, LeavesGroupFoldCosets)
{
    constexpr size_t log_arity = 2;
    const std::vector<fr> codeword = random_codeword(64);
    MerkleTree<TypeParam> tree(codeword, log_arity);
    ASSERT_EQ(tree.num_leaves(), 16U);
    for (size_t j = 0; j < tree.num_leaves(); ++j) {
        const auto opening = tree.open(j);
        ASSERT_EQ(opening.values.size(), 4U);
        for (size_t t = 0; t < 4; ++t) {
            EXPECT_EQ(opening.values[t], codeword[j + t * tree.num_leaves()]);
        }
    }
}

TYPED_TEST(WhirMerkleTreeTest, OpenVerifyRoundTrip)
{
    const std::vector<fr> codeword = random_codeword(128);
    MerkleTree<TypeParam> tree(codeword, 3);
    for (size_t j = 0; j < tree.num_leaves(); ++j) {
        EXPECT_TRUE(MerkleTree<TypeParam>::verify(tree.root(), j, tree.open(j)));
    }
}

TYPED_TEST(WhirMerkleTreeTest, TamperingRejected)
{
    const std::vector<fr> codeword = random_codeword(64);
    MerkleTree<TypeParam> tree(codeword, 2);

    // Wrong leaf index
    EXPECT_FALSE(MerkleTree<TypeParam>::verify(tree.root(), 1, tree.open(0)));

    // Tampered value
    auto opening = tree.open(5);
    opening.values[2] += fr(1);
    EXPECT_FALSE(MerkleTree<TypeParam>::verify(tree.root(), 5, opening));

    // Tampered path
    opening = tree.open(5);
    opening.path[1] = TypeParam::hash_node(opening.path[1], opening.path[1]);
    EXPECT_FALSE(MerkleTree<TypeParam>::verify(tree.root(), 5, opening));

    // Wrong root
    const std::vector<fr> other = random_codeword(64);
    MerkleTree<TypeParam> other_tree(other, 2);
    EXPECT_FALSE(MerkleTree<TypeParam>::verify(other_tree.root(), 5, tree.open(5)));
}

TYPED_TEST(WhirMerkleTreeTest, SaltedCommitmentIsHidingAndVerifies)
{
    const std::vector<fr> codeword = random_codeword(64);
    MerkleTree<TypeParam> salted_a(codeword, 2, /*salted=*/true);
    MerkleTree<TypeParam> salted_b(codeword, 2, /*salted=*/true);

    // Two salted commitments to the same data differ (statistically hiding).
    EXPECT_NE(salted_a.root(), salted_b.root());

    const auto opening = salted_a.open(3);
    ASSERT_TRUE(opening.salt.has_value());
    EXPECT_TRUE(MerkleTree<TypeParam>::verify(salted_a.root(), 3, opening));

    // Dropping the salt invalidates the opening.
    auto stripped = opening;
    stripped.salt.reset();
    EXPECT_FALSE(MerkleTree<TypeParam>::verify(salted_a.root(), 3, stripped));
}

TYPED_TEST(WhirMerkleTreeTest, DigestFieldRoundTrip)
{
    const std::vector<fr> codeword = random_codeword(16);
    MerkleTree<TypeParam> tree(codeword, 2);
    const auto fields = TypeParam::digest_to_fields(tree.root());
    EXPECT_EQ(TypeParam::digest_from_fields(fields), tree.root());
}

} // namespace bb::whir
