#include "barretenberg/commitment_schemes/whir/merkle_tree.hpp"

#include <gtest/gtest.h>
#include <numeric>

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

TYPED_TEST(WhirMerkleTreeTest, BatchOpeningAuthenticatesEveryLeaf)
{
    const std::vector<fr> codeword = random_codeword(1024);
    MerkleTree<TypeParam> tree(codeword, 2);
    using Tree = MerkleTree<TypeParam>;

    // Repeated and unsorted indices address the same distinct leaf set.
    const std::vector<size_t> indices = { 200, 7, 6, 201, 7, 0, 255 };
    const std::vector<size_t> leaves = Tree::batch_leaves(indices);
    EXPECT_EQ(leaves, (std::vector<size_t>{ 0, 6, 7, 200, 201, 255 }));

    const auto batch = tree.open_batch(indices);
    EXPECT_EQ(batch.values.size(), leaves.size());
    EXPECT_EQ(batch.siblings.size(), Tree::batch_num_siblings(leaves, tree.depth()));
    EXPECT_TRUE(Tree::verify_batch(tree.root(), leaves, tree.depth(), batch));

    // Each batched leaf carries the same coset the single-leaf opening does.
    for (size_t i = 0; i < leaves.size(); ++i) {
        EXPECT_EQ(batch.values[i], tree.open(leaves[i]).values);
    }

    // Sibling pairs (6,7) and (200,201) are derivable from each other, so a batch must be strictly
    // cheaper than the independent paths it replaces.
    EXPECT_LT(batch.siblings.size(), leaves.size() * tree.depth());
}

TYPED_TEST(WhirMerkleTreeTest, BatchOpeningTamperingRejected)
{
    const std::vector<fr> codeword = random_codeword(256);
    MerkleTree<TypeParam> tree(codeword, 2);
    using Tree = MerkleTree<TypeParam>;

    const std::vector<size_t> indices = { 3, 17, 40, 41 };
    const std::vector<size_t> leaves = Tree::batch_leaves(indices);
    ASSERT_TRUE(Tree::verify_batch(tree.root(), leaves, tree.depth(), tree.open_batch(indices)));

    auto tampered = tree.open_batch(indices);
    tampered.values[1][0] += fr(1);
    EXPECT_FALSE(Tree::verify_batch(tree.root(), leaves, tree.depth(), tampered));

    tampered = tree.open_batch(indices);
    tampered.siblings[0] = TypeParam::hash_node(tampered.siblings[0], tampered.siblings[0]);
    EXPECT_FALSE(Tree::verify_batch(tree.root(), leaves, tree.depth(), tampered));

    // Claiming a different leaf set for the same data must fail.
    EXPECT_FALSE(
        Tree::verify_batch(tree.root(), std::vector<size_t>{ 3, 17, 40, 42 }, tree.depth(), tree.open_batch(indices)));

    // A batch that omits or adds a sibling is rejected rather than silently accepted.
    tampered = tree.open_batch(indices);
    tampered.siblings.pop_back();
    EXPECT_FALSE(Tree::verify_batch(tree.root(), leaves, tree.depth(), tampered));
    tampered = tree.open_batch(indices);
    tampered.siblings.push_back(tree.root());
    EXPECT_FALSE(Tree::verify_batch(tree.root(), leaves, tree.depth(), tampered));
}

// Every leaf of the tree opened at once needs no siblings at all: the verifier rebuilds the whole
// tree from the leaves. This pins the walk's accounting at the extreme where sharing is total.
TYPED_TEST(WhirMerkleTreeTest, BatchOpeningOfEveryLeafSendsNoSiblings)
{
    const std::vector<fr> codeword = random_codeword(64);
    MerkleTree<TypeParam> tree(codeword, 2);
    using Tree = MerkleTree<TypeParam>;

    std::vector<size_t> all(tree.num_leaves());
    std::iota(all.begin(), all.end(), size_t(0));
    const auto batch = tree.open_batch(all);
    EXPECT_TRUE(batch.siblings.empty());
    EXPECT_TRUE(Tree::verify_batch(tree.root(), all, tree.depth(), batch));
}

// Leaves that all sit in one subtree still have to be climbed to the actual root: the walk is
// driven by the tree's depth, not by "one node left". Leaf 0 alone is the extreme case — its
// digest is already index 0 at every level, so a walk that stopped there would authenticate
// nothing at all.
TYPED_TEST(WhirMerkleTreeTest, BatchOpeningClimbsToTheRootFromAnySubtree)
{
    const std::vector<fr> codeword = random_codeword(256);
    MerkleTree<TypeParam> tree(codeword, 2);
    using Tree = MerkleTree<TypeParam>;

    for (const std::vector<size_t>& indices :
         std::vector<std::vector<size_t>>{ { 0 }, { 0, 1 }, { 0, 1, 2, 3 }, { 6, 7 }, { 63 }, { 0, 63 } }) {
        const std::vector<size_t> leaves = Tree::batch_leaves(indices);
        EXPECT_TRUE(Tree::verify_batch(tree.root(), leaves, tree.depth(), tree.open_batch(indices)))
            << "leaves starting at " << indices.front();
        // The same opening must not authenticate against a foreign root.
        MerkleTree<TypeParam> other(random_codeword(256), 2);
        EXPECT_FALSE(Tree::verify_batch(other.root(), leaves, tree.depth(), tree.open_batch(indices)));
    }
}

TYPED_TEST(WhirMerkleTreeTest, DigestFieldRoundTrip)
{
    const std::vector<fr> codeword = random_codeword(16);
    MerkleTree<TypeParam> tree(codeword, 2);
    const auto fields = TypeParam::digest_to_fields(tree.root());
    EXPECT_EQ(TypeParam::digest_from_fields(fields), tree.root());
}

} // namespace bb::whir
