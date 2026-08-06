#include "barretenberg/crypto/skyscraper/skyscraper.hpp"

#include "barretenberg/numeric/uint256/uint256.hpp"

#include <gtest/gtest.h>

// Vectors from ProveKit's reference implementation (worldfnd/ProveKit
// skyscraper/core/src/reference.rs), which in turn come from the paper's sage notebook.

namespace bb::crypto::skyscraper {
namespace {

fr from_hex(const char* hex)
{
    return fr(uint256_t(hex));
}

TEST(Skyscraper, SboxTable3)
{
    EXPECT_EQ(sbox(0xcd), 0xd3);
    EXPECT_EQ(sbox(0x17), 0x0e);
    EXPECT_EQ(sbox(0x83), 0x17);
    EXPECT_EQ(sbox(0x14), 0x28);
    EXPECT_EQ(sbox(0x2b), 0x46);
    EXPECT_EQ(sbox(0x1e), 0xbc);
}

TEST(Skyscraper, PermuteZero)
{
    const auto [l, r] = permute(fr(0), fr(0));
    EXPECT_EQ(l, from_hex("0x0ccee0e750cacbe110ab2b912d9cd38f0a4a74dbc4fa4bbcc2d3218600b3f9ea"));
    EXPECT_EQ(r, from_hex("0x1b2f71d974b15a2eccf059f57022bca6ffae279d81831a0884d26a76d2307925"));
}

TEST(Skyscraper, PermuteRandomVector)
{
    // The input left value exceeds r; ark's parser reduces it mod r, as does fr(uint256_t).
    const fr in_l = from_hex("0x6f7721ff66a1725a6647d22c3a9032b91f2d82e3bf61a6f5a88ac1c1df0de2f4");
    const fr in_r = from_hex("0x205325dcd29fb570ae478e12273840597b0d9adf8b76f6c8ed4ac3d9f1d8db4e");
    const auto [l, r] = permute(in_l, in_r);
    EXPECT_EQ(l, from_hex("0x12998f99c09d1c18162041642fd35a0b31cfdf560bc6ee14fa841165cb51664e"));
    EXPECT_EQ(r, from_hex("0x1a3d2642c9398e9bef8a84e5ede238a1fd395f9351be64ab377ecb11a0660fef"));
}

TEST(Skyscraper, CompressIsLeftFeedForward)
{
    const fr l = fr(5);
    const fr r = fr(7);
    EXPECT_EQ(compress(l, r), permute(l, r).first + l);
}

TEST(Skyscraper, FoldMatchesReduceSemantics)
{
    const std::vector<fr> values = { fr(1), fr(2), fr(3), fr(4) };
    const fr expected = compress(compress(compress(fr(1), fr(2)), fr(3)), fr(4));
    EXPECT_EQ(fold_compress(values), expected);
    EXPECT_EQ(fold_compress(std::span<const fr>(values.data(), 1)), fr(1));
}

} // namespace
} // namespace bb::crypto::skyscraper
