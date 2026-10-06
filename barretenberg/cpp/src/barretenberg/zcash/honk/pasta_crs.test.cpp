#include "pasta_crs.hpp"

#include <gtest/gtest.h>

using namespace bb;
using namespace bb::zcash;

namespace {
// halo2 / pasta_curves compressed encoding: x little-endian, sign of y in the top bit.
std::string compress(const vesta::g1::affine_element& p)
{
    uint256_t x(p.x);
    if (uint256_t(p.y).get_bit(0)) {
        x |= uint256_t(1) << 255;
    }
    std::string out;
    static constexpr char HEX[] = "0123456789abcdef";
    for (size_t i = 0; i < 32; ++i) {
        const auto byte = static_cast<uint8_t>(x.data[i / 8] >> (8 * (i % 8)));
        out += HEX[byte >> 4];
        out += HEX[byte & 0xf];
    }
    return out;
}
} // namespace

TEST(ZcashPastaCrs, MatchesHalo2Params)
{
    // halo2_proofs 0.4 Params::<vesta::Affine>::new(4), see scripts/zcash_orchard_reference/src/bin/params.rs
    const auto g = halo2_vesta_generators(16);
    EXPECT_EQ(compress(g[0]), "45065ed079bf389758f591131095ef419310e8c708a805852b9b77bed8c7ecbd");
    EXPECT_EQ(compress(g[1]), "e0c0802686d3ed571f7f3399526b24460b16ace461ebda9dcfe6e5b7b298c18c");
    EXPECT_EQ(compress(g[15]), "a370615e5b383d342f401237f6c2b72573bf10fafbd4f41c48857e4aef63cb81");
    EXPECT_EQ(compress(halo2_vesta_w()), "7520d96f3e5cd41760367151608b54821883c10c4b9a4ff2beae227bef94bcab");
    EXPECT_EQ(compress(halo2_vesta_u()), "379dc4dcfdbf61ccc7d5a0bb9759acf611694f24d0c040f249bad30a83b1a897");
}
