#include "group_hash.hpp"
#include "barretenberg/zcash/test_vectors/orchard_vectors.hpp"
#include "blake2b.hpp"

#include <gtest/gtest.h>

using namespace bb;
using namespace bb::zcash;

namespace {
std::string hex(std::span<const uint8_t> bytes)
{
    static const char* digits = "0123456789abcdef";
    std::string s;
    for (uint8_t b : bytes) {
        s.push_back(digits[b >> 4]);
        s.push_back(digits[b & 15]);
    }
    return s;
}
} // namespace

TEST(ZcashBlake2b, Rfc7693Vectors)
{
    auto abc = Blake2b(64).update("abc").finalize();
    EXPECT_EQ(
        hex(abc),
        "ba80a53f981c4d0d6a2797b69f12f6e94c212f14685ac4b74b12bb6fdbffa2d17d87c5392aab792dc252d5de4533cc9518d38aa8dbf"
        "1925ab92386edd4009923");
    auto empty = Blake2b(64).finalize();
    EXPECT_EQ(
        hex(empty),
        "786a02f742015903c6c6fd852552d272912f4740e15847618a86e217f71f5419d25e1031afee585313896444934eb04b903a685b1448"
        "b755d56f701afe9be2ce");
    // Incremental updates across the 128-byte block boundary match a single update.
    std::vector<uint8_t> msg(300);
    for (size_t i = 0; i < msg.size(); ++i) {
        msg[i] = static_cast<uint8_t>(i * 7);
    }
    auto one_shot = Blake2b(64).update(msg).finalize();
    Blake2b h(64);
    h.update(std::span<const uint8_t>(msg).subspan(0, 128));
    h.update(std::span<const uint8_t>(msg).subspan(128, 1));
    h.update(std::span<const uint8_t>(msg).subspan(129));
    EXPECT_EQ(hex(h.finalize()), hex(one_shot));
}

TEST(ZcashGroupHash, MatchesPastaCurvesAndOrchard)
{
    // Every vector was produced by pasta_curves' `pallas::Point::hash_to_curve(domain)(message)`; they cover the
    // pasta_curves test vector and every Orchard generator (SpendAuthG, NullifierK, ValueCommit V/R, commitment
    // randomness bases, Sinsemilla Q points, sample Sinsemilla S points).
    for (const auto& v : test_vectors::GROUP_HASH_VECTORS) {
        auto p = pallas_hash_to_curve(v.domain, v.message);
        EXPECT_EQ(uint256_t(p.x), uint256_t(v.point.x)) << v.domain;
        EXPECT_EQ(uint256_t(p.y), uint256_t(v.point.y)) << v.domain;
    }
}

TEST(ZcashGroupHash, PastaCurvesJacobianTestVector)
{
    // pasta_curves pallas.rs test_hash_to_curve, given in Jacobian coordinates (X, Y, Z) -> (X / Z^2, Y / Z^3).
    const pallas::fq X(uint256_t("0x36a6e3a9c50b7b6540cb002c977c82f37f8a875fb51eb35327ee1452e6ce7947"));
    const pallas::fq Y(uint256_t("0x01da3b4403d73252f2d7e9c19bc23dc6a080f2d02f8262fca4f7e3d756ac6a7c"));
    const pallas::fq Z(uint256_t("0x1d48103df8fcbb70d1809c1806c95651dd884a559fec0549658537ce9d94bed9"));
    auto p = pallas_hash_to_curve("z.cash:test", "Trans rights now!");
    const auto zinv = Z.invert();
    EXPECT_EQ(p.x, X * zinv.sqr());
    EXPECT_EQ(p.y, Y * zinv.sqr() * zinv);
}

TEST(ZcashGroupHash, VestaMatchesPastaCurves)
{
    const std::string_view msg = "hello";
    auto p = vesta_hash_to_curve("z.cash:test",
                                 std::span<const uint8_t>(reinterpret_cast<const uint8_t*>(msg.data()), msg.size()));
    EXPECT_EQ(uint256_t(p.x), uint256_t(test_vectors::VESTA_HASH_TO_CURVE_HELLO.x));
    EXPECT_EQ(uint256_t(p.y), uint256_t(test_vectors::VESTA_HASH_TO_CURVE_HELLO.y));
}
