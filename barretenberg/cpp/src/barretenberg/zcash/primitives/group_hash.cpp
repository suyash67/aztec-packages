#include "group_hash.hpp"
#include "barretenberg/common/assert.hpp"
#include "blake2b.hpp"

#include <array>
#include <optional>

namespace bb::zcash {

namespace {

template <typename Fq> struct SswuParams {
    Fq iso_a;
    Fq iso_b;
    Fq z;
    std::array<Fq, 13> iso;
    const char* curve_id;
};

// Constants from pasta_curves 0.6 (curves.rs): IsoEp / IsoEq coefficients, the SSWU Z, and the isogeny map.
SswuParams<pallas::fq> pallas_params()
{
    using F = pallas::fq;
    return {
        F(uint256_t("0x18354a2eb0ea8c9c49be2d7258370742b74134581a27a59f92bb4b0b657a014b")),
        F(1265),
        F(uint256_t("0x40000000000000000000000000000000224698fc094cf91b992d30ecfffffff4")),
        { F(uint256_t("0x0e38e38e38e38e38e38e38e38e38e38e4081775473d8375b775f6034aaaaaaab")),
          F(uint256_t("0x3509afd51872d88e267c7ffa51cf412a0f93b82ee4b994958cf863b02814fb76")),
          F(uint256_t("0x17329b9ec525375398c7d7ac3d98fd13380af066cfeb6d690eb64faef37ea4f7")),
          F(uint256_t("0x1c71c71c71c71c71c71c71c71c71c71c8102eea8e7b06eb6eebec06955555580")),
          F(uint256_t("0x1d572e7ddc099cff5a607fcce0494a799c434ac1c96b6980c47f2ab668bcd71f")),
          F(uint256_t("0x325669becaecd5d11d13bf2a7f22b105b4abf9fb9a1fc81c2aa3af1eae5b6604")),
          F(uint256_t("0x1a12f684bda12f684bda12f684bda12f7642b01ad461bad25ad985b5e38e38e4")),
          F(uint256_t("0x1a84d7ea8c396c47133e3ffd28e7a09507c9dc17725cca4ac67c31d8140a7dbb")),
          F(uint256_t("0x3fb98ff0d2ddcadd303216cce1db9ff11765e924f745937802e2be87d225b234")),
          F(uint256_t("0x025ed097b425ed097b425ed097b425ed0ac03e8e134eb3e493e53ab371c71c4f")),
          F(uint256_t("0x0c02c5bcca0e6b7f0790bfb3506defb65941a3a4a97aa1b35a28279b1d1b42ae")),
          F(uint256_t("0x17033d3c60c68173573b3d7f7d681310d976bbfabbc5661d4d90ab820b12320a")),
          F(uint256_t("0x40000000000000000000000000000000224698fc094cf91b992d30ecfffffde5")) },
        "pallas",
    };
}

SswuParams<vesta::fq> vesta_params()
{
    using F = vesta::fq;
    return {
        F(uint256_t("0x267f9b2ee592271a81639c4d96f787739673928c7d01b212c515ad7242eaa6b1")),
        F(1265),
        F(uint256_t("0x40000000000000000000000000000000224698fc0994a8dd8c46eb20fffffff4")),
        { F(uint256_t("0x38e38e38e38e38e38e38e38e38e38e390205dd51cfa0961a43cd42c800000001")),
          F(uint256_t("0x1d935247b4473d17acecf10f5f7c09a2216b8861ec72bd5d8b95c6aaf703bcc5")),
          F(uint256_t("0x18760c7f7a9ad20ded7ee4a9cdf78f8fd59d03d23b39cb11aeac67bbeb586a3d")),
          F(uint256_t("0x31c71c71c71c71c71c71c71c71c71c71e1c521a795ac8356fb539a6f0000002b")),
          F(uint256_t("0x0a2de485568125d51454798a5b5c56b2a3ad678129b604d3b7284f7eaf21a2e9")),
          F(uint256_t("0x14735171ee5427780c621de8b91c242a30cd6d53df49d235f169c187d2533465")),
          F(uint256_t("0x12f684bda12f684bda12f684bda12f685601f4709a8adcb36bef1642aaaaaaab")),
          F(uint256_t("0x2ec9a923da239e8bd6767887afbe04d121d910aefb03b31d8bee58e5fb81de63")),
          F(uint256_t("0x19b0d87e16e2578866d1466e9de10e6497a3ca5c24e9ea634986913ab4443034")),
          F(uint256_t("0x1ed097b425ed097b425ed097b425ed098bc32d36fb21a6a38f64842c55555533")),
          F(uint256_t("0x2f44d6c801c1b8bf9e7eb64f890a820c06a767bfc35b5bac58dfecce86b2745e")),
          F(uint256_t("0x3d59f455cafc7668252659ba2b546c7e926847fb9ddd76a1d43d449776f99d2f")),
          F(uint256_t("0x40000000000000000000000000000000224698fc0994a8dd8c46eb20fffffde5")) },
        "vesta",
    };
}

template <typename F> bool is_odd(const F& x)
{
    return uint256_t(x).get_bit(0);
}

template <typename F> F from_uniform_bytes_be(const std::array<uint8_t, 64>& be)
{
    // pasta_curves reverses the big-endian hash output and reduces the resulting 512-bit little-endian integer
    uint512_t acc = 0;
    for (uint8_t b : be) {
        acc = (acc << 8) + uint512_t(b);
    }
    return F((acc % uint512_t(F::modulus)).lo);
}

template <typename F>
std::array<F, 2> hash_to_field(const char* curve_id, std::string_view domain_prefix, std::span<const uint8_t> message)
{
    constexpr size_t CHUNKLEN = 64;
    const std::string_view cid(curve_id);
    BB_ASSERT_LT(22 + cid.size() + domain_prefix.size(), size_t{ 256 });
    const auto dst_len = static_cast<uint8_t>(22 + cid.size() + domain_prefix.size());
    auto suffix = [&](Blake2b& h) {
        h.update(domain_prefix).update("-").update(cid).update("_XMD:BLAKE2b_SSWU_RO_").update(dst_len);
    };

    std::array<uint8_t, 128> zeros{};
    Blake2b h0(CHUNKLEN);
    h0.update(zeros).update(message);
    const std::array<uint8_t, 3> lib_str = { 0, CHUNKLEN * 2, 0 };
    h0.update(lib_str);
    suffix(h0);
    const auto b0 = h0.finalize();

    Blake2b h1(CHUNKLEN);
    h1.update(b0).update(uint8_t{ 1 });
    suffix(h1);
    const auto b1 = h1.finalize();

    std::array<uint8_t, 64> x{};
    for (size_t i = 0; i < 64; ++i) {
        x[i] = b0[i] ^ b1[i];
    }
    Blake2b h2(CHUNKLEN);
    h2.update(x).update(uint8_t{ 2 });
    suffix(h2);
    const auto b2 = h2.finalize();

    return { from_uniform_bytes_be<F>(b1), from_uniform_bytes_be<F>(b2) };
}

// Affine points on the isogenous curve y^2 = x^3 + a x + b; nullopt is the identity.
template <typename F> using IsoPoint = std::optional<std::pair<F, F>>;

template <typename F> IsoPoint<F> map_to_curve_simple_swu(const F& u, const SswuParams<F>& p)
{
    const F z_u2 = p.z * u.sqr();
    const F ta = z_u2.sqr() + z_u2;
    const F num_x1 = p.iso_b * (ta + F::one());
    const F div = p.iso_a * (ta.is_zero() ? p.z : -ta);
    const F x1 = num_x1 * div.invert();
    const F gx1 = (x1.sqr() + p.iso_a) * x1 + p.iso_b;
    auto [gx1_square, y1] = gx1.sqrt();
    F x = x1;
    F y = y1;
    if (!gx1_square) {
        x = z_u2 * x1;
        const F gx2 = (x.sqr() + p.iso_a) * x + p.iso_b;
        auto [gx2_square, y2] = gx2.sqrt();
        BB_ASSERT(gx2_square, "SSWU: g(x2) must be a square when g(x1) is not");
        y = y2;
    }
    if (is_odd(u) != is_odd(y)) {
        y = -y;
    }
    return std::make_pair(x, y);
}

template <typename F> IsoPoint<F> iso_add(const IsoPoint<F>& P, const IsoPoint<F>& Q, const F& a)
{
    if (!P) {
        return Q;
    }
    if (!Q) {
        return P;
    }
    const auto& [x1, y1] = *P;
    const auto& [x2, y2] = *Q;
    F lambda;
    if (x1 == x2) {
        if ((y1 + y2).is_zero()) {
            return std::nullopt;
        }
        lambda = (x1.sqr() * F(3) + a) * (y1 + y1).invert();
    } else {
        lambda = (y2 - y1) * (x2 - x1).invert();
    }
    const F x3 = lambda.sqr() - x1 - x2;
    return std::make_pair(x3, lambda * (x1 - x3) - y1);
}

template <typename Group, typename F>
typename Group::affine_element hash_to_curve_impl(const SswuParams<F>& p,
                                                  std::string_view domain_prefix,
                                                  std::span<const uint8_t> message)
{
    const auto us = hash_to_field<F>(p.curve_id, domain_prefix, message);
    const auto r = iso_add(map_to_curve_simple_swu(us[0], p), map_to_curve_simple_swu(us[1], p), p.iso_a);
    if (!r) {
        return Group::affine_point_at_infinity;
    }
    const auto& [x, y] = *r;
    const auto& iso = p.iso;
    const F num_x = ((iso[0] * x + iso[1]) * x + iso[2]) * x + iso[3];
    const F div_x = (x + iso[4]) * x + iso[5];
    const F num_y = (((iso[6] * x + iso[7]) * x + iso[8]) * x + iso[9]) * y;
    const F div_y = ((x + iso[10]) * x + iso[11]) * x + iso[12];
    if (div_x.is_zero() || div_y.is_zero()) {
        return Group::affine_point_at_infinity;
    }
    typename Group::affine_element out(num_x * div_x.invert(), num_y * div_y.invert());
    BB_ASSERT(out.on_curve());
    return out;
}

} // namespace

pallas::g1::affine_element pallas_hash_to_curve(std::string_view domain_prefix, std::span<const uint8_t> message)
{
    static const auto params = pallas_params();
    return hash_to_curve_impl<pallas::g1>(params, domain_prefix, message);
}

vesta::g1::affine_element vesta_hash_to_curve(std::string_view domain_prefix, std::span<const uint8_t> message)
{
    static const auto params = vesta_params();
    return hash_to_curve_impl<vesta::g1>(params, domain_prefix, message);
}

} // namespace bb::zcash
