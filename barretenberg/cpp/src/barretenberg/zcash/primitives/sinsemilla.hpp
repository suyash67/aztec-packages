#pragma once

#include "barretenberg/common/assert.hpp"
#include "cycle.hpp"

#include <array>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace bb::zcash {

/**
 * @brief Native Sinsemilla (Zcash protocol spec §5.4.1.9) over the embedded curve of `Cycle`.
 * @details
 *   S(m) = GroupHash("z.cash:SinsemillaS", I2LEOSP_32(m)) for m in [0, 2^10),
 *   Q(D) = GroupHash("z.cash:SinsemillaQ", D),
 *   HashToPoint(D, M): A = Q(D); for each 10-bit word m_i of M (zero-padded): A = (A + S(m_i)) + A,
 *   with incomplete additions that return ⊥ on any exceptional case (equal x-coordinates or the identity),
 *   Commit_r(D, M) = HashToPoint(D || "-M", M) + [r] GroupHash(D || "-r", "").
 */
template <typename Cycle> class Sinsemilla {
  public:
    using FF = typename Cycle::FF;
    using AffineElement = typename Cycle::AffineElement;
    using Element = typename Cycle::Element;
    using Scalar = typename Cycle::EmbeddedScalar;

    static constexpr size_t K = 10;
    static constexpr size_t TABLE_SIZE = size_t{ 1 } << K;
    // Maximum number of words in one message (Sinsemilla's C parameter).
    static constexpr size_t C = 253;

    static const std::vector<AffineElement>& S_table()
    {
        static const std::vector<AffineElement> table = []() {
            std::vector<AffineElement> t(TABLE_SIZE);
            for (size_t m = 0; m < TABLE_SIZE; ++m) {
                const std::array<uint8_t, 4> le = { static_cast<uint8_t>(m),
                                                    static_cast<uint8_t>(m >> 8),
                                                    static_cast<uint8_t>(m >> 16),
                                                    static_cast<uint8_t>(m >> 24) };
                t[m] = Cycle::group_hash("z.cash:SinsemillaS", le);
            }
            return t;
        }();
        return table;
    }

    static AffineElement Q(std::string_view domain) { return group_hash<Cycle>("z.cash:SinsemillaQ", domain); }

    static AffineElement R(std::string_view commit_domain)
    {
        return group_hash<Cycle>(std::string(commit_domain) + "-r", "");
    }

    // Incomplete affine addition; nullopt is ⊥.
    static std::optional<AffineElement> add_incomplete(const AffineElement& a, const AffineElement& b)
    {
        if (a.is_point_at_infinity() || b.is_point_at_infinity() || a.x == b.x) {
            return std::nullopt;
        }
        const FF lambda = (b.y - a.y) * (b.x - a.x).invert();
        const FF x3 = lambda.sqr() - a.x - b.x;
        return AffineElement(x3, lambda * (a.x - x3) - a.y);
    }

    // Words of a little-endian bit string, zero-padded to a multiple of K bits.
    static std::vector<uint32_t> words(const std::vector<bool>& bits)
    {
        std::vector<uint32_t> out((bits.size() + K - 1) / K, 0);
        for (size_t i = 0; i < bits.size(); ++i) {
            if (bits[i]) {
                out[i / K] |= uint32_t{ 1 } << (i % K);
            }
        }
        return out;
    }

    static std::optional<AffineElement> hash_to_point(const AffineElement& q, const std::vector<bool>& bits)
    {
        const auto ws = words(bits);
        BB_ASSERT_LTE(ws.size(), C);
        const auto& table = S_table();
        std::optional<AffineElement> acc = q;
        for (uint32_t w : ws) {
            auto tmp = add_incomplete(*acc, table[w]);
            if (!tmp) {
                return std::nullopt;
            }
            acc = add_incomplete(*tmp, *acc);
            if (!acc) {
                return std::nullopt;
            }
        }
        return acc;
    }

    static std::optional<AffineElement> commit(std::string_view domain, const std::vector<bool>& bits, const Scalar& r)
    {
        auto h = hash_to_point(Q(std::string(domain) + "-M"), bits);
        if (!h) {
            return std::nullopt;
        }
        return AffineElement(Element(*h) + Element(R(domain)) * r);
    }

    static std::optional<FF> short_commit(std::string_view domain, const std::vector<bool>& bits, const Scalar& r)
    {
        auto c = commit(domain, bits, r);
        if (!c || c->is_point_at_infinity()) {
            return std::nullopt;
        }
        return c->x;
    }
};

// I2LEBSP_l(x): the l least-significant bits of the canonical integer x, little-endian.
template <typename T> void append_bits(std::vector<bool>& out, const T& value, size_t num_bits)
{
    const uint256_t v(value);
    for (size_t i = 0; i < num_bits; ++i) {
        out.push_back(v.get_bit(i));
    }
}

} // namespace bb::zcash
