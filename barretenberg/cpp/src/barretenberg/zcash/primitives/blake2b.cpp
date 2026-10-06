#include "blake2b.hpp"
#include "barretenberg/common/assert.hpp"

#include <algorithm>
#include <cstring>

namespace bb::zcash {

namespace {
constexpr std::array<uint64_t, 8> IV = { 0x6a09e667f3bcc908ULL, 0xbb67ae8584caa73bULL, 0x3c6ef372fe94f82bULL,
                                         0xa54ff53a5f1d36f1ULL, 0x510e527fade682d1ULL, 0x9b05688c2b3e6c1fULL,
                                         0x1f83d9abfb41bd6bULL, 0x5be0cd19137e2179ULL };

constexpr std::array<std::array<uint8_t, 16>, 12> SIGMA = { {
    { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 },
    { 14, 10, 4, 8, 9, 15, 13, 6, 1, 12, 0, 2, 11, 7, 5, 3 },
    { 11, 8, 12, 0, 5, 2, 15, 13, 10, 14, 3, 6, 7, 1, 9, 4 },
    { 7, 9, 3, 1, 13, 12, 11, 14, 2, 6, 5, 10, 4, 0, 15, 8 },
    { 9, 0, 5, 7, 2, 4, 10, 15, 14, 1, 11, 12, 6, 8, 3, 13 },
    { 2, 12, 6, 10, 0, 11, 8, 3, 4, 13, 7, 5, 15, 14, 1, 9 },
    { 12, 5, 1, 15, 14, 13, 4, 10, 0, 7, 6, 3, 9, 2, 8, 11 },
    { 13, 11, 7, 14, 12, 1, 3, 9, 5, 0, 15, 4, 8, 6, 2, 10 },
    { 6, 15, 14, 9, 11, 3, 0, 8, 12, 2, 13, 7, 1, 4, 10, 5 },
    { 10, 2, 8, 4, 7, 6, 1, 5, 15, 11, 9, 14, 3, 12, 13, 0 },
    { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 },
    { 14, 10, 4, 8, 9, 15, 13, 6, 1, 12, 0, 2, 11, 7, 5, 3 },
} };

constexpr uint64_t rotr(uint64_t x, int n)
{
    return (x >> n) | (x << (64 - n));
}

uint64_t load64(const uint8_t* p)
{
    uint64_t r = 0;
    for (size_t i = 0; i < 8; ++i) {
        r |= static_cast<uint64_t>(p[i]) << (8 * i);
    }
    return r;
}
} // namespace

Blake2b::Blake2b(size_t out_len, std::span<const uint8_t> personal)
    : out_len_(out_len)
{
    BB_ASSERT(out_len >= 1 && out_len <= 64);
    BB_ASSERT(personal.empty() || personal.size() == 16);
    h_ = IV;
    // parameter block: digest length, key length 0, fanout 1, depth 1; personalization in bytes 48..63
    h_[0] ^= 0x01010000ULL ^ static_cast<uint64_t>(out_len);
    if (!personal.empty()) {
        h_[6] ^= load64(personal.data());
        h_[7] ^= load64(personal.data() + 8);
    }
}

void Blake2b::compress(bool last)
{
    std::array<uint64_t, 16> m{};
    for (size_t i = 0; i < 16; ++i) {
        m[i] = load64(&buf_[8 * i]);
    }
    std::array<uint64_t, 16> v{};
    for (size_t i = 0; i < 8; ++i) {
        v[i] = h_[i];
        v[i + 8] = IV[i];
    }
    v[12] ^= t0_;
    v[13] ^= t1_;
    if (last) {
        v[14] = ~v[14];
    }
    auto g = [&](size_t a, size_t b, size_t c, size_t d, uint64_t x, uint64_t y) {
        v[a] = v[a] + v[b] + x;
        v[d] = rotr(v[d] ^ v[a], 32);
        v[c] = v[c] + v[d];
        v[b] = rotr(v[b] ^ v[c], 24);
        v[a] = v[a] + v[b] + y;
        v[d] = rotr(v[d] ^ v[a], 16);
        v[c] = v[c] + v[d];
        v[b] = rotr(v[b] ^ v[c], 63);
    };
    for (const auto& s : SIGMA) {
        g(0, 4, 8, 12, m[s[0]], m[s[1]]);
        g(1, 5, 9, 13, m[s[2]], m[s[3]]);
        g(2, 6, 10, 14, m[s[4]], m[s[5]]);
        g(3, 7, 11, 15, m[s[6]], m[s[7]]);
        g(0, 5, 10, 15, m[s[8]], m[s[9]]);
        g(1, 6, 11, 12, m[s[10]], m[s[11]]);
        g(2, 7, 8, 13, m[s[12]], m[s[13]]);
        g(3, 4, 9, 14, m[s[14]], m[s[15]]);
    }
    for (size_t i = 0; i < 8; ++i) {
        h_[i] ^= v[i] ^ v[i + 8];
    }
}

Blake2b& Blake2b::update(std::span<const uint8_t> data)
{
    size_t offset = 0;
    while (offset < data.size()) {
        // A full buffer is only compressed once more input arrives: the final block needs the `last` flag.
        if (buf_len_ == buf_.size()) {
            t0_ += buf_.size();
            if (t0_ < buf_.size()) {
                ++t1_;
            }
            compress(false);
            buf_len_ = 0;
        }
        const size_t take = std::min(buf_.size() - buf_len_, data.size() - offset);
        std::memcpy(&buf_[buf_len_], data.data() + offset, take);
        buf_len_ += take;
        offset += take;
    }
    return *this;
}

std::array<uint8_t, 64> Blake2b::finalize()
{
    t0_ += buf_len_;
    if (t0_ < buf_len_) {
        ++t1_;
    }
    std::fill(buf_.begin() + static_cast<std::ptrdiff_t>(buf_len_), buf_.end(), 0);
    compress(true);
    std::array<uint8_t, 64> out{};
    for (size_t i = 0; i < 64; ++i) {
        out[i] = static_cast<uint8_t>(h_[i / 8] >> (8 * (i % 8)));
    }
    for (size_t i = out_len_; i < 64; ++i) {
        out[i] = 0;
    }
    return out;
}

} // namespace bb::zcash
