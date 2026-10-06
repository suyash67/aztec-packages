#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <string_view>

namespace bb::zcash {

/**
 * @brief Incremental BLAKE2b (RFC 7693) with optional 16-byte personalization.
 * @details Used by the Pasta hash-to-curve (hash_to_field with BLAKE2b-512 XMD) that defines every Orchard
 * generator, and by halo2's parameter generation. Not constant-time with respect to message length.
 */
class Blake2b {
  public:
    explicit Blake2b(size_t out_len = 64, std::span<const uint8_t> personal = {});

    Blake2b& update(std::span<const uint8_t> data);
    Blake2b& update(std::string_view data)
    {
        return update(std::span<const uint8_t>(reinterpret_cast<const uint8_t*>(data.data()), data.size()));
    }
    Blake2b& update(uint8_t byte) { return update(std::span<const uint8_t>(&byte, 1)); }

    // Returns the digest in the first `out_len` bytes.
    std::array<uint8_t, 64> finalize();

  private:
    void compress(bool last);

    std::array<uint64_t, 8> h_{};
    std::array<uint8_t, 128> buf_{};
    size_t buf_len_ = 0;
    uint64_t t0_ = 0;
    uint64_t t1_ = 0;
    size_t out_len_;
};

} // namespace bb::zcash
