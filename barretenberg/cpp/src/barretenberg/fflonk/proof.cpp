#include "barretenberg/fflonk/proof.hpp"

namespace bb::fflonk_plonk {

namespace {

void append_word(std::vector<uint8_t>& buffer, const uint256_t& word)
{
    for (size_t limb = 4; limb-- > 0;) {
        for (size_t byte = 8; byte-- > 0;) {
            buffer.push_back(static_cast<uint8_t>((word.data[limb] >> (8 * byte)) & 0xff));
        }
    }
}

uint256_t read_word(std::span<const uint8_t> buffer, size_t offset)
{
    uint256_t word{ 0 };
    for (size_t limb = 4; limb-- > 0;) {
        uint64_t value = 0;
        for (size_t byte = 0; byte < 8; ++byte) {
            value = (value << 8) | buffer[offset++];
        }
        word.data[limb] = value;
    }
    return word;
}

void append_point(std::vector<uint8_t>& buffer, const Commitment& point)
{
    if (point.is_point_at_infinity()) {
        append_word(buffer, uint256_t(0));
        append_word(buffer, uint256_t(0));
    } else {
        append_word(buffer, static_cast<uint256_t>(point.x));
        append_word(buffer, static_cast<uint256_t>(point.y));
    }
}

bool read_point(std::span<const uint8_t> buffer, size_t& offset, Commitment& point)
{
    const uint256_t x = read_word(buffer, offset);
    const uint256_t y = read_word(buffer, offset + 32);
    offset += 64;

    // Reject non-canonical coordinates before building the field elements: `fq`'s constructor would
    // silently reduce them, which would let one point have many encodings and break the binding
    // between a proof's bytes and its transcript.
    if (x >= Curve::BaseField::modulus || y >= Curve::BaseField::modulus) {
        return false;
    }

    point = Commitment(Curve::BaseField(x), Curve::BaseField(y));
    return point.on_curve() && !point.is_point_at_infinity();
}

} // namespace

std::vector<uint8_t> Proof::to_buffer() const
{
    std::vector<uint8_t> buffer;
    buffer.reserve(SIZE_IN_BYTES);
    append_point(buffer, c1);
    append_point(buffer, c2);
    append_point(buffer, c3);
    append_point(buffer, w);
    append_point(buffer, w_prime);
    for (const FF& evaluation : evaluations) {
        append_word(buffer, static_cast<uint256_t>(evaluation));
    }
    return buffer;
}

bool Proof::from_buffer(std::span<const uint8_t> buffer, Proof& proof)
{
    if (buffer.size() != SIZE_IN_BYTES) {
        return false;
    }

    size_t offset = 0;
    for (Commitment* point : { &proof.c1, &proof.c2, &proof.c3, &proof.w, &proof.w_prime }) {
        if (!read_point(buffer, offset, *point)) {
            return false;
        }
    }

    for (FF& evaluation : proof.evaluations) {
        const uint256_t word = read_word(buffer, offset);
        offset += 32;
        if (word >= FF::modulus) {
            return false;
        }
        evaluation = FF(word);
    }

    return true;
}

} // namespace bb::fflonk_plonk
