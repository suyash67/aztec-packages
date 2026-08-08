#include "barretenberg/fflonk/proof.hpp"

#include "barretenberg/fflonk/encoding.hpp"

namespace bb::fflonk_plonk {

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
        if (!read_scalar(buffer, offset, evaluation)) {
            return false;
        }
    }

    return true;
}

} // namespace bb::fflonk_plonk
