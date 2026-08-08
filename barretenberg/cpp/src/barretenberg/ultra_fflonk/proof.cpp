#include "barretenberg/ultra_fflonk/proof.hpp"

#include "barretenberg/fflonk/encoding.hpp"

namespace bb::ultra_fflonk {

std::vector<uint8_t> Proof::to_buffer() const
{
    std::vector<uint8_t> buffer;
    buffer.reserve(SIZE_IN_BYTES);
    for (const Commitment& point : { wires, memory, grand_product, quotient, w, w_prime }) {
        fflonk_plonk::append_point(buffer, point);
    }
    for (const FF& evaluation : evaluations) {
        fflonk_plonk::append_word(buffer, static_cast<uint256_t>(evaluation));
    }
    return buffer;
}

bool Proof::from_buffer(std::span<const uint8_t> buffer, Proof& proof)
{
    if (buffer.size() != SIZE_IN_BYTES) {
        return false;
    }

    size_t offset = 0;
    for (Commitment* point :
         { &proof.wires, &proof.memory, &proof.grand_product, &proof.quotient, &proof.w, &proof.w_prime }) {
        if (!fflonk_plonk::read_point(buffer, offset, *point)) {
            return false;
        }
    }

    for (FF& evaluation : proof.evaluations) {
        if (!fflonk_plonk::read_scalar(buffer, offset, evaluation)) {
            return false;
        }
    }

    return true;
}

} // namespace bb::ultra_fflonk
