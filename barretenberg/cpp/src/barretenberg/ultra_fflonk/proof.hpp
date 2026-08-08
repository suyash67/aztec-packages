#pragma once

#include "barretenberg/ultra_fflonk/layout.hpp"

#include <array>
#include <cstdint>
#include <span>
#include <vector>

namespace bb::ultra_fflonk {

/**
 * @brief Six group elements and fifty-four scalars, whatever the circuit contains.
 *
 * @details `wires` packs `w_l, w_r, w_o`; `memory` packs the lookup read counts and tags with the
 * finalised `w_4`; `grand_product` packs the lookup inverses, `z_perm` and the two running sums;
 * `quotient` packs the six quotient chunks. `w` and `w_prime` are the two halves of the batched
 * opening. The four preprocessed groups' commitments live in the verification key.
 */
struct Proof {
    Commitment wires = Commitment::infinity();
    Commitment memory = Commitment::infinity();
    Commitment grand_product = Commitment::infinity();
    Commitment quotient = Commitment::infinity();
    Commitment w = Commitment::infinity();
    Commitment w_prime = Commitment::infinity();
    std::array<FF, NUM_EVALUATIONS> evaluations{};

    static constexpr size_t SIZE_IN_BYTES = (NUM_PROOF_COMMITMENTS * 64) + (NUM_EVALUATIONS * 32);

    /** @brief The four committed groups, in the order the `nu` powers apply after the preprocessed ones. */
    [[nodiscard]] std::array<Commitment, 4> committed_groups() const
    {
        return { wires, memory, grand_product, quotient };
    }

    /** @brief Big-endian words, in transcript order: the six points then the fifty-four scalars. */
    [[nodiscard]] std::vector<uint8_t> to_buffer() const;

    /**
     * @brief Parse a proof, rejecting anything a verifier must not accept.
     * @details Enforces the exact length, that every scalar is a canonical residue below the field
     * modulus, and that every point is a canonical coordinate pair on the curve.
     */
    [[nodiscard]] static bool from_buffer(std::span<const uint8_t> buffer, Proof& proof);
};

} // namespace bb::ultra_fflonk
