#pragma once

#include "barretenberg/fflonk/keys.hpp"
#include "barretenberg/fflonk/polynomial_utils.hpp"

#include <array>
#include <cstdint>
#include <span>
#include <vector>

namespace bb::fflonk_plonk {

/**
 * @brief Indices into `Proof::evaluations`, in the order they are sent and absorbed.
 * @details The order is protocol: the Solidity verifier reads the same offsets out of calldata.
 */
enum EvaluationIndex : size_t {
    EVAL_Q_L = 0,
    EVAL_Q_R = 1,
    EVAL_Q_O = 2,
    EVAL_Q_M = 3,
    EVAL_Q_C = 4,
    EVAL_S_1 = 5,
    EVAL_S_2 = 6,
    EVAL_S_3 = 7,
    EVAL_A = 8,
    EVAL_B = 9,
    EVAL_C = 10,
    EVAL_T0_LO = 11,
    EVAL_T0_HI = 12,
    EVAL_Z = 13,
    EVAL_Z_OMEGA = 14,
    EVAL_T1 = 15,
    EVAL_T2_LO = 16,
    EVAL_T2_MID = 17,
    EVAL_T2_HI = 18,
};

/**
 * @brief Five group elements and nineteen scalars: 928 bytes, whatever the circuit contains.
 *
 * @details `c1` packs the wires and the gate quotient, `c2` the grand product, `c3` the permutation
 * quotients; `w` and `w_prime` are the two halves of the batched opening. The preprocessed group's
 * commitment lives in the verification key rather than the proof.
 */
struct Proof {
    Commitment c1 = Commitment::infinity();
    Commitment c2 = Commitment::infinity();
    Commitment c3 = Commitment::infinity();
    Commitment w = Commitment::infinity();
    Commitment w_prime = Commitment::infinity();
    std::array<FF, NUM_EVALUATIONS> evaluations{};

    static constexpr size_t NUM_COMMITMENTS = 5;
    static constexpr size_t SIZE_IN_BYTES = NUM_COMMITMENTS * 64 + NUM_EVALUATIONS * 32;

    /** @brief Big-endian words, in transcript order: the five points then the nineteen scalars. */
    [[nodiscard]] std::vector<uint8_t> to_buffer() const;

    /**
     * @brief Parse a proof, rejecting anything a verifier must not accept.
     *
     * @details Enforces the exact length, that every scalar is a canonical residue below the field
     * modulus, and that every point is a canonical coordinate pair on the curve. BN254's G1 has
     * prime order, so on-curve implies in-subgroup and no separate cofactor check is needed. The
     * encoding of the point at infinity is `(0, 0)`, matching the EVM.
     */
    [[nodiscard]] static bool from_buffer(std::span<const uint8_t> buffer, Proof& proof);
};

} // namespace bb::fflonk_plonk
