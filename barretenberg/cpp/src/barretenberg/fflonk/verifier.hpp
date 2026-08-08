#pragma once

#include "barretenberg/fflonk/keys.hpp"
#include "barretenberg/fflonk/proof.hpp"

#include <span>

namespace bb::fflonk_plonk {

/**
 * @brief The outcome of every check the verifier makes, rather than a single bit.
 *
 * @details Two things a proof has to establish, and they are independent: that the claimed
 * evaluations really belong to the committed polynomials (`opening`), and that polynomials with
 * those evaluations describe a satisfied circuit (the three identities). A single boolean cannot
 * tell those apart, and "the proof was rejected" is not evidence that the right check rejected it -
 * so the checks are reported separately and none of them short-circuits.
 */
struct VerificationReport {
    bool well_formed = false;
    bool gate_identity = false;
    bool grand_product_start = false;
    bool permutation_identity = false;
    bool opening = false;

    [[nodiscard]] bool accepted() const
    {
        return well_formed && gate_identity && grand_product_start && permutation_identity && opening;
    }
};

/** @brief Run every check and report each outcome. Never throws; trusts no input. */
[[nodiscard]] VerificationReport verify_detailed(const VerificationKey& key,
                                                 const Proof& proof,
                                                 std::span<const FF> public_inputs);

/**
 * @brief Verify a proof against a verification key and its public inputs.
 *
 * @details Group work is seven scalar multiplications and one pairing, independent of the circuit.
 */
[[nodiscard]] bool verify(const VerificationKey& key, const Proof& proof, std::span<const FF> public_inputs);

/**
 * @brief Verify a serialised proof.
 *
 * @details The overload to prefer when the proof came from outside: it enforces the encoding as
 * well - exact length, canonical scalars, on-curve points - which the struct overload has to assume.
 */
[[nodiscard]] bool verify(const VerificationKey& key,
                          std::span<const uint8_t> proof_bytes,
                          std::span<const FF> public_inputs);

} // namespace bb::fflonk_plonk
