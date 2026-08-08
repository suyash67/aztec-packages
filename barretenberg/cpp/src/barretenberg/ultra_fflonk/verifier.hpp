#pragma once

#include "barretenberg/ultra_fflonk/keys.hpp"
#include "barretenberg/ultra_fflonk/proof.hpp"

#include <span>

namespace bb::ultra_fflonk {

/**
 * @brief The outcome of every check the verifier makes, rather than a single bit.
 *
 * @details Three things a proof has to establish, and they are independent: that it is encoded and
 * shaped the way the protocol says (`well_formed`), that the claimed evaluations really belong to
 * the committed polynomials (`opening`), and that polynomials with those evaluations describe a
 * satisfied circuit (`quotient`). A single boolean cannot tell those apart, and "the proof was
 * rejected" is not evidence that the right check rejected it - every proof element except `W'` is
 * absorbed into the transcript, so tampering with any of them breaks the pairing whether or not the
 * relations are checked at all. So the checks are reported separately and none of them
 * short-circuits.
 *
 * Unlike the three-wire system there is one relation check rather than three, because Sumcheck's
 * subrelations are batched by powers of `alpha` into a single quotient identity - individual
 * relations do not vanish at `xi`, only the batch does. What replaces the per-identity test is
 * per-evaluation sensitivity: every claimed evaluation the relations actually read must change
 * `quotient` when perturbed, and the three that no relation reads must not.
 */
struct VerificationReport {
    bool well_formed = false;
    bool quotient = false;
    bool opening = false;

    [[nodiscard]] bool accepted() const { return well_formed && quotient && opening; }
};

/** @brief Run every check and report each outcome. Never throws; trusts no input. */
[[nodiscard]] VerificationReport verify_detailed(const VerificationKey& key,
                                                 const Proof& proof,
                                                 std::span<const FF> public_inputs);

/**
 * @brief Verify a proof against a verification key and its public inputs.
 *
 * @details Group work is eleven scalar multiplications and one pairing, independent of the circuit;
 * field work is one evaluation of the flavor's relations.
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

} // namespace bb::ultra_fflonk
