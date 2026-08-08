#pragma once

#include "barretenberg/fflonk/keys.hpp"
#include "barretenberg/fflonk/proof.hpp"
#include "barretenberg/fflonk/transcript.hpp"

#include <array>

namespace bb::fflonk_plonk {

/**
 * @brief Prove the circuit `key` was preprocessed from.
 *
 * @details Aborts rather than returning an unsatisfiable proof: every quotient division is checked
 * for an exact remainder, the grand product is checked to close, and the three constraint identities
 * are re-checked against the prover's own claimed evaluations before the proof is returned. A
 * mistake in the circuit or in this file therefore surfaces here, not as a pairing failure with no
 * explanation on chain.
 */
Proof prove(const ProvingKey& key);

namespace prover_detail {

/**
 * @brief Rounds 3 to 5: open every group at the challenge and batch the openings into `W`, `W'`.
 *
 * @details Takes the four packed group polynomials, with the transcript positioned immediately
 * before `xi` is squeezed - so the caller has already absorbed the group commitments and written
 * them into `proof`. Fills in `proof.evaluations`, `proof.w` and `proof.w_prime`, and returns `xi`.
 *
 * Split out from `prove` so the soundness tests can drive it with deliberately wrong polynomials.
 * The opening argument only certifies that the claimed evaluations belong to the committed
 * polynomials; it says nothing about whether those polynomials satisfy the circuit. Separating the
 * two makes it testable that the verifier's constraint checks are what reject a forgery.
 */
FF open_and_batch(const ProvingKey& key,
                  const std::array<const std::vector<FF>*, NUM_GROUPS>& packed_groups,
                  Transcript& transcript,
                  Proof& proof);

} // namespace prover_detail

} // namespace bb::fflonk_plonk
