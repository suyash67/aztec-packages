#pragma once

#include "barretenberg/fflonk/transcript.hpp"
#include "barretenberg/ultra_fflonk/keys.hpp"
#include "barretenberg/ultra_fflonk/proof.hpp"

namespace bb::ultra_fflonk {

using Transcript = fflonk_plonk::Transcript;

/**
 * @brief Prove the circuit `key` was preprocessed from.
 *
 * @details Aborts rather than returning an unsatisfiable proof: the quotient division is checked for
 * an exact remainder, each running sum is checked to close, and the batched identity is re-checked
 * against the prover's own claimed evaluations before the proof is returned. A mistake in the
 * circuit or in this file therefore surfaces here, not as a pairing failure with no explanation on
 * chain.
 *
 * Mutates `key.instance`: `w_4` picks up its memory records once `eta` is known, and the lookup
 * inverses and grand product are rebuilt once `beta` and `gamma` are. Each is a full overwrite, so
 * proving twice from one key is safe and produces two independently blinded proofs.
 */
Proof prove(ProvingKey& key);

namespace prover_detail {

/**
 * @brief The closing rounds: open every group at the challenge and batch the openings into `W`, `W'`.
 *
 * @details Takes the eight packed group polynomials with the transcript positioned immediately
 * before `xi` is squeezed - so the caller has already absorbed the group commitments and written
 * them into `proof`. Fills in `proof.evaluations`, `proof.w` and `proof.w_prime`, and returns `xi`.
 *
 * Split out from `prove` so the soundness tests can drive it with deliberately wrong polynomials.
 * The opening argument only certifies that the claimed evaluations belong to the committed
 * polynomials; it says nothing about whether those polynomials satisfy the circuit. Separating the
 * two makes it testable that the verifier's relation check is what rejects a forgery.
 */
FF open_and_batch(const ProvingKey& key,
                  std::span<const std::vector<FF>* const> packed_groups,
                  Transcript& transcript,
                  Proof& proof);

} // namespace prover_detail

} // namespace bb::ultra_fflonk
