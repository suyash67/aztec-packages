#pragma once

#include "barretenberg/commitment_schemes/commitment_key.hpp"
#include "barretenberg/polynomials/evaluation_domain.hpp"
#include "barretenberg/ultra_fflonk/layout.hpp"
#include "barretenberg/ultra_honk/prover_instance.hpp"

#include <memory>
#include <vector>

namespace bb::ultra_fflonk {

/**
 * @brief The smallest circuit this proof system accepts.
 *
 * @details Two bounds meet at seven rows: the quotient domain `8n` must exceed the numerator's
 * `6n+6`, and the six quotient chunks must hold its `5n+7` coefficients. Eight is the next power of
 * two, and it also leaves room for the four disabled rows at the top of an Ultra trace.
 */
static constexpr size_t MIN_CIRCUIT_SIZE = 8;

/**
 * @brief Everything a verifier needs besides the proof and the public inputs.
 *
 * @details The 28 precomputed columns are committed as four interleaved groups, so the key holds
 * four group elements however many selectors, sigmas, ids and tables the arithmetization has.
 */
struct VerificationKey {
    size_t circuit_size = 0;
    size_t num_public_inputs = 0;
    /** @brief Row of the trace the public inputs start at; feeds the permutation's `public_input_delta`. */
    size_t pub_inputs_offset = 0;
    FF omega = FF::zero();
    std::array<Commitment, NUM_PREPROCESSED_GROUPS> preprocessed{};

    /** @brief The single word that binds a proof to this circuit, absorbed first in the transcript. */
    [[nodiscard]] FF hash() const;

    /** @brief Four metadata words then the four group commitments, in the proof's wire format. */
    static constexpr size_t NUM_METADATA_WORDS = 4;
    static constexpr size_t SIZE_IN_BYTES = (NUM_METADATA_WORDS * 32) + (NUM_PREPROCESSED_GROUPS * 64);
    [[nodiscard]] std::vector<uint8_t> to_buffer() const;
    [[nodiscard]] static bool from_buffer(std::span<const uint8_t> buffer, VerificationKey& key);

    bool operator==(const VerificationKey& other) const = default;
};

/**
 * @brief The prover's view of a preprocessed circuit.
 *
 * @details Holds the Honk prover instance - which is what makes every ACIR constraint, every custom
 * gate and the whole gate-count profile come for free - alongside the packed preprocessed groups and
 * the two evaluation domains this proof system needs on top of it.
 *
 * The instance is mutated by proving: `w_4` picks up memory records once `eta` is known, and the
 * lookup inverses and grand product are built once `beta` and `gamma` are. Every one of those is a
 * full overwrite, so the same key can be proven from repeatedly.
 */
struct ProvingKey {
    std::shared_ptr<ProverInstance_<Flavor>> instance;

    /** @brief `g_0 .. g_3`, in coefficient form; kept so the prover need not re-pack per proof. */
    std::array<std::vector<FF>, NUM_PREPROCESSED_GROUPS> packed_preprocessed;

    std::shared_ptr<EvaluationDomain<FF>> small_domain; // size n
    std::shared_ptr<EvaluationDomain<FF>> large_domain; // size 8n
    std::shared_ptr<CommitmentKey<Curve>> commitment_key;

    VerificationKey verification_key;

    [[nodiscard]] size_t circuit_size() const { return verification_key.circuit_size; }
    [[nodiscard]] const std::vector<FF>& public_inputs() const { return instance->public_inputs; }
};

/**
 * @brief Preprocess a circuit: build the Honk prover instance, then commit to the precomputed
 * columns as four interleaved groups.
 *
 * @details Consumes `UltraCircuitBuilder` unchanged - the same trace UltraHonk proves - so the ACIR
 * frontend, the lookup tables and the custom gates are shared rather than reimplemented.
 */
ProvingKey preprocess(Flavor::CircuitBuilder& circuit);

/** @brief Column `index` of a packed group, recovered by striding the interleaving. */
std::vector<FF> unpack_column(const std::vector<FF>& packed, size_t pack, size_t index);

} // namespace bb::ultra_fflonk
