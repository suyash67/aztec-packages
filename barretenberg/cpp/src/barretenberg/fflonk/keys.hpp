#pragma once

#include "barretenberg/commitment_schemes/commitment_key.hpp"
#include "barretenberg/fflonk/batched_opening.hpp"
#include "barretenberg/fflonk/circuit_builder.hpp"
#include "barretenberg/fflonk/polynomial_utils.hpp"
#include "barretenberg/polynomials/evaluation_domain.hpp"

#include <array>
#include <memory>

namespace bb::fflonk_plonk {

/** @brief Packing factors, one per group. The column order inside each group is protocol. */
static constexpr size_t PACK_PREPROCESSED = 8; // q_l q_r q_o q_m q_c s_1 s_2 s_3
static constexpr size_t PACK_WIRES = 5;        // a b c t0_lo t0_hi
static constexpr size_t PACK_GRAND_PRODUCT = 1;
static constexpr size_t PACK_QUOTIENTS = 4; // t1 t2_lo t2_mid t2_hi

/** @brief Where S_sigma1 starts inside the preprocessed group's column order. */
static constexpr size_t PREPROCESSED_SIGMA_OFFSET = 5;

static constexpr size_t NUM_GROUPS = 4;

/**
 * @brief The four groups, in the order the `nu` powers of the batched opening apply to them.
 * @details Only the grand product is opened at two points; it is alone in its group for exactly that
 * reason - see PROTOCOL.md section 7.
 */
static constexpr std::array<GroupShape, NUM_GROUPS> GROUP_SHAPES = {
    GroupShape{ .pack = PACK_PREPROCESSED, .two_point = false },
    GroupShape{ .pack = PACK_WIRES, .two_point = false },
    GroupShape{ .pack = PACK_GRAND_PRODUCT, .two_point = true },
    GroupShape{ .pack = PACK_QUOTIENTS, .two_point = false },
};

/** @brief The quotient FFT runs at this multiple of the circuit size; see PROTOCOL.md section 9. */
static constexpr size_t QUOTIENT_DOMAIN_FACTOR = 8;

/** @brief The SRS must hold this many points: the widest packed group is the preprocessed one. */
static constexpr size_t SRS_SIZE_FACTOR = PACK_PREPROCESSED;

/** @brief Number of scalars the prover sends: 8 + 5 + 2 + 4. */
static constexpr size_t NUM_EVALUATIONS = 19;

/**
 * @brief Everything a verifier needs besides the proof and the public inputs.
 *
 * @details `c0` commits to the whole preprocessed group at once, so the verification key holds one
 * group element however many selectors the arithmetization has.
 */
struct VerificationKey {
    size_t circuit_size = 0;
    size_t num_public_inputs = 0;
    FF omega = FF::zero();
    FF k1 = FF::zero();
    FF k2 = FF::zero();
    Commitment c0 = Commitment::infinity();

    /** @brief The single word that binds a proof to this circuit, absorbed first in the transcript. */
    [[nodiscard]] FF hash() const;

    /** @brief Five metadata words then `C0`, in the proof's wire format. */
    static constexpr size_t NUM_METADATA_WORDS = 5;
    static constexpr size_t SIZE_IN_BYTES = (NUM_METADATA_WORDS * 32) + 64;
    [[nodiscard]] std::vector<uint8_t> to_buffer() const;
    [[nodiscard]] static bool from_buffer(std::span<const uint8_t> buffer, VerificationKey& key);

    bool operator==(const VerificationKey& other) const = default;
};

/** @brief The preprocessed polynomials, in the group's column order, in coefficient form. */
struct PreprocessedPolynomials {
    std::array<std::vector<FF>, PACK_PREPROCESSED> columns;

    [[nodiscard]] const std::vector<FF>& q_l() const { return columns[0]; }
    [[nodiscard]] const std::vector<FF>& q_r() const { return columns[1]; }
    [[nodiscard]] const std::vector<FF>& q_o() const { return columns[2]; }
    [[nodiscard]] const std::vector<FF>& q_m() const { return columns[3]; }
    [[nodiscard]] const std::vector<FF>& q_c() const { return columns[4]; }
    [[nodiscard]] const std::vector<FF>& s_1() const { return columns[5]; }
    [[nodiscard]] const std::vector<FF>& s_2() const { return columns[6]; }
    [[nodiscard]] const std::vector<FF>& s_3() const { return columns[7]; }
};

/**
 * @brief The prover's view of a preprocessed circuit.
 *
 * @details Holds the trace, the preprocessed polynomials in both bases the prover needs (coefficient
 * form to pack and open, Lagrange form for the grand product), the two evaluation domains, and the
 * commitment key. The verification key is carried alongside so a prover cannot accidentally prove
 * against a key that was derived from a different circuit.
 */
struct ProvingKey {
    Trace trace;
    PreprocessedPolynomials preprocessed;
    std::array<std::vector<FF>, NUM_WIRES> sigma_lagrange;
    std::vector<FF> public_input_lagrange; // -x_i on the public-input rows, zero elsewhere
    std::vector<FF> public_input_poly;     // the same, in coefficient form

    std::vector<FF> packed_preprocessed; // g_0, kept so the prover need not re-pack per proof

    std::shared_ptr<EvaluationDomain<FF>> small_domain; // size n
    std::shared_ptr<EvaluationDomain<FF>> large_domain; // size 8n
    std::shared_ptr<CommitmentKey<Curve>> commitment_key;

    VerificationKey verification_key;
};

/**
 * @brief Preprocess a circuit: build the trace, derive the copy-constraint permutation, and commit
 * to the preprocessed group.
 *
 * @param minimum_size lower bound on the circuit size, rounded up to a power of two.
 */
ProvingKey preprocess(const CircuitBuilder& builder, size_t minimum_size = 0);

/** @brief Preprocess a trace that was already built (or loaded). */
ProvingKey preprocess(Trace trace);

} // namespace bb::fflonk_plonk
