#pragma once

#include "barretenberg/common/constexpr_utils.hpp"
#include "barretenberg/fflonk/batched_opening.hpp"
#include "barretenberg/flavor/ultra_flavor.hpp"

#include <array>
#include <cstddef>

namespace bb::ultra_fflonk {

using Flavor = UltraFlavor;
using Curve = fflonk_plonk::Curve;
using FF = fflonk_plonk::FF;
using Commitment = fflonk_plonk::Commitment;
using GroupElement = fflonk_plonk::GroupElement;
using fflonk_plonk::GroupEvaluations;
using fflonk_plonk::GroupShape;

using EntityId = Flavor::AllValues::EntityId;
using Relations = Flavor::Relations;

/** @brief Subrelations that must vanish at every row, batched by powers of `alpha`. */
static constexpr size_t NUM_SUBRELATIONS = Flavor::NUM_SUBRELATIONS;

namespace detail {

/** @brief How many of a relation's subrelations are only required to sum to zero over the trace. */
template <typename Relation> constexpr size_t num_dependent_subrelations()
{
    if constexpr (requires { Relation::SUBRELATION_LINEARLY_INDEPENDENT; }) {
        size_t count = 0;
        for (const bool independent : Relation::SUBRELATION_LINEARLY_INDEPENDENT) {
            if (!independent) {
                ++count;
            }
        }
        return count;
    } else {
        return 0;
    }
}

template <typename Tuple> constexpr size_t count_dependent_subrelations()
{
    return []<size_t... I>(std::index_sequence<I...>) {
        return (0 + ... + num_dependent_subrelations<std::tuple_element_t<I, Tuple>>());
    }(std::make_index_sequence<std::tuple_size_v<Tuple>>());
}

} // namespace detail

/**
 * @brief Subrelations that Sumcheck only requires to sum to zero across the whole trace.
 *
 * @details Sumcheck can enforce `sum_rows f = 0` directly; a quotient argument cannot - dividing by
 * `Z_H` says a polynomial vanishes at every row, which is a strictly stronger and different claim.
 * Each such subrelation therefore gets a committed running sum `S` with
 *
 * ```
 * S(omega X) - S(X) - f(X) = 0    for all X in H
 * ```
 *
 * Summing that identity over `H` telescopes around the cycle and leaves `sum_rows f = 0`, so the
 * mere existence of `S` is the statement Sumcheck was making. `S` is not pinned to start at zero
 * because it does not need to be: soundness comes from the sum, not from the starting value.
 *
 * `LogDerivLookupRelation`'s lookup subrelation and `MemoryRelation`'s ROM-LogUp sum subrelation are
 * the two in the Ultra set; the count is derived from the flavor rather than assumed, so that adding
 * a third breaks the build here instead of silently going unenforced.
 */
static constexpr size_t NUM_RUNNING_SUMS = detail::count_dependent_subrelations<Relations>();
static_assert(NUM_RUNNING_SUMS == 2, "the group layout below reserves exactly two running-sum columns");

/** @brief The 28 precomputed columns, packed seven at a time. */
static constexpr size_t NUM_PRECOMPUTED = Flavor::NUM_PRECOMPUTED_ENTITIES;
static constexpr size_t PACK_PREPROCESSED = 7;
static constexpr size_t NUM_PREPROCESSED_GROUPS = NUM_PRECOMPUTED / PACK_PREPROCESSED;
static_assert(NUM_PRECOMPUTED % PACK_PREPROCESSED == 0,
              "the preprocessed groups must tile the precomputed columns exactly");

/**
 * @brief The quotient is split into this many chunks of degree at most `n`.
 *
 * @details Derived in PROTOCOL.md: with witness columns of degree `n+1` and a maximum relation
 * degree of 6, the numerator reaches `6n+6`, so the quotient reaches `5n+6` and needs six chunks.
 * `prover.cpp` asserts the realised degree against this rather than trusting the derivation.
 */
static constexpr size_t NUM_QUOTIENT_CHUNKS = 6;

/** @brief The quotient FFT runs at this multiple of the circuit size; `6n+6 < 8n` for `n > 3`. */
static constexpr size_t QUOTIENT_DOMAIN_FACTOR = 8;

/**
 * @brief The SRS must cover the widest packed group: seven preprocessed columns of degree `n-1`
 * interleave to degree `7n-1`.
 */
static constexpr size_t SRS_SIZE_FACTOR = PACK_PREPROCESSED;

/** @brief Blinders per witness column, one per point the column is revealed at. */
static constexpr size_t NUM_WITNESS_BLINDERS = 2;

static constexpr size_t PACK_WIRES = 3;                            // w_l, w_r, w_o
static constexpr size_t PACK_MEMORY = 3;                           // lookup_read_counts, lookup_read_tags, w_4
static constexpr size_t PACK_GRAND_PRODUCT = 2 + NUM_RUNNING_SUMS; // lookup_inverses, z_perm, running sums

/** @brief Group indices, in the order the `nu` powers of the batched opening apply. */
static constexpr size_t GROUP_WIRES = NUM_PREPROCESSED_GROUPS;
static constexpr size_t GROUP_MEMORY = GROUP_WIRES + 1;
static constexpr size_t GROUP_GRAND_PRODUCT = GROUP_MEMORY + 1;
static constexpr size_t GROUP_QUOTIENT = GROUP_GRAND_PRODUCT + 1;
static constexpr size_t NUM_GROUPS = GROUP_QUOTIENT + 1;

/**
 * @brief Every group, in `nu`-power order.
 *
 * @details The three witness groups are opened at both `xi` and `xi*omega` because each contains at
 * least one column the relations read shifted; the columns that ride along (`lookup_read_counts`,
 * `lookup_read_tags`, `lookup_inverses`) are revealed at `xi*omega` as a consequence, which is what
 * their second blinder pays for. The preprocessed groups are public, and the quotient chunks are
 * never read shifted, so both are opened at `xi` alone.
 */
static constexpr std::array<GroupShape, NUM_GROUPS> GROUP_SHAPES = {
    GroupShape{ .pack = PACK_PREPROCESSED, .two_point = false },
    GroupShape{ .pack = PACK_PREPROCESSED, .two_point = false },
    GroupShape{ .pack = PACK_PREPROCESSED, .two_point = false },
    GroupShape{ .pack = PACK_PREPROCESSED, .two_point = false },
    GroupShape{ .pack = PACK_WIRES, .two_point = true },
    GroupShape{ .pack = PACK_MEMORY, .two_point = true },
    GroupShape{ .pack = PACK_GRAND_PRODUCT, .two_point = true },
    GroupShape{ .pack = NUM_QUOTIENT_CHUNKS, .two_point = false },
};

/** @brief Column positions inside the witness groups. The order is protocol. */
static constexpr size_t WIRES_W_L = 0;
static constexpr size_t WIRES_W_R = 1;
static constexpr size_t WIRES_W_O = 2;
static constexpr size_t MEMORY_READ_COUNTS = 0;
static constexpr size_t MEMORY_READ_TAGS = 1;
static constexpr size_t MEMORY_W_4 = 2;
static constexpr size_t GRAND_PRODUCT_LOOKUP_INVERSES = 0;
static constexpr size_t GRAND_PRODUCT_Z_PERM = 1;
static constexpr size_t GRAND_PRODUCT_RUNNING_SUM = 2; // the two running sums follow, in relation order

/** @brief How many field elements a proof carries: every group's opening set, flattened. */
consteval size_t num_evaluations()
{
    size_t total = 0;
    for (const GroupShape& shape : GROUP_SHAPES) {
        total += shape.num_evaluations();
    }
    return total;
}
static constexpr size_t NUM_EVALUATIONS = num_evaluations();

/** @brief Commitments in the proof: the four witness/quotient groups plus `W` and `W'`. */
static constexpr size_t NUM_PROOF_COMMITMENTS = 4 + 2;

/**
 * @brief Rebuild the flavor's per-row value container from a proof's claimed evaluations.
 *
 * @details The verifier evaluates the relations exactly as the prover does - through the same
 * accumulators in the `relations` module - so it needs the entities in the flavor's own container. The
 * shifted entities come from the `xi*omega` half of the two-point groups, which is the whole reason
 * those groups are opened twice.
 */
Flavor::AllValues to_all_values(std::span<const GroupEvaluations> evaluations);

/** @brief The claimed evaluations of the two running sums at `xi` and at `xi*omega`. */
struct RunningSumEvaluations {
    std::array<FF, NUM_RUNNING_SUMS> at_xi{};
    std::array<FF, NUM_RUNNING_SUMS> at_xi_omega{};
};

RunningSumEvaluations to_running_sum_evaluations(std::span<const GroupEvaluations> evaluations);

/**
 * @brief The proof's flat evaluation array: each group's `xi` openings, then its `xi*omega` ones.
 * @details The order is protocol - a Solidity verifier reads the same offsets out of calldata.
 */
std::array<FF, NUM_EVALUATIONS> flatten_evaluations(std::span<const GroupEvaluations> evaluations);

/** @brief The inverse of `flatten_evaluations`; the group shapes say where each boundary falls. */
std::vector<GroupEvaluations> unflatten_evaluations(std::span<const FF> flat);

} // namespace bb::ultra_fflonk
