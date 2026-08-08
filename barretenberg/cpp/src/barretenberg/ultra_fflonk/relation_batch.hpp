#pragma once

#include "barretenberg/relations/relation_parameters.hpp"
#include "barretenberg/ultra_fflonk/layout.hpp"

#include <array>
#include <span>

namespace bb::ultra_fflonk {

/** @brief What the relations say at one row, split by how each subrelation is enforced. */
struct RelationValues {
    /** @brief `sum_i alpha^i R_i` over the subrelations that must vanish at this row. */
    FF per_row = FF::zero();
    /** @brief The subrelations that only have to sum to zero across the trace, unbatched. */
    std::array<FF, NUM_RUNNING_SUMS> summed{};
    /** @brief `alpha^NUM_SUBRELATIONS`, where the running-sum identities pick up. */
    FF next_alpha_power = FF::one();
};

/**
 * @brief Evaluate every relation of the flavor at one row and batch the per-row subrelations.
 *
 * @details This is the whole arithmetization: the relations are barretenberg's own, called through
 * the same `accumulate` entry point Sumcheck uses, so a quotient built from this and a Sumcheck
 * proof are constraining the same circuit by construction. Prover and verifier both come through
 * here - the prover at every point of the quotient domain, the verifier once at `xi` - which is what
 * makes "the verifier checks the same relations the prover divided out" a property of the code
 * rather than of a review.
 *
 * The `alpha` power advances across *every* subrelation, including the ones that are held back for a
 * running sum, so a subrelation's power depends only on the flavor's relation order and not on which
 * subrelations happen to be linearly dependent.
 */
inline RelationValues accumulate_relations(const Flavor::AllValues& row,
                                           const RelationParameters<FF>& parameters,
                                           const FF& alpha)
{
    RelationValues result;
    FF alpha_power = FF::one();
    size_t summed_index = 0;

    bb::constexpr_for<0, std::tuple_size_v<Relations>, 1>([&]<size_t i>() {
        using Relation = std::tuple_element_t<i, Relations>;

        typename Relation::SumcheckArrayOfValuesOverSubrelations values;
        for (FF& value : values) {
            value = FF::zero();
        }
        Relation::accumulate(values, row, parameters, FF::one());

        for (size_t subrelation = 0; subrelation < values.size(); ++subrelation) {
            bool per_row = true;
            if constexpr (requires { Relation::SUBRELATION_LINEARLY_INDEPENDENT; }) {
                per_row = Relation::SUBRELATION_LINEARLY_INDEPENDENT[subrelation];
            }
            if (per_row) {
                result.per_row += alpha_power * values[subrelation];
            } else {
                result.summed[summed_index++] = values[subrelation];
            }
            alpha_power *= alpha;
        }
    });

    result.next_alpha_power = alpha_power;
    return result;
}

/**
 * @brief Every subrelation's raw value at one row, in the flavor's relation order, unbatched.
 *
 * @details What the Solidity verifier is checked against. `alpha` batches these into one number, and
 * a single wrong term inside that number is invisible; exporting them separately is what lets a
 * differential test say *which* relation a port got wrong.
 */
inline std::array<FF, NUM_SUBRELATIONS> all_subrelations(const Flavor::AllValues& row,
                                                         const RelationParameters<FF>& parameters)
{
    std::array<FF, NUM_SUBRELATIONS> values{};
    size_t index = 0;

    bb::constexpr_for<0, std::tuple_size_v<Relations>, 1>([&]<size_t i>() {
        using Relation = std::tuple_element_t<i, Relations>;

        typename Relation::SumcheckArrayOfValuesOverSubrelations accumulated;
        for (FF& value : accumulated) {
            value = FF::zero();
        }
        Relation::accumulate(accumulated, row, parameters, FF::one());

        for (const FF& value : accumulated) {
            values[index++] = value;
        }
    });

    BB_ASSERT_EQ(index, NUM_SUBRELATIONS, "the relation tuple must cover every subrelation");
    return values;
}

/**
 * @brief Just the subrelations that have to sum to zero across the trace, at one row.
 *
 * @details The prover needs these before it can build the running sums, and at that point the grand
 * product does not exist yet - so this skips every relation that has no such subrelation rather than
 * evaluating the whole set against a `z_perm` that is still zero.
 */
inline std::array<FF, NUM_RUNNING_SUMS> summed_subrelations(const Flavor::AllValues& row,
                                                            const RelationParameters<FF>& parameters)
{
    std::array<FF, NUM_RUNNING_SUMS> summed{};
    size_t summed_index = 0;

    bb::constexpr_for<0, std::tuple_size_v<Relations>, 1>([&]<size_t i>() {
        using Relation = std::tuple_element_t<i, Relations>;
        if constexpr (detail::num_dependent_subrelations<Relation>() > 0) {
            typename Relation::SumcheckArrayOfValuesOverSubrelations values;
            for (FF& value : values) {
                value = FF::zero();
            }
            Relation::accumulate(values, row, parameters, FF::one());

            for (size_t subrelation = 0; subrelation < values.size(); ++subrelation) {
                if (!Relation::SUBRELATION_LINEARLY_INDEPENDENT[subrelation]) {
                    summed[summed_index++] = values[subrelation];
                }
            }
        }
    });

    return summed;
}

/**
 * @brief The quotient numerator at one row: the per-row subrelations plus the running-sum identities.
 *
 * @details `sum_i alpha^i R_i + sum_j alpha^{NUM_SUBRELATIONS + j} (S_j(omega X) - S_j(X) - f_j(X))`.
 * Dividing this by `Z_H` is exactly the claim Sumcheck makes: every `R_i` vanishes at every row, and
 * every `f_j` sums to zero across the trace.
 */
inline FF quotient_numerator(const Flavor::AllValues& row,
                             const RelationParameters<FF>& parameters,
                             const FF& alpha,
                             std::span<const FF, NUM_RUNNING_SUMS> running_sum,
                             std::span<const FF, NUM_RUNNING_SUMS> running_sum_shifted)
{
    const RelationValues values = accumulate_relations(row, parameters, alpha);

    FF numerator = values.per_row;
    FF alpha_power = values.next_alpha_power;
    for (size_t j = 0; j < NUM_RUNNING_SUMS; ++j) {
        numerator += alpha_power * (running_sum_shifted[j] - running_sum[j] - values.summed[j]);
        alpha_power *= alpha;
    }
    return numerator;
}

} // namespace bb::ultra_fflonk
