#include "barretenberg/ultra_fflonk/layout.hpp"

namespace bb::ultra_fflonk {

Flavor::AllValues to_all_values(std::span<const GroupEvaluations> evaluations)
{
    BB_ASSERT_EQ(evaluations.size(), NUM_GROUPS, "one evaluation set per group");

    Flavor::AllValues row;

    // The preprocessed groups tile `get_precomputed()` in order, so the flavor's own layout decides
    // which selector lands in which group and nothing here has to name them.
    auto precomputed = row.get_precomputed();
    for (size_t g = 0; g < NUM_PREPROCESSED_GROUPS; ++g) {
        for (size_t i = 0; i < PACK_PREPROCESSED; ++i) {
            precomputed[(g * PACK_PREPROCESSED) + i] = evaluations[g].at_xi[i];
        }
    }

    const GroupEvaluations& wires = evaluations[GROUP_WIRES];
    const GroupEvaluations& memory = evaluations[GROUP_MEMORY];
    const GroupEvaluations& grand_product = evaluations[GROUP_GRAND_PRODUCT];

    row.w_l() = wires.at_xi[WIRES_W_L];
    row.w_r() = wires.at_xi[WIRES_W_R];
    row.w_o() = wires.at_xi[WIRES_W_O];
    row.lookup_read_counts() = memory.at_xi[MEMORY_READ_COUNTS];
    row.lookup_read_tags() = memory.at_xi[MEMORY_READ_TAGS];
    row.w_4() = memory.at_xi[MEMORY_W_4];
    row.lookup_inverses() = grand_product.at_xi[GRAND_PRODUCT_LOOKUP_INVERSES];
    row.z_perm() = grand_product.at_xi[GRAND_PRODUCT_Z_PERM];

    row.w_l_shift() = wires.at_xi_omega[WIRES_W_L];
    row.w_r_shift() = wires.at_xi_omega[WIRES_W_R];
    row.w_o_shift() = wires.at_xi_omega[WIRES_W_O];
    row.w_4_shift() = memory.at_xi_omega[MEMORY_W_4];
    row.z_perm_shift() = grand_product.at_xi_omega[GRAND_PRODUCT_Z_PERM];

    return row;
}

std::array<FF, NUM_EVALUATIONS> flatten_evaluations(std::span<const GroupEvaluations> evaluations)
{
    BB_ASSERT_EQ(evaluations.size(), NUM_GROUPS, "one evaluation set per group");

    std::array<FF, NUM_EVALUATIONS> flat{};
    size_t offset = 0;
    for (size_t g = 0; g < NUM_GROUPS; ++g) {
        BB_ASSERT_EQ(evaluations[g].at_xi.size(), GROUP_SHAPES[g].pack, "wrong number of evaluations at xi");
        for (const FF& evaluation : evaluations[g].at_xi) {
            flat[offset++] = evaluation;
        }
        if (GROUP_SHAPES[g].two_point) {
            BB_ASSERT_EQ(
                evaluations[g].at_xi_omega.size(), GROUP_SHAPES[g].pack, "wrong number of evaluations at xi*omega");
            for (const FF& evaluation : evaluations[g].at_xi_omega) {
                flat[offset++] = evaluation;
            }
        }
    }
    BB_ASSERT_EQ(offset, NUM_EVALUATIONS, "the flattening must consume every slot");
    return flat;
}

std::vector<GroupEvaluations> unflatten_evaluations(std::span<const FF> flat)
{
    BB_ASSERT_EQ(flat.size(), NUM_EVALUATIONS, "wrong number of evaluations");

    std::vector<GroupEvaluations> evaluations(NUM_GROUPS);
    size_t offset = 0;
    for (size_t g = 0; g < NUM_GROUPS; ++g) {
        const size_t pack = GROUP_SHAPES[g].pack;
        evaluations[g].at_xi.assign(flat.begin() + static_cast<std::ptrdiff_t>(offset),
                                    flat.begin() + static_cast<std::ptrdiff_t>(offset + pack));
        offset += pack;
        if (GROUP_SHAPES[g].two_point) {
            evaluations[g].at_xi_omega.assign(flat.begin() + static_cast<std::ptrdiff_t>(offset),
                                              flat.begin() + static_cast<std::ptrdiff_t>(offset + pack));
            offset += pack;
        }
    }
    return evaluations;
}

RunningSumEvaluations to_running_sum_evaluations(std::span<const GroupEvaluations> evaluations)
{
    const GroupEvaluations& grand_product = evaluations[GROUP_GRAND_PRODUCT];

    RunningSumEvaluations result;
    for (size_t j = 0; j < NUM_RUNNING_SUMS; ++j) {
        result.at_xi[j] = grand_product.at_xi[GRAND_PRODUCT_RUNNING_SUM + j];
        result.at_xi_omega[j] = grand_product.at_xi_omega[GRAND_PRODUCT_RUNNING_SUM + j];
    }
    return result;
}

} // namespace bb::ultra_fflonk
