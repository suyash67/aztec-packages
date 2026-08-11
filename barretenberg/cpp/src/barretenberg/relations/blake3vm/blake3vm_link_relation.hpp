#pragma once
/**
 * @brief Gathers the Blake3VM's claim data into one densely-indexed committed column.
 *
 * A consuming circuit binds to the VM by carrying the same claim vector in a committed column of
 * its own proof; since both proofs are KZG over BN254 against the same SRS, two columns holding the
 * same values at the same indices are the same group element, and comparing them is the whole link
 * (README.md, "The linking argument"). That only works if the VM's copy is dense — entry i at index
 * i — while the data itself is scattered through the trace: a block's eight message chunks sit on
 * its round-0 G rows and a digest's four chunks on the output rows of the call's final compression.
 *
 * `claim_value` is what a row contributes and is pinned row-locally by
 * `Blake3VMClaimRelation`; `link_value` is the dense copy, laid out so entry i sits at row
 * NUM_DISABLED_ROWS_IN_SUMCHECK + i — the same offset bb's databus polynomials use, so the section
 * aligns index for index with a consuming Mega circuit's calldata column. A log-derivative
 * permutation over the pairs (index, value) forces the two multisets to agree:
 *
 *   Σ_{claim rows} 1/(claim_index + OFFSET + β·claim_value + γ)
 *     − Σ_{link rows} 1/(row_index + β·link_value + γ) = 0
 *
 * Indices appear in both tuples, so a permutation of the *values* alone cannot satisfy it: entry i
 * on one side can only be matched by index OFFSET + i on the other.
 */

#include "barretenberg/constants.hpp"
#include "barretenberg/honk/proof_system/logderivative_library.hpp"
#include "barretenberg/relations/relation_types.hpp"

namespace bb {

template <typename FF_> class Blake3VMLinkRelationImpl {
  public:
    using FF = FF_;

    static constexpr size_t NUM_LOOKUP_TERMS = 1;
    static constexpr size_t NUM_TABLE_TERMS = 1;
    static constexpr size_t LENGTH = NUM_LOOKUP_TERMS + NUM_TABLE_TERMS + 3;

    static constexpr std::array<size_t, 2> SUBRELATION_PARTIAL_LENGTHS{ LENGTH, LENGTH };
    static constexpr std::array<bool, 2> SUBRELATION_LINEARLY_INDEPENDENT = { true, false };

    template <typename AllValues> static bool operation_exists_at_row(const AllValues& row)
    {
        return row.q_claim != FF(0) || row.q_link != FF(0);
    }

    template <typename AllEntities> static auto& get_inverse_polynomial(AllEntities& in) { return in.inv_link; }

    template <typename Accumulator, typename AllEntities>
    static Accumulator compute_inverse_exists(const AllEntities& in)
    {
        using View = typename Accumulator::View;
        const Accumulator claim(View(in.q_claim));
        const Accumulator link(View(in.q_link));
        return claim + link - claim * link;
    }

    /** @brief Each index appears once on each side, so the multiplicities are all 1. */
    template <typename Accumulator, size_t index, typename AllEntities>
    static Accumulator lookup_read_counts(const AllEntities&)
    {
        return Accumulator(FF(1));
    }

    template <typename Accumulator, size_t lookup_index, typename AllEntities>
    static Accumulator get_lookup_term_predicate(const AllEntities& in)
    {
        using View = typename Accumulator::View;
        static_assert(lookup_index == 0);
        return Accumulator(View(in.q_claim));
    }

    template <typename Accumulator, size_t table_index, typename AllEntities>
    static Accumulator get_table_term_predicate(const AllEntities& in)
    {
        using View = typename Accumulator::View;
        static_assert(table_index == 0);
        return Accumulator(View(in.q_link));
    }

    template <typename Accumulator, size_t lookup_index, typename AllEntities, typename Parameters>
    static Accumulator compute_lookup_term(const AllEntities& in, const Parameters& params)
    {
        using View = typename Accumulator::View;
        static_assert(lookup_index == 0);
        // The dense link section starts at NUM_DISABLED_ROWS_IN_SUMCHECK — the offset bb's databus
        // polynomials use — so that a consuming Mega circuit's calldata column carries the claims at
        // the same indices and the two commitments can be compared directly.
        return Accumulator(View(in.claim_index)) + FF(NUM_DISABLED_ROWS_IN_SUMCHECK) +
               Accumulator(View(in.claim_value)) * params.beta + params.gamma;
    }

    template <typename Accumulator, size_t table_index, typename AllEntities, typename Parameters>
    static Accumulator compute_table_term(const AllEntities& in, const Parameters& params)
    {
        using View = typename Accumulator::View;
        static_assert(table_index == 0);
        return Accumulator(View(in.row_index)) + Accumulator(View(in.link_value)) * params.beta + params.gamma;
    }

    template <typename ContainerOverSubrelations, typename AllEntities, typename Parameters>
    static void accumulate(ContainerOverSubrelations& accumulator,
                           const AllEntities& in,
                           const Parameters& params,
                           const FF& scaling_factor)
    {
        _accumulate_logderivative_subrelation_contributions<FF,
                                                            Blake3VMLinkRelationImpl<FF>,
                                                            ContainerOverSubrelations,
                                                            AllEntities,
                                                            Parameters,
                                                            false>(accumulator, in, params, scaling_factor);
    }
};

/**
 * @brief Pins each row's claim contribution to the trace data that row already carries.
 *
 * A round-0 G row claims its message pair `mx + 2³²·my`; an output row of a call's final
 * compression claims its output pair `out_{2t} + 2³²·out_{2t+1}`. Everywhere else the claim is
 * zero and `q_claim` is off. `claim_index` counts the claims seen so far, so the gathered order is
 * trace order: a call's message chunks in byte order, then its digest.
 */
template <typename FF_> class Blake3VMClaimRelationImpl {
  public:
    using FF = FF_;

    static constexpr size_t NUM_SUBRELATIONS = 6;
    static constexpr std::array<size_t, NUM_SUBRELATIONS> SUBRELATION_PARTIAL_LENGTHS = [] {
        std::array<size_t, NUM_SUBRELATIONS> lengths{};
        lengths.fill(5);
        return lengths;
    }();

    template <typename ContainerOverSubrelations, typename AllEntities, typename Parameters>
    static void accumulate(ContainerOverSubrelations& evals,
                           const AllEntities& in,
                           [[maybe_unused]] const Parameters& params,
                           const FF& scaling_factor)
    {
        using Accumulator = typename std::tuple_element_t<0, ContainerOverSubrelations>;
        using View = typename Accumulator::View;

        const Accumulator q_round0(View(in.q_round0));
        const Accumulator is_final(View(in.is_final));
        const Accumulator q_out =
            Accumulator(View(in.q_out_0)) + View(in.q_out_1) + View(in.q_out_2) + View(in.q_out_3);
        const Accumulator q_claim(View(in.q_claim));
        const Accumulator one(FF(1));
        const FF two_pow_32 = FF(1UL << 32);

        const auto word = [](const auto& b0, const auto& b1, const auto& b2, const auto& b3) {
            return Accumulator(View(b0)) + Accumulator(View(b1)) * FF(1 << 8) + Accumulator(View(b2)) * FF(1 << 16) +
                   Accumulator(View(b3)) * FF(1UL << 24);
        };
        // On a round-0 row the message pair is the two words the G function reads there.
        const Accumulator message_pair = word(in.mx_b_0, in.mx_b_1, in.mx_b_2, in.mx_b_3) +
                                         word(in.my_b_0, in.my_b_1, in.my_b_2, in.my_b_3) * two_pow_32;
        // On output row t the digest pair is the two output words that row computes, which the
        // wiring relation has already tied to z1/z2.
        const Accumulator digest_pair = word(in.z1_b_0, in.z1_b_1, in.z1_b_2, in.z1_b_3) +
                                        word(in.z2_b_0, in.z2_b_1, in.z2_b_2, in.z2_b_3) * two_pow_32;

        // 0: a row claims exactly when it is a round-0 row or a final compression's output row.
        const Accumulator claimed_out = q_out * is_final;
        std::get<0>(evals) += (q_claim - q_round0 - claimed_out) * scaling_factor;
        // 1: and the claim is that row's pair.
        std::get<1>(evals) += (q_round0 * (Accumulator(View(in.claim_value)) - message_pair) +
                               claimed_out * (Accumulator(View(in.claim_value)) - digest_pair)) *
                              scaling_factor;
        // 2: no other row claims anything.
        std::get<2>(evals) += ((one - q_claim) * Accumulator(View(in.claim_value))) * scaling_factor;
        // 3: the index advances by one per claim and never otherwise. The last row is exempt: its
        // shift wraps off the end of the trace.
        std::get<3>(evals) += (one - Accumulator(View(in.lagrange_last))) *
                              (Accumulator(View(in.claim_index_shift)) - Accumulator(View(in.claim_index)) - q_claim) *
                              scaling_factor;
        // 4/5: both selectors are boolean.
        std::get<4>(evals) += (q_claim * (q_claim - one)) * scaling_factor;
        std::get<5>(evals) += (is_final * (is_final - one)) * scaling_factor;
    }
};

template <typename FF> using Blake3VMLinkRelation = Relation<Blake3VMLinkRelationImpl<FF>>;
template <typename FF> using Blake3VMClaimRelation = Relation<Blake3VMClaimRelationImpl<FF>>;

} // namespace bb
