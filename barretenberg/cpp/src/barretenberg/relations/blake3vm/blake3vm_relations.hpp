#pragma once
/**
 * @brief Arithmetic and wiring relations of the Blake3VM.
 *
 * Notation (all little-endian byte limbs, ⌊w⌋ᵢ the i-th byte of the 32-bit word w):
 * - A G row at position p reads state words a = v_{A[p]}, b = v_{B[p]}, c = v_{C[p]}, d = v_{D[p]}
 *   and message words mx = m_{2p}, my = m_{2p+1}, and computes
 *       a1 = a + b + mx,  d1 = (d ⊕ a1) ⋙ 16,  c1 = c + d1,  b1 = (b ⊕ c1) ⋙ 12,
 *       a2 = a1 + b1 + my, d2 = (d1 ⊕ a2) ⋙ 8, c2 = c1 + d2, b2 = (b1 ⊕ c2) ⋙ 7.
 * - Additions mod 2³² are proven as x + y (+ z) = result + 2³²·carry with carry ∈ {0,1} for two
 *   terms and carry ∈ {0,1,2} for three; the result's byte decomposition (range-checked by the XOR
 *   lookups of blake3vm_lookup_relation.hpp) bounds it below 2³².
 * - Rotations by 16/8 permute byte limbs. Rotation by 12 splits z2's byte 1 into nibbles
 *   (z2_b_1 = z2_nib_lo + 16·z2_nib_hi); rotation by 7 splits z4's byte 0 as
 *   z4_b_0 = z4_lo7 + 128·z4_msb. The split pieces are range-checked by lookup set 4.
 *
 * Blake3VMGRelation carries the per-G-row algebra. Blake3VMWiringRelation carries state flow:
 * v_j updates through its fixed G role (j < 4 ↦ a2, j ∈ [4,8) ↦ b2, j ∈ [8,12) ↦ c2,
 * j ∈ [12,16) ↦ d2), message columns copy within a round and permute by MSG_PERM across rounds,
 * output rows compute outᵢ = vᵢ ⊕ vᵢ₊₈ two words at a time reusing the set-0/1 lookup columns, and
 * the boundary row pins the next row's state to IV / chained CV / block_len / flags.
 * Blake3VMZeroRowRelation forces every to-be-shifted column to 0 at row 0, which the shifted
 * polynomial views require.
 */

#include "barretenberg/blake3vm/blake3vm_circuit_builder.hpp"
#include "barretenberg/common/constexpr_utils.hpp"
#include "barretenberg/relations/relation_types.hpp"

namespace bb {

template <typename FF_> class Blake3VMGRelationImpl {
  public:
    using FF = FF_;
    using Layout = Blake3VMTraceLayout;

    static constexpr size_t NUM_SUBRELATIONS = 22;
    // Uniform length: the carry range checks q_g·c·(c−1)·(c−2) have degree 4.
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

        const std::array<View, 8> q_pos{ View(in.q_pos_0), View(in.q_pos_1), View(in.q_pos_2), View(in.q_pos_3),
                                         View(in.q_pos_4), View(in.q_pos_5), View(in.q_pos_6), View(in.q_pos_7) };
        Accumulator q_g(q_pos[0]);
        for (size_t p = 1; p < 8; ++p) {
            q_g += q_pos[p];
        }

        const std::array<View, 16> v{ View(in.v_0), View(in.v_1), View(in.v_2),  View(in.v_3),  View(in.v_4),  View(in.v_5),  View(in.v_6),  View(in.v_7),
                                      View(in.v_8), View(in.v_9), View(in.v_10), View(in.v_11), View(in.v_12), View(in.v_13), View(in.v_14), View(in.v_15) };
        const std::array<View, 16> m{ View(in.m_0), View(in.m_1), View(in.m_2),  View(in.m_3),  View(in.m_4),  View(in.m_5),  View(in.m_6),  View(in.m_7),
                                      View(in.m_8), View(in.m_9), View(in.m_10), View(in.m_11), View(in.m_12), View(in.m_13), View(in.m_14), View(in.m_15) };

        const auto word = [](const auto& b0, const auto& b1, const auto& b2, const auto& b3) {
            return Accumulator(View(b0)) + Accumulator(View(b1)) * FF(1 << 8) +
                   Accumulator(View(b2)) * FF(1 << 16) + Accumulator(View(b3)) * FF(1UL << 24);
        };
        // Σₚ q_pos_p · v_{TABLE[p]}
        const auto select_state = [&](const std::array<size_t, 8>& table) {
            Accumulator result = Accumulator(q_pos[0]) * v[table[0]];
            for (size_t p = 1; p < 8; ++p) {
                result += Accumulator(q_pos[p]) * v[table[p]];
            }
            return result;
        };
        const auto select_msg = [&](size_t parity) {
            Accumulator result = Accumulator(q_pos[0]) * m[parity];
            for (size_t p = 1; p < 8; ++p) {
                result += Accumulator(q_pos[p]) * m[2 * p + parity];
            }
            return result;
        };

        const Accumulator din_word = word(in.din_b_0, in.din_b_1, in.din_b_2, in.din_b_3);
        const Accumulator bin_word = word(in.bin_b_0, in.bin_b_1, in.bin_b_2, in.bin_b_3);
        const Accumulator mx_word = word(in.mx_b_0, in.mx_b_1, in.mx_b_2, in.mx_b_3);
        const Accumulator my_word = word(in.my_b_0, in.my_b_1, in.my_b_2, in.my_b_3);
        const Accumulator a1_word = word(in.a1_b_0, in.a1_b_1, in.a1_b_2, in.a1_b_3);
        const Accumulator c1_word = word(in.c1_b_0, in.c1_b_1, in.c1_b_2, in.c1_b_3);
        const Accumulator b1_word = word(in.b1_b_0, in.b1_b_1, in.b1_b_2, in.b1_b_3);
        const Accumulator a2_word = word(in.a2_b_0, in.a2_b_1, in.a2_b_2, in.a2_b_3);
        const Accumulator c2_word = word(in.c2_b_0, in.c2_b_1, in.c2_b_2, in.c2_b_3);
        // d1 = z1 ⋙ 16 and d2 = z3 ⋙ 8, as byte permutations.
        const Accumulator d1_word = word(in.z1_b_2, in.z1_b_3, in.z1_b_0, in.z1_b_1);
        const Accumulator d2_word = word(in.z3_b_1, in.z3_b_2, in.z3_b_3, in.z3_b_0);

        const Accumulator a1(View(in.a1));
        const Accumulator c1(View(in.c1));
        const Accumulator b1(View(in.b1));
        const Accumulator a2(View(in.a2));
        const Accumulator c2(View(in.c2));
        const Accumulator b2(View(in.b2));
        const Accumulator ca1(View(in.ca1));
        const Accumulator cc1(View(in.cc1));
        const Accumulator ca2(View(in.ca2));
        const Accumulator cc2(View(in.cc2));

        const FF two_pow_32 = FF(1UL << 32);
        const Accumulator one(FF(1));
        const Accumulator two(FF(2));

        // 0/1: operand byte decompositions of d and b.
        std::get<0>(evals) += (din_word * q_g - select_state(Layout::POS_D)) * scaling_factor;
        std::get<1>(evals) += (bin_word * q_g - select_state(Layout::POS_B)) * scaling_factor;
        // 2/3: message word selection; ties the byte-decomposed mx/my to the scheduled m columns.
        std::get<2>(evals) += (mx_word * q_g - select_msg(0)) * scaling_factor;
        std::get<3>(evals) += (my_word * q_g - select_msg(1)) * scaling_factor;
        // 4-6: a1 = a + b + mx (mod 2³²).
        std::get<4>(evals) +=
            (select_state(Layout::POS_A) + select_state(Layout::POS_B) + (mx_word - a1 - ca1 * two_pow_32) * q_g) *
            scaling_factor;
        std::get<5>(evals) += q_g * ca1 * (ca1 - one) * (ca1 - two) * scaling_factor;
        std::get<6>(evals) += q_g * (a1_word - a1) * scaling_factor;
        // 7-9: c1 = c + d1 (mod 2³²).
        std::get<7>(evals) += (select_state(Layout::POS_C) + (d1_word - c1 - cc1 * two_pow_32) * q_g) * scaling_factor;
        std::get<8>(evals) += q_g * cc1 * (cc1 - one) * scaling_factor;
        std::get<9>(evals) += q_g * (c1_word - c1) * scaling_factor;
        // 10-12: b1 = (b ⊕ c1) ⋙ 12 via the nibble split of z2's byte 1.
        std::get<10>(evals) += q_g *
                               (Accumulator(View(in.z2_b_1)) - Accumulator(View(in.z2_nib_lo)) -
                                Accumulator(View(in.z2_nib_hi)) * FF(16)) *
                               scaling_factor;
        const Accumulator rot12 = Accumulator(View(in.z2_nib_hi)) + Accumulator(View(in.z2_b_2)) * FF(1 << 4) +
                                  Accumulator(View(in.z2_b_3)) * FF(1 << 12) +
                                  Accumulator(View(in.z2_b_0)) * FF(1UL << 20) +
                                  Accumulator(View(in.z2_nib_lo)) * FF(1UL << 28);
        std::get<11>(evals) += q_g * (b1 - rot12) * scaling_factor;
        std::get<12>(evals) += q_g * (b1_word - b1) * scaling_factor;
        // 13-15: a2 = a1 + b1 + my (mod 2³²).
        std::get<13>(evals) += q_g * (a1 + b1 + my_word - a2 - ca2 * two_pow_32) * scaling_factor;
        std::get<14>(evals) += q_g * ca2 * (ca2 - one) * (ca2 - two) * scaling_factor;
        std::get<15>(evals) += q_g * (a2_word - a2) * scaling_factor;
        // 16-18: c2 = c1 + d2 (mod 2³²).
        std::get<16>(evals) += q_g * (c1 + d2_word - c2 - cc2 * two_pow_32) * scaling_factor;
        std::get<17>(evals) += q_g * cc2 * (cc2 - one) * scaling_factor;
        std::get<18>(evals) += q_g * (c2_word - c2) * scaling_factor;
        // 19-21: b2 = (b1 ⊕ c2) ⋙ 7 via the 7-bit/msb split of z4's byte 0.
        std::get<19>(evals) +=
            q_g *
            (Accumulator(View(in.z4_b_0)) - Accumulator(View(in.z4_lo7)) - Accumulator(View(in.z4_msb)) * FF(128)) *
            scaling_factor;
        std::get<20>(evals) +=
            q_g * Accumulator(View(in.z4_msb)) * (Accumulator(View(in.z4_msb)) - one) * scaling_factor;
        const Accumulator rot7 = Accumulator(View(in.z4_msb)) + Accumulator(View(in.z4_b_1)) * FF(2) +
                                 Accumulator(View(in.z4_b_2)) * FF(1 << 9) +
                                 Accumulator(View(in.z4_b_3)) * FF(1UL << 17) +
                                 Accumulator(View(in.z4_lo7)) * FF(1UL << 25);
        std::get<21>(evals) += q_g * (b2 - rot7) * scaling_factor;
    }
};

template <typename FF_> class Blake3VMWiringRelationImpl {
  public:
    using FF = FF_;
    using Layout = Blake3VMTraceLayout;

    static constexpr size_t NUM_SUBRELATIONS = 48;
    // Uniform length: the boundary terms q_out_3·chain·out_j have degree 3.
    static constexpr std::array<size_t, NUM_SUBRELATIONS> SUBRELATION_PARTIAL_LENGTHS = [] {
        std::array<size_t, NUM_SUBRELATIONS> lengths{};
        lengths.fill(4);
        return lengths;
    }();

    // For each state word j, the two G positions that update it (through its fixed role).
    static constexpr std::array<std::array<size_t, 2>, 16> UPDATE_POS = { {
        { 0, 4 },
        { 1, 5 },
        { 2, 6 },
        { 3, 7 }, // j ∈ [0,4): updated as a2
        { 0, 7 },
        { 1, 4 },
        { 2, 5 },
        { 3, 6 }, // j ∈ [4,8): updated as b2
        { 0, 6 },
        { 1, 7 },
        { 2, 4 },
        { 3, 5 }, // j ∈ [8,12): updated as c2
        { 0, 5 },
        { 1, 6 },
        { 2, 7 },
        { 3, 4 }, // j ∈ [12,16): updated as d2
    } };

    template <typename ContainerOverSubrelations, typename AllEntities, typename Parameters>
    static void accumulate(ContainerOverSubrelations& evals,
                           const AllEntities& in,
                           [[maybe_unused]] const Parameters& params,
                           const FF& scaling_factor)
    {
        using Accumulator = typename std::tuple_element_t<0, ContainerOverSubrelations>;
        using View = typename Accumulator::View;

        const std::array<View, 8> q_pos{ View(in.q_pos_0), View(in.q_pos_1), View(in.q_pos_2), View(in.q_pos_3),
                                         View(in.q_pos_4), View(in.q_pos_5), View(in.q_pos_6), View(in.q_pos_7) };
        Accumulator q_g(q_pos[0]);
        for (size_t p = 1; p < 8; ++p) {
            q_g += q_pos[p];
        }
        // Message columns copy on positions 0..6; on position 7 they permute into the next round
        // where q_mperm is set (all pos-7 rows except the final round's).
        Accumulator q_mcopy(q_pos[0]);
        for (size_t p = 1; p < 7; ++p) {
            q_mcopy += q_pos[p];
        }
        const std::array<View, 4> q_out{ View(in.q_out_0), View(in.q_out_1), View(in.q_out_2), View(in.q_out_3) };
        const Accumulator q_vcopy = Accumulator(q_out[0]) + q_out[1] + q_out[2];
        const Accumulator q_boundary(q_out[3]);
        const Accumulator q_mperm(View(in.q_mperm));
        const Accumulator chain(View(in.chain));
        const Accumulator lagrange_first(View(in.lagrange_first));

        const std::array<View, 16> v{ View(in.v_0), View(in.v_1), View(in.v_2),  View(in.v_3),  View(in.v_4),  View(in.v_5),  View(in.v_6),  View(in.v_7),
                                      View(in.v_8), View(in.v_9), View(in.v_10), View(in.v_11), View(in.v_12), View(in.v_13), View(in.v_14), View(in.v_15) };
        const std::array<View, 16> v_shift{ View(in.v_0_shift),  View(in.v_1_shift),  View(in.v_2_shift),  View(in.v_3_shift),
                                            View(in.v_4_shift),  View(in.v_5_shift),  View(in.v_6_shift),  View(in.v_7_shift),
                                            View(in.v_8_shift),  View(in.v_9_shift),  View(in.v_10_shift), View(in.v_11_shift),
                                            View(in.v_12_shift), View(in.v_13_shift), View(in.v_14_shift), View(in.v_15_shift) };
        const std::array<View, 16> m{ View(in.m_0), View(in.m_1), View(in.m_2),  View(in.m_3),  View(in.m_4),  View(in.m_5),  View(in.m_6),  View(in.m_7),
                                      View(in.m_8), View(in.m_9), View(in.m_10), View(in.m_11), View(in.m_12), View(in.m_13), View(in.m_14), View(in.m_15) };
        const std::array<View, 16> m_shift{ View(in.m_0_shift),  View(in.m_1_shift),  View(in.m_2_shift),  View(in.m_3_shift),
                                            View(in.m_4_shift),  View(in.m_5_shift),  View(in.m_6_shift),  View(in.m_7_shift),
                                            View(in.m_8_shift),  View(in.m_9_shift),  View(in.m_10_shift), View(in.m_11_shift),
                                            View(in.m_12_shift), View(in.m_13_shift), View(in.m_14_shift), View(in.m_15_shift) };
        const std::array<View, 8> out{ View(in.out_0), View(in.out_1), View(in.out_2), View(in.out_3),
                                       View(in.out_4), View(in.out_5), View(in.out_6), View(in.out_7) };
        const std::array<View, 8> out_shift{ View(in.out_0_shift), View(in.out_1_shift), View(in.out_2_shift), View(in.out_3_shift),
                                             View(in.out_4_shift), View(in.out_5_shift), View(in.out_6_shift), View(in.out_7_shift) };

        const auto word = [](const auto& b0, const auto& b1, const auto& b2, const auto& b3) {
            return Accumulator(View(b0)) + Accumulator(View(b1)) * FF(1 << 8) +
                   Accumulator(View(b2)) * FF(1 << 16) + Accumulator(View(b3)) * FF(1UL << 24);
        };
        const Accumulator d2_word = word(in.z3_b_1, in.z3_b_2, in.z3_b_3, in.z3_b_0);
        const std::array<Accumulator, 4> role_out{
            Accumulator(View(in.a2)), Accumulator(View(in.b2)), Accumulator(View(in.c2)), d2_word
        };

        // 0-15: state flow. On a G row, v_j either takes its role output (at its two update
        // positions) or copies; output rows 0-2 copy; the boundary row pins the next compression's
        // initial state.
        bb::constexpr_for<0, 16, 1>([&]<size_t J>() {
            const Accumulator q_update = Accumulator(q_pos[UPDATE_POS[J][0]]) + q_pos[UPDATE_POS[J][1]];
            const Accumulator q_copy = q_g - q_update + q_vcopy;
            const Accumulator shift(v_shift[J]);

            Accumulator boundary_value;
            if constexpr (J < 8) {
                // Initial CV: IV for a fresh hash, previous output words when chaining.
                boundary_value = chain * out[J] + (Accumulator(FF(1)) - chain) * FF(blake3::IV[J]);
            } else if constexpr (J < 12) {
                boundary_value = Accumulator(FF(blake3::IV[J - 8]));
            } else if constexpr (J < 14) {
                boundary_value = Accumulator(FF(0)); // single-chunk counter is always 0
            } else if constexpr (J == 14) {
                boundary_value = Accumulator(View(in.blk_len));
            } else {
                boundary_value = Accumulator(View(in.blk_flags));
            }

            const Accumulator constraint =
                q_update * (shift - role_out[J / 4]) + q_copy * (shift - v[J]) + q_boundary * (shift - boundary_value);
            std::get<J>(evals) += constraint * scaling_factor;
        });

        // 16-31: message flow (copy within a round, MSG_PERM across rounds).
        bb::constexpr_for<0, 16, 1>([&]<size_t J>() {
            const Accumulator constraint = q_mcopy * (Accumulator(m_shift[J]) - m[J]) +
                                           q_mperm * (Accumulator(m_shift[J]) - m[Layout::MSG_PERM[J]]);
            std::get<16 + J>(evals) += constraint * scaling_factor;
        });

        // 32-39: output words stay constant across the four output rows.
        bb::constexpr_for<0, 8, 1>([&]<size_t J>() {
            std::get<32 + J>(evals) += q_vcopy * (Accumulator(out_shift[J]) - out[J]) * scaling_factor;
        });

        // 40-45: output-row byte ties. Row t computes out_{2t} = v_{2t} ⊕ v_{2t+8} in the set-0
        // lookup columns (din/a1/z1) and out_{2t+1} = v_{2t+1} ⊕ v_{2t+9} in the set-1 columns.
        const Accumulator din_word = word(in.din_b_0, in.din_b_1, in.din_b_2, in.din_b_3);
        const Accumulator bin_word = word(in.bin_b_0, in.bin_b_1, in.bin_b_2, in.bin_b_3);
        const Accumulator a1_word = word(in.a1_b_0, in.a1_b_1, in.a1_b_2, in.a1_b_3);
        const Accumulator c1_word = word(in.c1_b_0, in.c1_b_1, in.c1_b_2, in.c1_b_3);
        const Accumulator z1_word = word(in.z1_b_0, in.z1_b_1, in.z1_b_2, in.z1_b_3);
        const Accumulator z2_word = word(in.z2_b_0, in.z2_b_1, in.z2_b_2, in.z2_b_3);
        Accumulator tie_out_even(FF(0));
        Accumulator tie_out_odd(FF(0));
        Accumulator tie_v_even(FF(0));
        Accumulator tie_v_odd(FF(0));
        Accumulator tie_vh_even(FF(0));
        Accumulator tie_vh_odd(FF(0));
        for (size_t t = 0; t < 4; ++t) {
            const Accumulator q_t(q_out[t]);
            tie_out_even += q_t * (Accumulator(out[2 * t]) - z1_word);
            tie_out_odd += q_t * (Accumulator(out[2 * t + 1]) - z2_word);
            tie_v_even += q_t * (din_word - v[2 * t]);
            tie_v_odd += q_t * (bin_word - v[2 * t + 1]);
            tie_vh_even += q_t * (a1_word - v[2 * t + 8]);
            tie_vh_odd += q_t * (c1_word - v[2 * t + 9]);
        }
        std::get<40>(evals) += tie_out_even * scaling_factor;
        std::get<41>(evals) += tie_out_odd * scaling_factor;
        std::get<42>(evals) += tie_v_even * scaling_factor;
        std::get<43>(evals) += tie_v_odd * scaling_factor;
        std::get<44>(evals) += tie_vh_even * scaling_factor;
        std::get<45>(evals) += tie_vh_odd * scaling_factor;

        // 46: chain is boolean on boundary rows. 47: the genesis compression cannot chain.
        std::get<46>(evals) += q_boundary * chain * (chain - Accumulator(FF(1))) * scaling_factor;
        std::get<47>(evals) += lagrange_first * chain * scaling_factor;
    }
};

/** @brief Forces every to-be-shifted column to zero at row 0 (required by the shifted views). */
template <typename FF_> class Blake3VMZeroRowRelationImpl {
  public:
    using FF = FF_;

    static constexpr size_t NUM_SUBRELATIONS = 41;
    static constexpr std::array<size_t, NUM_SUBRELATIONS> SUBRELATION_PARTIAL_LENGTHS = [] {
        std::array<size_t, NUM_SUBRELATIONS> lengths{};
        lengths.fill(3);
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

        const Accumulator lagrange_first(View(in.lagrange_first));
        const std::array<View, 41> to_be_shifted{ View(in.v_0),   View(in.v_1),   View(in.v_2),   View(in.v_3),   View(in.v_4),   View(in.v_5),   View(in.v_6),
                                                  View(in.v_7),   View(in.v_8),   View(in.v_9),   View(in.v_10),  View(in.v_11),  View(in.v_12),  View(in.v_13),
                                                  View(in.v_14),  View(in.v_15),  View(in.m_0),   View(in.m_1),   View(in.m_2),   View(in.m_3),   View(in.m_4),
                                                  View(in.m_5),   View(in.m_6),   View(in.m_7),   View(in.m_8),   View(in.m_9),   View(in.m_10),  View(in.m_11),
                                                  View(in.m_12),  View(in.m_13),  View(in.m_14),  View(in.m_15),  View(in.out_0), View(in.out_1), View(in.out_2),
                                                  View(in.out_3), View(in.out_4), View(in.out_5), View(in.out_6), View(in.out_7),
                                                  View(in.claim_index) };
        bb::constexpr_for<0, 41, 1>(
            [&]<size_t J>() { std::get<J>(evals) += lagrange_first * Accumulator(to_be_shifted[J]) * scaling_factor; });
    }
};

template <typename FF> using Blake3VMGRelation = Relation<Blake3VMGRelationImpl<FF>>;
template <typename FF> using Blake3VMWiringRelation = Relation<Blake3VMWiringRelationImpl<FF>>;
template <typename FF> using Blake3VMZeroRowRelation = Relation<Blake3VMZeroRowRelationImpl<FF>>;

} // namespace bb
