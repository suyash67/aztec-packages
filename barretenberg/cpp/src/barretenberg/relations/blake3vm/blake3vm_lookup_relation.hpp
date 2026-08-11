#pragma once
/**
 * @brief Log-derivative lookup relations binding the Blake3VM byte columns to the 8-bit XOR table.
 *
 * The table occupies rows [0, 2^16) of the precomputed columns (table_x, table_y, table_z) with
 * table_z = table_x ⊕ table_y, table_x, table_y ∈ [0, 256). A tuple (x, y, z) read from it
 * simultaneously proves z = x ⊕ y and the 8-bit range of all three components. Two derived uses:
 * - k-bit range check (k ≤ 8): read (v·2⁸⁻ᵏ, 0, v·2⁸⁻ᵏ), which is a table row iff v < 2ᵏ.
 * - pure byte-pair range check: read (x, y, j) with witness j = x ⊕ y, discarding j.
 *
 * The 23 reads per G row are split across six lookup sets, each with its own inverse column
 * (inv_0..inv_5) and read-count column (counts_0..counts_5), all writing against the same table:
 *
 *   set 0: z1ᵢ = dinᵢ ⊕ a1ᵢ, i ∈ [0,4)            (also reused on output rows: out_{2t} = v_{2t} ⊕ v_{2t+8})
 *   set 1: z2ᵢ = binᵢ ⊕ c1ᵢ                        (also reused on output rows: out_{2t+1} = v_{2t+1} ⊕ v_{2t+9})
 *   set 2: z3ᵢ = z1_{(i+2) mod 4} ⊕ a2ᵢ            (the shifted index is the rotation of d by 16)
 *   set 3: z4ᵢ = b1ᵢ ⊕ c2ᵢ
 *   set 4: nibble ranges for rot-12, 7-bit range for rot-7, and (mx₀, mx₁) byte range
 *   set 5: (mx₂, mx₃), (my₀, my₁), (my₂, my₃) byte ranges
 *
 * Subrelation 0 (linearly independent) proves each row's inverse is the inverse of the product of
 * all term denominators when any term is active. Subrelation 1 (linearly dependent) is the logup
 * balance Σ reads − Σ counts·writes = 0 over the whole trace. Both are implemented by
 * `_accumulate_logderivative_subrelation_contributions`.
 */

#include "barretenberg/honk/proof_system/logderivative_library.hpp"
#include "barretenberg/relations/relation_types.hpp"

namespace bb {

template <typename FF_, size_t SET> class Blake3VMLookupRelationImpl {
  public:
    using FF = FF_;
    static_assert(SET < 6);

    static constexpr size_t NUM_LOOKUP_TERMS = (SET < 5) ? 4 : 3;
    static constexpr size_t NUM_TABLE_TERMS = 1;
    // 1 + polynomial degree of this relation (matches the ECCVM lookup relation formula).
    static constexpr size_t LENGTH = NUM_LOOKUP_TERMS + NUM_TABLE_TERMS + 3;

    static constexpr std::array<size_t, 2> SUBRELATION_PARTIAL_LENGTHS{ LENGTH, LENGTH };
    static constexpr std::array<bool, 2> SUBRELATION_LINEARLY_INDEPENDENT = { true, false };

    // Sets 0 and 1 read on G rows and output rows; sets 2..5 read on G rows only.
    static constexpr bool READS_ON_OUT_ROWS = (SET < 2);

    template <typename AllValues> static bool operation_exists_at_row(const AllValues& row)
    {
        FF active = row.q_table;
        active += row.q_pos_0 + row.q_pos_1 + row.q_pos_2 + row.q_pos_3 + row.q_pos_4 + row.q_pos_5 + row.q_pos_6 +
                  row.q_pos_7;
        if constexpr (READS_ON_OUT_ROWS) {
            active += row.q_out_0 + row.q_out_1 + row.q_out_2 + row.q_out_3;
        }
        return active != FF(0);
    }

    template <typename AllEntities> static auto& get_inverse_polynomial(AllEntities& in)
    {
        if constexpr (SET == 0) {
            return in.inv_0;
        } else if constexpr (SET == 1) {
            return in.inv_1;
        } else if constexpr (SET == 2) {
            return in.inv_2;
        } else if constexpr (SET == 3) {
            return in.inv_3;
        } else if constexpr (SET == 4) {
            return in.inv_4;
        } else {
            return in.inv_5;
        }
    }

    // 0/1-valued: 1 on rows where this set has an active read.
    template <typename Accumulator, typename AllEntities> static Accumulator read_active(const AllEntities& in)
    {
        using View = typename Accumulator::View;
        Accumulator active = Accumulator(View(in.q_pos_0)) + View(in.q_pos_1) + View(in.q_pos_2) + View(in.q_pos_3) +
                             View(in.q_pos_4) + View(in.q_pos_5) + View(in.q_pos_6) + View(in.q_pos_7);
        if constexpr (READS_ON_OUT_ROWS) {
            active += Accumulator(View(in.q_out_0)) + View(in.q_out_1) + View(in.q_out_2) + View(in.q_out_3);
        }
        return active;
    }

    template <typename Accumulator, typename AllEntities>
    static Accumulator compute_inverse_exists(const AllEntities& in)
    {
        using View = typename Accumulator::View;
        const Accumulator active = read_active<Accumulator>(in);
        const Accumulator table(View(in.q_table));
        // Logical OR of two 0/1 selectors (G/output rows overlap the table region).
        return active + table - active * table;
    }

    template <typename Accumulator, size_t index, typename AllEntities>
    static Accumulator lookup_read_counts(const AllEntities& in)
    {
        using View = typename Accumulator::View;
        static_assert(index == 0);
        if constexpr (SET == 0) {
            return Accumulator(View(in.counts_0));
        } else if constexpr (SET == 1) {
            return Accumulator(View(in.counts_1));
        } else if constexpr (SET == 2) {
            return Accumulator(View(in.counts_2));
        } else if constexpr (SET == 3) {
            return Accumulator(View(in.counts_3));
        } else if constexpr (SET == 4) {
            return Accumulator(View(in.counts_4));
        } else {
            return Accumulator(View(in.counts_5));
        }
    }

    template <typename Accumulator, size_t lookup_index, typename AllEntities>
    static Accumulator get_lookup_term_predicate(const AllEntities& in)
    {
        static_assert(lookup_index < NUM_LOOKUP_TERMS);
        return read_active<Accumulator>(in);
    }

    template <typename Accumulator, size_t table_index, typename AllEntities>
    static Accumulator get_table_term_predicate(const AllEntities& in)
    {
        using View = typename Accumulator::View;
        static_assert(table_index == 0);
        return Accumulator(View(in.q_table));
    }

    // Fingerprint of the tuple (x, y, z): x + β·y + β²·z + γ.
    template <typename Accumulator>
    static Accumulator tuple_term(Accumulator x, Accumulator y, Accumulator z, const auto& params)
    {
        return x + y * params.beta + z * params.beta_sqr + params.gamma;
    }

    template <typename Accumulator, size_t lookup_index, typename AllEntities, typename Parameters>
    static Accumulator compute_lookup_term(const AllEntities& in, const Parameters& params)
    {
        using View = typename Accumulator::View;
        static_assert(lookup_index < NUM_LOOKUP_TERMS);
        constexpr size_t i = lookup_index;

        if constexpr (SET == 0) {
            const std::array<View, 4> x{ View(in.din_b_0), View(in.din_b_1), View(in.din_b_2), View(in.din_b_3) };
            const std::array<View, 4> y{ View(in.a1_b_0), View(in.a1_b_1), View(in.a1_b_2), View(in.a1_b_3) };
            const std::array<View, 4> z{ View(in.z1_b_0), View(in.z1_b_1), View(in.z1_b_2), View(in.z1_b_3) };
            return tuple_term(Accumulator(x[i]), Accumulator(y[i]), Accumulator(z[i]), params);
        } else if constexpr (SET == 1) {
            const std::array<View, 4> x{ View(in.bin_b_0), View(in.bin_b_1), View(in.bin_b_2), View(in.bin_b_3) };
            const std::array<View, 4> y{ View(in.c1_b_0), View(in.c1_b_1), View(in.c1_b_2), View(in.c1_b_3) };
            const std::array<View, 4> z{ View(in.z2_b_0), View(in.z2_b_1), View(in.z2_b_2), View(in.z2_b_3) };
            return tuple_term(Accumulator(x[i]), Accumulator(y[i]), Accumulator(z[i]), params);
        } else if constexpr (SET == 2) {
            // x component is z1 rotated by 16 bits: byte (i+2) mod 4.
            const std::array<View, 4> x{ View(in.z1_b_2), View(in.z1_b_3), View(in.z1_b_0), View(in.z1_b_1) };
            const std::array<View, 4> y{ View(in.a2_b_0), View(in.a2_b_1), View(in.a2_b_2), View(in.a2_b_3) };
            const std::array<View, 4> z{ View(in.z3_b_0), View(in.z3_b_1), View(in.z3_b_2), View(in.z3_b_3) };
            return tuple_term(Accumulator(x[i]), Accumulator(y[i]), Accumulator(z[i]), params);
        } else if constexpr (SET == 3) {
            const std::array<View, 4> x{ View(in.b1_b_0), View(in.b1_b_1), View(in.b1_b_2), View(in.b1_b_3) };
            const std::array<View, 4> y{ View(in.c2_b_0), View(in.c2_b_1), View(in.c2_b_2), View(in.c2_b_3) };
            const std::array<View, 4> z{ View(in.z4_b_0), View(in.z4_b_1), View(in.z4_b_2), View(in.z4_b_3) };
            return tuple_term(Accumulator(x[i]), Accumulator(y[i]), Accumulator(z[i]), params);
        } else if constexpr (SET == 4) {
            if constexpr (i == 0) {
                const Accumulator scaled = Accumulator(View(in.z2_nib_lo)) * FF(16);
                return tuple_term(scaled, Accumulator(FF(0)), scaled, params);
            } else if constexpr (i == 1) {
                const Accumulator scaled = Accumulator(View(in.z2_nib_hi)) * FF(16);
                return tuple_term(scaled, Accumulator(FF(0)), scaled, params);
            } else if constexpr (i == 2) {
                const Accumulator scaled = Accumulator(View(in.z4_lo7)) * FF(2);
                return tuple_term(scaled, Accumulator(FF(0)), scaled, params);
            } else {
                return tuple_term(Accumulator(View(in.mx_b_0)),
                                  Accumulator(View(in.mx_b_1)),
                                  Accumulator(View(in.junk_mx01)),
                                  params);
            }
        } else {
            if constexpr (i == 0) {
                return tuple_term(Accumulator(View(in.mx_b_2)),
                                  Accumulator(View(in.mx_b_3)),
                                  Accumulator(View(in.junk_mx23)),
                                  params);
            } else if constexpr (i == 1) {
                return tuple_term(Accumulator(View(in.my_b_0)),
                                  Accumulator(View(in.my_b_1)),
                                  Accumulator(View(in.junk_my01)),
                                  params);
            } else {
                return tuple_term(Accumulator(View(in.my_b_2)),
                                  Accumulator(View(in.my_b_3)),
                                  Accumulator(View(in.junk_my23)),
                                  params);
            }
        }
    }

    template <typename Accumulator, size_t table_index, typename AllEntities, typename Parameters>
    static Accumulator compute_table_term(const AllEntities& in, const Parameters& params)
    {
        using View = typename Accumulator::View;
        static_assert(table_index == 0);
        return tuple_term(
            Accumulator(View(in.table_x)), Accumulator(View(in.table_y)), Accumulator(View(in.table_z)), params);
    }

    template <typename ContainerOverSubrelations, typename AllEntities, typename Parameters>
    static void accumulate(ContainerOverSubrelations& accumulator,
                           const AllEntities& in,
                           const Parameters& params,
                           const FF& scaling_factor)
    {
        _accumulate_logderivative_subrelation_contributions<FF,
                                                            Blake3VMLookupRelationImpl<FF, SET>,
                                                            ContainerOverSubrelations,
                                                            AllEntities,
                                                            Parameters,
                                                            false>(accumulator, in, params, scaling_factor);
    }
};

template <typename FF, size_t SET> using Blake3VMLookupRelation = Relation<Blake3VMLookupRelationImpl<FF, SET>>;

} // namespace bb
