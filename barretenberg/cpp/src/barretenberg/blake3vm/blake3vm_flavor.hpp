#pragma once
/**
 * @brief Honk flavor for the Blake3VM: BN254, KZG, non-ZK sumcheck.
 *
 * Entities: 19 precomputed (row-type selectors + the 8-bit XOR table), 125 witness (40 to-be-shifted
 * state/message/output words, 79 per-row working columns, 6 logup inverses), 40 shifts. Relations:
 * the G-function algebra, the state/message/output wiring, the zero-row constraints, and six XOR
 * lookup relations. See blake3vm_circuit_builder.hpp for the trace layout and README.md for the
 * design rationale.
 */

#include "barretenberg/blake3vm/blake3vm_circuit_builder.hpp"
#include "barretenberg/commitment_schemes/commitment_key.hpp"
#include "barretenberg/commitment_schemes/kzg/kzg.hpp"
#include "barretenberg/constants.hpp"
#include "barretenberg/ecc/curves/bn254/bn254.hpp"
#include "barretenberg/flavor/flavor.hpp"
#include "barretenberg/flavor/flavor_macros.hpp"
#include "barretenberg/flavor/partially_evaluated_multivariates.hpp"
#include "barretenberg/numeric/bitop/get_msb.hpp"
#include "barretenberg/polynomials/polynomial.hpp"
#include "barretenberg/polynomials/univariate.hpp"
#include "barretenberg/relations/blake3vm/blake3vm_link_relation.hpp"
#include "barretenberg/relations/blake3vm/blake3vm_lookup_relation.hpp"
#include "barretenberg/relations/blake3vm/blake3vm_relations.hpp"
#include "barretenberg/relations/relation_parameters.hpp"
#include "barretenberg/relations/relation_tuple_helpers.hpp"
#include "barretenberg/transcript/transcript.hpp"

namespace bb {

class Blake3VMFlavor {
  public:
    using CircuitBuilder = Blake3VMCircuitBuilder;
    using Layout = Blake3VMTraceLayout;
    using Curve = curve::BN254;
    using PCS = KZG<Curve>;
    using GroupElement = Curve::Element;
    using Commitment = Curve::AffineElement;
    using CommitmentKey = bb::CommitmentKey<Curve>;
    using VerifierCommitmentKey = bb::VerifierCommitmentKey<Curve>;
    using FF = Curve::ScalarField;
    using BF = Curve::BaseField;
    using Polynomial = bb::Polynomial<FF>;
    using Codec = FrCodec;
    using HashFunction = crypto::Poseidon2<crypto::Poseidon2Bn254ScalarFieldParams>;
    using Transcript = BaseTranscript<Codec, HashFunction>;
    using Proof = HonkProof;

    static constexpr bool USE_SHORT_MONOMIALS = false;
    static constexpr bool HasZK = false;
    static constexpr size_t TRACE_OFFSET = 0;
    static constexpr bool USE_PADDING = false;

    template <typename FF>
    using Relations_ = std::tuple<Blake3VMGRelation<FF>,
                                  Blake3VMWiringRelation<FF>,
                                  Blake3VMZeroRowRelation<FF>,
                                  Blake3VMClaimRelation<FF>,
                                  Blake3VMLinkRelation<FF>,
                                  Blake3VMLookupRelation<FF, 0>,
                                  Blake3VMLookupRelation<FF, 1>,
                                  Blake3VMLookupRelation<FF, 2>,
                                  Blake3VMLookupRelation<FF, 3>,
                                  Blake3VMLookupRelation<FF, 4>,
                                  Blake3VMLookupRelation<FF, 5>>;
    using Relations = Relations_<FF>;

    static constexpr size_t NUM_SUBRELATIONS = compute_number_of_subrelations<Relations>();
    using SubrelationSeparators = std::array<FF, NUM_SUBRELATIONS - 1>;

    static constexpr size_t MAX_PARTIAL_RELATION_LENGTH = compute_max_partial_relation_length<Relations>();
    // Non-ZK: one extra degree for the pow (gate-separator) factor.
    static constexpr size_t BATCHED_RELATION_PARTIAL_LENGTH = MAX_PARTIAL_RELATION_LENGTH + 1;
    static constexpr size_t NUM_RELATIONS = std::tuple_size_v<Relations>;

    template <typename DataType_> class PrecomputedEntities {
      public:
        bool operator==(const PrecomputedEntities& other) const = default;
        using DataType = DataType_;
        DEFINE_FLAVOR_MEMBERS(DataType,
                              lagrange_first, // 1 at row 0
                              lagrange_last,  // 1 at the last row
                              q_pos_0,        // G-row selectors by position within a round
                              q_pos_1,
                              q_pos_2,
                              q_pos_3,
                              q_pos_4,
                              q_pos_5,
                              q_pos_6,
                              q_pos_7,
                              q_out_0, // output-row selectors; q_out_3 rows are boundary rows
                              q_out_1,
                              q_out_2,
                              q_out_3,
                              q_mperm,   // pos-7 G rows of rounds 0..5 (message permutation applies)
                              q_round0,  // round-0 G rows, whose message pair is a claim
                              row_index, // 0, 1, 2, ... — the dense link section's index
                              q_table,   // 1 on the 2^16 XOR-table rows
                              table_x,
                              table_y,
                              table_z)
    };

    /** @brief Witness columns read through shifts: state, message, and output words. */
    template <typename DataType> class WireToBeShiftedEntities {
      public:
        DEFINE_FLAVOR_MEMBERS(DataType,
                              v_0,
                              v_1,
                              v_2,
                              v_3,
                              v_4,
                              v_5,
                              v_6,
                              v_7,
                              v_8,
                              v_9,
                              v_10,
                              v_11,
                              v_12,
                              v_13,
                              v_14,
                              v_15,
                              m_0,
                              m_1,
                              m_2,
                              m_3,
                              m_4,
                              m_5,
                              m_6,
                              m_7,
                              m_8,
                              m_9,
                              m_10,
                              m_11,
                              m_12,
                              m_13,
                              m_14,
                              m_15,
                              out_0,
                              out_1,
                              out_2,
                              out_3,
                              out_4,
                              out_5,
                              out_6,
                              out_7,
                              claim_index)
    };

    /** @brief Per-row working columns of the G rows and output rows. */
    template <typename DataType> class WireNonShiftedEntities {
      public:
        DEFINE_FLAVOR_MEMBERS(DataType,
                              din_b_0,
                              din_b_1,
                              din_b_2,
                              din_b_3,
                              bin_b_0,
                              bin_b_1,
                              bin_b_2,
                              bin_b_3,
                              a1,
                              a1_b_0,
                              a1_b_1,
                              a1_b_2,
                              a1_b_3,
                              ca1,
                              z1_b_0,
                              z1_b_1,
                              z1_b_2,
                              z1_b_3,
                              c1,
                              c1_b_0,
                              c1_b_1,
                              c1_b_2,
                              c1_b_3,
                              cc1,
                              z2_b_0,
                              z2_b_1,
                              z2_b_2,
                              z2_b_3,
                              z2_nib_lo,
                              z2_nib_hi,
                              b1,
                              b1_b_0,
                              b1_b_1,
                              b1_b_2,
                              b1_b_3,
                              a2,
                              a2_b_0,
                              a2_b_1,
                              a2_b_2,
                              a2_b_3,
                              ca2,
                              z3_b_0,
                              z3_b_1,
                              z3_b_2,
                              z3_b_3,
                              c2,
                              c2_b_0,
                              c2_b_1,
                              c2_b_2,
                              c2_b_3,
                              cc2,
                              z4_b_0,
                              z4_b_1,
                              z4_b_2,
                              z4_b_3,
                              z4_lo7,
                              z4_msb,
                              b2,
                              mx_b_0,
                              mx_b_1,
                              mx_b_2,
                              mx_b_3,
                              my_b_0,
                              my_b_1,
                              my_b_2,
                              my_b_3,
                              junk_mx01,
                              junk_mx23,
                              junk_my01,
                              junk_my23,
                              blk_len,
                              blk_flags,
                              chain,
                              is_final,    // this compression ends its call, so its digest is claimed
                              claim_value, // the claim this row contributes, scattered through the trace
                              q_claim,     // 1 on rows that contribute a claim
                              link_value,  // the same claims, gathered densely at rows 0..N-1
                              q_link,      // 1 on rows [0, N)
                              counts_0,
                              counts_1,
                              counts_2,
                              counts_3,
                              counts_4,
                              counts_5)
    };

    /** @brief Logup inverse columns, committed after the beta/gamma challenges. */
    template <typename DataType> class DerivedWitnessEntities {
      public:
        DEFINE_FLAVOR_MEMBERS(DataType, inv_0, inv_1, inv_2, inv_3, inv_4, inv_5, inv_link)
    };

    template <typename DataType>
    class WitnessEntities : public WireToBeShiftedEntities<DataType>,
                            public WireNonShiftedEntities<DataType>,
                            public DerivedWitnessEntities<DataType> {
      public:
        DEFINE_COMPOUND_GET_ALL(WireToBeShiftedEntities<DataType>,
                                WireNonShiftedEntities<DataType>,
                                DerivedWitnessEntities<DataType>)

        auto get_wires()
        {
            return concatenate(WireToBeShiftedEntities<DataType>::get_all(),
                               WireNonShiftedEntities<DataType>::get_all());
        }
        static std::vector<std::string> get_wire_labels()
        {
            std::vector<std::string> labels = WireToBeShiftedEntities<DataType>::get_labels();
            const auto& non_shifted = WireNonShiftedEntities<DataType>::get_labels();
            labels.insert(labels.end(), non_shifted.begin(), non_shifted.end());
            return labels;
        }
        auto get_inverses() { return DerivedWitnessEntities<DataType>::get_all(); }
    };

    template <typename DataType> class ShiftedEntities {
      public:
        DEFINE_FLAVOR_MEMBERS(DataType,
                              v_0_shift,
                              v_1_shift,
                              v_2_shift,
                              v_3_shift,
                              v_4_shift,
                              v_5_shift,
                              v_6_shift,
                              v_7_shift,
                              v_8_shift,
                              v_9_shift,
                              v_10_shift,
                              v_11_shift,
                              v_12_shift,
                              v_13_shift,
                              v_14_shift,
                              v_15_shift,
                              m_0_shift,
                              m_1_shift,
                              m_2_shift,
                              m_3_shift,
                              m_4_shift,
                              m_5_shift,
                              m_6_shift,
                              m_7_shift,
                              m_8_shift,
                              m_9_shift,
                              m_10_shift,
                              m_11_shift,
                              m_12_shift,
                              m_13_shift,
                              m_14_shift,
                              m_15_shift,
                              out_0_shift,
                              out_1_shift,
                              out_2_shift,
                              out_3_shift,
                              out_4_shift,
                              out_5_shift,
                              out_6_shift,
                              out_7_shift,
                              claim_index_shift)
    };

    template <typename DataType>
    class AllEntities : public PrecomputedEntities<DataType>,
                        public WitnessEntities<DataType>,
                        public ShiftedEntities<DataType> {
      public:
        DEFINE_COMPOUND_GET_ALL(PrecomputedEntities<DataType>, WitnessEntities<DataType>, ShiftedEntities<DataType>)

        auto get_unshifted()
        {
            return concatenate(PrecomputedEntities<DataType>::get_all(), WitnessEntities<DataType>::get_all());
        }
        auto get_to_be_shifted() { return WireToBeShiftedEntities<DataType>::get_all(); }
        auto get_shifted() { return ShiftedEntities<DataType>::get_all(); }
        auto get_precomputed() { return PrecomputedEntities<DataType>::get_all(); }
    };

    static constexpr size_t NUM_PRECOMPUTED_ENTITIES = PrecomputedEntities<FF>::_members_size;
    static constexpr size_t NUM_WITNESS_ENTITIES = WireToBeShiftedEntities<FF>::_members_size +
                                                   WireNonShiftedEntities<FF>::_members_size +
                                                   DerivedWitnessEntities<FF>::_members_size;
    static constexpr size_t NUM_SHIFTED_ENTITIES = ShiftedEntities<FF>::_members_size;
    static constexpr size_t NUM_TO_BE_SHIFTED = WireToBeShiftedEntities<FF>::_members_size;
    static constexpr size_t NUM_ALL_ENTITIES = NUM_PRECOMPUTED_ENTITIES + WireToBeShiftedEntities<FF>::_members_size +
                                               WireNonShiftedEntities<FF>::_members_size +
                                               DerivedWitnessEntities<FF>::_members_size + NUM_SHIFTED_ENTITIES;
    static constexpr size_t NUM_WIRES_TO_COMMIT =
        WireToBeShiftedEntities<FF>::_members_size + WireNonShiftedEntities<FF>::_members_size;
    static constexpr size_t NUM_INVERSES = DerivedWitnessEntities<FF>::_members_size;

    class AllValues : public AllEntities<FF> {
      public:
        using Base = AllEntities<FF>;
        using Base::Base;
    };

    template <size_t LENGTH> using ProverUnivariates = AllEntities<bb::Univariate<FF, LENGTH>>;
    using ExtendedEdges = ProverUnivariates<MAX_PARTIAL_RELATION_LENGTH>;

    class ProverPolynomials : public AllEntities<Polynomial> {
      public:
        ProverPolynomials() = default;
        ProverPolynomials& operator=(const ProverPolynomials&) = delete;
        ProverPolynomials(const ProverPolynomials& o) = delete;
        ProverPolynomials(ProverPolynomials&& o) noexcept = default;
        ProverPolynomials& operator=(ProverPolynomials&& o) noexcept = default;
        ~ProverPolynomials() = default;

        [[nodiscard]] size_t get_polynomial_size() const { return this->q_table.virtual_size(); }

        AllValues get_row(const size_t row_idx) const
        {
            AllValues result;
            for (auto [result_field, polynomial] : zip_view(result.get_all(), this->get_all())) {
                result_field = polynomial.get(row_idx);
            }
            return result;
        }

        void set_shifted()
        {
            for (auto [shifted, to_be_shifted] : zip_view(get_shifted(), get_to_be_shifted())) {
                shifted = to_be_shifted.shifted();
            }
        }

        ProverPolynomials(const CircuitBuilder& builder)
        {
            using L = Layout;
            const size_t active_rows = builder.num_active_rows();
            const size_t min_size = std::max(active_rows, L::TABLE_ROWS);
            const size_t dyadic_size = 1UL << numeric::get_msb(2 * min_size - 1);
            const size_t num_compressions = builder.num_compressions();

            // Allocation: witness values live on [0, active_rows); the table and read counts on
            // [0, TABLE_ROWS); inverses on the union. Everything beyond is virtual zero.
            for (auto& poly : get_to_be_shifted()) {
                poly = Polynomial::shiftable(active_rows, dyadic_size);
            }
            // The running claim count is defined on every row, not just the active ones, so that its
            // increment relation holds past the end of the trace.
            this->claim_index = Polynomial::shiftable(dyadic_size, dyadic_size);
            for (auto& poly : WireNonShiftedEntities<Polynomial>::get_all()) {
                poly = Polynomial(active_rows, dyadic_size);
            }
            for (auto& poly : RefArray{
                     this->counts_0, this->counts_1, this->counts_2, this->counts_3, this->counts_4, this->counts_5 }) {
                poly = Polynomial(L::TABLE_ROWS, dyadic_size);
            }
            for (auto& poly : DerivedWitnessEntities<Polynomial>::get_all()) {
                poly = Polynomial(min_size, dyadic_size);
            }
            for (auto& poly : RefArray{ this->q_pos_0,
                                        this->q_pos_1,
                                        this->q_pos_2,
                                        this->q_pos_3,
                                        this->q_pos_4,
                                        this->q_pos_5,
                                        this->q_pos_6,
                                        this->q_pos_7,
                                        this->q_out_0,
                                        this->q_out_1,
                                        this->q_out_2,
                                        this->q_out_3,
                                        this->q_mperm }) {
                poly = Polynomial(active_rows, dyadic_size);
            }
            for (auto& poly : RefArray{ this->q_table, this->table_x, this->table_y, this->table_z }) {
                poly = Polynomial(L::TABLE_ROWS, dyadic_size);
            }
            this->row_index = Polynomial(dyadic_size);
            this->q_round0 = Polynomial(active_rows, dyadic_size);
            // The dense link section lives at rows [NUM_DISABLED_ROWS_IN_SUMCHECK,
            // NUM_DISABLED_ROWS_IN_SUMCHECK + num_claims): the offset bb's databus polynomials use,
            // so a consuming Mega circuit's calldata commitment can be compared with link_value's.
            const size_t num_claims = builder.claims.size();
            const size_t link_offset = NUM_DISABLED_ROWS_IN_SUMCHECK;
            BB_ASSERT_LTE(link_offset + num_claims, dyadic_size, "the link section must fit the trace");
            this->link_value = Polynomial(num_claims, dyadic_size, link_offset);
            this->q_link = Polynomial(num_claims, dyadic_size, link_offset);
            this->lagrange_first = Polynomial(1, dyadic_size);
            this->lagrange_first.at(0) = 1;
            this->lagrange_last = Polynomial(1, dyadic_size, dyadic_size - 1);
            this->lagrange_last.at(dyadic_size - 1) = 1;

            for (size_t i = 0; i < dyadic_size; ++i) {
                this->row_index.at(i) = i;
            }
            for (size_t i = 0; i < num_claims; ++i) {
                this->link_value.at(link_offset + i) = builder.claims[i];
                this->q_link.at(link_offset + i) = 1;
            }

            // XOR table and read counts.
            for (size_t i = 0; i < L::TABLE_ROWS; ++i) {
                const uint32_t x = static_cast<uint32_t>(i >> 8);
                const uint32_t y = static_cast<uint32_t>(i & 0xFF);
                this->q_table.at(i) = 1;
                this->table_x.at(i) = x;
                this->table_y.at(i) = y;
                this->table_z.at(i) = x ^ y;
            }
            auto counts_refs = RefArray{ this->counts_0, this->counts_1, this->counts_2,
                                         this->counts_3, this->counts_4, this->counts_5 };
            for (size_t set = 0; set < L::NUM_LOOKUP_SETS; ++set) {
                for (size_t i = 0; i < L::TABLE_ROWS; ++i) {
                    if (builder.read_counts[set][i] != 0) {
                        counts_refs[set].at(i) = builder.read_counts[set][i];
                    }
                }
            }

            // Selectors.
            for (size_t c = 0; c <= num_compressions; ++c) {
                this->q_out_3.at(L::boundary_row_of(c)) = 1;
            }
            RefArray q_pos{ this->q_pos_0, this->q_pos_1, this->q_pos_2, this->q_pos_3,
                            this->q_pos_4, this->q_pos_5, this->q_pos_6, this->q_pos_7 };
            RefArray q_out{ this->q_out_0, this->q_out_1, this->q_out_2, this->q_out_3 };
            for (size_t c = 0; c < num_compressions; ++c) {
                for (size_t r = 0; r < L::NUM_ROUNDS; ++r) {
                    for (size_t p = 0; p < L::NUM_G_PER_ROUND; ++p) {
                        q_pos[p].at(L::row_of_g(c, r, p)) = 1;
                    }
                    if (r + 1 < L::NUM_ROUNDS) {
                        this->q_mperm.at(L::row_of_g(c, r, 7)) = 1;
                    }
                }
                for (size_t t = 0; t + 1 < L::NUM_OUT_ROWS; ++t) {
                    q_out[t].at(L::row_of_out(c, t)) = 1;
                }
                for (size_t p = 0; p < L::NUM_G_PER_ROUND; ++p) {
                    this->q_round0.at(L::row_of_g(c, 0, p)) = 1;
                }
                // q_out_3 rows were set above as boundary rows (row_of_out(c, 3) == boundary_row_of(c + 1)).
            }

            // G-row witness values.
            RefArray v_cols{ this->v_0,  this->v_1,  this->v_2,  this->v_3, this->v_4,  this->v_5,
                             this->v_6,  this->v_7,  this->v_8,  this->v_9, this->v_10, this->v_11,
                             this->v_12, this->v_13, this->v_14, this->v_15 };
            RefArray m_cols{ this->m_0,  this->m_1,  this->m_2,  this->m_3, this->m_4,  this->m_5,
                             this->m_6,  this->m_7,  this->m_8,  this->m_9, this->m_10, this->m_11,
                             this->m_12, this->m_13, this->m_14, this->m_15 };
            RefArray out_cols{ this->out_0, this->out_1, this->out_2, this->out_3,
                               this->out_4, this->out_5, this->out_6, this->out_7 };

            const auto set_bytes = [](RefArray<Polynomial, 4> cols, size_t row, const std::array<uint8_t, 4>& bytes) {
                for (size_t i = 0; i < 4; ++i) {
                    cols[i].at(row) = bytes[i];
                }
            };

            for (size_t c = 0; c < num_compressions; ++c) {
                for (size_t r = 0; r < L::NUM_ROUNDS; ++r) {
                    for (size_t p = 0; p < L::NUM_G_PER_ROUND; ++p) {
                        const size_t row = L::row_of_g(c, r, p);
                        const auto& g = builder.g_rows[c * L::NUM_G_ROWS + r * L::NUM_G_PER_ROUND + p];
                        for (size_t j = 0; j < 16; ++j) {
                            v_cols[j].at(row) = g.v[j];
                            m_cols[j].at(row) = g.m[j];
                        }
                        set_bytes({ this->din_b_0, this->din_b_1, this->din_b_2, this->din_b_3 }, row, g.din_b);
                        set_bytes({ this->bin_b_0, this->bin_b_1, this->bin_b_2, this->bin_b_3 }, row, g.bin_b);
                        set_bytes({ this->a1_b_0, this->a1_b_1, this->a1_b_2, this->a1_b_3 }, row, g.a1_b);
                        set_bytes({ this->z1_b_0, this->z1_b_1, this->z1_b_2, this->z1_b_3 }, row, g.z1_b);
                        set_bytes({ this->c1_b_0, this->c1_b_1, this->c1_b_2, this->c1_b_3 }, row, g.c1_b);
                        set_bytes({ this->z2_b_0, this->z2_b_1, this->z2_b_2, this->z2_b_3 }, row, g.z2_b);
                        set_bytes({ this->b1_b_0, this->b1_b_1, this->b1_b_2, this->b1_b_3 }, row, g.b1_b);
                        set_bytes({ this->a2_b_0, this->a2_b_1, this->a2_b_2, this->a2_b_3 }, row, g.a2_b);
                        set_bytes({ this->z3_b_0, this->z3_b_1, this->z3_b_2, this->z3_b_3 }, row, g.z3_b);
                        set_bytes({ this->c2_b_0, this->c2_b_1, this->c2_b_2, this->c2_b_3 }, row, g.c2_b);
                        set_bytes({ this->z4_b_0, this->z4_b_1, this->z4_b_2, this->z4_b_3 }, row, g.z4_b);
                        set_bytes({ this->mx_b_0, this->mx_b_1, this->mx_b_2, this->mx_b_3 }, row, g.mx_b);
                        set_bytes({ this->my_b_0, this->my_b_1, this->my_b_2, this->my_b_3 }, row, g.my_b);
                        this->a1.at(row) = g.a1;
                        this->c1.at(row) = g.c1;
                        this->b1.at(row) = g.b1;
                        this->a2.at(row) = g.a2;
                        this->c2.at(row) = g.c2;
                        this->b2.at(row) = g.b2;
                        this->ca1.at(row) = g.ca1;
                        this->cc1.at(row) = g.cc1;
                        this->ca2.at(row) = g.ca2;
                        this->cc2.at(row) = g.cc2;
                        this->z2_nib_lo.at(row) = g.z2_nib_lo;
                        this->z2_nib_hi.at(row) = g.z2_nib_hi;
                        this->z4_lo7.at(row) = g.z4_lo7;
                        this->z4_msb.at(row) = g.z4_msb;
                        this->junk_mx01.at(row) = g.mx_b[0] ^ g.mx_b[1];
                        this->junk_mx23.at(row) = g.mx_b[2] ^ g.mx_b[3];
                        this->junk_my01.at(row) = g.my_b[0] ^ g.my_b[1];
                        this->junk_my23.at(row) = g.my_b[2] ^ g.my_b[3];
                    }
                }

                // Output rows: final state, output words, and the reused set-0/1 byte columns.
                const auto& blockdata = builder.out_blocks[c];
                for (size_t t = 0; t < L::NUM_OUT_ROWS; ++t) {
                    const size_t row = L::row_of_out(c, t);
                    for (size_t j = 0; j < 16; ++j) {
                        v_cols[j].at(row) = blockdata.v[j];
                    }
                    for (size_t j = 0; j < 8; ++j) {
                        out_cols[j].at(row) = blockdata.out[j];
                    }
                    const auto bytes_of = [](uint32_t w) {
                        return std::array<uint8_t, 4>{ static_cast<uint8_t>(w),
                                                       static_cast<uint8_t>(w >> 8),
                                                       static_cast<uint8_t>(w >> 16),
                                                       static_cast<uint8_t>(w >> 24) };
                    };
                    set_bytes({ this->din_b_0, this->din_b_1, this->din_b_2, this->din_b_3 },
                              row,
                              bytes_of(blockdata.v[2 * t]));
                    set_bytes({ this->bin_b_0, this->bin_b_1, this->bin_b_2, this->bin_b_3 },
                              row,
                              bytes_of(blockdata.v[2 * t + 1]));
                    set_bytes({ this->a1_b_0, this->a1_b_1, this->a1_b_2, this->a1_b_3 },
                              row,
                              bytes_of(blockdata.v[2 * t + 8]));
                    set_bytes({ this->c1_b_0, this->c1_b_1, this->c1_b_2, this->c1_b_3 },
                              row,
                              bytes_of(blockdata.v[2 * t + 9]));
                    set_bytes({ this->z1_b_0, this->z1_b_1, this->z1_b_2, this->z1_b_3 },
                              row,
                              bytes_of(blockdata.out[2 * t]));
                    set_bytes({ this->z2_b_0, this->z2_b_1, this->z2_b_2, this->z2_b_3 },
                              row,
                              bytes_of(blockdata.out[2 * t + 1]));
                }
            }

            // Claims, in trace order: each compression's eight message chunks, then the four digest
            // chunks of the compressions that end a call.
            {
                size_t claim_index = 0;
                for (size_t c = 0; c < num_compressions; ++c) {
                    for (size_t p = 0; p < L::NUM_G_PER_ROUND; ++p) {
                        const size_t row = L::row_of_g(c, 0, p);
                        const auto& g = builder.g_rows[c * L::NUM_G_ROWS + p];
                        this->q_claim.at(row) = 1;
                        this->claim_value.at(row) =
                            fr(static_cast<uint64_t>(g.m[2 * p]) + (static_cast<uint64_t>(g.m[2 * p + 1]) << 32));
                        this->claim_index.at(row) = claim_index++;
                    }
                    const bool final_compression = builder.is_final[c];
                    for (size_t t = 0; t < L::NUM_OUT_ROWS; ++t) {
                        const size_t row = L::row_of_out(c, t);
                        this->is_final.at(row) = final_compression ? 1 : 0;
                        if (final_compression) {
                            const auto& block = builder.out_blocks[c];
                            this->q_claim.at(row) = 1;
                            this->claim_value.at(row) = fr(static_cast<uint64_t>(block.out[2 * t]) +
                                                           (static_cast<uint64_t>(block.out[2 * t + 1]) << 32));
                            this->claim_index.at(row) = claim_index++;
                        }
                    }
                }
                BB_ASSERT_EQ(claim_index, num_claims, "the trace's claims must match the builder's claim vector");
                // claim_index holds the running count on every row, not just claiming ones.
                size_t running = 0;
                for (size_t row = 1; row < dyadic_size; ++row) {
                    running += (this->q_claim.get(row - 1) == fr(1)) ? 1 : 0;
                    this->claim_index.at(row) = running;
                }
            }

            // Boundary-row block parameters (the genesis row carries compression 0's).
            for (size_t c = 0; c < num_compressions; ++c) {
                const size_t row = L::boundary_row_of(c);
                this->blk_len.at(row) = builder.boundaries[c].block_len;
                this->blk_flags.at(row) = builder.boundaries[c].flags;
                this->chain.at(row) = builder.boundaries[c].chain ? 1 : 0;
            }

            // Ghost row: the last boundary row constrains the next row's state to a fresh IV state
            // with zero block parameters.
            const size_t ghost_row = L::boundary_row_of(num_compressions) + 1;
            for (size_t j = 0; j < 8; ++j) {
                v_cols[j].at(ghost_row) = blake3::IV[j];
            }
            for (size_t j = 8; j < 12; ++j) {
                v_cols[j].at(ghost_row) = blake3::IV[j - 8];
            }

            set_shifted();
        }
    };

    using PartiallyEvaluatedMultivariates =
        PartiallyEvaluatedMultivariatesBase<AllEntities<Polynomial>, ProverPolynomials, Polynomial>;

    class ProvingKey {
      public:
        size_t circuit_size;
        size_t log_circuit_size;
        ProverPolynomials polynomials;
        CommitmentKey commitment_key;

        ProvingKey(const CircuitBuilder& builder)
            : polynomials(builder)
        {
            circuit_size = polynomials.get_polynomial_size();
            log_circuit_size = numeric::get_msb(circuit_size);
            commitment_key = CommitmentKey(circuit_size);
        }
    };

    /** @brief Commitments to the precomputed polynomials plus the circuit size. */
    class VerificationKey : public PrecomputedEntities<Commitment> {
      public:
        size_t circuit_size = 0;
        size_t log_circuit_size = 0;

        VerificationKey() = default;
        VerificationKey(ProvingKey& proving_key)
            : circuit_size(proving_key.circuit_size)
            , log_circuit_size(proving_key.log_circuit_size)
        {
            for (auto [commitment, poly] : zip_view(this->get_all(), proving_key.polynomials.get_precomputed())) {
                commitment = proving_key.commitment_key.commit(poly);
            }
        }
    };

    class VerifierCommitments : public AllEntities<Commitment> {
      public:
        VerifierCommitments(const std::shared_ptr<VerificationKey>& verification_key)
        {
            for (auto [commitment, vk_commitment] : zip_view(this->get_precomputed(), verification_key->get_all())) {
                commitment = vk_commitment;
            }
        }
    };
};

} // namespace bb
