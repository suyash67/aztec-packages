#pragma once

#include "barretenberg/commitment_schemes/commitment_key.hpp"
#include "barretenberg/ecc/curves/pasta/pasta.hpp"
#include "barretenberg/flavor/flavor.hpp"
#include "barretenberg/flavor/flavor_macros.hpp"
#include "barretenberg/flavor/partially_evaluated_multivariates.hpp"
#include "barretenberg/polynomials/polynomial.hpp"
#include "barretenberg/polynomials/univariate.hpp"
#include "barretenberg/relations/relation_parameters.hpp"
#include "barretenberg/relations/relation_tuple_helpers.hpp"
#include "barretenberg/zcash/halo2/layout.hpp"
#include "barretenberg/zcash/honk/orchard_relations.hpp"
#include "barretenberg/zcash/honk/pasta_crs.hpp"
#include "barretenberg/zcash/honk/pasta_transcript.hpp"
#include "barretenberg/zcash/primitives/cycle.hpp"

// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)

namespace bb::zcash {

/**
 * @brief Honk flavor for the halo2 arithmetization of the Orchard Action circuit over the Pasta cycle.
 *
 * @details The execution trace is the halo2 table: 10 advice columns, 11 fixed columns, 45 selectors, the 3-column
 * Sinsemilla table, 14 equality-enabled columns. The relations are the halo2 gates (OrchardGateRelation), the copy
 * constraints (OrchardPermutationRelation) and the two lookup arguments. Gates that use halo2's Rotation::prev() are
 * anchored one row earlier and read advice cells shifted by two, so the advice columns are opened at their shifts by
 * one and by two.
 *
 * Zero knowledge: sumcheck is masked by Libra, the PCS by a Gemini masking polynomial, and every witness column holds
 * random values in rows 2..7. Rows 0..7 are disabled in sumcheck (TRACE_OFFSET = 8); rows 0 and 1 stay zero so that
 * the witness polynomials are divisible by X^2. The final opening is the halo2 IPA on Vesta.
 */
class OrchardFlavor {
  public:
    using Cycle = PastaCycle;
    using Curve = curve::Vesta;
    using FF = Curve::ScalarField;
    using BF = Curve::BaseField;
    using G1 = Curve::Group;
    using GroupElement = Curve::Element;
    using Commitment = Curve::AffineElement;
    using Polynomial = bb::Polynomial<FF>;
    using CommitmentKey = bb::CommitmentKey<Curve>;
    using Codec = PastaCodec;
    using Transcript = PastaTranscript;
    using Proof = std::vector<uint256_t>;
    static_assert(std::is_same_v<FF, Cycle::FF>);

    static constexpr bool HasZK = true;
    static constexpr bool HasGeminiMasking = true;
    static constexpr bool USE_SHORT_MONOMIALS = false;
    static constexpr bool USE_PADDING = false;
    static constexpr size_t LOG_NUM_DISABLED_ROWS = 3;
    static constexpr size_t TRACE_OFFSET = size_t{ 1 } << LOG_NUM_DISABLED_ROWS;
    // Leading zero rows of the witness polynomials (they are shifted by up to 2).
    static constexpr size_t NUM_ZERO_ROWS = 2;
    static constexpr size_t NUM_PERMUTATION_COLUMNS = halo2::NUM_PERMUTATION_COLUMNS;

    template <typename FF_>
    using Relations_ = std::tuple<OrchardGateRelation<OrchardFlavor, FF_>,
                                  OrchardPermutationRelation<OrchardFlavor, FF_>,
                                  OrchardSinsemillaLookupRelation<OrchardFlavor, FF_>,
                                  OrchardRangeLookupRelation<OrchardFlavor, FF_>>;
    using Relations = Relations_<FF>;

    static constexpr size_t NUM_SUBRELATIONS = compute_number_of_subrelations<Relations>();
    using SubrelationSeparators = std::array<FF, NUM_SUBRELATIONS - 1>;
    static constexpr size_t MAX_PARTIAL_RELATION_LENGTH = compute_max_partial_relation_length<Relations>();
    // +1 for the gate separator, +1 for the row-disabling polynomial.
    static constexpr size_t BATCHED_RELATION_PARTIAL_LENGTH = MAX_PARTIAL_RELATION_LENGTH + 2;
    static constexpr size_t NUM_RELATIONS = std::tuple_size_v<Relations>;
    static_assert(Curve::LIBRA_UNIVARIATES_LENGTH <= BATCHED_RELATION_PARTIAL_LENGTH);

    template <typename DataType> class FixedEntities {
      public:
        DEFINE_FLAVOR_MEMBERS(
            DataType, f_0, f_1, f_2, f_3, f_4, f_5, f_6, f_7, fixed_z, q_sinsemilla2_1, q_sinsemilla2_2)
    };
    template <typename DataType> class SelectorEntities {
      public:
        DEFINE_FLAVOR_MEMBERS(DataType,
                              q_orchard,
                              q_add_field,
                              q_point,
                              q_point_non_id,
                              q_add_incomplete,
                              q_add,
                              q_mul_hi_1,
                              q_mul_hi_2,
                              q_mul_hi_3,
                              q_mul_lo_1,
                              q_mul_lo_2,
                              q_mul_lo_3,
                              q_mul_decompose_var,
                              q_mul_lsb,
                              q_mul_overflow,
                              q_mul_fixed_running_sum,
                              q_mul_fixed_full,
                              q_mul_fixed_short,
                              q_mul_fixed_base_field,
                              q_poseidon_full,
                              q_poseidon_partial,
                              q_poseidon_pad_and_add,
                              q_sinsemilla1_1,
                              q_sinsemilla4_1,
                              q_sinsemilla1_2,
                              q_sinsemilla4_2,
                              q_merkle_decompose_1,
                              q_merkle_decompose_2,
                              q_swap_1,
                              q_swap_2,
                              q_lookup,
                              q_running,
                              q_bitshift,
                              q_commit_ivk,
                              q_notecommit_b,
                              q_notecommit_d,
                              q_notecommit_e,
                              q_notecommit_g,
                              q_notecommit_h,
                              q_notecommit_g_d,
                              q_notecommit_pk_d,
                              q_notecommit_value,
                              q_notecommit_rho,
                              q_notecommit_psi,
                              q_y_canon)
    };
    template <typename DataType> class TableEntities {
      public:
        DEFINE_FLAVOR_MEMBERS(DataType, table_idx, table_x, table_y, q_table)
    };
    template <typename DataType> class SigmaEntities {
      public:
        DEFINE_FLAVOR_MEMBERS(DataType,
                              sigma_0,
                              sigma_1,
                              sigma_2,
                              sigma_3,
                              sigma_4,
                              sigma_5,
                              sigma_6,
                              sigma_7,
                              sigma_8,
                              sigma_9,
                              sigma_10,
                              sigma_11,
                              sigma_12,
                              sigma_13)
    };
    template <typename DataType> class IdEntities {
      public:
        DEFINE_FLAVOR_MEMBERS(
            DataType, id_0, id_1, id_2, id_3, id_4, id_5, id_6, id_7, id_8, id_9, id_10, id_11, id_12, id_13)
    };
    template <typename DataType> class LagrangeEntities {
      public:
        DEFINE_FLAVOR_MEMBERS(DataType, lagrange_first, lagrange_last)
    };

    template <typename DataType>
    class PrecomputedEntities : public FixedEntities<DataType>,
                                public SelectorEntities<DataType>,
                                public TableEntities<DataType>,
                                public SigmaEntities<DataType>,
                                public IdEntities<DataType>,
                                public LagrangeEntities<DataType> {
      public:
        DEFINE_COMPOUND_GET_ALL(FixedEntities<DataType>,
                                SelectorEntities<DataType>,
                                TableEntities<DataType>,
                                SigmaEntities<DataType>,
                                IdEntities<DataType>,
                                LagrangeEntities<DataType>)
    };

    template <typename DataType> class AdviceEntities {
      public:
        DEFINE_FLAVOR_MEMBERS(DataType, a_0, a_1, a_2, a_3, a_4, a_5, a_6, a_7, a_8, a_9)
    };
    template <typename DataType> class DerivedWitnessEntities {
      public:
        DEFINE_FLAVOR_MEMBERS(DataType,
                              lookup_read_counts_sinsemilla,
                              lookup_read_counts_range,
                              lookup_inverses_sinsemilla,
                              lookup_inverses_range,
                              z_perm_mid,
                              z_perm)
    };
    template <typename DataType>
    class WitnessEntities : public AdviceEntities<DataType>, public DerivedWitnessEntities<DataType> {
      public:
        DEFINE_COMPOUND_GET_ALL(AdviceEntities<DataType>, DerivedWitnessEntities<DataType>)
    };
    template <typename DataType> class MaskingEntities {
      public:
        DEFINE_FLAVOR_MEMBERS(DataType, gemini_masking_poly)
    };
    template <typename DataType> class AdviceShiftEntities {
      public:
        DEFINE_FLAVOR_MEMBERS(DataType,
                              a_0_shift,
                              a_1_shift,
                              a_2_shift,
                              a_3_shift,
                              a_4_shift,
                              a_5_shift,
                              a_6_shift,
                              a_7_shift,
                              a_8_shift,
                              a_9_shift)
    };
    template <typename DataType> class ShiftedEntities : public AdviceShiftEntities<DataType> {
      public:
        DataType z_perm_shift;
        [[nodiscard]] auto get_all()
        {
            return concatenate(AdviceShiftEntities<DataType>::get_all(), RefArray{ z_perm_shift });
        }
        [[nodiscard]] auto get_all() const
        {
            return concatenate(AdviceShiftEntities<DataType>::get_all(), RefArray{ z_perm_shift });
        }
        static const std::vector<std::string>& get_labels()
        {
            static const auto labels =
                concatenate(AdviceShiftEntities<DataType>::get_labels(), std::vector<std::string>{ "z_perm_shift" });
            return labels;
        }
        static constexpr size_t size() { return AdviceShiftEntities<DataType>::size() + 1; }
    };
    template <typename DataType> class ShiftedByTwoEntities {
      public:
        DEFINE_FLAVOR_MEMBERS(DataType,
                              a_0_shift2,
                              a_1_shift2,
                              a_2_shift2,
                              a_3_shift2,
                              a_4_shift2,
                              a_5_shift2,
                              a_6_shift2,
                              a_7_shift2,
                              a_8_shift2,
                              a_9_shift2)
    };

    template <typename DataType>
    class AllEntities : public MaskingEntities<DataType>,
                        public PrecomputedEntities<DataType>,
                        public WitnessEntities<DataType>,
                        public ShiftedEntities<DataType>,
                        public ShiftedByTwoEntities<DataType> {
      public:
        DEFINE_COMPOUND_GET_ALL(MaskingEntities<DataType>,
                                PrecomputedEntities<DataType>,
                                WitnessEntities<DataType>,
                                ShiftedEntities<DataType>,
                                ShiftedByTwoEntities<DataType>)
        auto get_unshifted()
        {
            return concatenate(MaskingEntities<DataType>::get_all(),
                               PrecomputedEntities<DataType>::get_all(),
                               WitnessEntities<DataType>::get_all());
        }
        auto get_precomputed() { return PrecomputedEntities<DataType>::get_all(); }
        auto get_witness() { return WitnessEntities<DataType>::get_all(); }
        auto get_witness() const { return WitnessEntities<DataType>::get_all(); }
        auto get_to_be_shifted() { return concatenate(AdviceEntities<DataType>::get_all(), RefArray{ this->z_perm }); }
        auto get_shifted() { return ShiftedEntities<DataType>::get_all(); }
        auto get_to_be_shifted_by_two() { return AdviceEntities<DataType>::get_all(); }
        auto get_shifted_by_two() { return ShiftedByTwoEntities<DataType>::get_all(); }
    };

    static constexpr size_t NUM_PRECOMPUTED_ENTITIES =
        size_t{ halo2::NUM_FIXED } + size_t{ halo2::NUM_SELECTORS } + 4 + (2 * halo2::NUM_PERMUTATION_COLUMNS) + 2;
    static constexpr size_t NUM_WITNESS_ENTITIES = halo2::NUM_ADVICE + 6;
    static constexpr size_t NUM_ALL_ENTITIES =
        1 + NUM_PRECOMPUTED_ENTITIES + NUM_WITNESS_ENTITIES + (halo2::NUM_ADVICE + 1) + halo2::NUM_ADVICE;

    // Accessors used by the relations.
    template <typename E> static auto advice(const E& in)
    {
        return static_cast<const AdviceEntities<std::remove_cvref_t<decltype(in.a_0)>>&>(in).get_all();
    }
    template <typename E> static auto advice_shift(const E& in)
    {
        return static_cast<const AdviceShiftEntities<std::remove_cvref_t<decltype(in.a_0)>>&>(in).get_all();
    }
    template <typename E> static auto advice_shift2(const E& in)
    {
        return static_cast<const ShiftedByTwoEntities<std::remove_cvref_t<decltype(in.a_0)>>&>(in).get_all();
    }
    template <typename E> static auto fixed_columns(const E& in)
    {
        return static_cast<const FixedEntities<std::remove_cvref_t<decltype(in.a_0)>>&>(in).get_all();
    }
    template <typename E> static auto gate_selectors(const E& in)
    {
        return static_cast<const SelectorEntities<std::remove_cvref_t<decltype(in.a_0)>>&>(in).get_all();
    }
    template <typename E> static auto sigmas(const E& in)
    {
        return static_cast<const SigmaEntities<std::remove_cvref_t<decltype(in.a_0)>>&>(in).get_all();
    }
    template <typename E> static auto ids(const E& in)
    {
        return static_cast<const IdEntities<std::remove_cvref_t<decltype(in.a_0)>>&>(in).get_all();
    }
    // The 14 equality-enabled columns in permutation order: a0..a9, f0, f5, f6, f7.
    template <typename E> static auto permutation_columns(const E& in)
    {
        return RefArray{ in.a_0, in.a_1, in.a_2, in.a_3, in.a_4, in.a_5, in.a_6,
                         in.a_7, in.a_8, in.a_9, in.f_0, in.f_5, in.f_6, in.f_7 };
    }

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

        /**
         * @brief Allocates every polynomial with the given virtual size: precomputed polynomials from row 0, the
         * to-be-shifted ones from row NUM_ZERO_ROWS (advice, z_perm) and the other witnesses from row 0.
         */
        explicit ProverPolynomials(size_t circuit_size)
        {
            for (auto& poly : this->get_precomputed()) {
                poly = Polynomial(circuit_size);
            }
            for (auto& poly : this->get_witness()) {
                poly = Polynomial(circuit_size - NUM_ZERO_ROWS, circuit_size, NUM_ZERO_ROWS);
            }
            this->gemini_masking_poly = Polynomial(circuit_size);
        }

        [[nodiscard]] size_t get_polynomial_size() const { return this->lagrange_first.size(); }
        AllValues get_row(size_t row_idx) const
        {
            AllValues result;
            for (auto [result_field, polynomial] : zip_view(result.get_all(), this->get_all())) {
                result_field = polynomial.get(row_idx);
            }
            return result;
        }
        void set_shifted()
        {
            for (auto [shifted, to_be_shifted] : zip_view(this->get_shifted(), this->get_to_be_shifted())) {
                shifted = to_be_shifted.shifted();
            }
            for (auto [shifted, to_be_shifted] :
                 zip_view(this->get_shifted_by_two(), this->get_to_be_shifted_by_two())) {
                shifted = to_be_shifted.shifted().shifted();
            }
        }
    };

    using PartiallyEvaluatedMultivariates =
        PartiallyEvaluatedMultivariatesBase<AllEntities<Polynomial>, ProverPolynomials, Polynomial>;

    class CommitmentLabels : public AllEntities<std::string> {
      public:
        CommitmentLabels()
        {
            for (auto [label, name] : zip_view(this->get_all(), AllEntities<std::string>::get_labels())) {
                label = name;
            }
        }
    };

    /**
     * @brief Verification key: commitments to the precomputed polynomials and the public-input layout.
     * @details For public input i the VK stores the permutation label that the cell holding it would point to in
     * its copy cycle (`public_input_next_labels`); its sigma points to the unique label `special_label(i)` instead,
     * which makes the grand product equal to the public-input delta.
     */
    class VerificationKey : public PrecomputedEntities<Commitment> {
      public:
        size_t log_circuit_size = 0;
        size_t num_public_inputs = 0;
        std::vector<FF> public_input_next_labels;

        size_t circuit_size() const { return size_t{ 1 } << log_circuit_size; }
        FF special_label(size_t i) const { return FF(NUM_PERMUTATION_COLUMNS * circuit_size() + i); }

        FF hash() const
        {
            std::vector<uint256_t> words{ uint256_t(log_circuit_size), uint256_t(num_public_inputs) };
            for (const auto& label : public_input_next_labels) {
                words.push_back(uint256_t(label));
            }
            for (const auto& c : this->get_all()) {
                for (const auto& w : Codec::serialize_to_fields(c)) {
                    words.push_back(w);
                }
            }
            return FF(Blake2bTranscriptHash::hash(words));
        }

        // Public-input delta = prod_i (pi_i + beta * next_label_i + gamma) / (pi_i + beta * special_label_i + gamma).
        FF compute_public_input_delta(std::span<const FF> public_inputs, const FF& beta, const FF& gamma) const
        {
            BB_ASSERT_EQ(public_inputs.size(), num_public_inputs);
            FF num(1);
            FF den(1);
            for (size_t i = 0; i < public_inputs.size(); ++i) {
                num *= public_inputs[i] + beta * public_input_next_labels[i] + gamma;
                den *= public_inputs[i] + beta * special_label(i) + gamma;
            }
            return num / den;
        }
    };

    class VerifierCommitments : public AllEntities<Commitment> {
      public:
        explicit VerifierCommitments(const VerificationKey& vk)
        {
            for (auto [commitment, vk_commitment] : zip_view(this->get_precomputed(), vk.get_all())) {
                commitment = vk_commitment;
            }
        }
    };
};

static_assert(OrchardFlavor::PrecomputedEntities<int>{}.size() == OrchardFlavor::NUM_PRECOMPUTED_ENTITIES);
static_assert(OrchardFlavor::WitnessEntities<int>{}.size() == OrchardFlavor::NUM_WITNESS_ENTITIES);
static_assert(OrchardFlavor::AllEntities<int>{}.size() == OrchardFlavor::NUM_ALL_ENTITIES);

} // namespace bb::zcash

// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)
