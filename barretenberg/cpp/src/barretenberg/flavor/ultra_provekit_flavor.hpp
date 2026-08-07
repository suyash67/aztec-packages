#pragma once

#include "barretenberg/commitment_schemes/kzg/kzg.hpp"
#include "barretenberg/common/throw_or_abort.hpp"
#include "barretenberg/flavor/flavor.hpp"
#include "barretenberg/flavor/flavor_concepts.hpp"
#include "barretenberg/flavor/generated/ultra_provekit_flavor_generated.hpp"
#include "barretenberg/flavor/partially_evaluated_multivariates.hpp"
#include "barretenberg/flavor/prover_polynomials.hpp"
#include "barretenberg/polynomials/univariate.hpp"
#include "barretenberg/relations/relation_tuple_helpers.hpp"
#include "barretenberg/stdlib_circuit_builders/ultra_circuit_builder.hpp"
#include "barretenberg/transcript/transcript.hpp"

namespace bb {

/**
 * @brief UltraFlavor restricted to the gate kinds a ProveKit-style Noir circuit actually uses.
 *
 * @details Same arithmetization, same builder (`UltraCircuitBuilder`) and same trace layout as
 * `UltraFlavor`; only the entity list and the sumcheck relation tuple shrink. The elliptic,
 * non-native-field and both Poseidon2 relations are dropped, along with their four gate selectors
 * (`q_elliptic`, `q_nnf`, `q_poseidon2_external`, `q_poseidon2_internal`). This is sound exactly
 * when the corresponding trace blocks are empty: an empty block contributes no rows, so its
 * selector is identically zero and its subrelations vanish on every row of the trace. Circuits
 * that do in-circuit curve arithmetic, bignum arithmetic or Poseidon2 hashing must keep using
 * `UltraFlavor` — `assert_relations_are_dead()` checks the precondition on a finalized builder.
 *
 * Two savings follow. The precomputed entity count drops 28 → 24, so four fewer columns are
 * committed and opened. And because the Poseidon2 subrelations are the only degree-6 ones in
 * Ultra, `MAX_PARTIAL_RELATION_LENGTH` drops 7 → 6: sumcheck round polynomials carry one fewer
 * evaluation and the per-edge extension in the sumcheck hot loop is one degree shorter.
 *
 * There is no ZK variant and no Solidity verifier for this flavor; it exists to drive
 * `TransparentHonk` over hash-based PCS backends.
 */
class UltraProveKitFlavor : public UltraProveKitFlavor_Generated {
  public:
    using CircuitBuilder = UltraCircuitBuilder;
    using Curve = curve::BN254;
    using FF = Curve::ScalarField;
    using GroupElement = Curve::Element;
    using Commitment = Curve::AffineElement;
    using PCS = KZG<Curve>;
    using Polynomial = bb::Polynomial<FF>;
    using CommitmentKey = bb::CommitmentKey<Curve>;
    using Codec = FrCodec;
    using HashFunction = crypto::Poseidon2<crypto::Poseidon2Bn254ScalarFieldParams>;
    using Transcript = BaseTranscript<Codec, HashFunction>;

    static constexpr size_t VIRTUAL_LOG_N = CONST_PROOF_SIZE_LOG_N;
    static constexpr bool USE_SHORT_MONOMIALS = true;
    static constexpr bool USE_SIMD_SUMCHECK = true;
    static constexpr bool HasZK = false;
    // TransparentHonk runs sumcheck unpadded (`virtual_log_n = log_n`); nothing here needs the
    // constant-proof-size padding that the recursive UltraHonk verifier requires.
    static constexpr bool USE_PADDING = false;
    static constexpr size_t NUM_WIRES = CircuitBuilder::NUM_WIRES;

    using Generated = UltraProveKitFlavor_Generated;

    using Relations = Relations_<FF>;

    static constexpr size_t MAX_PARTIAL_RELATION_LENGTH = compute_max_partial_relation_length<Relations>();
    static_assert(MAX_PARTIAL_RELATION_LENGTH == 6,
                  "dropping the Poseidon2 relations should leave 6 as the longest subrelation");
    static constexpr size_t NUM_SUBRELATIONS = compute_number_of_subrelations<Relations>();
    using SubrelationSeparator = FF;

    static constexpr size_t BATCHED_RELATION_PARTIAL_LENGTH = MAX_PARTIAL_RELATION_LENGTH + 1;
    static constexpr size_t NUM_RELATIONS = std::tuple_size_v<Relations>;

    static constexpr size_t TRACE_OFFSET = NUM_DISABLED_ROWS_IN_SUMCHECK;

    static constexpr size_t FINAL_PCS_MSM_SIZE(size_t log_n = VIRTUAL_LOG_N)
    {
        return NUM_UNSHIFTED_ENTITIES + log_n + 2;
    }

    class AllValues : public AllEntities<FF> {
      public:
        using Base = AllEntities<FF>;
        using Base::Base;
    };

    static_assert(gemini_masking_layout_consistent<UltraProveKitFlavor>(),
                  "UltraProveKitFlavor gemini masking flag must match its entity layout");

    using ProverPolynomials = ProverPolynomialsBase<AllEntities<Polynomial>, AllValues, Polynomial>;

    using PrecomputedData = PrecomputedData_<Polynomial, NUM_PRECOMPUTED_ENTITIES>;

    using VerificationKey = NativeVerificationKey_<PrecomputedEntities<Commitment>, Codec, HashFunction, CommitmentKey>;

    using VKAndHash = VKAndHash_<FF, VerificationKey>;

    using PartiallyEvaluatedMultivariates =
        PartiallyEvaluatedMultivariatesBase<AllEntities<Polynomial>, ProverPolynomials, Polynomial>;

    template <size_t LENGTH> using ProverUnivariates = AllEntities<bb::Univariate<FF, LENGTH>>;

    using ExtendedEdges = ProverUnivariates<MAX_PARTIAL_RELATION_LENGTH>;

    using WitnessCommitments = WitnessEntities<Commitment>;

    using CommitmentLabels = AllEntities<std::string>;
    static const CommitmentLabels& commitment_labels()
    {
        static const CommitmentLabels instance = []() {
            CommitmentLabels result;
            const auto& src = AllEntities<std::string>::get_labels();
            std::copy(src.begin(), src.end(), result.data.begin());
            return result;
        }();
        return instance;
    }

    /**
     * @brief Aborts unless every trace block this flavor drops is empty.
     * @details The dropped relations are only vacuous when their gate blocks contribute no rows.
     * Call on a finalized builder, before handing it to a prover.
     */
    static void assert_relations_are_dead(const CircuitBuilder& circuit)
    {
        // Unconditional (not BB_ASSERT): a release build that skipped this would prove a statement
        // whose dropped relations are not vacuous, and the proof would still verify.
        const auto require_empty = [](size_t size, const char* gate_kind) {
            if (size != 0) {
                throw_or_abort(std::string("UltraProveKitFlavor cannot prove a circuit with ") + gate_kind + " gates");
            }
        };
        require_empty(circuit.blocks.elliptic.size(), "elliptic");
        require_empty(circuit.blocks.nnf.size(), "non-native field");
        require_empty(circuit.blocks.poseidon2_external.size(), "poseidon2 external");
        require_empty(circuit.blocks.poseidon2_internal.size(), "poseidon2 internal");
    }
};

} // namespace bb
