#pragma once

#include "barretenberg/commitment_schemes/commitment_key.hpp"
#include "barretenberg/ecc/curves/pasta/pasta.hpp"
#include "barretenberg/flavor/ultra_zk_flavor.hpp"
#include "barretenberg/zcash/honk/pasta_crs.hpp"
#include "barretenberg/zcash/honk/pasta_transcript.hpp"
#include "barretenberg/zcash/ultra_pasta/ultra_pasta_builder.hpp"

/**
 * @file ultra_pasta_flavor.hpp
 * @brief UltraHonk with ZK over the Pallas base field: the UltraZK entity layout and relations, with commitments on
 * Vesta, halo2's BLAKE2b transcript and the halo2 IPA as the PCS.
 */
namespace bb::zcash {

/**
 * @brief The UltraZK flavor over F_p (Pallas base field = Vesta scalar field).
 * @details The entity layout (selectors, wires, lookup and permutation polynomials, the Gemini masking polynomial) is
 * UltraZK's. The relation set is Ultra's without the two Poseidon2 relations, which are defined over BN254 constants;
 * their selectors stay in the layout and are identically zero. Proofs are not padded to a constant size: Sumcheck
 * and Gemini run over the actual log circuit size, as halo2 does.
 */
class UltraPastaZKFlavor : public UltraZKFlavor {
  public:
    using CircuitBuilder = UltraPastaCircuitBuilder;
    using Curve = curve::Vesta;
    using FF = Curve::ScalarField;
    using GroupElement = Curve::Element;
    using Commitment = Curve::AffineElement;
    using Polynomial = bb::Polynomial<FF>;
    using CommitmentKey = bb::CommitmentKey<Curve>;
    using Codec = PastaCodec;
    using HashFunction = Blake2bTranscriptHash;
    using Transcript = PastaTranscript;
    using Proof = std::vector<Codec::DataType>;
    static_assert(std::is_same_v<FF, pasta::fp>);

    static constexpr bool USE_PADDING = false;
    static constexpr bool USE_SIMD_SUMCHECK = false;

    template <typename FF_>
    using Relations_ = std::tuple<bb::UltraPermutationRelation<FF_>,
                                  bb::LogDerivLookupRelation<FF_>,
                                  bb::ArithmeticRelation<FF_>,
                                  bb::DeltaRangeConstraintRelation<FF_>,
                                  bb::EllipticRelation<FF_>,
                                  bb::MemoryRelation<FF_>,
                                  bb::NonNativeFieldRelation<FF_>>;
    using Relations = Relations_<FF>;

    static constexpr size_t MAX_PARTIAL_RELATION_LENGTH = compute_max_partial_relation_length<Relations>();
    static constexpr size_t NUM_SUBRELATIONS = compute_number_of_subrelations<Relations>();
    using SubrelationSeparator = FF;
    static constexpr size_t NUM_RELATIONS = std::tuple_size_v<Relations>;
    // +1 for the gate separator and +1 for the row-disabling polynomial; the Libra univariates (whose length is fixed
    // by the curve) must fit in a round univariate.
    static constexpr size_t BATCHED_RELATION_PARTIAL_LENGTH =
        std::max<size_t>(MAX_PARTIAL_RELATION_LENGTH + 2, Curve::LIBRA_UNIVARIATES_LENGTH);

    class AllValues : public AllEntities<FF> {
      public:
        using Base = AllEntities<FF>;
        using Base::Base;
    };

    using ProverPolynomials = ProverPolynomialsBase<AllEntities<Polynomial>, AllValues, Polynomial>;
    using PrecomputedData = PrecomputedData_<Polynomial, NUM_PRECOMPUTED_ENTITIES>;
    using VerificationKey = NativeVerificationKey_<PrecomputedEntities<Commitment>, Codec, HashFunction, CommitmentKey>;
    using VKAndHash = VKAndHash_<FF, VerificationKey>;
    using PartiallyEvaluatedMultivariates =
        PartiallyEvaluatedMultivariatesBase<AllEntities<Polynomial>, ProverPolynomials, Polynomial>;

    template <size_t LENGTH> using ProverUnivariates = AllEntities<bb::Univariate<FF, LENGTH>>;
    using ExtendedEdges = ProverUnivariates<MAX_PARTIAL_RELATION_LENGTH>;
    using WitnessCommitments = WitnessEntities<Commitment>;
};

} // namespace bb::zcash
