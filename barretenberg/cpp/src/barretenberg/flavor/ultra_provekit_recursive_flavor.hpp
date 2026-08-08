#pragma once

#include "barretenberg/flavor/ultra_provekit_flavor.hpp"
#include "barretenberg/stdlib/primitives/field/field.hpp"
#include "barretenberg/transcript/transcript.hpp"

namespace bb {

/**
 * @brief The recursive counterpart of `UltraProveKitFlavor`, for verifying one of its proofs in a
 * circuit.
 *
 * @details Same shape as `UltraRecursiveFlavor_`: the field type becomes `stdlib::field_t` and the
 * prover-side types are dropped, since a verifier circuit never emulates prover work. The relation
 * set is the flavor's own, re-instantiated over the circuit field, so `SumcheckVerifier<Flavor>`
 * evaluates exactly the subrelations the native verifier does.
 *
 * There is no `Commitment` type worth naming here. `TransparentHonk` commits with a hash-based
 * scheme, so a "commitment" is a Merkle root — an `FF`, not a curve point — which is the reason a
 * recursive verifier for this flavor needs no non-native arithmetic at all.
 */
template <typename BuilderType> class UltraProveKitRecursiveFlavor_ {
  public:
    using CircuitBuilder = BuilderType;
    using FF = stdlib::field_t<CircuitBuilder>;
    using NativeFlavor = UltraProveKitFlavor;
    using NativeVerificationKey = NativeFlavor::VerificationKey;
    using Codec = stdlib::StdlibCodec<FF>;
    using Transcript = StdlibTranscript<CircuitBuilder>;

    // A hash-based commitment is a root, so the "commitment" type is the field itself. Sumcheck
    // only needs the typedef to exist.
    using Commitment = FF;
    using GroupElement = FF;

    static constexpr size_t VIRTUAL_LOG_N = NativeFlavor::VIRTUAL_LOG_N;
    static constexpr bool USE_SHORT_MONOMIALS = NativeFlavor::USE_SHORT_MONOMIALS;
    static constexpr bool USE_PADDING = NativeFlavor::USE_PADDING;
    static constexpr bool HasZK = NativeFlavor::HasZK;
    static constexpr bool HasLogDerivLookup = NativeFlavor::HasLogDerivLookup;
    static constexpr bool HasElliptic = NativeFlavor::HasElliptic;
    static constexpr bool HasMemory = NativeFlavor::HasMemory;
    static constexpr bool HasNonNativeField = NativeFlavor::HasNonNativeField;
    static constexpr bool HasEccOpQueue = NativeFlavor::HasEccOpQueue;
    static constexpr bool HasDataBus = NativeFlavor::HasDataBus;

    static constexpr size_t NUM_WIRES = NativeFlavor::NUM_WIRES;
    static constexpr size_t NUM_ALL_ENTITIES = NativeFlavor::NUM_ALL_ENTITIES;
    static constexpr size_t NUM_PRECOMPUTED_ENTITIES = NativeFlavor::NUM_PRECOMPUTED_ENTITIES;
    static constexpr size_t NUM_WITNESS_ENTITIES = NativeFlavor::NUM_WITNESS_ENTITIES;
    static constexpr size_t NUM_UNSHIFTED_ENTITIES = NativeFlavor::NUM_UNSHIFTED_ENTITIES;
    static constexpr size_t NUM_SHIFTED_ENTITIES = NativeFlavor::NUM_SHIFTED_ENTITIES;
    static constexpr size_t TRACE_OFFSET = NativeFlavor::TRACE_OFFSET;

    using Relations = NativeFlavor::Relations_<FF>;

    static constexpr size_t MAX_PARTIAL_RELATION_LENGTH = compute_max_partial_relation_length<Relations>();
    static constexpr size_t BATCHED_RELATION_PARTIAL_LENGTH = MAX_PARTIAL_RELATION_LENGTH + 1;
    static constexpr size_t NUM_RELATIONS = std::tuple_size_v<Relations>;
    static constexpr size_t NUM_SUBRELATIONS = NativeFlavor::NUM_SUBRELATIONS;
    using SubrelationSeparator = FF;

    template <typename DataType> using AllEntities = NativeFlavor::AllEntities<DataType>;

    class AllValues : public AllEntities<FF> {
      public:
        using Base = AllEntities<FF>;
        using Base::Base;
    };

    using CommitmentLabels = NativeFlavor::CommitmentLabels;
    static const CommitmentLabels& commitment_labels() { return NativeFlavor::commitment_labels(); }
};

} // namespace bb
