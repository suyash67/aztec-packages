#pragma once
#include "barretenberg/blake3vm/blake3vm_flavor.hpp"
#include "barretenberg/honk/proof_system/types/proof.hpp"
#include "barretenberg/relations/relation_parameters.hpp"

namespace bb {

/**
 * @brief Proves a Blake3VM trace.
 *
 * Transcript schedule: the VK (circuit size + precomputed commitments) enters the hash buffer;
 * the 119 wire commitments are sent; challenges "beta"/"gamma" produce the six logup inverse
 * columns, whose commitments follow; "Sumcheck:alpha" and "Sumcheck:gate_challenge" drive the
 * non-ZK sumcheck; Shplemini and a KZG opening close the proof.
 */
class Blake3VMProver {
  public:
    using Flavor = Blake3VMFlavor;
    using FF = Flavor::FF;
    using CircuitBuilder = Flavor::CircuitBuilder;
    using ProvingKey = Flavor::ProvingKey;
    using VerificationKey = Flavor::VerificationKey;
    using Transcript = Flavor::Transcript;

    explicit Blake3VMProver(const CircuitBuilder& builder,
                            const std::shared_ptr<Transcript>& transcript = std::make_shared<Transcript>());

    HonkProof construct_proof();

    std::shared_ptr<ProvingKey> key;
    std::shared_ptr<VerificationKey> verification_key;
    std::shared_ptr<Transcript> transcript;
    RelationParameters<FF> relation_parameters;
};

} // namespace bb
