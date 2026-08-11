#pragma once
#include "barretenberg/blake3vm/blake3vm_flavor.hpp"
#include "barretenberg/honk/proof_system/types/proof.hpp"
#include "barretenberg/relations/relation_parameters.hpp"

namespace bb {

/** @brief Verifies a Blake3VM proof against a verification key. Mirrors Blake3VMProver's transcript schedule. */
class Blake3VMVerifier {
  public:
    using Flavor = Blake3VMFlavor;
    using FF = Flavor::FF;
    using VerificationKey = Flavor::VerificationKey;
    using Transcript = Flavor::Transcript;

    explicit Blake3VMVerifier(const std::shared_ptr<VerificationKey>& verification_key,
                              const std::shared_ptr<Transcript>& transcript = std::make_shared<Transcript>());

    bool verify_proof(const HonkProof& proof);

    std::shared_ptr<VerificationKey> verification_key;
    std::shared_ptr<Transcript> transcript;
    RelationParameters<FF> relation_parameters;
};

} // namespace bb
