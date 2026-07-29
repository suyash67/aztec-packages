#pragma once

#include "barretenberg/cq/cq_key.hpp"
#include "barretenberg/transcript/transcript.hpp"

namespace bb::cq {

/**
 * @brief Verifier for the cq lookup argument. Verification is O(1) group operations and 7 pairings, independent of
 * both the table size and the number of lookups (up to a log-size exponentiation for Z_V(gamma)).
 */
class CqVerifier {
  public:
    using Transcript = NativeTranscript;
    using Proof = Transcript::Proof;

    explicit CqVerifier(const CqVerificationKey& verification_key)
        : vk_(verification_key)
    {}

    bool verify_proof(const Proof& proof) const;

  private:
    const CqVerificationKey& vk_;
};

} // namespace bb::cq
