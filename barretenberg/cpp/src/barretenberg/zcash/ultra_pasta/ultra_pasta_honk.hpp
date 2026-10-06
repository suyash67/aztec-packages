#pragma once

#include "barretenberg/ultra_honk/prover_instance.hpp"
#include "barretenberg/zcash/ultra_pasta/ultra_pasta_flavor.hpp"

#include <memory>

/**
 * @file ultra_pasta_honk.hpp
 * @brief Proving and verification of UltraPastaZKFlavor circuits.
 * @details The prover is UltraHonk's: Oink (wire, lookup, memory and permutation commitments), ZK Sumcheck,
 * SmallSubgroupIPA for the Libra masking, Shplemini, and the halo2 inner-product argument as the PCS. The proof holds
 * the public inputs (sent by Oink).
 */
namespace bb::zcash {

using UltraPastaProverInstance = ProverInstance_<UltraPastaZKFlavor>;

UltraPastaZKFlavor::Proof ultra_pasta_prove(const std::shared_ptr<UltraPastaProverInstance>& instance,
                                            const std::shared_ptr<UltraPastaZKFlavor::VerificationKey>& vk);

/**
 * @brief Verifies a proof against `vk`; on success the proof's public inputs are written to `public_inputs`.
 */
bool ultra_pasta_verify(const std::shared_ptr<UltraPastaZKFlavor::VerificationKey>& vk,
                        const UltraPastaZKFlavor::Proof& proof,
                        std::vector<pasta::fp>* public_inputs = nullptr);

} // namespace bb::zcash
