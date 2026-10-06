#pragma once

#include "barretenberg/zcash/halo2/trace.hpp"
#include "barretenberg/zcash/honk/orchard_flavor.hpp"

#include <memory>
#include <vector>

/**
 * @file orchard_honk.hpp
 * @brief Key generation, proving and verification for the halo2-arithmetized Orchard Action circuit (OrchardFlavor).
 */
namespace bb::zcash {

/**
 * @brief Proving key: the precomputed polynomials of one circuit shape and its verification key.
 * @details The circuit shape (number of Actions, layout) fixes every precomputed polynomial: fixed columns, selectors,
 * table, permutation and Lagrange polynomials. The witness-independent part of `AnchoredTrace` is all that is used.
 */
template <typename Cycle> struct OrchardProvingKey_ {
    using Flavor = OrchardFlavor_<Cycle>;
    using FF = typename Flavor::FF;
    size_t circuit_size = 0;
    size_t log_circuit_size = 0;
    typename Flavor::ProverPolynomials precomputed; // only the precomputed polynomials are populated
    std::shared_ptr<typename Flavor::VerificationKey> vk;
    std::vector<std::pair<size_t, size_t>> public_input_cells;

    explicit OrchardProvingKey_(const halo2::AnchoredTrace<Cycle>& trace);
};
using OrchardProvingKey = OrchardProvingKey_<PastaCycle>;

/**
 * @brief Proves one execution of the circuit: the witness is read from the advice columns of `trace`, which must have
 * the shape the proving key was generated from.
 */
template <typename Cycle>
typename OrchardFlavor_<Cycle>::Proof orchard_prove(const OrchardProvingKey_<Cycle>& pk,
                                                    const halo2::AnchoredTrace<Cycle>& trace);

template <typename Cycle>
bool orchard_verify(const typename OrchardFlavor_<Cycle>::VerificationKey& vk,
                    const std::vector<typename Cycle::FF>& public_inputs,
                    const typename OrchardFlavor_<Cycle>::Proof& proof);

inline bool orchard_verify(const OrchardFlavor::VerificationKey& vk,
                           const std::vector<PastaCycle::FF>& public_inputs,
                           const OrchardFlavor::Proof& proof)
{
    return orchard_verify<PastaCycle>(vk, public_inputs, proof);
}
inline bool orchard_verify(const OrchardBn254Flavor::VerificationKey& vk,
                           const std::vector<Bn254Cycle::FF>& public_inputs,
                           const OrchardBn254Flavor::Proof& proof)
{
    return orchard_verify<Bn254Cycle>(vk, public_inputs, proof);
}

} // namespace bb::zcash
