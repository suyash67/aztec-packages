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
struct OrchardProvingKey {
    using Flavor = OrchardFlavor;
    using FF = Flavor::FF;
    size_t circuit_size = 0;
    size_t log_circuit_size = 0;
    Flavor::ProverPolynomials precomputed; // only the precomputed polynomials are populated
    std::shared_ptr<Flavor::VerificationKey> vk;
    std::vector<std::pair<size_t, size_t>> public_input_cells;

    explicit OrchardProvingKey(const halo2::AnchoredTrace<PastaCycle>& trace);
};

/**
 * @brief Proves one execution of the circuit: the witness is read from the advice columns of `trace`, which must have
 * the shape the proving key was generated from.
 */
OrchardFlavor::Proof orchard_prove(const OrchardProvingKey& pk, const halo2::AnchoredTrace<PastaCycle>& trace);

bool orchard_verify(const OrchardFlavor::VerificationKey& vk,
                    const std::vector<OrchardFlavor::FF>& public_inputs,
                    const OrchardFlavor::Proof& proof);

} // namespace bb::zcash
