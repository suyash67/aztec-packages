#pragma once

#include "barretenberg/ecc/curves/pasta/pasta.hpp"

#include <cstdint>
#include <span>
#include <string_view>

namespace bb::zcash {

/**
 * @brief GroupHash^P for the Pasta curves: the hash_to_curve construction of the pasta_curves crate
 * (draft-irtf-cfrg-hash-to-curve "XMD:BLAKE2b_SSWU_RO" with the 3-isogenous curves iso-Pallas / iso-Vesta).
 * @details Every Orchard generator (Sinsemilla S/Q/R points, SpendAuthG, NullifierK, ValueCommit V/R) is defined as
 * `pallas_hash_to_curve(domain)(message)`, and halo2's IPA parameters use `vesta_hash_to_curve("Halo2-Parameters")`.
 */
pallas::g1::affine_element pallas_hash_to_curve(std::string_view domain_prefix, std::span<const uint8_t> message);
vesta::g1::affine_element vesta_hash_to_curve(std::string_view domain_prefix, std::span<const uint8_t> message);

inline pallas::g1::affine_element pallas_hash_to_curve(std::string_view domain_prefix, std::string_view message)
{
    return pallas_hash_to_curve(
        domain_prefix, std::span<const uint8_t>(reinterpret_cast<const uint8_t*>(message.data()), message.size()));
}

} // namespace bb::zcash
