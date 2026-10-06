#pragma once

#include "barretenberg/ecc/curves/bn254/bn254.hpp"
#include "barretenberg/ecc/curves/grumpkin/grumpkin.hpp"
#include "barretenberg/ecc/curves/pasta/pasta.hpp"
#include "group_hash.hpp"

#include <span>
#include <string>
#include <string_view>

namespace bb::zcash {

/**
 * @brief The Pasta instantiation of the Action circuit, as deployed by Zcash: circuit field F_p (Pallas base), embedded
 * curve Pallas, polynomial commitments on Vesta (whose scalar field is F_p).
 */
struct PastaCycle {
    using FF = pallas::fq;
    using Group = pallas::g1;
    using AffineElement = Group::affine_element;
    using Element = Group::element;
    using EmbeddedScalar = pallas::fr;
    using CommitmentCurve = curve::Vesta;

    static constexpr const char* NAME = "pasta";

    static AffineElement group_hash(std::string_view domain, std::span<const uint8_t> message)
    {
        return pallas_hash_to_curve(domain, message);
    }
};

/**
 * @brief The BN254 instantiation of the Action circuit: circuit field F_r (BN254 scalar field = Grumpkin base field),
 * embedded curve Grumpkin, polynomial commitments on BN254 G1.
 * @details Grumpkin has no Orchard-defined generators, so GroupHash is barretenberg's Grumpkin hash-to-curve applied
 * to the same domain strings (Orchard's personalisation and message are concatenated into the hash seed).
 */
struct Bn254Cycle {
    using FF = bb::fr;
    using Group = grumpkin::g1;
    using AffineElement = Group::affine_element;
    using Element = Group::element;
    using EmbeddedScalar = grumpkin::fr;
    using CommitmentCurve = curve::BN254;

    static constexpr const char* NAME = "bn254";

    static AffineElement group_hash(std::string_view domain, std::span<const uint8_t> message)
    {
        std::vector<uint8_t> seed(domain.begin(), domain.end());
        seed.push_back(0);
        seed.insert(seed.end(), message.begin(), message.end());
        return AffineElement::hash_to_curve(seed);
    }
};

template <typename Cycle> typename Cycle::AffineElement group_hash(std::string_view domain, std::string_view message)
{
    return Cycle::group_hash(
        domain, std::span<const uint8_t>(reinterpret_cast<const uint8_t*>(message.data()), message.size()));
}

} // namespace bb::zcash
