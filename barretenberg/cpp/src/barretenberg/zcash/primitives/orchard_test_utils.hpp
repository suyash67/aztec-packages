#pragma once

#include "barretenberg/zcash/test_vectors/orchard_vectors.hpp"
#include "orchard.hpp"

#include <map>
#include <string>

namespace bb::zcash {

/**
 * @brief Converts an exported Action witness (name -> hex/decimal strings) into the native witness struct.
 */
template <typename Cycle>
typename Orchard<Cycle>::ActionWitness parse_action_witness(const test_vectors::ActionVector& vector)
{
    using O = Orchard<Cycle>;
    using FF = typename O::FF;
    using Scalar = typename O::Scalar;
    using AffineElement = typename O::AffineElement;
    std::map<std::string, std::string> kv(vector.witness.begin(), vector.witness.end());
    auto fe = [&](const std::string& k) { return FF(uint256_t(kv.at(k))); };
    auto sc = [&](const std::string& k) { return Scalar(uint256_t(kv.at(k))); };
    auto pt = [&](const std::string& k) {
        const uint256_t x(kv.at(k + "_x"));
        const uint256_t y(kv.at(k + "_y"));
        if (x == 0 && y == 0) {
            return AffineElement::infinity();
        }
        return AffineElement(FF(x), FF(y));
    };
    typename O::ActionWitness w;
    for (size_t i = 0; i < O::MERKLE_DEPTH; ++i) {
        w.path[i] = fe("path_" + std::to_string(i));
    }
    w.pos = static_cast<uint32_t>(std::stoul(kv.at("pos")));
    w.g_d_old = pt("g_d_old");
    w.pk_d_old = pt("pk_d_old");
    w.v_old = std::stoull(kv.at("v_old"));
    w.rho_old = fe("rho_old");
    w.psi_old = fe("psi_old");
    w.rcm_old = sc("rcm_old");
    w.cm_old = pt("cm_old");
    w.alpha = sc("alpha");
    w.ak = pt("ak");
    w.nk = fe("nk");
    w.rivk = sc("rivk");
    w.g_d_new = pt("g_d_new");
    w.pk_d_new = pt("pk_d_new");
    w.v_new = std::stoull(kv.at("v_new"));
    w.psi_new = fe("psi_new");
    w.rcm_new = sc("rcm_new");
    w.rcv = sc("rcv");
    return w;
}

} // namespace bb::zcash
