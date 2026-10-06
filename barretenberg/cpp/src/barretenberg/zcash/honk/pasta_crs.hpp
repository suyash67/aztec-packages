#pragma once

#include "barretenberg/ecc/curves/pasta/pasta.hpp"
#include "barretenberg/srs/factories/crs_factory.hpp"
#include "barretenberg/srs/global_crs.hpp"

#include <memory>
#include <mutex>
#include <vector>

/**
 * @file pasta_crs.hpp
 * @brief Transparent commitment parameters on Vesta, generated exactly like halo2's IPA `Params<vesta::Affine>`:
 * G_i = vesta_hash_to_curve("Halo2-Parameters")(0x00 || LE32(i)), plus the blinding generator W = hash(0x01) and the
 * IPA generator U = hash(0x02).
 */
namespace bb::zcash {

/// The first `n` halo2 Vesta generators G_0..G_{n-1}.
std::vector<vesta::g1::affine_element> halo2_vesta_generators(size_t n);
/// halo2's blinding generator W.
vesta::g1::affine_element halo2_vesta_w();
/// halo2's IPA generator U.
vesta::g1::affine_element halo2_vesta_u();

/**
 * @brief A CRS factory serving the halo2 Vesta generators; the generators are derived on demand and cached, and the
 * set only ever grows.
 */
class PastaCrsFactory : public srs::factories::CrsFactory<curve::Vesta> {
  public:
    std::shared_ptr<srs::factories::Crs<curve::Vesta>> get_crs(size_t degree) override;

  private:
    std::mutex mutex_;
    std::shared_ptr<srs::factories::Crs<curve::Vesta>> crs_;
};

} // namespace bb::zcash

namespace bb::srs {
template <> std::shared_ptr<factories::CrsFactory<curve::Vesta>> get_crs_factory<curve::Vesta>();
} // namespace bb::srs
