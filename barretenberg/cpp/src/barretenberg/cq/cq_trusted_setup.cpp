#include "cq_trusted_setup.hpp"

#include "barretenberg/common/assert.hpp"
#include "barretenberg/common/thread.hpp"
#include "barretenberg/ecc/curves/bn254/fr.hpp"
#include "barretenberg/srs/factories/mem_bn254_crs_factory.hpp"

namespace bb::cq {

namespace {

// CommitmentKey's only constructor pulls from the init-once global CRS factory; this subclass sets the protected
// srs member from a local mem factory instead, so each TestSrs gets a self-contained key.
class MemCommitmentKey : public CommitmentKey<curve::BN254> {
  public:
    MemCommitmentKey(const std::vector<g1::affine_element>& points, const g2::affine_element& g2_point)
    {
        srs::factories::MemBn254CrsFactory factory(points, g2_point);
        srs = factory.get_crs(points.size());
        srs_size = points.size();
    }
};

} // namespace

TestSrs TestSrs::create(size_t num_g1_points, size_t num_g2_points)
{
    const fr tau = fr::random_element();
    const size_t max_power = std::max(num_g1_points, num_g2_points);
    std::vector<fr> tau_powers(max_power);
    tau_powers[0] = fr::one();
    for (size_t i = 1; i < max_power; ++i) {
        tau_powers[i] = tau_powers[i - 1] * tau;
    }

    TestSrs srs;

    std::vector<g1::element> g1_jacobian(num_g1_points);
    parallel_for_range(num_g1_points, [&](size_t start, size_t end) {
        for (size_t i = start; i < end; ++i) {
            g1_jacobian[i] = g1::one * tau_powers[i];
        }
    });
    g1::element::batch_normalize(g1_jacobian.data(), num_g1_points);
    srs.g1_powers.resize(num_g1_points);
    for (size_t i = 0; i < num_g1_points; ++i) {
        srs.g1_powers[i] = { g1_jacobian[i].x, g1_jacobian[i].y };
    }

    std::vector<g2::element> g2_jacobian(num_g2_points);
    parallel_for_range(num_g2_points, [&](size_t start, size_t end) {
        for (size_t i = start; i < end; ++i) {
            g2_jacobian[i] = g2::one * tau_powers[i];
        }
    });
    g2::element::batch_normalize(g2_jacobian.data(), num_g2_points);
    srs.g2_powers.resize(num_g2_points);
    for (size_t i = 0; i < num_g2_points; ++i) {
        srs.g2_powers[i] = { g2_jacobian[i].x, g2_jacobian[i].y };
    }

    return srs;
}

CommitmentKey<curve::BN254> TestSrs::create_commitment_key() const
{
    BB_ASSERT_GTE(g2_powers.size(), 2U, "TestSrs: need at least [1]_2 and [tau]_2");
    return MemCommitmentKey(g1_powers, g2_powers[1]);
}

} // namespace bb::cq
