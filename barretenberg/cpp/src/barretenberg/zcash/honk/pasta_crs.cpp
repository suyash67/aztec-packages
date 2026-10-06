#include "pasta_crs.hpp"
#include "barretenberg/common/thread.hpp"
#include "barretenberg/zcash/primitives/group_hash.hpp"

namespace bb::zcash {

namespace {

constexpr const char* PARAMETERS_DOMAIN = "Halo2-Parameters";

class PastaCrs : public srs::factories::Crs<curve::Vesta> {
  public:
    explicit PastaCrs(std::vector<vesta::g1::affine_element> points)
        : points_(std::move(points))
    {}
    std::span<vesta::g1::affine_element> get_monomial_points() override { return points_; }
    size_t get_monomial_size() const override { return points_.size(); }
    vesta::g1::affine_element get_g1_identity() const override { return points_[0]; }

  private:
    std::vector<vesta::g1::affine_element> points_;
};

} // namespace

std::vector<vesta::g1::affine_element> halo2_vesta_generators(size_t n)
{
    std::vector<vesta::g1::affine_element> out(n);
    parallel_for(n, [&](size_t i) {
        const auto idx = static_cast<uint32_t>(i);
        const std::array<uint8_t, 5> message = { 0,
                                                 static_cast<uint8_t>(idx),
                                                 static_cast<uint8_t>(idx >> 8),
                                                 static_cast<uint8_t>(idx >> 16),
                                                 static_cast<uint8_t>(idx >> 24) };
        out[i] = vesta_hash_to_curve(PARAMETERS_DOMAIN, message);
    });
    return out;
}

vesta::g1::affine_element halo2_vesta_w()
{
    const std::array<uint8_t, 1> msg = { 1 };
    return vesta_hash_to_curve(PARAMETERS_DOMAIN, msg);
}

vesta::g1::affine_element halo2_vesta_u()
{
    const std::array<uint8_t, 1> msg = { 2 };
    return vesta_hash_to_curve(PARAMETERS_DOMAIN, msg);
}

std::shared_ptr<srs::factories::Crs<curve::Vesta>> PastaCrsFactory::get_crs(size_t degree)
{
    std::lock_guard<std::mutex> lock(mutex_);
    if (!crs_ || crs_->get_monomial_size() < degree) {
        size_t size = 1;
        while (size < degree) {
            size <<= 1;
        }
        crs_ = std::make_shared<PastaCrs>(halo2_vesta_generators(size));
    }
    return crs_;
}

} // namespace bb::zcash

namespace bb::srs {
template <> std::shared_ptr<factories::CrsFactory<curve::Vesta>> get_crs_factory<curve::Vesta>()
{
    static auto factory = std::make_shared<zcash::PastaCrsFactory>();
    return factory;
}
} // namespace bb::srs
