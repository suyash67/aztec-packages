#include "pasta.hpp"
#include "barretenberg/ecc/scalar_multiplication/scalar_multiplication.hpp"
#include <gtest/gtest.h>

using namespace bb;

namespace {

// Double-and-add without any curve-specific shortcuts, used as the reference for scalar multiplication.
template <typename Group> typename Group::affine_element naive_mul(const typename Group::element& base, uint256_t k)
{
    typename Group::element acc = Group::point_at_infinity;
    typename Group::element runner = base;
    for (size_t i = 0; i < 256; ++i) {
        if (k.get_bit(i)) {
            acc = acc + runner;
        }
        runner = runner.dbl();
    }
    return typename Group::affine_element(acc);
}

} // namespace

TEST(Pasta, ModuliAreTheCycle)
{
    // The scalar field of each curve is the base field of the other.
    EXPECT_EQ(pallas::fr::modulus, vesta::fq::modulus);
    EXPECT_EQ(pallas::fq::modulus, vesta::fr::modulus);
    EXPECT_EQ(pallas::fq::modulus, uint256_t("0x40000000000000000000000000000000224698fc094cf91b992d30ed00000001"));
    EXPECT_EQ(pallas::fr::modulus, uint256_t("0x40000000000000000000000000000000224698fc0994a8dd8c46eb2100000001"));
}

TEST(Pasta, GeneratorsMatchPastaCurves)
{
    // pasta_curves uses (-1, 2) as the generator of both curves.
    EXPECT_EQ(pallas::g1::affine_one.x, -pallas::fq::one());
    EXPECT_EQ(pallas::g1::affine_one.y, pallas::fq(2));
    EXPECT_EQ(vesta::g1::affine_one.x, -vesta::fq::one());
    EXPECT_EQ(vesta::g1::affine_one.y, vesta::fq(2));
    EXPECT_TRUE(pallas::g1::affine_one.on_curve());
    EXPECT_TRUE(vesta::g1::affine_one.on_curve());
}

TEST(Pasta, GroupOrders)
{
    // [r]G = O where r is the scalar field modulus; [r - 1]G = -G.
    auto pallas_minus_one = naive_mul<pallas::g1>(pallas::g1::one, pallas::fr::modulus - 1);
    EXPECT_EQ(pallas_minus_one, -pallas::g1::affine_one);
    auto vesta_minus_one = naive_mul<vesta::g1>(vesta::g1::one, vesta::fr::modulus - 1);
    EXPECT_EQ(vesta_minus_one, -vesta::g1::affine_one);
}

TEST(Pasta, CubeRootsMatchPastaCurvesZeta)
{
    // pasta_curves Fp::ZETA
    EXPECT_EQ(uint256_t(pasta::fp::cube_root_of_unity()),
              uint256_t("0x12ccca834acdba712caad5dc57aab1b01d1f8bd237ad31491dad5ebdfdfe4ab9"));
    EXPECT_EQ(pasta::fp::cube_root_of_unity().pow(3), pasta::fp::one());
    EXPECT_EQ(pasta::fq::cube_root_of_unity().pow(3), pasta::fq::one());
    EXPECT_NE(pasta::fp::cube_root_of_unity(), pasta::fp::one());
    EXPECT_NE(pasta::fq::cube_root_of_unity(), pasta::fq::one());
}

TEST(Pasta, EndomorphismIsConsistentOnBothCurves)
{
    // (x, y) -> (beta * x, y) equals [lambda] for the stored cube roots, on both curves.
    for (size_t i = 0; i < 4; ++i) {
        auto P = pallas::g1::element::random_element();
        pallas::g1::affine_element P_aff(P);
        auto lhs = naive_mul<pallas::g1>(P, uint256_t(pallas::fr::cube_root_of_unity()));
        EXPECT_EQ(lhs, pallas::g1::affine_element(P_aff.x * pallas::fq::cube_root_of_unity(), P_aff.y));

        auto Q = vesta::g1::element::random_element();
        vesta::g1::affine_element Q_aff(Q);
        auto rhs = naive_mul<vesta::g1>(Q, uint256_t(vesta::fr::cube_root_of_unity()));
        EXPECT_EQ(rhs, vesta::g1::affine_element(Q_aff.x * vesta::fq::cube_root_of_unity(), Q_aff.y));
    }
}

TEST(Pasta, ScalarMulMatchesNaive)
{
    for (size_t i = 0; i < 8; ++i) {
        auto k = pallas::fr::random_element();
        auto P = pallas::g1::element::random_element();
        EXPECT_EQ(pallas::g1::affine_element(P * k), naive_mul<pallas::g1>(P, uint256_t(k)));
        auto s = vesta::fr::random_element();
        auto Q = vesta::g1::element::random_element();
        EXPECT_EQ(vesta::g1::affine_element(Q * s), naive_mul<vesta::g1>(Q, uint256_t(s)));
    }
}

TEST(Pasta, SubgroupGenerators)
{
    auto check = [](auto gen, auto gen_inv, auto one) {
        EXPECT_EQ(gen * gen_inv, one);
        EXPECT_EQ(gen.pow(256), one);
        EXPECT_NE(gen.pow(128), one);
    };
    check(curve::Vesta::subgroup_generator, curve::Vesta::subgroup_generator_inverse, vesta::fr::one());
    check(curve::Pallas::subgroup_generator, curve::Pallas::subgroup_generator_inverse, pallas::fr::one());
}

TEST(Pasta, FieldArithmeticAgainstUint512)
{
    // Montgomery multiplication on the large-modulus code path, checked against schoolbook 512-bit arithmetic.
    auto check = [](auto a, auto b) {
        using F = decltype(a);
        uint512_t prod = uint512_t(uint256_t(a)) * uint512_t(uint256_t(b));
        uint256_t expected = (prod % uint512_t(F::modulus)).lo;
        EXPECT_EQ(uint256_t(a * b), expected);
        uint512_t sum = uint512_t(uint256_t(a)) + uint512_t(uint256_t(b));
        EXPECT_EQ(uint256_t(a + b), (sum % uint512_t(F::modulus)).lo);
        EXPECT_EQ(a * a.invert(), F::one());
    };
    for (size_t i = 0; i < 64; ++i) {
        check(pasta::fp::random_element(), pasta::fp::random_element());
        check(pasta::fq::random_element(), pasta::fq::random_element());
    }
    // values close to the modulus
    check(pasta::fp(pasta::fp::modulus - 1), pasta::fp(pasta::fp::modulus - 2));
    check(pasta::fq(pasta::fq::modulus - 1), pasta::fq(pasta::fq::modulus - 1));
}

TEST(Pasta, PippengerMatchesNaive)
{
    // The public MSM facade on Vesta (the commitment curve for circuits over F_p), across thread/window regimes.
    for (size_t n : { size_t{ 1 }, size_t{ 7 }, size_t{ 64 }, size_t{ 1000 }, size_t{ 1 } << 13 }) {
        std::vector<vesta::g1::affine_element> points(n);
        std::vector<vesta::fr> scalars(n);
        vesta::g1::element expected = vesta::g1::point_at_infinity;
        for (size_t i = 0; i < n; ++i) {
            points[i] = vesta::g1::affine_element(vesta::g1::element::random_element());
            scalars[i] = vesta::fr::random_element();
            expected += vesta::g1::element(points[i]) * scalars[i];
        }
        PolynomialSpan<const vesta::fr> span(0, scalars);
        auto result = scalar_multiplication::pippenger<curve::Vesta>(span, points);
        EXPECT_EQ(vesta::g1::affine_element(result), vesta::g1::affine_element(expected)) << "n = " << n;
    }
}
