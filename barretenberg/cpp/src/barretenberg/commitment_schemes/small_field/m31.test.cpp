#include "barretenberg/commitment_schemes/small_field/m31.hpp"

#include <gtest/gtest.h>

namespace bb::small_field {

namespace {
constexpr uint64_t P = 0x7fffffffULL;

uint32_t reference_mul(uint32_t a, uint32_t b)
{
    return static_cast<uint32_t>((static_cast<uint64_t>(a) * b) % P);
}
} // namespace

TEST(M31, ArithmeticMatchesReference)
{
    for (size_t trial = 0; trial < 1000; ++trial) {
        const m31 a = m31::random_element();
        const m31 b = m31::random_element();
        EXPECT_EQ((a + b).value(), (static_cast<uint64_t>(a.value()) + b.value()) % P);
        EXPECT_EQ((a - b).value(), (static_cast<uint64_t>(a.value()) + P - b.value()) % P);
        EXPECT_EQ((a * b).value(), reference_mul(a.value(), b.value()));
    }
}

TEST(M31, EdgeCases)
{
    const m31 max(P - 1);
    EXPECT_EQ((max + m31(1)).value(), 0U);
    EXPECT_EQ((max * max).value(), reference_mul(static_cast<uint32_t>(P - 1), static_cast<uint32_t>(P - 1)));
    EXPECT_EQ((-m31(0)).value(), 0U);
    EXPECT_EQ((m31(0) - m31(1)).value(), P - 1);
    // Unreduced constructor inputs
    EXPECT_EQ(m31(P).value(), 0U);
    EXPECT_EQ(m31(P + 5).value(), 5U);
    EXPECT_EQ(m31(0xffffffffffffffffULL).value(), 0xffffffffffffffffULL % P);
}

TEST(M31, InversionAndPow)
{
    for (size_t trial = 0; trial < 100; ++trial) {
        const m31 a = m31::random_element();
        if (a == m31::zero()) {
            continue;
        }
        EXPECT_EQ(a * a.invert(), m31::one());
    }
    // Fermat: a^{p-1} = 1
    EXPECT_EQ(m31(12345).pow(P - 1), m31::one());
}

TEST(CM31, ExtensionStructure)
{
    // i^2 = -1
    const cm31 i(m31::zero(), m31::one());
    EXPECT_EQ(i * i, cm31(-m31::one(), m31::zero()));
    for (size_t trial = 0; trial < 100; ++trial) {
        const cm31 a = cm31::random_element();
        if (a == cm31::zero()) {
            continue;
        }
        EXPECT_EQ(a * a.invert(), cm31::one());
        EXPECT_EQ(a * a.conjugate(), cm31(a.norm(), m31::zero()));
    }
}

TEST(QM31, ExtensionStructure)
{
    // u^2 = 2 + i
    const qm31 u(cm31::zero(), cm31::one());
    EXPECT_EQ(u * u, qm31(cm31(m31(2), m31(1)), cm31::zero()));
    for (size_t trial = 0; trial < 100; ++trial) {
        const qm31 a = qm31::random_element();
        const qm31 b = qm31::random_element();
        const qm31 c = qm31::random_element();
        if (a == qm31::zero()) {
            continue;
        }
        EXPECT_EQ(a * a.invert(), qm31::one());
        // Distributivity and commutativity spot checks
        EXPECT_EQ(a * (b + c), a * b + a * c);
        EXPECT_EQ(a * b, b * a);
    }
}

TEST(QM31, BaseFieldScaling)
{
    const qm31 a = qm31::random_element();
    const m31 s = m31::random_element();
    const qm31 s_embedded(cm31(s, m31::zero()), cm31::zero());
    EXPECT_EQ(a.scale(s), a * s_embedded);
}

} // namespace bb::small_field
