#include "barretenberg/zcash/ultra_pasta/ultra_pasta_honk.hpp"
#include "barretenberg/numeric/random/engine.hpp"

#include <gtest/gtest.h>

using namespace bb;
using namespace bb::zcash;

namespace {
using FF = pasta::fp;
using Point = pallas::g1::affine_element;
using Builder = UltraPastaCircuitBuilder;

// A circuit using every gate family the Orchard circuit uses: arithmetic, Pallas addition and doubling, ROM and
// range constraints, and a public input.
Builder small_circuit(bool corrupt = false)
{
    Builder builder;
    const FF a = FF::random_element();
    const FF b = FF::random_element();
    const auto a_idx = builder.add_variable(a);
    const auto b_idx = builder.add_variable(b);
    const auto c_idx = builder.add_variable(a * b);
    builder.create_big_mul_add_gate(
        { a_idx, b_idx, c_idx, builder.zero_idx(), FF(1), FF(0), FF(0), FF(-1), FF(0), FF(0) });
    const auto pub = builder.add_public_variable(a * b);
    builder.assert_equal(pub, c_idx);

    const Point p = pallas::g1::one;
    const Point q = Point(pallas::g1::element(p) * pallas::fr(7));
    const Point r = Point(pallas::g1::element(p) + pallas::g1::element(q));
    const Point d = Point(pallas::g1::element(r).dbl());
    std::array<uint32_t, 8> idx{};
    const std::array<FF, 8> coords{ p.x, p.y, q.x, q.y, r.x, r.y, d.x, d.y + (corrupt ? FF(1) : FF(0)) };
    for (size_t i = 0; i < 8; ++i) {
        idx[i] = builder.add_variable(coords[i]);
    }
    builder.create_ecc_add_gate({ idx[0], idx[1], idx[2], idx[3], idx[4], idx[5], true });
    builder.create_ecc_dbl_gate({ idx[4], idx[5], idx[6], idx[7] });

    const size_t rom = builder.create_ROM_array(4);
    for (size_t i = 0; i < 4; ++i) {
        builder.set_ROM_element(rom, i, idx[i]);
    }
    const auto read = builder.read_ROM_array(rom, builder.add_variable(FF(2)));
    builder.assert_equal(read, idx[2]);

    builder.create_dyadic_range_constraint(builder.add_variable(FF((1 << 20) - 1)), 20, "range");
    builder.create_dyadic_range_constraint(builder.add_variable(FF(uint256_t(1) << 100)), 128, "range");
    return builder;
}
} // namespace

TEST(UltraPastaHonk, ProveAndVerify)
{
    auto builder = small_circuit();
    EXPECT_FALSE(builder.failed());
    auto instance = std::make_shared<UltraPastaProverInstance>(builder);
    auto vk = std::make_shared<UltraPastaZKFlavor::VerificationKey>(instance->get_precomputed());
    auto proof = ultra_pasta_prove(instance, vk);
    std::vector<FF> public_inputs;
    EXPECT_TRUE(ultra_pasta_verify(vk, proof, &public_inputs));
    ASSERT_EQ(public_inputs.size(), 1U);
}

TEST(UltraPastaHonk, UnsatisfiedCircuitFails)
{
    auto builder = small_circuit(/*corrupt=*/true);
    auto instance = std::make_shared<UltraPastaProverInstance>(builder);
    auto vk = std::make_shared<UltraPastaZKFlavor::VerificationKey>(instance->get_precomputed());
    auto proof = ultra_pasta_prove(instance, vk);
    EXPECT_FALSE(ultra_pasta_verify(vk, proof));
}

// The VK hash checked by the verifier is the full F_p value that Oink hashes into the transcript (not reduced modulo
// the BN254 scalar field).
TEST(UltraPastaHonk, VerificationKeyHashIsTranscriptHash)
{
    auto builder = small_circuit();
    auto instance = std::make_shared<UltraPastaProverInstance>(builder);
    auto vk = std::make_shared<UltraPastaZKFlavor::VerificationKey>(instance->get_precomputed());
    const UltraPastaZKFlavor::VKAndHash vk_and_hash(vk);
    EXPECT_EQ(vk_and_hash.hash, FF(vk->hash_with_origin_tagging(OriginTag())));
    EXPECT_EQ(uint256_t(vk_and_hash.hash), UltraPastaZKFlavor::HashFunction::hash(vk->to_field_elements()));
}
