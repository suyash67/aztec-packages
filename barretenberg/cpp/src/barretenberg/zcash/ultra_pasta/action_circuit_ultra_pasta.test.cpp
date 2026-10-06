#include "barretenberg/zcash/ultra_pasta/action_circuit_ultra_pasta.hpp"
#include "barretenberg/numeric/random/engine.hpp"
#include "barretenberg/zcash/primitives/orchard_test_utils.hpp"
#include "barretenberg/zcash/ultra_pasta/ultra_pasta_honk.hpp"

#include <gtest/gtest.h>

using namespace bb;
using namespace bb::zcash;
using namespace bb::zcash::ultra_pasta;

namespace {
using O = Orchard<PastaCycle>;
using FF = PastaCycle::FF;

std::pair<std::vector<O::ActionWitness>, std::vector<FF>> random_bundle(size_t n)
{
    auto& engine = numeric::get_debug_randomness();
    std::vector<O::ActionWitness> ws;
    std::vector<FF> pis;
    for (size_t i = 0; i < n; ++i) {
        ws.push_back(O::random_witness(engine));
        auto pi = O::evaluate(ws.back());
        EXPECT_TRUE(pi.has_value());
        for (const auto& f : pi->to_field_elements()) {
            pis.push_back(f);
        }
    }
    return { ws, pis };
}

std::vector<FF> vector_public_inputs(const test_vectors::ActionVector& v)
{
    std::vector<FF> out;
    for (const auto& s : v.public_inputs) {
        out.emplace_back(uint256_t(s));
    }
    return out;
}

// Proves and verifies; returns whether the proof verified with the given public inputs.
bool prove_and_verify(Builder& builder, const std::vector<FF>& pis)
{
    auto instance = std::make_shared<UltraPastaProverInstance>(builder);
    auto vk = std::make_shared<UltraPastaZKFlavor::VerificationKey>(instance->get_precomputed());
    auto proof = ultra_pasta_prove(instance, vk);
    info("UltraZK Pasta Action: log_n = ", vk->log_circuit_size, ", proof size: ", proof.size() * 32, " bytes");
    std::vector<FF> proof_pis;
    const bool verified = ultra_pasta_verify(vk, proof, &proof_pis);
    return verified && proof_pis == pis;
}
} // namespace

TEST(ZcashActionUltraPasta, RandomActionProveAndVerify)
{
    auto [ws, pis] = random_bundle(1);
    Builder builder;
    ActionCircuitUltraPasta::build(builder, ws, pis);
    info("UltraHonk Pasta Action: ", builder.get_num_finalized_gates_inefficient(), " gates");
    EXPECT_FALSE(builder.failed()) << builder.err();
    EXPECT_TRUE(prove_and_verify(builder, pis));
}

TEST(ZcashActionUltraPasta, ProductionVectorsProveAndVerify)
{
    for (const auto* v : { &test_vectors::ACTION_REAL_SPEND,
                           &test_vectors::ACTION_DUMMY_SPEND_0,
                           &test_vectors::ACTION_DUMMY_SPEND_1 }) {
        const auto w = parse_action_witness<PastaCycle>(*v);
        const auto pis = vector_public_inputs(*v);
        Builder builder;
        ActionCircuitUltraPasta::build(builder, { w }, pis);
        EXPECT_FALSE(builder.failed()) << builder.err();
        EXPECT_TRUE(prove_and_verify(builder, pis));
    }
}

TEST(ZcashActionUltraPasta, WrongPublicInputFails)
{
    auto [ws, pis] = random_bundle(1);
    pis[halo2::layout::RK_X] += FF(1);
    Builder builder;
    ActionCircuitUltraPasta::build(builder, ws, pis);
    EXPECT_TRUE(builder.failed());
    EXPECT_FALSE(prove_and_verify(builder, pis));
}

TEST(ZcashActionUltraPasta, TwoActions)
{
    auto [ws, pis] = random_bundle(2);
    Builder builder;
    ActionCircuitUltraPasta::build(builder, ws, pis);
    EXPECT_FALSE(builder.failed()) << builder.err();
    EXPECT_TRUE(prove_and_verify(builder, pis));
}

// The canonicity gadget accepts p - 1 and rejects the 255-bit encoding p + 5 of 5.
TEST(ZcashActionUltraPasta, CanonicityCheck)
{
    const uint256_t p(FF::modulus);
    for (const auto& [encoding, canonical] : { std::make_pair(p - 1, true), std::make_pair(p + 5, false) }) {
        Builder builder;
        ActionCircuitUltraPasta circuit(builder);
        const uint32_t lo = circuit.witness(FF(encoding.slice(0, 130)));
        const uint32_t mid = circuit.witness(FF(encoding.slice(130, 254)));
        const uint32_t top = circuit.witness(FF(encoding.slice(254, 255)));
        circuit.range_constrain(lo, 130);
        circuit.range_constrain(mid, 124);
        circuit.assert_canonical(lo, mid, top, 130);
        EXPECT_EQ(builder.failed(), !canonical);
        EXPECT_EQ(prove_and_verify(builder, {}), canonical);
    }
}

TEST(ZcashActionUltraPasta, DISABLED_GateBreakdown)
{
    for (const size_t n : { size_t{ 1 }, size_t{ 2 } }) {
        auto [ws, pis] = random_bundle(n);
        Builder builder;
        ActionCircuitUltraPasta::build(builder, ws, pis);
        builder.finalize_circuit();
        info(n, " actions: ", builder.get_num_finalized_gates(), " gates");
        builder.blocks.summarize();
    }
}
