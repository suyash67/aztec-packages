#include "orchard_honk.hpp"
#include "barretenberg/zcash/halo2/action_circuit.hpp"
#include "barretenberg/zcash/primitives/orchard_test_utils.hpp"

#include <gtest/gtest.h>

using namespace bb;
using namespace bb::zcash;
using namespace bb::zcash::halo2;

namespace {
using Fp = PastaCycle::FF;

AnchoredTrace<PastaCycle> production_trace()
{
    const auto& v = test_vectors::ACTION_REAL_SPEND;
    std::vector<Fp> pis;
    for (const auto& s : v.public_inputs) {
        pis.emplace_back(uint256_t(s));
    }
    const auto w = parse_action_witness<PastaCycle>(v);
    auto table = ActionCircuit<PastaCycle>::build({ w }, pis);
    return AnchoredTrace<PastaCycle>::build(table, OrchardFlavor::TRACE_OFFSET);
}
} // namespace

TEST(ZcashOrchardHonk, ProveAndVerifyProductionVector)
{
    auto trace = production_trace();
    OrchardProvingKey pk(trace);
    auto proof = orchard_prove(pk, trace);
    info("proof size: ", proof.size() * 32, " bytes");
    EXPECT_TRUE(orchard_verify(*pk.vk, trace.public_inputs, proof));

    auto wrong_inputs = trace.public_inputs;
    wrong_inputs[RK_X] += Fp(1);
    EXPECT_FALSE(orchard_verify(*pk.vk, wrong_inputs, proof));
}

TEST(ZcashOrchardHonk, CorruptedGateWitnessIsRejected)
{
    auto trace = production_trace();
    OrchardProvingKey pk(trace);
    // lambda (a4) of a complete addition is not copy-constrained and not looked up: corrupting it breaks only q_add.
    size_t row = 0;
    for (size_t r = trace.row_offset; r < trace.num_rows; ++r) {
        if (!trace.selectors[Q_ADD][r].is_zero()) {
            row = r;
            break;
        }
    }
    ASSERT_NE(row, 0UL);
    trace.advice[4][row] += Fp(1);
    ASSERT_FALSE(trace.check().empty());
    auto proof = orchard_prove(pk, trace);
    EXPECT_FALSE(orchard_verify(*pk.vk, trace.public_inputs, proof));
}

TEST(ZcashOrchardHonk, TamperedOpeningProofIsRejected)
{
    auto trace = production_trace();
    OrchardProvingKey pk(trace);
    auto proof = orchard_prove(pk, trace);
    // The last two words are the IPA's final (c, f).
    auto tampered = proof;
    tampered[tampered.size() - 2] = uint256_t(Fp(tampered[tampered.size() - 2]) + Fp(1));
    EXPECT_FALSE(orchard_verify(*pk.vk, trace.public_inputs, tampered));
    tampered = proof;
    tampered.back() = uint256_t(Fp(tampered.back()) + Fp(1));
    EXPECT_FALSE(orchard_verify(*pk.vk, trace.public_inputs, tampered));
}

TEST(ZcashOrchardHonk, ProveAndVerifyMultipleActions)
{
    auto& engine = numeric::get_debug_randomness();
    std::vector<Orchard<PastaCycle>::ActionWitness> ws;
    std::vector<Fp> pis;
    for (size_t i = 0; i < 2; ++i) {
        ws.push_back(Orchard<PastaCycle>::random_witness(engine));
        for (const auto& f : Orchard<PastaCycle>::evaluate(ws.back())->to_field_elements()) {
            pis.push_back(f);
        }
    }
    auto table = ActionCircuit<PastaCycle>::build(ws, pis);
    auto trace = AnchoredTrace<PastaCycle>::build(table, OrchardFlavor::TRACE_OFFSET);
    OrchardProvingKey pk(trace);
    auto proof = orchard_prove(pk, trace);
    info("2 actions: n = ", trace.num_rows, ", proof size: ", proof.size() * 32, " bytes");
    EXPECT_TRUE(orchard_verify(*pk.vk, pis, proof));
}
