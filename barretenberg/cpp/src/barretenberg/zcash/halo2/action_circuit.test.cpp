#include "action_circuit.hpp"
#include "barretenberg/numeric/random/engine.hpp"
#include "barretenberg/zcash/primitives/orchard_test_utils.hpp"
#include "trace.hpp"

#include <gtest/gtest.h>

using namespace bb;
using namespace bb::zcash;
using namespace bb::zcash::halo2;

namespace {
using PastaAction = ActionCircuit<PastaCycle>;
using Fp = PastaCycle::FF;

std::vector<Fp> vector_public_inputs(const test_vectors::ActionVector& v)
{
    std::vector<Fp> out;
    for (const auto& s : v.public_inputs) {
        out.emplace_back(uint256_t(s));
    }
    return out;
}
} // namespace

TEST(ZcashHalo2Action, SynthesizesProductionVectors)
{
    for (const auto* v : { &test_vectors::ACTION_REAL_SPEND,
                           &test_vectors::ACTION_DUMMY_SPEND_0,
                           &test_vectors::ACTION_DUMMY_SPEND_1 }) {
        const auto w = parse_action_witness<PastaCycle>(*v);
        auto table = PastaAction::build({ w }, vector_public_inputs(*v));
        info("rows used by one Action: ", table.num_rows, ", copies: ", table.copies.size());
        EXPECT_LE(table.num_rows, 2048UL);
    }
}

TEST(ZcashHalo2Action, ProductionVectorSatisfiesAllConstraints)
{
    for (const auto* v : { &test_vectors::ACTION_REAL_SPEND, &test_vectors::ACTION_DUMMY_SPEND_0 }) {
        const auto w = parse_action_witness<PastaCycle>(*v);
        auto table = PastaAction::build({ w }, vector_public_inputs(*v));
        auto trace = AnchoredTrace<PastaCycle>::build(table, 8);
        EXPECT_EQ(trace.num_rows, 2048UL);
        auto failures = trace.check();
        for (const auto& f : failures) {
            info(f);
        }
        EXPECT_TRUE(failures.empty());
    }
}

TEST(ZcashHalo2Action, CorruptedTraceIsRejected)
{
    const auto w = parse_action_witness<PastaCycle>(test_vectors::ACTION_REAL_SPEND);
    auto table = PastaAction::build({ w }, vector_public_inputs(test_vectors::ACTION_REAL_SPEND));
    auto trace = AnchoredTrace<PastaCycle>::build(table, 8);
    ASSERT_TRUE(trace.check().empty());
    // Flip one advice cell in every column at a row used by some gate.
    for (size_t col = 0; col < NUM_ADVICE; ++col) {
        auto corrupted = trace;
        for (size_t row = trace.row_offset; row < trace.num_rows; ++row) {
            if (!corrupted.advice[col][row].is_zero()) {
                corrupted.advice[col][row] += Fp(1);
                break;
            }
        }
        EXPECT_FALSE(corrupted.check().empty()) << "column " << col;
    }
}

TEST(ZcashHalo2Action, MultipleRandomActions)
{
    auto& engine = numeric::get_debug_randomness();
    for (size_t num_actions : { 2UL, 4UL }) {
        std::vector<Orchard<PastaCycle>::ActionWitness> ws;
        std::vector<Fp> pis;
        for (size_t i = 0; i < num_actions; ++i) {
            ws.push_back(Orchard<PastaCycle>::random_witness(engine));
            auto pi = Orchard<PastaCycle>::evaluate(ws.back());
            ASSERT_TRUE(pi.has_value());
            for (const auto& f : pi->to_field_elements()) {
                pis.push_back(f);
            }
        }
        auto table = PastaAction::build(ws, pis);
        auto trace = AnchoredTrace<PastaCycle>::build(table, 8);
        info(num_actions, " actions: ", table.num_rows, " rows, trace size ", trace.num_rows);
        EXPECT_TRUE(trace.check().empty());
    }
}

TEST(ZcashHalo2Action, Bn254PortSatisfiesAllConstraints)
{
    using O = Orchard<Bn254Cycle>;
    auto& engine = numeric::get_debug_randomness();
    for (size_t num_actions : { 1UL, 2UL }) {
        std::vector<O::ActionWitness> ws;
        std::vector<fr> pis;
        for (size_t i = 0; i < num_actions; ++i) {
            ws.push_back(O::random_witness(engine));
            for (const auto& f : O::evaluate(ws.back())->to_field_elements()) {
                pis.push_back(f);
            }
        }
        auto table = ActionCircuit<Bn254Cycle>::build(ws, pis);
        auto trace = AnchoredTrace<Bn254Cycle>::build(table, 8);
        info("BN254 port, ", num_actions, " actions: ", table.num_rows, " rows, trace size ", trace.num_rows);
        auto failures = trace.check();
        for (const auto& f : failures) {
            info(f);
        }
        EXPECT_TRUE(failures.empty());
    }
}

TEST(ZcashHalo2Action, Bn254PortRejectsNonCanonicalDecomposition)
{
    // A corrupted canonicity witness (the t_hi cell of a CommitIvk / NoteCommit check) must break the gates.
    using O = Orchard<Bn254Cycle>;
    auto& engine = numeric::get_debug_randomness();
    auto w = O::random_witness(engine);
    std::vector<fr> pis;
    for (const auto& f : O::evaluate(w)->to_field_elements()) {
        pis.push_back(f);
    }
    auto table = ActionCircuit<Bn254Cycle>::build({ w }, pis);
    auto trace = AnchoredTrace<Bn254Cycle>::build(table, 8);
    ASSERT_TRUE(trace.check().empty());
    for (size_t s : { Q_NOTECOMMIT_G_D, Q_Y_CANON, Q_MUL_OVERFLOW }) {
        auto corrupted = trace;
        bool done = false;
        for (size_t row = 0; row < trace.num_rows && !done; ++row) {
            if (!corrupted.selectors[s][row].is_zero()) {
                // advice cells read by the gate at offsets 1..2 (anchor-relative)
                corrupted.advice[Q_Y_CANON == s ? 7 : 9][row + 1] += fr(1);
                corrupted.advice[8][row + 1] += fr(1);
                done = true;
            }
        }
        EXPECT_FALSE(corrupted.check().empty()) << "selector " << s;
    }
}
