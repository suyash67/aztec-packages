#include "action_circuit_ultra.hpp"
#include "barretenberg/circuit_checker/circuit_checker.hpp"
#include "barretenberg/flavor/ultra_zk_flavor.hpp"
#include "barretenberg/numeric/random/engine.hpp"
#include "barretenberg/srs/global_crs.hpp"
#include "barretenberg/ultra_honk/prover_instance.hpp"
#include "barretenberg/ultra_honk/ultra_prover.hpp"
#include "barretenberg/ultra_honk/ultra_verifier.hpp"

#include <gtest/gtest.h>

using namespace bb;
using namespace bb::zcash;
using namespace bb::zcash::ultra;

namespace {
using O = Orchard<Bn254Cycle>;

std::pair<std::vector<O::ActionWitness>, std::vector<fr>> random_bundle(size_t n)
{
    auto& engine = numeric::get_debug_randomness();
    std::vector<O::ActionWitness> ws;
    std::vector<fr> pis;
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
} // namespace

TEST(ZcashActionUltra, RandomActionSatisfiesCircuit)
{
    auto [ws, pis] = random_bundle(1);
    Builder builder;
    ActionCircuitUltra::build(builder, ws, pis);
    info("UltraHonk BN254 Action: ", builder.get_num_finalized_gates_inefficient(), " gates");
    EXPECT_TRUE(CircuitChecker::check(builder));
}

TEST(ZcashActionUltra, WrongPublicInputFails)
{
    auto [ws, pis] = random_bundle(1);
    pis[halo2::layout::RK_X] += fr(1);
    Builder builder;
    ActionCircuitUltra::build(builder, ws, pis);
    EXPECT_TRUE(builder.failed() || !CircuitChecker::check(builder));
}

TEST(ZcashActionUltra, ProveAndVerifyUltraZK)
{
    using Flavor = UltraZKFlavor;
    bb::srs::init_bn254_file_crs_factory(bb::srs::bb_crs_path());
    auto [ws, pis] = random_bundle(1);
    Builder builder;
    ActionCircuitUltra::build(builder, ws, pis);
    DefaultIO::add_default(builder);
    auto prover_instance = std::make_shared<ProverInstance_<Flavor>>(builder);
    auto vk = std::make_shared<Flavor::VerificationKey>(prover_instance->get_precomputed());
    auto vk_and_hash = std::make_shared<Flavor::VKAndHash>(vk);
    UltraProver_<Flavor> prover(prover_instance, vk);
    auto proof = prover.construct_proof();
    info("UltraZK BN254 Action: log_n = ", vk->log_circuit_size, ", proof size: ", proof.size() * 32, " bytes");
    // The Action's public inputs lead the proof.
    for (size_t i = 0; i < pis.size(); ++i) {
        EXPECT_EQ(proof[i], pis[i]);
    }
    UltraVerifier_<Flavor, DefaultIO> verifier(vk_and_hash);
    EXPECT_TRUE(verifier.verify_proof(proof).result);
}
