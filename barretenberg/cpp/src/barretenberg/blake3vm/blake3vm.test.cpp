#include <chrono>
#include <gtest/gtest.h>

#include "barretenberg/blake3vm/blake3vm_circuit_builder.hpp"
#include "barretenberg/blake3vm/blake3vm_flavor.hpp"
#include "barretenberg/blake3vm/blake3vm_prover.hpp"
#include "barretenberg/blake3vm/blake3vm_verifier.hpp"
#include "barretenberg/honk/proof_system/logderivative_library.hpp"
#include "barretenberg/honk/relation_checker.hpp"
#include "barretenberg/srs/global_crs.hpp"

namespace bb {

class Blake3VMTests : public ::testing::Test {
  public:
    using Flavor = Blake3VMFlavor;
    using FF = Flavor::FF;

    static void SetUpTestSuite() { bb::srs::init_file_crs_factory(bb::srs::bb_crs_path()); }

    // A batch covering the WHIR profile (64-byte Merkle nodes) plus single-, two-, and
    // many-block messages to exercise CV chaining.
    static Blake3VMCircuitBuilder typical_builder()
    {
        Blake3VMCircuitBuilder builder;
        for (const size_t size : { 64UL, 64UL, 64UL, 32UL, 65UL, 769UL }) {
            std::vector<uint8_t> message(size);
            for (size_t j = 0; j < size; ++j) {
                message[j] = static_cast<uint8_t>(size * 31 + j * 7);
            }
            builder.add_hash(message);
        }
        return builder;
    }

    static RelationParameters<FF> random_parameters()
    {
        RelationParameters<FF> params;
        params.beta = FF::random_element();
        params.beta_sqr = params.beta * params.beta;
        params.gamma = FF::random_element();
        return params;
    }

    static void compute_all_inverses(Flavor::ProverPolynomials& polynomials, RelationParameters<FF>& params)
    {
        bb::constexpr_for<0, Blake3VMTraceLayout::NUM_LOOKUP_SETS, 1>([&]<size_t SET>() {
            compute_logderivative_inverse<FF, Blake3VMLookupRelation<FF, SET>, Flavor::ProverPolynomials, false>(
                polynomials, params, 0);
        });
    }
};

// The trace builder checks each digest against blake3::blake3s internally; here we pin the row
// accounting: one compression per 64-byte block, 56 G rows and one output block per compression.
TEST_F(Blake3VMTests, TraceGeneration)
{
    Blake3VMCircuitBuilder builder = typical_builder();
    // 64,64,64,32 bytes: 1 block each; 65 bytes: 2 blocks; 769 bytes: 13 blocks.
    constexpr size_t expected_compressions = 4 + 2 + 13;
    EXPECT_EQ(builder.num_compressions(), expected_compressions);
    EXPECT_EQ(builder.g_rows.size(), expected_compressions * Blake3VMTraceLayout::NUM_G_ROWS);
    EXPECT_EQ(builder.out_blocks.size(), expected_compressions);
    EXPECT_EQ(builder.boundaries.size(), expected_compressions);
}

// Every relation must vanish on an honestly generated trace.
TEST_F(Blake3VMTests, RelationCorrectness)
{
    Blake3VMCircuitBuilder builder = typical_builder();
    Flavor::ProverPolynomials polynomials(builder);
    RelationParameters<FF> params = random_parameters();
    compute_all_inverses(polynomials, params);

    using Base = RelationChecker<void>;
    const auto expect_no_failures = [](const auto& failures, const std::string& label) {
        for (const auto& [subrelation, row] : failures) {
            ADD_FAILURE() << label << ": subrelation " << subrelation << " first fails at row " << row;
        }
    };

    expect_no_failures(Base::check<Blake3VMGRelation<FF>>(polynomials, params, "G"), "G");
    expect_no_failures(Base::check<Blake3VMWiringRelation<FF>>(polynomials, params, "Wiring"), "Wiring");
    expect_no_failures(Base::check<Blake3VMZeroRowRelation<FF>>(polynomials, params, "ZeroRow"), "ZeroRow");
    bb::constexpr_for<0, Blake3VMTraceLayout::NUM_LOOKUP_SETS, 1>([&]<size_t SET>() {
        const std::string label = "Lookup" + std::to_string(SET);
        expect_no_failures(
            Base::check<Blake3VMLookupRelation<FF, SET>, /*has_linearly_dependent=*/true>(polynomials, params, label),
            label);
    });
}

TEST_F(Blake3VMTests, ProveAndVerify)
{
    Blake3VMCircuitBuilder builder = typical_builder();
    Blake3VMProver prover(builder);
    const HonkProof proof = prover.construct_proof();

    Blake3VMVerifier verifier(prover.verification_key);
    EXPECT_TRUE(verifier.verify_proof(proof));
}

// Manual perf check at a realistic batch size (~1,000 Merkle nodes -> 2^17 rows). Run with
// --gtest_also_run_disabled_tests.
TEST_F(Blake3VMTests, DISABLED_LargeBatchProveAndVerify)
{
    Blake3VMCircuitBuilder builder;
    std::array<uint8_t, 64> node;
    for (size_t i = 0; i < 1000; ++i) {
        for (size_t j = 0; j < node.size(); ++j) {
            node[j] = static_cast<uint8_t>(i * 13 + j);
        }
        builder.add_hash(node);
    }
    const auto start_prove = std::chrono::steady_clock::now();
    Blake3VMProver prover(builder);
    const HonkProof proof = prover.construct_proof();
    const auto end_prove = std::chrono::steady_clock::now();

    Blake3VMVerifier verifier(prover.verification_key);
    EXPECT_TRUE(verifier.verify_proof(proof));
    const auto end_verify = std::chrono::steady_clock::now();

    info("Blake3VM: 1000 compressions, circuit size ",
         prover.key->circuit_size,
         ", prove (incl. trace/keys) ",
         std::chrono::duration_cast<std::chrono::milliseconds>(end_prove - start_prove).count(),
         " ms, verify ",
         std::chrono::duration_cast<std::chrono::milliseconds>(end_verify - end_prove).count(),
         " ms, proof size ",
         proof.size() * 32,
         " bytes");
}

// Corrupting a single output word must make the proof unverifiable.
TEST_F(Blake3VMTests, TamperedWitnessFails)
{
    Blake3VMCircuitBuilder builder = typical_builder();
    Blake3VMProver prover(builder);
    const size_t tamper_row = Blake3VMTraceLayout::row_of_out(0, 0);
    prover.key->polynomials.out_0.at(tamper_row) += 1;
    const HonkProof proof = prover.construct_proof();

    Blake3VMVerifier verifier(prover.verification_key);
    EXPECT_FALSE(verifier.verify_proof(proof));
}

} // namespace bb
