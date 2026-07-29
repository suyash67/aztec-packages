#include "barretenberg/cq/cq_key.hpp"
#include "barretenberg/cq/cq_prover.hpp"
#include "barretenberg/cq/cq_trusted_setup.hpp"
#include "barretenberg/cq/cq_verifier.hpp"
#include "barretenberg/stdlib_circuit_builders/plookup_tables/sha256.hpp"
#include "barretenberg/stdlib_circuit_builders/plookup_tables/uint.hpp"

#include <gtest/gtest.h>

namespace bb::cq {

/**
 * @brief cq proofs against real UltraHonk plookup tables (BasicTable), demonstrating that the tables UltraHonk
 * currently materializes in the execution trace can instead be preprocessed once and looked up at cost independent
 * of the table size.
 */
TEST(CqPlookup, UintXorBasicTable)
{
    // The 6-bit XOR BasicTable used by UltraHonk uint32 operations: 4096 rows of (a, b, a ^ b).
    const plookup::BasicTable xor_table = plookup::uint_tables::generate_xor_rotate_table<6, 0>(
        plookup::BasicTableId::UINT_XOR_SLICE_6_ROTATE_0, /*table_index=*/0);
    const size_t table_size = xor_table.size();
    ASSERT_EQ(table_size, 4096);

    const TestSrs srs = TestSrs::create(table_size, table_size + 1);
    const CqProvingKey pk = CqProvingKey::create(
        { xor_table.column_1, xor_table.column_2, xor_table.column_3 }, srs.g1_powers, srs.g2_powers);
    const CommitmentKey<curve::BN254> ck = srs.create_commitment_key();

    // Witness: 32 XOR operations, i.e. 32 rows (a, b, a ^ b).
    const size_t n = 32;
    std::vector<std::vector<fr>> lookups(3, std::vector<fr>(n));
    for (size_t j = 0; j < n; ++j) {
        const uint64_t a = (j * 13 + 5) % 64;
        const uint64_t b = (j * 29 + 1) % 64;
        lookups[0][j] = a;
        lookups[1][j] = b;
        lookups[2][j] = a ^ b;
    }
    CqProver prover(pk, ck);
    const auto proof = prover.construct_proof(lookups);
    CqVerifier verifier(pk.verification_key);
    EXPECT_TRUE(verifier.verify_proof(proof));

    // A wrong XOR output is not a table row and must be rejected.
    auto bad_lookups = lookups;
    bad_lookups[2][7] = bad_lookups[2][7] + fr::one();
    const auto bad_proof = prover.construct_proof(bad_lookups, /*map_missing_rows_to_first_entry_for_testing=*/true);
    EXPECT_FALSE(verifier.verify_proof(bad_proof));
}

TEST(CqPlookup, Sha256MajorityNormalizationBasicTable)
{
    // The SHA-256 majority normalization BasicTable: 16^3 = 4096 rows mapping a base-16 sparse value to its
    // normalized binary form (columns 1 and 2; column 3 is unused by this table).
    const plookup::BasicTable maj_table =
        plookup::sha256_tables::generate_majority_normalization_table(plookup::BasicTableId::SHA256_MAJ_NORMALIZE,
                                                                      /*table_index=*/0);
    const size_t table_size = maj_table.size();
    ASSERT_EQ(table_size, 4096);

    const TestSrs srs = TestSrs::create(table_size, table_size + 1);
    const CqProvingKey pk =
        CqProvingKey::create({ maj_table.column_1, maj_table.column_2 }, srs.g1_powers, srs.g2_powers);
    const CommitmentKey<curve::BN254> ck = srs.create_commitment_key();

    // Witness: 64 normalization lookups sampled from the table rows.
    const size_t n = 64;
    std::vector<std::vector<fr>> lookups(2, std::vector<fr>(n));
    for (size_t j = 0; j < n; ++j) {
        const size_t row = (j * 389 + 17) % table_size;
        lookups[0][j] = maj_table.column_1[row];
        lookups[1][j] = maj_table.column_2[row];
    }
    CqProver prover(pk, ck);
    const auto proof = prover.construct_proof(lookups);
    CqVerifier verifier(pk.verification_key);
    EXPECT_TRUE(verifier.verify_proof(proof));

    // A wrong normalization result must be rejected.
    auto bad_lookups = lookups;
    bad_lookups[1][3] = bad_lookups[1][3] + fr::one();
    const auto bad_proof = prover.construct_proof(bad_lookups, /*map_missing_rows_to_first_entry_for_testing=*/true);
    EXPECT_FALSE(verifier.verify_proof(bad_proof));
}

} // namespace bb::cq
