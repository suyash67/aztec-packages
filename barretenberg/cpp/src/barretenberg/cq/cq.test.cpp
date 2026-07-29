#include "barretenberg/common/log.hpp"
#include "barretenberg/cq/cq_key.hpp"
#include "barretenberg/cq/cq_prover.hpp"
#include "barretenberg/cq/cq_trusted_setup.hpp"
#include "barretenberg/cq/cq_verifier.hpp"
#include "barretenberg/ecc/curves/bn254/pairing.hpp"
#include "barretenberg/ecc/scalar_multiplication/scalar_multiplication.hpp"
#include "barretenberg/polynomials/polynomial_arithmetic.hpp"

#include <chrono>
#include <gtest/gtest.h>

namespace bb::cq {

namespace {

std::vector<fr> random_table(size_t size)
{
    std::vector<fr> table(size);
    for (auto& t : table) {
        t = fr::random_element();
    }
    return table;
}

std::vector<fr> sample_lookups(const std::vector<fr>& table, size_t num_lookups)
{
    std::vector<fr> lookups(num_lookups);
    for (size_t j = 0; j < num_lookups; ++j) {
        // Deterministic but scattered picks, with repeats.
        lookups[j] = table[(j * 7919 + 13) % table.size()];
    }
    return lookups;
}

} // namespace

class CqTest : public ::testing::Test {
  public:
    static constexpr size_t N = 256;
    static TestSrs srs;
    static std::vector<fr> table;
    static CqProvingKey pk;
    static CommitmentKey<curve::BN254> ck;

    static void SetUpTestSuite()
    {
        srs = TestSrs::create(N, N + 1);
        table = random_table(N);
        pk = CqProvingKey::create(table, srs.g1_powers, srs.g2_powers);
        ck = srs.create_commitment_key();
    }
};

TestSrs CqTest::srs;
std::vector<fr> CqTest::table;
CqProvingKey CqTest::pk;
CommitmentKey<curve::BN254> CqTest::ck;

TEST_F(CqTest, PreprocessedTableCommitmentsAreConsistent)
{
    // [T]_1 and [T]_2 must commit to the same polynomial: e([T]_1, [1]_2) == e([1]_1, [T]_2).
    const g1::affine_element g1_one(g1::one);
    const std::array<g1::affine_element, 2> p1{ pk.verification_key.table_commitments_g1[0],
                                                g1::affine_element(-g1::element(g1_one)) };
    const std::array<g2::affine_element, 2> p2{ srs.g2_powers[0], pk.verification_key.table_commitments_g2[0] };
    EXPECT_EQ(pairing::reduced_ate_pairing_batch(p1.data(), p2.data(), 2), fq12::one());

    // T interpolates the table over H.
    const SubgroupDomain domain(N);
    const std::vector<fr> roots = domain.root_powers();
    for (size_t i = 0; i < N; i += 17) {
        EXPECT_EQ(polynomial_arithmetic::evaluate<fr>(pk.column_monomials[0], roots[i]), table[i]);
    }
}

TEST_F(CqTest, CachedQuotientsSatisfyDefiningIdentity)
{
    // Q_i = L_i (T - t_i) / Z_H: recompute a few naively and compare commitments.
    const SubgroupDomain domain(N);
    const SubgroupDomain double_domain(2 * N);
    for (size_t i = 0; i < N; i += 41) {
        std::vector<fr> lagrange_coeffs(N, fr::zero());
        lagrange_coeffs[i] = fr::one();
        domain.ifft<fr>(lagrange_coeffs);
        std::vector<fr> shifted_table = pk.column_monomials[0];
        shifted_table[0] -= table[i];
        // product of degree-(N-1) polynomials via a 2N FFT
        std::vector<fr> a_padded(2 * N, fr::zero());
        std::vector<fr> b_padded(2 * N, fr::zero());
        std::copy(lagrange_coeffs.begin(), lagrange_coeffs.end(), a_padded.begin());
        std::copy(shifted_table.begin(), shifted_table.end(), b_padded.begin());
        double_domain.fft<fr>(a_padded);
        double_domain.fft<fr>(b_padded);
        for (size_t k = 0; k < 2 * N; ++k) {
            a_padded[k] *= b_padded[k];
        }
        double_domain.ifft<fr>(a_padded);
        // divide by Z_H = X^N - 1: quotient is the top block
        std::vector<fr> quotient(a_padded.begin() + N, a_padded.begin() + static_cast<std::ptrdiff_t>(2 * N - 1));
        for (size_t k = 0; k + 1 < N; ++k) {
            ASSERT_EQ(a_padded[k] + a_padded[k + N], fr::zero());
        }
        ASSERT_EQ(a_padded[N - 1], fr::zero());
        const g1::affine_element expected(scalar_multiplication::pippenger<curve::BN254>(
            PolynomialSpan<const fr>{ 0, quotient }, srs.g1_powers, /*handle_edge_cases=*/true));
        EXPECT_EQ(pk.cached_quotients[0][i], expected) << "cached quotient " << i;
    }
}

TEST_F(CqTest, ProveAndVerify)
{
    const std::vector<fr> lookups = sample_lookups(table, 16);
    CqProver prover(pk, ck);
    const auto proof = prover.construct_proof(lookups);
    CqVerifier verifier(pk.verification_key);
    EXPECT_TRUE(verifier.verify_proof(proof));
}

TEST_F(CqTest, ProveAndVerifyMaximalLookups)
{
    // n = N: every table entry looked up (and then some repeats).
    std::vector<fr> lookups(N);
    for (size_t j = 0; j < N; ++j) {
        lookups[j] = table[(j * 3) % N];
    }
    CqProver prover(pk, ck);
    const auto proof = prover.construct_proof(lookups);
    CqVerifier verifier(pk.verification_key);
    EXPECT_TRUE(verifier.verify_proof(proof));
}

TEST_F(CqTest, ProveAndVerifySingleRepeatedValue)
{
    const std::vector<fr> lookups(8, table[5]);
    CqProver prover(pk, ck);
    const auto proof = prover.construct_proof(lookups);
    CqVerifier verifier(pk.verification_key);
    EXPECT_TRUE(verifier.verify_proof(proof));
}

TEST_F(CqTest, ValueOutsideTableRejected)
{
    std::vector<fr> lookups = sample_lookups(table, 16);
    lookups[3] = fr::random_element(); // not in the table (w.h.p.)
    CqProver prover(pk, ck);
    const auto proof = prover.construct_proof(lookups, /*map_missing_values_to_first_entry_for_testing=*/true);
    CqVerifier verifier(pk.verification_key);
    EXPECT_FALSE(verifier.verify_proof(proof));
}

TEST_F(CqTest, TamperedEvaluationsRejected)
{
    const std::vector<fr> lookups = sample_lookups(table, 16);
    CqProver prover(pk, ck);
    const auto proof = prover.construct_proof(lookups);
    CqVerifier verifier(pk.verification_key);
    ASSERT_TRUE(verifier.verify_proof(proof));

    // The three evaluations (b0_gamma, f_gamma, a0) sit right before the two final opening proofs, each of which
    // serializes to 4 field elements.
    const size_t first_evaluation_index = proof.size() - 11;
    for (size_t idx : { first_evaluation_index, first_evaluation_index + 1, first_evaluation_index + 2 }) {
        auto tampered = proof;
        tampered[idx] += fr::one();
        EXPECT_FALSE(verifier.verify_proof(tampered)) << "tampered index " << idx;
    }
}

TEST_F(CqTest, ProofForDifferentTableRejected)
{
    // A valid proof for a different table of the same size must be rejected via the table commitment binding.
    std::vector<fr> other_table = random_table(N);
    CqProvingKey other_pk = CqProvingKey::create(other_table, srs.g1_powers, srs.g2_powers);
    CqProver prover(other_pk, ck);
    const auto proof = prover.construct_proof(sample_lookups(other_table, 8));
    EXPECT_TRUE(CqVerifier(other_pk.verification_key).verify_proof(proof));
    EXPECT_FALSE(CqVerifier(pk.verification_key).verify_proof(proof));
}

TEST(CqScaling, ProverTimeIndependentOfTableSize)
{
    // Same number of lookups against tables of increasing size: proving time must stay flat (the whole point of
    // cq). Not asserted (timing asserts are flaky); logged for inspection alongside a correctness check.
    const size_t num_lookups = 64;
    for (const size_t table_size : { 1UL << 10, 1UL << 12 }) {
        TestSrs srs = TestSrs::create(table_size, table_size + 1);
        std::vector<fr> table = random_table(table_size);
        const auto preprocess_start = std::chrono::steady_clock::now();
        CqProvingKey pk = CqProvingKey::create(table, srs.g1_powers, srs.g2_powers);
        const auto preprocess_end = std::chrono::steady_clock::now();
        CommitmentKey<curve::BN254> ck = srs.create_commitment_key();

        const std::vector<fr> lookups = sample_lookups(table, num_lookups);
        CqProver prover(pk, ck);
        const auto prove_start = std::chrono::steady_clock::now();
        const auto proof = prover.construct_proof(lookups);
        const auto prove_end = std::chrono::steady_clock::now();
        EXPECT_TRUE(CqVerifier(pk.verification_key).verify_proof(proof));

        const auto to_ms = [](auto duration) {
            return std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
        };
        info("cq: table size 2^",
             numeric::get_msb(table_size),
             ": preprocessing ",
             to_ms(preprocess_end - preprocess_start),
             " ms, proving ",
             to_ms(prove_end - prove_start),
             " ms (",
             num_lookups,
             " lookups)");
    }
}

} // namespace bb::cq
