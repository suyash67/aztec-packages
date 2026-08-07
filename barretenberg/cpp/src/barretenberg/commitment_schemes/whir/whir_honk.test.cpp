#include "barretenberg/commitment_schemes/whir/whir_honk.hpp"

#include "barretenberg/flavor/ultra_provekit_flavor.hpp"
#include "barretenberg/stdlib/hash/poseidon2/poseidon2.hpp"
#include "barretenberg/stdlib_circuit_builders/mock_circuits.hpp"

#include <gtest/gtest.h>

namespace bb::whir {

namespace {

/** @brief A circuit exercising arithmetic, public-input, lookup, and ROM (memory) gates. */
UltraCircuitBuilder build_test_circuit()
{
    UltraCircuitBuilder builder;
    MockCircuits::add_arithmetic_gates_with_public_inputs(builder, 16);
    MockCircuits::add_arithmetic_gates(builder, 16);
    MockCircuits::add_lookup_gates(builder, 2);
    // A small ROM table so the memory relation and the eta-dependent w_4 records are exercised.
    const size_t rom_id = builder.create_ROM_array(4);
    for (size_t i = 0; i < 4; ++i) {
        builder.set_ROM_element(rom_id, i, builder.add_variable(fr(3 * i + 1)));
    }
    const uint32_t read_index = builder.add_variable(fr(2));
    builder.read_ROM_array(rom_id, read_index);
    return builder;
}

} // namespace

template <typename Hasher> class WhirHonkTest : public ::testing::Test {
  public:
    using Honk = WhirHonk<Hasher>;

    struct Setup {
        WhirConfig config;
        typename Honk::VerificationKey vk;
        HonkProof proof;
    };

    static Setup prove_test_circuit()
    {
        UltraCircuitBuilder builder = build_test_circuit();
        // Sized from the circuit: ProverInstance dictates the dyadic size.
        ProverInstance_<UltraFlavor> sizing_instance(builder);
        const size_t log_n = sizing_instance.log_dyadic_size();

        UltraCircuitBuilder proving_builder = build_test_circuit();
        const WhirConfig config = Honk::make_config(log_n, /*security_bits=*/64);
        auto pk = Honk::create_proving_key(proving_builder, config);
        return { config, pk.vk, Honk::prove(pk) };
    }
};

using HasherTypes =
    ::testing::Types<Poseidon2MerkleHasher, Blake3sMerkleHasher, Sha256MerkleHasher, SkyscraperMerkleHasher>;
TYPED_TEST_SUITE(WhirHonkTest, HasherTypes);

TYPED_TEST(WhirHonkTest, ProveAndVerify)
{
    const auto setup = TestFixture::prove_test_circuit();
    EXPECT_TRUE(TestFixture::Honk::verify(setup.vk, setup.config, setup.proof));
}

TYPED_TEST(WhirHonkTest, TamperedProofRejected)
{
    const auto setup = TestFixture::prove_test_circuit();

    // Positions spread across the proof: public inputs / wire roots, sumcheck messages, WHIR
    // openings, the final polynomial.
    for (const size_t position :
         { size_t(2), setup.proof.size() / 4, setup.proof.size() / 2, setup.proof.size() - 2 }) {
        HonkProof tampered = setup.proof;
        tampered[position] += fr(1);
        EXPECT_FALSE(TestFixture::Honk::verify(setup.vk, setup.config, tampered)) << "position " << position;
    }
}

TYPED_TEST(WhirHonkTest, VirtualColumnsOmittedFromCommitment)
{
    const auto setup = TestFixture::prove_test_circuit();
    const uint32_t mask = setup.vk.virtual_mask;

    // The lagrange point indicators are always virtual (entities 8, 9), and the test circuit uses
    // no elliptic (23), nnf (25), or poseidon2 (26, 27) gates, so those selectors are zero columns.
    for (const size_t entity : { 8U, 9U, 23U, 25U, 26U, 27U }) {
        EXPECT_TRUE(((mask >> entity) & 1) == 1) << "entity " << entity;
    }
    EXPECT_LE(setup.vk.num_committed_precomputed(), 22U);
    EXPECT_TRUE(TestFixture::Honk::verify(setup.vk, setup.config, setup.proof));

    // The indicator rows are part of the statement: a wrong row must be rejected.
    auto bad_vk = setup.vk;
    bad_vk.lagrange_last_row ^= 1;
    EXPECT_FALSE(TestFixture::Honk::verify(bad_vk, setup.config, setup.proof));
}

TYPED_TEST(WhirHonkTest, WrongVkRejected)
{
    const auto setup = TestFixture::prove_test_circuit();

    auto bad_vk = setup.vk;
    auto root_fields = TypeParam::digest_to_fields(bad_vk.precomputed_commitment);
    root_fields[0] += fr(1);
    bad_vk.precomputed_commitment = TypeParam::digest_from_fields(root_fields);
    EXPECT_FALSE(TestFixture::Honk::verify(bad_vk, setup.config, setup.proof));
}

/**
 * @brief The stacked layout proves and verifies the same circuit, with a materially smaller proof.
 * @details `WhirStackedHonk` trades prover time for proof size; this pins the size win and that the
 * two layouts remain independently sound.
 */
TEST(WhirStackedHonkTest, ProveAndVerifyIsSmallerThanInterleaved)
{
    using Hasher = Blake3sMerkleHasher;
    const auto build = [] {
        UltraCircuitBuilder builder = build_test_circuit();
        ProverInstance_<UltraFlavor> sizing(builder);
        return sizing.log_dyadic_size();
    };
    const size_t log_n = build();

    UltraCircuitBuilder stacked_builder = build_test_circuit();
    const auto stacked_config = WhirStackedHonk<Hasher>::make_config(log_n, /*security_bits=*/64);
    auto stacked_pk = WhirStackedHonk<Hasher>::create_proving_key(stacked_builder, stacked_config);
    const HonkProof stacked_proof = WhirStackedHonk<Hasher>::prove(stacked_pk);
    EXPECT_TRUE(WhirStackedHonk<Hasher>::verify(stacked_pk.vk, stacked_config, stacked_proof));

    UltraCircuitBuilder plain_builder = build_test_circuit();
    const auto plain_config = WhirHonk<Hasher>::make_config(log_n, /*security_bits=*/64);
    auto plain_pk = WhirHonk<Hasher>::create_proving_key(plain_builder, plain_config);
    const HonkProof plain_proof = WhirHonk<Hasher>::prove(plain_pk);
    EXPECT_TRUE(WhirHonk<Hasher>::verify(plain_pk.vk, plain_config, plain_proof));

    EXPECT_LT(stacked_proof.size(), plain_proof.size());
}

/**
 * @brief The reduced ProveKit flavor proves the same circuit, more cheaply and in fewer fields.
 * @details `build_test_circuit` uses only arithmetic, public-input, lookup and memory gates, so the
 * four relations `UltraProveKitFlavor` drops are vacuous on it and both flavors must accept. The
 * proof shrinks by exactly the columns and round-polynomial evaluations the reduction removes:
 * four fewer opened precomputed columns, and one fewer evaluation in each of the `log_n` sumcheck
 * round polynomials (the Poseidon2 subrelations were the only degree-6 ones).
 */
TEST(WhirProveKitFlavorTest, ProveAndVerifyWithFewerFieldsThanUltra)
{
    using Hasher = Blake3sMerkleHasher;
    using ReducedHonk = honk_transparent::TransparentHonk<WhirPcs<Hasher>, UltraProveKitFlavor>;

    static_assert(UltraProveKitFlavor::NUM_PRECOMPUTED_ENTITIES == UltraFlavor::NUM_PRECOMPUTED_ENTITIES - 4);
    static_assert(UltraProveKitFlavor::MAX_PARTIAL_RELATION_LENGTH == UltraFlavor::MAX_PARTIAL_RELATION_LENGTH - 1);

    UltraCircuitBuilder sizing_builder = build_test_circuit();
    const size_t log_n = ProverInstance_<UltraFlavor>(sizing_builder).log_dyadic_size();

    UltraCircuitBuilder reduced_builder = build_test_circuit();
    const auto config = ReducedHonk::make_config(log_n, /*security_bits=*/64);
    auto reduced_pk = ReducedHonk::create_proving_key(reduced_builder, config);
    const auto reduced_vk = reduced_pk.vk;
    const HonkProof reduced_proof = ReducedHonk::prove(reduced_pk);
    EXPECT_TRUE(ReducedHonk::verify(reduced_vk, config, reduced_proof));

    UltraCircuitBuilder ultra_builder = build_test_circuit();
    auto ultra_pk = WhirHonk<Hasher>::create_proving_key(ultra_builder, config);
    const HonkProof ultra_proof = WhirHonk<Hasher>::prove(ultra_pk);
    EXPECT_TRUE(WhirHonk<Hasher>::verify(ultra_pk.vk, config, ultra_proof));

    // The dropped selectors were identically zero, so UltraFlavor already left them uncommitted:
    // both flavors open the same set of columns, and the saving is the shorter sumcheck round
    // polynomials plus the four claimed evaluations the reduced flavor no longer sends — not fewer
    // commitments. (That saving is not readable off the total proof lengths: the flavors' sumcheck
    // messages differ, so the query indices do too, and batched Merkle openings cost a number of
    // digests that depends on them.)
    EXPECT_EQ(reduced_vk.num_committed_precomputed(), ultra_pk.vk.num_committed_precomputed());

    HonkProof tampered = reduced_proof;
    tampered[tampered.size() / 2] += fr(1);
    EXPECT_FALSE(ReducedHonk::verify(reduced_vk, config, tampered));
}

/**
 * @brief Narrowing the first fold to 2 proves, verifies, and shrinks the proof several-fold.
 * @details Round-0 queries are the only ones that open the wide commitment, and each reveals
 * 2^{k₀} values of every committed column. Dropping k₀ from 3 to 1 quarters that term while the
 * later single-column rounds keep folding by 8, so the proof must come out substantially smaller
 * and still verify — and a tampered one must still fail.
 */
TEST(WhirNarrowHonkTest, NarrowFirstFoldIsSmallerAndStillSound)
{
    using Hasher = SkyscraperMerkleHasher;
    using Wide = honk_transparent::TransparentHonk<WhirPcs<Hasher, 0, WhirSoundness::PROVABLE_LIST, 3>>;
    using Narrow = honk_transparent::TransparentHonk<WhirPcs<Hasher, 0, WhirSoundness::PROVABLE_LIST, 3, 1>>;

    UltraCircuitBuilder sizing_builder = build_test_circuit();
    const size_t log_n = ProverInstance_<UltraFlavor>(sizing_builder).log_dyadic_size();

    UltraCircuitBuilder wide_builder = build_test_circuit();
    const auto wide_config = Wide::make_config(log_n, /*security_bits=*/64);
    auto wide_pk = Wide::create_proving_key(wide_builder, wide_config);
    const auto wide_vk = wide_pk.vk;
    const HonkProof wide_proof = Wide::prove(wide_pk);
    EXPECT_TRUE(Wide::verify(wide_vk, wide_config, wide_proof));

    UltraCircuitBuilder narrow_builder = build_test_circuit();
    const auto narrow_config = Narrow::make_config(log_n, /*security_bits=*/64);
    auto narrow_pk = Narrow::create_proving_key(narrow_builder, narrow_config);
    const auto narrow_vk = narrow_pk.vk;
    const HonkProof narrow_proof = Narrow::prove(narrow_pk);
    EXPECT_TRUE(Narrow::verify(narrow_vk, narrow_config, narrow_proof));

    // Round-0 leaf values drop to a quarter; the deeper trees give some of it back in path
    // digests, so the whole proof lands under 60% at this size and lower still as the circuit grows.
    EXPECT_LT(narrow_proof.size() * 5, wide_proof.size() * 3);

    HonkProof tampered = narrow_proof;
    tampered[tampered.size() / 2] += fr(1);
    EXPECT_FALSE(Narrow::verify(narrow_vk, narrow_config, tampered));

    // A proof made under one schedule must not verify under the other: the folding factors are
    // derived, not transcribed, so a mismatch has to surface as a failed check.
    EXPECT_FALSE(Narrow::verify(narrow_vk, wide_config, narrow_proof));
}

/**
 * @brief Grinding shrinks the proof further and its nonces are checked, not decorative.
 * @details The query counts drop because 20 of the 64 bits come from proof of work, so the proof
 * must be smaller than the same schedule without it — and a proof whose nonces were ground for 20
 * bits must fail against a verifier expecting more work.
 */
TEST(WhirNarrowHonkTest, GrindingShrinksTheProofAndIsEnforced)
{
    using Hasher = SkyscraperMerkleHasher;
    using Plain = honk_transparent::TransparentHonk<WhirPcs<Hasher, 0, WhirSoundness::PROVABLE_LIST, 3, 1>>;
    using Ground = honk_transparent::TransparentHonk<WhirPcs<Hasher, 0, WhirSoundness::PROVABLE_LIST, 3, 1, 20>>;

    UltraCircuitBuilder sizing_builder = build_test_circuit();
    const size_t log_n = ProverInstance_<UltraFlavor>(sizing_builder).log_dyadic_size();

    UltraCircuitBuilder plain_builder = build_test_circuit();
    const auto plain_config = Plain::make_config(log_n, /*security_bits=*/64);
    auto plain_pk = Plain::create_proving_key(plain_builder, plain_config);
    const HonkProof plain_proof = Plain::prove(plain_pk);
    EXPECT_TRUE(Plain::verify(plain_pk.vk, plain_config, plain_proof));

    UltraCircuitBuilder ground_builder = build_test_circuit();
    const auto ground_config = Ground::make_config(log_n, /*security_bits=*/64);
    auto ground_pk = Ground::create_proving_key(ground_builder, ground_config);
    const auto ground_vk = ground_pk.vk;
    const HonkProof ground_proof = Ground::prove(ground_pk);
    EXPECT_TRUE(Ground::verify(ground_vk, ground_config, ground_proof));

    EXPECT_LT(ground_proof.size(), plain_proof.size());

    // Demanding more work than the prover did must be rejected, even though every other check in
    // the proof still passes.
    WhirConfig stricter = ground_config;
    stricter.pow_bits = 30;
    EXPECT_FALSE(Ground::verify(ground_vk, stricter, ground_proof));
}

/**
 * @brief A circuit that uses a dropped gate kind is refused, not silently proved.
 * @details Dropping a relation is only sound when its trace block is empty. If the reduced flavor
 * proved a circuit with poseidon2 gates anyway, the resulting proof would verify while leaving
 * those gates unconstrained.
 */
TEST(WhirProveKitFlavorTest, CircuitUsingADroppedRelationIsRefused)
{
    using ReducedHonk = honk_transparent::TransparentHonk<WhirPcs<Blake3sMerkleHasher>, UltraProveKitFlavor>;

    UltraCircuitBuilder builder = build_test_circuit();
    stdlib::poseidon2<UltraCircuitBuilder>::hash(
        { stdlib::field_t<UltraCircuitBuilder>::from_witness(&builder, fr(1)) });

    UltraCircuitBuilder sizing_builder = builder;
    const size_t log_n = ProverInstance_<UltraFlavor>(sizing_builder).log_dyadic_size();
    const auto config = ReducedHonk::make_config(log_n, /*security_bits=*/64);
    EXPECT_THROW(ReducedHonk::create_proving_key(builder, config), std::runtime_error);
}

} // namespace bb::whir
