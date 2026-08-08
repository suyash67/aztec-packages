#include "barretenberg/srs/global_crs.hpp"
#include "barretenberg/stdlib/hash/poseidon2/poseidon2.hpp"
#include "barretenberg/stdlib/primitives/bigfield/bigfield.hpp"
#include "barretenberg/stdlib_circuit_builders/plookup_tables/fixed_base/fixed_base.hpp"
#include "barretenberg/ultra_fflonk/prover.hpp"
#include "barretenberg/ultra_fflonk/verifier.hpp"

#include <gtest/gtest.h>

using namespace bb;
using namespace bb::ultra_fflonk;

namespace {

using Builder = Flavor::CircuitBuilder;

/** @brief A small circuit: public inputs, a chain of arithmetic gates, and a copy constraint. */
Builder arithmetic_circuit(size_t num_public_inputs, size_t num_gates)
{
    Builder builder;

    std::vector<uint32_t> publics;
    publics.reserve(num_public_inputs);
    for (size_t i = 0; i < num_public_inputs; ++i) {
        publics.push_back(builder.add_public_variable(FF(static_cast<uint64_t>(i) + 3)));
    }

    const uint32_t a_idx = publics.empty() ? builder.add_variable(FF(7)) : publics[0];
    const FF a = builder.get_variable(a_idx);
    const FF b = FF::random_element();
    const uint32_t b_idx = builder.add_variable(b);
    const uint32_t c_idx = builder.add_variable(a + b);
    for (size_t i = 0; i < num_gates; ++i) {
        builder.create_add_gate({ a_idx, b_idx, c_idx, 1, 1, -1, 0 });
    }
    return builder;
}

/**
 * @brief A circuit that lights up every gate selector the flavor has.
 *
 * @details The point of consuming `UltraCircuitBuilder` unchanged is that lookups, RAM/ROM, elliptic
 * additions, range constraints, non-native field arithmetic and Poseidon2 all come along with it -
 * so the tests have to actually exercise them, or "it proves" would only be a statement about
 * addition gates.
 */
Builder assorted_circuit()
{
    Builder builder;

    // Arithmetic, including a big add gate that reaches into the next row (q_arith = 2).
    const FF a = FF::random_element();
    const uint32_t a_idx = builder.add_public_variable(a);
    const FF b = FF::random_element();
    const FF c = a + b;
    const FF d = a + c;
    const uint32_t b_idx = builder.add_variable(b);
    const uint32_t c_idx = builder.add_variable(c);
    const uint32_t d_idx = builder.add_variable(d);
    for (size_t i = 0; i < 8; ++i) {
        builder.create_add_gate({ a_idx, b_idx, c_idx, 1, 1, -1, 0 });
        builder.create_add_gate({ d_idx, c_idx, a_idx, 1, -1, -1, 0 });
    }
    const uint32_t e_idx = builder.add_variable(a + b + c + d);
    builder.create_big_add_gate({ a_idx, b_idx, c_idx, d_idx, -1, -1, -1, -1, 0 }, true);
    builder.create_big_add_gate({ builder.zero_idx(), builder.zero_idx(), builder.zero_idx(), e_idx, 0, 0, 0, 0, 0 },
                                false);

    // Lookups, via the fixed-base pedersen tables.
    const FF pedersen_input = FF::random_element();
    const auto input_hi =
        uint256_t(pedersen_input)
            .slice(plookup::fixed_base::table::BITS_PER_LO_SCALAR,
                   plookup::fixed_base::table::BITS_PER_LO_SCALAR + plookup::fixed_base::table::BITS_PER_HI_SCALAR);
    const auto input_lo = uint256_t(pedersen_input).slice(0, plookup::fixed_base::table::BITS_PER_LO_SCALAR);
    const uint32_t hi_idx = builder.add_variable(FF(input_hi));
    const uint32_t lo_idx = builder.add_variable(FF(input_lo));
    builder.create_gates_from_plookup_accumulators(
        plookup::MultiTableId::FIXED_BASE_LEFT_HI,
        plookup::get_lookup_accumulators(plookup::MultiTableId::FIXED_BASE_LEFT_HI, input_hi),
        hi_idx);
    builder.create_gates_from_plookup_accumulators(
        plookup::MultiTableId::FIXED_BASE_LEFT_LO,
        plookup::get_lookup_accumulators(plookup::MultiTableId::FIXED_BASE_LEFT_LO, input_lo),
        lo_idx);

    // Delta range.
    builder.enforce_small_deltas({ builder.add_variable(FF(0)),
                                   builder.add_variable(FF(1)),
                                   builder.add_variable(FF(2)),
                                   builder.add_variable(FF(3)) });

    // Elliptic curve addition.
    const grumpkin::g1::affine_element p1 = grumpkin::g1::affine_element::random_element();
    const grumpkin::g1::affine_element p2 = grumpkin::g1::affine_element::random_element();
    const grumpkin::g1::affine_element p3(grumpkin::g1::element(p1) - grumpkin::g1::element(p2));
    builder.create_ecc_add_gate({ builder.add_variable(p1.x),
                                  builder.add_variable(p1.y),
                                  builder.add_variable(p2.x),
                                  builder.add_variable(p2.y),
                                  builder.add_variable(p3.x),
                                  builder.add_variable(p3.y),
                                  /*is_addition=*/false });

    // RAM, which also exercises the ROM-LogUp running sum.
    std::array<uint32_t, 8> ram_values{};
    for (uint32_t& value : ram_values) {
        value = builder.add_variable(FF::random_element());
    }
    const size_t ram_id = builder.create_RAM_array(ram_values.size());
    for (size_t i = 0; i < ram_values.size(); ++i) {
        builder.init_RAM_element(ram_id, i, ram_values[i]);
    }
    const uint32_t read_a = builder.read_RAM_array(ram_id, builder.add_variable(FF(5)));
    const uint32_t read_b = builder.read_RAM_array(ram_id, builder.add_variable(FF(1)));
    builder.write_RAM_array(ram_id, builder.add_variable(FF(4)), builder.add_variable(FF(500)));
    const uint32_t read_c = builder.read_RAM_array(ram_id, builder.add_variable(FF(4)));
    const uint32_t ram_sum = builder.add_variable(builder.get_variable(read_a) + builder.get_variable(read_b) +
                                                  builder.get_variable(read_c));
    builder.create_big_add_gate({ read_a, read_b, read_c, ram_sum, 1, 1, 1, -1, 0 }, false);

    // Non-native field arithmetic.
    {
        using fq_ct = stdlib::bigfield<Builder, Bn254FqParams>;
        const fq_ct x = fq_ct::from_witness(&builder, fq::random_element());
        const fq_ct y = fq_ct::from_witness(&builder, fq::random_element());
        [[maybe_unused]] const fq_ct product = x * y;
    }

    // Poseidon2, both round types.
    {
        using field_ct = stdlib::field_t<Builder>;
        using witness_ct = stdlib::witness_t<Builder>;
        std::vector<field_ct> inputs;
        inputs.reserve(4);
        for (size_t i = 0; i < 4; ++i) {
            inputs.emplace_back(witness_ct(&builder, FF::random_element()));
        }
        stdlib::poseidon2<Builder>::hash(inputs);
    }

    return builder;
}

std::vector<FF> random_polynomial(size_t length)
{
    std::vector<FF> poly(length);
    for (FF& coefficient : poly) {
        coefficient = FF::random_element();
    }
    return poly;
}

Commitment perturb(const Commitment& point)
{
    return Commitment(GroupElement(point) + GroupElement(Commitment::one()));
}

/**
 * @brief A prover that opens whatever polynomials it is handed, honestly, without checking that they
 * describe anything.
 *
 * @details Everything the opening argument covers is correct in the result: the commitments are real
 * commitments to the given columns, the claimed evaluations are the real evaluations, and `W`, `W'`
 * are the real batched quotients. The only thing missing is any reason to believe the polynomials
 * satisfy the circuit - which is exactly the gap the verifier's quotient identity closes.
 *
 * The preprocessed groups are not forgeable: the verification key pins their commitments.
 */
Proof forge(const ProvingKey& key,
            std::span<const std::vector<FF>> wires,
            std::span<const std::vector<FF>> memory,
            std::span<const std::vector<FF>> grand_product,
            std::span<const std::vector<FF>> quotient)
{
    const std::vector<FF> packed_wires = fflonk_plonk::pack_columns(wires);
    const std::vector<FF> packed_memory = fflonk_plonk::pack_columns(memory);
    const std::vector<FF> packed_grand_product = fflonk_plonk::pack_columns(grand_product);
    const std::vector<FF> packed_quotient = fflonk_plonk::pack_columns(quotient);

    const CommitmentKey<Curve>& commitment_key = *key.commitment_key;
    Proof proof;
    proof.wires = commitment_key.commit(Polynomial<FF>(std::span<const FF>(packed_wires)));
    proof.memory = commitment_key.commit(Polynomial<FF>(std::span<const FF>(packed_memory)));
    proof.grand_product = commitment_key.commit(Polynomial<FF>(std::span<const FF>(packed_grand_product)));
    proof.quotient = commitment_key.commit(Polynomial<FF>(std::span<const FF>(packed_quotient)));

    Transcript transcript;
    transcript.absorb(key.verification_key.hash());
    for (const FF& public_input : key.public_inputs()) {
        transcript.absorb(public_input);
    }
    transcript.absorb(proof.wires);
    static_cast<void>(transcript.squeeze()); // eta
    static_cast<void>(transcript.squeeze()); // rom_logup_gamma
    transcript.absorb(proof.memory);
    static_cast<void>(transcript.squeeze()); // beta
    static_cast<void>(transcript.squeeze()); // gamma
    transcript.absorb(proof.grand_product);
    static_cast<void>(transcript.squeeze()); // alpha
    transcript.absorb(proof.quotient);

    std::array<const std::vector<FF>*, NUM_GROUPS> packed_groups{};
    for (size_t g = 0; g < NUM_PREPROCESSED_GROUPS; ++g) {
        packed_groups[g] = &key.packed_preprocessed[g];
    }
    packed_groups[GROUP_WIRES] = &packed_wires;
    packed_groups[GROUP_MEMORY] = &packed_memory;
    packed_groups[GROUP_GRAND_PRODUCT] = &packed_grand_product;
    packed_groups[GROUP_QUOTIENT] = &packed_quotient;

    static_cast<void>(prover_detail::open_and_batch(key, packed_groups, transcript, proof));
    return proof;
}

Proof forge_random(const ProvingKey& key)
{
    const size_t n = key.circuit_size();
    std::vector<std::vector<FF>> wires;
    for (size_t i = 0; i < PACK_WIRES; ++i) {
        wires.push_back(random_polynomial(n + NUM_WITNESS_BLINDERS));
    }
    std::vector<std::vector<FF>> memory;
    for (size_t i = 0; i < PACK_MEMORY; ++i) {
        memory.push_back(random_polynomial(n + NUM_WITNESS_BLINDERS));
    }
    std::vector<std::vector<FF>> grand_product;
    for (size_t i = 0; i < PACK_GRAND_PRODUCT; ++i) {
        grand_product.push_back(random_polynomial(n + NUM_WITNESS_BLINDERS));
    }
    std::vector<std::vector<FF>> quotient;
    for (size_t i = 0; i < NUM_QUOTIENT_CHUNKS; ++i) {
        quotient.push_back(random_polynomial(n + 1));
    }
    return forge(key, wires, memory, grand_product, quotient);
}

} // namespace

class UltraFflonkTest : public ::testing::Test {
  protected:
    static void SetUpTestSuite() { srs::init_file_crs_factory(srs::bb_crs_path()); }
};

// ---------------------------------------------------------------------------------------------
// Completeness
// ---------------------------------------------------------------------------------------------

TEST_F(UltraFflonkTest, ProvesAndVerifies)
{
    Builder builder = arithmetic_circuit(/*num_public_inputs=*/4, /*num_gates=*/16);
    ProvingKey key = preprocess(builder);
    const Proof proof = prove(key);
    EXPECT_TRUE(verify(key.verification_key, proof, key.public_inputs()));
}

TEST_F(UltraFflonkTest, ProvesAndVerifiesEveryGateType)
{
    Builder builder = assorted_circuit();
    ProvingKey key = preprocess(builder);

    // Every gate selector must be non-trivial, or "it proves" says nothing about that relation.
    for (const auto& selector : key.instance->polynomials.get_gate_selectors()) {
        bool non_zero = false;
        for (const FF& coefficient : selector.coeffs()) {
            non_zero = non_zero || !coefficient.is_zero();
        }
        EXPECT_TRUE(non_zero);
    }

    const Proof proof = prove(key);
    EXPECT_TRUE(verify(key.verification_key, proof, key.public_inputs()));
}

TEST_F(UltraFflonkTest, ProvesWithNoPublicInputs)
{
    Builder builder = arithmetic_circuit(/*num_public_inputs=*/0, /*num_gates=*/8);
    ProvingKey key = preprocess(builder);
    const Proof proof = prove(key);
    EXPECT_TRUE(key.public_inputs().empty());
    EXPECT_TRUE(verify(key.verification_key, proof, key.public_inputs()));
}

TEST_F(UltraFflonkTest, ProofsAreRandomisedButBothVerify)
{
    Builder builder = arithmetic_circuit(2, 8);
    ProvingKey key = preprocess(builder);

    const Proof first = prove(key);
    const Proof second = prove(key);

    EXPECT_NE(first.to_buffer(), second.to_buffer());
    EXPECT_TRUE(verify(key.verification_key, first, key.public_inputs()));
    EXPECT_TRUE(verify(key.verification_key, second, key.public_inputs()));
}

TEST_F(UltraFflonkTest, ProofSizeIsFixed)
{
    Builder small = arithmetic_circuit(1, 4);
    Builder large = arithmetic_circuit(1, 4000);
    ProvingKey small_key = preprocess(small);
    ProvingKey large_key = preprocess(large);

    EXPECT_NE(small_key.circuit_size(), large_key.circuit_size());
    EXPECT_EQ(prove(small_key).to_buffer().size(), Proof::SIZE_IN_BYTES);
    EXPECT_EQ(prove(large_key).to_buffer().size(), Proof::SIZE_IN_BYTES);
}

TEST_F(UltraFflonkTest, SerialisationRoundTrips)
{
    Builder builder = arithmetic_circuit(2, 6);
    ProvingKey key = preprocess(builder);
    const std::vector<uint8_t> bytes = prove(key).to_buffer();

    Proof parsed;
    ASSERT_TRUE(Proof::from_buffer(bytes, parsed));
    EXPECT_EQ(parsed.to_buffer(), bytes);
    EXPECT_TRUE(verify(key.verification_key, bytes, key.public_inputs()));
}

// ---------------------------------------------------------------------------------------------
// The verifier trusts nothing
// ---------------------------------------------------------------------------------------------

TEST_F(UltraFflonkTest, RejectsMalformedEncodings)
{
    Builder builder = arithmetic_circuit(2, 6);
    ProvingKey key = preprocess(builder);
    const std::vector<uint8_t> bytes = prove(key).to_buffer();
    const std::vector<FF> publics = key.public_inputs();
    ASSERT_TRUE(verify(key.verification_key, bytes, publics));

    {
        const std::vector<uint8_t> truncated(bytes.begin(), bytes.end() - 1);
        EXPECT_FALSE(verify(key.verification_key, truncated, publics));
    }
    {
        std::vector<uint8_t> extended = bytes;
        extended.push_back(0);
        EXPECT_FALSE(verify(key.verification_key, extended, publics));
    }
    {
        // A scalar at or above the field modulus must be rejected rather than silently reduced, or a
        // proof would have many encodings and the transcript binding would be broken.
        std::vector<uint8_t> non_canonical = bytes;
        const size_t scalar_offset = NUM_PROOF_COMMITMENTS * 64;
        for (size_t byte = 0; byte < 32; ++byte) {
            non_canonical[scalar_offset + byte] = 0xff;
        }
        Proof parsed;
        EXPECT_FALSE(Proof::from_buffer(non_canonical, parsed));
        EXPECT_FALSE(verify(key.verification_key, non_canonical, publics));
    }
    {
        std::vector<uint8_t> off_curve = bytes;
        off_curve[31] ^= 1;
        Proof parsed;
        EXPECT_FALSE(Proof::from_buffer(off_curve, parsed));
    }
    {
        // The point at infinity encodes as (0, 0) and is rejected: an honest prover never emits one.
        std::vector<uint8_t> infinity = bytes;
        for (size_t byte = 0; byte < 64; ++byte) {
            infinity[byte] = 0;
        }
        Proof parsed;
        EXPECT_FALSE(Proof::from_buffer(infinity, parsed));
    }
}

TEST_F(UltraFflonkTest, RejectsAnAllZeroProof)
{
    Builder builder = arithmetic_circuit(2, 6);
    ProvingKey key = preprocess(builder);

    const Proof empty;
    EXPECT_FALSE(verify(key.verification_key, empty, key.public_inputs()));

    const std::vector<uint8_t> zeros(Proof::SIZE_IN_BYTES, 0);
    EXPECT_FALSE(verify(key.verification_key, zeros, key.public_inputs()));
}

TEST_F(UltraFflonkTest, RejectsATamperedCommitment)
{
    Builder builder = arithmetic_circuit(2, 6);
    ProvingKey key = preprocess(builder);
    const Proof proof = prove(key);
    const std::vector<FF> publics = key.public_inputs();

    for (Commitment Proof::* field :
         { &Proof::wires, &Proof::memory, &Proof::grand_product, &Proof::quotient, &Proof::w, &Proof::w_prime }) {
        Proof tampered = proof;
        tampered.*field = perturb(tampered.*field);
        EXPECT_FALSE(verify(key.verification_key, tampered, publics));
    }
}

TEST_F(UltraFflonkTest, RejectsATamperedVerificationKey)
{
    Builder builder = arithmetic_circuit(2, 6);
    ProvingKey key = preprocess(builder);
    const Proof proof = prove(key);
    const std::vector<FF> publics = key.public_inputs();

    for (size_t g = 0; g < NUM_PREPROCESSED_GROUPS; ++g) {
        VerificationKey tampered = key.verification_key;
        tampered.preprocessed[g] = perturb(tampered.preprocessed[g]);
        EXPECT_FALSE(verify(tampered, proof, publics));
    }
    {
        VerificationKey tampered = key.verification_key;
        tampered.num_public_inputs += 1;
        EXPECT_FALSE(verify(tampered, proof, publics));
    }
    {
        VerificationKey tampered = key.verification_key;
        tampered.pub_inputs_offset += 1;
        EXPECT_FALSE(verify(tampered, proof, publics));
    }
}

TEST_F(UltraFflonkTest, RejectsAProofRecycledAcrossPublicInputs)
{
    // Same circuit shape, different public input values: the proof for one must not verify for the
    // other, which is what binds the statement rather than merely the circuit.
    Builder first = arithmetic_circuit(2, 6);
    ProvingKey first_key = preprocess(first);
    const Proof proof = prove(first_key);

    std::vector<FF> wrong = first_key.public_inputs();
    wrong[0] += FF::one();

    EXPECT_TRUE(verify(first_key.verification_key, proof, first_key.public_inputs()));
    EXPECT_FALSE(verify(first_key.verification_key, proof, wrong));
}

// ---------------------------------------------------------------------------------------------
// The two halves of soundness, separated: the opening argument, and the quotient identity
// ---------------------------------------------------------------------------------------------

TEST_F(UltraFflonkTest, AnHonestProofPassesEveryCheck)
{
    Builder builder = assorted_circuit();
    ProvingKey key = preprocess(builder);
    const VerificationReport report = verify_detailed(key.verification_key, prove(key), key.public_inputs());

    EXPECT_TRUE(report.well_formed);
    EXPECT_TRUE(report.quotient);
    EXPECT_TRUE(report.opening);
}

TEST_F(UltraFflonkTest, ArbitraryPolynomialsOpenCorrectlyAndAreRejectedByTheQuotientIdentity)
{
    // The test that says what the quotient identity is for. These polynomials are random, so they
    // describe no circuit at all - but they are committed and opened honestly, so the batched opening
    // and the pairing accept them. Only the relations stand in the way.
    Builder builder = arithmetic_circuit(3, 9);
    ProvingKey key = preprocess(builder);

    const Proof forged = forge_random(key);
    const VerificationReport report = verify_detailed(key.verification_key, forged, key.public_inputs());

    EXPECT_TRUE(report.well_formed);
    EXPECT_TRUE(report.opening) << "the forgery's openings should be valid - that is the point of the test";
    EXPECT_FALSE(report.quotient);
    EXPECT_FALSE(verify(key.verification_key, forged, key.public_inputs()));
}

TEST_F(UltraFflonkTest, ForgedOpeningsCannotSubstituteThePreprocessedGroups)
{
    // The preprocessed groups are the ones a forger cannot choose, because the verification key pins
    // their commitments. Confirm the opening really does bind them: opening a different preprocessed
    // polynomial makes the opening itself fail, not merely the quotient identity.
    Builder builder = arithmetic_circuit(3, 9);
    ProvingKey key = preprocess(builder);

    for (size_t g = 0; g < NUM_PREPROCESSED_GROUPS; ++g) {
        Builder other = arithmetic_circuit(3, 9);
        ProvingKey spliced = preprocess(other);
        spliced.packed_preprocessed[g][0] += FF::one();

        const Proof forged = forge_random(spliced);
        EXPECT_FALSE(verify_detailed(key.verification_key, forged, key.public_inputs()).opening)
            << "substituting preprocessed group " << g;
    }
}

TEST_F(UltraFflonkTest, EveryEvaluationTheRelationsReadIsConstrained)
{
    // Perturb one claimed evaluation at a time and read off whether the quotient identity notices.
    // Sumcheck's subrelations are batched into a single identity here, so the three-wire suite's
    // per-identity table becomes a per-evaluation one - and it is just as sharp, because it pins the
    // exact set of evaluations no relation reads.
    //
    // Those are the three columns that ride along in a two-point group without being shifted by any
    // relation: the lookup read counts, the read tags, and the lookup inverses, each at xi*omega.
    // Their second blinder exists precisely because they are revealed there for nothing.
    Builder builder = assorted_circuit();
    ProvingKey key = preprocess(builder);
    const Proof proof = prove(key);
    const std::vector<FF> publics = key.public_inputs();
    ASSERT_TRUE(verify(key.verification_key, proof, publics));

    // Offsets of the evaluations the relations never read.
    size_t offset = 0;
    std::vector<size_t> unread;
    for (size_t g = 0; g < NUM_GROUPS; ++g) {
        offset += GROUP_SHAPES[g].pack;
        if (!GROUP_SHAPES[g].two_point) {
            continue;
        }
        if (g == GROUP_MEMORY) {
            unread.push_back(offset + MEMORY_READ_COUNTS);
            unread.push_back(offset + MEMORY_READ_TAGS);
        }
        if (g == GROUP_GRAND_PRODUCT) {
            unread.push_back(offset + GRAND_PRODUCT_LOOKUP_INVERSES);
        }
        offset += GROUP_SHAPES[g].pack;
    }
    ASSERT_EQ(unread.size(), 3);

    for (size_t index = 0; index < NUM_EVALUATIONS; ++index) {
        Proof tampered = proof;
        tampered.evaluations[index] += FF::one();
        const VerificationReport report = verify_detailed(key.verification_key, tampered, publics);

        const bool is_unread = std::find(unread.begin(), unread.end(), index) != unread.end();
        EXPECT_EQ(report.quotient, is_unread) << "quotient identity, perturbing evaluation " << index;
        // Evaluations are absorbed before nu, so a perturbation also breaks the opening. The two
        // mechanisms are independent and both have to hold.
        EXPECT_FALSE(report.opening) << "opening, perturbing evaluation " << index;
        EXPECT_FALSE(verify(key.verification_key, tampered, publics));
    }
}

TEST_F(UltraFflonkTest, TheQuotientIdentityReadsThePublicInputsThroughTheirDelta)
{
    // A proof that is internally perfect - every opening honest, the pairing satisfied, the same
    // transcript on both sides - but where the verifier derives a different `public_input_delta` from
    // the prover. Nothing but the quotient identity, which is the only place the public inputs reach
    // the relations, can reject it.
    //
    // Unlike the three-wire system, the public inputs do not enter through an additive `PI(xi)` term
    // that the prover's own polynomials are free to ignore: they enter the permutation argument,
    // which the prover's quotient has to satisfy. That is why the lie has to be planted on the
    // verifier's side of the delta rather than in the prover's claim - see the prover-side test
    // below for the other direction.
    Builder builder = arithmetic_circuit(2, 6);
    ProvingKey key = preprocess(builder);
    const Proof proof = prove(key);
    const std::vector<FF> publics = key.public_inputs();
    ASSERT_TRUE(verify(key.verification_key, proof, publics));

    // `pub_inputs_offset` shifts the delta without touching anything else the verifier replays -
    // except the key hash, so put the mutation on both sides of the transcript.
    VerificationKey shifted = key.verification_key;
    shifted.pub_inputs_offset += 1;

    ProvingKey lying = preprocess(builder);
    lying.verification_key.pub_inputs_offset += 1;
    const Proof lie = prove(lying);

    const VerificationReport report = verify_detailed(shifted, lie, publics);
    EXPECT_TRUE(report.well_formed);
    EXPECT_TRUE(report.opening) << "the lie's openings are honest - that is what makes this a test of the identity";
    EXPECT_FALSE(report.quotient);
    EXPECT_FALSE(verify(shifted, lie, publics));
}

// ---------------------------------------------------------------------------------------------
// The prover refuses to prove a false statement
// ---------------------------------------------------------------------------------------------

TEST_F(UltraFflonkTest, ProverRejectsAWitnessThatViolatesAGate)
{
    Builder builder = arithmetic_circuit(1, 4);
    ProvingKey key = preprocess(builder);
    key.instance->polynomials.w_l().at(Flavor::TRACE_OFFSET + 1) += FF::one();
    EXPECT_ANY_THROW(prove(key));
}

TEST_F(UltraFflonkTest, ProverRejectsAClaimedPublicInputItsWitnessDoesNotContain)
{
    // The public inputs reach the relations through `public_input_delta`, so claiming one the wires
    // do not carry makes the permutation relation stop vanishing on the trace - and the prover finds
    // out when its quotient does not divide, rather than emitting a proof that fails on chain.
    Builder builder = arithmetic_circuit(2, 6);
    ProvingKey key = preprocess(builder);
    key.instance->public_inputs[0] += FF::one();
    EXPECT_ANY_THROW(prove(key));
}
