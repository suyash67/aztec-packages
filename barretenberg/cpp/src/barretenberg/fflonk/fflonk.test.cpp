#include "barretenberg/fflonk/prover.hpp"
#include "barretenberg/fflonk/transcript.hpp"
#include "barretenberg/fflonk/verifier.hpp"
#include "barretenberg/srs/global_crs.hpp"

#include <gtest/gtest.h>

using namespace bb;
using namespace bb::fflonk_plonk;

namespace {

/**
 * @brief A circuit exercising every construct the builder offers: public inputs, multiplication,
 * addition, constants, a boolean constraint, and a copy constraint added after the fact.
 */
CircuitBuilder build_circuit(size_t num_public_inputs, size_t num_rounds)
{
    CircuitBuilder builder;

    std::vector<uint32_t> publics;
    publics.reserve(num_public_inputs);
    for (size_t i = 0; i < num_public_inputs; ++i) {
        publics.push_back(builder.add_public_input(FF(static_cast<uint64_t>(i) + 3)));
    }

    uint32_t accumulator = publics.empty() ? builder.add_variable(FF(7)) : publics[0];
    for (size_t i = 0; i < num_rounds; ++i) {
        const uint32_t operand = builder.add_variable(FF(static_cast<uint64_t>(i) + 1));
        accumulator = builder.create_mul(accumulator, operand);
        accumulator = builder.create_add(accumulator, operand);
    }

    // A second variable holding the same value, merged after the fact: the copy constraint has to
    // carry it, since no gate relates the two.
    const uint32_t duplicate = builder.add_variable(builder.get_variable(accumulator));
    builder.assert_equal(duplicate, accumulator);
    builder.create_add_gate(
        FF::one(), duplicate, -FF::one(), accumulator, FF::zero(), builder.zero_variable(), FF::zero());

    builder.create_bool_gate(builder.one_variable());
    builder.fix_variable(accumulator, builder.get_variable(accumulator));

    return builder;
}

std::vector<FF> public_inputs_of(const ProvingKey& key)
{
    return key.trace.public_inputs;
}

Commitment perturb(const Commitment& point)
{
    return Commitment(GroupElement(point) + GroupElement(Commitment::one()));
}

std::vector<FF> random_polynomial(size_t length)
{
    std::vector<FF> poly(length);
    for (FF& coefficient : poly) {
        coefficient = FF::random_element();
    }
    return poly;
}

/**
 * @brief A prover that opens whatever polynomials it is handed, honestly, without checking that
 * they describe anything.
 *
 * @details Everything the opening argument covers is correct in the result: the commitments are
 * real commitments to `columns`, the claimed evaluations are the real evaluations, and `W`, `W'`
 * are the real batched quotients. The only thing missing is any reason to believe the polynomials
 * satisfy the circuit - which is exactly the gap the verifier's three constraint identities close.
 *
 * Group 0 is not forgeable: it is the preprocessed group, and the verification key pins its
 * commitment.
 */
Proof forge(const ProvingKey& key,
            std::span<const std::vector<FF>> wire_columns,
            std::span<const std::vector<FF>> grand_product_columns,
            std::span<const std::vector<FF>> quotient_columns)
{
    const std::vector<FF> packed_wires = pack_columns(wire_columns);
    const std::vector<FF> packed_grand_product = pack_columns(grand_product_columns);
    const std::vector<FF> packed_quotients = pack_columns(quotient_columns);

    Proof proof;
    proof.c1 = key.commitment_key->commit(Polynomial<FF>(std::span<const FF>(packed_wires)));
    proof.c2 = key.commitment_key->commit(Polynomial<FF>(std::span<const FF>(packed_grand_product)));
    proof.c3 = key.commitment_key->commit(Polynomial<FF>(std::span<const FF>(packed_quotients)));

    Transcript transcript;
    transcript.absorb(key.verification_key.hash());
    for (const FF& public_input : key.trace.public_inputs) {
        transcript.absorb(public_input);
    }
    transcript.absorb(proof.c1);
    static_cast<void>(transcript.squeeze()); // beta
    static_cast<void>(transcript.squeeze()); // gamma
    transcript.absorb(proof.c2);
    transcript.absorb(proof.c3);

    const std::array<const std::vector<FF>*, NUM_GROUPS> packed_groups = {
        &key.packed_preprocessed, &packed_wires, &packed_grand_product, &packed_quotients
    };
    static_cast<void>(prover_detail::open_and_batch(key, packed_groups, transcript, proof));
    return proof;
}

} // namespace

class FflonkTest : public ::testing::Test {
  protected:
    static void SetUpTestSuite() { srs::init_file_crs_factory(srs::bb_crs_path()); }
};

// ---------------------------------------------------------------------------------------------
// Completeness
// ---------------------------------------------------------------------------------------------

TEST_F(FflonkTest, ProvesAndVerifies)
{
    const CircuitBuilder builder = build_circuit(/*num_public_inputs=*/4, /*num_rounds=*/10);
    const ProvingKey key = preprocess(builder);

    std::string failure;
    ASSERT_TRUE(key.trace.check(failure)) << failure;

    const Proof proof = prove(key);
    EXPECT_TRUE(verify(key.verification_key, proof, public_inputs_of(key)));
}

TEST_F(FflonkTest, ProvesAndVerifiesAcrossSizes)
{
    for (const size_t rounds : { size_t(1), size_t(3), size_t(10), size_t(40), size_t(200) }) {
        const CircuitBuilder builder = build_circuit(/*num_public_inputs=*/2, rounds);
        const ProvingKey key = preprocess(builder);
        const Proof proof = prove(key);
        EXPECT_TRUE(verify(key.verification_key, proof, public_inputs_of(key)))
            << "failed at " << rounds << " rounds, circuit size " << key.verification_key.circuit_size;
    }
}

TEST_F(FflonkTest, ProvesWithNoPublicInputs)
{
    const CircuitBuilder builder = build_circuit(/*num_public_inputs=*/0, /*num_rounds=*/6);
    const ProvingKey key = preprocess(builder);
    ASSERT_EQ(key.verification_key.num_public_inputs, 0u);

    const Proof proof = prove(key);
    EXPECT_TRUE(verify(key.verification_key, proof, {}));
}

TEST_F(FflonkTest, ProvesWithManyPublicInputs)
{
    const CircuitBuilder builder = build_circuit(/*num_public_inputs=*/37, /*num_rounds=*/5);
    const ProvingKey key = preprocess(builder);
    const Proof proof = prove(key);
    EXPECT_TRUE(verify(key.verification_key, proof, public_inputs_of(key)));
}

TEST_F(FflonkTest, ProofsAreRandomisedButBothVerify)
{
    const CircuitBuilder builder = build_circuit(/*num_public_inputs=*/2, /*num_rounds=*/8);
    const ProvingKey key = preprocess(builder);

    const Proof first = prove(key);
    const Proof second = prove(key);

    // Blinding is live: the same statement proven twice must not produce the same bytes.
    EXPECT_NE(first.to_buffer(), second.to_buffer());
    EXPECT_TRUE(verify(key.verification_key, first, public_inputs_of(key)));
    EXPECT_TRUE(verify(key.verification_key, second, public_inputs_of(key)));
}

TEST_F(FflonkTest, ProofSizeIsFixed)
{
    for (const size_t rounds : { size_t(2), size_t(50) }) {
        const CircuitBuilder builder = build_circuit(/*num_public_inputs=*/1, rounds);
        const ProvingKey key = preprocess(builder);
        EXPECT_EQ(prove(key).to_buffer().size(), Proof::SIZE_IN_BYTES);
    }
    EXPECT_EQ(Proof::SIZE_IN_BYTES, 928u);
}

// ---------------------------------------------------------------------------------------------
// Internal consistency: the packing identity, and the trace checker
// ---------------------------------------------------------------------------------------------

TEST_F(FflonkTest, ResiduesAreTheColumnEvaluations)
{
    // g mod (X^t - z) = sum_i f_i(z) X^i is the identity the whole scheme rests on. Check it
    // against direct evaluation, on the real preprocessed group.
    const CircuitBuilder builder = build_circuit(/*num_public_inputs=*/3, /*num_rounds=*/7);
    const ProvingKey key = preprocess(builder);

    const FF point = FF::random_element();
    std::vector<FF> quotient;
    std::vector<FF> residue;
    divide_by_power_minus(key.packed_preprocessed, PACK_PREPROCESSED, point, quotient, residue);

    ASSERT_EQ(residue.size(), PACK_PREPROCESSED);
    for (size_t i = 0; i < PACK_PREPROCESSED; ++i) {
        EXPECT_EQ(residue[i], evaluate(key.preprocessed.columns[i], point)) << "column " << i;
    }

    // ... and the quotient really is the quotient.
    std::vector<FF> reconstructed(key.packed_preprocessed.size(), FF::zero());
    for (size_t i = 0; i < quotient.size(); ++i) {
        if (i + PACK_PREPROCESSED < reconstructed.size()) {
            reconstructed[i + PACK_PREPROCESSED] += quotient[i];
        }
        reconstructed[i] -= point * quotient[i];
    }
    for (size_t i = 0; i < residue.size(); ++i) {
        reconstructed[i] += residue[i];
    }
    EXPECT_EQ(reconstructed, key.packed_preprocessed);
}

TEST_F(FflonkTest, TraceCheckerCatchesABrokenGate)
{
    CircuitBuilder builder = build_circuit(/*num_public_inputs=*/1, /*num_rounds=*/3);
    Trace trace = builder.build_trace();

    std::string failure;
    ASSERT_TRUE(trace.check(failure)) << failure;

    trace.q_c[trace.num_public_inputs + 1] += FF::one();
    EXPECT_FALSE(trace.check(failure));
    EXPECT_NE(failure.find("gate identity"), std::string::npos);
}

TEST_F(FflonkTest, TraceCheckerCatchesABrokenCopyConstraint)
{
    CircuitBuilder builder = build_circuit(/*num_public_inputs=*/1, /*num_rounds=*/3);
    Trace trace = builder.build_trace();

    // Two slots hold the zero variable; break one of them without touching the other.
    trace.wires[1][trace.circuit_size - 1] = FF(5);
    trace.wires[2][trace.circuit_size - 1] = FF(5);
    std::string failure;
    EXPECT_FALSE(trace.check(failure));
    EXPECT_NE(failure.find("copy constraint"), std::string::npos);
}

TEST_F(FflonkTest, PreprocessingIsDeterministic)
{
    const ProvingKey first = preprocess(build_circuit(2, 9));
    const ProvingKey second = preprocess(build_circuit(2, 9));
    EXPECT_EQ(first.verification_key, second.verification_key);
    EXPECT_EQ(first.packed_preprocessed, second.packed_preprocessed);
}

// ---------------------------------------------------------------------------------------------
// Serialisation
// ---------------------------------------------------------------------------------------------

TEST_F(FflonkTest, SerialisationRoundTrips)
{
    const ProvingKey key = preprocess(build_circuit(3, 12));
    const Proof proof = prove(key);
    const std::vector<uint8_t> bytes = proof.to_buffer();

    Proof parsed;
    ASSERT_TRUE(Proof::from_buffer(bytes, parsed));
    EXPECT_EQ(parsed.to_buffer(), bytes);
    EXPECT_TRUE(verify(key.verification_key, bytes, public_inputs_of(key)));
}

TEST_F(FflonkTest, RejectsMalformedEncodings)
{
    const ProvingKey key = preprocess(build_circuit(2, 6));
    const std::vector<uint8_t> bytes = prove(key).to_buffer();
    const std::vector<FF> publics = public_inputs_of(key);

    {
        std::vector<uint8_t> truncated(bytes.begin(), bytes.end() - 1);
        EXPECT_FALSE(verify(key.verification_key, truncated, publics));
    }
    {
        std::vector<uint8_t> extended = bytes;
        extended.push_back(0);
        EXPECT_FALSE(verify(key.verification_key, extended, publics));
    }
    {
        // A scalar at or above the field modulus must be rejected rather than silently reduced,
        // or a proof would have many encodings and the transcript binding would be broken.
        std::vector<uint8_t> non_canonical = bytes;
        const size_t scalar_offset = Proof::NUM_COMMITMENTS * 64;
        for (size_t byte = 0; byte < 32; ++byte) {
            non_canonical[scalar_offset + byte] = 0xff;
        }
        Proof parsed;
        EXPECT_FALSE(Proof::from_buffer(non_canonical, parsed));
        EXPECT_FALSE(verify(key.verification_key, non_canonical, publics));
    }
    {
        // A point off the curve.
        std::vector<uint8_t> off_curve = bytes;
        off_curve[63] ^= 0x01;
        Proof parsed;
        EXPECT_FALSE(Proof::from_buffer(off_curve, parsed));
        EXPECT_FALSE(verify(key.verification_key, off_curve, publics));
    }
    {
        // A coordinate at or above the base-field modulus.
        std::vector<uint8_t> non_canonical_point = bytes;
        for (size_t byte = 0; byte < 32; ++byte) {
            non_canonical_point[byte] = 0xff;
        }
        Proof parsed;
        EXPECT_FALSE(Proof::from_buffer(non_canonical_point, parsed));
    }
    {
        // The point at infinity, encoded as (0, 0). Rejected so that this verifier and the Solidity
        // one - whose y^2 = x^3 + 3 test has no natural encoding for it - agree on every input.
        std::vector<uint8_t> infinity = bytes;
        for (size_t byte = 0; byte < 64; ++byte) {
            infinity[byte] = 0;
        }
        Proof parsed;
        EXPECT_FALSE(Proof::from_buffer(infinity, parsed));
        EXPECT_FALSE(verify(key.verification_key, infinity, publics));
    }
}

// ---------------------------------------------------------------------------------------------
// Soundness: every part of the proof and its context must be binding
// ---------------------------------------------------------------------------------------------

TEST_F(FflonkTest, RejectsATamperedEvaluation)
{
    const ProvingKey key = preprocess(build_circuit(3, 9));
    const Proof proof = prove(key);
    const std::vector<FF> publics = public_inputs_of(key);
    ASSERT_TRUE(verify(key.verification_key, proof, publics));

    for (size_t i = 0; i < NUM_EVALUATIONS; ++i) {
        Proof tampered = proof;
        tampered.evaluations[i] += FF::one();
        EXPECT_FALSE(verify(key.verification_key, tampered, publics)) << "evaluation " << i << " is not binding";
    }
}

TEST_F(FflonkTest, RejectsATamperedCommitment)
{
    const ProvingKey key = preprocess(build_circuit(3, 9));
    const Proof proof = prove(key);
    const std::vector<FF> publics = public_inputs_of(key);

    for (size_t i = 0; i < Proof::NUM_COMMITMENTS; ++i) {
        Proof tampered = proof;
        Commitment* points[] = { &tampered.c1, &tampered.c2, &tampered.c3, &tampered.w, &tampered.w_prime };
        *points[i] = perturb(*points[i]);
        EXPECT_FALSE(verify(key.verification_key, tampered, publics)) << "commitment " << i << " is not binding";
    }
}

TEST_F(FflonkTest, RejectsSwappedCommitments)
{
    const ProvingKey key = preprocess(build_circuit(2, 9));
    const Proof proof = prove(key);
    const std::vector<FF> publics = public_inputs_of(key);

    Proof swapped = proof;
    std::swap(swapped.c2, swapped.c3);
    EXPECT_FALSE(verify(key.verification_key, swapped, publics));

    Proof swapped_openings = proof;
    std::swap(swapped_openings.w, swapped_openings.w_prime);
    EXPECT_FALSE(verify(key.verification_key, swapped_openings, publics));
}

TEST_F(FflonkTest, RejectsWrongPublicInputs)
{
    const ProvingKey key = preprocess(build_circuit(4, 9));
    const Proof proof = prove(key);
    std::vector<FF> publics = public_inputs_of(key);
    ASSERT_TRUE(verify(key.verification_key, proof, publics));

    for (size_t i = 0; i < publics.size(); ++i) {
        std::vector<FF> tampered = publics;
        tampered[i] += FF::one();
        EXPECT_FALSE(verify(key.verification_key, proof, tampered)) << "public input " << i << " is not bound";
    }

    std::vector<FF> too_few(publics.begin(), publics.end() - 1);
    EXPECT_FALSE(verify(key.verification_key, proof, too_few));

    std::vector<FF> too_many = publics;
    too_many.push_back(FF::one());
    EXPECT_FALSE(verify(key.verification_key, proof, too_many));
}

TEST_F(FflonkTest, RejectsAProofForAnotherCircuit)
{
    const ProvingKey first = preprocess(build_circuit(2, 9));
    const ProvingKey second = preprocess(build_circuit(2, 11));
    ASSERT_NE(first.verification_key.c0, second.verification_key.c0);

    const Proof proof = prove(first);
    EXPECT_TRUE(verify(first.verification_key, proof, public_inputs_of(first)));
    EXPECT_FALSE(verify(second.verification_key, proof, public_inputs_of(first)));
}

TEST_F(FflonkTest, RejectsATamperedVerificationKey)
{
    const ProvingKey key = preprocess(build_circuit(2, 9));
    const Proof proof = prove(key);
    const std::vector<FF> publics = public_inputs_of(key);

    {
        VerificationKey tampered = key.verification_key;
        tampered.c0 = perturb(tampered.c0);
        EXPECT_FALSE(verify(tampered, proof, publics));
    }
    {
        VerificationKey tampered = key.verification_key;
        tampered.k1 += FF::one();
        EXPECT_FALSE(verify(tampered, proof, publics));
    }
    {
        VerificationKey tampered = key.verification_key;
        tampered.num_public_inputs += 1;
        EXPECT_FALSE(verify(tampered, proof, publics));
    }
}

TEST_F(FflonkTest, RejectsAnAllZeroProof)
{
    const ProvingKey key = preprocess(build_circuit(2, 9));
    const std::vector<FF> publics = public_inputs_of(key);

    Proof empty;
    EXPECT_FALSE(verify(key.verification_key, empty, publics));

    const std::vector<uint8_t> zeros(Proof::SIZE_IN_BYTES, 0);
    EXPECT_FALSE(verify(key.verification_key, zeros, publics));
}

TEST_F(FflonkTest, RejectsAProofRecycledAcrossPublicInputs)
{
    // Same circuit shape, different public input values: the proof for one must not verify for the
    // other, which is what binds the statement rather than merely the circuit.
    CircuitBuilder first;
    const uint32_t a = first.add_public_input(FF(11));
    first.fix_variable(a, FF(11));
    const ProvingKey first_key = preprocess(first, 64);

    CircuitBuilder second;
    const uint32_t b = second.add_public_input(FF(12));
    second.fix_variable(b, FF(12));
    const ProvingKey second_key = preprocess(second, 64);

    const std::vector<FF> eleven = { FF(11) };
    const std::vector<FF> twelve = { FF(12) };
    const Proof proof = prove(first_key);
    EXPECT_TRUE(verify(first_key.verification_key, proof, eleven));
    EXPECT_FALSE(verify(first_key.verification_key, proof, twelve));
    EXPECT_FALSE(verify(second_key.verification_key, proof, twelve));
}

// ---------------------------------------------------------------------------------------------
// The prover refuses to prove a false statement
// ---------------------------------------------------------------------------------------------

TEST_F(FflonkTest, ProverRejectsAnUnsatisfiedGate)
{
    CircuitBuilder builder = build_circuit(1, 4);
    Trace trace = builder.build_trace();
    trace.q_c[trace.num_public_inputs] += FF::one();

    const ProvingKey key = preprocess(std::move(trace));
    EXPECT_ANY_THROW(prove(key));
}

TEST_F(FflonkTest, ProverRejectsAViolatedCopyConstraint)
{
    CircuitBuilder builder;
    const uint32_t x = builder.add_variable(FF(3));
    const uint32_t y = builder.add_variable(FF(3));
    builder.assert_equal(x, y);
    // Two gates that read the merged variable, so the copy constraint has real work to do.
    builder.create_add_gate(FF::one(), x, -FF::one(), y, FF::zero(), builder.zero_variable(), FF::zero());

    Trace trace = builder.build_trace();
    // Break one slot of the cycle without touching the rest: the grand product cannot close.
    trace.wires[1][trace.num_public_inputs] = FF(4);

    const ProvingKey key = preprocess(std::move(trace));
    EXPECT_ANY_THROW(prove(key));
}

// ---------------------------------------------------------------------------------------------
// The two halves of soundness, separated: the opening argument, and the constraint identities
// ---------------------------------------------------------------------------------------------

TEST_F(FflonkTest, AnHonestProofPassesEveryCheck)
{
    const ProvingKey key = preprocess(build_circuit(3, 9));
    const VerificationReport report = verify_detailed(key.verification_key, prove(key), public_inputs_of(key));

    EXPECT_TRUE(report.well_formed);
    EXPECT_TRUE(report.gate_identity);
    EXPECT_TRUE(report.grand_product_start);
    EXPECT_TRUE(report.permutation_identity);
    EXPECT_TRUE(report.opening);
}

TEST_F(FflonkTest, ArbitraryPolynomialsOpenCorrectlyAndAreRejectedByTheConstraints)
{
    // The test that says what the constraint identities are for. These polynomials are random, so
    // they describe no circuit at all - but they are committed and opened honestly, so the batched
    // opening and the pairing accept them. Only the three identities stand in the way.
    const ProvingKey key = preprocess(build_circuit(3, 9));
    const size_t n = key.verification_key.circuit_size;

    std::vector<std::vector<FF>> wire_columns;
    for (size_t i = 0; i < PACK_WIRES; ++i) {
        wire_columns.push_back(random_polynomial(n + 2));
    }
    std::vector<std::vector<FF>> grand_product_columns = { random_polynomial(n + 3) };
    std::vector<std::vector<FF>> quotient_columns;
    for (size_t i = 0; i < PACK_QUOTIENTS; ++i) {
        quotient_columns.push_back(random_polynomial(n + 6));
    }

    const Proof forged = forge(key, wire_columns, grand_product_columns, quotient_columns);
    const VerificationReport report = verify_detailed(key.verification_key, forged, public_inputs_of(key));

    EXPECT_TRUE(report.well_formed);
    EXPECT_TRUE(report.opening) << "the forgery's openings should be valid - that is the point of the test";
    EXPECT_FALSE(report.gate_identity);
    EXPECT_FALSE(report.grand_product_start);
    EXPECT_FALSE(report.permutation_identity);
    EXPECT_FALSE(verify(key.verification_key, forged, public_inputs_of(key)));
}

TEST_F(FflonkTest, ForgedOpeningsCannotSubstituteThePreprocessedGroup)
{
    // The preprocessed group is the one a forger cannot choose, because the verification key pins
    // its commitment. Confirm the opening really does bind it: opening a different preprocessed
    // polynomial makes the opening itself fail, not merely the identities.
    const ProvingKey key = preprocess(build_circuit(3, 9));

    ProvingKey spliced = preprocess(build_circuit(3, 9));
    spliced.packed_preprocessed[0] += FF::one();

    const size_t n = key.verification_key.circuit_size;
    std::vector<std::vector<FF>> wire_columns;
    for (size_t i = 0; i < PACK_WIRES; ++i) {
        wire_columns.push_back(random_polynomial(n + 2));
    }
    std::vector<std::vector<FF>> grand_product_columns = { random_polynomial(n + 3) };
    std::vector<std::vector<FF>> quotient_columns;
    for (size_t i = 0; i < PACK_QUOTIENTS; ++i) {
        quotient_columns.push_back(random_polynomial(n + 6));
    }

    const Proof forged = forge(spliced, wire_columns, grand_product_columns, quotient_columns);
    EXPECT_FALSE(verify_detailed(key.verification_key, forged, public_inputs_of(key)).opening);
}

TEST_F(FflonkTest, EachIdentityConstrainsExactlyTheEvaluationsItShould)
{
    // Perturb one evaluation at a time and read off which identities notice. If an identity were
    // silently vacuous - a check that always passes - it would show up here as a false negative.
    struct Expectation {
        size_t index;
        bool gate;
        bool start;
        bool permutation;
        const char* name;
    };
    constexpr std::array<Expectation, NUM_EVALUATIONS> expectations = { {
        { EVAL_Q_L, false, true, true, "q_l" },
        { EVAL_Q_R, false, true, true, "q_r" },
        { EVAL_Q_O, false, true, true, "q_o" },
        { EVAL_Q_M, false, true, true, "q_m" },
        { EVAL_Q_C, false, true, true, "q_c" },
        { EVAL_S_1, true, true, false, "s_1" },
        { EVAL_S_2, true, true, false, "s_2" },
        { EVAL_S_3, true, true, false, "s_3" },
        { EVAL_A, false, true, false, "a" },
        { EVAL_B, false, true, false, "b" },
        { EVAL_C, false, true, false, "c" },
        { EVAL_T0_LO, false, true, true, "t0_lo" },
        { EVAL_T0_HI, false, true, true, "t0_hi" },
        { EVAL_Z, true, false, false, "z" },
        { EVAL_Z_OMEGA, true, true, false, "z_omega" },
        { EVAL_T1, true, false, true, "t1" },
        { EVAL_T2_LO, true, true, false, "t2_lo" },
        { EVAL_T2_MID, true, true, false, "t2_mid" },
        { EVAL_T2_HI, true, true, false, "t2_hi" },
    } };

    const ProvingKey key = preprocess(build_circuit(3, 9));
    const Proof proof = prove(key);
    const std::vector<FF> publics = public_inputs_of(key);
    ASSERT_TRUE(verify(key.verification_key, proof, publics));

    for (const Expectation& expectation : expectations) {
        Proof tampered = proof;
        tampered.evaluations[expectation.index] += FF::one();
        const VerificationReport report = verify_detailed(key.verification_key, tampered, publics);

        EXPECT_EQ(report.gate_identity, expectation.gate) << "gate identity, perturbing " << expectation.name;
        EXPECT_EQ(report.grand_product_start, expectation.start) << "start identity, perturbing " << expectation.name;
        EXPECT_EQ(report.permutation_identity, expectation.permutation)
            << "permutation identity, perturbing " << expectation.name;
        // Evaluations are absorbed before nu, so a perturbation also breaks the opening. The two
        // mechanisms are independent and both have to hold.
        EXPECT_FALSE(report.opening) << "opening, perturbing " << expectation.name;
    }
}

TEST_F(FflonkTest, LyingAboutAPublicInputIsCaughtByTheGateIdentity)
{
    // A prover whose proof is internally perfect - every opening honest, the pairing satisfied -
    // but who claims a public input its witness does not contain. Nothing but the gate identity,
    // which is the only place the verifier's own PI(xi) enters, can reject this.
    CircuitBuilder builder;
    const uint32_t input = builder.add_public_input(FF(11));
    builder.fix_variable(input, FF(11));

    ProvingKey key = preprocess(builder, 64);
    // The witness, the selectors and the PI polynomial all still say 11, so the prover's own
    // quotient divisions and self-check are satisfied; only the claim it broadcasts changes.
    key.trace.public_inputs[0] = FF(12);

    const Proof proof = prove(key);
    const std::vector<FF> claimed = { FF(12) };
    const VerificationReport report = verify_detailed(key.verification_key, proof, claimed);

    EXPECT_TRUE(report.well_formed);
    EXPECT_TRUE(report.opening);
    EXPECT_TRUE(report.grand_product_start);
    EXPECT_TRUE(report.permutation_identity);
    EXPECT_FALSE(report.gate_identity);
    EXPECT_FALSE(verify(key.verification_key, proof, claimed));
}

// ---------------------------------------------------------------------------------------------
// The transcript
// ---------------------------------------------------------------------------------------------

TEST_F(FflonkTest, TranscriptChainsState)
{
    Transcript transcript;
    const FF first = transcript.squeeze();
    const FF second = transcript.squeeze();
    EXPECT_NE(first, second);

    Transcript replay;
    EXPECT_EQ(replay.squeeze(), first);
    EXPECT_EQ(replay.squeeze(), second);

    Transcript diverged;
    diverged.absorb(FF::one());
    EXPECT_NE(diverged.squeeze(), first);
}

TEST_F(FflonkTest, TranscriptAbsorbsInfinityAsZero)
{
    Transcript with_infinity;
    with_infinity.absorb(Commitment::infinity());

    Transcript with_zeros;
    with_zeros.absorb_word(uint256_t(0));
    with_zeros.absorb_word(uint256_t(0));

    EXPECT_EQ(with_infinity.squeeze(), with_zeros.squeeze());
}
