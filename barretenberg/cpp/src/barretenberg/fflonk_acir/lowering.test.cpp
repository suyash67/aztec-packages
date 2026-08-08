#include "barretenberg/fflonk_acir/lowering.hpp"

#include "barretenberg/fflonk/prover.hpp"
#include "barretenberg/fflonk/verifier.hpp"
#include "barretenberg/srs/global_crs.hpp"
#include "barretenberg/stdlib/hash/poseidon2/poseidon2.hpp"

#include <gtest/gtest.h>

using namespace bb;
using namespace bb::fflonk_acir;

namespace {

using FF = fflonk_plonk::FF;

/** @brief An Ultra circuit of arithmetic gates only, with public inputs and copy constraints. */
UltraCircuitBuilder arithmetic_circuit()
{
    UltraCircuitBuilder builder;

    const FF a = FF(7);
    const FF b = FF(11);
    const uint32_t a_idx = builder.add_public_variable(a);
    const uint32_t b_idx = builder.add_public_variable(b);

    // A multiplication, so q_m is non-trivial.
    const uint32_t product = builder.add_variable(a * b);
    builder.create_big_mul_add_gate({ a_idx, b_idx, product, builder.zero_idx(), 1, 0, 0, -1, 0, 0 }, false);

    // A four-wire add, which is where the lowering has to split a row.
    const FF sum = a + b + (a * b);
    const uint32_t sum_idx = builder.add_variable(sum);
    builder.create_big_add_gate({ a_idx, b_idx, product, sum_idx, 1, 1, 1, -1, 0 }, false);

    // A second variable holding the same value, merged after the fact: only the copy constraint
    // relates the two, so the lowering has to carry it across.
    const uint32_t duplicate = builder.add_variable(sum);
    builder.assert_equal(duplicate, sum_idx);
    builder.create_big_add_gate({ duplicate, sum_idx, builder.zero_idx(), builder.zero_idx(), 1, -1, 0, 0, 0 }, false);

    return builder;
}

} // namespace

class FflonkLoweringTest : public ::testing::Test {
  protected:
    static void SetUpTestSuite() { srs::init_file_crs_factory(srs::bb_crs_path()); }
};

TEST_F(FflonkLoweringTest, LoweredCircuitSatisfiesTheThreeWireTrace)
{
    UltraCircuitBuilder circuit = arithmetic_circuit();
    const fflonk_plonk::CircuitBuilder lowered = lower(circuit);

    // The trace checker is the safety net the lowering is entitled to: a mistranslated selector or a
    // dropped copy constraint shows up here, as a named violation, rather than as a proof that fails
    // somewhere with no explanation.
    std::string failure;
    EXPECT_TRUE(lowered.build_trace().check(failure)) << failure;
}

TEST_F(FflonkLoweringTest, LoweredCircuitProvesAndVerifies)
{
    UltraCircuitBuilder circuit = arithmetic_circuit();
    const fflonk_plonk::CircuitBuilder lowered = lower(circuit);

    const fflonk_plonk::ProvingKey key = fflonk_plonk::preprocess(lowered);
    const fflonk_plonk::Proof proof = fflonk_plonk::prove(key);
    EXPECT_TRUE(fflonk_plonk::verify(key.verification_key, proof, key.trace.public_inputs));
}

TEST_F(FflonkLoweringTest, PublicInputsSurviveInOrder)
{
    UltraCircuitBuilder circuit = arithmetic_circuit();
    const fflonk_plonk::CircuitBuilder lowered = lower(circuit);
    const fflonk_plonk::Trace trace = lowered.build_trace();

    ASSERT_EQ(trace.public_inputs.size(), circuit.public_inputs().size());
    for (size_t i = 0; i < trace.public_inputs.size(); ++i) {
        EXPECT_EQ(trace.public_inputs[i], circuit.get_variable(circuit.public_inputs()[i]));
    }
}

TEST_F(FflonkLoweringTest, AWrongWitnessSurvivesLoweringAsAFailingTrace)
{
    // The lowering is a rewriting of the circuit, not a re-derivation of the witness: an Ultra
    // circuit whose witness does not satisfy it must lower to a three-wire trace that does not
    // satisfy it either, or the lowering would be laundering a false statement into a true one.
    UltraCircuitBuilder circuit;
    const uint32_t a_idx = circuit.add_public_variable(FF(7));
    const uint32_t b_idx = circuit.add_public_variable(FF(11));
    const uint32_t wrong_product = circuit.add_variable(FF(7) * FF(11) + FF::one());
    circuit.create_big_mul_add_gate({ a_idx, b_idx, wrong_product, circuit.zero_idx(), 1, 0, 0, -1, 0, 0 }, false);

    const fflonk_plonk::CircuitBuilder lowered = lower(circuit);
    std::string failure;
    EXPECT_FALSE(lowered.build_trace().check(failure));
}

TEST_F(FflonkLoweringTest, RejectsRangeConstraintsBecauseTheTagArgumentCannotBeCarried)
{
    // The delta-range quartic expands into three wires without trouble. What does not is the other
    // half of Ultra's range argument: that the sorted list the deltas run over is a permutation of the
    // constrained variables, which Ultra enforces with a tag multiset layered onto sigma rather than
    // with copy constraints. Lowering only the quartic would leave a prover free to choose the sorted
    // list, so the gate is refused instead.
    UltraCircuitBuilder circuit = arithmetic_circuit();
    circuit.create_dyadic_range_constraint(circuit.add_variable(FF(200)), 8, "range");
    EXPECT_ANY_THROW(static_cast<void>(lower(circuit)));
}

TEST_F(FflonkLoweringTest, LowersPoseidon2)
{
    // Noir's native hash. Each round becomes a chain of three-wire rows: an S-box per lane, then the
    // matrix.
    UltraCircuitBuilder circuit;
    const uint32_t a_idx = circuit.add_public_variable(FF(7));
    {
        using field_ct = stdlib::field_t<UltraCircuitBuilder>;
        using witness_ct = stdlib::witness_t<UltraCircuitBuilder>;
        std::vector<field_ct> inputs;
        inputs.emplace_back(field_ct::from_witness_index(&circuit, a_idx));
        for (size_t i = 1; i < 4; ++i) {
            inputs.emplace_back(witness_ct(&circuit, FF::random_element()));
        }
        stdlib::poseidon2<UltraCircuitBuilder>::hash(inputs);
    }

    const fflonk_plonk::CircuitBuilder lowered = lower(circuit);
    std::string failure;
    ASSERT_TRUE(lowered.build_trace().check(failure)) << failure;

    const fflonk_plonk::ProvingKey key = fflonk_plonk::preprocess(lowered);
    EXPECT_TRUE(fflonk_plonk::verify(key.verification_key, fflonk_plonk::prove(key), key.trace.public_inputs));
}

TEST_F(FflonkLoweringTest, RejectsGatesWithNoThreeWireEncoding)
{
    {
        UltraCircuitBuilder circuit = arithmetic_circuit();
        const grumpkin::g1::affine_element p1 = grumpkin::g1::affine_element::random_element();
        const grumpkin::g1::affine_element p2 = grumpkin::g1::affine_element::random_element();
        const grumpkin::g1::affine_element p3(grumpkin::g1::element(p1) + grumpkin::g1::element(p2));
        circuit.create_ecc_add_gate({ circuit.add_variable(p1.x),
                                      circuit.add_variable(p1.y),
                                      circuit.add_variable(p2.x),
                                      circuit.add_variable(p2.y),
                                      circuit.add_variable(p3.x),
                                      circuit.add_variable(p3.y),
                                      /*is_addition=*/true });
        EXPECT_ANY_THROW(static_cast<void>(lower(circuit)));
    }
    {
        // A dynamically indexed array needs the RAM/ROM argument, which this arithmetization has no
        // mechanism for at all - not merely no cheap encoding of.
        UltraCircuitBuilder circuit = arithmetic_circuit();
        const size_t rom_id = circuit.create_ROM_array(4);
        for (size_t i = 0; i < 4; ++i) {
            circuit.set_ROM_element(rom_id, i, circuit.add_variable(FF(static_cast<uint64_t>(i) + 1)));
        }
        static_cast<void>(circuit.read_ROM_array(rom_id, circuit.add_variable(FF(2))));
        EXPECT_ANY_THROW(static_cast<void>(lower(circuit)));
    }
}
