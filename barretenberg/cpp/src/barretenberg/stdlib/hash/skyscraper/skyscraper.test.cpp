#include "skyscraper.hpp"

#include "barretenberg/circuit_checker/circuit_checker.hpp"
#include "barretenberg/crypto/skyscraper/skyscraper.hpp"
#include "barretenberg/stdlib/hash/poseidon2/poseidon2.hpp"
#include "barretenberg/stdlib_circuit_builders/ultra_circuit_builder.hpp"

#include <gtest/gtest.h>

using namespace bb;

namespace {
using Builder = UltraCircuitBuilder;
using field_ct = stdlib::field_t<Builder>;
using witness_ct = stdlib::witness_t<Builder>;
using Skyscraper = stdlib::skyscraper::Skyscraper<Builder>;
} // namespace

// The in-circuit hash is only useful if it agrees with the native one bit for bit: a Merkle root
// computed by the prover has to be the root the circuit reconstructs.
TEST(StdlibSkyscraper, MatchesTheNativeImplementation)
{
    Builder builder;
    for (size_t i = 0; i < 8; ++i) {
        const fr left = fr::random_element();
        const fr right = fr::random_element();
        const field_ct l = witness_ct(&builder, left);
        const field_ct r = witness_ct(&builder, right);

        EXPECT_EQ(Skyscraper::compress(l, r).get_value(), crypto::skyscraper::compress(left, right))
            << "compression " << i;

        const auto [native_l, native_r] = crypto::skyscraper::permute(left, right);
        const auto [circuit_l, circuit_r] = Skyscraper::permute(l, r);
        EXPECT_EQ(circuit_l.get_value(), native_l) << "permutation left " << i;
        EXPECT_EQ(circuit_r.get_value(), native_r) << "permutation right " << i;
    }
    EXPECT_TRUE(CircuitChecker::check(builder));
}

// Edge values: the bar reads the canonical byte string, so a value whose representation is all
// zeros or all ones exercises the ends of the decomposition.
TEST(StdlibSkyscraper, MatchesTheNativeImplementationOnEdgeValues)
{
    Builder builder;
    const std::vector<fr> values{ fr(0), fr(1), fr(-1), fr(uint256_t(1) << 128) };
    for (const fr& left : values) {
        for (const fr& right : values) {
            const field_ct l = witness_ct(&builder, left);
            const field_ct r = witness_ct(&builder, right);
            EXPECT_EQ(Skyscraper::compress(l, r).get_value(), crypto::skyscraper::compress(left, right));
        }
    }
    EXPECT_TRUE(CircuitChecker::check(builder));
}

TEST(StdlibSkyscraper, FoldMatchesTheNativeImplementation)
{
    Builder builder;
    for (const size_t width : { size_t(1), size_t(2), size_t(5), size_t(24) }) {
        std::vector<fr> native(width);
        std::vector<field_ct> circuit;
        for (fr& value : native) {
            value = fr::random_element();
            circuit.push_back(witness_ct(&builder, value));
        }
        EXPECT_EQ(Skyscraper::fold_compress(circuit).get_value(), crypto::skyscraper::fold_compress(native))
            << "fold of width " << width;
    }
    EXPECT_TRUE(CircuitChecker::check(builder));
}

// Reports the marginal cost of a compression, which is what prices Skyscraper as a WHIR Merkle hash
// in a recursive verifier. The first compression also pays for the shared range lists and the
// 256-row S-box table, so the marginal figure is the one that scales.
// The other half of the Fiat-Shamir trade. Skyscraper's sponge absorbs one field element per
// permutation; bb's Poseidon2 transcript absorbs three, and Ultra gives Poseidon2 custom gates. So
// the hash that is cheaper natively is the dearer one in a recursive verifier, and this is the
// number that says by how much.
TEST(StdlibSkyscraper, TranscriptHashGateCount)
{
    for (const size_t width : { size_t(3), size_t(30), size_t(300) }) {
        Builder sky_builder;
        Builder p2_builder;
        std::vector<field_ct> sky_input;
        std::vector<stdlib::field_t<Builder>> p2_input;
        for (size_t i = 0; i < width; ++i) {
            sky_input.push_back(witness_ct(&sky_builder, fr::random_element()));
            p2_input.push_back(witness_ct(&p2_builder, fr::random_element()));
        }
        const size_t sky_before = sky_builder.get_num_finalized_gates_inefficient();
        const size_t p2_before = p2_builder.get_num_finalized_gates_inefficient();
        Skyscraper::hash(sky_input);
        stdlib::poseidon2<Builder>::hash(p2_input);
        const double sky = double(sky_builder.get_num_finalized_gates_inefficient() - sky_before) / double(width);
        const double p2 = double(p2_builder.get_num_finalized_gates_inefficient() - p2_before) / double(width);
        info("transcript hash of ", width, " elements: skyscraper ", sky, " gates/elt, poseidon2 ", p2, " gates/elt");
        EXPECT_TRUE(CircuitChecker::check(sky_builder));
    }
}

TEST(StdlibSkyscraper, GateCount)
{
    size_t previous = 0;
    for (const size_t count : { size_t(1), size_t(2), size_t(10) }) {
        Builder builder;
        field_ct l = witness_ct(&builder, fr::random_element());
        field_ct r = witness_ct(&builder, fr::random_element());
        const size_t before = builder.get_num_finalized_gates_inefficient();
        for (size_t i = 0; i < count; ++i) {
            l = Skyscraper::compress(l, r);
        }
        const size_t total = builder.get_num_finalized_gates_inefficient() - before;
        info("skyscraper x",
             count,
             ": ",
             total,
             " gates",
             count > 1 ? " | marginal per compression = " + std::to_string((total - previous) / (count - 1)) : "");
        if (count == 1) {
            previous = total;
        }
        EXPECT_TRUE(CircuitChecker::check(builder));
    }
}
