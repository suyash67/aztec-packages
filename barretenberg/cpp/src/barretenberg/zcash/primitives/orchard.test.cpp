#include "orchard.hpp"
#include "barretenberg/numeric/random/engine.hpp"
#include "barretenberg/zcash/test_vectors/orchard_vectors.hpp"
#include "orchard_test_utils.hpp"

#include <gtest/gtest.h>

using namespace bb;
using namespace bb::zcash;

using PastaOrchard = Orchard<PastaCycle>;
using Fp = PastaCycle::FF;

TEST(ZcashPoseidon, MatchesHalo2PoseidonPallasConstants)
{
    const auto& p = PoseidonP128Pow5T3<Fp>::params();
    // halo2_poseidon fp.rs ROUND_CONSTANTS[0] and MDS[0][0]
    EXPECT_EQ(uint256_t(p.round_constants[0][0]),
              uint256_t("0x360d7470611e473d353f628f76d110f34e71162f31003b7057538c2596426303"));
    EXPECT_EQ(uint256_t(p.round_constants[0][1]),
              uint256_t("0x2bab94d7ae222d135dc3c6c5febfaa314908ac2f12ebe06fbdb74213bf63188b"));
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            Fp acc = 0;
            for (size_t k = 0; k < 3; ++k) {
                acc += p.mds[i][k] * p.mds_inv[k][j];
            }
            EXPECT_EQ(acc, Fp(i == j ? 1 : 0));
        }
    }
}

TEST(ZcashPoseidon, PermutationMatchesHalo2Reference)
{
    // halo2_poseidon p128pow5t3.rs test_against_reference (Fp)
    PoseidonP128Pow5T3<Fp>::State s{ Fp(0), Fp(1), Fp(2) };
    PoseidonP128Pow5T3<Fp>::permute(s);
    EXPECT_EQ(uint256_t(s[0]), uint256_t("0x2a526acd0b64b45394efb364f966240ff7e69a71d0b642a0aeb1bc024aeca456"));
    EXPECT_EQ(uint256_t(s[1]), uint256_t("0x13c5d1568b4aa43076ff7dae343d5512dcd42e7fbed9dafe012a3e9628e5b82a"));
    EXPECT_EQ(uint256_t(s[2]), uint256_t("0x0a49c868c6976544256fcd597984561af7cfdfe1bda42c7b359029a1d34e9ddd"));
}

TEST(ZcashOrchard, GeneratorsMatchOrchard)
{
    const auto& c = PastaOrchard::constants();
    // orchard constants/fixed_bases/spend_auth_g.rs GENERATOR (little-endian bytes -> big-endian hex)
    EXPECT_EQ(uint256_t(c.spend_auth_g.x),
              uint256_t("0x375523b328f1d6063b8d187c3e5f445f0c7f0ce37b70a10c8d1a7284b875c963"));
    // Recomputing the published z-values for a few windows reproduces them (find_zs_and_us).
    auto zs = FixedBase<PastaCycle>::find_zs(c.value_commit_v, PastaOrchard::NUM_WINDOWS_SHORT);
    for (size_t w = 0; w < 3; ++w) {
        EXPECT_EQ(zs[w], orchard_constants::VALUE_COMMIT_V_Z[w]) << "window " << w;
    }
}

TEST(ZcashOrchard, FixedBaseLagrangeCoefficientsInterpolateWindowTable)
{
    const auto& fb = PastaOrchard::fixed_bases().fb_spend_auth_g;
    for (size_t w : { size_t{ 0 }, size_t{ 1 }, size_t{ 84 } }) {
        for (size_t k = 0; k < 8; ++k) {
            Fp x = 0;
            Fp kp = 1;
            for (size_t d = 0; d < 8; ++d) {
                x += fb.lagrange_coeffs[w][d] * kp;
                kp *= Fp(k);
            }
            EXPECT_EQ(x, fb.window_table[w][k].x);
            EXPECT_EQ(fb.u[w][k].sqr(), fb.window_table[w][k].y + Fp(fb.z[w]));
        }
    }
}

TEST(ZcashOrchard, NativeActionMatchesProductionVectors)
{
    // The witness and public inputs were exported from orchard 0.16.0's Action circuit for real bundles; evaluating the
    // Action statement natively must reproduce every public input exactly.
    for (const auto* v : { &test_vectors::ACTION_REAL_SPEND,
                           &test_vectors::ACTION_DUMMY_SPEND_0,
                           &test_vectors::ACTION_DUMMY_SPEND_1 }) {
        const auto w = parse_action_witness<PastaCycle>(*v);
        auto pi = PastaOrchard::evaluate(w);
        ASSERT_TRUE(pi.has_value());
        auto fields = pi->to_field_elements();
        // Dummy spends (v_old = 0) are not bound to the anchor: the circuit leaves the root unchecked.
        const size_t first = (w.v_old == 0) ? 1 : 0;
        for (size_t i = first; i < 10; ++i) {
            EXPECT_EQ(uint256_t(fields[i]), uint256_t(v->public_inputs[i])) << "public input " << i;
        }
    }
}

TEST(ZcashOrchard, RandomWitnessSatisfiesStatement)
{
    auto& engine = numeric::get_debug_randomness();
    auto w = PastaOrchard::random_witness(engine);
    EXPECT_TRUE(PastaOrchard::evaluate(w).has_value());
    w.pk_d_old = PastaOrchard::AffineElement(PastaOrchard::Element(w.pk_d_old) + PastaOrchard::Element(w.g_d_old));
    EXPECT_FALSE(PastaOrchard::evaluate(w).has_value());
}
