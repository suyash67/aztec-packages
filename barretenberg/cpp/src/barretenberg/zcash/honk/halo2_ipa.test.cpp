#include "halo2_ipa.hpp"
#include "barretenberg/zcash/honk/pasta_crs.hpp"
#include "barretenberg/zcash/honk/pasta_transcript.hpp"

#include <gtest/gtest.h>

using namespace bb;
using namespace bb::zcash;

namespace {
using Curve = curve::Vesta;
using Fr = Curve::ScalarField;
using IPA = Halo2IPA<Curve>;
} // namespace

TEST(ZcashHalo2IPA, OpensRandomPolynomial)
{
    for (size_t log_n : { 2UL, 5UL, 10UL }) {
        const size_t n = size_t{ 1 } << log_n;
        const auto g = halo2_vesta_generators(n);
        const IPA::Generators gens{ g, halo2_vesta_w(), halo2_vesta_u() };
        auto poly = Polynomial<Fr>::random(n);
        const Fr x = Fr::random_element();
        const Fr v = poly.evaluate(x);
        const Fr blind = Fr::random_element();
        const auto commitment =
            Curve::AffineElement(IPA::msm(std::span<const Fr>(poly.data(), n), g) + Curve::Element(gens.w) * blind);

        auto prover_transcript = std::make_shared<PastaTranscript>();
        IPA::prove(gens, ProverOpeningClaim<Curve>{ poly.share(), { x, v } }, blind, prover_transcript);
        const auto proof = prover_transcript->export_proof();
        EXPECT_EQ(proof.size(), 2 + (4 * log_n) + 2) << "S, L/R per round, c, f";

        auto verifier_transcript = std::make_shared<PastaTranscript>(proof);
        EXPECT_TRUE(IPA::verify(gens, OpeningClaim<Curve>{ { x, v }, commitment }, verifier_transcript));

        auto bad_transcript = std::make_shared<PastaTranscript>(proof);
        EXPECT_FALSE(IPA::verify(gens, OpeningClaim<Curve>{ { x, v + Fr(1) }, commitment }, bad_transcript));
    }
}
