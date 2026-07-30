#include "barretenberg/commitment_schemes/whir/whir.hpp"
#include "barretenberg/commitment_schemes/commitment_key.hpp"
#include "barretenberg/commitment_schemes/kzg/kzg.hpp"
#include "barretenberg/commitment_schemes/shplonk/shplemini.hpp"
#include "barretenberg/commitment_schemes/utils/mock_witness_generator.hpp"
#include "barretenberg/commitment_schemes/verification_key.hpp"
#include "barretenberg/srs/global_crs.hpp"
#include "barretenberg/transcript/transcript.hpp"

#include <benchmark/benchmark.h>

// Benchmarks WHIR against the Gemini+Shplonk+KZG (Shplemini) opening stack on the UltraHonk claim
// shape: 36 committed polynomials of which 5 are to-be-shifted (41 evaluation claims at a common
// point). Commit benchmarks measure one polynomial; open/verify benchmarks measure the full batched
// opening. Proof sizes are reported as counters (32 bytes per transcript field element).

namespace bb::whir {
namespace {

constexpr size_t NUM_POLYNOMIALS = 36;
constexpr size_t NUM_TO_BE_SHIFTED = 5;
constexpr size_t SECURITY_BITS = 100;

std::vector<fr> random_array(size_t size)
{
    std::vector<fr> array(size);
    for (fr& value : array) {
        value = fr::random_element();
    }
    return array;
}

std::vector<fr> random_point(size_t size)
{
    return random_array(size);
}

WhirConfig bench_config(size_t log_n, size_t log_inv_rate)
{
    return WhirConfig::create(log_n, SECURITY_BITS, log_inv_rate, /*folding_factor_bits=*/4, /*final_poly_bits=*/4);
}

/** @brief The Ultra-shaped batch: 31 unshifted-only + 5 to-be-shifted polynomials, 41 claims. */
template <typename Hasher> struct WhirBenchInstance {
    WhirCommitmentKey<Hasher> ck;
    std::vector<std::vector<fr>> arrays;
    std::vector<WhirProverData<Hasher>> data;
    std::vector<fr> u;
    typename WhirProver<Hasher>::Claims claims;
    typename WhirVerifier<Hasher>::Claims verifier_claims;

    WhirBenchInstance(size_t log_n,
                      size_t log_inv_rate,
                      size_t num_polynomials = NUM_POLYNOMIALS,
                      size_t num_to_be_shifted = NUM_TO_BE_SHIFTED)
        : ck(bench_config(log_n, log_inv_rate))
        , u(random_point(log_n))
    {
        const size_t n = size_t(1) << log_n;
        for (size_t p = 0; p < num_polynomials; ++p) {
            const bool to_be_shifted = p >= num_polynomials - num_to_be_shifted;
            std::vector<fr> array = random_array(n);
            if (to_be_shifted) {
                array[0] = fr::zero();
            }
            arrays.push_back(std::move(array));
            data.push_back(ck.commit(std::span<const fr>(arrays.back()), to_be_shifted));
        }
        for (size_t p = 0; p < num_polynomials; ++p) {
            claims.unshifted.push_back(&data[p]);
            claims.unshifted_evaluations.push_back(Polynomial<fr>(std::span<const fr>(arrays[p])).evaluate_mle(u));
        }
        for (size_t p = num_polynomials - num_to_be_shifted; p < num_polynomials; ++p) {
            std::vector<fr> shifted(arrays[p].begin() + 1, arrays[p].end());
            shifted.push_back(fr::zero());
            claims.to_be_shifted.push_back(&data[p]);
            claims.shifted_evaluations.push_back(Polynomial<fr>(std::span<const fr>(shifted)).evaluate_mle(u));
        }
        verifier_claims = { claims.unshifted_evaluations, claims.shifted_evaluations, {}, {} };
    }

    HonkProof prove() const
    {
        auto transcript = NativeTranscript::test_prover_init_empty();
        WhirProver<Hasher>::prove(ck, claims, u, transcript);
        return transcript->export_proof();
    }

    bool verify(const HonkProof& proof) const
    {
        auto transcript = std::make_shared<NativeTranscript>(proof);
        [[maybe_unused]] auto init = transcript->template receive_from_prover<fr>("Init");
        return WhirVerifier<Hasher>::verify(ck.config, verifier_claims, u, transcript);
    }
};

template <typename Hasher> void whir_commit(benchmark::State& state)
{
    const size_t log_n = static_cast<size_t>(state.range(0));
    const size_t log_inv_rate = static_cast<size_t>(state.range(1));
    WhirCommitmentKey<Hasher> ck(bench_config(log_n, log_inv_rate));
    const std::vector<fr> array = random_array(size_t(1) << log_n);
    for (auto _ : state) {
        auto data = ck.commit(std::span<const fr>(array));
        benchmark::DoNotOptimize(data.tree.root());
    }
}

template <typename Hasher> void whir_open(benchmark::State& state)
{
    const size_t log_n = static_cast<size_t>(state.range(0));
    const size_t log_inv_rate = static_cast<size_t>(state.range(1));
    WhirBenchInstance<Hasher> instance(log_n, log_inv_rate);
    HonkProof proof;
    for (auto _ : state) {
        proof = instance.prove();
    }
    state.counters["proof_KiB"] = static_cast<double>(proof.size() * 32) / 1024.0;
}

template <typename Hasher> void whir_verify(benchmark::State& state)
{
    const size_t log_n = static_cast<size_t>(state.range(0));
    const size_t log_inv_rate = static_cast<size_t>(state.range(1));
    WhirBenchInstance<Hasher> instance(log_n, log_inv_rate);
    const HonkProof proof = instance.prove();
    bool ok = true;
    for (auto _ : state) {
        ok = ok && instance.verify(proof);
    }
    if (!ok) {
        state.SkipWithError("WHIR verification failed");
    }
}

// Single-polynomial configuration matching the WHIR reference implementation's PCS benchmarks
// (WizardOfMenlo/whir defaults: one polynomial, one evaluation claim, k = 4, λ = 100 conjectured).

// Timed region: commit + claimed-evaluation computation + opening proof, as in the paper's
// "time to commit and open".
template <typename Hasher> void whir_single_commit_open(benchmark::State& state)
{
    const size_t log_n = static_cast<size_t>(state.range(0));
    const size_t log_inv_rate = static_cast<size_t>(state.range(1));
    WhirCommitmentKey<Hasher> ck(bench_config(log_n, log_inv_rate));
    const std::vector<fr> array = random_array(size_t(1) << log_n);
    const std::vector<fr> u = random_point(log_n);
    HonkProof proof;
    for (auto _ : state) {
        WhirProverData<Hasher> data = ck.commit(std::span<const fr>(array));
        typename WhirProver<Hasher>::Claims claims;
        claims.unshifted = { &data };
        claims.unshifted_evaluations = { Polynomial<fr>(std::span<const fr>(array)).evaluate_mle(u) };
        auto transcript = NativeTranscript::test_prover_init_empty();
        WhirProver<Hasher>::prove(ck, claims, u, transcript);
        proof = transcript->export_proof();
    }
    state.counters["proof_KiB"] = static_cast<double>(proof.size() * 32) / 1024.0;
}

template <typename Hasher> void whir_single_verify(benchmark::State& state)
{
    const size_t log_n = static_cast<size_t>(state.range(0));
    const size_t log_inv_rate = static_cast<size_t>(state.range(1));
    WhirCommitmentKey<Hasher> ck(bench_config(log_n, log_inv_rate));
    const std::vector<fr> array = random_array(size_t(1) << log_n);
    const std::vector<fr> u = random_point(log_n);
    WhirProverData<Hasher> data = ck.commit(std::span<const fr>(array));
    typename WhirProver<Hasher>::Claims claims;
    claims.unshifted = { &data };
    claims.unshifted_evaluations = { Polynomial<fr>(std::span<const fr>(array)).evaluate_mle(u) };
    auto prover_transcript = NativeTranscript::test_prover_init_empty();
    WhirProver<Hasher>::prove(ck, claims, u, prover_transcript);
    const HonkProof proof = prover_transcript->export_proof();

    const typename WhirVerifier<Hasher>::Claims verifier_claims{ claims.unshifted_evaluations, {}, {}, {} };
    bool ok = true;
    for (auto _ : state) {
        auto transcript = std::make_shared<NativeTranscript>(proof);
        [[maybe_unused]] auto init = transcript->template receive_from_prover<fr>("Init");
        ok = ok && WhirVerifier<Hasher>::verify(ck.config, verifier_claims, u, transcript);
    }
    if (!ok) {
        state.SkipWithError("WHIR verification failed");
    }
}

BENCHMARK_TEMPLATE(whir_single_commit_open, Blake3sMerkleHasher)
    ->ArgsProduct({ { 20, 22 }, { 1, 4 } })
    ->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(whir_single_commit_open, Poseidon2MerkleHasher)
    ->ArgsProduct({ { 22 }, { 1, 4 } })
    ->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(whir_single_verify, Blake3sMerkleHasher)
    ->ArgsProduct({ { 20, 22 }, { 1, 4 } })
    ->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(whir_single_verify, Poseidon2MerkleHasher)
    ->ArgsProduct({ { 22 }, { 1, 4 } })
    ->Unit(benchmark::kMillisecond);

// --- Gemini+Shplonk+KZG baseline on the identical claim shape ---

using Curve = curve::BN254;
using KzgCK = CommitmentKey<Curve>;
using ShpleminiProver = ShpleminiProver_<Curve>;
using ShpleminiVerifier = ShpleminiVerifier_<Curve, /*HasZK=*/false, /*has_gemini_masking=*/false>;

KzgCK create_kzg_ck(size_t n)
{
    srs::init_file_crs_factory(srs::bb_crs_path());
    return KzgCK(n);
}

void kzg_commit(benchmark::State& state)
{
    const size_t log_n = static_cast<size_t>(state.range(0));
    const size_t n = size_t(1) << log_n;
    KzgCK ck = create_kzg_ck(n);
    Polynomial<fr> poly = Polynomial<fr>::random(n);
    for (auto _ : state) {
        benchmark::DoNotOptimize(ck.commit(poly));
    }
}

void shplemini_open(benchmark::State& state)
{
    const size_t log_n = static_cast<size_t>(state.range(0));
    const size_t n = size_t(1) << log_n;
    KzgCK ck = create_kzg_ck(n);
    std::vector<fr> u = random_point(log_n);
    MockClaimGenerator<Curve> mock_claims(n, NUM_POLYNOMIALS, NUM_TO_BE_SHIFTED, u, ck);
    HonkProof proof;
    for (auto _ : state) {
        auto transcript = NativeTranscript::test_prover_init_empty();
        const auto opening_claim = ShpleminiProver::prove(n, mock_claims.polynomial_batcher, u, ck, transcript);
        KZG<Curve>::compute_opening_proof(ck, opening_claim, transcript);
        proof = transcript->export_proof();
    }
    state.counters["proof_KiB"] = static_cast<double>(proof.size() * 32) / 1024.0;
}

void shplemini_verify(benchmark::State& state)
{
    const size_t log_n = static_cast<size_t>(state.range(0));
    const size_t n = size_t(1) << log_n;
    KzgCK ck = create_kzg_ck(n);
    std::vector<fr> u = random_point(log_n);
    MockClaimGenerator<Curve> mock_claims(n, NUM_POLYNOMIALS, NUM_TO_BE_SHIFTED, u, ck);

    auto prover_transcript = NativeTranscript::test_prover_init_empty();
    const auto opening_claim = ShpleminiProver::prove(n, mock_claims.polynomial_batcher, u, ck, prover_transcript);
    KZG<Curve>::compute_opening_proof(ck, opening_claim, prover_transcript);
    const HonkProof proof = prover_transcript->export_proof();

    VerifierCommitmentKey<Curve> vk;
    bool ok = true;
    for (auto _ : state) {
        auto transcript = std::make_shared<NativeTranscript>(proof);
        [[maybe_unused]] auto init = transcript->template receive_from_prover<fr>("Init");
        auto batch_opening_claim = ShpleminiVerifier::compute_batch_opening_claim(
                                       mock_claims.claim_batcher, u, vk.get_g1_identity(), transcript)
                                       .batch_opening_claim;
        const auto pairing_points =
            KZG<Curve>::reduce_verify_batch_opening_claim(std::move(batch_opening_claim), transcript);
        ok = ok && pairing_points.check();
    }
    if (!ok) {
        state.SkipWithError("Shplemini verification failed");
    }
}

BENCHMARK_TEMPLATE(whir_commit, Blake3sMerkleHasher)
    ->ArgsProduct({ { 12, 16, 20 }, { 1, 2, 3 } })
    ->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(whir_commit, Poseidon2MerkleHasher)
    ->ArgsProduct({ { 12, 16, 20 }, { 1, 2, 3 } })
    ->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(whir_open, Blake3sMerkleHasher)
    ->ArgsProduct({ { 12, 16, 20 }, { 1, 2 } })
    ->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(whir_open, Poseidon2MerkleHasher)
    ->ArgsProduct({ { 12, 16, 20 }, { 2 } })
    ->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(whir_verify, Blake3sMerkleHasher)
    ->ArgsProduct({ { 12, 16, 20 }, { 1, 2 } })
    ->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(whir_verify, Poseidon2MerkleHasher)
    ->ArgsProduct({ { 12, 16, 20 }, { 2 } })
    ->Unit(benchmark::kMillisecond);
BENCHMARK(kzg_commit)->Args({ 12 })->Args({ 16 })->Args({ 20 })->Unit(benchmark::kMillisecond);
BENCHMARK(shplemini_open)->Args({ 12 })->Args({ 16 })->Args({ 20 })->Unit(benchmark::kMillisecond);
BENCHMARK(shplemini_verify)->Args({ 12 })->Args({ 16 })->Args({ 20 })->Unit(benchmark::kMillisecond);

} // namespace
} // namespace bb::whir

BENCHMARK_MAIN();
