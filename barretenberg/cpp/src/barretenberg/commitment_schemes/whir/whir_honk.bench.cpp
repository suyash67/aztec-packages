#include "barretenberg/commitment_schemes/whir/whir_honk.hpp"
#include "barretenberg/commitment_schemes/ligero/ligero_honk.hpp"
#include "barretenberg/commitment_schemes/mercury/mercury_honk.hpp"
#include "barretenberg/special_public_inputs/special_public_inputs.hpp"
#include "barretenberg/srs/global_crs.hpp"
#include "barretenberg/stdlib_circuit_builders/mock_circuits.hpp"
#include "barretenberg/ultra_honk/ultra_prover.hpp"
#include "barretenberg/ultra_honk/ultra_verifier.hpp"

#include <benchmark/benchmark.h>

// End-to-end UltraHonk benchmarks on identical circuits: WhirHonk (hash-based WHIR PCS) versus the
// standard UltraFlavor prover/verifier (Gemini+Shplonk+KZG). Timed regions are full proof
// construction and full verification; circuit building and proving/verification key setup are
// excluded on both sides. Proof sizes are reported as counters (32 bytes per field element).

namespace bb::whir {
namespace {

constexpr size_t WHIR_SECURITY_BITS = 100;
constexpr size_t WHIR_LOG_INV_RATE = 2;

/**
 * @brief Arithmetic + public-input + lookup + ROM circuit padded so the dyadic trace size is
 * 2^target_log_n (asserted by the callers via the sizing instance).
 */
UltraCircuitBuilder build_circuit(size_t target_log_n, bool add_default_io)
{
    UltraCircuitBuilder builder;
    MockCircuits::add_arithmetic_gates_with_public_inputs(builder, 16);
    MockCircuits::add_lookup_gates(builder, 2);
    const size_t rom_id = builder.create_ROM_array(4);
    for (size_t i = 0; i < 4; ++i) {
        builder.set_ROM_element(rom_id, i, builder.add_variable(fr(3 * i + 1)));
    }
    builder.read_ROM_array(rom_id, builder.add_variable(fr(2)));
    MockCircuits::add_arithmetic_gates(builder, (size_t(3) << (target_log_n - 2)) - (size_t(1) << (target_log_n - 3)));
    if (add_default_io) {
        DefaultIO::add_default(builder);
    }
    return builder;
}

size_t actual_log_n(size_t target_log_n, bool add_default_io)
{
    UltraCircuitBuilder builder = build_circuit(target_log_n, add_default_io);
    return ProverInstance_<UltraFlavor>(builder).log_dyadic_size();
}

template <typename Hasher> void whir_honk_prove(benchmark::State& state)
{
    using Honk = WhirHonk<Hasher>;
    const size_t target_log_n = static_cast<size_t>(state.range(0));
    const size_t log_n = actual_log_n(target_log_n, false);
    const WhirConfig config = Honk::make_config(log_n, WHIR_SECURITY_BITS, WHIR_LOG_INV_RATE);
    HonkProof proof;
    for (auto _ : state) {
        state.PauseTiming();
        UltraCircuitBuilder builder = build_circuit(target_log_n, false);
        auto pk = Honk::create_proving_key(builder, config);
        state.ResumeTiming();
        proof = Honk::prove(pk);
    }
    state.counters["proof_KiB"] = static_cast<double>(proof.size() * 32) / 1024.0;
    state.counters["log_n"] = static_cast<double>(log_n);
}

template <typename Hasher> void whir_honk_verify(benchmark::State& state)
{
    using Honk = WhirHonk<Hasher>;
    const size_t target_log_n = static_cast<size_t>(state.range(0));
    const size_t log_n = actual_log_n(target_log_n, false);
    const WhirConfig config = Honk::make_config(log_n, WHIR_SECURITY_BITS, WHIR_LOG_INV_RATE);
    UltraCircuitBuilder builder = build_circuit(target_log_n, false);
    auto pk = Honk::create_proving_key(builder, config);
    const auto vk = pk.vk;
    const HonkProof proof = Honk::prove(pk);
    bool ok = true;
    for (auto _ : state) {
        ok = ok && Honk::verify(vk, config, proof);
    }
    if (!ok) {
        state.SkipWithError("WhirHonk verification failed");
    }
}

template <typename Hasher> void ligero_honk_prove(benchmark::State& state)
{
    using Honk = ligero::LigeroHonk<Hasher>;
    const size_t target_log_n = static_cast<size_t>(state.range(0));
    const size_t log_n = actual_log_n(target_log_n, false);
    const ligero::LigeroConfig config = Honk::make_config(log_n, WHIR_SECURITY_BITS, WHIR_LOG_INV_RATE);
    HonkProof proof;
    for (auto _ : state) {
        state.PauseTiming();
        UltraCircuitBuilder builder = build_circuit(target_log_n, false);
        auto pk = Honk::create_proving_key(builder, config);
        state.ResumeTiming();
        proof = Honk::prove(pk);
    }
    state.counters["proof_KiB"] = static_cast<double>(proof.size() * 32) / 1024.0;
    state.counters["log_n"] = static_cast<double>(log_n);
}

template <typename Hasher> void ligero_honk_verify(benchmark::State& state)
{
    using Honk = ligero::LigeroHonk<Hasher>;
    const size_t target_log_n = static_cast<size_t>(state.range(0));
    const size_t log_n = actual_log_n(target_log_n, false);
    const ligero::LigeroConfig config = Honk::make_config(log_n, WHIR_SECURITY_BITS, WHIR_LOG_INV_RATE);
    UltraCircuitBuilder builder = build_circuit(target_log_n, false);
    auto pk = Honk::create_proving_key(builder, config);
    const auto vk = pk.vk;
    const HonkProof proof = Honk::prove(pk);
    bool ok = true;
    for (auto _ : state) {
        ok = ok && Honk::verify(vk, config, proof);
    }
    if (!ok) {
        state.SkipWithError("LigeroHonk verification failed");
    }
}

void mercury_honk_prove(benchmark::State& state)
{
    srs::init_file_crs_factory(srs::bb_crs_path());
    using Honk = mercury::MercuryHonk;
    const size_t target_log_n = static_cast<size_t>(state.range(0));
    const size_t log_n = actual_log_n(target_log_n, false);
    const mercury::MercuryConfig config = Honk::make_config(log_n);
    HonkProof proof;
    for (auto _ : state) {
        state.PauseTiming();
        UltraCircuitBuilder builder = build_circuit(target_log_n, false);
        auto pk = Honk::create_proving_key(builder, config);
        state.ResumeTiming();
        proof = Honk::prove(pk);
    }
    state.counters["proof_KiB"] = static_cast<double>(proof.size() * 32) / 1024.0;
    state.counters["log_n"] = static_cast<double>(log_n);
}

void mercury_honk_verify(benchmark::State& state)
{
    srs::init_file_crs_factory(srs::bb_crs_path());
    using Honk = mercury::MercuryHonk;
    const size_t target_log_n = static_cast<size_t>(state.range(0));
    const size_t log_n = actual_log_n(target_log_n, false);
    const mercury::MercuryConfig config = Honk::make_config(log_n);
    UltraCircuitBuilder builder = build_circuit(target_log_n, false);
    auto pk = Honk::create_proving_key(builder, config);
    const auto vk = pk.vk;
    const HonkProof proof = Honk::prove(pk);
    bool ok = true;
    for (auto _ : state) {
        ok = ok && Honk::verify(vk, config, proof);
    }
    if (!ok) {
        state.SkipWithError("MercuryHonk verification failed");
    }
}

void ultra_honk_kzg_prove(benchmark::State& state)
{
    srs::init_file_crs_factory(srs::bb_crs_path());
    const size_t log_n = static_cast<size_t>(state.range(0));
    HonkProof proof;
    for (auto _ : state) {
        state.PauseTiming();
        UltraCircuitBuilder builder = build_circuit(log_n, true);
        auto prover_instance = std::make_shared<ProverInstance_<UltraFlavor>>(builder);
        auto verification_key = std::make_shared<UltraFlavor::VerificationKey>(prover_instance->get_precomputed());
        UltraProver_<UltraFlavor> prover(prover_instance, verification_key);
        state.ResumeTiming();
        proof = prover.construct_proof();
    }
    state.counters["proof_KiB"] = static_cast<double>(proof.size() * 32) / 1024.0;
    state.counters["log_n"] = static_cast<double>(actual_log_n(log_n, true));
}

void ultra_honk_kzg_verify(benchmark::State& state)
{
    srs::init_file_crs_factory(srs::bb_crs_path());
    const size_t log_n = static_cast<size_t>(state.range(0));
    UltraCircuitBuilder builder = build_circuit(log_n, true);
    auto prover_instance = std::make_shared<ProverInstance_<UltraFlavor>>(builder);
    auto verification_key = std::make_shared<UltraFlavor::VerificationKey>(prover_instance->get_precomputed());
    auto vk_and_hash = std::make_shared<UltraFlavor::VKAndHash>(verification_key);
    UltraProver_<UltraFlavor> prover(prover_instance, verification_key);
    const HonkProof proof = prover.construct_proof();
    bool ok = true;
    for (auto _ : state) {
        UltraVerifier_<UltraFlavor, DefaultIO> verifier(vk_and_hash);
        ok = ok && verifier.verify_proof(proof).result;
    }
    if (!ok) {
        state.SkipWithError("UltraHonk verification failed");
    }
}

BENCHMARK_TEMPLATE(whir_honk_prove, Blake3sMerkleHasher)
    ->Arg(14)
    ->Arg(16)
    ->Arg(18)
    ->Arg(20)
    ->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(whir_honk_prove, Poseidon2MerkleHasher)
    ->Arg(14)
    ->Arg(16)
    ->Arg(18)
    ->Arg(20)
    ->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(whir_honk_verify, Blake3sMerkleHasher)
    ->Arg(14)
    ->Arg(16)
    ->Arg(18)
    ->Arg(20)
    ->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(whir_honk_verify, Poseidon2MerkleHasher)
    ->Arg(14)
    ->Arg(16)
    ->Arg(18)
    ->Arg(20)
    ->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(ligero_honk_prove, Blake3sMerkleHasher)
    ->Arg(14)
    ->Arg(16)
    ->Arg(18)
    ->Arg(20)
    ->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(ligero_honk_prove, Poseidon2MerkleHasher)->Arg(14)->Arg(16)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(ligero_honk_verify, Blake3sMerkleHasher)
    ->Arg(14)
    ->Arg(16)
    ->Arg(18)
    ->Arg(20)
    ->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(ligero_honk_verify, Poseidon2MerkleHasher)->Arg(14)->Arg(16)->Unit(benchmark::kMillisecond);
BENCHMARK(mercury_honk_prove)->Arg(14)->Arg(16)->Arg(18)->Arg(20)->Unit(benchmark::kMillisecond);
BENCHMARK(mercury_honk_verify)->Arg(14)->Arg(16)->Arg(18)->Arg(20)->Unit(benchmark::kMillisecond);
BENCHMARK(ultra_honk_kzg_prove)->Arg(14)->Arg(16)->Arg(18)->Arg(20)->Unit(benchmark::kMillisecond);
BENCHMARK(ultra_honk_kzg_verify)->Arg(14)->Arg(16)->Arg(18)->Arg(20)->Unit(benchmark::kMillisecond);

} // namespace
} // namespace bb::whir

BENCHMARK_MAIN();
