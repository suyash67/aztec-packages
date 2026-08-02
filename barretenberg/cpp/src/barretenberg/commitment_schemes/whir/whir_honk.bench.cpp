#include "barretenberg/commitment_schemes/whir/whir_honk.hpp"
#include "barretenberg/commitment_schemes/dory/dory_honk.hpp"
#include "barretenberg/commitment_schemes/hyrax/hyrax_honk.hpp"
#include "barretenberg/commitment_schemes/kzh/kzh3_honk.hpp"
#include "barretenberg/commitment_schemes/kzh/kzh_honk.hpp"
#include "barretenberg/commitment_schemes/ligero/ligero_honk.hpp"
#include "barretenberg/commitment_schemes/mercury/mercury_honk.hpp"
#include "barretenberg/commitment_schemes/pedersen_ipa/ipa_honk.hpp"
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
    if (target_log_n >= 14) {
        MockCircuits::add_lookup_gates(builder, 2);
    } else {
        // UINT32_XOR's 4096-row table cannot fit a 2^12 trace; use the 8-row dummy multitable so
        // the lookup argument stays active at small sizes.
        const fr a_value(1);
        const auto a_idx = builder.add_variable(a_value);
        const auto accumulators =
            plookup::get_lookup_accumulators(plookup::MultiTableId::HONK_DUMMY_MULTI, a_value, fr(0), true);
        builder.create_gates_from_plookup_accumulators(
            plookup::MultiTableId::HONK_DUMMY_MULTI, accumulators, a_idx, std::nullopt);
    }
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

// WHIR under the repaired up-to-capacity conjectures (Crites-Stewart, eprint 2025/2046) instead of
// the disproved capacity regime: same protocol, one extra query per aggressive-rate round.
template <typename Hasher> void whir_honk_repaired_prove(benchmark::State& state)
{
    using Honk = WhirHonk<Hasher>;
    const size_t target_log_n = static_cast<size_t>(state.range(0));
    const size_t log_n = actual_log_n(target_log_n, false);
    const WhirConfig config =
        WhirConfig::create(log_n, WHIR_SECURITY_BITS, WHIR_LOG_INV_RATE, 4, 4, WhirSoundness::REPAIRED_LIST);
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

template <typename Hasher> void whir_honk_repaired_verify(benchmark::State& state)
{
    using Honk = WhirHonk<Hasher>;
    const size_t target_log_n = static_cast<size_t>(state.range(0));
    const size_t log_n = actual_log_n(target_log_n, false);
    const WhirConfig config =
        WhirConfig::create(log_n, WHIR_SECURITY_BITS, WHIR_LOG_INV_RATE, 4, 4, WhirSoundness::REPAIRED_LIST);
    UltraCircuitBuilder builder = build_circuit(target_log_n, false);
    auto pk = Honk::create_proving_key(builder, config);
    const auto vk = pk.vk;
    const HonkProof proof = Honk::prove(pk);
    bool ok = true;
    for (auto _ : state) {
        ok = ok && Honk::verify(vk, config, proof);
    }
    if (!ok) {
        state.SkipWithError("WhirHonk repaired-preset verification failed");
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

void ipa_honk_prove(benchmark::State& state)
{
    using Honk = pedersen_ipa::IpaHonk;
    const size_t target_log_n = static_cast<size_t>(state.range(0));
    const size_t log_n = actual_log_n(target_log_n, false);
    const pedersen_ipa::IpaConfig config = Honk::make_config(log_n);
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

void ipa_honk_verify(benchmark::State& state)
{
    using Honk = pedersen_ipa::IpaHonk;
    const size_t target_log_n = static_cast<size_t>(state.range(0));
    const size_t log_n = actual_log_n(target_log_n, false);
    const pedersen_ipa::IpaConfig config = Honk::make_config(log_n);
    UltraCircuitBuilder builder = build_circuit(target_log_n, false);
    auto pk = Honk::create_proving_key(builder, config);
    const auto vk = pk.vk;
    const HonkProof proof = Honk::prove(pk);
    bool ok = true;
    for (auto _ : state) {
        ok = ok && Honk::verify(vk, config, proof);
    }
    if (!ok) {
        state.SkipWithError("IpaHonk verification failed");
    }
}

void hyrax_honk_prove(benchmark::State& state)
{
    using Honk = hyrax::HyraxHonk;
    const size_t target_log_n = static_cast<size_t>(state.range(0));
    const size_t log_n = actual_log_n(target_log_n, false);
    const hyrax::HyraxConfig config = Honk::make_config(log_n);
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

void hyrax_honk_verify(benchmark::State& state)
{
    using Honk = hyrax::HyraxHonk;
    const size_t target_log_n = static_cast<size_t>(state.range(0));
    const size_t log_n = actual_log_n(target_log_n, false);
    const hyrax::HyraxConfig config = Honk::make_config(log_n);
    UltraCircuitBuilder builder = build_circuit(target_log_n, false);
    auto pk = Honk::create_proving_key(builder, config);
    const auto vk = pk.vk;
    const HonkProof proof = Honk::prove(pk);
    bool ok = true;
    for (auto _ : state) {
        ok = ok && Honk::verify(vk, config, proof);
    }
    if (!ok) {
        state.SkipWithError("HyraxHonk verification failed");
    }
}

void kzh_honk_prove(benchmark::State& state)
{
    using Honk = kzh::KzhHonk;
    const size_t target_log_n = static_cast<size_t>(state.range(0));
    const size_t log_n = actual_log_n(target_log_n, false);
    const kzh::KzhConfig config = Honk::make_config(log_n);
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

void kzh_honk_verify(benchmark::State& state)
{
    using Honk = kzh::KzhHonk;
    const size_t target_log_n = static_cast<size_t>(state.range(0));
    const size_t log_n = actual_log_n(target_log_n, false);
    const kzh::KzhConfig config = Honk::make_config(log_n);
    UltraCircuitBuilder builder = build_circuit(target_log_n, false);
    auto pk = Honk::create_proving_key(builder, config);
    const auto vk = pk.vk;
    const HonkProof proof = Honk::prove(pk);
    bool ok = true;
    for (auto _ : state) {
        ok = ok && Honk::verify(vk, config, proof);
    }
    if (!ok) {
        state.SkipWithError("KzhHonk verification failed");
    }
}

void kzh3_honk_prove(benchmark::State& state)
{
    using Honk = kzh3::Kzh3Honk;
    const size_t target_log_n = static_cast<size_t>(state.range(0));
    const size_t log_n = actual_log_n(target_log_n, false);
    const kzh3::Kzh3Config config = Honk::make_config(log_n);
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

void kzh3_honk_verify(benchmark::State& state)
{
    using Honk = kzh3::Kzh3Honk;
    const size_t target_log_n = static_cast<size_t>(state.range(0));
    const size_t log_n = actual_log_n(target_log_n, false);
    const kzh3::Kzh3Config config = Honk::make_config(log_n);
    UltraCircuitBuilder builder = build_circuit(target_log_n, false);
    auto pk = Honk::create_proving_key(builder, config);
    const auto vk = pk.vk;
    const HonkProof proof = Honk::prove(pk);
    bool ok = true;
    for (auto _ : state) {
        ok = ok && Honk::verify(vk, config, proof);
    }
    if (!ok) {
        state.SkipWithError("Kzh3Honk verification failed");
    }
}

void dory_honk_prove(benchmark::State& state)
{
    using Honk = dory::DoryHonk;
    const size_t target_log_n = static_cast<size_t>(state.range(0));
    const size_t log_n = actual_log_n(target_log_n, false);
    const dory::DoryConfig config = Honk::make_config(log_n);
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

void dory_honk_verify(benchmark::State& state)
{
    using Honk = dory::DoryHonk;
    const size_t target_log_n = static_cast<size_t>(state.range(0));
    const size_t log_n = actual_log_n(target_log_n, false);
    const dory::DoryConfig config = Honk::make_config(log_n);
    UltraCircuitBuilder builder = build_circuit(target_log_n, false);
    auto pk = Honk::create_proving_key(builder, config);
    const auto vk = pk.vk;
    const HonkProof proof = Honk::prove(pk);
    bool ok = true;
    for (auto _ : state) {
        ok = ok && Honk::verify(vk, config, proof);
    }
    if (!ok) {
        state.SkipWithError("DoryHonk verification failed");
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

BENCHMARK_TEMPLATE(whir_honk_prove, Blake3sMerkleHasher)->DenseRange(12, 20, 2)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(whir_honk_prove, Poseidon2MerkleHasher)->DenseRange(12, 20, 2)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(whir_honk_verify, Blake3sMerkleHasher)->DenseRange(12, 20, 2)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(whir_honk_verify, Poseidon2MerkleHasher)->DenseRange(12, 20, 2)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(whir_honk_repaired_prove, Blake3sMerkleHasher)->DenseRange(12, 20, 2)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(whir_honk_repaired_prove, Poseidon2MerkleHasher)
    ->DenseRange(12, 20, 2)
    ->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(whir_honk_repaired_verify, Blake3sMerkleHasher)
    ->DenseRange(12, 20, 2)
    ->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(whir_honk_repaired_verify, Poseidon2MerkleHasher)
    ->DenseRange(12, 20, 2)
    ->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(ligero_honk_prove, Blake3sMerkleHasher)->DenseRange(12, 20, 2)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(ligero_honk_prove, Poseidon2MerkleHasher)->DenseRange(12, 20, 2)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(ligero_honk_verify, Blake3sMerkleHasher)->DenseRange(12, 20, 2)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(ligero_honk_verify, Poseidon2MerkleHasher)->DenseRange(12, 20, 2)->Unit(benchmark::kMillisecond);
BENCHMARK(mercury_honk_prove)->DenseRange(12, 20, 2)->Unit(benchmark::kMillisecond);
BENCHMARK(mercury_honk_verify)->DenseRange(12, 20, 2)->Unit(benchmark::kMillisecond);
BENCHMARK(ipa_honk_prove)->DenseRange(12, 20, 2)->Unit(benchmark::kMillisecond);
BENCHMARK(ipa_honk_verify)->DenseRange(12, 20, 2)->Unit(benchmark::kMillisecond);
BENCHMARK(hyrax_honk_prove)->DenseRange(12, 20, 2)->Unit(benchmark::kMillisecond);
BENCHMARK(hyrax_honk_verify)->DenseRange(12, 20, 2)->Unit(benchmark::kMillisecond);
BENCHMARK(kzh_honk_prove)->DenseRange(12, 20, 2)->Unit(benchmark::kMillisecond);
BENCHMARK(kzh_honk_verify)->DenseRange(12, 20, 2)->Unit(benchmark::kMillisecond);
BENCHMARK(kzh3_honk_prove)->DenseRange(12, 20, 2)->Unit(benchmark::kMillisecond);
BENCHMARK(kzh3_honk_verify)->DenseRange(12, 20, 2)->Unit(benchmark::kMillisecond);
BENCHMARK(dory_honk_prove)->DenseRange(12, 20, 2)->Unit(benchmark::kMillisecond);
BENCHMARK(dory_honk_verify)->DenseRange(12, 20, 2)->Unit(benchmark::kMillisecond);
BENCHMARK(ultra_honk_kzg_prove)->DenseRange(12, 20, 2)->Unit(benchmark::kMillisecond);
BENCHMARK(ultra_honk_kzg_verify)->DenseRange(12, 20, 2)->Unit(benchmark::kMillisecond);

} // namespace
} // namespace bb::whir

BENCHMARK_MAIN();
