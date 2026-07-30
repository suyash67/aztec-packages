#include "barretenberg/commitment_schemes/small_field/m31.hpp"

#include "barretenberg/crypto/blake3s/blake3s.hpp"
#include "barretenberg/ecc/curves/bn254/fr.hpp"
#include "barretenberg/numeric/uint256/uint256.hpp"

#include <benchmark/benchmark.h>
#include <cstring>

// Small-field vs BN254 Fr micro-benchmarks on the three prover-side workloads that dominate
// hash-based PCS backends: field multiplication streams, the sumcheck round kernel, and Blake3
// commitment hashing (where an M31 value is 4 bytes to Fr's 32). Single-threaded on purpose:
// these measure per-element work, not parallel scaling. See SMALL_FIELDS.md for the analysis.

namespace bb::small_field {
namespace {

constexpr size_t LOG_SIZE = 20;
constexpr size_t SIZE = size_t(1) << LOG_SIZE;

template <typename F> std::vector<F> random_vector(size_t size)
{
    std::vector<F> values(size);
    for (auto& value : values) {
        value = F::random_element();
    }
    return values;
}

template <> std::vector<fr> random_vector<fr>(size_t size)
{
    std::vector<fr> values(size);
    for (auto& value : values) {
        value = fr::random_element();
    }
    return values;
}

template <typename F> void mul_throughput(benchmark::State& state)
{
    const auto a = random_vector<F>(SIZE);
    const auto b = random_vector<F>(SIZE);
    for (auto _ : state) {
        F accumulator = a[0];
        for (size_t i = 0; i < SIZE; ++i) {
            accumulator += a[i] * b[i];
        }
        benchmark::DoNotOptimize(accumulator);
    }
    state.counters["Mmul/s"] =
        benchmark::Counter(static_cast<double>(SIZE) / 1e6, benchmark::Counter::kIsIterationInvariantRate);
}

/** @brief The degree-2 sumcheck round kernel: h(0), h(1), h(2) accumulation over paired tables. */
template <typename F> void sumcheck_round_kernel(benchmark::State& state)
{
    const auto f = random_vector<F>(SIZE);
    const auto w = random_vector<F>(SIZE);
    for (auto _ : state) {
        F h0{};
        F h1{};
        F h2{};
        for (size_t t = 0; t < SIZE / 2; ++t) {
            const F& f0 = f[2 * t];
            const F& f1 = f[2 * t + 1];
            const F& w0 = w[2 * t];
            const F& w1 = w[2 * t + 1];
            h0 += f0 * w0;
            h1 += f1 * w1;
            h2 += (f1 + f1 - f0) * (w1 + w1 - w0);
        }
        benchmark::DoNotOptimize(h0);
        benchmark::DoNotOptimize(h1);
        benchmark::DoNotOptimize(h2);
    }
}

/** @brief The mixed round: M31 witness table against a QM31 (challenge-field) weight table. */
void sumcheck_round_kernel_mixed(benchmark::State& state)
{
    const auto f = random_vector<m31>(SIZE);
    const auto w = random_vector<qm31>(SIZE);
    for (auto _ : state) {
        qm31 h0;
        qm31 h1;
        qm31 h2;
        for (size_t t = 0; t < SIZE / 2; ++t) {
            const m31& f0 = f[2 * t];
            const m31& f1 = f[2 * t + 1];
            const qm31& w0 = w[2 * t];
            const qm31& w1 = w[2 * t + 1];
            h0 += w0.scale(f0);
            h1 += w1.scale(f1);
            h2 += (w1 + w1 - w0).scale(f1 + f1 - f0);
        }
        benchmark::DoNotOptimize(h0);
        benchmark::DoNotOptimize(h1);
        benchmark::DoNotOptimize(h2);
    }
}

/** @brief Blake3 over the byte serialization of SIZE field values, in 768-byte chunks. */
void blake3_hash_values(benchmark::State& state, size_t bytes_per_value)
{
    std::vector<uint8_t> buffer(SIZE * bytes_per_value);
    for (size_t i = 0; i < buffer.size(); ++i) {
        buffer[i] = static_cast<uint8_t>(i * 0x9e3779b9U >> 24);
    }
    for (auto _ : state) {
        std::array<uint8_t, 32> digest{};
        std::vector<uint8_t> chunk(768 + 32);
        for (size_t offset = 0; offset < buffer.size(); offset += 768) {
            const size_t len = std::min<size_t>(768, buffer.size() - offset);
            chunk.assign(digest.begin(), digest.end());
            chunk.insert(chunk.end(),
                         buffer.begin() + static_cast<std::ptrdiff_t>(offset),
                         buffer.begin() + static_cast<std::ptrdiff_t>(offset + len));
            const auto out = blake3::blake3s(chunk);
            std::memcpy(digest.data(), out.data(), 32);
        }
        benchmark::DoNotOptimize(digest);
    }
    state.counters["MiB"] = static_cast<double>(buffer.size()) / (1024.0 * 1024.0);
}

void blake3_hash_m31_values(benchmark::State& state)
{
    blake3_hash_values(state, sizeof(uint32_t));
}
void blake3_hash_fr_values(benchmark::State& state)
{
    blake3_hash_values(state, 32);
}

BENCHMARK_TEMPLATE(mul_throughput, m31)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(mul_throughput, qm31)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(mul_throughput, fr)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(sumcheck_round_kernel, m31)->Unit(benchmark::kMillisecond);
BENCHMARK(sumcheck_round_kernel_mixed)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(sumcheck_round_kernel, fr)->Unit(benchmark::kMillisecond);
BENCHMARK(blake3_hash_m31_values)->Unit(benchmark::kMillisecond);
BENCHMARK(blake3_hash_fr_values)->Unit(benchmark::kMillisecond);

} // namespace
} // namespace bb::small_field

BENCHMARK_MAIN();
