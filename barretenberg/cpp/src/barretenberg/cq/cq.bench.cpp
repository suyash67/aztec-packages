#include "barretenberg/cq/cq_key.hpp"
#include "barretenberg/cq/cq_prover.hpp"
#include "barretenberg/cq/cq_trusted_setup.hpp"
#include "barretenberg/cq/cq_verifier.hpp"

#include <benchmark/benchmark.h>

namespace bb::cq {

namespace {

struct BenchSetup {
    TestSrs srs;
    std::vector<fr> table;
    CqProvingKey pk;
    CommitmentKey<curve::BN254> ck;
    std::vector<fr> lookups;
};

// One preprocessed table per size, shared across benchmark iterations.
BenchSetup& get_setup(size_t log_table_size, size_t num_lookups)
{
    static std::map<std::pair<size_t, size_t>, std::unique_ptr<BenchSetup>> setups;
    const auto key = std::make_pair(log_table_size, num_lookups);
    auto it = setups.find(key);
    if (it == setups.end()) {
        const size_t table_size = 1UL << log_table_size;
        auto setup = std::make_unique<BenchSetup>();
        setup->srs = TestSrs::create(table_size, table_size + 1);
        setup->table.resize(table_size);
        for (size_t i = 0; i < table_size; ++i) {
            setup->table[i] = fr(i * i + 1);
        }
        setup->pk = CqProvingKey::create(setup->table, setup->srs.g1_powers, setup->srs.g2_powers);
        setup->ck = setup->srs.create_commitment_key();
        setup->lookups.resize(num_lookups);
        for (size_t j = 0; j < num_lookups; ++j) {
            setup->lookups[j] = setup->table[(j * 7919 + 13) % table_size];
        }
        it = setups.emplace(key, std::move(setup)).first;
    }
    return *it->second;
}

/**
 * @brief Proving time at fixed lookup count against tables of growing size. The point of cq: this must stay flat
 * as the table grows.
 */
void bench_cq_prove(benchmark::State& state)
{
    const size_t log_table_size = static_cast<size_t>(state.range(0));
    const size_t num_lookups = static_cast<size_t>(state.range(1));
    BenchSetup& setup = get_setup(log_table_size, num_lookups);
    CqProver prover(setup.pk, setup.ck);
    for (auto _ : state) {
        benchmark::DoNotOptimize(prover.construct_proof(setup.lookups));
    }
}

/**
 * @brief One-time table preprocessing cost (Feist-Khovratovich cached quotients, Lagrange caches, G2 MSM).
 */
void bench_cq_preprocess(benchmark::State& state)
{
    const size_t table_size = 1UL << static_cast<size_t>(state.range(0));
    TestSrs srs = TestSrs::create(table_size, table_size + 1);
    std::vector<fr> table(table_size);
    for (size_t i = 0; i < table_size; ++i) {
        table[i] = fr(i * i + 1);
    }
    for (auto _ : state) {
        benchmark::DoNotOptimize(CqProvingKey::create(table, srs.g1_powers, srs.g2_powers));
    }
}

void bench_cq_verify(benchmark::State& state)
{
    const size_t log_table_size = static_cast<size_t>(state.range(0));
    const size_t num_lookups = static_cast<size_t>(state.range(1));
    BenchSetup& setup = get_setup(log_table_size, num_lookups);
    CqProver prover(setup.pk, setup.ck);
    const auto proof = prover.construct_proof(setup.lookups);
    CqVerifier verifier(setup.pk.verification_key);
    for (auto _ : state) {
        benchmark::DoNotOptimize(verifier.verify_proof(proof));
    }
}

} // namespace

BENCHMARK(bench_cq_prove)
    ->ArgsProduct({ { 10, 12, 14, 16 }, { 256, 1024 } })
    ->Unit(benchmark::kMillisecond)
    ->Iterations(3);
BENCHMARK(bench_cq_preprocess)->Arg(10)->Arg(12)->Arg(14)->Arg(16)->Unit(benchmark::kMillisecond)->Iterations(1);
BENCHMARK(bench_cq_verify)->Args({ 14, 1024 })->Unit(benchmark::kMillisecond)->Iterations(5);

} // namespace bb::cq

BENCHMARK_MAIN();
