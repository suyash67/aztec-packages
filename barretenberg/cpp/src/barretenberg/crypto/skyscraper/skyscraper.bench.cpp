#include "barretenberg/crypto/skyscraper/skyscraper.hpp"
#include "barretenberg/crypto/poseidon2/poseidon2.hpp"

#include "barretenberg/common/log.hpp"
#include "barretenberg/ecc/curves/bn254/fr.hpp"

#include <chrono>
#include <iomanip>
#include <iostream>
#include <vector>

/**
 * @brief Throughput of the Skyscraper compression, which is what a WHIR prover spends its time on.
 * @details A round-0 Merkle tree hashes millions of these, so the figure that matters is
 * compressions per second in a tight loop over data already in cache.
 */
int main()
{
    using namespace bb;

    constexpr size_t NUM_INPUTS = 1 << 12;
    constexpr size_t NUM_PASSES = 64;

    std::vector<fr> values(NUM_INPUTS);
    for (fr& value : values) {
        value = fr::random_element();
    }

    // Warm the cache and keep the optimiser honest by consuming the result.
    fr sink = fr::zero();
    for (size_t i = 1; i < NUM_INPUTS; ++i) {
        sink += crypto::skyscraper::compress(values[i - 1], values[i]);
    }

    const auto start = std::chrono::steady_clock::now();
    for (size_t pass = 0; pass < NUM_PASSES; ++pass) {
        for (size_t i = 1; i < NUM_INPUTS; ++i) {
            sink += crypto::skyscraper::compress(values[i - 1], values[i]);
        }
    }
    const auto ns =
        std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count();

    const size_t compressions = NUM_PASSES * (NUM_INPUTS - 1);
    std::cout << "skyscraper compress: " << std::fixed << std::setprecision(1)
              << static_cast<double>(ns) / static_cast<double>(compressions) << " ns each, "
              << static_cast<double>(compressions) * 1e9 / static_cast<double>(ns) / 1e6 << " M/s\n";
    // bb's field stream operator leaves the fill character set to '0' for its hex digits.
    std::cout << std::setfill(' ');

    // Where the time goes. A compression is 14 squaring half-rounds (a square and a multiply by
    // sigma^-1 each) and 4 bars (a Montgomery round trip and 32 S-box lookups each).
    const auto time_loop = [&](const char* label, size_t ops_per_compression, auto&& body) {
        fr acc = values[0];
        const auto begin = std::chrono::steady_clock::now();
        for (size_t pass = 0; pass < NUM_PASSES; ++pass) {
            for (size_t i = 1; i < NUM_INPUTS; ++i) {
                body(acc, values[i]);
            }
        }
        const auto elapsed =
            std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - begin).count();
        const double each = static_cast<double>(elapsed) / static_cast<double>(compressions);
        std::cout << "  " << std::left << std::setw(28) << label << std::right << std::setw(8) << std::fixed
                  << std::setprecision(2) << each << " ns/op   x" << ops_per_compression << " = " << std::setw(7)
                  << each * static_cast<double>(ops_per_compression) << " ns of a compression\n";
        sink += acc;
    };

    static const fr sigma_inv = crypto::skyscraper::sigma_inv();
    std::cout << "\ncomponents\n";
    time_loop("square", 14, [](fr& acc, const fr& v) { acc = acc.sqr() + v; });
    time_loop("square + mul by sigma_inv", 14, [](fr& acc, const fr& v) { acc = acc.sqr() * sigma_inv + v; });
    time_loop("montgomery round trip", 4, [](fr& acc, const fr& v) { acc = fr(uint256_t(acc)) + v; });
    time_loop("bar (round trip + sbox)", 4, [](fr& acc, const fr& v) {
        const uint256_t c(acc);
        acc =
            fr(uint256_t(crypto::skyscraper::sbox(static_cast<uint8_t>(c.data[2])), c.data[3], c.data[0], c.data[1])) +
            v;
    });
    // Fiat-Shamir: what a transcript hash costs per absorbed field element. bb's transcript is
    // Poseidon2 (rate 3, one permutation per three elements); Skyscraper's rate-1 duplex sponge is
    // one permutation each. Per element is therefore the fair unit, not per permutation.
    std::cout << "\nFiat-Shamir transcript hash, per absorbed field element\n";
    for (const size_t width : { size_t(8), size_t(64), size_t(512) }) {
        std::vector<fr> buffer(values.begin(), values.begin() + static_cast<long>(width));
        const size_t iterations = 1u << 14;

        auto begin = std::chrono::steady_clock::now();
        for (size_t i = 0; i < iterations; ++i) {
            sink += crypto::Poseidon2<crypto::Poseidon2Bn254ScalarFieldParams>::hash(buffer);
        }
        const double poseidon2_ns = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(
                                                            std::chrono::steady_clock::now() - begin)
                                                            .count()) /
                                    static_cast<double>(iterations * width);

        begin = std::chrono::steady_clock::now();
        for (size_t i = 0; i < iterations; ++i) {
            sink += crypto::skyscraper::hash(buffer);
        }
        const double skyscraper_ns = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(
                                                             std::chrono::steady_clock::now() - begin)
                                                             .count()) /
                                     static_cast<double>(iterations * width);

        std::cout << "  " << std::setw(4) << width << " elements: poseidon2 " << std::setw(7) << std::fixed
                  << std::setprecision(1) << poseidon2_ns << " ns/elt   skyscraper " << std::setw(7) << skyscraper_ns
                  << " ns/elt   ratio " << std::setprecision(2) << poseidon2_ns / skyscraper_ns << "x\n";
    }

    std::cout << std::setfill(' ') << "\n(sink " << sink << ")\n";
    return 0;
}
