#pragma once

#include "barretenberg/commitment_schemes/whir/rs_code.hpp"
#include "barretenberg/common/assert.hpp"
#include "barretenberg/common/thread.hpp"
#include "barretenberg/crypto/blake3s/blake3s.hpp"
#include "barretenberg/ecc/curves/bn254/fr.hpp"
#include "barretenberg/numeric/uint256/uint256.hpp"

#include <cmath>
#include <cstring>
#include <span>
#include <vector>

namespace bb::bolt {

/**
 * @brief Parameters of Bolt's sketched code (eprint 2026/310 §3, §6.2).
 *
 * @details The sketched code is `C_H(x) = (x, C(Hx))`: a *systematic* code that applies the
 * expensive base code `C` only to a compressed sketch `Hx`, where `H` is a sparse parity-check
 * matrix. `H` comes from the random `(j, k, q)`-LDPC ensemble — each of the `n` variable nodes has
 * degree `j`, each of the `nj/k` check nodes degree `k`, so `H` is `(j/k)n x n` and `alpha = j/k`.
 *
 * `gamma` is the certified minimum relative distance of that LDPC ensemble, which is exactly the
 * parity-check property the sketched code needs (`Hx != 0` for every `x` of weight below
 * `gamma * n`). It follows from Yang et al.'s bound as quoted in the paper's Theorem 6.5: `gamma` is
 * the root of `omega_{q,j,k}`, computed by `ldpc_distance_bound` below.
 *
 * ### The field size matters, and it is why the presets differ
 *
 * `omega` contains `H_q(x)`, whose leading term is `x * ln(q-1)`. That term grows with the field, so
 * a fixed column degree certifies less and less distance as `q` grows: the column degree has to
 * scale with `ln q`. Concretely, the paper's own instantiation `(j=16, k=128)` at `q = 2^32` gives
 * `gamma = 0.0941`, but the identical parameters at BN254 (`ln q ~ 176`) give `gamma = 0.00006`,
 * which is useless. Recovering a comparable distance at BN254 needs `j = 64`.
 *
 * This is the concrete reason Bolt's headline cost — `(3+eps)N` field *additions* — does not carry
 * to BN254. That figure is the Boolean-`H` instantiation of the paper's §6.1, where `j = 3` and the
 * entries are 0/1 so the sketch costs only additions. Over a 254-bit prime the same theorem requires
 * random field weights at `j = 64`, i.e. 64 *multiplications* per message symbol. What does survive
 * is the distance: `gamma = 0.15` against Brakedown's 0.07, and distance is the binding constraint
 * on proof size (see `../PCS_CANDIDATES.md` §7 and `../brakedown/README.md` §5).
 */
struct BoltParams {
    size_t column_degree; // j
    size_t check_degree;  // k
    double gamma;         // certified relative distance of the LDPC ensemble

    double alpha() const { return static_cast<double>(column_degree) / static_cast<double>(check_degree); }

    /** @brief BN254: the smallest even column degree that recovers a distance above Brakedown's. */
    static BoltParams bn254() { return { 64, 256, 0.1497 }; }

    /** @brief The paper's own instantiation at q = 2^32, kept as the validation reference. */
    static BoltParams reference_2p32() { return { 16, 128, 0.0941 }; }
};

/**
 * @brief The root of `omega_{q,j,k}`, i.e. the certified relative distance of the ensemble.
 * @details Theorem 6.5 of the paper (from Yang et al.):
 *
 *     omega(x)      = H_q(x) + (j/k) * (delta_qk(x) - ln q)
 *     H_q(x)        = x ln(1/x) + (1-x) ln(1/(1-x)) + x ln(q-1)
 *     delta_qk(x)   = inf_xhat { k D(x||xhat) + rho_qk(xhat) }
 *     rho_qk(xhat)  = ln[ 1 + (q-1) (1 - q xhat/(q-1))^k ]
 *     D(x||xhat)    = x ln(x/xhat) + (1-x) ln((1-x)/(1-xhat))
 *
 * `omega` is negative on `(0, x0)` and positive above, so the root is found by bisection. Returns 0
 * when no root exists below `max_distance`, which is the "this column degree certifies nothing at
 * this field size" case.
 *
 * @param log_field_size natural logarithm of `q`
 */
inline double ldpc_distance_bound(double log_field_size, size_t j, size_t k, double max_distance = 0.49)
{
    const double ln_q = log_field_size;
    // (q-1)/q and ln(q-1) are indistinguishable from 1 and ln q at every field size of interest.
    const auto entropy = [&](double x) {
        return (x * std::log(1.0 / x)) + ((1.0 - x) * std::log(1.0 / (1.0 - x))) + (x * ln_q);
    };
    const auto rho = [&](double xhat) {
        // ln[1 + (q-1)(1-xhat)^k] = ln[1 + exp(ln q + k ln(1-xhat))], computed stably.
        const double exponent = ln_q + (static_cast<double>(k) * std::log(1.0 - xhat));
        return exponent > 30.0 ? exponent : std::log1p(std::exp(exponent));
    };
    const auto divergence = [](double x, double xhat) {
        return (x * std::log(x / xhat)) + ((1.0 - x) * std::log((1.0 - x) / (1.0 - xhat)));
    };
    const auto delta_qk = [&](double x) {
        double lo = 1e-12;
        double hi = 1.0 - 1e-12;
        double best = std::numeric_limits<double>::infinity();
        for (size_t pass = 0; pass < 5; ++pass) {
            constexpr size_t SAMPLES = 1500;
            size_t best_index = 0;
            best = std::numeric_limits<double>::infinity();
            for (size_t i = 0; i <= SAMPLES; ++i) {
                const double xhat = lo + ((hi - lo) * static_cast<double>(i) / static_cast<double>(SAMPLES));
                const double value = (static_cast<double>(k) * divergence(x, xhat)) + rho(xhat);
                if (value < best) {
                    best = value;
                    best_index = i;
                }
            }
            const double step = (hi - lo) / static_cast<double>(SAMPLES);
            const double centre = lo + (step * static_cast<double>(best_index));
            lo = std::max(1e-12, centre - (2.0 * step));
            hi = std::min(1.0 - 1e-12, centre + (2.0 * step));
        }
        return best;
    };
    const auto omega = [&](double x) {
        return entropy(x) + ((static_cast<double>(j) / static_cast<double>(k)) * (delta_qk(x) - ln_q));
    };

    double lo = 1e-7;
    double hi = max_distance;
    if (omega(lo) > 0.0) {
        return 0.0; // this column degree certifies nothing at this field size
    }
    for (size_t iteration = 0; iteration < 50; ++iteration) {
        const double mid = (lo + hi) / 2.0;
        if (omega(mid) < 0.0) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return (lo + hi) / 2.0;
}

namespace detail {

/** @brief Deterministic Blake3 stream, so both parties derive the identical sketch matrix. */
class SeededPrng {
  public:
    SeededPrng(uint64_t seed, uint64_t domain)
        : seed_(seed)
        , domain_(domain)
    {}

    uint64_t next()
    {
        if (cursor_ == 0) {
            std::vector<uint8_t> input(24);
            std::memcpy(input.data(), &seed_, 8);
            std::memcpy(input.data() + 8, &domain_, 8);
            std::memcpy(input.data() + 16, &counter_, 8);
            block_ = blake3::blake3s(input);
            ++counter_;
            cursor_ = 4;
        }
        uint64_t value = 0;
        std::memcpy(&value, block_.data() + (8 * (4 - cursor_)), 8);
        --cursor_;
        return value;
    }

    size_t next_index(size_t bound) { return static_cast<size_t>(next() % bound); }

    fr next_nonzero_field()
    {
        // Locals first: constructor argument evaluation order is unspecified, and a
        // compiler-dependent matrix would break prover/verifier agreement.
        const uint64_t low = next();
        const uint64_t high = next();
        const fr candidate{ uint256_t(low, high, 0, 0) };
        return candidate.is_zero() ? fr::one() : candidate;
    }

  private:
    uint64_t seed_;
    uint64_t domain_;
    uint64_t counter_ = 0;
    size_t cursor_ = 0;
    std::vector<uint8_t> block_;
};

} // namespace detail

/**
 * @brief Bolt's sketched code over BN254: `C_H(x) = (x, C(Hx))` with `C` Reed-Solomon.
 *
 * @details `H` is sampled from the `(j, k)`-LDPC ensemble exactly as the paper describes it: a
 * uniformly random bijection between the `nj` variable sockets and the `nj` check sockets, with an
 * independent random nonzero weight per check socket. Stored compressed by row (one entry per check
 * node), because the product needed is `H x` — a gather per output row, which parallelizes.
 *
 * The base code is Reed-Solomon on the length-`alpha n` sketch. Bolt itself allows any base code;
 * RS is the natural pick here because the sketch is short, so its `O(m log m)` cost is a small term,
 * and because it brings a large distance to the second piece of the piecewise guarantee.
 *
 * The codeword is `(x, C(Hx))`: systematic in `[0, n)`, sketch in `[n, n + alpha n / rho)`. Those two
 * stretches carry *separate* distance guarantees — `gamma` and the base code's `delta` — which is
 * what `bolt_honk.hpp` exploits by querying them independently.
 */
class BoltCode {
  public:
    BoltCode(size_t message_length, const BoltParams& params, size_t log_inv_rate = 2, uint64_t seed = 0x626f6c74U)
        : params_(params)
        , message_length_(message_length)
    {
        BB_ASSERT_GT(message_length, size_t(0));
        BB_ASSERT_GT(params.gamma, 0.0, "this column degree certifies no distance at this field size");

        // One check node per k sockets; alpha = j/k.
        sketch_length_ = (message_length * params.column_degree) / params.check_degree;
        BB_ASSERT_GT(sketch_length_, size_t(0), "message too short for this sketch ratio");
        sketch_codeword_length_ = sketch_length_ << log_inv_rate;
        codeword_length_ = message_length + sketch_codeword_length_;

        sample_sketch_matrix(seed);
        domains_.get(sketch_codeword_length_);
    }

    size_t message_length() const { return message_length_; }
    size_t codeword_length() const { return codeword_length_; }
    size_t sketch_begin() const { return message_length_; }
    size_t sketch_codeword_length() const { return sketch_codeword_length_; }
    const BoltParams& params() const { return params_; }

    /** @brief Field multiplications one encoding costs: the sketch plus the base code's FFT. */
    size_t multiplications() const
    {
        const size_t sketch = message_length_ * params_.column_degree;
        const size_t fft = sketch_codeword_length_ * numeric::get_msb(sketch_codeword_length_);
        return sketch + fft;
    }

    std::vector<fr> encode(std::span<const fr> message) const
    {
        BB_ASSERT_EQ(message.size(), message_length_, "message length does not match the code");
        std::vector<fr> codeword(codeword_length_, fr::zero());
        std::copy(message.begin(), message.end(), codeword.begin());

        // sketch = H * message, one independent gather per check node.
        std::vector<fr> sketch(sketch_length_, fr::zero());
        parallel_for_range(sketch_length_, [&](size_t start, size_t end) {
            for (size_t row = start; row < end; ++row) {
                fr acc = fr::zero();
                for (size_t e = row_offsets_[row]; e < row_offsets_[row + 1]; ++e) {
                    acc += message[column_indices_[e]] * weights_[e];
                }
                sketch[row] = acc;
            }
        });

        const std::vector<fr> encoded_sketch =
            whir::rs_encode(sketch, domains_.get(sketch_codeword_length_), domains_.round_roots());
        std::copy(encoded_sketch.begin(),
                  encoded_sketch.end(),
                  codeword.begin() + static_cast<std::ptrdiff_t>(message_length_));
        return codeword;
    }

  private:
    /** @brief The random socket bijection of the (j, k) ensemble, stored row-major over check nodes. */
    void sample_sketch_matrix(uint64_t seed)
    {
        const size_t j = params_.column_degree;
        const size_t total_sockets = message_length_ * j;
        detail::SeededPrng prng(seed, 0x424f4c5400000000ULL + message_length_);

        // A uniformly random bijection variable-socket -> check-socket, by shuffling the check
        // sockets and reading them off in variable order.
        std::vector<uint32_t> permutation(total_sockets);
        for (size_t i = 0; i < total_sockets; ++i) {
            permutation[i] = static_cast<uint32_t>(i);
        }
        for (size_t i = total_sockets; i-- > 1;) {
            std::swap(permutation[i], permutation[prng.next_index(i + 1)]);
        }

        // Check socket s belongs to check node s / k.
        std::vector<size_t> counts(sketch_length_ + 1, 0);
        std::vector<uint32_t> node_of_socket(total_sockets);
        for (size_t variable = 0; variable < message_length_; ++variable) {
            for (size_t d = 0; d < j; ++d) {
                const size_t socket = permutation[(variable * j) + d];
                const size_t node = std::min(socket / params_.check_degree, sketch_length_ - 1);
                node_of_socket[(variable * j) + d] = static_cast<uint32_t>(node);
                ++counts[node + 1];
            }
        }
        row_offsets_.assign(sketch_length_ + 1, 0);
        for (size_t row = 0; row < sketch_length_; ++row) {
            row_offsets_[row + 1] = row_offsets_[row] + counts[row + 1];
        }
        column_indices_.resize(total_sockets);
        weights_.resize(total_sockets);
        std::vector<size_t> cursor(row_offsets_.begin(), row_offsets_.end() - 1);
        for (size_t variable = 0; variable < message_length_; ++variable) {
            for (size_t d = 0; d < j; ++d) {
                const size_t node = node_of_socket[(variable * j) + d];
                const size_t slot = cursor[node]++;
                column_indices_[slot] = static_cast<uint32_t>(variable);
                weights_[slot] = prng.next_nonzero_field();
            }
        }
    }

    BoltParams params_;
    size_t message_length_;
    size_t sketch_length_ = 0;
    size_t sketch_codeword_length_ = 0;
    size_t codeword_length_ = 0;
    std::vector<size_t> row_offsets_;
    std::vector<uint32_t> column_indices_;
    std::vector<fr> weights_;
    mutable whir::RSDomains domains_;
};

} // namespace bb::bolt
