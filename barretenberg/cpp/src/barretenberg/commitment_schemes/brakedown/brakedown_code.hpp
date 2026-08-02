#pragma once

#include "barretenberg/common/assert.hpp"
#include "barretenberg/common/thread.hpp"
#include "barretenberg/crypto/blake3s/blake3s.hpp"
#include "barretenberg/ecc/curves/bn254/fr.hpp"
#include "barretenberg/numeric/uint256/uint256.hpp"

#include <cmath>
#include <cstring>
#include <span>
#include <vector>

namespace bb::brakedown {

/** @brief Binary entropy in bits. */
inline double binary_entropy(double p)
{
    if (p <= 0.0 || p >= 1.0) {
        return 0.0;
    }
    return -p * std::log2(p) - (1.0 - p) * std::log2(1.0 - p);
}

/**
 * @brief Parameters of the Spielman-style code of Brakedown (eprint 2021/1043, Algorithm 1).
 *
 * @details The code is systematic, linear, of rate 1/r and relative distance beta/r. `alpha` sets
 * the recursion's contraction, `beta` the guaranteed relative weight of a nonzero codeword. The
 * admissible region is `0 < beta < alpha/1.28` and `r > (1 + 2 beta)/(1 - alpha)`.
 *
 * The named presets are the rows of the paper's Figure 2: each was chosen so that the code fails to
 * meet its distance with probability at most 2^-100 for messages up to 2^30 over a field of at
 * least 2^127 elements. BN254 Fr is larger than that, so the bound carries over.
 */
struct BrakedownParams {
    double alpha;
    double beta;
    double r;
    size_t base_length = 32; // recursion bottoms out at or below this message length

    double distance() const { return beta / r; }
    double rate() const { return 1.0 / r; }

    void validate() const
    {
        BB_ASSERT(beta > 0.0 && beta < alpha / 1.28, "Brakedown requires 0 < beta < alpha/1.28");
        BB_ASSERT(r > (1.0 + 2.0 * beta) / (1.0 - alpha), "Brakedown requires r > (1+2beta)/(1-alpha)");
    }

    /** @brief Figure 2 rows, indexed by the relative distance they achieve. */
    static BrakedownParams distance_2pct() { return { 0.1195, 0.0284, 1.42 }; }
    static BrakedownParams distance_3pct() { return { 0.138, 0.0444, 1.47 }; }
    static BrakedownParams distance_4pct() { return { 0.178, 0.061, 1.521 }; }
    static BrakedownParams distance_5pct() { return { 0.2, 0.082, 1.64 }; }
    static BrakedownParams distance_6pct() { return { 0.211, 0.097, 1.616 }; }
    static BrakedownParams distance_7pct() { return { 0.238, 0.1205, 1.72 }; }
};

/**
 * @brief Row sparsity of the first matrix, paper Equation (7).
 * @details c_n = min( max(1.28*beta*n, beta*n+4),
 *                     (110/n + H(beta) + alpha*H(1.28*beta/alpha)) / (beta*log2(alpha/(1.28*beta))) )
 * The 110 is the paper's failure-probability slack (2^-100 plus a margin).
 */
inline size_t row_sparsity_c(size_t n, const BrakedownParams& p)
{
    const double dn = static_cast<double>(n);
    const double cap = std::max(1.28 * p.beta * dn, p.beta * dn + 4.0);
    const double numerator = 110.0 / dn + binary_entropy(p.beta) + p.alpha * binary_entropy(1.28 * p.beta / p.alpha);
    const double denominator = p.beta * std::log2(p.alpha / (1.28 * p.beta));
    const double bound = numerator / denominator;
    return static_cast<size_t>(std::max(3.0, std::ceil(std::min(cap, bound))));
}

/**
 * @brief The real-valued bound underlying Equation (8), before the ceiling.
 * @details Exposed because Figure 2's published alpha/beta/r are rounded to three or four
 * significant figures, so the ceiling is not a stable thing to compare against; the bound is.
 */
inline double row_sparsity_d_bound(size_t n, const BrakedownParams& p, double log2_field_size)
{
    const double dn = static_cast<double>(n);
    const double mu = p.r - 1.0 - (p.r * p.alpha);
    const double nu = p.beta + (p.alpha * p.beta) + 0.03;
    const double r_alpha = p.r * p.alpha;
    const double two_beta = (2.0 * p.beta) + 0.03;

    const double first = ((r_alpha * binary_entropy(p.beta / p.r)) + (mu * binary_entropy(nu / mu)) + (110.0 / dn)) /
                         (p.alpha * p.beta * std::log2(mu / nu));
    const double second =
        ((r_alpha * binary_entropy(p.beta / r_alpha)) + (mu * binary_entropy(two_beta / mu)) + (110.0 / dn)) /
        (p.beta * std::log2(mu / two_beta));
    const double third =
        (two_beta * ((1.0 / (r_alpha - p.beta)) + (1.0 / (p.alpha * p.beta)) + (1.0 / (mu - two_beta)))) + 1.0;

    const double cap = ((2.0 * p.beta) + (((p.r - 1.0) + (110.0 / dn)) / log2_field_size)) * dn;
    return std::min(cap, std::max({ first, second, third }));
}

/**
 * @brief Row sparsity of the second matrix, paper Equation (8).
 * @details d_n = min( (2*beta + ((r-1) + 110/n)/log2(q)) * n, D ), with D the maximum of the three
 * expressions the paper gives; mu = r-1-r*alpha and nu = beta + alpha*beta + 0.03.
 */
inline size_t row_sparsity_d(size_t n, const BrakedownParams& p, double log2_field_size)
{
    return static_cast<size_t>(std::max(3.0, std::ceil(row_sparsity_d_bound(n, p, log2_field_size))));
}

namespace detail {

/** @brief Deterministic Blake3-based stream, so both parties derive the identical matrices. */
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
        std::memcpy(&value, block_.data() + 8 * (4 - cursor_), 8);
        --cursor_;
        return value;
    }

    /** @brief Uniform in [0, bound); the modulo bias is below 2^-40 for the sizes used here. */
    size_t next_index(size_t bound) { return static_cast<size_t>(next() % bound); }

    fr next_nonzero_field()
    {
        // Two 64-bit words are enough entropy for a nonzero sparse-matrix entry, and keeping the
        // draw small keeps preprocessing well under the encode cost it serves. The words are drawn
        // into locals first: argument evaluation order is unspecified, which would otherwise make
        // the generated matrices compiler-dependent and break prover/verifier agreement.
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

/**
 * @brief A sparse matrix with a fixed number of nonzeros per row, stored compressed by *column*.
 *
 * @details The construction specifies the matrices row-wise, but the product this code needs is
 * `vector * matrix`, which row-major storage turns into a scatter (`out[index] += ...`). Scatters
 * race, so they cannot be parallelised, and their random writes miss cache. Transposing to
 * column-major at build time makes each output entry an independent gather, which parallelises
 * cleanly and reads the input vector in a cache-friendly order. Preprocessing pays for it once.
 */
struct SparseMatrix {
    size_t num_rows = 0;
    size_t num_cols = 0;
    size_t per_row = 0;
    std::vector<size_t> col_offsets; // num_cols + 1 entries
    std::vector<uint32_t> row_indices;
    std::vector<fr> values;

    /** @brief out += row_vector * this. */
    void accumulate_product(std::span<const fr> row_vector, std::span<fr> out) const
    {
        BB_ASSERT_EQ(row_vector.size(), num_rows);
        BB_ASSERT_EQ(out.size(), num_cols);
        parallel_for_range(num_cols, [&](size_t start, size_t end) {
            for (size_t col = start; col < end; ++col) {
                fr acc = fr::zero();
                for (size_t k = col_offsets[col]; k < col_offsets[col + 1]; ++k) {
                    acc += row_vector[row_indices[k]] * values[k];
                }
                out[col] += acc;
            }
        });
    }
};

inline SparseMatrix sample_matrix(size_t num_rows, size_t num_cols, size_t per_row, SeededPrng& prng)
{
    SparseMatrix matrix;
    matrix.num_rows = num_rows;
    matrix.num_cols = num_cols;
    matrix.per_row = std::min(per_row, num_cols);

    // Sample row-wise as the construction specifies, then transpose into column-major.
    const size_t nnz = num_rows * matrix.per_row;
    std::vector<uint32_t> cols(nnz);
    std::vector<fr> values(nnz);
    std::vector<size_t> chosen;
    for (size_t row = 0; row < num_rows; ++row) {
        chosen.clear();
        while (chosen.size() < matrix.per_row) {
            const size_t candidate = prng.next_index(num_cols);
            if (std::find(chosen.begin(), chosen.end(), candidate) == chosen.end()) {
                chosen.push_back(candidate);
            }
        }
        for (size_t k = 0; k < matrix.per_row; ++k) {
            cols[(row * matrix.per_row) + k] = static_cast<uint32_t>(chosen[k]);
            values[(row * matrix.per_row) + k] = prng.next_nonzero_field();
        }
    }

    matrix.col_offsets.assign(num_cols + 1, 0);
    for (const uint32_t col : cols) {
        ++matrix.col_offsets[col + 1];
    }
    for (size_t col = 0; col < num_cols; ++col) {
        matrix.col_offsets[col + 1] += matrix.col_offsets[col];
    }
    matrix.row_indices.resize(nnz);
    matrix.values.resize(nnz);
    std::vector<size_t> cursor(matrix.col_offsets.begin(), matrix.col_offsets.end() - 1);
    for (size_t row = 0; row < num_rows; ++row) {
        for (size_t k = 0; k < matrix.per_row; ++k) {
            const size_t index = (row * matrix.per_row) + k;
            const size_t slot = cursor[cols[index]]++;
            matrix.row_indices[slot] = static_cast<uint32_t>(row);
            matrix.values[slot] = values[index];
        }
    }
    return matrix;
}

} // namespace detail

/**
 * @brief Brakedown's linear-time systematic code (eprint 2021/1043, Algorithm 1).
 *
 * @details Encoding a message x of length n:
 *
 *     y = x · A          A is n x ceil(alpha n),        c_n nonzeros per row
 *     z = Enc(y)         recursive, on the shorter message
 *     v = z · B          B is |z| x (rn - n - |z|),     d_n nonzeros per row
 *     Enc(x) = (x, z, v)
 *
 * Every level costs `n·c_n + |z|·d_n` multiplications and the message length contracts by `alpha`,
 * so the total is `O(n)` — no FFT, no smooth subgroup, no field structure beyond size. That is the
 * whole point: it is the linear-time, field-agnostic alternative to the Reed-Solomon encoder that
 * `whir/rs_code.hpp` provides, and it is the shared prerequisite for the Lightning and Bolt
 * constructions (see ../PCS_CANDIDATES.md §7).
 *
 * The matrices are sampled once per length in preprocessing from a seed, so prover and verifier
 * derive bit-identical codes from the same `(params, seed, message_length)`.
 */
class BrakedownCode {
  public:
    BrakedownCode(size_t message_length, const BrakedownParams& params, uint64_t seed = 0x62726b64u)
        : params_(params)
        , message_length_(message_length)
    {
        params.validate();
        BB_ASSERT_GT(message_length, size_t(0));
        build(message_length, seed);
    }

    size_t message_length() const { return message_length_; }
    size_t codeword_length() const { return codeword_length_; }
    const BrakedownParams& params() const { return params_; }

    /** @brief Number of field multiplications one encoding costs, summed over the recursion. */
    size_t multiplications() const
    {
        size_t total = 0;
        for (const Level& level : levels_) {
            total += level.a.num_rows * level.a.per_row + level.b.num_rows * level.b.per_row;
        }
        for (const auto& row : base_parity_) {
            total += row.size();
        }
        return total;
    }

    std::vector<fr> encode(std::span<const fr> message) const
    {
        BB_ASSERT_EQ(message.size(), message_length_, "message length does not match the code");
        std::vector<fr> codeword(codeword_length_, fr::zero());
        encode_into(message, codeword, 0);
        return codeword;
    }

  private:
    struct Level {
        size_t message_length;
        size_t codeword_length;
        size_t inner_message_length;  // |y|
        size_t inner_codeword_length; // |z|
        detail::SparseMatrix a;
        detail::SparseMatrix b;
    };

    static size_t scaled(double factor, size_t n)
    {
        return static_cast<size_t>(std::ceil(factor * static_cast<double>(n)));
    }

    void build(size_t n, uint64_t seed)
    {
        // Walk the recursion down, sizing every level, then fill the base case.
        std::vector<size_t> lengths;
        size_t current = n;
        while (current > params_.base_length) {
            lengths.push_back(current);
            current = scaled(params_.alpha, current);
        }
        base_length_ = current;
        base_codeword_length_ = std::max(scaled(params_.r, current), current + 1);
        build_base(seed);

        // Sizes resolve bottom-up because |z| is the child's codeword length.
        size_t child_codeword = base_codeword_length_;
        const double log2_field = 254.0; // BN254 Fr; the paper's tables assume at least 2^127
        for (size_t i = lengths.size(); i-- > 0;) {
            const size_t length = lengths[i];
            Level level;
            level.message_length = length;
            level.inner_message_length = scaled(params_.alpha, length);
            level.inner_codeword_length = child_codeword;
            level.codeword_length = std::max(scaled(params_.r, length), length + child_codeword + 1);
            const size_t tail = level.codeword_length - length - child_codeword;
            BB_ASSERT_GT(tail, size_t(0), "Brakedown parameters leave no room for the second block");

            detail::SeededPrng prng_a(seed, 0x4141410000000000ull + length);
            detail::SeededPrng prng_b(seed, 0x4242420000000000ull + length);
            level.a =
                detail::sample_matrix(length, level.inner_message_length, row_sparsity_c(length, params_), prng_a);
            level.b = detail::sample_matrix(child_codeword, tail, row_sparsity_d(length, params_, log2_field), prng_b);

            child_codeword = level.codeword_length;
            levels_.push_back(std::move(level));
        }
        // levels_ is bottom-up; encoding walks it in that order too.
        codeword_length_ = child_codeword;
    }

    /**
     * @brief Base case: a systematic Reed-Solomon-style code over distinct evaluation points.
     * @details At the base the message is tiny (<= `base_length`), so a dense parity block costs
     * nothing asymptotically and is MDS, whose relative distance (r-1)/r comfortably exceeds the
     * beta/r the recursion needs from it.
     */
    void build_base(uint64_t seed)
    {
        const size_t parity = base_codeword_length_ - base_length_;
        base_parity_.assign(parity, {});
        detail::SeededPrng prng(seed, 0x4343430000000000ull + base_length_);
        std::vector<fr> points(parity);
        for (size_t j = 0; j < parity; ++j) {
            points[j] = fr(static_cast<uint64_t>(j + 2));
        }
        static_cast<void>(prng);
        for (size_t j = 0; j < parity; ++j) {
            base_parity_[j].resize(base_length_);
            fr power = fr::one();
            for (size_t i = 0; i < base_length_; ++i) {
                base_parity_[j][i] = power;
                power *= points[j];
            }
        }
    }

    void encode_base(std::span<const fr> message, std::span<fr> out) const
    {
        std::copy(message.begin(), message.end(), out.begin());
        for (size_t j = 0; j < base_parity_.size(); ++j) {
            fr acc = fr::zero();
            for (size_t i = 0; i < base_length_; ++i) {
                acc += base_parity_[j][i] * message[i];
            }
            out[base_length_ + j] = acc;
        }
    }

    /** @brief Encode `message` (of the level's length) into `out`, `level` counted bottom-up. */
    void encode_into(std::span<const fr> message, std::span<fr> out, size_t level) const
    {
        if (level == levels_.size()) {
            std::vector<fr> padded(base_length_, fr::zero());
            std::copy(message.begin(), message.end(), padded.begin());
            encode_base(padded, out.subspan(0, base_codeword_length_));
            return;
        }
        // levels_ is bottom-up, so the outermost level is the last entry.
        const Level& current = levels_[levels_.size() - 1 - level];
        BB_ASSERT_EQ(message.size(), current.message_length);

        std::copy(message.begin(), message.end(), out.begin());
        std::vector<fr> y(current.inner_message_length, fr::zero());
        current.a.accumulate_product(message, y);

        const std::span<fr> z = out.subspan(current.message_length, current.inner_codeword_length);
        encode_into(y, z, level + 1);

        const std::span<fr> v = out.subspan(current.message_length + current.inner_codeword_length);
        std::fill(v.begin(), v.end(), fr::zero());
        current.b.accumulate_product(std::span<const fr>(z), v);
    }

    BrakedownParams params_;
    size_t message_length_;
    size_t codeword_length_ = 0;
    size_t base_length_ = 0;
    size_t base_codeword_length_ = 0;
    std::vector<Level> levels_; // bottom-up
    std::vector<std::vector<fr>> base_parity_;
};

} // namespace bb::brakedown
