#pragma once

#include "barretenberg/common/assert.hpp"
#include "barretenberg/numeric/uint256/uint256.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace bb::zcash {

/**
 * @brief Poseidon with width 3, rate 2, x^5 S-box, 8 full and 56 partial rounds ("P128Pow5T3"), the nullifier PRF of
 * Orchard.
 * @details Round constants and the MDS matrix are derived with the Grain LFSR exactly as halo2_poseidon's
 * `generate_constants` (the reference construction of the Poseidon paper, with field size `FF::modulus` bit length and
 * secure-MDS index 0). For the Pallas base field this reproduces halo2_poseidon's hardcoded `fp.rs` constants, and the
 * same derivation gives a P128Pow5T3 instance for any other prime field (e.g. BN254's scalar field).
 */
template <typename FF> class PoseidonP128Pow5T3 {
  public:
    static constexpr size_t WIDTH = 3;
    static constexpr size_t RATE = 2;
    static constexpr size_t FULL_ROUNDS = 8;
    static constexpr size_t PARTIAL_ROUNDS = 56;
    static constexpr size_t NUM_ROUNDS = FULL_ROUNDS + PARTIAL_ROUNDS;

    using State = std::array<FF, WIDTH>;
    using Matrix = std::array<std::array<FF, WIDTH>, WIDTH>;

    struct Params {
        std::vector<State> round_constants; // NUM_ROUNDS rows
        Matrix mds;
        Matrix mds_inv;
    };

    static const Params& params()
    {
        static const Params p = generate();
        return p;
    }

    static FF sbox(const FF& x)
    {
        const FF x2 = x.sqr();
        return x2.sqr() * x;
    }

    static void permute(State& state)
    {
        const auto& p = params();
        auto apply_mds = [&](State& s) {
            State out{};
            for (size_t i = 0; i < WIDTH; ++i) {
                FF acc = 0;
                for (size_t j = 0; j < WIDTH; ++j) {
                    acc += p.mds[i][j] * s[j];
                }
                out[i] = acc;
            }
            s = out;
        };
        size_t r = 0;
        for (size_t i = 0; i < FULL_ROUNDS / 2; ++i, ++r) {
            for (size_t j = 0; j < WIDTH; ++j) {
                state[j] = sbox(state[j] + p.round_constants[r][j]);
            }
            apply_mds(state);
        }
        for (size_t i = 0; i < PARTIAL_ROUNDS; ++i, ++r) {
            for (size_t j = 0; j < WIDTH; ++j) {
                state[j] += p.round_constants[r][j];
            }
            state[0] = sbox(state[0]);
            apply_mds(state);
        }
        for (size_t i = 0; i < FULL_ROUNDS / 2; ++i, ++r) {
            for (size_t j = 0; j < WIDTH; ++j) {
                state[j] = sbox(state[j] + p.round_constants[r][j]);
            }
            apply_mds(state);
        }
    }

    // Capacity element of the ConstantLength<L> domain: L * 2^64.
    static FF constant_length_capacity(size_t length) { return FF(uint256_t(length) << 64); }

    /**
     * @brief Poseidon::Hash<ConstantLength<2>> as used for Orchard's PRF^nf: one absorption of (a, b), one
     * permutation, output state[0].
     */
    static FF hash(const FF& a, const FF& b)
    {
        State s{ a, b, constant_length_capacity(2) };
        permute(s);
        return s[0];
    }

  private:
    // Grain LFSR in self-shrinking mode, as specified by the Poseidon reference implementation.
    class Grain {
      public:
        Grain(size_t field_bits, size_t width, size_t r_f, size_t r_p)
        {
            state_.fill(true);
            auto set_bits = [&](size_t offset, size_t len, uint64_t value) {
                for (size_t i = 0; i < len; ++i) {
                    state_[offset + len - 1 - i] = ((value >> i) & 1) != 0;
                }
            };
            set_bits(0, 2, 1); // prime field
            set_bits(2, 4, 0); // x^alpha S-box
            set_bits(6, 12, field_bits);
            set_bits(18, 12, width);
            set_bits(30, 10, r_f);
            set_bits(40, 10, r_p);
            next_bit_ = STATE;
            for (size_t i = 0; i < 20; ++i) {
                load_next_8_bits();
                next_bit_ = STATE;
            }
        }

        bool next()
        {
            while (!get_next_bit()) {
                get_next_bit();
            }
            return get_next_bit();
        }

      private:
        static constexpr size_t STATE = 80;
        std::array<bool, STATE> state_{};
        size_t next_bit_ = STATE;

        void load_next_8_bits()
        {
            std::array<bool, 8> new_bits{};
            for (size_t i = 0; i < 8; ++i) {
                new_bits[i] =
                    state_[i + 62] ^ state_[i + 51] ^ state_[i + 38] ^ state_[i + 23] ^ state_[i + 13] ^ state_[i];
            }
            std::rotate(state_.begin(), state_.begin() + 8, state_.end());
            next_bit_ -= 8;
            for (size_t i = 0; i < 8; ++i) {
                state_[next_bit_ + i] = new_bits[i];
            }
        }

        bool get_next_bit()
        {
            if (next_bit_ == STATE) {
                load_next_8_bits();
            }
            return state_[next_bit_++];
        }
    };

    static size_t num_bits() { return static_cast<size_t>(uint256_t(FF::modulus).get_msb()) + 1; }

    // Big-endian bit string of length NUM_BITS read from the LFSR.
    static uint256_t next_bits(Grain& grain)
    {
        const size_t nb = num_bits();
        uint256_t v = 0;
        for (size_t i = 0; i < nb; ++i) {
            v = (v << 1) + uint256_t(grain.next() ? 1 : 0);
        }
        return v;
    }

    static FF next_field_element(Grain& grain)
    {
        while (true) {
            const uint256_t v = next_bits(grain);
            if (v < FF::modulus) {
                return FF(v);
            }
        }
    }

    static FF next_field_element_without_rejection(Grain& grain)
    {
        const uint256_t v = next_bits(grain);
        return FF((uint512_t(v) % uint512_t(FF::modulus)).lo);
    }

    static Params generate()
    {
        Grain grain(num_bits(), WIDTH, FULL_ROUNDS, PARTIAL_ROUNDS);
        Params p;
        p.round_constants.resize(NUM_ROUNDS);
        for (auto& row : p.round_constants) {
            for (auto& c : row) {
                c = next_field_element(grain);
            }
        }
        // Cauchy MDS from the first set of 2 * WIDTH distinct elements (secure-MDS index 0).
        std::array<FF, 2 * WIDTH> vals{};
        while (true) {
            for (auto& v : vals) {
                v = next_field_element_without_rejection(grain);
            }
            bool distinct = true;
            for (size_t i = 0; i < vals.size(); ++i) {
                for (size_t j = i + 1; j < vals.size(); ++j) {
                    distinct &= (vals[i] != vals[j]);
                }
            }
            if (distinct) {
                break;
            }
        }
        std::array<FF, WIDTH> xs{};
        std::array<FF, WIDTH> ys{};
        for (size_t i = 0; i < WIDTH; ++i) {
            xs[i] = vals[i];
            ys[i] = vals[WIDTH + i];
        }
        for (size_t i = 0; i < WIDTH; ++i) {
            for (size_t j = 0; j < WIDTH; ++j) {
                const FF sum = xs[i] + ys[j];
                BB_ASSERT(!sum.is_zero());
                p.mds[i][j] = sum.invert();
            }
        }
        // Closed-form inverse of the Cauchy matrix (halo2_poseidon mds.rs).
        auto lagrange = [](const std::array<FF, WIDTH>& pts, size_t j, const FF& x) {
            FF acc = 1;
            for (size_t m = 0; m < WIDTH; ++m) {
                if (m != j) {
                    acc *= (x - pts[m]) * (pts[j] - pts[m]).invert();
                }
            }
            return acc;
        };
        std::array<FF, WIDTH> neg_ys{};
        for (size_t i = 0; i < WIDTH; ++i) {
            neg_ys[i] = -ys[i];
        }
        for (size_t i = 0; i < WIDTH; ++i) {
            for (size_t j = 0; j < WIDTH; ++j) {
                p.mds_inv[i][j] = (xs[j] - neg_ys[i]) * lagrange(xs, j, neg_ys[i]) * lagrange(neg_ys, i, xs[j]);
            }
        }
        return p;
    }
};

} // namespace bb::zcash
