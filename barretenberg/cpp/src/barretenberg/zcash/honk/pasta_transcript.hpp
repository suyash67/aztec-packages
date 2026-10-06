#pragma once

#include "barretenberg/common/assert.hpp"
#include "barretenberg/common/throw_or_abort.hpp"
#include "barretenberg/common/type_traits.hpp"
#include "barretenberg/ecc/curves/pasta/pasta.hpp"
#include "barretenberg/numeric/uint256/uint256.hpp"
#include "barretenberg/numeric/uintx/uintx.hpp"
#include "barretenberg/polynomials/univariate.hpp"
#include "barretenberg/transcript/transcript.hpp"
#include "barretenberg/zcash/primitives/blake2b.hpp"

#include <array>
#include <vector>

/**
 * @file pasta_transcript.hpp
 * @brief Fiat-Shamir transcript for proofs over the Pasta cycle, hashing with BLAKE2b as halo2's transcript does.
 */
namespace bb::zcash {

/**
 * @brief Transcript codec for Pasta proofs: every transcript element is one or more 256-bit words.
 * @details Field elements of F_p and F_q are one word each; Pallas and Vesta points are two words (x, y), the point at
 * infinity being (0, 0). Challenges are full-width elements of F_p (which are also canonical elements of F_q, as
 * p < q), or 127-bit halves of one.
 */
class PastaCodec {
  public:
    using DataType = uint256_t;
    using fp = pasta::fp;
    using fq = pasta::fq;
    using pallas_point = pallas::g1::affine_element;
    using vesta_point = vesta::g1::affine_element;

    template <typename T> static constexpr size_t calc_num_fields()
    {
        if constexpr (IsAnyOf<T, uint32_t, uint64_t, uint256_t, bool, fp, fq>) {
            return 1;
        } else if constexpr (IsAnyOf<T, pallas_point, vesta_point>) {
            return 2;
        } else {
            return calc_num_fields<typename T::value_type>() * (std::tuple_size<T>::value);
        }
    }

    template <typename T> static T deserialize_from_fields(std::span<const uint256_t> vec)
    {
        BB_ASSERT_EQ(vec.size(), calc_num_fields<T>());
        if constexpr (IsAnyOf<T, bool>) {
            return static_cast<bool>(vec[0]);
        } else if constexpr (IsAnyOf<T, fp, fq>) {
            BB_ASSERT_LT(vec[0], uint256_t(T::modulus), "Non-canonical field element");
            return T(vec[0]);
        } else if constexpr (IsAnyOf<T, uint32_t, uint64_t, uint256_t>) {
            return static_cast<T>(vec[0]);
        } else if constexpr (IsAnyOf<T, pallas_point, vesta_point>) {
            using BaseField = typename T::Fq;
            T val;
            val.x = deserialize_from_fields<BaseField>(vec.subspan(0, 1));
            val.y = deserialize_from_fields<BaseField>(vec.subspan(1, 1));
            if (val.x.is_zero() && val.y.is_zero()) {
                val.self_set_infinity();
            }
            if (!val.on_curve()) {
                throw_or_abort("Deserialized point is not on the curve");
            }
            return val;
        } else {
            T val;
            constexpr size_t SZ = calc_num_fields<typename T::value_type>();
            size_t i = 0;
            for (auto& x : val) {
                x = deserialize_from_fields<typename T::value_type>(vec.subspan(SZ * i, SZ));
                ++i;
            }
            return val;
        }
    }

    template <typename T> static std::vector<uint256_t> serialize_to_fields(const T& val)
    {
        if constexpr (IsAnyOf<T, bool, uint32_t, uint64_t, uint256_t>) {
            return { uint256_t(val) };
        } else if constexpr (IsAnyOf<T, fp, fq>) {
            return { uint256_t(val) };
        } else if constexpr (IsAnyOf<T, pallas_point, vesta_point>) {
            if (val.is_point_at_infinity()) {
                return { uint256_t(0), uint256_t(0) };
            }
            return { uint256_t(val.x), uint256_t(val.y) };
        } else {
            std::vector<uint256_t> out;
            for (auto& e : val) {
                auto tmp = serialize_to_fields(e);
                out.insert(out.end(), tmp.begin(), tmp.end());
            }
            return out;
        }
    }

    static std::array<uint256_t, 2> split_challenge(const uint256_t& challenge)
    {
        static constexpr size_t LO_BITS = 127;
        return { challenge.slice(0, LO_BITS), challenge.slice(LO_BITS, 2 * LO_BITS) };
    }

    template <typename T> static T convert_short_challenge(const uint256_t& challenge) { return T(challenge); }
    template <typename T> static T convert_full_challenge(const uint256_t& challenge) { return T(challenge); }
};

/**
 * @brief BLAKE2b-512 with halo2's transcript personalization; the 64-byte digest is reduced modulo p.
 */
struct Blake2bTranscriptHash {
    static uint256_t hash(const std::vector<uint256_t>& input)
    {
        static constexpr std::array<uint8_t, 16> PERSONAL = { 'H', 'a', 'l', 'o', '2', '-', 'T', 'r',
                                                              'a', 'n', 's', 'c', 'r', 'i', 'p', 't' };
        Blake2b hasher(64, PERSONAL);
        for (const auto& word : input) {
            std::array<uint8_t, 32> bytes{};
            for (size_t i = 0; i < 4; ++i) {
                for (size_t j = 0; j < 8; ++j) {
                    bytes[(8 * i) + j] = static_cast<uint8_t>(word.data[i] >> (8 * j));
                }
            }
            hasher.update(bytes);
        }
        const auto digest = hasher.finalize();
        // Interpret the digest as a little-endian 512-bit integer and reduce it modulo p.
        uint512_t wide(0);
        for (size_t i = 0; i < 64; ++i) {
            wide += uint512_t(uint256_t(digest[i])) << (8 * i);
        }
        return (wide % uint512_t(uint256_t(pasta::fp::modulus))).lo;
    }
};

using PastaTranscript = BaseTranscript<PastaCodec, Blake2bTranscriptHash>;

} // namespace bb::zcash
