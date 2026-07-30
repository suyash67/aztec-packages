#pragma once

#include "barretenberg/common/assert.hpp"
#include "barretenberg/numeric/random/engine.hpp"

#include <cstdint>
#include <vector>

namespace bb::small_field {

/**
 * @brief The Mersenne-31 prime field, p = 2^31 - 1: the small-field design point used by Circle
 * STARKs / stwo / Plonky3. One value is 4 bytes (an eighth of a BN254 Fr element); reduction is
 * shift-and-add, no Montgomery form. See SMALL_FIELDS.md for how (and whether) this field can back
 * a Honk-style stack.
 */
class m31 {
  public:
    static constexpr uint32_t modulus = 0x7fffffffU;

    constexpr m31() = default;
    constexpr explicit m31(uint64_t value)
        : value_(reduce64(value))
    {}

    static constexpr m31 zero() { return m31(0); }
    static constexpr m31 one() { return m31(1); }
    constexpr uint32_t value() const { return value_; }

    constexpr m31 operator+(const m31& other) const
    {
        uint32_t sum = value_ + other.value_; // <= 2p - 2 < 2^32
        sum = (sum & modulus) + (sum >> 31);
        return from_reduced(sum >= modulus ? sum - modulus : sum);
    }
    constexpr m31 operator-(const m31& other) const
    {
        return *this + from_reduced(other.value_ == 0 ? 0 : modulus - other.value_);
    }
    constexpr m31 operator-() const { return from_reduced(value_ == 0 ? 0 : modulus - value_); }
    constexpr m31 operator*(const m31& other) const
    {
        const uint64_t product = static_cast<uint64_t>(value_) * other.value_; // < 2^62
        uint32_t folded = reduce64(product);
        return from_reduced(folded);
    }
    constexpr m31& operator+=(const m31& other) { return *this = *this + other; }
    constexpr m31& operator-=(const m31& other) { return *this = *this - other; }
    constexpr m31& operator*=(const m31& other) { return *this = *this * other; }
    constexpr bool operator==(const m31& other) const = default;

    constexpr m31 sqr() const { return *this * *this; }

    constexpr m31 pow(uint64_t exponent) const
    {
        m31 result = one();
        m31 base = *this;
        while (exponent != 0) {
            if ((exponent & 1) != 0) {
                result *= base;
            }
            base = base.sqr();
            exponent >>= 1;
        }
        return result;
    }

    /** @brief Multiplicative inverse; maps zero to zero (Fermat exponentiation). */
    constexpr m31 invert() const { return pow(modulus - 2); }

    static m31 random_element()
    {
        auto& engine = numeric::get_randomness();
        return m31(static_cast<uint64_t>(engine.get_random_uint32()));
    }

  private:
    static constexpr uint32_t reduce64(uint64_t value)
    {
        // Fold twice: 2^31 = 1 (mod p), then a final conditional subtraction.
        uint64_t folded = (value & modulus) + (value >> 31); // < 2^33
        folded = (folded & modulus) + (folded >> 31);        // <= p + 3
        return static_cast<uint32_t>(folded >= modulus ? folded - modulus : folded);
    }
    static constexpr m31 from_reduced(uint32_t value)
    {
        m31 result;
        result.value_ = value;
        return result;
    }

    uint32_t value_ = 0;
};

/** @brief The complex extension CM31 = F_p[i]/(i^2 + 1); p = 3 (mod 4) makes -1 a non-residue. */
class cm31 {
  public:
    constexpr cm31() = default;
    constexpr cm31(m31 real, m31 imag)
        : real_(real)
        , imag_(imag)
    {}
    constexpr explicit cm31(uint64_t value)
        : real_(m31(value))
    {}

    static constexpr cm31 zero() { return {}; }
    static constexpr cm31 one() { return cm31(m31::one(), m31::zero()); }
    constexpr m31 real() const { return real_; }
    constexpr m31 imag() const { return imag_; }

    constexpr cm31 operator+(const cm31& other) const { return { real_ + other.real_, imag_ + other.imag_ }; }
    constexpr cm31 operator-(const cm31& other) const { return { real_ - other.real_, imag_ - other.imag_ }; }
    constexpr cm31 operator*(const cm31& other) const
    {
        // (a + bi)(c + di) = (ac - bd) + (ad + bc)i
        return { real_ * other.real_ - imag_ * other.imag_, real_ * other.imag_ + imag_ * other.real_ };
    }
    constexpr cm31& operator+=(const cm31& other) { return *this = *this + other; }
    constexpr cm31& operator*=(const cm31& other) { return *this = *this * other; }
    constexpr bool operator==(const cm31& other) const = default;

    constexpr cm31 conjugate() const { return { real_, -imag_ }; }
    constexpr m31 norm() const { return real_ * real_ + imag_ * imag_; }
    constexpr cm31 invert() const
    {
        // 1/(a+bi) = (a-bi)/(a^2+b^2)
        const m31 norm_inverse = norm().invert();
        return { real_ * norm_inverse, -imag_ * norm_inverse };
    }

    static cm31 random_element() { return { m31::random_element(), m31::random_element() }; }

  private:
    m31 real_;
    m31 imag_;
};

/**
 * @brief The degree-4 extension QM31 = CM31[u]/(u^2 - (2 + i)), the "secure field" of the M31
 * ecosystem (~124 bits): the field Fiat-Shamir challenges must live in when witnesses are M31.
 */
class qm31 {
  public:
    constexpr qm31() = default;
    constexpr qm31(cm31 first, cm31 second)
        : first_(first)
        , second_(second)
    {}
    constexpr explicit qm31(uint64_t value)
        : first_(cm31(value))
    {}

    static constexpr qm31 zero() { return {}; }
    static constexpr qm31 one() { return qm31(cm31::one(), cm31::zero()); }
    constexpr cm31 first() const { return first_; }
    constexpr cm31 second() const { return second_; }

    constexpr qm31 operator+(const qm31& other) const { return { first_ + other.first_, second_ + other.second_ }; }
    constexpr qm31 operator-(const qm31& other) const { return { first_ - other.first_, second_ - other.second_ }; }
    // u^2 in CM31
    static constexpr cm31 non_residue() { return { m31(2), m31(1) }; }

    constexpr qm31 operator*(const qm31& other) const
    {
        // (x1 + y1 u)(x2 + y2 u) = x1 x2 + (2 + i) y1 y2 + (x1 y2 + y1 x2) u
        return { first_ * other.first_ + non_residue() * (second_ * other.second_),
                 first_ * other.second_ + second_ * other.first_ };
    }
    constexpr qm31& operator+=(const qm31& other) { return *this = *this + other; }
    constexpr qm31& operator*=(const qm31& other) { return *this = *this * other; }
    constexpr bool operator==(const qm31& other) const = default;

    constexpr qm31 invert() const
    {
        // 1/(x + yu) = (x - yu)/(x^2 - (2+i) y^2), the denominator landing in CM31.
        const cm31 denominator = first_ * first_ - non_residue() * (second_ * second_);
        const cm31 denominator_inverse = denominator.invert();
        return { first_ * denominator_inverse, (cm31::zero() - second_) * denominator_inverse };
    }

    /** @brief Multiply a QM31 accumulator by an M31 base-field value (the mixed op sumcheck needs). */
    constexpr qm31 scale(const m31& scalar) const
    {
        const cm31 scalar_ext{ scalar, m31::zero() };
        return { first_ * scalar_ext, second_ * scalar_ext };
    }

    static qm31 random_element() { return { cm31::random_element(), cm31::random_element() }; }

  private:
    cm31 first_;
    cm31 second_;
};

} // namespace bb::small_field
