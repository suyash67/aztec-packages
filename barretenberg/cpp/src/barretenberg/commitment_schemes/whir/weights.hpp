#pragma once

#include "barretenberg/common/assert.hpp"
#include "barretenberg/ecc/curves/bn254/fr.hpp"

#include <span>
#include <vector>

namespace bb::whir {

/** @brief eq tensor table over the given challenge slice, LSB variable first: out[j] = eq_j(u). */
inline std::vector<fr> eq_tensor(std::span<const fr> u)
{
    std::vector<fr> table(size_t(1) << u.size());
    table[0] = fr::one();
    for (size_t i = 0; i < u.size(); ++i) {
        const size_t built = size_t(1) << i;
        for (size_t j = 0; j < built; ++j) {
            table[j + built] = table[j] * u[i];
            table[j] -= table[j + built];
        }
    }
    return table;
}

/**
 * @brief A weight polynomial that is a scalar times a product of per-variable affine factors,
 * W(X) = coeff · ∏ⱼ (aⱼ + bⱼ·Xⱼ). Covers both claim types of README.md eq. (3.2):
 *  - eq(u, X): aⱼ = 1-uⱼ, bⱼ = 2uⱼ-1 (an MLE evaluation claim at u);
 *  - pow_y(X): aⱼ = 1, bⱼ = y^{2^j}-1 (a univariate evaluation claim at y).
 * Variable j = 0 is the low (first-folded) variable. `bind_variable` consumes variables in order,
 * mirroring the prover-side stride-2 table fold.
 */
class WeightTerm {
  public:
    static WeightTerm eq_weight(std::span<const fr> u, const fr& coeff)
    {
        WeightTerm term;
        term.coeff_ = coeff;
        term.a_.reserve(u.size());
        term.b_.reserve(u.size());
        for (const fr& u_j : u) {
            term.a_.push_back(fr(1) - u_j);
            term.b_.push_back(u_j + u_j - fr(1));
        }
        return term;
    }

    static WeightTerm pow_weight(const fr& y, size_t num_variables, const fr& coeff)
    {
        WeightTerm term;
        term.coeff_ = coeff;
        term.a_.assign(num_variables, fr(1));
        term.b_.reserve(num_variables);
        fr power = y; // y^{2^j}
        for (size_t j = 0; j < num_variables; ++j) {
            term.b_.push_back(power - fr(1));
            power = power.sqr();
        }
        return term;
    }

    size_t num_remaining_variables() const { return a_.size() - bound_; }

    /** @brief Consume the next variable at value alpha: coeff *= (aⱼ + bⱼ·α). */
    void bind_variable(const fr& alpha)
    {
        BB_ASSERT_LT(bound_, a_.size(), "no variables left to bind");
        coeff_ *= a_[bound_] + b_[bound_] * alpha;
        ++bound_;
    }

    /** @brief W at the hypercube point with bits `cube_index` (bit j = remaining variable j). */
    fr evaluate(size_t cube_index) const
    {
        fr result = coeff_;
        for (size_t j = bound_; j < a_.size(); ++j) {
            result *= ((cube_index >> (j - bound_)) & 1) ? (a_[j] + b_[j]) : a_[j];
        }
        return result;
    }

    /** @brief table[b] += W(b) over the remaining-variable hypercube (tensor expansion). */
    void accumulate_table(std::span<fr> table) const
    {
        const size_t remaining = num_remaining_variables();
        BB_ASSERT_EQ(table.size(), size_t(1) << remaining, "table size mismatch");
        std::vector<fr> tensor(table.size());
        tensor[0] = coeff_;
        for (size_t j = 0; j < remaining; ++j) {
            const size_t built = size_t(1) << j;
            const fr& a = a_[bound_ + j];
            const fr& b = b_[bound_ + j];
            for (size_t t = 0; t < built; ++t) {
                tensor[t + built] = tensor[t] * (a + b);
                tensor[t] *= a;
            }
        }
        for (size_t t = 0; t < table.size(); ++t) {
            table[t] += tensor[t];
        }
    }

    /** @brief ∑_b W(b)·values[b] over the remaining-variable hypercube. */
    fr weighted_sum(std::span<const fr> values) const
    {
        BB_ASSERT_EQ(values.size(), size_t(1) << num_remaining_variables(), "value table size mismatch");
        fr total = fr::zero();
        for (size_t b = 0; b < values.size(); ++b) {
            total += evaluate(b) * values[b];
        }
        return total;
    }

  private:
    fr coeff_ = fr::zero();
    std::vector<fr> a_;
    std::vector<fr> b_;
    size_t bound_ = 0;
};

} // namespace bb::whir
