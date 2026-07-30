#include "barretenberg/commitment_schemes/whir/weights.hpp"

#include "barretenberg/commitment_schemes/whir/rs_code.hpp"
#include "barretenberg/polynomials/polynomial.hpp"
#include "barretenberg/polynomials/polynomial_arithmetic.hpp"

#include <gtest/gtest.h>

namespace bb::whir {

namespace {
std::vector<fr> random_array(size_t size)
{
    std::vector<fr> array(size);
    for (fr& value : array) {
        value = fr::random_element();
    }
    return array;
}

std::vector<fr> random_point(size_t size)
{
    return random_array(size);
}
} // namespace

// ∑_b eq(u, b)·a_b = MLE of a at u.
TEST(WhirWeights, EqWeightedSumIsMleEvaluation)
{
    constexpr size_t m = 5;
    const std::vector<fr> array = random_array(size_t(1) << m);
    const std::vector<fr> u = random_point(m);

    const WeightTerm term = WeightTerm::eq_weight(u, fr::one());
    const Polynomial<fr> reference{ std::span<const fr>(array) };
    EXPECT_EQ(term.weighted_sum(array), reference.evaluate_mle(u));
}

// ∑_b pow_y(b)·a_b = A(y), the univariate evaluation: README.md eq. (3.2).
TEST(WhirWeights, PowWeightedSumIsUnivariateEvaluation)
{
    constexpr size_t m = 5;
    const std::vector<fr> array = random_array(size_t(1) << m);
    const fr y = fr::random_element();

    const WeightTerm term = WeightTerm::pow_weight(y, m, fr::one());
    EXPECT_EQ(term.weighted_sum(array), polynomial_arithmetic::evaluate(array.data(), y, array.size()));
}

// accumulate_table produces exactly the per-point evaluations, scaled by coeff.
TEST(WhirWeights, TableMatchesEvaluate)
{
    constexpr size_t m = 4;
    const std::vector<fr> u = random_point(m);
    const fr coeff = fr::random_element();
    const WeightTerm term = WeightTerm::eq_weight(u, coeff);

    std::vector<fr> table(size_t(1) << m, fr::zero());
    term.accumulate_table(table);
    for (size_t b = 0; b < table.size(); ++b) {
        EXPECT_EQ(table[b], term.evaluate(b));
    }
}

// Binding a variable in the term equals folding the term's table: the prover-side table fold and the
// verifier-side symbolic bind stay consistent.
TEST(WhirWeights, BindMatchesTableFold)
{
    constexpr size_t m = 4;
    const fr y = fr::random_element();
    WeightTerm term = WeightTerm::pow_weight(y, m, fr::random_element());

    std::vector<fr> table(size_t(1) << m, fr::zero());
    term.accumulate_table(table);

    const fr alpha = fr::random_element();
    fold_array_in_place(table, alpha);
    term.bind_variable(alpha);

    std::vector<fr> bound_table(size_t(1) << (m - 1), fr::zero());
    term.accumulate_table(bound_table);
    EXPECT_EQ(bound_table, table);
}

// Weighted sums are additive across terms, matching the combined-table prover representation.
TEST(WhirWeights, CombinedTableMatchesTermSum)
{
    constexpr size_t m = 3;
    const std::vector<fr> array = random_array(size_t(1) << m);
    const WeightTerm eq_term = WeightTerm::eq_weight(random_point(m), fr::random_element());
    const WeightTerm pow_term = WeightTerm::pow_weight(fr::random_element(), m, fr::random_element());

    std::vector<fr> table(size_t(1) << m, fr::zero());
    eq_term.accumulate_table(table);
    pow_term.accumulate_table(table);

    fr combined = fr::zero();
    for (size_t b = 0; b < array.size(); ++b) {
        combined += table[b] * array[b];
    }
    EXPECT_EQ(combined, eq_term.weighted_sum(array) + pow_term.weighted_sum(array));
}

} // namespace bb::whir
