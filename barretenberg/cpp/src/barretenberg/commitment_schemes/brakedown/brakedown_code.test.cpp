#include "barretenberg/commitment_schemes/brakedown/brakedown_code.hpp"

#include <gtest/gtest.h>

namespace bb::brakedown {

namespace {

std::vector<fr> random_message(size_t size)
{
    std::vector<fr> message(size);
    for (fr& value : message) {
        value = fr::random_element();
    }
    return message;
}

size_t hamming_weight(std::span<const fr> word)
{
    size_t weight = 0;
    for (const fr& value : word) {
        if (!value.is_zero()) {
            ++weight;
        }
    }
    return weight;
}

} // namespace

// The paper's Figure 2 lists c_n and d_n for each parameter row "for all large enough n". Deriving
// them from Equations (7) and (8) and reproducing the published integers is the strongest available
// check that both formulas are transcribed correctly.
//
// c_n reproduces every row exactly. For d_n the check is on the real-valued bound rather than its
// ceiling, because Figure 2 quotes alpha, beta and r to three or four significant figures: at those
// rounded parameters two rows evaluate to 26.0033 and 22.0014 against published 26 and 22, so the
// true optima sit a hair under integers the ceiling rounds up. Asserting that the bound lands
// within rounding slack of the published value pins the formula more tightly than comparing
// ceilings would, and the implementation keeps the conservative ceiling — one extra nonzero per row
// only improves distance.
TEST(BrakedownParamsTest, SparsitiesMatchPaperFigure2)
{
    struct Row {
        BrakedownParams params;
        double distance;
        size_t c_n;
        size_t d_n;
    };
    const std::vector<Row> rows = {
        { BrakedownParams::distance_2pct(), 0.02, 6, 33 }, { BrakedownParams::distance_3pct(), 0.03, 7, 26 },
        { BrakedownParams::distance_4pct(), 0.04, 7, 22 }, { BrakedownParams::distance_5pct(), 0.05, 8, 19 },
        { BrakedownParams::distance_6pct(), 0.06, 9, 21 }, { BrakedownParams::distance_7pct(), 0.07, 10, 23 },
    };
    // "Large enough n": the 110/n slack has to be negligible, which it is well before 2^20.
    const size_t n = size_t(1) << 20;
    for (const Row& row : rows) {
        row.params.validate();
        EXPECT_NEAR(row.params.distance(), row.distance, 5e-4) << "distance for r=" << row.params.r;
        EXPECT_EQ(row_sparsity_c(n, row.params), row.c_n) << "c_n for r=" << row.params.r;

        const double bound = row_sparsity_d_bound(n, row.params, 127.0);
        EXPECT_LT(bound, static_cast<double>(row.d_n) + 0.01) << "d_n bound for r=" << row.params.r;
        EXPECT_GT(bound, static_cast<double>(row.d_n) - 1.0) << "d_n bound for r=" << row.params.r;
        // The shipped integer is the paper's value or one above it, never below.
        const size_t sparsity = row_sparsity_d(n, row.params, 127.0);
        EXPECT_GE(sparsity, row.d_n) << "d_n for r=" << row.params.r;
        EXPECT_LE(sparsity, row.d_n + 1) << "d_n for r=" << row.params.r;
    }
}

TEST(BrakedownCodeTest, IsSystematicAndHasTheExpectedRate)
{
    const BrakedownParams params = BrakedownParams::distance_7pct();
    const size_t n = 1024;
    const BrakedownCode code(n, params);

    const std::vector<fr> message = random_message(n);
    const std::vector<fr> codeword = code.encode(message);

    ASSERT_EQ(codeword.size(), code.codeword_length());
    for (size_t i = 0; i < n; ++i) {
        EXPECT_EQ(codeword[i], message[i]) << "systematic prefix broken at " << i;
    }
    // Rate 1/r, up to the ceilings the recursion introduces at each level.
    const double rate = static_cast<double>(n) / static_cast<double>(code.codeword_length());
    EXPECT_GT(rate, params.rate() * 0.85);
    EXPECT_LE(rate, params.rate() * 1.02);
}

TEST(BrakedownCodeTest, EncodingIsLinear)
{
    const BrakedownCode code(512, BrakedownParams::distance_5pct());
    const std::vector<fr> a = random_message(512);
    const std::vector<fr> b = random_message(512);
    const fr scalar = fr::random_element();

    std::vector<fr> combination(512);
    for (size_t i = 0; i < 512; ++i) {
        combination[i] = a[i] + scalar * b[i];
    }

    const std::vector<fr> encoded_a = code.encode(a);
    const std::vector<fr> encoded_b = code.encode(b);
    const std::vector<fr> encoded_combination = code.encode(combination);

    for (size_t i = 0; i < encoded_a.size(); ++i) {
        EXPECT_EQ(encoded_combination[i], encoded_a[i] + scalar * encoded_b[i]) << "nonlinear at " << i;
    }
}

TEST(BrakedownCodeTest, SameSeedGivesTheSameCode)
{
    const std::vector<fr> message = random_message(256);
    const BrakedownCode first(256, BrakedownParams::distance_7pct(), /*seed=*/7);
    const BrakedownCode second(256, BrakedownParams::distance_7pct(), /*seed=*/7);
    const BrakedownCode other(256, BrakedownParams::distance_7pct(), /*seed=*/8);

    EXPECT_EQ(first.encode(message), second.encode(message));
    EXPECT_NE(first.encode(message), other.encode(message));
}

// The construction's guarantee is that a nonzero message encodes to a word of weight at least
// beta*n. Sparse and adversarially-shaped messages are the interesting cases: a single nonzero
// coordinate is exactly what the sparsity-preserving matrices exist to spread out.
TEST(BrakedownCodeTest, NonzeroMessagesMeetTheDistanceBound)
{
    const BrakedownParams params = BrakedownParams::distance_7pct();
    const size_t n = 1024;
    const BrakedownCode code(n, params);
    const auto floor = static_cast<size_t>(params.beta * static_cast<double>(n));

    std::vector<std::vector<fr>> messages;
    for (const size_t position : { size_t(0), size_t(1), n / 2, n - 1 }) {
        std::vector<fr> unit(n, fr::zero());
        unit[position] = fr::one();
        messages.push_back(std::move(unit));
    }
    std::vector<fr> pair(n, fr::zero());
    pair[3] = fr::one();
    pair[n - 4] = -fr::one();
    messages.push_back(std::move(pair));
    messages.push_back(random_message(n));

    for (size_t m = 0; m < messages.size(); ++m) {
        const std::vector<fr> codeword = code.encode(messages[m]);
        EXPECT_GE(hamming_weight(codeword), floor)
            << "message " << m << " encodes below the beta*n = " << floor << " floor";
    }
}

TEST(BrakedownCodeTest, EncodingCostIsLinearInTheMessage)
{
    const BrakedownParams params = BrakedownParams::distance_7pct();
    const BrakedownCode small(1024, params);
    const BrakedownCode large(16384, params);

    const double small_ratio = static_cast<double>(small.multiplications()) / 1024.0;
    const double large_ratio = static_cast<double>(large.multiplications()) / 16384.0;
    // A 16x longer message must not raise the per-symbol cost: that is the linear-time claim, and
    // it is what separates this code from the O(n log n) RS encoder it stands in for.
    EXPECT_LT(large_ratio, small_ratio * 1.2);
    EXPECT_LT(large_ratio, 40.0);
}

} // namespace bb::brakedown
