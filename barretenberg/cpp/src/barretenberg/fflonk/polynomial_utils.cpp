#include "barretenberg/fflonk/polynomial_utils.hpp"

#include "barretenberg/common/thread.hpp"

#include <algorithm>
#include <ranges>

namespace bb::fflonk_plonk {

std::vector<FF> pack_columns(std::span<const std::vector<FF>> columns)
{
    const size_t t = columns.size();
    BB_ASSERT_GT(t, size_t(0), "a packed group needs at least one column");

    size_t longest = 0;
    for (const std::vector<FF>& column : columns) {
        longest = std::max(longest, column.size());
    }

    std::vector<FF> packed(t * longest, FF::zero());
    for (size_t i = 0; i < t; ++i) {
        const std::vector<FF>& column = columns[i];
        for (size_t j = 0; j < column.size(); ++j) {
            packed[j * t + i] = column[j];
        }
    }
    return packed;
}

void divide_by_power_minus(
    std::span<const FF> poly, const size_t t, const FF& c, std::vector<FF>& quotient, std::vector<FF>& residue)
{
    BB_ASSERT_GT(t, size_t(0), "cannot divide by X^0 - c");
    const size_t size = poly.size();

    quotient.assign(size, FF::zero());
    for (size_t j = size; j-- > 0;) {
        if (j + t < size) {
            quotient[j] = poly[j + t] + c * quotient[j + t];
        }
    }

    residue.assign(t, FF::zero());
    for (size_t i = 0; i < t && i < size; ++i) {
        residue[i] = poly[i] + c * quotient[i];
    }
}

std::vector<FF> divide_by_linear(std::span<const FF> poly, const FF& r, FF& remainder)
{
    const size_t size = poly.size();
    if (size == 0) {
        remainder = FF::zero();
        return {};
    }

    std::vector<FF> quotient(size - 1, FF::zero());
    FF accumulator = FF::zero();
    for (size_t i = size; i-- > 0;) {
        const FF next = poly[i] + r * accumulator;
        if (i > 0) {
            quotient[i - 1] = next;
        } else {
            remainder = next;
        }
        accumulator = next;
    }
    return quotient;
}

bool divide_by_vanishing(std::span<const FF> poly, const size_t n, std::vector<FF>& quotient)
{
    std::vector<FF> residue;
    divide_by_power_minus(poly, n, FF::one(), quotient, residue);

    for (const FF& coefficient : residue) {
        if (!coefficient.is_zero()) {
            return false;
        }
    }
    // The quotient of a degree-d polynomial by X^n - 1 has degree d - n; the top n slots of the
    // buffer are structurally zero and only exist because the recurrence is written in place.
    quotient.resize(quotient.size() > n ? quotient.size() - n : static_cast<size_t>(0));
    return true;
}

void trim(std::vector<FF>& poly)
{
    while (poly.size() > 1 && poly.back().is_zero()) {
        poly.pop_back();
    }
}

FF evaluate(std::span<const FF> poly, const FF& x)
{
    FF accumulator = FF::zero();
    for (size_t i = poly.size(); i-- > 0;) {
        accumulator = accumulator * x + poly[i];
    }
    return accumulator;
}

std::vector<FF> coefficients_to_evaluations(std::span<const FF> coefficients, const EvaluationDomain<FF>& domain)
{
    BB_ASSERT_LTE(coefficients.size(), domain.size, "polynomial does not fit the evaluation domain");

    std::vector<FF> padded(domain.size, FF::zero());
    std::ranges::copy(coefficients, padded.begin());

    std::vector<FF> evaluations(domain.size, FF::zero());
    polynomial_arithmetic::fft_inner_parallel(
        padded.data(), evaluations.data(), domain, domain.root, domain.get_round_roots());
    return evaluations;
}

std::vector<FF> evaluations_to_coefficients(std::span<const FF> evaluations, const EvaluationDomain<FF>& domain)
{
    BB_ASSERT_EQ(evaluations.size(), domain.size, "inverse FFT needs one evaluation per domain point");

    std::vector<FF> source(evaluations.begin(), evaluations.end());
    std::vector<FF> coefficients(domain.size, FF::zero());
    polynomial_arithmetic::ifft(source.data(), coefficients.data(), domain);
    return coefficients;
}

std::vector<FF> compute_power_table(const FF& root, const size_t count)
{
    std::vector<FF> powers(count, FF::zero());
    if (count == 0) {
        return powers;
    }

    // Chunked so the running product stays sequential inside a chunk while chunks run in parallel;
    // each chunk pays one exponentiation to find its own starting power.
    const size_t num_chunks = std::min(count, static_cast<size_t>(get_num_cpus()));
    const size_t chunk_size = (count + num_chunks - 1) / num_chunks;
    parallel_for(num_chunks, [&](size_t chunk) {
        const size_t start = chunk * chunk_size;
        const size_t end = std::min(start + chunk_size, count);
        if (start >= end) {
            return;
        }
        FF accumulator = root.pow(static_cast<uint64_t>(start));
        for (size_t i = start; i < end; ++i) {
            powers[i] = accumulator;
            accumulator *= root;
        }
    });
    return powers;
}

} // namespace bb::fflonk_plonk
