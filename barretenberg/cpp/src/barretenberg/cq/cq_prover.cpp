#include "cq_prover.hpp"

#include "barretenberg/common/assert.hpp"
#include "barretenberg/common/throw_or_abort.hpp"
#include "barretenberg/cq/cq_domain.hpp"
#include "barretenberg/ecc/scalar_multiplication/scalar_multiplication.hpp"
#include "barretenberg/polynomials/polynomial_arithmetic.hpp"

#include <map>

namespace bb::cq {

namespace {

g1::affine_element sparse_msm(std::span<const fr> scalars, std::span<const g1::affine_element> points)
{
    BB_ASSERT_EQ(scalars.size(), points.size(), "sparse_msm: scalars and points must have equal length");
    return g1::affine_element(scalar_multiplication::pippenger<curve::BN254>(
        PolynomialSpan<const fr>{ 0, scalars }, points, /*handle_edge_cases=*/true));
}

/**
 * @brief Product of two polynomials of degree < n, via FFTs over a domain of size 2n.
 */
std::vector<fr> polynomial_product(std::span<const fr> a, std::span<const fr> b, size_t n)
{
    const SubgroupDomain domain(2 * n);
    std::vector<fr> a_padded(2 * n, fr::zero());
    std::vector<fr> b_padded(2 * n, fr::zero());
    std::copy(a.begin(), a.end(), a_padded.begin());
    std::copy(b.begin(), b.end(), b_padded.begin());
    domain.fft<fr>(a_padded);
    domain.fft<fr>(b_padded);
    for (size_t k = 0; k < 2 * n; ++k) {
        a_padded[k] *= b_padded[k];
    }
    domain.ifft<fr>(a_padded);
    return a_padded;
}

/**
 * @brief Exact division of a polynomial of degree <= 2n-2 by X^n - 1. The quotient is simply the upper coefficient
 * block; exactness requires p_k + p_{k+n} = 0 for k < n-1 and p_{n-1} = 0.
 */
std::vector<fr> divide_by_vanishing(std::span<const fr> p, size_t n)
{
    std::vector<fr> quotient(p.begin() + static_cast<std::ptrdiff_t>(n),
                             p.begin() + static_cast<std::ptrdiff_t>(2 * n - 1));
    for (size_t k = 0; k + 1 < n; ++k) {
        BB_ASSERT(p[k] + p[k + n] == fr::zero(), "divide_by_vanishing: input not divisible by X^n - 1");
    }
    BB_ASSERT(p[n - 1] == fr::zero(), "divide_by_vanishing: input not divisible by X^n - 1");
    return quotient;
}

} // namespace

CqProver::Proof CqProver::construct_proof(const std::vector<std::vector<fr>>& lookup_columns,
                                          bool map_missing_rows_to_first_entry_for_testing)
{
    const size_t big_n = key_.table_size;
    const size_t num_columns = key_.num_columns();
    BB_ASSERT_EQ(lookup_columns.size(), num_columns, "CqProver: lookup columns must match the table columns");
    const size_t n = lookup_columns[0].size();
    for (const auto& column : lookup_columns) {
        BB_ASSERT_EQ(column.size(), n, "CqProver: all lookup columns must have the same size");
    }
    BB_ASSERT(n >= 2 && (n & (n - 1)) == 0, "CqProver: number of lookups must be a power of 2 and >= 2");
    BB_ASSERT_LTE(n, big_n, "CqProver: number of lookups must not exceed the table size");

    auto transcript = std::make_shared<Transcript>();
    const SubgroupDomain domain_v(n);

    // Round 0: bind the statement (table identity and sizes) into the transcript.
    transcript->send_to_verifier("CQ:table_size", static_cast<uint32_t>(big_n));
    transcript->send_to_verifier("CQ:num_lookups", static_cast<uint32_t>(n));
    transcript->send_to_verifier("CQ:num_columns", static_cast<uint32_t>(num_columns));
    for (size_t k = 0; k < num_columns; ++k) {
        transcript->send_to_verifier("CQ:T_" + std::to_string(k), key_.verification_key.table_commitments_g1[k]);
    }

    // Round 1a: multiplicity polynomial m (sparse over H) from the witness rows, and per-column witness
    // polynomials f^k (dense over V).
    std::map<uint32_t, fr> multiplicities;
    std::vector<fr> row(num_columns);
    for (size_t j = 0; j < n; ++j) {
        for (size_t k = 0; k < num_columns; ++k) {
            row[k] = lookup_columns[k][j];
        }
        std::optional<uint32_t> index = key_.find_row(row);
        if (!index.has_value()) {
            if (!map_missing_rows_to_first_entry_for_testing) {
                throw_or_abort("CqProver: lookup row not present in the table");
            }
            index = 0;
        }
        multiplicities[index.value()] += fr::one();
    }
    const size_t support_size = multiplicities.size();
    std::vector<uint32_t> support_indices;
    std::vector<fr> m_values;
    support_indices.reserve(support_size);
    m_values.reserve(support_size);
    for (const auto& [index, count] : multiplicities) {
        support_indices.push_back(index);
        m_values.push_back(count);
    }
    std::vector<g1::affine_element> support_lagrange_points(support_size);
    std::vector<g1::affine_element> support_zero_quotient_points(support_size);
    for (size_t s = 0; s < support_size; ++s) {
        support_lagrange_points[s] = key_.lagrange_commitments[support_indices[s]];
        support_zero_quotient_points[s] = key_.lagrange_zero_quotients[support_indices[s]];
    }
    const g1::affine_element m_commitment = sparse_msm(m_values, support_lagrange_points);
    transcript->send_to_verifier("CQ:m", m_commitment);

    std::vector<std::vector<fr>> f_column_monomials(num_columns);
    for (size_t k = 0; k < num_columns; ++k) {
        f_column_monomials[k].assign(lookup_columns[k].begin(), lookup_columns[k].end());
        domain_v.ifft<fr>(f_column_monomials[k]);
        const g1::affine_element f_commitment = ck_.commit(PolynomialSpan<const fr>{ 0, f_column_monomials[k] });
        transcript->send_to_verifier("CQ:f_" + std::to_string(k), f_commitment);
    }

    // Round 1b: the column-combination challenge theta reduces rows to scalars.
    const fr theta = transcript->template get_challenge<fr>("CQ:theta");
    std::vector<fr> theta_powers(num_columns);
    theta_powers[0] = fr::one();
    for (size_t k = 1; k < num_columns; ++k) {
        theta_powers[k] = theta_powers[k - 1] * theta;
    }
    std::vector<fr> combined_lookups(n, fr::zero());
    for (size_t k = 0; k < num_columns; ++k) {
        for (size_t j = 0; j < n; ++j) {
            combined_lookups[j] += theta_powers[k] * lookup_columns[k][j];
        }
    }
    std::vector<fr> f_monomial(n, fr::zero());
    for (size_t k = 0; k < num_columns; ++k) {
        for (size_t d = 0; d < n; ++d) {
            f_monomial[d] += theta_powers[k] * f_column_monomials[k][d];
        }
    }
    std::vector<fr> combined_table_on_support(support_size, fr::zero());
    for (size_t k = 0; k < num_columns; ++k) {
        for (size_t s = 0; s < support_size; ++s) {
            combined_table_on_support[s] += theta_powers[k] * key_.columns[k][support_indices[s]];
        }
    }

    const fr beta = transcript->template get_challenge<fr>("CQ:beta");

    // Round 2, table side: A_i = m_i / (beta + t_i(theta)) on the support; commitments via sparse MSMs over cached
    // points. The cached quotients combine homomorphically: [Q_A] = sum_k theta^k sum_i A_i [Q^k_i].
    std::vector<fr> a_values(support_size);
    for (size_t s = 0; s < support_size; ++s) {
        a_values[s] = beta + combined_table_on_support[s];
    }
    fr::batch_invert(a_values.data(), support_size);
    fr a_sum = fr::zero();
    for (size_t s = 0; s < support_size; ++s) {
        a_values[s] *= m_values[s];
        a_sum += a_values[s];
    }
    const fr a_zero = a_sum * fr(big_n).invert();
    const g1::affine_element a_commitment = sparse_msm(a_values, support_lagrange_points);
    g1::element q_a_accumulator;
    q_a_accumulator.self_set_infinity();
    std::vector<g1::affine_element> support_quotient_points(support_size);
    for (size_t k = 0; k < num_columns; ++k) {
        for (size_t s = 0; s < support_size; ++s) {
            support_quotient_points[s] = key_.cached_quotients[k][support_indices[s]];
        }
        q_a_accumulator += g1::element(sparse_msm(a_values, support_quotient_points)) * theta_powers[k];
    }
    const g1::affine_element q_a_commitment(q_a_accumulator);

    // Round 2, witness side: B_j = 1/(beta + f_j(theta)) dense over V, plus the quotient Q_B and the degree-check
    // shift.
    std::vector<fr> b_values(n);
    for (size_t j = 0; j < n; ++j) {
        b_values[j] = beta + combined_lookups[j];
    }
    fr::batch_invert(b_values.data(), n);
    std::vector<fr> b_monomial = b_values;
    domain_v.ifft<fr>(b_monomial);
    const fr b_zero = b_monomial[0];
    if (!map_missing_rows_to_first_entry_for_testing) {
        BB_ASSERT(fr(big_n) * a_zero == fr(n) * b_zero, "CqProver: log-derivative sums do not match");
    }

    const std::vector<fr> b_shifted(b_monomial.begin() + 1, b_monomial.end()); // B_0(X) = (B(X) - B(0))/X

    std::vector<fr> product = polynomial_product(b_monomial, f_monomial, n);
    for (size_t k = 0; k < n; ++k) {
        product[k] += beta * b_monomial[k];
    }
    product[0] -= fr::one();
    const std::vector<fr> q_b = divide_by_vanishing(product, n);

    const g1::affine_element b_shifted_commitment = ck_.commit(PolynomialSpan<const fr>{ 0, b_shifted });
    const g1::affine_element q_b_commitment = ck_.commit(PolynomialSpan<const fr>{ 0, q_b });
    // P_B = B_0 * X^{N-n+1}: committable iff deg(B_0) <= n-2, which is the degree bound the verifier checks.
    const g1::affine_element degree_shift_commitment = ck_.commit(PolynomialSpan<const fr>{ big_n - n + 1, b_shifted });

    transcript->send_to_verifier("CQ:A", a_commitment);
    transcript->send_to_verifier("CQ:Q_A", q_a_commitment);
    transcript->send_to_verifier("CQ:B_0", b_shifted_commitment);
    transcript->send_to_verifier("CQ:Q_B", q_b_commitment);
    transcript->send_to_verifier("CQ:P_B", degree_shift_commitment);
    const fr gamma = transcript->template get_challenge<fr>("CQ:gamma");

    // Round 3: evaluations at gamma, then the batched KZG opening at gamma and the opening of A at zero.
    const fr b_shifted_at_gamma = polynomial_arithmetic::evaluate<fr>(b_shifted, gamma);
    const fr f_at_gamma = polynomial_arithmetic::evaluate<fr>(f_monomial, gamma);
    transcript->send_to_verifier("CQ:b0_gamma", b_shifted_at_gamma);
    transcript->send_to_verifier("CQ:f_gamma", f_at_gamma);
    transcript->send_to_verifier("CQ:a0", a_zero);
    const fr eta = transcript->template get_challenge<fr>("CQ:eta");

    const fr q_b_at_gamma = polynomial_arithmetic::evaluate<fr>(q_b, gamma);
    std::vector<fr> batched(n, fr::zero());
    for (size_t k = 0; k + 1 < n; ++k) {
        batched[k] += b_shifted[k];
    }
    for (size_t k = 0; k < n; ++k) {
        batched[k] += eta * f_monomial[k];
    }
    const fr eta_sqr = eta.sqr();
    for (size_t k = 0; k + 1 < n; ++k) {
        batched[k] += eta_sqr * q_b[k];
    }
    batched[0] -= b_shifted_at_gamma + eta * f_at_gamma + eta_sqr * q_b_at_gamma;
    polynomial_arithmetic::factor_roots<fr>(batched, gamma);
    const g1::affine_element opening_at_gamma = ck_.commit(PolynomialSpan<const fr>{ 0, batched });
    const g1::affine_element opening_at_zero = sparse_msm(a_values, support_zero_quotient_points);

    transcript->send_to_verifier("CQ:pi_gamma", opening_at_gamma);
    transcript->send_to_verifier("CQ:pi_zero", opening_at_zero);

    return transcript->export_proof();
}

} // namespace bb::cq
