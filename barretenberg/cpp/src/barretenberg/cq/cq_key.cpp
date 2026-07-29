#include "cq_key.hpp"

#include "barretenberg/common/assert.hpp"
#include "barretenberg/common/thread.hpp"
#include "barretenberg/cq/cq_fk.hpp"
#include "barretenberg/ecc/scalar_multiplication/scalar_multiplication.hpp"

namespace bb::cq {

namespace {

/**
 * @brief G2 MSM via per-thread accumulation. There is no pippenger over G2 in barretenberg; this is a one-time
 * preprocessing cost so a parallel naive MSM suffices.
 */
g2::affine_element g2_msm(std::span<const fr> scalars, std::span<const g2::affine_element> points)
{
    const size_t num_chunks = get_num_cpus();
    const size_t n = scalars.size();
    const size_t chunk_size = (n + num_chunks - 1) / num_chunks;
    std::vector<g2::element> partial_sums(num_chunks);
    parallel_for(num_chunks, [&](size_t chunk) {
        g2::element acc;
        acc.self_set_infinity();
        const size_t start = chunk * chunk_size;
        const size_t end = std::min(start + chunk_size, n);
        for (size_t i = start; i < end; ++i) {
            if (!scalars[i].is_zero()) {
                acc += g2::element(points[i]) * scalars[i];
            }
        }
        partial_sums[chunk] = acc;
    });
    g2::element total;
    total.self_set_infinity();
    for (const auto& partial : partial_sums) {
        total += partial;
    }
    return g2::affine_element(total);
}

std::vector<g1::affine_element> to_affine_batch(std::vector<g1::element>&& points)
{
    g1::element::batch_normalize(points.data(), points.size());
    std::vector<g1::affine_element> result(points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        result[i] = points[i].is_point_at_infinity() ? g1::affine_element::infinity()
                                                     : g1::affine_element(points[i].x, points[i].y);
    }
    return result;
}

} // namespace

CqProvingKey CqProvingKey::create(std::vector<std::vector<fr>> table_columns,
                                  std::span<const g1::affine_element> g1_powers,
                                  std::span<const g2::affine_element> g2_powers)
{
    BB_ASSERT(!table_columns.empty(), "CqProvingKey: at least one table column required");
    const size_t n = table_columns[0].size();
    const size_t num_columns = table_columns.size();
    for (const auto& column : table_columns) {
        BB_ASSERT_EQ(column.size(), n, "CqProvingKey: all table columns must have the same size");
    }
    BB_ASSERT(n >= 4 && (n & (n - 1)) == 0, "CqProvingKey: table size must be a power of 2 and >= 4");
    BB_ASSERT_GTE(g1_powers.size(), n, "CqProvingKey: G1 SRS must have at least N points");
    BB_ASSERT_GTE(g2_powers.size(), n + 1, "CqProvingKey: G2 SRS must have at least N+1 points");

    CqProvingKey key;
    key.table_size = n;
    key.columns = std::move(table_columns);

    const SubgroupDomain domain(n);
    const std::vector<fr> root_powers = domain.root_powers();

    // Column-independent cached KZG data: Lagrange commitments and Lagrange zero-opening quotients.
    key.lagrange_commitments = lagrange_basis_commitments(g1_powers, n);
    key.lagrange_zero_quotients = lagrange_zero_quotient_commitments(g1_powers, n);

    key.verification_key.table_size = n;
    key.verification_key.vanishing_commitment_g2 =
        g2::affine_element(g2::element(g2_powers[n]) - g2::element(g2_powers[0]));
    key.verification_key.g2_powers.assign(g2_powers.begin(), g2_powers.end());

    // Per-column data: the interpolating polynomial T^k, its G1/G2 commitments, and all cached quotients
    // [Q^k_i] = (omega^i / N) * [(T^k - t^k_i)/(X - omega^i)], each O(N log N) via group FFTs.
    for (size_t k = 0; k < num_columns; ++k) {
        std::vector<fr> monomial = key.columns[k];
        domain.ifft<fr>(monomial);

        std::vector<g1::affine_element> opening_proofs = fk_all_kzg_opening_proofs(monomial, g1_powers);
        std::vector<g1::element> scaled(n);
        parallel_for_range(n, [&](size_t start, size_t end) {
            for (size_t i = start; i < end; ++i) {
                scaled[i] = g1::element(opening_proofs[i]) * (root_powers[i] * domain.size_inverse());
            }
        });
        key.cached_quotients.emplace_back(to_affine_batch(std::move(scaled)));

        key.verification_key.table_commitments_g1.emplace_back(scalar_multiplication::pippenger<curve::BN254>(
            PolynomialSpan<const fr>{ 0, monomial }, g1_powers, /*handle_edge_cases=*/true));
        key.verification_key.table_commitments_g2.emplace_back(g2_msm(monomial, g2_powers.subspan(0, n)));
        key.column_monomials.emplace_back(std::move(monomial));
    }

    std::vector<fr> row(num_columns);
    for (uint32_t i = 0; i < n; ++i) {
        for (size_t k = 0; k < num_columns; ++k) {
            row[k] = key.columns[k][i];
        }
        key.row_index_.emplace(row_key(row), i);
    }

    return key;
}

} // namespace bb::cq
