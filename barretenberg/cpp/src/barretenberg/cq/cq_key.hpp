#pragma once

#include "barretenberg/cq/cq_domain.hpp"
#include "barretenberg/ecc/curves/bn254/fr.hpp"
#include "barretenberg/ecc/curves/bn254/g1.hpp"
#include "barretenberg/ecc/curves/bn254/g2.hpp"
#include "barretenberg/numeric/uint256/uint256.hpp"

#include <map>
#include <optional>
#include <span>
#include <vector>

namespace bb::cq {

/**
 * @brief Verifier-side preprocessed data for a cq lookup table of one or more columns.
 *
 * @details A row of a K-column table is reduced to the single field element t_i(theta) = sum_k theta^k t^k_i for a
 * challenge theta drawn after the witness commitments; commitments to the combined table polynomial are obtained
 * homomorphically from the per-column commitments stored here. Verification needs O(1) group elements per proof:
 * the [T^k]_2, [tau^N - 1]_2, [1]_2, [tau]_2 and the degree-check power [tau^{N-n+1}]_2 (which depends on the
 * number of lookups n). The full G2 power vector is kept so any n can be served; a production verification key for
 * a fixed circuit size would store the single required power.
 */
struct CqVerificationKey {
    size_t table_size = 0;
    std::vector<g1::affine_element> table_commitments_g1; // [T^k]_1 per column, binds the table in the transcript
    std::vector<g2::affine_element> table_commitments_g2; // [T^k(tau)]_2 per column
    g2::affine_element vanishing_commitment_g2;           // [tau^N - 1]_2
    std::vector<g2::affine_element> g2_powers;            // [tau^k]_2, k = 0..N

    size_t num_columns() const { return table_commitments_g2.size(); }
    g2::affine_element g2_identity() const { return g2_powers[0]; }
    g2::affine_element g2_x() const { return g2_powers[1]; }
    /** @brief [tau^{N-n+1}]_2, pairing counterpart of the shifted commitment [B_0 * X^{N-n+1}]_1. */
    g2::affine_element degree_check_power(size_t num_lookups) const { return g2_powers[table_size - num_lookups + 1]; }
};

/**
 * @brief Prover-side preprocessed data for a cq lookup table: O(N) cached commitments per column enabling proofs
 * whose cost depends only on the number of lookups, never on the table size N.
 *
 * @details Preprocessing costs O(N log N) group operations per column (via Feist-Khovratovich batched KZG openings
 * and group FFTs) and is done once per table. Afterwards each proof uses only sparse MSMs over the cached points.
 */
class CqProvingKey {
  public:
    size_t table_size = 0;
    std::vector<std::vector<fr>> columns;                          // t^k_i, k = 0..K-1, i = 0..N-1
    std::vector<std::vector<fr>> column_monomials;                 // coefficients of T^k(X) interpolating t^k over H
    std::vector<g1::affine_element> lagrange_commitments;          // [L_i(tau)]_1 (column-independent)
    std::vector<g1::affine_element> lagrange_zero_quotients;       // [(L_i(tau) - L_i(0)) / tau]_1 (column-independent)
    std::vector<std::vector<g1::affine_element>> cached_quotients; // [Q^k_i(tau)]_1, Q^k_i = L_i (T^k - t^k_i) / Z_H
    CqVerificationKey verification_key;

    size_t num_columns() const { return columns.size(); }

    /**
     * @param table_columns the table, one vector per column, all of the same power-of-2 size (pad with a repeated
     *        row if needed)
     * @param g1_powers monomial G1 SRS [tau^d]_1, d = 0..N-1. Exactly N points should exist for a production
     *        deployment (the truncation is what degree-bounds the prover's committed polynomials).
     * @param g2_powers G2 SRS [tau^k]_2, k = 0..N
     */
    static CqProvingKey create(std::vector<std::vector<fr>> table_columns,
                               std::span<const g1::affine_element> g1_powers,
                               std::span<const g2::affine_element> g2_powers);

    /** @brief Single-column convenience overload. */
    static CqProvingKey create(std::vector<fr> table,
                               std::span<const g1::affine_element> g1_powers,
                               std::span<const g2::affine_element> g2_powers)
    {
        std::vector<std::vector<fr>> table_columns;
        table_columns.emplace_back(std::move(table));
        return create(std::move(table_columns), g1_powers, g2_powers);
    }

    /** @brief Index of a row in the table, if present. */
    std::optional<uint32_t> find_row(std::span<const fr> row) const
    {
        auto it = row_index_.find(row_key(row));
        return it == row_index_.end() ? std::nullopt : std::optional<uint32_t>(it->second);
    }

  private:
    static std::vector<uint256_t> row_key(std::span<const fr> row)
    {
        std::vector<uint256_t> key;
        key.reserve(row.size());
        for (const fr& value : row) {
            key.emplace_back(static_cast<uint256_t>(value));
        }
        return key;
    }

    std::map<std::vector<uint256_t>, uint32_t> row_index_;
};

} // namespace bb::cq
