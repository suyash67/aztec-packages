#pragma once

#include "barretenberg/commitment_schemes/commitment_key.hpp"
#include "barretenberg/cq/cq_key.hpp"
#include "barretenberg/transcript/transcript.hpp"

#include <span>
#include <vector>

namespace bb::cq {

/**
 * @brief Prover for the cq lookup argument (Eagen-Fiore-Gabizon, https://eprint.iacr.org/2022/1763), extended to
 * multi-column (vector) lookups.
 *
 * @details Proves that every row of a witness of n K-tuples appears among the N rows of the preprocessed table.
 * All proving work is O(K n log n) field operations plus MSMs of size O(n); nothing depends on N, which is the
 * point: tables far larger than the circuit become affordable because their cost is paid once, at preprocessing.
 *
 * Rows are reduced to single field elements with a challenge theta drawn after the witness column commitments:
 * t_i(theta) = sum_k theta^k t^k_i and f_j(theta) = sum_k theta^k f^k_j. The scalar protocol is the log-derivative
 * identity sum_i m_i/(beta + t_i) = sum_j 1/(beta + f_j) checked via committed polynomials: the left side over the
 * table domain H (|H| = N) using cached quotient commitments for the product identity
 * A(X)(T(X) + beta) - m(X) = Q_A(X) Z_H(X), and the right side over the witness domain V (|V| = n) using an
 * ordinary dense quotient B(X)(f(X) + beta) - 1 = Q_B(X) Z_V(X). The two sides are linked through the sums
 * N * A(0) = n * B(0). Non-zero-knowledge variant.
 */
class CqProver {
  public:
    using Transcript = NativeTranscript;
    using Proof = Transcript::Proof;

    CqProver(const CqProvingKey& proving_key, const CommitmentKey<curve::BN254>& commitment_key)
        : key_(proving_key)
        , ck_(commitment_key)
    {}

    /**
     * @param lookup_columns witness values, one vector per table column, all of the same power-of-2 size in [2, N];
     *        row j (the tuple of the j-th entries) must appear as a row of the table
     * @param map_missing_rows_to_first_entry_for_testing soundness-test hook: instead of aborting on a witness row
     *        that is not in the table, count it toward the multiplicity of table row 0, producing a proof an honest
     *        verifier must reject
     */
    Proof construct_proof(const std::vector<std::vector<fr>>& lookup_columns,
                          bool map_missing_rows_to_first_entry_for_testing = false);

    /** @brief Single-column convenience overload. */
    Proof construct_proof(std::span<const fr> lookups, bool map_missing_rows_to_first_entry_for_testing = false)
    {
        std::vector<std::vector<fr>> columns;
        columns.emplace_back(lookups.begin(), lookups.end());
        return construct_proof(columns, map_missing_rows_to_first_entry_for_testing);
    }

  private:
    const CqProvingKey& key_;
    const CommitmentKey<curve::BN254>& ck_;
};

} // namespace bb::cq
