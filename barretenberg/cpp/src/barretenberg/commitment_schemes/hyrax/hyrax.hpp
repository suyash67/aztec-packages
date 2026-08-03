#pragma once

#include "barretenberg/commitment_schemes/whir/weights.hpp"
#include "barretenberg/common/assert.hpp"
#include "barretenberg/common/thread.hpp"
#include "barretenberg/ecc/curves/bn254/bn254.hpp"
#include "barretenberg/polynomials/polynomial.hpp"

#include <map>
#include <memory>
#include <span>
#include <string>
#include <vector>

namespace bb::hyrax {

using Curve = curve::BN254;
using Commitment = Curve::AffineElement;
using GroupElement = Curve::Element;

/** @brief Reference to one polynomial of one committed group (duck-type shared with the suite). */
struct HyraxColumnRef {
    size_t group;
    size_t column;
};

/** @brief Hyrax parameters: the square-ish matrix split (columns = low variables). */
struct HyraxConfig {
    size_t num_variables;
    size_t log_num_cols;

    size_t num_rows() const { return size_t(1) << (num_variables - log_num_cols); }
    size_t num_cols() const { return size_t(1) << log_num_cols; }

    static HyraxConfig create(size_t num_variables, size_t /*security_bits*/ = 100, size_t /*log_inv_rate*/ = 0)
    {
        BB_ASSERT_GT(num_variables, size_t(1));
        return { num_variables, (num_variables + 1) / 2 };
    }
};

/** @brief Prover-side commitment to a group: dense arrays and per-row Pedersen commitments. */
struct HyraxGroupData {
    std::vector<std::vector<fr>> coefficients;
    // commitments[c][r]: row r of polynomial c
    std::vector<std::vector<Commitment>> commitments;

    size_t num_columns() const { return coefficients.size(); }
};

/** @brief Transparent commitment key: hash-derived generators for one matrix row. */
class HyraxCommitmentKey {
  public:
    explicit HyraxCommitmentKey(const HyraxConfig& config)
        : config(config)
        , generators(Curve::Group::derive_generators(std::vector<uint8_t>{ 'b', 'b', '_', 'h', 'y', 'r', 'a', 'x' },
                                                     config.num_cols()))
    {}

    HyraxGroupData commit_group(std::vector<std::vector<fr>> payload_columns,
                                const std::vector<bool>& /*to_be_shifted*/ = {}) const
    {
        const size_t n = size_t(1) << config.num_variables;
        const size_t num_rows = config.num_rows();
        const size_t num_cols = config.num_cols();
        HyraxGroupData data;
        for (auto& column : payload_columns) {
            BB_ASSERT_LTE(column.size(), n, "polynomial too large for the configured size");
            column.resize(n, fr::zero());
            std::vector<Commitment> row_commitments(num_rows);
            for (size_t r = 0; r < num_rows; ++r) {
                std::vector<fr> row(column.begin() + static_cast<std::ptrdiff_t>(r * num_cols),
                                    column.begin() + static_cast<std::ptrdiff_t>((r + 1) * num_cols));
                row_commitments[r] = Commitment::batch_mul(std::span<const Commitment>(generators), std::span<fr>(row));
            }
            data.commitments.push_back(std::move(row_commitments));
            data.coefficients.push_back(std::move(column));
        }
        return data;
    }

    HyraxGroupData commit_group(std::span<const Polynomial<fr>* const> polynomials,
                                const std::vector<bool>& to_be_shifted = {}) const
    {
        std::vector<std::vector<fr>> payload_columns;
        for (const Polynomial<fr>* polynomial : polynomials) {
            std::vector<fr> column(polynomial->end_index(), fr::zero());
            for (size_t i = polynomial->start_index(); i < polynomial->end_index(); ++i) {
                column[i] = (*polynomial)[i];
            }
            payload_columns.push_back(std::move(column));
        }
        return commit_group(std::move(payload_columns), to_be_shifted);
    }

    HyraxConfig config;
    std::vector<Commitment> generators;
};

namespace detail {

inline std::string hyrax_label(const std::string& name, size_t i)
{
    return "HYRAX:" + name + "_" + std::to_string(i);
}

} // namespace detail

/**
 * @brief Hyrax opening prover: the matrix tensor sandwich v = b^T M a with Pedersen row
 * commitments. The prover sends three combined rows (unshifted; shifted part 1; shifted part 2 of
 * the rank-2 shift decomposition, exactly as in the Ligero backend); the verifier checks each
 * against the homomorphically folded row commitments and evaluates the claim equation. No queries
 * and no extra challenges: binding is the Pedersen homomorphism itself, which is what makes the
 * proof sqrt(N) field elements with a deterministic check.
 */
class HyraxProver {
  public:
    struct Claims {
        std::vector<const HyraxGroupData*> groups;
        std::vector<HyraxColumnRef> unshifted;
        std::vector<fr> unshifted_evaluations;
        std::vector<HyraxColumnRef> to_be_shifted;
        std::vector<fr> shifted_evaluations;
        bool send_roots = true;
    };

    template <typename Transcript>
    static void prove(const HyraxCommitmentKey& ck,
                      const Claims& claims,
                      std::span<const fr> u,
                      const std::shared_ptr<Transcript>& transcript)
    {
        const HyraxConfig& config = ck.config;
        const size_t num_rows = config.num_rows();
        const size_t num_cols = config.num_cols();
        BB_ASSERT_EQ(u.size(), config.num_variables, "opening point size mismatch");

        if (claims.send_roots) {
            for (size_t g = 0; g < claims.groups.size(); ++g) {
                for (size_t c = 0; c < claims.groups[g]->num_columns(); ++c) {
                    for (size_t r = 0; r < num_rows; ++r) {
                        transcript->send_to_verifier(detail::hyrax_label("root", g) + "_" + std::to_string(c) + "_" +
                                                         std::to_string(r),
                                                     claims.groups[g]->commitments[c][r]);
                    }
                }
            }
        }
        const fr rho = transcript->template get_challenge<fr>("HYRAX:rho");

        const std::vector<fr> b = whir::eq_tensor(u.subspan(config.log_num_cols));

        // Per-(group,column) scalars, then the three combined rows in one pass (Ligero layout).
        struct RowTask {
            const std::vector<fr>* dense;
            fr scalar_u;
            fr scalar_s;
        };
        std::map<std::pair<size_t, size_t>, RowTask> tasks;
        fr rho_power = fr::one();
        for (const HyraxColumnRef& ref : claims.unshifted) {
            auto& task = tasks[{ ref.group, ref.column }];
            task.dense = &claims.groups[ref.group]->coefficients[ref.column];
            task.scalar_u += rho_power;
            rho_power *= rho;
        }
        for (const HyraxColumnRef& ref : claims.to_be_shifted) {
            auto& task = tasks[{ ref.group, ref.column }];
            task.dense = &claims.groups[ref.group]->coefficients[ref.column];
            task.scalar_s += rho_power;
            rho_power *= rho;
        }

        std::vector<fr> w_u(num_cols, fr::zero());
        std::vector<fr> w_s(num_cols, fr::zero());
        std::vector<fr> w_s2(num_cols, fr::zero());
        // Threads split the column range, so each owns a disjoint slice of all three combined rows.
        std::vector<RowTask> task_list;
        task_list.reserve(tasks.size());
        for (const auto& [key, task] : tasks) {
            task_list.push_back(task);
        }
        parallel_for_range(num_cols, [&](size_t start, size_t end) {
            for (const RowTask& task : task_list) {
                for (size_t r = 0; r < num_rows; ++r) {
                    const fr* row = task.dense->data() + r * num_cols;
                    const fr cu = task.scalar_u * b[r];
                    const fr cs = task.scalar_s * b[r];
                    const fr cs2 = (r > 0) ? task.scalar_s * b[r - 1] : fr::zero();
                    for (size_t c = start; c < end; ++c) {
                        w_u[c] += cu * row[c];
                        w_s[c] += cs * row[c];
                        w_s2[c] += cs2 * row[c];
                    }
                }
            }
        });
        for (size_t c = 0; c < num_cols; ++c) {
            transcript->send_to_verifier(detail::hyrax_label("w_u", c), w_u[c]);
            transcript->send_to_verifier(detail::hyrax_label("w_s", c), w_s[c]);
            transcript->send_to_verifier(detail::hyrax_label("w_s2", c), w_s2[c]);
        }
    }
};

/** @brief The mirror of `HyraxProver`: two MSM families and the rank-2 claim equation. */
class HyraxVerifier {
  public:
    struct Claims {
        std::vector<size_t> group_num_columns;
        // commitments[g][c][r]; when non-empty, transcript-bound
        std::vector<std::vector<std::vector<Commitment>>> group_commitments;
        std::vector<HyraxColumnRef> unshifted;
        std::vector<fr> unshifted_evaluations;
        std::vector<HyraxColumnRef> to_be_shifted;
        std::vector<fr> shifted_evaluations;
    };

    template <typename Transcript>
    static bool verify(const HyraxConfig& config,
                       const Claims& claims,
                       std::span<const fr> u,
                       const std::shared_ptr<Transcript>& transcript,
                       const HyraxCommitmentKey& ck)
    {
        const size_t num_rows = config.num_rows();
        const size_t num_cols = config.num_cols();
        BB_ASSERT_EQ(u.size(), config.num_variables, "opening point size mismatch");

        std::vector<std::vector<std::vector<Commitment>>> group_commitments = claims.group_commitments;
        if (group_commitments.empty()) {
            for (size_t g = 0; g < claims.group_num_columns.size(); ++g) {
                std::vector<std::vector<Commitment>> group;
                for (size_t c = 0; c < claims.group_num_columns[g]; ++c) {
                    std::vector<Commitment> rows;
                    for (size_t r = 0; r < num_rows; ++r) {
                        rows.push_back(transcript->template receive_from_prover<Commitment>(
                            detail::hyrax_label("root", g) + "_" + std::to_string(c) + "_" + std::to_string(r)));
                    }
                    group.push_back(std::move(rows));
                }
                group_commitments.push_back(std::move(group));
            }
        }
        const fr rho = transcript->template get_challenge<fr>("HYRAX:rho");

        const std::vector<fr> a = whir::eq_tensor(u.subspan(0, config.log_num_cols));
        const std::vector<fr> b = whir::eq_tensor(u.subspan(config.log_num_cols));

        // Expected combined-row commitments, as one MSM per w-vector over all claimed rows.
        std::vector<Commitment> points_u;
        std::vector<fr> scalars_u;
        std::vector<Commitment> points_s;
        std::vector<fr> scalars_s;
        std::vector<Commitment> points_s2;
        std::vector<fr> scalars_s2;
        fr expected = fr::zero();
        fr rho_power = fr::one();
        for (size_t i = 0; i < claims.unshifted.size(); ++i) {
            const auto& rows = group_commitments[claims.unshifted[i].group][claims.unshifted[i].column];
            for (size_t r = 0; r < num_rows; ++r) {
                points_u.push_back(rows[r]);
                scalars_u.push_back(rho_power * b[r]);
            }
            expected += rho_power * claims.unshifted_evaluations[i];
            rho_power *= rho;
        }
        for (size_t l = 0; l < claims.to_be_shifted.size(); ++l) {
            const auto& rows = group_commitments[claims.to_be_shifted[l].group][claims.to_be_shifted[l].column];
            for (size_t r = 0; r < num_rows; ++r) {
                points_s.push_back(rows[r]);
                scalars_s.push_back(rho_power * b[r]);
                if (r > 0) {
                    points_s2.push_back(rows[r]);
                    scalars_s2.push_back(rho_power * b[r - 1]);
                }
            }
            expected += rho_power * claims.shifted_evaluations[l];
            rho_power *= rho;
        }

        std::vector<fr> w_u(num_cols);
        std::vector<fr> w_s(num_cols);
        std::vector<fr> w_s2(num_cols);
        for (size_t c = 0; c < num_cols; ++c) {
            w_u[c] = transcript->template receive_from_prover<fr>(detail::hyrax_label("w_u", c));
            w_s[c] = transcript->template receive_from_prover<fr>(detail::hyrax_label("w_s", c));
            w_s2[c] = transcript->template receive_from_prover<fr>(detail::hyrax_label("w_s2", c));
        }

        // Commitment consistency: Commit(w_x) must equal the folded row commitments.
        if (!check_combined_row(ck, w_u, points_u, scalars_u)) {
            return false;
        }
        if (!claims.to_be_shifted.empty()) {
            if (!check_combined_row(ck, w_s, points_s, scalars_s)) {
                return false;
            }
            if (!check_combined_row(ck, w_s2, points_s2, scalars_s2)) {
                return false;
            }
        }

        // Claim equation with the rank-2 shift decomposition (Ligero README §2).
        fr total = fr::zero();
        for (size_t c = 0; c < num_cols; ++c) {
            total += w_u[c] * a[c];
            if (c >= 1) {
                total += w_s[c] * a[c - 1];
            }
        }
        total += a[num_cols - 1] * w_s2[0];
        return total == expected;
    }

  private:
    static bool check_combined_row(const HyraxCommitmentKey& ck,
                                   const std::vector<fr>& w,
                                   const std::vector<Commitment>& points,
                                   const std::vector<fr>& scalars)
    {
        if (points.empty()) {
            // No contributing rows: the combined row must be identically zero.
            for (const fr& value : w) {
                if (value != fr::zero()) {
                    return false;
                }
            }
            return true;
        }
        // batch_mul may clobber its scalar buffer; work on copies.
        std::vector<fr> scalars_copy = scalars;
        const Commitment expected =
            Commitment::batch_mul(std::span<const Commitment>(points), std::span<fr>(scalars_copy));
        std::vector<fr> w_copy = w;
        const Commitment computed =
            Commitment::batch_mul(std::span<const Commitment>(ck.generators), std::span<fr>(w_copy));
        return computed == expected;
    }
};

} // namespace bb::hyrax
