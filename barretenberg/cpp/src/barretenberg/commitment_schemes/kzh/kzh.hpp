#pragma once

#include "barretenberg/commitment_schemes/utils/batch_accumulate.hpp"
#include "barretenberg/commitment_schemes/whir/weights.hpp"
#include "barretenberg/common/assert.hpp"
#include "barretenberg/ecc/curves/bn254/bn254.hpp"
#include "barretenberg/crypto/blake3s/blake3s.hpp"
#include "barretenberg/ecc/curves/bn254/pairing.hpp"
#include "barretenberg/polynomials/polynomial.hpp"

#include <cstring>
#include <map>
#include <memory>
#include <span>
#include <string>
#include <vector>

namespace bb::kzh {

using Curve = curve::BN254;
using Commitment = Curve::AffineElement;
using GroupElement = Curve::Element;
using G2Affine = Curve::G2AffineElement;

/** @brief Reference to one polynomial of one committed group (duck-type shared with the suite). */
struct KzhColumnRef {
    size_t group;
    size_t column;
};

/** @brief KZH2 parameters: the matrix split (columns = low variables). */
struct KzhConfig {
    size_t num_variables;
    size_t log_num_cols;

    size_t num_rows() const { return size_t(1) << (num_variables - log_num_cols); }
    size_t num_cols() const { return size_t(1) << log_num_cols; }

    static KzhConfig create(size_t num_variables, size_t /*security_bits*/ = 100, size_t /*log_inv_rate*/ = 0)
    {
        BB_ASSERT_GT(num_variables, size_t(1));
        return { num_variables, (num_variables + 1) / 2 };
    }
};

/**
 * @brief Prover-side commitment to a group: dense arrays, one KZH commitment per column, and the
 * per-column row commitments D_i (the opening's auxiliary data, computed at commit time).
 */
struct KzhGroupData {
    std::vector<std::vector<fr>> coefficients;
    std::vector<Commitment> commitments;
    std::vector<std::vector<Commitment>> row_commitments; // [column][row]

    size_t num_columns() const { return coefficients.size(); }
};

/**
 * @brief KZH2 structured setup (alinush.github.io/kzh; the k=2 scheme of the KZH/KZH-Fold line):
 * G1 generators A_j for one matrix row, row trapdoors tau_i with V_i = tau_i V2 in G2. The
 * commitment C = sum_{i,j} f_{i,j} (tau_i A_j) binds through the pairing identity
 * e(C, V2) = prod_i e(D_i, V_i) with D_i = sum_j f_{i,j} A_j.
 *
 * @warning Test-only setup: this key SAMPLES the trapdoors tau (a real deployment runs a ceremony
 * and discards them) and uses them to commit via the two-level MSM shortcut C = sum_i tau_i D_i,
 * which computes the identical group elements a production prover derives from the published
 * H_{i,j} = tau_i A_j grid at the same asymptotic cost.
 */
class KzhCommitmentKey {
  public:
    explicit KzhCommitmentKey(const KzhConfig& config)
        : config(config)
        , generators(
              Curve::Group::derive_generators(std::vector<uint8_t>{ 'b', 'b', '_', 'k', 'z', 'h' }, config.num_cols()))
    {
        const size_t num_rows = config.num_rows();
        row_trapdoors.resize(num_rows);
        v2_powers.resize(num_rows);
        for (size_t i = 0; i < num_rows; ++i) {
            // Deterministic pseudo-trapdoors so prover and verifier derive the same test-only
            // setup; a real deployment samples these in a ceremony and publishes only the group
            // elements.
            std::vector<uint8_t> seed = { 'k', 'z', 'h', '_', 't', 'a', 'u' };
            for (size_t byte = 0; byte < 8; ++byte) {
                seed.push_back(static_cast<uint8_t>((i >> (8 * byte)) & 0xff));
            }
            const auto digest = blake3::blake3s(seed);
            uint256_t value(0);
            std::memcpy(value.data, digest.data(), 32);
            row_trapdoors[i] = fr(value);
            v2_powers[i] = G2Affine(bb::g2::element(bb::g2::affine_one) * row_trapdoors[i]);
        }
    }

    KzhGroupData commit_group(std::vector<std::vector<fr>> payload_columns,
                              const std::vector<bool>& /*to_be_shifted*/ = {}) const
    {
        const size_t n = size_t(1) << config.num_variables;
        const size_t num_rows = config.num_rows();
        const size_t num_cols = config.num_cols();
        KzhGroupData data;
        for (auto& column : payload_columns) {
            BB_ASSERT_LTE(column.size(), n, "polynomial too large for the configured size");
            column.resize(n, fr::zero());
            std::vector<Commitment> rows(num_rows);
            for (size_t r = 0; r < num_rows; ++r) {
                std::vector<fr> row(column.begin() + static_cast<std::ptrdiff_t>(r * num_cols),
                                    column.begin() + static_cast<std::ptrdiff_t>((r + 1) * num_cols));
                rows[r] = Commitment::batch_mul(std::span<const Commitment>(generators), std::span<fr>(row));
            }
            std::vector<fr> trapdoor_copy = row_trapdoors;
            data.commitments.push_back(
                Commitment::batch_mul(std::span<const Commitment>(rows), std::span<fr>(trapdoor_copy)));
            data.row_commitments.push_back(std::move(rows));
            data.coefficients.push_back(std::move(column));
        }
        return data;
    }

    KzhGroupData commit_group(std::span<const Polynomial<fr>* const> polynomials,
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

    KzhConfig config;
    std::vector<Commitment> generators; // A_j
    std::vector<fr> row_trapdoors;      // tau_i (test-only)
    std::vector<G2Affine> v2_powers;    // V_i = tau_i V2
};

namespace detail {

inline std::string kzh_label(const std::string& name, size_t i)
{
    return "KZH:" + name + "_" + std::to_string(i);
}

} // namespace detail

/**
 * @brief KZH2 opening prover for the batched Honk claim set: chain A carries the unshifted
 * rho-combination, chain B the (unshifted!) to-be-shifted combination whose shift lives in the
 * verifier-side rank-2 query decomposition, as in the Hyrax/Ligero backends. Each chain sends its
 * combined row commitments D_i and combined row vectors; the verifier binds the D_i to the
 * homomorphically-combined commitment with one multipairing per chain.
 */
class KzhProver {
  public:
    struct Claims {
        std::vector<const KzhGroupData*> groups;
        std::vector<KzhColumnRef> unshifted;
        std::vector<fr> unshifted_evaluations;
        std::vector<KzhColumnRef> to_be_shifted;
        std::vector<fr> shifted_evaluations;
        bool send_roots = true;
    };

    template <typename Transcript>
    static void prove(const KzhCommitmentKey& ck,
                      const Claims& claims,
                      std::span<const fr> u,
                      const std::shared_ptr<Transcript>& transcript)
    {
        const KzhConfig& config = ck.config;
        const size_t n = size_t(1) << config.num_variables;
        const size_t num_rows = config.num_rows();
        const size_t num_cols = config.num_cols();
        BB_ASSERT_EQ(u.size(), config.num_variables, "opening point size mismatch");

        if (claims.send_roots) {
            for (size_t g = 0; g < claims.groups.size(); ++g) {
                for (size_t c = 0; c < claims.groups[g]->num_columns(); ++c) {
                    transcript->send_to_verifier(detail::kzh_label("root", g) + "_" + std::to_string(c),
                                                 claims.groups[g]->commitments[c]);
                }
            }
        }
        const fr rho = transcript->template get_challenge<fr>("KZH:rho");

        // Combined arrays for the two chains.
        std::vector<std::vector<fr>> chains;
        chains.emplace_back(n, fr::zero());
        std::vector<pcs_utils::ScaledTerm> terms_a;
        std::vector<pcs_utils::ScaledTerm> terms_b;
        fr rho_power = fr::one();
        for (size_t i = 0; i < claims.unshifted.size(); ++i) {
            const auto& column = claims.groups[claims.unshifted[i].group]->coefficients[claims.unshifted[i].column];
            terms_a.push_back({ column.data(), rho_power, false });
            rho_power *= rho;
        }
        if (!claims.to_be_shifted.empty()) {
            chains.emplace_back(n, fr::zero());
            for (size_t l = 0; l < claims.to_be_shifted.size(); ++l) {
                const auto& column =
                    claims.groups[claims.to_be_shifted[l].group]->coefficients[claims.to_be_shifted[l].column];
                terms_b.push_back({ column.data(), rho_power, false });
                rho_power *= rho;
            }
        }
        pcs_utils::accumulate_scaled(chains[0], terms_a);
        if (chains.size() > 1) {
            pcs_utils::accumulate_scaled(chains[1], terms_b);
        }

        const std::vector<fr> b = whir::eq_tensor(u.subspan(config.log_num_cols));
        for (size_t chain = 0; chain < chains.size(); ++chain) {
            // Combined row commitments D_i of the chain array (R MSMs of size C).
            for (size_t r = 0; r < num_rows; ++r) {
                std::vector<fr> row(chains[chain].begin() + static_cast<std::ptrdiff_t>(r * num_cols),
                                    chains[chain].begin() + static_cast<std::ptrdiff_t>((r + 1) * num_cols));
                const Commitment d_r =
                    Commitment::batch_mul(std::span<const Commitment>(ck.generators), std::span<fr>(row));
                transcript->send_to_verifier(detail::kzh_label("D", chain) + "_" + std::to_string(r), d_r);
            }
            // Combined rows: w_u (chain 0, weights b) or w_s, w_s2 (chain 1, weights b and shr(b)).
            std::vector<fr> w(num_cols, fr::zero());
            std::vector<fr> w2(num_cols, fr::zero());
            for (size_t r = 0; r < num_rows; ++r) {
                const fr* row = chains[chain].data() + r * num_cols;
                for (size_t c = 0; c < num_cols; ++c) {
                    w[c] += b[r] * row[c];
                    if (chain == 1 && r > 0) {
                        w2[c] += b[r - 1] * row[c];
                    }
                }
            }
            for (size_t c = 0; c < num_cols; ++c) {
                transcript->send_to_verifier(detail::kzh_label("w", chain) + "_" + std::to_string(c), w[c]);
            }
            if (chain == 1) {
                for (size_t c = 0; c < num_cols; ++c) {
                    transcript->send_to_verifier(detail::kzh_label("w2", chain) + "_" + std::to_string(c), w2[c]);
                }
            }
        }
    }
};

/** @brief The mirror of `KzhProver`: one multipairing and three MSM identities per opening. */
class KzhVerifier {
  public:
    struct Claims {
        std::vector<size_t> group_num_columns;
        std::vector<std::vector<Commitment>> group_commitments; // when non-empty, transcript-bound
        std::vector<KzhColumnRef> unshifted;
        std::vector<fr> unshifted_evaluations;
        std::vector<KzhColumnRef> to_be_shifted;
        std::vector<fr> shifted_evaluations;
    };

    template <typename Transcript>
    static bool verify(const KzhConfig& config,
                       const Claims& claims,
                       std::span<const fr> u,
                       const std::shared_ptr<Transcript>& transcript,
                       const KzhCommitmentKey& ck)
    {
        const size_t num_rows = config.num_rows();
        const size_t num_cols = config.num_cols();
        BB_ASSERT_EQ(u.size(), config.num_variables, "opening point size mismatch");

        std::vector<std::vector<Commitment>> group_commitments = claims.group_commitments;
        if (group_commitments.empty()) {
            for (size_t g = 0; g < claims.group_num_columns.size(); ++g) {
                std::vector<Commitment> commitments;
                for (size_t c = 0; c < claims.group_num_columns[g]; ++c) {
                    commitments.push_back(transcript->template receive_from_prover<Commitment>(
                        detail::kzh_label("root", g) + "_" + std::to_string(c)));
                }
                group_commitments.push_back(std::move(commitments));
            }
        }
        const fr rho = transcript->template get_challenge<fr>("KZH:rho");

        struct ChainView {
            GroupElement commitment = GroupElement::infinity();
            fr claimed_evaluation = fr::zero();
        };
        std::vector<ChainView> chains(claims.to_be_shifted.empty() ? 1 : 2);
        fr rho_power = fr::one();
        for (size_t i = 0; i < claims.unshifted.size(); ++i) {
            chains[0].commitment +=
                GroupElement(group_commitments[claims.unshifted[i].group][claims.unshifted[i].column]) * rho_power;
            chains[0].claimed_evaluation += rho_power * claims.unshifted_evaluations[i];
            rho_power *= rho;
        }
        if (chains.size() > 1) {
            for (size_t l = 0; l < claims.to_be_shifted.size(); ++l) {
                chains[1].commitment +=
                    GroupElement(group_commitments[claims.to_be_shifted[l].group][claims.to_be_shifted[l].column]) *
                    rho_power;
                chains[1].claimed_evaluation += rho_power * claims.shifted_evaluations[l];
                rho_power *= rho;
            }
        }

        const std::vector<fr> a = whir::eq_tensor(u.subspan(0, config.log_num_cols));
        const std::vector<fr> b = whir::eq_tensor(u.subspan(config.log_num_cols));

        for (size_t chain = 0; chain < chains.size(); ++chain) {
            std::vector<Commitment> row_commitments(num_rows);
            for (size_t r = 0; r < num_rows; ++r) {
                row_commitments[r] = transcript->template receive_from_prover<Commitment>(
                    detail::kzh_label("D", chain) + "_" + std::to_string(r));
            }
            std::vector<fr> w(num_cols);
            for (size_t c = 0; c < num_cols; ++c) {
                w[c] = transcript->template receive_from_prover<fr>(detail::kzh_label("w", chain) + "_" +
                                                                    std::to_string(c));
            }
            std::vector<fr> w2(num_cols, fr::zero());
            if (chain == 1) {
                for (size_t c = 0; c < num_cols; ++c) {
                    w2[c] = transcript->template receive_from_prover<fr>(detail::kzh_label("w2", chain) + "_" +
                                                                         std::to_string(c));
                }
            }

            // Check 1 (binding): e(-C, V2) * prod_i e(D_i, V_i) == 1.
            std::vector<Commitment> pairing_g1;
            std::vector<G2Affine> pairing_g2;
            pairing_g1.push_back(Commitment(-chains[chain].commitment));
            pairing_g2.push_back(G2Affine(bb::g2::affine_one));
            for (size_t r = 0; r < num_rows; ++r) {
                pairing_g1.push_back(row_commitments[r]);
                pairing_g2.push_back(ck.v2_powers[r]);
            }
            if (pairing::reduced_ate_pairing_batch(pairing_g1.data(), pairing_g2.data(), pairing_g1.size()) !=
                fq12::one()) {
                return false;
            }

            // Check 2 (row combination): Commit_A(w) == sum_i weight_i D_i.
            if (!check_combined_row(ck, w, row_commitments, b, /*shift_weights=*/false)) {
                return false;
            }
            if (chain == 1 && !check_combined_row(ck, w2, row_commitments, b, /*shift_weights=*/true)) {
                return false;
            }

            // Check 3 (claim): the rank-2 evaluation identity.
            fr total = fr::zero();
            if (chain == 0) {
                for (size_t c = 0; c < num_cols; ++c) {
                    total += w[c] * a[c];
                }
            } else {
                for (size_t c = 1; c < num_cols; ++c) {
                    total += w[c] * a[c - 1];
                }
                total += a[num_cols - 1] * w2[0];
            }
            if (total != chains[chain].claimed_evaluation) {
                return false;
            }
        }
        return true;
    }

  private:
    static bool check_combined_row(const KzhCommitmentKey& ck,
                                   const std::vector<fr>& w,
                                   const std::vector<Commitment>& row_commitments,
                                   const std::vector<fr>& b,
                                   bool shift_weights)
    {
        std::vector<fr> weights(row_commitments.size(), fr::zero());
        for (size_t r = 0; r < row_commitments.size(); ++r) {
            if (shift_weights) {
                if (r > 0) {
                    weights[r] = b[r - 1];
                }
            } else {
                weights[r] = b[r];
            }
        }
        const Commitment expected =
            Commitment::batch_mul(std::span<const Commitment>(row_commitments), std::span<fr>(weights));
        std::vector<fr> w_copy = w;
        const Commitment computed =
            Commitment::batch_mul(std::span<const Commitment>(ck.generators), std::span<fr>(w_copy));
        return computed == expected;
    }
};

} // namespace bb::kzh
