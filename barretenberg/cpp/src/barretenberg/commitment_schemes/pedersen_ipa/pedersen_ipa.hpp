#pragma once

#include "barretenberg/commitment_schemes/ipa/ipa.hpp"
#include "barretenberg/common/assert.hpp"
#include "barretenberg/ecc/curves/bn254/bn254.hpp"
#include "barretenberg/ecc/scalar_multiplication/scalar_multiplication.hpp"
#include "barretenberg/polynomials/eq_polynomial.hpp"
#include "barretenberg/polynomials/polynomial.hpp"
#include "barretenberg/polynomials/shifted_eq_polynomial.hpp"

#include <map>
#include <memory>
#include <mutex>
#include <span>
#include <string>
#include <vector>

namespace bb::pedersen_ipa {

using Curve = curve::BN254;
using Commitment = Curve::AffineElement;
using GroupElement = Curve::Element;

/** @brief Reference to one polynomial of one committed group (duck-type shared with the suite). */
struct IpaColumnRef {
    size_t group;
    size_t column;
};

/** @brief Parameters: only the size; the scheme is transparent and unstructured. */
struct IpaConfig {
    size_t num_variables;

    static IpaConfig create(size_t num_variables, size_t /*security_bits*/ = 100, size_t /*log_inv_rate*/ = 0)
    {
        return { num_variables };
    }
};

/** @brief Prover-side commitment to a group: dense arrays and one Pedersen commitment per column. */
struct IpaGroupData {
    std::vector<std::vector<fr>> coefficients;
    std::vector<Commitment> commitments;

    size_t num_columns() const { return coefficients.size(); }
};

/**
 * @brief Transparent commitment key: hash-derived independent generators (the powers-of-tau SRS
 * cannot serve as Pedersen/IPA generators — its points are DL-related). The auxiliary generator of
 * the inner-product relation is `Commitment::one()`, DL-independent of the derived vector. The
 * derived vector is a public parameter shared by every key of the same size, so it is memoized:
 * constructing a key (e.g. the verifier recreating one per proof) does not re-run the O(n)
 * hash-to-curve derivation.
 */
class IpaCommitmentKey {
  public:
    explicit IpaCommitmentKey(const IpaConfig& config)
        : config(config)
        , generators_ptr(shared_generators(config.num_variables))
        , generators(*generators_ptr)
    {}

  private:
    static std::shared_ptr<const std::vector<Commitment>> shared_generators(size_t num_variables)
    {
        static std::mutex mutex;
        static std::map<size_t, std::shared_ptr<const std::vector<Commitment>>> cache;
        const std::lock_guard<std::mutex> lock(mutex);
        auto& entry = cache[num_variables];
        if (!entry) {
            entry = std::make_shared<const std::vector<Commitment>>(Curve::Group::derive_generators(
                std::vector<uint8_t>{ 'b', 'b', '_', 'i', 'p', 'a' }, size_t(1) << num_variables));
        }
        return entry;
    }

  public:
    Commitment commit(std::span<const fr> coefficients) const
    {
        std::vector<fr> scalars(coefficients.begin(), coefficients.end());
        scalars.resize(size_t(1) << config.num_variables, fr::zero());
        return Commitment::batch_mul(std::span<const Commitment>(generators), std::span<fr>(scalars));
    }

    IpaGroupData commit_group(std::vector<std::vector<fr>> payload_columns,
                              const std::vector<bool>& /*to_be_shifted*/ = {}) const
    {
        const size_t n = size_t(1) << config.num_variables;
        IpaGroupData data;
        for (auto& column : payload_columns) {
            BB_ASSERT_LTE(column.size(), n, "polynomial too large for the configured size");
            data.commitments.push_back(commit(column));
            column.resize(n, fr::zero());
            data.coefficients.push_back(std::move(column));
        }
        return data;
    }

    IpaGroupData commit_group(std::span<const Polynomial<fr>* const> polynomials,
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

    IpaConfig config;

  private:
    std::shared_ptr<const std::vector<Commitment>> generators_ptr;

  public:
    const std::vector<Commitment>& generators;
};

namespace detail {

inline std::string root_label(size_t group, size_t column)
{
    return "IPA:root_" + std::to_string(group) + "_" + std::to_string(column);
}

/**
 * @brief Instantiate `fn` at the compile-time log-size matching the runtime `log_n`.
 * @details The shared IPA core (`bb::IPA`) is compile-time sized; the transparent suite configures
 * sizes at runtime, so prover and verifier dispatch through this switch.
 */
template <typename Fn> decltype(auto) with_log_n(size_t log_n, Fn&& fn)
{
    switch (log_n) {
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define BB_PEDERSEN_IPA_CASE(L)                                                                                        \
    case (L):                                                                                                          \
        return fn.template operator()<(L)>();
        BB_PEDERSEN_IPA_CASE(1)
        BB_PEDERSEN_IPA_CASE(2)
        BB_PEDERSEN_IPA_CASE(3)
        BB_PEDERSEN_IPA_CASE(4)
        BB_PEDERSEN_IPA_CASE(5)
        BB_PEDERSEN_IPA_CASE(6)
        BB_PEDERSEN_IPA_CASE(7)
        BB_PEDERSEN_IPA_CASE(8)
        BB_PEDERSEN_IPA_CASE(9)
        BB_PEDERSEN_IPA_CASE(10)
        BB_PEDERSEN_IPA_CASE(11)
        BB_PEDERSEN_IPA_CASE(12)
        BB_PEDERSEN_IPA_CASE(13)
        BB_PEDERSEN_IPA_CASE(14)
        BB_PEDERSEN_IPA_CASE(15)
        BB_PEDERSEN_IPA_CASE(16)
        BB_PEDERSEN_IPA_CASE(17)
        BB_PEDERSEN_IPA_CASE(18)
        BB_PEDERSEN_IPA_CASE(19)
        BB_PEDERSEN_IPA_CASE(20)
        BB_PEDERSEN_IPA_CASE(21)
        BB_PEDERSEN_IPA_CASE(22)
        BB_PEDERSEN_IPA_CASE(23)
        BB_PEDERSEN_IPA_CASE(24)
#undef BB_PEDERSEN_IPA_CASE
    default:
        throw_or_abort("pedersen_ipa: unsupported num_variables " + std::to_string(log_n));
    }
}

} // namespace detail

/**
 * @brief Multilinear Bulletproofs inner-product argument over hash-derived BN254 generators
 * (non-hiding): the {transparent, small-proof, linear-verifier} corner of the PCS suite.
 *
 * The claims are reduced to a single run of the shared IPA core (`IPA<Curve,
 * log_n>::compute_inner_product_proof_internal`), inheriting its prover optimizations: 127-bit
 * short round challenges, the rescaled generator fold G' = u·G_lo + G_hi (short raw challenge
 * instead of the full-width inverse, unscaled once at the end), and the fused two-round
 * `batch_two_round_fold` (see `ipa/ELEMENT_IMPL_FOLD.md`).
 *
 * Reduction (the TripleIPA pattern, with the pow tensor omitted):
 *   1. `IPA:rho` batches the unshifted claims into F = Σ ρⁱ·fᵢ with b-vector eq(u), and the
 *      to-be-shifted claims into F' = Σ ρˡ·gˡ with b-vector b_sh (b_sh[0] = 0, b_sh[i] = eq[i-1],
 *      so ⟨F', b_sh⟩ is the batched shifted evaluation — no shifted commitments needed).
 *   2. The prover sends the cross-sum `IPA:cross_F_shift` = ⟨F, b_sh⟩ + ⟨F', eq⟩; the challenges
 *      `IPA:zeta_F`, `IPA:zeta_shift` merge the two claims into one:
 *        witness  A = ζ_F·F + ζ_sh·F',    tensor  b = ζ_F·eq + ζ_sh·b_sh,
 *        ⟨A, b⟩ = ζ_F²·v_F + ζ_sh²·v_sh + ζ_F·ζ_sh·cross.
 *   3. One IPA core run opens ⟨A, b⟩ over the derived generators.
 *
 * The verifier's b₀ contraction is O(log n) via the closed forms
 * `ShiftedEqPolynomial::evaluate_eq_folded` / `evaluate_folded`; its linear work is the single
 * size-n MSM certifying the prover-claimed G₀ (the Halo amortization point — deferrable in
 * principle, checked inline here).
 */
class IpaProver {
  public:
    struct Claims {
        std::vector<const IpaGroupData*> groups;
        std::vector<IpaColumnRef> unshifted;
        std::vector<fr> unshifted_evaluations;
        std::vector<IpaColumnRef> to_be_shifted;
        std::vector<fr> shifted_evaluations;
        bool send_roots = true;
    };

    template <typename Transcript>
    static void prove(const IpaCommitmentKey& ck,
                      const Claims& claims,
                      std::span<const fr> u,
                      const std::shared_ptr<Transcript>& transcript)
    {
        BB_ASSERT_EQ(u.size(), ck.config.num_variables, "opening point size mismatch");
        detail::with_log_n(ck.config.num_variables,
                           [&]<size_t LOG_N>() { prove_impl<LOG_N>(ck, claims, u, transcript); });
    }

  private:
    template <size_t LOG_N, typename Transcript>
    static void prove_impl(const IpaCommitmentKey& ck,
                           const Claims& claims,
                           std::span<const fr> u,
                           const std::shared_ptr<Transcript>& transcript)
    {
        using Ipa = IPA<Curve, LOG_N>;
        using ShiftedEq = ShiftedEqPolynomial<Curve, LOG_N>;
        constexpr size_t n = size_t(1) << LOG_N;

        if (claims.send_roots) {
            for (size_t g = 0; g < claims.groups.size(); ++g) {
                for (size_t c = 0; c < claims.groups[g]->num_columns(); ++c) {
                    transcript->send_to_verifier(detail::root_label(g, c), claims.groups[g]->commitments[c]);
                }
            }
        }
        const fr rho = transcript->template get_challenge<fr>("IPA:rho");

        // rho-batch the unshifted claims (F) and, continuing the same rho powers, the to-be-shifted
        // claims (F'). The shift lives entirely in F''s b-vector.
        fr rho_power = fr::one();
        const auto batch_chain = [&](std::span<const IpaColumnRef> refs, std::span<const fr> evaluations) {
            std::vector<PolynomialSpan<const fr>> sources;
            std::vector<fr> weights;
            GroupElement commitment = GroupElement::infinity();
            fr evaluation = fr::zero();
            for (size_t i = 0; i < refs.size(); ++i) {
                const auto& group = *claims.groups[refs[i].group];
                sources.push_back({ 0, std::span<const fr>(group.coefficients[refs[i].column]) });
                weights.push_back(rho_power);
                commitment += GroupElement(group.commitments[refs[i].column]) * rho_power;
                evaluation += rho_power * evaluations[i];
                rho_power *= rho;
            }
            Polynomial<fr> witness(n);
            add_scaled_batch(witness, std::span<const PolynomialSpan<const fr>>(sources), std::span<const fr>(weights));
            return std::make_tuple(std::move(witness), commitment, evaluation);
        };
        auto [unshifted_witness, unshifted_commitment, unshifted_evaluation] =
            batch_chain(claims.unshifted, claims.unshifted_evaluations);

        Polynomial<fr> eq = ProverEqPolynomial<fr>::construct(u, LOG_N);

        Polynomial<fr> witness;
        Polynomial<fr> b_vec;
        GroupElement combined_commitment = GroupElement::infinity();
        fr evaluation = fr::zero();
        if (claims.to_be_shifted.empty()) {
            witness = std::move(unshifted_witness);
            b_vec = std::move(eq);
            combined_commitment = unshifted_commitment;
            evaluation = unshifted_evaluation;
        } else {
            auto [shifted_witness, shifted_commitment, shifted_evaluation] =
                batch_chain(claims.to_be_shifted, claims.shifted_evaluations);

            // cross = ⟨F, b_sh⟩ + ⟨F', eq⟩; sent before the zetas that weight it.
            const fr cross = ShiftedEq::evaluate_from_eq(eq, unshifted_witness) + shifted_witness.evaluate_mle(u);
            transcript->send_to_verifier("IPA:cross_F_shift", cross);
            const auto zeta =
                transcript->template get_challenges<fr>(std::array<std::string, 2>{ "IPA:zeta_F", "IPA:zeta_shift" });

            witness = Polynomial<fr>(n);
            witness.add_scaled(unshifted_witness, zeta[0]);
            witness.add_scaled(shifted_witness, zeta[1]);
            b_vec = Polynomial<fr>(n);
            b_vec.add_scaled(eq, zeta[0]);
            ShiftedEq::add_scaled(b_vec, eq, zeta[1]);
            combined_commitment = unshifted_commitment * zeta[0] + shifted_commitment * zeta[1];
            evaluation =
                zeta[0].sqr() * unshifted_evaluation + zeta[1].sqr() * shifted_evaluation + zeta[0] * zeta[1] * cross;
        }

        add_combined_claim_to_hash_buffer(transcript, u, Commitment(combined_commitment), evaluation);
        Ipa::compute_inner_product_proof_internal(
            std::span<const Commitment>(ck.generators), witness, std::move(b_vec), transcript);
    }

    template <typename Transcript>
    static void add_combined_claim_to_hash_buffer(const std::shared_ptr<Transcript>& transcript,
                                                  std::span<const fr> u,
                                                  const Commitment& commitment,
                                                  const fr& evaluation)
    {
        for (size_t coordinate_idx = 0; coordinate_idx < u.size(); ++coordinate_idx) {
            transcript->add_to_hash_buffer("IPA:u_" + std::to_string(coordinate_idx), u[coordinate_idx]);
        }
        transcript->add_to_hash_buffer("IPA:combined_commitment", commitment);
        transcript->add_to_hash_buffer("IPA:combined_evaluation", evaluation);
    }

    friend class IpaVerifier;
};

/** @brief The mirror of `IpaProver`; the SRS-style MSM certifying G₀ makes verification O(n). */
class IpaVerifier {
  public:
    struct Claims {
        std::vector<size_t> group_num_columns;
        std::vector<std::vector<Commitment>> group_commitments; // when non-empty, transcript-bound
        std::vector<IpaColumnRef> unshifted;
        std::vector<fr> unshifted_evaluations;
        std::vector<IpaColumnRef> to_be_shifted;
        std::vector<fr> shifted_evaluations;
    };

    template <typename Transcript>
    static bool verify(const IpaConfig& config,
                       const Claims& claims,
                       std::span<const fr> u,
                       const std::shared_ptr<Transcript>& transcript,
                       const IpaCommitmentKey& ck)
    {
        BB_ASSERT_EQ(u.size(), config.num_variables, "opening point size mismatch");
        return detail::with_log_n(config.num_variables,
                                  [&]<size_t LOG_N>() { return verify_impl<LOG_N>(claims, u, transcript, ck); });
    }

  private:
    template <size_t LOG_N, typename Transcript>
    static bool verify_impl(const Claims& claims,
                            std::span<const fr> u,
                            const std::shared_ptr<Transcript>& transcript,
                            const IpaCommitmentKey& ck)
    {
        using Ipa = IPA<Curve, LOG_N>;
        using ShiftedEq = ShiftedEqPolynomial<Curve, LOG_N>;
        constexpr size_t n = size_t(1) << LOG_N;

        std::vector<std::vector<Commitment>> group_commitments = claims.group_commitments;
        if (group_commitments.empty()) {
            for (size_t g = 0; g < claims.group_num_columns.size(); ++g) {
                std::vector<Commitment> commitments;
                for (size_t c = 0; c < claims.group_num_columns[g]; ++c) {
                    commitments.push_back(
                        transcript->template receive_from_prover<Commitment>(detail::root_label(g, c)));
                }
                group_commitments.push_back(std::move(commitments));
            }
        }
        const fr rho = transcript->template get_challenge<fr>("IPA:rho");

        fr rho_power = fr::one();
        const auto batch_chain = [&](std::span<const IpaColumnRef> refs, std::span<const fr> evaluations) {
            GroupElement commitment = GroupElement::infinity();
            fr evaluation = fr::zero();
            for (size_t i = 0; i < refs.size(); ++i) {
                commitment += GroupElement(group_commitments[refs[i].group][refs[i].column]) * rho_power;
                evaluation += rho_power * evaluations[i];
                rho_power *= rho;
            }
            return std::make_pair(commitment, evaluation);
        };
        auto [unshifted_commitment, unshifted_evaluation] = batch_chain(claims.unshifted, claims.unshifted_evaluations);

        const bool has_shifted = !claims.to_be_shifted.empty();
        GroupElement combined_commitment = unshifted_commitment;
        fr evaluation = unshifted_evaluation;
        fr zeta_unshifted = fr::one();
        fr zeta_shifted = fr::zero();
        if (has_shifted) {
            auto [shifted_commitment, shifted_evaluation] =
                batch_chain(claims.to_be_shifted, claims.shifted_evaluations);
            const fr cross = transcript->template receive_from_prover<fr>("IPA:cross_F_shift");
            const auto zeta =
                transcript->template get_challenges<fr>(std::array<std::string, 2>{ "IPA:zeta_F", "IPA:zeta_shift" });
            zeta_unshifted = zeta[0];
            zeta_shifted = zeta[1];
            combined_commitment = unshifted_commitment * zeta_unshifted + shifted_commitment * zeta_shifted;
            evaluation = zeta_unshifted.sqr() * unshifted_evaluation + zeta_shifted.sqr() * shifted_evaluation +
                         zeta_unshifted * zeta_shifted * cross;
        }

        const Commitment combined_commitment_affine(combined_commitment);
        IpaProver::add_combined_claim_to_hash_buffer(transcript, u, combined_commitment_affine, evaluation);

        // b₀ = ζ_F·⟨eq, s⟩ + ζ_sh·⟨b_sh, s⟩, both contractions O(log n) via the closed forms.
        const auto b_zero_from_round_challenges = [&](std::span<const fr> round_challenges_inv) {
            fr result = zeta_unshifted * ShiftedEq::evaluate_eq_folded(u, round_challenges_inv);
            if (has_shifted) {
                result += zeta_shifted * ShiftedEq::evaluate_folded(u, round_challenges_inv);
            }
            return result;
        };
        const auto data = Ipa::read_inner_product_transcript_data(
            combined_commitment_affine, evaluation, b_zero_from_round_challenges, transcript);

        // Certify the prover-claimed G₀ = ⟨s, G⟩ with the single size-n MSM over the derived
        // generators, then check the IPA group relation.
        BB_ASSERT_EQ(ck.generators.size(), n, "generator count mismatch");
        Commitment G_zero;
        G_zero = scalar_multiplication::pippenger_unsafe<Curve>(data.s_vec, { ck.generators.data(), n });
        if (G_zero != data.G_zero_from_prover) {
            return false;
        }
        const Commitment aux_generator = Commitment::one() * data.gen_challenge;
        const GroupElement right_hand_side = G_zero * data.a_zero + aux_generator * (data.a_zero * data.b_zero);
        return data.C_zero.normalize() == right_hand_side.normalize();
    }
};

} // namespace bb::pedersen_ipa
