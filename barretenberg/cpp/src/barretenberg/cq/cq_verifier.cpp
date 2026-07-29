#include "cq_verifier.hpp"

#include "barretenberg/common/log.hpp"
#include "barretenberg/ecc/curves/bn254/fq12.hpp"
#include "barretenberg/ecc/curves/bn254/pairing.hpp"

#include <array>

namespace bb::cq {

bool CqVerifier::verify_proof(const Proof& proof) const
{
    auto transcript = std::make_shared<Transcript>(proof);

    const auto table_size = transcript->template receive_from_prover<uint32_t>("CQ:table_size");
    const auto num_lookups = transcript->template receive_from_prover<uint32_t>("CQ:num_lookups");
    const auto num_columns = transcript->template receive_from_prover<uint32_t>("CQ:num_columns");
    if (table_size != vk_.table_size || num_columns != vk_.num_columns()) {
        vinfo("cq: proof was constructed for a different table shape");
        return false;
    }
    const size_t n = num_lookups;
    if (n < 2 || n > table_size || (n & (n - 1)) != 0) {
        vinfo("cq: invalid number of lookups");
        return false;
    }
    for (size_t k = 0; k < num_columns; ++k) {
        const auto table_commitment =
            transcript->template receive_from_prover<g1::affine_element>("CQ:T_" + std::to_string(k));
        if (table_commitment != vk_.table_commitments_g1[k]) {
            vinfo("cq: proof was constructed for a different table");
            return false;
        }
    }

    const auto m_commitment = transcript->template receive_from_prover<g1::affine_element>("CQ:m");
    std::vector<g1::affine_element> f_commitments(num_columns);
    for (size_t k = 0; k < num_columns; ++k) {
        f_commitments[k] = transcript->template receive_from_prover<g1::affine_element>("CQ:f_" + std::to_string(k));
    }
    const fr theta = transcript->template get_challenge<fr>("CQ:theta");
    const fr beta = transcript->template get_challenge<fr>("CQ:beta");

    // Homomorphic column combination: [f] = sum_k theta^k [f^k] and [T]_2 = sum_k theta^k [T^k]_2.
    g1::element f_combined = g1::element(f_commitments[0]);
    g2::element table_combined_g2 = g2::element(vk_.table_commitments_g2[0]);
    fr theta_power = fr::one();
    for (size_t k = 1; k < num_columns; ++k) {
        theta_power *= theta;
        f_combined += g1::element(f_commitments[k]) * theta_power;
        table_combined_g2 += g2::element(vk_.table_commitments_g2[k]) * theta_power;
    }
    const g1::affine_element f_commitment(f_combined);
    const g2::affine_element table_commitment_g2(table_combined_g2);

    const auto a_commitment = transcript->template receive_from_prover<g1::affine_element>("CQ:A");
    const auto q_a_commitment = transcript->template receive_from_prover<g1::affine_element>("CQ:Q_A");
    const auto b_shifted_commitment = transcript->template receive_from_prover<g1::affine_element>("CQ:B_0");
    const auto q_b_commitment = transcript->template receive_from_prover<g1::affine_element>("CQ:Q_B");
    const auto degree_shift_commitment = transcript->template receive_from_prover<g1::affine_element>("CQ:P_B");
    const fr gamma = transcript->template get_challenge<fr>("CQ:gamma");

    const fr b_shifted_at_gamma = transcript->template receive_from_prover<fr>("CQ:b0_gamma");
    const fr f_at_gamma = transcript->template receive_from_prover<fr>("CQ:f_gamma");
    const fr a_zero = transcript->template receive_from_prover<fr>("CQ:a0");
    const fr eta = transcript->template get_challenge<fr>("CQ:eta");

    const auto opening_at_gamma = transcript->template receive_from_prover<g1::affine_element>("CQ:pi_gamma");
    const auto opening_at_zero = transcript->template receive_from_prover<g1::affine_element>("CQ:pi_zero");
    const fr rho = transcript->template get_challenge<fr>("CQ:rho");

    const g2::affine_element g2_one = vk_.g2_identity();

    // Check 1 (table side): A(X)(T(X) + beta) - m(X) = Q_A(X) Z_H(X) at tau, i.e.
    //   e([A], [T]_2) * e(beta*[A] - [m], [1]_2) * e(-[Q_A], [Z_H]_2) == 1.
    {
        const g1::affine_element scaled_a(g1::element(a_commitment) * beta - g1::element(m_commitment));
        const std::array<g1::affine_element, 3> p1{ a_commitment, scaled_a, g1::affine_element(-q_a_commitment) };
        const std::array<g2::affine_element, 3> p2{ table_commitment_g2, g2_one, vk_.vanishing_commitment_g2 };
        if (pairing::reduced_ate_pairing_batch(p1.data(), p2.data(), 3) != fq12::one()) {
            vinfo("cq: table-side identity check failed");
            return false;
        }
    }

    // Check 2 (degree bound): P_B = B_0 * X^{N-n+1}, i.e. e([B_0], [tau^{N-n+1}]_2) * e(-[P_B], [1]_2) == 1.
    // Together with the SRS being truncated at degree N-1 this enforces deg(B_0) <= n-2.
    {
        const std::array<g1::affine_element, 2> p1{ b_shifted_commitment,
                                                    g1::affine_element(-degree_shift_commitment) };
        const std::array<g2::affine_element, 2> p2{ vk_.degree_check_power(n), g2_one };
        if (pairing::reduced_ate_pairing_batch(p1.data(), p2.data(), 2) != fq12::one()) {
            vinfo("cq: degree check failed");
            return false;
        }
    }

    // Checks 3 and 4, batched by rho: KZG openings of (B_0, f, Q_B) at gamma (batched by eta) and of A at zero.
    // The claimed evaluation of Q_B is derived from the witness-side identity B(X)(f(X) + beta) - 1 = Q_B(X) Z_V(X)
    // with B(gamma) = gamma*B_0(gamma) + b_0 and the sum condition b_0 = N*a_0/n linking the two sides.
    {
        const fr b_zero = a_zero * fr(table_size) / fr(n);
        const fr b_at_gamma = gamma * b_shifted_at_gamma + b_zero;
        const fr vanishing_v_at_gamma = gamma.pow(n) - fr::one();
        if (vanishing_v_at_gamma.is_zero()) {
            vinfo("cq: evaluation challenge landed in the witness domain");
            return false;
        }
        const fr q_b_at_gamma = (b_at_gamma * (f_at_gamma + beta) - fr::one()) / vanishing_v_at_gamma;

        const fr eta_sqr = eta.sqr();
        const g1::element batched_commitment =
            g1::element(b_shifted_commitment) + g1::element(f_commitment) * eta + g1::element(q_b_commitment) * eta_sqr;
        const fr batched_evaluation = b_shifted_at_gamma + eta * f_at_gamma + eta_sqr * q_b_at_gamma;

        // e(C - v*[1] + gamma*pi_gamma + rho*([A] - a_0*[1]), [1]_2) * e(-(pi_gamma + rho*pi_zero), [tau]_2) == 1.
        const g1::element p0 = batched_commitment - g1::one * batched_evaluation +
                               g1::element(opening_at_gamma) * gamma +
                               (g1::element(a_commitment) - g1::one * a_zero) * rho;
        const g1::element p1 = -(g1::element(opening_at_gamma) + g1::element(opening_at_zero) * rho);
        const std::array<g1::affine_element, 2> lhs{ g1::affine_element(p0), g1::affine_element(p1) };
        const std::array<g2::affine_element, 2> rhs{ g2_one, vk_.g2_x() };
        if (pairing::reduced_ate_pairing_batch(lhs.data(), rhs.data(), 2) != fq12::one()) {
            vinfo("cq: batched KZG opening check failed");
            return false;
        }
    }

    return true;
}

} // namespace bb::cq
