#pragma once

#include "barretenberg/commitment_schemes/claim.hpp"
#include "barretenberg/commitment_schemes/commitment_key.hpp"
#include "barretenberg/common/bb_bench.hpp"
#include "barretenberg/common/thread.hpp"
#include "barretenberg/ecc/scalar_multiplication/scalar_multiplication.hpp"
#include "barretenberg/numeric/bitop/get_msb.hpp"
#include "barretenberg/numeric/random/engine.hpp"
#include "barretenberg/polynomials/polynomial.hpp"
#include "barretenberg/zcash/honk/pasta_crs.hpp"

#include <span>
#include <vector>

/**
 * @file halo2_ipa.hpp
 * @brief The inner-product-argument opening proof of halo2 (halo2_proofs poly/commitment/{prover,verifier}.rs), as
 * used by Zcash's Orchard prover, over the barretenberg commitment key.
 *
 * @details Proves that the polynomial p committed in P = <p, G> + [blind] W evaluates to v at x:
 *   1. The prover samples a random s with s(x) = 0, sends S = <s, G> + [s_blind] W, and receives xi, z.
 *   2. p' = p + xi s - v (as constant term) has p'(x) = 0 and commitment P' = P - [v] G_0 + [xi] S.
 *   3. k = log2(n) rounds halve p', b = (1, x, x^2, ...) and G; round j sends
 *        L_j = <p'_hi, G_lo> + [z <p'_hi, b_lo>] U + [l_j] W,  R_j = <p'_lo, G_hi> + [z <p'_lo, b_hi>] U + [r_j] W,
 *      receives u_j and folds p' <- p'_lo + u_j^{-1} p'_hi, b <- b_lo + u_j b_hi, G <- G_lo + u_j G_hi.
 *   4. The prover sends c = p'_final and the blind f.
 * The verifier accepts iff
 *   P - [v] G_0 + [xi] S + sum_j ([u_j^{-1}] L_j + [u_j] R_j) - [c] (<s(u), G> + [z b(x, u)] U) - [f] W = 0,
 * where s(u) is the folding vector of the challenges (an MSM of size n, as in halo2).
 * The size n is the runtime power of two given by the claim's polynomial (zero-padded).
 */
namespace bb::zcash {

template <typename Curve> class Halo2IPA {
  public:
    using Fr = typename Curve::ScalarField;
    using Commitment = typename Curve::AffineElement;
    using GroupElement = typename Curve::Element;

    struct Generators {
        std::span<const Commitment> g; // G_0..G_{n-1}
        Commitment w;                  // blinding generator
        Commitment u;                  // inner-product generator
    };

    /**
     * @brief Multi-scalar multiplications sharing one point vector: result i is sum_k scalars[i][k] * points[k].
     * @param handle_edge_cases false only if the points are known to be independent generators
     */
    static std::vector<Commitment> batch_msm(std::span<const Commitment> points,
                                             std::vector<std::vector<Fr>>& scalars,
                                             bool handle_edge_cases)
    {
        std::vector<PolynomialSpan<Fr>> spans;
        for (auto& s : scalars) {
            BB_ASSERT_LTE(s.size(), points.size());
            spans.emplace_back(0, std::span<Fr>(s));
        }
        return scalar_multiplication::MSM<Curve>::batch_multi_scalar_mul(points, spans, handle_edge_cases);
    }

    static GroupElement msm(std::span<const Fr> scalars, std::span<const Commitment> points)
    {
        if (scalars.empty()) {
            return GroupElement::infinity();
        }
        std::vector<std::vector<Fr>> s{ std::vector<Fr>(scalars.begin(), scalars.end()) };
        return GroupElement(batch_msm(points.subspan(0, scalars.size()), s, /*handle_edge_cases=*/true)[0]);
    }

    static Fr inner_product(std::span<const Fr> a, std::span<const Fr> b)
    {
        Fr acc(0);
        for (size_t i = 0; i < a.size(); ++i) {
            acc += a[i] * b[i];
        }
        return acc;
    }

    // halo2 compute_b: prod_j (1 + u_{k-1-j} x^{2^j})
    static Fr compute_b(const Fr& x, std::span<const Fr> u)
    {
        Fr tmp(1);
        Fr cur = x;
        for (auto it = u.rbegin(); it != u.rend(); ++it) {
            tmp *= Fr(1) + *it * cur;
            cur = cur.sqr();
        }
        return tmp;
    }

    // halo2 compute_s: the coefficients of prod_j (1 + u_{k-1-j} X^{2^j}), scaled by init.
    static std::vector<Fr> compute_s(std::span<const Fr> u, const Fr& init)
    {
        std::vector<Fr> v(size_t{ 1 } << u.size(), Fr(0));
        v[0] = init;
        size_t len = 1;
        for (auto it = u.rbegin(); it != u.rend(); ++it) {
            for (size_t i = 0; i < len; ++i) {
                v[len + i] = v[i] * *it;
            }
            len <<= 1;
        }
        return v;
    }

    // Signed 4-bit windows of a scalar: k = sum_i d_i 16^i with d_i in [-8, 8).
    static constexpr size_t NUM_SIGNED_WINDOWS = 65;
    using SignedDigits = std::array<int8_t, NUM_SIGNED_WINDOWS>;
    static SignedDigits signed_digits(const Fr& scalar)
    {
        const uint256_t k(scalar);
        SignedDigits d{};
        int carry = 0;
        for (size_t i = 0; i < NUM_SIGNED_WINDOWS; ++i) {
            int digit = (i < 64 ? static_cast<int>(k.slice(4 * i, (4 * i) + 4).data[0]) : 0) + carry;
            carry = 0;
            if (digit >= 8) {
                digit -= 16;
                carry = 1;
            }
            d[i] = static_cast<int8_t>(digit);
        }
        return d;
    }

    /**
     * @brief sum_t scalars_t * points_t (scalars given by their signed digits) by Straus' method: the terms share the
     * doublings, each adds from its table of [1..8] P_t.
     */
    static GroupElement straus(std::span<const Commitment> points, std::span<const SignedDigits> digits)
    {
        const size_t num = points.size();
        std::vector<std::array<GroupElement, 8>> tables(num);
        for (size_t t = 0; t < num; ++t) {
            auto& table = tables[t];
            table[0] = GroupElement(points[t]);
            table[1] = table[0].dbl();
            for (size_t j = 2; j < 8; ++j) {
                table[j] = table[j - 1] + points[t];
            }
        }
        GroupElement acc = GroupElement::infinity();
        for (size_t w = NUM_SIGNED_WINDOWS; w-- > 0;) {
            if (w + 1 < NUM_SIGNED_WINDOWS) {
                for (size_t i = 0; i < 4; ++i) {
                    acc.self_dbl();
                }
            }
            for (size_t t = 0; t < num; ++t) {
                const int d = digits[t][w];
                if (d > 0) {
                    acc += tables[t][static_cast<size_t>(d - 1)];
                } else if (d < 0) {
                    acc -= tables[t][static_cast<size_t>(-d - 1)];
                }
            }
        }
        return acc;
    }

    static std::vector<Commitment> normalize(std::vector<GroupElement>& elements)
    {
        GroupElement::batch_normalize(elements.data(), elements.size());
        std::vector<Commitment> out(elements.size());
        parallel_for_range(elements.size(), [&](size_t start, size_t end) {
            for (size_t i = start; i < end; ++i) {
                out[i] = Commitment(elements[i].x, elements[i].y);
            }
        });
        return out;
    }

    /**
     * @brief Prove p(x) = v for the claim's polynomial, where P = <p, G> + [p_blind] W.
     */
    template <typename Transcript>
    static void prove(const Generators& gens,
                      const ProverOpeningClaim<Curve>& claim,
                      const Fr& p_blind,
                      const std::shared_ptr<Transcript>& transcript,
                      numeric::RNG& rng = numeric::get_randomness())
    {
        const size_t n = gens.g.size();
        BB_ASSERT((n & (n - 1)) == 0, "IPA size must be a power of two");
        const size_t k = numeric::get_msb(n);
        const Fr& x = claim.opening_pair.challenge;
        const Fr& v = claim.opening_pair.evaluation;

        std::vector<Fr> p(n, Fr(0));
        const auto& poly = claim.polynomial;
        BB_ASSERT_LTE(poly.end_index(), n);
        for (size_t i = poly.start_index(); i < poly.end_index(); ++i) {
            p[i] = poly[i];
        }

        // Random s with s(x) = 0.
        std::vector<Fr> s(n);
        for (auto& c : s) {
            c = Fr::random_element(&rng);
        }
        Fr s_at_x(0);
        for (size_t i = n; i-- > 0;) {
            s_at_x = s_at_x * x + s[i];
        }
        s[0] -= s_at_x;
        const Fr s_blind = Fr::random_element(&rng);
        BB_BENCH_NAME("halo2_ipa/prove");
        const Commitment s_commitment = GroupElement(msm(s, gens.g)) + GroupElement(gens.w) * s_blind;
        transcript->send_to_verifier("IPA:S", s_commitment);
        const Fr xi = transcript->template get_challenge<Fr>("IPA:xi");
        const Fr z = transcript->template get_challenge<Fr>("IPA:z");

        // p' = p + xi s, with constant term reduced by v so that p'(x) = 0.
        for (size_t i = 0; i < n; ++i) {
            p[i] += xi * s[i];
        }
        p[0] -= v;
        Fr f = s_blind * xi + p_blind;

        std::vector<Fr> b(n);
        b[0] = Fr(1);
        for (size_t i = 1; i < n; ++i) {
            b[i] = b[i - 1] * x;
        }
        // For the first rounds the folded generators are not materialized: after j rounds, G'[i] = sum_t w_t G[t * m +
        // i] with m = n / 2^j and w the tensor of the challenges (w_t = prod_r u_r^{bit_r(t)}, most significant bit
        // first), so L_j and R_j are MSMs over the original generators. After TENSOR_ROUNDS rounds the folded
        // generators are computed (one Straus MSM of 2^TENSOR_ROUNDS terms per generator) and folded explicitly, as
        // halo2 does, from then on.
        constexpr size_t TENSOR_ROUNDS = 4;
        const size_t tensor_rounds = std::min(k, TENSOR_ROUNDS);
        std::vector<Fr> tensor{ Fr(1) };
        std::vector<Commitment> g_folded;
        for (size_t j = 0; j < k; ++j) {
            const size_t m = n >> j;
            const size_t half = m >> 1;
            std::span<const Fr> p_lo(p.data(), half);
            std::span<const Fr> p_hi(p.data() + half, half);
            const Fr value_l = inner_product(p_hi, std::span<const Fr>(b.data(), half));
            const Fr value_r = inner_product(p_lo, std::span<const Fr>(b.data() + half, half));
            const Fr l_rand = Fr::random_element(&rng);
            const Fr r_rand = Fr::random_element(&rng);
            std::vector<Commitment> lr;
            if (j < tensor_rounds) {
                BB_BENCH_NAME("halo2_ipa/tensor_round");
                std::vector<std::vector<Fr>> scalars(2, std::vector<Fr>(n, Fr(0)));
                parallel_for(tensor.size(), [&](size_t t) {
                    for (size_t i = 0; i < half; ++i) {
                        scalars[0][(t * m) + i] = p_hi[i] * tensor[t];
                        scalars[1][(t * m) + half + i] = p_lo[i] * tensor[t];
                    }
                });
                lr = batch_msm(gens.g.subspan(0, n), scalars, /*handle_edge_cases=*/false);
            } else {
                if (j == tensor_rounds) {
                    BB_BENCH_NAME("halo2_ipa/materialize");
                    std::vector<SignedDigits> digits(tensor.size());
                    for (size_t t = 0; t < tensor.size(); ++t) {
                        digits[t] = signed_digits(tensor[t]);
                    }
                    std::vector<GroupElement> folded(m);
                    parallel_for_range(m, [&](size_t start, size_t end) {
                        std::vector<Commitment> terms(tensor.size());
                        for (size_t i = start; i < end; ++i) {
                            for (size_t t = 0; t < tensor.size(); ++t) {
                                terms[t] = gens.g[(t * m) + i];
                            }
                            folded[i] = straus(terms, digits);
                        }
                    });
                    g_folded = normalize(folded);
                }
                std::vector<std::vector<Fr>> lo{ std::vector<Fr>(p_hi.begin(), p_hi.end()) };
                std::vector<std::vector<Fr>> hi{ std::vector<Fr>(p_lo.begin(), p_lo.end()) };
                const std::span<const Commitment> g(g_folded);
                lr.push_back(batch_msm(g.subspan(0, half), lo, /*handle_edge_cases=*/false)[0]);
                lr.push_back(batch_msm(g.subspan(half, half), hi, /*handle_edge_cases=*/false)[0]);
            }
            const Commitment l_j =
                GroupElement(lr[0]) + GroupElement(gens.u) * (value_l * z) + GroupElement(gens.w) * l_rand;
            const Commitment r_j =
                GroupElement(lr[1]) + GroupElement(gens.u) * (value_r * z) + GroupElement(gens.w) * r_rand;
            const std::string index = std::to_string(j);
            transcript->send_to_verifier("IPA:L_" + index, l_j);
            transcript->send_to_verifier("IPA:R_" + index, r_j);
            const Fr u_j = transcript->template get_challenge<Fr>("IPA:u_" + index);
            const Fr u_j_inv = u_j.invert();

            for (size_t i = 0; i < half; ++i) {
                p[i] += p[i + half] * u_j_inv;
                b[i] += b[i + half] * u_j;
            }
            p.resize(half);
            b.resize(half);
            if (j < tensor_rounds) {
                // G'_lo + u_j G'_hi: the tensor doubles, with u_j on the entries whose new bit is 1.
                std::vector<Fr> next(2 * tensor.size());
                for (size_t t = 0; t < tensor.size(); ++t) {
                    next[2 * t] = tensor[t];
                    next[(2 * t) + 1] = tensor[t] * u_j;
                }
                tensor = std::move(next);
            } else if (half > 1) {
                BB_BENCH_NAME("halo2_ipa/fold");
                const std::array<SignedDigits, 2> digits{ signed_digits(Fr(1)), signed_digits(u_j) };
                std::vector<GroupElement> folded(half);
                parallel_for_range(half, [&](size_t start, size_t end) {
                    for (size_t i = start; i < end; ++i) {
                        const std::array<Commitment, 2> terms{ g_folded[i], g_folded[i + half] };
                        folded[i] = straus(terms, digits);
                    }
                });
                g_folded = normalize(folded);
            }

            f += l_rand * u_j_inv + r_rand * u_j;
        }
        BB_ASSERT_EQ(p.size(), size_t{ 1 });
        transcript->send_to_verifier("IPA:c", p[0]);
        transcript->send_to_verifier("IPA:f", f);
    }

    /**
     * @brief Verify an opening claim (commitment P, point x, value v) with a halo2 IPA proof of size n = |G|.
     */
    template <typename Transcript>
    static bool verify(const Generators& gens,
                       const OpeningClaim<Curve>& claim,
                       const std::shared_ptr<Transcript>& transcript)
    {
        const size_t n = gens.g.size();
        const size_t k = numeric::get_msb(n);
        const Fr& x = claim.opening_pair.challenge;
        const Fr& v = claim.opening_pair.evaluation;

        const auto s_commitment = transcript->template receive_from_prover<Commitment>("IPA:S");
        const Fr xi = transcript->template get_challenge<Fr>("IPA:xi");
        const Fr z = transcript->template get_challenge<Fr>("IPA:z");

        std::vector<Commitment> ls(k);
        std::vector<Commitment> rs(k);
        std::vector<Fr> u(k);
        for (size_t j = 0; j < k; ++j) {
            const std::string index = std::to_string(j);
            ls[j] = transcript->template receive_from_prover<Commitment>("IPA:L_" + index);
            rs[j] = transcript->template receive_from_prover<Commitment>("IPA:R_" + index);
            u[j] = transcript->template get_challenge<Fr>("IPA:u_" + index);
        }
        const Fr c = transcript->template receive_from_prover<Fr>("IPA:c");
        const Fr f = transcript->template receive_from_prover<Fr>("IPA:f");

        // One MSM over G (scalars -c s(u), with -v added to G_0) and the 2k + 4 other points.
        std::vector<Fr> scalars = compute_s(u, -c);
        scalars[0] -= v;
        std::vector<Commitment> points(gens.g.begin(), gens.g.begin() + static_cast<std::ptrdiff_t>(n));
        auto append = [&](const Fr& scalar, const Commitment& point) {
            scalars.push_back(scalar);
            points.push_back(point);
        };
        append(Fr(1), claim.commitment);
        append(xi, s_commitment);
        for (size_t j = 0; j < k; ++j) {
            append(u[j].invert(), ls[j]);
            append(u[j], rs[j]);
        }
        append(-(c * compute_b(x, u) * z), gens.u);
        append(-f, gens.w);
        // Points at infinity (e.g. a zero commitment) are skipped by the MSM.
        std::vector<Fr> msm_scalars;
        std::vector<Commitment> msm_points;
        for (size_t i = 0; i < scalars.size(); ++i) {
            if (!points[i].is_point_at_infinity()) {
                msm_scalars.push_back(scalars[i]);
                msm_points.push_back(points[i]);
            }
        }
        return GroupElement(msm(msm_scalars, msm_points)).is_point_at_infinity();
    }
};

/**
 * @brief The halo2 IPA generators of size n on Vesta: G from the commitment key (halo2's parameters), and halo2's W
 * and U.
 */
inline Halo2IPA<curve::Vesta>::Generators halo2_vesta_ipa_generators(const CommitmentKey<curve::Vesta>& ck, size_t n)
{
    using Commitment = curve::Vesta::AffineElement;
    static const Commitment w = halo2_vesta_w();
    static const Commitment u = halo2_vesta_u();
    return { std::span<const Commitment>(ck.get_monomial_points().data(), n), w, u };
}

} // namespace bb::zcash
