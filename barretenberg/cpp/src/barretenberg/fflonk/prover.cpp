#include "barretenberg/fflonk/prover.hpp"

#include "barretenberg/common/thread.hpp"
#include "barretenberg/common/throw_or_abort.hpp"

#include <algorithm>

namespace bb::fflonk_plonk {

namespace {

/** @brief Ranges shorter than this are not worth handing to the thread pool. */
constexpr size_t PARALLEL_THRESHOLD = 1 << 12;

/** @brief `poly += scalar * X^power`, growing the polynomial if it has to. */
void add_monomial(std::vector<FF>& poly, size_t power, const FF& scalar)
{
    if (poly.size() <= power) {
        poly.resize(power + 1, FF::zero());
    }
    poly[power] += scalar;
}

/** @brief Coefficients of `poly(wX)`, given the powers of `w`. */
std::vector<FF> shift_by_root(std::span<const FF> poly, std::span<const FF> root_powers)
{
    std::vector<FF> shifted(poly.size());
    for (size_t i = 0; i < poly.size(); ++i) {
        shifted[i] = poly[i] * root_powers[i % root_powers.size()];
    }
    return shifted;
}

/** @brief A group's columns packed into one polynomial, and its commitment. */
struct Group {
    std::vector<FF> packed;
    Commitment commitment = Commitment::infinity();
};

Group build_group(std::span<const std::vector<FF>> columns, const CommitmentKey<Curve>& commitment_key)
{
    Group group;
    group.packed = pack_columns(columns);
    group.commitment = commitment_key.commit(Polynomial<FF>(std::span<const FF>(group.packed)));
    return group;
}

} // namespace

FF prover_detail::open_and_batch(const ProvingKey& key,
                                 const std::array<const std::vector<FF>*, NUM_GROUPS>& packed_groups,
                                 Transcript& transcript,
                                 Proof& proof)
{
    const VerificationKey& vk = key.verification_key;
    const size_t n = vk.circuit_size;
    const CommitmentKey<Curve>& commitment_key = *key.commitment_key;

    const FF xi = transcript.squeeze();
    if (xi.is_zero() || xi.pow(static_cast<uint64_t>(n)) == FF::one()) {
        throw_or_abort("fflonk: degenerate evaluation challenge");
    }
    const FF xi_omega = xi * vk.omega;

    // Every group's evaluations are the residue of its packed polynomial modulo the group's
    // vanishing polynomial, so the division that builds W also produces what is sent.
    BatchedOpeningProver opening(packed_groups, GROUP_SHAPES, xi, xi_omega);
    const std::vector<GroupEvaluations>& evaluations = opening.evaluations();

    for (size_t i = 0; i < PACK_PREPROCESSED; ++i) {
        proof.evaluations[EVAL_Q_L + i] = evaluations[0].at_xi[i];
    }
    for (size_t i = 0; i < PACK_WIRES; ++i) {
        proof.evaluations[EVAL_A + i] = evaluations[1].at_xi[i];
    }
    proof.evaluations[EVAL_Z] = evaluations[2].at_xi[0];
    proof.evaluations[EVAL_Z_OMEGA] = evaluations[2].at_xi_omega[0];
    for (size_t i = 0; i < PACK_QUOTIENTS; ++i) {
        proof.evaluations[EVAL_T1 + i] = evaluations[3].at_xi[i];
    }

    for (const FF& evaluation : proof.evaluations) {
        transcript.absorb(evaluation);
    }
    const FF nu = transcript.squeeze();

    const std::vector<FF> w_polynomial = opening.compute_w(nu);
    proof.w = commitment_key.commit(Polynomial<FF>(std::span<const FF>(w_polynomial)));

    transcript.absorb(proof.w);
    const FF y = transcript.squeeze();

    const std::vector<FF> w_prime_polynomial = opening.compute_w_prime(nu, y, w_polynomial);
    proof.w_prime = commitment_key.commit(Polynomial<FF>(std::span<const FF>(w_prime_polynomial)));

    return xi;
}

Proof prove(const ProvingKey& key)
{
    const VerificationKey& vk = key.verification_key;
    const size_t n = vk.circuit_size;
    const EvaluationDomain<FF>& small_domain = *key.small_domain;
    const EvaluationDomain<FF>& large_domain = *key.large_domain;
    const size_t large_size = large_domain.size;
    const CommitmentKey<Curve>& commitment_key = *key.commitment_key;

    const FF omega = vk.omega;
    const FF n_inverse = FF(static_cast<uint64_t>(n)).invert();
    const std::array<FF, NUM_WIRES> coset_shifts = { FF::one(), vk.k1, vk.k2 };

    Transcript transcript;
    transcript.absorb(vk.hash());
    for (const FF& public_input : key.trace.public_inputs) {
        transcript.absorb(public_input);
    }

    // ---------------------------------------------------------------------------------------------
    // Round 1: the wires, blinded, and the gate quotient.
    // ---------------------------------------------------------------------------------------------
    std::array<std::vector<FF>, NUM_WIRES> wires;
    for (size_t wire = 0; wire < NUM_WIRES; ++wire) {
        wires[wire] = evaluations_to_coefficients(key.trace.wires[wire], small_domain);
        wires[wire].resize(n + 2, FF::zero());

        // + (b_0 X + b_1)(X^n - 1), which vanishes on H and so changes nothing the circuit says.
        const FF b_0 = FF::random_element();
        const FF b_1 = FF::random_element();
        wires[wire][0] -= b_1;
        wires[wire][1] -= b_0;
        wires[wire][n] += b_1;
        wires[wire][n + 1] += b_0;
    }

    std::array<std::vector<FF>, NUM_WIRES> wire_evaluations;
    for (size_t wire = 0; wire < NUM_WIRES; ++wire) {
        wire_evaluations[wire] = coefficients_to_evaluations(wires[wire], large_domain);
    }

    std::vector<FF> t0;
    {
        std::vector<FF> gate_numerator(large_size, FF::zero());
        const std::vector<FF>& a = wire_evaluations[0];
        const std::vector<FF>& b = wire_evaluations[1];
        const std::vector<FF>& c = wire_evaluations[2];

        // Accumulated one selector at a time so only one large-domain selector buffer is ever alive.
        const auto accumulate = [&](const std::vector<FF>& selector_coefficients, const auto& term) {
            const std::vector<FF> selector = coefficients_to_evaluations(selector_coefficients, large_domain);
            parallel_for_range(
                large_size,
                [&](size_t start, size_t end) {
                    for (size_t i = start; i < end; ++i) {
                        gate_numerator[i] += selector[i] * term(i);
                    }
                },
                PARALLEL_THRESHOLD);
        };

        accumulate(key.preprocessed.q_m(), [&](size_t i) { return a[i] * b[i]; });
        accumulate(key.preprocessed.q_l(), [&](size_t i) { return a[i]; });
        accumulate(key.preprocessed.q_r(), [&](size_t i) { return b[i]; });
        accumulate(key.preprocessed.q_o(), [&](size_t i) { return c[i]; });
        accumulate(key.preprocessed.q_c(), [](size_t) { return FF::one(); });
        accumulate(key.public_input_poly, [](size_t) { return FF::one(); });

        const std::vector<FF> numerator_coefficients = evaluations_to_coefficients(gate_numerator, large_domain);
        if (!divide_by_vanishing(numerator_coefficients, n, t0)) {
            throw_or_abort("fflonk: the gate identity does not vanish on the domain - the witness is invalid");
        }
        trim(t0);
    }

    const size_t t0_split_point = std::min(n, t0.size());
    std::vector<FF> t0_lo(t0.begin(), t0.begin() + static_cast<std::ptrdiff_t>(t0_split_point));
    std::vector<FF> t0_hi(t0.begin() + static_cast<std::ptrdiff_t>(t0_split_point), t0.end());
    {
        // Randomise inside the kernel of the split map: t0_lo + X^n t0_hi is unchanged, but the two
        // revealed evaluations are no longer jointly determined by the gate identity.
        const FF b_split = FF::random_element();
        add_monomial(t0_lo, n, b_split);
        add_monomial(t0_hi, 0, -b_split);
    }
    t0 = {};

    const std::array<std::vector<FF>, PACK_WIRES> wire_columns = { wires[0], wires[1], wires[2], t0_lo, t0_hi };
    const Group group_wires = build_group(wire_columns, commitment_key);

    transcript.absorb(group_wires.commitment);
    const FF beta = transcript.squeeze();
    const FF gamma = transcript.squeeze();

    // ---------------------------------------------------------------------------------------------
    // Round 2: the grand product, and the two permutation quotients.
    // ---------------------------------------------------------------------------------------------
    const std::vector<FF> domain_powers = compute_power_table(omega, n);

    std::vector<FF> grand_product;
    {
        std::vector<FF> numerators(n);
        std::vector<FF> denominators(n);
        for (size_t row = 0; row < n; ++row) {
            FF numerator = FF::one();
            FF denominator = FF::one();
            for (size_t wire = 0; wire < NUM_WIRES; ++wire) {
                const FF value = key.trace.wires[wire][row];
                numerator *= value + beta * coset_shifts[wire] * domain_powers[row] + gamma;
                denominator *= value + beta * key.sigma_lagrange[wire][row] + gamma;
            }
            numerators[row] = numerator;
            denominators[row] = denominator;
        }
        FF::batch_invert(denominators);

        std::vector<FF> lagrange(n);
        lagrange[0] = FF::one();
        for (size_t row = 0; row + 1 < n; ++row) {
            lagrange[row + 1] = lagrange[row] * numerators[row] * denominators[row];
        }
        // The running product over every row must return to 1, which happens exactly when sigma is
        // a permutation and the witness respects the copy constraints.
        if (lagrange[n - 1] * numerators[n - 1] * denominators[n - 1] != FF::one()) {
            throw_or_abort("fflonk: the grand product does not close - the copy constraints are violated");
        }

        grand_product = evaluations_to_coefficients(lagrange, small_domain);
        grand_product.resize(n + 3, FF::zero());

        // + (b_0 X^2 + b_1 X + b_2)(X^n - 1): three blinders for two revealed evaluations.
        const FF b_0 = FF::random_element();
        const FF b_1 = FF::random_element();
        const FF b_2 = FF::random_element();
        grand_product[0] -= b_2;
        grand_product[1] -= b_1;
        grand_product[2] -= b_0;
        grand_product[n] += b_2;
        grand_product[n + 1] += b_1;
        grand_product[n + 2] += b_0;
    }

    // T1 = (z - 1) L_0 / Z_H = (z - 1) / (n (X - 1)), since L_0 = (X^n - 1)/(n(X - 1)).
    std::vector<FF> t1;
    {
        std::vector<FF> shifted = grand_product;
        shifted[0] -= FF::one();
        FF remainder = FF::zero();
        t1 = divide_by_linear(shifted, FF::one(), remainder);
        if (!remainder.is_zero()) {
            throw_or_abort("fflonk: the grand product does not start at one");
        }
        for (FF& coefficient : t1) {
            coefficient *= n_inverse;
        }
        trim(t1);
    }

    std::vector<FF> t2;
    {
        const std::vector<FF> large_domain_powers = compute_power_table(large_domain.root, large_size);
        const std::vector<FF> z_evaluations = coefficients_to_evaluations(grand_product, large_domain);
        const std::vector<FF> z_shifted_evaluations =
            coefficients_to_evaluations(shift_by_root(grand_product, domain_powers), large_domain);

        // The sigma side, folded one wire at a time so only one sigma buffer is alive at once.
        std::vector<FF> sigma_product(large_size, FF::one());
        for (size_t wire = 0; wire < NUM_WIRES; ++wire) {
            const std::vector<FF> sigma =
                coefficients_to_evaluations(key.preprocessed.columns[PREPROCESSED_SIGMA_OFFSET + wire], large_domain);
            const std::vector<FF>& value = wire_evaluations[wire];
            parallel_for_range(
                large_size,
                [&](size_t start, size_t end) {
                    for (size_t i = start; i < end; ++i) {
                        sigma_product[i] *= value[i] + beta * sigma[i] + gamma;
                    }
                },
                PARALLEL_THRESHOLD);
        }

        std::vector<FF> permutation_numerator(large_size, FF::zero());
        parallel_for_range(
            large_size,
            [&](size_t start, size_t end) {
                for (size_t i = start; i < end; ++i) {
                    FF identity_product = FF::one();
                    for (size_t wire = 0; wire < NUM_WIRES; ++wire) {
                        identity_product *=
                            wire_evaluations[wire][i] + beta * coset_shifts[wire] * large_domain_powers[i] + gamma;
                    }
                    permutation_numerator[i] =
                        identity_product * z_evaluations[i] - sigma_product[i] * z_shifted_evaluations[i];
                }
            },
            PARALLEL_THRESHOLD);

        const std::vector<FF> numerator_coefficients = evaluations_to_coefficients(permutation_numerator, large_domain);
        if (!divide_by_vanishing(numerator_coefficients, n, t2)) {
            throw_or_abort("fflonk: the permutation identity does not vanish on the domain");
        }
        trim(t2);
    }

    for (std::vector<FF>& evaluations : wire_evaluations) {
        evaluations = {};
    }

    std::array<std::vector<FF>, 3> t2_split;
    {
        const auto slice = [&](size_t from, size_t to) {
            const size_t lo = std::min(from, t2.size());
            const size_t hi = std::min(to, t2.size());
            return std::vector<FF>(t2.begin() + static_cast<std::ptrdiff_t>(lo),
                                   t2.begin() + static_cast<std::ptrdiff_t>(hi));
        };
        t2_split[0] = slice(0, n);
        t2_split[1] = slice(n, 2 * n);
        t2_split[2] = slice(2 * n, t2.size());

        // Kernel randomisation again: t2_lo + X^n t2_mid + X^2n t2_hi is unchanged.
        const FF b_0 = FF::random_element();
        const FF b_1 = FF::random_element();
        add_monomial(t2_split[0], n, b_0);
        add_monomial(t2_split[1], 0, -b_0);
        add_monomial(t2_split[1], n, b_1);
        add_monomial(t2_split[2], 0, -b_1);
    }
    t2 = {};

    const std::array<std::vector<FF>, PACK_GRAND_PRODUCT> grand_product_columns = { grand_product };
    const std::array<std::vector<FF>, PACK_QUOTIENTS> quotient_columns = { t1, t2_split[0], t2_split[1], t2_split[2] };
    const Group group_grand_product = build_group(grand_product_columns, commitment_key);
    const Group group_quotients = build_group(quotient_columns, commitment_key);

    Proof proof;
    proof.c1 = group_wires.commitment;
    proof.c2 = group_grand_product.commitment;
    proof.c3 = group_quotients.commitment;

    transcript.absorb(proof.c2);
    transcript.absorb(proof.c3);

    const std::array<const std::vector<FF>*, NUM_GROUPS> packed_groups = {
        &key.packed_preprocessed, &group_wires.packed, &group_grand_product.packed, &group_quotients.packed
    };
    const FF xi = prover_detail::open_and_batch(key, packed_groups, transcript, proof);

    // ---------------------------------------------------------------------------------------------
    // The verifier's three identities, checked here against the prover's own claimed evaluations.
    // ---------------------------------------------------------------------------------------------
    {
        const auto& e = proof.evaluations;
        const FF xi_pow_n = xi.pow(static_cast<uint64_t>(n));
        const FF vanishing = xi_pow_n - FF::one();
        const FF lagrange_first = vanishing * n_inverse * (xi - FF::one()).invert();
        const FF public_input_evaluation = evaluate(key.public_input_poly, xi);

        const FF gate = e[EVAL_Q_L] * e[EVAL_A] + e[EVAL_Q_R] * e[EVAL_B] + e[EVAL_Q_O] * e[EVAL_C] +
                        e[EVAL_Q_M] * e[EVAL_A] * e[EVAL_B] + e[EVAL_Q_C] + public_input_evaluation;
        if (gate != (e[EVAL_T0_LO] + xi_pow_n * e[EVAL_T0_HI]) * vanishing) {
            throw_or_abort("fflonk: the gate identity does not hold at the challenge");
        }
        if ((e[EVAL_Z] - FF::one()) * lagrange_first != e[EVAL_T1] * vanishing) {
            throw_or_abort("fflonk: the grand-product start identity does not hold at the challenge");
        }
        const FF identity_side = (e[EVAL_A] + beta * xi + gamma) * (e[EVAL_B] + beta * vk.k1 * xi + gamma) *
                                 (e[EVAL_C] + beta * vk.k2 * xi + gamma) * e[EVAL_Z];
        const FF sigma_side = (e[EVAL_A] + beta * e[EVAL_S_1] + gamma) * (e[EVAL_B] + beta * e[EVAL_S_2] + gamma) *
                              (e[EVAL_C] + beta * e[EVAL_S_3] + gamma) * e[EVAL_Z_OMEGA];
        const FF quotient = e[EVAL_T2_LO] + xi_pow_n * e[EVAL_T2_MID] + xi_pow_n * xi_pow_n * e[EVAL_T2_HI];
        if (identity_side - sigma_side != quotient * vanishing) {
            throw_or_abort("fflonk: the permutation identity does not hold at the challenge");
        }
    }

    return proof;
}

} // namespace bb::fflonk_plonk
