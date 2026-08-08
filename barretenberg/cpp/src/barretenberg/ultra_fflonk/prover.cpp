#include "barretenberg/ultra_fflonk/prover.hpp"

#include "barretenberg/common/thread.hpp"
#include "barretenberg/common/throw_or_abort.hpp"
#include "barretenberg/ultra_fflonk/relation_batch.hpp"
#include "barretenberg/ultra_honk/oink_prover.hpp"

#include <algorithm>

namespace bb::ultra_fflonk {

namespace {

using fflonk_plonk::coefficients_to_evaluations;
using fflonk_plonk::divide_by_vanishing;
using fflonk_plonk::evaluations_to_coefficients;
using fflonk_plonk::pack_columns;

/** @brief Ranges shorter than this are not worth handing to the thread pool. */
constexpr size_t PARALLEL_THRESHOLD = 1 << 10;

/** @brief The unshifted entities of the flavor, whose evaluations the relations read directly. */
constexpr size_t NUM_UNSHIFTED = Flavor::NUM_UNSHIFTED_ENTITIES;

/**
 * @brief A witness column in coefficient form, blinded so its two revealed evaluations say nothing.
 *
 * @details Adds `(b_0 X + b_1)(X^n - 1)`, which vanishes on `H` - so every row of the trace, and
 * therefore every relation, is untouched - while making the column's evaluations off `H` uniformly
 * random. Two blinders because every witness group is opened at two points.
 */
std::vector<FF> blinded_coefficients(const std::vector<FF>& lagrange,
                                     const size_t n,
                                     const EvaluationDomain<FF>& domain)
{
    std::vector<FF> coefficients = evaluations_to_coefficients(lagrange, domain);
    coefficients.resize(n + NUM_WITNESS_BLINDERS, FF::zero());

    const FF b_0 = FF::random_element();
    const FF b_1 = FF::random_element();
    coefficients[0] -= b_1;
    coefficients[1] -= b_0;
    coefficients[n] += b_1;
    coefficients[n + 1] += b_0;
    return coefficients;
}

/** @brief A Honk prover polynomial read out as `n` Lagrange values over the trace. */
std::vector<FF> lagrange_values(const Flavor::Polynomial& polynomial, const size_t n)
{
    std::vector<FF> lagrange(n, FF::zero());
    for (size_t row = 0; row < n; ++row) {
        lagrange[row] = polynomial[row];
    }
    return lagrange;
}

std::vector<FF> blinded_coefficients(const Flavor::Polynomial& polynomial,
                                     const size_t n,
                                     const EvaluationDomain<FF>& domain)
{
    return blinded_coefficients(lagrange_values(polynomial, n), n, domain);
}

/**
 * @brief `poly` evaluated on the coset `generator * H`.
 *
 * @details The quotient domain of size `8n` is eight cosets of `H`, indexed so that the `8n`-th root
 * `w_{8n}^{c + 8j}` is the `j`-th point of coset `c`. Working one coset at a time keeps only `n`
 * evaluations of each entity alive instead of `8n`, and it makes the shift free: multiplying a coset
 * point by `omega` moves it one position along the *same* coset, so `p(omega X)` on a coset is that
 * coset's evaluations rotated by one.
 *
 * Reducing modulo `X^n - 1` after scaling is exact for this purpose because the coset points satisfy
 * no relation other than being `n`-th roots of `generator^n`.
 */
std::vector<FF> evaluate_on_coset(std::span<const FF> poly, const FF& generator, const EvaluationDomain<FF>& domain)
{
    const size_t n = domain.size;
    std::vector<FF> folded(n, FF::zero());
    FF power = FF::one();
    for (size_t k = 0; k < poly.size(); ++k) {
        folded[k % n] += poly[k] * power;
        power *= generator;
    }
    return coefficients_to_evaluations(folded, domain);
}

/** @brief `poly += scalar * X^power`, growing the polynomial if it has to. */
void add_monomial(std::vector<FF>& poly, const size_t power, const FF& scalar)
{
    if (poly.size() <= power) {
        poly.resize(power + 1, FF::zero());
    }
    poly[power] += scalar;
}

} // namespace

FF prover_detail::open_and_batch(const ProvingKey& key,
                                 std::span<const std::vector<FF>* const> packed_groups,
                                 Transcript& transcript,
                                 Proof& proof)
{
    const VerificationKey& vk = key.verification_key;
    const size_t n = vk.circuit_size;
    const CommitmentKey<Curve>& commitment_key = *key.commitment_key;

    const FF xi = transcript.squeeze();
    if (xi.is_zero() || xi.pow(static_cast<uint64_t>(n)) == FF::one()) {
        throw_or_abort("ultra_fflonk: degenerate evaluation challenge");
    }
    const FF xi_omega = xi * vk.omega;

    fflonk_plonk::BatchedOpeningProver opening(packed_groups, GROUP_SHAPES, xi, xi_omega);
    proof.evaluations = flatten_evaluations(opening.evaluations());

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

Proof prove(ProvingKey& key)
{
    ProverInstance_<Flavor>& instance = *key.instance;
    const VerificationKey& vk = key.verification_key;
    const size_t n = vk.circuit_size;
    const EvaluationDomain<FF>& small_domain = *key.small_domain;
    const EvaluationDomain<FF>& large_domain = *key.large_domain;
    const size_t large_size = large_domain.size;
    const CommitmentKey<Curve>& commitment_key = *key.commitment_key;

    Proof proof;
    Transcript transcript;
    transcript.absorb(vk.hash());
    for (const FF& public_input : instance.public_inputs) {
        transcript.absorb(public_input);
    }

    const auto commit_group = [&](std::span<const std::vector<FF>> columns, std::vector<FF>& packed) {
        packed = pack_columns(columns);
        return commitment_key.commit(Polynomial<FF>(std::span<const FF>(packed)));
    };

    // ---------------------------------------------------------------------------------------------
    // Round 1: the first three wires. Mirrors OinkProver::commit_to_wires.
    // ---------------------------------------------------------------------------------------------
    std::vector<FF> packed_wires;
    std::array<std::vector<FF>, PACK_WIRES> wire_columns;
    {
        auto wires = instance.polynomials.get_wires();
        wire_columns[WIRES_W_L] = blinded_coefficients(wires[0], n, small_domain);
        wire_columns[WIRES_W_R] = blinded_coefficients(wires[1], n, small_domain);
        wire_columns[WIRES_W_O] = blinded_coefficients(wires[2], n, small_domain);
        proof.wires = commit_group(wire_columns, packed_wires);
    }
    transcript.absorb(proof.wires);

    instance.relation_parameters.compute_eta_powers(transcript.squeeze());
    instance.relation_parameters.rom_logup_gamma = transcript.squeeze();

    // ---------------------------------------------------------------------------------------------
    // Round 2: the memory records land in w_4, then the lookup counts and tags join it.
    // ---------------------------------------------------------------------------------------------
    OinkProver<Flavor>::add_ram_rom_memory_records_to_wire_4(instance);
    OinkProver<Flavor>::add_rom_logup_inverses_to_wire_4(instance);

    std::vector<FF> packed_memory;
    std::array<std::vector<FF>, PACK_MEMORY> memory_columns;
    {
        memory_columns[MEMORY_READ_COUNTS] =
            blinded_coefficients(instance.polynomials.lookup_read_counts(), n, small_domain);
        memory_columns[MEMORY_READ_TAGS] =
            blinded_coefficients(instance.polynomials.lookup_read_tags(), n, small_domain);
        memory_columns[MEMORY_W_4] = blinded_coefficients(instance.polynomials.w_4(), n, small_domain);
        proof.memory = commit_group(memory_columns, packed_memory);
    }
    transcript.absorb(proof.memory);

    {
        const FF beta = transcript.squeeze();
        const FF gamma = transcript.squeeze();
        instance.relation_parameters.compute_beta_powers(beta);
        instance.relation_parameters.gamma = gamma;
    }

    // ---------------------------------------------------------------------------------------------
    // Round 3: the lookup inverses, the permutation grand product, and one running sum for each
    // subrelation that Sumcheck would have enforced as a sum rather than row by row.
    // ---------------------------------------------------------------------------------------------
    OinkProver<Flavor>::compute_logderivative_inverses(instance);
    uint32_t z_perm_dup_count = 0;
    OinkProver<Flavor>::compute_grand_product_polynomial(instance, z_perm_dup_count);

    std::array<std::vector<FF>, NUM_RUNNING_SUMS> running_sum_lagrange;
    {
        for (std::vector<FF>& sum : running_sum_lagrange) {
            sum.assign(n, FF::zero());
        }
        std::vector<std::array<FF, NUM_RUNNING_SUMS>> per_row(n);
        parallel_for_range(
            n,
            [&](size_t start, size_t end) {
                for (size_t row = start; row < end; ++row) {
                    per_row[row] = summed_subrelations(instance.polynomials.get_row(row), instance.relation_parameters);
                }
            },
            PARALLEL_THRESHOLD);

        for (size_t j = 0; j < NUM_RUNNING_SUMS; ++j) {
            FF accumulator = FF::zero();
            for (size_t row = 0; row + 1 < n; ++row) {
                accumulator += per_row[row][j];
                running_sum_lagrange[j][row + 1] = accumulator;
            }
            // The identity is cyclic, so closing the loop is exactly the statement that the
            // subrelation sums to zero across the trace - which is what Sumcheck was enforcing.
            if (accumulator + per_row[n - 1][j] != FF::zero()) {
                throw_or_abort("ultra_fflonk: a trace-sum subrelation does not sum to zero");
            }
        }
    }

    std::vector<FF> packed_grand_product;
    std::array<std::vector<FF>, PACK_GRAND_PRODUCT> grand_product_columns;
    {
        grand_product_columns[GRAND_PRODUCT_LOOKUP_INVERSES] =
            blinded_coefficients(instance.polynomials.lookup_inverses(), n, small_domain);
        grand_product_columns[GRAND_PRODUCT_Z_PERM] =
            blinded_coefficients(instance.polynomials.z_perm(), n, small_domain);
        for (size_t j = 0; j < NUM_RUNNING_SUMS; ++j) {
            grand_product_columns[GRAND_PRODUCT_RUNNING_SUM + j] =
                blinded_coefficients(running_sum_lagrange[j], n, small_domain);
        }
        proof.grand_product = commit_group(grand_product_columns, packed_grand_product);
    }
    transcript.absorb(proof.grand_product);

    const FF alpha = transcript.squeeze();

    // ---------------------------------------------------------------------------------------------
    // Round 4: the quotient. Every entity is evaluated on the quotient domain one coset at a time,
    // the flavor's own relation accumulators are called at each point, and the result is divided by
    // the vanishing polynomial.
    // ---------------------------------------------------------------------------------------------
    std::array<std::vector<FF>, NUM_QUOTIENT_CHUNKS> quotient_columns;
    {
        std::array<std::vector<FF>, NUM_UNSHIFTED> entity_coefficients;
        for (size_t g = 0; g < NUM_PREPROCESSED_GROUPS; ++g) {
            for (size_t i = 0; i < PACK_PREPROCESSED; ++i) {
                entity_coefficients[(g * PACK_PREPROCESSED) + i] =
                    unpack_column(key.packed_preprocessed[g], PACK_PREPROCESSED, i);
            }
        }
        entity_coefficients[static_cast<size_t>(EntityId::w_l)] = wire_columns[WIRES_W_L];
        entity_coefficients[static_cast<size_t>(EntityId::w_r)] = wire_columns[WIRES_W_R];
        entity_coefficients[static_cast<size_t>(EntityId::w_o)] = wire_columns[WIRES_W_O];
        entity_coefficients[static_cast<size_t>(EntityId::w_4)] = memory_columns[MEMORY_W_4];
        entity_coefficients[static_cast<size_t>(EntityId::lookup_read_counts)] = memory_columns[MEMORY_READ_COUNTS];
        entity_coefficients[static_cast<size_t>(EntityId::lookup_read_tags)] = memory_columns[MEMORY_READ_TAGS];
        entity_coefficients[static_cast<size_t>(EntityId::z_perm)] = grand_product_columns[GRAND_PRODUCT_Z_PERM];
        entity_coefficients[static_cast<size_t>(EntityId::lookup_inverses)] =
            grand_product_columns[GRAND_PRODUCT_LOOKUP_INVERSES];

        std::vector<FF> numerator(large_size, FF::zero());
        std::array<std::vector<FF>, NUM_UNSHIFTED> coset;
        std::array<std::vector<FF>, NUM_RUNNING_SUMS> running_sum_coset;

        for (size_t c = 0; c < QUOTIENT_DOMAIN_FACTOR; ++c) {
            const FF generator = large_domain.root.pow(static_cast<uint64_t>(c));
            for (size_t e = 0; e < NUM_UNSHIFTED; ++e) {
                coset[e] = evaluate_on_coset(entity_coefficients[e], generator, small_domain);
            }
            for (size_t j = 0; j < NUM_RUNNING_SUMS; ++j) {
                running_sum_coset[j] =
                    evaluate_on_coset(grand_product_columns[GRAND_PRODUCT_RUNNING_SUM + j], generator, small_domain);
            }

            parallel_for_range(
                n,
                [&](size_t start, size_t end) {
                    for (size_t p = start; p < end; ++p) {
                        const size_t next = (p + 1) % n;

                        Flavor::AllValues row;
                        for (size_t e = 0; e < NUM_UNSHIFTED; ++e) {
                            row.data[e] = coset[e][p];
                        }
                        for (size_t s = 0; s < Flavor::NUM_SHIFTED_ENTITIES; ++s) {
                            row.data[NUM_UNSHIFTED + s] = coset[static_cast<size_t>(EntityId::w_l) + s][next];
                        }

                        std::array<FF, NUM_RUNNING_SUMS> sum_here{};
                        std::array<FF, NUM_RUNNING_SUMS> sum_next{};
                        for (size_t j = 0; j < NUM_RUNNING_SUMS; ++j) {
                            sum_here[j] = running_sum_coset[j][p];
                            sum_next[j] = running_sum_coset[j][next];
                        }

                        numerator[c + (QUOTIENT_DOMAIN_FACTOR * p)] =
                            quotient_numerator(row, instance.relation_parameters, alpha, sum_here, sum_next);
                    }
                },
                PARALLEL_THRESHOLD);
        }

        const std::vector<FF> numerator_coefficients = evaluations_to_coefficients(numerator, large_domain);
        numerator = {};

        std::vector<FF> quotient;
        if (!divide_by_vanishing(numerator_coefficients, n, quotient)) {
            throw_or_abort("ultra_fflonk: the relations do not vanish on the domain - the witness is invalid");
        }
        // The derivation in PROTOCOL.md bounds the quotient at `5n+6`; check the realised degree
        // rather than trusting it, because overflowing the split would silently drop coefficients.
        fflonk_plonk::trim(quotient);
        if (quotient.size() > NUM_QUOTIENT_CHUNKS * n) {
            throw_or_abort("ultra_fflonk: the quotient is wider than the group layout reserves for it");
        }

        for (size_t chunk = 0; chunk < NUM_QUOTIENT_CHUNKS; ++chunk) {
            const size_t lo = std::min(chunk * n, quotient.size());
            const size_t hi = std::min((chunk + 1) * n, quotient.size());
            quotient_columns[chunk] = std::vector<FF>(quotient.begin() + static_cast<std::ptrdiff_t>(lo),
                                                      quotient.begin() + static_cast<std::ptrdiff_t>(hi));
        }
        // Randomise inside the kernel of the split map: sum_c X^{cn} t_c is unchanged, but the
        // revealed chunk evaluations are no longer jointly determined by the quotient identity.
        for (size_t chunk = 0; chunk + 1 < NUM_QUOTIENT_CHUNKS; ++chunk) {
            const FF blinder = FF::random_element();
            add_monomial(quotient_columns[chunk], n, blinder);
            add_monomial(quotient_columns[chunk + 1], 0, -blinder);
        }
    }

    std::vector<FF> packed_quotient;
    proof.quotient = commit_group(quotient_columns, packed_quotient);
    transcript.absorb(proof.quotient);

    std::array<const std::vector<FF>*, NUM_GROUPS> packed_groups{};
    for (size_t g = 0; g < NUM_PREPROCESSED_GROUPS; ++g) {
        packed_groups[g] = &key.packed_preprocessed[g];
    }
    packed_groups[GROUP_WIRES] = &packed_wires;
    packed_groups[GROUP_MEMORY] = &packed_memory;
    packed_groups[GROUP_GRAND_PRODUCT] = &packed_grand_product;
    packed_groups[GROUP_QUOTIENT] = &packed_quotient;

    const FF xi = prover_detail::open_and_batch(key, packed_groups, transcript, proof);

    // ---------------------------------------------------------------------------------------------
    // The verifier's identity, checked here against the prover's own claimed evaluations.
    // ---------------------------------------------------------------------------------------------
    {
        const std::vector<GroupEvaluations> evaluations = unflatten_evaluations(proof.evaluations);
        const RunningSumEvaluations running_sums = to_running_sum_evaluations(evaluations);
        const FF numerator = quotient_numerator(to_all_values(evaluations),
                                                instance.relation_parameters,
                                                alpha,
                                                running_sums.at_xi,
                                                running_sums.at_xi_omega);

        const FF xi_pow_n = xi.pow(static_cast<uint64_t>(n));
        FF quotient_at_xi = FF::zero();
        for (size_t chunk = NUM_QUOTIENT_CHUNKS; chunk-- > 0;) {
            quotient_at_xi = (quotient_at_xi * xi_pow_n) + evaluations[GROUP_QUOTIENT].at_xi[chunk];
        }
        if (numerator != quotient_at_xi * (xi_pow_n - FF::one())) {
            throw_or_abort("ultra_fflonk: the batched identity does not hold at the challenge");
        }
    }

    return proof;
}

} // namespace bb::ultra_fflonk
