#include "barretenberg/fflonk/keys.hpp"

#include "barretenberg/common/throw_or_abort.hpp"
#include "barretenberg/crypto/keccak/keccak.hpp"
#include "barretenberg/numeric/bitop/get_msb.hpp"

#include <algorithm>

namespace bb::fflonk_plonk {

FF VerificationKey::hash() const
{
    return crypto::Keccak::hash({ uint256_t(circuit_size),
                                  uint256_t(num_public_inputs),
                                  static_cast<uint256_t>(k1),
                                  static_cast<uint256_t>(k2),
                                  c0.is_point_at_infinity() ? uint256_t(0) : static_cast<uint256_t>(c0.x),
                                  c0.is_point_at_infinity() ? uint256_t(0) : static_cast<uint256_t>(c0.y) });
}

namespace {

/**
 * @brief The copy-constraint permutation, in Lagrange form.
 *
 * @details Slot `(wire, row)` is identified with the field element `k_wire * w^row`. The three
 * cosets `H`, `k1 H`, `k2 H` are disjoint - asserted by the caller - so that identification is
 * injective over all `3n` slots, which is exactly what makes the grand product a permutation check
 * rather than a multiset check over colliding labels.
 *
 * Slots holding the same variable form a cycle, and `sigma` sends each slot to the next one in its
 * cycle; a variable used once is a fixed point.
 */
std::array<std::vector<FF>, NUM_WIRES> compute_sigma_polynomials(const Trace& trace,
                                                                 const std::array<FF, NUM_WIRES>& coset_shifts,
                                                                 const std::vector<FF>& domain_powers)
{
    const size_t n = trace.circuit_size;
    const size_t num_slots = NUM_WIRES * n;

    uint32_t max_variable = 0;
    for (size_t wire = 0; wire < NUM_WIRES; ++wire) {
        for (const uint32_t variable : trace.wire_variables[wire]) {
            max_variable = std::max(max_variable, variable);
        }
    }

    constexpr uint32_t NO_SLOT = std::numeric_limits<uint32_t>::max();
    std::vector<uint32_t> first_slot(static_cast<size_t>(max_variable) + 1, NO_SLOT);
    std::vector<uint32_t> last_slot(static_cast<size_t>(max_variable) + 1, NO_SLOT);
    std::vector<uint32_t> next_slot(num_slots, NO_SLOT);

    for (size_t wire = 0; wire < NUM_WIRES; ++wire) {
        for (size_t row = 0; row < n; ++row) {
            const auto slot = static_cast<uint32_t>(wire * n + row);
            const uint32_t variable = trace.wire_variables[wire][row];
            if (first_slot[variable] == NO_SLOT) {
                first_slot[variable] = slot;
            } else {
                next_slot[last_slot[variable]] = slot;
            }
            last_slot[variable] = slot;
        }
    }
    for (uint32_t variable = 0; variable <= max_variable; ++variable) {
        if (first_slot[variable] != NO_SLOT) {
            next_slot[last_slot[variable]] = first_slot[variable];
        }
    }

    std::array<std::vector<FF>, NUM_WIRES> sigma;
    for (size_t wire = 0; wire < NUM_WIRES; ++wire) {
        sigma[wire].resize(n);
        for (size_t row = 0; row < n; ++row) {
            const uint32_t target = next_slot[wire * n + row];
            BB_ASSERT_NEQ(target, NO_SLOT, "every wire slot belongs to a cycle");
            sigma[wire][row] = coset_shifts[target / n] * domain_powers[target % n];
        }
    }
    return sigma;
}

} // namespace

ProvingKey preprocess(const CircuitBuilder& builder, const size_t minimum_size)
{
    return preprocess(builder.build_trace(minimum_size));
}

ProvingKey preprocess(Trace trace)
{
    const size_t n = trace.circuit_size;
    // The lower bound matches the Solidity verifier's, so a key one accepts the other accepts too.
    if (n < MIN_CIRCUIT_SIZE || (n & (n - 1)) != 0) {
        throw_or_abort("fflonk: circuit size must be a power of two of at least 8");
    }
    if (trace.num_public_inputs >= n) {
        throw_or_abort("fflonk: more public inputs than rows");
    }

    ProvingKey key;
    key.small_domain = std::make_shared<EvaluationDomain<FF>>(n);
    key.small_domain->compute_lookup_table();
    key.large_domain = std::make_shared<EvaluationDomain<FF>>(QUOTIENT_DOMAIN_FACTOR * n);
    key.large_domain->compute_lookup_table();

    const FF omega = key.small_domain->root;
    const FF k1 = FF::coset_generator();
    const FF k2 = k1 * k1;

    // The permutation labels slots as k_j * w^i; if two of the three cosets met, distinct slots
    // would share a label and the grand product would accept a non-permutation.
    const auto in_domain = [n](const FF& value) { return value.pow(static_cast<uint64_t>(n)) == FF::one(); };
    if (in_domain(k1) || in_domain(k2) || in_domain(k1 * k2.invert())) {
        throw_or_abort("fflonk: coset shifts collide with the evaluation domain at this circuit size");
    }

    const std::vector<FF> domain_powers = compute_power_table(omega, n);
    key.sigma_lagrange = compute_sigma_polynomials(trace, { FF::one(), k1, k2 }, domain_powers);

    key.public_input_lagrange.assign(n, FF::zero());
    for (size_t i = 0; i < trace.num_public_inputs; ++i) {
        key.public_input_lagrange[i] = -trace.public_inputs[i];
    }
    key.public_input_poly = evaluations_to_coefficients(key.public_input_lagrange, *key.small_domain);

    const std::array<const std::vector<FF>*, PACK_PREPROCESSED> lagrange_columns = {
        &trace.q_l,
        &trace.q_r,
        &trace.q_o,
        &trace.q_m,
        &trace.q_c,
        &key.sigma_lagrange[0],
        &key.sigma_lagrange[1],
        &key.sigma_lagrange[2],
    };
    for (size_t i = 0; i < PACK_PREPROCESSED; ++i) {
        key.preprocessed.columns[i] = evaluations_to_coefficients(*lagrange_columns[i], *key.small_domain);
    }

    key.packed_preprocessed = pack_columns(key.preprocessed.columns);

    key.commitment_key = std::make_shared<CommitmentKey<Curve>>(SRS_SIZE_FACTOR * n);
    const Commitment c0 = key.commitment_key->commit(Polynomial<FF>(std::span<const FF>(key.packed_preprocessed)));

    key.verification_key = VerificationKey{
        .circuit_size = n, .num_public_inputs = trace.num_public_inputs, .omega = omega, .k1 = k1, .k2 = k2, .c0 = c0
    };
    key.trace = std::move(trace);
    return key;
}

} // namespace bb::fflonk_plonk
