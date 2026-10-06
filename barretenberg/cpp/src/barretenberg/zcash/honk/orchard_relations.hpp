#pragma once

#include "barretenberg/relations/relation_types.hpp"
#include "barretenberg/zcash/halo2/gates.hpp"

#include <array>
#include <tuple>

/**
 * @file orchard_relations.hpp
 * @brief Honk relations of the Orchard flavor: the halo2 custom gates, the copy-constraint permutation and the two
 * lookup arguments of the Action circuit.
 *
 * @details All relations read entities through the static accessors of the flavor (`Flavor::advice(in, c)` etc.), so
 * they are written once for the prover (univariates) and the verifier (field elements).
 */
namespace bb::zcash {

namespace detail {
// Index a tuple of identical subrelation accumulators at runtime.
template <typename Tuple> auto accumulator_pointers(Tuple& tuple)
{
    constexpr size_t N = std::tuple_size_v<Tuple>;
    using T = std::remove_reference_t<decltype(std::get<0>(tuple))>;
    std::array<T*, N> out{};
    [&]<size_t... I>(std::index_sequence<I...>) {
        ((out[I] = &std::get<I>(tuple)), ...);
    }(std::make_index_sequence<N>{});
    return out;
}

template <typename T> bool is_zero_value(const T& v)
{
    if constexpr (requires { v.evaluations; }) {
        for (const auto& e : v.evaluations) {
            if (!e.is_zero()) {
                return false;
            }
        }
        return true;
    } else {
        return v.is_zero();
    }
}
} // namespace detail

/**
 * @brief The halo2 custom gates: one subrelation per polynomial identity of `OrchardGates`.
 */
template <typename Flavor_, typename FF_> class OrchardGateRelationImpl {
  public:
    using FF = FF_;
    using Cycle = typename Flavor_::Cycle;
    static constexpr size_t NUM = halo2::OrchardGates<Cycle>::NUM_GATE_CONSTRAINTS;
    static constexpr size_t LENGTH = halo2::MAX_GATE_DEGREE + 1;

    static constexpr std::array<size_t, NUM> SUBRELATION_PARTIAL_LENGTHS = [] {
        std::array<size_t, NUM> out{};
        out.fill(LENGTH);
        return out;
    }();

    template <typename AllEntities> static bool skip(const AllEntities& in)
    {
        for (const auto& s : Flavor_::gate_selectors(in)) {
            if (!detail::is_zero_value(s)) {
                return false;
            }
        }
        return true;
    }

    template <typename ContainerOverSubrelations, typename AllEntities, typename Parameters>
    static void accumulate(ContainerOverSubrelations& evals,
                           const AllEntities& in,
                           const Parameters& /*params*/,
                           const FF& scaling_factor)
    {
        using Accumulator = std::tuple_element_t<0, ContainerOverSubrelations>;
        using View = typename Accumulator::View;
        auto acc = detail::accumulator_pointers(evals);
        const auto advice = Flavor_::advice(in);
        const auto advice_shift = Flavor_::advice_shift(in);
        const auto advice_shift2 = Flavor_::advice_shift2(in);
        const auto fixed = Flavor_::fixed_columns(in);
        const auto selectors = Flavor_::gate_selectors(in);
        struct Input {
            decltype(advice)& a;
            decltype(advice_shift)& a1;
            decltype(advice_shift2)& a2;
            decltype(fixed)& f;
            decltype(selectors)& q;
            Accumulator adv(size_t col, size_t off) const
            {
                switch (off) {
                case 0:
                    return Accumulator(View(a[col]));
                case 1:
                    return Accumulator(View(a1[col]));
                default:
                    return Accumulator(View(a2[col]));
                }
            }
            Accumulator fix(size_t col) const { return Accumulator(View(f[col])); }
            Accumulator sel(size_t s) const { return Accumulator(View(q[s])); }
            bool active(size_t s) const { return !detail::is_zero_value(q[s]); }
        } input{ advice, advice_shift, advice_shift2, fixed, selectors };
        halo2::OrchardGates<Cycle>::template evaluate<Accumulator>(
            input, [&](size_t idx, const Accumulator& value) { *acc[idx] += value * scaling_factor; });
    }
};

/**
 * @brief Copy constraints over the 14 equality-enabled columns, split into two chunks of 7 columns linked by the
 * intermediate product z_perm_mid (halo2 also chunks its permutation argument to bound the degree):
 *
 *   (z_perm + L_first) * prod_{j<7} (w_j + beta id_j + gamma) = z_perm_mid * prod_{j<7} (w_j + beta sigma_j + gamma)
 *   z_perm_mid * prod_{j>=7} (w_j + beta id_j + gamma) = (z_perm_shift + L_last delta) * prod_{j>=7} (w_j + beta
 * sigma_j + gamma)
 *
 * where delta is the public-input correction computed by the verifier.
 */
template <typename Flavor_, typename FF_> class OrchardPermutationRelationImpl {
  public:
    using FF = FF_;
    static constexpr size_t CHUNK = 7;
    static constexpr std::array<size_t, 2> SUBRELATION_PARTIAL_LENGTHS{ CHUNK + 2, CHUNK + 2 };

    template <typename ContainerOverSubrelations, typename AllEntities, typename Parameters>
    static void accumulate(ContainerOverSubrelations& evals,
                           const AllEntities& in,
                           const Parameters& params,
                           const FF& scaling_factor)
    {
        using Accumulator = std::tuple_element_t<0, ContainerOverSubrelations>;
        using View = typename Accumulator::View;
        const auto& beta = params.beta;
        const auto& gamma = params.gamma;
        const auto columns = Flavor_::permutation_columns(in);
        const auto sigmas = Flavor_::sigmas(in);
        const auto ids = Flavor_::ids(in);

        auto chunk_products = [&](size_t start) {
            Accumulator num(View(columns[start]) + View(ids[start]) * beta + gamma);
            Accumulator den(View(columns[start]) + View(sigmas[start]) * beta + gamma);
            for (size_t j = start + 1; j < start + CHUNK; ++j) {
                num *= Accumulator(View(columns[j]) + View(ids[j]) * beta + gamma);
                den *= Accumulator(View(columns[j]) + View(sigmas[j]) * beta + gamma);
            }
            return std::make_pair(num, den);
        };
        const auto z_perm = Accumulator(View(in.z_perm));
        const auto z_perm_mid = Accumulator(View(in.z_perm_mid));
        const auto z_perm_shift = Accumulator(View(in.z_perm_shift));
        const auto lagrange_first = Accumulator(View(in.lagrange_first));
        const auto lagrange_last = Accumulator(View(in.lagrange_last));
        {
            auto [num, den] = chunk_products(0);
            std::get<0>(evals) += ((z_perm + lagrange_first) * num - z_perm_mid * den) * scaling_factor;
        }
        {
            auto [num, den] = chunk_products(CHUNK);
            std::get<1>(evals) +=
                (z_perm_mid * num - (z_perm_shift + lagrange_last * params.public_input_delta) * den) * scaling_factor;
        }
    }
};

/**
 * @brief Log-derivative lookup of the Sinsemilla generator table, read by both Sinsemilla chips.
 * @details With R_c = gamma + m_c + eta x_c + eta^2 y_c the read tuple of chip c (active iff q_sinsemilla1 of the
 * chip), T = gamma + idx + eta x + eta^2 y (active iff q_table) and I the inverse column:
 *   I * R_1 * R_2 * T - exists = 0            on every row, exists = 1 - (1 - q_1)(1 - q_2)(1 - q_table),
 *   sum over rows of I * (q_1 R_2 T + q_2 R_1 T - counts R_1 R_2) = 0.
 */
template <typename Flavor_, typename FF_> class OrchardSinsemillaLookupRelationImpl {
  public:
    using FF = FF_;
    using Cycle = typename Flavor_::Cycle;
    static constexpr std::array<size_t, 2> SUBRELATION_PARTIAL_LENGTHS{ 9, 9 };
    static constexpr std::array<bool, 2> SUBRELATION_LINEARLY_INDEPENDENT{ true, false };

    template <typename AllEntities> static bool skip(const AllEntities& in)
    {
        return detail::is_zero_value(in.q_sinsemilla1_1) && detail::is_zero_value(in.q_sinsemilla1_2) &&
               detail::is_zero_value(in.q_table) && detail::is_zero_value(in.lookup_read_counts_sinsemilla);
    }

    template <typename Accumulator, typename AllEntities, typename Parameters>
    static Accumulator read_term(const AllEntities& in, const Parameters& params, size_t chip)
    {
        using View = typename Accumulator::View;
        const auto advice = Flavor_::advice(in);
        const auto advice_shift = Flavor_::advice_shift(in);
        const auto fixed = Flavor_::fixed_columns(in);
        struct Input {
            decltype(advice)& a;
            decltype(advice_shift)& a1;
            decltype(fixed)& f;
            Accumulator adv(size_t col, size_t off) const
            {
                return off == 0 ? Accumulator(View(a[col])) : Accumulator(View(a1[col]));
            }
            Accumulator fix(size_t col) const { return Accumulator(View(f[col])); }
        } input{ advice, advice_shift, fixed };
        const auto [m, x, y] = halo2::OrchardGates<Cycle>::template sinsemilla_lookup_value<Accumulator>(input, chip);
        return m + x * params.eta + y * params.eta_two + params.gamma;
    }

    template <typename Accumulator, typename AllEntities, typename Parameters>
    static Accumulator write_term(const AllEntities& in, const Parameters& params)
    {
        using View = typename Accumulator::View;
        return Accumulator(View(in.table_idx) + View(in.table_x) * params.eta + View(in.table_y) * params.eta_two +
                           params.gamma);
    }

    template <typename ContainerOverSubrelations, typename AllEntities, typename Parameters>
    static void accumulate(ContainerOverSubrelations& evals,
                           const AllEntities& in,
                           const Parameters& params,
                           const FF& scaling_factor)
    {
        using Accumulator = std::tuple_element_t<0, ContainerOverSubrelations>;
        using View = typename Accumulator::View;
        const auto q_1 = Accumulator(View(in.q_sinsemilla1_1));
        const auto q_2 = Accumulator(View(in.q_sinsemilla1_2));
        const auto q_table = Accumulator(View(in.q_table));
        const auto inverse = Accumulator(View(in.lookup_inverses_sinsemilla));
        const auto counts = Accumulator(View(in.lookup_read_counts_sinsemilla));
        const auto r_1 = read_term<Accumulator>(in, params, 0);
        const auto r_2 = read_term<Accumulator>(in, params, 1);
        const auto t = write_term<Accumulator>(in, params);
        const auto exists = -((-q_1 + FF(1)) * (-q_2 + FF(1)) * (-q_table + FF(1))) + FF(1);
        std::get<0>(evals) += (inverse * r_1 * r_2 * t - exists) * scaling_factor;
        std::get<1>(evals) += inverse * (q_1 * r_2 * t + q_2 * r_1 * t - counts * r_1 * r_2);
    }
};

/**
 * @brief Log-derivative lookup of the 10-bit range table (halo2 LookupRangeCheckConfig).
 * @details Read R = gamma + v (active iff q_lookup, v = q_running (z - 2^K z_next) + (1 - q_running) z), write
 * T = gamma + table_idx (active iff q_table):
 *   I * R * T - exists = 0,  exists = 1 - (1 - q_lookup)(1 - q_table),
 *   sum over rows of I * (q_lookup T - counts R) = 0.
 */
template <typename Flavor_, typename FF_> class OrchardRangeLookupRelationImpl {
  public:
    using FF = FF_;
    using Cycle = typename Flavor_::Cycle;
    static constexpr std::array<size_t, 2> SUBRELATION_PARTIAL_LENGTHS{ 6, 6 };
    static constexpr std::array<bool, 2> SUBRELATION_LINEARLY_INDEPENDENT{ true, false };

    template <typename AllEntities> static bool skip(const AllEntities& in)
    {
        return detail::is_zero_value(in.q_lookup) && detail::is_zero_value(in.q_table) &&
               detail::is_zero_value(in.lookup_read_counts_range);
    }

    template <typename Accumulator, typename AllEntities, typename Parameters>
    static Accumulator read_term(const AllEntities& in, const Parameters& params)
    {
        using View = typename Accumulator::View;
        const auto z = Accumulator(View(Flavor_::advice(in)[halo2::RANGE_CHECK_COLUMN]));
        const auto z_next = Accumulator(View(Flavor_::advice_shift(in)[halo2::RANGE_CHECK_COLUMN]));
        const auto q_running = Accumulator(View(in.q_running));
        return z - q_running * z_next * FF(uint64_t{ 1 } << 10) + params.gamma;
    }

    template <typename ContainerOverSubrelations, typename AllEntities, typename Parameters>
    static void accumulate(ContainerOverSubrelations& evals,
                           const AllEntities& in,
                           const Parameters& params,
                           const FF& scaling_factor)
    {
        using Accumulator = std::tuple_element_t<0, ContainerOverSubrelations>;
        using View = typename Accumulator::View;
        const auto q_lookup = Accumulator(View(in.q_lookup));
        const auto q_table = Accumulator(View(in.q_table));
        const auto inverse = Accumulator(View(in.lookup_inverses_range));
        const auto counts = Accumulator(View(in.lookup_read_counts_range));
        const auto r = read_term<Accumulator>(in, params);
        const auto t = Accumulator(View(in.table_idx) + params.gamma);
        const auto exists = -((-q_lookup + FF(1)) * (-q_table + FF(1))) + FF(1);
        std::get<0>(evals) += (inverse * r * t - exists) * scaling_factor;
        std::get<1>(evals) += inverse * (q_lookup * t - counts * r);
    }
};

template <typename Flavor, typename FF> using OrchardGateRelation = Relation<OrchardGateRelationImpl<Flavor, FF>>;
template <typename Flavor, typename FF>
using OrchardPermutationRelation = Relation<OrchardPermutationRelationImpl<Flavor, FF>>;
template <typename Flavor, typename FF>
using OrchardSinsemillaLookupRelation = Relation<OrchardSinsemillaLookupRelationImpl<Flavor, FF>>;
template <typename Flavor, typename FF>
using OrchardRangeLookupRelation = Relation<OrchardRangeLookupRelationImpl<Flavor, FF>>;

} // namespace bb::zcash
