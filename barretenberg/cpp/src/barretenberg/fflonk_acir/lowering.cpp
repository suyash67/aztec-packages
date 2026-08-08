#include "barretenberg/fflonk_acir/lowering.hpp"

#include "barretenberg/common/throw_or_abort.hpp"
#include "barretenberg/relations/poseidon2_internal_relation.hpp"

#include <unordered_map>

namespace bb::fflonk_acir {

namespace {

using FF = fflonk_plonk::FF;
using Gate = fflonk_plonk::CircuitBuilder::Gate;

/** @brief Gate selectors that have no three-wire encoding worth having. */
constexpr std::array<GateKind, 3> UNSUPPORTED_GATES = { GateKind::Lookup, GateKind::Elliptic, GateKind::Nnf };

const char* name_of(const GateKind kind)
{
    switch (kind) {
    case GateKind::Lookup:
        return "lookup";
    case GateKind::Elliptic:
        return "elliptic curve addition";
    case GateKind::Nnf:
        return "non-native field";
    case GateKind::Memory:
        return "RAM/ROM";
    default:
        return "unsupported";
    }
}

/**
 * @brief Rewrites an Ultra trace row by row into three-wire gates.
 *
 * @details Every Ultra gate spans four wires and, for several gate types, reaches into the next row
 * as well. Three wires hold none of those directly, so each is re-expressed as a chain of three-wire
 * rows threaded by fresh intermediate variables - which is sound because each link is itself a
 * constraint, so the chain enforces exactly the composite the Ultra relation did.
 */
class Lowering {
  public:
    explicit Lowering(UltraCircuitBuilder& circuit)
        : circuit_(circuit)
    {}

    fflonk_plonk::CircuitBuilder run();

  private:
    /**
     * @brief The three-wire variable standing for an Ultra one.
     * @details Keyed on the representative after Ultra's equality merging, so two slots that Ultra's
     * copy constraints identify land on one variable here and need no constraint of their own.
     */
    uint32_t translate(uint32_t ultra_index)
    {
        const uint32_t representative = circuit_.real_variable_index[ultra_index];
        const auto found = translated_.find(representative);
        if (found != translated_.end()) {
            return found->second;
        }
        const uint32_t created = out_.add_variable(circuit_.get_variable(ultra_index));
        translated_.emplace(representative, created);
        return created;
    }

    [[nodiscard]] const FF& value(const uint32_t index) const { return out_.get_variable(index); }

    /** @brief `q_m a b + q_l a + q_r b + q_o c + q_c = 0`. */
    void gate(const FF& q_m,
              const FF& q_l,
              const FF& q_r,
              const FF& q_o,
              const FF& q_c,
              const uint32_t a,
              const uint32_t b,
              const uint32_t c)
    {
        out_.create_gate(Gate{ .q_m = q_m, .q_l = q_l, .q_r = q_r, .q_o = q_o, .q_c = q_c, .a = a, .b = b, .c = c });
    }

    /** @brief A fresh variable constrained to equal `q_l a + q_r b + q_c`. */
    uint32_t linear(const FF& q_l, const uint32_t a, const FF& q_r, const uint32_t b, const FF& q_c)
    {
        const uint32_t result = out_.add_variable((q_l * value(a)) + (q_r * value(b)) + q_c);
        gate(FF::zero(), q_l, q_r, -FF::one(), q_c, a, b, result);
        return result;
    }

    /** @brief A fresh variable constrained to equal `q_m a b + q_l a + q_r b`. */
    uint32_t linear_from_product(const FF& q_m, const uint32_t a, const uint32_t b, const FF& q_l, const FF& q_r)
    {
        const uint32_t result = out_.add_variable((q_m * value(a) * value(b)) + (q_l * value(a)) + (q_r * value(b)));
        gate(q_m, q_l, q_r, -FF::one(), FF::zero(), a, b, result);
        return result;
    }

    /** @brief A fresh variable constrained to equal `a * b`. */
    uint32_t product(const uint32_t a, const uint32_t b)
    {
        const uint32_t result = out_.add_variable(value(a) * value(b));
        gate(FF::one(), FF::zero(), FF::zero(), -FF::one(), FF::zero(), a, b, result);
        return result;
    }

    /** @brief A fresh variable constrained to equal `x^5`, the Poseidon2 S-box. */
    uint32_t pow5(const uint32_t x)
    {
        const uint32_t squared = product(x, x);
        return product(product(squared, squared), x);
    }

    /** @brief The wire of the row after `idx`, refusing to guess across a block boundary. */
    uint32_t next_row_wire(auto& block, size_t idx, size_t wire, const char* gate_name);

    void expand_arithmetic(auto& block, size_t idx, const FF& q_arith);
    void expand_poseidon2_external(auto& block, size_t idx);
    void expand_poseidon2_internal(auto& block, size_t idx);

    UltraCircuitBuilder& circuit_;
    fflonk_plonk::CircuitBuilder out_;
    std::unordered_map<uint32_t, uint32_t> translated_;
};

uint32_t Lowering::next_row_wire(auto& block, const size_t idx, const size_t wire, const char* gate_name)
{
    // Rows within a block are contiguous in the final trace, but the row after a block's last one
    // belongs to whichever block the trace layout puts next - which is not this pass's to know, so a
    // gate in that position is refused rather than guessed at.
    if (idx + 1 >= block.size()) {
        throw_or_abort(std::string("fflonk: a ") + gate_name +
                       " gate reads the next trace row across a block boundary");
    }
    return translate(block.wires[wire][idx + 1]);
}

/**
 * @brief `q_m w_l w_r + q_l w_l + q_r w_r + q_o w_o + q_4 w_4 + q_c [+ w_4(next)] = 0`.
 *
 * @details Five wire slots, or six at `q_arith == 2`, against three - so two or three rows. The
 * product term's coefficient is `(q_arith - 3) q_m (-1/2)`: `q_m` at mode 1, half of it at mode 2.
 */
void Lowering::expand_arithmetic(auto& block, const size_t idx, const FF& q_arith)
{
    const bool uses_next_row = q_arith == FF(2);

    const FF q_m = block.q_m()[idx];
    const FF q_l = block.q_1()[idx];
    const FF q_r = block.q_2()[idx];
    const FF q_o = block.q_3()[idx];
    const FF q_4 = block.q_4()[idx];
    const FF q_c = block.q_c()[idx];

    const uint32_t a = translate(block.w_l()[idx]);
    const uint32_t b = translate(block.w_r()[idx]);
    const uint32_t c = translate(block.w_o()[idx]);
    const uint32_t d = translate(block.w_4()[idx]);

    const FF product_coefficient = uses_next_row ? q_m * FF(2).invert() : q_m;
    const uint32_t partial = linear_from_product(product_coefficient, a, b, q_l, q_r);

    if (!uses_next_row) {
        gate(FF::zero(), FF::one(), q_o, q_4, q_c, partial, c, d);
        return;
    }

    const uint32_t next = next_row_wire(block, idx, 3, "arithmetic");
    const uint32_t running = linear(FF::one(), partial, q_o, c, FF::zero());
    gate(FF::zero(), FF::one(), q_4, FF::one(), q_c, running, d, next);
}

/**
 * @brief The full external Poseidon2 round: an S-box on each lane, then the MDS matrix.
 *
 * @details The matrix rows are the ones barretenberg's relation computes with its fourteen
 * additions; written out they are
 * `v1 = 5u1+7u2+u3+3u4`, `v2 = 4u1+6u2+u3+u4`, `v3 = u1+3u2+5u3+7u4`, `v4 = u1+u2+4u3+6u4`.
 */
void Lowering::expand_poseidon2_external(auto& block, const size_t idx)
{
    const std::array<FF, 4> constants = { block.q_1()[idx], block.q_2()[idx], block.q_3()[idx], block.q_4()[idx] };

    std::array<uint32_t, 4> u{};
    for (size_t lane = 0; lane < 4; ++lane) {
        const uint32_t wire = translate(block.wires[lane][idx]);
        u[lane] = pow5(linear(FF::one(), wire, FF::zero(), out_.zero_variable(), constants[lane]));
    }

    static constexpr std::array<std::array<uint64_t, 4>, 4> MATRIX = {
        { { 5, 7, 1, 3 }, { 4, 6, 1, 1 }, { 1, 3, 5, 7 }, { 1, 1, 4, 6 } }
    };

    for (size_t row = 0; row < 4; ++row) {
        const uint32_t next = next_row_wire(block, idx, row, "Poseidon2");
        const uint32_t partial = linear(FF(MATRIX[row][0]), u[0], FF(MATRIX[row][1]), u[1], FF::zero());
        const uint32_t running = linear(FF::one(), partial, FF(MATRIX[row][2]), u[2], FF::zero());
        gate(FF::zero(), FF::one(), FF(MATRIX[row][3]), -FF::one(), FF::zero(), running, u[3], next);
    }
}

/**
 * @brief The internal Poseidon2 round: one S-box, then the diagonal matrix.
 *
 * @details `v_i = u_i + sum_j u_j` with `u_j = w_j` for `j > 1`, so every row shares the partial sum
 * `w_2 + w_3 + w_4`. The diagonal constants come from the relation itself rather than being copied,
 * so a change to the Poseidon2 parameters cannot leave this pass behind.
 */
void Lowering::expand_poseidon2_internal(auto& block, const size_t idx)
{
    using Internal = Poseidon2InternalRelationImpl<FF>;

    const uint32_t w_1 = translate(block.w_l()[idx]);
    const uint32_t w_2 = translate(block.w_r()[idx]);
    const uint32_t w_3 = translate(block.w_o()[idx]);
    const uint32_t w_4 = translate(block.w_4()[idx]);

    const uint32_t u_1 = pow5(linear(FF::one(), w_1, FF::zero(), out_.zero_variable(), block.q_1()[idx]));
    const uint32_t partial_sum =
        linear(FF::one(), linear(FF::one(), w_2, FF::one(), w_3, FF::zero()), FF::one(), w_4, FF::zero());

    // v_1 = D_1 u_1 + (w_2 + w_3 + w_4)
    gate(FF::zero(),
         Internal::D1,
         FF::one(),
         -FF::one(),
         FF::zero(),
         u_1,
         partial_sum,
         next_row_wire(block, idx, 0, "Poseidon2"));

    // v_i = u_1 + (D_i - 1) w_i + (w_2 + w_3 + w_4), for i = 2, 3, 4
    const uint32_t base = linear(FF::one(), u_1, FF::one(), partial_sum, FF::zero());
    const std::array<FF, 3> diagonal = { Internal::D2_minus_1, Internal::D3_minus_1, Internal::D4_minus_1 };
    const std::array<uint32_t, 3> wires = { w_2, w_3, w_4 };
    for (size_t lane = 0; lane < 3; ++lane) {
        gate(FF::zero(),
             FF::one(),
             diagonal[lane],
             -FF::one(),
             FF::zero(),
             base,
             wires[lane],
             next_row_wire(block, idx, lane + 1, "Poseidon2"));
    }
}

fflonk_plonk::CircuitBuilder Lowering::run()
{
    circuit_.finalize_circuit();

    for (const uint32_t public_input : circuit_.public_inputs()) {
        // `add_public_input` allocates its own variable on a public-input row; tying it to the
        // circuit's variable is what makes the claimed value the one the gates actually used.
        const uint32_t exposed = out_.add_public_input(circuit_.get_variable(public_input));
        out_.assert_equal(exposed, translate(public_input));
    }

    for (auto& block : circuit_.blocks.get()) {
        for (size_t idx = 0; idx < block.size(); ++idx) {
            for (const GateKind kind : UNSUPPORTED_GATES) {
                if (!read_gate_selector(block, kind, idx).is_zero()) {
                    throw_or_abort(std::string("fflonk: this circuit uses a ") + name_of(kind) +
                                   " gate, which the three-wire arithmetization does not have. Prove it with "
                                   "--scheme ultra_fflonk instead.");
                }
            }
            if (!read_gate_selector(block, GateKind::Memory, idx).is_zero()) {
                throw_or_abort("fflonk: this circuit indexes an array dynamically, which needs the RAM/ROM argument "
                               "the three-wire arithmetization does not have. Prove it with --scheme ultra_fflonk "
                               "instead.");
            }

            const FF q_arith = read_gate_selector(block, GateKind::Arith, idx);
            if (!q_arith.is_zero()) {
                if (q_arith != FF::one() && q_arith != FF(2)) {
                    throw_or_abort("fflonk: unsupported arithmetic gate mode");
                }
                expand_arithmetic(block, idx, q_arith);
            }
            if (!read_gate_selector(block, GateKind::DeltaRange, idx).is_zero()) {
                // The delta-range quartic itself expands into three wires without trouble; it is only
                // half of Ultra's range argument. The other half is that the sorted
                // list the deltas run over is a *permutation of the constrained variables*, and Ultra
                // enforces that with `real_variable_tags` and `tau`: a multiset equality layered onto
                // sigma, separate from the copy constraints in `real_variable_index`.
                //
                // This pass carries `real_variable_index` and nothing else, so lowering only the
                // quartic would leave a prover free to put any sorted list it liked in those wires.
                // `ABrokenRangeConstraintSurvivesLoweringAsAFailingTrace` is the test that caught it.
                // Carrying the tags needs the three-wire builder to express a permutation cycle that
                // is not an equality class, which its union-find cannot - so this is an
                // arithmetization change to `fflonk/`, not a gap in this file.
                throw_or_abort("fflonk: this circuit uses a range constraint, whose soundness rests on Ultra's tag "
                               "multiset argument. The three-wire arithmetization has no mechanism for it. Prove it "
                               "with --scheme ultra_fflonk instead.");
            }
            if (!read_gate_selector(block, GateKind::Poseidon2Ext, idx).is_zero()) {
                expand_poseidon2_external(block, idx);
            }
            if (!read_gate_selector(block, GateKind::Poseidon2Int, idx).is_zero()) {
                expand_poseidon2_internal(block, idx);
            }
        }
    }

    return std::move(out_);
}

} // namespace

fflonk_plonk::CircuitBuilder lower(UltraCircuitBuilder& circuit)
{
    return Lowering(circuit).run();
}

} // namespace bb::fflonk_acir
