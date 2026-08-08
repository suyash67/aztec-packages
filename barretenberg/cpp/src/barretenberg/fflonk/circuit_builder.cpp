#include "barretenberg/fflonk/circuit_builder.hpp"

#include "barretenberg/common/log.hpp"
#include "barretenberg/common/throw_or_abort.hpp"

#include <algorithm>

namespace bb::fflonk_plonk {

namespace {
size_t round_up_to_power_of_two(size_t value)
{
    size_t power = 1;
    while (power < value) {
        power <<= 1;
    }
    return power;
}
} // namespace

CircuitBuilder::CircuitBuilder()
{
    // Both constants are constrained by a gate rather than merely stored, so a circuit that uses
    // them as operands is relying on something the proof actually enforces.
    zero_variable_ = add_variable(FF::zero());
    fix_variable(zero_variable_, FF::zero());
    one_variable_ = add_variable(FF::one());
    fix_variable(one_variable_, FF::one());
}

uint32_t CircuitBuilder::add_variable(const FF& value)
{
    const auto index = static_cast<uint32_t>(variables_.size());
    variables_.push_back(value);
    equality_class_.push_back(index);
    return index;
}

uint32_t CircuitBuilder::add_public_input(const FF& value)
{
    const uint32_t index = add_variable(value);
    public_inputs_.push_back(index);
    return index;
}

uint32_t CircuitBuilder::find(uint32_t index) const
{
    BB_ASSERT_LT(static_cast<size_t>(index), equality_class_.size(), "unknown variable index");
    uint32_t root = index;
    while (equality_class_[root] != root) {
        root = equality_class_[root];
    }
    while (equality_class_[index] != root) {
        const uint32_t next = equality_class_[index];
        equality_class_[index] = root;
        index = next;
    }
    return root;
}

void CircuitBuilder::assert_equal(uint32_t lhs, uint32_t rhs)
{
    const uint32_t left_root = find(lhs);
    const uint32_t right_root = find(rhs);
    if (left_root == right_root) {
        return;
    }
    BB_ASSERT(variables_[left_root] == variables_[right_root],
              "assert_equal on variables holding different values: the circuit is unsatisfiable");

    // Lowest index wins, so the representative of a class does not depend on the order the merges
    // happened to be issued in - the permutation, and therefore the verification key, stays stable.
    const uint32_t root = std::min(left_root, right_root);
    const uint32_t child = std::max(left_root, right_root);
    equality_class_[child] = root;
}

void CircuitBuilder::create_gate(const Gate& gate)
{
    BB_ASSERT_LT(static_cast<size_t>(gate.a), variables_.size(), "gate references an unknown variable");
    BB_ASSERT_LT(static_cast<size_t>(gate.b), variables_.size(), "gate references an unknown variable");
    BB_ASSERT_LT(static_cast<size_t>(gate.c), variables_.size(), "gate references an unknown variable");
    gates_.push_back(gate);
}

void CircuitBuilder::create_add_gate(
    const FF& q_l, uint32_t a, const FF& q_r, uint32_t b, const FF& q_o, uint32_t c, const FF& q_c)
{
    create_gate({ .q_m = FF::zero(), .q_l = q_l, .q_r = q_r, .q_o = q_o, .q_c = q_c, .a = a, .b = b, .c = c });
}

uint32_t CircuitBuilder::create_mul(uint32_t a, uint32_t b)
{
    const uint32_t out = add_variable(get_variable(a) * get_variable(b));
    create_gate({ .q_m = FF::one(),
                  .q_l = FF::zero(),
                  .q_r = FF::zero(),
                  .q_o = -FF::one(),
                  .q_c = FF::zero(),
                  .a = a,
                  .b = b,
                  .c = out });
    return out;
}

uint32_t CircuitBuilder::create_add(uint32_t a, uint32_t b)
{
    const uint32_t out = add_variable(get_variable(a) + get_variable(b));
    create_add_gate(FF::one(), a, FF::one(), b, -FF::one(), out, FF::zero());
    return out;
}

void CircuitBuilder::fix_variable(uint32_t variable, const FF& value)
{
    create_add_gate(FF::one(), variable, FF::zero(), zero_variable_, FF::zero(), zero_variable_, -value);
}

void CircuitBuilder::create_bool_gate(uint32_t variable)
{
    // v*v - v = 0
    create_gate({ .q_m = FF::one(),
                  .q_l = -FF::one(),
                  .q_r = FF::zero(),
                  .q_o = FF::zero(),
                  .q_c = FF::zero(),
                  .a = variable,
                  .b = variable,
                  .c = zero_variable_ });
}

const FF& CircuitBuilder::get_variable(uint32_t index) const
{
    return variables_[find(index)];
}

Trace CircuitBuilder::build_trace(const size_t minimum_size) const
{
    const size_t num_public = public_inputs_.size();
    const size_t rows = num_public + gates_.size();
    const size_t circuit_size = round_up_to_power_of_two(std::max({ rows, minimum_size, MIN_CIRCUIT_SIZE }));

    if (rows > circuit_size) {
        throw_or_abort("fflonk: circuit does not fit the requested size");
    }

    Trace trace;
    trace.circuit_size = circuit_size;
    trace.num_public_inputs = num_public;
    trace.q_m.assign(circuit_size, FF::zero());
    trace.q_l.assign(circuit_size, FF::zero());
    trace.q_r.assign(circuit_size, FF::zero());
    trace.q_o.assign(circuit_size, FF::zero());
    trace.q_c.assign(circuit_size, FF::zero());
    for (size_t wire = 0; wire < NUM_WIRES; ++wire) {
        trace.wire_variables[wire].assign(circuit_size, find(zero_variable_));
        trace.wires[wire].assign(circuit_size, FF::zero());
    }

    const uint32_t zero_root = find(zero_variable_);
    const FF zero_value = variables_[zero_root];
    for (size_t wire = 0; wire < NUM_WIRES; ++wire) {
        std::fill(trace.wires[wire].begin(), trace.wires[wire].end(), zero_value);
    }

    auto place = [&](size_t row, size_t wire, uint32_t variable) {
        const uint32_t root = find(variable);
        trace.wire_variables[wire][row] = root;
        trace.wires[wire][row] = variables_[root];
    };

    // Public inputs first, one row each: q_L = 1 on wire a, so the row reads x_i + PI(w^i) = 0.
    for (size_t i = 0; i < num_public; ++i) {
        trace.q_l[i] = FF::one();
        place(i, 0, public_inputs_[i]);
        trace.public_inputs.push_back(variables_[find(public_inputs_[i])]);
    }

    for (size_t g = 0; g < gates_.size(); ++g) {
        const size_t row = num_public + g;
        const Gate& gate = gates_[g];
        trace.q_m[row] = gate.q_m;
        trace.q_l[row] = gate.q_l;
        trace.q_r[row] = gate.q_r;
        trace.q_o[row] = gate.q_o;
        trace.q_c[row] = gate.q_c;
        place(row, 0, gate.a);
        place(row, 1, gate.b);
        place(row, 2, gate.c);
    }

    return trace;
}

bool Trace::check(std::string& failure) const
{
    for (size_t row = 0; row < circuit_size; ++row) {
        const FF a = wires[0][row];
        const FF b = wires[1][row];
        const FF c = wires[2][row];
        const FF public_input = row < num_public_inputs ? -public_inputs[row] : FF::zero();
        const FF value = q_m[row] * a * b + q_l[row] * a + q_r[row] * b + q_o[row] * c + q_c[row] + public_input;
        if (!value.is_zero()) {
            failure = "gate identity violated at row " + std::to_string(row);
            return false;
        }
    }

    // Copy constraints: every slot holding the same variable must hold the same value.
    std::vector<std::pair<uint32_t, FF>> seen;
    seen.reserve(circuit_size * NUM_WIRES);
    for (size_t wire = 0; wire < NUM_WIRES; ++wire) {
        for (size_t row = 0; row < circuit_size; ++row) {
            seen.emplace_back(wire_variables[wire][row], wires[wire][row]);
        }
    }
    std::ranges::sort(seen, [](const auto& lhs, const auto& rhs) { return lhs.first < rhs.first; });
    for (size_t i = 1; i < seen.size(); ++i) {
        if (seen[i].first == seen[i - 1].first && seen[i].second != seen[i - 1].second) {
            failure = "copy constraint violated for variable " + std::to_string(seen[i].first);
            return false;
        }
    }

    return true;
}

} // namespace bb::fflonk_plonk
