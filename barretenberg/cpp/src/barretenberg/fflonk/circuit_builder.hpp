#pragma once

#include "barretenberg/fflonk/polynomial_utils.hpp"

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace bb::fflonk_plonk {

/** @brief The three wire columns, in the order the permutation argument indexes them. */
static constexpr size_t NUM_WIRES = 3;

/**
 * @brief The smallest circuit this proof system accepts.
 * @details Small enough not to matter, large enough that no degree bound in the prover degenerates.
 * The Solidity verifier's constructor enforces the same bound, so a key one accepts the other does.
 */
static constexpr size_t MIN_CIRCUIT_SIZE = 8;

/**
 * @brief One circuit row, laid out over `n` rows of the domain and padded to a power of two.
 *
 * @details `wire_variables` carries, per wire and row, the *representative* variable index after
 * equality merging; the copy-constraint permutation is read off it directly. The selector and
 * witness columns are in Lagrange form over `H`.
 */
struct Trace {
    size_t circuit_size = 0;
    size_t num_public_inputs = 0;

    std::vector<FF> q_m;
    std::vector<FF> q_l;
    std::vector<FF> q_r;
    std::vector<FF> q_o;
    std::vector<FF> q_c;

    std::array<std::vector<uint32_t>, NUM_WIRES> wire_variables;
    std::array<std::vector<FF>, NUM_WIRES> wires;

    std::vector<FF> public_inputs;

    /**
     * @brief Every gate identity and every copy constraint, checked natively.
     * @details A failing trace can still be "proven" by a prover that shares the builder's mistake,
     * so this is the check that says the circuit means what it was written to mean.
     * @param failure set to a description of the first violation found.
     */
    [[nodiscard]] bool check(std::string& failure) const;
};

/**
 * @brief A plonkish circuit: `q_M·a·b + q_L·a + q_R·b + q_O·c + q_C + PI = 0` per row, plus copy
 * constraints between wire slots that hold the same variable.
 *
 * @details Variables are values, wires are slots. Two slots holding the same variable index are
 * constrained equal by the permutation argument, so a circuit is written by allocating variables and
 * referring to them from as many gates as needed. `assert_equal` merges two variables after the
 * fact, which is how a gate output gets identified with a later gate's input without an extra row.
 *
 * Public inputs are not gates the caller writes: they are materialised by `build_trace` into the
 * first `ℓ` rows, each carrying `q_L = 1` on a wire holding the public variable, so the row reads
 * `x_i + PI(ω^i) = 0`. That both binds the value into the transcript and puts the variable into the
 * copy-constraint system, so the rest of the circuit can refer to it like any other.
 */
class CircuitBuilder {
  public:
    struct Gate {
        FF q_m = FF::zero();
        FF q_l = FF::zero();
        FF q_r = FF::zero();
        FF q_o = FF::zero();
        FF q_c = FF::zero();
        uint32_t a = 0;
        uint32_t b = 0;
        uint32_t c = 0;
    };

    CircuitBuilder();

    /** @brief Allocate a variable holding `value`. */
    uint32_t add_variable(const FF& value);

    /** @brief Allocate a variable holding `value` and expose it as the next public input. */
    uint32_t add_public_input(const FF& value);

    /** @brief A variable fixed to zero, used for unconstrained wire slots and padding. */
    [[nodiscard]] uint32_t zero_variable() const { return zero_variable_; }

    /** @brief A variable fixed to one. */
    [[nodiscard]] uint32_t one_variable() const { return one_variable_; }

    /** @brief Constrain two variables to be equal by merging their copy-constraint classes. */
    void assert_equal(uint32_t lhs, uint32_t rhs);

    void create_gate(const Gate& gate);

    /** @brief `q_l·a + q_r·b + q_o·c + q_c = 0`. */
    void create_add_gate(
        const FF& q_l, uint32_t a, const FF& q_r, uint32_t b, const FF& q_o, uint32_t c, const FF& q_c);

    /** @brief Allocate and return `a * b`, constrained. */
    uint32_t create_mul(uint32_t a, uint32_t b);

    /** @brief Allocate and return `a + b`, constrained. */
    uint32_t create_add(uint32_t a, uint32_t b);

    /** @brief Constrain `variable` to equal the constant `value`. */
    void fix_variable(uint32_t variable, const FF& value);

    /** @brief Constrain `variable` to be 0 or 1. */
    void create_bool_gate(uint32_t variable);

    [[nodiscard]] const FF& get_variable(uint32_t index) const;
    [[nodiscard]] size_t num_gates() const { return gates_.size(); }
    [[nodiscard]] size_t num_variables() const { return variables_.size(); }
    [[nodiscard]] size_t num_public_inputs() const { return public_inputs_.size(); }

    /** @brief The number of rows the circuit needs, before rounding up to a power of two. */
    [[nodiscard]] size_t num_rows() const { return gates_.size() + public_inputs_.size(); }

    /**
     * @brief Materialise the padded trace.
     * @param minimum_size a lower bound on the row count, rounded up to a power of two along with
     * the circuit's own requirement. Useful for pinning a verification key to a fixed size.
     */
    [[nodiscard]] Trace build_trace(size_t minimum_size = 0) const;

  private:
    /** @brief The representative of `index`'s equality class, with path compression. */
    uint32_t find(uint32_t index) const;

    std::vector<FF> variables_;
    std::vector<Gate> gates_;
    std::vector<uint32_t> public_inputs_;
    mutable std::vector<uint32_t> equality_class_;

    uint32_t zero_variable_ = 0;
    uint32_t one_variable_ = 0;
};

} // namespace bb::fflonk_plonk
