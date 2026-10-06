#pragma once

#include "barretenberg/stdlib/primitives/bool/bool.hpp"
#include "barretenberg/stdlib/primitives/field/field.hpp"
#include "barretenberg/stdlib/primitives/group/cycle_group.hpp"
#include "barretenberg/stdlib/primitives/memory/twin_rom_table.hpp"
#include "barretenberg/stdlib_circuit_builders/ultra_circuit_builder.hpp"
#include "barretenberg/zcash/halo2/layout.hpp"
#include "barretenberg/zcash/primitives/orchard.hpp"

#include <optional>
#include <vector>

/**
 * @file action_circuit_ultra.hpp
 * @brief The Orchard Action statement over BN254 / Grumpkin as an UltraHonk circuit built from barretenberg's stdlib.
 *
 * @details The statement is identical to Orchard<Bn254Cycle> (the native reference), i.e. the Orchard Action with the
 * Pallas primitives replaced by their Grumpkin counterparts: Sinsemilla over Grumpkin with the same domain strings,
 * Poseidon P128Pow5T3 over F_r (parameters from the Poseidon reference generator), Grumpkin fixed bases. The gadgets
 * are barretenberg's:
 *   - Sinsemilla S-table: a twin ROM table of the 1024 (x, y) pairs; every hash step is two `cycle_group` incomplete
 *     additions (Ultra elliptic gates), as in halo2's Sinsemilla chip.
 *   - message decomposition: 10-bit words are slices of the message fields, range-constrained with Ultra's range
 *     constraints; full field elements in NoteCommit / CommitIvk are decomposed canonically
 * (`validate_split_in_field`).
 *   - fixed-base scalar multiplications: `cycle_group::fixed_batch_mul` (ROM-based Straus tables);
 *     variable-base: `cycle_group::batch_mul`.
 *   - Poseidon: the P128Pow5T3 permutation with `field_t` arithmetic.
 * Public inputs (10 per Action, in instance order) are the circuit's public inputs.
 */
namespace bb::zcash::ultra {

using Builder = UltraCircuitBuilder;
using field_ct = stdlib::field_t<Builder>;
using bool_ct = stdlib::bool_t<Builder>;
using cycle_group_ct = stdlib::cycle_group<Builder>;
using cycle_scalar_ct = stdlib::cycle_scalar<Builder>;
using twin_rom_table_ct = stdlib::twin_rom_table<Builder>;

class ActionCircuitUltra {
  public:
    using Cycle = Bn254Cycle;
    using O = Orchard<Cycle>;
    using FF = Cycle::FF;
    using Scalar = Cycle::EmbeddedScalar;
    using AffineElement = Cycle::AffineElement;
    using Witness = O::ActionWitness;

    /**
     * @brief Adds every Action to the builder; the public inputs are added in order (10 per Action).
     */
    static void build(Builder& builder, const std::vector<Witness>& witnesses, const std::vector<FF>& public_inputs);

    explicit ActionCircuitUltra(Builder& builder);

    void synthesize_action(const Witness& w, std::span<const FF> public_inputs);

    // A field of a Sinsemilla message: `num_bits` low bits of `value`; `canonical` requires value < r as an integer
    // (for full field elements encoded on 255 bits).
    struct Segment {
        field_ct value;
        size_t num_bits;
        bool canonical = false;
    };
    std::vector<field_ct> message_words(const std::vector<Segment>& segments);
    cycle_group_ct hash_to_point(const AffineElement& q, const std::vector<field_ct>& words);
    field_ct poseidon_hash(const field_ct& a, const field_ct& b);
    cycle_group_ct fixed_base_mul(const AffineElement& base, const cycle_scalar_ct& scalar);
    cycle_scalar_ct base_field_scalar(const field_ct& x);
    // (least-significant bit of y, canonically)
    field_ct y_lsb(const field_ct& y);
    cycle_group_ct note_commit(const cycle_group_ct& g_d,
                               const cycle_group_ct& pk_d,
                               const field_ct& value,
                               const field_ct& rho,
                               const field_ct& psi,
                               const cycle_scalar_ct& rcm);
    field_ct commit_ivk(const field_ct& ak, const field_ct& nk, const cycle_scalar_ct& rivk);
    field_ct merkle_root(const field_ct& leaf, const std::array<FF, O::MERKLE_DEPTH>& path, uint32_t pos);

  private:
    Builder& builder_;
    twin_rom_table_ct s_table_;
};

} // namespace bb::zcash::ultra
