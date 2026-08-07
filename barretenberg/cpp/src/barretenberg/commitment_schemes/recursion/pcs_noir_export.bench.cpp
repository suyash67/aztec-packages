#include "barretenberg/commitment_schemes/ligero/ligero_honk.hpp"
#include "barretenberg/commitment_schemes/recursion/noir_export.hpp"
#include "barretenberg/commitment_schemes/whir/whir_honk.hpp"
#include "barretenberg/flavor/ultra_provekit_flavor.hpp"
#include "barretenberg/honk/library/grand_product_delta.hpp"
#include "barretenberg/stdlib_circuit_builders/mock_circuits.hpp"
#include "barretenberg/ultra_honk/prover_instance.hpp"

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <vector>

// Exports a WHIR-Honk or Ligero-Honk proof in the structured form a Noir recursive verifier
// consumes: `params.nr` (every compile-time constant of the proof family) and `Prover.toml` (this
// proof's field elements). The verifier's read order is reproduced here by a mirror of
// `TransparentHonk::verify` that segments the flat proof as it walks the Fiat-Shamir chain; the
// mirror re-checks the proof from the segmented data alone, so a successful export is evidence that
// the exported data is exactly what a verifier needs and nothing more.
//
// Usage: pcs_noir_export_bench --pcs <whir|ligero> --log-n <m> --out <dir>
//        [--security <bits>] [--rate <log_inv_rate>] [--k <bits>] [--k0 <bits>] [--soundness <s>]
//
// Both backends run the Poseidon2 Merkle hasher: an Fr digest and an algebraic compression are the
// only sensible choices inside a BN254 circuit.

namespace {

using namespace bb;
using namespace bb::pcs_recursion;

using Hasher = whir::Poseidon2MerkleHasher;
using Flavor = UltraProveKitFlavor;
using FF = fr;

constexpr size_t NUM_PRECOMPUTED = Flavor::NUM_PRECOMPUTED_ENTITIES;       // 24
constexpr size_t NUM_UNSHIFTED = Flavor::NUM_UNSHIFTED_ENTITIES;           // 32
constexpr size_t NUM_ALL_ENTITIES = Flavor::NUM_ALL_ENTITIES;              // 37
constexpr size_t NUM_SUBRELATIONS = Flavor::NUM_SUBRELATIONS;              // 20
constexpr size_t BATCHED_LENGTH = Flavor::BATCHED_RELATION_PARTIAL_LENGTH; // 7

// ---------------------------------------------------------------------------------------------
// A Fiat-Shamir transcript that also segments the proof stream.
// ---------------------------------------------------------------------------------------------

/**
 * @brief Mirror of `bb::NativeTranscript`'s verifier side that records what it reads.
 * @details Reproduces `BaseTranscript`'s duplex exactly: a challenge is
 * `Poseidon2(previous_challenge ++ round_buffer)`, the buffer holding everything received through
 * `receive` since the last challenge. `receive_unhashed` reads without absorbing - those values are
 * bound to the transcript by a Merkle root instead, and re-absorbing them would cost the in-circuit
 * verifier a permutation per three field elements for no soundness gain.
 */
class MirrorTranscript {
  public:
    explicit MirrorTranscript(const HonkProof& proof)
        : proof_(proof)
    {}

    std::vector<FF> receive(size_t count)
    {
        std::vector<FF> values = read(count);
        buffer_.insert(buffer_.end(), values.begin(), values.end());
        return values;
    }

    FF receive_one() { return receive(1)[0]; }

    std::vector<FF> receive_unhashed(size_t count) { return read(count); }

    void absorb(const FF& value) { buffer_.push_back(value); }
    void absorb(std::span<const FF> values) { buffer_.insert(buffer_.end(), values.begin(), values.end()); }

    FF challenge()
    {
        std::vector<FF> input;
        input.reserve(buffer_.size() + 1);
        if (!first_) {
            input.push_back(previous_);
        }
        first_ = false;
        input.insert(input.end(), buffer_.begin(), buffer_.end());
        buffer_.clear();
        previous_ = crypto::Poseidon2<crypto::Poseidon2Bn254ScalarFieldParams>::hash(input);
        return previous_;
    }

    size_t cursor() const { return cursor_; }
    bool exhausted() const { return cursor_ == proof_.size(); }

  private:
    std::vector<FF> read(size_t count)
    {
        if (cursor_ + count > proof_.size()) {
            throw_or_abort("pcs_noir_export: proof stream underrun");
        }
        std::vector<FF> values(proof_.begin() + static_cast<std::ptrdiff_t>(cursor_),
                               proof_.begin() + static_cast<std::ptrdiff_t>(cursor_ + count));
        cursor_ += count;
        return values;
    }

    const HonkProof& proof_;
    size_t cursor_ = 0;
    std::vector<FF> buffer_;
    FF previous_ = FF::zero();
    bool first_ = true;
};

// ---------------------------------------------------------------------------------------------
// Claim layout shared by both backends (mirrors TransparentHonk's private helpers).
// ---------------------------------------------------------------------------------------------

struct ColumnRef {
    size_t group;
    size_t column;
};

struct ClaimLayout {
    std::vector<ColumnRef> unshifted;
    std::vector<size_t> unshifted_entity; // index into get_all()
    std::vector<ColumnRef> shifted;
    std::vector<size_t> shifted_entity;
    std::vector<size_t> group_columns;
    std::vector<bool> is_virtual;
    std::vector<size_t> precomputed_dense_column; // entity -> column within group 0
};

ClaimLayout build_claim_layout(uint32_t virtual_mask, size_t num_committed_precomputed)
{
    ClaimLayout layout;
    layout.is_virtual.resize(NUM_PRECOMPUTED);
    layout.precomputed_dense_column.assign(NUM_PRECOMPUTED, 0);
    size_t column = 0;
    for (size_t e = 0; e < NUM_PRECOMPUTED; ++e) {
        const bool hidden = ((virtual_mask >> e) & 1) == 1;
        layout.is_virtual[e] = hidden;
        if (!hidden) {
            layout.precomputed_dense_column[e] = column;
            layout.unshifted.push_back({ 0, column++ });
            layout.unshifted_entity.push_back(e);
        }
    }
    BB_ASSERT_EQ(column, num_committed_precomputed);
    // Witness entities in get_unshifted() order after the precomputed block, with the (group,
    // column) each was committed at: wires | counts_w4 | inverses_z_perm.
    const std::vector<std::pair<size_t, ColumnRef>> witness = {
        { 24, { 1, 0 } }, // w_l
        { 25, { 1, 1 } }, // w_r
        { 26, { 1, 2 } }, // w_o
        { 27, { 2, 2 } }, // w_4
        { 28, { 3, 1 } }, // z_perm
        { 29, { 3, 0 } }, // lookup_inverses
        { 30, { 2, 0 } }, // lookup_read_counts
        { 31, { 2, 1 } }, // lookup_read_tags
    };
    for (const auto& [entity, ref] : witness) {
        layout.unshifted.push_back(ref);
        layout.unshifted_entity.push_back(entity);
    }
    layout.shifted = { { 1, 0 }, { 1, 1 }, { 1, 2 }, { 2, 2 }, { 3, 1 } };
    layout.shifted_entity = { 32, 33, 34, 35, 36 };
    layout.group_columns = { num_committed_precomputed, 3, 3, 2 };
    return layout;
}

// ---------------------------------------------------------------------------------------------
// Honk shell: everything the two backends share, up to the opening argument.
// ---------------------------------------------------------------------------------------------

struct ShellData {
    size_t log_n = 0;
    size_t num_public_inputs = 0;
    size_t pub_inputs_offset = 0;
    uint32_t virtual_mask = 0;
    size_t lagrange_first_row = 0;
    size_t lagrange_last_row = 0;
    FF vk_root;

    FF init;
    std::vector<FF> public_inputs;
    std::array<FF, 3> witness_roots{};                       // wires, counts_w4, inverses_z_perm
    std::vector<std::array<FF, BATCHED_LENGTH>> univariates; // one per sumcheck round
    std::array<FF, NUM_ALL_ENTITIES> evaluations{};

    // Derived (recomputed by the Noir circuit; kept here for the golden trace).
    std::vector<FF> challenge;
    FF alpha;
};

/** @brief Walk the shell's Fiat-Shamir rounds, filling `ShellData` and leaving the transcript at
 * the start of the opening argument. */
ShellData read_shell(MirrorTranscript& transcript,
                     size_t log_n,
                     size_t num_public_inputs,
                     size_t pub_inputs_offset,
                     uint32_t virtual_mask,
                     size_t lagrange_first_row,
                     size_t lagrange_last_row,
                     const FF& vk_root)
{
    ShellData shell;
    shell.log_n = log_n;
    shell.num_public_inputs = num_public_inputs;
    shell.pub_inputs_offset = pub_inputs_offset;
    shell.virtual_mask = virtual_mask;
    shell.lagrange_first_row = lagrange_first_row;
    shell.lagrange_last_row = lagrange_last_row;
    shell.vk_root = vk_root;

    shell.init = transcript.receive_one();
    // absorb_vk
    transcript.absorb(FF(log_n));
    transcript.absorb(FF(num_public_inputs));
    transcript.absorb(FF(pub_inputs_offset));
    transcript.absorb(FF(virtual_mask));
    transcript.absorb(FF(lagrange_first_row));
    transcript.absorb(FF(lagrange_last_row));
    transcript.absorb(vk_root);

    shell.public_inputs = transcript.receive(num_public_inputs);
    shell.witness_roots[0] = transcript.receive_one(); // HONK:wires
    transcript.challenge();                            // eta
    transcript.challenge();                            // rom_logup_gamma
    shell.witness_roots[1] = transcript.receive_one(); // HONK:counts_w4
    transcript.challenge();                            // beta
    transcript.challenge();                            // gamma
    shell.witness_roots[2] = transcript.receive_one(); // HONK:inverses_z_perm
    shell.alpha = transcript.challenge();              // alpha
    transcript.challenge();                            // Sumcheck:gate_challenge

    shell.univariates.resize(log_n);
    for (size_t round = 0; round < log_n; ++round) {
        const std::vector<FF> values = transcript.receive(BATCHED_LENGTH);
        std::copy(values.begin(), values.end(), shell.univariates[round].begin());
        shell.challenge.push_back(transcript.challenge()); // Sumcheck:u_i
    }
    const std::vector<FF> evaluations = transcript.receive(NUM_ALL_ENTITIES);
    std::copy(evaluations.begin(), evaluations.end(), shell.evaluations.begin());
    return shell;
}

void emit_shell_params(NoirWriter& params, const ShellData& shell, const ClaimLayout& layout)
{
    params.comment("Honk shell (UltraProveKitFlavor)");
    params.global_u32("LOG_N", shell.log_n);
    params.global_u32("NUM_PUBLIC_INPUTS", shell.num_public_inputs);
    params.global_u32("PUB_INPUTS_OFFSET", shell.pub_inputs_offset);
    params.global_field("VK_VIRTUAL_MASK", FF(shell.virtual_mask));
    params.global_u32("LAGRANGE_FIRST_ROW", shell.lagrange_first_row);
    params.global_u32("LAGRANGE_LAST_ROW", shell.lagrange_last_row);
    params.global_u32("NUM_PRECOMPUTED", NUM_PRECOMPUTED);
    params.global_u32("NUM_UNSHIFTED", NUM_UNSHIFTED);
    params.global_u32("NUM_ALL_ENTITIES", NUM_ALL_ENTITIES);
    params.global_u32("NUM_SUBRELATIONS", NUM_SUBRELATIONS);
    params.global_u32("BATCHED_LENGTH", BATCHED_LENGTH);
    params.global_bool_array("IS_VIRTUAL", layout.is_virtual);
    params.blank();

    params.comment("Batched-opening claim layout");
    params.global_u32("NUM_GROUPS", layout.group_columns.size());
    params.global_u32_array("GROUP_COLUMNS", layout.group_columns);
    std::vector<size_t> group_column_offsets;
    size_t total_columns = 0;
    for (const size_t columns : layout.group_columns) {
        group_column_offsets.push_back(total_columns);
        total_columns += columns;
    }
    params.global_u32_array("GROUP_COLUMN_OFFSET", group_column_offsets);
    params.global_u32("TOTAL_COLUMNS", total_columns);
    params.global_u32("NUM_UNSHIFTED_CLAIMS", layout.unshifted.size());
    std::vector<size_t> groups;
    std::vector<size_t> columns;
    for (const ColumnRef& ref : layout.unshifted) {
        groups.push_back(ref.group);
        columns.push_back(ref.column);
    }
    params.global_u32_array("UNSHIFTED_GROUP", groups);
    params.global_u32_array("UNSHIFTED_COLUMN", columns);
    params.global_u32_array("UNSHIFTED_ENTITY", layout.unshifted_entity);
    groups.clear();
    columns.clear();
    for (const ColumnRef& ref : layout.shifted) {
        groups.push_back(ref.group);
        columns.push_back(ref.column);
    }
    params.global_u32("NUM_SHIFTED_CLAIMS", layout.shifted.size());
    params.global_u32_array("SHIFTED_GROUP", groups);
    params.global_u32_array("SHIFTED_COLUMN", columns);
    params.global_u32_array("SHIFTED_ENTITY", layout.shifted_entity);
    params.blank();
}

void emit_shell_toml(TomlWriter& toml, const ShellData& shell)
{
    toml.field("init", shell.init);
    toml.field("vk_root", shell.vk_root);
    toml.field_array("public_inputs", shell.public_inputs);
    toml.field_array("witness_roots", { shell.witness_roots[0], shell.witness_roots[1], shell.witness_roots[2] });
    std::vector<FF> flat;
    for (const auto& univariate : shell.univariates) {
        flat.insert(flat.end(), univariate.begin(), univariate.end());
    }
    toml.field_array_2d("sumcheck_univariates", flat, BATCHED_LENGTH);
    toml.field_array("sumcheck_evaluations", std::vector<FF>(shell.evaluations.begin(), shell.evaluations.end()));
}

// ---------------------------------------------------------------------------------------------
// Merkle batch bookkeeping, mirrored from MerkleTree so the hints match the prover's walk.
// ---------------------------------------------------------------------------------------------

struct BatchPlan {
    std::vector<size_t> leaves;    // distinct, ascending
    std::vector<size_t> positions; // per query: index into `leaves`
    size_t num_siblings = 0;
};

BatchPlan plan_batch(const std::vector<size_t>& indices, size_t depth)
{
    BatchPlan plan;
    plan.leaves = whir::MerkleTree<Hasher>::batch_leaves(indices);
    for (const size_t index : indices) {
        const auto it = std::lower_bound(plan.leaves.begin(), plan.leaves.end(), index);
        plan.positions.push_back(static_cast<size_t>(it - plan.leaves.begin()));
    }
    plan.num_siblings = whir::MerkleTree<Hasher>::batch_num_siblings(plan.leaves, depth);
    return plan;
}

/** @brief One authenticated batch, padded to the worst case the circuit is compiled for. */
struct BatchOpening {
    BatchPlan plan;
    std::vector<FF> values;   // padded: max_leaves * leaf_width
    std::vector<FF> siblings; // padded: max_leaves * depth
};

BatchOpening read_batch(MirrorTranscript& transcript,
                        const std::vector<size_t>& indices,
                        size_t leaf_width,
                        size_t depth,
                        size_t max_leaves)
{
    BatchOpening opening;
    opening.plan = plan_batch(indices, depth);
    const size_t num_leaves = opening.plan.leaves.size();
    const std::vector<FF> flat = transcript.receive_unhashed(num_leaves * leaf_width + opening.plan.num_siblings);
    opening.values.assign(flat.begin(), flat.begin() + static_cast<std::ptrdiff_t>(num_leaves * leaf_width));
    opening.values.resize(max_leaves * leaf_width, FF::zero());
    opening.siblings.assign(flat.begin() + static_cast<std::ptrdiff_t>(num_leaves * leaf_width), flat.end());
    opening.siblings.resize(max_leaves * depth, FF::zero());
    return opening;
}

/** @brief Recompute a batch's root exactly as `MerkleTree::verify_batch` does. */
FF verify_batch_root(const BatchOpening& opening, size_t leaf_width, size_t depth)
{
    const std::vector<size_t>& leaves = opening.plan.leaves;
    std::vector<size_t> level(leaves.begin(), leaves.end());
    std::vector<FF> digests(leaves.size());
    for (size_t i = 0; i < leaves.size(); ++i) {
        std::vector<FF> input;
        input.push_back(FF::zero());
        input.insert(input.end(),
                     opening.values.begin() + static_cast<std::ptrdiff_t>(i * leaf_width),
                     opening.values.begin() + static_cast<std::ptrdiff_t>((i + 1) * leaf_width));
        digests[i] = crypto::Poseidon2<crypto::Poseidon2Bn254ScalarFieldParams>::hash(input);
    }
    size_t cursor = 0;
    for (size_t l = 0; l < depth; ++l) {
        std::vector<size_t> parents;
        std::vector<FF> parent_digests;
        for (size_t p = 0; p < level.size();) {
            const bool paired = p + 1 < level.size() && level[p + 1] == (level[p] ^ 1);
            parents.push_back(level[p] >> 1);
            if (paired) {
                parent_digests.push_back(Hasher::hash_node(digests[p], digests[p + 1]));
                p += 2;
                continue;
            }
            const FF sibling = opening.siblings[cursor++];
            parent_digests.push_back((level[p] & 1) ? Hasher::hash_node(sibling, digests[p])
                                                    : Hasher::hash_node(digests[p], sibling));
            p += 1;
        }
        level = std::move(parents);
        digests = std::move(parent_digests);
    }
    return digests[0];
}

} // namespace

#include "pcs_noir_export_ligero.inc"
#include "pcs_noir_export_whir.inc"

namespace {

std::string arg_value(int argc, char** argv, const std::string& name, const std::string& fallback)
{
    for (int i = 1; i + 1 < argc; ++i) {
        if (name == argv[i]) {
            return argv[i + 1];
        }
    }
    return fallback;
}

size_t arg_size(int argc, char** argv, const std::string& name, size_t fallback)
{
    const std::string value = arg_value(argc, argv, name, "");
    return value.empty() ? fallback : static_cast<size_t>(std::stoul(value));
}

/**
 * @brief The circuit whose proof is recursively verified: arithmetic, public inputs, lookups and a
 * small ROM table, padded with arithmetic gates to a chosen dyadic size.
 * @details Every relation of `UltraProveKitFlavor` is exercised, so the recursive verifier really
 * evaluates all twenty subrelations rather than a degenerate subset.
 */
UltraCircuitBuilder build_inner_circuit(size_t target_log_n)
{
    UltraCircuitBuilder builder;
    MockCircuits::add_arithmetic_gates_with_public_inputs(builder, 8);
    if (target_log_n >= 14) {
        MockCircuits::add_lookup_gates(builder, 2);
    } else {
        // UINT32_XOR's 4096-row table does not fit a small trace; the 8-row dummy multitable keeps
        // the lookup argument active without dictating the circuit size.
        const FF a_value(1);
        const auto a_idx = builder.add_variable(a_value);
        const auto accumulators =
            plookup::get_lookup_accumulators(plookup::MultiTableId::HONK_DUMMY_MULTI, a_value, FF(0), true);
        builder.create_gates_from_plookup_accumulators(
            plookup::MultiTableId::HONK_DUMMY_MULTI, accumulators, a_idx, std::nullopt);
    }
    const size_t rom_id = builder.create_ROM_array(4);
    for (size_t i = 0; i < 4; ++i) {
        builder.set_ROM_element(rom_id, i, builder.add_variable(FF(3 * i + 1)));
    }
    const uint32_t read_index = builder.add_variable(FF(2));
    builder.read_ROM_array(rom_id, read_index);

    // Pad to the requested size. `ProverInstance` rounds the finalized trace up to a power of two,
    // so grow in blocks and stop as soon as the target is reached.
    for (size_t attempt = 0; attempt < 64; ++attempt) {
        UltraCircuitBuilder probe = builder;
        ProverInstance_<Flavor> instance(probe);
        if (instance.log_dyadic_size() >= target_log_n) {
            return builder;
        }
        const size_t deficit = (size_t(1) << target_log_n) - instance.dyadic_size();
        MockCircuits::add_arithmetic_gates(builder, std::max<size_t>(deficit / 2, 16));
    }
    return builder;
}

} // namespace

int main(int argc, char** argv)
{
    const std::string pcs = arg_value(argc, argv, "--pcs", "");
    const std::string out = arg_value(argc, argv, "--out", "");
    const size_t log_n = arg_size(argc, argv, "--log-n", 12);
    const size_t security_bits = arg_size(argc, argv, "--security", 100);
    const size_t log_inv_rate = arg_size(argc, argv, "--rate", 2);
    if (pcs.empty() || out.empty()) {
        std::cerr << "usage: pcs_noir_export_bench --pcs <whir|ligero> --out <dir> [--log-n m] "
                     "[--security bits] [--rate r] [--k bits] [--k0 bits] [--soundness s]\n";
        return 1;
    }

    UltraCircuitBuilder builder = build_inner_circuit(log_n);
    UltraCircuitBuilder sizing = builder;
    ProverInstance_<Flavor> sizing_instance(sizing);
    const size_t actual_log_n = sizing_instance.log_dyadic_size();
    std::cerr << "inner circuit: log_n = " << actual_log_n
              << ", public inputs = " << sizing_instance.num_public_inputs() << "\n";

    if (pcs == "ligero") {
        return export_ligero(builder, actual_log_n, security_bits, log_inv_rate, out);
    }
    if (pcs == "whir") {
        return export_whir(builder,
                           actual_log_n,
                           security_bits,
                           log_inv_rate,
                           arg_size(argc, argv, "--k", 4),
                           arg_size(argc, argv, "--k0", 0),
                           arg_value(argc, argv, "--soundness", "provable"),
                           out);
    }
    std::cerr << "unknown pcs: " << pcs << "\n";
    return 1;
}
