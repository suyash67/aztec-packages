#include "barretenberg/commitment_schemes/fflonk/fflonk_honk.hpp"
#include "barretenberg/commitment_schemes/mercury/vela_honk.hpp"
#include "barretenberg/flavor/ultra_keccak_flavor.hpp"
#include "barretenberg/common/get_bytecode.hpp"
#include "barretenberg/crypto/keccak/keccak.hpp"
#include "barretenberg/dsl/acir_format/acir_format.hpp"
#include "barretenberg/dsl/acir_format/acir_to_constraint_buf.hpp"
#include "barretenberg/dsl/acir_format/serde/witness_stack.hpp"
#include "barretenberg/srs/global_crs.hpp"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// Proves an ACIR circuit with UltraHonk + Vela under a Keccak Fiat-Shamir transcript - the hash an
// EVM verifier can afford - and writes the opening argument's inputs for the Solidity gas harness in
// `recursion/sol`.
//
// The output is a flat file of 32-byte big-endian words, so the Solidity side needs no JSON parser:
//   [0] log_n | [1] num_unshifted claims | [2] num_shifted claims | [3] num_all_entities
//   [4] the transcript's previous challenge entering the opening phase
//   then: all claimed sumcheck evaluations (the round buffer the first opening challenge absorbs),
//         the claim -> entity index map, the sumcheck challenge u, the claim commitments in claim
//         order, and the opening argument itself (C_h, v0_a, v1_a, v0_b, v1_b, w0, C_q, pi_L).
//
// Usage: vela_evm_export_bench -b <bytecode> -w <witness.gz> -o <file>

namespace {

using namespace bb;

using Flavor = UltraKeccakFlavor;
using VelaHonkEvm = honk_transparent::TransparentHonk<vela::VelaPcs, Flavor, typename Flavor::Transcript>;
using FflonkHonkEvm = honk_transparent::TransparentHonk<fflonk::FflonkPcs, Flavor, typename Flavor::Transcript>;

constexpr size_t NUM_PRECOMPUTED = Flavor::NUM_PRECOMPUTED_ENTITIES;
constexpr size_t NUM_ALL_ENTITIES = Flavor::NUM_ALL_ENTITIES;
constexpr size_t BATCHED_LENGTH = Flavor::BATCHED_RELATION_PARTIAL_LENGTH;

/** @brief Verifier-side mirror of `KeccakTranscript` that also segments the proof stream. */
class KeccakMirror {
  public:
    explicit KeccakMirror(const std::vector<uint256_t>& proof)
        : proof_(proof)
    {}

    std::vector<uint256_t> receive(size_t count)
    {
        std::vector<uint256_t> values = read(count);
        buffer_.insert(buffer_.end(), values.begin(), values.end());
        return values;
    }
    uint256_t receive_one() { return receive(1)[0]; }
    void absorb(const uint256_t& value) { buffer_.push_back(value); }

    fr challenge()
    {
        std::vector<uint256_t> input;
        if (!first_) {
            input.push_back(uint256_t(previous_));
        }
        first_ = false;
        input.insert(input.end(), buffer_.begin(), buffer_.end());
        buffer_.clear();
        previous_ = crypto::Keccak::hash(input);
        return previous_;
    }

    fr previous() const { return previous_; }
    size_t cursor() const { return cursor_; }

  private:
    std::vector<uint256_t> read(size_t count)
    {
        if (cursor_ + count > proof_.size()) {
            throw_or_abort("vela_evm_export: proof stream underrun");
        }
        std::vector<uint256_t> values(proof_.begin() + static_cast<std::ptrdiff_t>(cursor_),
                                      proof_.begin() + static_cast<std::ptrdiff_t>(cursor_ + count));
        cursor_ += count;
        return values;
    }

    const std::vector<uint256_t>& proof_;
    size_t cursor_ = 0;
    std::vector<uint256_t> buffer_;
    fr previous_ = fr::zero();
    bool first_ = true;
};

void push_word(std::vector<uint256_t>& out, const uint256_t& value)
{
    out.push_back(value);
}

} // namespace

int main(int argc, char** argv)
{
    std::string bytecode_path;
    std::string witness_path;
    std::string output_path;
    for (int i = 1; i + 1 < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "-b") {
            bytecode_path = argv[++i];
        } else if (arg == "-w") {
            witness_path = argv[++i];
        } else if (arg == "-o") {
            output_path = argv[++i];
        }
    }
    if (bytecode_path.empty() || witness_path.empty() || output_path.empty()) {
        std::cerr << "usage: vela_evm_export_bench -b <bytecode> -w <witness.gz> -o <file>\n";
        return 1;
    }

    srs::init_file_crs_factory(srs::bb_crs_path());

    acir_format::AcirProgram program{
        acir_format::circuit_buf_to_acir_format(get_bytecode(bytecode_path), /*is_mega=*/false), {}
    };
    program.witness = acir_format::witness_buf_to_witness_vector(get_bytecode(witness_path));
    UltraCircuitBuilder builder = acir_format::create_circuit<UltraCircuitBuilder>(program);

    UltraCircuitBuilder sizing = builder;
    const size_t log_n = ProverInstance_<Flavor>(sizing).log_dyadic_size();
    const auto config = VelaHonkEvm::make_config(log_n);
    auto pk = VelaHonkEvm::create_proving_key(builder, config);
    const auto vk = pk.vk;
    const std::vector<uint256_t> proof = VelaHonkEvm::prove(pk);
    if (!VelaHonkEvm::verify(vk, config, proof)) {
        std::cerr << "native Vela-Honk (keccak) verification failed\n";
        return 1;
    }

    // Walk the shell to the opening phase, recording what the Solidity verifier needs.
    KeccakMirror transcript(proof);
    transcript.receive_one(); // Init
    transcript.absorb(uint256_t(log_n));
    transcript.absorb(uint256_t(vk.num_public_inputs));
    transcript.absorb(uint256_t(vk.pub_inputs_offset));
    transcript.absorb(uint256_t(vk.virtual_mask));
    transcript.absorb(uint256_t(vk.lagrange_first_row));
    transcript.absorb(uint256_t(vk.lagrange_last_row));
    for (const auto& commitment : vk.precomputed_commitment) {
        transcript.absorb(uint256_t(commitment.x));
        transcript.absorb(uint256_t(commitment.y));
    }
    transcript.receive(vk.num_public_inputs);

    std::vector<std::vector<uint256_t>> witness_commitments;
    const std::array<size_t, 3> group_columns = VelaHonkEvm::WITNESS_GROUP_COLUMNS;
    for (size_t g = 0; g < 3; ++g) {
        witness_commitments.push_back(transcript.receive(2 * group_columns[g]));
        if (g == 0) {
            transcript.challenge(); // eta
            transcript.challenge(); // rom_logup_gamma
        } else if (g == 1) {
            transcript.challenge(); // beta
            transcript.challenge(); // gamma
        }
    }
    transcript.challenge(); // alpha
    transcript.challenge(); // Sumcheck:gate_challenge

    std::vector<uint256_t> sumcheck_challenge;
    for (size_t round = 0; round < log_n; ++round) {
        transcript.receive(BATCHED_LENGTH);
        sumcheck_challenge.push_back(uint256_t(transcript.challenge()));
    }
    const std::vector<uint256_t> evaluations = transcript.receive(NUM_ALL_ENTITIES);
    const fr previous_challenge = transcript.previous();

    // Vela opening: C_h, the five transmitted values, C_q, pi_L, with the challenges between them.
    transcript.challenge(); // rho (absorbs the evaluations)
    transcript.challenge(); // lambda
    const std::vector<uint256_t> c_h = transcript.receive(2);
    transcript.challenge(); // z
    const std::vector<uint256_t> values = transcript.receive(5);
    transcript.challenge(); // alpha
    const std::vector<uint256_t> c_q = transcript.receive(2);
    transcript.challenge(); // zeta
    const std::vector<uint256_t> pi_l = transcript.receive(2);
    if (transcript.cursor() != proof.size()) {
        std::cerr << "mirror: " << proof.size() - transcript.cursor() << " unread proof fields\n";
        return 1;
    }

    // Claim order: the committed precomputed columns, then the witness columns, then the shifted
    // ones - `TransparentHonk`'s `append_unshifted_refs` / `append_shifted_refs`.
    std::vector<size_t> claim_entity;
    std::vector<std::pair<uint256_t, uint256_t>> claim_commitment;
    size_t column = 0;
    for (size_t e = 0; e < NUM_PRECOMPUTED; ++e) {
        if (((vk.virtual_mask >> e) & 1) == 0) {
            claim_entity.push_back(e);
            claim_commitment.emplace_back(uint256_t(vk.precomputed_commitment[column].x),
                                          uint256_t(vk.precomputed_commitment[column].y));
            ++column;
        }
    }
    const auto witness_point = [&](size_t group, size_t index) {
        return std::make_pair(witness_commitments[group][2 * index], witness_commitments[group][2 * index + 1]);
    };
    const std::vector<std::pair<size_t, std::pair<size_t, size_t>>> witness_layout = {
        { NUM_PRECOMPUTED + 0, { 0, 0 } }, // w_l
        { NUM_PRECOMPUTED + 1, { 0, 1 } }, // w_r
        { NUM_PRECOMPUTED + 2, { 0, 2 } }, // w_o
        { NUM_PRECOMPUTED + 3, { 1, 2 } }, // w_4
        { NUM_PRECOMPUTED + 4, { 2, 1 } }, // z_perm
        { NUM_PRECOMPUTED + 5, { 2, 0 } }, // lookup_inverses
        { NUM_PRECOMPUTED + 6, { 1, 0 } }, // lookup_read_counts
        { NUM_PRECOMPUTED + 7, { 1, 1 } }, // lookup_read_tags
    };
    for (const auto& [entity, ref] : witness_layout) {
        claim_entity.push_back(entity);
        claim_commitment.push_back(witness_point(ref.first, ref.second));
    }
    const size_t num_unshifted = claim_entity.size();
    const std::vector<std::pair<size_t, std::pair<size_t, size_t>>> shifted_layout = {
        { NUM_PRECOMPUTED + 8, { 0, 0 } },  // w_l_shift
        { NUM_PRECOMPUTED + 9, { 0, 1 } },  // w_r_shift
        { NUM_PRECOMPUTED + 10, { 0, 2 } }, // w_o_shift
        { NUM_PRECOMPUTED + 11, { 1, 2 } }, // w_4_shift
        { NUM_PRECOMPUTED + 12, { 2, 1 } }, // z_perm_shift
    };
    for (const auto& [entity, ref] : shifted_layout) {
        claim_entity.push_back(entity);
        claim_commitment.push_back(witness_point(ref.first, ref.second));
    }
    const size_t num_shifted = claim_entity.size() - num_unshifted;

    std::vector<uint256_t> out;
    push_word(out, uint256_t(log_n));
    push_word(out, uint256_t(num_unshifted));
    push_word(out, uint256_t(num_shifted));
    push_word(out, uint256_t(NUM_ALL_ENTITIES));
    push_word(out, uint256_t(previous_challenge));
    for (const uint256_t& value : evaluations) {
        push_word(out, value);
    }
    for (const size_t entity : claim_entity) {
        push_word(out, uint256_t(entity));
    }
    for (const uint256_t& value : sumcheck_challenge) {
        push_word(out, value);
    }
    for (const auto& [x, y] : claim_commitment) {
        push_word(out, x);
        push_word(out, y);
    }
    for (const uint256_t& word : c_h) {
        push_word(out, word);
    }
    for (const uint256_t& word : values) {
        push_word(out, word);
    }
    for (const uint256_t& word : c_q) {
        push_word(out, word);
    }
    for (const uint256_t& word : pi_l) {
        push_word(out, word);
    }

    const auto write_words = [](const std::string& path, const std::vector<uint256_t>& words) {
        std::ofstream stream(path, std::ios::binary);
        for (const uint256_t& word : words) {
            for (size_t limb = 4; limb-- > 0;) {
                for (size_t byte = 8; byte-- > 0;) {
                    const char value = static_cast<char>((word.data[limb] >> (8 * byte)) & 0xff);
                    stream.write(&value, 1);
                }
            }
        }
    };
    // The proof itself, so the harness can price its calldata against the KZG proof's.
    write_words(output_path + ".proof", proof);

    std::ofstream stream(output_path, std::ios::binary);
    for (const uint256_t& word : out) {
        for (size_t limb = 4; limb-- > 0;) {
            for (size_t byte = 8; byte-- > 0;) {
                const char value = static_cast<char>((word.data[limb] >> (8 * byte)) & 0xff);
                stream.write(&value, 1);
            }
        }
    }
    std::cerr << "vela-evm: log_n " << log_n << ", claims " << num_unshifted << " + " << num_shifted << ", proof "
              << proof.size() << " words (" << proof.size() * 32 / 1024 << " KiB), export " << out.size() << " words\n";
    return 0;
}
