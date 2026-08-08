#include "barretenberg/commitment_schemes/fflonk/fflonk_honk.hpp"
#include "barretenberg/common/get_bytecode.hpp"
#include "barretenberg/crypto/keccak/keccak.hpp"
#include "barretenberg/dsl/acir_format/acir_format.hpp"
#include "barretenberg/dsl/acir_format/acir_to_constraint_buf.hpp"
#include "barretenberg/dsl/acir_format/serde/witness_stack.hpp"
#include "barretenberg/flavor/ultra_keccak_flavor.hpp"
#include "barretenberg/srs/global_crs.hpp"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// Proves an ACIR circuit with UltraHonk + the fflonk opening under a Keccak transcript and writes the
// opening argument for the Solidity gas harness in `recursion/sol`, in the same flat 32-byte-word
// form `vela_evm_export` uses.
//
// Usage: fflonk_evm_export_bench -b <bytecode> -w <witness.gz> -o <file>

namespace {

using namespace bb;

using Flavor = UltraKeccakFlavor;
using FflonkHonkEvm = honk_transparent::TransparentHonk<fflonk::FflonkPcs, Flavor, typename Flavor::Transcript>;

constexpr size_t NUM_PRECOMPUTED = Flavor::NUM_PRECOMPUTED_ENTITIES;
constexpr size_t NUM_ALL_ENTITIES = Flavor::NUM_ALL_ENTITIES;
constexpr size_t BATCHED_LENGTH = Flavor::BATCHED_RELATION_PARTIAL_LENGTH;

/** @brief Verifier-side mirror of `KeccakTranscript` that also segments the proof stream. */
class KeccakMirror {
  public:
    explicit KeccakMirror(const std::vector<uint256_t>& proof)
        : proof_(&proof)
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
        if (cursor_ + count > proof_->size()) {
            throw_or_abort("fflonk_evm_export: proof stream underrun");
        }
        std::vector<uint256_t> values(proof_->begin() + static_cast<std::ptrdiff_t>(cursor_),
                                      proof_->begin() + static_cast<std::ptrdiff_t>(cursor_ + count));
        cursor_ += count;
        return values;
    }

    const std::vector<uint256_t>* proof_;
    size_t cursor_ = 0;
    std::vector<uint256_t> buffer_;
    fr previous_ = fr::zero();
    bool first_ = true;
};

void write_words(const std::string& path, const std::vector<uint256_t>& words)
{
    std::ofstream stream(path, std::ios::binary);
    for (const uint256_t& word : words) {
        for (size_t limb = 4; limb-- > 0;) {
            for (size_t byte = 8; byte-- > 0;) {
                const char value = static_cast<char>((word.data[limb] >> (8 * byte)) & 0xff);
                stream.write(&value, 1);
            }
        }
    }
}

size_t next_power_of_two(size_t value)
{
    size_t power = 1;
    while (power < value) {
        power <<= 1;
    }
    return power;
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
        std::cerr << "usage: fflonk_evm_export_bench -b <bytecode> -w <witness.gz> -o <file>\n";
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
    const auto config = FflonkHonkEvm::make_config(log_n);
    auto pk = FflonkHonkEvm::create_proving_key(builder, config);
    const auto vk = pk.vk;
    const std::vector<uint256_t> proof = FflonkHonkEvm::prove(pk);
    if (!FflonkHonkEvm::verify(vk, config, proof)) {
        std::cerr << "native fflonk-Honk (keccak) verification failed\n";
        return 1;
    }

    // Walk the shell to the opening phase.
    KeccakMirror transcript(proof);
    transcript.receive_one(); // Init
    transcript.absorb(uint256_t(log_n));
    transcript.absorb(uint256_t(vk.num_public_inputs));
    transcript.absorb(uint256_t(vk.pub_inputs_offset));
    transcript.absorb(uint256_t(vk.virtual_mask));
    transcript.absorb(uint256_t(vk.lagrange_first_row));
    transcript.absorb(uint256_t(vk.lagrange_last_row));
    transcript.absorb(uint256_t(vk.precomputed_commitment.x));
    transcript.absorb(uint256_t(vk.precomputed_commitment.y));
    transcript.receive(vk.num_public_inputs);

    std::vector<std::pair<uint256_t, uint256_t>> group_commitments;
    group_commitments.emplace_back(uint256_t(vk.precomputed_commitment.x), uint256_t(vk.precomputed_commitment.y));
    for (size_t g = 0; g < 3; ++g) {
        const std::vector<uint256_t> point = transcript.receive(2);
        group_commitments.emplace_back(point[0], point[1]);
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

    // Claim layout: the committed precomputed columns, the witness columns, then the shifted ones.
    struct Claim {
        size_t group;
        size_t column;
        size_t entity;
    };
    std::vector<Claim> claims;
    std::vector<size_t> group_columns;
    size_t column = 0;
    for (size_t e = 0; e < NUM_PRECOMPUTED; ++e) {
        if (((vk.virtual_mask >> e) & 1) == 0) {
            claims.push_back({ 0, column++, e });
        }
    }
    group_columns = { column, 3, 3, 2 };
    const std::vector<Claim> witness_layout = {
        { 1, 0, NUM_PRECOMPUTED + 0 }, { 1, 1, NUM_PRECOMPUTED + 1 }, { 1, 2, NUM_PRECOMPUTED + 2 },
        { 2, 2, NUM_PRECOMPUTED + 3 }, { 3, 1, NUM_PRECOMPUTED + 4 }, { 3, 0, NUM_PRECOMPUTED + 5 },
        { 2, 0, NUM_PRECOMPUTED + 6 }, { 2, 1, NUM_PRECOMPUTED + 7 },
    };
    claims.insert(claims.end(), witness_layout.begin(), witness_layout.end());
    const size_t num_unshifted = claims.size();
    const std::vector<Claim> shifted_layout = {
        { 1, 0, NUM_PRECOMPUTED + 8 },  { 1, 1, NUM_PRECOMPUTED + 9 },  { 1, 2, NUM_PRECOMPUTED + 10 },
        { 2, 2, NUM_PRECOMPUTED + 11 }, { 3, 1, NUM_PRECOMPUTED + 12 },
    };
    claims.insert(claims.end(), shifted_layout.begin(), shifted_layout.end());
    const size_t num_shifted = claims.size() - num_unshifted;

    std::vector<size_t> packs;
    for (const size_t columns : group_columns) {
        packs.push_back(next_power_of_two(columns));
    }
    packs.push_back(1); // h

    // fflonk opening.
    transcript.challenge(); // rho (absorbs the evaluations)
    transcript.challenge(); // lambda
    const std::vector<uint256_t> c_h = transcript.receive(2);
    transcript.challenge(); // z
    std::vector<uint256_t> evaluation_pairs;
    for (size_t g = 0; g + 1 < packs.size(); ++g) {
        for (size_t c = 0; c < packs[g]; ++c) {
            const std::vector<uint256_t> pair = transcript.receive(2);
            evaluation_pairs.push_back(pair[0]);
            evaluation_pairs.push_back(pair[1]);
        }
    }
    const std::vector<uint256_t> h_evals = transcript.receive(2);
    evaluation_pairs.push_back(h_evals[0]);
    evaluation_pairs.push_back(h_evals[1]);
    transcript.challenge(); // nu
    const std::vector<uint256_t> c_w = transcript.receive(2);
    transcript.challenge(); // y
    const std::vector<uint256_t> c_wp = transcript.receive(2);
    if (transcript.cursor() != proof.size()) {
        std::cerr << "mirror: " << proof.size() - transcript.cursor() << " unread proof fields\n";
        return 1;
    }

    std::vector<uint256_t> out;
    out.emplace_back(log_n);
    out.emplace_back(group_columns.size());
    out.emplace_back(NUM_ALL_ENTITIES);
    out.emplace_back(num_unshifted);
    out.emplace_back(num_shifted);
    out.emplace_back(uint256_t(previous_challenge));
    for (const size_t pack : packs) {
        out.emplace_back(pack);
    }
    out.insert(out.end(), evaluations.begin(), evaluations.end());
    for (const Claim& claim : claims) {
        out.emplace_back(claim.group);
    }
    for (const Claim& claim : claims) {
        out.emplace_back(claim.column);
    }
    for (const Claim& claim : claims) {
        out.emplace_back(claim.entity);
    }
    out.insert(out.end(), sumcheck_challenge.begin(), sumcheck_challenge.end());
    for (const auto& [x, y] : group_commitments) {
        out.push_back(x);
        out.push_back(y);
    }
    out.insert(out.end(), c_h.begin(), c_h.end());
    out.insert(out.end(), evaluation_pairs.begin(), evaluation_pairs.end());
    out.insert(out.end(), c_w.begin(), c_w.end());
    out.insert(out.end(), c_wp.begin(), c_wp.end());

    write_words(output_path, out);
    write_words(output_path + ".proof", proof);
    std::cerr << "fflonk-evm: log_n " << log_n << ", claims " << num_unshifted << " + " << num_shifted << ", packs";
    for (const size_t pack : packs) {
        std::cerr << " " << pack;
    }
    std::cerr << ", proof " << proof.size() << " words (" << proof.size() * 32 / 1024 << " KiB), export " << out.size()
              << " words\n";
    return 0;
}
