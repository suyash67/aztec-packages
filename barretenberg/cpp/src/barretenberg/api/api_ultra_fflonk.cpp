#include "barretenberg/api/api_ultra_fflonk.hpp"

#include "barretenberg/api/file_io.hpp"
#include "barretenberg/common/get_bytecode.hpp"
#include "barretenberg/common/throw_or_abort.hpp"
#include "barretenberg/dsl/acir_format/acir_to_constraint_buf.hpp"
#include "barretenberg/dsl/acir_format/serde/witness_stack.hpp"
#include "barretenberg/dsl/acir_proofs/ultra_fflonk_contract.hpp"
#include "barretenberg/ultra_fflonk/prover.hpp"
#include "barretenberg/ultra_fflonk/verifier.hpp"

#include <iostream>

namespace bb {

namespace {

using Builder = ultra_fflonk::Flavor::CircuitBuilder;

/**
 * @brief Build the circuit from ACIR, through the frontend the Honk path uses.
 * @details `witness` may be empty, which yields a circuit suitable for preprocessing but not proving.
 */
Builder build_circuit(std::vector<uint8_t>&& bytecode, std::vector<uint8_t>&& witness)
{
    acir_format::AcirProgram program{ acir_format::circuit_buf_to_acir_format(std::move(bytecode),
                                                                              /*is_mega=*/false),
                                      {} };
    if (!witness.empty()) {
        program.witness = acir_format::witness_buf_to_witness_vector(std::move(witness));
    }
    return acir_format::create_circuit<Builder>(program);
}

/** @brief A full-width `0x`-prefixed literal, which is what a `uint256` constructor argument wants. */
std::string to_hex(const uint256_t& value)
{
    static constexpr std::string_view DIGITS = "0123456789abcdef";
    std::string out = "0x";
    for (size_t limb = 4; limb-- > 0;) {
        for (size_t nibble = 16; nibble-- > 0;) {
            out += DIGITS[(value.data[limb] >> (4 * nibble)) & 0xf];
        }
    }
    return out;
}

std::vector<uint8_t> public_inputs_to_buffer(const std::vector<ultra_fflonk::FF>& public_inputs)
{
    std::vector<uint256_t> words;
    words.reserve(public_inputs.size());
    for (const ultra_fflonk::FF& public_input : public_inputs) {
        words.push_back(static_cast<uint256_t>(public_input));
    }
    return to_buffer(words);
}

} // namespace

bool UltraFflonkAPI::check(const Flags&,
                           const std::filesystem::path& bytecode_path,
                           const std::filesystem::path& witness_path)
{
    // Building the prover instance runs every witness computation the proof depends on, and
    // preprocessing rejects a trace this proof system cannot accept, so a circuit that survives both
    // is one `prove` will not abort on for structural reasons.
    Builder builder = build_circuit(get_bytecode(bytecode_path), get_bytecode(witness_path));
    ultra_fflonk::ProvingKey key = ultra_fflonk::preprocess(builder);
    info("Circuit accepted: ", key.circuit_size(), " rows, ", key.public_inputs().size(), " public inputs");
    return true;
}

void UltraFflonkAPI::prove(const Flags& flags,
                           const std::filesystem::path& bytecode_path,
                           const std::filesystem::path& witness_path,
                           const std::filesystem::path& vk_path,
                           const std::filesystem::path& output_dir)
{
    if (output_dir == "-") {
        throw_or_abort("Stdout output is not supported. Please specify an output directory.");
    }

    Builder builder = build_circuit(get_bytecode(bytecode_path), get_bytecode(witness_path));
    ultra_fflonk::ProvingKey key = ultra_fflonk::preprocess(builder);

    // A supplied key must be the one this circuit preprocesses to; proving against a different one
    // would produce a proof nothing accepts, which is a worse failure than saying so here.
    if (!vk_path.empty() && !flags.write_vk) {
        ultra_fflonk::VerificationKey supplied;
        if (!ultra_fflonk::VerificationKey::from_buffer(read_file(vk_path), supplied)) {
            throw_or_abort("ultra_fflonk: could not parse the supplied verification key");
        }
        if (!(supplied == key.verification_key)) {
            throw_or_abort("ultra_fflonk: the supplied verification key does not match this circuit");
        }
    }

    const ultra_fflonk::Proof proof = ultra_fflonk::prove(key);

    write_file(output_dir / "proof", proof.to_buffer());
    write_file(output_dir / "public_inputs", public_inputs_to_buffer(key.public_inputs()));
    info("Proof saved to ", output_dir / "proof");
    info("Public inputs saved to ", output_dir / "public_inputs");

    if (flags.write_vk) {
        write_file(output_dir / "vk", key.verification_key.to_buffer());
        info("VK saved to ", output_dir / "vk");
    }
}

bool UltraFflonkAPI::verify(const Flags&,
                            const std::filesystem::path& public_inputs_path,
                            const std::filesystem::path& proof_path,
                            const std::filesystem::path& vk_path)
{
    ultra_fflonk::VerificationKey key;
    if (!ultra_fflonk::VerificationKey::from_buffer(read_file(vk_path), key)) {
        info("Proof verification failed: the verification key is malformed");
        return false;
    }

    const std::vector<uint256_t> words =
        many_from_buffer_exact<uint256_t>(read_file(public_inputs_path), "ultra_fflonk public inputs file");
    std::vector<ultra_fflonk::FF> public_inputs;
    public_inputs.reserve(words.size());
    for (const uint256_t& word : words) {
        if (word >= ultra_fflonk::FF::modulus) {
            info("Proof verification failed: a public input is not a canonical field element");
            return false;
        }
        public_inputs.emplace_back(word);
    }

    return ultra_fflonk::verify(key, read_file(proof_path), public_inputs);
}

void UltraFflonkAPI::write_vk(const Flags&,
                              const std::filesystem::path& bytecode_path,
                              const std::filesystem::path& output_dir)
{
    if (output_dir == "-") {
        throw_or_abort("Stdout output is not supported. Please specify an output directory.");
    }

    Builder builder = build_circuit(get_bytecode(bytecode_path), {});
    ultra_fflonk::ProvingKey key = ultra_fflonk::preprocess(builder);

    write_file(output_dir / "vk", key.verification_key.to_buffer());
    info("VK saved to ", output_dir / "vk");
}

void UltraFflonkAPI::gates(const Flags& flags, const std::filesystem::path& bytecode_path)
{
    acir_format::AcirProgram program{ acir_format::circuit_buf_to_acir_format(get_bytecode(bytecode_path),
                                                                              /*is_mega=*/false),
                                      {} };
    const size_t num_acir_opcodes = program.constraints.num_acir_opcodes;

    acir_format::ProgramMetadata metadata;
    metadata.collect_gates_per_opcode = flags.include_gates_per_opcode;
    Builder builder = acir_format::create_circuit<Builder>(program, metadata);
    builder.finalize_circuit();

    const size_t num_gates = builder.get_finalized_total_circuit_size();

    std::string gates_per_opcode;
    if (flags.include_gates_per_opcode) {
        for (size_t i = 0; i < program.constraints.gates_per_opcode.size(); ++i) {
            if (i != 0) {
                gates_per_opcode += ",";
            }
            gates_per_opcode += std::to_string(program.constraints.gates_per_opcode[i]);
        }
    }

    std::cout << format(
        "{\"functions\": [\n  {\n        \"acir_opcodes\": ",
        num_acir_opcodes,
        ",\n        \"circuit_size\": ",
        num_gates,
        (flags.include_gates_per_opcode ? format(",\n        \"gates_per_opcode\": [", gates_per_opcode, "]") : ""),
        "\n  }\n]}");
}

void UltraFflonkAPI::write_solidity_verifier(const Flags&,
                                             const std::filesystem::path& output_path,
                                             const std::filesystem::path& vk_path)
{
    ultra_fflonk::VerificationKey key;
    if (!ultra_fflonk::VerificationKey::from_buffer(read_vk_file(vk_path), key)) {
        throw_or_abort("ultra_fflonk: could not parse the verification key");
    }

    // The verifier takes its key as constructor arguments, so one contract source serves every
    // circuit. Emitting a derived contract that supplies them keeps the audited source verbatim and
    // still hands the user a file they can deploy without knowing what the arguments mean.
    std::string contract = std::string(ULTRA_FFLONK_VERIFIER_SOURCE);
    contract += "\n/// @notice `UltraFflonkVerifier` with this circuit's verification key baked in.\n";
    contract += "contract UltraVerifier is UltraFflonkVerifier {\n";
    contract += "    constructor()\n        UltraFflonkVerifier(\n";
    contract += "            " + std::to_string(key.circuit_size) + ",\n";
    contract += "            " + std::to_string(key.num_public_inputs) + ",\n";
    contract += "            " + std::to_string(key.pub_inputs_offset) + ",\n";
    contract += "            " + to_hex(static_cast<uint256_t>(key.omega)) + ",\n";
    contract += "            [\n";
    for (size_t g = 0; g < ultra_fflonk::NUM_PREPROCESSED_GROUPS; ++g) {
        const std::string prefix = g == 0 ? "                uint256(" : "                ";
        const std::string suffix = g == 0 ? ")" : "";
        contract += prefix + to_hex(static_cast<uint256_t>(key.preprocessed[g].x)) + suffix + ",\n";
        contract += "                " + to_hex(static_cast<uint256_t>(key.preprocessed[g].y)) +
                    (g + 1 == ultra_fflonk::NUM_PREPROCESSED_GROUPS ? "\n" : ",\n");
    }
    contract += "            ]\n        )\n    {}\n}\n";

    if (output_path == "-") {
        std::cout << contract;
    } else {
        write_file(output_path,
                   std::span<const uint8_t>(reinterpret_cast<const uint8_t*>(contract.c_str()), contract.size()));
        info("ultra_fflonk solidity verifier saved to ", output_path);
    }
}

} // namespace bb
