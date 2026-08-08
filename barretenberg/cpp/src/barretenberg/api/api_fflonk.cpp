#include "barretenberg/api/api_fflonk.hpp"

#include "barretenberg/api/file_io.hpp"
#include "barretenberg/common/get_bytecode.hpp"
#include "barretenberg/common/throw_or_abort.hpp"
#include "barretenberg/dsl/acir_format/acir_to_constraint_buf.hpp"
#include "barretenberg/dsl/acir_format/serde/witness_stack.hpp"
#include "barretenberg/dsl/acir_proofs/fflonk_contract.hpp"
#include "barretenberg/fflonk/prover.hpp"
#include "barretenberg/fflonk/verifier.hpp"
#include "barretenberg/fflonk_acir/lowering.hpp"

#include <iostream>

namespace bb {

namespace {

using fflonk_plonk::FF;

/**
 * @brief Build the ACIR circuit, then rewrite it into the three-wire arithmetization.
 * @details The trace is checked before it is proven, so a lowering bug surfaces as a named circuit
 * violation rather than as a proof that fails with no explanation.
 */
fflonk_plonk::CircuitBuilder build_circuit(std::vector<uint8_t>&& bytecode, std::vector<uint8_t>&& witness)
{
    acir_format::AcirProgram program{ acir_format::circuit_buf_to_acir_format(std::move(bytecode), /*is_mega=*/false),
                                      {} };
    const bool has_witness = !witness.empty();
    if (has_witness) {
        program.witness = acir_format::witness_buf_to_witness_vector(std::move(witness));
    }

    UltraCircuitBuilder ultra = acir_format::create_circuit<UltraCircuitBuilder>(program);
    fflonk_plonk::CircuitBuilder lowered = fflonk_acir::lower(ultra);

    if (has_witness) {
        std::string failure;
        if (!lowered.build_trace().check(failure)) {
            throw_or_abort("fflonk: the lowered circuit is not satisfied: " + failure);
        }
    }
    return lowered;
}

std::vector<uint8_t> public_inputs_to_buffer(const std::vector<FF>& public_inputs)
{
    std::vector<uint256_t> words;
    words.reserve(public_inputs.size());
    for (const FF& public_input : public_inputs) {
        words.push_back(static_cast<uint256_t>(public_input));
    }
    return to_buffer(words);
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

} // namespace

bool FflonkAPI::check(const Flags&,
                      const std::filesystem::path& bytecode_path,
                      const std::filesystem::path& witness_path)
{
    const fflonk_plonk::CircuitBuilder circuit = build_circuit(get_bytecode(bytecode_path), get_bytecode(witness_path));
    info("Circuit accepted: ",
         circuit.num_rows(),
         " rows before padding, ",
         circuit.num_public_inputs(),
         " public inputs");
    return true;
}

void FflonkAPI::prove(const Flags& flags,
                      const std::filesystem::path& bytecode_path,
                      const std::filesystem::path& witness_path,
                      const std::filesystem::path& vk_path,
                      const std::filesystem::path& output_dir)
{
    if (output_dir == "-") {
        throw_or_abort("Stdout output is not supported. Please specify an output directory.");
    }

    const fflonk_plonk::CircuitBuilder circuit = build_circuit(get_bytecode(bytecode_path), get_bytecode(witness_path));
    const fflonk_plonk::ProvingKey key = fflonk_plonk::preprocess(circuit);

    if (!vk_path.empty() && !flags.write_vk) {
        fflonk_plonk::VerificationKey supplied;
        if (!fflonk_plonk::VerificationKey::from_buffer(read_file(vk_path), supplied)) {
            throw_or_abort("fflonk: could not parse the supplied verification key");
        }
        if (!(supplied == key.verification_key)) {
            throw_or_abort("fflonk: the supplied verification key does not match this circuit");
        }
    }

    const fflonk_plonk::Proof proof = fflonk_plonk::prove(key);

    write_file(output_dir / "proof", proof.to_buffer());
    write_file(output_dir / "public_inputs", public_inputs_to_buffer(key.trace.public_inputs));
    info("Proof saved to ", output_dir / "proof");
    info("Public inputs saved to ", output_dir / "public_inputs");

    if (flags.write_vk) {
        write_file(output_dir / "vk", key.verification_key.to_buffer());
        info("VK saved to ", output_dir / "vk");
    }
}

bool FflonkAPI::verify(const Flags&,
                       const std::filesystem::path& public_inputs_path,
                       const std::filesystem::path& proof_path,
                       const std::filesystem::path& vk_path)
{
    fflonk_plonk::VerificationKey key;
    if (!fflonk_plonk::VerificationKey::from_buffer(read_file(vk_path), key)) {
        info("Proof verification failed: the verification key is malformed");
        return false;
    }

    const std::vector<uint256_t> words =
        many_from_buffer_exact<uint256_t>(read_file(public_inputs_path), "fflonk public inputs file");
    std::vector<FF> public_inputs;
    public_inputs.reserve(words.size());
    for (const uint256_t& word : words) {
        if (word >= FF::modulus) {
            info("Proof verification failed: a public input is not a canonical field element");
            return false;
        }
        public_inputs.emplace_back(word);
    }

    return fflonk_plonk::verify(key, read_file(proof_path), public_inputs);
}

void FflonkAPI::write_vk(const Flags&,
                         const std::filesystem::path& bytecode_path,
                         const std::filesystem::path& output_dir)
{
    if (output_dir == "-") {
        throw_or_abort("Stdout output is not supported. Please specify an output directory.");
    }

    const fflonk_plonk::CircuitBuilder circuit = build_circuit(get_bytecode(bytecode_path), {});
    const fflonk_plonk::ProvingKey key = fflonk_plonk::preprocess(circuit);

    write_file(output_dir / "vk", key.verification_key.to_buffer());
    info("VK saved to ", output_dir / "vk");
}

void FflonkAPI::gates(const Flags&, const std::filesystem::path& bytecode_path)
{
    acir_format::AcirProgram program{ acir_format::circuit_buf_to_acir_format(get_bytecode(bytecode_path),
                                                                              /*is_mega=*/false),
                                      {} };
    const size_t num_acir_opcodes = program.constraints.num_acir_opcodes;

    UltraCircuitBuilder ultra = acir_format::create_circuit<UltraCircuitBuilder>(program);
    const fflonk_plonk::CircuitBuilder lowered = fflonk_acir::lower(ultra);

    std::cout << format("{\"functions\": [\n  {\n        \"acir_opcodes\": ",
                        num_acir_opcodes,
                        ",\n        \"circuit_size\": ",
                        lowered.num_rows(),
                        "\n  }\n]}");
}

void FflonkAPI::write_solidity_verifier(const Flags&,
                                        const std::filesystem::path& output_path,
                                        const std::filesystem::path& vk_path)
{
    fflonk_plonk::VerificationKey key;
    if (!fflonk_plonk::VerificationKey::from_buffer(read_vk_file(vk_path), key)) {
        throw_or_abort("fflonk: could not parse the verification key");
    }

    // The verifier takes its key as constructor arguments, so one contract source serves every
    // circuit. Emitting a derived contract that supplies them keeps the audited source verbatim and
    // still hands the user a file they can deploy without knowing what the arguments mean.
    const std::string contract = std::string(FFLONK_VERIFIER_SOURCE) + R"CONTRACT(
/// @notice `FflonkVerifier` with this circuit's verification key baked in.
contract UltraVerifier is FflonkVerifier {
    constructor()
        FflonkVerifier()CONTRACT" +
                                 "\n            " + std::to_string(key.circuit_size) + "," + "\n            " +
                                 std::to_string(key.num_public_inputs) + "," + "\n            " +
                                 to_hex(static_cast<uint256_t>(key.omega)) + "," + "\n            " +
                                 to_hex(static_cast<uint256_t>(key.k1)) + "," + "\n            " +
                                 to_hex(static_cast<uint256_t>(key.k2)) + "," + "\n            " +
                                 to_hex(static_cast<uint256_t>(key.c0.x)) + "," + "\n            " +
                                 to_hex(static_cast<uint256_t>(key.c0.y)) + "\n        )\n    {}\n}\n";

    if (output_path == "-") {
        std::cout << contract;
    } else {
        write_file(output_path,
                   std::span<const uint8_t>(reinterpret_cast<const uint8_t*>(contract.c_str()), contract.size()));
        info("fflonk solidity verifier saved to ", output_path);
    }
}

} // namespace bb
