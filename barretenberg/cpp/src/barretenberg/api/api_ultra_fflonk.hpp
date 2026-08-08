#pragma once

#include "barretenberg/api/api.hpp"

#include <filesystem>

namespace bb {

/**
 * @brief `--scheme ultra_fflonk`: fflonk over the arithmetization UltraHonk already proves.
 *
 * @details The circuit is built by the same ACIR frontend and the same `UltraCircuitBuilder` the
 * Honk path uses, so every constraint type, lookup table and custom gate is shared rather than
 * reimplemented; only the proof system on top of the trace is different. What that buys is a
 * verifier whose group work is a fixed handful of scalar multiplications and one pairing, and whose
 * proof is a fixed 2,112 bytes - the shape that is cheap to verify on Ethereum.
 */
class UltraFflonkAPI : public API {
  public:
    bool check(const Flags& flags,
               const std::filesystem::path& bytecode_path,
               const std::filesystem::path& witness_path) override;

    void prove(const Flags& flags,
               const std::filesystem::path& bytecode_path,
               const std::filesystem::path& witness_path,
               const std::filesystem::path& vk_path,
               const std::filesystem::path& output_dir);

    bool verify(const Flags& flags,
                const std::filesystem::path& public_inputs_path,
                const std::filesystem::path& proof_path,
                const std::filesystem::path& vk_path) override;

    void write_vk(const Flags& flags,
                  const std::filesystem::path& bytecode_path,
                  const std::filesystem::path& output_path) override;

    void gates(const Flags& flags, const std::filesystem::path& bytecode_path) override;

    void write_solidity_verifier(const Flags& flags,
                                 const std::filesystem::path& output_path,
                                 const std::filesystem::path& vk_path) override;
};

} // namespace bb
