#pragma once

#include "barretenberg/api/api.hpp"

#include <filesystem>

namespace bb {

/**
 * @brief `--scheme fflonk`: the three-wire system, reached from ACIR by lowering an Ultra circuit.
 *
 * @details The cheapest verifier in this repository - 237,795 gas for a whole transaction, flat in
 * circuit size - at the cost of an arithmetization with three wires, five selectors and no lookups.
 * The lowering in `fflonk_acir/` covers arithmetic gates only and rejects everything else, so this
 * backend serves arithmetic-only programs that want the smallest possible on-chain verifier and the
 * smallest possible trusted code base. Anything with lookups, range constraints, RAM/ROM, elliptic
 * additions or Poseidon2 belongs on `--scheme ultra_fflonk`.
 */
class FflonkAPI : public API {
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
