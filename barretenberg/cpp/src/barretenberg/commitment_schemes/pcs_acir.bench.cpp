#include "barretenberg/commitment_schemes/dory/dory_honk.hpp"
#include "barretenberg/commitment_schemes/hyrax/hyrax_honk.hpp"
#include "barretenberg/commitment_schemes/kzh/kzh_honk.hpp"
#include "barretenberg/commitment_schemes/ligero/ligero_honk.hpp"
#include "barretenberg/commitment_schemes/mercury/mercury_honk.hpp"
#include "barretenberg/commitment_schemes/pedersen_ipa/ipa_honk.hpp"
#include "barretenberg/commitment_schemes/whir/whir_honk.hpp"
#include "barretenberg/common/get_bytecode.hpp"
#include "barretenberg/dsl/acir_format/acir_format.hpp"
#include "barretenberg/dsl/acir_format/acir_to_constraint_buf.hpp"
#include "barretenberg/dsl/acir_format/serde/witness_stack.hpp"
#include "barretenberg/special_public_inputs/special_public_inputs.hpp"
#include "barretenberg/srs/global_crs.hpp"
#include "barretenberg/ultra_honk/ultra_prover.hpp"
#include "barretenberg/ultra_honk/ultra_verifier.hpp"

#include "barretenberg/env/logstr.hpp"

#include <chrono>
#include <cstdint>
#include <iostream>
#include <string>

// Proves and verifies one ACIR circuit (bytecode + witness, as produced by nargo) under a
// user-selected polynomial commitment scheme, all through the shared UltraHonk arithmetization:
// the KZG baseline runs the production UltraFlavor prover, every other backend runs
// TransparentHonk<Pcs> on the identical UltraCircuitBuilder. Emits a single JSON object on
// stdout so a driver can aggregate medians across fresh-process repetitions (peak_rss_bytes is
// only meaningful when each run is its own process).
//
// Usage: pcs_acir_bench -b <bytecode> -w <witness.gz> --pcs <kzg|mercury|whir|whir-p2|whir-sky|ligero|hyrax|kzh2|ipa|dory>

namespace {

using namespace bb;

constexpr size_t SECURITY_BITS = 100;
constexpr size_t LOG_INV_RATE = 2;

struct Timings {
    double circuit_ms = 0;
    double pk_ms = 0;
    double prove_ms = 0;
    double verify_ms = 0;
    size_t proof_fields = 0;
    size_t log_n = 0;
    size_t num_gates = 0;
    bool verified = false;
};

double ms_since(const std::chrono::steady_clock::time_point& start)
{
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
}

UltraCircuitBuilder build_circuit(const std::string& bytecode_path, const std::string& witness_path)
{
    acir_format::AcirProgram program{
        acir_format::circuit_buf_to_acir_format(get_bytecode(bytecode_path), /*is_mega=*/false), {}
    };
    program.witness = acir_format::witness_buf_to_witness_vector(get_bytecode(witness_path));
    return acir_format::create_circuit<UltraCircuitBuilder>(program);
}

Timings run_kzg(UltraCircuitBuilder& builder, Timings timings)
{
    auto start = std::chrono::steady_clock::now();
    auto prover_instance = std::make_shared<ProverInstance_<UltraFlavor>>(builder);
    auto verification_key = std::make_shared<UltraFlavor::VerificationKey>(prover_instance->get_precomputed());
    timings.pk_ms = ms_since(start);
    timings.log_n = prover_instance->log_dyadic_size();
    timings.num_gates = prover_instance->dyadic_size();

    UltraProver_<UltraFlavor> prover(prover_instance, verification_key);
    start = std::chrono::steady_clock::now();
    HonkProof proof = prover.construct_proof();
    timings.prove_ms = ms_since(start);
    timings.proof_fields = proof.size();

    auto vk_and_hash = std::make_shared<UltraFlavor::VKAndHash>(verification_key);
    start = std::chrono::steady_clock::now();
    UltraVerifier_<UltraFlavor, DefaultIO> verifier(vk_and_hash);
    timings.verified = verifier.verify_proof(proof).result;
    timings.verify_ms = ms_since(start);
    return timings;
}

template <typename Honk> Timings run_transparent(UltraCircuitBuilder& builder, Timings timings)
{
    // The backend config must know the dyadic size up front; size a throwaway copy since
    // ProverInstance finalizes (mutates) the builder it consumes.
    UltraCircuitBuilder sizing_copy = builder;
    const size_t log_n = ProverInstance_<UltraFlavor>(sizing_copy).log_dyadic_size();
    const auto config = Honk::make_config(log_n, SECURITY_BITS, LOG_INV_RATE);

    auto start = std::chrono::steady_clock::now();
    auto pk = Honk::create_proving_key(builder, config);
    timings.pk_ms = ms_since(start);
    timings.log_n = pk.instance->log_dyadic_size();
    timings.num_gates = pk.instance->dyadic_size();
    const auto vk = pk.vk;

    start = std::chrono::steady_clock::now();
    HonkProof proof = Honk::prove(pk);
    timings.prove_ms = ms_since(start);
    timings.proof_fields = proof.size();

    start = std::chrono::steady_clock::now();
    timings.verified = Honk::verify(vk, config, proof);
    timings.verify_ms = ms_since(start);
    return timings;
}

} // namespace

int main(int argc, char** argv)
{
    std::string bytecode_path;
    std::string witness_path;
    std::string pcs;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "-b" && i + 1 < argc) {
            bytecode_path = argv[++i];
        } else if (arg == "-w" && i + 1 < argc) {
            witness_path = argv[++i];
        } else if (arg == "--pcs" && i + 1 < argc) {
            pcs = argv[++i];
        } else {
            std::cerr << "unknown argument: " << arg << "\n";
            return 1;
        }
    }
    if (bytecode_path.empty() || witness_path.empty() || pcs.empty()) {
        std::cerr << "usage: pcs_acir_bench -b <bytecode> -w <witness.gz> --pcs "
                     "<kzg|mercury|whir|whir-p2|whir-sky|ligero|hyrax|kzh2|ipa|dory>\n";
        return 1;
    }

    bb::srs::init_file_crs_factory(bb::srs::bb_crs_path());

    auto start = std::chrono::steady_clock::now();
    bb::UltraCircuitBuilder builder = build_circuit(bytecode_path, witness_path);
    Timings timings;
    timings.circuit_ms = ms_since(start);

    if (pcs == "kzg") {
        timings = run_kzg(builder, timings);
    } else if (pcs == "mercury") {
        timings = run_transparent<bb::mercury::MercuryHonk>(builder, timings);
    } else if (pcs == "whir") {
        timings = run_transparent<bb::whir::WhirHonk<bb::whir::Blake3sMerkleHasher>>(builder, timings);
    } else if (pcs == "whir-p2") {
        timings = run_transparent<bb::whir::WhirHonk<bb::whir::Poseidon2MerkleHasher>>(builder, timings);
    } else if (pcs == "whir-sky") {
        timings = run_transparent<bb::whir::WhirHonk<bb::whir::SkyscraperMerkleHasher>>(builder, timings);
    } else if (pcs == "ligero") {
        timings = run_transparent<bb::ligero::LigeroHonk<bb::whir::Blake3sMerkleHasher>>(builder, timings);
    } else if (pcs == "hyrax") {
        timings = run_transparent<bb::hyrax::HyraxHonk>(builder, timings);
    } else if (pcs == "kzh2") {
        timings = run_transparent<bb::kzh::KzhHonk>(builder, timings);
    } else if (pcs == "ipa") {
        timings = run_transparent<bb::pedersen_ipa::IpaHonk>(builder, timings);
    } else if (pcs == "dory") {
        timings = run_transparent<bb::dory::DoryHonk>(builder, timings);
    } else {
        std::cerr << "unknown pcs: " << pcs << "\n";
        return 1;
    }

    std::cout << "{\"pcs\":\"" << pcs << "\""
              << ",\"log_n\":" << timings.log_n << ",\"num_gates\":" << timings.num_gates
              << ",\"circuit_ms\":" << timings.circuit_ms << ",\"pk_ms\":" << timings.pk_ms
              << ",\"prove_ms\":" << timings.prove_ms << ",\"verify_ms\":" << timings.verify_ms
              << ",\"proof_bytes\":" << timings.proof_fields * 32 << ",\"peak_rss_bytes\":" << peak_rss_bytes()
              << ",\"verified\":" << (timings.verified ? "true" : "false") << "}\n";
    return timings.verified ? 0 : 2;
}
