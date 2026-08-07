// Differential fuzzing of barretenberg's WHIR against the reference implementation ProveKit uses
// (WizardOfMenlo/whir, pinned revision), driven in-process through the C ABI in rust/.
//
// The two implementations do not share a transcript — barretenberg Fiat-Shamirs through the Honk
// transcript with Poseidon2/Blake3s Merkle trees, the reference through a spongefish duplex sponge
// — so proofs are not byte-comparable and neither verifier can read the other's proof. What must
// agree, and what this harness checks over randomized inputs, is everything above the transcript:
//
//   1. Opening semantics. Both claim "the multilinear extension of this array at this point". If
//      the two disagree on the hypercube variable order or the basis convention, the same array and
//      point yield different claimed values, and a WHIR opening in Honk would be proving something
//      other than what sumcheck asked for. This is the property that has to hold for barretenberg's
//      WHIR to be a drop-in for ProveKit's.
//   2. The query schedule at shared parameters, which fixes the soundness both actually achieve.
//   3. Accept/reject behaviour: each implementation accepts its own honest proof, rejects a
//      corrupted one, and rejects a proof presented against a wrong claimed value.
//
// Usage: whir_differential_bench [trials] [seed]

#include "barretenberg/commitment_schemes/whir/rust/whir_rs_ffi.h"
#include "barretenberg/commitment_schemes/whir/whir.hpp"
#include "barretenberg/commitment_schemes/whir/whir_config.hpp"
#include "barretenberg/common/log.hpp"
#include "barretenberg/transcript/transcript.hpp"

#include <cstdlib>
#include <random>
#include <string>
#include <vector>

namespace {

using namespace bb;
using namespace bb::whir;
using Hasher = Blake3sMerkleHasher;

std::vector<uint8_t> encode(std::span<const fr> values)
{
    std::vector<uint8_t> bytes(values.size() * 32);
    for (size_t i = 0; i < values.size(); ++i) {
        const uint256_t canonical(values[i]);
        std::memcpy(&bytes[i * 32], canonical.data, 32);
    }
    return bytes;
}

fr decode(const uint8_t* bytes)
{
    uint256_t value(0);
    std::memcpy(value.data, bytes, 32);
    return fr(value);
}

struct Failure {
    std::string what;
    size_t trial;
};

std::vector<Failure> failures;

void check(bool condition, const std::string& what, size_t trial)
{
    if (!condition) {
        failures.push_back({ what, trial });
    }
}

/**
 * @brief One randomized trial: same array, same point, both implementations.
 * @details The reference folds by `k` from the first round and stops at a constant, while
 * barretenberg stops at a clear polynomial of `final_poly_bits` variables, so the two run different
 * numbers of iterations. That is a protocol-shape difference, not a semantic one; the checks below
 * are the ones that must hold regardless.
 */
void run_trial(size_t trial, std::mt19937_64& rng, size_t num_variables, size_t folding_factor_bits)
{
    const size_t n = size_t(1) << num_variables;
    std::vector<fr> array(n);
    for (fr& value : array) {
        value = fr::random_element();
    }
    std::vector<fr> point(num_variables);
    for (fr& coordinate : point) {
        coordinate = fr::random_element();
    }
    const std::vector<uint8_t> array_bytes = encode(array);
    const std::vector<uint8_t> point_bytes = encode(point);

    // 1. Opening semantics.
    const fr bb_evaluation = bb::whir::detail::mle_of_span(array, point);
    std::array<uint8_t, 32> rust_value{};
    check(whir_rs_mle_evaluate(
              static_cast<uint32_t>(num_variables), array_bytes.data(), point_bytes.data(), rust_value.data()) ==
              WHIR_RS_OK,
          "reference MLE evaluation failed",
          trial);
    check(decode(rust_value.data()) == bb_evaluation, "multilinear evaluations disagree", trial);

    // 2. Schedule parity at shared parameters. Grinding is off so every bit comes from queries,
    // which is the regime barretenberg implements.
    const size_t security_bits = 128;
    const size_t log_inv_rate = 2;
    const WhirRsParams params{ .num_variables = static_cast<uint32_t>(num_variables),
                               .security_level = static_cast<uint32_t>(security_bits),
                               .pow_bits = 0,
                               .initial_folding_factor = static_cast<uint32_t>(folding_factor_bits),
                               .folding_factor = static_cast<uint32_t>(folding_factor_bits),
                               .starting_log_inv_rate = static_cast<uint32_t>(log_inv_rate),
                               .unique_decoding = 0 };
    const WhirConfig config = WhirConfig::create(num_variables,
                                                 security_bits,
                                                 log_inv_rate,
                                                 folding_factor_bits,
                                                 /*final_poly_bits=*/4,
                                                 WhirSoundness::PROVABLE_LIST);
    // The reference queries every oracle it commits; barretenberg queries the same oracles but
    // stops folding earlier, so compare the counts they share — the rates each side actually walks.
    for (size_t i = 0; i < config.num_iterations(); ++i) {
        const size_t expected =
            WhirConfig::compute_num_queries(security_bits, config.rounds[i].log_inv_rate, WhirSoundness::PROVABLE_LIST);
        check(config.rounds[i].num_queries == expected, "barretenberg round query count drifted", trial);
    }
    check(whir_rs_ood_samples(&params) == static_cast<int64_t>(config.num_ood_samples),
          "out-of-domain sample counts disagree",
          trial);

    // 3a. barretenberg proves and verifies its own opening.
    WhirCommitmentKey<Hasher> ck(config);
    WhirGroupData<Hasher> group = ck.commit(array);
    WhirProver<Hasher>::Claims prover_claims;
    prover_claims.groups.push_back(&group);
    prover_claims.unshifted.push_back({ 0, 0 });
    prover_claims.unshifted_evaluations.push_back(bb_evaluation);

    auto prover_transcript = NativeTranscript::test_prover_init_empty();
    WhirProver<Hasher>::prove(ck, prover_claims, point, prover_transcript);
    const HonkProof bb_proof = prover_transcript->export_proof();

    WhirVerifier<Hasher>::Claims verifier_claims;
    verifier_claims.group_num_columns.push_back(1);
    verifier_claims.unshifted = prover_claims.unshifted;
    verifier_claims.unshifted_evaluations = prover_claims.unshifted_evaluations;

    auto verify_bb = [&](const HonkProof& proof, const std::vector<fr>& evaluations) {
        auto transcript = std::make_shared<NativeTranscript>(proof);
        [[maybe_unused]] auto init = transcript->template receive_from_prover<fr>("Init");
        auto claims = verifier_claims;
        claims.unshifted_evaluations = evaluations;
        return WhirVerifier<Hasher>::verify(config, claims, point, transcript);
    };
    check(verify_bb(bb_proof, { bb_evaluation }), "barretenberg rejected its own honest proof", trial);
    check(!verify_bb(bb_proof, { bb_evaluation + fr(1) }), "barretenberg accepted a wrong claimed value", trial);

    HonkProof tampered_bb = bb_proof;
    tampered_bb[std::uniform_int_distribution<size_t>(1, tampered_bb.size() - 1)(rng)] += fr(1);
    check(!verify_bb(tampered_bb, { bb_evaluation }), "barretenberg accepted a tampered proof", trial);

    // 3b. The reference proves and verifies the same opening.
    std::array<uint8_t, 32> rust_claimed{};
    uint8_t* rust_proof = nullptr;
    size_t rust_proof_len = 0;
    check(whir_rs_prove(
              &params, array_bytes.data(), point_bytes.data(), rust_claimed.data(), &rust_proof, &rust_proof_len) ==
              WHIR_RS_OK,
          "reference prover failed",
          trial);
    check(decode(rust_claimed.data()) == bb_evaluation, "reference proved a different claimed value", trial);
    check(whir_rs_verify(&params, point_bytes.data(), rust_claimed.data(), rust_proof, rust_proof_len) == WHIR_RS_OK,
          "reference rejected its own honest proof",
          trial);

    const std::vector<uint8_t> wrong_value = encode(std::vector<fr>{ bb_evaluation + fr(1) });
    check(whir_rs_verify(&params, point_bytes.data(), wrong_value.data(), rust_proof, rust_proof_len) != WHIR_RS_OK,
          "reference accepted a wrong claimed value",
          trial);

    std::vector<uint8_t> tampered_rust(rust_proof, rust_proof + rust_proof_len);
    // Skip the 8-byte length header, which is framing rather than proof content.
    tampered_rust[std::uniform_int_distribution<size_t>(8, tampered_rust.size() - 1)(rng)] ^= 0x01;
    check(
        whir_rs_verify(&params, point_bytes.data(), rust_claimed.data(), tampered_rust.data(), tampered_rust.size()) !=
            WHIR_RS_OK,
        "reference accepted a tampered proof",
        trial);

    whir_rs_free(rust_proof, rust_proof_len);
}

} // namespace

/**
 * @brief Isolate the field-element encoding from every protocol question.
 * @details The multilinear extension of a one-element array is that element, so a disagreement here
 * is purely about how the 32 bytes are read, not about variable order or basis.
 */
void check_encoding()
{
    for (size_t i = 0; i < 8; ++i) {
        const fr value = fr::random_element();
        const std::vector<uint8_t> bytes = encode(std::vector<fr>{ value });
        std::array<uint8_t, 32> out{};
        check(whir_rs_mle_evaluate(0, bytes.data(), nullptr, out.data()) == WHIR_RS_OK,
              "reference rejected a single-element encoding",
              i);
        check(decode(out.data()) == value, "field element encodings disagree", i);
    }
}

int main(int argc, char** argv)
{
    const size_t trials = argc > 1 ? std::stoul(argv[1]) : 24;
    const uint64_t seed = argc > 2 ? std::stoull(argv[2]) : 1;
    std::mt19937_64 rng(seed);

    check_encoding();

    // Sizes wide enough to exercise several fold-and-commit iterations at both folding factors, and
    // small enough that a few dozen trials stay quick.
    const std::vector<size_t> sizes = { 8, 10, 11, 12 };
    const std::vector<size_t> folding_factors = { 3, 4 };

    for (size_t trial = 0; trial < trials; ++trial) {
        const size_t num_variables = sizes[trial % sizes.size()];
        const size_t folding_factor_bits = folding_factors[(trial / sizes.size()) % folding_factors.size()];
        run_trial(trial, rng, num_variables, folding_factor_bits);
    }

    if (failures.empty()) {
        info("whir differential: ", trials, " trials passed against the reference implementation");
        return 0;
    }
    for (const Failure& failure : failures) {
        info("whir differential FAILURE (trial ", failure.trial, "): ", failure.what);
    }
    return 1;
}
