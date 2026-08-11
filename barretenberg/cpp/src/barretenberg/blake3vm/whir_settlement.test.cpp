#include <chrono>
#include <gtest/gtest.h>
#include <iomanip>
#include <sstream>

#include "barretenberg/blake3vm/blake3vm_circuit_builder.hpp"
#include "barretenberg/blake3vm/blake3vm_prover.hpp"
#include "barretenberg/blake3vm/blake3vm_verifier.hpp"
#include "barretenberg/commitment_schemes/whir/stdlib/transparent_honk_recursive_verifier.hpp"
#include "barretenberg/special_public_inputs/special_public_inputs.hpp"
#include "barretenberg/srs/global_crs.hpp"
#include "barretenberg/stdlib_circuit_builders/mock_circuits.hpp"
#include "barretenberg/ultra_honk/prover_instance.hpp"
#include "barretenberg/ultra_honk/ultra_prover.hpp"
#include "barretenberg/ultra_honk/ultra_verifier.hpp"

/**
 * The settlement pipeline a sequencer would run for client-side WHIR-Honk proofs, measured end to
 * end, with the Merkle hashing of the recursion step priced both ways: performed inside the
 * UltraHonk circuit, and delegated to the Blake3VM.
 *
 * 1. A client proves a statement with transparent UltraHonk over WHIR, Blake3s Merkle trees.
 * 2. The sequencer verifies it natively - milliseconds - and soft-confirms.
 * 3. The sequencer verifies it again inside an UltraHonk circuit. Every Blake3s compression on that
 *    path is what the Blake3VM exists to remove.
 * 4. That circuit is proved with UltraHonk + KZG: the proof that settles.
 * 5. Several inner proofs share one outer circuit, so the settled proof size is constant.
 */
namespace bb {
namespace {

using namespace bb::whir;
using namespace bb::whir::recursion;

using Builder = UltraCircuitBuilder;
using FF = stdlib::field_t<Builder>;
using Bool = stdlib::bool_t<Builder>;

/**
 * @brief The hash workload of one in-circuit WHIR verification, and the fingerprint binding it.
 *
 * @details Records every Blake3s call the verifier makes, in the form the VM consumes (the exact
 * byte string), and accumulates the same calls into a Horner fingerprint
 * `acc ← acc·γ + x` over the field elements the circuit holds: for each call, the values absorbed,
 * the chained-in digest where there is one, and the resulting digest. A linking argument would
 * make the VM recompute this fingerprint from its own trace and assert the two agree; `γ` is a
 * witness here so the accumulator costs what a real challenge-driven one would.
 */
class HashOracle {
  public:
    struct Call {
        std::vector<uint8_t> input;
        std::array<uint8_t, 32> digest;
    };

    void reset(Builder& builder)
    {
        calls_.clear();
        gamma_ = FF::from_witness(&builder, fr::random_element());
        accumulator_ = FF::from_witness(&builder, fr(0));
        // Both stand in for values a real link would draw from the transcript, so neither may enter
        // the circuit as a free witness: absorbing transcript-derived data into one would read as an
        // unconstrained witness meeting provenanced data.
        gamma_.set_origin_tag(OriginTag::constant());
        accumulator_.set_origin_tag(OriginTag::constant());
        link_gates_ = 0;
        absorbed_ = 0;
        builder_ = &builder;
    }

    /** @brief Run one Blake3s call, record it for the VM, and return its digest. */
    std::array<uint8_t, 32> hash(std::span<const uint8_t> input)
    {
        std::array<uint8_t, 32> digest{};
        blake3::blake3s(input, digest);
        calls_.push_back({ std::vector<uint8_t>(input.begin(), input.end()), digest });
        return digest;
    }

    /**
     * @brief Fold one field element into the fingerprint.
     * @details The accumulator's provenance is reset each step. A real link draws its challenge
     * after every hash claim has been committed, so the accumulator would carry that final round's
     * tag; carrying the merged provenance of everything absorbed instead makes the round-ordering
     * check read absorptions from a later round as data folded under an earlier challenge.
     */
    void absorb(const FF& element)
    {
        const size_t before = builder_->num_gates();
        accumulator_ = accumulator_.madd(gamma_, element);
        accumulator_.set_origin_tag(OriginTag::constant());
        link_gates_ += builder_->num_gates() - before;
        ++absorbed_;
    }

    const std::vector<Call>& calls() const { return calls_; }
    size_t link_gates() const { return link_gates_; }
    size_t absorbed() const { return absorbed_; }
    const FF& fingerprint() const { return accumulator_; }

  private:
    std::vector<Call> calls_;
    FF gamma_;
    FF accumulator_;
    size_t link_gates_ = 0;
    size_t absorbed_ = 0;
    Builder* builder_ = nullptr;
};

// The stdlib hasher interface is entirely static, so the active oracle is process-global. One
// circuit is built at a time here.
HashOracle& oracle()
{
    static HashOracle instance;
    return instance;
}

/**
 * @brief `StdlibBlake3sHasher` with the hashing delegated instead of constrained.
 *
 * @details Same interface, same digests, but a digest enters the circuit as a witness bound into
 * the link fingerprint rather than as the output of an in-circuit compression. Two consequences,
 * both of which the measurement below sees: the compressions leave the circuit, and a digest
 * becomes the two field elements it travels the transcript as rather than 32 bytes, so the Merkle
 * path's conditional swap costs two selects per level instead of thirty-two.
 */
template <typename Builder_> class DelegatingBlake3sHasher {
  public:
    using FF = stdlib::field_t<Builder_>;
    using Bool = stdlib::bool_t<Builder_>;
    using NativeHasher = Blake3sMerkleHasher;
    using Digest = std::array<FF, 2>;
    static constexpr size_t DIGEST_NUM_FIELDS = NativeHasher::DIGEST_NUM_FIELDS;
    static constexpr size_t LEAF_CHUNK_VALUES = NativeHasher::LEAF_CHUNK_VALUES;

    static Digest hash_leaf(Builder_& builder, std::span<const FF> values)
    {
        BB_ASSERT(!values.empty(), "a WHIR leaf is never empty");
        std::vector<uint8_t> buffer;
        buffer.push_back(uint8_t(0)); // leaf tag
        std::array<uint8_t, 32> digest{};
        Digest result{};
        OriginTag tag = OriginTag::constant();

        size_t absorbed = 0;
        bool chained = false;
        while (absorbed < values.size()) {
            const size_t chunk = std::min(LEAF_CHUNK_VALUES, values.size() - absorbed);
            for (size_t t = 0; t < chunk; ++t) {
                const FF& value = values[absorbed + t];
                tag = OriginTag(tag, value.get_origin_tag());
                whir::detail::append_fr_bytes(buffer, value.get_value());
            }
            // One VM hash per call: absorb its inputs, then its output.
            if (chained) {
                oracle().absorb(result[0]);
                oracle().absorb(result[1]);
            }
            for (size_t t = 0; t < chunk; ++t) {
                oracle().absorb(values[absorbed + t]);
            }
            absorbed += chunk;
            digest = oracle().hash(buffer);
            result = witness_digest(builder, digest, tag);
            oracle().absorb(result[0]);
            oracle().absorb(result[1]);
            // The next chunk chains on the digest, replacing the tag and everything after it.
            buffer.assign(digest.begin(), digest.end());
            chained = true;
        }
        return result;
    }

    static Digest hash_node(Builder_& builder, const Digest& left, const Digest& right)
    {
        std::array<uint8_t, 65> input{};
        input[0] = uint8_t(1); // node tag
        const auto left_bytes =
            whir::detail::byte_digest_from_fields(std::array<fr, 2>{ left[0].get_value(), left[1].get_value() });
        const auto right_bytes =
            whir::detail::byte_digest_from_fields(std::array<fr, 2>{ right[0].get_value(), right[1].get_value() });
        std::memcpy(input.data() + 1, left_bytes.data(), 32);
        std::memcpy(input.data() + 33, right_bytes.data(), 32);

        const std::array<uint8_t, 32> digest = oracle().hash(input);
        const OriginTag tag(digest_tag(left), digest_tag(right));
        const Digest result = witness_digest(builder, digest, tag);

        for (const FF& element : { left[0], left[1], right[0], right[1], result[0], result[1] }) {
            oracle().absorb(element);
        }
        return result;
    }

    static Digest from_fields(Builder_&, std::span<const FF> fields) { return { fields[0], fields[1] }; }
    static std::vector<FF> to_fields(const Digest& digest) { return { digest[0], digest[1] }; }

    static Digest conditional_assign(const Bool& predicate, const Digest& lhs, const Digest& rhs)
    {
        return { FF::conditional_assign(predicate, lhs[0], rhs[0]), FF::conditional_assign(predicate, lhs[1], rhs[1]) };
    }

    static void assert_equal(const Digest& lhs, const Digest& rhs, const std::string& message)
    {
        lhs[0].assert_equal(rhs[0], message);
        lhs[1].assert_equal(rhs[1], message);
    }

  private:
    /** @brief The oracle's answer, as witnesses carrying the provenance of what produced them. */
    static Digest witness_digest(Builder_& builder, const std::array<uint8_t, 32>& digest, const OriginTag& tag)
    {
        const std::array<fr, 2> fields = whir::detail::byte_digest_to_fields(digest);
        Digest result{ FF::from_witness(&builder, fields[0]), FF::from_witness(&builder, fields[1]) };
        result[0].set_origin_tag(tag);
        result[1].set_origin_tag(tag);
        return result;
    }

    static OriginTag digest_tag(const Digest& digest)
    {
        return OriginTag(digest[0].get_origin_tag(), digest[1].get_origin_tag());
    }
};

using InCircuitVerifier =
    TransparentHonkRecursiveVerifier<Builder, RecursionBlake3sWhirPcs, StdlibBlake3sHasher<Builder>>;
using DelegatedVerifier =
    TransparentHonkRecursiveVerifier<Builder, RecursionBlake3sWhirPcs, DelegatingBlake3sHasher<Builder>>;
using InnerHonk = InCircuitVerifier::NativeHonk;

size_t elapsed_ms(std::chrono::steady_clock::time_point start)
{
    return size_t(
        std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count());
}

/** @brief A small application circuit standing in for the client's statement. */
Builder build_client_circuit(size_t num_gates)
{
    Builder builder;
    MockCircuits::add_arithmetic_gates_with_public_inputs(builder, 16);
    MockCircuits::add_arithmetic_gates(builder, num_gates);
    MockCircuits::add_lookup_gates(builder, 2);
    const size_t rom_id = builder.create_ROM_array(4);
    for (size_t i = 0; i < 4; ++i) {
        builder.set_ROM_element(rom_id, i, builder.add_variable(fr(3 * i + 1)));
    }
    builder.read_ROM_array(rom_id, builder.add_variable(fr(2)));
    return builder;
}

WhirConfig recursion_config(size_t log_n, bool zk)
{
    WhirConfig config = WhirConfig::create(log_n,
                                           /*security_bits=*/100,
                                           /*log_inv_rate=*/4,
                                           /*folding_factor_bits=*/4,
                                           /*final_poly_bits=*/4,
                                           WhirSoundness::REPAIRED_LIST,
                                           zk,
                                           /*max_stack_bits=*/0,
                                           /*initial_folding_factor_bits=*/1,
                                           /*pow_bits=*/20);
    // Grinding is verified with Poseidon2 so the only Blake3s work on the recursion path is Merkle
    // hashing, which is what the VM takes over.
    config.enable_recursion_profile();
    return config;
}

/** @brief One proved client statement, ready for the sequencer. */
struct ClientProof {
    WhirConfig config;
    InnerHonk::ProvingKey pk;
    HonkProof proof;
    size_t log_n = 0;
    size_t prove_ms = 0;
};

ClientProof prove_on_client(size_t num_gates, bool zk)
{
    Builder sizing = build_client_circuit(num_gates);
    const size_t log_n = ProverInstance_<UltraFlavor>(sizing).log_dyadic_size();
    WhirConfig config = recursion_config(log_n, zk);

    Builder circuit = build_client_circuit(num_gates);
    auto pk = InnerHonk::create_proving_key(circuit, config);
    const auto start = std::chrono::steady_clock::now();
    HonkProof proof = InnerHonk::prove(pk);
    return { config, std::move(pk), std::move(proof), log_n, elapsed_ms(start) };
}

/** @brief Convert a native proof into circuit witnesses. */
StdlibTranscript<Builder>::Proof to_stdlib_proof(Builder& builder, const HonkProof& proof)
{
    StdlibTranscript<Builder>::Proof stdlib_proof;
    stdlib_proof.reserve(proof.size());
    for (const fr& element : proof) {
        stdlib_proof.push_back(FF::from_witness(&builder, element));
    }
    return stdlib_proof;
}

std::string with_unit(size_t value, const char* unit)
{
    return std::to_string(value) + " " + unit;
}

/**
 * @brief One titled block of the report, emitted as a single log line.
 * @details `info` annotates every call with the process's resident size, so a row-per-call table
 * would carry an annotation per row.
 */
class Section {
  public:
    explicit Section(const std::string& title) { text_ << "\n" << title << "\n"; }
    Section& row(const char* name, const std::string& value)
    {
        text_ << "    " << std::left << std::setw(34) << name << std::right << std::setw(20) << value << "\n";
        return *this;
    }
    ~Section() { info(text_.str()); }

  private:
    std::ostringstream text_;
};

} // namespace

class WhirSettlementTests : public ::testing::Test {
  public:
    static void SetUpTestSuite() { bb::srs::init_file_crs_factory(bb::srs::bb_crs_path()); }
};

/**
 * @brief The whole flow, with the recursion step's Merkle hashing priced both ways.
 * @details Every number printed is measured on this run: nothing is modelled. The one part of the
 * delegated design that is not yet built is the argument linking the circuit's fingerprint to the
 * VM's trace, so its cost is reported as the fingerprint accumulator alone.
 */
TEST_F(WhirSettlementTests, FullPipelineWithAndWithoutTheVM)
{
    constexpr size_t CLIENT_GATES = 100;
    constexpr size_t NUM_AGGREGATED = 2;

    info("\nWHIR-Honk settlement pipeline, Blake3s Merkle trees\n", std::string(74, '='));

    // ---- 1. Client: prove the statement on the user's device -------------------------------
    // The client proof is zk: WHIR salts its Merkle leaves and blinds the committed polynomials.
    const ClientProof zk_client = prove_on_client(CLIENT_GATES, /*zk=*/true);
    auto start = std::chrono::steady_clock::now();
    const bool zk_ok = InnerHonk::verify(zk_client.pk.vk, zk_client.config, zk_client.proof);
    const size_t zk_verify_ms = elapsed_ms(start);
    EXPECT_TRUE(zk_ok);

    Section("[1] client proves on device (zk WHIR)")
        .row("inner circuit", "2^" + std::to_string(zk_client.log_n))
        .row("prove", with_unit(zk_client.prove_ms, "ms"))
        .row("proof", with_unit(zk_client.proof.size() * sizeof(fr), "B"));

    // ---- 2. Sequencer: native verification, then a soft confirmation ------------------------
    Section("[2] sequencer verifies natively -> soft confirmation")
        .row("verify", with_unit(zk_verify_ms, "ms"))
        .row("soft confirmation", zk_ok ? "sent" : "WITHHELD");

    // The in-circuit verifier does not accept salted leaves yet (`whir_recursive_verifier.hpp`
    // asserts `!config.zk`), so the recursion stage below runs on the same statement proved
    // without zk. Everything else about the two proofs is identical.
    const ClientProof client = prove_on_client(CLIENT_GATES, /*zk=*/false);
    ASSERT_TRUE(InnerHonk::verify(client.pk.vk, client.config, client.proof));

    // ---- 3a. Recursion with the Blake3s compressions inside the circuit ---------------------
    GateReport baseline_report;
    size_t baseline_gates = 0;
    {
        // Scoped: this circuit is several million gates, and only its measurements outlive it.
        Builder baseline;
        const auto stdlib_proof = to_stdlib_proof(baseline, client.proof);
        const std::vector<FF> public_inputs =
            InCircuitVerifier::verify(baseline, client.pk.vk, client.config, stdlib_proof, &baseline_report);
        for (const FF& input : public_inputs) {
            baseline.set_public_input(input.get_witness_index());
        }
        baseline_gates = baseline.get_num_finalized_gates_inefficient();
    }
    size_t hashing_gates = 0;
    for (const std::string& phase : baseline_report.phases()) {
        if (phase == "merkle: leaf hashing" || phase == "merkle: path walk" || phase == "merkle: cap fold") {
            hashing_gates += baseline_report.gates(phase);
        }
    }

    // ---- 3b. Recursion with the hashing delegated to the VM ---------------------------------
    Builder delegated;
    oracle().reset(delegated);
    {
        const auto stdlib_proof = to_stdlib_proof(delegated, client.proof);
        const std::vector<FF> public_inputs =
            DelegatedVerifier::verify(delegated, client.pk.vk, client.config, stdlib_proof);
        for (const FF& input : public_inputs) {
            delegated.set_public_input(input.get_witness_index());
        }
    }
    const size_t link_gates = oracle().link_gates();
    const size_t absorbed = oracle().absorbed();
    const std::vector<HashOracle::Call> hash_workload = oracle().calls();
    // The fingerprint is what a linking argument would expose to the VM side.
    delegated.set_public_input(oracle().fingerprint().normalize().get_witness_index());
    const size_t delegated_gates = delegated.get_num_finalized_gates_inefficient();

    Section("[3] sequencer verifies the proof inside an UltraHonk circuit")
        .row("(a) Blake3s in circuit", with_unit(baseline_gates, "gates"))
        .row("      of which Merkle hashing", with_unit(hashing_gates, "gates"))
        .row("(b) hashing delegated to the VM", with_unit(delegated_gates, "gates"))
        .row("      of which link fingerprint", with_unit(link_gates, "gates"))
        .row("      field elements bound", std::to_string(absorbed))
        .row("recursion circuit shrinks by",
             std::to_string(baseline_gates / std::max<size_t>(delegated_gates, 1)) + "x");

    // ---- 3c. The Blake3VM proves exactly that hash workload ---------------------------------
    Blake3VMCircuitBuilder vm_builder;
    size_t oracle_calls_replayed = 0;
    for (const HashOracle::Call& call : hash_workload) {
        const std::array<uint8_t, 32> digest = vm_builder.add_hash(call.input);
        ASSERT_EQ(digest, call.digest) << "the VM disagrees with the circuit's hash oracle";
        ++oracle_calls_replayed;
    }

    start = std::chrono::steady_clock::now();
    Blake3VMProver vm_prover(vm_builder);
    const HonkProof vm_proof = vm_prover.construct_proof();
    const size_t vm_prove_ms = elapsed_ms(start);

    start = std::chrono::steady_clock::now();
    Blake3VMVerifier vm_verifier(vm_prover.verification_key);
    const bool vm_ok = vm_verifier.verify_proof(vm_proof);
    const size_t vm_verify_ms = elapsed_ms(start);
    EXPECT_TRUE(vm_ok);

    Section("    Blake3VM proof over the same hash workload")
        .row("blake3s calls", std::to_string(hash_workload.size()))
        .row("compressions", std::to_string(vm_builder.num_compressions()))
        .row("VM trace", "2^" + std::to_string(vm_prover.key->log_circuit_size))
        .row("prove", with_unit(vm_prove_ms, "ms"))
        .row("verify", with_unit(vm_verify_ms, "ms"))
        .row("proof", with_unit(vm_proof.size() * sizeof(fr), "B"));

    // ---- 4. Settle the delegated recursion circuit with UltraHonk + KZG ---------------------
    DefaultIO::add_default(delegated);
    auto prover_instance = std::make_shared<ProverInstance_<UltraFlavor>>(delegated);
    auto verification_key = std::make_shared<UltraFlavor::VerificationKey>(prover_instance->get_precomputed());
    UltraProver_<UltraFlavor> outer_prover(prover_instance, verification_key);
    start = std::chrono::steady_clock::now();
    const HonkProof outer_proof = outer_prover.construct_proof();
    const size_t outer_prove_ms = elapsed_ms(start);

    auto vk_and_hash = std::make_shared<UltraFlavor::VKAndHash>(verification_key);
    start = std::chrono::steady_clock::now();
    UltraVerifier_<UltraFlavor, DefaultIO> outer_verifier(vk_and_hash);
    const bool outer_ok = outer_verifier.verify_proof(outer_proof).result;
    const size_t outer_verify_ms = elapsed_ms(start);
    EXPECT_TRUE(outer_ok);

    Section("[4] sequencer settles it with UltraHonk + KZG")
        .row("outer circuit", "2^" + std::to_string(prover_instance->log_dyadic_size()))
        .row("prove", with_unit(outer_prove_ms, "ms"))
        .row("verify", with_unit(outer_verify_ms, "ms"))
        .row("proof", with_unit(outer_proof.size() * sizeof(fr), "B"));

    // ---- 5. Aggregate several client proofs into one settled proof --------------------------
    Builder aggregator;
    oracle().reset(aggregator);
    for (size_t i = 0; i < NUM_AGGREGATED; ++i) {
        const auto stdlib_proof = to_stdlib_proof(aggregator, client.proof);
        const std::vector<FF> public_inputs =
            DelegatedVerifier::verify(aggregator, client.pk.vk, client.config, stdlib_proof);
        for (const FF& input : public_inputs) {
            aggregator.set_public_input(input.get_witness_index());
        }
    }
    aggregator.set_public_input(oracle().fingerprint().normalize().get_witness_index());
    const size_t aggregator_gates = aggregator.get_num_finalized_gates_inefficient();
    DefaultIO::add_default(aggregator);

    auto aggregator_instance = std::make_shared<ProverInstance_<UltraFlavor>>(aggregator);
    auto aggregator_vk = std::make_shared<UltraFlavor::VerificationKey>(aggregator_instance->get_precomputed());
    UltraProver_<UltraFlavor> aggregator_prover(aggregator_instance, aggregator_vk);
    start = std::chrono::steady_clock::now();
    const HonkProof aggregated_proof = aggregator_prover.construct_proof();
    const size_t aggregated_prove_ms = elapsed_ms(start);

    auto aggregator_vk_and_hash = std::make_shared<UltraFlavor::VKAndHash>(aggregator_vk);
    UltraVerifier_<UltraFlavor, DefaultIO> aggregated_verifier(aggregator_vk_and_hash);
    const bool aggregated_ok = aggregated_verifier.verify_proof(aggregated_proof).result;
    EXPECT_TRUE(aggregated_ok);

    Section("[5] " + std::to_string(NUM_AGGREGATED) + " client proofs aggregated into one settled proof")
        .row("aggregator circuit", "2^" + std::to_string(aggregator_instance->log_dyadic_size()))
        .row("gates per client proof", with_unit(aggregator_gates / NUM_AGGREGATED, "gates"))
        .row("prove", with_unit(aggregated_prove_ms, "ms"))
        .row("proof", with_unit(aggregated_proof.size() * sizeof(fr), "B"));
    info("\n    The settled proof is the same size whatever the batch holds; an outer circuit of\n"
         "    2^k rows carries 2^k / (gates per client proof) of them.\n");

    // Delegation is the whole point: the recursion circuit must drop by orders of magnitude, and
    // the hashing must be what left it.
    EXPECT_LT(delegated_gates * 10, baseline_gates);
    EXPECT_GT(hashing_gates * 10, baseline_gates * 9);
    // Every hash the circuit consumed from the oracle is one the VM proved.
    EXPECT_EQ(hash_workload.size(), oracle_calls_replayed);
}

} // namespace bb
