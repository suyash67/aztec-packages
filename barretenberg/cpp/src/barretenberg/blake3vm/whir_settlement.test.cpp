#include <algorithm>
#include <chrono>
#include <gtest/gtest.h>
#include <iomanip>
#include <optional>
#include <sstream>

#include "barretenberg/blake3vm/blake3vm_circuit_builder.hpp"
#include "barretenberg/blake3vm/blake3vm_prover.hpp"
#include "barretenberg/blake3vm/blake3vm_verifier.hpp"
#include "barretenberg/circuit_checker/circuit_checker.hpp"
#include "barretenberg/commitment_schemes/commitment_key.hpp"
#include "barretenberg/commitment_schemes/whir/stdlib/transparent_honk_recursive_verifier.hpp"
#include "barretenberg/flavor/mega_flavor.hpp"
#include "barretenberg/numeric/bitop/get_msb.hpp"
#include "barretenberg/polynomials/polynomial.hpp"
#include "barretenberg/special_public_inputs/special_public_inputs.hpp"
#include "barretenberg/srs/global_crs.hpp"
#include "barretenberg/stdlib_circuit_builders/mega_circuit_builder.hpp"
#include "barretenberg/stdlib_circuit_builders/mock_circuits.hpp"
#include "barretenberg/ultra_honk/prover_instance.hpp"
#include "barretenberg/ultra_honk/ultra_prover.hpp"
#include "barretenberg/ultra_honk/ultra_verifier.hpp"

/**
 * The settlement pipeline a sequencer would run for client-side WHIR-Honk proofs, measured end to
 * end, with the Merkle hashing of the recursion step priced both ways: performed inside the
 * UltraHonk circuit, and delegated to the Blake3VM under the linking argument.
 *
 * 1. A client proves a statement with transparent UltraHonk over WHIR, Blake3s Merkle trees.
 * 2. The sequencer verifies it natively - milliseconds - and soft-confirms.
 * 3. The sequencer verifies it again inside a Honk circuit. Every Blake3s compression on that path
 *    is delegated to the Blake3VM; the circuit witnesses the hashed bytes eight at a time and
 *    commits them in its databus calldata column, whose commitment must equal the VM proof's link
 *    column commitment (README.md, "The linking argument").
 * 4. That circuit is proved with MegaHonk + KZG: the proof that settles. The settlement verifier
 *    checks both proofs and the commitment equality between them.
 * 5. Several inner proofs share one outer circuit, so the settled proof size is constant.
 */
namespace bb {
namespace {

using namespace bb::whir;
using namespace bb::whir::recursion;

using Builder = UltraCircuitBuilder;
using MegaBuilder = MegaCircuitBuilder;
using FF = stdlib::field_t<Builder>;
using MFF = stdlib::field_t<MegaBuilder>;
using Commitment = curve::BN254::AffineElement;

/**
 * @brief The hash workload of one delegated WHIR verification.
 *
 * @details Records every Blake3s call the verifier makes, in the form the VM consumes (the exact
 * byte string), together with bookkeeping the tests use: which proof-stream witnesses reached the
 * oracle as digests (the tamper candidates), how many 8-byte chunks the circuit committed, and what
 * the chunk constraints cost. With `replay` set the oracle answers with given digests rather than
 * real hashes - precisely the freedom a prover has over an oracle nothing binds to the VM.
 */
class HashOracle {
  public:
    struct Call {
        std::vector<uint8_t> input;
        std::array<uint8_t, 32> digest;
    };

    void reset()
    {
        calls_.clear();
        replay_.clear();
        received_digest_indices_.clear();
        link_gates_ = 0;
        chunks_ = 0;
    }

    /** @brief Answer one Blake3s call, recording the input the circuit actually supplied. */
    std::array<uint8_t, 32> hash(std::span<const uint8_t> input)
    {
        std::array<uint8_t, 32> digest{};
        if (replay_.empty()) {
            blake3::blake3s(input, digest);
        } else {
            digest = replay_.at(calls_.size()).digest;
        }
        calls_.push_back({ std::vector<uint8_t>(input.begin(), input.end()), digest });
        return digest;
    }

    void set_replay(std::vector<Call> answers) { replay_ = std::move(answers); }

    /** @brief A digest the circuit read off the proof stream rather than got from the oracle. */
    void record_received_digest(uint32_t index) { received_digest_indices_.push_back(index); }
    const std::vector<uint32_t>& received_digests() const { return received_digest_indices_; }

    void add_link_gates(size_t gates) { link_gates_ += gates; }
    void add_chunk() { ++chunks_; }

    const std::vector<Call>& calls() const { return calls_; }
    size_t link_gates() const { return link_gates_; }
    size_t chunks() const { return chunks_; }

  private:
    std::vector<Call> calls_;
    std::vector<Call> replay_;
    size_t link_gates_ = 0;
    size_t chunks_ = 0;
    std::vector<uint32_t> received_digest_indices_;
};

// The stdlib hasher interface is entirely static, so the active oracle is process-global. One
// circuit is built at a time here.
HashOracle& oracle()
{
    static HashOracle instance;
    return instance;
}

/**
 * @brief `StdlibBlake3sHasher` with the hashing delegated to the Blake3VM.
 *
 * @details Same interface, same digests, but a digest enters the circuit as a witness rather than
 * as the output of an in-circuit compression, and every byte the oracle hashes is witnessed eight
 * at a time and committed in the circuit's databus calldata column - the circuit's half of the
 * linking argument. The VM proof carries the same chunk sequence in its link column, and equality
 * of the two commitments is what binds the oracle's digests to hashes the VM actually proved.
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
            // One VM hash per iteration: witness the chunks of its input, then of its output. A
            // chained call's buffer opens with the previous digest instead of the tag byte, which
            // leaves every piece 8-byte aligned.
            ChunkEmitter emitter(builder, chained ? std::nullopt : std::optional<uint8_t>(uint8_t(0)));
            if (chained) {
                emitter.add_piece(result[0], 16);
                emitter.add_piece(result[1], 16);
            }
            for (size_t t = 0; t < chunk; ++t) {
                const FF& value = values[absorbed + t];
                tag = OriginTag(tag, value.get_origin_tag());
                whir::detail::append_fr_bytes(buffer, value.get_value());
                emitter.add_piece(value, 32);
            }
            absorbed += chunk;
            digest = oracle().hash(buffer);
            result = witness_digest(builder, digest, tag);
            emitter.finish_message_chunks();
            emitter.add_digest(result);
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

        ChunkEmitter emitter(builder, uint8_t(1));
        for (const FF& half : { left[0], left[1], right[0], right[1] }) {
            emitter.add_piece(half, 16);
        }
        emitter.finish_message_chunks();
        emitter.add_digest(result);
        return result;
    }

    static Digest from_fields(Builder_&, std::span<const FF> fields)
    {
        for (size_t half = 0; half < DIGEST_NUM_FIELDS; ++half) {
            if (!fields[half].is_constant()) {
                oracle().record_received_digest(fields[half].get_witness_index());
            }
        }
        return { fields[0], fields[1] };
    }
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
    /**
     * @brief Witnesses one oracle call's 8-byte chunks, proves they recombine to the field elements
     * the circuit holds, and appends each chunk to the databus calldata column - in exactly the
     * order the VM's claim vector lists them: the (zero-padded) message blocks, then the digest.
     *
     * @details A tagged buffer offsets every piece one byte against the chunk grid, so each piece
     * splits as (7 bytes, 8-byte middles, 1 byte); the leading 7 bytes complete the pending chunk
     * and the top byte waits for the next piece. Both split parts are range-constrained: the chunks
     * are pinned to the VM's byte-built claims by the commitment equality, and the ranges make the
     * piece's byte decomposition unique, so the recombination pins the field element itself.
     * Without them a prover could shift one piece and a neighbour's split in tandem, claiming
     * different field elements recombine to the same honest chunks. Untagged (chained) buffers are
     * 8-byte aligned: the split parts are the chunks, pinned directly, and need no range checks.
     */
    class ChunkEmitter {
      public:
        ChunkEmitter(Builder_& builder, std::optional<uint8_t> tag)
            : builder_(builder)
            , gates_before_(builder.num_gates())
        {
            if (tag.has_value()) {
                carry_ = FF(&builder, fr(uint64_t(*tag)));
                data_bytes_ = 1;
            }
        }
        ~ChunkEmitter() { oracle().add_link_gates(builder_.num_gates() - gates_before_); }
        ChunkEmitter(const ChunkEmitter&) = delete;
        ChunkEmitter& operator=(const ChunkEmitter&) = delete;

        /** @brief One 16- or 32-byte piece of the hashed buffer. */
        void add_piece(const FF& element, size_t bytes)
        {
            const uint256_t value(element.get_value());
            const OriginTag tag = element.get_origin_tag();
            FF recombined(fr(0));
            if (carry_.has_value()) {
                FF lo7 = witness_slice(value, 0, 56, tag);
                lo7.create_range_constraint(56, "link: piece split low 7 bytes");
                push_chunk((*carry_ + lo7 * FF(fr(256))).normalize());
                recombined = lo7;
                size_t bit = 56;
                for (; bit + 64 <= 8 * bytes - 8; bit += 64) {
                    const FF middle = witness_slice(value, bit, 64, tag);
                    push_chunk(middle);
                    recombined += middle * FF(fr(uint256_t(1) << bit));
                }
                FF top = witness_slice(value, bit, 8, tag);
                top.create_range_constraint(8, "link: piece split top byte");
                recombined += top * FF(fr(uint256_t(1) << bit));
                carry_ = top;
            } else {
                for (size_t bit = 0; bit < 8 * bytes; bit += 64) {
                    const FF chunk = witness_slice(value, bit, 64, tag);
                    push_chunk(chunk);
                    recombined = (bit == 0) ? chunk : recombined + chunk * FF(fr(uint256_t(1) << bit));
                }
            }
            element.assert_equal(recombined, "link: piece must recombine to its chunks");
            data_bytes_ += bytes;
        }

        /** @brief Flush the trailing byte's chunk and the final block's zero padding. */
        void finish_message_chunks()
        {
            if (carry_.has_value()) {
                push_chunk(carry_->normalize());
                carry_.reset();
            }
            const size_t num_blocks = data_bytes_ <= 64 ? 1 : (data_bytes_ + 63) / 64;
            while (chunks_emitted_ < 8 * num_blocks) {
                push_zero_chunk();
            }
        }

        /** @brief The call's digest, four aligned chunks over its two 16-byte halves. */
        void add_digest(const Digest& digest)
        {
            for (const FF& half : digest) {
                const uint256_t value(half.get_value());
                const OriginTag tag = half.get_origin_tag();
                const FF lo = witness_slice(value, 0, 64, tag);
                const FF hi = witness_slice(value, 64, 64, tag);
                half.assert_equal(lo + hi * FF(fr(uint256_t(1) << 64)), "link: digest half must recombine");
                push_chunk(lo);
                push_chunk(hi);
            }
        }

      private:
        FF witness_slice(const uint256_t& value, size_t bit, size_t width, const OriginTag& tag)
        {
            const uint256_t mask = (uint256_t(1) << width) - 1;
            FF slice = FF::from_witness(&builder_, fr((value >> bit) & mask));
            slice.set_origin_tag(tag);
            return slice;
        }

        void push_chunk(const FF& chunk)
        {
            builder_.add_public_calldata(BusId::KERNEL_CALLDATA, chunk.get_witness_index());
            ++chunks_emitted_;
            oracle().add_chunk();
        }

        void push_zero_chunk()
        {
            builder_.add_public_calldata(BusId::KERNEL_CALLDATA, builder_.zero_idx());
            ++chunks_emitted_;
            oracle().add_chunk();
        }

        Builder_& builder_;
        size_t gates_before_;
        std::optional<FF> carry_; // top byte of the previous piece, or the tag constant
        size_t data_bytes_ = 0;
        size_t chunks_emitted_ = 0;
    };

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
    TransparentHonkRecursiveVerifier<MegaBuilder, RecursionBlake3sWhirPcs, DelegatingBlake3sHasher<MegaBuilder>>;
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

/** @brief Convert a native proof into circuit witnesses, optionally reporting where each one landed. */
template <typename Builder_>
typename StdlibTranscript<Builder_>::Proof to_stdlib_proof(Builder_& builder,
                                                           const HonkProof& proof,
                                                           std::vector<uint32_t>* witness_indices = nullptr)
{
    typename StdlibTranscript<Builder_>::Proof stdlib_proof;
    stdlib_proof.reserve(proof.size());
    for (const fr& element : proof) {
        stdlib_proof.push_back(stdlib::field_t<Builder_>::from_witness(&builder, element));
        if (witness_indices != nullptr) {
            witness_indices->push_back(stdlib_proof.back().get_witness_index());
        }
    }
    return stdlib_proof;
}

/**
 * @brief The commitment a Blake3VM proof would carry for this claim vector: the claims at rows
 * [NUM_DISABLED_ROWS_IN_SUMCHECK, +N) of the link column. The VM's `link_value` polynomial is
 * exactly this layout, so committing it directly gives the same group element without proving.
 */
Commitment commit_link_column(const std::vector<fr>& claims)
{
    const size_t offset = NUM_DISABLED_ROWS_IN_SUMCHECK;
    const size_t dyadic = 1UL << numeric::get_msb(2 * (offset + claims.size()) - 1);
    Polynomial<fr> link(claims.size(), dyadic, offset);
    for (size_t i = 0; i < claims.size(); ++i) {
        link.at(offset + i) = claims[i];
    }
    CommitmentKey<curve::BN254> commitment_key(dyadic);
    return commitment_key.commit(link);
}

/**
 * @brief The commitment a Mega proof of this circuit carries for its calldata column, computed
 * directly: the databus polynomial the prover commits is built from the circuit's bus vector, so
 * committing it against the same SRS gives the same group element without proving.
 */
Commitment commit_calldata(MegaBuilder& circuit)
{
    ProverInstance_<MegaFlavor> instance(circuit);
    CommitmentKey<curve::BN254> commitment_key(instance.dyadic_size());
    return commitment_key.commit(instance.polynomials.kernel_calldata());
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
 * @brief The linking argument catches the forgery an unlinked oracle allows.
 * @details The recursion circuit reads Merkle siblings from the proof stream and feeds them to the
 * hash oracle. With the hashing delegated, nothing *inside the circuit* pins a sibling: tampering
 * with one makes native verification fail while the circuit stays satisfiable when the oracle
 * replays the honest digests. The link closes this outside the circuit: the tampered sibling's
 * bytes flow into the circuit's committed chunk column, so whatever workload the VM proves - the
 * honest calls, or the calls the tampered circuit actually made - its link commitment diverges
 * from the circuit's calldata commitment and the settlement verifier rejects.
 */
TEST_F(WhirSettlementTests, TamperedSiblingDivergesTheLinkCommitments)
{
    const ClientProof client = prove_on_client(/*num_gates=*/100, /*zk=*/false);
    ASSERT_TRUE(InnerHonk::verify(client.pk.vk, client.config, client.proof));

    // An honest build: learn the digests, find the tamper candidates, and show the link accepts.
    std::vector<uint32_t> proof_witnesses;
    std::vector<size_t> candidates;
    std::vector<HashOracle::Call> honest_calls;
    Commitment honest_link;
    {
        MegaBuilder honest;
        oracle().reset();
        const auto stdlib_proof = to_stdlib_proof(honest, client.proof, &proof_witnesses);
        std::ignore = DelegatedVerifier::verify(honest, client.pk.vk, client.config, stdlib_proof);
        ASSERT_TRUE(CircuitChecker::check(honest)) << "the honest delegated circuit must be satisfiable";
        honest_calls = oracle().calls();

        // Merkle path siblings are read off the proof stream as digests. The query openings sit at
        // the end of the proof and are read unhashed, so tampering there leaves every Fiat-Shamir
        // challenge untouched; take the deepest candidates first.
        for (const uint32_t index : oracle().received_digests()) {
            const auto found = std::find(proof_witnesses.begin(), proof_witnesses.end(), index);
            if (found != proof_witnesses.end()) {
                candidates.push_back(size_t(std::distance(proof_witnesses.begin(), found)));
            }
        }
        std::sort(candidates.begin(), candidates.end(), std::greater<>());
        ASSERT_FALSE(candidates.empty()) << "expected digests reaching the oracle from the proof stream";

        // Completeness: the honest circuit's chunk column is the honest VM's claim vector.
        Blake3VMCircuitBuilder vm_builder;
        for (const HashOracle::Call& call : honest_calls) {
            vm_builder.add_hash(call.input);
        }
        ASSERT_EQ(vm_builder.claims.size(), oracle().chunks()) << "circuit and VM disagree on the chunk count";
        honest_link = commit_link_column(vm_builder.claims);
        EXPECT_EQ(commit_calldata(honest), honest_link) << "the honest run must pass the link check";
    }

    // Tamper a sibling so the proof fails natively, and replay the honest digests - the freedom an
    // unlinked oracle would leave a prover.
    bool forgery_caught = false;
    size_t forged_index = 0;
    for (size_t attempt = 0; attempt < 8 && attempt < candidates.size(); ++attempt) {
        forged_index = candidates[attempt];
        HonkProof forged = client.proof;
        forged[forged_index] += 1;
        if (InnerHonk::verify(client.pk.vk, client.config, forged)) {
            continue; // this element is not actually binding for the statement
        }

        MegaBuilder attack;
        oracle().reset();
        oracle().set_replay(honest_calls);
        const auto stdlib_proof = to_stdlib_proof(attack, forged);
        std::ignore = DelegatedVerifier::verify(attack, client.pk.vk, client.config, stdlib_proof);
        // The circuit itself still cannot tell: the oracle's answers are unconstrained witnesses.
        ASSERT_TRUE(CircuitChecker::check(attack)) << "the forged circuit itself stays satisfiable";

        // But the tampered bytes are now committed. Whichever workload the VM proves, its link
        // commitment cannot match: the honest workload differs in the tampered message chunks, and
        // the circuit's actual calls hash to digests that differ from the replayed honest ones.
        const Commitment attack_calldata = commit_calldata(attack);
        EXPECT_NE(attack_calldata, honest_link) << "a VM proof of the honest workload must not match";

        Blake3VMCircuitBuilder tampered_vm;
        for (const HashOracle::Call& call : oracle().calls()) {
            tampered_vm.add_hash(call.input);
        }
        EXPECT_NE(attack_calldata, commit_link_column(tampered_vm.claims))
            << "a VM proof of the tampered workload must not match either";

        forgery_caught = true;
        break;
    }

    EXPECT_TRUE(forgery_caught) << "no binding sibling found to demonstrate the link's rejection";
    info("\n    Forgery caught: proof element ",
         forged_index,
         " is a Merkle sibling the circuit only hands to the hash oracle. Tampering with it\n"
         "    leaves the delegated circuit satisfiable, but its committed chunk column now\n"
         "    diverges from the link column of every VM proof a prover could produce, so the\n"
         "    settlement verifier's commitment-equality check rejects the forgery.\n");
}

/**
 * @brief The whole flow, with the recursion step's Merkle hashing priced both ways.
 * @details Every number printed is measured on this run: nothing is modelled. The delegated
 * circuit's proof carries its chunk column commitment, the VM's proof carries its link column
 * commitment, and the settlement stage checks their equality - the complete linking argument.
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
    MegaBuilder delegated;
    oracle().reset();
    {
        const auto stdlib_proof = to_stdlib_proof(delegated, client.proof);
        const std::vector<MFF> public_inputs =
            DelegatedVerifier::verify(delegated, client.pk.vk, client.config, stdlib_proof);
        for (const MFF& input : public_inputs) {
            delegated.set_public_input(input.get_witness_index());
        }
    }
    const size_t link_gates = oracle().link_gates();
    const size_t chunks_committed = oracle().chunks();
    const std::vector<HashOracle::Call> hash_workload = oracle().calls();
    const size_t delegated_gates = delegated.get_num_finalized_gates_inefficient();

    Section("[3] sequencer verifies the proof inside a Honk circuit")
        .row("(a) Blake3s in circuit", with_unit(baseline_gates, "gates"))
        .row("      of which Merkle hashing", with_unit(hashing_gates, "gates"))
        .row("(b) hashing delegated to the VM", with_unit(delegated_gates, "gates"))
        .row("      of which the link chunks", with_unit(link_gates, "gates"))
        .row("      chunks committed", std::to_string(chunks_committed))
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
    ASSERT_EQ(vm_builder.claims.size(), chunks_committed) << "circuit and VM disagree on the chunk count";

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

    // ---- 4. Settle the delegated recursion circuit with MegaHonk + KZG ----------------------
    DefaultIO::add_default(delegated);
    auto prover_instance = std::make_shared<ProverInstance_<MegaFlavor>>(delegated);
    auto verification_key = std::make_shared<MegaFlavor::VerificationKey>(prover_instance->get_precomputed());
    UltraProver_<MegaFlavor> outer_prover(prover_instance, verification_key);
    start = std::chrono::steady_clock::now();
    const HonkProof outer_proof = outer_prover.construct_proof();
    const size_t outer_prove_ms = elapsed_ms(start);

    auto vk_and_hash = std::make_shared<MegaFlavor::VKAndHash>(verification_key);
    start = std::chrono::steady_clock::now();
    UltraVerifier_<MegaFlavor, DefaultIO> outer_verifier(vk_and_hash);
    const bool outer_ok = outer_verifier.verify_proof(outer_proof).result;
    const size_t outer_verify_ms = elapsed_ms(start);
    EXPECT_TRUE(outer_ok);

    // The linking argument's check: the recursion proof's calldata commitment and the VM proof's
    // link commitment are the same group element, because both columns carry the claim sequence at
    // the same indices over the same SRS. The settlement verifier performs this equality alongside
    // the two proof verifications; nothing about the claims themselves crosses the transcript.
    const Commitment calldata_commitment = outer_verifier.get_witness_commitments().kernel_calldata();
    EXPECT_EQ(calldata_commitment, vm_verifier.link_commitment) << "the link must accept the honest run";

    Section("[4] sequencer settles it with MegaHonk + KZG")
        .row("outer circuit", "2^" + std::to_string(prover_instance->log_dyadic_size()))
        .row("prove", with_unit(outer_prove_ms, "ms"))
        .row("verify", with_unit(outer_verify_ms, "ms"))
        .row("proof", with_unit(outer_proof.size() * sizeof(fr), "B"))
        .row("link check C_calldata == C_link",
             calldata_commitment == vm_verifier.link_commitment ? "accepted" : "REJECTED");

    // ---- 5. Aggregate several client proofs into one settled proof --------------------------
    MegaBuilder aggregator;
    oracle().reset();
    for (size_t i = 0; i < NUM_AGGREGATED; ++i) {
        const auto stdlib_proof = to_stdlib_proof(aggregator, client.proof);
        const std::vector<MFF> public_inputs =
            DelegatedVerifier::verify(aggregator, client.pk.vk, client.config, stdlib_proof);
        for (const MFF& input : public_inputs) {
            aggregator.set_public_input(input.get_witness_index());
        }
    }
    const size_t aggregator_gates = aggregator.get_num_finalized_gates_inefficient();
    // One VM batch proves both verifications' hashing; its claim vector is the aggregator's
    // concatenated chunk column, so the same single commitment equality links the whole batch.
    Blake3VMCircuitBuilder aggregate_vm;
    for (const HashOracle::Call& call : oracle().calls()) {
        aggregate_vm.add_hash(call.input);
    }
    ASSERT_EQ(aggregate_vm.claims.size(), oracle().chunks());
    DefaultIO::add_default(aggregator);

    auto aggregator_instance = std::make_shared<ProverInstance_<MegaFlavor>>(aggregator);
    auto aggregator_vk = std::make_shared<MegaFlavor::VerificationKey>(aggregator_instance->get_precomputed());
    UltraProver_<MegaFlavor> aggregator_prover(aggregator_instance, aggregator_vk);
    start = std::chrono::steady_clock::now();
    const HonkProof aggregated_proof = aggregator_prover.construct_proof();
    const size_t aggregated_prove_ms = elapsed_ms(start);

    auto aggregator_vk_and_hash = std::make_shared<MegaFlavor::VKAndHash>(aggregator_vk);
    UltraVerifier_<MegaFlavor, DefaultIO> aggregated_verifier(aggregator_vk_and_hash);
    const bool aggregated_ok = aggregated_verifier.verify_proof(aggregated_proof).result;
    EXPECT_TRUE(aggregated_ok);
    EXPECT_EQ(aggregated_verifier.get_witness_commitments().kernel_calldata(), commit_link_column(aggregate_vm.claims))
        << "one commitment equality links the aggregated batch to its VM workload";

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
