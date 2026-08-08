#include "barretenberg/ultra_fflonk/keys.hpp"

#include "barretenberg/common/throw_or_abort.hpp"
#include "barretenberg/crypto/keccak/keccak.hpp"
#include "barretenberg/fflonk/encoding.hpp"

#include <limits>

namespace bb::ultra_fflonk {

FF VerificationKey::hash() const
{
    std::vector<uint256_t> words;
    words.reserve(3 + (2 * NUM_PREPROCESSED_GROUPS));
    words.emplace_back(circuit_size);
    words.emplace_back(num_public_inputs);
    words.emplace_back(pub_inputs_offset);
    for (const Commitment& commitment : preprocessed) {
        if (commitment.is_point_at_infinity()) {
            words.emplace_back(0);
            words.emplace_back(0);
        } else {
            words.push_back(static_cast<uint256_t>(commitment.x));
            words.push_back(static_cast<uint256_t>(commitment.y));
        }
    }
    return crypto::Keccak::hash(words);
}

std::vector<uint8_t> VerificationKey::to_buffer() const
{
    std::vector<uint8_t> buffer;
    buffer.reserve(SIZE_IN_BYTES);
    fflonk_plonk::append_word(buffer, uint256_t(circuit_size));
    fflonk_plonk::append_word(buffer, uint256_t(num_public_inputs));
    fflonk_plonk::append_word(buffer, uint256_t(pub_inputs_offset));
    fflonk_plonk::append_word(buffer, static_cast<uint256_t>(omega));
    for (const Commitment& commitment : preprocessed) {
        fflonk_plonk::append_point(buffer, commitment);
    }
    return buffer;
}

bool VerificationKey::from_buffer(std::span<const uint8_t> buffer, VerificationKey& key)
{
    if (buffer.size() != SIZE_IN_BYTES) {
        return false;
    }

    size_t offset = 0;
    const uint256_t circuit_size = fflonk_plonk::read_word(buffer, offset);
    const uint256_t num_public_inputs = fflonk_plonk::read_word(buffer, offset + 32);
    const uint256_t pub_inputs_offset = fflonk_plonk::read_word(buffer, offset + 64);
    offset += 96;

    // The sizes index into buffers and drive loop bounds, so anything that would not round-trip
    // through `size_t` is rejected here rather than truncated into something plausible.
    if (circuit_size > uint256_t(std::numeric_limits<uint32_t>::max()) ||
        num_public_inputs > uint256_t(std::numeric_limits<uint32_t>::max()) ||
        pub_inputs_offset > uint256_t(std::numeric_limits<uint32_t>::max())) {
        return false;
    }
    key.circuit_size = static_cast<size_t>(static_cast<uint64_t>(circuit_size));
    key.num_public_inputs = static_cast<size_t>(static_cast<uint64_t>(num_public_inputs));
    key.pub_inputs_offset = static_cast<size_t>(static_cast<uint64_t>(pub_inputs_offset));

    if (!fflonk_plonk::read_scalar(buffer, offset, key.omega)) {
        return false;
    }
    for (Commitment& commitment : key.preprocessed) {
        if (!fflonk_plonk::read_point(buffer, offset, commitment)) {
            return false;
        }
    }
    return true;
}

std::vector<FF> unpack_column(const std::vector<FF>& packed, const size_t pack, const size_t index)
{
    BB_ASSERT_GT(pack, index, "column index outside the group");
    std::vector<FF> column;
    column.reserve((packed.size() + pack - 1) / pack);
    for (size_t j = index; j < packed.size(); j += pack) {
        column.push_back(packed[j]);
    }
    return column;
}

ProvingKey preprocess(Flavor::CircuitBuilder& circuit)
{
    ProvingKey key;
    key.instance = std::make_shared<ProverInstance_<Flavor>>(circuit);

    const size_t n = key.instance->dyadic_size();
    if (n < MIN_CIRCUIT_SIZE || (n & (n - 1)) != 0) {
        throw_or_abort("ultra_fflonk: circuit size must be a power of two of at least 8");
    }
    if (key.instance->num_public_inputs() >= n) {
        throw_or_abort("ultra_fflonk: more public inputs than rows");
    }

    key.small_domain = std::make_shared<EvaluationDomain<FF>>(n);
    key.small_domain->compute_lookup_table();
    key.large_domain = std::make_shared<EvaluationDomain<FF>>(QUOTIENT_DOMAIN_FACTOR * n);
    key.large_domain->compute_lookup_table();
    key.commitment_key = std::make_shared<CommitmentKey<Curve>>(SRS_SIZE_FACTOR * n);

    // The precomputed columns arrive in Lagrange form over the trace; the packing that fflonk
    // commits needs them in coefficient form, so each one costs a single inverse FFT here and never
    // again.
    auto precomputed = key.instance->polynomials.get_precomputed();
    std::array<std::vector<FF>, PACK_PREPROCESSED> columns;
    for (size_t g = 0; g < NUM_PREPROCESSED_GROUPS; ++g) {
        for (size_t i = 0; i < PACK_PREPROCESSED; ++i) {
            const auto& polynomial = precomputed[(g * PACK_PREPROCESSED) + i];
            std::vector<FF> lagrange(n, FF::zero());
            for (size_t row = 0; row < n; ++row) {
                lagrange[row] = polynomial[row];
            }
            columns[i] = fflonk_plonk::evaluations_to_coefficients(lagrange, *key.small_domain);
        }
        key.packed_preprocessed[g] = fflonk_plonk::pack_columns(columns);
        key.verification_key.preprocessed[g] =
            key.commitment_key->commit(Polynomial<FF>(std::span<const FF>(key.packed_preprocessed[g])));
    }

    key.verification_key.circuit_size = n;
    key.verification_key.num_public_inputs = key.instance->num_public_inputs();
    key.verification_key.pub_inputs_offset = key.instance->pub_inputs_offset();
    key.verification_key.omega = key.small_domain->root;

    return key;
}

} // namespace bb::ultra_fflonk
