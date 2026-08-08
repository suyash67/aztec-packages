#pragma once

#include "barretenberg/crypto/keccak/keccak.hpp"
#include "barretenberg/fflonk/polynomial_utils.hpp"

#include <vector>

namespace bb::fflonk_plonk {

/**
 * @brief Keccak Fiat-Shamir, in the one shape a Solidity verifier can reproduce for free.
 *
 * @details A challenge is `keccak256(state || absorbed words) mod r`, and the digest becomes the new
 * state, so the whole history is chained through a single 32-byte value. Every absorbed item is one
 * 32-byte big-endian word; a group element absorbs as `(x, y)` with the point at infinity absorbing
 * as `(0, 0)`, matching how the EVM's precompiles represent it.
 *
 * The initial state is zero rather than absent, so the first challenge has the same shape as every
 * other one - there is no special case for a Solidity mirror to get wrong.
 */
class Transcript {
  public:
    void absorb(const FF& value) { buffer_.push_back(static_cast<uint256_t>(value)); }

    void absorb_word(const uint256_t& value) { buffer_.push_back(value); }

    void absorb(const Commitment& point)
    {
        if (point.is_point_at_infinity()) {
            buffer_.emplace_back(0);
            buffer_.emplace_back(0);
        } else {
            buffer_.push_back(static_cast<uint256_t>(point.x));
            buffer_.push_back(static_cast<uint256_t>(point.y));
        }
    }

    FF squeeze()
    {
        std::vector<uint256_t> input;
        input.reserve(buffer_.size() + 1);
        input.push_back(static_cast<uint256_t>(state_));
        input.insert(input.end(), buffer_.begin(), buffer_.end());
        buffer_.clear();

        state_ = crypto::Keccak::hash(input);
        return state_;
    }

    const FF& state() const { return state_; }

  private:
    std::vector<uint256_t> buffer_;
    FF state_ = FF::zero();
};

} // namespace bb::fflonk_plonk
