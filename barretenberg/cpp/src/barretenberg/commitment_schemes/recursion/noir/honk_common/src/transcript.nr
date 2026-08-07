//! Fiat-Shamir transcript matching `bb::NativeTranscript` (`BaseTranscript<FrCodec, Poseidon2>`).
//!
//! A challenge is `Poseidon2(previous_challenge || round_buffer)`, where the round buffer holds
//! everything absorbed since the last challenge and the previous challenge is omitted only for the
//! very first one. Values the prover sends "unhashed" never enter the buffer - they are bound to the
//! transcript by a Merkle root instead - so the recursive verifier takes them as separate witnesses
//! and pays nothing to re-absorb them.

use crate::poseidon;

pub struct Transcript {
    pub previous: Field,
}

impl Transcript {
    /// Open a transcript on its first challenge: `Poseidon2(round_buffer)`, no previous challenge.
    pub fn start<let N: u32>(data: [Field; N]) -> (Self, Field) {
        let challenge = poseidon::hash(data);
        (Transcript { previous: challenge }, challenge)
    }

    /// The next challenge over a round buffer of `data`.
    pub fn challenge<let N: u32>(&mut self, data: [Field; N]) -> Field {
        self.previous = poseidon::hash_prefixed(self.previous, data);
        self.previous
    }

    /// The next challenge over an empty round buffer (nothing was sent since the last one).
    pub fn challenge_empty(&mut self) -> Field {
        self.previous = poseidon::hash([self.previous]);
        self.previous
    }

    /// The next challenge over a single absorbed element.
    pub fn challenge_1(&mut self, a: Field) -> Field {
        self.previous = poseidon::hash([self.previous, a]);
        self.previous
    }
}

/// Low `bits` bits of a challenge's canonical integer, the query-index derivation both WHIR and
/// Ligero use (`uint256(challenge).data[0] & mask`; for a power-of-two bound the mask and the
/// modulo agree). `to_le_bytes` is canonical, so the decomposition is unique.
pub fn index_from_challenge(challenge: Field, bits: u32) -> u32 {
    let bytes: [u8; 32] = challenge.to_le_bytes();
    let mut index: u32 = 0;
    // 32 bits is more than any schedule here needs, and the mask below discards the rest.
    for i in 0..4 {
        index += (bytes[i] as u32) << (8 * i);
    }
    index & ((1 << bits) - 1)
}
