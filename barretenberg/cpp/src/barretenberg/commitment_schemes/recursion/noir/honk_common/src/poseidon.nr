//! Poseidon2 over BN254 Fr, bit-identical to `bb::crypto::Poseidon2<Poseidon2Bn254ScalarFieldParams>`.
//!
//! bb's sponge has rate 3, capacity 1 and domain IV = `input_length << 64` placed in `state[3]`. It
//! absorbs into a 3-element cache and permutes whenever the cache is full or the hash is squeezed,
//! which is the same thing as absorbing the zero-padded input three elements at a time. Writing it
//! that way keeps every index a compile-time constant.

global IV_SHIFT: Field = 18446744073709551616; // 2^64

/// Poseidon2 hash of `input`.
pub fn hash<let N: u32>(input: [Field; N]) -> Field {
    let num_chunks: u32 = (N + 2) / 3;
    let mut state: [Field; 4] = [0, 0, 0, (N as Field) * IV_SHIFT];
    for c in 0..num_chunks {
        for j in 0..3 {
            let idx = c * 3 + j;
            if idx < N {
                state[j] += input[idx];
            }
        }
        state = std::hash::poseidon2_permutation(state);
    }
    state[0]
}

/// Poseidon2 hash of `prefix` followed by `values`, without materializing the concatenation.
pub fn hash_prefixed<let N: u32>(prefix: Field, values: [Field; N]) -> Field {
    let total: u32 = N + 1;
    let num_chunks: u32 = (total + 2) / 3;
    let mut state: [Field; 4] = [0, 0, 0, (total as Field) * IV_SHIFT];
    for c in 0..num_chunks {
        for j in 0..3 {
            let idx = c * 3 + j;
            if idx == 0 {
                state[j] += prefix;
            } else if idx < total {
                state[j] += values[idx - 1];
            }
        }
        state = std::hash::poseidon2_permutation(state);
    }
    state[0]
}

/// `Poseidon2MerkleHasher::hash_leaf` with no salt: the leaf tag 0 followed by the leaf's values.
pub fn hash_leaf<let N: u32>(values: [Field; N]) -> Field {
    hash_prefixed(0, values)
}

/// `Poseidon2MerkleHasher::hash_node`: the node tag 1 followed by the two child digests.
pub fn hash_node(left: Field, right: Field) -> Field {
    hash([1, left, right])
}

/// Leaf hash of the first `width` entries of `values`, for arrays sized to a worst case above the
/// width a given tree actually uses. `width` is a compile-time constant at every call site, so the
/// permutation count is exactly the one the real leaf needs.
pub fn hash_leaf_prefix<let N: u32>(values: [Field; N], width: u32) -> Field {
    let total: u32 = width + 1;
    let max_chunks: u32 = (N + 3) / 3;
    let mut state: [Field; 4] = [0, 0, 0, (total as Field) * IV_SHIFT];
    for c in 0..max_chunks {
        if c * 3 < total {
            for j in 0..3 {
                let idx = c * 3 + j;
                if (idx != 0) & (idx < total) {
                    state[j] += values[idx - 1];
                }
            }
            state = std::hash::poseidon2_permutation(state);
        }
    }
    state[0]
}
