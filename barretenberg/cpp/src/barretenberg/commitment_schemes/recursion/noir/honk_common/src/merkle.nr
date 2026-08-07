//! Merkle authentication against `bb::whir::MerkleTree<Poseidon2MerkleHasher>`.

use crate::poseidon;

/// Root implied by one leaf and its authentication path (`MerkleTree::verify`).
pub fn root_from_path<let W: u32, let D: u32>(
    index: u32,
    values: [Field; W],
    path: [Field; D],
) -> Field {
    let mut digest = poseidon::hash_leaf(values);
    let mut walk = index;
    for level in 0..D {
        let sibling = path[level];
        digest = if (walk & 1) == 1 {
            poseidon::hash_node(sibling, digest)
        } else {
            poseidon::hash_node(digest, sibling)
        };
        walk >>= 1;
    }
    digest
}

/// Root implied by a batch of leaves that share their common path nodes
/// (`MerkleTree::verify_batch`).
///
/// `leaves` must be strictly ascending over its first `num_leaves` entries - the caller checks that
/// and that every query index it needs appears there. The walk climbs the tree level by level,
/// pairing adjacent held nodes when they are siblings and otherwise consuming one transmitted
/// sibling, which is exactly the order the prover emitted them in. Every array is sized for the
/// worst case (all queries distinct) and the per-level node count is carried as a witness, so a
/// shorter batch simply leaves the tail of each level inactive.
/// `batch_root` over leaves whose real width and tree depth are compile-time constants below the
/// array bounds `W` and `D` the circuit is compiled for. One code path then serves batches of
/// different shapes - the four Honk commitment trees and the folded oracles of every WHIR round.
pub fn batch_root<let T: u32, let W: u32, let D: u32, let S: u32>(
    leaves: [u32; T],
    num_leaves: u32,
    values: [[Field; W]; T],
    width: u32,
    siblings: [Field; S],
    depth: u32,
) -> Field {
    let mut level: [u32; T] = leaves;
    let mut digests: [Field; T] = [0; T];
    for i in 0..T {
        digests[i] = poseidon::hash_leaf_prefix(values[i], width);
    }

    let mut size = num_leaves;
    let mut cursor: u32 = 0;
    for depth_index in 0..D {
        if depth_index < depth {
        let mut parents: [u32; T] = [0; T];
        let mut parent_digests: [Field; T] = [0; T];
        let mut out: u32 = 0;
        let mut p: u32 = 0;
        for _step in 0..T {
            // A paired step advances `p` by two, so it can run one past the last held node on the
            // final iteration; the reads are clamped because an inactive step's values are unused.
            let active = p < size;
            let slot = if p < T {
                p
            } else {
                0
            };
            let next_slot = if p + 1 < T {
                p + 1
            } else {
                0
            };
            let current = level[slot];
            let next = level[next_slot];
            let paired = (p + 1 < size) & (next == (current ^ 1));

            let current_digest = digests[slot];
            let sibling = siblings[if cursor < S {
                cursor
            } else {
                0
            }];
            let other = if paired {
                digests[next_slot]
            } else {
                sibling
            };
            // A paired node is the left child (the list ascends, so `current` is even); an unpaired
            // one takes the transmitted sibling on the side its own index bit dictates.
            let swap = !paired & ((current & 1) == 1);
            let left = if swap {
                other
            } else {
                current_digest
            };
            let right = if swap {
                current_digest
            } else {
                other
            };
            let parent_digest = poseidon::hash_node(left, right);

            let write_slot = if out < T {
                out
            } else {
                0
            };
            parents[write_slot] = if active {
                current >> 1
            } else {
                parents[write_slot]
            };
            parent_digests[write_slot] = if active {
                parent_digest
            } else {
                parent_digests[write_slot]
            };
            out += if active {
                1
            } else {
                0
            };
            p += if active {
                if paired {
                    2
                } else {
                    1
                }
            } else {
                0
            };
            cursor += if active & !paired {
                1
            } else {
                0
            };
        }
        level = parents;
        digests = parent_digests;
        size = out;
        }
    }
    assert(size == 1, "merkle batch did not collapse to the root");
    assert(level[0] == 0, "merkle batch ended off the root path");
    digests[0]
}
