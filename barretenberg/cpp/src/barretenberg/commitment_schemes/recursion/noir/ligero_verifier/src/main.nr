//! Recursive verifier for a Ligero-Honk proof (transparent UltraHonk with the Ligero tensor PCS).
//!
//! Mirrors `bb::ligero::LigeroVerifier::verify`. The prover interleaves each committed group's row
//! codewords so that one Merkle leaf is one column of the whole group; a query opens that column,
//! authenticates it, and checks the three claimed row combinations against the encodings of the
//! combined rows at that column. The verifier never needs a whole encoding, only its value at the
//! queried index, so the circuit evaluates `enc(w)[idx] = sum_c w[c] * (omega^idx)^c` directly.
//!
//! The public statement is the inner circuit's public inputs together with the verification key's
//! precomputed-column root; the proof itself is private witness.

mod params;

use honk_common::poseidon;
use honk_common::shell;
use honk_common::transcript::index_from_challenge;

use params::{
    GROUP_COLUMN_OFFSET, GROUP_COLUMNS, IS_VIRTUAL, LG_CODEWORD_LENGTH, LG_GROUP_ROW_OFFSET,
    LG_GROUP_ROWS, LG_INDEX_BITS, LG_LOG_NUM_COLS, LG_NUM_COLS, LG_NUM_QUERIES, LG_NUM_ROWS,
    LG_OMEGA, LG_TOTAL_GROUP_ROWS, LOG_N, NUM_ALL_ENTITIES, NUM_GROUPS, NUM_PUBLIC_INPUTS,
    NUM_SHIFTED_CLAIMS, NUM_UNSHIFTED_CLAIMS, PUB_INPUTS_OFFSET, SHIFTED_COLUMN, SHIFTED_ENTITY,
    SHIFTED_GROUP, TOTAL_COLUMNS, UNSHIFTED_COLUMN, UNSHIFTED_ENTITY, UNSHIFTED_GROUP,
    VK_VIRTUAL_MASK, LAGRANGE_FIRST_ROW, LAGRANGE_LAST_ROW, BATCHED_LENGTH,
};

global PATH_TOTAL: u32 = LG_NUM_QUERIES * NUM_GROUPS * LG_INDEX_BITS;
global VALUE_TOTAL: u32 = LG_NUM_QUERIES * LG_TOTAL_GROUP_ROWS;

/// The eq tensor of `u`: `table[i] = prod_j (i_j ? u_j : 1 - u_j)`, `N = 2^K`.
fn eq_tensor<let K: u32, let N: u32>(u: [Field; K]) -> [Field; N] {
    let mut table: [Field; N] = [0; N];
    table[0] = 1;
    for i in 0..K {
        let built: u32 = 1 << (i);
        for j in 0..N {
            if j < built {
                let high = table[j] * u[i];
                table[j + built] = high;
                table[j] = table[j] - high;
            }
        }
    }
    table
}

/// `enc(row)[idx]`: the row's univariate evaluated at the codeword domain point `omega^idx`.
fn encode_at<let N: u32>(row: [Field; N], point: Field) -> Field {
    let mut accumulator: Field = 0;
    for c in 0..N {
        accumulator = accumulator * point + row[N - 1 - c];
    }
    accumulator
}

fn main(
    init: Field,
    vk_root: pub Field,
    public_inputs: pub [Field; NUM_PUBLIC_INPUTS],
    witness_roots: [Field; 3],
    sumcheck_univariates: [[Field; BATCHED_LENGTH]; LOG_N],
    sumcheck_evaluations: [Field; NUM_ALL_ENTITIES],
    combined_rows: [Field; 3 * LG_NUM_COLS],
    combined_digest: Field,
    opened_values: [Field; VALUE_TOTAL],
    opened_paths: [Field; PATH_TOTAL],
) {
    let vk = shell::Vk {
        log_n: LOG_N as Field,
        num_public_inputs: NUM_PUBLIC_INPUTS as Field,
        pub_inputs_offset: PUB_INPUTS_OFFSET as Field,
        virtual_mask: VK_VIRTUAL_MASK,
        lagrange_first_row: LAGRANGE_FIRST_ROW,
        lagrange_last_row: LAGRANGE_LAST_ROW,
        root: vk_root,
    };
    let out = shell::verify(
        init,
        vk,
        public_inputs,
        witness_roots,
        sumcheck_univariates,
        sumcheck_evaluations,
        IS_VIRTUAL,
    );
    let mut transcript = out.transcript;
    let u = out.challenge;

    // The opening argument's first challenge is what absorbs the claimed evaluations.
    let rho = transcript.challenge(sumcheck_evaluations);

    assert(combined_digest == poseidon::hash_leaf(combined_rows), "combined-row digest");

    // Split the sumcheck point across the matrix: low variables index columns, high ones rows.
    let mut u_cols: [Field; LG_LOG_NUM_COLS] = [0; LG_LOG_NUM_COLS];
    for i in 0..LG_LOG_NUM_COLS {
        u_cols[i] = u[i];
    }
    let mut u_rows: [Field; LOG_N - LG_LOG_NUM_COLS] = [0; LOG_N - LG_LOG_NUM_COLS];
    for i in 0..(LOG_N - LG_LOG_NUM_COLS) {
        u_rows[i] = u[LG_LOG_NUM_COLS + i];
    }
    let a: [Field; LG_NUM_COLS] = eq_tensor(u_cols);
    let b: [Field; LG_NUM_ROWS] = eq_tensor(u_rows);

    // Per-(group, column) batching scalars and the batched claim value.
    let mut scalar_u: [Field; TOTAL_COLUMNS] = [0; TOTAL_COLUMNS];
    let mut scalar_s: [Field; TOTAL_COLUMNS] = [0; TOTAL_COLUMNS];
    let mut expected: Field = 0;
    let mut rho_power: Field = 1;
    for i in 0..NUM_UNSHIFTED_CLAIMS {
        let slot = GROUP_COLUMN_OFFSET[UNSHIFTED_GROUP[i]] + UNSHIFTED_COLUMN[i];
        scalar_u[slot] += rho_power;
        expected += rho_power * sumcheck_evaluations[UNSHIFTED_ENTITY[i]];
        rho_power *= rho;
    }
    for l in 0..NUM_SHIFTED_CLAIMS {
        let slot = GROUP_COLUMN_OFFSET[SHIFTED_GROUP[l]] + SHIFTED_COLUMN[l];
        scalar_s[slot] += rho_power;
        expected += rho_power * sumcheck_evaluations[SHIFTED_ENTITY[l]];
        rho_power *= rho;
    }

    let mut w_u: [Field; LG_NUM_COLS] = [0; LG_NUM_COLS];
    let mut w_s: [Field; LG_NUM_COLS] = [0; LG_NUM_COLS];
    let mut w_s2: [Field; LG_NUM_COLS] = [0; LG_NUM_COLS];
    for c in 0..LG_NUM_COLS {
        w_u[c] = combined_rows[c];
        w_s[c] = combined_rows[LG_NUM_COLS + c];
        w_s2[c] = combined_rows[2 * LG_NUM_COLS + c];
    }

    // Claim equation: the rank-2 shift decomposition of the batched multilinear claim.
    let mut total: Field = 0;
    for c in 0..LG_NUM_COLS {
        total += w_u[c] * a[c];
        if c >= 1 {
            total += w_s[c] * a[c - 1];
        }
    }
    total += a[LG_NUM_COLS - 1] * w_s2[0];
    assert(total == expected, "ligero claim equation");

    // Column queries.
    let mut value_cursor: u32 = 0;
    let mut path_cursor: u32 = 0;
    for s in 0..LG_NUM_QUERIES {
        let challenge = if s == 0 {
            transcript.challenge_1(combined_digest)
        } else {
            transcript.challenge_empty()
        };
        // The codeword length is a power of two, so the bounded reduction is a bit mask.
        let idx = index_from_challenge(challenge, LG_INDEX_BITS);

        let mut acc_u: Field = 0;
        let mut acc_s: Field = 0;
        let mut acc_s2: Field = 0;
        for g in 0..NUM_GROUPS {
            let mut column: [Field; LG_TOTAL_GROUP_ROWS] = [0; LG_TOTAL_GROUP_ROWS];
            for r in 0..LG_GROUP_ROWS[g] {
                column[r] = opened_values[value_cursor + r];
            }
            let mut path: [Field; LG_INDEX_BITS] = [0; LG_INDEX_BITS];
            for l in 0..LG_INDEX_BITS {
                path[l] = opened_paths[path_cursor + l];
            }
            // Only the group's own rows are hashed; the tail of `column` is padding.
            let root = merkle_root_prefix(idx, column, LG_GROUP_ROWS[g], path);
            let expected_root = if g == 0 {
                vk_root
            } else {
                witness_roots[g - 1]
            };
            assert(root == expected_root, "ligero column authentication");

            for p in 0..GROUP_COLUMNS[g] {
                let slot = GROUP_COLUMN_OFFSET[g] + p;
                for r in 0..LG_NUM_ROWS {
                    let value = column[p * LG_NUM_ROWS + r];
                    acc_u += scalar_u[slot] * b[r] * value;
                    acc_s += scalar_s[slot] * b[r] * value;
                    if r > 0 {
                        acc_s2 += scalar_s[slot] * b[r - 1] * value;
                    }
                }
            }
            value_cursor += LG_GROUP_ROWS[g];
            path_cursor += LG_INDEX_BITS;
        }

        let point = pow_u32(LG_OMEGA, idx, LG_INDEX_BITS);
        assert(acc_u == encode_at(w_u, point), "ligero column check (unshifted)");
        assert(acc_s == encode_at(w_s, point), "ligero column check (shifted)");
        assert(acc_s2 == encode_at(w_s2, point), "ligero column check (shift kernel)");
    }
    assert(value_cursor == VALUE_TOTAL, "opened values consumed");
    assert(path_cursor == PATH_TOTAL, "opened paths consumed");
    // `LG_CODEWORD_LENGTH` documents the query domain; the index mask above enforces it.
    assert(LG_CODEWORD_LENGTH == 1 << (LG_INDEX_BITS));
}

/// Merkle root from a leaf whose values are the first `width` entries of `column`.
fn merkle_root_prefix<let N: u32, let D: u32>(
    index: u32,
    column: [Field; N],
    width: u32,
    path: [Field; D],
) -> Field {
    // `width` is a compile-time constant at every call site, so the leaf hash is a fixed shape.
    let mut digest = poseidon::hash_leaf_prefix(column, width);
    let mut walk = index;
    for level in 0..D {
        digest = if (walk & 1) == 1 {
            poseidon::hash_node(path[level], digest)
        } else {
            poseidon::hash_node(digest, path[level])
        };
        walk >>= 1;
    }
    digest
}

/// `base^exponent` for an exponent below `2^bits`.
fn pow_u32(base: Field, exponent: u32, bits: u32) -> Field {
    let mut result: Field = 1;
    let mut square = base;
    for i in 0..bits {
        result *= if ((exponent >> (i)) & 1) == 1 {
            square
        } else {
            1
        };
        square *= square;
    }
    result
}
