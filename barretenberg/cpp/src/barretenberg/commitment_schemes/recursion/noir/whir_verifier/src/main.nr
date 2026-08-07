//! Recursive verifier for a WHIR-Honk proof (transparent UltraHonk with the WHIR PCS).
//!
//! Mirrors `bb::whir::WhirVerifier::verify` in the interleaved, non-zk layout `WhirHonk` uses: every
//! claimed column is its own constituent of the rho-batched oracle and all of them are claimed at the
//! single sumcheck point, so there is one weight term to start with, no column-stacking challenges
//! and no cross evaluations. Each iteration runs `k` sumcheck rounds against the running claim,
//! receives the folded oracle's root, samples it out of domain, and queries it in domain; the queries
//! of iteration 0 hit the four Honk commitment trees and every later query hits the single folded
//! oracle. The final polynomial arrives in the clear and closes both the oracle-consistency check and
//! the accumulated weighted-sum claim.
//!
//! The public statement is the inner circuit's public inputs together with the verification key's
//! precomputed-column root; the proof itself is private witness.

mod params;

use honk_common::merkle;
use honk_common::shell;
use honk_common::transcript::index_from_challenge;

use params::{
    BATCHED_LENGTH, GROUP_COLUMNS, IS_VIRTUAL, LAGRANGE_FIRST_ROW, LAGRANGE_LAST_ROW, LOG_N,
    NUM_ALL_ENTITIES, NUM_GROUPS, NUM_PUBLIC_INPUTS, PUB_INPUTS_OFFSET, VK_VIRTUAL_MASK,
    WH_BATCH_DEPTH, WH_BATCH_LEAF_OFFSET, WH_BATCH_LEAF_WIDTH, WH_BATCH_MAX_LEAVES,
    WH_BATCH_POSITION_OFFSET, WH_BATCH_SIBLING_OFFSET, WH_BATCH_VALUE_OFFSET,
    WH_CONSTITUENT_COLUMN, WH_CONSTITUENT_ENTITY, WH_CONSTITUENT_SHIFTED, WH_CONSTITUENT_TREE,
    WH_FINAL_SIZE, WH_FINAL_VARS, WH_MAX_BATCH_LEAVES, WH_MAX_DEPTH, WH_MAX_LEAF_WIDTH,
    WH_NUM_BATCHES, WH_NUM_CONSTITUENTS, WH_NUM_ITERATIONS, WH_NUM_OOD, WH_NUM_TERMS,
    WH_MAX_ARITY, WH_MAX_K, WH_MAX_QUERIES, WH_NUM_VARIABLES, WH_ROUND0_ARITY,
    WH_ROUND_ETA_INV, WH_ROUND_INDEX_BITS, WH_ROUND_K,
    WH_ROUND_NUM_VARIABLES, WH_ROUND_OMEGA, WH_ROUND_QUERIES, WH_SUMCHECK_OFFSET,
    WH_TOTAL_LEAVES, WH_TOTAL_POSITIONS, WH_TOTAL_SIBLINGS, WH_TOTAL_SUMCHECK_ROUNDS,
    WH_TOTAL_VALUES,
};

/// 2^-1 over BN254 Fr, the constant `fold_coset` divides the even/odd split by.
global TWO_INV: Field =
    10944121435919637611123202872628637544274182200208017171849102093287904247809;

/// A weight term of the running claim, in the two forms the schedule produces.
///
/// `WeightTerm` is `coeff * prod_j (a_j + b_j X_j)`, and every term here ends the protocol with the
/// final round's variables unbound. That makes the two forms collapse to a scalar plus one point:
/// `eq(u, X)` contributes the multilinear extension of the final polynomial at the unbound tail of
/// `u`, and `pow_y(X)` contributes the final polynomial - read as coefficients - evaluated at
/// `y^(2^bound)`. Binding a variable multiplies `coeff` and squares the running power.
struct Term {
    coeff: Field,
    power: Field, // y^(2^bound); unused by the eq term, which is always term 0
}

fn main(
    init: Field,
    vk_root: pub Field,
    public_inputs: pub [Field; NUM_PUBLIC_INPUTS],
    witness_roots: [Field; 3],
    sumcheck_univariates: [[Field; BATCHED_LENGTH]; LOG_N],
    sumcheck_evaluations: [Field; NUM_ALL_ENTITIES],
    ood_values: [[Field; WH_NUM_CONSTITUENTS]; WH_NUM_OOD],
    round_sumcheck: [[Field; 3]; WH_TOTAL_SUMCHECK_ROUNDS],
    round_roots: [Field; WH_NUM_ITERATIONS],
    round_y_ood: [Field; WH_NUM_ITERATIONS],
    round_nonces: [Field; WH_NUM_ITERATIONS + 1],
    final_poly: [Field; WH_FINAL_SIZE],
    opened_values: [Field; WH_TOTAL_VALUES],
    opened_siblings: [Field; WH_TOTAL_SIBLINGS],
    hint_leaves: [u32; WH_TOTAL_LEAVES],
    hint_num_leaves: [u32; WH_NUM_BATCHES],
    hint_positions: [u32; WH_TOTAL_POSITIONS],
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

    let mut roots: [Field; NUM_GROUPS] = [0; NUM_GROUPS];
    roots[0] = vk_root;
    for g in 1..NUM_GROUPS {
        roots[g] = witness_roots[g - 1];
    }

    // Out-of-domain samples against the initial oracle. The first challenge absorbs the claimed
    // sumcheck evaluations; each later one absorbs the previous sample's answers.
    let mut z_ood_initial: [Field; WH_NUM_OOD] = [0; WH_NUM_OOD];
    for s in 0..WH_NUM_OOD {
        z_ood_initial[s] = if s == 0 {
            transcript.challenge(sumcheck_evaluations)
        } else {
            transcript.challenge(ood_values[s - 1])
        };
    }
    let rho = transcript.challenge(ood_values[WH_NUM_OOD - 1]);
    let gamma_claims = transcript.challenge_empty();

    // sigma_0 and the initial weight terms: the batched claim at u, then one per OOD sample.
    let mut rho_powers: [Field; WH_NUM_CONSTITUENTS] = [0; WH_NUM_CONSTITUENTS];
    let mut rho_power: Field = 1;
    for c in 0..WH_NUM_CONSTITUENTS {
        rho_powers[c] = rho_power;
        rho_power *= rho;
    }

    let mut terms: [Term; WH_NUM_TERMS] = [Term { coeff: 0, power: 0 }; WH_NUM_TERMS];
    let mut num_terms: u32 = 0;
    let mut sigma: Field = 0;
    let mut gamma_power: Field = 1;
    {
        let mut oracle_value: Field = 0;
        for c in 0..WH_NUM_CONSTITUENTS {
            oracle_value += rho_powers[c] * sumcheck_evaluations[WH_CONSTITUENT_ENTITY[c]];
        }
        sigma += gamma_power * oracle_value;
        terms[0] = Term { coeff: gamma_power, power: 0 };
        num_terms = 1;
        gamma_power *= gamma_claims;
    }
    for s in 0..WH_NUM_OOD {
        let mut oracle_value: Field = 0;
        for c in 0..WH_NUM_CONSTITUENTS {
            oracle_value += rho_powers[c] * ood_values[s][c];
        }
        sigma += gamma_power * oracle_value;
        terms[num_terms] = Term { coeff: gamma_power, power: z_ood_initial[s] };
        num_terms += 1;
        gamma_power *= gamma_claims;
    }

    let mut folded_root: Field = 0;
    let mut bound: u32 = 0; // sumcheck variables consumed so far, i.e. the eq term's bound index

    for iteration in 0..WH_NUM_ITERATIONS {
        let k = WH_ROUND_K[iteration];
        let queries = WH_ROUND_QUERIES[iteration];
        let index_bits = WH_ROUND_INDEX_BITS[iteration];
        let omega = WH_ROUND_OMEGA[iteration];
        let eta_inv = WH_ROUND_ETA_INV[iteration];

        // 1. Sumcheck rounds against the running claim.
        let mut alphas: [Field; WH_MAX_K] = [0; WH_MAX_K];
        for j in 0..k {
            let h = round_sumcheck[WH_SUMCHECK_OFFSET[iteration] + j];
            let alpha = transcript.challenge(h);
            assert(h[0] + h[1] == sigma, "whir sumcheck round sum");
            // h given by its values at 0, 1, 2; evaluate at alpha.
            let a2 = (h[2] - h[1] - h[1] + h[0]) * TWO_INV;
            let a1 = h[1] - h[0] - a2;
            sigma = h[0] + alpha * (a1 + alpha * a2);
            alphas[j] = alpha;
            for t in 0..WH_NUM_TERMS {
                if t < num_terms {
                    let factor = if t == 0 {
                        (1 - u[bound]) + (u[bound] + u[bound] - 1) * alpha
                    } else {
                        1 + (terms[t].power - 1) * alpha
                    };
                    terms[t] = Term {
                        coeff: terms[t].coeff * factor,
                        power: terms[t].power * terms[t].power,
                    };
                }
            }
            bound += 1;
        }

        // 2/3. The folded oracle's root and its out-of-domain sample.
        let next_root = round_roots[iteration];
        let z_ood = transcript.challenge_1(next_root);
        let y_ood = round_y_ood[iteration];
        let _pow_seed = transcript.challenge_1(y_ood);
        let nonce = round_nonces[iteration];

        // 4. In-domain queries.
        let mut indices: [u32; WH_MAX_QUERIES] = [0; WH_MAX_QUERIES];
        for s in 0..queries {
            let challenge = if s == 0 {
                transcript.challenge_1(nonce)
            } else {
                transcript.challenge_empty()
            };
            indices[s] = index_from_challenge(challenge, index_bits);
        }

        let arity: u32 = 1 << k;
        let mut cosets: [[Field; WH_MAX_ARITY]; WH_MAX_QUERIES] =
            [[0; WH_MAX_ARITY]; WH_MAX_QUERIES];
        if iteration == 0 {
            // Round 0 opens all four Honk trees; the batched oracle's coset values are the
            // rho-combination of the opened columns, with shifted constituents scaled by x^-1.
            let mut positions: [[u32; WH_MAX_QUERIES]; NUM_GROUPS] =
                [[0; WH_MAX_QUERIES]; NUM_GROUPS];
            for g in 0..NUM_GROUPS {
                positions[g] = authenticate_batch(
                    g,
                    indices,
                    queries,
                    roots[g],
                    opened_values,
                    opened_siblings,
                    hint_leaves,
                    hint_num_leaves,
                    hint_positions,
                );
            }
            for s in 0..queries {
                let x_base = pow_index(omega, indices[s], index_bits);
                let x_base_inv = 1 / x_base;
                let mut x_inverses: [Field; WH_ROUND0_ARITY] = [0; WH_ROUND0_ARITY];
                let mut running = x_base_inv;
                for t in 0..WH_ROUND0_ARITY {
                    x_inverses[t] = running;
                    running *= eta_inv;
                }
                for c in 0..WH_NUM_CONSTITUENTS {
                    let tree = WH_CONSTITUENT_TREE[c];
                    let base = WH_BATCH_VALUE_OFFSET[tree]
                        + positions[tree][s] * WH_BATCH_LEAF_WIDTH[tree]
                        + WH_CONSTITUENT_COLUMN[c] * WH_ROUND0_ARITY;
                    for t in 0..WH_ROUND0_ARITY {
                        let mut value = rho_powers[c] * opened_values[base + t];
                        if WH_CONSTITUENT_SHIFTED[c] {
                            value *= x_inverses[t];
                        }
                        cosets[s][t] += value;
                    }
                }
            }
        } else {
            let batch = NUM_GROUPS + iteration - 1;
            let positions = authenticate_batch(
                batch,
                indices,
                queries,
                folded_root,
                opened_values,
                opened_siblings,
                hint_leaves,
                hint_num_leaves,
                hint_positions,
            );
            for s in 0..queries {
                let base = WH_BATCH_VALUE_OFFSET[batch] + positions[s] * WH_BATCH_LEAF_WIDTH[batch];
                for t in 0..WH_BATCH_LEAF_WIDTH[batch] {
                    cosets[s][t] = opened_values[base + t];
                }
            }
        }

        // 5. Fold each queried coset and combine into the next claim.
        let gamma = transcript.challenge_empty();
        let next_vars = WH_ROUND_NUM_VARIABLES[iteration] - k;
        let mut gamma_power_round = gamma;
        sigma += gamma_power_round * y_ood;
        terms[num_terms] = Term { coeff: gamma_power_round, power: z_ood };
        num_terms += 1;
        for s in 0..queries {
            let x_base = pow_index(omega, indices[s], index_bits);
            let folded = fold_coset(cosets[s], 1 / x_base, eta_inv, alphas, k, arity);
            gamma_power_round *= gamma;
            sigma += gamma_power_round * folded;
            // The folded oracle's domain point: omega^(idx * 2^k).
            let folded_point = pow_index(
                WH_ROUND_OMEGA[iteration],
                indices[s] << k,
                WH_ROUND_INDEX_BITS[iteration] + k,
            );
            terms[num_terms] = Term { coeff: gamma_power_round, power: folded_point };
            num_terms += 1;
        }
        assert(next_vars == WH_ROUND_NUM_VARIABLES[iteration + 1], "whir schedule");
        folded_root = next_root;
    }

    // Final phase: the clear polynomial, its agreement with the last oracle, and the claim.
    let _final_pow_seed = transcript.challenge(final_poly);
    let final_nonce = round_nonces[WH_NUM_ITERATIONS];
    let final_k = WH_ROUND_K[WH_NUM_ITERATIONS];
    let final_queries = WH_ROUND_QUERIES[WH_NUM_ITERATIONS];
    let final_index_bits = WH_ROUND_INDEX_BITS[WH_NUM_ITERATIONS];
    let final_omega = WH_ROUND_OMEGA[WH_NUM_ITERATIONS];
    let final_eta = 1 / WH_ROUND_ETA_INV[WH_NUM_ITERATIONS];

    let mut final_indices: [u32; WH_MAX_QUERIES] = [0; WH_MAX_QUERIES];
    for s in 0..final_queries {
        let challenge = if s == 0 {
            transcript.challenge_1(final_nonce)
        } else {
            transcript.challenge_empty()
        };
        final_indices[s] = index_from_challenge(challenge, final_index_bits);
    }
    let final_batch = NUM_GROUPS + WH_NUM_ITERATIONS - 1;
    let final_positions = authenticate_batch(
        final_batch,
        final_indices,
        final_queries,
        folded_root,
        opened_values,
        opened_siblings,
        hint_leaves,
        hint_num_leaves,
        hint_positions,
    );
    for s in 0..final_queries {
        let base = WH_BATCH_VALUE_OFFSET[final_batch]
            + final_positions[s] * WH_BATCH_LEAF_WIDTH[final_batch];
        let mut point = pow_index(final_omega, final_indices[s], final_index_bits);
        for t in 0..WH_BATCH_LEAF_WIDTH[final_batch] {
            assert(
                opened_values[base + t] == evaluate_coefficients(final_poly, point),
                "whir final oracle consistency",
            );
            point *= final_eta;
        }
    }

    // The accumulated claim: every weight term against the final polynomial.
    let mut total: Field = 0;
    let mut u_tail: [Field; WH_FINAL_VARS] = [0; WH_FINAL_VARS];
    for i in 0..WH_FINAL_VARS {
        u_tail[i] = u[WH_NUM_VARIABLES - WH_FINAL_VARS + i];
    }
    total += terms[0].coeff * multilinear_at(final_poly, u_tail);
    for t in 1..WH_NUM_TERMS {
        total += terms[t].coeff * evaluate_coefficients(final_poly, terms[t].power);
    }
    assert(total == sigma, "whir final claim");
    assert(num_terms == WH_NUM_TERMS, "term count");
    assert(bound == WH_NUM_VARIABLES - WH_FINAL_VARS, "bound variable count");
}

/// Authenticate one batch of Merkle openings and return each query's slot in the batch.
///
/// The distinct, ascending leaf set is a hint: the circuit checks it is strictly ascending and that
/// every query index appears at its claimed slot, which is all the batch walk needs. Extra entries
/// would only add authentication work, so no upper bound has to be enforced.
fn authenticate_batch(
    batch: u32,
    indices: [u32; WH_MAX_QUERIES],
    num_queries: u32,
    root: Field,
    opened_values: [Field; WH_TOTAL_VALUES],
    opened_siblings: [Field; WH_TOTAL_SIBLINGS],
    hint_leaves: [u32; WH_TOTAL_LEAVES],
    hint_num_leaves: [u32; WH_NUM_BATCHES],
    hint_positions: [u32; WH_TOTAL_POSITIONS],
) -> [u32; WH_MAX_QUERIES] {
    let max_leaves = WH_BATCH_MAX_LEAVES[batch];
    let width = WH_BATCH_LEAF_WIDTH[batch];
    let depth = WH_BATCH_DEPTH[batch];
    let num_leaves = hint_num_leaves[batch];

    let mut leaves: [u32; WH_MAX_BATCH_LEAVES] = [0; WH_MAX_BATCH_LEAVES];
    for i in 0..max_leaves {
        leaves[i] = hint_leaves[WH_BATCH_LEAF_OFFSET[batch] + i];
    }
    for i in 0..max_leaves {
        if i + 1 < num_leaves {
            assert(leaves[i] < leaves[i + 1], "batch leaves must ascend");
        }
    }

    let mut positions: [u32; WH_MAX_QUERIES] = [0; WH_MAX_QUERIES];
    for s in 0..num_queries {
        let slot = hint_positions[WH_BATCH_POSITION_OFFSET[batch] + s];
        assert(slot < num_leaves, "query slot out of range");
        assert(leaves[slot] == indices[s], "query slot mismatch");
        positions[s] = slot;
    }

    let mut values: [[Field; WH_MAX_LEAF_WIDTH]; WH_MAX_BATCH_LEAVES] =
        [[0; WH_MAX_LEAF_WIDTH]; WH_MAX_BATCH_LEAVES];
    for i in 0..max_leaves {
        for j in 0..width {
            values[i][j] = opened_values[WH_BATCH_VALUE_OFFSET[batch] + i * width + j];
        }
    }
    let mut siblings: [Field; WH_MAX_BATCH_LEAVES * WH_MAX_DEPTH] =
        [0; WH_MAX_BATCH_LEAVES * WH_MAX_DEPTH];
    for i in 0..(max_leaves * depth) {
        siblings[i] = opened_siblings[WH_BATCH_SIBLING_OFFSET[batch] + i];
    }

    let computed = merkle::batch_root::<WH_MAX_BATCH_LEAVES, WH_MAX_LEAF_WIDTH, WH_MAX_DEPTH, _>(
        leaves,
        num_leaves,
        values,
        width,
        siblings,
        depth,
    );
    assert(computed == root, "whir merkle batch");
    positions
}

/// `Fold(A, alpha)(x_base^(2^k))` from the 2^k coset values, README eq. (3.1) applied k times.
fn fold_coset<let A: u32, let K: u32>(
    values: [Field; A],
    x_base_inv: Field,
    eta_inv: Field,
    alphas: [Field; K],
    k: u32,
    arity: u32,
) -> Field {
    let mut current: [Field; A] = values;
    let mut size = arity;
    let mut base_inv = x_base_inv;
    let mut level_eta_inv = eta_inv;
    for level in 0..K {
        if level < k {
            let half = size / 2;
            let mut point_inv = base_inv;
            let mut next: [Field; A] = current;
            for t in 0..A {
                if t < half {
                    let v_pos = current[t];
                    let v_neg = current[t + half];
                    let even = (v_pos + v_neg) * TWO_INV;
                    let odd = (v_pos - v_neg) * TWO_INV * point_inv;
                    next[t] = even + alphas[level] * (odd - even);
                    point_inv *= level_eta_inv;
                }
            }
            current = next;
            base_inv *= base_inv;
            level_eta_inv *= level_eta_inv;
            size = half;
        }
    }
    current[0]
}

/// `base^exponent` for an exponent below `2^bits`.
fn pow_index(base: Field, exponent: u32, bits: u32) -> Field {
    let mut result: Field = 1;
    let mut square = base;
    for i in 0..bits {
        result *= if ((exponent >> i) & 1) == 1 {
            square
        } else {
            1
        };
        square *= square;
    }
    result
}

/// The univariate with coefficients `poly` at `point`.
fn evaluate_coefficients<let N: u32>(poly: [Field; N], point: Field) -> Field {
    let mut accumulator: Field = 0;
    for i in 0..N {
        accumulator = accumulator * point + poly[N - 1 - i];
    }
    accumulator
}

/// The multilinear extension of `poly` (evaluation form, low variable first) at `point`.
fn multilinear_at<let N: u32, let K: u32>(poly: [Field; N], point: [Field; K]) -> Field {
    let mut current: [Field; N] = poly;
    let mut size = N;
    for i in 0..K {
        let half = size / 2;
        let mut next: [Field; N] = current;
        for j in 0..N {
            if j < half {
                next[j] = current[2 * j] + point[i] * (current[2 * j + 1] - current[2 * j]);
            }
        }
        current = next;
        size = half;
    }
    current[0]
}
