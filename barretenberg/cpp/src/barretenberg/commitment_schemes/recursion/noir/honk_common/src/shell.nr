//! The transparent-UltraHonk shell: Fiat-Shamir schedule, sumcheck, and the claims handed to the
//! polynomial commitment scheme. Mirrors `honk_transparent::TransparentHonk::verify`.

use crate::poseidon;
use crate::relations;
use crate::relations::RelationParameters;
use crate::transcript::Transcript;

global IV_SHIFT: Field = 18446744073709551616; // 2^64
global PERMUTATION_ARGUMENT_VALUE_SEPARATOR: Field = 268435456; // 1 << 28
global BATCHED_LENGTH: u32 = 7;

/// The transparent verification key: circuit metadata plus the Merkle root of the committed
/// precomputed columns.
pub struct Vk {
    pub log_n: Field,
    pub num_public_inputs: Field,
    pub pub_inputs_offset: Field,
    pub virtual_mask: Field,
    pub lagrange_first_row: u32,
    pub lagrange_last_row: u32,
    pub root: Field,
}

/// Barycentric denominators d_i = prod_{j != i} (i - j) over the domain {0, ..., 6}, as
/// `BarycentricData<Fr, 7, 7>::lagrange_denominators` computes them.
global LAGRANGE_DENOMINATORS: [Field; BATCHED_LENGTH] = [720, -120, 48, -36, 48, -120, 720];

/// A degree-6 univariate given by its evaluations at 0..6, evaluated at `u`.
fn evaluate_univariate(values: [Field; BATCHED_LENGTH], u: Field) -> Field {
    let mut numerator: Field = 1;
    for i in 0..BATCHED_LENGTH {
        numerator *= u - i as Field;
    }
    let mut result: Field = 0;
    for i in 0..BATCHED_LENGTH {
        result += values[i] / (LAGRANGE_DENOMINATORS[i] * (u - i as Field));
    }
    result * numerator
}

/// The first Fiat-Shamir challenge, over the whole opening round buffer:
/// `Init || vk metadata || vk root || public inputs || wires root`. Spelled out as an index walk so
/// no concatenated array has to be materialized.
fn first_challenge<let NUM_PI: u32>(
    init: Field,
    vk: Vk,
    public_inputs: [Field; NUM_PI],
    wires_root: Field,
) -> Field {
    let total: u32 = NUM_PI + 9;
    let num_chunks: u32 = (total + 2) / 3;
    let mut state: [Field; 4] = [0, 0, 0, (total as Field) * IV_SHIFT];
    for c in 0..num_chunks {
        for j in 0..3 {
            let idx = c * 3 + j;
            if idx == 0 {
                state[j] += init;
            } else if idx == 1 {
                state[j] += vk.log_n;
            } else if idx == 2 {
                state[j] += vk.num_public_inputs;
            } else if idx == 3 {
                state[j] += vk.pub_inputs_offset;
            } else if idx == 4 {
                state[j] += vk.virtual_mask;
            } else if idx == 5 {
                state[j] += vk.lagrange_first_row as Field;
            } else if idx == 6 {
                state[j] += vk.lagrange_last_row as Field;
            } else if idx == 7 {
                state[j] += vk.root;
            } else if idx < total - 1 {
                state[j] += public_inputs[idx - 8];
            } else if idx == total - 1 {
                state[j] += wires_root;
            }
        }
        state = std::hash::poseidon2_permutation(state);
    }
    state[0]
}

/// The grand-product public-input correction (`compute_public_input_delta`).
fn public_input_delta<let NUM_PI: u32>(
    public_inputs: [Field; NUM_PI],
    beta: Field,
    gamma: Field,
    offset: Field,
) -> Field {
    let mut numerator: Field = 1;
    let mut denominator: Field = 1;
    let mut numerator_acc = gamma + beta * (PERMUTATION_ARGUMENT_VALUE_SEPARATOR + offset);
    let mut denominator_acc = gamma - beta * (offset + 1);
    for i in 0..NUM_PI {
        numerator *= numerator_acc + public_inputs[i];
        denominator *= denominator_acc + public_inputs[i];
        if i < NUM_PI - 1 {
            numerator_acc += beta;
            denominator_acc -= beta;
        }
    }
    numerator / denominator
}

/// eq(bits(row), u): the multilinear extension of a point indicator at the sumcheck challenge.
fn eq_at_row<let LOG_N: u32>(row: u32, u: [Field; LOG_N]) -> Field {
    let mut eq: Field = 1;
    for i in 0..LOG_N {
        eq *= if ((row >> i) & 1) == 1 {
            u[i]
        } else {
            1 - u[i]
        };
    }
    eq
}

pub struct ShellOutput<let LOG_N: u32> {
    pub challenge: [Field; LOG_N],
    pub transcript: Transcript,
}

/// Run the shell: absorb the commitments, derive the relation parameters, verify sumcheck, and
/// check the claimed evaluations of the columns that are never committed.
///
/// Leaves the transcript positioned exactly where the opening argument begins.
pub fn verify<let LOG_N: u32, let NUM_PI: u32>(
    init: Field,
    vk: Vk,
    public_inputs: [Field; NUM_PI],
    witness_roots: [Field; 3],
    univariates: [[Field; BATCHED_LENGTH]; LOG_N],
    evaluations: [Field; relations::NUM_ALL_ENTITIES],
    is_virtual: [bool; 24],
) -> ShellOutput<LOG_N> {
    let mut transcript =
        Transcript { previous: first_challenge(init, vk, public_inputs, witness_roots[0]) };
    let eta = transcript.previous;
    let rom_logup_gamma = transcript.challenge_empty();
    let beta = transcript.challenge_1(witness_roots[1]);
    let gamma = transcript.challenge_empty();
    let alpha = transcript.challenge_1(witness_roots[2]);
    let gate_challenge = transcript.challenge_empty();

    let params = RelationParameters {
        eta,
        eta_two: eta * eta,
        eta_three: eta * eta * eta,
        rom_logup_gamma,
        beta,
        beta_sqr: beta * beta,
        beta_cube: beta * beta * beta,
        gamma,
        public_input_delta: public_input_delta(public_inputs, beta, gamma, vk.pub_inputs_offset),
    };

    // Sumcheck: check each round univariate against the running target, then fold it.
    let mut target: Field = 0;
    let mut challenge: [Field; LOG_N] = [0; LOG_N];
    for round in 0..LOG_N {
        let univariate = univariates[round];
        let u = transcript.challenge(univariate);
        assert(univariate[0] + univariate[1] == target, "sumcheck round sum");
        target = evaluate_univariate(univariate, u);
        challenge[round] = u;
    }
    // The claimed evaluations stay in the round buffer: the opening argument's first challenge is
    // what absorbs them, exactly as the native verifier does.

    // pow_beta at the challenge: prod_k ((1 - u_k) + u_k * beta_k), beta_k the dyadic powers.
    let mut pow: Field = 1;
    let mut beta_k = gate_challenge;
    for k in 0..LOG_N {
        pow *= (1 - challenge[k]) + challenge[k] * beta_k;
        beta_k *= beta_k;
    }

    let subrelations = relations::accumulate(evaluations, params, pow);
    assert(relations::batch(subrelations, alpha) == target, "sumcheck final check");

    // Columns that are never committed: their claimed evaluation must equal the value the verifier
    // computes itself - zero for an identically-zero column, eq(row, u) for a point indicator.
    for e in 0..24 {
        if is_virtual[e] {
            let expected = if e == relations::LAGRANGE_FIRST {
                eq_at_row(vk.lagrange_first_row, challenge)
            } else if e == relations::LAGRANGE_LAST {
                eq_at_row(vk.lagrange_last_row, challenge)
            } else {
                0
            };
            assert(evaluations[e] == expected, "virtual column evaluation");
        }
    }

    ShellOutput { challenge, transcript }
}
