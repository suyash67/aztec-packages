//! The twenty subrelations of `UltraProveKitFlavor`, evaluated at the sumcheck challenge.
//!
//! Ported one-for-one from `barretenberg/relations/{permutation, logderiv_lookup, ultra_arithmetic,
//! delta_range_constraint, memory}_relation.hpp`. The verifier evaluates each relation at a single
//! point rather than as a univariate, so every accumulator is one field element and the prover-only
//! fast paths (which branch on witness values) are absent.
//!
//! Subrelation 1 of the lookup relation and subrelation 7 of the memory relation are linearly
//! dependent: they assert a sum over the whole trace rather than a per-row identity, so they are not
//! multiplied by the `pow_beta` scaling factor.

pub global NUM_ALL_ENTITIES: u32 = 37;
pub global NUM_SUBRELATIONS: u32 = 20;

// `AllEntities` indices, in EntityId order.
pub global SIGMA_1: u32 = 0;
pub global SIGMA_2: u32 = 1;
pub global SIGMA_3: u32 = 2;
pub global SIGMA_4: u32 = 3;
pub global ID_1: u32 = 4;
pub global ID_2: u32 = 5;
pub global ID_3: u32 = 6;
pub global ID_4: u32 = 7;
pub global LAGRANGE_FIRST: u32 = 8;
pub global LAGRANGE_LAST: u32 = 9;
pub global Q_LOOKUP: u32 = 10;
pub global TABLE_1: u32 = 11;
pub global TABLE_2: u32 = 12;
pub global TABLE_3: u32 = 13;
pub global TABLE_4: u32 = 14;
pub global Q_M: u32 = 15;
pub global Q_R: u32 = 16;
pub global Q_O: u32 = 17;
pub global Q_C: u32 = 18;
pub global Q_L: u32 = 19;
pub global Q_4: u32 = 20;
pub global Q_ARITH: u32 = 21;
pub global Q_DELTA_RANGE: u32 = 22;
pub global Q_MEMORY: u32 = 23;
pub global W_L: u32 = 24;
pub global W_R: u32 = 25;
pub global W_O: u32 = 26;
pub global W_4: u32 = 27;
pub global Z_PERM: u32 = 28;
pub global LOOKUP_INVERSES: u32 = 29;
pub global LOOKUP_READ_COUNTS: u32 = 30;
pub global LOOKUP_READ_TAGS: u32 = 31;
pub global W_L_SHIFT: u32 = 32;
pub global W_R_SHIFT: u32 = 33;
pub global W_O_SHIFT: u32 = 34;
pub global W_4_SHIFT: u32 = 35;
pub global Z_PERM_SHIFT: u32 = 36;

/// (-2)^-1 over BN254 Fr, the `neg_half` constant of the primary arithmetic subrelation.
global NEG_HALF: Field =
    10944121435919637611123202872628637544274182200208017171849102093287904247808;

pub struct RelationParameters {
    pub eta: Field,
    pub eta_two: Field,
    pub eta_three: Field,
    pub rom_logup_gamma: Field,
    pub beta: Field,
    pub beta_sqr: Field,
    pub beta_cube: Field,
    pub gamma: Field,
    pub public_input_delta: Field,
}

/// Every subrelation's value at the sumcheck challenge, in flavor tuple order:
/// permutation (3), lookup (3), arithmetic (2), delta range (4), memory (8).
pub fn accumulate(
    p: [Field; NUM_ALL_ENTITIES],
    params: RelationParameters,
    pow: Field,
) -> [Field; NUM_SUBRELATIONS] {
    let mut out: [Field; NUM_SUBRELATIONS] = [0; NUM_SUBRELATIONS];

    // --- UltraPermutationRelation ---------------------------------------------------------
    let w_1 = p[W_L];
    let w_2 = p[W_R];
    let w_3 = p[W_O];
    let w_4 = p[W_4];
    let gamma = params.gamma;
    let beta = params.beta;

    let w_1_plus_gamma = w_1 + gamma;
    let w_2_plus_gamma = w_2 + gamma;
    let w_3_plus_gamma = w_3 + gamma;
    let w_4_plus_gamma = w_4 + gamma;

    let t1 = (p[ID_1] * beta + w_1_plus_gamma) * pow;
    let t2 = p[ID_2] * beta + w_2_plus_gamma;
    let t3 = p[ID_3] * beta + w_3_plus_gamma;
    let t4 = p[ID_4] * beta + w_4_plus_gamma;
    let t5 = (p[SIGMA_1] * beta + w_1_plus_gamma) * pow;
    let t6 = p[SIGMA_2] * beta + w_2_plus_gamma;
    let t7 = p[SIGMA_3] * beta + w_3_plus_gamma;
    let t8 = p[SIGMA_4] * beta + w_4_plus_gamma;

    let numerator = t1 * t2 * t3 * t4;
    let denominator = t5 * t6 * t7 * t8;
    let lagrange_first = p[LAGRANGE_FIRST];
    let lagrange_last = p[LAGRANGE_LAST];
    let public_input_term = lagrange_last * params.public_input_delta + p[Z_PERM_SHIFT];
    out[0] = (p[Z_PERM] + lagrange_first) * numerator - public_input_term * denominator;
    out[1] = lagrange_last * p[Z_PERM_SHIFT] * pow;
    out[2] = lagrange_first * p[Z_PERM] * pow;

    // --- LogDerivLookupRelation -----------------------------------------------------------
    let inverses = p[LOOKUP_INVERSES];
    let read_counts = p[LOOKUP_READ_COUNTS];
    let read_selector = p[Q_LOOKUP];
    let read_tag = p[LOOKUP_READ_TAGS];
    let beta_sqr = params.beta_sqr;
    let beta_cube = params.beta_cube;

    let inverse_exists = -(read_tag * read_selector) + read_tag + read_selector;
    let table_term = p[TABLE_2] * beta + p[TABLE_3] * beta_sqr + p[TABLE_4] * beta_cube + p[TABLE_1]
        + gamma;
    let derived_1 = p[Q_R] * p[W_L_SHIFT] + (w_1 + gamma);
    let derived_2 = p[Q_M] * p[W_R_SHIFT] + w_2;
    let derived_3 = p[Q_C] * p[W_O_SHIFT] + w_3;
    let lookup_term = derived_2 * beta + derived_3 * beta_sqr + derived_1 + p[Q_O] * beta_cube;

    out[3] = (lookup_term * table_term * inverses - inverse_exists) * pow;
    out[4] = (read_selector * table_term - read_counts * lookup_term) * inverses;
    out[5] = (read_tag * read_tag - read_tag) * pow;

    // --- ArithmeticRelation ---------------------------------------------------------------
    let q_arith = p[Q_ARITH];
    let q_m = p[Q_M];
    let q_arith_sub_1 = q_arith - 1;
    let scaled_q_arith = q_arith * pow;

    let tmp0 = (w_2 * w_1 * NEG_HALF) * ((q_arith - 3) * q_m);
    let mut tmp1 = p[Q_L] * w_1 + p[Q_R] * w_2 + p[Q_O] * w_3 + p[Q_4] * w_4 + p[Q_C];
    tmp1 += q_arith_sub_1 * p[W_4_SHIFT];
    out[6] = (tmp0 + tmp1) * scaled_q_arith;

    let secondary = w_1 + w_4 - p[W_L_SHIFT] + q_m;
    out[7] = (secondary * (q_arith - 2)) * (q_arith_sub_1 * scaled_q_arith);

    // --- DeltaRangeConstraintRelation -----------------------------------------------------
    let q_delta_range_scaled = p[Q_DELTA_RANGE] * pow;
    let delta_1 = w_2 - w_1;
    let delta_2 = w_3 - w_2;
    let delta_3 = w_4 - w_3;
    let delta_4 = p[W_L_SHIFT] - w_4;
    let deltas = [delta_1, delta_2, delta_3, delta_4];
    for i in 0..4 {
        let t = (deltas[i] - 3) * deltas[i];
        out[8 + i] = t * (t + 2) * q_delta_range_scaled;
    }

    // --- MemoryRelation -------------------------------------------------------------------
    let eta = params.eta;
    let eta_two = params.eta_two;
    let eta_three = params.eta_three;
    let w_1_shift = p[W_L_SHIFT];
    let w_2_shift = p[W_R_SHIFT];
    let w_3_shift = p[W_O_SHIFT];
    let w_4_shift = p[W_4_SHIFT];
    let q_1 = p[Q_L];
    let q_2 = p[Q_R];
    let q_3 = p[Q_O];
    let q_4_sel = p[Q_4];
    let q_c = p[Q_C];
    let q_memory = p[Q_MEMORY];

    let partial_record_check = w_3 * eta_three + w_2 * eta_two + w_1 * eta + q_c;
    let memory_record_check = partial_record_check - w_4;

    let neg_index_delta = w_1 - w_1_shift;
    let index_delta_is_zero = neg_index_delta + 1;
    let record_delta = w_4_shift - w_4;
    let index_increases_by_zero_or_one =
        neg_index_delta * neg_index_delta + neg_index_delta;
    let adjacent_values_match = index_delta_is_zero * record_delta;

    let q_memory_by_scaling = q_memory * pow;
    let q_one_by_two = q_1 * q_2;
    let q_one_by_two_by_memory_by_scaling = q_one_by_two * q_memory_by_scaling;
    out[13] = adjacent_values_match * q_one_by_two_by_memory_by_scaling;
    out[14] = index_increases_by_zero_or_one * q_one_by_two_by_memory_by_scaling;

    let rom_consistency_check_identity = memory_record_check * q_one_by_two;

    let neg_access_type = partial_record_check - w_4;
    let access_check = neg_access_type * neg_access_type + neg_access_type;
    let neg_next_gate_access_type =
        w_3_shift * eta_three + w_2_shift * eta_two + w_1_shift * eta - w_4_shift;
    let value_delta = w_3_shift - w_3;
    let adjacent_values_match_and_next_is_read =
        (index_delta_is_zero * value_delta) * (neg_next_gate_access_type + 1);
    let next_gate_access_type_is_boolean =
        neg_next_gate_access_type * neg_next_gate_access_type + neg_next_gate_access_type;

    let q_3_by_memory_and_scaling = q_3 * q_memory_by_scaling;
    out[15] = adjacent_values_match_and_next_is_read * q_3_by_memory_and_scaling;
    out[16] = index_increases_by_zero_or_one * q_3_by_memory_and_scaling;
    out[17] = next_gate_access_type_is_boolean * q_3_by_memory_and_scaling;

    let ram_consistency_check_identity = access_check * q_3_by_memory_and_scaling;

    let timestamp_delta = w_2_shift - w_2;
    let ram_timestamp_check_identity = index_delta_is_zero * timestamp_delta - w_3;

    let mut memory_identity = rom_consistency_check_identity;
    memory_identity += ram_timestamp_check_identity * (q_4_sel * q_1);
    memory_identity += memory_record_check * (q_m * q_1);
    memory_identity *= q_memory_by_scaling;
    memory_identity += ram_consistency_check_identity;
    out[12] = memory_identity;

    let one_minus_q_1 = -q_1 + 1;
    let q_logup_table = q_2 * one_minus_q_1;
    let q_logup_read = q_4_sel * one_minus_q_1;
    let q_logup_any = q_logup_table + q_logup_read;
    let denom = w_1 + w_2 * eta + q_c * eta_two + params.rom_logup_gamma;

    out[18] = q_memory_by_scaling * q_logup_any * (w_4 * denom - 1);
    out[19] = (q_logup_read - q_logup_table * w_3) * w_4 * q_memory;

    out
}

/// The batched relation value: subrelation 0 carries coefficient 1, subrelation i coefficient
/// `alpha^i` (`initialize_relation_separator` followed by `scale_and_batch_elements`).
pub fn batch(evaluations: [Field; NUM_SUBRELATIONS], alpha: Field) -> Field {
    let mut result = evaluations[0];
    let mut power = alpha;
    for i in 1..NUM_SUBRELATIONS {
        result += evaluations[i] * power;
        power *= alpha;
    }
    result
}
