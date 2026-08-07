//! C ABI over the reference WHIR implementation (WizardOfMenlo/whir at the revision ProveKit
//! pins), so barretenberg's C++ WHIR can be differentially tested against it in-process and so a
//! Rust-backed PCS variant can be prototyped without reimplementing the protocol.
//!
//! Every field element crosses the boundary as 32 canonical little-endian bytes of BN254 Fr, the
//! same representation `bb::fr`'s `uint256_t` conversion produces. Arrays are in the little-endian
//! hypercube order both implementations use: element i is the multilinear extension's value at
//! `bits(i)` with bit 0 the first variable.
//!
//! # Opening-point variable order
//!
//! The two implementations order the *point* oppositely, and this is the one convention that does
//! not line up on its own. The reference's `MultilinearExtension` binds `point[0]` to the **most
//! significant** index bit — `eval_eq` splits the accumulator's high half on `point[0]`, and
//! `mixed_multilinear_extend` interpolates the last coordinate over adjacent entries. Barretenberg
//! binds `point[0]` to index bit 0, the variable its stride-2 fold consumes first. Points therefore
//! cross this boundary **reversed**.
//!
//! Every entry point here takes and returns points in barretenberg's order and reverses internally,
//! so C++ callers never have to think about it. Nothing else needs adjusting: the committed array
//! itself has the same meaning on both sides.

use std::{borrow::Cow, panic, ptr, slice, sync::Arc};

use ark_bn254::Fr;
use ark_ff::{BigInteger, PrimeField};
use whir::{
    algebra::{
        embedding::Identity,
        linear_form::{Evaluate, LinearForm, MultilinearExtension},
        ntt::{NttEngine, ReedSolomon, NTT},
    },
    hash,
    parameters::ProtocolParameters,
    protocols::whir::Config as WhirConfig,
    transcript::{codecs::Empty, DomainSeparator, Proof, ProverState, VerifierState},
};

/// Fiat-Shamir session label; prover and verifier must agree, nothing else depends on it.
const SESSION: &str = "barretenberg whir differential harness";

pub const WHIR_RS_OK: i32 = 0;
pub const WHIR_RS_ERR_ARGS: i32 = -1;
pub const WHIR_RS_ERR_ENCODING: i32 = -2;
pub const WHIR_RS_ERR_VERIFY: i32 = -3;
pub const WHIR_RS_ERR_PANIC: i32 = -4;

/// Protocol parameters, mirroring `whir::parameters::ProtocolParameters`. The defaults ProveKit
/// uses are security_level 128, pow_bits 10, both folding factors 3, log_inv_rate 2 and
/// `unique_decoding` false (the Johnson bound).
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct WhirRsParams {
    pub num_variables: u32,
    pub security_level: u32,
    pub pow_bits: u32,
    pub initial_folding_factor: u32,
    pub folding_factor: u32,
    pub starting_log_inv_rate: u32,
    /// Non-zero selects unique decoding; zero selects list decoding at the Johnson bound.
    pub unique_decoding: u32,
}

/// BN254 is not in the reference crate's built-in NTT registry (ProveKit registers its own
/// hand-tuned engine). The crate's generic Cooley-Tukey engine works for any `FftField`, so
/// registering that keeps this harness self-contained.
fn ensure_bn254_registered() {
    use std::sync::Once;
    static ONCE: Once = Once::new();
    ONCE.call_once(|| {
        let engine: Arc<dyn ReedSolomon<Fr>> = Arc::new(NttEngine::<Fr>::new_from_fftfield());
        NTT.insert(engine);
    });
}

fn build_config(params: &WhirRsParams) -> WhirConfig<Identity<Fr>> {
    ensure_bn254_registered();
    let proto = ProtocolParameters {
        unique_decoding: params.unique_decoding != 0,
        security_level: params.security_level as usize,
        pow_bits: params.pow_bits as usize,
        initial_folding_factor: params.initial_folding_factor as usize,
        folding_factor: params.folding_factor as usize,
        starting_log_inv_rate: params.starting_log_inv_rate as usize,
        batch_size: 1,
        hash_id: hash::BLAKE3,
    };
    WhirConfig::<Identity<Fr>>::new(1 << params.num_variables, &proto)
}

/// Decode `count` field elements from 32-byte little-endian canonical limbs.
fn read_field_elements(ptr: *const u8, count: usize) -> Option<Vec<Fr>> {
    if count == 0 {
        return Some(Vec::new());
    }
    if ptr.is_null() {
        return None;
    }
    let bytes = unsafe { slice::from_raw_parts(ptr, count * 32) };
    bytes
        .chunks_exact(32)
        .map(|chunk| {
            let value = Fr::from_le_bytes_mod_order(chunk);
            // Reject non-canonical encodings rather than silently reducing them, so a mismatch in
            // the C++ encoder surfaces here instead of as a wrong evaluation.
            let mut round_trip = value.into_bigint().to_bytes_le();
            round_trip.resize(32, 0);
            (round_trip == chunk).then_some(value)
        })
        .collect()
}

/// Read a point given in barretenberg's variable order and return it in the reference's.
/// See the module docs: `point[0]` is barretenberg's first-folded variable but the reference's
/// most significant one.
fn read_point(ptr: *const u8, count: usize) -> Option<Vec<Fr>> {
    let mut point = read_field_elements(ptr, count)?;
    point.reverse();
    Some(point)
}

fn write_field_element(value: &Fr, out: *mut u8) {
    let mut bytes = value.into_bigint().to_bytes_le();
    bytes.resize(32, 0);
    unsafe { ptr::copy_nonoverlapping(bytes.as_ptr(), out, 32) };
}

/// Wire format for a `Proof`: an 8-byte little-endian length for the narg string, then the narg
/// string, then the hints.
fn proof_to_bytes(proof: &Proof) -> Vec<u8> {
    let mut out = Vec::with_capacity(8 + proof.narg_string.len() + proof.hints.len());
    out.extend_from_slice(&(proof.narg_string.len() as u64).to_le_bytes());
    out.extend_from_slice(&proof.narg_string);
    out.extend_from_slice(&proof.hints);
    out
}

fn proof_from_bytes(bytes: &[u8]) -> Option<Proof> {
    if bytes.len() < 8 {
        return None;
    }
    let narg_len = u64::from_le_bytes(bytes[0..8].try_into().ok()?) as usize;
    if bytes.len() < 8 + narg_len {
        return None;
    }
    Some(Proof {
        narg_string: bytes[8..8 + narg_len].to_vec(),
        hints: bytes[8 + narg_len..].to_vec(),
        #[cfg(debug_assertions)]
        pattern: Vec::new(),
    })
}

/// Multilinear extension of `coeffs` at `point`, in the reference implementation's convention.
///
/// # Safety
/// `coeffs` must point to `2^num_variables` 32-byte elements, `point` to `num_variables` of them,
/// and `out_value` to 32 writable bytes.
#[no_mangle]
pub unsafe extern "C" fn whir_rs_mle_evaluate(
    num_variables: u32,
    coeffs: *const u8,
    point: *const u8,
    out_value: *mut u8,
) -> i32 {
    let result = panic::catch_unwind(|| {
        let Some(vector) = read_field_elements(coeffs, 1 << num_variables) else {
            return WHIR_RS_ERR_ENCODING;
        };
        let Some(point) = read_point(point, num_variables as usize) else {
            return WHIR_RS_ERR_ENCODING;
        };
        if out_value.is_null() {
            return WHIR_RS_ERR_ARGS;
        }
        let form = MultilinearExtension { point };
        write_field_element(&form.evaluate(&Identity::<Fr>::new(), &vector), out_value);
        WHIR_RS_OK
    });
    result.unwrap_or(WHIR_RS_ERR_PANIC)
}

/// Commit to `coeffs` and prove its multilinear opening at `point`.
///
/// Writes the claimed evaluation to `out_value` and allocates the proof, whose pointer and length
/// go to `out_proof` / `out_proof_len`. Release it with `whir_rs_free`.
///
/// # Safety
/// Pointer sizes are as in `whir_rs_mle_evaluate`; `out_proof` and `out_proof_len` must be
/// writable.
#[no_mangle]
pub unsafe extern "C" fn whir_rs_prove(
    params: *const WhirRsParams,
    coeffs: *const u8,
    point: *const u8,
    out_value: *mut u8,
    out_proof: *mut *mut u8,
    out_proof_len: *mut usize,
) -> i32 {
    let result = panic::catch_unwind(|| {
        if params.is_null() || out_value.is_null() || out_proof.is_null() || out_proof_len.is_null()
        {
            return WHIR_RS_ERR_ARGS;
        }
        let params = unsafe { *params };
        let Some(vector) = read_field_elements(coeffs, 1 << params.num_variables) else {
            return WHIR_RS_ERR_ENCODING;
        };
        let Some(point) = read_point(point, params.num_variables as usize) else {
            return WHIR_RS_ERR_ENCODING;
        };

        let config = build_config(&params);
        let evaluation = MultilinearExtension {
            point: point.clone(),
        }
        .evaluate(config.embedding(), &vector);

        let ds = DomainSeparator::protocol(&config)
            .session(&SESSION)
            .instance(&Empty);
        let mut prover_state = ProverState::new_std(&ds);
        let witness = config.commit(&mut prover_state, &[&vector]);
        let forms: Vec<Box<dyn LinearForm<Fr>>> = vec![Box::new(MultilinearExtension { point })];
        let _ = config.prove(
            &mut prover_state,
            vec![Cow::Owned(vector)],
            vec![Cow::Owned(witness)],
            forms,
            Cow::Owned(vec![evaluation]),
        );

        let bytes = proof_to_bytes(&prover_state.proof());
        write_field_element(&evaluation, out_value);
        let mut boxed = bytes.into_boxed_slice();
        unsafe {
            *out_proof_len = boxed.len();
            *out_proof = boxed.as_mut_ptr();
        }
        std::mem::forget(boxed);
        WHIR_RS_OK
    });
    result.unwrap_or(WHIR_RS_ERR_PANIC)
}

/// Verify a proof produced by `whir_rs_prove`. Returns `WHIR_RS_OK` on acceptance and
/// `WHIR_RS_ERR_VERIFY` on rejection.
///
/// # Safety
/// `proof` must point to `proof_len` readable bytes; the other pointers are as above.
#[no_mangle]
pub unsafe extern "C" fn whir_rs_verify(
    params: *const WhirRsParams,
    point: *const u8,
    value: *const u8,
    proof: *const u8,
    proof_len: usize,
) -> i32 {
    let result = panic::catch_unwind(|| {
        if params.is_null() || proof.is_null() {
            return WHIR_RS_ERR_ARGS;
        }
        let params = unsafe { *params };
        let Some(point) = read_point(point, params.num_variables as usize) else {
            return WHIR_RS_ERR_ENCODING;
        };
        let Some(value) = read_field_elements(value, 1) else {
            return WHIR_RS_ERR_ENCODING;
        };
        let Some(proof) = proof_from_bytes(unsafe { slice::from_raw_parts(proof, proof_len) })
        else {
            return WHIR_RS_ERR_ENCODING;
        };

        let config = build_config(&params);
        let ds = DomainSeparator::protocol(&config)
            .session(&SESSION)
            .instance(&Empty);
        let mut verifier_state = VerifierState::new_std(&ds, &proof);
        let Ok(commitment) = config.receive_commitment(&mut verifier_state) else {
            return WHIR_RS_ERR_VERIFY;
        };
        let Ok(final_claim) = config.verify(&mut verifier_state, &[&commitment], &value) else {
            return WHIR_RS_ERR_VERIFY;
        };
        let form = MultilinearExtension { point };
        let forms: Vec<&dyn LinearForm<Fr>> = vec![&form];
        if final_claim.verify(forms).is_err() {
            return WHIR_RS_ERR_VERIFY;
        }
        WHIR_RS_OK
    });
    result.unwrap_or(WHIR_RS_ERR_PANIC)
}

/// Release a proof buffer returned by `whir_rs_prove`.
///
/// # Safety
/// `ptr`/`len` must be exactly what `whir_rs_prove` wrote, and freed at most once.
#[no_mangle]
pub unsafe extern "C" fn whir_rs_free(ptr: *mut u8, len: usize) {
    if ptr.is_null() {
        return;
    }
    drop(unsafe { Box::from_raw(slice::from_raw_parts_mut(ptr, len)) });
}

/// Total in-domain query count of the schedule these parameters produce, so the C++ side can
/// assert its own schedule matches without duplicating the derivation.
#[no_mangle]
pub unsafe extern "C" fn whir_rs_total_queries(params: *const WhirRsParams) -> i64 {
    let result = panic::catch_unwind(|| {
        if params.is_null() {
            return i64::from(WHIR_RS_ERR_ARGS);
        }
        let config = build_config(unsafe { &*params });
        let mut total = config.initial_committer.in_domain_samples;
        for round in &config.round_configs {
            total += round.irs_committer.in_domain_samples;
        }
        total as i64
    });
    result.unwrap_or(i64::from(WHIR_RS_ERR_PANIC))
}

/// Out-of-domain samples the reference takes against each committed oracle.
#[no_mangle]
pub unsafe extern "C" fn whir_rs_ood_samples(params: *const WhirRsParams) -> i64 {
    let result = panic::catch_unwind(|| {
        if params.is_null() {
            return i64::from(WHIR_RS_ERR_ARGS);
        }
        let config = build_config(unsafe { &*params });
        config.initial_committer.out_domain_samples as i64
    });
    result.unwrap_or(i64::from(WHIR_RS_ERR_PANIC))
}

#[cfg(test)]
mod tests {
    use ark_ff::{AdditiveGroup, Field};
    use ark_std::{
        rand::{rngs::StdRng, SeedableRng},
        UniformRand,
    };

    use super::*;

    fn params(num_variables: u32) -> WhirRsParams {
        WhirRsParams {
            num_variables,
            security_level: 32,
            pow_bits: 0,
            initial_folding_factor: 3,
            folding_factor: 3,
            starting_log_inv_rate: 2,
            unique_decoding: 0,
        }
    }

    fn encode(values: &[Fr]) -> Vec<u8> {
        let mut out = vec![0u8; values.len() * 32];
        for (i, value) in values.iter().enumerate() {
            write_field_element(value, unsafe { out.as_mut_ptr().add(i * 32) });
        }
        out
    }

    #[test]
    fn prove_verify_round_trip() {
        let mut rng = StdRng::seed_from_u64(7);
        let params = params(10);
        let vector: Vec<Fr> = (0..1 << params.num_variables)
            .map(|_| Fr::rand(&mut rng))
            .collect();
        let point: Vec<Fr> = (0..params.num_variables)
            .map(|_| Fr::rand(&mut rng))
            .collect();
        let coeffs = encode(&vector);
        let point_bytes = encode(&point);

        let mut value = [0u8; 32];
        let mut proof: *mut u8 = ptr::null_mut();
        let mut proof_len = 0usize;
        let status = unsafe {
            whir_rs_prove(
                &params,
                coeffs.as_ptr(),
                point_bytes.as_ptr(),
                value.as_mut_ptr(),
                &mut proof,
                &mut proof_len,
            )
        };
        assert_eq!(status, WHIR_RS_OK);

        // The claimed evaluation must be the multilinear extension at the point.
        let mut expected = [0u8; 32];
        assert_eq!(
            unsafe {
                whir_rs_mle_evaluate(
                    params.num_variables,
                    coeffs.as_ptr(),
                    point_bytes.as_ptr(),
                    expected.as_mut_ptr(),
                )
            },
            WHIR_RS_OK
        );
        assert_eq!(value, expected);

        assert_eq!(
            unsafe {
                whir_rs_verify(
                    &params,
                    point_bytes.as_ptr(),
                    value.as_ptr(),
                    proof,
                    proof_len,
                )
            },
            WHIR_RS_OK
        );

        // A wrong claimed value must be rejected.
        let mut wrong = value;
        wrong[0] ^= 1;
        assert_eq!(
            unsafe {
                whir_rs_verify(
                    &params,
                    point_bytes.as_ptr(),
                    wrong.as_ptr(),
                    proof,
                    proof_len,
                )
            },
            WHIR_RS_ERR_VERIFY
        );

        unsafe { whir_rs_free(proof, proof_len) };
    }

    // Pins the point reversal on a case that can be read off by hand. With the array holding
    // f(bits(i)) little-endian, the corner (x0, x1) = (1, 0) is index 1, so a point that is 1 in
    // barretenberg's first variable and 0 in its second must select array[1].
    #[test]
    fn point_order_follows_barretenberg() {
        let array = [Fr::from(10u64), Fr::from(20u64), Fr::from(30u64), Fr::from(40u64)];
        let coeffs = encode(&array);
        let point = encode(&[Fr::ONE, Fr::ZERO]);
        let mut value = [0u8; 32];
        assert_eq!(
            unsafe { whir_rs_mle_evaluate(2, coeffs.as_ptr(), point.as_ptr(), value.as_mut_ptr()) },
            WHIR_RS_OK
        );
        assert_eq!(Fr::from_le_bytes_mod_order(&value), array[1]);

        // And the opposite corner selects array[2], which is what would break if the reversal were
        // dropped: without it this same call would return array[2] and the one above array[1].
        let point = encode(&[Fr::ZERO, Fr::ONE]);
        assert_eq!(
            unsafe { whir_rs_mle_evaluate(2, coeffs.as_ptr(), point.as_ptr(), value.as_mut_ptr()) },
            WHIR_RS_OK
        );
        assert_eq!(Fr::from_le_bytes_mod_order(&value), array[2]);
    }

    #[test]
    fn non_canonical_encodings_are_rejected() {
        // The BN254 modulus itself is not a canonical element encoding.
        let modulus_bytes = {
            let mut bytes = Fr::MODULUS.to_bytes_le();
            bytes.resize(32, 0);
            bytes
        };
        assert!(read_field_elements(modulus_bytes.as_ptr(), 1).is_none());
    }

    #[test]
    fn schedule_matches_provekit_parameters() {
        // ProveKit's own parameters at 2^19: 128-bit security with 10 bits of grinding.
        let mut provekit = params(19);
        provekit.security_level = 128;
        provekit.pow_bits = 10;
        assert_eq!(unsafe { whir_rs_total_queries(&provekit) }, 305);
        assert_eq!(unsafe { whir_rs_ood_samples(&provekit) }, 1);

        // With grinding disabled every bit comes from queries, which is barretenberg's situation.
        provekit.pow_bits = 0;
        assert_eq!(unsafe { whir_rs_total_queries(&provekit) }, 330);
    }
}
