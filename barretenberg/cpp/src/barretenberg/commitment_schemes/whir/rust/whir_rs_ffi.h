#pragma once

// C ABI exported by the `whir-rs-ffi` crate in this directory, which wraps the reference WHIR
// implementation (WizardOfMenlo/whir at the revision ProveKit pins). Field elements are 32
// canonical little-endian bytes of BN254 Fr — the same layout `uint256_t(bb::fr)` produces.
//
// Build the library with `cargo build --release` here; it lands at
// `rust/target/release/libwhir_rs_ffi.a`.

#include <cstddef>
#include <cstdint>

extern "C" {

enum WhirRsStatus : int32_t {
    WHIR_RS_OK = 0,
    WHIR_RS_ERR_ARGS = -1,
    WHIR_RS_ERR_ENCODING = -2,
    WHIR_RS_ERR_VERIFY = -3,
    WHIR_RS_ERR_PANIC = -4,
};

struct WhirRsParams {
    uint32_t num_variables;
    uint32_t security_level;
    uint32_t pow_bits;
    uint32_t initial_folding_factor;
    uint32_t folding_factor;
    uint32_t starting_log_inv_rate;
    uint32_t unique_decoding; // non-zero selects unique decoding, zero the Johnson bound
};

int32_t whir_rs_mle_evaluate(uint32_t num_variables, const uint8_t* coeffs, const uint8_t* point, uint8_t* out_value);

int32_t whir_rs_prove(const WhirRsParams* params,
                      const uint8_t* coeffs,
                      const uint8_t* point,
                      uint8_t* out_value,
                      uint8_t** out_proof,
                      size_t* out_proof_len);

int32_t whir_rs_verify(
    const WhirRsParams* params, const uint8_t* point, const uint8_t* value, const uint8_t* proof, size_t proof_len);

void whir_rs_free(uint8_t* ptr, size_t len);

int64_t whir_rs_total_queries(const WhirRsParams* params);
int64_t whir_rs_ood_samples(const WhirRsParams* params);
}
