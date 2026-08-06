barretenberg_module(commitment_schemes common transcript polynomials ecc numeric srs crypto_blake3s crypto_poseidon2 crypto_skyscraper)

# pcs_acir_bench (from pcs_acir.bench.cpp) drives an ACIR circuit through every PCS backend; it
# needs the dsl frontend (and its vm2 stub), which the commitment_schemes library itself must not
# depend on.
if(TARGET pcs_acir_bench)
    target_link_libraries(pcs_acir_bench PRIVATE dsl ultra_honk stdlib_honk_verifier vm2_stub)
endif()