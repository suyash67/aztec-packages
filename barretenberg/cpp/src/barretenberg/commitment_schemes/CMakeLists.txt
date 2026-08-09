barretenberg_module(commitment_schemes common transcript polynomials ecc numeric srs crypto_blake3s crypto_poseidon2 crypto_sha256 crypto_skyscraper stdlib_primitives stdlib_poseidon2 stdlib_blake3s circuit_checker)

# pcs_acir_bench (from pcs_acir.bench.cpp) drives an ACIR circuit through every PCS backend; it
# needs the dsl frontend (and its vm2 stub), which the commitment_schemes library itself must not
# depend on.
if(TARGET pcs_acir_bench)
    target_link_libraries(pcs_acir_bench PRIVATE dsl ultra_honk stdlib_honk_verifier vm2_stub)
endif()

# pcs_noir_export_bench (from recursion/pcs_noir_export.bench.cpp) writes the inputs of the Noir
# recursive verifiers in recursion/noir/; it drives the transparent Honk prover directly.
if(TARGET pcs_noir_export_bench)
    target_link_libraries(pcs_noir_export_bench PRIVATE ultra_honk)
endif()

# vela_evm_export_bench proves a Noir circuit with UltraHonk + Vela under a Keccak transcript and
# writes the opening argument for the Solidity gas harness; it loads ACIR, so it needs the frontend.
if(TARGET vela_evm_export_bench)
    target_link_libraries(vela_evm_export_bench PRIVATE dsl ultra_honk vm2_stub)
endif()

# whir_differential_bench compares barretenberg's WHIR against the reference Rust implementation
# ProveKit builds on, through the C ABI in whir/rust/. That library is not part of the C++ build:
# build it with `cargo build --release` in whir/rust/, then reconfigure. Without it the target is
# skipped rather than failing the build, since cargo is not a barretenberg build dependency.
set(WHIR_RS_FFI_LIB
    ${CMAKE_CURRENT_SOURCE_DIR}/whir/rust/target/release/${CMAKE_STATIC_LIBRARY_PREFIX}whir_rs_ffi${CMAKE_STATIC_LIBRARY_SUFFIX})
if(TARGET whir_differential_bench)
    if(EXISTS ${WHIR_RS_FFI_LIB})
        target_link_libraries(whir_differential_bench PRIVATE ${WHIR_RS_FFI_LIB})
        if(APPLE)
            # Rust's std pulls these in on macOS.
            target_link_libraries(whir_differential_bench PRIVATE "-framework CoreFoundation" "-framework Security")
        endif()
    else()
        message(STATUS "whir_rs_ffi not built; excluding whir_differential_bench")
        set_target_properties(whir_differential_bench PROPERTIES EXCLUDE_FROM_ALL TRUE)
    endif()
endif()
# whir_recursion_bench measures the in-circuit WHIR verifier; its `outer` mode additionally proves
# the aggregation circuit with UltraHonk + KZG, which needs the prover.
if(TARGET whir_recursion_bench)
    target_link_libraries(whir_recursion_bench PRIVATE ultra_honk)
endif()
