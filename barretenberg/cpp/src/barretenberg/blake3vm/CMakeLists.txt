barretenberg_module(blake3vm sumcheck commitment_schemes crypto_blake3s honk)

# whir_settlement.test.cpp drives the sequencer pipeline: a WHIR-Honk proof verified inside an
# UltraHonk circuit, and that circuit proved. The library itself must not depend on the prover.
if(TARGET blake3vm_tests)
    target_link_libraries(blake3vm_tests PRIVATE ultra_honk stdlib_blake3s)
endif()
