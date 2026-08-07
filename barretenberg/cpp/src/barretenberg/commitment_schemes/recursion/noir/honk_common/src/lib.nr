//! In-circuit verification of a transparent UltraHonk proof: the pieces shared by the WHIR and
//! Ligero recursive verifiers.
//!
//! Everything here is native BN254 field arithmetic plus Poseidon2 - no elliptic-curve or
//! non-native arithmetic appears anywhere in a hash-based verifier.

pub mod merkle;
pub mod poseidon;
pub mod relations;
pub mod shell;
pub mod transcript;
