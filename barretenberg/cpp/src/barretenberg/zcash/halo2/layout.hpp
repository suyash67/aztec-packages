#pragma once

#include <cstddef>

/**
 * @file layout.hpp
 * @brief Column layout of the halo2 Orchard Action circuit (orchard 0.16.0 `Config::configure`), shared by the chip
 * ports, the Honk relations and the flavor.
 *
 * @details
 *   advice  a0..a9   all equality-enabled
 *   fixed   f0..f7   Lagrange coefficients of fixed-base mul; f2..f4 = Poseidon rc_a, f5..f7 = rc_b;
 *                    f0 = constants column and Sinsemilla #1 y_Q; f1 = Sinsemilla #2 y_Q
 *           fixed_z  fixed-base mul z-values (private to the ECC chip)
 *           q_s2_1,  Sinsemilla q_sinsemilla2 fixed columns (one per Sinsemilla chip)
 *           q_s2_2
 *   lookup  table_idx, table_x, table_y  (Sinsemilla S table; table_idx alone is the 10-bit range table)
 *   equality-enabled columns: a0..a9, f0 (constants), f5..f7 (Poseidon rc_b)
 *
 * Selectors are listed in the order their gates are configured. halo2 configures the NoteCommit gates once per
 * NoteCommit chip (old/new note) with identical columns; the two copies define identical polynomials, so a single set
 * of NoteCommit selectors is used here. Public inputs (halo2's instance column) are held in otherwise-unused rows of
 * advice column a0 and bound by the permutation argument's public-input delta, as in Honk.
 */
namespace bb::zcash::halo2::layout {

inline constexpr size_t NUM_ADVICE = 10;

enum Fixed : size_t {
    F0 = 0,
    F1,
    F2,
    F3,
    F4,
    F5,
    F6,
    F7,
    FIXED_Z,
    Q_SINSEMILLA2_1,
    Q_SINSEMILLA2_2,
    NUM_FIXED,
};

enum Selector : size_t {
    Q_ORCHARD = 0,
    Q_ADD_FIELD,
    // ECC
    Q_POINT,
    Q_POINT_NON_ID,
    Q_ADD_INCOMPLETE,
    Q_ADD,
    Q_MUL_HI_1,
    Q_MUL_HI_2,
    Q_MUL_HI_3,
    Q_MUL_LO_1,
    Q_MUL_LO_2,
    Q_MUL_LO_3,
    Q_MUL_DECOMPOSE_VAR,
    Q_MUL_LSB,
    Q_MUL_OVERFLOW,
    Q_MUL_FIXED_RUNNING_SUM,
    Q_MUL_FIXED_FULL,
    Q_MUL_FIXED_SHORT,
    Q_MUL_FIXED_BASE_FIELD,
    // Poseidon
    Q_POSEIDON_FULL,
    Q_POSEIDON_PARTIAL,
    Q_POSEIDON_PAD_AND_ADD,
    // Sinsemilla / Merkle (one set per chip instance)
    Q_SINSEMILLA1_1,
    Q_SINSEMILLA4_1,
    Q_SINSEMILLA1_2,
    Q_SINSEMILLA4_2,
    Q_MERKLE_DECOMPOSE_1,
    Q_MERKLE_DECOMPOSE_2,
    Q_SWAP_1,
    Q_SWAP_2,
    // Lookup range check
    Q_LOOKUP,
    Q_RUNNING,
    Q_BITSHIFT,
    // CommitIvk
    Q_COMMIT_IVK,
    // NoteCommit
    Q_NOTECOMMIT_B,
    Q_NOTECOMMIT_D,
    Q_NOTECOMMIT_E,
    Q_NOTECOMMIT_G,
    Q_NOTECOMMIT_H,
    Q_NOTECOMMIT_G_D,
    Q_NOTECOMMIT_PK_D,
    Q_NOTECOMMIT_VALUE,
    Q_NOTECOMMIT_RHO,
    Q_NOTECOMMIT_PSI,
    Q_Y_CANON,
    NUM_SELECTORS,
};

// Number of equality-enabled columns in the permutation argument: a0..a9, f0, f5, f6, f7.
inline constexpr size_t NUM_PERMUTATION_COLUMNS = 14;

// Number of public inputs per Action (instance rows 0..9).
inline constexpr size_t NUM_PUBLIC_INPUTS_PER_ACTION = 10;

enum PublicInput : size_t {
    ANCHOR = 0,
    CV_NET_X,
    CV_NET_Y,
    NF_OLD,
    RK_X,
    RK_Y,
    CMX,
    ENABLE_SPEND,
    ENABLE_OUTPUT,
    DISABLE_CROSS_ADDRESS,
};

} // namespace bb::zcash::halo2::layout
