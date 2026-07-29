#pragma once

#include "barretenberg/commitment_schemes/commitment_key.hpp"
#include "barretenberg/ecc/curves/bn254/g1.hpp"
#include "barretenberg/ecc/curves/bn254/g2.hpp"

#include <vector>

namespace bb::cq {

/**
 * @brief TEST-ONLY structured reference string with G2 powers of tau.
 *
 * @details cq requires pairing checks against [T(tau)]_2 and [tau^N]_2, i.e. a setup with N+1 powers of tau in G2.
 * Powers-of-tau ceremonies (e.g. perpetual powers of tau) publish such G2 powers, but the SRS files shipped with
 * barretenberg only contain [tau]_2, so for tests we synthesize a fresh SRS from a locally sampled tau. The toxic
 * waste is known to the process; NEVER use this outside tests.
 *
 * The G1 side is deliberately truncated to exactly the requested number of points: cq's soundness relies on the
 * prover being unable to commit to polynomials of degree >= N (this is what bounds deg(A) < N and, via the shifted
 * commitment [B_0 * X^{N-n+1}], deg(B_0) <= n-2).
 */
struct TestSrs {
    std::vector<g1::affine_element> g1_powers; // [tau^0]_1, ..., [tau^{num_g1-1}]_1
    std::vector<g2::affine_element> g2_powers; // [tau^0]_2, ..., [tau^{num_g2-1}]_2

    static TestSrs create(size_t num_g1_points, size_t num_g2_points);

    /**
     * @brief A CommitmentKey over this SRS's G1 powers. The global CRS factory is init-once, so tests that need
     * several independent SRSs (e.g. different table sizes, each with its own tau) construct keys directly.
     */
    CommitmentKey<curve::BN254> create_commitment_key() const;
};

} // namespace bb::cq
