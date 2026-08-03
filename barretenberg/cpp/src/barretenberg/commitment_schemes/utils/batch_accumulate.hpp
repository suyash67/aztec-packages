#pragma once

#include "barretenberg/common/thread.hpp"
#include "barretenberg/ecc/curves/bn254/fr.hpp"

#include <span>
#include <vector>

namespace bb::pcs_utils {

/**
 * @brief One term of a random linear combination: `scalar * source`, optionally left-shifted by one.
 * @details `shifted` selects the Honk shift contract's virtual polynomial `A(X)/X`, whose dense array
 * is `source[i + 1]`; the term stops one index short of the output length because `source[n]` does not
 * exist (the shift contract's zero constant coefficient makes the omitted slot irrelevant).
 */
struct ScaledTerm {
    const fr* source;
    fr scalar;
    bool shifted = false;
};

/**
 * @brief `out[i] += Σ_k terms[k].scalar · terms[k].source[i + shifted]`, in one parallel pass over `i`.
 *
 * @details The batched claim sets every backend opens are random linear combinations of 40-ish dense
 * arrays. Accumulating them one term at a time is a serial `num_terms · N` pass that dominates the
 * opening phase; splitting the *output* range across threads and looping terms inside each chunk keeps
 * writes disjoint (so no locking) and touches each output cache line once per chunk rather than once
 * per term. This is the `std::vector`-backed counterpart of `bb::add_scaled_batch`, which requires
 * `Polynomial` operands.
 */
inline void accumulate_scaled(std::span<fr> out, std::span<const ScaledTerm> terms)
{
    if (terms.empty() || out.empty()) {
        return;
    }
    const size_t n = out.size();
    parallel_for_range(n, [&](size_t start, size_t end) {
        for (const ScaledTerm& term : terms) {
            const size_t offset = term.shifted ? 1 : 0;
            const size_t stop = term.shifted ? std::min(end, n - 1) : end;
            for (size_t i = start; i < stop; ++i) {
                out[i] += term.scalar * term.source[i + offset];
            }
        }
    });
}

} // namespace bb::pcs_utils
