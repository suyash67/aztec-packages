#include "orchard_honk.hpp"
#include "barretenberg/commitment_schemes/claim_batcher.hpp"
#include "barretenberg/commitment_schemes/shplonk/shplemini.hpp"
#include "barretenberg/commitment_schemes/small_subgroup_ipa/small_subgroup_ipa_impl.hpp"
#include "barretenberg/common/bb_bench.hpp"
#include "barretenberg/sumcheck/sumcheck.hpp"
#include "barretenberg/zcash/honk/halo2_ipa.hpp"

#include <numeric>

namespace bb {
template class SmallSubgroupIPAProver<zcash::OrchardFlavor>;
} // namespace bb

namespace bb::zcash {

namespace {

using Flavor = OrchardFlavor;
using FF = Flavor::FF;
using Curve = Flavor::Curve;
using Commitment = Flavor::Commitment;
using Polynomial = Flavor::Polynomial;
using Trace = halo2::AnchoredTrace<PastaCycle>;
using Gates = halo2::OrchardGates<PastaCycle>;
using IPA = Halo2IPA<Curve>;

constexpr size_t CHUNK = 7;

// First and last rows of the masked region of every witness polynomial.
constexpr size_t MASK_BEGIN = Flavor::NUM_ZERO_ROWS;
constexpr size_t MASK_END = Flavor::TRACE_OFFSET;

void mask(Polynomial& poly)
{
    for (size_t row = MASK_BEGIN; row < MASK_END; ++row) {
        poly.at(row) = FF::random_element();
    }
}

IPA::Generators ipa_generators(const Flavor::CommitmentKey& ck, size_t n)
{
    static const Commitment w = halo2_vesta_w();
    static const Commitment u = halo2_vesta_u();
    return { std::span<const Commitment>(ck.get_monomial_points().data(), n), w, u };
}

// halo2 label of a permutation cell (column j, row r).
FF permutation_label(size_t column, size_t row, size_t n)
{
    return FF((column * n) + row);
}

void hash_preamble(Flavor::Transcript& transcript,
                   const Flavor::VerificationKey& vk,
                   const std::vector<FF>& public_inputs)
{
    transcript.add_to_hash_buffer("vk_hash", vk.hash());
    for (size_t i = 0; i < public_inputs.size(); ++i) {
        transcript.add_to_hash_buffer("public_input_" + std::to_string(i), public_inputs[i]);
    }
}

Commitment batch_mul(std::span<const Commitment> commitments, std::span<const FF> scalars)
{
    std::vector<FF> s;
    std::vector<Commitment> p;
    for (size_t i = 0; i < commitments.size(); ++i) {
        if (!commitments[i].is_point_at_infinity()) {
            s.push_back(scalars[i]);
            p.push_back(commitments[i]);
        }
    }
    return Commitment(IPA::msm(s, p));
}

} // namespace

OrchardProvingKey::OrchardProvingKey(const Trace& t)
    : circuit_size(t.num_rows)
    , log_circuit_size(numeric::get_msb(t.num_rows))
    , precomputed(t.num_rows)
    , vk(std::make_shared<Flavor::VerificationKey>())
    , public_input_cells(t.public_input_cells)
{
    BB_BENCH_NAME("OrchardProvingKey");
    BB_ASSERT_EQ(t.row_offset, Flavor::TRACE_OFFSET);
    const size_t n = circuit_size;
    auto& p = precomputed;
    auto fixed = static_cast<Flavor::FixedEntities<Polynomial>&>(p).get_all();
    for (size_t c = 0; c < halo2::NUM_FIXED; ++c) {
        for (size_t r = 0; r < n; ++r) {
            fixed[c].at(r) = t.fixed[c][r];
        }
    }
    auto selectors = static_cast<Flavor::SelectorEntities<Polynomial>&>(p).get_all();
    for (size_t s = 0; s < halo2::NUM_SELECTORS; ++s) {
        for (size_t r = 0; r < n; ++r) {
            selectors[s].at(r) = t.selectors[s][r];
        }
    }
    for (size_t r = 0; r < n; ++r) {
        p.table_idx.at(r) = t.table[0][r];
        p.table_x.at(r) = t.table[1][r];
        p.table_y.at(r) = t.table[2][r];
        p.q_table.at(r) = t.q_table[r];
    }
    p.lagrange_first.at(Flavor::TRACE_OFFSET) = FF(1);
    p.lagrange_last.at(n - 1) = FF(1);

    // Copy cycles over the cells (column j, row r) -> j * n + r.
    const size_t num_cells = Flavor::NUM_PERMUTATION_COLUMNS * n;
    std::vector<uint32_t> parent(num_cells);
    std::iota(parent.begin(), parent.end(), 0);
    auto find = [&](uint32_t x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    };
    for (const auto& [a, b] : t.copies) {
        const auto ra = find(static_cast<uint32_t>((a.first * n) + a.second));
        const auto rb = find(static_cast<uint32_t>((b.first * n) + b.second));
        if (ra != rb) {
            parent[ra] = rb;
        }
    }
    // next[c] = following cell in c's cycle (cycles in increasing cell order).
    std::vector<uint32_t> next(num_cells);
    std::iota(next.begin(), next.end(), 0);
    {
        std::vector<uint32_t> first_of(num_cells, UINT32_MAX);
        std::vector<uint32_t> last_of(num_cells, UINT32_MAX);
        for (uint32_t c = 0; c < num_cells; ++c) {
            const auto root = find(c);
            if (first_of[root] == UINT32_MAX) {
                first_of[root] = c;
            } else {
                next[last_of[root]] = c;
            }
            last_of[root] = c;
        }
        for (uint32_t root = 0; root < num_cells; ++root) {
            if (last_of[root] != UINT32_MAX) {
                next[last_of[root]] = first_of[root];
            }
        }
    }
    auto sigmas = static_cast<Flavor::SigmaEntities<Polynomial>&>(p).get_all();
    auto ids = static_cast<Flavor::IdEntities<Polynomial>&>(p).get_all();
    for (size_t j = 0; j < Flavor::NUM_PERMUTATION_COLUMNS; ++j) {
        for (size_t r = 0; r < n; ++r) {
            const uint32_t cell = static_cast<uint32_t>((j * n) + r);
            ids[j].at(r) = permutation_label(j, r, n);
            sigmas[j].at(r) = FF(next[cell]);
        }
    }
    vk->log_circuit_size = log_circuit_size;
    vk->num_public_inputs = public_input_cells.size();
    for (size_t i = 0; i < public_input_cells.size(); ++i) {
        const auto [col, row] = public_input_cells[i];
        vk->public_input_next_labels.push_back(sigmas[col][row]);
        sigmas[col].at(row) = vk->special_label(i);
    }

    Flavor::CommitmentKey ck(n);
    for (auto [commitment, poly] : zip_view(vk->get_all(), p.get_precomputed())) {
        commitment = ck.commit(poly);
    }
}

Flavor::Proof orchard_prove(const OrchardProvingKey& pk, const Trace& trace)
{
    BB_BENCH_NAME("orchard_prove");
    const size_t n = pk.circuit_size;
    const size_t log_n = pk.log_circuit_size;
    BB_ASSERT_EQ(trace.num_rows, n);
    auto transcript = std::make_shared<Flavor::Transcript>();
    Flavor::CommitmentKey ck(n);
    Flavor::CommitmentLabels labels;

    hash_preamble(*transcript, *pk.vk, trace.public_inputs);

    Flavor::ProverPolynomials polys;
    {
        auto& pre = const_cast<Flavor::ProverPolynomials&>(pk.precomputed);
        for (auto [dst, src] : zip_view(polys.get_precomputed(), pre.get_precomputed())) {
            dst = src.share();
        }
    }
    for (auto& poly : polys.get_witness()) {
        poly = Polynomial(n - Flavor::NUM_ZERO_ROWS, n, Flavor::NUM_ZERO_ROWS);
    }

    // Gemini masking polynomial
    polys.gemini_masking_poly = Polynomial::random(n);
    transcript->send_to_verifier("Gemini:masking_poly_comm", ck.commit(polys.gemini_masking_poly));

    // Advice columns and lookup read counts
    auto advice = static_cast<Flavor::AdviceEntities<Polynomial>&>(polys).get_all();
    parallel_for(halo2::NUM_ADVICE, [&](size_t c) {
        for (size_t r = MASK_END; r < n; ++r) {
            advice[c].at(r) = trace.advice[c][r];
        }
    });
    for (auto& a : advice) {
        mask(a);
    }
    {
        for (size_t r = MASK_END; r + 2 < n; ++r) {
            typename Trace::RowView view{ trace, r };
            for (size_t chip = 0; chip < 2; ++chip) {
                const size_t q = (chip == 0) ? halo2::Q_SINSEMILLA1_1 : halo2::Q_SINSEMILLA1_2;
                if (!trace.selectors[q][r].is_zero()) {
                    const auto [m, x, y] = Gates::sinsemilla_lookup_value<FF>(view, chip);
                    const uint256_t idx(m);
                    BB_ASSERT_LT(idx, uint256_t(halo2::SINSEMILLA_TABLE_SIZE));
                    polys.lookup_read_counts_sinsemilla.at(trace.row_offset + idx.data[0]) += FF(1);
                }
            }
            if (!trace.selectors[halo2::Q_LOOKUP][r].is_zero()) {
                const uint256_t v(Gates::range_lookup_value<FF>(view));
                BB_ASSERT_LT(v, uint256_t(halo2::SINSEMILLA_TABLE_SIZE));
                polys.lookup_read_counts_range.at(trace.row_offset + v.data[0]) += FF(1);
            }
        }
        mask(polys.lookup_read_counts_sinsemilla);
        mask(polys.lookup_read_counts_range);
    }
    {
        auto batch = ck.start_batch();
        for (auto [poly, label] :
             zip_view(advice, static_cast<Flavor::AdviceEntities<std::string>&>(labels).get_all())) {
            batch.add_to_batch(poly, label);
        }
        batch.add_to_batch(polys.lookup_read_counts_sinsemilla, labels.lookup_read_counts_sinsemilla);
        batch.add_to_batch(polys.lookup_read_counts_range, labels.lookup_read_counts_range);
        batch.commit_and_send_to_verifier(transcript);
    }

    RelationParameters<FF> params;
    params.eta = transcript->template get_challenge<FF>("eta");
    params.eta_two = params.eta.sqr();
    auto [beta, gamma] = transcript->template get_challenges<FF>(std::array<std::string, 2>{ "beta", "gamma" });
    params.beta = beta;
    params.gamma = gamma;
    params.public_input_delta = pk.vk->compute_public_input_delta(trace.public_inputs, beta, gamma);

    // Lookup inverses
    {
        BB_BENCH_NAME("orchard_prove/lookup_inverses");
        const size_t count = n - MASK_END;
        std::vector<FF> sinsemilla(count, FF(0));
        std::vector<FF> range(count, FF(0));
        parallel_for_range(count, [&](size_t start, size_t end) {
            for (size_t i = start; i < end; ++i) {
                const size_t r = MASK_END + i;
                typename Trace::RowView view{ trace, r };
                // Rows past the end read zeros for the next two rows.
                const bool in_bounds = r + 2 < n;
                const bool q1 = in_bounds && !trace.selectors[halo2::Q_SINSEMILLA1_1][r].is_zero();
                const bool q2 = in_bounds && !trace.selectors[halo2::Q_SINSEMILLA1_2][r].is_zero();
                const bool qt = !trace.q_table[r].is_zero();
                const bool ql = in_bounds && !trace.selectors[halo2::Q_LOOKUP][r].is_zero();
                const FF t = trace.table[0][r] + trace.table[1][r] * params.eta + trace.table[2][r] * params.eta_two +
                             params.gamma;
                if (q1 || q2 || qt) {
                    FF product = t;
                    for (size_t chip = 0; chip < 2; ++chip) {
                        FF read = params.gamma;
                        if (in_bounds) {
                            const auto [m, x, y] = Gates::sinsemilla_lookup_value<FF>(view, chip);
                            read += m + x * params.eta + y * params.eta_two;
                        } else {
                            BB_ASSERT(false, "lookup table row in the last two rows");
                        }
                        product *= read;
                    }
                    sinsemilla[i] = product;
                }
                if (ql || qt) {
                    FF read = params.gamma;
                    if (in_bounds) {
                        read += Gates::range_lookup_value<FF>(view);
                    }
                    range[i] = read * (trace.table[0][r] + params.gamma);
                }
            }
        });
        FF::batch_invert(std::span<FF>(sinsemilla));
        FF::batch_invert(std::span<FF>(range));
        for (size_t i = 0; i < count; ++i) {
            polys.lookup_inverses_sinsemilla.at(MASK_END + i) = sinsemilla[i];
            polys.lookup_inverses_range.at(MASK_END + i) = range[i];
        }
        mask(polys.lookup_inverses_sinsemilla);
        mask(polys.lookup_inverses_range);
    }

    // Permutation grand product (two chunks of 7 columns)
    {
        BB_BENCH_NAME("orchard_prove/grand_product");
        const size_t first = Flavor::TRACE_OFFSET;
        const size_t count = n - first;
        auto sigmas = static_cast<Flavor::SigmaEntities<Polynomial>&>(polys).get_all();
        auto ids = static_cast<Flavor::IdEntities<Polynomial>&>(polys).get_all();
        std::vector<FF> num_a(count);
        std::vector<FF> den_a(count);
        std::vector<FF> num_b(count);
        std::vector<FF> den_b(count);
        parallel_for_range(count, [&](size_t start, size_t end) {
            for (size_t i = start; i < end; ++i) {
                const size_t r = first + i;
                FF na(1);
                FF da(1);
                FF nb(1);
                FF db(1);
                for (size_t j = 0; j < Flavor::NUM_PERMUTATION_COLUMNS; ++j) {
                    const FF& w = trace.permutation_value(j, r);
                    const FF num = w + beta * ids[j][r] + gamma;
                    const FF den = w + beta * sigmas[j][r] + gamma;
                    if (j < CHUNK) {
                        na *= num;
                        da *= den;
                    } else {
                        nb *= num;
                        db *= den;
                    }
                }
                num_a[i] = na;
                den_a[i] = da;
                num_b[i] = nb;
                den_b[i] = db;
            }
        });
        FF::batch_invert(std::span<FF>(den_a));
        FF::batch_invert(std::span<FF>(den_b));
        FF z(1);
        for (size_t i = 0; i < count; ++i) {
            const size_t r = first + i;
            if (i > 0) {
                polys.z_perm.at(r) = z;
            }
            const FF mid = z * num_a[i] * den_a[i];
            polys.z_perm_mid.at(r) = mid;
            z = mid * num_b[i] * den_b[i];
        }
        BB_ASSERT(z == params.public_input_delta, "permutation grand product does not match the public input delta");
        mask(polys.z_perm);
        mask(polys.z_perm_mid);
    }
    {
        auto batch = ck.start_batch();
        batch.add_to_batch(polys.lookup_inverses_sinsemilla, labels.lookup_inverses_sinsemilla);
        batch.add_to_batch(polys.lookup_inverses_range, labels.lookup_inverses_range);
        batch.add_to_batch(polys.z_perm_mid, labels.z_perm_mid);
        batch.add_to_batch(polys.z_perm, labels.z_perm);
        batch.commit_and_send_to_verifier(transcript);
    }
    polys.set_shifted();

    // Sumcheck
    const FF alpha = transcript->template get_challenge<FF>("Sumcheck:alpha");
    std::vector<FF> gate_challenges =
        transcript->template get_dyadic_powers_of_challenge<FF>("Sumcheck:gate_challenge", log_n);
    SumcheckProver<Flavor> sumcheck(n, polys, transcript, alpha, gate_challenges, params, log_n);
    auto zk_sumcheck_data = [&]() {
        BB_BENCH_NAME("orchard_prove/zk_sumcheck_data");
        return ZKSumcheckData<Flavor>(log_n, transcript, ck);
    }();
    auto sumcheck_output = [&]() {
        BB_BENCH_NAME("orchard_prove/sumcheck");
        return sumcheck.prove(zk_sumcheck_data);
    }();

    // Shplemini + halo2 IPA
    SmallSubgroupIPAProver<Flavor> small_subgroup_ipa(
        zk_sumcheck_data, sumcheck_output.challenge, sumcheck_output.claimed_libra_evaluation, transcript, ck);
    {
        BB_BENCH_NAME("orchard_prove/small_subgroup_ipa");
        small_subgroup_ipa.prove();
    }
    using PolynomialBatcher = GeminiProver_<Curve>::PolynomialBatcher;
    PolynomialBatcher batcher(n, n);
    batcher.set_unshifted(polys.get_unshifted());
    batcher.set_to_be_shifted_by_one(polys.get_to_be_shifted());
    batcher.set_to_be_shifted_by_two(polys.get_to_be_shifted_by_two());
    auto opening_claim = ShpleminiProver_<Curve>::prove(
        n, batcher, sumcheck_output.challenge, ck, transcript, small_subgroup_ipa.get_witness_polynomials());
    {
        BB_BENCH_NAME("orchard_prove/halo2_ipa");
        IPA::prove(ipa_generators(ck, n), opening_claim, FF(0), transcript);
    }
    return transcript->export_proof();
}

bool orchard_verify(const Flavor::VerificationKey& vk, const std::vector<FF>& public_inputs, const Flavor::Proof& proof)
{
    BB_BENCH_NAME("orchard_verify");
    if (public_inputs.size() != vk.num_public_inputs) {
        return false;
    }
    const size_t n = vk.circuit_size();
    const size_t log_n = vk.log_circuit_size;
    auto transcript = std::make_shared<Flavor::Transcript>(proof);
    Flavor::CommitmentLabels labels;
    hash_preamble(*transcript, vk, public_inputs);

    Flavor::VerifierCommitments comms(vk);
    comms.gemini_masking_poly = transcript->template receive_from_prover<Commitment>("Gemini:masking_poly_comm");
    for (auto [comm, label] : zip_view(static_cast<Flavor::AdviceEntities<Commitment>&>(comms).get_all(),
                                       static_cast<Flavor::AdviceEntities<std::string>&>(labels).get_all())) {
        comm = transcript->template receive_from_prover<Commitment>(label);
    }
    comms.lookup_read_counts_sinsemilla =
        transcript->template receive_from_prover<Commitment>(labels.lookup_read_counts_sinsemilla);
    comms.lookup_read_counts_range =
        transcript->template receive_from_prover<Commitment>(labels.lookup_read_counts_range);

    RelationParameters<FF> params;
    params.eta = transcript->template get_challenge<FF>("eta");
    params.eta_two = params.eta.sqr();
    auto [beta, gamma] = transcript->template get_challenges<FF>(std::array<std::string, 2>{ "beta", "gamma" });
    params.beta = beta;
    params.gamma = gamma;
    params.public_input_delta = vk.compute_public_input_delta(public_inputs, beta, gamma);

    comms.lookup_inverses_sinsemilla =
        transcript->template receive_from_prover<Commitment>(labels.lookup_inverses_sinsemilla);
    comms.lookup_inverses_range = transcript->template receive_from_prover<Commitment>(labels.lookup_inverses_range);
    comms.z_perm_mid = transcript->template receive_from_prover<Commitment>(labels.z_perm_mid);
    comms.z_perm = transcript->template receive_from_prover<Commitment>(labels.z_perm);

    const FF alpha = transcript->template get_challenge<FF>("Sumcheck:alpha");
    std::vector<FF> gate_challenges =
        transcript->template get_dyadic_powers_of_challenge<FF>("Sumcheck:gate_challenge", log_n);
    SumcheckVerifier<Flavor> sumcheck(transcript, alpha, log_n);
    std::array<Commitment, NUM_SMALL_IPA_COMMITMENTS> libra_commitments = {};
    libra_commitments[0] = transcript->template receive_from_prover<Commitment>("Libra:concatenation_commitment");
    auto sumcheck_output = sumcheck.verify(params, gate_challenges);
    libra_commitments[1] = transcript->template receive_from_prover<Commitment>("Libra:grand_sum_commitment");
    libra_commitments[2] = transcript->template receive_from_prover<Commitment>("Libra:quotient_commitment");

    using ClaimBatcher = ClaimBatcher_<Curve>;
    using Batch = ClaimBatcher::Batch;
    auto& evals = sumcheck_output.claimed_evaluations;
    ClaimBatcher claim_batcher{
        .unshifted = Batch{ comms.get_unshifted(), evals.get_unshifted() },
        .shifted = Batch{ comms.get_to_be_shifted(), evals.get_shifted() },
        .shifted_by_two = Batch{ comms.get_to_be_shifted_by_two(), evals.get_shifted_by_two() },
    };
    Flavor::CommitmentKey ck(n);
    const auto gens = ipa_generators(ck, n);
    auto shplemini_output =
        ShpleminiVerifier_<Curve, true, true>::compute_batch_opening_claim(claim_batcher,
                                                                           sumcheck_output.challenge,
                                                                           gens.g[0],
                                                                           transcript,
                                                                           {},
                                                                           libra_commitments,
                                                                           sumcheck_output.claimed_libra_evaluation);
    const auto& batch_claim = shplemini_output.batch_opening_claim;
    const OpeningClaim<Curve> opening_claim{ { batch_claim.evaluation_point, FF(0) },
                                             batch_mul(batch_claim.commitments, batch_claim.scalars) };
    const bool ipa_verified = IPA::verify(gens, opening_claim, transcript);
    vinfo("orchard verifier: sumcheck ",
          sumcheck_output.verified,
          ", libra consistency ",
          shplemini_output.consistency_checked,
          ", ipa ",
          ipa_verified);
    return sumcheck_output.verified && shplemini_output.consistency_checked && ipa_verified;
}

} // namespace bb::zcash
