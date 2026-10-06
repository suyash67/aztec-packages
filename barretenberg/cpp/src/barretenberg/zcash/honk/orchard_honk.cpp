#include "orchard_honk.hpp"
#include "barretenberg/commitment_schemes/claim_batcher.hpp"
#include "barretenberg/commitment_schemes/kzg/kzg.hpp"
#include "barretenberg/commitment_schemes/shplonk/shplemini.hpp"
#include "barretenberg/commitment_schemes/small_subgroup_ipa/small_subgroup_ipa_impl.hpp"
#include "barretenberg/common/bb_bench.hpp"
#include "barretenberg/srs/global_crs.hpp"
#include "barretenberg/sumcheck/sumcheck.hpp"
#include "barretenberg/zcash/honk/halo2_ipa.hpp"

#include <numeric>

namespace bb {
template class SmallSubgroupIPAProver<zcash::OrchardFlavor>;
template class SmallSubgroupIPAProver<zcash::OrchardBn254Flavor>;
} // namespace bb

namespace bb::zcash {

namespace {

constexpr size_t CHUNK = 7;

#define ORCHARD_TYPES(Cycle)                                                                                           \
    using Flavor [[maybe_unused]] = OrchardFlavor_<Cycle>;                                                             \
    using FF [[maybe_unused]] = typename Flavor::FF;                                                                   \
    using Curve [[maybe_unused]] = typename Flavor::Curve;                                                             \
    using Commitment [[maybe_unused]] = typename Flavor::Commitment;                                                   \
    using Polynomial [[maybe_unused]] = typename Flavor::Polynomial;                                                   \
    using Trace [[maybe_unused]] = halo2::AnchoredTrace<Cycle>;                                                        \
    using Gates [[maybe_unused]] = halo2::OrchardGates<Cycle>;                                                         \
    using ProverPolynomials [[maybe_unused]] = typename Flavor::ProverPolynomials;                                     \
    using CommitmentKey [[maybe_unused]] = typename Flavor::CommitmentKey;                                             \
    using Transcript [[maybe_unused]] = typename Flavor::Transcript;                                                   \
    [[maybe_unused]] constexpr size_t MASK_BEGIN = Flavor::NUM_ZERO_ROWS;                                              \
    [[maybe_unused]] constexpr size_t MASK_END = Flavor::TRACE_OFFSET;                                                 \
    static_assert(MASK_BEGIN < MASK_END)

// Random values in the masked rows NUM_ZERO_ROWS..TRACE_OFFSET-1 (shared by both flavors).
template <typename FF> void mask(Polynomial<FF>& poly)
{
    for (size_t row = OrchardFlavor::NUM_ZERO_ROWS; row < OrchardFlavor::TRACE_OFFSET; ++row) {
        poly.at(row) = FF::random_element();
    }
}

// halo2 label of a permutation cell (column j, row r).
template <typename FF> FF permutation_label(size_t column, size_t row, size_t n)
{
    return FF((column * n) + row);
}

template <typename Flavor>
void hash_preamble(typename Flavor::Transcript& transcript,
                   const typename Flavor::VerificationKey& vk,
                   const std::vector<typename Flavor::FF>& public_inputs)
{
    transcript.add_to_hash_buffer("vk_hash", vk.hash());
    for (size_t i = 0; i < public_inputs.size(); ++i) {
        transcript.add_to_hash_buffer("public_input_" + std::to_string(i), public_inputs[i]);
    }
}

template <typename Curve>
typename Curve::AffineElement batch_mul(std::span<const typename Curve::AffineElement> commitments,
                                        std::span<const typename Curve::ScalarField> scalars)
{
    std::vector<typename Curve::ScalarField> s;
    std::vector<typename Curve::AffineElement> p;
    for (size_t i = 0; i < commitments.size(); ++i) {
        if (!commitments[i].is_point_at_infinity()) {
            s.push_back(scalars[i]);
            p.push_back(commitments[i]);
        }
    }
    return typename Curve::AffineElement(Halo2IPA<Curve>::msm(s, p));
}

template <typename Cycle> void init_crs()
{
    if constexpr (std::is_same_v<Cycle, Bn254Cycle>) {
        srs::init_bn254_file_crs_factory(srs::bb_crs_path());
    }
}

} // namespace

template <typename Cycle>
OrchardProvingKey_<Cycle>::OrchardProvingKey_(const halo2::AnchoredTrace<Cycle>& t)
    : circuit_size(t.num_rows)
    , log_circuit_size(numeric::get_msb(t.num_rows))
    , precomputed(t.num_rows)
    , vk(std::make_shared<typename OrchardFlavor_<Cycle>::VerificationKey>())
    , public_input_cells(t.public_input_cells)
{
    BB_BENCH_NAME("OrchardProvingKey");
    ORCHARD_TYPES(Cycle);
    init_crs<Cycle>();
    BB_ASSERT_EQ(t.row_offset, Flavor::TRACE_OFFSET);
    const size_t n = circuit_size;
    auto& p = precomputed;
    auto fixed = static_cast<typename Flavor::template FixedEntities<Polynomial>&>(p).get_all();
    for (size_t c = 0; c < halo2::NUM_FIXED; ++c) {
        for (size_t r = 0; r < n; ++r) {
            fixed[c].at(r) = t.fixed[c][r];
        }
    }
    auto selectors = static_cast<typename Flavor::template SelectorEntities<Polynomial>&>(p).get_all();
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
    auto sigmas = static_cast<typename Flavor::template SigmaEntities<Polynomial>&>(p).get_all();
    auto ids = static_cast<typename Flavor::template IdEntities<Polynomial>&>(p).get_all();
    for (size_t j = 0; j < Flavor::NUM_PERMUTATION_COLUMNS; ++j) {
        for (size_t r = 0; r < n; ++r) {
            const uint32_t cell = static_cast<uint32_t>((j * n) + r);
            ids[j].at(r) = permutation_label<FF>(j, r, n);
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

    CommitmentKey ck(n);
    for (auto [commitment, poly] : zip_view(vk->get_all(), p.get_precomputed())) {
        commitment = ck.commit(poly);
    }
}

template <typename Cycle>
typename OrchardFlavor_<Cycle>::Proof orchard_prove(const OrchardProvingKey_<Cycle>& pk,
                                                    const halo2::AnchoredTrace<Cycle>& trace)
{
    BB_BENCH_NAME("orchard_prove");
    ORCHARD_TYPES(Cycle);
    init_crs<Cycle>();
    const size_t n = pk.circuit_size;
    const size_t log_n = pk.log_circuit_size;
    BB_ASSERT_EQ(trace.num_rows, n);
    auto transcript = std::make_shared<Transcript>();
    CommitmentKey ck(n);
    typename Flavor::CommitmentLabels labels;

    hash_preamble<Flavor>(*transcript, *pk.vk, trace.public_inputs);

    ProverPolynomials polys;
    {
        auto& pre = const_cast<ProverPolynomials&>(pk.precomputed);
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
    auto advice = static_cast<typename Flavor::template AdviceEntities<Polynomial>&>(polys).get_all();
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
                    const auto [m, x, y] = Gates::template sinsemilla_lookup_value<FF>(view, chip);
                    const uint256_t idx(m);
                    BB_ASSERT_LT(idx, uint256_t(halo2::SINSEMILLA_TABLE_SIZE));
                    polys.lookup_read_counts_sinsemilla.at(trace.row_offset + idx.data[0]) += FF(1);
                }
            }
            if (!trace.selectors[halo2::Q_LOOKUP][r].is_zero()) {
                const uint256_t v(Gates::template range_lookup_value<FF>(view));
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
             zip_view(advice, static_cast<typename Flavor::template AdviceEntities<std::string>&>(labels).get_all())) {
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
                            const auto [m, x, y] = Gates::template sinsemilla_lookup_value<FF>(view, chip);
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
                        read += Gates::template range_lookup_value<FF>(view);
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
        auto sigmas = static_cast<typename Flavor::template SigmaEntities<Polynomial>&>(polys).get_all();
        auto ids = static_cast<typename Flavor::template IdEntities<Polynomial>&>(polys).get_all();
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
    using PolynomialBatcher = typename GeminiProver_<Curve>::PolynomialBatcher;
    PolynomialBatcher batcher(n, n);
    batcher.set_unshifted(polys.get_unshifted());
    batcher.set_to_be_shifted_by_one(polys.get_to_be_shifted());
    batcher.set_to_be_shifted_by_two(polys.get_to_be_shifted_by_two());
    auto opening_claim = ShpleminiProver_<Curve>::prove(
        n, batcher, sumcheck_output.challenge, ck, transcript, small_subgroup_ipa.get_witness_polynomials());
    if constexpr (Flavor::IS_PASTA) {
        BB_BENCH_NAME("orchard_prove/halo2_ipa");
        Halo2IPA<Curve>::prove(halo2_vesta_ipa_generators(ck, n), opening_claim, FF(0), transcript);
    } else {
        KZG<Curve>::compute_opening_proof(ck, opening_claim, transcript);
    }
    return transcript->export_proof();
}

template <typename Cycle>
bool orchard_verify(const typename OrchardFlavor_<Cycle>::VerificationKey& vk,
                    const std::vector<typename Cycle::FF>& public_inputs,
                    const typename OrchardFlavor_<Cycle>::Proof& proof)
{
    BB_BENCH_NAME("orchard_verify");
    ORCHARD_TYPES(Cycle);
    init_crs<Cycle>();
    if (public_inputs.size() != vk.num_public_inputs) {
        return false;
    }
    const size_t n = vk.circuit_size();
    const size_t log_n = vk.log_circuit_size;
    auto transcript = std::make_shared<Transcript>(proof);
    typename Flavor::CommitmentLabels labels;
    hash_preamble<Flavor>(*transcript, vk, public_inputs);

    typename Flavor::VerifierCommitments comms(vk);
    comms.gemini_masking_poly = transcript->template receive_from_prover<Commitment>("Gemini:masking_poly_comm");
    for (auto [comm, label] :
         zip_view(static_cast<typename Flavor::template AdviceEntities<Commitment>&>(comms).get_all(),
                  static_cast<typename Flavor::template AdviceEntities<std::string>&>(labels).get_all())) {
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
    using Batch = typename ClaimBatcher::Batch;
    auto& evals = sumcheck_output.claimed_evaluations;
    ClaimBatcher claim_batcher{
        .unshifted = Batch{ comms.get_unshifted(), evals.get_unshifted() },
        .shifted = Batch{ comms.get_to_be_shifted(), evals.get_shifted() },
        .shifted_by_two = Batch{ comms.get_to_be_shifted_by_two(), evals.get_shifted_by_two() },
    };
    CommitmentKey ck(n);
    const Commitment g1_identity = ck.get_monomial_points()[0];
    auto shplemini_output =
        ShpleminiVerifier_<Curve, true, true>::compute_batch_opening_claim(claim_batcher,
                                                                           sumcheck_output.challenge,
                                                                           g1_identity,
                                                                           transcript,
                                                                           {},
                                                                           libra_commitments,
                                                                           sumcheck_output.claimed_libra_evaluation);
    bool pcs_verified = false;
    if constexpr (Flavor::IS_PASTA) {
        const auto& batch_claim = shplemini_output.batch_opening_claim;
        const OpeningClaim<Curve> opening_claim{ { batch_claim.evaluation_point, FF(0) },
                                                 batch_mul<Curve>(batch_claim.commitments, batch_claim.scalars) };
        pcs_verified = Halo2IPA<Curve>::verify(halo2_vesta_ipa_generators(ck, n), opening_claim, transcript);
    } else {
        auto pairing_points =
            KZG<Curve>::reduce_verify_batch_opening_claim(std::move(shplemini_output.batch_opening_claim), transcript);
        pcs_verified = pairing_points.check();
    }
    vinfo("orchard verifier: sumcheck ",
          sumcheck_output.verified,
          ", libra consistency ",
          shplemini_output.consistency_checked,
          ", pcs ",
          pcs_verified);
    return sumcheck_output.verified && shplemini_output.consistency_checked && pcs_verified;
}

template struct OrchardProvingKey_<PastaCycle>;
template struct OrchardProvingKey_<Bn254Cycle>;
template OrchardFlavor::Proof orchard_prove<PastaCycle>(const OrchardProvingKey_<PastaCycle>&,
                                                        const halo2::AnchoredTrace<PastaCycle>&);
template OrchardBn254Flavor::Proof orchard_prove<Bn254Cycle>(const OrchardProvingKey_<Bn254Cycle>&,
                                                             const halo2::AnchoredTrace<Bn254Cycle>&);
template bool orchard_verify<PastaCycle>(const OrchardFlavor::VerificationKey&,
                                         const std::vector<PastaCycle::FF>&,
                                         const OrchardFlavor::Proof&);
template bool orchard_verify<Bn254Cycle>(const OrchardBn254Flavor::VerificationKey&,
                                         const std::vector<Bn254Cycle::FF>&,
                                         const OrchardBn254Flavor::Proof&);

} // namespace bb::zcash
