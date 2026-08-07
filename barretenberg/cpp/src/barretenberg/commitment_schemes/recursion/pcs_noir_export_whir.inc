// WHIR-Honk export: included by pcs_noir_export.bench.cpp.

namespace {

using WhirHonkPcs = whir::WhirPcs<Hasher>;
using WhirHonkExport = honk_transparent::TransparentHonk<WhirHonkPcs, Flavor>;

whir::WhirSoundness parse_soundness(const std::string& name)
{
    if (name == "conjectured") {
        return whir::WhirSoundness::CONJECTURED_LIST;
    }
    if (name == "repaired") {
        return whir::WhirSoundness::REPAIRED_LIST;
    }
    if (name == "unique") {
        return whir::WhirSoundness::UNIQUE_DECODING;
    }
    return whir::WhirSoundness::PROVABLE_LIST;
}

/**
 * @brief Segment a WHIR opening argument and re-verify it from the segmented data alone.
 *
 * @details Mirrors `WhirVerifier::verify` in the interleaved, non-zk layout that `WhirHonk` uses:
 * every claimed column is its own constituent of the batched oracle, all claimed at the single
 * sumcheck point, so there are no column-stacking challenges and no cross evaluations. Each
 * iteration runs `k` sumcheck rounds against the running claim, commits the folded oracle, samples
 * it out of domain, and queries it in domain; the queries of iteration 0 hit the four Honk
 * commitment trees and every later query hits the single folded oracle.
 */
int export_whir(UltraCircuitBuilder& builder,
                size_t log_n,
                size_t security_bits,
                size_t log_inv_rate,
                size_t folding_factor_bits,
                size_t initial_folding_factor_bits,
                const std::string& soundness_name,
                const std::string& out)
{
    const whir::WhirConfig config = whir::WhirConfig::create(log_n,
                                                             security_bits,
                                                             log_inv_rate,
                                                             folding_factor_bits,
                                                             /*final_poly_bits=*/4,
                                                             parse_soundness(soundness_name),
                                                             /*zk=*/false,
                                                             /*max_stack_bits=*/0,
                                                             initial_folding_factor_bits,
                                                             /*pow_bits=*/0);
    WhirHonkExport::ProvingKey pk = WhirHonkExport::create_proving_key(builder, config);
    const auto vk = pk.vk;
    const HonkProof proof = WhirHonkExport::prove(pk);
    if (!WhirHonkExport::verify(vk, config, proof)) {
        std::cerr << "native WHIR-Honk verification failed\n";
        return 1;
    }

    const ClaimLayout layout = build_claim_layout(vk.virtual_mask, vk.num_committed_precomputed());
    MirrorTranscript transcript(proof);
    const ShellData shell = read_shell(transcript,
                                       vk.log_dyadic_size,
                                       vk.num_public_inputs,
                                       vk.pub_inputs_offset,
                                       vk.virtual_mask,
                                       vk.lagrange_first_row,
                                       vk.lagrange_last_row,
                                       vk.precomputed_commitment);

    const std::vector<FF> roots = {
        vk.precomputed_commitment, shell.witness_roots[0], shell.witness_roots[1], shell.witness_roots[2]
    };
    const size_t k0 = config.initial_folding_factor_bits;
    const size_t round0_arity = size_t(1) << k0;
    const size_t round0_index_bits = config.num_variables + config.log_inv_rate - k0;

    // Constituents: the unshifted claims then the shifted ones, in claim order.
    struct Constituent {
        size_t tree;
        size_t leaf_column;
        bool shifted;
        FF value;
    };
    std::vector<Constituent> constituents;
    for (size_t i = 0; i < layout.unshifted.size(); ++i) {
        constituents.push_back({ layout.unshifted[i].group,
                                 layout.unshifted[i].column,
                                 false,
                                 shell.evaluations[layout.unshifted_entity[i]] });
    }
    for (size_t l = 0; l < layout.shifted.size(); ++l) {
        constituents.push_back({ layout.shifted[l].group,
                                 layout.shifted[l].column,
                                 true,
                                 shell.evaluations[layout.shifted_entity[l]] });
    }
    const size_t num_constituents = constituents.size();

    // Out-of-domain samples against the initial oracle.
    std::vector<FF> z_ood_initial(config.num_ood_samples);
    std::vector<std::vector<FF>> ood_values(config.num_ood_samples, std::vector<FF>(num_constituents));
    for (size_t s = 0; s < config.num_ood_samples; ++s) {
        z_ood_initial[s] = transcript.challenge();
        for (size_t c = 0; c < num_constituents; ++c) {
            ood_values[s][c] = transcript.receive_one();
        }
    }

    const FF rho = transcript.challenge();
    const FF gamma_claims = transcript.challenge();

    std::vector<FF> rho_powers(num_constituents);
    FF rho_power = FF::one();
    for (size_t c = 0; c < num_constituents; ++c) {
        rho_powers[c] = rho_power;
        rho_power *= rho;
    }
    FF sigma = FF::zero();
    std::vector<whir::WeightTerm> terms;
    FF gamma_power = FF::one();
    {
        FF oracle_value = FF::zero();
        for (size_t c = 0; c < num_constituents; ++c) {
            oracle_value += rho_powers[c] * constituents[c].value;
        }
        sigma += gamma_power * oracle_value;
        terms.push_back(whir::WeightTerm::eq_weight(shell.challenge, gamma_power));
        gamma_power *= gamma_claims;
    }
    for (size_t s = 0; s < config.num_ood_samples; ++s) {
        FF oracle_value = FF::zero();
        for (size_t c = 0; c < num_constituents; ++c) {
            oracle_value += rho_powers[c] * ood_values[s][c];
        }
        sigma += gamma_power * oracle_value;
        terms.push_back(whir::WeightTerm::pow_weight(z_ood_initial[s], config.num_variables, gamma_power));
        gamma_power *= gamma_claims;
    }

    // Flat, per-round exports. Round r of `rounds` is followed by the final round.
    std::vector<FF> round_sumcheck;  // k triples per round
    std::vector<FF> round_roots_g;   // one per iteration
    std::vector<FF> round_y_ood;     // one per iteration
    std::vector<FF> round_nonces;    // one per iteration and one for the final round
    std::vector<FF> opened_values;   // padded per batch
    std::vector<FF> opened_siblings; // padded per batch
    std::vector<size_t> hint_leaves;
    std::vector<size_t> hint_num_leaves;
    std::vector<size_t> hint_positions;
    std::vector<size_t> hint_num_siblings;
    std::vector<FF> sigma_trace;

    std::optional<FF> folded_root;
    const auto record_batch = [&](const BatchOpening& opening, size_t leaf_width, size_t depth, size_t max_leaves) {
        opened_values.insert(opened_values.end(), opening.values.begin(), opening.values.end());
        opened_siblings.insert(opened_siblings.end(), opening.siblings.begin(), opening.siblings.end());
        std::vector<size_t> padded_leaves = opening.plan.leaves;
        padded_leaves.resize(max_leaves, 0);
        hint_leaves.insert(hint_leaves.end(), padded_leaves.begin(), padded_leaves.end());
        hint_num_leaves.push_back(opening.plan.leaves.size());
        hint_positions.insert(
            hint_positions.end(), opening.plan.positions.begin(), opening.plan.positions.end());
        hint_num_siblings.push_back(opening.plan.num_siblings);
        static_cast<void>(leaf_width);
        static_cast<void>(depth);
    };

    const size_t num_iterations = config.rounds.size();
    for (size_t i = 0; i < num_iterations + 1; ++i) {
        const bool is_final = i == num_iterations;
        const whir::WhirRound& round = is_final ? config.final_round : config.rounds[i];
        const size_t k = round.folding_factor_bits;

        if (!is_final) {
            std::vector<FF> alphas(k);
            for (size_t j = 0; j < k; ++j) {
                const std::vector<FF> h = transcript.receive(3);
                round_sumcheck.insert(round_sumcheck.end(), h.begin(), h.end());
                if (h[0] + h[1] != sigma) {
                    std::cerr << "mirror: WHIR sumcheck check failed at round " << i << " step " << j << "\n";
                    return 1;
                }
                alphas[j] = transcript.challenge();
                sigma = whir::detail::evaluate_quadratic({ h[0], h[1], h[2] }, alphas[j]);
                for (whir::WeightTerm& term : terms) {
                    term.bind_variable(alphas[j]);
                }
            }
            sigma_trace.push_back(sigma);
            round_roots_g.push_back(transcript.receive_one());
            const FF z_ood = transcript.challenge();
            round_y_ood.push_back(transcript.receive_one());

            const size_t index_bits = round.log_domain_size - k;
            transcript.challenge(); // pow seed; pow_bits = 0, so no grinding check
            round_nonces.push_back(transcript.receive_one());
            std::vector<size_t> indices(round.num_queries);
            for (size_t s = 0; s < round.num_queries; ++s) {
                indices[s] = whir::detail::index_from_challenge(transcript.challenge(), index_bits);
            }
            const FF omega = fr::get_root_of_unity(round.log_domain_size);
            const FF eta_inv = omega.pow(uint256_t(uint64_t(1)) << index_bits).invert();

            std::vector<std::vector<FF>> cosets(round.num_queries);
            if (i == 0) {
                std::vector<BatchOpening> batches;
                for (size_t g = 0; g < roots.size(); ++g) {
                    const size_t leaf_width = layout.group_columns[g] * round0_arity;
                    BatchOpening opening =
                        read_batch(transcript, indices, leaf_width, round0_index_bits, round.num_queries);
                    if (verify_batch_root(opening, leaf_width, round0_index_bits) != roots[g]) {
                        std::cerr << "mirror: round-0 Merkle batch failed for group " << g << "\n";
                        return 1;
                    }
                    record_batch(opening, leaf_width, round0_index_bits, round.num_queries);
                    batches.push_back(std::move(opening));
                }
                for (size_t s = 0; s < round.num_queries; ++s) {
                    std::vector<FF> x_inverses(round0_arity);
                    FF x_inv = omega.pow(uint256_t(uint64_t(indices[s]))).invert();
                    for (size_t t = 0; t < round0_arity; ++t) {
                        x_inverses[t] = x_inv;
                        x_inv *= eta_inv;
                    }
                    std::vector<FF> coset(round0_arity, FF::zero());
                    for (size_t c = 0; c < num_constituents; ++c) {
                        const Constituent& constituent = constituents[c];
                        const BatchOpening& batch = batches[constituent.tree];
                        const size_t leaf_width = layout.group_columns[constituent.tree] * round0_arity;
                        const size_t slot = batch.plan.positions[s];
                        for (size_t t = 0; t < round0_arity; ++t) {
                            FF value =
                                rho_powers[c] *
                                batch.values[slot * leaf_width + constituent.leaf_column * round0_arity + t];
                            if (constituent.shifted) {
                                value *= x_inverses[t];
                            }
                            coset[t] += value;
                        }
                    }
                    cosets[s] = std::move(coset);
                }
            } else {
                BatchOpening opening =
                    read_batch(transcript, indices, size_t(1) << k, index_bits, round.num_queries);
                if (verify_batch_root(opening, size_t(1) << k, index_bits) != *folded_root) {
                    std::cerr << "mirror: folded Merkle batch failed at round " << i << "\n";
                    return 1;
                }
                record_batch(opening, size_t(1) << k, index_bits, round.num_queries);
                for (size_t s = 0; s < round.num_queries; ++s) {
                    const size_t slot = opening.plan.positions[s];
                    cosets[s].assign(opening.values.begin() + static_cast<std::ptrdiff_t>(slot << k),
                                     opening.values.begin() + static_cast<std::ptrdiff_t>((slot + 1) << k));
                }
            }

            const FF gamma = transcript.challenge();
            const size_t next_vars = round.num_variables - k;
            FF gamma_power_round = gamma;
            sigma += gamma_power_round * round_y_ood.back();
            terms.push_back(whir::WeightTerm::pow_weight(z_ood, next_vars, gamma_power_round));
            for (size_t s = 0; s < round.num_queries; ++s) {
                const FF x_base_inv = omega.pow(uint256_t(uint64_t(indices[s]))).invert();
                const FF folded = whir::fold_coset(cosets[s], x_base_inv, eta_inv, alphas);
                gamma_power_round *= gamma;
                sigma += gamma_power_round * folded;
                terms.push_back(whir::WeightTerm::pow_weight(
                    omega.pow(uint256_t(uint64_t(indices[s])) << k), next_vars, gamma_power_round));
            }
            folded_root = round_roots_g.back();
            continue;
        }

        // Final phase: the clear polynomial, its consistency with the last oracle, and the claim.
        const size_t final_size = size_t(1) << config.final_round.num_variables;
        const std::vector<FF> final_poly = transcript.receive(final_size);
        const size_t index_bits = config.final_round.log_domain_size - k;
        const FF omega = fr::get_root_of_unity(config.final_round.log_domain_size);
        const FF eta = omega.pow(uint256_t(uint64_t(1)) << index_bits);
        const FF eta_inv = eta.invert();
        transcript.challenge(); // final pow seed
        round_nonces.push_back(transcript.receive_one());
        std::vector<size_t> final_indices(config.final_round.num_queries);
        for (size_t s = 0; s < config.final_round.num_queries; ++s) {
            final_indices[s] = whir::detail::index_from_challenge(transcript.challenge(), index_bits);
        }
        std::vector<std::vector<FF>> cosets(config.final_round.num_queries);
        if (num_iterations == 0) {
            std::cerr << "mirror: schedules with no folding iteration are not exported\n";
            return 1;
        }
        BatchOpening opening =
            read_batch(transcript, final_indices, size_t(1) << k, index_bits, config.final_round.num_queries);
        if (verify_batch_root(opening, size_t(1) << k, index_bits) != *folded_root) {
            std::cerr << "mirror: final Merkle batch failed\n";
            return 1;
        }
        record_batch(opening, size_t(1) << k, index_bits, config.final_round.num_queries);
        for (size_t s = 0; s < config.final_round.num_queries; ++s) {
            const size_t slot = opening.plan.positions[s];
            cosets[s].assign(opening.values.begin() + static_cast<std::ptrdiff_t>(slot << k),
                             opening.values.begin() + static_cast<std::ptrdiff_t>((slot + 1) << k));
            FF point = omega.pow(uint256_t(uint64_t(final_indices[s])));
            for (const FF& value : cosets[s]) {
                if (value != polynomial_arithmetic::evaluate(final_poly.data(), point, final_size)) {
                    std::cerr << "mirror: final oracle/polynomial mismatch at query " << s << "\n";
                    return 1;
                }
                point *= eta;
            }
        }
        FF total = FF::zero();
        for (const whir::WeightTerm& term : terms) {
            total += term.weighted_sum(final_poly);
        }
        if (total != sigma) {
            std::cerr << "mirror: final WHIR claim failed\n";
            return 1;
        }
        static_cast<void>(eta_inv);

        if (!transcript.exhausted()) {
            std::cerr << "mirror: " << proof.size() - transcript.cursor() << " unread proof fields\n";
            return 1;
        }

        // --- emit ---------------------------------------------------------------------------
        NoirWriter params(out + "/params.nr");
        params.comment("WHIR-Honk recursive verifier parameters");
        params.comment("inner circuit 2^" + std::to_string(log_n) + " rows, lambda = " +
                       std::to_string(security_bits) + ", rate 2^-" + std::to_string(log_inv_rate) +
                       ", k = " + std::to_string(folding_factor_bits) + ", k0 = " + std::to_string(k0) +
                       ", soundness = " + soundness_name);
        params.blank();
        emit_shell_params(params, shell, layout);

        params.comment("WHIR schedule");
        params.global_u32("WH_NUM_VARIABLES", config.num_variables);
        params.global_u32("WH_LOG_INV_RATE", config.log_inv_rate);
        params.global_u32("WH_K0", k0);
        params.global_u32("WH_NUM_OOD", config.num_ood_samples);
        params.global_u32("WH_NUM_CONSTITUENTS", num_constituents);
        params.global_u32("WH_NUM_ITERATIONS", num_iterations);
        params.global_u32("WH_ROUND0_ARITY", round0_arity);
        params.global_u32("WH_ROUND0_INDEX_BITS", round0_index_bits);
        std::vector<size_t> constituent_tree;
        std::vector<size_t> constituent_column;
        std::vector<bool> constituent_shifted;
        std::vector<size_t> constituent_entity;
        for (size_t c = 0; c < num_constituents; ++c) {
            constituent_tree.push_back(constituents[c].tree);
            constituent_column.push_back(constituents[c].leaf_column);
            constituent_shifted.push_back(constituents[c].shifted);
            constituent_entity.push_back(c < layout.unshifted.size()
                                              ? layout.unshifted_entity[c]
                                              : layout.shifted_entity[c - layout.unshifted.size()]);
        }
        params.global_u32_array("WH_CONSTITUENT_TREE", constituent_tree);
        params.global_u32_array("WH_CONSTITUENT_COLUMN", constituent_column);
        params.global_bool_array("WH_CONSTITUENT_SHIFTED", constituent_shifted);
        params.global_u32_array("WH_CONSTITUENT_ENTITY", constituent_entity);

        std::vector<size_t> round_k;
        std::vector<size_t> round_queries;
        std::vector<size_t> round_log_domain;
        std::vector<size_t> round_index_bits;
        std::vector<size_t> round_num_variables;
        std::vector<FF> round_omega;
        std::vector<FF> round_eta_inv;
        std::vector<size_t> sumcheck_offsets;
        size_t sumcheck_offset = 0;
        for (size_t r = 0; r < num_iterations + 1; ++r) {
            const whir::WhirRound& entry = r == num_iterations ? config.final_round : config.rounds[r];
            const size_t bits = entry.log_domain_size - entry.folding_factor_bits;
            const FF root = fr::get_root_of_unity(entry.log_domain_size);
            round_k.push_back(entry.folding_factor_bits);
            round_queries.push_back(entry.num_queries);
            round_log_domain.push_back(entry.log_domain_size);
            round_index_bits.push_back(bits);
            round_num_variables.push_back(entry.num_variables);
            round_omega.push_back(root);
            round_eta_inv.push_back(root.pow(uint256_t(uint64_t(1)) << bits).invert());
            sumcheck_offsets.push_back(sumcheck_offset);
            sumcheck_offset += r == num_iterations ? 0 : entry.folding_factor_bits;
        }
        params.global_u32_array("WH_ROUND_K", round_k);
        params.global_u32_array("WH_ROUND_QUERIES", round_queries);
        params.global_u32_array("WH_ROUND_LOG_DOMAIN", round_log_domain);
        params.global_u32_array("WH_ROUND_INDEX_BITS", round_index_bits);
        params.global_u32_array("WH_ROUND_NUM_VARIABLES", round_num_variables);
        params.global_field_array("WH_ROUND_OMEGA", round_omega);
        params.global_field_array("WH_ROUND_ETA_INV", round_eta_inv);
        params.global_u32_array("WH_SUMCHECK_OFFSET", sumcheck_offsets);
        params.global_u32("WH_TOTAL_SUMCHECK_ROUNDS", sumcheck_offset);
        params.global_u32("WH_FINAL_VARS", config.final_round.num_variables);
        params.global_u32("WH_FINAL_SIZE", final_size);
        params.global_u32("WH_MAX_QUERIES", *std::max_element(round_queries.begin(), round_queries.end()));
        params.global_u32("WH_MAX_K", *std::max_element(round_k.begin(), round_k.end()));
        params.global_u32("WH_MAX_ARITY", size_t(1) << *std::max_element(round_k.begin(), round_k.end()));
        params.global_u32("WH_MAX_INDEX_BITS",
                          *std::max_element(round_index_bits.begin(), round_index_bits.end()));
        params.global_u32("WH_NUM_TERMS", terms.size());

        // Batch layout: round 0 opens the four Honk trees, every later round the folded oracle.
        std::vector<size_t> batch_leaf_width;
        std::vector<size_t> batch_depth;
        std::vector<size_t> batch_max_leaves;
        std::vector<size_t> batch_value_offset;
        std::vector<size_t> batch_sibling_offset;
        std::vector<size_t> batch_leaf_offset;
        std::vector<size_t> batch_position_offset;
        size_t value_offset = 0;
        size_t sibling_offset = 0;
        size_t leaf_offset = 0;
        size_t position_offset = 0;
        const auto add_batch = [&](size_t leaf_width, size_t depth, size_t max_leaves) {
            batch_leaf_width.push_back(leaf_width);
            batch_depth.push_back(depth);
            batch_max_leaves.push_back(max_leaves);
            batch_value_offset.push_back(value_offset);
            batch_sibling_offset.push_back(sibling_offset);
            batch_leaf_offset.push_back(leaf_offset);
            batch_position_offset.push_back(position_offset);
            value_offset += max_leaves * leaf_width;
            sibling_offset += max_leaves * depth;
            leaf_offset += max_leaves;
            position_offset += max_leaves;
        };
        for (size_t g = 0; g < roots.size(); ++g) {
            add_batch(layout.group_columns[g] * round0_arity, round0_index_bits, config.rounds[0].num_queries);
        }
        for (size_t r = 1; r < num_iterations + 1; ++r) {
            const whir::WhirRound& entry = r == num_iterations ? config.final_round : config.rounds[r];
            add_batch(size_t(1) << entry.folding_factor_bits,
                      entry.log_domain_size - entry.folding_factor_bits,
                      entry.num_queries);
        }
        params.global_u32("WH_NUM_BATCHES", batch_leaf_width.size());
        params.global_u32_array("WH_BATCH_LEAF_WIDTH", batch_leaf_width);
        params.global_u32_array("WH_BATCH_DEPTH", batch_depth);
        params.global_u32_array("WH_BATCH_MAX_LEAVES", batch_max_leaves);
        params.global_u32_array("WH_BATCH_VALUE_OFFSET", batch_value_offset);
        params.global_u32_array("WH_BATCH_SIBLING_OFFSET", batch_sibling_offset);
        params.global_u32_array("WH_BATCH_LEAF_OFFSET", batch_leaf_offset);
        params.global_u32_array("WH_BATCH_POSITION_OFFSET", batch_position_offset);
        params.global_u32("WH_MAX_LEAF_WIDTH",
                          *std::max_element(batch_leaf_width.begin(), batch_leaf_width.end()));
        params.global_u32("WH_MAX_DEPTH", *std::max_element(batch_depth.begin(), batch_depth.end()));
        params.global_u32("WH_MAX_BATCH_LEAVES",
                          *std::max_element(batch_max_leaves.begin(), batch_max_leaves.end()));
        params.global_u32("WH_TOTAL_VALUES", value_offset);
        params.global_u32("WH_TOTAL_SIBLINGS", sibling_offset);
        params.global_u32("WH_TOTAL_LEAVES", leaf_offset);
        params.global_u32("WH_TOTAL_POSITIONS", position_offset);
        params.write();

        TomlWriter toml(out + "/Prover.toml");
        emit_shell_toml(toml, shell);
        std::vector<FF> ood_flat;
        for (const auto& row : ood_values) {
            ood_flat.insert(ood_flat.end(), row.begin(), row.end());
        }
        toml.field_array_2d("ood_values", ood_flat, num_constituents);
        toml.field_array_2d("round_sumcheck", round_sumcheck, 3);
        toml.field_array("round_roots", round_roots_g);
        toml.field_array("round_y_ood", round_y_ood);
        toml.field_array("round_nonces", round_nonces);
        toml.field_array("final_poly", final_poly);
        toml.field_array("opened_values", opened_values);
        toml.field_array("opened_siblings", opened_siblings);
        toml.u32_array("hint_leaves", hint_leaves);
        toml.u32_array("hint_num_leaves", hint_num_leaves);
        toml.u32_array("hint_positions", hint_positions);
        toml.write();

        std::cerr << "whir: iterations " << num_iterations << ", constituents " << num_constituents
                  << ", round-0 queries " << config.rounds[0].num_queries << ", arity " << round0_arity
                  << ", depth " << round0_index_bits << ", proof " << proof.size() << " fields ("
                  << (proof.size() * 32) / 1024 << " KiB)\n";
        return 0;
    }
    return 1;
}

} // namespace
