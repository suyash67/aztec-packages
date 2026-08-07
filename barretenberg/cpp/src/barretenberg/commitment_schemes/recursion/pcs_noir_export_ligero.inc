// Ligero-Honk export: included by pcs_noir_export.bench.cpp.

namespace {

using LigeroHonk = honk_transparent::TransparentHonk<ligero::LigeroPcs<Hasher>, Flavor>;

/**
 * @brief Segment a Ligero opening argument and re-verify it from the segmented data alone.
 *
 * @details Mirrors `LigeroVerifier::verify`. The verifier's per-query work is: draw a column index,
 * read that column of every committed group, authenticate it against the group's Merkle root, and
 * check that the three claimed row combinations agree with the encodings of the combined rows at
 * that column. The circuit never needs the whole encoding - only its value at the queried index -
 * so the export records the combined rows and the circuit evaluates
 * `enc(w)[idx] = sum_c w[c] * (omega^idx)^c` directly.
 */
int export_ligero(UltraCircuitBuilder& builder,
                  size_t log_n,
                  size_t security_bits,
                  size_t log_inv_rate,
                  const std::string& out)
{
    const ligero::LigeroConfig config = LigeroHonk::make_config(log_n, security_bits, log_inv_rate);
    LigeroHonk::ProvingKey pk = LigeroHonk::create_proving_key(builder, config);
    const auto vk = pk.vk;
    const HonkProof proof = LigeroHonk::prove(pk);
    if (!LigeroHonk::verify(vk, config, proof)) {
        std::cerr << "native Ligero-Honk verification failed\n";
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

    const size_t num_rows = config.num_rows();
    const size_t num_cols = config.num_cols();
    const size_t codeword_length = config.codeword_length();
    const size_t index_bits = numeric::get_msb(codeword_length);
    const std::vector<FF> roots = {
        vk.precomputed_commitment, shell.witness_roots[0], shell.witness_roots[1], shell.witness_roots[2]
    };

    const FF rho = transcript.challenge();
    const std::vector<FF> combined = transcript.receive_unhashed(3 * num_cols);
    const FF w_digest = transcript.receive_one();
    if (w_digest != Hasher::hash_leaf(combined, std::nullopt)) {
        std::cerr << "mirror: combined-row digest mismatch\n";
        return 1;
    }

    // Per-(group, column) rho scalars, in the same order the claims are batched.
    std::map<std::pair<size_t, size_t>, std::pair<FF, FF>> scalars;
    FF expected = FF::zero();
    FF rho_power = FF::one();
    for (size_t i = 0; i < layout.unshifted.size(); ++i) {
        scalars[{ layout.unshifted[i].group, layout.unshifted[i].column }].first += rho_power;
        expected += rho_power * shell.evaluations[layout.unshifted_entity[i]];
        rho_power *= rho;
    }
    for (size_t l = 0; l < layout.shifted.size(); ++l) {
        scalars[{ layout.shifted[l].group, layout.shifted[l].column }].second += rho_power;
        expected += rho_power * shell.evaluations[layout.shifted_entity[l]];
        rho_power *= rho;
    }

    const std::vector<FF> a = whir::eq_tensor(std::span<const FF>(shell.challenge).subspan(0, config.log_num_cols));
    const std::vector<FF> b = whir::eq_tensor(std::span<const FF>(shell.challenge).subspan(config.log_num_cols));
    const std::span<const FF> w_u(combined.data(), num_cols);
    const std::span<const FF> w_s(combined.data() + num_cols, num_cols);
    const std::span<const FF> w_s2(combined.data() + 2 * num_cols, num_cols);

    FF total = FF::zero();
    for (size_t c = 0; c < num_cols; ++c) {
        total += w_u[c] * a[c];
        if (c >= 1) {
            total += w_s[c] * a[c - 1];
        }
    }
    total += a[num_cols - 1] * w_s2[0];
    if (total != expected) {
        std::cerr << "mirror: Ligero claim equation failed\n";
        return 1;
    }

    // The circuit evaluates the row encodings pointwise; check that against the real encoder once.
    const FF omega = fr::get_root_of_unity(index_bits);
    ligero::RSCodePolicy code(config);
    const std::vector<FF> enc_u = code.encode(w_u);
    const auto encode_at = [&](std::span<const FF> row, size_t idx) {
        const FF point = omega.pow(uint256_t(uint64_t(idx)));
        FF accumulator = FF::zero();
        for (size_t c = num_cols; c-- > 0;) {
            accumulator = accumulator * point + row[c];
        }
        return accumulator;
    };
    for (size_t idx = 0; idx < codeword_length; idx += std::max<size_t>(codeword_length / 8, 1)) {
        if (encode_at(w_u, idx) != enc_u[idx]) {
            std::cerr << "mirror: pointwise RS evaluation disagrees with the encoder\n";
            return 1;
        }
    }
    const std::vector<FF> enc_s = code.encode(w_s);
    const std::vector<FF> enc_s2 = code.encode(w_s2);

    // Column queries.
    std::vector<size_t> query_indices;
    std::vector<FF> opened_values; // [query][group][column * num_rows + row]
    std::vector<FF> opened_paths;  // [query][group][level]
    for (size_t s = 0; s < config.num_queries; ++s) {
        const FF challenge = transcript.challenge();
        const size_t idx = static_cast<size_t>(uint256_t(challenge).data[0] % codeword_length);
        query_indices.push_back(idx);
        FF acc_u = FF::zero();
        FF acc_s = FF::zero();
        FF acc_s2 = FF::zero();
        for (size_t g = 0; g < layout.group_columns.size(); ++g) {
            const size_t group_rows = layout.group_columns[g] * num_rows;
            const std::vector<FF> flat = transcript.receive_unhashed(group_rows + index_bits);
            const std::vector<FF> values(flat.begin(), flat.begin() + static_cast<std::ptrdiff_t>(group_rows));
            const std::vector<FF> path(flat.begin() + static_cast<std::ptrdiff_t>(group_rows), flat.end());
            opened_values.insert(opened_values.end(), values.begin(), values.end());
            opened_paths.insert(opened_paths.end(), path.begin(), path.end());

            FF digest = Hasher::hash_leaf(values, std::nullopt);
            size_t index = idx;
            for (const FF& sibling : path) {
                digest = (index & 1) != 0 ? Hasher::hash_node(sibling, digest) : Hasher::hash_node(digest, sibling);
                index >>= 1;
            }
            if (digest != roots[g]) {
                std::cerr << "mirror: Merkle authentication failed at query " << s << " group " << g << "\n";
                return 1;
            }
            for (size_t p = 0; p < layout.group_columns[g]; ++p) {
                const auto it = scalars.find({ g, p });
                if (it == scalars.end()) {
                    continue;
                }
                const auto& [scalar_u, scalar_s] = it->second;
                for (size_t r = 0; r < num_rows; ++r) {
                    const FF& value = values[p * num_rows + r];
                    acc_u += scalar_u * b[r] * value;
                    acc_s += scalar_s * b[r] * value;
                    if (r > 0) {
                        acc_s2 += scalar_s * b[r - 1] * value;
                    }
                }
            }
        }
        if (acc_u != enc_u[idx] || acc_s != enc_s[idx] || acc_s2 != enc_s2[idx]) {
            std::cerr << "mirror: Ligero column check failed at query " << s << "\n";
            return 1;
        }
    }
    if (!transcript.exhausted()) {
        std::cerr << "mirror: " << proof.size() - transcript.cursor() << " unread proof fields\n";
        return 1;
    }

    // --- emit ---------------------------------------------------------------------------------
    NoirWriter params(out + "/params.nr");
    params.comment("Ligero-Honk recursive verifier parameters");
    params.comment("inner circuit 2^" + std::to_string(log_n) + " rows, lambda = " +
                   std::to_string(security_bits) + ", rate 2^-" + std::to_string(log_inv_rate));
    params.blank();
    emit_shell_params(params, shell, layout);
    params.comment("Ligero");
    params.global_u32("LG_LOG_NUM_COLS", config.log_num_cols);
    params.global_u32("LG_NUM_COLS", num_cols);
    params.global_u32("LG_NUM_ROWS", num_rows);
    params.global_u32("LG_CODEWORD_LENGTH", codeword_length);
    params.global_u32("LG_INDEX_BITS", index_bits);
    params.global_u32("LG_NUM_QUERIES", config.num_queries);
    params.global_field("LG_OMEGA", omega);
    std::vector<size_t> group_rows;
    std::vector<size_t> group_row_offsets;
    size_t offset = 0;
    for (const size_t columns : layout.group_columns) {
        group_rows.push_back(columns * num_rows);
        group_row_offsets.push_back(offset);
        offset += columns * num_rows;
    }
    params.global_u32_array("LG_GROUP_ROWS", group_rows);
    params.global_u32_array("LG_GROUP_ROW_OFFSET", group_row_offsets);
    params.global_u32("LG_TOTAL_GROUP_ROWS", offset);
    params.write();

    TomlWriter toml(out + "/Prover.toml");
    emit_shell_toml(toml, shell);
    toml.field_array("combined_rows", combined);
    toml.field("combined_digest", w_digest);
    toml.field_array("opened_values", opened_values);
    toml.field_array("opened_paths", opened_paths);
    toml.write();

    std::cerr << "ligero: queries " << config.num_queries << ", rows " << num_rows << ", cols " << num_cols
              << ", index bits " << index_bits << ", proof " << proof.size() << " fields ("
              << (proof.size() * 32) / 1024 << " KiB)\n";
    return 0;
}

} // namespace
