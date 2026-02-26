#ifndef HOMOLOGIZER_RELABEL_HPP
#define HOMOLOGIZER_RELABEL_HPP

#include "utils.hpp"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace homologizer {

struct RelabeledPhaseTable {
    std::vector<std::int64_t> iterations;                // row keys
    std::vector<std::string> column_names;               // e.g., L9c1, L9c2, ...
    std::vector<std::vector<std::string>> columns;       // empty string => missing
};

struct PhaseLogWideData {
    std::set<std::int64_t> iterations;
    std::set<int> copies;
    std::map<std::int64_t, std::map<int, std::string>> iter_copy_to_phase;
};

inline bool parse_copy_number(const std::string& token, int& copy_num) {
    const std::size_t cpos = token.rfind('c');
    if (cpos == std::string::npos || cpos + 1 >= token.size()) return false;
    for (std::size_t i = cpos + 1; i < token.size(); ++i) {
        if (std::isdigit(static_cast<unsigned char>(token[i])) == 0) return false;
    }
    copy_num = std::stoi(token.substr(cpos + 1));
    return true;
}

inline PhaseLogWideData read_phase_log_for_wide_table(const std::string& filepath) {
    std::ifstream in(filepath);
    if (!in) throw std::runtime_error("Failed to open phase file: " + filepath);

    std::string line;
    if (!std::getline(in, line)) throw std::runtime_error("Phase file is empty: " + filepath);

    const std::vector<std::string> header = split_tsv_line(line);
    if (header.size() < 2) {
        throw std::runtime_error("Phase file must have Iteration + phase columns: " + filepath);
    }
    if (normalize_header(header[0]) != "iteration") {
        throw std::runtime_error("First column must be Iteration: " + filepath);
    }

    std::vector<std::string> phase_labels;
    phase_labels.reserve(header.size() - 1);
    for (std::size_t j = 1; j < header.size(); ++j) {
        const std::string phase = trim(header[j]);
        if (phase.empty()) {
            throw std::runtime_error("Empty phase-column name in file: " + filepath);
        }
        phase_labels.push_back(phase);
    }

    PhaseLogWideData out;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> row = split_tsv_line(line);
        if (row.empty()) continue;

        const std::int64_t iter = parse_int64(row[0]);
        out.iterations.insert(iter);

        for (std::size_t j = 0; j < phase_labels.size(); ++j) {
            if (j + 1 >= row.size()) continue;
            const std::string value = trim(row[j + 1]);
            if (is_blank_value(value)) continue;

            int copy_num = -1;
            if (!parse_copy_number(value, copy_num)) continue;

            auto& copy_map = out.iter_copy_to_phase[iter];
            auto it_existing = copy_map.find(copy_num);
            if (it_existing != copy_map.end() && it_existing->second != phase_labels[j]) {
                throw std::runtime_error(
                    "Copy c" + std::to_string(copy_num) + " appears in multiple phase columns at iteration " +
                    std::to_string(iter) + " in file: " + filepath
                );
            }
            copy_map[copy_num] = phase_labels[j];
            out.copies.insert(copy_num);
        }
    }

    return out;
}

inline RelabeledPhaseTable relabel_phase_logs_to_wide_table(
    const std::vector<std::string>& tags,
    const std::vector<std::string>& filepaths
) {
    if (tags.size() != filepaths.size()) {
        throw std::runtime_error("tags and filepaths must have the same length.");
    }
    if (tags.empty()) {
        throw std::runtime_error("At least one phase file must be provided.");
    }

    for (std::size_t i = 0; i < tags.size(); ++i) {
        if (tags[i].empty()) throw std::runtime_error("Encountered empty tag.");
        if (filepaths[i].empty()) throw std::runtime_error("Encountered empty phase file path.");
    }

    std::vector<PhaseLogWideData> per_file;
    per_file.reserve(filepaths.size());
    std::set<std::int64_t> all_iterations;
    for (const auto& path : filepaths) {
        PhaseLogWideData dat = read_phase_log_for_wide_table(path);
        all_iterations.insert(dat.iterations.begin(), dat.iterations.end());
        per_file.push_back(std::move(dat));
    }

    RelabeledPhaseTable out;
    out.iterations.assign(all_iterations.begin(), all_iterations.end());
    const std::size_t n_rows = out.iterations.size();

    for (std::size_t f = 0; f < per_file.size(); ++f) {
        const std::string& tag = tags[f];
        const PhaseLogWideData& dat = per_file[f];

        for (int copy_num : dat.copies) {
            const std::string col_name = tag + "c" + std::to_string(copy_num);
            std::vector<std::string> col(n_rows, "");

            for (std::size_t r = 0; r < n_rows; ++r) {
                const std::int64_t iter = out.iterations[r];
                auto it_iter = dat.iter_copy_to_phase.find(iter);
                if (it_iter == dat.iter_copy_to_phase.end()) continue;
                auto it_copy = it_iter->second.find(copy_num);
                if (it_copy == it_iter->second.end()) continue;
                col[r] = it_copy->second;
            }

            out.column_names.push_back(col_name);
            out.columns.push_back(std::move(col));
        }
    }

    return out;
}

struct LabelSwapResult {
    RelabeledPhaseTable relabeled_table;
    std::vector<std::string> labels;                    // source label universe
    std::vector<std::vector<int>> relabel_map_index;    // iteration -> old_label_idx -> new_label_idx
};

inline std::vector<std::vector<int>> group_set_indices_by_size(
    const std::vector<std::vector<int>>& sets_as_idx
) {
    std::map<std::size_t, std::vector<int>> by_size;
    for (std::size_t i = 0; i < sets_as_idx.size(); ++i) {
        by_size[sets_as_idx[i].size()].push_back(static_cast<int>(i));
    }
    std::vector<std::vector<int>> out;
    out.reserve(by_size.size());
    for (auto& kv : by_size) out.push_back(kv.second);
    return out;
}

inline void enumerate_allowed_label_maps(
    const std::vector<std::vector<int>>& sets_as_idx,
    std::size_t n_labels,
    std::vector<std::vector<int>>& out_maps
) {
    const std::vector<std::vector<int>> buckets = group_set_indices_by_size(sets_as_idx);
    std::vector<int> cur_map(n_labels, -1);

    std::function<void(std::size_t)> recurse_bucket;
    recurse_bucket = [&](std::size_t bidx) {
        if (bidx == buckets.size()) {
            out_maps.push_back(cur_map);
            return;
        }

        std::vector<int> src_sets = buckets[bidx];
        std::vector<int> dst_sets = buckets[bidx];
        std::sort(dst_sets.begin(), dst_sets.end());

        do {
            std::function<void(std::size_t)> recurse_set_in_bucket;
            recurse_set_in_bucket = [&](std::size_t sidx) {
                if (sidx == src_sets.size()) {
                    recurse_bucket(bidx + 1);
                    return;
                }

                const int src_set_id = src_sets[sidx];
                const int dst_set_id = dst_sets[sidx];
                const std::vector<int>& src_labels = sets_as_idx[(std::size_t)src_set_id];
                std::vector<int> dst_labels = sets_as_idx[(std::size_t)dst_set_id];
                std::sort(dst_labels.begin(), dst_labels.end());

                do {
                    std::vector<int> touched;
                    touched.reserve(src_labels.size());
                    for (std::size_t k = 0; k < src_labels.size(); ++k) {
                        const int src_lab = src_labels[k];
                        const int dst_lab = dst_labels[k];
                        cur_map[(std::size_t)src_lab] = dst_lab;
                        touched.push_back(src_lab);
                    }

                    recurse_set_in_bucket(sidx + 1);

                    for (int idx : touched) cur_map[(std::size_t)idx] = -1;
                } while (std::next_permutation(dst_labels.begin(), dst_labels.end()));
            };

            recurse_set_in_bucket(0);
        } while (std::next_permutation(dst_sets.begin(), dst_sets.end()));
    };

    recurse_bucket(0);
}

inline int mismatch_cost_between_rows(
    const std::vector<int>& row_prev,
    const std::vector<int>& row_cur,
    const std::vector<int>& map_prev,
    const std::vector<int>& map_cur
) {
    int cost = 0;
    for (std::size_t c = 0; c < row_prev.size(); ++c) {
        const int a = row_prev[c];
        const int b = row_cur[c];
        if (a < 0 || b < 0) continue;
        const int ma = map_prev[(std::size_t)a];
        const int mb = map_cur[(std::size_t)b];
        if (ma != mb) cost += 1;
    }
    return cost;
}

inline LabelSwapResult relabel_wide_table_by_allelic_sets(
    const RelabeledPhaseTable& input,
    const std::vector<std::vector<std::string>>& allelic_sets
) {
    if (input.iterations.empty()) {
        throw std::runtime_error("Input table has no rows.");
    }
    if (input.columns.size() != input.column_names.size()) {
        throw std::runtime_error("Input table has inconsistent column names/columns.");
    }
    for (const auto& col : input.columns) {
        if (col.size() != input.iterations.size()) {
            throw std::runtime_error("Input table has columns with inconsistent row counts.");
        }
    }
    if (allelic_sets.empty()) {
        throw std::runtime_error("allelic_tips must contain at least one set.");
    }

    std::set<std::string> observed_labels;
    for (const auto& col : input.columns) {
        for (const auto& v : col) {
            if (!v.empty()) observed_labels.insert(v);
        }
    }
    if (observed_labels.empty()) {
        throw std::runtime_error("Input table has no non-missing phase labels.");
    }

    std::set<std::string> partition_labels;
    for (const auto& set_labels : allelic_sets) {
        if (set_labels.empty()) {
            throw std::runtime_error("allelic_tips contains an empty label set.");
        }
        for (const auto& lab : set_labels) {
            if (lab.empty()) throw std::runtime_error("allelic_tips contains an empty label.");
            if (!partition_labels.insert(lab).second) {
                throw std::runtime_error("Label appears more than once in allelic_tips: " + lab);
            }
        }
    }
    if (partition_labels != observed_labels) {
        throw std::runtime_error(
            "allelic_tips labels must exactly match labels present in the input table."
        );
    }

    std::vector<std::string> labels;
    labels.reserve(partition_labels.size());
    for (const auto& set_labels : allelic_sets) {
        for (const auto& lab : set_labels) labels.push_back(lab);
    }

    std::unordered_map<std::string, int> label_to_idx;
    label_to_idx.reserve(labels.size());
    for (std::size_t i = 0; i < labels.size(); ++i) {
        label_to_idx[labels[i]] = static_cast<int>(i);
    }

    std::vector<std::vector<int>> sets_as_idx;
    sets_as_idx.reserve(allelic_sets.size());
    for (const auto& set_labels : allelic_sets) {
        std::vector<int> idx_set;
        idx_set.reserve(set_labels.size());
        for (const auto& lab : set_labels) idx_set.push_back(label_to_idx.at(lab));
        sets_as_idx.push_back(std::move(idx_set));
    }

    std::vector<std::vector<int>> allowed_maps;
    enumerate_allowed_label_maps(sets_as_idx, labels.size(), allowed_maps);
    if (allowed_maps.empty()) {
        throw std::runtime_error("No allowed relabeling maps were generated.");
    }

    int identity_state = -1;
    for (std::size_t s = 0; s < allowed_maps.size(); ++s) {
        bool is_identity = true;
        for (std::size_t i = 0; i < labels.size(); ++i) {
            if (allowed_maps[s][i] != static_cast<int>(i)) {
                is_identity = false;
                break;
            }
        }
        if (is_identity) {
            identity_state = static_cast<int>(s);
            break;
        }
    }
    if (identity_state < 0) {
        throw std::runtime_error("Identity relabeling not found in allowed map set.");
    }

    const std::size_t T = input.iterations.size();
    const std::size_t S = allowed_maps.size();
    const std::size_t C = input.columns.size();

    std::vector<std::vector<int>> rows(T, std::vector<int>(C, -1));
    for (std::size_t c = 0; c < C; ++c) {
        for (std::size_t t = 0; t < T; ++t) {
            const std::string& v = input.columns[c][t];
            if (v.empty()) continue;
            rows[t][c] = label_to_idx.at(v);
        }
    }

    const int INF = std::numeric_limits<int>::max() / 4;
    std::vector<std::vector<int>> dp(T, std::vector<int>(S, INF));
    std::vector<std::vector<int>> back(T, std::vector<int>(S, -1));

    dp[0][(std::size_t)identity_state] = 0;
    for (std::size_t t = 1; t < T; ++t) {
        for (std::size_t s_cur = 0; s_cur < S; ++s_cur) {
            int best = INF;
            int best_prev = -1;
            for (std::size_t s_prev = 0; s_prev < S; ++s_prev) {
                if (dp[t - 1][s_prev] >= INF) continue;
                const int trans =
                    mismatch_cost_between_rows(rows[t - 1], rows[t], allowed_maps[s_prev], allowed_maps[s_cur]);
                const int cand = dp[t - 1][s_prev] + trans;
                if (cand < best) {
                    best = cand;
                    best_prev = static_cast<int>(s_prev);
                }
            }
            dp[t][s_cur] = best;
            back[t][s_cur] = best_prev;
        }
    }

    int best_last_state = 0;
    int best_last_cost = dp[T - 1][0];
    for (std::size_t s = 1; s < S; ++s) {
        if (dp[T - 1][s] < best_last_cost) {
            best_last_cost = dp[T - 1][s];
            best_last_state = static_cast<int>(s);
        }
    }

    std::vector<int> chosen_state(T, identity_state);
    chosen_state[T - 1] = best_last_state;
    for (std::size_t t = T - 1; t > 0; --t) {
        chosen_state[t - 1] = back[t][(std::size_t)chosen_state[t]];
    }

    RelabeledPhaseTable out = input;
    for (std::size_t t = 0; t < T; ++t) {
        const std::vector<int>& map_idx = allowed_maps[(std::size_t)chosen_state[t]];
        for (std::size_t c = 0; c < C; ++c) {
            const int raw_idx = rows[t][c];
            if (raw_idx < 0) {
                out.columns[c][t] = "";
            } else {
                out.columns[c][t] = labels[(std::size_t)map_idx[(std::size_t)raw_idx]];
            }
        }
    }

    LabelSwapResult result;
    result.relabeled_table = std::move(out);
    result.labels = std::move(labels);
    result.relabel_map_index.reserve(T);
    for (std::size_t t = 0; t < T; ++t) {
        result.relabel_map_index.push_back(allowed_maps[(std::size_t)chosen_state[t]]);
    }
    return result;
}

} // namespace homologizer

#endif // HOMOLOGIZER_RELABEL_HPP
