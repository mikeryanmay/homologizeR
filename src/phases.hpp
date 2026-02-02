#ifndef HOMOLOGIZER_PHASE_HPP
#define HOMOLOGIZER_PHASE_HPP

#include "utils.hpp"

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace homologizer {

using PhaseRow   = std::vector<std::string>;           // values aligned to labels
using PhaseTable = std::map<std::int64_t, PhaseRow>;   // iteration -> row

struct PhaseFile {
    std::string tag;               // e.g. "L10"
    std::vector<std::string> labels;
    PhaseTable table;
};

inline PhaseFile read_phase_file_tsv(const std::string& tag,
                                     const std::string& filepath) {
    std::ifstream in(filepath);
    if (!in) throw std::runtime_error("Failed to open phase file: " + filepath);

    std::string line;

    if (!std::getline(in, line))
        throw std::runtime_error("Phase file is empty: " + filepath);

    auto header = split_tsv_line(line);
    if (header.size() < 2)
        throw std::runtime_error("Phase file header must have Iteration + >=1 label columns: " + filepath);

    if (normalize_header(header[0]) != "iteration")
        throw std::runtime_error("First column must be 'Iteration' in phase file: " + filepath);

    PhaseFile pf;
    pf.tag = tag;
    pf.labels.reserve(header.size() - 1);

    for (std::size_t j = 1; j < header.size(); ++j) {
        std::string lab = trim(header[j]);
        if (lab.empty())
            throw std::runtime_error("Empty label in header col " + std::to_string(j + 1) + " in " + filepath);
        pf.labels.push_back(lab);
    }

    while (std::getline(in, line)) {
        if (line.empty()) continue;
        auto row = split_tsv_line(line);
        if (row.empty()) continue;

        std::int64_t iter = parse_int64(row[0]);

        PhaseRow vals(pf.labels.size(), "");
        for (std::size_t j = 0; j < pf.labels.size(); ++j) {
            if (j + 1 < row.size()) vals[j] = trim(row[j + 1]);
        }

        pf.table[iter] = std::move(vals);
    }

    return pf;
}

inline std::vector<PhaseFile> read_phase_files_checked(
        const std::vector<std::string>& tags,
        const std::vector<std::string>& filepaths
) {
    if (tags.size() != filepaths.size()) {
        throw std::runtime_error("tags and filepaths must have the same length");
    }
    if (filepaths.empty()) return {};

    std::vector<PhaseFile> files;
    files.reserve(filepaths.size());

    PhaseFile first = read_phase_file_tsv(tags[0], filepaths[0]);
    files.push_back(std::move(first));
    const auto& ref_labels = files[0].labels;

    for (std::size_t i = 1; i < filepaths.size(); ++i) {
        PhaseFile pf = read_phase_file_tsv(tags[i], filepaths[i]);
        if (pf.labels != ref_labels) {
            throw std::runtime_error("Phase file headers do not match across files. Mismatch in: " + filepaths[i]);
        }
        files.push_back(std::move(pf));
    }

    return files;
}

inline std::set<std::int64_t> intersect_iterations(const std::vector<PhaseFile>& phase_files) {
    std::set<std::int64_t> inter;
    if (phase_files.empty()) return inter;

    for (const auto& kv : phase_files[0].table) inter.insert(kv.first);

    for (std::size_t i = 1; i < phase_files.size(); ++i) {
        std::set<std::int64_t> next;
        for (const auto& kv : phase_files[i].table) next.insert(kv.first);

        std::set<std::int64_t> tmp;
        std::set_intersection(inter.begin(), inter.end(),
                              next.begin(), next.end(),
                              std::inserter(tmp, tmp.begin()));
        inter.swap(tmp);

        if (inter.empty()) break;
    }
    return inter;
}

// Core worker: takes already-read PhaseFile objects
inline std::map<std::int64_t, std::map<std::string, std::vector<std::string>>>
    collect_phase_labels_all_iterations_intersection(
        const std::vector<PhaseFile>& phase_files,
        bool unique_values = false
    ) {
        using LabelToValues = std::map<std::string, std::vector<std::string>>;
        std::map<std::int64_t, LabelToValues> out;
        if (phase_files.empty()) return out;

        const auto& labels = phase_files[0].labels;
        std::set<std::int64_t> iters = intersect_iterations(phase_files);

        for (std::int64_t iter : iters) {
            LabelToValues by_label;
            for (const auto& lab : labels) by_label[lab] = {};

            std::map<std::string, std::unordered_set<std::string>> seen;
            if (unique_values) {
                for (const auto& lab : labels) seen[lab] = {};
            }

            for (const auto& pf : phase_files) {
                const PhaseRow& vals = pf.table.at(iter);

                for (std::size_t j = 0; j < labels.size(); ++j) {
                    const std::string& raw = vals[j];
                    if (is_blank_value(raw)) continue;

                    // ✅ append tag to value
                    std::string value = append_tag(raw, pf.tag);

                    if (unique_values) {
                        auto& s = seen[labels[j]];
                        if (s.insert(value).second) by_label[labels[j]].push_back(std::move(value));
                    } else {
                        by_label[labels[j]].push_back(std::move(value));
                    }
                }
            }

            out.emplace(iter, std::move(by_label));
        }

        return out;
    }

// Convenience overload: take tags + filepaths directly (no intermediate objects)
inline std::map<std::int64_t, std::map<std::string, std::vector<std::string>>>
    collect_phase_labels_all_iterations_intersection(
        const std::vector<std::string>& tags,
        const std::vector<std::string>& filepaths,
        bool unique_values = false
    ) {
        std::vector<PhaseFile> files = read_phase_files_checked(tags, filepaths);
        return collect_phase_labels_all_iterations_intersection(files, unique_values);
    }

} // namespace homologizer

#endif // HOMOLOGIZER_PHASE_HPP
