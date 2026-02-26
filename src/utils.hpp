#ifndef HOMOLOGIZER_UTILS_HPP
#define HOMOLOGIZER_UTILS_HPP

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>
#include <fstream>

namespace homologizer {

inline std::string trim(std::string s) {
    auto not_space = [](unsigned char c){ return !std::isspace(c); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), not_space));
    s.erase(std::find_if(s.rbegin(), s.rend(), not_space).base(), s.end());
    return s;
}

inline std::string to_lower(std::string s) {
    for (char& c : s) c = (char)std::tolower((unsigned char)c);
    return s;
}

inline std::vector<std::string> split_tsv_line(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    for (char c : line) {
        if (c == '\t') { out.push_back(cur); cur.clear(); }
        else cur.push_back(c);
    }
    out.push_back(cur);
    return out;
}

inline std::int64_t parse_int64(const std::string& s) {
    std::string t = trim(s);
    if (t.empty()) throw std::runtime_error("Empty integer field");
    char* endptr = nullptr;
    long long v = std::strtoll(t.c_str(), &endptr, 10);
    if (endptr == t.c_str() || *endptr != '\0') {
        throw std::runtime_error("Invalid integer token: '" + t + "'");
    }
    return (std::int64_t)v;
}

// Normalize header for comparisons: lowercase + alnum only
inline std::string normalize_header(std::string x) {
    x = to_lower(trim(x));
    std::string out;
    out.reserve(x.size());
    for (char c : x) {
        unsigned char uc = (unsigned char)c;
        if (std::isalnum(uc)) out.push_back((char)uc);
    }
    return out;
}

// Blank predicate for phase values (customize if needed)
inline bool is_blank_value(const std::string& raw) {
    std::string v = trim(raw);
    if (v.empty()) return true;

    std::string vl = to_lower(v);
    if (vl == "na" || vl == "nan" || vl == "null") return true;

    // Your example uses IDs containing "BLANK"
    if (vl.find("blank") != std::string::npos) return true;

    return false;
}

inline std::size_t find_column_by_name(
        const std::vector<std::string>& header,
        const std::vector<std::string>& candidates
) {
    std::vector<std::string> cand_norm;
    cand_norm.reserve(candidates.size());
    for (const auto& c : candidates) cand_norm.push_back(normalize_header(c));

    for (std::size_t j = 0; j < header.size(); ++j) {
        std::string hj = normalize_header(header[j]);
        for (const auto& cn : cand_norm) {
            if (hj == cn) return j;
        }
    }
    throw std::runtime_error("Required column not found in header");
}

// How we append the phase-file tag to a sample label.
// You said “append those names”; this does that as: sample + "|" + tag.
// Change the delimiter if you prefer.
inline std::string append_tag(const std::string& sample, const std::string& tag) {
    if (tag.empty()) return sample;
    return sample + "_" + tag;
}

// Read first nonempty line for MAP Newick
static std::string read_single_newick_line(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open MAP tree file: " + path);
    std::string line;
    while (std::getline(in, line)) {
        line = trim(line);
        if (!line.empty()) return line;
    }
    throw std::runtime_error("MAP tree file is empty: " + path);
}

/// Strip NEXUS-style bracket comments: [ ... ]
/// Non-nested (sufficient for typical TREE blocks).
inline std::string strip_nexus_bracket_comments(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    bool in_comment = false;
    for (char c : s) {
        if (!in_comment) {
            if (c == '[') {
                in_comment = true;
            } else {
                out.push_back(c);
            }
        } else {
            if (c == ']') {
                in_comment = false;
            }
        }
    }
    return out;
}

/// Read exactly one Newick tree from a NEXUS file.
/// Assumes the file contains exactly one `tree ... = <newick>;` entry.
/// Returns the Newick string including the trailing `;`.
inline std::string read_single_newick_from_nexus_file(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Cannot open NEXUS file: " + path);
    }

    std::string line;
    bool found_tree = false;
    std::ostringstream newick;

    while (std::getline(in, line)) {
        // Remove comments early
        std::string no_comments = strip_nexus_bracket_comments(line);

        // Use your existing trim helper
        std::string t = trim(no_comments);
        if (t.empty()) continue;

        if (!found_tree) {
            // Case-insensitive check for leading "tree"
            std::string lower = t;
            std::transform(lower.begin(), lower.end(), lower.begin(),
                           [](unsigned char c) { return std::tolower(c); });

            if (lower.rfind("tree", 0) == 0) {
                auto eq = t.find('=');
                if (eq == std::string::npos) continue;

                found_tree = true;
                std::string after = trim(t.substr(eq + 1));
                newick << after;

                // Stop immediately if we already hit ';'
                if (after.find(';') != std::string::npos) break;
            }
        } else {
            // Continuation lines until ';'
            newick << ' ' << t;
            if (t.find(';') != std::string::npos) break;
        }
    }

    if (!found_tree) {
        throw std::runtime_error(
            "No `tree ... = <newick>;` entry found in NEXUS file: " + path
        );
    }

    std::string nwk = strip_nexus_bracket_comments(newick.str());
    nwk = trim(nwk);

    auto semi = nwk.find(';');
    if (semi == std::string::npos) {
        throw std::runtime_error(
            "Tree entry found but missing terminating ';' in NEXUS file: " + path
        );
    }

    nwk = trim(nwk.substr(0, semi + 1));

    // Very cheap sanity check
    if (nwk.find('(') == std::string::npos) {
        throw std::runtime_error(
            "Extracted tree does not look like Newick (no '(' found): " + path
        );
    }

    return nwk;
}

} // namespace homologizer

#endif // HOMOLOGIZER_UTILS_HPP

