#ifndef HOMOLOGIZER_UTILS_HPP
#define HOMOLOGIZER_UTILS_HPP

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

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

} // namespace homologizer

#endif // HOMOLOGIZER_UTILS_HPP
