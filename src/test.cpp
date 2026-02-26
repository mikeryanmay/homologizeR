#include <Rcpp.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <regex>
#include <unordered_map>
#include <unordered_set>
#include "trees.hpp"
#include "phases.hpp"
#include "core_set.hpp"
#include "relabel.hpp"

using namespace Rcpp;
using namespace homologizer;

typedef std::map<std::int64_t, Tree> treeMap;
typedef std::map<std::int64_t, std::map<std::string, std::vector<std::string>>> phaseMap;

namespace {

std::string join_tsv_line(const std::vector<std::string>& fields) {
    if (fields.empty()) return "";
    std::string out = fields[0];
    for (std::size_t i = 1; i < fields.size(); ++i) {
        out.push_back('\t');
        out += fields[i];
    }
    return out;
}

std::string regex_escape(const std::string& s) {
    std::string out;
    out.reserve(s.size() * 2);
    for (char ch : s) {
        switch (ch) {
            case '\\':
            case '^':
            case '$':
            case '.':
            case '|':
            case '?':
            case '*':
            case '+':
            case '(':
            case ')':
            case '[':
            case ']':
            case '{':
            case '}':
                out.push_back('\\');
                out.push_back(ch);
                break;
            default:
                out.push_back(ch);
                break;
        }
    }
    return out;
}

void replace_tip_label_token(std::string& newick,
                             const std::string& from,
                             const std::string& to) {
    if (from == to) return;
    const std::string pat = "([\\(,])" + regex_escape(from) + "(?=\\[|:|,|\\)|;)";
    const std::regex re(pat);
    const std::string rep = "$1" + to;
    newick = std::regex_replace(newick, re, rep);
}

std::vector<std::vector<std::string>> parse_allelic_tips_list(const Rcpp::List& allelic_tips) {
    if (allelic_tips.size() == 0) {
        throw std::runtime_error("allelic_tips must contain at least one set.");
    }
    Rcpp::CharacterVector set_names = allelic_tips.names();
    if (set_names.size() != allelic_tips.size()) {
        throw std::runtime_error("allelic_tips must be a named list.");
    }

    std::vector<std::vector<std::string>> sets;
    sets.reserve((std::size_t)allelic_tips.size());
    for (R_xlen_t i = 0; i < allelic_tips.size(); ++i) {
        if (set_names[i] == NA_STRING) {
            throw std::runtime_error("allelic_tips has an NA set name.");
        }
        Rcpp::CharacterVector labs = allelic_tips[i];
        if (labs.size() == 0) {
            throw std::runtime_error("allelic_tips contains an empty set.");
        }
        std::vector<std::string> one_set;
        one_set.reserve((std::size_t)labs.size());
        for (R_xlen_t k = 0; k < labs.size(); ++k) {
            if (labs[k] == NA_STRING) throw std::runtime_error("allelic_tips contains an NA label.");
            std::string lab = Rcpp::as<std::string>(labs[k]);
            if (lab.empty()) throw std::runtime_error("allelic_tips contains an empty label.");
            one_set.push_back(lab);
        }
        sets.push_back(std::move(one_set));
    }
    return sets;
}

struct ParsedAllelicTips {
    std::vector<std::string> group_names;
    std::vector<std::vector<std::string>> groups;
};

ParsedAllelicTips parse_allelic_tips_with_names(const Rcpp::List& allelic_tips) {
    if (allelic_tips.size() == 0) {
        throw std::runtime_error("allelic_tips must contain at least one set.");
    }
    Rcpp::CharacterVector set_names = allelic_tips.names();
    if (set_names.size() != allelic_tips.size()) {
        throw std::runtime_error("allelic_tips must be a named list.");
    }

    ParsedAllelicTips out;
    out.group_names.reserve((std::size_t)allelic_tips.size());
    out.groups.reserve((std::size_t)allelic_tips.size());
    for (R_xlen_t i = 0; i < allelic_tips.size(); ++i) {
        if (set_names[i] == NA_STRING) {
            throw std::runtime_error("allelic_tips has an NA set name.");
        }
        std::string set_name = Rcpp::as<std::string>(set_names[i]);
        if (set_name.empty()) {
            throw std::runtime_error("allelic_tips has an empty set name.");
        }
        Rcpp::CharacterVector labs = allelic_tips[i];
        if (labs.size() == 0) {
            throw std::runtime_error("allelic_tips contains an empty set.");
        }
        std::vector<std::string> one_set;
        one_set.reserve((std::size_t)labs.size());
        for (R_xlen_t k = 0; k < labs.size(); ++k) {
            if (labs[k] == NA_STRING) throw std::runtime_error("allelic_tips contains an NA label.");
            std::string lab = Rcpp::as<std::string>(labs[k]);
            if (lab.empty()) throw std::runtime_error("allelic_tips contains an empty label.");
            one_set.push_back(lab);
        }
        out.group_names.push_back(set_name);
        out.groups.push_back(std::move(one_set));
    }
    return out;
}

void write_label_marginals_file(const RelabeledPhaseTable& table,
                                const ParsedAllelicTips& parsed,
                                const std::string& outpath) {
    std::ofstream out(outpath);
    if (!out) throw std::runtime_error("Failed to open output marginals file for writing.");
    out << std::fixed << std::setprecision(6);

    std::vector<std::string> labels;
    std::unordered_map<std::string, int> label_idx;
    for (const auto& g : parsed.groups) {
        for (const auto& lab : g) {
            if (label_idx.find(lab) != label_idx.end()) {
                throw std::runtime_error("Label appears more than once in allelic_tips: " + lab);
            }
            label_idx[lab] = static_cast<int>(labels.size());
            labels.push_back(lab);
        }
    }
    if (labels.empty()) throw std::runtime_error("No labels found in allelic_tips.");

    out << "locus_copy";
    for (const auto& lab : labels) out << "\t" << lab;
    for (const auto& gname : parsed.group_names) out << "\t" << gname;
    out << "\n";

    const std::size_t n_rows = table.iterations.size();
    for (std::size_t c = 0; c < table.column_names.size(); ++c) {
        std::vector<int> counts(labels.size(), 0);
        int observed = 0;

        for (std::size_t r = 0; r < n_rows; ++r) {
            const std::string& v = table.columns[c][r];
            if (v.empty()) continue;
            auto it = label_idx.find(v);
            if (it == label_idx.end()) {
                throw std::runtime_error("Found label not present in allelic_tips: " + v);
            }
            counts[(std::size_t)it->second] += 1;
            observed += 1;
        }

        out << table.column_names[c];
        if (observed == 0) {
            for (std::size_t i = 0; i < labels.size(); ++i) out << "\tNA";
            for (std::size_t g = 0; g < parsed.groups.size(); ++g) out << "\tNA";
            out << "\n";
            continue;
        }

        for (std::size_t i = 0; i < labels.size(); ++i) {
            const double p = static_cast<double>(counts[i]) / static_cast<double>(observed);
            out << "\t" << p;
        }
        for (const auto& grp : parsed.groups) {
            int sum = 0;
            for (const auto& lab : grp) sum += counts[(std::size_t)label_idx.at(lab)];
            const double p = static_cast<double>(sum) / static_cast<double>(observed);
            out << "\t" << p;
        }
        out << "\n";
    }
}

void write_locus_manifest_file(const RelabeledPhaseTable& table,
                               const std::string& outpath) {
    std::ofstream out(outpath);
    if (!out) throw std::runtime_error("Failed to open output locus manifest file for writing.");

    out << "locus_copy\tlocus\tcopy\n";
    for (const auto& col_name : table.column_names) {
        std::string locus = col_name;
        std::string copy = "";

        const std::size_t cpos = col_name.rfind('c');
        if (cpos != std::string::npos && cpos + 1 < col_name.size()) {
            bool numeric_suffix = true;
            for (std::size_t i = cpos + 1; i < col_name.size(); ++i) {
                if (std::isdigit(static_cast<unsigned char>(col_name[i])) == 0) {
                    numeric_suffix = false;
                    break;
                }
            }
            if (numeric_suffix) {
                locus = col_name.substr(0, cpos);
                copy = col_name.substr(cpos);
            }
        }

        out << col_name << "\t" << locus << "\t" << copy << "\n";
    }
}

RelabeledPhaseTable apply_burnin_fraction(const RelabeledPhaseTable& input, double burnin_prop) {
    if (!(burnin_prop >= 0.0 && burnin_prop < 1.0)) {
        throw std::runtime_error("burnin must be in [0, 1).");
    }
    const std::size_t n = input.iterations.size();
    if (n == 0) {
        throw std::runtime_error("Cannot apply burnin to an empty table.");
    }

    const std::size_t drop = static_cast<std::size_t>(std::floor(burnin_prop * static_cast<double>(n)));
    if (drop >= n) {
        throw std::runtime_error("burnin discards all iterations.");
    }
    if (drop == 0) return input;

    RelabeledPhaseTable out;
    out.iterations.assign(input.iterations.begin() + drop, input.iterations.end());
    out.column_names = input.column_names;
    out.columns.reserve(input.columns.size());

    for (const auto& col : input.columns) {
        out.columns.emplace_back(
            col.begin() + drop,
            col.end()
        );
    }

    return out;
}

} // namespace

// [[Rcpp::export]]
void phaseLabelTrees(std::string tree_file, CharacterVector phase_files, std::string output_file) {

    // read the tree file
    treeMap trees = read_tree_tsv_by_iteration(tree_file);

    // extract phase file paths + tags from named character vector
    std::vector<std::string> paths = as<std::vector<std::string>>(phase_files);
    CharacterVector nms = phase_files.names();
    if (nms.size() != phase_files.size()) {
        stop("phase_files must be a *named* character vector (names are tags like L10, L100, ...).");
    }

    std::vector<std::string> tags;
    tags.reserve((std::size_t)nms.size());
    for (R_xlen_t i = 0; i < nms.size(); ++i) {
        if (nms[i] == NA_STRING) stop("phase_files has an NA name at position %d", (int)i + 1);
        std::string tag = as<std::string>(nms[i]);
        if (tag.empty()) stop("phase_files has an empty name at position %d", (int)i + 1);
        tags.push_back(tag);
    }

    // read + aggregate phase files; values become "sample|TAG"
    phaseMap phasings = collect_phase_labels_all_iterations_intersection(tags, paths);

    // check that each iteration in trees exists in the phasings
    for (auto it = trees.begin(); it != trees.end(); ++it) {
        if (phasings.find(it->first) == phasings.end()) {
            stop("Found trees for iterations that are not present in the phasing files. "
                     "Make sure each sampled tree has a matching sampled phasing.");
        }
    }

    // for each sampled tree
    for (auto it = trees.begin(); it != trees.end(); ++it) {
        // relabel the trees
        Tree& this_tree = it->second;
        const auto& this_phasing = phasings.at(it->first);
        this_tree.apply_phasing(this_phasing);
    }

    // write to file
    write_tree_map_tsv(output_file, trees);

}

// [[Rcpp::export]]
Rcpp::List identifyCorePhasedSets(std::string phase_labeled_tree_file,
                                  double C = 0.5,
                                  bool aggregate_alleles = true) {
    
    if (phase_labeled_tree_file.empty()) {
        stop("phase_labeled_tree_file must be a non-empty path.");
    }
    
    if (!(C > 0.0 && C <= 1.0)) {
        stop("C must be in (0, 1].");
    }

    Rcpp::List result = posterior_core_sets_from_phase_labeled_trees(
        phase_labeled_tree_file,
        C,
        aggregate_alleles
    );

	    return Rcpp::List::create(
	        Rcpp::Named("core_sets") = result["core_sets"],
	        Rcpp::Named("core_set_marginal_probability") = result["core_set_marginal_probability"],
	        Rcpp::Named("co_phasing_probability") = result["co_phasing_probability"]
	    );

	}

// [[Rcpp::export]]
Rcpp::DataFrame relabelPhaseLogs(Rcpp::CharacterVector phase_files) {
    if (phase_files.size() == 0) {
        stop("phase_files must contain at least one file path.");
    }

    Rcpp::CharacterVector nms = phase_files.names();
    if (nms.size() != phase_files.size()) {
        stop("phase_files must be a named character vector.");
    }

    // Validate names/paths here so errors are immediate at the Rcpp boundary.
    for (R_xlen_t i = 0; i < phase_files.size(); ++i) {
        if (nms[i] == NA_STRING) stop("phase_files has an NA name at position %d.", (int)i + 1);
        std::string tag = as<std::string>(nms[i]);
        if (tag.empty()) stop("phase_files has an empty name at position %d.", (int)i + 1);

        if (phase_files[i] == NA_STRING) stop("phase_files has an NA path at position %d.", (int)i + 1);
        std::string path = as<std::string>(phase_files[i]);
        if (path.empty()) stop("phase_files has an empty path at position %d.", (int)i + 1);
    }

    // Convert named R vector into STL inputs expected by relabel.hpp.
    std::vector<std::string> tags;
    std::vector<std::string> paths;
    tags.reserve((std::size_t)phase_files.size());
    paths.reserve((std::size_t)phase_files.size());
    for (R_xlen_t i = 0; i < phase_files.size(); ++i) {
        tags.push_back(as<std::string>(nms[i]));
        paths.push_back(as<std::string>(phase_files[i]));
    }

    // Run pure C++ relabel-table construction.
    RelabeledPhaseTable table = relabel_phase_logs_to_wide_table(tags, paths);

    // Convert STL table to an R data.frame with Iteration + locus-copy columns.
    const int n_rows = static_cast<int>(table.iterations.size());
    Rcpp::List out_cols;
    Rcpp::CharacterVector out_names;

    Rcpp::NumericVector iter_col(n_rows);
    for (int r = 0; r < n_rows; ++r) {
        iter_col[r] = static_cast<double>(table.iterations[(std::size_t)r]);
    }
    out_cols.push_back(iter_col);
    out_names.push_back("Iteration");

    for (std::size_t c = 0; c < table.column_names.size(); ++c) {
        Rcpp::CharacterVector col(n_rows, NA_STRING);
        for (int r = 0; r < n_rows; ++r) {
            const std::string& v = table.columns[c][(std::size_t)r];
            if (!v.empty()) col[r] = v;
        }
        out_cols.push_back(col);
        out_names.push_back(table.column_names[c]);
    }

    out_cols.attr("names") = out_names;
    out_cols.attr("class") = "data.frame";
    out_cols.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, -n_rows);
    return Rcpp::DataFrame(out_cols);
    
}

// [[Rcpp::export]]
Rcpp::List alignPhaseLabels(Rcpp::DataFrame phase_table, Rcpp::List allelic_tips) {
    Rcpp::CharacterVector col_names = phase_table.names();
    if (col_names.size() < 2) {
        stop("phase_table must contain Iteration plus at least one locus-copy column.");
    }
    if (Rcpp::as<std::string>(col_names[0]) != "Iteration") {
        stop("phase_table first column must be named 'Iteration'.");
    }

    const int n_rows = phase_table.nrows();
    if (n_rows <= 0) {
        stop("phase_table has no rows.");
    }

    // Convert input data.frame to STL table used by core alignment code.
    RelabeledPhaseTable input;
    input.iterations.reserve((std::size_t)n_rows);

    Rcpp::NumericVector iter_col = phase_table[0];
    for (int r = 0; r < n_rows; ++r) {
        input.iterations.push_back((std::int64_t)iter_col[r]);
    }

    input.column_names.reserve((std::size_t)col_names.size() - 1);
    input.columns.reserve((std::size_t)col_names.size() - 1);
    for (R_xlen_t j = 1; j < col_names.size(); ++j) {
        input.column_names.push_back(Rcpp::as<std::string>(col_names[j]));
        Rcpp::CharacterVector rcol = phase_table[j];
        std::vector<std::string> col((std::size_t)n_rows, "");
        for (int r = 0; r < n_rows; ++r) {
            if (rcol[r] == NA_STRING) continue;
            col[(std::size_t)r] = Rcpp::as<std::string>(rcol[r]);
        }
        input.columns.push_back(std::move(col));
    }

    if (allelic_tips.size() == 0) {
        stop("allelic_tips must contain at least one set.");
    }
    Rcpp::CharacterVector set_names = allelic_tips.names();
    if (set_names.size() != allelic_tips.size()) {
        stop("allelic_tips must be a named list.");
    }

    // Convert allelic_tips list into vector-of-vectors for pure C++ alignment.
    std::vector<std::vector<std::string>> sets;
    sets.reserve((std::size_t)allelic_tips.size());
    for (R_xlen_t i = 0; i < allelic_tips.size(); ++i) {
        if (set_names[i] == NA_STRING) stop("allelic_tips has an NA set name at position %d.", (int)i + 1);
        Rcpp::CharacterVector labs = allelic_tips[i];
        if (labs.size() == 0) stop("allelic_tips contains an empty set at position %d.", (int)i + 1);
        std::vector<std::string> one_set;
        one_set.reserve((std::size_t)labs.size());
        for (R_xlen_t k = 0; k < labs.size(); ++k) {
            if (labs[k] == NA_STRING) stop("allelic_tips contains an NA label.");
            std::string lab = Rcpp::as<std::string>(labs[k]);
            if (lab.empty()) stop("allelic_tips contains an empty label.");
            one_set.push_back(lab);
        }
        sets.push_back(std::move(one_set));
    }

    // Compute iteration-wise constrained relabeling that minimizes cross-iteration mismatch.
    LabelSwapResult aligned = relabel_wide_table_by_allelic_sets(input, sets);

    // Convert relabeled table back to an R data.frame.
    Rcpp::List out_cols;
    Rcpp::CharacterVector out_names;

    Rcpp::NumericVector out_iter(n_rows);
    for (int r = 0; r < n_rows; ++r) out_iter[r] = static_cast<double>(aligned.relabeled_table.iterations[(std::size_t)r]);
    out_cols.push_back(out_iter);
    out_names.push_back("Iteration");

    for (std::size_t c = 0; c < aligned.relabeled_table.column_names.size(); ++c) {
        Rcpp::CharacterVector col(n_rows, NA_STRING);
        for (int r = 0; r < n_rows; ++r) {
            const std::string& v = aligned.relabeled_table.columns[c][(std::size_t)r];
            if (!v.empty()) col[r] = v;
        }
        out_cols.push_back(col);
        out_names.push_back(aligned.relabeled_table.column_names[c]);
    }
    out_cols.attr("names") = out_names;
    out_cols.attr("class") = "data.frame";
    out_cols.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, -n_rows);
    Rcpp::DataFrame relabeled_df(out_cols);

    // Build per-iteration old->new label map table for debugging/auditing.
    Rcpp::List map_cols;
    Rcpp::CharacterVector map_names;
    map_cols.push_back(out_iter);
    map_names.push_back("Iteration");
    for (std::size_t lab = 0; lab < aligned.labels.size(); ++lab) {
        Rcpp::CharacterVector col(n_rows, NA_STRING);
        for (int r = 0; r < n_rows; ++r) {
            const int mapped_idx = aligned.relabel_map_index[(std::size_t)r][lab];
            col[r] = aligned.labels[(std::size_t)mapped_idx];
        }
        map_cols.push_back(col);
        map_names.push_back(aligned.labels[lab]);
    }
    map_cols.attr("names") = map_names;
    map_cols.attr("class") = "data.frame";
    map_cols.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, -n_rows);
    Rcpp::DataFrame map_df(map_cols);

    return Rcpp::List::create(
        Rcpp::Named("relabeled_table") = relabeled_df,
        Rcpp::Named("relabel_map") = map_df
    );
}

// [[Rcpp::export]]
Rcpp::DataFrame relabeledLabelMarginals(Rcpp::DataFrame relabeled_table, Rcpp::List allelic_tips) {
    Rcpp::CharacterVector col_names = relabeled_table.names();
    if (col_names.size() < 2) {
        stop("relabeled_table must contain Iteration plus at least one locus-copy column.");
    }
    if (Rcpp::as<std::string>(col_names[0]) != "Iteration") {
        stop("relabeled_table first column must be named 'Iteration'.");
    }
    const int n_rows = relabeled_table.nrows();
    if (n_rows <= 0) {
        stop("relabeled_table has no rows.");
    }

    ParsedAllelicTips parsed = parse_allelic_tips_with_names(allelic_tips);
    const std::vector<std::vector<std::string>>& groups = parsed.groups;
    const std::vector<std::string>& group_names = parsed.group_names;

    std::vector<std::string> labels;
    std::unordered_set<std::string> seen;
    for (const auto& g : groups) {
        for (const auto& lab : g) {
            if (!seen.insert(lab).second) {
                stop("Label appears more than once in allelic_tips: " + lab);
            }
            labels.push_back(lab);
        }
    }
    if (labels.empty()) {
        stop("allelic_tips did not provide any labels.");
    }

    std::unordered_map<std::string, int> label_idx;
    label_idx.reserve(labels.size());
    for (std::size_t i = 0; i < labels.size(); ++i) {
        label_idx[labels[i]] = static_cast<int>(i);
    }

    std::unordered_set<std::string> labels_in_table;
    for (R_xlen_t j = 1; j < col_names.size(); ++j) {
        Rcpp::CharacterVector col = relabeled_table[j];
        for (int r = 0; r < n_rows; ++r) {
            if (col[r] == NA_STRING) continue;
            labels_in_table.insert(Rcpp::as<std::string>(col[r]));
        }
    }
    if (labels_in_table != seen) {
        stop("Labels in relabeled_table do not match allelic_tips labels.");
    }

    const int n_samples = static_cast<int>(col_names.size() - 1);
    const int n_labels = static_cast<int>(labels.size());
    const int n_groups = static_cast<int>(groups.size());

    std::vector<std::vector<double>> label_probs((std::size_t)n_labels, std::vector<double>((std::size_t)n_samples, NA_REAL));
    std::vector<std::vector<double>> group_probs((std::size_t)n_groups, std::vector<double>((std::size_t)n_samples, NA_REAL));
    std::vector<std::string> sample_names;
    sample_names.reserve((std::size_t)n_samples);

    for (int s = 0; s < n_samples; ++s) {
        Rcpp::CharacterVector col = relabeled_table[s + 1];
        sample_names.push_back(Rcpp::as<std::string>(col_names[s + 1]));

        std::vector<int> counts((std::size_t)n_labels, 0);
        int observed = 0;
        for (int r = 0; r < n_rows; ++r) {
            if (col[r] == NA_STRING) continue;
            std::string v = Rcpp::as<std::string>(col[r]);
            auto it = label_idx.find(v);
            if (it == label_idx.end()) continue;
            counts[(std::size_t)it->second] += 1;
            observed += 1;
        }
        if (observed == 0) continue;

        for (int l = 0; l < n_labels; ++l) {
            label_probs[(std::size_t)l][(std::size_t)s] =
                static_cast<double>(counts[(std::size_t)l]) / static_cast<double>(observed);
        }
        for (int g = 0; g < n_groups; ++g) {
            int sum = 0;
            for (const auto& lab : groups[(std::size_t)g]) {
                sum += counts[(std::size_t)label_idx.at(lab)];
            }
            group_probs[(std::size_t)g][(std::size_t)s] =
                static_cast<double>(sum) / static_cast<double>(observed);
        }
    }

    Rcpp::List out_cols;
    Rcpp::CharacterVector out_names;
    out_cols.push_back(Rcpp::CharacterVector(sample_names.begin(), sample_names.end()));
    out_names.push_back("locus_copy");

    for (int l = 0; l < n_labels; ++l) {
        Rcpp::NumericVector v(n_samples);
        for (int s = 0; s < n_samples; ++s) v[s] = label_probs[(std::size_t)l][(std::size_t)s];
        out_cols.push_back(v);
        out_names.push_back(labels[(std::size_t)l]);
    }

    for (int g = 0; g < n_groups; ++g) {
        Rcpp::NumericVector v(n_samples);
        for (int s = 0; s < n_samples; ++s) v[s] = group_probs[(std::size_t)g][(std::size_t)s];
        out_cols.push_back(v);
        out_names.push_back(group_names[(std::size_t)g]);
    }

    out_cols.attr("names") = out_names;
    out_cols.attr("class") = "data.frame";
    out_cols.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, -n_samples);
    return Rcpp::DataFrame(out_cols);
}

// [[Rcpp::export]]
void relabelPhaseAndTrees(Rcpp::CharacterVector phase_files,
                          Rcpp::List allelic_tips,
                          std::string tree_file = "scratch/22B02_r1/22B02_r1trees.trees",
                          std::string output_prefix = "scratch/22B02_r1/22B02_r1",
                          std::string output_partition_suffix = "_partition_relabeled.log",
                          std::string output_tree_suffix = "_trees_relabeled.trees",
                          std::string output_relabel_map_suffix = "_relabel_map.tsv",
                          std::string output_locus_manifest_suffix = "_locus_manifest.tsv",
                          std::string output_marginal_suffix = "_label_marginals.tsv",
                          double burnin = 0.25) {
                            
    if (tree_file.empty()) stop("tree_file must be non-empty.");
    if (output_prefix.empty()) stop("output_prefix must be non-empty.");
    if (output_partition_suffix.empty()) stop("output_partition_suffix must be non-empty.");
    if (output_tree_suffix.empty()) stop("output_tree_suffix must be non-empty.");
    if (output_relabel_map_suffix.empty()) stop("output_relabel_map_suffix must be non-empty.");
    if (output_locus_manifest_suffix.empty()) stop("output_locus_manifest_suffix must be non-empty.");
    if (output_marginal_suffix.empty()) stop("output_marginal_suffix must be non-empty.");
    if (!(burnin >= 0.0 && burnin < 1.0)) stop("burnin must be in [0, 1).");

    const std::string output_partition_file = output_prefix + output_partition_suffix;
    const std::string output_tree_file = output_prefix + output_tree_suffix;
    const std::string output_relabel_map_file = output_prefix + output_relabel_map_suffix;
    const std::string output_locus_manifest_file = output_prefix + output_locus_manifest_suffix;
    const std::string output_marginal_file = output_prefix + output_marginal_suffix;

    if (phase_files.size() == 0) {
        stop("phase_files must contain at least one file path.");
    }
    Rcpp::CharacterVector nms = phase_files.names();
    if (nms.size() != phase_files.size()) {
        stop("phase_files must be a named character vector.");
    }

    // Convert phase_files named vector into STL tags/paths.
    std::vector<std::string> tags;
    std::vector<std::string> paths;
    tags.reserve((std::size_t)phase_files.size());
    paths.reserve((std::size_t)phase_files.size());
    for (R_xlen_t i = 0; i < phase_files.size(); ++i) {
        if (nms[i] == NA_STRING) stop("phase_files has an NA name at position %d.", (int)i + 1);
        if (phase_files[i] == NA_STRING) stop("phase_files has an NA path at position %d.", (int)i + 1);
        std::string tag = as<std::string>(nms[i]);
        std::string path = as<std::string>(phase_files[i]);
        if (tag.empty()) stop("phase_files has an empty name at position %d.", (int)i + 1);
        if (path.empty()) stop("phase_files has an empty path at position %d.", (int)i + 1);
        tags.push_back(tag);
        paths.push_back(path);
    }

    // Build and align the partition table in STL.
    ParsedAllelicTips parsed_tips = parse_allelic_tips_with_names(allelic_tips);
    std::vector<std::vector<std::string>> sets = parsed_tips.groups;
    RelabeledPhaseTable raw_table = relabel_phase_logs_to_wide_table(tags, paths);
    raw_table = apply_burnin_fraction(raw_table, burnin);
    LabelSwapResult aligned = relabel_wide_table_by_allelic_sets(raw_table, sets);

    // Write aligned partition table.
    {
        std::ofstream out_part(output_partition_file);
        if (!out_part) {
            stop("Failed to open output partition file for writing.");
        }

        out_part << "Iteration";
        for (const auto& nm : aligned.relabeled_table.column_names) out_part << "\t" << nm;
        out_part << "\n";

        const std::size_t n_rows = aligned.relabeled_table.iterations.size();
        const std::size_t n_cols = aligned.relabeled_table.column_names.size();
        for (std::size_t r = 0; r < n_rows; ++r) {
            out_part << aligned.relabeled_table.iterations[r];
            for (std::size_t c = 0; c < n_cols; ++c) {
                out_part << "\t" << aligned.relabeled_table.columns[c][r];
            }
            out_part << "\n";
        }
    }

    // Index iteration -> relabel map for fast tree-string relabeling.
    std::unordered_map<std::int64_t, std::vector<int>> map_by_iter;
    map_by_iter.reserve(aligned.relabeled_table.iterations.size());
    for (std::size_t i = 0; i < aligned.relabeled_table.iterations.size(); ++i) {
        map_by_iter[aligned.relabeled_table.iterations[i]] = aligned.relabel_map_index[i];
    }

    // Write iteration-wise relabel map table.
    {
        std::ofstream out_map(output_relabel_map_file);
        if (!out_map) {
            stop("Failed to open output relabel map file for writing.");
        }

        out_map << "Iteration";
        for (const auto& lab : aligned.labels) out_map << "\t" << lab;
        out_map << "\n";

        for (std::size_t r = 0; r < aligned.relabeled_table.iterations.size(); ++r) {
            out_map << aligned.relabeled_table.iterations[r];
            const std::vector<int>& idx_map = aligned.relabel_map_index[r];
            for (std::size_t j = 0; j < aligned.labels.size(); ++j) {
                out_map << "\t" << aligned.labels[(std::size_t)idx_map[j]];
            }
            out_map << "\n";
        }
    }

    write_label_marginals_file(aligned.relabeled_table, parsed_tips, output_marginal_file);
    write_locus_manifest_file(aligned.relabeled_table, output_locus_manifest_file);

    // Stream tree TSV and relabel Newick strings by token replacement (no Newick parsing).
    std::ifstream in_tree(tree_file);
    if (!in_tree) stop("Failed to open tree_file.");
    std::ofstream out_tree(output_tree_file);
    if (!out_tree) stop("Failed to open output_tree_file for writing.");

    std::string line;
    if (!std::getline(in_tree, line)) stop("tree_file is empty.");
    std::vector<std::string> header = split_tsv_line(line);
    const std::size_t iter_col = find_column_by_name(header, {"Iteration", "iter", "iteration"});

    std::vector<std::vector<std::string>> sample_rows;
    std::vector<std::vector<std::string>> buffered_rows;
    while (sample_rows.size() < 5 && std::getline(in_tree, line)) {
        if (line.empty()) continue;
        std::vector<std::string> row = split_tsv_line(line);
        buffered_rows.push_back(row);
        sample_rows.push_back(row);
    }
    std::size_t newick_col = 0;
    if (!sample_rows.empty()) {
        newick_col = detect_newick_column(header, sample_rows);
    } else {
        newick_col = find_column_by_name(header, {"psi", "Tree", "tree", "newick"});
    }

    out_tree << join_tsv_line(header) << "\n";

    auto process_row = [&](std::vector<std::string> row) {
        if (row.empty()) return;
        if (iter_col >= row.size() || newick_col >= row.size()) {
            out_tree << join_tsv_line(row) << "\n";
            return;
        }
        const std::int64_t iter = parse_int64(row[iter_col]);
        auto it_map = map_by_iter.find(iter);
        if (it_map == map_by_iter.end()) {
            // Burn-in rows (or other unmatched rows) are dropped from output.
            return;
        }

        std::string nwk = row[newick_col];
        const std::vector<int>& idx_map = it_map->second;
        std::vector<std::string> placeholders(aligned.labels.size());
        for (std::size_t i = 0; i < aligned.labels.size(); ++i) {
            placeholders[i] = "__TMP_LABEL_" + std::to_string(i) + "__";
            replace_tip_label_token(nwk, aligned.labels[i], placeholders[i]);
        }
        for (std::size_t i = 0; i < aligned.labels.size(); ++i) {
            const int mapped_idx = idx_map[i];
            replace_tip_label_token(nwk, placeholders[i], aligned.labels[(std::size_t)mapped_idx]);
        }

        row[newick_col] = nwk;
        out_tree << join_tsv_line(row) << "\n";
    };

    for (auto& row : buffered_rows) process_row(row);
    while (std::getline(in_tree, line)) {
        if (line.empty()) continue;
        process_row(split_tsv_line(line));
    }

}
