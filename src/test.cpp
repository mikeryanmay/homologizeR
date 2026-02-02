#include <Rcpp.h>
#include "trees.hpp"
#include "phases.hpp"

using namespace Rcpp;
using namespace homologizer;

typedef std::map<std::int64_t, Tree> treeMap;
typedef std::map<std::int64_t, std::map<std::string, std::vector<std::string>>> phaseMap;

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
