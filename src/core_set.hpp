#ifndef HOMOLOGIZER_CORE_SET_HPP
#define HOMOLOGIZER_CORE_SET_HPP

#include <Rcpp.h>

#include <string>

namespace homologizer {

Rcpp::List posterior_core_sets_from_phase_labeled_trees(
    const std::string& tree_file,
    double C = 0.5,
    bool aggregate_alleles = true
);

} // namespace homologizer

#endif // HOMOLOGIZER_CORE_SET_HPP
