#include "core_set.hpp"
#include "trees.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace homologizer {

namespace {

// A phase group is the set of sample labels currently co-phased together.
using Group = std::vector<std::string>;

struct IterGroups {
    // Groups extracted from one posterior tree draw.
    std::vector<Group> groups;
    // Phased sample labels observed in that tree draw (used for pairwise denominators).
    std::unordered_set<std::string> present_samples;
};

// A phase polytomy is the synthetic node created by apply_phasing:
// internal node with >=2 leaf children, each child branch length exactly 0.
inline bool is_phase_polytomy(const Node* n) {
    if (!n) return false;
    const auto& kids = n->children();
    if (kids.size() < 2) return false;

    for (const Node* child : kids) {
        if (!child) return false;
        if (!child->is_leaf()) return false;
        if (child->name().empty()) return false;
        if (child->branch_length() != 0.0) return false;
    }

    return true;
}

void collect_all_nodes_dfs(const Node* n, std::vector<const Node*>& out) {
    if (!n) return;
    out.push_back(n);
    for (const Node* c : n->children()) collect_all_nodes_dfs(c, out);
}

Group dedup_preserve_order(const Group& in) {
    Group out;
    out.reserve(in.size());
    std::unordered_set<std::string> seen;
    for (const auto& s : in) {
        if (seen.insert(s).second) out.push_back(s);
    }
    return out;
}

bool natural_less(const std::string& a, const std::string& b) {
    std::size_t i = 0;
    std::size_t j = 0;

    while (i < a.size() && j < b.size()) {
        const unsigned char ca = static_cast<unsigned char>(a[i]);
        const unsigned char cb = static_cast<unsigned char>(b[j]);
        const bool da = std::isdigit(ca) != 0;
        const bool db = std::isdigit(cb) != 0;

        if (da && db) {
            std::size_t i0 = i;
            std::size_t j0 = j;
            while (i < a.size() && std::isdigit(static_cast<unsigned char>(a[i])) != 0) ++i;
            while (j < b.size() && std::isdigit(static_cast<unsigned char>(b[j])) != 0) ++j;

            std::size_t ai = i0;
            std::size_t bj = j0;
            while (ai < i && a[ai] == '0') ++ai;
            while (bj < j && b[bj] == '0') ++bj;

            const std::size_t len_a = i - ai;
            const std::size_t len_b = j - bj;
            if (len_a != len_b) return len_a < len_b;

            for (std::size_t k = 0; k < len_a; ++k) {
                const char xa = a[ai + k];
                const char xb = b[bj + k];
                if (xa != xb) return xa < xb;
            }

            const std::size_t full_len_a = i - i0;
            const std::size_t full_len_b = j - j0;
            if (full_len_a != full_len_b) return full_len_a < full_len_b;
            continue;
        }

        if (ca != cb) return ca < cb;
        ++i;
        ++j;
    }

    return a.size() < b.size();
}

struct ParsedPhaseLabel {
    bool ok = false;
    int copy = -1;
    int locus = -1;
};

ParsedPhaseLabel parse_phase_label(const std::string& s) {
    ParsedPhaseLabel out;
    const std::size_t cpos = s.rfind('c');
    const std::size_t lpos = s.rfind("_L");
    if (cpos == std::string::npos || lpos == std::string::npos) return out;
    if (cpos + 1 >= lpos || lpos + 2 >= s.size()) return out;

    for (std::size_t i = cpos + 1; i < lpos; ++i) {
        if (std::isdigit(static_cast<unsigned char>(s[i])) == 0) return out;
    }
    for (std::size_t i = lpos + 2; i < s.size(); ++i) {
        if (std::isdigit(static_cast<unsigned char>(s[i])) == 0) return out;
    }

    out.ok = true;
    out.copy = std::stoi(s.substr(cpos + 1, lpos - (cpos + 1)));
    out.locus = std::stoi(s.substr(lpos + 2));
    return out;
}

bool phase_label_less(const std::string& a, const std::string& b) {
    const ParsedPhaseLabel pa = parse_phase_label(a);
    const ParsedPhaseLabel pb = parse_phase_label(b);
    if (pa.ok && pb.ok) {
        if (pa.locus != pb.locus) return pa.locus < pb.locus; // L first
        if (pa.copy != pb.copy) return pa.copy < pb.copy;     // then c
    } else if (pa.ok != pb.ok) {
        // Keep recognizable phased labels grouped before unmatched labels.
        return pa.ok;
    }
    return natural_less(a, b);
}

IterGroups extract_iteration_groups(const Tree& t, bool aggregate_alleles) {
    IterGroups out;
    const Node* root = t.root();
    if (!root) return out;

    std::vector<const Node*> nodes;
    collect_all_nodes_dfs(root, nodes);

    std::vector<const Node*> phase_nodes;
    phase_nodes.reserve(nodes.size());
    // Collect candidate phase polytomies independent of global topology.
    for (const Node* n : nodes) {
        if (is_phase_polytomy(n)) phase_nodes.push_back(n);
    }

    if (!aggregate_alleles) {
        // Allele-specific mode: keep each phase polytomy as its own group.
        for (const Node* pn : phase_nodes) {
            Group g;
            g.reserve(pn->children().size());
            for (const Node* child : pn->children()) g.push_back(child->name());
            for (const auto& sample : g) out.present_samples.insert(sample);
            out.groups.push_back(std::move(g));
        }
        return out;
    }

    // Population-level mode: merge sister phase polytomies (2-allele representation).
    std::unordered_map<const Node*, std::size_t> idx;
    idx.reserve(phase_nodes.size());
    for (std::size_t i = 0; i < phase_nodes.size(); ++i) idx[phase_nodes[i]] = i;

    std::vector<char> used(phase_nodes.size(), 0);
    for (std::size_t i = 0; i < phase_nodes.size(); ++i) {
        if (used[i]) continue;
        const Node* a = phase_nodes[i];
        const Node* p = a->parent();

        if (!p || p->children().size() != 2) continue;
        const Node* c0 = p->children()[0];
        const Node* c1 = p->children()[1];
        auto it0 = idx.find(c0);
        auto it1 = idx.find(c1);
        if (it0 == idx.end() || it1 == idx.end()) continue;

        const std::size_t j0 = it0->second;
        const std::size_t j1 = it1->second;
        if (used[j0] || used[j1]) continue;

        Group merged;
        merged.reserve(c0->children().size() + c1->children().size());
        for (const Node* child : c0->children()) merged.push_back(child->name());
        for (const Node* child : c1->children()) merged.push_back(child->name());
        Group merged_dedup = dedup_preserve_order(merged);
        for (const auto& sample : merged_dedup) out.present_samples.insert(sample);
        out.groups.push_back(std::move(merged_dedup));
        used[j0] = 1;
        used[j1] = 1;
    }

    // Keep unpaired phase polytomies as standalone groups.
    for (std::size_t i = 0; i < phase_nodes.size(); ++i) {
        if (used[i]) continue;
        Group g;
        g.reserve(phase_nodes[i]->children().size());
        for (const Node* child : phase_nodes[i]->children()) g.push_back(child->name());
        for (const auto& sample : g) out.present_samples.insert(sample);
        out.groups.push_back(std::move(g));
    }

    return out;
}

std::vector<std::vector<int>> threshold_components(
    const std::vector<std::vector<double>>& probs,
    const std::vector<std::vector<std::int64_t>>& co_observed,
    double C
) {
    // Build graph edges for sample pairs with posterior co-phasing >= C.
    const int n = static_cast<int>(probs.size());
    std::vector<std::vector<int>> adj(n);
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (co_observed[i][j] <= 0) continue;
            if (probs[i][j] >= C) {
                adj[i].push_back(j);
                adj[j].push_back(i);
            }
        }
    }

    // Connected components define output core sets.
    std::vector<char> seen(n, 0);
    std::vector<std::vector<int>> comps;
    for (int i = 0; i < n; ++i) {
        if (seen[i]) continue;
        std::queue<int> q;
        q.push(i);
        seen[i] = 1;
        std::vector<int> comp;
        while (!q.empty()) {
            int u = q.front();
            q.pop();
            comp.push_back(u);
            for (int v : adj[u]) {
                if (!seen[v]) {
                    seen[v] = 1;
                    q.push(v);
                }
            }
        }
        comps.push_back(std::move(comp));
    }
    return comps;
}

void enumerate_cliques_from_threshold_graph(
    const std::vector<std::vector<double>>& probs,
    const std::vector<std::vector<std::int64_t>>& co_observed,
    double C,
    std::vector<std::vector<int>>& out_cliques
) {
    const int n = static_cast<int>(probs.size());
    std::vector<std::vector<char>> adj(n, std::vector<char>(n, 0));
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (co_observed[i][j] > 0 && probs[i][j] >= C) {
                adj[i][j] = 1;
                adj[j][i] = 1;
            }
        }
    }

    std::vector<int> current;
    current.reserve(n);

    const auto is_compatible = [&](int v, const std::vector<int>& clique) {
        for (int u : clique) {
            if (!adj[u][v]) return false;
        }
        return true;
    };

    std::function<void(int)> dfs = [&](int start) {
        for (int v = start; v < n; ++v) {
            if (!is_compatible(v, current)) continue;
            current.push_back(v);
            if (current.size() >= 2) out_cliques.push_back(current);
            dfs(v + 1);
            current.pop_back();
        }
    };

    dfs(0);
}

} // namespace

Rcpp::List posterior_core_sets_from_phase_labeled_trees(
    const std::string& tree_file,
    double C,
    bool aggregate_alleles
) {
    
    if (!(C > 0.0 && C <= 1.0)) {
        throw std::runtime_error("C must be in (0, 1].");
    }

    // Read posterior trees and extract per-iteration phase groups.
    const std::map<std::int64_t, Tree> trees = read_tree_tsv_by_iteration(tree_file);
    if (trees.empty()) {
        throw std::runtime_error("No trees were read from input file.");
    }

    std::vector<IterGroups> by_iter;
    by_iter.reserve(trees.size());
    std::unordered_set<std::string> sample_set;
    for (const auto& kv : trees) {
        IterGroups ig = extract_iteration_groups(kv.second, aggregate_alleles);
        for (const auto& s : ig.present_samples) sample_set.insert(s);
        by_iter.push_back(std::move(ig));
    }

    // Stable index over all samples that appear in any iteration.
    std::vector<std::string> samples(sample_set.begin(), sample_set.end());
    std::sort(samples.begin(), samples.end(), phase_label_less);
    const std::size_t n = samples.size();
    std::unordered_map<std::string, std::size_t> sample_idx;
    sample_idx.reserve(n);
    for (std::size_t i = 0; i < n; ++i) sample_idx[samples[i]] = i;

    std::vector<std::vector<std::int64_t>> co_observed(n, std::vector<std::int64_t>(n, 0));
    std::vector<std::vector<std::int64_t>> co_phased(n, std::vector<std::int64_t>(n, 0));
    std::vector<std::unordered_map<std::size_t, int>> sample_group_by_iter;
    sample_group_by_iter.reserve(by_iter.size());

    // Count pairwise co-observation (denominator) and co-phasing (numerator).
    for (const auto& ig : by_iter) {
        std::vector<std::size_t> present;
        present.reserve(ig.present_samples.size());
        for (const auto& s : ig.present_samples) {
            auto it = sample_idx.find(s);
            if (it != sample_idx.end()) present.push_back(it->second);
        }
        std::sort(present.begin(), present.end());
        for (std::size_t a = 0; a < present.size(); ++a) {
            const std::size_t i = present[a];
            co_observed[i][i] += 1;
            for (std::size_t b = a + 1; b < present.size(); ++b) {
                const std::size_t j = present[b];
                co_observed[i][j] += 1;
                co_observed[j][i] += 1;
            }
        }

        // Guard against accidental duplicate assignment in one iteration.
        std::unordered_map<std::size_t, int> assigned_group;
        for (std::size_t g = 0; g < ig.groups.size(); ++g) {
            for (const auto& s : ig.groups[g]) {
                auto it = sample_idx.find(s);
                if (it == sample_idx.end()) continue;
                const std::size_t idx_s = it->second;
                auto ag = assigned_group.find(idx_s);
                if (ag != assigned_group.end() && ag->second != static_cast<int>(g)) {
                    throw std::runtime_error("Sample assigned to multiple phase groups in one tree: " + s);
                }
                assigned_group[idx_s] = static_cast<int>(g);
            }
        }
        sample_group_by_iter.push_back(std::move(assigned_group));

        for (const auto& g : ig.groups) {
            std::vector<std::size_t> members;
            members.reserve(g.size());
            for (const auto& s : g) {
                auto it = sample_idx.find(s);
                if (it != sample_idx.end()) members.push_back(it->second);
            }
            std::sort(members.begin(), members.end());
            members.erase(std::unique(members.begin(), members.end()), members.end());
            for (std::size_t a = 0; a < members.size(); ++a) {
                const std::size_t i = members[a];
                co_phased[i][i] += 1;
                for (std::size_t b = a + 1; b < members.size(); ++b) {
                    const std::size_t j = members[b];
                    co_phased[i][j] += 1;
                    co_phased[j][i] += 1;
                }
            }
        }
    }

    // Convert counts to posterior probabilities; NaN if never co-observed.
    std::vector<std::vector<double>> probs(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            if (i == j) {
                probs[i][j] = 1.0;
                continue;
            }
            if (co_observed[i][j] == 0) {
                probs[i][j] = std::numeric_limits<double>::quiet_NaN();
            } else {
                probs[i][j] = static_cast<double>(co_phased[i][j]) /
                              static_cast<double>(co_observed[i][j]);
            }
        }
    }

    std::vector<std::vector<int>> candidate_cliques;
    enumerate_cliques_from_threshold_graph(probs, co_observed, C, candidate_cliques);

    // Prepare standard-library outputs first; convert to Rcpp types only at return.
    std::vector<double> prob_flat(n * n, 0.0);
    std::vector<int> phased_flat(n * n, 0);
    std::vector<int> observed_flat(n * n, 0);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            const std::size_t k = i + (n * j);
            prob_flat[k] = probs[i][j];
            phased_flat[k] = static_cast<int>(co_phased[i][j]);
            observed_flat[k] = static_cast<int>(co_observed[i][j]);
        }
    }

    std::vector<std::vector<std::string>> core_sets_std;
    core_sets_std.reserve(candidate_cliques.size());
    for (std::size_t k = 0; k < candidate_cliques.size(); ++k) {
        std::vector<std::string> members;
        members.reserve(candidate_cliques[k].size());
        for (int idx_i : candidate_cliques[k]) members.push_back(samples[static_cast<std::size_t>(idx_i)]);
        std::sort(members.begin(), members.end(), phase_label_less);
        core_sets_std.push_back(std::move(members));
    }
    std::sort(core_sets_std.begin(), core_sets_std.end(),
              [](const std::vector<std::string>& a, const std::vector<std::string>& b) {
                  if (a.empty() || b.empty()) return a.size() < b.size();
                  return phase_label_less(a.front(), b.front());
              });

    std::vector<double> core_set_marginal_probability(
        core_sets_std.size(), std::numeric_limits<double>::quiet_NaN()
    );
    for (std::size_t k = 0; k < core_sets_std.size(); ++k) {
        std::int64_t observed = 0;
        std::int64_t together = 0;
        const auto& members = core_sets_std[k];

        for (std::size_t it = 0; it < sample_group_by_iter.size(); ++it) {
            const auto& group_map = sample_group_by_iter[it];
            int group_id = -1;
            bool all_present = true;
            bool same_group = true;

            for (const auto& label : members) {
                auto sit = sample_idx.find(label);
                if (sit == sample_idx.end()) {
                    all_present = false;
                    break;
                }
                auto git = group_map.find(sit->second);
                if (git == group_map.end()) {
                    all_present = false;
                    break;
                }
                if (group_id < 0) {
                    group_id = git->second;
                } else if (git->second != group_id) {
                    same_group = false;
                }
            }

            if (!all_present) continue;
            observed += 1;
            if (same_group) together += 1;
        }

        if (observed > 0) {
            core_set_marginal_probability[k] =
                static_cast<double>(together) / static_cast<double>(observed);
        }
    }

    Rcpp::NumericMatrix prob_mat(static_cast<int>(n), static_cast<int>(n));
    Rcpp::IntegerMatrix phased_mat(static_cast<int>(n), static_cast<int>(n));
    Rcpp::IntegerMatrix observed_mat(static_cast<int>(n), static_cast<int>(n));
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            const std::size_t k = i + (n * j);
            prob_mat(static_cast<int>(i), static_cast<int>(j)) = prob_flat[k];
            phased_mat(static_cast<int>(i), static_cast<int>(j)) = phased_flat[k];
            observed_mat(static_cast<int>(i), static_cast<int>(j)) = observed_flat[k];
        }
    }
    Rcpp::CharacterVector sample_vec(samples.begin(), samples.end());
    Rcpp::List dimnames = Rcpp::List::create(sample_vec, sample_vec);
    prob_mat.attr("dimnames") = dimnames;
    phased_mat.attr("dimnames") = dimnames;
    observed_mat.attr("dimnames") = dimnames;

    struct Candidate {
        std::size_t idx = 0;
        double marginal = std::numeric_limits<double>::quiet_NaN();
    };
    std::vector<Candidate> candidates;
    candidates.reserve(core_sets_std.size());
    for (std::size_t k = 0; k < core_sets_std.size(); ++k) {
        const double p = core_set_marginal_probability[k];
        if (!std::isnan(p) && p >= C) {
            candidates.push_back({k, p});
        }
    }

    // Enforce non-overlapping returned core sets:
    // greedily keep highest-marginal candidates, then larger sets, then deterministic label order.
    std::sort(candidates.begin(), candidates.end(),
              [&](const Candidate& a, const Candidate& b) {
                  if (a.marginal != b.marginal) return a.marginal > b.marginal;
                  const std::size_t asz = core_sets_std[a.idx].size();
                  const std::size_t bsz = core_sets_std[b.idx].size();
                  if (asz != bsz) return asz > bsz;
                  return std::lexicographical_compare(
                      core_sets_std[a.idx].begin(), core_sets_std[a.idx].end(),
                      core_sets_std[b.idx].begin(), core_sets_std[b.idx].end(),
                      phase_label_less
                  );
              });

    std::vector<char> used_sample(n, 0);
    std::vector<std::vector<std::string>> core_sets_filtered;
    std::vector<double> core_set_marginal_probability_filtered;
    core_sets_filtered.reserve(candidates.size());
    core_set_marginal_probability_filtered.reserve(candidates.size());
    for (const auto& cand : candidates) {
        const auto& members = core_sets_std[cand.idx];
        bool overlaps = false;
        for (const auto& s : members) {
            const auto it = sample_idx.find(s);
            if (it != sample_idx.end() && used_sample[it->second]) {
                overlaps = true;
                break;
            }
        }
        if (overlaps) continue;

        for (const auto& s : members) {
            const auto it = sample_idx.find(s);
            if (it != sample_idx.end()) used_sample[it->second] = 1;
        }
        core_sets_filtered.push_back(members);
        core_set_marginal_probability_filtered.push_back(cand.marginal);
    }

    return Rcpp::List::create(
        Rcpp::Named("samples") = Rcpp::wrap(samples),                       // Canonical sample order for matrix rows/columns.
        Rcpp::Named("co_phasing_probability") = prob_mat,                   // Posterior pairwise co-phasing probability matrix P(i, j).
        Rcpp::Named("co_phased_count") = phased_mat,                        // Numerator: iterations where i and j are co-phased.
        Rcpp::Named("co_observed_count") = observed_mat,                    // Denominator: iterations where i and j are both observed.
        Rcpp::Named("core_sets") = Rcpp::wrap(core_sets_filtered),          // Core sets retained under marginal-rule threshold: Pr(set together | observed) >= C.
        Rcpp::Named("core_set_marginal_probability") = Rcpp::wrap(core_set_marginal_probability_filtered), // Marginal probabilities aligned to returned core_sets.
        Rcpp::Named("threshold_C") = C,                                     // Threshold used for pairwise edge/clique candidacy and marginal filtering.
        Rcpp::Named("aggregate_alleles") = aggregate_alleles,               // Whether sister allelic phase tips were collapsed.
        Rcpp::Named("n_iterations") = static_cast<int>(trees.size())        // Number of posterior trees/iterations summarized.
    );
}

} // namespace homologizer
