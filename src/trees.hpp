#ifndef HOMOLOGIZER_TREES_HPP
#define HOMOLOGIZER_TREES_HPP

#include "utils.hpp"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace homologizer {

class Tree; // forward

class Node {
public:
    const std::string& name() const { return name_; }
    double branch_length() const { return branch_length_; }   // parent -> this
    Node* parent() const { return parent_; }
    const std::vector<Node*>& children() const { return children_; }

    bool is_leaf() const { return children_.empty(); }
    std::size_t degree() const { return children_.size(); }

private:
    friend class Tree;
    friend class NewickParser;

    std::string name_;
    double branch_length_ = 0.0;
    Node* parent_ = nullptr;
    std::vector<Node*> children_;
};

class Tree {
public:
    Node* root() const { return root_; }

    // Note: this counts allocated nodes (including detached ones after edits).
    // Use your own traversal if you want "active" node count.
    std::size_t allocated_nodes() const { return storage_.size(); }

    // Apply phasing map:
    // old_tip_label -> list of replacement tip labels (possibly empty).
    void apply_phasing(const std::map<std::string, std::vector<std::string>>& phasing) {
        for (const auto& kv : phasing) {
            const std::string& old_label = kv.first;
            const std::vector<std::string>& new_labels = kv.second;

            Node* tip = find_tip_by_name(old_label);
            if (!tip) {
                throw std::runtime_error(
                        "Phasing refers to tip '" + old_label + "' but no such tip exists in the tree."
                );
            }

            if (new_labels.empty()) {
                delete_leaf_and_cleanup(tip);
                continue;
            }

            if (new_labels.size() == 1) {
                // rename
                tip->name_ = new_labels[0];
                continue;
            }

            // expand: reuse existing tip as attachment/internal node
            tip->name_.clear();
            tip->children_.clear();

            for (const auto& nm : new_labels) {
                Node* child = make_node();
                child->name_ = nm;
                child->branch_length_ = 0.0;   // new tips get length 0 unless you want otherwise
                child->parent_ = tip;
                tip->children_.push_back(child);
            }
        }

        // Optional: if root becomes unary after deletions, clean it
        suppress_unary_root();
    }

private:
    friend class NewickParser;

    // ---- storage / allocation ----
    Node* make_node() {
        storage_.emplace_back(std::make_unique<Node>());
        return storage_.back().get();
    }

    std::vector<std::unique_ptr<Node>> storage_;
    Node* root_ = nullptr;

    // ---- editing helpers (all private, so they can touch Node internals) ----

    Node* find_tip_by_name(const std::string& name) {
        for (const auto& up : storage_) {
            Node* n = up.get();
            if (!n) continue;
            if (n->children_.empty() && n->name_ == name) return n;
        }
        return nullptr;
    }

    static void remove_child_link(Node* parent, Node* child) {
        auto& kids = parent->children_;
        kids.erase(std::remove(kids.begin(), kids.end(), child), kids.end());
    }

    // Suppress unary internal nodes upward, starting at `start`.
    // Preserve path length: child's branch += removed node's branch.
    void suppress_unary_upward(Node* start) {
        Node* cur = start;

        while (cur) {
            if (cur == root_) break; // handle root separately

            if (cur->children_.size() == 1) {
                Node* child  = cur->children_[0];
                Node* parent = cur->parent_;
                if (!parent) break;

                // preserve path length: parent->cur + cur->child becomes parent->child
                child->branch_length_ += cur->branch_length_;

                // replace cur with child in parent's children list
                for (auto& k : parent->children_) {
                    if (k == cur) { k = child; break; }
                }
                child->parent_ = parent;

                // detach cur
                cur->children_.clear();
                cur->parent_ = nullptr;

                // continue upwards: parent may now be unary
                cur = parent;
                continue;
            }

            // stop if not unary
            break;
        }

        // root cleanup handled separately
        suppress_unary_root();
    }

    void suppress_unary_root() {
        // If root has exactly one child, splice it out.
        while (root_ && root_->children_.size() == 1) {
            Node* child = root_->children_[0];
            child->parent_ = nullptr;

            // Root has no incoming edge in our representation.
            // Best convention: set branch length at new root to 0.
            child->branch_length_ = 0.0;

            // detach old root
            root_->children_.clear();
            root_->parent_ = nullptr;

            root_ = child;
        }
    }

    void delete_leaf_and_cleanup(Node* leaf) {
        if (!leaf) return;
        if (!leaf->children_.empty()) {
            throw std::runtime_error("delete_leaf_and_cleanup called on non-leaf node");
        }

        Node* parent = leaf->parent_;
        if (!parent) {
            // leaf is root: tree becomes empty
            root_ = nullptr;
            return;
        }

        remove_child_link(parent, leaf);

        // detach leaf
        leaf->parent_ = nullptr;

        // If parent became unary, suppress upward. If parent became degree 0 (rare), suppress upward won’t help,
        // but that would mean the tree is collapsing toward emptiness; we handle root cleanup anyway.
        suppress_unary_upward(parent);
    }
};

// -------------------- Newick parser --------------------
// Ignores [ ... ] blocks; supports polytomies.
class NewickParser {
public:
    explicit NewickParser(std::string_view input) : s_(input) {}

    Tree parse_tree() {
        Tree t;
        skip_ws_and_brackets();
        Node* r = parse_subtree(t, nullptr);

        skip_ws_and_brackets();
        if (!consume_char(';')) throw error("Expected ';' at end of Newick string");
        skip_ws_and_brackets();
        if (!eof()) throw error("Unexpected trailing characters after ';'");

        t.root_ = r;
        return t;
    }

private:
    std::string_view s_;
    std::size_t i_ = 0;

    Node* parse_subtree(Tree& t, Node* parent) {
        skip_ws_and_brackets();

        Node* n = t.make_node();
        n->parent_ = parent;
        if (parent) parent->children_.push_back(n);

        if (peek_raw() == '(') {
            ++i_; // '('
            skip_ws_and_brackets();

            // first child
            parse_subtree(t, n);
            skip_ws_and_brackets();

            // remaining children
            while (consume_char(',')) {
                skip_ws_and_brackets();
                parse_subtree(t, n);
                skip_ws_and_brackets();
            }

            if (!consume_char(')')) throw error("Expected ')'");

            skip_ws_and_brackets();
            parse_optional_label_and_length(*n);
        } else {
            parse_optional_label_and_length(*n);
        }

        return n;
    }

    void parse_optional_label_and_length(Node& n) {
        skip_ws_and_brackets();

        std::string label = parse_label();
        if (!label.empty()) n.name_ = std::move(label);

        skip_ws_and_brackets();
        if (consume_char(':')) {
            skip_ws_and_brackets();
            n.branch_length_ = parse_number();
            skip_ws_and_brackets();
        }
        skip_ws_and_brackets();
    }

    std::string parse_label() {
        skip_ws_and_brackets();
        std::string out;
        while (!eof()) {
            char c = s_[i_];
            if (c == ':' || c == ',' || c == ')' || c == ';' ||
                c == '[' || std::isspace((unsigned char)c)) {
                break;
            }
            out.push_back(c);
            ++i_;
        }
        return out;
    }

    double parse_number() {
        if (eof()) throw error("Expected number, got end-of-input");

        std::size_t start = i_;
        while (!eof()) {
            char c = s_[i_];
            if (c == ',' || c == ')' || c == ';' || c == '[' ||
                std::isspace((unsigned char)c)) {
                break;
            }
            ++i_;
        }

        std::string token(std::string(s_.substr(start, i_ - start)));
        if (token.empty()) throw error("Expected number after ':'");

        char* endptr = nullptr;
        double val = std::strtod(token.c_str(), &endptr);
        if (endptr == token.c_str() || *endptr != '\0') {
            throw error("Invalid number token: '" + token + "'");
        }
        return val;
    }

    void skip_ws_and_brackets() {
        for (;;) {
            while (!eof() && std::isspace((unsigned char)s_[i_])) ++i_;
            if (!eof() && s_[i_] == '[') {
                ++i_;
                while (!eof() && s_[i_] != ']') ++i_;
                if (eof()) throw error("Unterminated bracket block '[' ...");
                ++i_;
                continue;
            }
            break;
        }
    }

    bool consume_char(char expected) {
        skip_ws_and_brackets();
        if (!eof() && s_[i_] == expected) {
            ++i_;
            return true;
        }
        return false;
    }

    char peek_raw() const { return eof() ? '\0' : s_[i_]; }
    bool eof() const { return i_ >= s_.size(); }

    std::runtime_error error(const std::string& msg) const {
        return std::runtime_error("Newick parse error at pos " + std::to_string(i_) + ": " + msg);
    }
};

// -------------------- TSV reading for trees --------------------

inline std::string strip_brackets(const std::string& s) {
    std::string out;
    bool in_bracket = false;
    for (char c : s) {
        if (!in_bracket) {
            if (c == '[') in_bracket = true;
            else out.push_back(c);
        } else {
            if (c == ']') in_bracket = false;
        }
    }
    return out;
}

inline bool parens_balanced(const std::string& s) {
    int bal = 0;
    for (char c : s) {
        if (c == '(') ++bal;
        else if (c == ')') {
            --bal;
            if (bal < 0) return false;
        }
    }
    return bal == 0;
}

inline int newick_score(const std::string& field_raw) {
    std::string f = trim(field_raw);
    if (f.empty()) return 0;

    std::string g = strip_brackets(f);

    bool has_parens = g.find('(') != std::string::npos && g.find(')') != std::string::npos;
    bool has_semicolon = g.find(';') != std::string::npos;
    bool ends_semicolon = !g.empty() && g.back() == ';';
    bool balanced = parens_balanced(g);
    bool has_colon = g.find(':') != std::string::npos;

    int score = 0;
    if (has_parens) score += 3;
    if (has_semicolon) score += 3;
    if (ends_semicolon) score += 2;
    if (balanced) score += 3;
    if (has_colon) score += 1;

    return score;
}

inline std::size_t detect_newick_column(
        const std::vector<std::string>& header,
        const std::vector<std::vector<std::string>>& sample_rows
) {
    if (header.empty()) throw std::runtime_error("Empty header");

    std::vector<int> scores(header.size(), 0);
    for (std::size_t j = 0; j < header.size(); ++j) {
        for (const auto& row : sample_rows) {
            if (j < row.size()) scores[j] += newick_score(row[j]);
        }
    }

    auto it = std::max_element(scores.begin(), scores.end());
    if (it == scores.end() || *it == 0) {
        throw std::runtime_error("Could not detect Newick column (no column scored > 0)");
    }
    return (std::size_t)std::distance(scores.begin(), it);
}

// Return: iteration -> Tree
inline std::map<std::int64_t, Tree>
    read_tree_tsv_by_iteration(const std::string& filepath,
                               std::size_t sample_rows_to_detect = 5) {
    std::ifstream in(filepath);
    if (!in) throw std::runtime_error("Failed to open file: " + filepath);

    std::string line;

    if (!std::getline(in, line))
        throw std::runtime_error("File is empty: " + filepath);

    auto header = split_tsv_line(line);
    std::size_t iter_col = find_column_by_name(header, {"Iteration", "iter", "iteration"});

    std::vector<std::vector<std::string>> samples;
    std::vector<std::vector<std::string>> buffered_rows;

    while (samples.size() < sample_rows_to_detect && std::getline(in, line)) {
        if (line.empty()) continue;
        auto row = split_tsv_line(line);
        buffered_rows.push_back(row);
        samples.push_back(row);
    }

    std::size_t newick_col = detect_newick_column(header, samples);

    std::map<std::int64_t, Tree> trees;

    auto parse_row = [&](const std::vector<std::string>& row) {
        if (iter_col >= row.size() || newick_col >= row.size()) return;

        std::int64_t iter = parse_int64(row[iter_col]);
        std::string nwk = trim(row[newick_col]);
        if (nwk.empty()) return;

        trees[iter] = NewickParser(nwk).parse_tree();
    };

    for (const auto& r : buffered_rows) parse_row(r);

    while (std::getline(in, line)) {
        if (line.empty()) continue;
        auto row = split_tsv_line(line);
        parse_row(row);
    }

    return trees;
}

inline void write_newick_node(std::ostream& os,
                              const Node* n,
                              int precision = 6) {
    if (!n) return;

    // Internal: write children list
    if (!n->children().empty()) {
        os << "(";
        const auto& kids = n->children();
        for (std::size_t i = 0; i < kids.size(); ++i) {
            if (i) os << ",";
            write_newick_node(os, kids[i], precision);
        }
        os << ")";
        // Optional internal-node name (usually empty in your edited trees)
        if (!n->name().empty()) os << n->name();
    } else {
        // Leaf: name required (if empty, Newick leaf is unnamed)
        os << n->name();
    }

    // Branch length: in our representation, this is the edge from parent -> this.
    // Root has no parent edge; by convention we omit root branch length.
    if (n->parent() != nullptr) {
        os << ":" << std::fixed << std::setprecision(precision) << n->branch_length();
    }
}

inline std::string to_newick_string(const Tree& t, int precision = 6) {
    std::ostringstream oss;
    write_newick_node(oss, t.root(), precision);
    oss << ";";
    return oss.str();
}

// --- Write relabeled trees to a TSV file ---
// Map: iteration -> Tree
inline void write_tree_map_tsv(const std::string& outpath,
                               const std::map<std::int64_t, Tree>& trees,
                               int precision = 6,
                               const std::string& iteration_header = "Iteration",
                               const std::string& tree_header = "Tree") {
    std::ofstream out(outpath);
    if (!out) throw std::runtime_error("Failed to open output file for writing: " + outpath);

    // header
    out << iteration_header << "\t" << tree_header << "\n";

    for (const auto& kv : trees) {
        std::int64_t iter = kv.first;
        const Tree& tree = kv.second;
        out << iter << "\t" << to_newick_string(tree, precision) << "\n";
    }
}

} // namespace homologizer

#endif // HOMOLOGIZER_TREES_HPP
