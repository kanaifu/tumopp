/*! @file tree.cpp
    @brief Implementation of Tree class
*/
#include "tree.hpp"

namespace tumopp {

void Tree::add_child(unsigned id_, unsigned id_child_, double dist) {
    if (tree.size() <= id_) {
        // Extend tree until it gets to id
        std::vector<std::vector<Edge>> extension(id_ - tree.size() + 1);
        tree.insert(tree.end(), extension.begin(), extension.end());
    }
    tree[id_].push_back({id_child_, dist});
}

std::string Tree::get_subtree_(unsigned id_, double par_dist) const {
    std::string newick = "(";
    if (id_ < tree.size()) {
        // Node might have children.
        for (auto &[child_id_, dist] : tree[id_]) {
            newick += get_subtree_(child_id_, dist);
        }
    }
    if (newick.size() > 1) {
        newick.push_back(')');
    } else {
        newick = "";
    }
    newick += std::to_string(id_);
    if (par_dist != 0.0) {
        newick += ":" + std::to_string(par_dist);
    }
    return newick;
}

}