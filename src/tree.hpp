/*! @file tree.hpp
    @brief Interface of Tree class.
*/
#pragma once
#ifndef TUMOPP_TREE_HPP_
#define TUMOPP_TREE_HPP_

#include <vector>

namespace tumopp {

struct Edge {
    unsigned id_;
    double dist;
};

class Tree {
public:
    ~Tree() = default;
    void add_child(unsigned id_, unsigned id_child_, double dist);

    /*
        @id_ root id.
        @par_dist dist to this node from parent (0 if no dist).
    */
    std::string get_subtree(unsigned id_, double par_dist) const {
        return get_subtree_(id_, par_dist) + ";";
    }

private:
    std::vector<std::vector<Edge>> tree;

    std::string get_subtree_(unsigned id_, double par_dist) const;
};

} // namespace tumopp

#endif // TUMOPP_TREE_HPP_