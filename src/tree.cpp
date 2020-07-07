#include "tree.hpp"

namespace strom
{
    Tree::Tree()
    {
        std::cout << "Constructing a Tree" << std::endl;
        clear();
    }

    Tree::~Tree()
    {
        std::cout << "Destroying a Tree" << std::endl;
    }

    void Tree::clear()
    {
        _is_rooted = false;
        _root = 0;
        _nodes.clear();
        _preorder.clear();
        _levelorder.clear();
    }

    bool Tree::isRooted() const
    {
        return _is_rooted;
    }

    unsigned Tree::numLeaves() const
    {
        return _nleaves;
    }

    unsigned Tree::numInternals() const
    {
        return _ninternals;
    }

    unsigned Tree::numNodes() const
    {
        return (unsigned)_nodes.size();
    }

} // namespace strom