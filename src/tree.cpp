#include "tree.hpp"

namespace strom
{
    Tree::Tree()
    {
        //std::cout << "Constructing a Tree" << std::endl;
        clear();
    }

    Tree::~Tree()
    {
        //std::cout << "Destroying a Tree" << std::endl;
    }

    void Tree::clear()
    {
        _is_rooted = false;
        _is_height_built = false;
        _root = 0;
        _nodes.clear();
        _preorder.clear();
        _levelorder.clear();
    }

    Node::PtrVector Tree::getPreOrder()
    {
        return _preorder;
    }

    Node::PtrVector Tree::getAllInternals()
    {
        if (_all_internals.empty())
        {
            Node *real_root = _root->_left_child;
            _all_internals.reserve(_ninternals);
            for (auto nd : _preorder)
            {
                if (nd != real_root && nd->_left_child != NULL)
                {
                    _all_internals.push_back(nd);
                }
            }
        }
        return _all_internals;
    }

    void Tree::buildNodesHeightInfo()
    {
        if (_is_height_built)
            return;
        updateNodesHeightInfo();
        _is_height_built = true;
    }

    void Tree::updateNodesHeightInfo()
    {
        for (auto nd : _preorder)
        {
            if (nd->_parent == _root) {
                nd->_height = 0.0;
            } else {
                nd->_height = nd->_parent->_height + nd->_edge_length;
            }
        }
    }

    double Tree::getTreeMaxHeight()
    {
        double max_height = 0.0;
        for (auto nd : _preorder)
        {
            max_height = std::max(max_height, nd->_height);
        }
        return max_height;
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