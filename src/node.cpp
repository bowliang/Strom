#include "node.hpp"

namespace strom
{
    Node::Node()
    {
        clear();
    }

    Node::~Node()
    {
    }

    void Node::clear()
    {
        _left_child = 0;
        _right_sib = 0;
        _parent = 0;
        _number = -1;
        _name = "";
        _edge_length = _smallest_edge_length;
        _added_t = false;
    }

    void Node::setEdgeLength(double v)
    {
        _edge_length = (v < _smallest_edge_length ? _smallest_edge_length : v);
    }

    Node *Node::getRightChild()
    {
        if (_left_child != NULL)
            return _left_child->getRightSib();
        else
            return NULL;
    }

    void Node::setName(std::string name)
    {
        if (!_added_t) {
            _name = name;
            _added_t = true;
        }
    }
} // namespace strom