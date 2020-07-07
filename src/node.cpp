#include "node.hpp"

namespace strom
{
    Node::Node()
    {
        //std::cout << "Creating Node object" << std::endl;
        clear();
    }

    Node::~Node()
    {
        //std::cout << "Destroying Node object" << std::endl;
    }

    void Node::clear()
    {
        _left_child = 0;
        _right_sib = 0;
        _parent = 0;
        _number = -1;
        _name = "";
        _edge_length = _smallest_edge_length;
    }

    void Node::setEdgeLength(double v)
    {
        _edge_length = (v < _smallest_edge_length ? _smallest_edge_length : v);
    }
} // namespace strom