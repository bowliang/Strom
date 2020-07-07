#include <memory>
#include <iostream>
#include "node.hpp"

namespace strom
{
    class TreeManip;

    class Tree
    {

        friend class TreeManip;

    public:
        Tree();
        ~Tree();

        bool isRooted() const;
        unsigned numLeaves() const;
        unsigned numInternals() const;
        unsigned numNodes() const;

    private:
        void clear();

        bool _is_rooted;
        Node *_root;
        unsigned _nleaves;
        unsigned _ninternals;
        Node::PtrVector _preorder;
        Node::PtrVector _levelorder;
        Node::Vector _nodes;

    public:
        typedef std::shared_ptr<Tree> SharedPtr;
    };

} // namespace strom
