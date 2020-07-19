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
        Node::PtrVector getPreOrder();
        Node* getRoot() { return _root; };
        Node::PtrVector getAllInternals();
        void buildNodesHeightInfo();
        void updateNodesHeightInfo();
        double getTreeMaxHeight();

    private:
        void clear();

        bool _is_rooted;
        bool _is_height_built;
        Node *_root;
        unsigned _nleaves;
        unsigned _ninternals;
        Node::PtrVector _preorder;
        Node::PtrVector _levelorder;
        Node::PtrVector _all_internals;
        Node::Vector _nodes;

    public:
        typedef std::shared_ptr<Tree> SharedPtr;
    };

} // namespace strom
