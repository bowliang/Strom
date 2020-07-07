#include <cassert>
#include <memory>
#include <stack>
#include <queue>
#include <set>
#include <boost/regex.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/format.hpp>
#include "tree.hpp"
#include "xstrom.hpp"

namespace strom
{

    class TreeManip
    {
    public:
        TreeManip();
        TreeManip(Tree::SharedPtr t);
        ~TreeManip();

        void setTree(Tree::SharedPtr t);
        Tree::SharedPtr getTree();

        double calcTreeLength() const;
        unsigned countEdges() const;
        void scaleAllEdgeLengths(double scaler);

        void createTestTree();
        std::string makeNewick(unsigned precision, bool use_names = false) const;

        void buildFromNewick(const std::string newick, bool rooted, bool allow_polytomies);
        void rerootAtNodeNumber(int node_number);

        void clear();

    private:
        Tree::SharedPtr _tree;

        Node *findNextPreorder(Node *nd);
        void refreshPreorder();
        void refreshLevelorder();
        void renumberInternals();
        void rerootAtNode(Node *prospective_root);
        void extractNodeNumberFromName(Node *nd, std::set<unsigned> &used);
        void extractEdgeLen(Node *nd, std::string edge_length_string);
        unsigned countNewickLeaves(std::string newick);
        void stripOutNexusComments(std::string &newick);
        bool canHaveSibling(Node *nd, bool rooted, bool allow_polytomies);

    public:
        typedef std::shared_ptr<TreeManip> SharedPtr;
    };

} // namespace strom
