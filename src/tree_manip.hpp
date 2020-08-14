#include <cassert>
#include <memory>
#include <stack>
#include <queue>
#include <set>
#include <map>
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
        TreeManip(Tree t);
        TreeManip(const TreeManip &tm1);
        ~TreeManip();

        void setTree(Tree t);
        Tree getTree();

        double calcTreeLength() const;
        unsigned countEdges() const;
        void scaleAllEdgeLengths(double scaler);

        void createTestTree();
        std::string makeNewick(unsigned precision, bool use_names = false) const;

        void buildFromNewick(const std::string newick, bool rooted, bool allow_polytomies);
        void buildFromNewickFile(const std::string filename, bool rooted = true, bool allow_polytomies = false);
        void rerootAtNodeNumber(int node_number);

        void buildNodeNameAndNumberMap();
        int getNodeNumberByName(std::string name);
        void addTToName();

        void narrowExchangeFrom(TreeManip &original_tm);
        void wideExchangeFrom(TreeManip &original_tm);
        void FNPRExchangeFrom(TreeManip &original_tm);
        void totalLengthChange(double &proposed_rootage, double &rate_prior_ratio, double &rate_proposal_ratio, double mean, double stddev);
        void randomLengthChange(double delta_time, double &time_proposal_ratio);
        void allInternalLengthChange(double delta_time, double &time_proposal_ratio);

        void buildNodesPossibilitiesInfo(double m10, double m01);
        Node::PtrVector getPreOrder() { return _tree.getPreOrder(); };

        void clear();

        int getNumLeaves() { return _tree.numLeaves(); };
        int getNumNodes() { return _tree.numNodes(); };
        Node *getRootNode() { return _tree.getRoot(); };
        void updateNodesHeightInfo() { _tree.updateNodesHeightInfo(); };
        double getTreeMaxHeight() { return _tree.getTreeMaxHeight(); };
        void updateSubstitutionMatrix(double br_len, double m10, double m01);

    private:
        Tree _tree;
        std::map<std::string, int> _node_name_and_number_map;
        double _lambda_edge = 2 * log(1.2);
        double _substitution_matrix[2][2];

        Node *findNextPreorder(Node *nd);
        void refreshPreorder();
        void addNodeToPreOrder(Node *nd);
        void refreshLevelorder();
        void renumberInternals();
        void rerootAtNode(Node *prospective_root);
        void extractNodeNumberFromName(Node *nd, std::set<unsigned> &used);
        void extractEdgeLen(Node *nd, std::string edge_length_string);
        unsigned countNewickLeaves(std::string newick);
        void stripOutNexusComments(std::string &newick);
        bool canHaveSibling(Node *nd, bool rooted, bool allow_polytomies);
        void exchangeNodes(Node *child, Node *uncle, Node *parent, Node *grand_parent);
        void replace(Node* parent, Node* child, Node* replacement);
        void pruneChild(Node* parent, Node* child);
        void graftChild(Node* parent, Node* child);
        void changeInternalNode(Node* internal, double delta_time, double &time_proposal_ratio);
    };

} // namespace strom
