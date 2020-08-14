#include <string>
#include <vector>
#include <iostream>

namespace strom
{
    class Tree;
    class TreeManip;

    class Node
    {
        friend class Tree;
        friend class TreeManip;
        //friend class Likelihood;
        //friend class Updater;

    public:
        Node();
        ~Node();

        Node *getParent() { return _parent; }
        Node *getLeftChild() { return _left_child; }
        Node *getRightChild();
        Node *getRightSib() { return _right_sib; }
        int getNumber() { return _number; }
        std::string getName() { return _name; }
        void setName(std::string name);
        //Split               getSplit()      {return _split;}

        double getEdgeLength() { return _edge_length; }
        void setEdgeLength(double v);
        double getHeight() { return _height; }

        void setP00(double p00) { _p00 = p00; }
        void setP01(double p01) { _p01 = p01; }
        void setP10(double p10) { _p10 = p10; }
        void setP11(double p11) { _p11 = p11; }
        double getP00() { return _p00; }
        double getP01() { return _p01; }
        double getP10() { return _p10; }
        double getP11() { return _p11; }

        static const double _smallest_edge_length;

        typedef std::vector<Node> Vector;
        typedef std::vector<Node *> PtrVector;

    private:
        void clear();

        Node *_left_child;
        Node *_right_child;
        Node *_right_sib;
        Node *_parent;
        int _number;
        std::string _name;
        double _edge_length;
        double _height;
        double _p00;
        double _p01;
        double _p10;
        double _p11;
        bool _added_t;
        //Split               _split;
    };

} // namespace strom
