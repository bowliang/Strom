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

        bool isComputedState() { return _computed_state; }
        void setComputed(bool is_computed) { _computed_state = is_computed; }
        double getL0() { return _l0; }
        double getL1() { return _l1; }
        void setL0(double l0) { _l0 = l0; }
        void setL1(double l1) { _l1 = l1; }
        void setState0(int state_0) { _state_0 = state_0;}
        void setState1(int state_1) { _state_1 = state_1;}
        int getState0() { return _state_0;}
        int getState1() { return _state_1;}
        void setObservedState(int observed_state) { _observed_state = observed_state;}
        int getObservedState() { return _observed_state;}
        void setFinalState(int final_state) { _final_state = final_state;}
        int getFinalState() { return _final_state;}

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
        int _observed_state;
        int _final_state;
        int _state_0;
        int _state_1;
        double _l0;
        double _l1;
        bool _computed_state;
        //Split               _split;
    };

} // namespace strom
