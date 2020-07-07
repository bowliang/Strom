#include <iostream> 
#include <string>
#include "tree_manip.hpp"
#include "xstrom.hpp"

using namespace strom;

const double Node::_smallest_edge_length = 1.0e-12;

int main(int argc, const char * argv[]) {
    std::cout << "Starting..." << std::endl;
    TreeManip tm;
    std::string newick = std::string("(1:1,(4:0.676357614,(5:0.4311142368,(2:0.4203683334,3:0.4203683334):0.01074590333):0.2452433772):0.323642386);");
    std::cout << "Input: " << newick << std::endl;
    try {
        tm.buildFromNewick(newick, true, false);
        std::cout << "Output: " << tm.makeNewick(3) << std::endl;
        tm.rerootAtNodeNumber(4);                                   
        std::cout << "Output: " << tm.makeNewick(3) << std::endl;   
    }
    catch (XStrom x) {
        std::cout << "Error: " << x.what() << std::endl;
    }
    std::cout << "\nFinished!" << std::endl;

    return 0;
}   