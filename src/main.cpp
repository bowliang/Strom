#include <iostream>
#include <string>
#include <random>
#include <vector>
#include <math.h>
#include <unistd.h>  
#include <chrono> 
#include <ctime>

#include "tree_manip.hpp"
#include "utility.hpp"
#include "xstrom.hpp"

using namespace strom;

const double Node::_smallest_edge_length = 1.0e-12;

int main(int argc, char *argv[])
{
    // Get starting timepoint 
    auto start = std::chrono::high_resolution_clock::now(); 

    generator.seed(time(0));

    // 0. read arguments
    // Read the arguments passing from command line and validate
    int opt, niter = 1000000;
    std::string output_filename;
    while ((opt = getopt(argc, argv, "n:f:")) != -1)
    {
        switch (opt)
        {
        case 'n':
            niter = atoi(optarg);
            break;
        case 'f':
            output_filename = optarg;
            std::cout<<"output_filename: " <<output_filename<<"\n";
            break;
        default:
            fprintf(stderr, "Wrong arguments!\n");
            exit(1);
        }
    }

    // 1. Start with input random tree file and input observed data
    std::string true_tree_filename = "./input/tree/RandomTreeScale_3.tre";
    std::string true_time_tree_filename = "./input/tree/RandomTree_3.tre";
    std::string tree_data_filename = "./Long_Obs_binary_obs_0_1_tip_mu_005_alpha_001_beta_001_matrix_2.csv";

    // 2. Build tree and get tree info
    TreeManip true_tree_tm, true_time_tree_tm;
    true_tree_tm.buildFromNewickFile(true_tree_filename);
    true_time_tree_tm.buildFromNewickFile(true_time_tree_filename);

    std::cout << "true_tree_tm: " << true_tree_tm.makeNewick(3) << std::endl;
    std::cout << "true_time_tree_tm: " << true_time_tree_tm.makeNewick(3) << std::endl;

    std::vector<std::string> data_label;
    std::vector<std::vector<int>> tree_data_original;
    readObservedDataFromCSVFile(5000, tree_data_filename, data_label, tree_data_original);

    for (auto label : data_label) 
    {
        std::cout<<"label "<<label<<"\n";
    }

    for (auto data : tree_data_original[0]) 
    {
        std::cout<<"data "<<data<<"\n";
    }

    // 3. initialize value based on some distribution (normal, dirichlet)
    alpha = 0.01; // getUniformDistribution(0, 0.45);
    beta = 0.01; // getUniformDistribution(0, 0.45);
    m01 = 1, m00 = -m01, rootage = 1.045;
    m10 = getUniformDistribution(0, 1.0);
    m11 = -m10;
    pi_error_alpha = 0.0, pi_error_beta = 0.0, pi_mutation_M10 = 0.3, delta_time = 0.20;
    lambda_edge = 2 * log(1.2), lambda_root = 2 * log(1.2);
    alpha_vector.reserve(niter);
    alpha_vector.push_back(alpha);
    beta_vector.reserve(niter);
    beta_vector.push_back(beta);
    m10_vector.reserve(niter);
    m10_vector.push_back(m10);
    rootage_vector.reserve(niter);
    rootage_vector.push_back(rootage);

    // Two values may need to generate from R
    //double x[4] = {0.344875, 0.1136095, 0.3938868, 0.1476288}; // rdirichlet
    std::string br_tree_string = "((2:0.4994737979,3:0.4994737979):0.5293488252,((4:0.4210867618,1:0.4210867618):0.1298891809,5:0.5509759427):0.4778466805);";
    TreeManip br_tree;
    br_tree.buildFromNewick(br_tree_string, true, false);
    std::cout << "br_tree: " << br_tree.makeNewick(3) << std::endl;

    Tmax = 1.045;
    //br_tree.scaleAllEdgeLengths(Tmax);
    std::cout << "br_tree scaled: " << br_tree.makeNewick(3) << std::endl;
    // find map between tip node name and node number
    br_tree.buildNodeNameAndNumberMap();

    //TreeManip br_tree_exchange;
    // br_tree_exchange.buildFromNewick(br_tree.makeNewick(8), true, false);
    //double time_proposal_ratio = 1.0;
    //br_tree_exchange.narrowExchangeFrom(br_tree);
    //std::cout << "br_tree_exchange: " << br_tree_exchange.makeNewick(3) << std::endl;

    // 4. Run iterations with many iterations and four differnet conditions
    runTreeAndOrderSearchDP(tree_data_original, data_label, true_tree_tm, br_tree, niter);

    std::cout << "average of alpha_vector: " << getAverageOfVector(alpha_vector) << "\n";
    std::cout << "average of beta_vector: " << getAverageOfVector(beta_vector) << "\n";
    std::cout << "average of m10_vector: " << getAverageOfVector(m10_vector) << "\n";
    std::cout << "max of m10_vector: " << getMaxOfVector(m10_vector) << "\n";
    std::cout << "min of m10_vector: " << getMinOfVector(m10_vector) << "\n";
    std::cout << "average of rootage_vector: " << getAverageOfVector(rootage_vector) << "\n";

    // 5. Output results to files
    outputToCSVFile(output_filename);

    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start); 
    std::cout << "Time taken by function: " << duration.count() << " seconds\n"; 

    return 0;
}