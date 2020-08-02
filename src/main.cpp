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

/**
int main(int argc, char *argv[])
{
    // Get starting timepoint 
    auto start = std::chrono::high_resolution_clock::now(); 

    generator.seed(time(0));

    // 0. read arguments
    // Read the arguments passing from command line and validate
    int opt, niter = 1000000;
    std::string output_filename, br_tree_string;
    while ((opt = getopt(argc, argv, "n:f:t:")) != -1)
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
        case 't':
            br_tree_string = optarg;
            std::cout<<"start tree: " <<br_tree_string<<"\n";
            break;
        default:
            fprintf(stderr, "Wrong arguments!\n");
            exit(1);
        }
    }

    // 1. Start with input random tree file and input observed data
    std::string true_tree_filename = "./input/tree/RandomTreeScale_3.tre";
    std::string true_time_tree_filename = "./input/tree/RandomTree_3.tre";
    std::string tree_data_filename = "./input/data/Long_Obs_binary_obs_0_1_tip_mu_005_alpha_001_beta_001_matrix_1.csv";

    // 2. Build tree and get tree info
    TreeManip true_tree_tm, true_time_tree_tm;
    true_tree_tm.buildFromNewickFile(true_tree_filename);
    true_time_tree_tm.buildFromNewickFile(true_time_tree_filename);

    std::cout << "true_tree_tm: " << true_tree_tm.makeNewick(3) << std::endl;
    std::cout << "true_time_tree_tm: " << true_time_tree_tm.makeNewick(3) << std::endl;

    std::vector<std::string> data_label, mutation_data_label;
    std::vector<std::vector<int>> tree_data_original, mutation_tree_data_original;
    readObservedDataFromCSVFile(5000, tree_data_filename, data_label, tree_data_original);

    // 3. initialize value based on some distribution (normal, dirichlet)
    alpha = getUniformDistribution(0, 0.45);
    beta = getUniformDistribution(0, 0.45);
    m01 = 1, m00 = -m01, rootage = 1.045;
    m10 = getUniformDistribution(0, 1.0);
    m11 = -m10;
    pi_error_alpha = 0.1, pi_error_beta = 0.2, pi_mutation_M10 = 0.3, delta_time = 0.20;
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
    //std::string br_tree_string = "((4:0.1533530955,2:0.1533530955):0.8916469045,(5:0.6116008526,(3:0.338719725,1:0.338719725):0.2728811275):0.4333991474);";
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
    bool run_with_mutation = false;
    runTreeAndOrderSearchDP(tree_data_original, mutation_tree_data_original, data_label, mutation_data_label, br_tree, run_with_mutation, niter);

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
*/

int main(int argc, char *argv[])
{

    // Get starting timepoint
    auto start = std::chrono::high_resolution_clock::now();

    generator.seed(time(0));

    // 0. read arguments
    // Read the arguments passing from command line and validate
    int opt, niter = 1000000;
    std::string output_filename, br_tree_string;
    while ((opt = getopt(argc, argv, "n:f:t:")) != -1)
    {
        switch (opt)
        {
        case 'n':
            niter = atoi(optarg);
            break;
        case 'f':
            output_filename = optarg;
            std::cout << "output_filename: " << output_filename << "\n";
            break;
        case 't':
            br_tree_string = optarg;
            std::cout << "start tree: " << br_tree_string << "\n";
            break;
        default:
            fprintf(stderr, "Wrong arguments!\n");
            exit(1);
        }
    }

    // 1. Start with input random tree file and input observed data
    std::string true_tree_filename = "./input/tree/RandomTreeScale_3.tre";
    std::string true_time_tree_filename = "./input/tree/RandomTree_3.tre";
    std::string tree_data_filename = "./input/data/Long_Obs_binary_obs_0_1_tip_mu_005_alpha_001_beta_001_matrix_1.csv";

    // 2. Build tree and get tree info
    readObservedDataFromCSVFile(5000, tree_data_filename);

    for (auto l : data_label)
    {
        std::cout << "l: " << l << std::endl;
    }
    for (auto l : tree_data_original[0])
    {
        std::cout << "data: " << l << std::endl;
    }
    for (auto l : mutation_tree_data[0])
    {
        std::cout << "mu data: " << l << std::endl;
    }
    std::cout << "tree_data_original: " << tree_data_original.size() << std::endl;
    std::cout << "mutation_tree_data: " << mutation_tree_data.size() << std::endl;

    // 3. initialize value based on some distribution (normal, dirichlet)
    alpha = getUniformDistribution(0, 0.45);
    beta = getUniformDistribution(0, 0.45);
    m01 = 1, m00 = -m01, rootage = 1.045;
    m10 = getUniformDistribution(0, 1.0);
    m11 = -m10;
    pi_error_alpha = 0.1, pi_error_beta = 0.2, pi_mutation_M10 = 0.3, delta_time = 0.20;
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
    //std::string br_tree_string1 = "((4:0.1533530955,2:0.1533530955):0.8916469045,(5:0.6116008526,(3:0.338719725,1:0.338719725):0.2728811275):0.4333991474);";
    //"((2:0.9883721203,3:0.9883721203):0.05670532564,(5:0.8869449045,(4:0.7599123994,1:0.7599123994):0.1270325051):0.1581325414);"
    TreeManip br_tree_tm;
    br_tree_tm.buildFromNewick(br_tree_string, true, false);
    std::cout << "br_tree: " << br_tree_tm.makeNewick(3) << std::endl;

    // TreeManip br_tree_tm1;
    // br_tree_tm1.buildFromNewick(br_tree_string1, true, false);
    // for (auto nd : br_tree_tm.getTree()->getPreOrder())
    // {
    //     if (nd->getName().size() == 0) 
    //     {
    //         nd->setName(std::to_string(nd->getNumber() + 1));
    //     }
    //     std::cout << "nd: " << nd->getNumber() <<", name " << nd->getName() << std::endl;
    // }
    // std::cout << "br_tree_tm: " << br_tree_tm.makeNewick(3, true) << std::endl;

    Tmax = 0.9125;
    //br_tree_tm.scaleAllEdgeLengths(Tmax);
    //std::cout << "br_tree_tm scaled: " << br_tree_tm.makeNewick(3) << std::endl;
    // find map between tip node name and node number
    br_tree_tm.buildNodeNameAndNumberMap();

    // 4. Run iterations with many iterations and four differnet conditions
    bool run_with_mutation = true;
    runTreeAndOrderSearchDP(br_tree_tm, run_with_mutation, niter);

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

