#include <iostream>
#include <string>
#include <random>
#include <vector>
#include <math.h>
#include <unistd.h>
#include <chrono>
#include <ctime>
#include <fstream>

#include "tree_manip.hpp"
#include "utility.hpp"
#include "xstrom.hpp"

using namespace strom;

const double Node::_smallest_edge_length = 1.0e-12;
double substitution_matrix[2][2], error_matrix[2][2];
double const STD_DEV = 0.005;
double const STD_DEV_M10 = 0.02;
double bd_alpha_alpha, bd_beta_alpha, bd_alpha_beta, bd_beta_beta, rootage_mean, rootage_stddev;
std::discrete_distribution<int> discrete_distribution{1, 1, 2, 2, 2, 1};
bool verbose = false;

void mapObservedStateToRightLabel(TreeManip &start_tree_tm, std::vector<int> &observed_state, std::vector<int> &mapped_observed_state, std::vector<std::string> &data_label)
{
    int id = 0;
    for (auto state : observed_state)
    {
        std::string label = data_label[id];
        mapped_observed_state[start_tree_tm.getNodeNumberByName(label)] = state;
        id++;
    }
}

void readObservedDataFromCSVFile(const int num_lines, const std::string tree_data_filename, std::vector<std::string> &data_label, std::vector<std::vector<int>> &tree_data,
                                 std::vector<std::vector<int>> &mutation_tree_data_original)
{
    // read input tree file
    std::ifstream fin;
    fin.open(tree_data_filename);
    std::string label, line;
    getline(fin, label);
    // erase " and "t" from data
    label.erase(std::remove(label.begin(), label.end(), '"'), label.end());
    label.erase(std::remove(label.begin(), label.end(), 't'), label.end());
    label.erase(std::remove(label.begin(), label.end(), '\r'), label.end());
    label.erase(std::remove(label.begin(), label.end(), '\n'), label.end());

    // get label
    std::string delimiter = ",";
    //label.erase(0, 1); // remove first comma
    size_t pos = 0;
    int num;
    while ((pos = label.find(delimiter)) != std::string::npos)
    {
        data_label.push_back(label.substr(0, pos));
        label.erase(0, pos + delimiter.length());
    }
    // data_label.push_back(label);

    int line_i = 0;
    std::vector<int> line_data;
    while (getline(fin, line))
    {
        // erase " from data
        line.erase(std::remove(line.begin(), line.end(), '"'), line.end());
        //line.erase(0, line.find(delimiter) + delimiter.length()); // remove first number and comma

        while ((pos = line.find(delimiter)) != std::string::npos)
        {
            num = (line.substr(0, pos))[0] - '0';
            line.erase(0, pos + delimiter.length());
            line_data.push_back(num);
        }
        int label = line[0] - '0';
        if (label == 0 && line_i < num_lines)
        {
            tree_data.push_back(line_data);
        }
        else if (label == 1)
        {
            mutation_tree_data_original.push_back(line_data);
        }

        line_data.clear();
        line_i++;
    }

    fin.close();
}

void outputToCSVFile(const std::string output_filename, std::vector<double> &alpha_vector, std::vector<double> &beta_vector, std::vector<double> &m10_vector,
                     std::vector<double> &rootage_vector, std::vector<double> &loglikelihood_vector, std::vector<std::string> &tree_string_vector)
{
    std::ofstream myfile;
    myfile.open(output_filename + ".csv");
    myfile << "alpha, beta, m10, rootage, loglikelihood\n";
    for (int i = 0; i < alpha_vector.size(); i++)
    {
        myfile << alpha_vector[i] << ", ";
        myfile << beta_vector[i] << ", ";
        myfile << m10_vector[i] << ", ";
        myfile << rootage_vector[i] << ", ";
        myfile << loglikelihood_vector[i] << "\n";
    }

    myfile.close();

    myfile.open(output_filename + "_trees.txt");
    for (int i = 0; i < tree_string_vector.size(); i++)
    {
        myfile << tree_string_vector[i];
    }
    myfile.close();
}

double getAverageOfVector(std::vector<double> &vector, int ignore_first_n = 0)
{
    double result = 0.0;
    for (int i = 0; i < vector.size(); i++)
    {
        if (i < ignore_first_n)
        {
            continue;
        }
        result += vector[i];
    }
    return (result / (vector.size() - ignore_first_n));
}

double getMaxOfVector(std::vector<double> &vector)
{
    double max = -std::numeric_limits<double>::max();
    for (auto v : vector)
    {
        max = std::max(max, v);
    }
    return max;
}

double getMinOfVector(std::vector<double> &vector)
{
    double min = std::numeric_limits<double>::max();
    for (auto v : vector)
    {
        min = std::min(min, v);
    }
    return min;
}

double adjustInfinity(double ratio)
{
    if (ratio == std::numeric_limits<double>::infinity())
    {
        return std::numeric_limits<double>::max();
    }
    else if (ratio == -std::numeric_limits<float>::infinity())
    {
        return -std::numeric_limits<double>::max();
    }
    else
    {
        return ratio;
    }
}

int encodeStateToValue(std::vector<int> state)
{
    int value = 0;
    for (auto s : state)
    {
        value <<= 1;
        value += s;
    }
    return value;
}

void encodeValueToState(int value, int num_tips, std::vector<int> &state)
{
    int i = 0;
    while (i < num_tips)
    {
        state.insert(state.begin(), value & 1);
        value >>= 1;
        i++;
    }
}

void buildCountsTableForTreeData(std::vector<std::vector<int>> &tree_data_original, std::map<int, int> &state_code_and_count_map)
{
    // loop over all tree data, and increment corresponding encode value's count for each line
    for (auto state : tree_data_original)
    {
        int encode_value = encodeStateToValue(state);
        auto it = state_code_and_count_map.find(encode_value);
        if (it == state_code_and_count_map.end())
        {
            state_code_and_count_map[encode_value] = 1;
        }
        else
        {
            (it->second)++;
        }
    }
}

void updateSubstitutionMatrix(double br_len, double m10)
{
    // Note to change m01
    double m01 = 1.0;
    double divisor = m01 + m10;
    double exp_len = exp((-m01 - m10) * br_len);

    substitution_matrix[0][0] = (m10 + m01 * exp_len) / divisor;
    substitution_matrix[0][1] = (m01 - m01 * exp_len) / divisor;
    substitution_matrix[1][0] = (m10 - m10 * exp_len) / divisor;
    substitution_matrix[1][1] = (m01 + m10 * exp_len) / divisor;
}

void updateErrorMatrix(double alpha, double beta)
{
    error_matrix[0][0] = 1 - alpha;
    error_matrix[0][1] = alpha;
    error_matrix[1][0] = beta;
    error_matrix[1][1] = 1 - beta;
}

void buildLeaveLikelihoodMatrix(std::vector<int> &observed_state, std::vector<std::vector<double>> &leaves_likelihood)
{
    int leave_num = 0;
    for (auto state : observed_state)
    {
        if (state == 1)
        {
            leaves_likelihood[0][leave_num] = error_matrix[0][1];
            leaves_likelihood[1][leave_num] = error_matrix[1][1];
        }
        else
        {
            leaves_likelihood[0][leave_num] = error_matrix[0][0];
            leaves_likelihood[1][leave_num] = error_matrix[1][0];
        }
        leave_num++;
    }
}

double felsensteinBinaryDP(std::vector<int> &observed_state, Node *cur_node, int cur_state, std::vector<std::vector<double>> &likelihood_table,
                           std::vector<std::vector<double>> &leaves_likelihood, double m10)
{
    if (likelihood_table[cur_state][cur_node->getNumber()] != 0.0)
    {
        return likelihood_table[cur_state][cur_node->getNumber()];
    }

    if (cur_node->getLeftChild() == NULL)
    {
        // if it's a tip, return the likelihood based on leaves_likelihood
        likelihood_table[cur_state][cur_node->getNumber()] = leaves_likelihood[cur_state][cur_node->getNumber()];
        return likelihood_table[cur_state][cur_node->getNumber()];
    }
    else
    {
        // if not tip, get its children
        Node *left_child = cur_node->getLeftChild();
        Node *right_child = cur_node->getRightChild();
        updateSubstitutionMatrix(left_child->getEdgeLength(), m10);
        double pi0_left_child = substitution_matrix[cur_state][0];
        double pi1_left_child = substitution_matrix[cur_state][1];
        updateSubstitutionMatrix(right_child->getEdgeLength(), m10);
        double pi0_right_child = substitution_matrix[cur_state][0];
        double pi1_right_child = substitution_matrix[cur_state][1];

        double likelihood_left_child_0 = felsensteinBinaryDP(observed_state, left_child, 0, likelihood_table, leaves_likelihood, m10);
        double likelihood_left_child_1 = felsensteinBinaryDP(observed_state, left_child, 1, likelihood_table, leaves_likelihood, m10);
        double likelihood_right_child_0 = felsensteinBinaryDP(observed_state, right_child, 0, likelihood_table, leaves_likelihood, m10);
        double likelihood_right_child_1 = felsensteinBinaryDP(observed_state, right_child, 1, likelihood_table, leaves_likelihood, m10);

        double likelihood = (pi0_left_child * likelihood_left_child_0 + pi1_left_child * likelihood_left_child_1) *
                            (pi0_right_child * likelihood_right_child_0 + pi1_right_child * likelihood_right_child_1);
        likelihood_table[cur_state][cur_node->getNumber()] = likelihood;
        return likelihood;
    }
}

double treeCompareFelsensteinDP(std::vector<std::vector<int>> &tree_data_original, std::vector<std::string> &data_label, TreeManip &start_tree_tm, std::map<int, int> &state_code_and_count_map,
                                std::vector<std::vector<double>> &likelihood_table, std::vector<std::vector<double>> &leaves_likelihood, double m10)
{
    // initialize values
    int num_tips = start_tree_tm.getNumLeaves();

    // std::cout << "num_tips: " << num_tips << "\n";
    // std::cout << "tree_data_original size: " << tree_data_original.size() << "\n";
    // std::cout << "data_label size: " << data_label.size() << "\n";
    // std::cout << "state_code_and_count_map size: " << state_code_and_count_map.size() << "\n";
    // std::cout << "likelihood_table size: " << likelihood_table.size() << "\n";
    // std::cout << "leaves_likelihood size: " << leaves_likelihood.size() << "\n";
    // std::cout << "m10: " << m10 << "\n";

    Node *root = start_tree_tm.getRootNode()->getLeftChild();
    double total_loglikelihood = 0.0;
    for (auto it = state_code_and_count_map.begin(); it != state_code_and_count_map.end(); it++)
    {
        // transfer value to state vector
        int value = it->first;
        int count = it->second;
        // std::cout << "value: " << value << ", count: " << count << "\n";
        std::vector<int> observed_state;
        observed_state.reserve(num_tips);
        encodeValueToState(value, num_tips, observed_state);
        std::vector<int> mapped_observed_state(observed_state);
        mapObservedStateToRightLabel(start_tree_tm, observed_state, mapped_observed_state, data_label);

        buildLeaveLikelihoodMatrix(mapped_observed_state, leaves_likelihood);
        double likelihood = felsensteinBinaryDP(mapped_observed_state, root, 0, likelihood_table, leaves_likelihood, m10);
        total_loglikelihood += log(likelihood) * count;
        std::fill(likelihood_table[0].begin(), likelihood_table[0].end(), 0);
        std::fill(likelihood_table[1].begin(), likelihood_table[1].end(), 0);
    }
    return total_loglikelihood;
}

static void update_tree(int random_sample, TreeManip &cur_tree_tm, TreeManip &proposed_tree_tm, double &proposed_rootage, double &topology_prior_ratio,
                        double &time_prior_ratio, double &rate_prior_ratio, double &topology_proposal_ratio, double &time_proposal_ratio, double &rate_proposal_ratio, double delta_time)
{
    switch (random_sample)
    {
    case 0:
    {
        std::cout << "narrow exchange.\n";
        proposed_tree_tm.narrowExchangeFrom(cur_tree_tm);
        break;
    }
    case 1:
        std::cout << "wide exchange.\n";
        proposed_tree_tm.wideExchangeFrom(cur_tree_tm);
        break;
    case 2:
        std::cout << "FNPR exchange.\n";
        proposed_tree_tm.FNPRExchangeFrom(cur_tree_tm);
        break;
    case 3:
        std::cout << "Total length change.\n";
        proposed_tree_tm.buildFromNewick(cur_tree_tm.makeNewick(8), true, false);
        proposed_tree_tm.totalLengthChange(proposed_rootage, rate_prior_ratio, rate_proposal_ratio, rootage_mean, rootage_stddev);
        proposed_tree_tm.updateNodesHeightInfo();
        break;
    case 4:
        std::cout << "Random length change.\n";
        proposed_tree_tm.buildFromNewick(cur_tree_tm.makeNewick(8), true, false);
        proposed_tree_tm.updateNodesHeightInfo();
        proposed_tree_tm.randomLengthChange(delta_time, time_proposal_ratio);
        proposed_tree_tm.updateNodesHeightInfo();
        break;
    case 5:
        std::cout << "All internal length change.\n";
        proposed_tree_tm.buildFromNewick(cur_tree_tm.makeNewick(8), true, false);
        proposed_tree_tm.updateNodesHeightInfo();
        proposed_tree_tm.allInternalLengthChange(delta_time, time_proposal_ratio);
        proposed_tree_tm.updateNodesHeightInfo();
        break;
    default:
        std::cout << "Error: Now we only support 7 different cases.\n";
    }
    proposed_tree_tm.buildNodeNameAndNumberMap();
    proposed_rootage = proposed_tree_tm.getTreeMaxHeight();
    std::cout << "proposed_rootage end " << proposed_rootage << "\n";
}

void assignOneToAllChildNodesAndSelf(Node *nd, std::vector<int> &state)
{
    state[nd->getNumber()] = 1;
    Node *left_child = nd->getLeftChild();
    if (left_child != NULL)
    {
        assignOneToAllChildNodesAndSelf(left_child, state);
        assignOneToAllChildNodesAndSelf(left_child->getRightSib(), state);
    }
}

void buildPossibilityMatrix(TreeManip &br_tree_tm, std::vector<std::vector<int>> &possibility_matrix)
{
    possibility_matrix.clear();
    Node *root = br_tree_tm.getRootNode()->getLeftChild();
    int num_nodes = br_tree_tm.getNumNodes() - 1;
    for (auto nd : br_tree_tm.getPreOrder())
    {
        if (nd == root)
        {
            continue;
        }

        std::vector<int> state(num_nodes, 0);
        assignOneToAllChildNodesAndSelf(nd, state);
        possibility_matrix.push_back(state);
    }
}

double getCorrespondingPossibilityFromState(Node *nd, int nd_state, int nd_parent_state)
{
    if (nd_parent_state == 0 && nd_state == 0)
        return nd->getP00();
    else if (nd_parent_state == 0 && nd_state == 1)
        return nd->getP01();
    else if (nd_parent_state == 1 && nd_state == 0)
        return nd->getP10();
    else if (nd_parent_state == 1 && nd_state == 1)
        return nd->getP11();
    std::cout<<"Error: should not reach here\n";
    return 0.0;
}

void normalizeVector(std::vector<double> &value_vector)
{
    double total_value = 0.0;
    for (auto &p : value_vector)
        total_value += p;

    for (auto &p : value_vector)
        p /= total_value;
    for (auto v : value_vector) {
        if (isnan(v)) {
            std::cout<<"Error: find NAN\n";
            break;
        }
    }
}

void buildPossibilityValueVector(TreeManip &br_tree_tm, std::vector<std::vector<int>> &possibility_matrix, std::vector<double> &possibility_value_vector)
{
    int i = 0;
    Node *root = br_tree_tm.getRootNode()->getLeftChild();
    for (auto state : possibility_matrix)
    {
        double p = 1.0;
        for (auto nd : br_tree_tm.getPreOrder())
        {
            if (nd == root)
            {
                continue;
            }

            p *= getCorrespondingPossibilityFromState(nd, state[nd->getNumber()], state[nd->getParent()->getNumber()]);
        }
        possibility_value_vector[i] = p;
        i++;
    }

    normalizeVector(possibility_value_vector);
}

void computeErrorMatrixOnObservedDate(TreeManip &br_tree_tm, std::vector<std::vector<int>> &tree_data_original, std::vector<std::string> &data_label, std::vector<std::vector<double>> &leaves_likelihood,
                                      std::vector<double> &possibility_value_vector, std::vector<std::vector<int>> &possibility_matrix, std::vector<std::vector<double>> &error_matrix_on_observed_data)
{
    int num_tips = br_tree_tm.getNumLeaves();
    leaves_likelihood.resize(2); // 2 * num_tips
    for (auto &it : leaves_likelihood)
    {
        it.resize(num_tips);
    }

    for (auto observed_state : tree_data_original)
    {
        std::vector<int> mapped_observed_state(observed_state);
        mapObservedStateToRightLabel(br_tree_tm, observed_state, mapped_observed_state, data_label);

        buildLeaveLikelihoodMatrix(mapped_observed_state, leaves_likelihood);

        std::vector<double> error_vector(possibility_value_vector);
        for (int i = 0; i < possibility_matrix.size(); i++)
        {
            auto state = possibility_matrix[i];

            for (auto nd : br_tree_tm.getPreOrder())
            {
                if (nd->getLeftChild() == NULL)
                {
                    error_vector[i] *= leaves_likelihood[state[nd->getNumber()]][nd->getNumber()];
                    //std::cout << "nd name: " << nd->getName() << ", error: "<< leaves_likelihood[state[nd->getNumber()]][nd->getNumber()] <<". ";
                }
            }
        }
        //std::cout << "\n";
        error_matrix_on_observed_data.push_back(error_vector);
    }
}

std::vector<double> computeMutationPossibilitiesForCurrentTree(TreeManip &start_tree_tm, std::vector<std::vector<int>> &mutation_tree_data, std::vector<std::string> &data_label, 
    std::vector<int> &mutation_branch_vector, std::vector<std::vector<int>> &possibility_matrix, std::vector<double> &possibility_value_vector, std::vector<std::vector<double>> &leaves_likelihood)
{
    std::vector<double> loglikelihoods(2, 0.0);
    // 0. compute normalized possibility
    std::fill(possibility_value_vector.begin(), possibility_value_vector.end(), 1.0);
    buildPossibilityValueVector(start_tree_tm, possibility_matrix, possibility_value_vector);

    // 1. compute error on observed data
    std::vector<std::vector<double>> error_matrix_on_observed_data;
    computeErrorMatrixOnObservedDate(start_tree_tm, mutation_tree_data, data_label, leaves_likelihood, possibility_value_vector, possibility_matrix, error_matrix_on_observed_data);

    // 2.save mutation branch on all observed data
    int i = 0;
    double mutation_loglikehood = 0.0, mutation_scaled_prior_loglikelihood = 0.0;

    for (auto error_vector : error_matrix_on_observed_data)
    {
        mutation_loglikehood += log(error_vector[mutation_branch_vector[i]]);

        // scale error_vector
        normalizeVector(error_vector);
        mutation_scaled_prior_loglikelihood += log(error_vector[mutation_branch_vector[i]]);

        i++;
    }
    loglikelihoods[0] = mutation_loglikehood;
    loglikelihoods[1] = mutation_scaled_prior_loglikelihood;
    return loglikelihoods;
}

std::vector<double> computeMutationPossibilities(TreeManip &start_tree_tm, std::vector<int> &mutation_branch_vector, std::vector<std::vector<int>> &possibility_matrix, std::vector<double> &possibility_value_vector, 
    std::vector<std::vector<int>> &mutation_tree_data, std::vector<std::string> &data_label, std::vector<std::vector<double>> &leaves_likelihood)
{
    //std::cout << "current error_matrix[0][1]: " << error_matrix[0][1] << "\n";
    mutation_branch_vector.clear();

    std::vector<double> loglikelihoods(2, 0.0);
    // 1. build all possibility matrix for current tree
    buildPossibilityMatrix(start_tree_tm, possibility_matrix);

    // 2. compute normalized possibility
    if (possibility_value_vector.empty())
    {
        possibility_value_vector.resize(possibility_matrix.size());
    }
    std::fill(possibility_value_vector.begin(), possibility_value_vector.end(), 1.0);
    buildPossibilityValueVector(start_tree_tm, possibility_matrix, possibility_value_vector);

    // 3. compute error on observed data
    std::vector<std::vector<double>> error_matrix_on_observed_data;
    computeErrorMatrixOnObservedDate(start_tree_tm, mutation_tree_data, data_label, leaves_likelihood, possibility_value_vector, possibility_matrix, error_matrix_on_observed_data);

    // 4. do sampling and save mutation branch on all observed data
    int num_nodes = start_tree_tm.getNumNodes() - 1, i = 0;
    double mutation_loglikehood = 0.0, mutation_scaled_prior_loglikelihood = 0.0;

    for (auto error_vector : error_matrix_on_observed_data)
    {
        //std::discrete_distribution<int> error_discrete_distribution(error_vector.begin(), error_vector.end());
        //int random_sample = error_discrete_distribution(generator);
        int random_sample = 0;
        mutation_branch_vector.push_back(random_sample);
        mutation_loglikehood += log(error_vector[mutation_branch_vector.back()]);

        // scale error_vector
        normalizeVector(error_vector);
        mutation_scaled_prior_loglikelihood += log(error_vector[mutation_branch_vector.back()]);

        i++;
    }
    loglikelihoods[0] = mutation_loglikehood;
    loglikelihoods[1] = mutation_scaled_prior_loglikelihood;
    return loglikelihoods;
}

void runTreeAndOrderSearchDP(std::vector<std::vector<int>> &tree_data_original, std::vector<std::vector<int>> &mutation_tree_data_original, std::vector<std::string> &data_label, TreeManip &start_tree_tm, int niter, std::vector<double> &alpha_vector, std::vector<double> &beta_vector,
                             std::vector<double> &m10_vector, std::vector<double> &rootage_vector, std::vector<double> &loglikelihood_vector, std::vector<std::string> &tree_string_vector)
{
    // initialize all values
    int num_tips = start_tree_tm.getNumLeaves();
    start_tree_tm.updateNodesHeightInfo();
    int num_nodes = start_tree_tm.getNumNodes();
    std::vector<double> possibility_value_vector;
    std::vector<std::vector<int>> possibility_matrix;

    std::vector<std::vector<double>> leaves_likelihood, likelihood_table;
    leaves_likelihood.resize(2); // 2 * num_tips
    likelihood_table.resize(2);  // 2 * num_nodes
    for (auto &it : leaves_likelihood)
    {
        it.resize(num_tips);
    }
    for (auto &it : likelihood_table)
    {
        it.resize(num_nodes);
    }
    std::cout << "start_tree_tm: " << start_tree_tm.makeNewick(8) << std::endl;

    // build counts table for all possible observed data
    std::map<int, int> state_code_and_count_map;
    buildCountsTableForTreeData(tree_data_original, state_code_and_count_map);

    double alpha = alpha_vector.front();
    double beta = beta_vector.front();
    updateErrorMatrix(alpha, beta);
    double m10 = m10_vector.front();
    double m01 = 1.0;

    // run felsensteinBinaryDP for tree
    double total_loglikelihood = treeCompareFelsensteinDP(tree_data_original, data_label, start_tree_tm, state_code_and_count_map, likelihood_table, leaves_likelihood, m10);
    std::cout << "total_loglikelihood: " << total_loglikelihood << "\n";

    // initialize vectors for iterations
    start_tree_tm.buildNodesPossibilitiesInfo(m10, m01);
    TreeManip start_tree_tm_copy, cur_tree_manip;
    start_tree_tm_copy.buildFromNewick(start_tree_tm.makeNewick(8), true, false);
    start_tree_tm_copy.updateNodesHeightInfo();
    start_tree_tm_copy.buildNodeNameAndNumberMap();
    cur_tree_manip.buildFromNewick(start_tree_tm.makeNewick(8), true, false);
    cur_tree_manip.updateNodesHeightInfo();
    cur_tree_manip.buildNodeNameAndNumberMap();
    cur_tree_manip.buildNodesPossibilitiesInfo(m10, m01);
    start_tree_tm_copy.addTToName();
    tree_string_vector.push_back(start_tree_tm_copy.makeNewick(8, true) + ";\n");
    loglikelihood_vector.reserve(niter);
    loglikelihood_vector.push_back(total_loglikelihood);
    double pi_error_alpha = 0.1, pi_error_beta = 0.2, pi_mutation_M10 = 0.3, delta_time = 0.2;

    // run with mutation, compute initial branch possibilities
    double cur_mutation_loglikehood = 0.0, cur_mutation_scaled_prior_loglikelihood = 0.0;
    std::vector<int> cur_mutation_branch_vector;
    std::vector<double> likelihoods = computeMutationPossibilities(start_tree_tm, cur_mutation_branch_vector, possibility_matrix, possibility_value_vector, mutation_tree_data_original, data_label, leaves_likelihood);
    cur_mutation_loglikehood = likelihoods[0];
    cur_mutation_scaled_prior_loglikelihood = likelihoods[1];
    std::cout << "mutation_loglikehood: " << cur_mutation_loglikehood << "\n";
    std::cout << "mutation_scaled_prior_loglikelihood: " << cur_mutation_scaled_prior_loglikelihood << "\n";

    // start component-wise update
    for (int i = 1; i < niter; i++)
    {
        double r = getUniformDistribution(0, 1);
        std::cout << "iter: " << i << "\n";
        if (r < pi_error_alpha)
        {
            std::cout << "update alpha.\n";

            double cur_alpha = alpha_vector[i - 1];
            double cur_beta = beta_vector[i - 1];
            double cur_m10 = m10_vector[i - 1];
            double cur_loglikelihood = loglikelihood_vector[i - 1];
            double cur_rootage = rootage_vector[i - 1];
            updateErrorMatrix(cur_alpha, cur_beta);
            m10 = cur_m10;

            if (verbose)
                std::cout << "initial tree before update alpha: " << cur_tree_manip.makeNewick(8) << "\n";

            rootage_vector.push_back(cur_rootage);
            m10_vector.push_back(cur_m10);
            TreeManip cur_tree_manip_copy;
            cur_tree_manip_copy.buildFromNewick(cur_tree_manip.makeNewick(8), true, false);
            cur_tree_manip_copy.addTToName();
            tree_string_vector.push_back(cur_tree_manip_copy.makeNewick(8, true) + ";\n");

            double proposed_alpha = getNormalDistribution(cur_alpha, STD_DEV);
            // adjust it between (0, 1)
            if (proposed_alpha > 1)
            {
                proposed_alpha = 2 - proposed_alpha;
            }
            else if (proposed_alpha < 0)
            {
                proposed_alpha = -proposed_alpha;
            }

            // proposed loglikelihood
            updateErrorMatrix(proposed_alpha, cur_beta);
            double proposed_loglikelihood = treeCompareFelsensteinDP(tree_data_original, data_label, cur_tree_manip, state_code_and_count_map, likelihood_table, leaves_likelihood, m10);

            // proposed mutation loglikelihood
            std::vector<double> likelihoods = computeMutationPossibilitiesForCurrentTree(cur_tree_manip, mutation_tree_data_original, data_label, cur_mutation_branch_vector, possibility_matrix, possibility_value_vector, leaves_likelihood);
            double proposed_mutation_loglikelihood = likelihoods[0];
            double proposed_mutation_scaled_prior_loglikelihood = likelihoods[1];
            std::cout << "proposed_mutation_loglikelihood: " << proposed_mutation_loglikelihood << "\n";
            std::cout << "proposed_mutation_scaled_prior_loglikelihood: " << proposed_mutation_scaled_prior_loglikelihood << "\n";
            std::cout << "cur_loglikelihood: " << cur_loglikelihood << "\n";
            std::cout << "cur_mutation_loglikelihood: " << cur_mutation_loglikehood << "\n";
            std::cout << "cur_mutationscaled_prior_loglikelihood: " << cur_mutation_scaled_prior_loglikelihood << "\n";

            // std::cout<<"cur_alpha: "<< cur_alpha<<"\n";
            // std::cout<<"proposed_alpha: "<< proposed_alpha<<"\n";
            if (verbose)
            {
                std::cout << "cur_loglikelihood: " << cur_loglikelihood << "\n";
                std::cout << "proposed_loglikelihood: " << proposed_loglikelihood << "\n";
            }
            double dnorm_on_cur_alpha = getNormalDistributionDensity(cur_alpha, proposed_alpha, STD_DEV);
            double dnorm_on_proposed_alpha = getNormalDistributionDensity(proposed_alpha, cur_alpha, STD_DEV);
            double alpha_error_proposal_ratio = exp(log(dnorm_on_cur_alpha) - log(dnorm_on_proposed_alpha));

            double dbeta_on_proposed_alpha = getBetaDistributionDensity(proposed_alpha, bd_alpha_alpha, bd_beta_alpha);
            double dbeta_on_cur_alpha = getBetaDistributionDensity(cur_alpha, bd_alpha_alpha, bd_beta_alpha);
            double alpha_error_Prior_ratio = exp(log(dbeta_on_proposed_alpha) - log(dbeta_on_cur_alpha));

            double target_error_ratio = exp(proposed_loglikelihood + proposed_mutation_loglikelihood - cur_loglikelihood - cur_mutation_loglikehood 
                + cur_mutation_scaled_prior_loglikelihood - proposed_mutation_scaled_prior_loglikelihood) * alpha_error_proposal_ratio * alpha_error_Prior_ratio;
            target_error_ratio = adjustInfinity(target_error_ratio);

            double r_error_ratio = getUniformDistribution(0, 1);
            if (r_error_ratio < target_error_ratio)
            {
                if (verbose)
                {
                    std::cout << "accept r_error_ratio: " << r_error_ratio << "\n";
                    std::cout << "accept target_error_ratio: " << target_error_ratio << "\n";
                }
                alpha_vector.push_back(proposed_alpha);
                beta_vector.push_back(cur_beta);
                loglikelihood_vector.push_back(proposed_loglikelihood);
                cur_mutation_loglikehood = proposed_mutation_loglikelihood;
                cur_mutation_scaled_prior_loglikelihood = proposed_mutation_scaled_prior_loglikelihood;
            }
            else
            {
                if (verbose)
                {
                    std::cout << "reject r_error_ratio: " << r_error_ratio << "\n";
                    std::cout << "reject target_error_ratio: " << target_error_ratio << "\n";
                }
                alpha_vector.push_back(cur_alpha);
                beta_vector.push_back(cur_beta);
                loglikelihood_vector.push_back(cur_loglikelihood);
            }
        }
        else if (r < pi_error_beta)
        {
            std::cout << "update beta.\n";

            double cur_alpha = alpha_vector[i - 1];
            double cur_beta = beta_vector[i - 1];
            double cur_m10 = m10_vector[i - 1];
            double cur_loglikelihood = loglikelihood_vector[i - 1];
            double cur_rootage = rootage_vector[i - 1];
            updateErrorMatrix(cur_alpha, cur_beta);
            m10 = cur_m10;

            if (verbose)
                std::cout << "initial tree before update beta: " << cur_tree_manip.makeNewick(8) << "\n";

            rootage_vector.push_back(cur_rootage);
            m10_vector.push_back(cur_m10);
            TreeManip cur_tree_manip_copy;
            cur_tree_manip_copy.buildFromNewick(cur_tree_manip.makeNewick(8), true, false);
            cur_tree_manip_copy.addTToName();
            tree_string_vector.push_back(cur_tree_manip_copy.makeNewick(8, true) + ";\n");

            double proposed_beta = getNormalDistribution(cur_beta, STD_DEV);
            // adjust it between (0, 1)
            if (proposed_beta > 1)
            {
                proposed_beta = 2 - proposed_beta;
            }
            else if (proposed_beta < 0)
            {
                proposed_beta = -proposed_beta;
            }

            if (verbose)
            {
                std::cout << "cur_beta: " << cur_beta << "\n";
                std::cout << "proposed_beta: " << proposed_beta << "\n";
            }

            // proposed loglikelihood
            updateErrorMatrix(cur_alpha, proposed_beta);
            double proposed_loglikelihood = treeCompareFelsensteinDP(tree_data_original, data_label, cur_tree_manip, state_code_and_count_map, likelihood_table, leaves_likelihood, m10);

            // proposed mutation loglikelihood
            std::vector<double> likelihoods = computeMutationPossibilitiesForCurrentTree(cur_tree_manip, mutation_tree_data_original, data_label, cur_mutation_branch_vector, possibility_matrix, possibility_value_vector, leaves_likelihood);
            double proposed_mutation_loglikelihood = likelihoods[0];
            double proposed_mutation_scaled_prior_loglikelihood = likelihoods[1];
            std::cout << "proposed_mutation_loglikelihood: " << proposed_mutation_loglikelihood << "\n";
            std::cout << "proposed_mutation_scaled_prior_loglikelihood: " << proposed_mutation_scaled_prior_loglikelihood << "\n";
            std::cout << "cur_loglikelihood: " << cur_loglikelihood << "\n";
            std::cout << "cur_mutation_loglikelihood: " << cur_mutation_loglikehood << "\n";
            std::cout << "cur_mutationscaled_prior_loglikelihood: " << cur_mutation_scaled_prior_loglikelihood << "\n";

            if (verbose)
            {
                std::cout << "cur_loglikelihood: " << cur_loglikelihood << "\n";
                std::cout << "proposed_loglikelihood: " << proposed_loglikelihood << "\n";
            }

            double dnorm_on_cur_beta = getNormalDistributionDensity(cur_beta, proposed_beta, STD_DEV);
            double dnorm_on_proposed_beta = getNormalDistributionDensity(proposed_beta, cur_beta, STD_DEV);
            double beta_error_proposal_ratio = exp(log(dnorm_on_cur_beta) - log(dnorm_on_proposed_beta));

            if (verbose)
            {
                std::cout << "dnorm_on_cur_beta: " << dnorm_on_cur_beta << "\n";
                std::cout << "dnorm_on_proposed_beta: " << dnorm_on_proposed_beta << "\n";
                std::cout << "beta_error_proposal_ratio: " << beta_error_proposal_ratio << "\n";
            }

            double dbeta_on_proposed_beta = getBetaDistributionDensity(proposed_beta, bd_alpha_beta, bd_beta_beta);
            double dbeta_on_cur_beta = getBetaDistributionDensity(cur_beta, bd_alpha_beta, bd_beta_beta);
            double beta_error_Prior_ratio = exp(log(dbeta_on_proposed_beta) - log(dbeta_on_cur_beta));

            if (verbose)
            {
                std::cout << "dbeta_on_proposed_beta: " << dbeta_on_proposed_beta << "\n";
                std::cout << "dbeta_on_cur_beta: " << dbeta_on_cur_beta << "\n";
                std::cout << "beta_error_Prior_ratio: " << beta_error_Prior_ratio << "\n";
            }
            
            double target_error_ratio = exp(proposed_loglikelihood + proposed_mutation_loglikelihood - cur_loglikelihood - cur_mutation_loglikehood 
                + cur_mutation_scaled_prior_loglikelihood - proposed_mutation_scaled_prior_loglikelihood) * beta_error_proposal_ratio * beta_error_Prior_ratio;
            target_error_ratio = adjustInfinity(target_error_ratio);

            double r_error_ratio = getUniformDistribution(0, 1);
            if (r_error_ratio < target_error_ratio)
            {
                if (verbose)
                {
                    std::cout << "accept r_error_ratio: " << r_error_ratio << "\n";
                    std::cout << "accept target_error_ratio: " << target_error_ratio << "\n";
                }
                alpha_vector.push_back(cur_alpha);
                beta_vector.push_back(proposed_beta);
                loglikelihood_vector.push_back(proposed_loglikelihood);
                cur_mutation_loglikehood = proposed_mutation_loglikelihood;
                cur_mutation_scaled_prior_loglikelihood = proposed_mutation_scaled_prior_loglikelihood;
            }
            else
            {
                if (verbose)
                {
                    std::cout << "reject r_error_ratio: " << r_error_ratio << "\n";
                    std::cout << "reject target_error_ratio: " << target_error_ratio << "\n";
                }
                alpha_vector.push_back(cur_alpha);
                beta_vector.push_back(cur_beta);
                loglikelihood_vector.push_back(cur_loglikelihood);
            }
        }
        else if (r < pi_mutation_M10)
        {
            std::cout << "update M10.\n";

            double cur_alpha = alpha_vector[i - 1];
            double cur_beta = beta_vector[i - 1];
            double cur_m10 = m10_vector[i - 1];
            double cur_loglikelihood = loglikelihood_vector[i - 1];
            double cur_rootage = rootage_vector[i - 1];
            updateErrorMatrix(cur_alpha, cur_beta);
            m10 = cur_m10;

            if (verbose)
                std::cout << "initial tree before update m10: " << cur_tree_manip.makeNewick(8) << "\n";

            rootage_vector.push_back(cur_rootage);
            alpha_vector.push_back(cur_alpha);
            beta_vector.push_back(cur_beta);
            TreeManip cur_tree_manip_copy;
            cur_tree_manip_copy.buildFromNewick(cur_tree_manip.makeNewick(8), true, false);
            cur_tree_manip_copy.addTToName();
            tree_string_vector.push_back(cur_tree_manip_copy.makeNewick(8, true) + ";\n");

            double proposed_m10 = std::abs(getNormalDistribution(cur_m10, STD_DEV_M10));
            if (verbose)
            {
                std::cout << "cur_m10: " << cur_m10 << "\n";
                std::cout << "proposed_m10: " << proposed_m10 << "\n";
            }

            // std::cout << "cur_tree_manip: " << cur_tree_manip.makeNewick(8) << std::endl;

            // proposed loglikelihood
            double proposed_loglikelihood = treeCompareFelsensteinDP(tree_data_original, data_label, cur_tree_manip, state_code_and_count_map, likelihood_table, leaves_likelihood, proposed_m10);

            if (verbose)
            {
                std::cout << "cur_loglikelihood: " << cur_loglikelihood << "\n";
                std::cout << "proposed_loglikelihood: " << proposed_loglikelihood << "\n";
            }

            // proposed mutation loglikelihood
            cur_tree_manip.buildNodesPossibilitiesInfo(proposed_m10, m01);
            std::vector<double> likelihoods = computeMutationPossibilitiesForCurrentTree(cur_tree_manip, mutation_tree_data_original, data_label, cur_mutation_branch_vector, possibility_matrix, possibility_value_vector, leaves_likelihood);
            double proposed_mutation_loglikelihood = likelihoods[0];
            double proposed_mutation_scaled_prior_loglikelihood = likelihoods[1];
            std::cout << "proposed_mutation_loglikelihood: " << proposed_mutation_loglikelihood << "\n";
            std::cout << "proposed_mutation_scaled_prior_loglikelihood: " << proposed_mutation_scaled_prior_loglikelihood << "\n";
            std::cout << "cur_loglikelihood: " << cur_loglikelihood << "\n";
            std::cout << "cur_mutation_loglikelihood: " << cur_mutation_loglikehood << "\n";
            std::cout << "cur_mutationscaled_prior_loglikelihood: " << cur_mutation_scaled_prior_loglikelihood << "\n";

            double dbeta_on_proposed_m10 = getUniformDistributionDensity(proposed_m10, 0, 1);
            double dbeta_on_cur_m10 = getUniformDistributionDensity(cur_m10, 0, 1);
            double m10_error_Prior_ratio = exp(log(dbeta_on_proposed_m10) - log(dbeta_on_cur_m10));

            if (verbose)
            {
                std::cout << "dbeta_on_proposed_m10: " << dbeta_on_proposed_m10 << "\n";
                std::cout << "dbeta_on_cur_m10: " << dbeta_on_cur_m10 << "\n";
                std::cout << "m10_error_Prior_ratio: " << m10_error_Prior_ratio << "\n";
            }

            double target_error_ratio = exp(proposed_loglikelihood + proposed_mutation_loglikelihood - cur_loglikelihood - cur_mutation_loglikehood 
                + cur_mutation_scaled_prior_loglikelihood - proposed_mutation_scaled_prior_loglikelihood) * m10_error_Prior_ratio;
            target_error_ratio = adjustInfinity(target_error_ratio);

            double r_error_ratio = getUniformDistribution(0, 1);
            if (r_error_ratio < target_error_ratio)
            {
                if (verbose)
                {
                    std::cout << "accept r_error_ratio: " << r_error_ratio << "\n";
                    std::cout << "accept target_error_ratio: " << target_error_ratio << "\n";
                }
                m10_vector.push_back(proposed_m10);
                loglikelihood_vector.push_back(proposed_loglikelihood);
                cur_mutation_loglikehood = proposed_mutation_loglikelihood;
                cur_mutation_scaled_prior_loglikelihood = proposed_mutation_scaled_prior_loglikelihood;
            }
            else
            {
                if (verbose)
                {
                    std::cout << "reject r_error_ratio: " << r_error_ratio << "\n";
                    std::cout << "reject target_error_ratio: " << target_error_ratio << "\n";
                }
                m10_vector.push_back(cur_m10);
                loglikelihood_vector.push_back(cur_loglikelihood);
                cur_tree_manip.buildNodesPossibilitiesInfo(cur_m10, m01);
            }
        }
        else
        {
            std::cout << "update tree.\n";

            double cur_alpha = alpha_vector[i - 1];
            double cur_beta = beta_vector[i - 1];
            double cur_m10 = m10_vector[i - 1];
            double cur_loglikelihood = loglikelihood_vector[i - 1];
            double cur_rootage = rootage_vector[i - 1]; //cur_tree->getTreeMaxHeight();
            updateErrorMatrix(cur_alpha, cur_beta);
            m10 = cur_m10;

            if (verbose)
                std::cout << "initial tree before update tree begin: " << cur_tree_manip.makeNewick(8) << "\n";

            alpha_vector.push_back(cur_alpha);
            beta_vector.push_back(cur_beta);
            m10_vector.push_back(cur_m10);

            double proposed_rootage = cur_rootage;
            TreeManip proposed_tree_tm;
            int random_sample = discrete_distribution(generator);
            double topology_proposal_ratio = 1.0, time_proposal_ratio = 1.0, rate_proposal_ratio = 1.0, topology_prior_ratio = 1.0, time_prior_ratio = 1.0, rate_prior_ratio = 1.0;
            update_tree(random_sample, cur_tree_manip, proposed_tree_tm, proposed_rootage, topology_prior_ratio, time_prior_ratio,
                        rate_prior_ratio, topology_proposal_ratio, time_proposal_ratio, rate_proposal_ratio, delta_time);

            if (verbose)
            {
                std::cout << "cur_tree_manip: " << cur_tree_manip.makeNewick(8) << std::endl;
                std::cout << "proposed_tree_tm: " << proposed_tree_tm.makeNewick(8) << std::endl;
                std::cout << "alpha: " << error_matrix[0][1] << std::endl;
                std::cout << "beta: " << error_matrix[1][0] << std::endl;
                std::cout << "tree m10: " << m10 << std::endl;
            }

            // proposed loglikelihood
            double proposed_loglikelihood = treeCompareFelsensteinDP(tree_data_original, data_label, proposed_tree_tm, state_code_and_count_map, likelihood_table, leaves_likelihood, m10);

            if (verbose)
            {
                std::cout << "cur_loglikelihood: " << cur_loglikelihood << "\n";
                std::cout << "proposed_loglikelihood: " << proposed_loglikelihood << "\n";

                std::cout << "topology_proposal_ratio: " << topology_proposal_ratio << "\n";
                std::cout << "time_proposal_ratio: " << time_proposal_ratio << "\n";
                std::cout << "rate_proposal_ratio: " << rate_proposal_ratio << "\n";
                std::cout << "topology_prior_ratio: " << topology_prior_ratio << "\n";
                std::cout << "time_prior_ratio: " << time_prior_ratio << "\n";
                std::cout << "rate_prior_ratio: " << rate_prior_ratio << "\n";
            }

            proposed_tree_tm.buildNodesPossibilitiesInfo(m10, m01);
            std::vector<int> proposed_mutation_branch_vector;
            std::vector<double> likelihoods = computeMutationPossibilities(proposed_tree_tm, proposed_mutation_branch_vector, possibility_matrix, possibility_value_vector, mutation_tree_data_original, data_label, leaves_likelihood);
            double proposed_mutation_loglikelihood = likelihoods[0];
            double proposed_mutation_scaled_prior_loglikelihood = likelihoods[1];
            std::cout << "proposed_mutation_loglikelihood: " << proposed_mutation_loglikelihood << "\n";
            std::cout << "proposed_mutation_scaled_prior_loglikelihood: " << proposed_mutation_scaled_prior_loglikelihood << "\n";
            std::cout << "cur_loglikelihood: " << cur_loglikelihood << "\n";
            std::cout << "cur_mutation_loglikelihood: " << cur_mutation_loglikehood << "\n";
            std::cout << "cur_mutationscaled_prior_loglikelihood: " << cur_mutation_scaled_prior_loglikelihood << "\n";

            double target_error_ratio = exp(proposed_loglikelihood + proposed_mutation_loglikelihood - cur_loglikelihood - cur_mutation_loglikehood + cur_mutation_scaled_prior_loglikelihood - proposed_mutation_scaled_prior_loglikelihood) *
                                        topology_proposal_ratio * time_proposal_ratio * rate_proposal_ratio * topology_prior_ratio * time_prior_ratio * rate_prior_ratio;
            target_error_ratio = adjustInfinity(target_error_ratio);

            double r_error_ratio = getUniformDistribution(0, 1);
            if (r_error_ratio < target_error_ratio)
            {
                if (verbose)
                {
                    std::cout << "accept r_error_ratio: " << r_error_ratio << "\n";
                    std::cout << "accept target_error_ratio: " << target_error_ratio << "\n";
                }
                loglikelihood_vector.push_back(proposed_loglikelihood);
                cur_tree_manip.clear();
                cur_tree_manip.buildFromNewick(proposed_tree_tm.makeNewick(8), true, false);
                cur_tree_manip.updateNodesHeightInfo();
                cur_tree_manip.buildNodeNameAndNumberMap();
                cur_tree_manip.buildNodesPossibilitiesInfo(m10, m01);

                TreeManip proposed_tree_tm_copy;
                proposed_tree_tm_copy.buildFromNewick(proposed_tree_tm.makeNewick(8), true, false);
                proposed_tree_tm_copy.addTToName();
                tree_string_vector.push_back(proposed_tree_tm_copy.makeNewick(8, true) + ";\n");

                rootage_vector.push_back(proposed_rootage);
                cur_mutation_branch_vector = proposed_mutation_branch_vector;
                cur_mutation_loglikehood = proposed_mutation_loglikelihood;
                cur_mutation_scaled_prior_loglikelihood = proposed_mutation_scaled_prior_loglikelihood;
            }
            else
            {
                if (verbose)
                {
                    std::cout << "reject r_error_ratio: " << r_error_ratio << "\n";
                    std::cout << "reject target_error_ratio: " << target_error_ratio << "\n";
                }
                loglikelihood_vector.push_back(cur_loglikelihood);

                TreeManip cur_tree_manip_copy;
                cur_tree_manip_copy.buildFromNewick(cur_tree_manip.makeNewick(8), true, false);
                cur_tree_manip_copy.addTToName();
                tree_string_vector.push_back(cur_tree_manip_copy.makeNewick(8, true) + ";\n");

                rootage_vector.push_back(cur_rootage);
            }
        }
    }
}

int main(int argc, char *argv[])
{
    // Get starting timepoint
    auto start = std::chrono::high_resolution_clock::now();

    generator.seed(time(0));

    // 0. read arguments
    // Read the arguments passing from command line and validate
    int opt, niter = 1000000;
    std::string output_filename, br_tree_string, tree_data_filename;

    while ((opt = getopt(argc, argv, "n:o:t:a:b:c:d:f:m:s:v")) != -1)
    {
        switch (opt)
        {
        case 'n':
            niter = atoi(optarg);
            std::cout << "num iter: " << niter << "\n";
            break;
        case 'o':
            output_filename = optarg;
            std::cout << "output_filename: " << output_filename << "\n";
            break;
        case 't':
            br_tree_string = optarg;
            std::cout << "start tree: " << br_tree_string << "\n";
            break;
        case 'a':
            bd_alpha_alpha = atof(optarg);
            std::cout << "bd_alpha for alpha: " << bd_alpha_alpha << "\n";
            break;
        case 'b':
            bd_beta_alpha = atof(optarg);
            std::cout << "bd_beta for alpha: " << bd_beta_alpha << "\n";
            break;
        case 'c':
            bd_alpha_beta = atof(optarg);
            std::cout << "bd_alpha for beta: " << bd_alpha_beta << "\n";
            break;
        case 'd':
            bd_beta_beta = atof(optarg);
            std::cout << "bd_beta for beta: " << bd_beta_beta << "\n";
            break;
        case 'f':
            tree_data_filename = optarg;
            std::cout << "tree_data_filename: " << tree_data_filename << "\n";
            break;
        case 'm':
            rootage_mean = atof(optarg);
            std::cout << "rootage mean: " << rootage_mean << "\n";
            break;
        case 's':
            rootage_stddev = atof(optarg);
            std::cout << "rootage stddev: " << rootage_stddev << "\n";
            break;
        case 'v':
            verbose = true;
            std::cout << "verbose on\n";
            break;
        default:
            fprintf(stderr, "Wrong arguments!\n");
            exit(1);
        }
    }

    std::cout << "done parsing cmd line \n";

    TreeManip br_tree;
    br_tree.buildFromNewick(br_tree_string, true, false);
    if (verbose)
        std::cout << "br_tree: " << br_tree.makeNewick(8) << std::endl;

    std::vector<std::string> data_label;
    std::vector<std::vector<int>> tree_data_original, mutation_tree_data_original;
    readObservedDataFromCSVFile(5000, tree_data_filename, data_label, tree_data_original, mutation_tree_data_original);

    for (auto l : data_label)
    {
        std::cout << "label: " << l << std::endl;
    }
    for (auto data : tree_data_original[0])
    {
        std::cout << "data: " << data << std::endl;
    }
    for (auto data : mutation_tree_data_original[0])
    {
        std::cout << "mu data: " << data << std::endl;
    }
    std::cout << "tree_data_original: " << tree_data_original.size() << std::endl;
    std::cout << "mutation_tree_data_original: " << mutation_tree_data_original.size() << std::endl;

    // 3. initialize value based on some distribution (normal, dirichlet)
    double alpha, beta, m01, m00, m10, m11, rootage;
    double lambda_edge, lambda_root;
    alpha = 0.1; //getUniformDistribution(0, 0.45);
    beta = 0.1; //getUniformDistribution(0, 0.45);
    m01 = 1, m00 = -m01, rootage = 1.045;
    m10 = 0.1; //getUniformDistribution(0, 1.0);
    m11 = -m10;
    lambda_edge = 2 * log(1.2), lambda_root = 2 * log(1.2);

    std::vector<double> alpha_vector, beta_vector, m10_vector, rootage_vector, loglikelihood_vector;
    std::vector<std::string> tree_string_vector;
    alpha_vector.reserve(niter);
    alpha_vector.push_back(alpha);
    beta_vector.reserve(niter);
    beta_vector.push_back(beta);
    m10_vector.reserve(niter);
    m10_vector.push_back(m10);
    rootage_vector.reserve(niter);
    rootage_vector.push_back(rootage);
    tree_string_vector.reserve(niter);

    // Two values may need to generate from R
    //std::string br_tree_string = "(3:0.9125,(5:0.8466827877,(1:0.01848638469,(4:0.01632424983,2:0.01632424983):0.002162134862):0.828196403):0.06581721228);";
    // TreeManip br_tree;
    // br_tree.buildFromNewick(br_tree_string, true, false);
    // std::cout << "br_tree: " << br_tree.makeNewick(8) << std::endl;

    //br_tree.scaleAllEdgeLengths(Tmax);
    //std::cout << "br_tree scaled: " << br_tree.makeNewick(8) << std::endl;
    // find map between tip node name and node number
    br_tree.buildNodeNameAndNumberMap();

    // 4. Run iterations with many iterations and four differnet conditions
    runTreeAndOrderSearchDP(tree_data_original, mutation_tree_data_original, data_label, br_tree, niter, alpha_vector, beta_vector, m10_vector, rootage_vector, loglikelihood_vector, tree_string_vector);

    std::cout << "average of alpha_vector: " << getAverageOfVector(alpha_vector) << "\n";
    std::cout << "average of beta_vector: " << getAverageOfVector(beta_vector) << "\n";
    std::cout << "average of m10_vector: " << getAverageOfVector(m10_vector) << "\n";
    std::cout << "max of m10_vector: " << getMaxOfVector(m10_vector) << "\n";
    std::cout << "min of m10_vector: " << getMinOfVector(m10_vector) << "\n";
    std::cout << "average of rootage_vector: " << getAverageOfVector(rootage_vector) << "\n";

    // 5. Output results to files
    outputToCSVFile(output_filename, alpha_vector, beta_vector, m10_vector, rootage_vector, loglikelihood_vector, tree_string_vector);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Time taken by function: " << duration.count() << " seconds\n";

    return 0;
}