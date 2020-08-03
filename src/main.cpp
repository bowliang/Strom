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
double substitution_matrix[2][2], error_matrix[2][2];

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
        //std::cout << "value: " << value << ", count: " << count << "\n";
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
        proposed_tree_tm.totalLengthChange(proposed_rootage, rate_prior_ratio, rate_proposal_ratio);
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

void runTreeAndOrderSearchDP(std::vector<std::vector<int>> &tree_data_original, std::vector<std::string> &data_label, TreeManip &start_tree_tm, int niter, std::vector<double> &alpha_vector, std::vector<double> &beta_vector,
                             std::vector<double> &m10_vector, std::vector<double> &rootage_vector, std::vector<double> &loglikelihood_vector, std::vector<TreeManip> &tree_manip_vector)
{
    // initialize all values
    int num_tips = start_tree_tm.getNumLeaves();
    start_tree_tm.updateNodesHeightInfo();
    int num_nodes = start_tree_tm.getNumNodes();

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
    std::cout << "start_tree_tm: " << start_tree_tm.makeNewick(3) << std::endl;

    // build counts table for all possible observed data
    std::map<int, int> state_code_and_count_map;
    buildCountsTableForTreeData(tree_data_original, state_code_and_count_map);

    double alpha = alpha_vector.front();
    double beta = beta_vector.front();
    updateErrorMatrix(alpha, beta);
    double m10 = m10_vector.front();

    // run felsensteinBinaryDP for tree
    double total_loglikelihood = treeCompareFelsensteinDP(tree_data_original, data_label, start_tree_tm, state_code_and_count_map, likelihood_table, leaves_likelihood, m10);
    std::cout << "total_loglikelihood: " << total_loglikelihood << "\n";

    // initialize vectors for iterations
    tree_manip_vector.push_back(start_tree_tm);
    loglikelihood_vector.reserve(niter);
    loglikelihood_vector.push_back(total_loglikelihood);
    double pi_error_alpha = 0.1, pi_error_beta = 0.2, pi_mutation_M10 = 0.3, delta_time = 0.2;

    // start component-wise update
    for (int i = 1; i < niter; i++)
    {
        double r = getUniformDistribution(0, 1);
        std::cout << "iter: " << i << "\n";
        if (r < pi_error_alpha)
        {
            std::cout << "update alpha.\n";

            TreeManip cur_tree_manip = tree_manip_vector[i - 1];
            double cur_alpha = alpha_vector[i - 1];
            double cur_beta = beta_vector[i - 1];
            double cur_m10 = m10_vector[i - 1];
            double cur_loglikelihood = loglikelihood_vector[i - 1];
            double cur_rootage = rootage_vector[i - 1];
            updateErrorMatrix(cur_alpha, cur_beta);
            m10 = cur_m10;

            std::cout << "initial tree before update alpha: " << cur_tree_manip.makeNewick(3) << "\n";

            rootage_vector.push_back(cur_rootage);
            m10_vector.push_back(cur_m10);
            tree_manip_vector.push_back(cur_tree_manip);

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

            // std::cout<<"cur_alpha: "<< cur_alpha<<"\n";
            // std::cout<<"proposed_alpha: "<< proposed_alpha<<"\n";
            std::cout<<"cur_loglikelihood: "<< cur_loglikelihood<<"\n";
            std::cout<<"proposed_loglikelihood: "<< proposed_loglikelihood<<"\n";
            double dnorm_on_cur_alpha = getNormalDistributionDensity(cur_alpha, proposed_alpha, STD_DEV);
            double dnorm_on_proposed_alpha = getNormalDistributionDensity(proposed_alpha, cur_alpha, STD_DEV);
            double alpha_error_proposal_ratio = exp(log(dnorm_on_cur_alpha) - log(dnorm_on_proposed_alpha));

            double dbeta_on_proposed_alpha = getBetaDistributionDensity(proposed_alpha, BD_ALPHA, BD_BETA);
            double dbeta_on_cur_alpha = getBetaDistributionDensity(cur_alpha, BD_ALPHA, BD_BETA);
            double alpha_error_Prior_ratio = exp(log(dbeta_on_proposed_alpha) - log(dbeta_on_cur_alpha));

            double target_error_ratio = exp(proposed_loglikelihood - cur_loglikelihood) * alpha_error_proposal_ratio * alpha_error_Prior_ratio;
            target_error_ratio = adjustInfinity(target_error_ratio);

            double r_error_ratio = getUniformDistribution(0, 1);
            if (r_error_ratio < target_error_ratio)
            {
                std::cout << "accept r_error_ratio: " << r_error_ratio << "\n";
                std::cout << "accept target_error_ratio: " << target_error_ratio << "\n";
                alpha_vector.push_back(proposed_alpha);
                beta_vector.push_back(cur_beta);
                loglikelihood_vector.push_back(proposed_loglikelihood);
            }
            else
            {
                std::cout << "reject r_error_ratio: " << r_error_ratio << "\n";
                std::cout << "reject target_error_ratio: " << target_error_ratio << "\n";
                alpha_vector.push_back(cur_alpha);
                beta_vector.push_back(cur_beta);
                loglikelihood_vector.push_back(cur_loglikelihood);
            }
        }
        else if (r < pi_error_beta)
        {
            std::cout << "update beta.\n";

            TreeManip cur_tree_manip = tree_manip_vector[i - 1];
            double cur_alpha = alpha_vector[i - 1];
            double cur_beta = beta_vector[i - 1];
            double cur_m10 = m10_vector[i - 1];
            double cur_loglikelihood = loglikelihood_vector[i - 1];
            double cur_rootage = rootage_vector[i - 1];
            updateErrorMatrix(cur_alpha, cur_beta);
            m10 = cur_m10;

            std::cout << "initial tree before update beta: " << cur_tree_manip.makeNewick(3) << "\n";

            rootage_vector.push_back(cur_rootage);
            m10_vector.push_back(cur_m10);
            tree_manip_vector.push_back(cur_tree_manip);

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

            // proposed loglikelihood
            updateErrorMatrix(cur_alpha, proposed_beta);
            double proposed_loglikelihood = treeCompareFelsensteinDP(tree_data_original, data_label, cur_tree_manip, state_code_and_count_map, likelihood_table, leaves_likelihood, m10);

            std::cout<<"cur_loglikelihood: "<< cur_loglikelihood<<"\n";
            std::cout<<"proposed_loglikelihood: "<< proposed_loglikelihood<<"\n";

            double dnorm_on_cur_beta = getNormalDistributionDensity(cur_beta, proposed_beta, STD_DEV);
            double dnorm_on_proposed_beta = getNormalDistributionDensity(proposed_beta, cur_beta, STD_DEV);
            double beta_error_proposal_ratio = exp(log(dnorm_on_cur_beta) - log(dnorm_on_proposed_beta));

            double dbeta_on_proposed_beta = getBetaDistributionDensity(proposed_beta, BD_ALPHA, BD_BETA);
            double dbeta_on_cur_beta = getBetaDistributionDensity(cur_beta, BD_ALPHA, BD_BETA);
            double beta_error_Prior_ratio = exp(log(dbeta_on_proposed_beta) - log(dbeta_on_cur_beta));

            double target_error_ratio = exp(proposed_loglikelihood - cur_loglikelihood) * beta_error_proposal_ratio * beta_error_Prior_ratio;
            target_error_ratio = adjustInfinity(target_error_ratio);

            double r_error_ratio = getUniformDistribution(0, 1);
            if (r_error_ratio < target_error_ratio)
            {
                std::cout << "accept r_error_ratio: " << r_error_ratio << "\n";
                std::cout << "accept target_error_ratio: " << target_error_ratio << "\n";
                alpha_vector.push_back(cur_alpha);
                beta_vector.push_back(proposed_beta);
                loglikelihood_vector.push_back(proposed_loglikelihood);
            }
            else
            {
                std::cout << "reject r_error_ratio: " << r_error_ratio << "\n";
                std::cout << "reject target_error_ratio: " << target_error_ratio << "\n";
                alpha_vector.push_back(cur_alpha);
                beta_vector.push_back(cur_beta);
                loglikelihood_vector.push_back(cur_loglikelihood);
            }
        }
        else if (r < pi_mutation_M10)
        {
            std::cout << "update M10.\n";

            TreeManip cur_tree_manip = tree_manip_vector[i - 1];
            double cur_alpha = alpha_vector[i - 1];
            double cur_beta = beta_vector[i - 1];
            double cur_m10 = m10_vector[i - 1];
            double cur_loglikelihood = loglikelihood_vector[i - 1];
            double cur_rootage = rootage_vector[i - 1];
            updateErrorMatrix(cur_alpha, cur_beta);
            m10 = cur_m10;

            std::cout << "initial tree before update m10: " << cur_tree_manip.makeNewick(3) << "\n";

            rootage_vector.push_back(cur_rootage);
            alpha_vector.push_back(cur_alpha);
            beta_vector.push_back(cur_beta);
            tree_manip_vector.push_back(cur_tree_manip);

            double proposed_m10 = std::abs(getNormalDistribution(cur_m10, STD_DEV_M10));
            std::cout<<"cur_m10: "<< cur_m10<<"\n";
            std::cout<<"proposed_m10: "<< proposed_m10<<"\n";
            // std::cout << "cur_tree_manip: " << cur_tree_manip.makeNewick(3) << std::endl;

            // proposed loglikelihood
            double proposed_loglikelihood = treeCompareFelsensteinDP(tree_data_original, data_label, cur_tree_manip, state_code_and_count_map, likelihood_table, leaves_likelihood, proposed_m10);

            std::cout<<"cur_loglikelihood: "<< cur_loglikelihood<<"\n";
            std::cout<<"proposed_loglikelihood: "<< proposed_loglikelihood<<"\n";

            double dbeta_on_proposed_m10 = getUniformDistributionDensity(proposed_m10, 0, 1);
            double dbeta_on_cur_m10 = getUniformDistributionDensity(cur_m10, 0, 1);
            double m10_error_Prior_ratio = exp(log(dbeta_on_proposed_m10) - log(dbeta_on_cur_m10));

            double target_error_ratio = exp(proposed_loglikelihood - cur_loglikelihood) * m10_error_Prior_ratio;
            target_error_ratio = adjustInfinity(target_error_ratio);

            double r_error_ratio = getUniformDistribution(0, 1);
            if (r_error_ratio < target_error_ratio)
            {
                std::cout << "accept r_error_ratio: " << r_error_ratio << "\n";
                std::cout << "accept target_error_ratio: " << target_error_ratio << "\n";
                m10_vector.push_back(proposed_m10);
                loglikelihood_vector.push_back(proposed_loglikelihood);
            }
            else
            {
                std::cout << "reject r_error_ratio: " << r_error_ratio << "\n";
                std::cout << "reject target_error_ratio: " << target_error_ratio << "\n";
                m10_vector.push_back(cur_m10);
                loglikelihood_vector.push_back(cur_loglikelihood);
            }
            std::cout << "m10_vector.back(): " << m10_vector.back() << "\n";
        }
        else
        {
            std::cout << "update tree.\n";

            TreeManip cur_tree_manip = tree_manip_vector[i - 1];
            double cur_alpha = alpha_vector[i - 1];
            double cur_beta = beta_vector[i - 1];
            double cur_m10 = m10_vector[i - 1];
            double cur_loglikelihood = loglikelihood_vector[i - 1];
            double cur_rootage = rootage_vector[i - 1]; //cur_tree->getTreeMaxHeight();
            updateErrorMatrix(cur_alpha, cur_beta);
            m10 = cur_m10;

            std::cout << "initial tree before update tree begin: " << cur_tree_manip.makeNewick(3) << "\n";

            alpha_vector.push_back(cur_alpha);
            beta_vector.push_back(cur_beta);
            m10_vector.push_back(cur_m10);

            double proposed_rootage = cur_rootage;
            TreeManip proposed_tree_tm;
            int random_sample = discrete_distribution(generator);
            double topology_proposal_ratio = 1.0, time_proposal_ratio = 1.0, rate_proposal_ratio = 1.0, topology_prior_ratio = 1.0, time_prior_ratio = 1.0, rate_prior_ratio = 1.0;
            update_tree(random_sample, cur_tree_manip, proposed_tree_tm, proposed_rootage, topology_prior_ratio, time_prior_ratio, 
                rate_prior_ratio, topology_proposal_ratio, time_proposal_ratio, rate_proposal_ratio, delta_time);

            std::cout << "cur_tree_manip: " << cur_tree_manip.makeNewick(3) << std::endl;
            std::cout << "proposed_tree_tm: " << proposed_tree_tm.makeNewick(3) << std::endl;
            std::cout << "alpha: " << error_matrix[0][1] << std::endl;
            std::cout << "beta: " << error_matrix[1][0] << std::endl;
            std::cout << "tree m10: " << m10 << std::endl;

            // proposed loglikelihood
            double proposed_loglikelihood = treeCompareFelsensteinDP(tree_data_original, data_label, proposed_tree_tm, state_code_and_count_map, likelihood_table, leaves_likelihood, m10);            
            
            std::cout<<"cur_loglikelihood: "<< cur_loglikelihood<<"\n";
            std::cout<<"proposed_loglikelihood: "<< proposed_loglikelihood<<"\n";

            std::cout<<"topology_proposal_ratio: "<< topology_proposal_ratio<<"\n";
            std::cout<<"time_proposal_ratio: "<< time_proposal_ratio<<"\n";
            std::cout<<"rate_proposal_ratio: "<< rate_proposal_ratio<<"\n";
            std::cout<<"topology_prior_ratio: "<< topology_prior_ratio<<"\n";
            std::cout<<"time_prior_ratio: "<< time_prior_ratio<<"\n";
            std::cout<<"rate_prior_ratio: "<< rate_prior_ratio<<"\n";

            double target_error_ratio = exp(proposed_loglikelihood - cur_loglikelihood) *
                                        topology_proposal_ratio * time_proposal_ratio * rate_proposal_ratio * topology_prior_ratio * time_prior_ratio * rate_prior_ratio;
            target_error_ratio = adjustInfinity(target_error_ratio);

            double r_error_ratio = getUniformDistribution(0, 1);
            if (r_error_ratio < target_error_ratio)
            {
                std::cout << "accept r_error_ratio: " << r_error_ratio << "\n";
                std::cout << "accept target_error_ratio: " << target_error_ratio << "\n";
                loglikelihood_vector.push_back(proposed_loglikelihood);
                tree_manip_vector.push_back(proposed_tree_tm);
                rootage_vector.push_back(proposed_rootage);
                std::cout << "use new tree.\n";
            }
            else
            {
                std::cout << "reject r_error_ratio: " << r_error_ratio << "\n";
                std::cout << "reject target_error_ratio: " << target_error_ratio << "\n";
                loglikelihood_vector.push_back(cur_loglikelihood);
                tree_manip_vector.push_back(cur_tree_manip);
                rootage_vector.push_back(cur_rootage);
                std::cout << "use old tree.\n";
            }
        }
        std::cout << "rootage_vector.back() " << rootage_vector.back() << "\n";
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
            std::cout << "output_filename: " << output_filename << "\n";
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

    // 3. initialize value based on some distribution (normal, dirichlet)
    double alpha, beta, m01, m00, m10, m11, rootage;
    double lambda_edge, lambda_root;
    alpha = getUniformDistribution(0, 0.45);
    beta = getUniformDistribution(0, 0.45);
    m01 = 1, m00 = -m01, rootage = 1.045;
    m10 = getUniformDistribution(0, 1.0);
    m11 = -m10;
    lambda_edge = 2 * log(1.2), lambda_root = 2 * log(1.2);

    std::vector<double> alpha_vector, beta_vector, m10_vector, rootage_vector, loglikelihood_vector;
    std::vector<TreeManip> tree_manip_vector;
    alpha_vector.reserve(niter);
    alpha_vector.push_back(alpha);
    beta_vector.reserve(niter);
    beta_vector.push_back(beta);
    m10_vector.reserve(niter);
    m10_vector.push_back(m10);
    rootage_vector.reserve(niter);
    rootage_vector.push_back(rootage);
    tree_manip_vector.reserve(niter);

    // Two values may need to generate from R
    std::string br_tree_string = "((2:0.4994737979,3:0.4994737979):0.5293488252,((4:0.4210867618,1:0.4210867618):0.1298891809,5:0.5509759427):0.4778466805);";
    TreeManip br_tree;
    br_tree.buildFromNewick(br_tree_string, true, false);
    std::cout << "br_tree: " << br_tree.makeNewick(3) << std::endl;

    //br_tree.scaleAllEdgeLengths(Tmax);
    std::cout << "br_tree scaled: " << br_tree.makeNewick(3) << std::endl;
    // find map between tip node name and node number
    br_tree.buildNodeNameAndNumberMap();

    // 4. Run iterations with many iterations and four differnet conditions
    runTreeAndOrderSearchDP(tree_data_original, data_label, br_tree, niter, alpha_vector, beta_vector, m10_vector, rootage_vector, loglikelihood_vector, tree_manip_vector);

    std::cout << "average of alpha_vector: " << getAverageOfVector(alpha_vector) << "\n";
    std::cout << "average of beta_vector: " << getAverageOfVector(beta_vector) << "\n";
    std::cout << "average of m10_vector: " << getAverageOfVector(m10_vector) << "\n";
    std::cout << "max of m10_vector: " << getMaxOfVector(m10_vector) << "\n";
    std::cout << "min of m10_vector: " << getMinOfVector(m10_vector) << "\n";
    std::cout << "average of rootage_vector: " << getAverageOfVector(rootage_vector) << "\n";

    // 5. Output results to files
    outputToCSVFile(output_filename, alpha_vector, beta_vector, m10_vector, rootage_vector, loglikelihood_vector, tree_manip_vector);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Time taken by function: " << duration.count() << " seconds\n";

    return 0;
}