#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <random>
#include <math.h>
#include <map>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/beta.hpp>
#include <boost/math/distributions/uniform.hpp>

namespace strom
{
    static double alpha, beta, m01, m00, m10, m11, rootage, Tmax;
    static double pi_error_alpha, pi_error_beta, pi_mutation_M10, delta_time;
    static double lambda_edge, lambda_root;
    static std::vector<double> alpha_vector, beta_vector, m10_vector, rootage_vector, tree_count;
    static std::vector<double> loglikelihood_vector;
    static std::vector<int> tree_id_vector;
    static std::vector<TreeManip> tree_manip_vector;
    static std::vector<std::vector<double>> leaves_likelihood, likelihood_table;
    static std::map<int, int> state_code_and_count_map;
    static double substitution_matrix[2][2], error_matrix[2][2];
    static double const STD_DEV = 0.005;
    static double const STD_DEV_M10 = 0.02;
    static double const BD_ALPHA = 1000 * 0.01;
    static double const BD_BETA = 1000 * (1 - 0.01);
    static std::discrete_distribution<int> discrete_distribution{0, 0, 0, 0, 1, 1, 1};
    static std::default_random_engine generator;

    static void readObservedDataFromCSVFile(const int num_lines, const std::string tree_data_filename, std::vector<std::string> &data_label, std::vector<std::vector<int>> &tree_data)
    {
        // read input tree file
        std::ifstream fin;
        fin.open(tree_data_filename);
        std::string label, line;
        getline(fin, label);
        // erase " and "t" from data
        label.erase(std::remove(label.begin(), label.end(), '"'), label.end());
        label.erase(std::remove(label.begin(), label.end(), 't'), label.end());

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
        data_label.push_back(label);

        int line_i = 0;
        std::vector<int> line_data;
        while (line_i < num_lines && getline(fin, line))
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
            line_data.push_back(line[0] - '0');
            tree_data.push_back(line_data);
            line_data.clear();
            line_i++;
        }

        fin.close();
    }

    static void outputToCSVFile(const std::string output_filename)
    {
        std::ofstream myfile;
        myfile.open(output_filename + ".csv");
        myfile << "alpha, beta, m10, rootage, loglikelihood, tree_id\n";
        for (int i = 0; i < alpha_vector.size(); i++)
        {
            myfile << alpha_vector[i] << ", ";
            myfile << beta_vector[i] << ", ";
            myfile << m10_vector[i] << ", ";
            myfile << rootage_vector[i] << ", ";
            myfile << loglikelihood_vector[i] << ", ";
            myfile << tree_id_vector[i] << "\n";
        }

        myfile.close();

        myfile.open(output_filename + "_trees.txt");
        for (int i = 0; i < tree_manip_vector.size(); i++)
        {
            tree_manip_vector[i].addTToName();
            myfile << tree_manip_vector[i].makeNewick(3, true) << ";\n";
        }
        myfile.close();
    }

    static double getNormalDistribution(double mean, double stddev)
    {
        std::normal_distribution<double> distribution(mean, stddev);
        return distribution(generator);
    }

    static double getNormalDistributionDensity(double x, double mean, double stddev)
    {
        boost::math::normal nd(mean, stddev);
        return boost::math::pdf(nd, x);
    }

    static double getBetaDistributionDensity(double x, double alpha, double beta)
    {
        boost::math::beta_distribution<> bd(alpha, beta);
        return boost::math::pdf(bd, x);
    }

    static double getUniformDistribution(double min, double max)
    {
        std::uniform_real_distribution<double> distribution(min, max);
        return distribution(generator);
    }

    static double getUniformDistributionDensity(double x, double min, double max)
    {
        boost::math::uniform_distribution<> ud(min, max);
        return boost::math::pdf(ud, x);
    }

    static void getSubstitutionMatrix(double br_len, double m10)
    {
        double divisor = m01 + m10;
        double exp_len = exp((-m01 - m10) * br_len);

        substitution_matrix[0][0] = (m10 + m01 * exp_len) / divisor;
        substitution_matrix[0][1] = (m01 - m01 * exp_len) / divisor;
        substitution_matrix[1][0] = (m10 - m10 * exp_len) / divisor;
        substitution_matrix[1][1] = (m01 + m10 * exp_len) / divisor;
    }

    static void getErrorMatrix(double alpha, double beta)
    {
        error_matrix[0][0] = 1 - alpha;
        error_matrix[0][1] = alpha;
        error_matrix[1][0] = beta;
        error_matrix[1][1] = 1 - beta;
    }

    static double getAverageOfVector(std::vector<double> &vector, int ignore_first_n = 0)
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

    static double getMaxOfVector(std::vector<double> &vector)
    {
        double max = -std::numeric_limits<double>::max();
        for (auto v : vector)
        {
            max = std::max(max, v);
        }
        return max;
    }

    static double getMinOfVector(std::vector<double> &vector)
    {
        double min = std::numeric_limits<double>::max();
        for (auto v : vector)
        {
            min = std::min(min, v);
        }
        return min;
    }

    static double adjustInfinity(double ratio)
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

    static int encodeStateToValue(std::vector<int> state)
    {
        int value = 0;
        for (auto s : state)
        {
            value <<= 1;
            value += s;
        }
        return value;
    }

    static void encodeValueToState(int value, int num_tips, std::vector<int> &state)
    {
        int i = 0;
        while (i < num_tips)
        {
            state.insert(state.begin(), value & 1);
            value >>= 1;
            i++;
        }
    }

    static void mapObservedStateToRightLabel(TreeManip &start_tree_tm, std::vector<int> &observed_state, std::vector<int> &mapped_observed_state, std::vector<std::string> &data_label)
    {
        int id = 0;
        for (auto state : observed_state)
        {
            mapped_observed_state[start_tree_tm.getNodeNumberByName(data_label[id])] = state;
            id++;
        }
    }

    static void buildCountsTableForTreeData(std::vector<std::vector<int>> &tree_data_original, std::map<int, int> &state_code_and_count_map)
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

    static void buildLeaveLikelihoodMatrix(std::vector<int> &observed_state)
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

    static double felsensteinBinaryDP(std::vector<int> &observed_state, Node *cur_node, int cur_state, double cur_m10)
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
            getSubstitutionMatrix(left_child->getEdgeLength(), cur_m10);
            double pi0_left_child = substitution_matrix[cur_state][0];
            double pi1_left_child = substitution_matrix[cur_state][1];
            getSubstitutionMatrix(right_child->getEdgeLength(), cur_m10);
            double pi0_right_child = substitution_matrix[cur_state][0];
            double pi1_right_child = substitution_matrix[cur_state][1];

            double likelihood_left_child_0 = felsensteinBinaryDP(observed_state, left_child, 0, cur_m10);
            double likelihood_left_child_1 = felsensteinBinaryDP(observed_state, left_child, 1, cur_m10);
            double likelihood_right_child_0 = felsensteinBinaryDP(observed_state, right_child, 0, cur_m10);
            double likelihood_right_child_1 = felsensteinBinaryDP(observed_state, right_child, 1, cur_m10);

            double likelihood = (pi0_left_child * likelihood_left_child_0 + pi1_left_child * likelihood_left_child_1) *
                                (pi0_right_child * likelihood_right_child_0 + pi1_right_child * likelihood_right_child_1);
            likelihood_table[cur_state][cur_node->getNumber()] = likelihood;
            return likelihood;
        }
    }

    static double treeCompareFelsensteinDP(std::vector<std::vector<int>> &tree_data_original, std::vector<std::string> &data_label, TreeManip &start_tree_tm, double cur_m10)
    {
        // initialize values
        Tree::SharedPtr start_tree = start_tree_tm.getTree();
        int num_tips = start_tree->numLeaves();
        //int num_edges = 2 * num_tips - 2;
        //int num_nodes = start_tree->numNodes();

        std::cout << "DP cur_m10: " << cur_m10 << "\n";

        Node *root = start_tree->getRoot()->getLeftChild();
        double total_loglikelihood = 0.0;
        for (auto it = state_code_and_count_map.begin(); it != state_code_and_count_map.end(); it++)
        {
            // transfer value to state vector
            int value = it->first;
            int count = it->second;
            // std::cout << "value " << value << "\n";
            // std::cout << "count " << count << "\n";
            std::vector<int> observed_state;
            observed_state.reserve(num_tips);
            encodeValueToState(value, num_tips, observed_state);
            std::vector<int> mapped_observed_state(observed_state);
            mapObservedStateToRightLabel(start_tree_tm, observed_state, mapped_observed_state, data_label);

            buildLeaveLikelihoodMatrix(mapped_observed_state);
            double likelihood = felsensteinBinaryDP(mapped_observed_state, root, 0, cur_m10);
            total_loglikelihood += log(likelihood) * count;
            std::fill(likelihood_table[0].begin(), likelihood_table[0].end(), 0);
            std::fill(likelihood_table[1].begin(), likelihood_table[1].end(), 0);
        }
        return total_loglikelihood;
    }

    static void update_tree(int random_sample, TreeManip &cur_tree_tm, TreeManip &proposed_tree_tm, double &proposed_rootage, double &topology_prior_ratio, double &time_prior_ratio,
                            double &rate_prior_ratio, double &topology_proposal_ratio, double &time_proposal_ratio, double &rate_proposal_ratio)
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
        {
            std::cout << "wide exchange.\n";
            proposed_tree_tm.wideExchangeFrom(cur_tree_tm);
            break;
        }
        case 2:
        {
            std::cout << "FNPR exchange.\n";
            proposed_tree_tm.FNPRExchangeFrom(cur_tree_tm);
            break;
        }
        case 3:
        {
            std::cout << "Total length change.\n";
            proposed_tree_tm.buildFromNewick(cur_tree_tm.makeNewick(8), true, false);
            std::cout << "initial propose tree: " << proposed_tree_tm.makeNewick(3) << "\n";
            proposed_tree_tm.totalLengthChange(proposed_rootage, rate_prior_ratio, rate_proposal_ratio);
            std::cout << "total length changed propose tree: " << proposed_tree_tm.makeNewick(3) << "\n";
            break;
        }
        case 4:
        {
            std::cout << "Random length change.\n";
            proposed_tree_tm.buildFromNewick(cur_tree_tm.makeNewick(8), true, false);
            proposed_tree_tm.getTree()->updateNodesHeightInfo();
            std::cout << "initial propose tree: " << proposed_tree_tm.makeNewick(3) << "\n";

            proposed_tree_tm.getTree()->updateNodesHeightInfo();
            Node::PtrVector internal_nodes = proposed_tree_tm.getTree()->getAllInternals();
            std::random_shuffle(internal_nodes.begin(), internal_nodes.end());
            Node *internal = internal_nodes.front();

            // get parent and child edge height
            double parent_edge_height = internal->getHeight();
            double parent_edge_length = internal->getEdgeLength();
            double child_edge_height = internal->getLeftChild()->getHeight();
            double child_edge_length = internal->getLeftChild()->getEdgeLength();
            bool find_right_child = false;
            if (child_edge_length > internal->getLeftChild()->getRightSib()->getEdgeLength())
            {
                child_edge_length = internal->getLeftChild()->getRightSib()->getEdgeLength();
                find_right_child = true;
            }
            double std_dev = std::min(parent_edge_length, child_edge_length) * delta_time / 2.0;
            double max_height = proposed_tree_tm.getTree()->getTreeMaxHeight();
            double proposed_edge_height = max_height - getNormalDistribution((max_height - parent_edge_height), std_dev);
            internal->setEdgeLength(proposed_edge_height - internal->getParent()->getHeight());
            internal->getLeftChild()->setEdgeLength(parent_edge_length + internal->getLeftChild()->getEdgeLength() - internal->getEdgeLength());
            internal->getLeftChild()->getRightSib()->setEdgeLength(parent_edge_length + internal->getLeftChild()->getRightSib()->getEdgeLength() - internal->getEdgeLength());
            proposed_tree_tm.getTree()->updateNodesHeightInfo();

            std::cout << "delta_time " << delta_time << "\n";
            std::cout << "proposed_edge_height " << proposed_edge_height << "\n";
            std::cout << "parent_edge_height " << parent_edge_height << "\n";

            double proposed_parent_br_len = internal->getEdgeLength();
            double proposed_child_edge_length;
            if (find_right_child)
            {
                proposed_child_edge_length = internal->getLeftChild()->getRightSib()->getEdgeLength();
            }
            else
            {
                proposed_child_edge_length = internal->getLeftChild()->getEdgeLength();
            }

            double cur_on_proposed_normal_density = getNormalDistributionDensity(max_height - parent_edge_height, max_height - proposed_edge_height,
                                                                                 std::min(proposed_parent_br_len, proposed_child_edge_length) * delta_time / 2.0);
            double proposed_on_cur_normal_density = getNormalDistributionDensity(max_height - proposed_edge_height, max_height - parent_edge_height, std_dev);
            time_proposal_ratio *= exp(log(cur_on_proposed_normal_density) - log(proposed_on_cur_normal_density));

            std::cout << "random length changed propose tree: " << proposed_tree_tm.makeNewick(3) << "\n";
            proposed_tree_tm.getTree()->updateNodesHeightInfo();
            break;
        }
        case 5:
        {
            std::cout << "All internal length change.\n";
            proposed_tree_tm.buildFromNewick(cur_tree_tm.makeNewick(8), true, false);
            proposed_tree_tm.getTree()->updateNodesHeightInfo();
            proposed_tree_tm.allInternalLengthChange(delta_time, time_proposal_ratio);
            proposed_tree_tm.getTree()->updateNodesHeightInfo();
            break;
        }
        case 6:
        {
            std::cout << "Root change only.\n";
            proposed_tree_tm.buildFromNewick(cur_tree_tm.makeNewick(8), true, false);
            proposed_tree_tm.getTree()->updateNodesHeightInfo();
            std::cout<< "root changed initial tree: " << proposed_tree_tm.makeNewick(3) << "\n";

            double edge_u = getUniformDistribution(0.0, 1.0);
            double edge_m = exp(lambda_root * (edge_u - 0.5));
            double cur_rootage = proposed_rootage;

            Node *root = proposed_tree_tm.getTree()->getRoot()->getLeftChild();
            double child_edge_length = root->getLeftChild()->getEdgeLength();
            bool find_right_child = false;
            if (child_edge_length > root->getLeftChild()->getRightSib()->getEdgeLength())
            {
                child_edge_length = root->getLeftChild()->getRightSib()->getEdgeLength();
                find_right_child = true;
            }
            double max_height = proposed_tree_tm.getTree()->getTreeMaxHeight();
            double proposed_root_height = child_edge_length * edge_m;
            double delta_root_height = proposed_root_height - child_edge_length;
            root->setEdgeLength(0.0);
            root->getLeftChild()->setEdgeLength(delta_root_height + root->getLeftChild()->getEdgeLength());
            root->getLeftChild()->getRightSib()->setEdgeLength(delta_root_height + root->getLeftChild()->getRightSib()->getEdgeLength());

            proposed_tree_tm.getTree()->updateNodesHeightInfo();
            proposed_rootage = proposed_tree_tm.getTree()->getTreeMaxHeight();
            std::cout<<"child_edge_length " << child_edge_length <<"\n";
            std::cout<<"proposed_root_height " << proposed_root_height <<"\n";
            std::cout<<"delta_root_height " << delta_root_height <<"\n";
            std::cout<< "root changed proposed tree: " << proposed_tree_tm.makeNewick(3) << "\n";
            for (auto nd : proposed_tree_tm.getTree()->getPreOrder())
            {
                std::cout<<"nd num " << nd->getNumber() << ", height " << nd->getHeight() <<"\n";
            }
            std::cout<<"proposed_rootage " << proposed_rootage <<"\n";
            std::cout<<"cur_rootage " << cur_rootage <<"\n";

            double proposed_rootage_dnorm = getNormalDistributionDensity(proposed_rootage, 0.9125, 0.2);
            double cur_rootage_dnorm = getNormalDistributionDensity(cur_rootage, 0.9125, 0.2);
            rate_prior_ratio = proposed_rootage_dnorm / cur_rootage_dnorm;
            rate_proposal_ratio = proposed_rootage / cur_rootage;

            break;
        }
        default:
            std::cout << "Error: Now we only support 7 different cases.\n";
        }
        proposed_tree_tm.buildNodeNameAndNumberMap();
        // proposed_rootage = proposed_tree_tm.getTree()->getTreeMaxHeight();
        std::cout << "proposed_rootage end " << proposed_rootage << "\n";
    }

    static void runTreeAndOrderSearchDP(std::vector<std::vector<int>> &tree_data_original, std::vector<std::string> &data_label, TreeManip &true_tree_tm, TreeManip &start_tree_tm, int niter)
    {
        // initialize all values
        Tree::SharedPtr start_tree = start_tree_tm.getTree();
        int num_tips = start_tree->numLeaves();
        start_tree->updateNodesHeightInfo();
        start_tree_tm.buildNodeNameAndNumberMap();
        //int num_edges = 2 * num_tips - 2;
        int num_nodes = start_tree->numNodes();
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
        buildCountsTableForTreeData(tree_data_original, state_code_and_count_map);

        getErrorMatrix(alpha, beta);

        // run felsensteinBinaryDP for tree
        double total_loglikelihood = treeCompareFelsensteinDP(tree_data_original, data_label, start_tree_tm, m10);
        std::cout << "total_loglikelihood: " << total_loglikelihood << "\n";

        // initialize vectors for iterations
        int tree_id = 0;
        tree_id_vector.reserve(niter);
        tree_id_vector.push_back(tree_id);
        loglikelihood_vector.reserve(niter);
        loglikelihood_vector.push_back(total_loglikelihood);
        tree_manip_vector.reserve(niter);
        tree_manip_vector.push_back(start_tree_tm);

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
                getErrorMatrix(cur_alpha, cur_beta);
                m10 = cur_m10;

                rootage_vector.push_back(cur_rootage);
                m10_vector.push_back(cur_m10);
                tree_id_vector.push_back(tree_id);
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
                getErrorMatrix(proposed_alpha, cur_beta);
                double proposed_loglikelihood = treeCompareFelsensteinDP(tree_data_original, data_label, cur_tree_manip, cur_m10);

                // std::cout<<"cur_alpha: "<< cur_alpha<<"\n";
                // std::cout<<"proposed_alpha: "<< proposed_alpha<<"\n";
                // std::cout<<"proposed_loglikelihood: "<< proposed_loglikelihood<<"\n";
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
                getErrorMatrix(cur_alpha, cur_beta);
                m10 = cur_m10;

                rootage_vector.push_back(cur_rootage);
                m10_vector.push_back(cur_m10);
                tree_id_vector.push_back(tree_id);
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
                getErrorMatrix(cur_alpha, proposed_beta);
                double proposed_loglikelihood = treeCompareFelsensteinDP(tree_data_original, data_label, cur_tree_manip, cur_m10);

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
                getErrorMatrix(cur_alpha, cur_beta);
                m10 = cur_m10;

                std::cout << "initial tree" << cur_tree_manip.makeNewick(3) << "\n";
                cur_tree_manip.getTree()->updateNodesHeightInfo();
                cur_tree_manip.buildNodeNameAndNumberMap();

                rootage_vector.push_back(cur_rootage);
                alpha_vector.push_back(cur_alpha);
                beta_vector.push_back(cur_beta);
                tree_id_vector.push_back(tree_id);
                tree_manip_vector.push_back(cur_tree_manip);

                double proposed_m10 = std::abs(getNormalDistribution(cur_m10, STD_DEV_M10));
                std::cout << "cur_m10: " << cur_m10 << "\n";
                std::cout << "proposed_m10: " << proposed_m10 << "\n";
                std::cout << "cur_tree_manip: " << cur_tree_manip.makeNewick(3) << std::endl;

                // proposed loglikelihood
                m10 = proposed_m10;
                double proposed_loglikelihood = treeCompareFelsensteinDP(tree_data_original, data_label, cur_tree_manip, proposed_m10);
                std::cout << "cur_loglikelihood: " << cur_loglikelihood << "\n";
                std::cout << "proposed_loglikelihood: " << proposed_loglikelihood << "\n";

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
                getErrorMatrix(cur_alpha, cur_beta);
                m10 = cur_m10;

                alpha_vector.push_back(cur_alpha);
                beta_vector.push_back(cur_beta);
                m10_vector.push_back(cur_m10);

                double proposed_rootage = cur_rootage;
                TreeManip proposed_tree_tm;
                int random_sample = discrete_distribution(generator);
                double topology_prior_ratio = 1.0, time_prior_ratio = 1.0, rate_prior_ratio = 1.0;
                double topology_proposal_ratio = 1.0, time_proposal_ratio = 1.0, rate_proposal_ratio = 1.0;
                update_tree(random_sample, cur_tree_manip, proposed_tree_tm, proposed_rootage, topology_prior_ratio, time_prior_ratio, rate_prior_ratio, topology_proposal_ratio, time_proposal_ratio, rate_proposal_ratio);

                std::cout << "cur_rootage: " << cur_rootage << "\n";
                std::cout << "proposed_rootage: " << proposed_rootage << "\n";

                // std::cout << "cur_tree_manip: " << cur_tree_manip.makeNewick(3) << std::endl;
                // std::cout << "proposed_tree_tm: " << proposed_tree_tm.makeNewick(3) << std::endl;
                // std::cout << "alpha: " << error_matrix[0][1] << std::endl;
                // std::cout << "beta: " << error_matrix[1][0] << std::endl;
                std::cout << "tree m10: " << m10 << std::endl;

                // proposed loglikelihood
                double proposed_loglikelihood = treeCompareFelsensteinDP(tree_data_original, data_label, proposed_tree_tm, cur_m10);
                std::cout << "cur_loglikelihood: " << cur_loglikelihood << "\n";
                std::cout << "proposed_loglikelihood: " << proposed_loglikelihood << "\n";

                std::cout << "topology_proposal_ratio: " << topology_proposal_ratio << "\n";
                std::cout << "time_proposal_ratio: " << time_proposal_ratio << "\n";
                std::cout << "rate_proposal_ratio: " << rate_proposal_ratio << "\n";
                std::cout << "topology_prior_ratio: " << topology_prior_ratio << "\n";
                std::cout << "time_prior_ratio: " << time_prior_ratio << "\n";
                std::cout << "rate_prior_ratio: " << rate_prior_ratio << "\n";

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
                    tree_id++;
                    tree_id_vector.push_back(tree_id);
                    rootage_vector.push_back(proposed_rootage);
                    std::cout << "use new tree.\n";
                }
                else
                {
                    std::cout << "reject r_error_ratio: " << r_error_ratio << "\n";
                    std::cout << "reject target_error_ratio: " << target_error_ratio << "\n";
                    loglikelihood_vector.push_back(cur_loglikelihood);
                    tree_manip_vector.push_back(cur_tree_manip);
                    tree_id_vector.push_back(tree_id);
                    rootage_vector.push_back(cur_rootage);
                    std::cout << "use old tree.\n";
                }
            }
            std::cout << "rootage_vector.back() " << rootage_vector.back() << "\n";
        }
    }

} // namespace strom