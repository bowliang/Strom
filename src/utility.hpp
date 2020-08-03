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
    // static double topology_prior_ratio, time_prior_ratio, rate_prior_ratio;
    // static double topology_proposal_ratio, time_proposal_ratio, rate_proposal_ratio;
    static double const STD_DEV = 0.005;
    static double const STD_DEV_M10 = 0.02;
    static double const BD_ALPHA = 1000 * 0.01;
    static double const BD_BETA = 1000 * (1 - 0.01);
    static std::discrete_distribution<int> discrete_distribution{1, 1, 2, 2, 2, 1};
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

    static void outputToCSVFile(const std::string output_filename, std::vector<double> &alpha_vector, std::vector<double> &beta_vector, std::vector<double> &m10_vector, 
        std::vector<double> &rootage_vector, std::vector<double> &loglikelihood_vector, std::vector<TreeManip> &tree_manip_vector)
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
            myfile << loglikelihood_vector[i] << ", ";
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

} // namespace strom