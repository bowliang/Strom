#include <string>
#include <vector>
#include <iostream>
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
    static std::default_random_engine generator;

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

} // namespace strom