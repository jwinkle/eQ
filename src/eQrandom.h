#ifndef EQRANDOM_H
#define EQRANDOM_H

#include <random>
#include <chrono>
#include <iostream>
#include <memory>
namespace eQ {

class uniformRandomNumber
{
public:
	uniformRandomNumber()
	{
//      use a random device generator to generate seed:
//            std::random_device rdev{};
//            generator.seed(rdev());
//            distribution =  std::make_shared<std::uniform_real_distribution<double>>(0.0,1.0);
		// construct a trivial random generator engine from a time-based seed:
		auto seed = std::chrono::system_clock::now().time_since_epoch().count();
		std::cout<<"Called eQ::uniformRandomNumber() default constructor, which generated seed="<<seed<<std::endl;
		init(size_t(seed));
	}
	uniformRandomNumber(size_t seed)
	{
		init(seed);
	}
	void init(size_t seed)
	{
		generator.seed(unsigned(seed));
		distribution =  std::make_shared<std::uniform_real_distribution<double>>(0.0,1.0);
		std::cout<<"\t Initializing uniformRandomNumber(size_t seed) with seed="<<seed<<std::endl;
	}
	double randomNumber(double scale)
	{
		return scale * (*distribution)(generator);
	}
	double randomNumber()
	{
		return (*distribution)(generator);
	}
private:
	std::default_random_engine generator;
	std::shared_ptr<std::uniform_real_distribution<double>> distribution;
};
//=======================================================================================
class logNormalRandomNumber
{
public:
	logNormalRandomNumber(std::vector<std::pair<double, double>> &meansStdevs)
	{
		//      use a random device generator to generate seed:
		//        std::random_device rdev{};
		//        generator.seed(rdev());
				// construct a trivial random generator engine from a time-based seed:
				auto seed = std::chrono::system_clock::now().time_since_epoch().count();
				generator.seed(static_cast<unsigned long>(seed));
				for(auto &init : meansStdevs)
				{
					//convert to (lognormal) mean and variance of the log of target distribution
					double lnmu = log(init.first/sqrt(1.0 + init.second/(init.first*init.first)));
					double lns2 = log(1.0 + init.second/(init.first*init.first));

					lndistributions.push_back(std::make_shared<std::lognormal_distribution<double>>(lnmu, lns2));
				}
	}
	double getLogNormalRandomNumber(size_t which)
	{
		return (lndistributions.size() > which) ? (*(lndistributions[which]))(generator) : 0.0;
	}
private:
	std::default_random_engine      generator;
	std::vector<std::shared_ptr<std::lognormal_distribution<double>>>
									lndistributions;
};

//    static double getLogNormalMean(size_t which)
//    {
//        return (lndistributions.size() > which) ? lndistributions[which]->m(): 0.0;
//    }
//    static double getLogNormalRandomNumber(size_t which, double mean, double var)
//    {
//        if(lndistributions.size() > which)
//        {
//            //update the mean, var
//            lndistributions[which]->param(std::lognormal_distribution<double>::param_type(mean, var));
//            return (*(lndistributions[which]))(generator);
//        }
//        return mean;
//    }



}

#endif // EQRANDOM_H
