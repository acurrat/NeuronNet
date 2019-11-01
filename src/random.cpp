#include "random.h"




RandomNumbers::RandomNumbers(unsigned long int s)
{
	if(s == 0) //on cherche à générer une valeur aléatoire
	{
		std::random_device rd;
		seed = rd(); //génère seed aléatoire
	} else {
		seed = s; //seed donnée par utilisateur non nulle
	}
	rng = std::mt19937(seed);
}


void RandomNumbers::uniform_double(std::vector<double>& vec, double lower, double upper)
{
	std::uniform_real_distribution<> distribution(lower, upper);
	///for(auto& val : vec) val = distribution(rng);
	for (auto i = vec.begin(); i != vec.end(); i++) *i = distribution(rng);
}


double RandomNumbers::uniform_double(double lower, double upper)
{
	std::uniform_real_distribution<> distribution(lower, upper);
	return distribution(rng);
}


void RandomNumbers::normal(std::vector<double>& vec, double mean, double sd)
{
	std::normal_distribution<> distribution(mean, sd);
	for(auto i = vec.begin(); i != vec.end(); i++) *i = distribution (rng);
}


double RandomNumbers::normal(double mean, double sd)
{
	std::normal_distribution<> distribution(mean, sd);
	return distribution(rng);
}


void RandomNumbers::poisson(std::vector<int>& vec, double mean)
{
	std::poisson_distribution<> distribution(mean);
	for(auto i = vec.begin(); i != vec.end(); i++) *i = distribution(rng);
}


int RandomNumbers::poisson(double mean)
{
	std::poisson_distribution<> distribution(mean);
	return distribution(rng);
}

