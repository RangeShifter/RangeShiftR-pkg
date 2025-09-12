/*----------------------------------------------------------------------------
 *
 *	Copyright (C) 2020 Greta Bocedi, Stephen C.F. Palmer, Justin M.J. Travis, Anne-Kathleen Malchow, Damaris Zurell
 *
 *	This file is part of RangeShifter.
 *
 *	RangeShifter is free software: you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *	the Free Software Foundation, either version 3 of the License, or
 *	(at your option) any later version.
 *
 *	RangeShifter is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *	GNU General Public License for more details.
 *
 *	You should have received a copy of the GNU General Public License
 *	along with RangeShifter. If not, see <https://www.gnu.org/licenses/>.
 *
 --------------------------------------------------------------------------*/

#include "RSrandom.h"
#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP
#ifndef NDEBUG
#include "Parameters.h"
extern paramSim* paramsSim;
#endif

#if RS_RCPP
std::uint32_t RS_random_seed = 0;
#else
int RS_random_seed = 0;
#endif

// C'tor
#if RS_RCPP
// if parameter seed is negative, a random seed will be generated, else it is used as seed
RSrandom::RSrandom(std::int64_t seed)
{
	// get seed
	std::vector<std::uint32_t> random_seed(3);
	random_seed[0] = 1967593562;
	random_seed[1] = 3271254416;
	if (seed < 0) {
		// random seed
#if RSWIN64
		random_seed[2] = std::time(NULL) + (seed * (-17));
#else
		std::random_device device;
		random_seed[2] = device();
#endif
	}
	else {
		// fixed seed
		random_seed[2] = seed;
	}

	RS_random_seed = random_seed[2];

	// set up Mersenne Twister random number generator with seed sequence
	std::seed_seq seq(random_seed.begin(), random_seed.end());

#ifdef _OPENMP
	int nb_generators = omp_get_max_threads();
	gens.reserve(nb_generators);
	for (int i = 0; i < nb_generators; i++)
		gens.emplace_back(seq);
#else
	gens.reserve(1);
	gens.emplace_back(seq);
#endif // _OPENMP


	// Set up standard uniform distribution
	pRandom01 = new uniform_real_distribution<double>(0.0, 1.0);
	// Set up standard normal distribution
	pNormal = new normal_distribution<double>(0.0, 1.0);
}
#else
RSrandom::RSrandom() {

#ifndef NDEBUG
	// fixed seed
	RS_random_seed = 666;
#else
	// random seed
#if LINUX_CLUSTER
	std::random_device device;
	RS_random_seed = device(); // old versions of g++ on Windows return a constant value within a given Windows
	// session; in this case better use time stamp
#else
	RS_random_seed = std::time(NULL);
#endif
#endif // NDEBUG

	// set up Mersenne Twister RNG
#ifdef _OPENMP
	int nb_generators = omp_get_max_threads();
	gens.reserve(nb_generators);
	for (int i = 0; i < nb_generators; i++)
		gens.emplace_back(RS_random_seed + i);
#else
	gens.reserve(1);
	gens.emplace_back(RS_random_seed);
#endif // _OPENMP	

	// Set up standard uniform distribution
	pRandom01 = new uniform_real_distribution<double>(0.0, 1.0);
	// Set up standard normal distribution
	pNormal = new normal_distribution<double>(0.0, 1.0);
}
#endif // RS_RCPP

RSrandom::~RSrandom(void) {
    gens.clear();
    if (pRandom01 != 0)
        delete pRandom01;
    if (pNormal != 0)
        delete pNormal;
}

mt19937 RSrandom::getRNG() {
#ifdef _OPENMP
	return gens[omp_get_thread_num() % gens.size()];
#else
	return gens[0];
#endif // _OPENMP
}

double RSrandom::Random() {
    // return random number between 0 and 1
#ifdef _OPENMP
	return pRandom01->operator()(gens[omp_get_thread_num() % gens.size()]);
#else
	return pRandom01->operator()(gens[0]);
#endif // _OPENMP
}

int RSrandom::IRandom(int min, int max) {
    // return random integer in the interval min <= x <= max
    if (min == max)
        return min;

    uniform_int_distribution<int> unif(min, max);
#ifdef _OPENMP
	return unif(gens[omp_get_thread_num() % gens.size()]);
#else
	return unif(gens[0]);
#endif // _OPENMP
}

float RSrandom::FRandom(float min, float max) {
    if (min == max) return min;
    // return random double in the interval min <= x <= max
    uniform_real_distribution<float> unif(min, max);

#ifdef _OPENMP
	return unif(gens[omp_get_thread_num() % gens.size()]);
#else
	return unif(gens[0]);
#endif
}

int RSrandom::Bernoulli(double p) {
	if (p < 0) {
		throw runtime_error("Bernoulli's p cannot be negative.\n");
	}
	if (p > 1) {
		throw runtime_error("Bernoulli's p cannot be above 1.\n");
	}
    return Random() < p;
}

int RSrandom::Binomial(const int& n, const double& p) {
	binomial_distribution<int> binom(n, p);
#ifdef _OPENMP
	return binom(gens[omp_get_thread_num() % gens.size()]);
#else
	return binom(gens[0]);
#endif
}

double RSrandom::Normal(double mean, double sd) {
#ifdef _OPENMP
	return mean + sd * pNormal->operator()(gens[omp_get_thread_num() % gens.size()]);
#else
	return mean + sd * pNormal->operator()(gens[0]);
#endif // _OPENMP
}

int RSrandom::Poisson(double mean) {
    poisson_distribution<int> poiss(mean);
#ifdef _OPENMP
	return poiss(gens[omp_get_thread_num() % gens.size()]);
#else
	return poiss(gens[0]);
#endif // _OPENMP
}

double RSrandom::Gamma(double shape, double scale) { //scale  = mean/shape, shape must be positive and scale can be positive or negative

    gamma_distribution<> gamma(shape, abs(scale));

#ifdef _OPENMP
	double x = poiss(gens[omp_get_thread_num() % gens.size()]);
#else
	double x = gamma(gens[0]);
#endif // _OPENMP
    if (scale < 0) x = -x;
    return x;
}

double RSrandom::NegExp(double mean) {
    double r1 = 0.0000001 + this->Random() * (1.0 - 0.0000001);
    double x = (-1.0 * mean) * log(r1);
    return x;
}

void RSrandom::fixNewSeed(int seed) {
#ifdef _OPENMP
	for (int i = 0; i < omp_get_max_threads(); i++)
		gens[i].seed(seed + i);
#else
	gens[0].seed(seed);
#endif // _OPENMP
}

//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------

#ifndef NDEBUG
#if !RS_RCPP
	void testRSrandom() {

		{
			// Bernoulli distribution
			// Abuse cases
			assert_error("Bernoulli's p cannot be negative.\n", []{
				RSrandom rsr;
				rsr.Bernoulli(-0.3);
				});
			assert_error("Bernoulli's p cannot be above 1.\n", [] {
				RSrandom rsr;
				rsr.Bernoulli(1.1);
				});
			// Use cases
			RSrandom rsr;
			assert(rsr.Bernoulli(0) == 0);
			assert(rsr.Bernoulli(1) == 1);
			int bern_trial = rsr.Bernoulli(0.5);
			assert(bern_trial == 0 || bern_trial == 1);
		}
	}
#endif
#endif // NDEBUG
//---------------------------------------------------------------------------
