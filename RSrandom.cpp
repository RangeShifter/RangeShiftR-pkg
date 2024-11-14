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
#ifndef NDEBUG
#include "Parameters.h"
extern paramSim* paramsSim;
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
	gen = new mt19937(seq);

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
	gen = new mt19937(RS_random_seed);
	// Set up standard uniform distribution
	pRandom01 = new uniform_real_distribution<double>(0.0, 1.0);
	// Set up standard normal distribution
	pNormal = new normal_distribution<double>(0.0, 1.0);
}
#endif // RS_RCPP

RSrandom::~RSrandom(void) {
    delete gen;
    if (pRandom01 != 0)
        delete pRandom01;
    if (pNormal != 0)
        delete pNormal;
}

mt19937 RSrandom::getRNG(void) {
    return *gen;
}

double RSrandom::Random(void) {
    // return random number between 0 and 1
    return pRandom01->operator()(*gen);
}

int RSrandom::IRandom(int min, int max) {
    // return random integer in the interval min <= x <= max
    if (min == max)
        return min;

    uniform_int_distribution<int> unif(min, max);
    return unif(*gen);
}

float RSrandom::FRandom(float min, float max) {
    if (min == max)
        return min;
    // return random double in the interval min <= x <= max
    uniform_real_distribution<float> unif(min, max);
    return unif(*gen);
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
	return binom(*gen);
}

double RSrandom::Normal(double mean, double sd) {
    return mean + sd * pNormal->operator()(*gen);
}

int RSrandom::Poisson(double mean) {
    poisson_distribution<int> poiss(mean);
    return poiss(*gen);
}

double RSrandom::Gamma(double shape, double scale) { //scale  = mean/shape, shape must be positive and scale can be positive or negative

    gamma_distribution<> gamma(shape, abs(scale));

    double x = gamma(*gen);
    if (scale < 0) x = -x;
    return x;
}

double RSrandom::NegExp(double mean) {
    double r1 = 0.0000001 + this->Random() * (1.0 - 0.0000001);
    double x = (-1.0 * mean) * log(r1);
    return x;
}

void RSrandom::fixNewSeed(int seed) {
    gen->seed(seed);
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
