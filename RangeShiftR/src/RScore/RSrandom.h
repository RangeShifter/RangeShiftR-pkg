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


 /*------------------------------------------------------------------------------

 RangeShifter v2.0 RSrandom

 Implements the RSrandom class

 Authors: Steve Palmer, University of Aberdeen
				  Anne-Kathleen Malchow, Potsdam University

 Last updated: 12 January 2021 by Steve Palmer

 ------------------------------------------------------------------------------*/

#ifndef RSrandomH
#define RSrandomH

#include <stdlib.h>
#include <fstream>
#include <cassert>
#include <cmath>
#include <random>
#include <set>
#include "Utils.h"
#if !LINUX_CLUSTER
#include <ctime>
#endif

using namespace std;

#if RS_RCPP
typedef uint32_t seed_t;
#else
typedef int seed_t;
#endif

class RSrandom
{

public:
#if !RS_RCPP
	RSrandom(void);
#else
	RSrandom(std::int64_t);       // if int is negative, a random seed will be generated, else it is used as seed
#endif
	~RSrandom(void);
	double Random(void);
	int IRandom(int, int);
	float FRandom(float, float);
	int Bernoulli(double);
	int Binomial(const int& n, const double& p);
	double Normal(double, double);
	double Gamma(double, double);
	double NegExp(double);
	int Poisson(double);
	mt19937 getRNG(void);
	void fixNewSeed(int);
	seed_t getSeed() const { return RS_random_seed; };

private:
	seed_t RS_random_seed;
	mt19937* gen;
	std::uniform_real_distribution<>* pRandom01;
	std::normal_distribution<>* pNormal;
};

#ifndef NDEBUG
	void testRSrandom();
#endif // NDEBUG

//---------------------------------------------------------------------------

#endif // RSrandomH



