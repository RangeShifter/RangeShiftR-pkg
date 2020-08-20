/*------------------------------------------------------------------------------

RangeShifter v2.0 RSrandom

Implements the RSrandom class

Author: Steve Palmer, University of Aberdeen

Last updated: 17 January 2020 by Steve Palmer

Modified: 19 February 2020 by Anne-Kathleen Malchow, Humboldt University Berlin

------------------------------------------------------------------------------*/

#ifndef RSrandomH
#define RSrandomH

#include <stdlib.h>
#include <fstream>
//#include <iostream>

#include "Version.h"

using namespace std;

#if RSDEBUG
extern ofstream DEBUGLOG;
#endif

#if !RS_RCPP

	#if CLUSTER
	//#include <random>
	//#include <tr1/random>
	#include "maths.h"
	#else
	#if RSWIN64
	#include <dinkumware64/random>
	#else
	#include <dinkumware/random>
	#endif
	#endif

	class RSrandom {

	public:
	#if RS_ABC
		RSrandom(int);
	#else
		RSrandom(void);
	#endif
		~RSrandom(void);
		double Random(void);
		int IRandom(int,int);
		int Bernoulli(double);
		double Normal(double,double);
		int Poisson(double);
	#if RS_ABC
		double Beta(double,double);
		double Gamma(double,double);
	#endif

	private:
		double normal_x2; int normal_x2_valid; // variables used by Normal distribution
	#if !CLUSTER
		tr1::uniform_real<> *pRandom01;
		tr1::normal_distribution<> *pNormal;
	#if RS_ABC
		tr1::mt19937 *gen;
	#endif
	#endif
	};

#else // if RS_RCPP

	#include <cmath>
	#include <random>
	#if RSWIN64
	#include <ctime>
	#endif

	class RSrandom {

	public:
		RSrandom(std::int64_t);              // if int is negative, a random seed will be generated, else it is used as seed
		~RSrandom(void);
		double Random(void);
		int IRandom(int,int);
		int Bernoulli(double);
		double Normal(double,double);
		int Poisson(double);
	/* ADDITIONAL DISTRIBUTIONS
		double Beta(double,double);
		double Gamma(double,double); // !! make sure coorect definition is used: using shape and scale (as defined here) OR using shape/alpha and rate/beta (=1/scale)
		double Cauchy(double,double);
	*/

	private:
		mt19937 *gen;
		std::uniform_real_distribution<> *pRandom01;
		std::normal_distribution<> *pNormal;
	};

#endif // RS_RCPP

//---------------------------------------------------------------------------
#endif // RSrandomH
