/*------------------------------------------------------------------------------

RangeShifter v2.0 RSrandom

Implements the RSrandom class

Author: Steve Palmer, University of Aberdeen

Last updated: 6 January 2020 by Steve Palmer

Modified: 14 January 2020 by Anne-Kathleen Malchow, Humboldt University Berlin


------------------------------------------------------------------------------*/

#ifndef RSrandomH
#define RSrandomH

#include <stdlib.h>
#include <fstream>
//#include <iostream>

#include "Version.h"

#include <cmath>

#include <random>

using namespace std;

class RSrandom {

public:
	RSrandom(int);              // if int is negative, a random seed will be generated, else it is used as seed
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

#if RSDEBUG
extern ofstream DEBUGLOG;
#endif

//---------------------------------------------------------------------------
#endif
