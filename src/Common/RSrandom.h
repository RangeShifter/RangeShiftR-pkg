/*------------------------------------------------------------------------------

RangeShifter v2.0 RSrandom

Implements the RSrandom class

Author: Steve Palmer, University of Aberdeen

Last updated: 12 June 2018 by Steve Palmer

------------------------------------------------------------------------------*/

#ifndef RSrandomH
#define RSrandomH

#include <stdlib.h>
#include <fstream>
//#include <iostream>

#include "Version.h"
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
using namespace std;

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

#if RSDEBUG
extern ofstream DEBUGLOG;
#endif

//---------------------------------------------------------------------------
#endif


