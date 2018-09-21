//---------------------------------------------------------------------------

#pragma hdrstop

#include "RSrandom.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)

#if RSDEBUG
#include "Parameters.h"
extern paramSim *paramsSim;
#include <fstream>
//ofstream RSRANDOMLOG;
#endif

#if !RS_ABC

#if !CLUSTER
// set up Mersenne Twister random generator
#if GROUPDISP
#if RSDEBUG
int RS_random_seed = 123;
tr1::mt19937 gen(RS_random_seed);
#else
int RS_random_seed = time(0);
tr1::mt19937 gen(RS_random_seed);
#endif // RSDEBUG 
#else
#if RSDEBUG
tr1::mt19937 gen(666);
#else
tr1::mt19937 gen(time(0));
#endif // RSDEBUG 
#endif // GROUPDISP 
#endif // !CLUSTER

#endif // !RS_ABC

#if RS_ABC
int RS_random_seed;
RSrandom::RSrandom(int seed)
#else
RSrandom::RSrandom(void)
#endif
{

#if RS_ABC

#if !CLUSTER
// set up Mersenne Twister random generator
#if RSDEBUG
RS_random_seed = 666;
gen = new tr1::mt19937 (RS_random_seed);
// RS random initialisation log
#else
RS_random_seed = time(0)+seed;
gen = new tr1::mt19937 (RS_random_seed);
#endif // RSDEBUG
#endif // !CLUSTER

#else

#if GROUPDISP || RS_ABC
#if RSDEBUG
DEBUGLOG << "RSrandom::RSrandom(): RS_random_seed=" << RS_random_seed
	<< endl;
#endif // RSDEBUG 
#endif // GROUPDISP || RS_ABC

#endif // RS_ABC 

#if RSDEBUG
// RS random initialisation log added by SCFP 25/8/16
//string name = paramsSim->getDir(2) + "RandomLog.txt";
//RSRANDOMLOG.open(name.c_str());
//RSRANDOMLOG << "RSrandom::RSrandom(): creating new random stream" << endl;
#endif
/*
#  if defined(BOOST_HAS_CPP_0X)
cout << "BOOST_HAS_CPP_0X" << endl;
#  else
cout << "***NOT*** BOOST_HAS_CPP_0X" << endl;
#  endif
#  if defined(BOOST_HAS_TR1_RANDOM)
cout << "BOOST_HAS_TR1_RANDOM" << endl;
#     ifdef BOOST_HAS_INCLUDE_NEXT
cout << "BOOST_HAS_INCLUDE_NEXT" << endl;
//#        include_next <random>
#     else
cout << "***NOT*** BOOST_HAS_INCLUDE_NEXT" << endl;
cout << "BOOST_TR1_STD_HEADER" << endl;
//#        include BOOST_TR1_STD_HEADER(random)
#     endif
#  else
cout << "***NOT*** BOOST_HAS_TR1_RANDOM" << endl;
#  endif
*/

normal_x2_valid = 0;
#if !CLUSTER
// Set up standard uniform distribution
pRandom01 = new	tr1::uniform_real<> (0.0,1.0);
// Set up standard normal distribution
pNormal = new	tr1::normal_distribution<> (0.0,1.0);
#endif

//cout << endl;
//for (int i = 0; i < 5; i++) {
//	cout << gen() << " ";
//}
//cout << endl << endl;

#if RSDEBUG
//RSRANDOMLOG.close(); RSRANDOMLOG.clear();
#endif
}

RSrandom::~RSrandom(void) {
#if !CLUSTER
if (pRandom01 != 0) delete pRandom01;
if (pNormal != 0) delete pNormal;
#endif
}

double RSrandom::Random(void) {
#if CLUSTER
return unif_rand();
#else
#if RS_ABC
return pRandom01->operator()(*gen);
#else
return pRandom01->operator()(gen);
#endif
//double result_Random;
//result_Random = pRandom01->operator()(*gen);
//#if RSDEBUG
//DEBUGLOG << "RSrandom::Random(): result_Random=" << result_Random
//	<< endl;
//#endif
//return result_Random;
#endif
}

int RSrandom::IRandom(int min,int max) {
//#if CLUSTER
//return irand(min,max);
//#else
//tr1::uniform_int<> dist(min,max);
//return dist(gen);
//#endif
// output random integer in the interval min <= x <= max (copied from mersenne.cpp)
int r;
r = int((max - min + 1) * Random()) + min; // multiply interval with random and truncate
if (r > max) r = max;
if (max < min) return 0x80000000;
return r;
}

int RSrandom::Bernoulli(double p) {
#if RSDEBUG
//DEBUGLOG << "RSrandom::Bernoulli(): p=" << p
//	<< endl;
#endif
#if CLUSTER
return unif_rand() < p;
#else
return Random() < p;
#endif
}

double RSrandom::Normal(double mean,double sd) {
#if RSDEBUG
//DEBUGLOG << "RSrandom::Normal(): mean=" << mean << " sd=" << sd
//	<< endl;
#endif
#if CLUSTER
// normal distribution derived from Agner Fog's method
double normal_x1;          // first random coordinate (normal_x2 is member of class)
double w;                  // radius
#if RSDEBUG
//DEBUGLOG << "RSrandom::Normal(): normal_x2_valid=" << normal_x2_valid
//	<< endl;
#endif
if (normal_x2_valid) {     // we have a valid result from last call
	normal_x2_valid = 0;
	return normal_x2 * sd + mean;
}
// make two normally distributed variates by Box-Muller transformation
do {
	normal_x1 = 2. * Random() - 1.;
	normal_x2 = 2. * Random() - 1.;
	w = normal_x1*normal_x1 + normal_x2*normal_x2;
} while (w >= 1. || w < 1E-30);
w = sqrt(log(w)*(-2./w));
normal_x1 *= w;  normal_x2 *= w;     // normal_x1 and normal_x2 are independent normally distributed variates
normal_x2_valid = 1;                 // save normal_x2 for next call
#if RSDEBUG
//DEBUGLOG << "RSrandom::Normal(): normal_x1=" << normal_x1
//	<< endl;
#endif
return normal_x1 * sd + mean;
#else
//double norm = pNormal->operator()(gen);
//#if RSDEBUG
//DEBUGLOG << "RSrandom::Normal(): norm=" << norm
//	<< endl;
//#endif
//return mean + sd * norm;
#if RS_ABC
return mean + sd * pNormal->operator()(*gen);
#else
return mean + sd * pNormal->operator()(gen);
#endif
#endif
}

int RSrandom::Poisson(double mean) {
#if CLUSTER
return rpois(mean);
#else
#if RS_ABC
if (mean > 50.0) {
	return this->Normal(mean,sqrt(mean));
}
else {
	tr1::poisson_distribution<> poiss(mean);
	return poiss(*gen);
}
#else
tr1::poisson_distribution<> poiss(mean);
return poiss(gen);
#endif
#endif
}

#if RS_ABC

// Beta distribution - sample from two gamma distributions
double RSrandom::Beta(double p0,double p1) {
double g0,g1,beta;
#if RSDEBUG
//DEBUGLOG << "RSrandom::Beta(): p0=" << p0 << " p1=" << p1
//	<< endl;
#endif
if (p0 > 0.0 && p1 > 0.0) { // valid beta parameters
#if CLUSTER
	g0 = rgamma(p0,1.0);
	g1 = rgamma(p1,1.0);
#else
	tr1::gamma_distribution<> gamma0(p0,1.0);
	tr1::gamma_distribution<> gamma1(p1,1.0);
	g0 = gamma0(*gen);
	g1 = gamma1(*gen);
#endif // CLUSTER
	beta = g0 / (g0 + g1);
#if RSDEBUG
//DEBUGLOG << "RSrandom::Beta(): g0=" << g0 << " g1=" << g1 << " beta=" << beta
//	<< endl;
#endif
}
else { // return invalid value
	beta = -666.0;
}
return  beta;
}

// Gamma distribution
double RSrandom::Gamma(double p0,double p1) {
double p2,gamma;
if (p0 > 0.0 && p1 > 0.0) { // valid gamma parameters
	p2 = 1.0 / p1;
#if CLUSTER
	gamma = rgamma(p0,p2);
#else
	tr1::gamma_distribution<> gamma0(p0,p2);
	gamma = gamma0(*gen);
#endif // CLUSTER
}
else { // return invalid value
	gamma = -666.0;
}
return  gamma;
}

#endif // RS_ABC

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

