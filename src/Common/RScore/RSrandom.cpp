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

//--------------- 2.) New version of RSrandom.cpp


#if RSDEBUG
#include "Parameters.h"
extern paramSim* paramsSim;
// ofstream RSRANDOMLOG;
#endif

int RS_random_seed = 0;

// C'tor
RSrandom::RSrandom()
{
#if RSDEBUG
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
#endif // RSDEBUG

#if BATCH && RSDEBUG
    DEBUGLOG << "RSrandom::RSrandom(): RS_random_seed=" << RS_random_seed << endl;
#endif // RSDEBUG

    // set up Mersenne Twister RNG
    gen = new mt19937(RS_random_seed);

    // Set up standard uniform distribution
    pRandom01 = new uniform_real_distribution<double>(0.0, 1.0);
    // Set up standard normal distribution
    pNormal = new normal_distribution<double>(0.0, 1.0);
}

RSrandom::~RSrandom(void)
{
    delete gen;
    if(pRandom01 != 0)
	delete pRandom01;
    if(pNormal != 0)
	delete pNormal;
}

mt19937 RSrandom::getRNG(void)
{
    return *gen;
}

double RSrandom::Random(void)
{
    // return random number between 0 and 1
    return pRandom01->operator()(*gen);
}

int RSrandom::IRandom(int min, int max)
{
    // return random integer in the interval min <= x <= max
    uniform_int_distribution<int> unif(min, max);
    return unif(*gen);
}

int RSrandom::Bernoulli(double p)
{
    return Random() < p;
}

double RSrandom::Normal(double mean, double sd)
{
    return mean + sd * pNormal->operator()(*gen);
}

int RSrandom::Poisson(double mean)
{
    poisson_distribution<int> poiss(mean);
    return poiss(*gen);
}


//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------
