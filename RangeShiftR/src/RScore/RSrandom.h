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
//#include <iostream>

//#if RS_RCPP && !R_CMD
#include "../Version.h"
//#endif

//#if !RS_RCPP && R_CMD
//#include "../../Batch/Version.h"
//#endif

using namespace std;

#if RSDEBUG
extern ofstream DEBUGLOG;
#endif



	#include <cmath>
	#include <random>
	#if RSWIN64
	#include <ctime>
	#endif

	class RSrandom {

	public:
		RSrandom(std::int64_t);       // if int is negative, a random seed will be generated, else it is used as seed
		~RSrandom(void);
		mt19937 getRNG(void);
		double Random(void);
		int IRandom(int,int);
		int Bernoulli(double);
		double Normal(double,double);
		int Poisson(double);
	/* ADDITIONAL DISTRIBUTIONS
		double Beta(double,double);
		double Gamma(double,double); // !! make sure correct definition is used: using shape and scale (as defined here) OR using shape/alpha and rate/beta (=1/scale)
		double Cauchy(double,double);
	*/

	private:
		mt19937 *gen;
		std::uniform_real_distribution<> *pRandom01;
		std::normal_distribution<> *pNormal;
	};



//---------------------------------------------------------------------------

#endif // RSrandomH



