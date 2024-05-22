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

#if LINUX_CLUSTER || R_CMD
#include <unistd.h>
#else
#include <tchar.h>
#endif

#include <stdlib.h>
#include <string>
#include <iostream>
#include <cassert>
#include "Individual.h"
#include "Community.h"
#include "RSrandom.h"
#include "Utils.h"
#include "Parameters.h"
#include "Landscape.h"
#include "Species.h"
#include "SubCommunity.h"

using namespace std;

void run_unit_tests() {
	cout << "******* Unit test output *******" << endl;
	testRSrandom();
	testIndividual();
	cout << endl << "************************" << endl;
}

// Global vars
string habmapname, patchmapname, distnmapname;
string costmapname, genfilename;
string landFile;
paramGrad* paramsGrad;
paramStoch* paramsStoch;
paramInit* paramsInit;
paramSim* paramsSim;
RSrandom* pRandom;
ofstream DEBUGLOG;
ofstream MUTNLOG;
vector <string> hfnames;
Species* pSpecies;
Community* pComm;

#if LINUX_CLUSTER || RS_RCPP
int main(int argc, char* argv[])
#else
int _tmain(int argc, _TCHAR* argv[])
#endif
{
#ifdef NDEBUG
	cout << "This code is only for running tests and not meant to run in release." << endl;
	return 1;
# else
	assert(0.1 > 0.0); // assert does run correctly
	try
	{
		run_unit_tests();
	}
	catch (const std::exception& e)
	{
		cerr << endl << "Error: " << e.what() << endl;
	}
	cout << "All tests have completed." << endl;
	return 0; // if tests succeed, we are happy
# endif // NDEBUG
}
