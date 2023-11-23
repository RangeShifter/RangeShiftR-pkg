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

#include <string>
#include <iostream>
#include <cassert>

using namespace std;

#if LINUX_CLUSTER || R_CMD
#include <unistd.h>
#else
#include <tchar.h>
#endif

void assert_error(const string& exptd_err_msg, void (*x)(void)) {
	string err_msg{ "No error.\n" };
	try { x(); }
	catch (exception& e) { err_msg = e.what(); }
	assert(err_msg == exptd_err_msg);
}

void run_unit_tests() {
	cout << "******* Unit test output *******" << endl;
	testRSrandom();
	testIndividual();
	cout << endl << "************************" << endl;
}

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
	run_unit_tests();
	cout << "All tests have completed." << endl;
	return 0; // if tests succeed, we are happy
# endif // NDEBUG
}