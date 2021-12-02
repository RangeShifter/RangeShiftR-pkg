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

RangeShifter v2.0 FractalGenerator

Implements the midpoint displacement algorithm for generating a fractal Landscape,
following:

Saupe, D. (1988). Algorithms for random fractals. In: The Science of Fractal Images
(eds. Pietgen, H.O. & Saupe, D.). Springer, New York, pp. 71–113.

and applying a diamond-square algorithm partially derived from a C++ implementation
published on-line by Nick O'Brien 10/8/2018 at:

https://medium.com/@nickobrien/
diamond-square-algorithm-explanation-and-c-implementation-5efa891e486f
(last accessed 17/10/2021)

For full details of RangeShifter, please see:
Bocedi G., Palmer S.C.F., Pe’er G., Heikkinen R.K., Matsinos Y.G., Watts K.
and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
eco-evolutionary dynamics and species’ responses to environmental changes.
Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

Authors: Greta Bocedi & Steve Palmer, University of Aberdeen

Last updated: 2 December 2021 by Steve Palmer

------------------------------------------------------------------------------*/

#ifndef FractalGeneratorH
#define FractalGeneratorH

#include <algorithm>
#include <vector>
//using namespace std;

#include "Parameters.h"

class land
{
 public:
	land();
	int x_coord;
	int y_coord;
	float value;
//	int avail; // if 0 the patch is not available as habitat, if 1 it is
 private:
};

// IMPORTANT NOTE: X AND Y ARE TRANSPOSED, i.e. X IS THE VERTICAL CO-ORDINATE
// ==========================================================================

vector<land>& fractal_landscape(
	int,		// X dimension (Y of LandScape)
	int,		// Y dimension (X of LandScape)
	double,	// Hurst exponent
	double,	// proportion of NON-suitable habitat
	double,	// maximum quality value
	double	// minimum quality value
);
bool compare(const land&, const land&);

double squareStep(double**,int,int,int,int,int,double);
double diamondStep(double**,int,int,int,int,int,double);
void diamondSquare(double**,int,int,int,double,double);

extern RSrandom *pRandom;
#if RSDEBUG
extern void DebugGUI(string);
extern ofstream MUTNLOG;
#if BATCH
extern ofstream DEBUGLOG;
#endif
#endif

//---------------------------------------------------------------------------
#endif
