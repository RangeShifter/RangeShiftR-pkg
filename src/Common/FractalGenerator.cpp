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
 
 
//---------------------------------------------------------------------------
#if RS_EMBARCADERO
#pragma hdrstop
#endif

#include "FractalGenerator.h"
//---------------------------------------------------------------------------
#if RS_EMBARCADERO 
#pragma package(smart_init)
#endif

vector<land> patches;

//----- Landscape creation --------------------------------------------------

land::land(): x_coord(0), y_coord(0), value(0.0) {}

bool compare(const land& z, const land& zz) //compares only the values of the cells
{
return z.value < zz.value;
}

vector<land>& fractal_landscape(int X,int Y,double Hurst,double prop,
	double maxValue,double minValue)
{
#if RSDEBUG
DEBUGLOG << "fractal_landscape(): X=" << X << " Y=" << Y 
	<< " Hurst=" << Hurst << " prop=" << prop 
	<< " maxValue=" << maxValue << " minValue=" << minValue  
	<< endl;
#endif

double range;	//range to draw random numbers at each iteration
int Nno; 			// number of cells NON suitable as habitat

double **arena = new double *[X];
for(int ii = 0; ii < X; ii++) {
	arena[ii] = new double[Y];
}

patches.clear();
// initialise all the landscape with zeroes
for (int jj = 0; jj < X; jj++) {
	for (int ii = 0; ii < Y; ii++) {
		arena[jj][ii]=0;
	}
}

#if RSDEBUG
MUTNLOG << "X,Y,Value" << endl;
#endif

// initialise corners of largest square(s)
int size = min(X,Y)-1;
for (int x = 0; x < X; x+=size) {
	for (int y = 0; y < Y; y+=size) {
		arena[x][y] = 1.0 + pRandom->Random() * (maxValue-1.0);
#if RSDEBUG
MUTNLOG << x << "," << y << "," << arena[x][y] << endl;
#endif		
	}	
}

// run the recursive diamond-square algorithm to populate the landscape array
diamondSquare(arena,(min(X,Y)-1),X,Y,(maxValue-1.0),Hurst);

#if RSDEBUG
double meanq = 0.0;
double minq = 10000000.0;
double maxq = 0.0;
int count0 = 0;
for (int x = 0; x < X; x++) {
	for (int y = 0; y < Y; y++) {
		meanq += arena[x][y];
		if (arena[x][y] < minq) minq = arena[x][y];
		if (arena[x][y] > maxq) maxq = arena[x][y];
		if (arena[x][y] <= 0.0) count0++;
	}
}	
meanq /= (X*Y);
int above = 0;
int below = 0;
for (int x = 0; x < X; x++) {
	for (int y = 0; y < Y; y++) {
		if (arena[x][y] < meanq) below++; 
		if (arena[x][y] > meanq) above++; 
	}
}	
DEBUGLOG << "fractal_landscape(): 3333 - count0=" << count0 
	<< " meanq=" << meanq << " minq=" << minq << " maxq=" << maxq  
	<< " below=" << below << " above=" << above  
	<< endl;
#endif

// sort the cells and set Nno cells to be matrix, i.e. with K = 0

land *patch;
for (int x = 0; x < X; x++) // put all the cells with their values in a vector
{
	for (int y = 0; y < Y; y++)
	{
		patch = new land;
		patch->x_coord = x;
		patch->y_coord = y;
		patch->value = (float)arena[x][y];

		patches.push_back(*patch);

		delete patch;
	}
}
sort(patches.begin(),patches.end(),compare);  // sorts the vector

// change by SCFP 17/10/2021
// remove half of the unsuitable squares from the bottom end of the generated distribution
// and half from the top end to prevent retained non-zero distribution being skewed
//Nno = (int)(prop*(double)X*(double)Y);
Nno = (int)(prop*(double)X*(double)Y/2.0);

// variables for the rescaling
double min = (double)patches[Nno].value;        
//double max = (double)patches[X*Y-1].value;
double max = (double)patches[X*Y-Nno-1].value;

double diff = max - min;
double diffK = maxValue-minValue;
double new_value;

#if RSDEBUG
DEBUGLOG << "fractal_landscape(): 4444 - Nno=" << Nno 
	<< " min=" << min << " max=" << max << " diff=" << diff << " diffK=" << diffK   
	<< endl;
#endif

int npatches = patches.size();
for (int i = 0; i < Nno; i++) {
	new_value = 0.0;
#if RSDEBUG
//DEBUGLOG << "fractal_landscape(): 5555 - i=" << i 
//	<< " patches[i].value=" << patches[i].value 
//	<< " new_value=" << new_value 
//	<< endl;
#endif
	patches[i].value = new_value;
}
for (int i = Nno; i < (npatches-Nno); i++) {
	new_value = minValue + diffK * (patches[i].value - min) / diff;
#if RSDEBUG
//DEBUGLOG << "fractal_landscape(): 6666 - i=" << i 
//	<< " patches[i].value=" << patches[i].value 
//	<< " new_value=" << new_value 
//	<< endl;
#endif
	patches[i].value = new_value;	
}
for (int i = (npatches-Nno); i < npatches; i++) {
	new_value = 0.0;
#if RSDEBUG
//DEBUGLOG << "fractal_landscape(): 7777 - i=" << i 
//	<< " patches[i].value=" << patches[i].value 
//	<< " new_value=" << new_value 
//	<< endl;
#endif
	patches[i].value = new_value;
}

if (arena != NULL) {
#if RSDEBUG
//DebugGUI(("fractal_landscape(): arena=" + Int2Str((int)arena)
//	+ " X=" + Int2Str(X) + " Y=" + Int2Str(Y)
//	).c_str());
#endif
	for(int ii = 0; ii < X; ii++) {
#if RSDEBUG
//DebugGUI(("fractal_landscape(): ii=" + Int2Str(ii)
//	+ " arena[ii]=" + Int2Str((int)arena[ii])
//	).c_str());
#endif
		delete[] arena[ii];
	}
	delete[] arena;
}

return patches;

}

void diamondSquare(double **Array,int size,int dimX,int dimY,double range,double Hurst)
{
int delta = size / 2;
#if RSDEBUG
DEBUGLOG << "diamondSquare(): AAAA xx=" " size=" << size << " delta=" << delta 
	<< " dimX=" << dimX << " dimY=" << dimY << " range=" << range   
	<< endl;
#endif
if (delta < 1) return ;
//square steps
for (int y = delta; y < dimY; y+=size) {
	for (int x = delta; x < dimX; x+=size) {
#if RSDEBUG
//DEBUGLOG << "diamondSquare(): BBBB x=" << x << " y=" << y 
//	<< endl;
#endif
		Array[x][y] = squareStep(Array,x%dimX,y%dimY,delta,dimX,dimY,(pRandom->Random()*2.0*range - range));
#if RSDEBUG
MUTNLOG << x << "," << y << "," << Array[x][y] << endl;
#endif
	}
}
// diamond steps
int col = 0;
for (int x = 0; x < dimX; x += delta) {
	col++;
	if (col % 2 == 1) // odd column
		for (int y = delta; y < dimY; y += size) {
#if RSDEBUG
//DEBUGLOG << "diamondSquare(): CCCC x=" << x << " y=" << y << " (x % dimX)=" << x % dimX << " (y % dimY)=" << y % dimY 
//	<< endl;
#endif
			Array[x][y] = diamondStep(Array,x%dimX,y%dimY,delta,dimX,dimY,(pRandom->Random()*2.0*range - range));
#if RSDEBUG
MUTNLOG << x << "," << y << "," << Array[x][y] << endl;
#endif
		}
	else // even column
		for (int y = 0; y < dimY; y += size) {
#if RSDEBUG
//DEBUGLOG << "diamondSquare(): DDDD x=" << x << " y=" << y << " (x % dimX)=" << x % dimX << " (y % dimY)=" << y % dimY 
//	<< endl;
#endif
			Array[x][y] = diamondStep(Array,x%dimX,y%dimY,delta,dimX,dimY,(pRandom->Random()*2.0*range - range));
#if RSDEBUG
MUTNLOG << x << "," << y << "," << Array[x][y] << endl;
#endif
		}
}
range *= pow(2,-Hurst); // reduce the random number range
diamondSquare(Array,delta,dimX,dimY,range,Hurst);
}

double squareStep(double **Array,int x,int y,int delta,int dimX,int dimY,double r) {
int count = 0;
double result = 0.0;
#if RSDEBUG
//DEBUGLOG << "squareStep(): x=" << x << " y=" << y << " delta=" << delta
//	<< " dimX=" << dimX << " dimY=" << dimY;
#endif
if (x-delta >= 0 && y-delta >= 0) {
#if RSDEBUG
//DEBUGLOG	<< " Array[x-delta][y-delta]=" << Array[x-delta][y-delta]; 
#endif
	result += Array[x-delta][y-delta]; count++;
}
if (x-delta >= 0 && y+delta < dimY) {
#if RSDEBUG
//DEBUGLOG	<< " Array[x-delta][y+delta]=" << Array[x-delta][y+delta]; 
#endif
	result += Array[x-delta][y+delta]; count++; 
}
if (x+delta < dimX && y-delta >= 0) {
#if RSDEBUG
//DEBUGLOG	<< " Array[x+delta][y-delta]=" << Array[x+delta][y-delta]; 
#endif
	result += Array[x+delta][y-delta]; count++; 
}
if (x+delta < dimX && y+delta < dimY) {
#if RSDEBUG
//DEBUGLOG	<< " Array[x+delta][y+delta]=" << Array[x+delta][y+delta]; 
#endif
	result += Array[x+delta][y+delta]; count++; 
}
if (count > 0) result /= (double)count; result += r;
#if RSDEBUG
//DEBUGLOG	<< " COUNT=" << count	<< " RESULT=" << result; 
//DEBUGLOG	<< endl;
#endif
return result;
}

double diamondStep(double **Array,int x,int y,int delta,int dimX,int dimY,double r) {
int count = 0;
double result = 0.0;
#if RSDEBUG
//DEBUGLOG << "diamondStep(): x=" << x << " y=" << y << " delta=" << delta 
//	<< " dimX=" << dimX << " dimY=" << dimY;
#endif
if (y-delta >= 0) {
#if RSDEBUG
//DEBUGLOG	<< " Array[x][y-delta]=" << Array[x][y-delta];
#endif
	result += Array[x][y-delta]; count++; 	
}
if (y+delta < dimY) {
#if RSDEBUG
//DEBUGLOG	<< " Array[x][y+delta]=" << Array[x][y+delta];
#endif
	result += Array[x][y+delta]; count++; 
}
if (x-delta >= 0) {
#if RSDEBUG
//DEBUGLOG	<< " Array[x-delta][y]=" << Array[x-delta][y];
#endif
	result += Array[x-delta][y]; count++; 
}
if (x+delta < dimX) {
#if RSDEBUG
//DEBUGLOG	<< " Array[x+delta][y]=" << Array[x+delta][y];
#endif
	result += Array[x+delta][y]; count++; 
}
if (count > 0) result /= (double)count; result += r;
#if RSDEBUG
//DEBUGLOG	<< " COUNT=" << count	<< " RESULT=" << result; 
//DEBUGLOG	<< endl;
#endif
return result;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

