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

RangeShifter v2.0 BatchMode

Functions for running in BATCH MODE

For full details of RangeShifter, please see:
Bocedi G., Palmer S.C.F., Pe’er G., Heikkinen R.K., Matsinos Y.G., Watts K.
and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
eco-evolutionary dynamics and species’ responses to environmental changes.
Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

Authors: Greta Bocedi & Steve Palmer, University of Aberdeen

Last updated: 26 October 2021 by Steve Palmer

------------------------------------------------------------------------------*/

#ifndef BatchModeH
#define BatchModeH

#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
using namespace std;

#include "../Common/RScore/Parameters.h"
#include "../Common/RScore/Landscape.h"
#include "../Common/RScore/Species.h"
#include "../Common/RScore/Model.h"

struct batchfiles {
	bool ok;
	int batchNum;
	int nSimuls;
	int nLandscapes;
	int patchmodel, resolution, landtype, maxNhab, speciesdist, distresolution;
	int reproductn;
	int repseasons;
	int stagestruct, stages, transfer;
	int sexesDem;		// no. of explicit sexes for demographic model
	int sexesDisp;	// no. of explicit sexes for dispersal model
	string parameterFile;
	string landFile;
	string stageStructFile;
	string emigrationFile;
	string transferFile;
	string settleFile;
	string geneticsFile;
	string initFile;
};

struct simCheck {
	bool newsimul;
	int simul,simlines,reqdsimlines,errors;
};

batchfiles ParseControlFile(string,string,string);
int ParseParameterFile(void);
int ParseLandFile(int,string);
int ParseDynamicFile(string,string);
int ParseStageFile(string);
int ParseTransitionFile(short,short);
int ParseWeightsFile(string);
int ParseEmigFile(void);
int ParseTransferFile(string);
int ParseSettleFile(void);
int ParseGeneticsFile(string);
int ParseArchFile(void);
int ParseInitFile(string);
int ParseInitIndsFile(void);
simCheck CheckStageSex(string,int,int,simCheck,int,int,int,int,int,bool,bool);

void BatchError(
	string,	// file name
	int,		// line number
	int,		// option
	string	// fieldname
);
/* Write error message to batch log file. Options are as follows:
0 - general message only, no reference to field name
1 - fieldname must be 0 or 1
2 - fieldname must be 0, 1 or 2
3 - fieldname must be 0, 1, 2 or 3
4 - fieldname must be from 0 to 4
5 - fieldname must be from 0 to 5
6 - fieldname must be from 0 to 6
7 - fieldname must be from 0 to 7
10 - fieldname must be >0
11 - fieldname must be >=1
12 - fieldname must be >=2
13 - fieldname must be >=3
18 - fieldname must be >1
19 - fieldname must be >=0
20 - fieldname must be between 0 and 1
21 - fieldname must be >1
33 - fieldname must be 1, 2 or 3
44 - fieldname must be from 1 to 4
55 - fieldname must be from 1 to 5
66 - fieldname must be from 1 to 6
100 - fieldname must be between 0 and 100
111 - simulation must match first simulation in ParameterFile
222 - simulations must be sequential integers
223 - seasons must be sequential integers
333 - columns must match no. of habitats
444 - columns must be one fewer than no. of stages
555 - columns must match no. of stages
666 - fieldname must be a unique positive integer
*/
void BatchError(
	string,	// file name
	int,		// line number
	int,		// option
	string,	// fieldname
	string	// fieldname2
);
/* Write error message to batch log file. Options are as follows:
1 - fieldname must be greater than fieldname2
2 - fieldname must be greater than or equal to fieldname2
3 - fieldname must be less than or equal to fieldname2
4 - fieldname must be less than fieldname2
*/

int power2check(int x);

void CtrlFormatError(void);
void ArchFormatError(void);
void FormatError(string,int);
void OpenError(string,string);
void EOFerror(string);
void FileOK(string,int,int);
void FileHeadersOK(string);
void SimulnCountError(string);

void RunBatch(int,int);
int ReadParameters(int,Landscape*);
int ReadLandFile(int);
int ReadLandFile(int,Landscape*);
int ReadDynLandFile(Landscape*);
int ReadStageStructure(int);
int ReadTransitionMatrix(
	short,	// no. of stages
	short,	// no. of sexes represented for demography 
	short,	// habitat index
	short		// season
);
int ReadStageWeights(int);
int ReadEmigration(int);
int ReadTransfer(int,Landscape*);
int ReadSettlement(int);
int ReadGenetics(int);
int ReadArchFile(string);
int ReadInitialisation(int,Landscape*);
int ReadInitIndsFile(int,Landscape*,string);

#if RSDEBUG
extern ofstream DEBUGLOG;
#endif

// external pointers to parameter sets
extern paramGrad *paramsGrad;
extern paramStoch *paramsStoch;
extern paramInit *paramsInit;
extern paramSim *paramsSim;

extern Species *pSpecies;
extern string costmapname;	// see FormMove.cpp (VCL) OR Main.cpp (batch)
extern string genfilename;	// see FormGenetics.cpp (VCL) OR Main.cpp (batch)
extern int RS_random_seed;

//---------------------------------------------------------------------------
#endif
