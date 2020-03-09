/*------------------------------------------------------------------------------

RangeShifter v2.0 Rinterface

Implements the interface to the R-package RangeshiftR.

Includes the following functions:

For full details of RangeShifter, please see:
Bocedi G., Palmer S.C.F., Pe'er G., Heikkinen R.K., Matsinos Y.G., Watts K.
and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
eco-evolutionary dynamics and species' responses to environmental changes.
Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

Author: Anne-Kathleen Malchow, Humboldt University Berlin
        large parts modified from 'Main.cpp' and 'BatchMode.cpp' created by
        Steve Palmer, University of Aberdeen

------------------------------------------------------------------------------*/

#ifndef RinterfaceH
#define RinterfaceH

//---------------------------------------------------------------------------

#include <string>
#include <stdio.h>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <locale>
#if !RSWIN64
#include <codecvt>
#endif

using namespace std;

#include <Rcpp.h>
#include "Version.h"
#include "Parameters.h"
#include "Landscape.h"
#include "Species.h"
#include "SubCommunity.h"
#include "RandomCheck.h"
#include "Model.h"


//---------------------------------------------------------------------------

//Rcpp::List run_from_R(Rcpp::S4, Rcpp::String); // entry functions from R
//Rcpp::List BatchMainFile(string, Rcpp::S4);	 // missing dependend functions for file reading/parsing, from BatchMode.h/.cpp
#if RSDEBUG
Rcpp::List BatchMainR(std::string, Rcpp::S4);
#endif // RSDEBUG
bool ReadLandParamsR(Landscape*, Rcpp::S4);
int ReadDynLandR(Landscape*, Rcpp::S4);
int ReadParametersR(Landscape*, Rcpp::S4);
int ReadStageStructureR(Rcpp::S4);
int ReadEmigrationR(Rcpp::S4);
int ReadTransferR(Landscape*, Rcpp::S4);
int ReadSettlementR(Rcpp::S4);
int ReadInitialisationR(Landscape*, Rcpp::S4);
int ReadGeneticsR(Rcpp::S4);

int ParseInitIndsFileR(wifstream&);
int ReadInitIndsFileR(int,Landscape*);
int ReadArchFileR(wifstream&);

Rcpp::List RunBatchR(int, int, Rcpp::S4);
void setglobalvarsR(Rcpp::S4);


//---------------------------------------------------------------------------

#if !RSWIN64
string check_bom(string); // check BOM of a text file for UTF-16 encoding
#endif
rasterdata ParseRasterHead(string); //parse, read and return ASCII raster head data

void BatchErrorR(string,int,int,string);
void BatchErrorR(string,int,int,string,string);

//void CtrlFormatError(void);
//void ArchFormatError(void);
#if VIRTUALECOLOGIST
//void SampleFormatError(void);
#endif // VIRTUALECOLOGIST
void FormatErrorR(string,int);
void OpenErrorR(string,string);
void EOFerrorR(string);
void StreamErrorR(string);
void ArchFormatErrorR(void);
//void FileOK(string,int,int);
//void FileHeadersOK(string);
//void SimulnCountError(string);



//---------------------------------------------------------------------------

// Dummy functions corresponding to those used in GUI version

/* Batch mode of v2.0 currently has no facility to save maps (unless initiated from GUI).
*/

const string Int2Str(const int);
const string Float2Str(const float);
const string Double2Str(const double);
void MemoLine(string);
#if RSDEBUG
void DebugGUI(string);
#endif

traitCanvas SetupTraitCanvas(void);


//---------------------------------------------------------------------------

// external pointers

#if RSDEBUG
extern ofstream DEBUGLOG;
#endif

extern paramGrad *paramsGrad;
extern paramStoch *paramsStoch;
extern paramInit *paramsInit;
extern paramSim *paramsSim;

extern Species *pSpecies;
extern string costmapname;	// see FormMove.cpp (VCL) OR Main.cpp (batch)
extern string genfilename;	// see FormGenetics.cpp (VCL) OR Main.cpp (batch)
#if VIRTUALECOLOGIST
extern string locfilename;		// see FormVirtEcol.cpp (VCL) OR Main.cpp (batch)
extern string patchfilename;	// see [NOT YET CODED FOR GUI] (VCL) OR Main.cpp (batch)
#endif // VIRTUALECOLOGIST 
#if TEMPMORT
extern string mortfilename;	// see [NOT YET CODED FOR GUI] (VCL) OR Main.cpp (batch)
#endif // TEMPMORT
#if !CLUSTER || RS_RCPP
extern int RS_random_seed;			// see RSrandom.cpp
#endif // !CLUSTER || RS_RCPP

//---------------------------------------------------------------------------
#endif
