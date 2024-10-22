/*----------------------------------------------------------------------------
 *
 *	Copyright (C) 2020 Anne-Kathleen Malchow, Greta Bocedi, Stephen C.F. Palmer, Justin M.J. Travis, Damaris Zurell
 *
 *	This file is part of RangeShiftR.
 *
 *	RangeShiftR is free software: you can redistribute it and/or modify
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
 *	along with RangeShiftR. If not, see <https://www.gnu.org/licenses/>.
 *
 --------------------------------------------------------------------------*/


/*------------------------------------------------------------------------------

RangeShifter v2.0 Rinterface

Implements the interface to the R-package RangeshiftR.

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

#include "RScore/Parameters.h"
#include "RScore/Landscape.h"
#include "RScore/Species.h"
#include "RScore/SubCommunity.h"
#include "RScore/RSrandom.h"
#include "RScore/Model.h"
#include "RScore/Management.h"
#include "RScore/Cell.h"



//---------------------------------------------------------------------------


//Rcpp::List run_from_R(Rcpp::S4, Rcpp::String); // entry functions from R
//Rcpp::List BatchMainFile(string, Rcpp::S4);	 // missing dependend functions for file reading/parsing, from BatchMode.h/.cpp
Rcpp::List BatchMainR(std::string, Rcpp::S4);

bool ReadLandParamsR(Landscape*, Rcpp::S4);
int ReadDynLandR(Landscape*, Rcpp::S4);
int ReadParametersR(Landscape*, Rcpp::S4);
int ReadStageStructureR(Rcpp::S4);
int ReadEmigrationR(Rcpp::S4);
int ReadTransferR(Landscape*, Rcpp::S4);
int ReadSettlementR(Rcpp::S4);
int ReadInitialisationR(Landscape*, Rcpp::S4);
int ReadGeneticsR(Rcpp::S4);
int ReadTraitsR(Rcpp::S4);
const sex_t intToSex(const int& i);
TraitType stringToTraitType(const std::string& str);
TraitType addSexDepToTrait(const TraitType& t, const sex_t& sex);
GenParamType strToGenParamType(const string& str);
ExpressionType stringToExpressionType(const std::string& str);
map<GenParamType, float> NumericToParameterMap(string parameterString, Rcpp::NumericVector parameter);
set<int> selectRandomLociPositions(int noLoci, const int& genomeSize);
DistributionType stringToDistributionType(const std::string& str);
void setUpSpeciesTrait(string TraitTypeR, set<int> positions, string ExpressionTypeR,
                       string initDistR, Rcpp::NumericVector initParamsR,
                       string DominanceDistR, Rcpp::NumericVector DominanceParamsR,
                       bool isInherited, string MutationDistR,
                       Rcpp::NumericVector MutationParamsR, float MutationRateR,
                       int sexdep, bool isOutputR);
int ReadTranslocationR(Landscape*,Rcpp::S4);

int ParseInitIndsFileR(wifstream&);
int ReadInitIndsFileR(int,Landscape*);


Rcpp::List RunBatchR(int, int, Rcpp::S4);
void setglobalvarsR(Rcpp::S4);

struct DispersalTraitInputOptions {
    bool isEmigIndVar = false;
    bool isEmigDensDep = false;
    bool isEmigSexDep = false;

    bool isSettIndVar = false;
    bool isSettSexDep = false;

    bool isKernTransfIndVar = false;
    bool isKernTransfSexDep = false;
    bool usesTwoKernels = false;

    bool isSMSTransfIndVar = false;
    bool usesSMSGoalBias = false;

    bool isCRWTransfIndVar = false;
};



//---------------------------------------------------------------------------

#if !RSWIN64
string check_bom(string); // check BOM of a text file for UTF-16 encoding
#endif
rasterdata ParseRasterHead(string); //parse, read and return ASCII raster head data

void BatchErrorR(string,int,int,string);
void BatchErrorR(string,int,int,string,string);

//void CtrlFormatError(void);
//void ArchFormatError(void);
void FormatErrorR(string,int);
void OpenErrorR(string,string);
void EOFerrorR(string);
void StreamErrorR(string);
void ArchFormatErrorR(void);
//void FileOK(string,int,int);
//void FileHeadersOK(string);
//void SimulnCountError(string);


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
extern std::uint32_t RS_random_seed;			// see RSrandom.cpp

//---------------------------------------------------------------------------
#endif
