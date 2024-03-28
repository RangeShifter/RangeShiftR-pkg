/*----------------------------------------------------------------------------
 *
 *	Copyright (C) 2024 Greta Bocedi, Stephen C.F. Palmer, Justin M.J. Travis, Anne-Kathleen Malchow, Damaris Zurell, Jette Reeg
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

 RangeShifter v2.0 Parameters

 Implements the following classes:

 paramManagement  - Management parameters
 paramTranslocation  - Translocation parameters


 Last updated: 12 March 2024 by Jette Reeg

 ------------------------------------------------------------------------------*/


#ifndef ManagementH
#define ManagementH

//#if LINUX_CLUSTER
//#include <string.h>
//#else
#include <string>
//#endif
#include <fstream>
#include <iostream>
//#include <io.h>
#include <iomanip>
#include <stdlib.h>
#include <vector>
#include <map>
using namespace std;

#include <Rcpp.h> // for Rcpp::Rcout
#include "Parameters.h"
#include "Species.h"
#include "Cell.h"
#include "Landscape.h"

#include "SubCommunity.h"
#include "Population.h"


#if RS_RCPP
typedef intptr_t intptr;
#else
#if RSWIN64
typedef unsigned long long intptr;
#else
typedef unsigned int intptr;
#endif
#endif



//---------------------------------------------------------------------------

/*
 * Management settings
 */

// Structure for management parameters
struct managementParams {
    bool translocation; // Translocation
};

// Structure for translocation parameters
struct translocationParams {
    double catching_rate; // Catching rate
    std::vector<int> translocation_years; // Number of years of translocation -> will be increased at the beginning of a simulation
    std::map< int, std::vector <locn> > source; // Source patch or cell: should be a vector of arrays
    std::map< int, std::vector <locn> > target; // Target patch or cell
    std::map< int, std::vector <int> > nb; // number of ttanslocated individuals
    std::map< int, std::vector <int> > min_age; // Minimum age of translocated individuals
    std::map< int, std::vector <int> > max_age; // Maximum age of translocated individuals
    std::map< int, std::vector <int> > stage; // Stage of translocated individuals
    std::map< int, std::vector <int> > sex; // Sex of translocated individuals
};


//---------------------------------------------------------------------------

class Management{
public:
    Management(void);
    ~Management(void);
    void setManagementParams( // function to set management parameters
            const managementParams	// structure holding general management parameters
    );
    managementParams getManagementParams(void); // get management parameters
    void setTranslocationParams( // function to set translocation parameters
            const translocationParams	// structure holding translocation parameters
    );
    translocationParams getTranslocationParams(void);
    void translocate(   // Translocation
            int  ,       // year of translocation
            Landscape* , // pointer to the landscape
            // Community*, // pointer to the community
            Species*   // pointer to the species
            );

    //
    bool translocation; // Translocation
    double catching_rate; // Catching rate
    std::vector<int> translocation_years; // Number of years of translocation -> should be a dynamic vector
    std::map< int, std::vector <locn> > source; // Source patch or cell: should be a vector of arrays
    std::map< int, std::vector <locn> > target; // Target patch or cell
    std::map< int, std::vector <int> > nb; // number of ttanslocated individuals
    std::map< int, std::vector <int> > min_age; // Minimum age of translocated individuals
    std::map< int, std::vector <int> > max_age; // Maximum age of translocated individuals
    std::map< int, std::vector <int> > stage; // Stage of translocated individuals
    std::map< int, std::vector <int> > sex; // Sex of translocated individuals

};

//---------------------------------------------------------------------------

extern paramSim *paramsSim;

//---------------------------------------------------------------------------
#endif
