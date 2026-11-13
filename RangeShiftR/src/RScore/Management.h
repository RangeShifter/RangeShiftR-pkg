/*----------------------------------------------------------------------------
 *
 *	Copyright (C) 2026 Greta Bocedi, Stephen C.F. Palmer, Justin M.J. Travis, Anne-Kathleen Malchow, Roslyn Henry, Théo Pannetier, Jette Wolff, Damaris Zurell
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
 * File Created by Jette Wolff
 --------------------------------------------------------------------------*/


/*------------------------------------------------------------------------------

 RangeShifter v2.0 Parameters

 Implements the following classes:

 paramManagement  - Management parameters
 paramTranslocation  - Translocation parameters
 paramHarvesting  - Harvesting parameters


 Last updated: 12 March 2024 by Jette Reeg

 ------------------------------------------------------------------------------*/


#ifndef ManagementH
#define ManagementH

#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <vector>
#include <map>
using namespace std;
#if RS_RCPP
#include <RcppArmadillo.h>
#endif
#include "Parameters.h"
#include "Species.h"
#include "Cell.h"
#include "Landscape.h"

#include "SubCommunity.h"
#include "Population.h"


#if RS_RCPP
typedef intptr_t intptr;
#else
typedef unsigned long long intptr;
#endif // RS_RCPP



//---------------------------------------------------------------------------

/*
 * Management settings
 */

// Structure for management parameters
struct managementParams {
    bool translocation; // Translocation
    bool harvesting; // Harvesting
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

// Structure for harvesting parameters
struct harvestingParams {
    double harvesting_success; // Harvesting success rate
    std::vector<int> harvesting_years; // Number of years of harvesting events -> will be increased at the beginning of a simulation
    std::map< int, std::vector <locn> > harvestLoc; // Patch or cell: should be a vector of arrays
    std::map< int, std::vector <int> > harvestThres; // Threshold for harvesting
    std::map< int, std::vector <int> > harvestNb; // number of harvested individuals
    std::map< int, std::vector <int> > harvestMin_age; // Minimum age of harvested individuals
    std::map< int, std::vector <int> > harvestMax_age; // Maximum age of harvested individuals
    std::map< int, std::vector <int> > harvestStage; // Stage of harvested individuals
    std::map< int, std::vector <int> > harvestSex; // Sex of harvested individuals
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
    void setHarvestingParams( // function to set harvesting parameters
            const harvestingParams	// structure holding harvesting parameters
    );
    translocationParams getTranslocationParams(void);
    harvestingParams getHarvestingParams(void);
    void translocate(   // Translocation
            int  ,       // year of translocation
            Landscape* , // pointer to the landscape
            // Community*, // pointer to the community
            Species*   // pointer to the species
            );

    void harvest(   // Harvesting
            int  ,       // year of harvesting
            Landscape* , // pointer to the landscape
            // Community*, // pointer to the community
            Species*   // pointer to the species
            );

    //
    bool translocation; // Translocation
    double catching_rate; // Catching rate
    bool non_dispersed; // whether non-dispersed individuals should be translocated
    std::vector<int> translocation_years; // Number of years of translocation -> should be a dynamic vector
    std::map< int, std::vector <locn> > source; // Source patch or cell: should be a vector of arrays
    std::map< int, std::vector <locn> > target; // Target patch or cell
    std::map< int, std::vector <int> > nb; // number of ttanslocated individuals
    std::map< int, std::vector <int> > min_age; // Minimum age of translocated individuals
    std::map< int, std::vector <int> > max_age; // Maximum age of translocated individuals
    std::map< int, std::vector <int> > stage; // Stage of translocated individuals
    std::map< int, std::vector <int> > sex; // Sex of translocated individuals

    bool harvesting; // Harvesting
    double harvesting_success; // Harvesting success rate
    std::vector<int> harvesting_years; // Number of years of harvesting events -> should be a dynamic vector
    std::map< int, std::vector <locn> > harvestLoc; // Source patch or cell: should be a vector of arrays
    std::map< int, std::vector <int> > harvestThres; // Threshold for harvesting
    std::map< int, std::vector <int> > harvestNb; // number of harvested individuals
    std::map< int, std::vector <int> > harvestMin_age; // Minimum age of harvested individuals
    std::map< int, std::vector <int> > harvestMax_age; // Maximum age of harvested individuals
    std::map< int, std::vector <int> > harvestStage; // Stage of harvested individuals
    std::map< int, std::vector <int> > harvestSex; // Sex of harvested individuals

};

//---------------------------------------------------------------------------

extern paramSim *paramsSim;

//---------------------------------------------------------------------------
#endif
