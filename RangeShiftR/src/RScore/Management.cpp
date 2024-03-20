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

#include "Management.h"

//---------------------------------------------------------------------------
/*
 * Initialize management class
 */

Management::Management(void) {
    translocation = false;
    catching_rate = 1.0; // Catching rate
    std::vector<int> translocation_years; // Number of years of translocation
    std::map< int, std::vector <locn> > source; // Source patch or cell: should be a vector of arrays
    std::map< int, std::vector <locn> > target; // Target patch or cell
    std::map< int, std::vector <int> > nb; // number of ttanslocated individuals
    std::map< int, std::vector <int> > min_age; // Minimum age of translocated individuals
    std::map< int, std::vector <int> > max_age; // Maximum age of translocated individuals
    std::map< int, std::vector <int> > stage; // Stage of translocated individuals
    std::map< int, std::vector <int> > sex; // Sex of translocated individuals

}

Management::~Management(void) {}

managementParams Management::getManagementParams(void) {
    managementParams m;
    m.translocation = translocation;
    return m;
}

void Management::setManagementParams(const managementParams m){
    translocation = m.translocation;
};

translocationParams Management::getTranslocationParams(void) {
    translocationParams t;
    t.translocation_years = translocation_years;
    t.source = source;
    t.target = target;
    t.nb = nb;
    t.min_age = min_age;
    t.max_age = max_age;
    t.stage = stage;
    t.sex = sex;
    return t;
}

// not sure if this is a good way, so won't use it for now
void Management::setTranslocationParams(const translocationParams t){
    translocation_years = t.translocation_years;
    source = t.source;
    target = t.target;
    nb = t.nb;
    min_age = t.min_age;
    max_age = t.max_age;
    stage = t.stage;
    sex = t.sex;

};

void Management::translocate(int yr
                                 , Landscape* pLandscape
                                 // , Community* pComm
                                 , Species* pSpecies
                                 ){
    landParams ppLand = pLandscape->getLandParams();
    auto it = nb.find(yr); // the number of translocation events is determined by the number of elements of the maps at year yr
    auto nb_it = nb.find(yr);
    auto source_it = source.find(yr);
    auto target_it = target.find(yr);
    auto min_age_it = min_age.find(yr);
    auto max_age_it = max_age.find(yr);
    auto stage_it = stage.find(yr);
    auto sex_it = sex.find(yr);

    // iterate over the number of events
    for (int e = 0; e < it->second.size(); e++) {
        // find the source patch
        Patch* patch;
        if(ppLand.patchModel){
            if(pLandscape->existsPatch(source_it->second[e].x)){
                patch = pLandscape->findPatch(source_it->second[e].x);
            }
       } else{
           Cell* pCell;
           pCell = pLandscape->findCell(source_it->second[e].x, source_it->second[e].y);
           if (pCell != 0) {
               intptr ppatch = pCell->getPatch();
               if (ppatch != 0) {
                   patch = (Patch*)ppatch;
               }
           }
        }
       const auto pPop = (Population*)patch->getPopn((intptr)pSpecies); // returns the population in that cell
        // get individuals with the given characteristics in that population
        int min_age = min_age_it->second[e];
        int max_age = max_age_it->second[e];
        int stage = stage_it->second[e];
        int sex = sex_it->second[e];
        std::set <Individual*> sampledInds;
        sampledInds = pPop->getIndsWithCharacteristics(min_age, max_age, stage, sex); // checking values was done when reading in the parameters
        // for the number of individuals to be translocated or the length of the vectorof individuals
        // sample randomly from sampledInds
        // try to catch the individual
        // if individual is caught: translocate to new location

    }

};