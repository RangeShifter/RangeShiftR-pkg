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
    non_dispersed = false; // do not consider non-dispersed individuals
    std::vector<int> translocation_years; // Number of years of translocation
    std::map< int, std::vector <locn> > source; // Source patch or cell: should be a vector of arrays
    std::map< int, std::vector <locn> > target; // Target patch or cell
    std::map< int, std::vector <int> > nb; // number of ttanslocated individuals
    std::map< int, std::vector <int> > min_age; // Minimum age of translocated individuals
    std::map< int, std::vector <int> > max_age; // Maximum age of translocated individuals
    std::map< int, std::vector <int> > stage; // Stage of translocated individuals
    std::map< int, std::vector <int> > sex; // Sex of translocated individuals

}

Management::~Management(void) {
    translocation_years.clear();
    source.clear();
    target.clear();
    nb.clear();
    min_age.clear();
    max_age.clear();
    stage.clear();
}

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
    t.catching_rate = catching_rate;
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
    catching_rate = t.catching_rate;
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
                                 , Species* pSpecies
                                 ){
    Rcpp::Rcout << "Start translocation events in year " << yr << endl;
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
#if RS_RCPP
        Rcpp::Rcout << "Translocation event " << e << " in year " << yr << endl;
#endif
        // find the source patch
        Patch* s_patch;
        Population* s_pPop;
        if(ppLand.patchModel){
            if(pLandscape->existsPatch(source_it->second[e].x)){
#if RS_RCPP
                Rcpp::Rcout << "Source patch exist." << endl;
#endif
                s_patch = pLandscape->findPatch(source_it->second[e].x);
                if (s_patch != 0) {
                    // test if population in patch is not zero
                    s_pPop = (Population*)s_patch->getPopn((intptr)pSpecies); // returns the population of the species in that cell
                    if (s_pPop != 0 && s_pPop->getNInds() > 0){
                    } else {
#if RS_RCPP
                        Rcpp::Rcout << "Population does not exist in source patch or is 0! skipping translocation event." << endl;
#endif
                        return;
                    }
                } else {
#if RS_RCPP
                    Rcpp::Rcout << "Source patch was found but NULL! skipping translocation event." << endl; // not sure if this ever happens
#endif
                    return;
                }
    //
            } else{
#if RS_RCPP
                Rcpp::Rcout << "Source patch was not found in landscape! skipping translocation event." << endl;
#endif
                return;
            }
       } else{
           Cell* pCell;
           pCell = pLandscape->findCell(source_it->second[e].x, source_it->second[e].y);
           if (pCell != 0) {
#if RS_RCPP
               Rcpp::Rcout << "Source cell was found" << endl;
#endif
               intptr s_ppatch = pCell->getPatch();
               if (s_ppatch != 0) {
                   s_patch = (Patch*)s_ppatch;
                   // test if population in patch is not zero
                   s_pPop = (Population*)s_patch->getPopn((intptr)pSpecies); // returns the population of the species in that cell
                   if (s_pPop != 0 && s_pPop->getNInds() > 0){
                   } else {
#if RS_RCPP
                       Rcpp::Rcout << "Population does not exist in source cell or is 0! skipping translocation event." << endl;
#endif
                       return;
                   }
               } else {
#if RS_RCPP
                   Rcpp::Rcout << "Source cell does not exist! skipping translocation event." << endl;
#endif
                   return;
               }
           } else {
#if RS_RCPP
               Rcpp::Rcout << "Cell does not belong to landscape! skipping translocation event." << endl;
#endif
               return;
           }
        }
       // find the target patch and check for existence
       Patch* t_patch;
       Population* t_pPop;
       if(ppLand.patchModel){
           if(pLandscape->existsPatch(target_it->second[e].x)){
#if RS_RCPP
                Rcpp::Rcout << "Target patch exist." << endl;
#endif
                t_patch = pLandscape->findPatch(target_it->second[e].x);
           } else{
#if RS_RCPP
               Rcpp::Rcout << "Target patch was not found in landscape! skipping translocation event." << endl;
#endif
               return;
           }
       } else{
           Cell* pCell;
           pCell = pLandscape->findCell(target_it->second[e].x, target_it->second[e].y);
           if (pCell != 0) {
#if RS_RCPP
               Rcpp::Rcout << "Target cell was found" << endl;
#endif
               intptr t_ppatch = pCell->getPatch();
               if (t_ppatch != 0) {
                   t_patch = (Patch*)t_ppatch;
               } else {
#if RS_RCPP
                   Rcpp::Rcout << "Target cell does not exist! skipping translocation event." << endl;
#endif
                   return;
               }
           } else {
#if RS_RCPP
               Rcpp::Rcout << "Target cell does not belong to landscape! skipping translocation event." << endl;
#endif
               return;
           }
       }

       // only if source and target cell/patch exist, we can translocate individuals:
        // get individuals with the given characteristics in that population
        int min_age = min_age_it->second[e];
        int max_age = max_age_it->second[e];
        int stage = stage_it->second[e];
        int sex = sex_it->second[e];
        int nb = nb_it->second[e];
        int nbSampledInds = 0;
        // We made already sure by now that in s_pPop at least some individuals exist
        nbSampledInds = s_pPop->sampleIndividuals(nb, min_age, max_age, stage, sex); // checking values was done when reading in the parameters
        popStats s_stats = s_pPop->getStats();
        Individual* catched_individual;
        int translocated = 0;
        // loop over all indsividuals, extract sampled individuals, try to catch individual + translocate them to new patch
        for (int j = 0; j < s_stats.nInds; j++) {
            // if there are individuals to catch
            if(s_pPop->getSizeSampledInds()){
                // if this individual is matching one of the sampled individuals
                catched_individual = s_pPop->catchIndividual(catching_rate, j); // catch individual in the source patch
                if (catched_individual !=NULL) { // translocated individual - has already been removed from natal population
                    // Check if a population of this species already exists in target patch t_patch
                    t_pPop = (Population*)t_patch->getPopn((intptr)pSpecies);
                    if (t_pPop == 0) { // translocated individual is the first in a previously uninhabited patch
#if RS_RCPP
                        Rcpp::Rcout << "Population does not exist in target patch. Creating new population." << endl;
#endif
                        // create a new population in the corresponding sub-community
                        SubCommunity* pSubComm = (SubCommunity*)t_patch->getSubComm();
                        t_pPop = pSubComm->newPopn(pLandscape, pSpecies, t_patch, 0);
                    }
                    catched_individual->setStatus(10); // make sure individual is not dispersing after the translocation
                    t_pPop->recruit(catched_individual); // recruit individual to target population TODO:  maybe use a specified function which also updates pCurrCell + pPrevCell to a random cell in target patch?
                    translocated ++;
                    // NOTE:
                    // the variables pCurrCell and pPrevCell are not updated! These are important for the dispersal process!
                    // currently, translocated individuals are not considered as potential emigrants, thus there is no problem in changing that
                    // however, if we want to consider dispersal events after translocation, we need to adapt that; but that would also mean, that we might loose the information
                    // about the natal patch of an individual?
                    simParams sim = paramsSim->getSim();
                    if (sim.outConnect) { // increment connectivity totals
                        int newpatch = t_patch->getSeqNum();
                        int prevpatch = s_patch->getSeqNum();
                        pLandscape->incrConnectMatrix(prevpatch, newpatch);
                    }

                }
            }
        }
#if RS_RCPP
        Rcpp::Rcout << "Successfully translocated " << translocated << " out of " << nb_it->second[e] << " individuals in translocation event " << e <<"." << endl;
#endif
        // remove pointers to sampled individuals
        s_pPop->clean();
    }
};