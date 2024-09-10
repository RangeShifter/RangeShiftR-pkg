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

#include "Model.h"

ofstream outPar;
using namespace std::chrono;
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
#if RS_RCPP && !R_CMD
Rcpp::List RunModel(Landscape* pLandscape, int seqsim)
#else
int RunModel(Landscape* pLandscape, int seqsim)
#endif
{
	int yr, totalInds;
	bool filesOK;

	landParams ppLand = pLandscape->getLandParams();
	envGradParams grad = paramsGrad->getGradient();
	envStochParams env = paramsStoch->getStoch();
	demogrParams dem = pSpecies->getDemogrParams();
	stageParams sstruct = pSpecies->getStageParams();
	//emigRules emig = pSpecies->getEmig();
	transferRules trfr = pSpecies->getTransferRules();
	initParams init = paramsInit->getInit();
	simParams sim = paramsSim->getSim();
	simView v = paramsSim->getViews();

	if (!ppLand.generated) {
		if (!ppLand.patchModel) { // cell-based landscape
			// create patches for suitable cells, adding unsuitable cells to the matrix
			// NB this is an overhead here, but is necessary in case the identity of
			// suitable habitats has been changed from one simulation to another (GUI or batch)
			// substantial time savings may result during simulation in certain landscapes
			// if using neutral markers, set up patches to sample from 
			pLandscape->allocatePatches(pSpecies);
		}
		pComm = new Community(pLandscape); // set up community
		// set up a sub-community associated with each patch (incl. the matrix)
		pLandscape->updateCarryingCapacity(pSpecies, 0, 0);
		patchData ppp;
		int npatches = pLandscape->patchCount();
		for (int i = 0; i < npatches; i++) {
			ppp = pLandscape->getPatchData(i);
			pComm->addSubComm(ppp.pPatch, ppp.patchNum); // SET UP ALL SUB-COMMUNITIES
		}
		if (init.seedType == 0 && init.freeType < 2 && init.initFrzYr > 0) {
			// restrict available landscape to initialised region
			pLandscape->setLandLimits(init.minSeedX, init.minSeedY,
				init.maxSeedX, init.maxSeedY);
		}
		else {
			pLandscape->resetLandLimits();
		}

		// Random patches are sampled once per landscape
		if (sim.patchSamplingOption == "random") {
			int nbToSample = pSpecies->getNbPatchesToSample();
			auto patchesToSample = pLandscape->samplePatches(sim.patchSamplingOption, nbToSample, pSpecies);
			pSpecies->setSamplePatchList(patchesToSample);
		}
	}

#if RS_RCPP && !R_CMD
	Rcpp::List list_outPop;
#endif

	// Loop through replicates
	for (int rep = 0; rep < sim.reps; rep++) {

		cout << "Running replicate " << rep + 1 << " / " << sim.reps << endl;

#if RS_RCPP && !R_CMD
		Rcpp::Rcout << endl << "starting replicate " << rep << endl;
#endif

		if (sim.saveVisits && !ppLand.generated) {
			pLandscape->resetVisits();
		}

		if (sim.fixReplicateSeed) {
			pRandom->fixNewSeed(rep);
		}
		patchChange patchchange;
		costChange costchange;
		int npatchchanges = pLandscape->numPatchChanges();
		int ncostchanges = pLandscape->numCostChanges();
		int ixpchchg = 0;
		int ixcostchg = 0;

		if (ppLand.generated) {
			// delete previous community (if any)
			// Note: this must be BEFORE the landscape is reset (as a sub-community accesses
			// its corresponding patch upon deletion)
			if (pComm != 0) delete pComm;
			// generate new cell-based landscape
			pLandscape->resetLand();
			pLandscape->generatePatches();
			pComm = new Community(pLandscape); // set up community
			// set up a sub-community associated with each patch (incl. the matrix)
			pLandscape->updateCarryingCapacity(pSpecies, 0, 0);
			patchData ppp;
			int npatches = pLandscape->patchCount();
			for (int i = 0; i < npatches; i++) {
				ppp = pLandscape->getPatchData(i);
#if RSWIN64
#if LINUX_CLUSTER
				pComm->addSubComm(ppp.pPatch, ppp.patchNum); // SET UP ALL SUB-COMMUNITIES
#else
				SubCommunity* pSubComm = pComm->addSubComm(ppp.pPatch, ppp.patchNum); // SET UP ALL SUB-COMMUNITIES
#endif
#else
				pComm->addSubComm(ppp.pPatch, ppp.patchNum); // SET UP ALL SUB-COMMUNITIES
#endif
			}
			if (sim.patchSamplingOption == "random") {
				// Then patches must be resampled for new landscape
				int nbToSample = pSpecies->getNbPatchesToSample();
				auto patchesToSample = pLandscape->samplePatches(sim.patchSamplingOption, nbToSample, pSpecies);
				pSpecies->setSamplePatchList(patchesToSample);
			}
		}
		if (init.seedType == 0 && init.freeType < 2 && init.initFrzYr > 0) {
			// restrict available landscape to initialised region
			pLandscape->setLandLimits(init.minSeedX, init.minSeedY,
				init.maxSeedX, init.maxSeedY);
		}
		else {
			pLandscape->resetLandLimits();
		}

		filesOK = true;
		if (rep == 0) {
			// open output files
			if (sim.outRange) { // open Range file
				if (!pComm->outRangeHeaders(pSpecies, ppLand.landNum)) {
					filesOK = false;
				}
			}
			if (sim.outOccup && sim.reps > 1)
				if (!pComm->outOccupancyHeaders(0)) {
					filesOK = false;
				}
			if (sim.outPop) {
				// open Population file
				if (!pComm->outPopHeaders(pSpecies, ppLand.landNum)) {
					filesOK = false;
				}
			}
			if (sim.outTraitsCells)
				if (!pComm->outTraitsHeaders(pSpecies, ppLand.landNum)) {
					filesOK = false;
				}
			if (sim.outTraitsRows)
				if (!pComm->outTraitsRowsHeaders(pSpecies, ppLand.landNum)) {
					filesOK = false;
				}
			if (sim.outConnect && ppLand.patchModel) // open Connectivity file
				if (!pLandscape->outConnectHeaders(0)) {
					filesOK = false;
				}
			if (sim.outputWeirCockerham || sim.outputWeirHill) { // open neutral genetics file
				if (!pComm->openNeutralOutputFile(pSpecies, ppLand.landNum)) {
					filesOK = false;
				}
			}
		}
		if (!filesOK) {
			// close any files which may be open
			if (sim.outRange) {
				pComm->outRangeHeaders(pSpecies, -999);
			}
			if (sim.outOccup && sim.reps > 1)
				pComm->outOccupancyHeaders(-999);
			if (sim.outPop) {
				pComm->outPopHeaders(pSpecies, -999);
			}
			if (sim.outTraitsCells)
				pComm->outTraitsHeaders(pSpecies, -999);
			if (sim.outTraitsRows)
				pComm->outTraitsRowsHeaders(pSpecies, -999);
			if (sim.outConnect && ppLand.patchModel)
				pLandscape->outConnectHeaders(-999);
			if (sim.outputWeirCockerham) {
				pComm->openNeutralOutputFile(pSpecies, -999);
			}
#if RS_RCPP && !R_CMD
			return Rcpp::List::create(Rcpp::Named("Errors") = 666);
#else
			return 666;
#endif
		}

		if (env.stoch && !env.local) {
			// create time series in case of global environmental stochasticity
			pLandscape->setGlobalStoch(sim.years + 1);
		}

		if (grad.gradient) { // set up environmental gradient
			pLandscape->setEnvGradient(pSpecies, true);
		}

		if (sim.outConnect && ppLand.patchModel)
			pLandscape->createConnectMatrix();

		// variables to control dynamic landscape
		landChange landChg; landChg.chgnum = 0; landChg.chgyear = 999999;
		if (!ppLand.generated && ppLand.dynamic) {
			landChg = pLandscape->getLandChange(0); // get first change year
		}

		// set up populations in the community
		pLandscape->updateCarryingCapacity(pSpecies, 0, 0);
		pComm->initialise(pSpecies, -1);
		bool updateland = false;
		int landIx = 0; // landscape change index

#if BATCH && RS_RCPP && !R_CMD
		Rcpp::Rcout << "RunModel(): completed initialisation " << endl;
#endif

		// open a new individuals file for each replicate
		if (sim.outInds)
			pComm->outInds(rep, 0, 0, ppLand.landNum);

		if (sim.outputGeneValues) {
			if (!pComm->openOutGenesFile(pSpecies->isDiploid(), ppLand.landNum, rep))
				throw logic_error("Output gene value file could not be initialised.");
		}

		// open a new genetics file for each replicate for per locus and pairwise stats
		if (sim.outputWeirCockerham) {
			pComm->openPerLocusFstFile(pSpecies, pLandscape, ppLand.landNum, rep);
		}
		if (sim.outputWeirHill) {
			pComm->openPairwiseFstFile(pSpecies, pLandscape, ppLand.landNum, rep);
		}
#if RS_RCPP
		// open a new movement paths file for each replicate
		if (sim.outPaths)
			pLandscape->outPathsHeaders(rep, 0);
#endif

		// years loop
		for (yr = 0; yr < sim.years; yr++) {
#if RS_RCPP && !R_CMD
			Rcpp::checkUserInterrupt();
#endif
			bool updateCC = false;
			if (yr < 4
				|| (yr < 31 && yr % 10 == 0)
				|| (yr < 301 && yr % 100 == 0)
				|| (yr < 3001 && yr % 1000 == 0)
				|| (yr < 30001 && yr % 10000 == 0)
				|| (yr < 300001 && yr % 100000 == 0)
				|| (yr < 3000001 && yr % 1000000 == 0)
				) {
#if RS_RCPP && !R_CMD
				Rcpp::Rcout << "Starting year " << yr << "..." << endl;
#else
				cout << "Starting year " << yr << endl;
#endif
			}
			if (init.seedType == 0 && init.freeType < 2) {
				// apply any range restrictions
				if (yr == init.initFrzYr) {
					// release initial frozen range - reset landscape to its full extent
					pLandscape->resetLandLimits();
					updateCC = true;
				}
				if (init.restrictRange) {
					if (yr > init.initFrzYr && yr < init.finalFrzYr) {
						if ((yr - init.initFrzYr) % init.restrictFreq == 0) {
							// apply dynamic range restriction
							commStats s = pComm->getStats();
							int minY = s.maxY - init.restrictRows;
							if (minY < 0) minY = 0;
							pLandscape->setLandLimits(ppLand.minX, minY, ppLand.maxX, ppLand.maxY);
							updateCC = true;
						}
					}
					if (yr == init.finalFrzYr) {
						// apply final range restriction
						commStats s = pComm->getStats();
						pLandscape->setLandLimits(ppLand.minX, s.minY, ppLand.maxX, s.maxY);
						updateCC = true;
					}
				}
			}
			// environmental gradient, stochasticity & local extinction
			// or dynamic landscape
			updateland = false;
			if (env.stoch || grad.gradient || ppLand.dynamic) {
				if (grad.shifting && yr > grad.shift_begin && yr < grad.shift_stop) {
					paramsGrad->incrOptY();
					pLandscape->setEnvGradient(pSpecies, false);
					updateCC = true;
				}
				if (env.stoch) {
					if (env.local) pLandscape->updateLocalStoch();
					updateCC = true;
				}
				if (ppLand.dynamic) {
					if (yr == landChg.chgyear) { // apply landscape change
						landIx = landChg.chgnum;
						updateland = updateCC = true;
						if (ppLand.patchModel) { // apply any patch changes
							Patch* pPatch;
							Cell* pCell;
							patchchange = pLandscape->getPatchChange(ixpchchg++);
							while (patchchange.chgnum <= landIx && ixpchchg <= npatchchanges) {
								// move cell from original patch to new patch
								pCell = pLandscape->findCell(patchchange.x, patchchange.y);
								if (patchchange.oldpatch != 0) { // not matrix
									pPatch = pLandscape->findPatch(patchchange.oldpatch);
									pPatch->removeCell(pCell);
								}
								if (patchchange.newpatch == 0) { // matrix
									pPatch = 0;
								}
								else {
									pPatch = pLandscape->findPatch(patchchange.newpatch);
									pPatch->addCell(pCell, patchchange.x, patchchange.y);
								}
								pCell->setPatch((intptr)pPatch);
								// get next patch change
								patchchange = pLandscape->getPatchChange(ixpchchg++);
							}
							ixpchchg--;
							pLandscape->resetPatches(); // reset patch limits
						}
						if (landChg.costfile != "NULL") { // apply any SMS cost changes
							Cell* pCell;
							costchange = pLandscape->getCostChange(ixcostchg++);
							while (costchange.chgnum <= landIx && ixcostchg <= ncostchanges) {
								pCell = pLandscape->findCell(costchange.x, costchange.y);
								if (pCell != 0) {
									pCell->setCost(costchange.newcost);
								}
								costchange = pLandscape->getCostChange(ixcostchg++);
							}
							ixcostchg--;
							pLandscape->resetEffCosts();
						}
						if (landIx < pLandscape->numLandChanges()) { // get next change
							landChg = pLandscape->getLandChange(landIx);
						}
						else {
							landChg.chgyear = 9999999;
						}
					}
				}
			} // end of environmental gradient, etc.

			if (updateCC) {
				pLandscape->updateCarryingCapacity(pSpecies, yr, landIx);
			}

			if (sim.outConnect && ppLand.patchModel)
				pLandscape->resetConnectMatrix();

			if (ppLand.dynamic && updateland) {
				if (trfr.usesMovtProc && trfr.moveType == 1) { // SMS
					if (!trfr.costMap) pLandscape->resetCosts(); // in case habitats have changed
				}
				// apply effects of landscape change to species present in changed patches
				pComm->patchChanges();
#if RS_RCPP
				pComm->dispersal(landIx, yr);
#else
				pComm->dispersal(landIx);
#endif // RS_RCPP
			}
			if (init.restrictRange) {
				// remove any population from region removed from restricted range
				if (yr > init.initFrzYr && yr < init.finalFrzYr) {
					if ((yr - init.initFrzYr) % init.restrictFreq == 0) {
						pComm->patchChanges();
					}
				}
			}

			if (init.seedType == 2) {
				// add any new initial individuals for the current year
				pComm->initialise(pSpecies, yr);
			}

			for (int gen = 0; gen < dem.repSeasons; gen++) // generation loop
			{
#ifndef NDEBUG
				// TEMPORARY RANDOM STREAM CHECK
				if (yr % 1 == 0)
				{
					DEBUGLOG << endl << "RunModel(): start of gen " << gen << " in year " << yr
						<< " for rep " << rep << " (";
					for (int i = 0; i < 5; i++) {
						int rrrr = pRandom->IRandom(1000, 2000);
						DEBUGLOG << " " << rrrr;
					}
					DEBUGLOG << " )" << endl;
				}
#endif

				// Output and pop. visualisation before reproduction
				if (v.viewPop || v.viewTraits || sim.outOccup
					|| sim.outTraitsCells || sim.outTraitsRows || sim.saveMaps)
					PreReproductionOutput(pLandscape, pComm, rep, yr, gen);
				// for non-structured population, also produce range and population output now
				if (!dem.stageStruct && (sim.outRange || sim.outPop))
					RangePopOutput(pComm, rep, yr, gen);
#if RS_RCPP && !R_CMD
				if (sim.ReturnPopRaster && sim.outPop && yr >= sim.outStartPop && yr % sim.outIntPop == 0) {
					list_outPop.push_back(pComm->addYearToPopList(rep, yr), "rep" + std::to_string(rep) + "_year" + std::to_string(yr));
				}
#endif
				// apply local extinction for generation 0 only
				// CHANGED TO *BEFORE* RANGE & POPN OUTPUT PRODUCTION IN v1.1,
				// SO THAT NOS. OF JUVENILES BORN CAN BE REPORTED
				if (!ppLand.patchModel && gen == 0) {
					if (env.localExt) pComm->localExtinction(0);
					if (grad.gradient && grad.gradType == 3) pComm->localExtinction(1);
				}

				// reproduction
				pComm->reproduction(yr);

				if (dem.stageStruct) {
					if (sstruct.survival == 0) { // at reproduction
						pComm->survival(0, 2, 1); // survival of all non-juvenile stages
					}
				}

				// Output and pop. visualisation AFTER reproduction
				if (dem.stageStruct && (sim.outRange || sim.outPop))
					RangePopOutput(pComm, rep, yr, gen);

				// Dispersal
				pComm->emigration();
#if RS_RCPP
				pComm->dispersal(landIx, yr);
#else
				pComm->dispersal(landIx);
#endif // RS_RCPP

				// survival part 0
				if (dem.stageStruct) {
					if (sstruct.survival == 0) { // at reproduction
						pComm->survival(0, 0, 1); // survival of juveniles only
					}
					if (sstruct.survival == 1) { // between reproduction events
						pComm->survival(0, 1, 1); // survival of all stages
					}
					if (sstruct.survival == 2) { // annually
						pComm->survival(0, 1, 0); // development only of all stages
					}
				}
				else { // non-structured population
					pComm->survival(0, 1, 1);
				}

				// output Individuals
				if (sim.outInds && yr >= sim.outStartInd && yr % sim.outIntInd == 0)
					pComm->outInds(rep, yr, gen, -1);

				if ((sim.outputGeneValues || sim.outputWeirCockerham || sim.outputWeirHill)
					&& yr >= sim.outStartGenetics
					&& yr % sim.outputGeneticInterval == 0) {

					simParams sim = paramsSim->getSim();
					if (sim.patchSamplingOption != "list" && sim.patchSamplingOption != "random") {
						// then patches must be re-sampled every gen
						int nbToSample = pSpecies->getNbPatchesToSample();
						auto patchesToSample = pLandscape->samplePatches(sim.patchSamplingOption, nbToSample, pSpecies);
						pSpecies->setSamplePatchList(patchesToSample);
					}
					// otherwise always use the user-specified list (even if patches are empty)
					pComm->sampleIndividuals(pSpecies);

					if (sim.outputGeneValues) {
						pComm->outputGeneValues(yr, gen, pSpecies);
					}
					if (sim.outputWeirCockerham || sim.outputWeirHill) {
						pComm->outNeutralGenetics(pSpecies, rep, yr, gen, sim.outputWeirCockerham, sim.outputWeirHill);
					}
				}
				if (dem.stageStruct) {
					pComm->survival(1, 0, 1);
				}
				else { // non-structured population
					pComm->survival(1, 0, 1);
				}

			} // end of the generation loop

			totalInds = pComm->totalInds();
			if (totalInds <= 0) { 
				cout << "All populations went extinct." << endl;
				yr++; 
				break; 
			}

			// Connectivity Matrix
			if (sim.outConnect && ppLand.patchModel
				&& yr >= sim.outStartConn && yr % sim.outIntConn == 0)
				pLandscape->outConnect(rep, yr);

			if (dem.stageStruct && sstruct.survival == 2) {  // annual survival - all stages
				pComm->survival(0, 1, 2);
				pComm->survival(1, 0, 1);
			}

			if (dem.stageStruct) {
				pComm->ageIncrement(); // increment age of all individuals
				if (sim.outInds && yr >= sim.outStartInd && yr % sim.outIntInd == 0)
					pComm->outInds(rep, yr, -1, -1); // list any individuals dying having reached maximum age
				pComm->survival(1, 0, 1);					// delete any such individuals
				totalInds = pComm->totalInds();
				if (totalInds <= 0) { 
					cout << "All populations went extinct." << endl;
					yr++; 
					break;

				}
			}

		} // end of the years loop

		// Final output
		// produce final summary output
		if (v.viewPop || v.viewTraits || sim.outOccup
			|| sim.outTraitsCells || sim.outTraitsRows || sim.saveMaps)
			PreReproductionOutput(pLandscape, pComm, rep, yr, 0);
		if (sim.outRange || sim.outPop)
			RangePopOutput(pComm, rep, yr, 0);

		pComm->resetPopns();

		//Reset the gradient optimum
		if (grad.gradient) paramsGrad->resetOptY();

		pLandscape->resetLandLimits();
		if (ppLand.patchModel && ppLand.dynamic && ixpchchg > 0) {
			// apply any patch changes to reset landscape to original configuration
			// (provided that at least one has already occurred)
			patchChange patchchange;
			Patch* pPatch;
			Cell* pCell;
			patchchange = pLandscape->getPatchChange(ixpchchg++);
			while (patchchange.chgnum <= 666666 && ixpchchg <= npatchchanges) {
				// move cell from original patch to new patch
				pCell = pLandscape->findCell(patchchange.x, patchchange.y);
				if (patchchange.oldpatch != 0) { // not matrix
					pPatch = pLandscape->findPatch(patchchange.oldpatch);
					pPatch->removeCell(pCell);
				}
				if (patchchange.newpatch == 0) { // matrix
					pPatch = 0;
				}
				else {
					pPatch = pLandscape->findPatch(patchchange.newpatch);
					pPatch->addCell(pCell, patchchange.x, patchchange.y);
				}
				pCell->setPatch((intptr)pPatch);
				// get next patch change
				patchchange = pLandscape->getPatchChange(ixpchchg++);
			}
			ixpchchg--;
			pLandscape->resetPatches();
		}
		if (ppLand.dynamic) {
			transferRules trfr = pSpecies->getTransferRules();
			if (trfr.usesMovtProc && trfr.moveType == 1) { // SMS
				if (ixcostchg > 0) {
					// apply any cost changes to reset landscape to original configuration
					// (provided that at least one has already occurred)
					Cell* pCell;
					costchange = pLandscape->getCostChange(ixcostchg++);
					while (costchange.chgnum <= 666666 && ixcostchg <= ncostchanges) {
						pCell = pLandscape->findCell(costchange.x, costchange.y);
						if (pCell != 0) {
							pCell->setCost(costchange.newcost);
						}
						costchange = pLandscape->getCostChange(ixcostchg++);
					}
					ixcostchg--;
					pLandscape->resetEffCosts();
				}
				if (!trfr.costMap) pLandscape->resetCosts(); // in case habitats have changed
			}
		}

		if (sim.outConnect && ppLand.patchModel)
			pLandscape->resetConnectMatrix(); // set connectivity matrix to zeroes

		if (sim.outInds) // close Individuals output file
			pComm->outInds(rep, 0, 0, -999);

		if (sim.outputGeneValues) { // close genetic values output file
			pComm->openOutGenesFile(false, -999, rep);
		}

		if (sim.outputWeirCockerham) //close per locus file 
			pComm->openPerLocusFstFile(pSpecies, pLandscape, -999, rep);
		if (sim.outputWeirHill) //close per locus file 
			pComm->openPairwiseFstFile(pSpecies, pLandscape, -999, rep);

		if (sim.saveVisits) {
			pLandscape->outVisits(rep, ppLand.landNum);
			pLandscape->resetVisits();
		}

#if RS_RCPP
		if (sim.outPaths)
			pLandscape->outPathsHeaders(rep, -999);
#endif

	} // end of the replicates loop

	if (sim.outConnect && ppLand.patchModel) {
		pLandscape->deleteConnectMatrix();
		pLandscape->outConnectHeaders(-999); // close Connectivity Matrix file
	}

	// Occupancy outputs
	if (sim.outOccup && sim.reps > 1) {
		pComm->outOccupancy();
		pComm->outOccSuit(v.viewGraph);
		pComm->deleteOccupancy((sim.years / sim.outIntOcc) + 1);
		pComm->outOccupancyHeaders(-999);
	}

	if (sim.outRange) {
		pComm->outRangeHeaders(pSpecies, -999); // close Range file
	}
	if (sim.outPop) {
		pComm->outPopHeaders(pSpecies, -999); // close Population file
	}
	if (sim.outTraitsCells)
		pComm->outTraitsHeaders(pSpecies, -999); // close Traits file
	if (sim.outTraitsRows)
		pComm->outTraitsRowsHeaders(pSpecies, -999); // close Traits rows file
	// close Individuals & Genetics output files if open
	// they can still be open if the simulation was stopped by the user
	if (sim.outInds) pComm->outInds(0, 0, 0, -999);
	if (sim.outputGeneValues) pComm->openOutGenesFile(0, -999, 0);
	if (sim.outputWeirCockerham || sim.outputWeirHill) {
		pComm->openNeutralOutputFile(pSpecies, -999);
	}
	if (sim.outputWeirCockerham) {
		pComm->openPerLocusFstFile(pSpecies, pLandscape, -999, 0);
	}
	if (sim.outputWeirHill) pComm->openPairwiseFstFile(pSpecies, pLandscape, -999, 0);

	delete pComm; 
	pComm = 0;

#if RS_RCPP && !R_CMD
	return list_outPop;
#else
	return 0;
#endif

}

#if LINUX_CLUSTER || RS_RCPP
// Check whether a specified directory path exists
bool is_directory(const char* pathname) {
	struct stat info;
	if (stat(pathname, &info) != 0) return false; // path does not exist
	if (S_ISDIR(info.st_mode)) return true;
	return false;
}
#endif

//---------------------------------------------------------------------------
bool CheckDirectory(void)
{
	bool errorfolder = false;

	string subfolder;

	subfolder = paramsSim->getDir(0) + "Inputs";
	const char* inputs = subfolder.c_str();
	if (!is_directory(inputs)) errorfolder = true;
	subfolder = paramsSim->getDir(0) + "Outputs";
	const char* outputs = subfolder.c_str();
	if (!is_directory(outputs)) errorfolder = true;
	subfolder = paramsSim->getDir(0) + "Output_Maps";
	const char* outputmaps = subfolder.c_str();
	if (!is_directory(outputmaps)) errorfolder = true;

	return errorfolder;
}

//---------------------------------------------------------------------------
//For outputs and population visualisations pre-reproduction
void PreReproductionOutput(Landscape* pLand, Community* pComm, int rep, int yr, int gen)
{
#ifndef NDEBUG
	landParams ppLand = pLand->getLandParams();
#endif
	simParams sim = paramsSim->getSim();
	simView v = paramsSim->getViews();

#ifndef NDEBUG
	DEBUGLOG << "PreReproductionOutput(): 11111 rep=" << rep << " yr=" << yr << " gen=" << gen
		<< " landNum=" << ppLand.landNum << " maxX=" << ppLand.maxX << " maxY=" << ppLand.maxY
		<< endl;
	DEBUGLOG << "PreReproductionOutput(): 11112 outRange=" << sim.outRange
		<< " outIntRange=" << sim.outIntRange
		<< " outPop=" << sim.outPop << " outIntPop=" << sim.outIntPop
		<< endl;
#endif

	// trait outputs and visualisation
	if (v.viewTraits
		|| ((sim.outTraitsCells && yr >= sim.outStartTraitCell && yr % sim.outIntTraitCell == 0) ||
			(sim.outTraitsRows && yr >= sim.outStartTraitRow && yr % sim.outIntTraitRow == 0)))
	{
		pComm->outTraits(pSpecies, rep, yr, gen);
	}
	if (sim.outOccup && yr % sim.outIntOcc == 0 && gen == 0)
		pComm->updateOccupancy(yr / sim.outIntOcc, rep);
}

//For outputs and population visualisations pre-reproduction
void RangePopOutput(Community* pComm, int rep, int yr, int gen)
{
	simParams sim = paramsSim->getSim();

	if (sim.outRange && (yr % sim.outIntRange == 0 || pComm->totalInds() <= 0))
		pComm->outRange(pSpecies, rep, yr, gen);

	if (sim.outPop && yr >= sim.outStartPop && yr % sim.outIntPop == 0)
		pComm->outPop(rep, yr, gen);

}

//---------------------------------------------------------------------------
void OutParameters(Landscape* pLandscape)
{
	double k;
	//int nrows,ncols,nsexes,nstages;
	int nsexes, nstages;

	landParams ppLand = pLandscape->getLandParams();
	genLandParams ppGenLand = pLandscape->getGenLandParams();
	envGradParams grad = paramsGrad->getGradient();
	envStochParams env = paramsStoch->getStoch();
	demogrParams dem = pSpecies->getDemogrParams();
	stageParams sstruct = pSpecies->getStageParams();
	emigRules emig = pSpecies->getEmigRules();
	transferRules trfr = pSpecies->getTransferRules();
	settleType sett = pSpecies->getSettle();
	settleRules srules;
	settleSteps ssteps;
	settleTraits settleDD;
	simParams sim = paramsSim->getSim();

	string name;
	if (sim.batchMode)
		name = paramsSim->getDir(2)
		+ "Batch" + to_string(sim.batchNum) + "_"
		+ "Sim" + to_string(sim.simulation)
		+ "_Land" + to_string(ppLand.landNum) + "_Parameters.txt";
	else
		name = paramsSim->getDir(2) + "Sim" + to_string(sim.simulation) + "_Parameters.txt";
	outPar.open(name.c_str());

	outPar << "RangeShifter 2.0 ";

#if !RS_RCPP
#if RSWIN64
	outPar << " - 64 bit implementation";
#else
	outPar << " - 32 bit implementation";
#endif
#endif
	outPar << endl;

	outPar << "================ ";

	outPar << "   =====================";
	outPar << endl << endl;

	outPar << "BATCH MODE \t";
	if (sim.batchMode) outPar << "yes" << endl; else outPar << "no" << endl;
#if RS_RCPP
	outPar << "SEED \t" << RS_random_seed << endl;
#endif
	outPar << "REPLICATES \t" << sim.reps << endl;
	outPar << "YEARS \t" << sim.years << endl;
	outPar << "REPRODUCTIVE SEASONS / YEAR\t" << dem.repSeasons << endl;
	if (ppLand.patchModel) {
		outPar << "PATCH-BASED MODEL" << endl;
		outPar << "No. PATCHES \t" << pLandscape->patchCount() - 1 << endl;
	}
	else
		outPar << "CELL-BASED MODEL" << endl;
	outPar << "BOUNDARIES \t";
	if (sim.absorbing) outPar << "absorbing" << endl;
	else outPar << "reflective" << endl;
	outPar << endl;

	outPar << "LANDSCAPE:\t";
	if (ppLand.generated) {
		outPar << "artificially generated map" << endl;
		outPar << "TYPE: \t";
		if (ppGenLand.continuous) outPar << "continuous \t";
		else outPar << "discrete \t";
		if (ppGenLand.fractal) outPar << "fractal";
		else outPar << "random";
		outPar << endl << "PROPORTION OF SUITABLE HABITAT (p)\t" << ppGenLand.propSuit << endl;
		if (ppGenLand.fractal) outPar << "HURST EXPONENT\t" << ppGenLand.hurst << endl;
	}
	else {
		outPar << "imported map" << endl;
		outPar << "TYPE: \t";
		switch (ppLand.rasterType) {
		case 0:
			outPar << "habitat codes" << endl;
			break;
		case 1:
			outPar << "habitat % cover" << endl;
			break;
		case 2:
			outPar << "habitat quality" << endl;
			break;
		}
		outPar << "FILE NAME: ";
#if RS_RCPP
		if (ppLand.dynamic) {
			outPar << name_landscape << endl;
		}
		else {
			outPar << name_landscape << endl;
		}
		if (ppLand.patchModel) {
			outPar << "PATCH FILE: " << name_patch << endl;
		}
		if (trfr.costMap) {
			outPar << "COSTS FILE: " << name_costfile << endl;
		}
#else
		if (sim.batchMode) outPar << " (see batch file) " << landFile << endl;
		else {
			outPar << habmapname << endl;
			if (ppLand.rasterType == 1) { // habitat % cover - list additional layers
				for (int i = 0; i < ppLand.nHab - 1; i++) {
					outPar << "           " << hfnames[i] << endl;
				}
			}
			if (ppLand.patchModel) {
				outPar << "PATCH FILE: " << patchmapname << endl;
			}
		}
#endif
		outPar << "No. HABITATS:\t" << ppLand.nHab << endl;
	}
	outPar << "RESOLUTION (m): \t" << ppLand.resol << endl;
	outPar << "DIMENSIONS:  X " << ppLand.dimX << "  Y " << ppLand.dimY << endl;
	outPar << "AVAILABLE:   min.X " << ppLand.minX << " min.Y " << ppLand.minY
		<< "  max.X " << ppLand.maxX << " max.Y " << ppLand.maxY << endl;
	if (!ppLand.generated && ppLand.dynamic) {
		landChange chg;
		outPar << "DYNAMIC LANDSCAPE: " << endl;
		int nchanges = pLandscape->numLandChanges();
		for (int i = 0; i < nchanges; i++) {
			chg = pLandscape->getLandChange(i);
			outPar << "Change no. " << chg.chgnum << " in year " << chg.chgyear << endl;
			outPar << "Landscape: " << chg.habfile << endl;
			if (ppLand.patchModel) {
				outPar << "Patches  : " << chg.pchfile << endl;
			}
			if (chg.costfile != "none" && chg.costfile != "NULL") {
				outPar << "Costs    : " << chg.costfile << endl;
			}
			//		outPar << "Change no. " << chg.chgnum << " in year " << chg.chgyear
			//			<< " habitat map: " << chg.habfile << endl;
		}
	}
	outPar << endl << "SPECIES DISTRIBUTION LOADED: \t";
	if (ppLand.spDist)
	{
		outPar << "yes" << endl;
		outPar << "RESOLUTION (m)\t" << ppLand.spResol << endl;
		outPar << "FILE NAME: ";
#if !RS_RCPP
		if (sim.batchMode) outPar << " (see batch file) " << landFile << endl;
		else {
			outPar << distnmapname << endl;
		}
#else
		outPar << name_sp_dist << endl;
#endif
	}
	else outPar << "no" << endl;

	outPar << endl << "ENVIRONMENTAL GRADIENT:\t ";
	if (grad.gradient)
	{
		switch (grad.gradType) {
		case 1:
			if (dem.stageStruct) outPar << "Density dependence strength (1/b)" << endl;
			else outPar << "Carrying capacity (K)" << endl;
			break;
		case 2:
			if (dem.stageStruct) outPar << "Fecundity" << endl;
			else outPar << "Intrinsic growth rate (r)" << endl;
			break;
		case 3:
			outPar << "Local extinction probability" << endl;
			break;
		default:
			outPar << "ERROR ERROR ERROR" << endl;
			;
		}
		outPar << "G:\t\t " << grad.grad_inc << endl;
		outPar << "optimum Y:\t " << grad.opt_y << endl;
		outPar << "f:\t\t " << grad.factor << endl;
		if (grad.gradType == 3) outPar << "Local extinction prob. at optimum:\t "
			<< grad.extProbOpt << endl;
		outPar << "GRADIENT SHIFTING:\t ";
		if (grad.shifting)
		{
			outPar << "yes" << endl;
			outPar << "SHIFTING RATE  (rows/year):\t " << grad.shift_rate << endl;
			outPar << "SHIFTING START (year):\t\t " << grad.shift_begin << endl;
			outPar << "SHIFTING STOP  (year):\t\t " << grad.shift_stop << endl;
		}
		else   outPar << "no" << endl;
	}
	else outPar << "no";
	outPar << endl;
	outPar << "ENVIRONMENTAL STOCHASTICITY:\t";
	if (env.stoch) {
		outPar << "yes" << endl;
		outPar << "TYPE\t in ";
		if (dem.stageStruct) {
			if (env.inK) outPar << "1/b" << endl;
			else outPar << "fecundity" << endl;
		}
		else {
			if (env.inK) outPar << "K" << endl;
			else outPar << "R" << endl;
		}
		outPar << "SPATIAL AUTOCORRELATION\t ";
		if (env.local) outPar << "local" << endl;
		else outPar << "global" << endl;
		outPar << "TEMPORAL AUTOCORRELATION (ac)\t" << env.ac << endl;
		outPar << "AMPLITUDE (std)\t" << env.std << endl;
		if (dem.stageStruct) {
			if (env.inK) {
				outPar << "MIN. 1/b\t" << pSpecies->getMinMax(0)
					* (10000.0 / (float)(ppLand.resol * ppLand.resol)) << endl;
				outPar << "MAX. 1/b\t" << pSpecies->getMinMax(1)
					* (10000.0 / (float)(ppLand.resol * ppLand.resol)) << endl;
			}
			else {
				outPar << "MIN. fecundity\t" << pSpecies->getMinMax(0) << endl;
				outPar << "MAX. fecundity\t" << pSpecies->getMinMax(1) << endl;
			}
		}
		else {
			if (env.inK) {
				outPar << "MIN. K\t" << pSpecies->getMinMax(0)
					* (10000.0 / (float)(ppLand.resol * ppLand.resol)) << endl;
				outPar << "MAX. K\t" << pSpecies->getMinMax(1)
					* (10000.0 / (float)(ppLand.resol * ppLand.resol)) << endl;
			}
			else {
				outPar << "MIN. r\t" << pSpecies->getMinMax(0) << endl;
				outPar << "MAX. r\t" << pSpecies->getMinMax(1) << endl;
			}
		}
	}
	else outPar << "no" << endl;
	outPar << "LOCAL EXTINCTION PROBABILITY:\t";
	if (env.localExt) outPar << env.locExtProb << endl;
	else outPar << "0.0" << endl;

	outPar << endl << "SPECIES' PARAMETERS." << endl;
	outPar << "REPRODUCTION:" << endl;
	outPar << "TYPE: ";
	switch (dem.repType) {
	case 0:
		outPar << "Asexual / Only female model" << endl;
		break;
	case 1:
		outPar << "Sexual model (simple)";
		outPar << endl;
		outPar << "PROP. of MALES\t" << dem.propMales << endl;
		break;
	case 2:
		outPar << "Sexual model (explicit mating system)" << endl;
		outPar << "PROP. of MALES\t" << dem.propMales << endl;
		outPar << "MAX. HAREM SIZE (h)\t" << dem.harem << endl;
		break;
	}
	outPar << "STAGE STRUCTURE:\t";
	if (dem.stageStruct) {
		outPar << "yes" << endl;
		outPar << "PROBABILITY OF REPRODUCING IN SUBSEQUENT SEASONS\t" << sstruct.probRep << endl;
		outPar << "No. OF REP. SEASONS BEFORE SUBSEQUENT REPRODUCTIONS\t" << sstruct.repInterval << endl;
		if (!ppLand.generated && ppLand.dynamic) {
			outPar << "ACTION AFTER POPULATION DESTRUCTION: all individuals ";
			if (sstruct.disperseOnLoss) outPar << "disperse" << endl;
			else outPar << "die" << endl;
		}
		outPar << "No. STAGES\t" << sstruct.nStages << endl;
		outPar << "MAX. AGE\t" << sstruct.maxAge << endl;
		// no sex-specific demographic parameters
		if (dem.repType != 2) {
			outPar << "MIN. AGES:" << endl;
			for (int i = 0; i < sstruct.nStages; i++) {
				outPar << "stage\t" << i << ":\t" << pSpecies->getMinAge(i, 0) << "\tyears" << endl;
			}
			outPar << "FECUNDITIES:" << endl;
			for (int i = 0; i < sstruct.nStages; i++) {
				outPar << "stage\t" << i << ":\t" << pSpecies->getFec(i, 0) << endl;
			}
			outPar << "DEVELOPMENT PROB.:" << endl;
			for (int i = 0; i < sstruct.nStages; i++) {
				outPar << "stage\t" << i << ":\t" << pSpecies->getDev(i, 0) << endl;
			}
			outPar << "SURVIVAL PROB.:" << endl;
			for (int i = 0; i < sstruct.nStages; i++) {
				outPar << "stage\t" << i << ":\t" << pSpecies->getSurv(i, 0) << endl;
			}
		}
		// sex-specific demographic parameters
		else {
			outPar << "MIN. AGES:" << endl;
			for (int i = 0; i < sstruct.nStages; i++) {
				outPar << "males " << i << ":\t" << pSpecies->getMinAge(i, 1) << " years;\t";
				outPar << "females " << i << ":\t" << pSpecies->getMinAge(i, 0) << " years" << endl;
			}
			outPar << "FECUNDITIES:" << endl;
			for (int i = 0; i < sstruct.nStages; i++) {
				outPar << "males   " << i << ":\t" << pSpecies->getFec(i, 1) << endl;
				outPar << "females " << i << ":\t" << pSpecies->getFec(i, 0) << endl;
			}
			outPar << "DEVELOPMENT PROB.:" << endl;
			for (int i = 0; i < sstruct.nStages; i++) {
				outPar << "males   " << i << ":\t" << pSpecies->getDev(i, 1) << endl;
				outPar << "females " << i << ":\t" << pSpecies->getDev(i, 0) << endl;
			}
			outPar << "SURVIVAL PROB.:" << endl;
			for (int i = 0; i < sstruct.nStages; i++) {
				outPar << "males   " << i << ":\t" << pSpecies->getSurv(i, 1) << endl;
				outPar << "females " << i << ":\t" << pSpecies->getSurv(i, 0) << endl;
			}
		}

		outPar << "SCHEDULING OF SURVIVAL: ";
		switch (sstruct.survival) {
		case 0:
			outPar << "At reproduction" << endl;
			break;
		case 1:
			outPar << "Between reproductive events" << endl;
			break;
		case 2:
			outPar << "Annually" << endl;
			break;
		}

		int mSize; // index for weights matrices
		if (dem.repType == 2) mSize = sstruct.nStages * gMaxNbSexes;
		else mSize = sstruct.nStages;

		outPar << "DENSITY-DEPENDENCE IN FECUNDITY:\t";
		if (sstruct.fecDens) {
			outPar << "yes" << endl;
			if (sstruct.fecStageDens) {
				outPar << "STAGE'S WEIGHTS:" << endl;
				for (int i = 0; i < mSize; i++) {
					if (dem.repType == 2) {
						outPar << "stage " << i / gMaxNbSexes << " ";
						if (i % gMaxNbSexes == 0) outPar << "males  : \t";
						else outPar << "females: \t";
					}
					else outPar << "stage " << i << ": \t";
					for (int j = 0; j < mSize; j++) outPar << pSpecies->getDDwtFec(j, i) << "\t";
					outPar << endl;
				}
			}
			else outPar << "not stage-dependent" << endl;
		}
		else outPar << "no" << endl;

		densDepParams ddparams = pSpecies->getDensDep();

		outPar << "DENSITY-DEPENDENCE IN DEVELOPMENT:\t";
		if (sstruct.devDens) {
			outPar << "yes - coefficient: " << ddparams.devCoeff << endl;
			if (sstruct.devStageDens) {
				outPar << "STAGE'S WEIGHTS:" << endl;
				for (int i = 0; i < mSize; i++) {
					if (dem.repType == 2) {
						outPar << "stage " << i / gMaxNbSexes << " ";
						if (i % gMaxNbSexes == 0) outPar << "males  : \t";
						else outPar << "females: \t";
					}
					else outPar << "stage " << i << ": \t";
					for (int j = 0; j < mSize; j++) outPar << pSpecies->getDDwtDev(j, i) << "\t";
					outPar << endl;
				}
			}
			else outPar << "not stage-dependent" << endl;
		}
		else outPar << "no" << endl;

		outPar << "DENSITY-DEPENDENCE IN SURVIVAL:\t\t";
		if (sstruct.survDens) {
			outPar << "yes - coefficient: " << ddparams.survCoeff << endl;
			if (sstruct.survStageDens) {
				outPar << "STAGE'S WEIGHTS:" << endl;
				for (int i = 0; i < mSize; i++) {
					if (dem.repType == 2) {
						outPar << "stage " << i / gMaxNbSexes << " ";
						if (i % gMaxNbSexes == 0) outPar << "males  : \t";
						else outPar << "females: \t";
					}
					else outPar << "stage " << i << ": \t";
					for (int j = 0; j < mSize; j++) outPar << pSpecies->getDDwtSurv(j, i) << "\t";
					outPar << endl;
				}
			}
			else outPar << "not stage-dependent" << endl;
		}
		else outPar << "no" << endl;
	} // end of if (dem.stageStruct)
	else { // not stage-strutured
		outPar << "no" << endl;
		outPar << "Rmax\t" << dem.lambda << endl;
		outPar << "bc\t" << dem.bc << endl;
	}

	if (dem.stageStruct) {
		outPar << endl << "HABITAT SPECIFIC 1/b:" << endl;
	}
	else {
		outPar << endl << "CARRYING CAPACITIES:" << endl;
	}
	int nhab = ppLand.nHab;
	if (ppLand.generated) {
		if (ppGenLand.continuous) nhab = 1;
	}
	for (int i = 0; i < nhab; i++) {
		k = pSpecies->getHabK(i) * (10000.0 / (float)(ppLand.resol * ppLand.resol));
		if (!ppLand.generated && ppLand.rasterType == 0) { // imported & habitat codes
			outPar << "Habitat " << pLandscape->getHabCode(i) << ": \t";
		}
		else {
			outPar << "Habitat " << i << ": ";
		}
		if (dem.stageStruct) outPar << "1/b ";
		else outPar << "K ";
		outPar << k << endl;
	}
	emigTraits ep0, ep1;
	string sexdept = "SEX-DEPENDENT:   ";
	string stgdept = "STAGE-DEPENDENT: ";
	string indvar = "INDIVIDUAL VARIABILITY: ";
	string emigstage = "EMIGRATION STAGE: ";

	outPar << endl << "DISPERSAL - EMIGRATION:\t";
	if (emig.densDep) {
		outPar << "density-dependent" << endl;
		if (emig.sexDep) {
			outPar << sexdept << "yes" << endl;
			if (emig.stgDep) {
				outPar << stgdept << "yes" << endl;
				outPar << indvar << "no" << endl;
				for (int i = 0; i < sstruct.nStages; i++) {
					outPar << "stage " << i << ":" << endl;
					ep0 = pSpecies->getSpEmigTraits(i, 0);
					ep1 = pSpecies->getSpEmigTraits(i, 1);
					outPar << "D0:    females " << ep0.d0 << "  males " << ep1.d0 << endl;
					outPar << "alpha: females " << ep0.alpha << "  males " << ep1.alpha << endl;
					outPar << "beta:  females " << ep0.beta << "  males " << ep1.beta << endl;
				}
			}
			else { // !emig.stgDep
				outPar << stgdept << "no" << endl;
				ep0 = pSpecies->getSpEmigTraits(0, 0);
				ep1 = pSpecies->getSpEmigTraits(0, 1);
				outPar << "D0:    females " << ep0.d0 << "  males " << ep1.d0 << endl;
				outPar << "alpha: females " << ep0.alpha << "  males " << ep1.alpha << endl;
				outPar << "beta:  females " << ep0.beta << "  males " << ep1.beta << endl;
		}
	}
		else { // !emig.sexDep
			outPar << sexdept << "no" << endl;
			if (emig.stgDep) {
				outPar << stgdept << "yes" << endl;
				outPar << indvar << "no" << endl;
				for (int i = 0; i < sstruct.nStages; i++) {
					ep0 = pSpecies->getSpEmigTraits(i, 0);
					outPar << "stage " << i << ": \t" << "D0: " << ep0.d0;
					outPar << " \talpha: " << ep0.alpha << " \tbeta: " << ep0.beta << endl;
				}
			}
			else { // !emig.stgDep
				outPar << stgdept << "no" << endl;
				ep0 = pSpecies->getSpEmigTraits(0, 0);
				outPar << "D0:    " << ep0.d0 << endl;
				outPar << "alpha: " << ep0.alpha << endl;
				outPar << "beta:  " << ep0.beta << endl;
			}
		}
	}
	else { // not density-dependent
		string initprob = "INITIAL EMIGRATION PROB. ";
		outPar << "density-independent" << endl;
		if (!trfr.usesMovtProc) { // transfer by kernel
			outPar << "USE FULL KERNEL TO DETERMINE EMIGRATION: ";
			if (pSpecies->useFullKernel()) outPar << "yes";
			else outPar << "no";
			outPar << endl;
		}

		if (emig.sexDep) {
			outPar << sexdept << "yes" << endl;
			if (emig.stgDep) {
				outPar << stgdept << "yes" << endl;
				outPar << indvar << "no" << endl;
				for (int i = 0; i < sstruct.nStages; i++) {
					outPar << "stage " << i << ": \t" << "EMIGRATION PROB.: \tfemales "
						<< pSpecies->getSpEmigD0(i, 0) << " \tmales " << pSpecies->getSpEmigD0(i, 1) << endl;
				}
			}
			else { // !emig.stgDep
				outPar << stgdept << "no" << endl;
				outPar << "EMIGRATION PROB.: \tfemales " << pSpecies->getSpEmigD0(0, 0)
					<< "\t males " << pSpecies->getSpEmigD0(0, 1) << endl;
			}
		}
		else { // !emig.sexDep
			outPar << sexdept << "no" << endl;
			if (emig.stgDep) {
				outPar << stgdept << "yes" << endl;
				outPar << indvar << "no" << endl;
				for (int i = 0; i < sstruct.nStages; i++) {
					outPar << "stage " << i << ": \t" << "EMIGRATION PROB.: "
						<< pSpecies->getSpEmigD0(i, 0) << endl;
				}
			}
			else { // !emig.stgDep
				outPar << stgdept << "no" << endl;
				outPar << "EMIGRATION PROB.:\t" << pSpecies->getSpEmigD0(0, 0) << endl;
			}
		}
	}

	// Transfer

	outPar << endl << "DISPERSAL - TRANSFER: \t";

	if (trfr.usesMovtProc) {
		bool straightenPath;
		if (trfr.moveType == 1) { // SMS
			trfrSMSTraits move = pSpecies->getSpSMSTraits();
			straightenPath = move.straightenPath;
			if (trfr.costMap) {
				outPar << "SMS\tcosts from imported cost map" << endl;
#if !RS_RCPP
				outPar << "FILE NAME: " << costmapname << endl;
#endif
			}
			else {
				outPar << "SMS\tcosts:" << endl;
				if (!ppLand.generated && ppLand.rasterType == 0) {
					for (int i = 0; i < ppLand.nHab; i++)
						outPar << "\thab. " << pLandscape->getHabCode(i) << "\t"
						<< pSpecies->getHabCost(i) << endl;
				}
				else {
					for (int i = 0; i < ppLand.nHab; i++)
						outPar << "\thab. " << i << "\t"
						<< pSpecies->getHabCost(i) << endl;
				}
			}
			string pr = "PERCEPTUAL RANGE";
			outPar << pr << ":        " << move.pr << endl;
			outPar << pr << " METHOD: " << move.prMethod << endl;
			if (!trfr.indVar) outPar << "DIRECTIONAL PERSISTENCE: " << move.dp << endl;
			outPar << "MEMORY SIZE: " << move.memSize << endl;
			outPar << "GOAL TYPE:   " << move.goalType << endl;
			if (!trfr.indVar) {
				if (move.goalType == 2) { //  dispersal bias
					outPar << "GOAL BIAS:   " << move.gb << endl;
					outPar << "ALPHA DB:    " << move.alphaDB << endl;
					outPar << "BETA DB:     " << move.betaDB << endl;
				}
			}
			outPar << indvar << "no " << endl;
		}
		else { // CRW
			trfrCRWTraits move = pSpecies->getSpCRWTraits();
			straightenPath = move.straightenPath;
			outPar << "CRW" << endl;
			string lgth = "STEP LENGTH (m) ";
			string corr = "STEP CORRELATION";
			outPar << lgth << ": " << move.stepLength << endl;
			outPar << corr << ": " << move.rho << endl;
		}
		outPar << "STRAIGHTEN PATH AFTER DECISION NOT TO SETTLE: ";
		if (straightenPath) outPar << "yes" << endl;
		else outPar << "no" << endl;
		outPar << "STEP MORTALITY:\t" << endl;
		if (trfr.habMort)
		{
			outPar << "habitat dependent:\t" << endl;
			if (!ppLand.generated && ppLand.rasterType == 0) {
				for (int i = 0; i < ppLand.nHab; i++)
					outPar << "\thab. " << pLandscape->getHabCode(i) << "\t"
					<< pSpecies->getHabMort(i) << endl;
			}
			else {
				for (int i = 0; i < ppLand.nHab; i++)
					outPar << "\thab. " << i << "\t"
					<< pSpecies->getHabMort(i) << endl;
			}
		}
		else
		{
			trfrCRWTraits move = pSpecies->getSpCRWTraits();
			outPar << "constant " << move.stepMort << endl;
		}
	} // end of movement process
	else { // kernel
		string meandist = "MEAN DISTANCE";
		string probkern = "PROB. KERNEL I";
		trfrKernelParams kern0, kern1;
		outPar << "dispersal kernel" << endl << "TYPE: \t";
		if (trfr.twinKern) outPar << "double ";
		outPar << "negative exponential" << endl;

		if (trfr.sexDep) {
			outPar << sexdept << "yes" << endl;
			if (trfr.stgDep) {
				outPar << stgdept << "yes" << endl;
				outPar << indvar << "no" << endl;
				for (int i = 0; i < sstruct.nStages; i++) {
					outPar << "stage " << i << ":" << endl;
					kern0 = pSpecies->getSpKernTraits(i, 0);
					kern1 = pSpecies->getSpKernTraits(i, 1);
					outPar << meandist << " I: \tfemales " << kern0.meanDist1 << " \tmales " << kern1.meanDist1 << endl;
					if (trfr.twinKern)
					{
						outPar << meandist << " II: \tfemales " << kern0.meanDist2 << " \tmales " << kern1.meanDist2 << endl;
						outPar << probkern << ": \tfemales " << kern0.probKern1 << " \tmales " << kern1.probKern1 << endl;
					}
				}
			}
			else { // !trfr.stgDep
				outPar << stgdept << "no" << endl;
				kern0 = pSpecies->getSpKernTraits(0, 0);
				kern1 = pSpecies->getSpKernTraits(0, 1);
				outPar << meandist << " I: \tfemales " << kern0.meanDist1 << " \tmales " << kern1.meanDist1 << endl;
				if (trfr.twinKern)
				{
					outPar << meandist << " II: \tfemales " << kern0.meanDist2 << " \tmales " << kern1.meanDist2 << endl;
					outPar << probkern << ": \tfemales " << kern0.probKern1 << " \tmales " << kern1.probKern1 << endl;
				}
			}
		}
		else { // !trfr.sexDep
			outPar << sexdept << "no" << endl;
			if (trfr.stgDep) {
				outPar << stgdept << "yes" << endl;
				outPar << indvar << "no" << endl;
				for (int i = 0; i < sstruct.nStages; i++) {
					kern0 = pSpecies->getSpKernTraits(i, 0);
					outPar << "stage " << i << ": \t" << meandist << " I: " << kern0.meanDist1;
					if (trfr.twinKern)
					{
						outPar << " \t" << meandist << " II: " << kern0.meanDist2;
						outPar << " \t" << probkern << ": " << kern0.probKern1;
					}
					outPar << endl;
				}
			}
			else { // !trfr.stgDep
				outPar << stgdept << "no" << endl;
				kern0 = pSpecies->getSpKernTraits(0, 0);
				outPar << meandist << " I: \t" << kern0.meanDist1 << endl;
				if (trfr.twinKern)
				{
					outPar << meandist << " II: \t" << kern0.meanDist2 << endl;
					outPar << probkern << ": \t" << kern0.probKern1 << endl;
				}
			}
		}

		outPar << "DISPERSAL MORTALITY:   ";
		trfrMortParams mort = pSpecies->getMortParams();
		if (trfr.distMort) {
			outPar << "distance-dependent" << endl;
			outPar << "SLOPE: " << mort.mortAlpha << " \tINFLECTION POINT: " << mort.mortBeta << endl;
		}
		else {
			outPar << "constant" << endl << "MORTALITY PROBABILITY: " << mort.fixedMort << endl;
		}
	} // end of kernel transfer

	// Settlement

	outPar << endl << "DISPERSAL - SETTLEMENT:" << endl;

	if (trfr.usesMovtProc) {
		string plusmating = "+ mating requirements";

		if (sett.sexDep) {
			nsexes = 2;
			outPar << sexdept << "yes" << endl;
			if (sett.stgDep) {
				nstages = sstruct.nStages;
				outPar << stgdept << "yes" << endl;
				for (int i = 0; i < nstages; i++) {
				    if (dem.stageStruct && nstages > 1) outPar << "stage " << i << ": " << endl;
				    for (int sx = 0; sx < nsexes; sx++) {
				        if (sx == 0) outPar << "FEMALES:" << endl;
				        else outPar << "MALES:" << endl;
				        ssteps = pSpecies->getSteps(i, sx);

				        outPar << "MIN. No. OF STEPS:\t " << ssteps.minSteps << endl;
				        outPar << "MAX. No. OF STEPS:\t ";
				        if (ssteps.maxSteps == 99999999) outPar << "not applied" << endl;
				        else outPar << ssteps.maxSteps << endl;
				    }
				}
			}
			else { // !sett.stgDep
				nstages = 1;
				outPar << stgdept << "no" << endl;
				for (int sx = 0; sx < nsexes; sx++) {
				    if (sx == 0) outPar << "FEMALES:" << endl;
				    else outPar << "MALES:" << endl;
				    ssteps = pSpecies->getSteps(0, sx);

				    outPar << "MIN. No. OF STEPS:\t " << ssteps.minSteps << endl;
				    outPar << "MAX. No. OF STEPS:\t ";
				    if (ssteps.maxSteps == 99999999) outPar << "not applied" << endl;
				    else outPar << ssteps.maxSteps << endl;
				}
			}
		}
		else { // !sett.sexDep
			nsexes = 1;
			outPar << sexdept << "no" << endl;
			if (sett.stgDep) {
				nstages = sstruct.nStages;
				outPar << stgdept << "yes" << endl;
				for (int i = 0; i < nstages; i++) {
				    if (dem.stageStruct && nstages > 1) outPar << "stage " << i << ": " << endl;
				    ssteps = pSpecies->getSteps(i, 0);

				    outPar << "MIN. No. OF STEPS:\t " << ssteps.minSteps << endl;
				    outPar << "MAX. No. OF STEPS:\t ";
				    if (ssteps.maxSteps == 99999999) outPar << "not applied" << endl;
				    else outPar << ssteps.maxSteps << endl;
				}
			}
			else { // !sett.stgDep
				nstages = 1;
				outPar << stgdept << "no" << endl;
				ssteps = pSpecies->getSteps(0, 0);

				outPar << "MIN. No. OF STEPS:\t " << ssteps.minSteps << endl;
				outPar << "MAX. No. OF STEPS:\t ";
				if (ssteps.maxSteps == 99999999) outPar << "not applied" << endl;
				else outPar << ssteps.maxSteps << endl;
			}
		}
		for (int sx = 0; sx < nsexes; sx++) {
			if (sett.sexDep) {
				if (sx == 0) outPar << "FEMALES:" << endl;
				else outPar << "MALES:" << endl;
			}
			outPar << "SETTLE IF: ";
			for (int i = 0; i < nstages; i++) {
				if (dem.stageStruct && nstages > 1) outPar << "stage " << i << ": " << endl;
				outPar << "find a suitable cell/patch ";
				srules = pSpecies->getSettRules(i, sx);
				if (srules.densDep) {
					settleDD = pSpecies->getSpSettTraits(i, sx);
					outPar << "+ density dependence ";
					if (srules.findMate) outPar << plusmating;
					outPar << endl;
					if (!sett.indVar) {
						outPar << "S0: " << settleDD.s0 << "  AlphaS: " << settleDD.alpha
							<< "  BetaS: " << settleDD.beta << endl;
					}
				}
				else {
					if (srules.findMate) outPar << plusmating << endl;
					else outPar << "(not the natal one)" << endl;
				}
				if (dem.stageStruct) {
					ssteps = pSpecies->getSteps(i, sx);
					outPar << "MAX. No. OF STEPS/YEAR:\t ";
					if (ssteps.maxStepsYr == 99999999) outPar << "not applied" << endl;
					else outPar << ssteps.maxStepsYr << endl;
				}
			}
		}
	}
	else { // kernel-based transfer
		string notsuit = "IF THE ARRIVAL CELL/PATCH IS UNSUITABLE: ";
		string rchoose = " randomly choose a suitable neighb. cell/patch or ";
		string matereq = "MATING REQUIREMENTS: ";
		if (sett.sexDep) {
			nsexes = 2;
			outPar << sexdept << "yes" << endl;
			if (sett.stgDep) {
				nstages = sstruct.nStages;
				outPar << stgdept << "yes" << endl;
				outPar << notsuit << endl;
			}
			else {
				nstages = 1;
				outPar << stgdept << "no" << endl;
			}
		}
		else {
			nsexes = 1;
			outPar << sexdept << "no" << endl;
			if (sett.stgDep) {
				nstages = sstruct.nStages;
				outPar << stgdept << "yes" << endl;
				outPar << notsuit << endl;
			}
			else {
				nstages = 1;
				outPar << stgdept << "no" << endl;
				outPar << notsuit;
			}
		}
		for (int i = 0; i < nstages; i++) {
			if (sett.stgDep) {
				outPar << "stage " << i << ":" << endl;
			}
			for (int sx = 0; sx < nsexes; sx++) {
				if (sett.sexDep) {
					if (sx == 0) outPar << "FEMALES: ";
					else outPar << "MALES:   ";
					if (!sett.stgDep) outPar << notsuit;
				}
				srules = pSpecies->getSettRules(i, sx);
				if (srules.go2nbrLocn) {
					outPar << rchoose;
					if (srules.wait) outPar << "wait" << endl;
					else outPar << "die" << endl;
				}
				else {
					if (srules.wait) outPar << "wait" << endl;
					else outPar << "die" << endl;
				}
				outPar << matereq;
				if (srules.findMate) outPar << "yes" << endl;
				else outPar << "no" << endl;
			}
		}
	}

	// Genetics
	outPar << endl << "GENETICS:" << endl;
	set<TraitType> traitList = pSpecies->getTraitTypes();

	if (pSpecies->isDiploid()) outPar << "DIPLOID" << endl; else outPar << "HAPLOID" << endl;
	outPar << "Genome size: " << pSpecies->getGenomeSize() << endl;
	outPar << "Chromosome breaks : ";

	for (auto end : pSpecies->getChromosomeEnds())
		outPar << end << " ";
	outPar << endl;
	outPar << "Recombination rate: " << pSpecies->getRecombinationRate() << endl;
	outPar << "Traits modelled:  " << endl;
	for (auto trait : traitList)
		outPar << trait << endl;

	// Initialisation

	initParams init = paramsInit->getInit();
	outPar << endl << "INITIALISATION CONDITIONS:" << endl;
	switch (init.seedType) {
	case 0:
		outPar << "Free initialisation: \t";
		switch (init.freeType) {
		case 0:
			outPar << "Random \t";
			outPar << "No. of cells/patches: " << init.nSeedPatches << endl;
			break;
		case 1:
			outPar << "all suitable cells/patches" << endl;
			break;
		case 2:
			outPar << "manually selected cells/patches" << endl;
			break;
		}
		break;
	case 1:
		outPar << "From species distribution: \t" << endl;
		switch (init.spDistType) {
		case 0:
			outPar << "all presence cells/patches" << endl;
			break;
		case 1:
			outPar << "some random presence cells/patches" << endl;
			break;
		case 2:
			outPar << "all cells/patches within selected distribution cells" << endl;
			break;
		}
		break;
	case 2:
		outPar << "From initial individuals file: " << paramsSim->getDir(1) + init.indsFile << endl;
		break;
	case 3:
		outPar << "From file" << endl;
		break;
	}
	if (init.seedType != 2) {
		outPar << "INITIAL NO. OF INDIVIDUALS: \t";
		switch (init.initDens) {
		case 0:
			outPar << "at carrying capacity" << endl;
			break;
		case 1:
			outPar << "at half carrying capacity" << endl;
			break;
		case 2:
			if (ppLand.patchModel) {
				outPar << init.indsHa << " individuals per ha" << endl;
			}
			else {
				outPar << init.indsCell << " individuals per cell" << endl;
			}
			break;
		}
		if (dem.stageStruct) {
			outPar << "INITIAL STAGE PROPORTIONS:" << endl;
			for (int i = 1; i < sstruct.nStages; i++) {
				outPar << "stage " << i << ": " << paramsInit->getProp(i) << " \t";
			}
			outPar << endl;
			outPar << "Initial age distribution: ";
			switch (init.initAge) {
			case 0:
				outPar << "lowest possible age";
				break;
			case 1:
				outPar << "randomised";
				break;
			case 2:
				outPar << "quasi-equilibrium";
				break;
			}
			outPar << endl;
		}
		outPar << "GEOGRAPHICAL CONSTRAINTS (cell numbers): " << endl;
		outPar << "min X: " << init.minSeedX << " max X: " << init.maxSeedX << endl;
		outPar << "min Y: " << init.minSeedY << " max Y: " << init.maxSeedY << endl;
		//	if (init.seedType != 1 && init.freeType < 2 && init.initFrzYr > 0) {
		//		outPar << "Freeze initial range until year " << init.initFrzYr << endl;
		//	}
		if (init.seedType == 0 && init.freeType < 2) {
			if (init.initFrzYr > 0) {
				outPar << "Freeze initial range until year " << init.initFrzYr << endl;
			}
			if (init.restrictRange) {
				outPar << "Restrict range to northern " << init.restrictRows
					<< " rows every " << init.restrictFreq << " years" << endl;
				if (init.finalFrzYr < sim.years) {
					outPar << "Freeze range at year " << init.finalFrzYr << endl;
				}
			}
		}
	}

	outPar << endl << "OUTPUTS:" << endl;
	if (sim.outRange) {
		outPar << "Range - every " << sim.outIntRange << " year";
		if (sim.outIntRange > 1) outPar << "s";
		//	if (sim.outStartRange > 0) outPar << " starting year " << sim.outStartRange;
		outPar << endl;
	}
	if (sim.outOccup) {
		outPar << "Occupancy - every " << sim.outIntOcc << " year";
		if (sim.outIntOcc > 1) outPar << "s";
		//	if (sim.outStartOcc > 0) outPar << " starting year " << sim.outStartOcc;
		outPar << endl;
	}
	if (sim.outPop) {
		outPar << "Populations - every " << sim.outIntPop << " year";
		if (sim.outIntPop > 1) outPar << "s";
		if (sim.outStartPop > 0) outPar << " starting year " << sim.outStartPop;
		outPar << endl;
	}
	if (sim.outInds) {
		outPar << "Individuals - every " << sim.outIntInd << " year";
		if (sim.outIntInd > 1) outPar << "s";
		if (sim.outStartInd > 0) outPar << " starting year " << sim.outStartInd;
		outPar << endl;
	}
	if (sim.outputWeirCockerham || sim.outputWeirHill) {
		outPar << "Neutral genetics - every " << sim.outputGeneticInterval << " year";
		if (sim.outputGeneticInterval > 1) outPar << "s";
		if (sim.outputWeirHill) outPar << " outputting pairwise patch fst";
		if (sim.outputWeirCockerham) outPar << " outputting per locus fst ";
		outPar << endl;
	}

	if (sim.outTraitsCells) {
		outPar << "Traits per ";
		if (ppLand.patchModel) outPar << "patch"; else outPar << "cell";
		outPar << " - every " << sim.outIntTraitCell << " year";
		if (sim.outIntTraitCell > 1) outPar << "s";
		if (sim.outStartTraitCell > 0) outPar << " starting year " << sim.outStartTraitCell;
		outPar << endl;
	}
	if (sim.outTraitsRows) {
		outPar << "Traits per row - every " << sim.outIntTraitRow << " year";
		if (sim.outIntTraitRow > 1) outPar << "s";
		if (sim.outStartTraitRow > 0) outPar << " starting year " << sim.outStartTraitRow;
		outPar << endl;
	}
	if (sim.outConnect) {
		outPar << "Connectivity matrix - every " << sim.outIntConn << " year";
		if (sim.outIntConn > 1) outPar << "s";
		if (sim.outStartConn > 0) outPar << " starting year " << sim.outStartConn;
		outPar << endl;
	}
#if RS_RCPP
	if (sim.outPaths) {
		outPar << "SMS paths - every " << sim.outIntPaths << " year";
		if (sim.outIntPaths > 1) outPar << "s";
		if (sim.outStartPaths > 0) outPar << " starting year " << sim.outStartPaths;
		outPar << endl;
	}
#endif
	outPar << "SAVE MAPS: ";
	if (sim.saveMaps) {
		outPar << "yes - every " << sim.mapInt << " year";
		if (sim.mapInt > 1) outPar << "s";
		outPar << endl;
	}
	else outPar << "no" << endl;
	outPar << "SAVE TRAITS MAPS: ";
	if (sim.saveTraitMaps) {
		outPar << "yes - every " << sim.traitInt << " year";
		if (sim.traitInt > 1) outPar << "s";
		outPar << endl;
	}
	else outPar << "no" << endl;
	if (trfr.usesMovtProc && trfr.moveType == 1) {
		outPar << "SMS HEAT MAPS: ";
		if (sim.saveVisits) outPar << "yes" << endl;
		else outPar << "no" << endl;
	}
	outPar.close(); outPar.clear();
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
