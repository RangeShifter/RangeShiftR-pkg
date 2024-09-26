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

#include "SubCommunity.h"
//---------------------------------------------------------------------------

ofstream outtraits;

//---------------------------------------------------------------------------

SubCommunity::SubCommunity(Patch* pPch, int num) {
	subCommNum = num;
	pPatch = pPch;
	// record the new sub-community no. in the patch
	pPatch->setSubComm((intptr)this);
	initial = false;
	occupancy = 0;
}

SubCommunity::~SubCommunity() {
	pPatch->setSubComm(0);
	int npops = (int)popns.size();
	for (int i = 0; i < npops; i++) { // all populations
		delete popns[i];
	}
	popns.clear();
	if (occupancy != 0) delete[] occupancy;
}

intptr SubCommunity::getNum(void) { return subCommNum; }

Patch* SubCommunity::getPatch(void) { return pPatch; }

locn SubCommunity::getLocn(void) {
	locn loc = pPatch->getCellLocn(0);
	return loc;
}

void SubCommunity::setInitial(bool b) { initial = b; }

void SubCommunity::initialise(Landscape* pLandscape, Species* pSpecies)
{
	int ncells;
	landParams ppLand = pLandscape->getLandParams();
	initParams init = paramsInit->getInit();

	// determine size of initial population
	int nInds = 0;
	if (subCommNum == 0 // matrix patch
		|| !initial)   		// not in initial region or distribution
		nInds = 0;
	else {
		float k = pPatch->getK();
		if (k > 0.0) { // patch is currently suitable for this species
			switch (init.initDens) {
			case 0: // at carrying capacity
				nInds = (int)k;
				break;
			case 1: // at half carrying capacity
				nInds = (int)(k / 2.0);
				break;
			case 2: // specified no. per cell or density
				ncells = pPatch->getNCells();
				if (ppLand.patchModel) {
					nInds = (int)(init.indsHa * (float)(ncells * ppLand.resol * ppLand.resol) / 10000.0);
				}
				else {
					nInds = init.indsCell * ncells;
				}
				break;
			}
		}
		else nInds = 0;
	}

	// create new population only if it is non-zero or the matrix popn
	if (subCommNum == 0 || nInds > 0) {
		newPopn(pLandscape, pSpecies, pPatch, nInds);
	}

}

// initialise a specified individual
void SubCommunity::initialInd(Landscape* pLandscape, Species* pSpecies,
	Patch* pPatch, Cell* pCell, int ix)
{
	demogrParams dem = pSpecies->getDemogrParams();
	stageParams sstruct = pSpecies->getStageParams();
	emigRules emig = pSpecies->getEmigRules();
	transferRules trfr = pSpecies->getTransferRules();
	settleType sett = pSpecies->getSettle();
	short stg, age, repInt;
	Individual* pInd;
	float probmale;

	// create new population if not already in existence
	int npopns = (int)popns.size();
	if (npopns < 1) {
		newPopn(pLandscape, pSpecies, pPatch, 0);
	}

	// create new individual
	initInd iind = paramsInit->getInitInd(ix);
	if (dem.stageStruct) {
		stg = iind.stage; age = iind.age;   repInt = sstruct.repInterval;
	}
	else {
		age = stg = 1;   repInt = 0;
	}
	if (dem.repType == 0) {
		probmale = 0.0;
	}
	else {
		if (iind.sex == 1) probmale = 1.0; else probmale = 0.0;
	}
	pInd = new Individual(pCell, pPatch, stg, age, repInt, probmale, trfr.usesMovtProc, trfr.moveType);

	// add new individual to the population
	// NB THIS WILL NEED TO BE CHANGED FOR MULTIPLE SPECIES...
	popns[0]->recruit(pInd);

	if (pSpecies->getNTraits() > 0)
	{
		// individual variation - set up genetics
		landData land = pLandscape->getLandData();
		pInd->setUpGenes(pSpecies, land.resol);
	}

}

// Create a new population, and return its address
Population* SubCommunity::newPopn(Landscape* pLandscape, Species* pSpecies,
	Patch* pPatch, int nInds)
{
	landParams land = pLandscape->getLandParams();
	int npopns = (int)popns.size();
	popns.push_back(new Population(pSpecies, pPatch, nInds, land.resol));
	return popns[npopns];
}

popStats SubCommunity::getPopStats(void) {
	popStats p, pop;
	p.pSpecies = 0; p.spNum = 0; p.nInds = p.nAdults = p.nNonJuvs = 0; p.breeding = false;
	p.pPatch = pPatch;
	// FOR SINGLE SPECIES IMPLEMENTATION, THERE IS ONLY ONE POPULATION IN THE PATCH
	int npops = (int)popns.size();
	for (int i = 0; i < npops; i++) { // all populations
		pop = popns[i]->getStats();
		p.pSpecies = pop.pSpecies;
		p.spNum = pop.spNum;
		p.nInds += pop.nInds;
		p.nNonJuvs += pop.nNonJuvs;
		p.nAdults += pop.nAdults;
		p.breeding = pop.breeding;
	}
	return p;
}

void SubCommunity::resetPopns(void) {
	int npops = (int)popns.size();
	for (int i = 0; i < npops; i++) { // all populations
		delete popns[i];
	}
	popns.clear();
	// clear the list of populations in the corresponding patch
	pPatch->resetPopn();
}

void SubCommunity::resetPossSettlers(void) {
	if (subCommNum == 0) return; // not applicable in the matrix
	pPatch->resetPossSettlers();
}

// Extirpate all populations according to
// option 0 - random local extinction probability
// option 1 - local extinction probability gradient
// NB only applied for cell-based model
void SubCommunity::localExtinction(int option) {
	double pExtinct = 0.0;
	if (option == 0) {
		envStochParams env = paramsStoch->getStoch();
		if (env.localExt) pExtinct = env.locExtProb;
	}
	else {
		envGradParams grad = paramsGrad->getGradient();
		Cell* pCell = pPatch->getRandomCell(); // get only cell in the patch
		// extinction prob is complement of cell gradient value plus any non-zero prob at the optimum
		pExtinct = 1.0 - pCell->getEnvVal() + grad.extProbOpt;
		if (pExtinct > 1.0) pExtinct = 1.0;
	}
	if (pRandom->Bernoulli(pExtinct)) {
		int npops = (int)popns.size();
		for (int i = 0; i < npops; i++) { // all populations
			popns[i]->extirpate();
		}
	}
}

// Action in event of patch becoming unsuitable owing to landscape change
void SubCommunity::patchChange(void) {
	if (subCommNum == 0) return; // no reproduction in the matrix
	Species* pSpecies;
	float localK = 0.0;
	int npops = (int)popns.size();
	// THE FOLLOWING MAY BE MORE EFFICIENT WHILST THERE IS ONLY ONE SPECIES ...
	if (npops < 1) return;
	localK = pPatch->getK();
	if (localK <= 0.0) { // patch in dynamic landscape has become unsuitable
		for (int i = 0; i < npops; i++) { // all populations
			pSpecies = popns[i]->getSpecies();
			demogrParams dem = pSpecies->getDemogrParams();
			if (dem.stageStruct) {
				stageParams sstruct = pSpecies->getStageParams();
				if (sstruct.disperseOnLoss) popns[i]->allEmigrate();
				else popns[i]->extirpate();
			}
			else { // non-stage-structured species is destroyed
				popns[i]->extirpate();
			}
		}
	}
}

void SubCommunity::reproduction(int resol, float epsGlobal, short rasterType, bool patchModel)
{
	if (subCommNum == 0) return; // no reproduction in the matrix
	float localK, envval;
	Cell* pCell;
	envGradParams grad = paramsGrad->getGradient();
	envStochParams env = paramsStoch->getStoch();

	int npops = (int)popns.size();
	// THE FOLLOWING MAY BE MORE EFFICIENT WHILST THERE IS ONLY ONE SPECIES ...
	if (npops < 1) return;

	localK = pPatch->getK();
	if (localK > 0.0) {
		if (patchModel) {
			envval = 1.0; // environmental gradient is currently not applied for patch-based model
		}
		else { // cell-based model
			if (grad.gradient && grad.gradType == 2)
			{ // gradient in fecundity
				Cell* pCell = pPatch->getRandomCell(); // locate the only cell in the patch
				envval = pCell->getEnvVal();
			}
			else envval = 1.0;
		}
		if (env.stoch && !env.inK) { // stochasticity in fecundity
			if (env.local) {
				if (!patchModel) { // only permitted for cell-based model
					pCell = pPatch->getRandomCell();
					if (pCell != 0) envval += pCell->getEps();
				}
			}
			else { // global stochasticity
				envval += epsGlobal;
			}
		}
		for (int i = 0; i < npops; i++) { // all populations
			popns[i]->reproduction(localK, envval, resol);
			popns[i]->fledge();
		}
	}
}

void SubCommunity::emigration(void)
{
	if (subCommNum == 0) return; // no emigration from the matrix
	int npops = static_cast<int>(popns.size());
	if (npops < 1) return;
	float localK = pPatch->getK();
	// NOTE that even if K is zero, it could have been >0 in previous time-step, and there
	// might be emigrants if there is non-juvenile emigration
	for (int i = 0; i < npops; i++) { // all populations
		popns[i]->emigration(localK);
	}
}

// Remove emigrants from their natal patch and add to patch 0 (matrix)
void SubCommunity::initiateDispersal(SubCommunity* matrix) {
	if (subCommNum == 0) return; // no dispersal initiation in the matrix
	popStats pop;
	disperser disp;

	int npops = (int)popns.size();
	for (int i = 0; i < npops; i++) { // all populations
		pop = popns[i]->getStats();
		for (int j = 0; j < pop.nInds; j++) {
			disp = popns[i]->extractDisperser(j);
			if (disp.yes) { // disperser - has already been removed from natal population
				// add to matrix population
				matrix->recruit(disp.pInd, pop.pSpecies);
			}
		}
		// remove pointers to emigrants
		popns[i]->clean();
	}
}

// Add an individual into the local population of its species in the patch
void SubCommunity::recruit(Individual* pInd, Species* pSpecies) {
	int npops = (int)popns.size();
	for (int i = 0; i < npops; i++) { // all populations
		if (pSpecies == popns[i]->getSpecies()) {
			popns[i]->recruit(pInd);
		}
	}
}

// Transfer through the matrix - run for the matrix sub-community only
#if RS_RCPP
int SubCommunity::transfer(Landscape* pLandscape, short landIx, short nextseason)
#else
int SubCommunity::transfer(Landscape* pLandscape, short landIx)
#endif // RS_RCPP
{
	int ndispersers = 0;
	int npops = (int)popns.size();
	for (int i = 0; i < npops; i++) { // all populations
#if RS_RCPP
		ndispersers += popns[i]->transfer(pLandscape, landIx, nextseason);
#else
		ndispersers += popns[i]->transfer(pLandscape, landIx);
#endif // RS_RCPP

	}
	return ndispersers;
}

//---------------------------------------------------------------------------

// Remove emigrants from patch 0 (matrix) and transfer to sub-community
// in which their destination co-ordinates fall
// This function is executed for the matrix patch only

void SubCommunity::completeDispersal(Landscape* pLandscape, bool connect)
{
	int popsize;
	disperser settler;
	Species* pSpecies;
	Population* pPop;
	Patch* pPrevPatch;
	Patch* pNewPatch;
	Cell* pPrevCell;
	SubCommunity* pSubComm;

	int npops = (int)popns.size();
	for (int i = 0; i < npops; i++) { // all populations
		pSpecies = popns[i]->getSpecies();
		popsize = popns[i]->getNInds();
		for (int j = 0; j < popsize; j++) {
			bool settled;
			settler = popns[i]->extractSettler(j);
			settled = settler.yes;
			if (settled) {
			// settler - has already been removed from matrix population
			// find new patch
				pNewPatch = (Patch*)settler.pCell->getPatch();
				// find population within the patch (if there is one)
				pPop = (Population*)pNewPatch->getPopn((intptr)pSpecies);
				if (pPop == 0) { // settler is the first in a previously uninhabited patch
					// create a new population in the corresponding sub-community
					pSubComm = (SubCommunity*)pNewPatch->getSubComm();
					pPop = pSubComm->newPopn(pLandscape, pSpecies, pNewPatch, 0);
				}
				pPop->recruit(settler.pInd);
				if (connect) { // increment connectivity totals
					int newpatch = pNewPatch->getSeqNum();
					pPrevCell = settler.pInd->getLocn(0); // previous cell
					intptr patch = pPrevCell->getPatch();
					if (patch != 0) {
						pPrevPatch = (Patch*)patch;
						int prevpatch = pPrevPatch->getSeqNum();
						pLandscape->incrConnectMatrix(prevpatch, newpatch);
					}
				}
			}
			else { // for group dispersal only
			}
		}
		// remove pointers in the matrix popn to settlers
		popns[i]->clean();
	}

}

//---------------------------------------------------------------------------

void SubCommunity::survival(short part, short option0, short option1)
{
	int npops = (int)popns.size();
	if (npops < 1) return;
	if (part == 0) {
		float localK = pPatch->getK();
		for (int i = 0; i < npops; i++) { // all populations
			popns[i]->survival0(localK, option0, option1);
		}
	}
	else {
		for (int i = 0; i < npops; i++) { // all populations
			popns[i]->survival1();
		}
	}
}

void SubCommunity::ageIncrement(void) {
	int npops = (int)popns.size();
	for (int i = 0; i < npops; i++) { // all populations
		popns[i]->ageIncrement();
	}
}

// Find the population of a given species in a given patch
Population* SubCommunity::findPop(Species* pSp, Patch* pPch) {

	Population* pPop = 0;
	popStats pop;
	int npops = (int)popns.size();

	for (int i = 0; i < npops; i++) { // all populations
		pop = popns[i]->getStats();
		if (pop.pSpecies == pSp && pop.pPatch == pPch) { // population located
			pPop = popns[i];
			break;
		}
		else pPop = 0;
	}
	return pPop;
}

//---------------------------------------------------------------------------

void SubCommunity::createOccupancy(int nrows) {
	if (occupancy != 0) deleteOccupancy();
	if (nrows > 0) {
		occupancy = new int[nrows];
		for (int i = 0; i < nrows; i++) occupancy[i] = 0;
	}
}

void SubCommunity::updateOccupancy(int row) {
	popStats pop;
	int npops = (int)popns.size();
	for (int i = 0; i < npops; i++) {
		pop = popns[i]->getStats();
		if (pop.nInds > 0 && pop.breeding) {
			occupancy[row]++;
			i = npops;
		}
	}
}

int SubCommunity::getOccupancy(int row) {
	if (row >= 0) return occupancy[row];
	else return 0;
}

void SubCommunity::deleteOccupancy(void) {
	delete[] occupancy;
	occupancy = 0;
}

//---------------------------------------------------------------------------
// Open population file and write header record
bool SubCommunity::outPopHeaders(Landscape* pLandscape, Species* pSpecies, int option)
{
	bool fileOK;
	Population* pPop;
	landParams land = pLandscape->getLandParams();

	if (option == -999) { // close the file
		// as all populations may have been deleted, set up a dummy one
		// species is not necessary
		pPop = new Population();
		fileOK = pPop->outPopHeaders(-999, land.patchModel);
		delete pPop;
	}
	else { // open the file
		// as no population has yet been created, set up a dummy one
		// species is necessary, as columns depend on stage and sex structure
		pPop = new Population(pSpecies, pPatch, 0, land.resol);
		fileOK = pPop->outPopHeaders(land.landNum, land.patchModel);
		delete pPop;
	}
	return fileOK;
}

// Write records to population file
void SubCommunity::outPop(Landscape* pLandscape, int rep, int yr, int gen)
{
	landParams land = pLandscape->getLandParams();
	envGradParams grad = paramsGrad->getGradient();
	envStochParams env = paramsStoch->getStoch();
	bool writeEnv = false;
	bool gradK = false;
	if (grad.gradient) {
		writeEnv = true;
		if (grad.gradType == 1) gradK = true; // ... carrying capacity
	}
	if (env.stoch) writeEnv = true;

	// generate output for each population within the sub-community (patch)
	// provided that the patch is suitable (i.e. non-zero carrying capacity)
	// or the population is above zero (possible if there is stochasticity or a moving gradient)
	// or it is the matrix patch in a patch-based model
	int npops = (int)popns.size();
	int patchnum;
	Cell* pCell;
	float localK;
	float eps = 0.0;
	if (env.stoch) {
		if (env.local) {
			pCell = pPatch->getRandomCell();
			if (pCell != 0) eps = pCell->getEps();
		}
		else {
			eps = pLandscape->getGlobalStoch(yr);
		}
	}

	patchnum = pPatch->getPatchNum();
	for (int i = 0; i < npops; i++) { // all populations
		localK = pPatch->getK();
		if (localK > 0.0 || (land.patchModel && patchnum == 0)) {
			popns[i]->outPopulation(rep, yr, gen, eps, land.patchModel, writeEnv, gradK);
		}
		else {
			if (popns[i]->totalPop() > 0) {
				popns[i]->outPopulation(rep, yr, gen, eps, land.patchModel, writeEnv, gradK);
			}
		}
	}
}

// Write records to individuals file
void SubCommunity::outInds(Landscape* pLandscape, int rep, int yr, int gen, int landNr) {
	landParams ppLand = pLandscape->getLandParams();
	Population* pPop;
	if (landNr >= 0) { // open the file
		popns[0]->outIndsHeaders(rep, landNr, ppLand.patchModel);
		return;
	}
	if (landNr == -999) { // close the file

		// as all populations may have been deleted, set up a dummy one
	// species is not necessary
		pPop = new Population();
		pPop->outIndsHeaders(rep, -999, ppLand.patchModel);
		delete pPop;
		return;
	}
	// generate output for each population within the sub-community (patch)
	int npops = (int)popns.size();
	for (int i = 0; i < npops; i++) { // all populations
		popns[i]->outIndividual(pLandscape, rep, yr, gen, pPatch->getPatchNum());
	}
}


// Population size of a specified stage
int SubCommunity::stagePop(int stage) {
	int popsize = 0;
	int npops = (int)popns.size();
	for (int i = 0; i < npops; i++) { // all populations
		popsize += popns[i]->stagePop(stage);
	}
	return popsize;
}

// Open traits file and write header record
bool SubCommunity::outTraitsHeaders(Landscape* pLandscape, Species* pSpecies, int landNr)
{
	landParams land = pLandscape->getLandParams();
	if (landNr == -999) { // close file
		if (outtraits.is_open()) outtraits.close();
		outtraits.clear();
		return true;
	}

	string name;
	emigRules emig = pSpecies->getEmigRules();
	transferRules trfr = pSpecies->getTransferRules();
	settleType sett = pSpecies->getSettle();
	simParams sim = paramsSim->getSim();

	string DirOut = paramsSim->getDir(2);
	if (sim.batchMode) {
		if (land.patchModel) {
			name = DirOut
				+ "Batch" + to_string(sim.batchNum) + "_"
				+ "Sim" + to_string(sim.simulation) + "_Land" + to_string(landNr)
				+ "_TraitsXpatch.txt";
		}
		else {
			name = DirOut
				+ "Batch" + to_string(sim.batchNum) + "_"
				+ "Sim" + to_string(sim.simulation) + "_Land" + to_string(landNr)
				+ "_TraitsXcell.txt";
		}
	}
	else {
		if (land.patchModel) {
			name = DirOut + "Sim" + to_string(sim.simulation) + "_TraitsXpatch.txt";
		}
		else {
			name = DirOut + "Sim" + to_string(sim.simulation) + "_TraitsXcell.txt";
		}
	}
	outtraits.open(name.c_str());

	outtraits << "Rep\tYear\tRepSeason";
	if (land.patchModel) outtraits << "\tPatchID";
	else
		outtraits << "\tx\ty";

	if (emig.indVar) {
		if (emig.sexDep) {
			if (emig.densDep) {
				outtraits << "\tF_meanD0\tF_stdD0\tM_meanD0\tM_stdD0";
				outtraits << "\tF_meanAlpha\tF_stdAlpha\tM_meanAlpha\tM_stdAlpha";
				outtraits << "\tF_meanBeta\tF_stdBeta\tM_meanBeta\tM_stdBeta";
			}
			else {
				outtraits << "\tF_meanEP\tF_stdEP\tM_meanEP\tM_stdEP";
			}
		}
		else {
			if (emig.densDep) {
				outtraits << "\tmeanD0\tstdD0\tmeanAlpha\tstdAlpha";
				outtraits << "\tmeanBeta\tstdBeta";
			}
			else {
				outtraits << "\tmeanEP\tstdEP";
			}
		}
	}
	if (trfr.indVar) {
		if (trfr.usesMovtProc) {
			if (trfr.moveType == 1) {
				outtraits << "\tmeanDP\tstdDP\tmeanGB\tstdGB";
				outtraits << "\tmeanAlphaDB\tstdAlphaDB\tmeanBetaDB\tstdBetaDB";
			}
			if (trfr.moveType == 2) {
				outtraits << "\tmeanStepLength\tstdStepLength\tmeanRho\tstdRho";
			}
		}
		else {
			if (trfr.sexDep) {
				outtraits << "\tF_mean_distI\tF_std_distI\tM_mean_distI\tM_std_distI";
				if (trfr.twinKern)
					outtraits << "\tF_mean_distII\tF_std_distII\tM_mean_distII\tM_std_distII"
					<< "\tF_meanPfirstKernel\tF_stdPfirstKernel"
					<< "\tM_meanPfirstKernel\tM_stdPfirstKernel";
			}
			else {
				outtraits << "\tmean_distI\tstd_distI";
				if (trfr.twinKern)
					outtraits << "\tmean_distII\tstd_distII\tmeanPfirstKernel\tstdPfirstKernel";
			}
		}
	}
	if (sett.indVar) {
		if (sett.sexDep) {
			outtraits << "\tF_meanS0\tF_stdS0\tM_meanS0\tM_stdS0";
			outtraits << "\tF_meanAlphaS\tF_stdAlphaS\tM_meanAlphaS\tM_stdAlphaS";
			outtraits << "\tF_meanBetaS\tF_stdBetaS\tM_meanBetaS\tM_stdBetaS";
		}
		else {
			outtraits << "\tmeanS0\tstdS0";
			outtraits << "\tmeanAlphaS\tstdAlphaS";
			outtraits << "\tmeanBetaS\tstdBetaS";
		}
	}
	if (pSpecies->getNbGenLoadTraits() > 0) {
		if (pSpecies->getDemogrParams().repType > 0) {
			outtraits << "\tF_meanGenFitness\tF_stdGenFitness\tM_meanGenFitness\tM_stdGenFitness";
		}
		else {
			outtraits << "\tmeanGenFitness\tstdGenFitness";
		}
	}

	outtraits << endl;

	return outtraits.is_open();
}

// Write records to traits file and return aggregated sums
traitsums SubCommunity::outTraits(Landscape* pLandscape, int rep, int yr, int gen, bool commlevel)
{
	int popsize, ploidy;
	landParams land = pLandscape->getLandParams();
	simParams sim = paramsSim->getSim();
	bool writefile = false;
	if (sim.outTraitsCells && yr % sim.outIntTraitCell == 0 && !commlevel)
		writefile = true;
	traitsums ts, indTraitsSums;
	for (int i = 0; i < gMaxNbSexes; i++) {
		ts.ninds[i] = 0;
		ts.sumD0[i] = ts.ssqD0[i] = 0.0;
		ts.sumAlpha[i] = ts.ssqAlpha[i] = 0.0; ts.sumBeta[i] = ts.ssqBeta[i] = 0.0;
		ts.sumDist1[i] = ts.ssqDist1[i] = 0.0; ts.sumDist2[i] = ts.ssqDist2[i] = 0.0;
		ts.sumProp1[i] = ts.ssqProp1[i] = 0.0;
		ts.sumDP[i] = ts.ssqDP[i] = 0.0;
		ts.sumGB[i] = ts.ssqGB[i] = 0.0;
		ts.sumAlphaDB[i] = ts.ssqAlphaDB[i] = 0.0;
		ts.sumBetaDB[i] = ts.ssqBetaDB[i] = 0.0;
		ts.sumStepL[i] = ts.ssqStepL[i] = 0.0; ts.sumRho[i] = ts.ssqRho[i] = 0.0;
		ts.sumS0[i] = ts.ssqS0[i] = 0.0;
		ts.sumAlphaS[i] = ts.ssqAlphaS[i] = 0.0; ts.sumBetaS[i] = ts.ssqBetaS[i] = 0.0;
		ts.sumGeneticFitness[i] = ts.ssqGeneticFitness[i] = 0.0;
	}

	// generate output for each population within the sub-community (patch)
	// provided that the patch is suitable (i.e. non-zero carrying capacity)
	int npops = (int)popns.size();
	Species* pSpecies;
	float localK;

	for (int iPop = 0; iPop < npops; iPop++) { // all populations
		localK = pPatch->getK();
		if (localK > 0.0 && popns[iPop]->getNInds() > 0) {
			pSpecies = popns[iPop]->getSpecies();
			demogrParams dem = pSpecies->getDemogrParams();
			emigRules emig = pSpecies->getEmigRules();
			transferRules trfr = pSpecies->getTransferRules();
			settleType sett = pSpecies->getSettle();
			indTraitsSums = popns[iPop]->getIndTraitsSums(pSpecies);

			if (writefile) {
				outtraits << rep << "\t" << yr << "\t" << gen;
				if (land.patchModel) {
					outtraits << "\t" << pPatch->getPatchNum();
				}
				else {
					locn loc = pPatch->getCellLocn(0);
					outtraits << "\t" << loc.x << "\t" << loc.y;
				}
			}

			if (emig.indVar) {
				if (emig.sexDep) { // must be a sexual species
					ploidy = 2;
				}
				else {
					if (dem.repType == 0) { // asexual reproduction
						ploidy = 1;
					}
					else { // sexual reproduction
						ploidy = 1;
					}
				}
				double mnD0[2], mnAlpha[2], mnBeta[2], sdD0[2], sdAlpha[2], sdBeta[2];
				for (int whichChromosome = 0; whichChromosome < ploidy; whichChromosome++) {
					mnD0[whichChromosome] = mnAlpha[whichChromosome] = mnBeta[whichChromosome] = sdD0[whichChromosome] = sdAlpha[whichChromosome] = sdBeta[whichChromosome] = 0.0;
					// individuals may have been counted by sex if there was
					// sex dependency in another dispersal phase
					if (ploidy == 2) popsize = indTraitsSums.ninds[whichChromosome];
					else popsize = indTraitsSums.ninds[0] + indTraitsSums.ninds[1];
					if (popsize > 0) {
						mnD0[whichChromosome] = indTraitsSums.sumD0[whichChromosome] / (double)popsize;
						mnAlpha[whichChromosome] = indTraitsSums.sumAlpha[whichChromosome] / (double)popsize;
						mnBeta[whichChromosome] = indTraitsSums.sumBeta[whichChromosome] / (double)popsize;
						if (popsize > 1) {
							sdD0[whichChromosome] = indTraitsSums.ssqD0[whichChromosome] / (double)popsize - mnD0[whichChromosome] * mnD0[whichChromosome];
							if (sdD0[whichChromosome] > 0.0) sdD0[whichChromosome] = sqrt(sdD0[whichChromosome]); else sdD0[whichChromosome] = 0.0;
							sdAlpha[whichChromosome] = indTraitsSums.ssqAlpha[whichChromosome] / (double)popsize - mnAlpha[whichChromosome] * mnAlpha[whichChromosome];
							if (sdAlpha[whichChromosome] > 0.0) sdAlpha[whichChromosome] = sqrt(sdAlpha[whichChromosome]); else sdAlpha[whichChromosome] = 0.0;
							sdBeta[whichChromosome] = indTraitsSums.ssqBeta[whichChromosome] / (double)popsize - mnBeta[whichChromosome] * mnBeta[whichChromosome];
							if (sdBeta[whichChromosome] > 0.0) sdBeta[whichChromosome] = sqrt(sdBeta[whichChromosome]); else sdBeta[whichChromosome] = 0.0;
						}
						else {
							sdD0[whichChromosome] = sdAlpha[whichChromosome] = sdBeta[whichChromosome] = 0.0;
						}
					}
				}
				if (writefile) {
					if (emig.sexDep) {
						outtraits << "\t" << mnD0[0] << "\t" << sdD0[0];
						outtraits << "\t" << mnD0[1] << "\t" << sdD0[1];
						if (emig.densDep) {
							outtraits << "\t" << mnAlpha[0] << "\t" << sdAlpha[0];
							outtraits << "\t" << mnAlpha[1] << "\t" << sdAlpha[1];
							outtraits << "\t" << mnBeta[0] << "\t" << sdBeta[0];
							outtraits << "\t" << mnBeta[1] << "\t" << sdBeta[1];
						}
					}
					else { // sex-independent
						outtraits << "\t" << mnD0[0] << "\t" << sdD0[0];
						if (emig.densDep) {
							outtraits << "\t" << mnAlpha[0] << "\t" << sdAlpha[0];
							outtraits << "\t" << mnBeta[0] << "\t" << sdBeta[0];
						}
					}
				}
			}

			if (trfr.indVar) {
				if (trfr.usesMovtProc) {
					// CURRENTLY INDIVIDUAL VARIATION CANNOT BE SEX-DEPENDENT
					ploidy = 1;
				}
				else {
					if (trfr.sexDep) { // must be a sexual species
						ploidy = 2;
					}
					else {
						ploidy = 1;
					}
				}
				double mnDist1[2], mnDist2[2], mnProp1[2], mnStepL[2], mnRho[2];
				double sdDist1[2], sdDist2[2], sdProp1[2], sdStepL[2], sdRho[2];
				double mnDP[2], mnGB[2], mnAlphaDB[2], mnBetaDB[2];
				double sdDP[2], sdGB[2], sdAlphaDB[2], sdBetaDB[2];
				for (int whichChromosome = 0; whichChromosome < ploidy; whichChromosome++) {
					mnDist1[whichChromosome] = mnDist2[whichChromosome] = mnProp1[whichChromosome] = mnStepL[whichChromosome] = mnRho[whichChromosome] = 0.0;
					sdDist1[whichChromosome] = sdDist2[whichChromosome] = sdProp1[whichChromosome] = sdStepL[whichChromosome] = sdRho[whichChromosome] = 0.0;
					mnDP[whichChromosome] = mnGB[whichChromosome] = mnAlphaDB[whichChromosome] = mnBetaDB[whichChromosome] = 0.0;
					sdDP[whichChromosome] = sdGB[whichChromosome] = sdAlphaDB[whichChromosome] = sdBetaDB[whichChromosome] = 0.0;
					// individuals may have been counted by sex if there was
					// sex dependency in another dispersal phase
					if (ploidy == 2) popsize = indTraitsSums.ninds[whichChromosome];
					else popsize = indTraitsSums.ninds[0] + indTraitsSums.ninds[1];
					if (popsize > 0) {
						mnDist1[whichChromosome] = indTraitsSums.sumDist1[whichChromosome] / (double)popsize;
						mnDist2[whichChromosome] = indTraitsSums.sumDist2[whichChromosome] / (double)popsize;
						mnProp1[whichChromosome] = indTraitsSums.sumProp1[whichChromosome] / (double)popsize;
						mnStepL[whichChromosome] = indTraitsSums.sumStepL[whichChromosome] / (double)popsize;
						mnRho[whichChromosome] = indTraitsSums.sumRho[whichChromosome] / (double)popsize;
						mnDP[whichChromosome] = indTraitsSums.sumDP[whichChromosome] / (double)popsize;
						mnGB[whichChromosome] = indTraitsSums.sumGB[whichChromosome] / (double)popsize;
						mnAlphaDB[whichChromosome] = indTraitsSums.sumAlphaDB[whichChromosome] / (double)popsize;
						mnBetaDB[whichChromosome] = indTraitsSums.sumBetaDB[whichChromosome] / (double)popsize;
						if (popsize > 1) {
							sdDist1[whichChromosome] = indTraitsSums.ssqDist1[whichChromosome] / (double)popsize - mnDist1[whichChromosome] * mnDist1[whichChromosome];
							if (sdDist1[whichChromosome] > 0.0) sdDist1[whichChromosome] = sqrt(sdDist1[whichChromosome]); else sdDist1[whichChromosome] = 0.0;
							sdDist2[whichChromosome] = indTraitsSums.ssqDist2[whichChromosome] / (double)popsize - mnDist2[whichChromosome] * mnDist2[whichChromosome];
							if (sdDist2[whichChromosome] > 0.0) sdDist2[whichChromosome] = sqrt(sdDist2[whichChromosome]); else sdDist2[whichChromosome] = 0.0;
							sdProp1[whichChromosome] = indTraitsSums.ssqProp1[whichChromosome] / (double)popsize - mnProp1[whichChromosome] * mnProp1[whichChromosome];
							if (sdProp1[whichChromosome] > 0.0) sdProp1[whichChromosome] = sqrt(sdProp1[whichChromosome]); else sdProp1[whichChromosome] = 0.0;
							sdStepL[whichChromosome] = indTraitsSums.ssqStepL[whichChromosome] / (double)popsize - mnStepL[whichChromosome] * mnStepL[whichChromosome];
							if (sdStepL[whichChromosome] > 0.0) sdStepL[whichChromosome] = sqrt(sdStepL[whichChromosome]); else sdStepL[whichChromosome] = 0.0;
							sdRho[whichChromosome] = indTraitsSums.ssqRho[whichChromosome] / (double)popsize - mnRho[whichChromosome] * mnRho[whichChromosome];
							if (sdRho[whichChromosome] > 0.0) sdRho[whichChromosome] = sqrt(sdRho[whichChromosome]); else sdRho[whichChromosome] = 0.0;
							sdDP[whichChromosome] = indTraitsSums.ssqDP[whichChromosome] / (double)popsize - mnDP[whichChromosome] * mnDP[whichChromosome];
							if (sdDP[whichChromosome] > 0.0) sdDP[whichChromosome] = sqrt(sdDP[whichChromosome]); else sdDP[whichChromosome] = 0.0;
							sdGB[whichChromosome] = indTraitsSums.ssqGB[whichChromosome] / (double)popsize - mnGB[whichChromosome] * mnGB[whichChromosome];
							if (sdGB[whichChromosome] > 0.0) sdGB[whichChromosome] = sqrt(sdGB[whichChromosome]); else sdGB[whichChromosome] = 0.0;
							sdAlphaDB[whichChromosome] = indTraitsSums.ssqAlphaDB[whichChromosome] / (double)popsize - mnAlphaDB[whichChromosome] * mnAlphaDB[whichChromosome];
							if (sdAlphaDB[whichChromosome] > 0.0) sdAlphaDB[whichChromosome] = sqrt(sdAlphaDB[whichChromosome]); else sdAlphaDB[whichChromosome] = 0.0;
							sdBetaDB[whichChromosome] = indTraitsSums.ssqBetaDB[whichChromosome] / (double)popsize - mnBetaDB[whichChromosome] * mnBetaDB[whichChromosome];
							if (sdBetaDB[whichChromosome] > 0.0) sdBetaDB[whichChromosome] = sqrt(sdBetaDB[whichChromosome]); else sdBetaDB[whichChromosome] = 0.0;
						}
					}
				}
				if (writefile) {
					if (trfr.usesMovtProc) {
						if (trfr.moveType == 1) {
							outtraits << "\t" << mnDP[0] << "\t" << sdDP[0];
							outtraits << "\t" << mnGB[0] << "\t" << sdGB[0];
							outtraits << "\t" << mnAlphaDB[0] << "\t" << sdAlphaDB[0];
							outtraits << "\t" << mnBetaDB[0] << "\t" << sdBetaDB[0];
						}
						if (trfr.moveType == 2) {
							outtraits << "\t" << mnStepL[0] << "\t" << sdStepL[0];
							outtraits << "\t" << mnRho[0] << "\t" << sdRho[0];
						}
					}
					else {
						if (trfr.sexDep) {
							outtraits << "\t" << mnDist1[0] << "\t" << sdDist1[0];
							outtraits << "\t" << mnDist1[1] << "\t" << sdDist1[1];
							if (trfr.twinKern)
							{
								outtraits << "\t" << mnDist2[0] << "\t" << sdDist2[0];
								outtraits << "\t" << mnDist2[1] << "\t" << sdDist2[1];
								outtraits << "\t" << mnProp1[0] << "\t" << sdProp1[0];
								outtraits << "\t" << mnProp1[1] << "\t" << sdProp1[1];
							}
						}
						else { // sex-independent
							outtraits << "\t" << mnDist1[0] << "\t" << sdDist1[0];
							if (trfr.twinKern)
							{
								outtraits << "\t" << mnDist2[0] << "\t" << sdDist2[0];
								outtraits << "\t" << mnProp1[0] << "\t" << sdProp1[0];
							}
						}
					}
				}
			}

			if (sett.indVar) {
				if (sett.sexDep) { // must be a sexual species
					ploidy = 2;
				}
				else {
					if (dem.repType == 0) { // asexual reproduction
						ploidy = 1;
					}
					else { // sexual reproduction
						ploidy = 1;
					}
				}
				// CURRENTLY INDIVIDUAL VARIATION CANNOT BE SEX-DEPENDENT
				double mnS0[2], mnAlpha[2], mnBeta[2], sdS0[2], sdAlpha[2], sdBeta[2];
				for (int whichChromosome = 0; whichChromosome < ploidy; whichChromosome++) {
					mnS0[whichChromosome] = mnAlpha[whichChromosome] = mnBeta[whichChromosome] = sdS0[whichChromosome] = sdAlpha[whichChromosome] = sdBeta[whichChromosome] = 0.0;
					// individuals may have been counted by sex if there was
					// sex dependency in another dispersal phase
					if (ploidy == 2) popsize = indTraitsSums.ninds[whichChromosome];
					else popsize = indTraitsSums.ninds[0] + indTraitsSums.ninds[1];
					
					if (popsize > 0) {
						
						mnS0[whichChromosome] = indTraitsSums.sumS0[whichChromosome] / (double)popsize;
						mnAlpha[whichChromosome] = indTraitsSums.sumAlphaS[whichChromosome] / (double)popsize;
						mnBeta[whichChromosome] = indTraitsSums.sumBetaS[whichChromosome] / (double)popsize;
						
						if (popsize > 1) {
							sdS0[whichChromosome] = indTraitsSums.ssqS0[whichChromosome] / (double)popsize - mnS0[whichChromosome] * mnS0[whichChromosome];
							if (sdS0[whichChromosome] > 0.0) sdS0[whichChromosome] = sqrt(sdS0[whichChromosome]); else sdS0[whichChromosome] = 0.0;
							sdAlpha[whichChromosome] = indTraitsSums.ssqAlphaS[whichChromosome] / (double)popsize - mnAlpha[whichChromosome] * mnAlpha[whichChromosome];
							if (sdAlpha[whichChromosome] > 0.0) sdAlpha[whichChromosome] = sqrt(sdAlpha[whichChromosome]); else sdAlpha[whichChromosome] = 0.0;
							sdBeta[whichChromosome] = indTraitsSums.ssqBetaS[whichChromosome] / (double)popsize - mnBeta[whichChromosome] * mnBeta[whichChromosome];
							if (sdBeta[whichChromosome] > 0.0) sdBeta[whichChromosome] = sqrt(sdBeta[whichChromosome]); else sdBeta[whichChromosome] = 0.0;
						}
						else {
							sdS0[whichChromosome] = sdAlpha[whichChromosome] = sdBeta[whichChromosome] = 0.0;
						}
					}
				}
				if (writefile) {
					if (sett.sexDep) {
						outtraits << "\t" << mnS0[0] << "\t" << sdS0[0];
						outtraits << "\t" << mnS0[1] << "\t" << sdS0[1];
						outtraits << "\t" << mnAlpha[0] << "\t" << sdAlpha[0];
						outtraits << "\t" << mnAlpha[1] << "\t" << sdAlpha[1];
						outtraits << "\t" << mnBeta[0] << "\t" << sdBeta[0];
						outtraits << "\t" << mnBeta[1] << "\t" << sdBeta[1];
					}
					else { // sex-independent
						outtraits << "\t" << mnS0[0] << "\t" << sdS0[0];
						outtraits << "\t" << mnAlpha[0] << "\t" << sdAlpha[0];
						outtraits << "\t" << mnBeta[0] << "\t" << sdBeta[0];
					}
				}
			}

			// Genetic load
			if (pSpecies->getNbGenLoadTraits() > 0) {

				ploidy = pSpecies->isDiploid() + 1;
				double mnGenFitness[2], sdGenFitness[2];

				for (int whichChromosome = 0; whichChromosome < ploidy; whichChromosome++) {
					mnGenFitness[whichChromosome] = sdGenFitness[whichChromosome] = 0.0;

					if (ploidy == 2) popsize = indTraitsSums.ninds[whichChromosome];
					else popsize = indTraitsSums.ninds[0] + indTraitsSums.ninds[1];

					if (popsize > 0) {

						mnGenFitness[whichChromosome] = indTraitsSums.sumGeneticFitness[whichChromosome] / (double)popsize;
						if (popsize > 1) {
							sdGenFitness[whichChromosome] = indTraitsSums.ssqGeneticFitness[whichChromosome] / (double)popsize - mnGenFitness[whichChromosome] * mnGenFitness[whichChromosome];
							if (sdGenFitness[whichChromosome] > 0.0) sdGenFitness[whichChromosome] = sqrt(sdGenFitness[whichChromosome]); else sdGenFitness[whichChromosome] = 0.0;
						}
						else {
							sdGenFitness[whichChromosome] = 0.0;
						}
					}
				}

				if (writefile) {
					if (pSpecies->getDemogrParams().repType > 0) {
						outtraits << "\t" << mnGenFitness[0] << "\t" << sdGenFitness[0];
						outtraits << "\t" << mnGenFitness[1] << "\t" << sdGenFitness[1];
					}
					else { // sex-independent
						outtraits << "\t" << mnGenFitness[0] << "\t" << sdGenFitness[0];
					}
				}
			}

			// CURRENTLY INDIVIDUAL VARIATION CANNOT BE SEX-DEPENDENT
			if (writefile) outtraits << endl;

			for (int iSex = 0; iSex < gMaxNbSexes; iSex++) {
				ts.ninds[iSex] += indTraitsSums.ninds[iSex];
				ts.sumD0[iSex] += indTraitsSums.sumD0[iSex];     
				ts.ssqD0[iSex] += indTraitsSums.ssqD0[iSex];
				ts.sumAlpha[iSex] += indTraitsSums.sumAlpha[iSex];  
				ts.ssqAlpha[iSex] += indTraitsSums.ssqAlpha[iSex];
				ts.sumBeta[iSex] += indTraitsSums.sumBeta[iSex];   
				ts.ssqBeta[iSex] += indTraitsSums.ssqBeta[iSex];
				ts.sumDist1[iSex] += indTraitsSums.sumDist1[iSex];  
				ts.ssqDist1[iSex] += indTraitsSums.ssqDist1[iSex];
				ts.sumDist2[iSex] += indTraitsSums.sumDist2[iSex];  
				ts.ssqDist2[iSex] += indTraitsSums.ssqDist2[iSex];
				ts.sumProp1[iSex] += indTraitsSums.sumProp1[iSex];  
				ts.ssqProp1[iSex] += indTraitsSums.ssqProp1[iSex];
				ts.sumDP[iSex] += indTraitsSums.sumDP[iSex];    
				ts.ssqDP[iSex] += indTraitsSums.ssqDP[iSex];
				ts.sumGB[iSex] += indTraitsSums.sumGB[iSex];    
				ts.ssqGB[iSex] += indTraitsSums.ssqGB[iSex];
				ts.sumAlphaDB[iSex] += indTraitsSums.sumAlphaDB[iSex];
				ts.ssqAlphaDB[iSex] += indTraitsSums.ssqAlphaDB[iSex];
				ts.sumBetaDB[iSex] += indTraitsSums.sumBetaDB[iSex]; 
				ts.ssqBetaDB[iSex] += indTraitsSums.ssqBetaDB[iSex];
				ts.sumStepL[iSex] += indTraitsSums.sumStepL[iSex]; 
				ts.ssqStepL[iSex] += indTraitsSums.ssqStepL[iSex];
				ts.sumRho[iSex] += indTraitsSums.sumRho[iSex];
				ts.ssqRho[iSex] += indTraitsSums.ssqRho[iSex];
				ts.sumS0[iSex] += indTraitsSums.sumS0[iSex];
				ts.ssqS0[iSex] += indTraitsSums.ssqS0[iSex];
				ts.sumAlphaS[iSex] += indTraitsSums.sumAlphaS[iSex]; 
				ts.ssqAlphaS[iSex] += indTraitsSums.ssqAlphaS[iSex];
				ts.sumBetaS[iSex] += indTraitsSums.sumBetaS[iSex]; 
				ts.ssqBetaS[iSex] += indTraitsSums.ssqBetaS[iSex];
				ts.sumGeneticFitness[iSex] += indTraitsSums.sumGeneticFitness[iSex];  
				ts.ssqGeneticFitness[iSex] += indTraitsSums.ssqGeneticFitness[iSex];
			}
		}
	}
	return ts;
}

//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


