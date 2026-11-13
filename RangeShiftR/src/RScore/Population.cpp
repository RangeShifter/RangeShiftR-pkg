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
 --------------------------------------------------------------------------*/


 //---------------------------------------------------------------------------

#include "Population.h"

// Draw the number of offspring of a breeding female.
// - fecSD < 0  : Poisson distribution with mean meanFecundity (default)
// - fecSD >= 0 : Normal distribution with mean meanFecundity and SD fecSD;
//                the draw is truncated towards zero (deliberately, NOT rounded)
//                and negative values are set to 0. Note that truncation lowers
//                the realised mean fecundity relative to meanFecundity.
// - maxOffspring >= 0 : upper limit; larger draws are set to maxOffspring
//                       (applies to both distributions; < 0: no limit)
static int drawNbOffspring(double meanFecundity, float fecSD, int maxOffspring) {
	if (meanFecundity <= 0.0) return 0;
	int n;
	if (fecSD < 0.0) n = pRandom->Poisson(meanFecundity);
	else n = std::max(0, static_cast<int>(pRandom->Normal(meanFecundity, fecSD)));
	if (maxOffspring >= 0 && n > maxOffspring) n = maxOffspring;
	return n;
}

#include <algorithm>
//---------------------------------------------------------------------------

// Population and Individuals output files: text (default) or binary, see OutputTable
OutputTable outPop;
OutputTable outInds;

//---------------------------------------------------------------------------
// The inds vector holds pointers to individuals that are scattered over the heap
// (and is shuffled after every clean()). Loops over all individuals are therefore
// dominated by cache misses on the first access to each Individual. Asking the CPU
// to start loading the individual a few iterations ahead hides most of this
// latency. Pure performance hint: no effect on results.
static inline void prefetchInd(const std::vector<Individual*>& v, int i) {
#if defined(__GNUC__) || defined(__clang__)
	constexpr int distance = 8;
	if (i + distance < (int)v.size()) __builtin_prefetch(v[i + distance]);
#endif
}

//---------------------------------------------------------------------------

Population::Population(void) {
	nSexes = nStages = 0;
	pPatch = NULL;
	pSpecies = NULL;
	return;
}

Population::Population(Species* pSp, Patch* pPch, int ninds, int resol)
{
	// constructor for a Population of a specified size

	int n, nindivs, age = 0, minage, maxage, nAges = 0;
	int cumtotal = 0;
	float probmale;
	double ageprob, ageprobsum;
	std::vector <double> ageProb; // for quasi-equilibrium initial age distribution
	Cell* pCell;

	if (ninds > 0) {
		inds.reserve(ninds);
		juvs.reserve(ninds);
	}

	pSpecies = pSp;
	pPatch = pPch;
	// record the new population in the patch
	patchPopn pp;
	pp.pSp = pSpecies; pp.pPop = this;
	pPatch->addPopn(pp);

	demogrParams dem = pSpecies->getDemogrParams();
	stageParams sstruct = pSpecies->getStageParams();
	emigRules emig = pSpecies->getEmigRules();
	transferRules trfr = pSpecies->getTransferRules();
	settleType sett = pSpecies->getSettle();
	initParams init = paramsInit->getInit();

	// determine no. of stages and sexes of species to initialise
	if (dem.stageStruct) {
		nStages = sstruct.nStages;
	}
	else // non-structured population has 2 stages, but user only ever sees stage 1
		nStages = 2;
	if (dem.repType == 0) { nSexes = 1; probmale = 0.0; }
	else { nSexes = 2; probmale = dem.propMales; }

	// set up population sub-totals
	for (int stg = 0; stg < gMaxNbStages; stg++) {
		for (int sex = 0; sex < gMaxNbSexes; sex++) {
			nInds[stg][sex] = 0;
		}
	}

	// set up local copy of minimum age table
	short minAge[gMaxNbStages][gMaxNbSexes];
	for (int stg = 0; stg < nStages; stg++) {
		for (int sex = 0; sex < nSexes; sex++) {
			if (dem.stageStruct) {
				if (dem.repType == 1) { // simple sexual model
					// both sexes use minimum ages recorded for females
					minAge[stg][sex] = pSpecies->getMinAge(stg, 0);
				}
				else {
					minAge[stg][sex] = pSpecies->getMinAge(stg, sex);
				}
			}
			else { // non-structured population
				minAge[stg][sex] = 0;
			}
		}
	}

	// individuals of new population must be >= stage 1
	for (int stg = 1; stg < nStages; stg++) {
		if (dem.stageStruct) { // allocate to stages according to initialisation conditions
			// final stage is treated separately to ensure that correct total
			// no. of individuals is created
			if (stg == nStages - 1) {
				n = ninds - cumtotal;
			}
			else {
				n = (int)(ninds * paramsInit->getProp(stg) + 0.5);
				cumtotal += n;
			}
		}
		else { // non-structured - all individuals go into stage 1
			n = ninds;
		}
		// establish initial age distribution
		minage = maxage = stg;
		if (dem.stageStruct) {
			// allow for stage-dependent minimum ages (use whichever sex is greater)
			if (minAge[stg][0] > 0 && minage < minAge[stg][0]) minage = minAge[stg][0];
			if (nSexes == 2 && minAge[stg][1] > 0 && minage < minAge[stg][1]) minage = minAge[stg][1];
			// allow for specified age distribution
			if (init.initAge != 0) { // not lowest age
				if (stg == nStages - 1) maxage = sstruct.maxAge; // final stage
				else { // all other stages - use female max age, as sex of individuals is not predetermined
					maxage = minAge[stg + 1][0] - 1;
				}
				if (maxage < minage) maxage = minage;
				nAges = maxage - minage + 1;
				if (init.initAge == 2) { // quasi-equilibrium distribution
					double psurv = (double)pSpecies->getSurv(stg, 0); // use female survival for the stage
					ageProb.clear();
					ageprobsum = 0.0;
					ageprob = 1.0;
					for (int i = 0; i < nAges; i++) {
						ageProb.push_back(ageprob); ageprobsum += ageprob; ageprob *= psurv;
					}
					for (int i = 0; i < nAges; i++) {
						ageProb[i] /= ageprobsum;
						if (i > 0) ageProb[i] += ageProb[i - 1]; // to give cumulative probability
					}
				}
			}
		}
		// create individuals
		int sex;
		nindivs = (int)inds.size();
		for (int i = 0; i < n; i++) {
			pCell = pPatch->getRandomCell();
			if (dem.stageStruct) {
				switch (init.initAge) {
				case 0: // lowest possible age
					age = minage;
					break;
				case 1: // randomised
					if (maxage > minage) age = pRandom->IRandom(minage, maxage);
					else age = minage;
					break;
				case 2: // quasi-equilibrium
					if (nAges > 1) {
						double rrr = pRandom->Random();
						int ageclass = 0;
						// guard against rounding: the last cumulative value may be < 1
						while (ageclass < nAges - 1 && rrr > ageProb[ageclass]) ageclass++;
						age = minage + ageclass;
					}
					else age = minage;
					break;
				}
			}
			else age = stg;

			Individual* newInd = new Individual(pSpecies, pCell, pPatch, stg, age, sstruct.repInterval,
				probmale, trfr.usesMovtProc, trfr.moveType);

			if (pSpecies->getNTraits() > 0) {
				newInd->setUpGenes(pSpecies, resol);
			}
			inds.push_back(newInd);
			nInds[stg][newInd->getSex()]++;
		}
	}
}

Population::~Population(void) {
	int ninds = (int)inds.size();
	for (int i = 0; i < ninds; i++) {
		if (inds[i] != NULL) delete inds[i];
	}
	inds.clear();
	int njuvs = (int)juvs.size();
	for (int i = 0; i < njuvs; i++) {
		if (juvs[i] != NULL) delete juvs[i];
	}
	juvs.clear();
	for (Individual* pInd : harvestedInds) delete pInd;
	harvestedInds.clear();
	int nsampledInds = (int)sampledInds.size();
	for (int i = 0; i < nsampledInds; i++) {
	    if (sampledInds[i] != NULL) sampledInds[i]=NULL;
}
	sampledInds.clear();
}

traitsums Population::getIndTraitsSums(Species* pSpecies) {

	traitsums ts = traitsums();
	for (int sex = 0; sex < gMaxNbSexes; sex++) {
		ts.ninds[sex] = 0;
		ts.sumD0[sex] = ts.ssqD0[sex] = 0.0;
		ts.sumAlpha[sex] = ts.ssqAlpha[sex] = 0.0;
		ts.sumBeta[sex] = ts.ssqBeta[sex] = 0.0;
		ts.sumDist1[sex] = ts.ssqDist1[sex] = 0.0;
		ts.sumDist2[sex] = ts.ssqDist2[sex] = 0.0;
		ts.sumProp1[sex] = ts.ssqProp1[sex] = 0.0;
		ts.sumDP[sex] = ts.ssqDP[sex] = 0.0;
		ts.sumGB[sex] = ts.ssqGB[sex] = 0.0;
		ts.sumAlphaDB[sex] = ts.ssqAlphaDB[sex] = 0.0;
		ts.sumBetaDB[sex] = ts.ssqBetaDB[sex] = 0.0;
		ts.sumStepL[sex] = ts.ssqStepL[sex] = 0.0;
		ts.sumRho[sex] = ts.ssqRho[sex] = 0.0;
		ts.sumS0[sex] = ts.ssqS0[sex] = 0.0;
		ts.sumAlphaS[sex] = ts.ssqAlphaS[sex] = 0.0;
		ts.sumBetaS[sex] = ts.ssqBetaS[sex] = 0.0;
		ts.sumGeneticFitness[sex] = ts.ssqGeneticFitness[sex] = 0.0;
	}

	emigRules emig = pSpecies->getEmigRules();
	transferRules trfr = pSpecies->getTransferRules();
	settleType sett = pSpecies->getSettle();

	for (auto& ind : inds) {

		int sex = ind->getSex();
		ts.ninds[sex] += 1;

		// emigration traits
		emigTraits e = ind->getIndEmigTraits();
		ts.sumD0[sex] += e.d0;    
		ts.ssqD0[sex] += e.d0 * e.d0;
		ts.sumAlpha[sex] += e.alpha; 
		ts.ssqAlpha[sex] += e.alpha * e.alpha;
		ts.sumBeta[sex] += e.beta;  
		ts.ssqBeta[sex] += e.beta * e.beta;

		// transfer traits
		if (trfr.usesMovtProc) {

			switch (trfr.moveType) {

			case 1: { // SMS
				trfrSMSTraits sms = ind->getIndSMSTraits();
				ts.sumDP[sex] += sms.dp; 
				ts.ssqDP[sex] += sms.dp * sms.dp;
				ts.sumGB[sex] += sms.gb;
				ts.ssqGB[sex] += sms.gb * sms.gb;
				ts.sumAlphaDB[sex] += sms.alphaDB;
				ts.ssqAlphaDB[sex] += sms.alphaDB * sms.alphaDB;
				ts.sumBetaDB[sex] += sms.betaDB; 
				ts.ssqBetaDB[sex] += sms.betaDB * sms.betaDB;
				break;
			}
			case 2: {
				trfrCRWTraits c = ind->getIndCRWTraits();
				ts.sumStepL[sex] += c.stepLength;
				ts.ssqStepL[sex] += c.stepLength * c.stepLength;
				ts.sumRho[sex] += c.rho;       
				ts.ssqRho[sex] += c.rho * c.rho;
				break;
			}
			default:
				throw runtime_error("usesMoveProcess is ON but moveType is neither 1 (SMS) or 2 (CRW).");
				break;
			}
		}
		else {
			trfrKernelParams k = ind->getIndKernTraits();
			ts.sumDist1[sex] += k.meanDist1; 
			ts.ssqDist1[sex] += k.meanDist1 * k.meanDist1;
			ts.sumDist2[sex] += k.meanDist2;
			ts.ssqDist2[sex] += k.meanDist2 * k.meanDist2;
			ts.sumProp1[sex] += k.probKern1; 
			ts.ssqProp1[sex] += k.probKern1 * k.probKern1;
		}
		// settlement traits
		settleTraits s = ind->getIndSettTraits();
		ts.sumS0[sex] += s.s0;     
		ts.ssqS0[sex] += s.s0 * s.s0;
		ts.sumAlphaS[sex] += s.alpha; 
		ts.ssqAlphaS[sex] += s.alpha * s.alpha;
		ts.sumBetaS[sex] += s.beta;   
		ts.ssqBetaS[sex] += s.beta * s.beta;

		double fitness = ind->getGeneticFitness();
		ts.sumGeneticFitness[sex] += fitness;
		ts.ssqGeneticFitness[sex] += fitness * fitness;
	}
	return ts;
}
//int Population::getNInds() { return static_cast<int>(inds.size()); }

// ----------------------------------------------------------------------------------------
// reset allele table
// ----------------------------------------------------------------------------------------
void Population::resetPopNeutralTables() {
	for (auto& entry : popNeutralCountTables) {
		entry.reset();
	}
}

// ----------------------------------------------------------------------------------------
// Populate population-level NEUTRAL count tables
// Update allele occurrence and heterozygosity counts, and allele frequencies
// ----------------------------------------------------------------------------------------
void Population::updatePopNeutralTables() {

	const int nLoci = pSpecies->getNPositionsForTrait(NEUTRAL);
	const int nAlleles = pSpecies->getSpTrait(NEUTRAL)->getNbNeutralAlleles();
	const auto& positions = pSpecies->getSpTrait(NEUTRAL)->getGenePositions();
	const int ploidy = pSpecies->isDiploid() ? 2 : 1;

	// Create /reset empty tables
	if (popNeutralCountTables.size() != 0)
		resetPopNeutralTables();
	else {
		popNeutralCountTables.reserve(nLoci);

		for (int l = 0; l < nLoci; l++) {
			popNeutralCountTables.push_back(NeutralCountsTable(nAlleles));
		}
	}

	// Fill tallies for each locus
	for (Individual* individual : sampledInds) {

		// direct access by locus index (positions are iterated in ascending order,
		// which is the order of the loci in the neutral genes array)
		const NeutralTrait* trait = static_cast<const NeutralTrait*>(individual->getTrait(NEUTRAL));
		for (int whichLocus = 0; whichLocus < (int)positions.size(); whichLocus++) {

			int alleleOnChromA = (int)trait->getAlleleAtIndex(0, whichLocus);
			popNeutralCountTables[whichLocus].incrementTally(alleleOnChromA);

			if (ploidy == 2) { // second allele and heterozygosity
				int alleleOnChromB = (int)trait->getAlleleAtIndex(1, whichLocus);
				popNeutralCountTables[whichLocus].incrementTally(alleleOnChromB);

				bool isHetero = alleleOnChromA != alleleOnChromB;
				if (isHetero) {
					popNeutralCountTables[whichLocus].incrementHeteroTally(alleleOnChromA);
					popNeutralCountTables[whichLocus].incrementHeteroTally(alleleOnChromB);
				}
			}
		}
	}

	// Fill frequencies
	if (sampledInds.size() > 0) {
		std::for_each(
			popNeutralCountTables.begin(),
			popNeutralCountTables.end(),
			[&](NeutralCountsTable& thisLocus) -> void {
				thisLocus.setFrequencies(static_cast<int>(sampledInds.size()) * ploidy);
			});
	}
}

double Population::getAlleleFrequency(int thisLocus, int whichAllele) {
	return popNeutralCountTables[thisLocus].getFrequency(whichAllele);
}

int Population::getAlleleTally(int thisLocus, int whichAllele) {
	return popNeutralCountTables[thisLocus].getTally(whichAllele);
}

int Population::getHeteroTally(int thisLocus, int whichAllele) {
	return popNeutralCountTables[thisLocus].getHeteroTally(whichAllele);
}

// ----------------------------------------------------------------------------------------
// Count number of heterozygotes loci in sampled individuals
// ----------------------------------------------------------------------------------------
int Population::countHeterozygoteLoci() {
	int nbHetero = 0;
	if (pSpecies->isDiploid()) {
		for (Individual* ind : sampledInds) {
			const NeutralTrait* trait = (NeutralTrait*)(ind->getTrait(NEUTRAL));
			nbHetero += trait->countHeterozygoteLoci();
		}
	}
	return nbHetero;
}

// ----------------------------------------------------------------------------------------
// Count number of heterozygotes among sampled individuals for each locus
// ----------------------------------------------------------------------------------------
vector<int> Population::countNbHeterozygotesEachLocus() {
	const auto& positions = pSpecies->getSpTrait(NEUTRAL)->getGenePositions();
	vector<int> hetero(positions.size(), 0);

	if (pSpecies->isDiploid()) {
		for (Individual* ind : sampledInds) {
			const NeutralTrait* trait = (NeutralTrait*)ind->getTrait(NEUTRAL);
			for (int counter = 0; counter < (int)positions.size(); counter++) {
				hetero[counter] += trait->isHeterozygoteAtIndex(counter);
			}
		}
	}
	return hetero;
}

// ----------------------------------------------------------------------------------------
//	Compute the expected heterozygosity for population
// ----------------------------------------------------------------------------------------
double Population::computeHs() {
	int nLoci = pSpecies->getNPositionsForTrait(NEUTRAL);
	int nAlleles = pSpecies->getSpTrait(NEUTRAL)->getNbNeutralAlleles();
	double hs = 0;
	double freq;
	vector<double> locihet(nLoci, 1);

	if (sampledInds.size() > 0) {
		for (int thisLocus = 0; thisLocus < nLoci; ++thisLocus) {
			for (int allele = 0; allele < nAlleles; ++allele) {
				freq = getAlleleFrequency(thisLocus, allele);
				freq *= freq; //squared frequencies (expected _homozygosity)
				locihet[thisLocus] -= freq; // 1 - sum of p2 = expected heterozygosity
			}
			hs += locihet[thisLocus];
		}
	}
	return hs;
}

popStats Population::getStats(const std::vector <float>& localDemoScaling)
{
	popStats p = popStats();
	int ninds;
	float fec;
	bool breeders[2] = { false, false };
	demogrParams dem = pSpecies->getDemogrParams();
	p.pSpecies = pSpecies;
	p.pPatch = pPatch;
	p.spNum = pSpecies->getSpNum();
	p.nInds = (int)inds.size();
	p.nNonJuvs = p.nAdults = 0;
	p.breeding = false;
	for (int stg = 1; stg < nStages; stg++) {
		for (int sex = 0; sex < nSexes; sex++) {
			ninds = nInds[stg][sex];
			p.nNonJuvs += ninds;
			if (ninds > 0) {
				if (pSpecies->stageStructured()) {
					if (dem.repType == 2) {
						if (pSpecies->getFecSpatial() && pSpecies->getFecLayer(stg,sex)>=0){
						fec = pSpecies->getFec(stg,sex)*localDemoScaling[pSpecies->getFecLayer(stg,sex)];
						}
						else fec = pSpecies->getFec(stg,sex);
					}
					else {
						if (pSpecies->getFecSpatial() && pSpecies->getFecLayer(stg,0)>=0){
							fec = pSpecies->getFec(stg,0)*localDemoScaling[pSpecies->getFecLayer(stg,0)];
						}
					else fec = pSpecies->getFec(stg, 0);
					}
					if (fec > 0.0) { breeders[sex] = true; p.nAdults += ninds; }
				}
				else breeders[sex] = true;
			}
		}
	}
	// is there a breeding population present?
	if (nSexes == 1) {
		p.breeding = breeders[0];
	}
	else {
		if (breeders[0] && breeders[1]) p.breeding = true;
	}
	return p;
}

Species* Population::getSpecies(void) { return pSpecies; }

int Population::getNbInds() const {
	return inds.size();
}

int Population::getNbInds(int stg) const {
	int t = 0;
	if (stg < 0 || stg >= nStages) throw runtime_error("Attempt to get nb individuals for stage " + to_string(stg) + ", no such stage.");
	for (int sex = 0; sex < nSexes; sex++) {
		t += nInds[stg][sex];
	}
	return t;
}

int Population::getNbInds(int stg, int sex) const {
	if (stg < 0 || stg >= nStages) throw runtime_error("Attempt to get nb individuals for stage " + to_string(stg) + ", no such stage.");
	return nInds[stg][sex];
}

//---------------------------------------------------------------------------
// Remove all Individuals
void Population::extirpate(void) {
	int ninds = (int)inds.size();
	for (int i = 0; i < ninds; i++) {
		if (inds[i] != NULL) delete inds[i];
	}
	inds.clear();
	int njuvs = (int)juvs.size();
	for (int i = 0; i < njuvs; i++) {
		if (juvs[i] != NULL) delete juvs[i];
	}
	juvs.clear();
	for (int sex = 0; sex < nSexes; sex++) {
		for (int stg = 0; stg < nStages; stg++) {
			nInds[stg][sex] = 0;
		}
	}
}

//---------------------------------------------------------------------------
// Produce juveniles and hold them in the juvs vector
void Population::reproduction(const float localK, const float envval, const int resol, const std::vector <float>& localDemoScaling)
{

	// get population size at start of reproduction
	int ninds = (int)inds.size();
	if (ninds == 0) return;

	int nsexes, stage, sex, njuvs, nj, nmales, nfemales;
	Cell* pCell;
	indStats ind;
	double meanFecundity;
	bool skipbreeding;

	envStochParams env = paramsStoch->getStoch();
	demogrParams dem = pSpecies->getDemogrParams();
	stageParams sstruct = pSpecies->getStageParams();
	emigRules emig = pSpecies->getEmigRules();
	transferRules trfr = pSpecies->getTransferRules();
	settleType sett = pSpecies->getSettle();

	float fecStdDev = dem.fecSD;

	if (dem.repType == 0)
		nsexes = 1;
	else nsexes = 2;


	// set up local copy of species fecundity table
	float fec[gMaxNbStages][gMaxNbSexes];
	for (int stg = 0; stg < sstruct.nStages; stg++) {
		for (int sex = 0; sex < nsexes; sex++) {
			if (dem.stageStruct) {
				if (dem.repType == 1) { // simple sexual model
					// both sexes use fecundity recorded for females
					if (pSpecies->getFecSpatial() && pSpecies->getFecLayer(stg,0)>=0){
						fec[stg][sex] = pSpecies->getFec(stg,0)*localDemoScaling[pSpecies->getFecLayer(stg,0)];
				}
					else fec[stg][sex] = pSpecies->getFec(stg,0);
				}
				else {
					if (pSpecies->getFecSpatial() && pSpecies->getFecLayer(stg,sex)>=0){
						fec[stg][sex] = pSpecies->getFec(stg,sex)*localDemoScaling[pSpecies->getFecLayer(stg,sex)];
					}
				else fec[stg][sex] = pSpecies->getFec(stg, sex);
			}
			}
			else { // non-structured population
				if (stg == 1) fec[stg][sex] = dem.lambda; // adults
				else fec[stg][sex] = 0.0; // juveniles
			}
		}
	}

	if (dem.stageStruct) {
		// apply environmental effects and density dependence
		// to all non-zero female non-juvenile stages
		for (int stg = 1; stg < nStages; stg++) {
			if (fec[stg][0] > 0.0) {
				// apply any effect of environmental gradient and/or stochasticty
				fec[stg][0] *= envval;
				if (env.stoch && !env.inK) {
					// fecundity (at low density) is constrained to lie between limits specified
					// for the species
					float limit;
					limit = pSpecies->getMinMax(0);
					if (fec[stg][0] < limit) fec[stg][0] = limit;
					limit = pSpecies->getMinMax(1);
					if (fec[stg][0] > limit) fec[stg][0] = limit;
				}
				if (sstruct.fecDens) { // apply density dependence
					float effect = 0.0;
					if (sstruct.fecStageDens) { // stage-specific density dependence
						// NOTE: matrix entry [row][col] is the effect of the COLUMN stage (acting) on the
						// ROW stage (affected), w_ij = effect of j on i, as in the R documentation;
						// in the complex sexual model the matrix index puts males before females
						// (2*stage = m, 2*stage+1 = f), although internally sex 0 = female, 1 = male
						float weight = 0.0;
						for (int effstg = 0; effstg < nStages; effstg++) {
							for (int effsex = 0; effsex < nSexes; effsex++) {
								if (dem.repType == 2) {
									if (effsex == 0) weight = pSpecies->getDDwtFec(2 * stg + 1, 2 * effstg + 1);
									else weight = pSpecies->getDDwtFec(2 * stg + 1, 2 * effstg);
								}
								else {
									weight = pSpecies->getDDwtFec(stg, effstg);
								}
								effect += (float)nInds[effstg][effsex] * weight;
							}
						}
					}
					else // not stage-specific
						effect = (float)getNbInds();
					if (localK > 0.0) fec[stg][0] *= exp(-effect / localK);
				}
			}
		}
	}
	else { // non-structured - set fecundity for adult females only
		// apply any effect of environmental gradient and/or stochasticty
		fec[1][0] *= envval;
		if (env.stoch && !env.inK) {
			// fecundity (at low density) is constrained to lie between limits specified
			// for the species
			float limit;
			limit = pSpecies->getMinMax(0);
			if (fec[1][0] < limit) fec[1][0] = limit;
			limit = pSpecies->getMinMax(1);
			if (fec[1][0] > limit) fec[1][0] = limit;
		}
		// apply density dependence
		if (localK > 0.0) {
			if (dem.repType == 1 || dem.repType == 2) { // sexual model
				// apply factor of 2 (as in manual, eqn. 6)
				fec[1][0] *= 2.0;
			}
			fec[1][0] /= (1.0f + fabs(dem.lambda - 1.0f) * pow(((float)ninds / localK), dem.bc));
		}
	}

	double propBreed;
	Individual* father = nullptr;
	std::vector <Individual*> fathers;

	switch (dem.repType) {

	case 0: // asexual model
		for (int i = 0; i < ninds; i++) {
			prefetchInd(inds, i);
			stage = inds[i]->breedingFem();
			if (stage > 0) { // female of breeding age
				if (dem.stageStruct) {
					// determine whether she must miss current breeding attempt
					ind = inds[i]->getStats();
					if (ind.fallow >= sstruct.repInterval) {
						if (pRandom->Bernoulli(sstruct.probRep)) skipbreeding = false;
						else skipbreeding = true;
					}
					else skipbreeding = true; // cannot breed this time
				}
				else skipbreeding = false; // not structured - always breed
				if (skipbreeding) {
					inds[i]->incFallow();
				}
				else { // attempt to breed
					inds[i]->resetFallow();

					meanFecundity = fec[stage][0];
					njuvs = drawNbOffspring(meanFecundity, fecStdDev, dem.maxOffspring);

					pCell = pPatch->getRandomCell();
					for (int j = 0; j < njuvs; j++) {

						Individual* newJuv;
						// asexual (female-only) model: offspring are always female
						newJuv = new Individual(pSpecies, pCell, pPatch, 0, 0, 0, 0.0, trfr.usesMovtProc, trfr.moveType);

						if (pSpecies->getNTraits() > 0) {
							newJuv->inheritTraits(pSpecies, inds[i], resol);
						}

						if (!newJuv->isViable()) {
							delete newJuv;
						}
						else {
							juvs.push_back(newJuv);
						nInds[0][0]++;
						}
					}
				}
			}
		}
		break;

	case 1: // simple sexual model
	case 2: // complex sexual model
		// count breeding females and males
		// add breeding males to list of potential fathers
		nfemales = nmales = 0;
		for (int i = 0; i < ninds; i++) {
			prefetchInd(inds, i);
			ind = inds[i]->getStats();
			if (ind.sex == 0 && fec[ind.stage][0] > 0.0) nfemales++;
			if (ind.sex == 1 && fec[ind.stage][1] > 0.0) {
				fathers.push_back(inds[i]);
				nmales++;
			}
		}
		if (nfemales > 0 && nmales > 0)
		{ // population can breed
			if (dem.repType == 2) { // complex sexual model
				// calculate proportion of eligible females which breed
				propBreed = (2.0 * dem.harem * nmales) / (nfemales + dem.harem * nmales);
				if (propBreed > 1.0) propBreed = 1.0;
			}
			else propBreed = 1.0;
			for (int i = 0; i < ninds; i++) {
				prefetchInd(inds, i);
				stage = inds[i]->breedingFem();
				if (stage > 0 && fec[stage][0] > 0.0) { // (potential) breeding female
					if (dem.stageStruct) {
						// determine whether she must miss current breeding attempt
						ind = inds[i]->getStats();
						if (ind.fallow >= sstruct.repInterval) {
							if (pRandom->Bernoulli(sstruct.probRep)) skipbreeding = false;
							else skipbreeding = true;
						}
						else skipbreeding = true; // cannot breed this time
					}
					else skipbreeding = false; // not structured - always breed
					if (skipbreeding) {
						inds[i]->incFallow();
					}
					else { // attempt to breed
						inds[i]->resetFallow();
						// NOTE: FOR COMPLEX SEXUAL MODEL, NO. OF FEMALES *ACTUALLY* BREEDING DOES NOT
						// NECESSARILY EQUAL THE EXPECTED NO. FROM EQN. 7 IN THE MANUAL...
						if (pRandom->Bernoulli(propBreed)) {
							meanFecundity = fec[stage][0]; // breeds
							njuvs = drawNbOffspring(meanFecundity, fecStdDev, dem.maxOffspring);
						}
						else njuvs = 0;

						if (njuvs > 0) {
							// select father at random from breeding males ...
							int rrr = 0;
							if (nmales > 1) rrr = pRandom->IRandom(0, nmales - 1);
							father = fathers[rrr];
							pCell = pPatch->getRandomCell();

							for (int j = 0; j < njuvs; j++) {
								Individual* newJuv;

								newJuv = new Individual(pSpecies, pCell, pPatch, 0, 0, 0, dem.propMales, trfr.usesMovtProc, trfr.moveType);

								if (pSpecies->getNTraits() > 0) {
									newJuv->inheritTraits(pSpecies, inds[i], father, resol);
								}

								if (!newJuv->isViable()) {
									delete newJuv;
								}
								else {
									juvs.push_back(newJuv);
									sex = newJuv->getSex();
								nInds[0][sex]++;
								}
							}
						}
					}
				}
			}
		}
		fathers.clear();
		break;

	} // end of switch (dem.repType)

// THIS MAY NOT BE CORRECT FOR MULTIPLE SPECIES IF THERE IS SOME FORM OF
// CROSS-SPECIES DENSITY-DEPENDENT FECUNDITY
}

// Following reproduction of ALL species, add juveniles to the population prior to dispersal
void Population::fledge(void)
{
	demogrParams dem = pSpecies->getDemogrParams();

	if (dem.stageStruct) { // juveniles are added to the individuals vector
		inds.insert(inds.end(), juvs.begin(), juvs.end());
	}
	else { // all adults die and juveniles replace adults
		int ninds = (int)inds.size();
		for (int i = 0; i < ninds; i++) {
			delete inds[i];
		}
		inds.clear();
		for (int sex = 0; sex < nSexes; sex++) {
			nInds[1][sex] = 0; // set count of adults to zero
		}
		inds = std::move(juvs);
	}
	juvs.clear();
}

Individual* Population::sampleInd() const {
	if (inds.empty()) return nullptr;
	int index = pRandom->IRandom(0, static_cast<int>(inds.size() - 1));
	return inds[index];
}

void Population::sampleIndsWithoutReplacement(string strNbToSample, const set<int>& sampleStages) {

	if (sampledInds.size() > 0) {
		sampledInds.clear();
	}
	auto& rng = pRandom->getRNG();// bugfix: auto& - plain auto would copy the generator again
	vector<Individual*> stagedInds;

	// Stage individuals in eligible stages
	for (int stage : sampleStages) {
		vector<Individual*> toAdd = getIndividualsInStage(stage);
		// exclude individuals already doomed this season (status 6-9; removed by survival1())
		erase_if(toAdd, [](Individual* ind) { const int st = ind->getStatus(); return st >= 6 && st != 10; });
		stagedInds.insert(stagedInds.begin(), toAdd.begin(), toAdd.end());
	}

	if (strNbToSample == "all") {
		// Sample all individuals in selected stages
		sampledInds = stagedInds;
	}
	else { // random
		int nbToSample = stoi(strNbToSample);
		if (stagedInds.size() <= nbToSample) {
			// Sample all individuals in selected stages
			sampledInds = stagedInds;
		}
		else {
			// Sample n individuals across selected stages
			sample(stagedInds.begin(), stagedInds.end(), std::back_inserter(sampledInds), nbToSample, rng);
		}
	}
}

int Population::sampleSize() const {
	return static_cast<int>(sampledInds.size());
}

vector<Individual*> Population::getIndividualsInStage(int stage) {
	vector<Individual*> indsInStage;
	for (auto ind : inds) {
		if (ind->getStats().stage == stage)
			indsInStage.push_back(ind);
	}
	return indsInStage;
}

// Determine which individuals will disperse
void Population::emigration(float localK)
{
	int nsexes;
	double disp, pbDisp, NK;
	demogrParams dem = pSpecies->getDemogrParams();
	stageParams sstruct = pSpecies->getStageParams();
	emigRules emig = pSpecies->getEmigRules();
	emigTraits eparams;
	transferRules trfr = pSpecies->getTransferRules();
	indStats ind;

// to avoid division by zero, assume carrying capacity is at least one individual
// localK can be zero if there is a moving gradient or stochasticity in K
	if (localK < 1.0) localK = 1.0;
	NK = static_cast<float>(getNbInds()) / localK;

	int ninds = static_cast<int>(inds.size());

	// set up local copy of emigration probability table
	// used when there is no individual variability
	// NB - IT IS DOUBTFUL THIS CONTRIBUTES ANY SUBSTANTIAL TIME SAVING
	if (dem.repType == 0) nsexes = 1;
	else nsexes = 2;
	double pbEmig[gMaxNbStages][gMaxNbSexes];

	for (int stg = 0; stg < sstruct.nStages; stg++) {
		for (int sex = 0; sex < nsexes; sex++) {
			if (emig.indVar) pbEmig[stg][sex] = 0.0;
			else {
				if (emig.densDep) {
					if (emig.sexDep) {
						if (emig.stgDep) {
							eparams = pSpecies->getSpEmigTraits(stg, sex);
						}
						else {
							eparams = pSpecies->getSpEmigTraits(0, sex);
						}
					}
					else { // !emig.sexDep
						if (emig.stgDep) {
							eparams = pSpecies->getSpEmigTraits(stg, 0);
						}
						else {
							eparams = pSpecies->getSpEmigTraits(0, 0);
						}
					}
					pbEmig[stg][sex] = eparams.d0 / (1.0 + exp(-(NK - eparams.beta) * eparams.alpha));
				}
				else { // density-independent
					if (emig.sexDep) {
						if (emig.stgDep) {
							pbEmig[stg][sex] = pSpecies->getSpEmigD0(stg, sex);
						}
						else { // !emig.stgDep
							pbEmig[stg][sex] = pSpecies->getSpEmigD0(0, sex);
						}
					}
					else { // !emig.sexDep
						if (emig.stgDep) {
							pbEmig[stg][sex] = pSpecies->getSpEmigD0(stg, 0);
						}
						else { // !emig.stgDep
							pbEmig[stg][sex] = pSpecies->getSpEmigD0(0, 0);
						}
					}
				}
			} // end of !emig.indVar
		}
	}

	for (int i = 0; i < ninds; i++) {
		prefetchInd(inds, i);
		ind = inds[i]->getStats();
		if (ind.status < 1) // ToDo: Maybe allow dispersal after translocation? If so, we need to update the pPrevCell and pCurrCell variables of the translocated individuals!
		{
			if (emig.indVar) { // individual variability in emigration
				if (dem.stageStruct && ind.stage != emig.emigStage) {
					// emigration may not occur
					pbDisp = 0.0;
				}
				else { // non-structured or individual is in emigration stage
					eparams = inds[i]->getIndEmigTraits();
					if (emig.densDep) { // density-dependent
						NK = (float)getNbInds() / localK;
						pbDisp = eparams.d0 / (1.0 + exp(-(NK - eparams.beta) * eparams.alpha));
					}
					else { // density-independent
						if (emig.sexDep) {
							pbDisp = pbEmig[0][ind.sex] + eparams.d0;
						}
						else {
							pbDisp = pbEmig[0][0] + eparams.d0;
						}
					}
				}
			} // end of individual variability
			else { // no individual variability

				if (emig.densDep) {
					if (emig.sexDep) {
						if (emig.stgDep) {
							pbDisp = pbEmig[ind.stage][ind.sex];
						}
						else {
							pbDisp = pbEmig[0][ind.sex];
						}
					}
					else { // !emig.sexDep
						if (emig.stgDep) {
							pbDisp = pbEmig[ind.stage][0];
						}
						else {
							pbDisp = pbEmig[0][0];
						}
					}
				}
				else { // density-independent
					if (emig.sexDep) {
						if (emig.stgDep) {
							pbDisp = pbEmig[ind.stage][ind.sex];
						}
						else { // !emig.stgDep
							pbDisp = pbEmig[0][ind.sex];
						}
					}
					else { // !emig.sexDep
						if (emig.stgDep) {
							pbDisp = pbEmig[ind.stage][0];
						}
						else { // !emig.stgDep
							pbDisp = pbEmig[0][0];
						}
					}
				}
			} // end of no individual variability
			disp = pRandom->Bernoulli(pbDisp);

			if (disp == 1) { // emigrant
				inds[i]->setStatus(1);
			}
		} // end of if (ind.status < 1) condition
	} // end of for loop
}

// All individuals emigrate after patch destruction
void Population::allEmigrate(void) {
	int ninds = (int)inds.size();
	for (int i = 0; i < ninds; i++) {
		inds[i]->setStatus(1);
	}
}

// Remove an Individual from the Population
Individual* Population::extractIndividual(int ix) {
	Individual* pInd = inds[ix];
	indStats ind = pInd->getStats();
	inds[ix] = nullptr;
	nInds[ind.stage][ind.sex]--;
	return pInd;
}

// If an Individual has been identified as an emigrant, remove it from the Population
disperser Population::extractDisperser(int ix) {
	disperser d = disperser();
	indStats ind = inds[ix]->getStats();
	if (ind.status == 1) { // emigrant
		d.pInd = inds[ix];
		d.yes = true;
		inds[ix] = nullptr;
		nInds[ind.stage][ind.sex]--;
	}
	else {
		d.pInd = NULL; d.yes = false;
	}
	return d;
}


// For an individual identified as being in the matrix population:
// if it is a settler, return its new location and remove it from the current population
// otherwise, leave it in the matrix population for possible reporting before deletion
disperser Population::extractSettler(int ix) {
	disperser d = disperser();
	Cell* pCell;

	indStats ind = inds[ix]->getStats();

	pCell = inds[ix]->getCurrCell();
	d.pInd = inds[ix];  d.pCell = pCell; d.yes = false;
	if (ind.status == 4 || ind.status == 5) { // settled
		d.yes = true;
		inds[ix] = 0;
		nInds[ind.stage][ind.sex]--;
	}
	return d;
}

// Add a specified individual to the new/current dispersal group
// Add a specified individual to the population
void Population::recruit(Individual* pInd) {
	indStats ind = pInd->getStats(); // potentially I need to add localscalings or so?
	nInds[ind.stage][ind.sex]++;
#ifdef _OPENMP
	const std::lock_guard<std::mutex> lock(inds_mutex);
#endif // _OPENMP
	inds.push_back(pInd);
}

// Add specified individuals to the population
void Population::recruitMany(std::vector<Individual*>& recruits) {
	if (recruits.empty()) return;
	for (Individual* pInd : recruits) {
		indStats ind = pInd->getStats();
		nInds[ind.stage][ind.sex]++;
	}
#ifdef _OPENMP
	const std::lock_guard<std::mutex> lock(inds_mutex);
#endif // _OPENMP
	inds.insert(inds.end(), recruits.begin(), recruits.end());
}

//---------------------------------------------------------------------------
// Determine survival and development and record in individual's status code
// Changes are NOT applied to the Population at this stage

// FOR MULTIPLE SPECIES, MAY NEED TO SEPARATE OUT THIS IDENTIFICATION STAGE,
// SO THAT IT CAN BE PERFORMED FOR ALL SPECIES BEFORE ANY UPDATING OF POPULATIONS

void Population::survival0(float localK, short option0, short option1, const std::vector <float>& localDemoScaling)
{
	// option0:	0 - stage 0 (juveniles) only
	//					1 - all stages
	//					2 - stage 1 and above (all non-juveniles)
	//
	// option1:	0 - development only (when survival is annual)
	//	  	 		1 - development and survival
	//	  	 		2 - survival only (when survival is annual)
	densDepParams ddparams = pSpecies->getDensDep();
	demogrParams dem = pSpecies->getDemogrParams();
	stageParams sstruct = pSpecies->getStageParams();

	// get current population size
	int ninds = inds.size();
	if (ninds == 0) return;

	// set up local copies of species development and survival tables
	int nsexes = dem.repType == 0 ? 1 : 2;
	float dev[gMaxNbStages][gMaxNbSexes];
	float surv[gMaxNbStages][gMaxNbSexes];
	short minAge[gMaxNbStages][gMaxNbSexes];

	for (int stg = 0; stg < sstruct.nStages; stg++) {
		for (int sex = 0; sex < nsexes; sex++) {
			if (dem.stageStruct) {
				if (dem.repType == 1) { // simple sexual model
					// both sexes use development and survival recorded for females
					if (pSpecies->getDevSpatial() && pSpecies->getDevLayer(stg,0)>=0){
						dev[stg][sex] = pSpecies->getDev(stg,0)*localDemoScaling[pSpecies->getDevLayer(stg,0)];
					} 
					else dev[stg][sex] = pSpecies->getDev(stg,0);
					if (pSpecies->getSurvSpatial() && pSpecies->getSurvLayer(stg,0)>=0){
						surv[stg][sex] = pSpecies->getSurv(stg,0)*localDemoScaling[pSpecies->getSurvLayer(stg,0)];
					}
					else surv[stg][sex] = pSpecies->getSurv(stg,0);
					minAge[stg][sex] = pSpecies->getMinAge(stg, 0);
				}
				else {
					if (pSpecies->getDevSpatial() && pSpecies->getDevLayer(stg,sex)>=0){
						dev[stg][sex] = pSpecies->getDev(stg,sex)*localDemoScaling[pSpecies->getDevLayer(stg,sex)];
					}
					else dev[stg][sex] = pSpecies->getDev(stg,sex);
					if (pSpecies->getSurvSpatial() && pSpecies->getSurvLayer(stg,sex)>=0){
										surv[stg][sex] = pSpecies->getSurv(stg,sex)*localDemoScaling[pSpecies->getSurvLayer(stg,sex)];
					}
					else surv[stg][sex] = pSpecies->getSurv(stg,sex);
					minAge[stg][sex] = pSpecies->getMinAge(stg, sex);
				}
				if (option1 == 0) surv[stg][sex] = 1.0; // development only - all survive
				if (option1 == 2) dev[stg][sex] = 0.0;  // survival only - none develops
			}
			else { // non-structured population
				if (stg == 1) { // adults
					dev[stg][sex] = 0.0; surv[stg][sex] = 0.0; minAge[stg][sex] = 0;
				}
				else { // juveniles
					dev[stg][sex] = 1.0; surv[stg][sex] = 1.0; minAge[stg][sex] = 0;
				}
			}
		}
	}
	if (dem.stageStruct) {
		for (int stg = 0; stg < nStages; stg++) {
			for (int sex = 0; sex < nsexes; sex++) {
				if (option1 != 2 && sstruct.devDens && stg > 0) {
				// NB DD in development does NOT apply to juveniles,
					float effect = 0.0;
					if (sstruct.devStageDens) { // stage-specific density dependence
						// NOTE: matrix entry [row][col] is the effect of the COLUMN stage (acting) on the
						// ROW stage (affected), w_ij = effect of j on i, as in the R documentation;
						// in the complex sexual model the matrix index puts males before females
						// (2*stage = m, 2*stage+1 = f), although internally sex 0 = female, 1 = male
						float weight = 0.0;
						for (int effstg = 0; effstg < nStages; effstg++) {
							for (int effsex = 0; effsex < nSexes; effsex++) {
								if (dem.repType == 2) {
									int rowincr, colincr;
									if (effsex == 0) rowincr = 1; else rowincr = 0;
									if (sex == 0) colincr = 1; else colincr = 0;
									weight = pSpecies->getDDwtDev(2 * stg + colincr, 2 * effstg + rowincr);
								}
								else {
									weight = pSpecies->getDDwtDev(stg, effstg);
								}
								effect += (float)nInds[effstg][effsex] * weight;
							}
						}
					}
					else // not stage-specific
						effect = (float)getNbInds();
					if (localK > 0.0)
						dev[stg][sex] *= exp(-(ddparams.devCoeff * effect) / localK);
				} // end of if (sstruct.devDens && stg > 0)
				if (option1 != 0 && sstruct.survDens) {

					float effect = 0.0;
					if (sstruct.survStageDens) { // stage-specific density dependence
						// NOTE: matrix entry [row][col] is the effect of the COLUMN stage (acting) on the
						// ROW stage (affected), w_ij = effect of j on i, as in the R documentation;
						// in the complex sexual model the matrix index puts males before females
						// (2*stage = m, 2*stage+1 = f), although internally sex 0 = female, 1 = male
						float weight = 0.0;
						for (int effstg = 0; effstg < nStages; effstg++) {
							for (int effsex = 0; effsex < nSexes; effsex++) {
								if (dem.repType == 2) {
									int rowincr, colincr;
									if (effsex == 0) rowincr = 1; else rowincr = 0;
									if (sex == 0) colincr = 1; else colincr = 0;
									weight = pSpecies->getDDwtSurv(2 * stg + colincr, 2 * effstg + rowincr);
								}
								else {
									weight = pSpecies->getDDwtSurv(stg, effstg);
								}
								effect += (float)nInds[effstg][effsex] * weight;
							}
						}
					}
					else // not stage-specific
						effect = (float)getNbInds();
					if (localK > 0.0)
						surv[stg][sex] *= exp(-(ddparams.survCoeff * effect) / localK);
				} // end of if (sstruct.survDens)
			}
		}
	}
	// identify which individuals die or develop
	for (int i = 0; i < ninds; i++) {
		prefetchInd(inds, i);
		indStats ind = inds[i]->getStats();

		if ((ind.stage == 0 && option0 < 2) || (ind.stage > 0 && option0 > 0)) {
			// condition for processing the stage is met...
			if (ind.status < 6 || ind.status == 10) { // not already doomed
				double probsurv = surv[ind.stage][ind.sex];
				// does the individual survive?
				if (pRandom->Bernoulli(probsurv)) { // survives
					// does the individual develop?
					double probdev = dev[ind.stage][ind.sex];
					if (ind.stage < nStages - 1) { // not final stage
						if (ind.age >= minAge[ind.stage + 1][ind.sex]) { // old enough to enter next stage
							if (pRandom->Bernoulli(probdev)) {
								inds[i]->setToDevelop();
							}
						}
					}
				}
				else { // doomed to die
					inds[i]->setStatus(8);
				}
			}
		}
	}
}

// Apply survival changes to the population
void Population::survival1(void)
{
	// individuals harvested in this season have been written to the individual output by now
	for (Individual* pInd : harvestedInds) delete pInd;
	harvestedInds.clear();

	int ninds = (int)inds.size();

	for (int i = 0; i < ninds; i++) {
		prefetchInd(inds, i);
		indStats ind = inds[i]->getStats();

		if (ind.status > 5 && ind.status != 10) { // doomed to die; status 10 is translocated
			if (ind.status != 10) //not going into cold storage 
			delete inds[i];
			inds[i] = nullptr;
			nInds[ind.stage][ind.sex]--;
		}
		else {
			if (ind.isDeveloping) { // develops to next stage
				nInds[ind.stage][ind.sex]--;
				inds[i]->develop();
				nInds[ind.stage + 1][ind.sex]++;
			}
		}
	}
	clean();
}

void Population::ageIncrement(void) {
	int ninds = (int)inds.size();
	stageParams sstruct = pSpecies->getStageParams();
	for (int i = 0; i < ninds; i++) {
		prefetchInd(inds, i);
		inds[i]->ageIncrement(sstruct.maxAge);
	}
}

//---------------------------------------------------------------------------
// Remove zero pointers to dead or dispersed individuals
void Population::clean(void)
{
	int ninds = (int)inds.size();
	if (ninds > 0) {
		inds.erase(std::remove(inds.begin(), inds.end(), (Individual *)NULL), inds.end());
#if RS_RCPP
		shuffle(inds.begin(), inds.end(), pRandom->getRNG());
#else

#ifdef NDEBUG
		// do not randomise individuals in DEBUG mode, as the function uses rand()
		// and therefore the randomisation will differ between identical runs of RS
		shuffle(inds.begin(), inds.end(), pRandom->getRNG());
#endif // NDEBUG

#endif // RS_RCPP
	}
}

//---------------------------------------------------------------------------
// Close population file
bool Population::outPopFinishLandscape() {
	outPop.close();
	return true;
}

//---------------------------------------------------------------------------
// Open population file and write header record
bool Population::outPopStartLandscape(int landNr, bool patchModel) {

	string name;
	simParams sim = paramsSim->getSim();
	envGradParams grad = paramsGrad->getGradient();

	// NEED TO REPLACE CONDITIONAL COLUMNS BASED ON ATTRIBUTES OF ONE SPECIES TO COVER
	// ATTRIBUTES OF *ALL* SPECIES AS DETECTED AT MODEL LEVEL
	demogrParams dem = pSpecies->getDemogrParams();
	stageParams sstruct = pSpecies->getStageParams();
	// file name without extension (.txt or .h5 is added by OutputTable)
	if (sim.batchMode) {
		name = paramsSim->getDir(2)
			+ "Batch" + to_string(sim.batchNum) + "_"
			+ "Sim" + to_string(sim.simulation) + "_Land" + to_string(landNr) + "_Pop";
	}
	else {
		name = paramsSim->getDir(2) + "Sim" + to_string(sim.simulation) + "_Pop";
	}
	typedef OutputTable T;
	outPop.clearColumns();
	outPop.addColumn("Rep", T::INT); outPop.addColumn("Year", T::INT); outPop.addColumn("RepSeason", T::INT);
	if (patchModel) { outPop.addColumn("PatchID", T::INT); outPop.addColumn("Ncells", T::INT); }
	else { outPop.addColumn("x", T::INT); outPop.addColumn("y", T::INT); }
	// determine whether environmental data need be written for populations
	bool writeEnv = false;
	if (grad.gradient) writeEnv = true;
	if (paramsStoch->envStoch()) writeEnv = true;
	if (writeEnv) {
		outPop.addColumn("Epsilon", T::FLOAT); outPop.addColumn("Gradient", T::FLOAT); outPop.addColumn("Local_K", T::FLOAT);
	}
	outPop.addColumn("Species", T::INT); outPop.addColumn("NInd", T::INT);
	if (dem.stageStruct) {
		if (dem.repType == 0)
		{
			for (int i = 1; i < sstruct.nStages; i++) outPop.addColumn("NInd_stage" + to_string(i), T::INT);
			outPop.addColumn("NJuvs", T::INT);
		}
		else {
			for (int i = 1; i < sstruct.nStages; i++) {
				outPop.addColumn("Nfemales_stage" + to_string(i), T::INT);
				outPop.addColumn("Nmales_stage" + to_string(i), T::INT);
		}
			outPop.addColumn("NJuvFemales", T::INT); outPop.addColumn("NJuvMales", T::INT);
	}
	}
	else {
		if (dem.repType != 0) { outPop.addColumn("Nfemales", T::INT); outPop.addColumn("Nmales", T::INT); }
	}
	return outPop.open(name, static_cast<OutputFormat>(sim.outputFormat));
}

//---------------------------------------------------------------------------
// Write record to population file
void Population::outPopulation(int rep, int yr, int gen, float eps,
	bool patchModel, bool writeEnv, bool gradK)
{
	Cell* pCell;

// NEED TO REPLACE CONDITIONAL COLUMNS BASED ON ATTRIBUTES OF ONE SPECIES TO COVER
	demogrParams dem = pSpecies->getDemogrParams();

	popStats p;
	outPop << rep << yr << gen;
	if (patchModel) {
		outPop << pPatch->getPatchNum();
		outPop << pPatch->getNCells();
	}
	else {
		locn loc = pPatch->getCellLocn(0);
		outPop << loc.x << loc.y;
	}
	if (writeEnv) {
		if (pPatch->getPatchNum() == 0) { // matrix
			outPop << 0 << 0 << 0;
		}
		else {
			float k = pPatch->getK();
			float envval = 0.0;
			pCell = pPatch->getRandomCell();
			if (pCell != 0) envval = pCell->getEnvVal();
			outPop << eps << envval << k;
		}
	}
	outPop << pSpecies->getSpNum();
	if (dem.stageStruct) {
		p = getStats(pPatch->getDemoScaling());
		outPop << p.nNonJuvs;
		// non-juvenile stage totals from permanent array
		for (int stg = 1; stg < nStages; stg++) {
			for (int sex = 0; sex < nSexes; sex++) {
				outPop << (int)nInds[stg][sex];
			}
		}
		// juveniles from permanent array
		for (int sex = 0; sex < nSexes; sex++) {
			outPop << (int)nInds[0][sex];
		}
	}
	else { // non-structured population
		outPop << getNbInds();
		if (dem.repType != 0)
		{ // sexual model
			outPop << (int)nInds[1][0] << (int)nInds[1][1];
		}
	}
	outPop.endRow();
}

//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// Close individuals file
void Population::outIndsFinishReplicate()
{
	outInds.close();
	return;
}

//---------------------------------------------------------------------------
// Open individuals file and write header record
void Population::outIndsStartReplicate(int rep, int landNr, bool patchModel)
{
	string name;
	demogrParams dem = pSpecies->getDemogrParams();
	emigRules emig = pSpecies->getEmigRules();
	transferRules trfr = pSpecies->getTransferRules();
	settleType sett = pSpecies->getSettle();
	simParams sim = paramsSim->getSim();

		// file name without extension (.txt or .h5 is added by OutputTable)
		name = paramsSim->getDir(2)
		+ (sim.batchMode ? "Batch" + to_string(sim.batchNum) + "_" : "")
			+ "Sim" + to_string(sim.simulation)
			+ "_Land" + to_string(landNr) + "_Rep" + to_string(rep) + "_Inds";

	typedef OutputTable T;
	auto col = [](const char* n, T::ColType t) { outInds.addColumn(n, t); };
	outInds.clearColumns();
	col("Rep", T::INT); col("Year", T::INT); col("RepSeason", T::INT); col("Species", T::INT); col("IndID", T::INT); col("Status", T::INT);
	if (patchModel) { col("Natal_patch", T::INT); col("PatchID", T::INT); }
	else { col("Natal_X", T::INT); col("Natal_Y", T::INT); col("X", T::INT); col("Y", T::INT); }
	if (dem.repType != 0) col("Sex", T::INT);
	if (dem.stageStruct) { col("Age", T::INT); col("Stage", T::INT); }
	if (pSpecies->getNbGenLoadTraits() > 0) col("ProbViable", T::FLOAT);
	if (emig.indVar) {
		if (emig.densDep) { col("D0", T::FLOAT); col("Alpha", T::FLOAT); col("Beta", T::FLOAT); }
		else col("EP", T::FLOAT);
	}
	if (trfr.indVar) {
		if (trfr.usesMovtProc) {
			if (trfr.moveType == 1) { // SMS
				col("DP", T::FLOAT); col("GB", T::FLOAT); col("AlphaDB", T::FLOAT); col("BetaDB", T::FLOAT);
			}
			if (trfr.moveType == 2) { // CRW
				col("StepLength", T::FLOAT); col("Rho", T::FLOAT);
			}
		}
		else { // kernel
			col("MeanDistI", T::FLOAT);
			if (trfr.twinKern) { col("MeanDistII", T::FLOAT); col("PKernelI", T::FLOAT); }
		}
	}
	if (sett.indVar) {
		col("S0", T::FLOAT); col("AlphaS", T::FLOAT); col("BetaS", T::FLOAT);
	}
	col("DistMoved", T::FLOAT);
#ifndef NDEBUG
	col("Nsteps", T::INT);
#else
	if (trfr.usesMovtProc) col("Nsteps", T::INT);
#endif
	outInds.open(name, static_cast<OutputFormat>(sim.outputFormat));
}

//---------------------------------------------------------------------------
// Write records to individuals file
void Population::outIndividual(Landscape* pLandscape, int rep, int yr, int gen,
	int patchNum)
{
	bool writeInd;
	pathSteps steps;
	Cell* pCell;
	landParams ppLand = pLandscape->getLandParams();
	demogrParams dem = pSpecies->getDemogrParams();
	emigRules emig = pSpecies->getEmigRules();
	transferRules trfr = pSpecies->getTransferRules();
	settleType sett = pSpecies->getSettle();
	short spNum = pSpecies->getSpNum();

	// living individuals plus those harvested in this season (status 11)
	std::vector<Individual*> toWrite(inds);
	toWrite.insert(toWrite.end(), harvestedInds.begin(), harvestedInds.end());
	int ninds = (int)toWrite.size();
	for (int i = 0; i < ninds; i++) {
		Individual* pInd = toWrite[i];
		indStats ind = pInd->getStats();
		if (yr == -1) { // write all initialised individuals
			writeInd = true;
			outInds << rep << yr << dem.repSeasons - 1;
		}
		else {
			if (dem.stageStruct && gen < 0) { // write status 9 individuals only
				if (ind.status == 9) {
					writeInd = true;
					outInds << rep << yr << dem.repSeasons - 1;
				}
				else writeInd = false;
			}
			else {
				writeInd = true;
				outInds << rep << yr << gen;
			}
		}
		if (writeInd) {
			outInds << "\t" << spNum << "\t" << pInd->getId();
			if (dem.stageStruct) outInds << "\t" << ind.status;
			else { // non-structured population
				outInds << "\t" << ind.status;
			}
			pCell = pInd->getCurrCell();
			locn loc;
			if (pCell == 0) loc.x = loc.y = -1; // beyond boundary or in no-data cell
			else loc = pCell->getLocn();
			pCell = pInd->getPrevCell();
			locn natalloc = pCell->getLocn();
			if (ppLand.patchModel) {
				outInds << "\t" << pInd->getNatalPatch()->getPatchNum();
				if (loc.x == -1) outInds << "\t-1";
				else outInds << "\t" << patchNum;
			}
			else { // cell-based model
				outInds << (float)natalloc.x << natalloc.y;
				outInds << (float)loc.x << (float)loc.y;
			}
			if (dem.repType != 0) outInds << (int)ind.sex;
			if (dem.stageStruct) outInds << ind.age << ind.stage;

			if (pSpecies->getNbGenLoadTraits() > 0) outInds << "\t" << pInd->getGeneticFitness();

			if (emig.indVar) {
				emigTraits e = pInd->getIndEmigTraits();
				if (emig.densDep) {
					outInds << e.d0 << e.alpha << e.beta;
				}
				else {
					outInds << e.d0;
				}
			} // end of if (emig.indVar)
			if (trfr.indVar) {
				if (trfr.usesMovtProc) {
					if (trfr.moveType == 1) { // SMS
						trfrSMSTraits s = pInd->getIndSMSTraits();
						outInds << s.dp << s.gb;
						outInds << s.alphaDB << s.betaDB;
					} // end of SMS
					if (trfr.moveType == 2) { // CRW
						trfrCRWTraits c = pInd->getIndCRWTraits();
						outInds << c.stepLength << c.rho;
					} // end of CRW
				}
				else { // kernel
					trfrKernelParams k = pInd->getIndKernTraits();
					if (trfr.twinKern)
					{
						outInds << k.meanDist1 << k.meanDist2 << k.probKern1;
					}
					else {
						outInds << k.meanDist1;
					}
				}
			}

			if (sett.indVar) {
				settleTraits s = pInd->getIndSettTraits();
				outInds << s.s0 << s.alpha << s.beta;
			}

			// distance moved (metres)
			if (loc.x == -1) outInds << -1;
			else {
				float d = ppLand.resol * sqrt((float)((natalloc.x - loc.x) * (natalloc.x - loc.x)
					+ (natalloc.y - loc.y) * (natalloc.y - loc.y)));
				outInds << d;
			}
#ifndef NDEBUG
			// ALWAYS WRITE NO. OF STEPS
			steps = pInd->getSteps();
			outInds << steps.year;
#else
			if (trfr.usesMovtProc) {
				steps = pInd->getSteps();
				outInds << steps.year;
			}
#endif
			outInds.endRow();
		} // end of writeInd condition

	}
}

void Population::outputGeneValues(OutputTable& ofsGenes, const int& yr, const int& gen) const {

	const bool isDiploid = pSpecies->isDiploid();
	int indID;
	float alleleOnChromA, alleleOnChromB;
	float domCoefA, domCoefB;

	// Subset traits that are selected to be output
	const set<TraitType>& traitTypes = pSpecies->getTraitTypes();
	set<TraitType> outputTraitTypes;
	for (auto trType : traitTypes) {
		if (pSpecies->getSpTrait(trType)->isOutput())
			outputTraitTypes.insert(trType);
}

	// Fetch map to positions for each trait
	// Presumably faster than fetching for every individual
	map<TraitType, const set<int>*> allGenePositions;
	for (auto trType : outputTraitTypes) {
		allGenePositions.emplace(trType, &pSpecies->getSpTrait(trType)->getGenePositions());
		}

	for (Individual* ind : sampledInds) {
		indID = ind->getId();
		for (auto trType : outputTraitTypes) {
			const set<int>& positions = *allGenePositions[trType]; // no copy per individual
			auto indTrait = ind->getTrait(trType);
			for (auto pos : positions) {
				alleleOnChromA = indTrait->getAlleleValueAtLocus(0, pos);
				if (trType == GENETIC_LOAD1 || trType == GENETIC_LOAD2 || trType == GENETIC_LOAD3 || trType == GENETIC_LOAD4 || trType == GENETIC_LOAD5) {
					domCoefA = indTrait->getDomCoefAtLocus(0, pos);
	}
	else {
					domCoefA = 0.0;
	}
				ofsGenes << yr << gen << indID << to_string(trType) << pos << alleleOnChromA << domCoefA;
				if (isDiploid) {
					alleleOnChromB = indTrait->getAlleleValueAtLocus(1, pos);
					if (trType == GENETIC_LOAD1 || trType == GENETIC_LOAD2 || trType == GENETIC_LOAD3 || trType == GENETIC_LOAD4 || trType == GENETIC_LOAD5) {
						domCoefB = indTrait->getDomCoefAtLocus(1, pos);
}
					else {
						domCoefB = 0.0;
	}
					ofsGenes << alleleOnChromB << domCoefB;
				}
				ofsGenes.endRow();
	}
		}
	}
}

// ---------------------------------------------------------------------------
// Extract all individuals of a population with certain characteristics based on age, stage and sex
// returns a set of pointers to the individuals
// ---------------------------------------------------------------------------
std::vector <Individual*> Population::getIndsWithCharacteristics( // Select a set of individuals with specified characteristics
        int min_age,	// min age (0 if not set)
        int max_age,    // max age (max age if not set)
        int stage,    // stage
        int sex,     // sex
        bool includeTranslocated // also accept translocated individuals (status 10)
){
    // get all suitable individuals based on settings
    std::vector <Individual*> filteredInds;
	int ninds = (int)inds.size();
#if RS_RCPP
    Rcpp::Rcout << "Number individuals in cell: " << ninds << endl;
#endif
    if (ninds > 0) {
        // copy ALL individuals to filteredInds
	for (int i = 0; i < ninds; i++) {
            filteredInds.push_back(inds[i]);
		}

        // check status of inividuals
        for (int i = 0; i < ninds; i++) {
            const short st = inds[i] != NULL ? inds[i]->getStats().status : 0;
            const bool accepted = st == 0 || st == 4 || st == 5 || (includeTranslocated && st == 10);
            if (inds[i] != NULL && !accepted){ // only accept individuals with status 0, 4 or 5 (not in transfer phase + not dead/harvested); status 10 (translocated) only if requested, i.e. for harvesting - a translocated individual is never translocated again
                // Rcpp::Rcout << "Status: " << inds[i]->getStats().status << endl;
                filteredInds[i] = NULL; // set it to NULL
	}
        }

        // Check minimal age
        if (min_age!=-9){
            // loop over all number of individuals in cell
            for (int i = 0; i < ninds; i++) {
                if (filteredInds[i] != NULL && inds[i]->getStats().age < min_age){ // if not already NULL + age too young
                    filteredInds[i] = NULL; // set it to NULL
}
            }
        }
        // check max age
        if (max_age!=-9){
            // loop over all number of individuals in cell
            for (int i = 0; i < ninds; i++) {
                if (filteredInds[i] != NULL && inds[i]->getStats().age > max_age){// if not already NULL + age too old; bugfix: check was inverted before 
                    if (filteredInds[i] != NULL) filteredInds[i] = NULL; // set it to NULL if not already NULL
                }
            }
        }
        // check stage
        if (stage!=-9){
            // loop over all number of individuals in cell
            for (int i = 0; i < ninds; i++) {
                if (filteredInds[i] != NULL && inds[i]->getStats().stage != stage){// if not already NULL + stage not correct
                    if (filteredInds[i] != NULL) filteredInds[i] = NULL; // set it to NULL if not already NULL
                }
            }
        }
        // check sex
        if (sex!=-9){
            // loop over all number of individuals in cell
            for (int i = 0; i < ninds; i++) {
                if (filteredInds[i] != NULL && inds[i]->getStats().sex != sex){// if not already NULL + sex not correct
                    if (filteredInds[i] != NULL) filteredInds[i] = NULL; // set it to NULL if not already NULL
                }
            }
        }
    } else {
#if RS_RCPP
        Rcpp::Rcout << "No individuals in source patch" << endl;
#endif
        return filteredInds;
        }
    int nfiltered = 0;
    for ( auto filtered : filteredInds){
        if (filtered != NULL) nfiltered++;
    }

    // loop over iterator of filteredInds and remove NULL values
    filteredInds.erase(std::remove(filteredInds.begin(), filteredInds.end(), nullptr), filteredInds.end());

    return filteredInds;
};
//---------------------------------------------------------------------------
// Clean the sampled individuals
//---------------------------------------------------------------------------
void Population::cleanSampledInds(Individual* pInd // Return a set of individuals with specified characteristics
){
    // find inds[j] and remove it from sampledInds
    sampledInds.erase(std::remove(sampledInds.begin(), sampledInds.end(), pInd), sampledInds.end());
};
//---------------------------------------------------------------------------
// Sample N individuals from the population with a given set of characteristics
// ---------------------------------------------------------------------------
int Population::sampleIndividuals( // Select a set of individuals with specified characteristics
// void Population::sampleIndividuals( // Select a set of individuals with specified characteristics
        int nb,	// number of individuals to sample
        int min_age,	// min age (0 if not set)
        int max_age,    // max age (max age if not set)
        int stage,    // stage
        int sex,     // sex
        bool includeTranslocated // also accept translocated individuals (status 10)
        ){
    if(sampledInds.size() > 0)  sampledInds.clear(); // clear old vector
    auto& rng = pRandom->getRNG(); // random number for sampling from suitable individuals bugfix: auto& - plain auto would copy the generator again

    // get individuals with the characteristics
    std::vector <Individual*> filtered;
    filtered = getIndsWithCharacteristics(min_age, max_age, stage, sex, includeTranslocated);
#if RS_RCPP
    Rcpp::Rcout << "Number of individuals with fitting characteristics: " << filtered.size() << endl;
#endif
    if (filtered.size() <= nb)
        // Sample all individuals in selected stages
        sampledInds = filtered;
    else {
        vector<Individual*> out;
        // Sample n individuals across filtered individuals
        std::sample(filtered.begin(), filtered.end(), std::back_inserter(out), nb, rng);
        std::copy(out.begin(), out.end(), std::inserter(sampledInds, sampledInds.end()));
    }

    int nb_sampled = 0;
    if (sampledInds.size() > 0) {
        for (int i = 0; i < (int)sampledInds.size(); i++) {
            if (sampledInds[i] != NULL) nb_sampled++;
        }
    }
    return nb_sampled;
}
// ---------------------------------------------------------------------------
// catch individuals according to catching rate
// ---------------------------------------------------------------------------
Individual* Population::catchIndividual( // Translocate a set of individuals with specified characteristics
        double catching_rate,
        int j
){
    Individual* catched;
    int id = inds[j]->getId();
    // If individual is part of the sampledInds vector:
    if (std::find(sampledInds.begin(), sampledInds.end(), inds[j]) != std::end(sampledInds)){
        // try to catch individual
#if RS_RCPP
        if(catching_rate > 1) Rcpp::Rcout << "Catching rate: " << catching_rate << std::endl;
#endif
        if (pRandom->Bernoulli(catching_rate)){
            indStats indstat = inds[j]->getStats();
            catched = inds[j];
            // remove individual from source patch
            inds[j] = 0;
            nInds[indstat.stage][indstat.sex]--;
            cleanSampledInds(catched); // clean vector of sampled individuals after the event
            return catched;
        }else {
            cleanSampledInds(inds[j]); // clean vector of sampled individuals after the event
            return NULL;
            }
    } else {
        return NULL;
    }
}

// ---------------------------------------------------------------------------
// harvest individuals according to harvest success
// ---------------------------------------------------------------------------
Individual* Population::harvestIndividual( // Harvest one sampled individual with probability harvest_success
        double harvest_success,
        int j
){
    // If individual is part of the sampledInds vector:
    if (std::find(sampledInds.begin(), sampledInds.end(), inds[j]) != std::end(sampledInds)){
        if (pRandom->Bernoulli(harvest_success)){
            indStats indstat = inds[j]->getStats();
            Individual* harvested = inds[j];
            // remove the individual from the living population immediately (as in catchIndividual()),
            // so that it no longer counts for density dependence, reproduction or population output;
            // it is kept in harvestedInds for the individual output (status 11) and deleted in survival1()
            inds[j] = nullptr;
            nInds[indstat.stage][indstat.sex]--;
            harvested->setStatus(11);
            harvestedInds.push_back(harvested);
            cleanSampledInds(harvested); // clean vector of sampled individuals after the event
            return harvested;
        } else {
            cleanSampledInds(inds[j]); // clean vector of sampled individuals after the event
            return NULL;
        }
    } else {
        return NULL;
    }
}

// ---------------------------------------------------------------------------
bool Population::getSizeSampledInds(
){
    bool size = false;
    if (sampledInds.size() > 0) size = true;
    return size;
};

//---------------------------------------------------------------------------


