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

#include "Population.h"
//---------------------------------------------------------------------------

ofstream outPop;
ofstream outInds;

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
	patchPopn pp = patchPopn();
	pp.pSp = (intptr)pSpecies; 
	pp.pPop = (intptr)this;
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
						while (rrr > ageProb[ageclass]) ageclass++;
						age = minage + ageclass;
					}
					else age = minage;
					break;
				}
			}
			else age = stg;

			inds.push_back(new Individual(pCell, pPatch, stg, age, sstruct.repInterval,
				probmale, trfr.usesMovtProc, trfr.moveType));

			sex = inds[nindivs + i]->getSex();
			if (pSpecies->getNTraits() > 0) {
				// individual variation - set up genetics
				inds[nindivs + i]->setUpGenes(pSpecies, resol);
			}
			nInds[stg][sex]++;
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
}

traitsums Population::getIndTraitsSums(Species* pSpecies) {
	int g;
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

	int ninds = (int)inds.size();
	for (int iInd = 0; iInd < ninds; iInd++) {
		int sex = inds[iInd]->getSex();
		if (emig.sexDep || trfr.sexDep || sett.sexDep) 
			g = sex; 
		else g = 0;
		ts.ninds[g] += 1;

		// emigration traits
		emigTraits e = inds[iInd]->getIndEmigTraits();
		if (emig.sexDep) g = sex; 
		else g = 0;
		ts.sumD0[g] += e.d0;    
		ts.ssqD0[g] += e.d0 * e.d0;
		ts.sumAlpha[g] += e.alpha; 
		ts.ssqAlpha[g] += e.alpha * e.alpha;
		ts.sumBeta[g] += e.beta;  
		ts.ssqBeta[g] += e.beta * e.beta;

		// transfer traits
		if (trfr.usesMovtProc) {

			switch (trfr.moveType) {

			case 1: // SMS
			{
				trfrSMSTraits sms = inds[iInd]->getIndSMSTraits();
				g = 0; // CURRENTLY INDIVIDUAL VARIATION CANNOT BE SEX-DEPENDENT
				ts.sumDP[g] += sms.dp; 
				ts.ssqDP[g] += sms.dp * sms.dp;
				ts.sumGB[g] += sms.gb;
				ts.ssqGB[g] += sms.gb * sms.gb;
				ts.sumAlphaDB[g] += sms.alphaDB;
				ts.ssqAlphaDB[g] += sms.alphaDB * sms.alphaDB;
				ts.sumBetaDB[g] += sms.betaDB; 
				ts.ssqBetaDB[g] += sms.betaDB * sms.betaDB;
				break;
			}
			case 2:
			{
				trfrCRWTraits c = inds[iInd]->getIndCRWTraits();
				g = 0; // CURRENTLY INDIVIDUAL VARIATION CANNOT BE SEX-DEPENDENT
				ts.sumStepL[g] += c.stepLength;
				ts.ssqStepL[g] += c.stepLength * c.stepLength;
				ts.sumRho[g] += c.rho;       
				ts.ssqRho[g] += c.rho * c.rho;
				break;
			}
			default:
				throw runtime_error("usesMoveProcess is ON but moveType is neither 1 (SMS) or 2 (CRW).");
				break;
			}
		}
		else {
			trfrKernelParams k = inds[iInd]->getIndKernTraits();
			if (trfr.sexDep) g = sex; 
			else g = 0;
			ts.sumDist1[g] += k.meanDist1; 
			ts.ssqDist1[g] += k.meanDist1 * k.meanDist1;
			ts.sumDist2[g] += k.meanDist2;
			ts.ssqDist2[g] += k.meanDist2 * k.meanDist2;
			ts.sumProp1[g] += k.probKern1; 
			ts.ssqProp1[g] += k.probKern1 * k.probKern1;
		}
		// settlement traits
		settleTraits s = inds[iInd]->getIndSettTraits();
		if (sett.sexDep) g = sex; 
		else g = 0;
		//	g = 0; // CURRENTLY INDIVIDUAL VARIATION CANNOT BE SEX-DEPENDENT
		ts.sumS0[g] += s.s0;     
		ts.ssqS0[g] += s.s0 * s.s0;
		ts.sumAlphaS[g] += s.alpha; 
		ts.ssqAlphaS[g] += s.alpha * s.alpha;
		ts.sumBetaS[g] += s.beta;   
		ts.ssqBetaS[g] += s.beta * s.beta;

		if (gMaxNbSexes > 1) g = sex; 
		else g = 0;

		ts.sumGeneticFitness[g] += inds[iInd]->getGeneticFitness();
		ts.ssqGeneticFitness[g] += inds[iInd]->getGeneticFitness() * inds[iInd]->getGeneticFitness();
	}
	return ts;
}

int Population::getNInds(void) { return (int)inds.size(); }

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

		const auto trait = individual->getTrait(NEUTRAL);
		int whichLocus = 0;
		for (auto position : positions) {

			int alleleOnChromA = (int)trait->getAlleleValueAtLocus(0, position);
			popNeutralCountTables[whichLocus].incrementTally(alleleOnChromA);

			if (ploidy == 2) { // second allele and heterozygosity
				int alleleOnChromB = (int)trait->getAlleleValueAtLocus(1, position);
				popNeutralCountTables[whichLocus].incrementTally(alleleOnChromB);

				bool isHetero = alleleOnChromA != alleleOnChromB;
				if (isHetero) {
					popNeutralCountTables[whichLocus].incrementHeteroTally(alleleOnChromA);
					popNeutralCountTables[whichLocus].incrementHeteroTally(alleleOnChromB);
				}
			}
			whichLocus++;
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
			int counter = 0;
			for (auto position : positions) {
				hetero[counter] += trait->isHeterozygoteAtLocus(position);
				counter++;
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

popStats Population::getStats(void)
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
					if (dem.repType == 2) fec = pSpecies->getFec(stg, sex);
					else fec = pSpecies->getFec(stg, 0);
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

int Population::totalPop(void) {
	int t = 0;
	for (int stg = 0; stg < nStages; stg++) {
		for (int sex = 0; sex < nSexes; sex++) {
			t += nInds[stg][sex];
		}
	}
	return t;
}

int Population::stagePop(int stg) {
	int t = 0;
	if (stg < 0 || stg >= nStages) return t;
	for (int sex = 0; sex < nSexes; sex++) {
		t += nInds[stg][sex];
	}
	return t;
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
void Population::reproduction(const float localK, const float envval, const int resol)
{

	// get population size at start of reproduction
	int ninds = (int)inds.size();
	if (ninds == 0) return;

	int nsexes, stage, sex, njuvs, nj, nmales, nfemales;
	Cell* pCell;
	indStats ind;
	double expected;
	bool skipbreeding;

	envStochParams env = paramsStoch->getStoch();
	demogrParams dem = pSpecies->getDemogrParams();
	stageParams sstruct = pSpecies->getStageParams();
	emigRules emig = pSpecies->getEmigRules();
	transferRules trfr = pSpecies->getTransferRules();
	settleType sett = pSpecies->getSettle();

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
					fec[stg][sex] = pSpecies->getFec(stg, 0);
				}
				else
					fec[stg][sex] = pSpecies->getFec(stg, sex);
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
						// NOTE: matrix entries represent effect of ROW on COLUMN 
						// AND males precede females
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
						effect = (float)totalPop();
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
					expected = fec[stage][0];
					if (expected <= 0.0) njuvs = 0;
					else njuvs = pRandom->Poisson(expected);
					nj = (int)juvs.size();
					pCell = pPatch->getRandomCell();
					for (int j = 0; j < njuvs; j++) {

						Individual* newJuv;
						newJuv = new Individual(pCell, pPatch, 0, 0, 0, dem.propMales, trfr.usesMovtProc, trfr.moveType);

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
							expected = fec[stage][0]; // breeds
						}
						else expected = 0.0; // fails to breed
						if (expected <= 0.0) njuvs = 0;
						else njuvs = pRandom->Poisson(expected);

						if (njuvs > 0)
						{
							nj = (int)juvs.size();
							// select father at random from breeding males ...
							int rrr = 0;
							if (nmales > 1) rrr = pRandom->IRandom(0, nmales - 1);
							father = fathers[rrr];
							pCell = pPatch->getRandomCell();
							for (int j = 0; j < njuvs; j++) {
								Individual* newJuv;

								newJuv = new Individual(pCell, pPatch, 0, 0, 0, dem.propMales, trfr.usesMovtProc, trfr.moveType);

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
		inds = juvs;
	}
	juvs.clear();

}

Individual* Population::sampleInd() const {
	int index = pRandom->IRandom(0, static_cast<int>(inds.size() - 1));
	return inds[index];
}

void Population::sampleIndsWithoutReplacement(string strNbToSample, const set<int>& sampleStages) {

	if (sampledInds.size() > 0) {
		sampledInds.clear();
	}
	auto rng = pRandom->getRNG();
	vector<Individual*> stagedInds;

	// Stage individuals in eligible stages
	for (int stage : sampleStages) {
		vector<Individual*> toAdd = getIndividualsInStage(stage);
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
	indStats ind;

	// to avoid division by zero, assume carrying capacity is at least one individual
	// localK can be zero if there is a moving gradient or stochasticity in K
	if (localK < 1.0) localK = 1.0;
	NK = static_cast<float>(totalPop()) / localK;

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
		ind = inds[i]->getStats();
		if (ind.status < 1) {
			if (emig.indVar) { // individual variability in emigration
				if (dem.stageStruct && ind.stage != emig.emigStage) {
					// emigration may not occur
					pbDisp = 0.0;
				}
				else { // non-structured or individual is in emigration stage
					eparams = inds[i]->getIndEmigTraits();
					if (emig.densDep) { // density-dependent
						NK = (float)totalPop() / localK;
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

// If an Individual has been identified as an emigrant, remove it from the Population
disperser Population::extractDisperser(int ix) {
	disperser d = disperser();
	indStats ind = inds[ix]->getStats();
	if (ind.status == 1) { // emigrant
		d.pInd = inds[ix]; 
		d.yes = true;
		inds[ix] = 0;
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
	pCell = inds[ix]->getLocn(1);
	d.pInd = inds[ix];  d.pCell = pCell; d.yes = false;
	if (ind.status == 4 || ind.status == 5) { // settled
		d.yes = true;
		inds[ix] = 0;
		nInds[ind.stage][ind.sex]--;
	}
	return d;
}

// Add a specified individual to the new/current dispersal group
void Population::recruit(Individual* pInd) {
	inds.push_back(pInd);
	indStats ind = pInd->getStats();
	nInds[ind.stage][ind.sex]++;
}

//---------------------------------------------------------------------------

// Transfer is run for populations in the matrix only
#if RS_RCPP // included also SEASONAL
int Population::transfer(Landscape* pLandscape, short landIx, short nextseason)
#else
int Population::transfer(Landscape* pLandscape, short landIx)
#endif
{
	int ndispersers = 0;
	int disperser;
	short othersex;
	bool mateOK, densdepOK;
	intptr patch, popn;
	int patchnum;
	double localK, popsize, settprob;
	Patch* pPatch = 0;
	Cell* pCell = 0;
	indStats ind;
	Population* pNewPopn = 0;
	locn newloc = locn();
	locn nbrloc = locn();

	landData ppLand = pLandscape->getLandData();
	short reptype = pSpecies->getRepType();
	transferRules trfr = pSpecies->getTransferRules();
	settleType settletype = pSpecies->getSettle();
	settleRules sett;
	settleTraits settDD;
	settlePatch settle;
	simParams sim = paramsSim->getSim();
	// each individual takes one step
	// for dispersal by kernel, this should be the only step taken
	int ninds = (int)inds.size();

	for (int i = 0; i < ninds; i++) {
		if (trfr.usesMovtProc) {

			disperser = inds[i]->moveStep(pLandscape, pSpecies, landIx, sim.absorbing);
		}
		else {
			disperser = inds[i]->moveKernel(pLandscape, pSpecies, sim.absorbing);
		}
		ndispersers += disperser;
		if (disperser) {
			if (reptype > 0)
			{ // sexual species - record as potential settler in new patch
				if (inds[i]->getStatus() == 2)
				{ // disperser has found a patch
					pCell = inds[i]->getLocn(1);
					patch = pCell->getPatch();
					if (patch != 0) { // not no-data area
						pPatch = (Patch*)patch;
						pPatch->incrPossSettler(pSpecies, inds[i]->getSex());
					}
				}
			}
		}
	}

	for (int i = 0; i < ninds; i++) {
		ind = inds[i]->getStats();
		if (ind.sex == 0) othersex = 1; else othersex = 0;
		if (settletype.stgDep) {
			if (settletype.sexDep) sett = pSpecies->getSettRules(ind.stage, ind.sex);
			else sett = pSpecies->getSettRules(ind.stage, 0);
		}
		else {
			if (settletype.sexDep) sett = pSpecies->getSettRules(0, ind.sex);
			else sett = pSpecies->getSettRules(0, 0);
		}
		if (ind.status == 2)
		{ // awaiting settlement
			pCell = inds[i]->getLocn(1);
			if (pCell == 0) {
				// this condition can occur in a patch-based model at the time of a dynamic landscape
				// change when there is a range restriction in place, since a patch can straddle the
				// range restriction and an individual forced to disperse upon patch removal could
				// start its trajectory beyond the boundary of the restrictyed range - such a model is 
				// not good practice, but the condition must be handled by killing the individual conceerned
				ind.status = 6;
			}
			else {
				mateOK = false;
				if (sett.findMate) {
					// determine whether at least one individual of the opposite sex is present in the
					// new population
					if (matePresent(pCell, othersex)) mateOK = true;
				}
				else { // no requirement to find a mate
					mateOK = true;
				}
				densdepOK = false;
				settle = inds[i]->getSettPatch();
				if (sett.densDep)
				{
					patch = pCell->getPatch();
					if (patch != 0) { // not no-data area
						pPatch = (Patch*)patch;
						if (settle.settleStatus == 0
							|| settle.pSettPatch != pPatch)
							// note: second condition allows for having moved from one patch to another
							// adjacent one
						{
							// determine whether settlement occurs in the (new) patch
							localK = (double)pPatch->getK();
							popn = pPatch->getPopn((intptr)pSpecies);

							if (popn == 0) { // population has not been set up in the new patch
								popsize = 0.0;
							}
							else {
								pNewPopn = (Population*)popn;
								popsize = (double)pNewPopn->totalPop();
							}
							if (localK > 0.0) {
								// make settlement decision
								if (settletype.indVar) settDD = inds[i]->getIndSettTraits();
								else settDD = pSpecies->getSpSettTraits(ind.stage, ind.sex);
								settprob = settDD.s0 /
									(1.0 + exp(-(popsize / localK - (double)settDD.beta) * (double)settDD.alpha));
								if (pRandom->Bernoulli(settprob)) { // settlement allowed
									densdepOK = true;
									settle.settleStatus = 2;
								}
								else { // settlement procluded
									settle.settleStatus = 1;
								}
								settle.pSettPatch = pPatch;
							}
							inds[i]->setSettPatch(settle);
						}
						else {
							if (settle.settleStatus == 2) { // previously allowed to settle
								densdepOK = true;
							}
						}
					}
				}
				else { // no density-dependent settlement
					densdepOK = true;
					settle.settleStatus = 2;
					settle.pSettPatch = pPatch;
					inds[i]->setSettPatch(settle);
				}

				if (mateOK && densdepOK) { // can recruit to patch
					ind.status = 4;
					ndispersers--;
				}
				else { // does not recruit
					if (trfr.usesMovtProc) {
						ind.status = 1; // continue dispersing, unless ...
						// ... maximum steps has been exceeded
						pathSteps steps = inds[i]->getSteps();
						settleSteps settsteps = pSpecies->getSteps(ind.stage, ind.sex);
						if (steps.year >= settsteps.maxStepsYr) {
							ind.status = 3; // waits until next year
						}
						if (steps.total >= settsteps.maxSteps) {
							ind.status = 6; // dies
						}
					}
					else { // dispersal kernel
						if (sett.wait) {
							ind.status = 3; // wait until next dispersal event
						}
						else {
							ind.status = 6; // (dies unless a neighbouring cell is suitable)
						}
						ndispersers--;
					}
				}
			}

			inds[i]->setStatus(ind.status);
		}
#if RS_RCPP
		// write each individuals current movement step and status to paths file
		if (trfr.usesMovtProc && sim.outPaths) {
			if (nextseason >= sim.outStartPaths && nextseason % sim.outIntPaths == 0) {
				inds[i]->outMovePath(nextseason);
			}
		}
#endif
		if (!trfr.usesMovtProc && sett.go2nbrLocn && (ind.status == 3 || ind.status == 6))
		{
			// for kernel-based transfer only ...
			// determine whether recruitment to a neighbouring cell is possible
			pCell = inds[i]->getLocn(1);
			newloc = pCell->getLocn();
			vector <Cell*> nbrlist;
			for (int dx = -1; dx < 2; dx++) {
				for (int dy = -1; dy < 2; dy++) {
					if (dx != 0 || dy != 0) { //cell is not the current cell
						nbrloc.x = newloc.x + dx; nbrloc.y = newloc.y + dy;
						if (nbrloc.x >= 0 && nbrloc.x <= ppLand.maxX
							&& nbrloc.y >= 0 && nbrloc.y <= ppLand.maxY) { // within landscape
							// add to list of potential neighbouring cells if suitable, etc.
							pCell = pLandscape->findCell(nbrloc.x, nbrloc.y);
							if (pCell != 0) { // not no-data area
								patch = pCell->getPatch();
								if (patch != 0) { // not no-data area
									pPatch = (Patch*)patch;
									patchnum = pPatch->getPatchNum();
									if (patchnum > 0 && pPatch != inds[i]->getNatalPatch())
									{ // not the matrix or natal patch
										if (pPatch->getK() > 0.0)
										{ // suitable
											if (sett.findMate) {
												if (matePresent(pCell, othersex)) nbrlist.push_back(pCell);
											}
											else
												nbrlist.push_back(pCell);
										}
									}
								}
							}
						}
					}
				}
			}
			int listsize = (int)nbrlist.size();
			if (listsize > 0) { // there is at least one suitable neighbouring cell
				if (listsize == 1) {
					inds[i]->moveto(nbrlist[0]);
				}
				else { // select at random from the list
					int rrr = pRandom->IRandom(0, listsize - 1);
					inds[i]->moveto(nbrlist[rrr]);
				}
			}
			// else list empty - do nothing - individual retains its current location and status
		}
	}
	return ndispersers;
}

// Determine whether there is a potential mate present in a patch which a potential
// settler has reached
bool Population::matePresent(Cell* pCell, short othersex)
{
	int patch;
	Patch* pPatch;
	Population* pNewPopn;
	int popsize = 0;
	bool matefound = false;

	patch = (int)pCell->getPatch();
	if (patch != 0) {
		pPatch = (Patch*)pCell->getPatch();
		if (pPatch->getPatchNum() > 0) { // not the matrix patch
			if (pPatch->getK() > 0.0)
			{ // suitable
				pNewPopn = (Population*)pPatch->getPopn((intptr)pSpecies);
				if (pNewPopn != 0) {
					// count members of other sex already resident in the patch
					for (int stg = 0; stg < nStages; stg++) {
						popsize += pNewPopn->nInds[stg][othersex];
					}
				}
				if (popsize < 1) {
					// add any potential settlers of the other sex
					popsize += pPatch->getPossSettlers(pSpecies, othersex);
				}
			}
		}
	}
	if (popsize > 0) matefound = true;
	return matefound;
}

//---------------------------------------------------------------------------
// Determine survival and development and record in individual's status code
// Changes are NOT applied to the Population at this stage

// FOR MULTIPLE SPECIES, MAY NEED TO SEPARATE OUT THIS IDENTIFICATION STAGE,
// SO THAT IT CAN BE PERFORMED FOR ALL SPECIES BEFORE ANY UPDATING OF POPULATIONS

void Population::survival0(float localK, short option0, short option1)
{
	// option0:	0 - stage 0 (juveniles) only
	//			1 - all stages
	//			2 - stage 1 and above (all non-juveniles)
	// 
	// option1:	0 - development only (when survival is annual)
	//	  	 	1 - development and survival
	//	  	 	2 - survival only (when survival is annual)
	densDepParams ddparams = pSpecies->getDensDep();
	demogrParams dem = pSpecies->getDemogrParams();
	stageParams sstruct = pSpecies->getStageParams();

	// get current population size
	int ninds = (int)inds.size();
	if (ninds == 0) return;
	// set up local copies of species development and survival tables
	int nsexes;
	if (dem.repType == 0) nsexes = 1; else nsexes = 2;
	float dev[gMaxNbStages][gMaxNbSexes];
	float surv[gMaxNbStages][gMaxNbSexes];
	short minAge[gMaxNbStages][gMaxNbSexes];
	for (int stg = 0; stg < sstruct.nStages; stg++) {
		for (int sex = 0; sex < nsexes; sex++) {
			if (dem.stageStruct) {
				if (dem.repType == 1) { // simple sexual model
					// both sexes use development and survival recorded for females
					dev[stg][sex] = pSpecies->getDev(stg, 0);
					surv[stg][sex] = pSpecies->getSurv(stg, 0);
					minAge[stg][sex] = pSpecies->getMinAge(stg, 0);
				}
				else {
					dev[stg][sex] = pSpecies->getDev(stg, sex);
					surv[stg][sex] = pSpecies->getSurv(stg, sex);
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
						// NOTE: matrix entries represent effect of ROW on COLUMN 
						// AND males precede females
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
						effect = (float)totalPop();
					if (localK > 0.0)
						dev[stg][sex] *= exp(-(ddparams.devCoeff * effect) / localK);
				} // end of if (sstruct.devDens && stg > 0)
				if (option1 != 0 && sstruct.survDens) {

					float effect = 0.0;
					if (sstruct.survStageDens) { // stage-specific density dependence
						// NOTE: matrix entries represent effect of ROW on COLUMN 
						// AND males precede females
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
						effect = (float)totalPop();
					if (localK > 0.0)
						surv[stg][sex] *= exp(-(ddparams.survCoeff * effect) / localK);
				} // end of if (sstruct.survDens)
			}
		}
	}
	// identify which individuals die or develop
	for (int i = 0; i < ninds; i++) {
		indStats ind = inds[i]->getStats();

		if ((ind.stage == 0 && option0 < 2) || (ind.stage > 0 && option0 > 0)) {
			// condition for processing the stage is met...
			if (ind.status < 6) { // not already doomed
				if (ind.sex < sex_t::FEM || ind.sex > sex_t::MAL)
					// ?? MSVC believes it's important to bound check ind.sex
					throw runtime_error("Individual sex is out of bounds");
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
	int ninds = (int)inds.size();

	for (int i = 0; i < ninds; i++) {
		indStats ind = inds[i]->getStats();

		if (ind.status > 5) { // doomed to die
			if (ind.status != 10) //not going into cold storage
				delete inds[i];
			inds[i] = NULL;
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
		inds[i]->ageIncrement(sstruct.maxAge);
	}
}

//---------------------------------------------------------------------------
// Remove zero pointers to dead or dispersed individuals
void Population::clean(void)
{
	int ninds = (int)inds.size();
	if (ninds > 0) {
			// ALTERNATIVE METHOD: AVOIDS SLOW SORTING OF POPULATION
		std::vector <Individual*> survivors; // all surviving individuals
		for (int i = 0; i < ninds; i++) {
			if (inds[i] != NULL) {
				survivors.push_back(inds[i]);
			}
		}
		inds.clear();
		inds = survivors;
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
// Open population file and write header record
bool Population::outPopHeaders(int landNr, bool patchModel) {

	if (landNr == -999) { // close file
		if (outPop.is_open()) outPop.close();
		outPop.clear();
		return true;
	}
	string name;
	simParams sim = paramsSim->getSim();
	envGradParams grad = paramsGrad->getGradient();

	// NEED TO REPLACE CONDITIONAL COLUMNS BASED ON ATTRIBUTES OF ONE SPECIES TO COVER
	// ATTRIBUTES OF *ALL* SPECIES AS DETECTED AT MODEL LEVEL
	demogrParams dem = pSpecies->getDemogrParams();
	stageParams sstruct = pSpecies->getStageParams();
	name = paramsSim->getDir(2)
		+ (sim.batchMode ? "Batch" + to_string(sim.batchNum) + "_" : "")
		+ "Batch" + to_string(sim.batchNum) + "_"
		+ "Sim" + to_string(sim.simulation) + "_Land" + to_string(landNr) + "_Pop.txt";

	outPop.open(name.c_str());
	outPop << "Rep\tYear\tRepSeason";
	if (patchModel) outPop << "\tPatchID\tNcells";
	else outPop << "\tx\ty";
	// determine whether environmental data need be written for populations
	bool writeEnv = false;
	if (grad.gradient) writeEnv = true;
	if (paramsStoch->envStoch()) writeEnv = true;
	if (writeEnv) outPop << "\tEpsilon\tGradient\tLocal_K";
	outPop << "\tSpecies\tNInd";
	if (dem.stageStruct) {
		if (dem.repType == 0)
		{
			for (int i = 1; i < sstruct.nStages; i++) outPop << "\tNInd_stage" << i;
			outPop << "\tNJuvs";
		}
		else {
			for (int i = 1; i < sstruct.nStages; i++)
				outPop << "\tNfemales_stage" << i << "\tNmales_stage" << i;
			outPop << "\tNJuvFemales\tNJuvMales";
		}
	}
	else {
		if (dem.repType != 0) outPop << "\tNfemales\tNmales";
	}
	outPop << endl;
	return outPop.is_open();
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
	outPop << rep << "\t" << yr << "\t" << gen;
	if (patchModel) {
		outPop << "\t" << pPatch->getPatchNum();
		outPop << "\t" << pPatch->getNCells();
	}
	else {
		locn loc = pPatch->getCellLocn(0);
		outPop << "\t" << loc.x << "\t" << loc.y;
	}
	if (writeEnv) {
		if (pPatch->getPatchNum() == 0) { // matrix
			outPop << "\t0\t0\t0";
		}
		else {
			float k = pPatch->getK();
			float envval = 0.0;
			pCell = pPatch->getRandomCell();
			if (pCell != 0) envval = pCell->getEnvVal();
			outPop << "\t" << eps << "\t" << envval << "\t" << k;
		}
	}
	outPop << "\t" << pSpecies->getSpNum();
	if (dem.stageStruct) {
		p = getStats();
		outPop << "\t" << p.nNonJuvs;
		// non-juvenile stage totals from permanent array
		for (int stg = 1; stg < nStages; stg++) {
			for (int sex = 0; sex < nSexes; sex++) {
				outPop << "\t" << nInds[stg][sex];
			}
		}
		// juveniles from permanent array
		for (int sex = 0; sex < nSexes; sex++) {
			outPop << "\t" << nInds[0][sex];
		}
	}
	else { // non-structured population
		outPop << "\t" << totalPop();
		if (dem.repType != 0)
		{ // sexual model
			outPop << "\t" << nInds[1][0] << "\t" << nInds[1][1];
		}
	}
	outPop << endl;
}

//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// Open individuals file and write header record
void Population::outIndsHeaders(int rep, int landNr, bool patchModel)
{
	if (landNr == -999) { // close file
		if (outInds.is_open()) {
			outInds.close(); outInds.clear();
		}
		return;
	}

	string name;
	demogrParams dem = pSpecies->getDemogrParams();
	emigRules emig = pSpecies->getEmigRules();
	transferRules trfr = pSpecies->getTransferRules();
	settleType sett = pSpecies->getSettle();
	simParams sim = paramsSim->getSim();

	name = paramsSim->getDir(2)
		+ (sim.batchMode ? "Batch" + to_string(sim.batchNum) + "_" : "")
		+ "Sim" + to_string(sim.simulation)
		+ "_Land" + to_string(landNr) + "_Rep" + to_string(rep) + "_Inds.txt";

	outInds.open(name.c_str());
	outInds << "Rep\tYear\tRepSeason\tSpecies\tIndID\tStatus";
	if (patchModel) outInds << "\tNatal_patch\tPatchID";
	else outInds << "\tNatal_X\tNatal_Y\tX\tY";
	if (dem.repType != 0) outInds << "\tSex";
	if (dem.stageStruct) outInds << "\tAge\tStage";
	if (pSpecies->getNbGenLoadTraits() > 0) outInds << "\tProbViable";
	if (emig.indVar) {
		if (emig.densDep) outInds << "\tD0\tAlpha\tBeta";
		else outInds << "\tEP";
	}
	if (trfr.indVar) {
		if (trfr.usesMovtProc) {
			if (trfr.moveType == 1) { // SMS
				outInds << "\tDP\tGB\tAlphaDB\tBetaDB";
			}
			if (trfr.moveType == 2) { // CRW
				outInds << "\tStepLength\tRho";
			}
		}
		else { // kernel
			outInds << "\tMeanDistI";
			if (trfr.twinKern) outInds << "\tMeanDistII\tPKernelI";
		}
	}
	if (sett.indVar) {
		outInds << "\tS0\tAlphaS\tBetaS";
	}
	outInds << "\tDistMoved";
#ifndef NDEBUG
	outInds << "\tNsteps";
#else
	if (trfr.usesMovtProc) outInds << "\tNsteps";
#endif
	outInds << endl;
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

	int ninds = (int)inds.size();
	for (int i = 0; i < ninds; i++) {
		indStats ind = inds[i]->getStats();
		if (yr == -1) { // write all initialised individuals
			writeInd = true;
			outInds << rep << "\t" << yr << "\t" << dem.repSeasons - 1;
		}
		else {
			if (dem.stageStruct && gen < 0) { // write status 9 individuals only
				if (ind.status == 9) {
					writeInd = true;
					outInds << rep << "\t" << yr << "\t" << dem.repSeasons - 1;
				}
				else writeInd = false;
			}
			else {
				writeInd = true;
				outInds << rep << "\t" << yr << "\t" << gen;
			}
		}
		if (writeInd) {
			outInds << "\t" << spNum << "\t" << inds[i]->getId();
			if (dem.stageStruct) outInds << "\t" << ind.status;
			else { // non-structured population
				outInds << "\t" << ind.status;
			}
			pCell = inds[i]->getLocn(1);
			locn loc = locn();
			if (pCell == 0) loc.x = loc.y = -1; // beyond boundary or in no-data cell
			else loc = pCell->getLocn();
			pCell = inds[i]->getLocn(0);
			locn natalloc = pCell->getLocn();
			if (ppLand.patchModel) {
				outInds << "\t" << inds[i]->getNatalPatch()->getPatchNum();
				if (loc.x == -1) outInds << "\t-1";
				else outInds << "\t" << patchNum;
			}
			else { // cell-based model
				outInds << "\t" << (float)natalloc.x << "\t" << natalloc.y;
				outInds << "\t" << (float)loc.x << "\t" << (float)loc.y;
			}
			if (dem.repType != 0) outInds << "\t" << ind.sex;
			if (dem.stageStruct) outInds << "\t" << ind.age << "\t" << ind.stage;

			if (pSpecies->getNbGenLoadTraits() > 0) outInds << "\t" << inds[i]->getGeneticFitness();
		
			if (emig.indVar) {
				emigTraits e = inds[i]->getIndEmigTraits();
				if (emig.densDep) {
					outInds << "\t" << e.d0 << "\t" << e.alpha << "\t" << e.beta;
				}
				else {
					outInds << "\t" << e.d0;
				}
			} // end of if (emig.indVar)
			if (trfr.indVar) {
				if (trfr.usesMovtProc) {
					if (trfr.moveType == 1) { // SMS
						trfrSMSTraits s = inds[i]->getIndSMSTraits();
						outInds << "\t" << s.dp << "\t" << s.gb;
						outInds << "\t" << s.alphaDB << "\t" << s.betaDB;
					} // end of SMS
					if (trfr.moveType == 2) { // CRW
						trfrCRWTraits c = inds[i]->getIndCRWTraits();
						outInds << "\t" << c.stepLength << "\t" << c.rho;
					} // end of CRW
				}
				else { // kernel
					trfrKernelParams k = inds[i]->getIndKernTraits();
					if (trfr.twinKern)
					{
						outInds << "\t" << k.meanDist1 << "\t" << k.meanDist2 << "\t" << k.probKern1;
					}
					else {
						outInds << "\t" << k.meanDist1;
					}
				}
			}

			if (sett.indVar) {
				settleTraits s = inds[i]->getIndSettTraits();
				outInds << "\t" << s.s0 << "\t" << s.alpha << "\t" << s.beta;
			}

			// distance moved (metres)
			if (loc.x == -1) outInds << "\t-1";
			else {
				float d = ppLand.resol * sqrt((float)((natalloc.x - loc.x) * (natalloc.x - loc.x)
					+ (natalloc.y - loc.y) * (natalloc.y - loc.y)));
				outInds << "\t" << d;
			}
#ifndef NDEBUG
			// ALWAYS WRITE NO. OF STEPS
			steps = inds[i]->getSteps();
			outInds << "\t" << steps.year;
#else
			if (trfr.usesMovtProc) {
				steps = inds[i]->getSteps();
				outInds << "\t" << steps.year;
			}
#endif
			outInds << endl;
		} // end of writeInd condition

	}
}

void Population::outputGeneValues(ofstream& ofsGenes, const int& yr, const int& gen) const {
	
	const bool isDiploid = pSpecies->isDiploid();
	int indID;
	float alleleOnChromA, alleleOnChromB;
	float domCoefA, domCoefB;

	// Subset traits that are selected to be output
	set<TraitType> traitTypes = pSpecies->getTraitTypes();
	set<TraitType> outputTraitTypes;
	for (auto trType : traitTypes) {
		if (pSpecies->getSpTrait(trType)->isOutput())
			outputTraitTypes.insert(trType);
	}

	// Fetch map to positions for each trait
	// Presumably faster than fetching for every individual
	map<TraitType, set<int>> allGenePositions;
	for (auto trType : outputTraitTypes) {
		set<int> traitPositions = pSpecies->getSpTrait(trType)->getGenePositions();
		allGenePositions.insert(make_pair(trType, traitPositions));
	}

	set<int> positions;
	for (Individual* ind : sampledInds) {
		indID = ind->getId();
		for (auto trType : outputTraitTypes) {
			positions = allGenePositions[trType];
			auto indTrait = ind->getTrait(trType);
			for (auto pos : positions) {
				alleleOnChromA = indTrait->getAlleleValueAtLocus(0, pos);
				if (trType == GENETIC_LOAD1 || trType == GENETIC_LOAD2 || trType == GENETIC_LOAD3 || trType == GENETIC_LOAD4 || trType == GENETIC_LOAD5) {
					domCoefA = indTrait->getDomCoefAtLocus(0, pos);
				}
				else {
					domCoefA = 0.0;
				}
				ofsGenes << yr << '\t' << gen << '\t' << indID << '\t' << to_string(trType) << '\t' << pos << '\t' << alleleOnChromA << '\t' << domCoefA;
				if (isDiploid) {
					alleleOnChromB = indTrait->getAlleleValueAtLocus(1, pos);
					if (trType == GENETIC_LOAD1 || trType == GENETIC_LOAD2 || trType == GENETIC_LOAD3 || trType == GENETIC_LOAD4 || trType == GENETIC_LOAD5) {
						domCoefB = indTrait->getDomCoefAtLocus(1, pos);
					}
					else {
						domCoefB = 0.0;
					}
					ofsGenes << '\t' << alleleOnChromB << '\t' << domCoefB;
				}
				ofsGenes << endl;
			}
		}
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
