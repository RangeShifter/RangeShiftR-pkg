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

#include "Individual.h"
//---------------------------------------------------------------------------

int Individual::indCounter = 0;
TraitFactory Individual::traitFactory = TraitFactory();

//---------------------------------------------------------------------------

// Individual constructor
Individual::Individual(Cell* pCell, Patch* pPatch, short stg, short a, short repInt,
	float probmale, bool movt, short moveType)
{
	indId = indCounter; indCounter++; // unique identifier for each individual
	geneticFitness = 1.0;
	stage = stg;
	if (probmale <= 0.0) sex = FEM;
	else sex = pRandom->Bernoulli(probmale) ? MAL : FEM;
	age = a;
	status = 0;

	if (sex == 0 && repInt > 0) { // set no. of fallow seasons for female
		fallow = pRandom->IRandom(0, repInt);
	}
	else fallow = 9999;
	isDeveloping = false;
	pPrevCell = pCurrCell = pCell;
	pNatalPatch = pPatch;
	pTrfrData = nullptr; //set to null as default
	if (movt) {
		locn loc = pCell->getLocn();
		path = new pathData;
		path->year = 0; path->total = 0; path->out = 0;
		path->pSettPatch = 0; path->settleStatus = 0;
		if (moveType == 1) { // SMS
			// set up location data for SMS
			pTrfrData = make_unique<smsData>(loc, loc);

		}
		if (moveType == 2) { // CRW
			// set up continuous co-ordinates etc. for CRW movement
			float xc = ((float)pRandom->Random() * 0.999f) + (float)loc.x;
			float yc = (float)(pRandom->Random() * 0.999f) + (float)loc.y;
			float prevdrn = (float)(pRandom->Random() * 2.0 * PI);
			pTrfrData = make_unique<crwData>(prevdrn, xc, yc);
		}
	}
	else {
		path = 0;
		pTrfrData = make_unique<kernelData>(0.0, 0.0, 0.0);
	}
}

Individual::~Individual(void) {
	if (path != 0) delete path;
}


Individual* Individual::traitClone(Cell* pCell, Patch* pPatch, float probmale, bool movt, short moveType) {

	Individual* myTraitClone = new Individual(pCell, pPatch, 0, 0, 0, probmale, movt, moveType);
	myTraitClone->pEmigTraits = make_unique<emigTraits>(*pEmigTraits);
	myTraitClone->pSettleTraits = make_unique<settleTraits>(*pSettleTraits);
	myTraitClone->pTrfrData->clone(*pTrfrData);

	return myTraitClone;
}

void Individual::setEmigTraits(const emigTraits& emig) {
	pEmigTraits = make_unique<emigTraits>(emig);
}

void Individual::setSettleTraits(const settleTraits& settle) {
	pSettleTraits = make_unique<settleTraits>(settle);
}

QuantitativeTrait* Individual::getTrait(TraitType trait) const {
	auto p = this->spTraitTable.find(trait);
	if (p == spTraitTable.end())
		throw runtime_error("Trait does not exist in trait table.");
	else return p->second.get();
}

set<TraitType> Individual::getTraitTypes() {
	auto kv = std::views::keys(this->spTraitTable);
	set< TraitType > keys{ kv.begin(), kv.end() };
	return keys;
}

//---------------------------------------------------------------------------
// Inheritance for diploid, sexual species
//---------------------------------------------------------------------------
void Individual::inherit(Species* pSpecies, const Individual* mother, const Individual* father) {

	int events = 0;
	const set<int> chromosomeEnds = pSpecies->getChromosomeEnds();
	const int genomeSize = pSpecies->getGenomeSize();

	int maternalStartingChromosome = pRandom->Bernoulli(0.5);
	int paternalStartingChromosome = pRandom->Bernoulli(0.5);

	set<unsigned int> maternalRecomPositions;
	set<unsigned int> paternalRecomPositions;

	// Determine which parental chromosomes are inherited
	for (int pos : chromosomeEnds) {
		if (pRandom->Bernoulli(0.5)) // switch strand for next chromosome
			maternalRecomPositions.insert(pos);
		if (pRandom->Bernoulli(0.5))
			paternalRecomPositions.insert(pos);
	}

	// Draw recombination events for maternal genome
	if (pSpecies->getRecombinationRate() > 0.0)
		events = pRandom->Binomial(genomeSize, pSpecies->getRecombinationRate());
	// if poisson exceeds genomeSize, bound to genomeSize
	int nbrCrossOvers = events + maternalRecomPositions.size();
	if (nbrCrossOvers > genomeSize) {
		nbrCrossOvers = genomeSize;
	}
	while (maternalRecomPositions.size() < nbrCrossOvers) {
		// Sample recombination sites
		maternalRecomPositions.insert(pRandom->IRandom(0, genomeSize));
	}

	// Draw recombination events for paternal genome
	if (pSpecies->getRecombinationRate() > 0.0)
		events = pRandom->Binomial(genomeSize, pSpecies->getRecombinationRate());
	nbrCrossOvers = events + paternalRecomPositions.size();
	if (nbrCrossOvers > genomeSize) {
		nbrCrossOvers = genomeSize;
	}
	while (paternalRecomPositions.size() < nbrCrossOvers) {
		paternalRecomPositions.insert(pRandom->IRandom(0, genomeSize));
	}

	// Inherit genes for each trait
	const auto& spTraits = pSpecies->getTraitTypes();
	for (auto const& trait : spTraits)
	{
		const auto motherTrait = mother->getTrait(trait);
		const auto fatherTrait = father->getTrait(trait);
		auto newTrait = motherTrait->clone(); // shallow copy pointer to species-level attributes

		// Inherit from mother first
		newTrait->inheritGenes(true, motherTrait, maternalRecomPositions, maternalStartingChromosome);
		if (newTrait->isInherited()) {
			// Inherit father trait values
			newTrait->inheritGenes(false, fatherTrait, paternalRecomPositions, paternalStartingChromosome);
			if (newTrait->getMutationRate() > 0 && pSpecies->areMutationsOn())
				newTrait->mutate();
		}
		if (trait == GENETIC_LOAD1 || trait == GENETIC_LOAD2 || trait == GENETIC_LOAD3 || trait == GENETIC_LOAD4 || trait == GENETIC_LOAD5)
			geneticFitness *= newTrait->express();

		// Add the inherited trait and genes to the newborn's list
		spTraitTable.insert(make_pair(trait, move(newTrait)));
	}
}

//---------------------------------------------------------------------------
// Inheritance for haploid, asexual species
//---------------------------------------------------------------------------
void Individual::inherit(Species* pSpecies, const Individual* mother) {
	set<unsigned int> recomPositions; //not used here cos haploid but need it for inherit function, not ideal 
	int startingChromosome = 0;

	const auto& motherTraits = getTraitTypes();

	for (auto const& trait : motherTraits)
	{
		const auto motherTrait = mother->getTrait(trait);
		auto newTrait = motherTrait->clone(); // shallow copy, pointer to species trait initialised and empty sequence

		newTrait->inheritGenes(true, motherTrait, recomPositions, startingChromosome);
		if (newTrait->isInherited()) {
			if (newTrait->getMutationRate() > 0 && pSpecies->areMutationsOn())
				newTrait->mutate();
		}
		if (trait == GENETIC_LOAD1 || trait == GENETIC_LOAD2 || trait == GENETIC_LOAD3 || trait == GENETIC_LOAD4 || trait == GENETIC_LOAD5)
			geneticFitness *= newTrait->express();

		// Add the inherited trait and genes to the newborn's list
		spTraitTable.insert(make_pair(trait, move(newTrait)));
	}
}

// Initialise individual trait genes from species-level traits
void Individual::setUpGenes(Species* pSpecies, int resol) {

	// this way to keep spp trait table immutable i.e. not able to call getTraitTable, 
	// could pass it back by value (copy) instead but could be heavy if large map
	const auto& traitTypes = pSpecies->getTraitTypes();
	for (auto const& traitType : traitTypes)
	{
		const auto spTrait = pSpecies->getSpTrait(traitType);
		this->spTraitTable.emplace(traitType, traitFactory.Create(traitType, spTrait));
	}
	setDispersalPhenotypes(pSpecies, resol);
}

void Individual::setDispersalPhenotypes(Species* pSpecies, int resol) {

	const emigRules emig = pSpecies->getEmigRules();
	const transferRules trfr = pSpecies->getTransferRules();
	const settleType sett = pSpecies->getSettle();

	// record phenotypic traits
	if (emig.indVar)
		this->setEmigTraits(pSpecies, emig.sexDep, emig.densDep);
	if (trfr.indVar)
		this->setTransferTraits(pSpecies, trfr, resol);
	if (sett.indVar)
		this->setSettlementTraits(pSpecies, sett.sexDep);
}

void Individual::setTransferTraits(Species* pSpecies, transferRules trfr, int resol) {
	if (trfr.usesMovtProc) {
		if (trfr.moveType == 1) {
			setIndSMSTraits(pSpecies);
		}
		else
			setIndCRWTraits(pSpecies);
	}
	else
		setIndKernelTraits(pSpecies, trfr.sexDep, trfr.twinKern, resol);
}

void Individual::setSettlementTraits(Species* pSpecies, bool sexDep) {

	settleTraits s; s.s0 = s.alpha = s.beta = 0.0;
	if (sexDep && this->getSex() == MAL) {
		s.s0 = getTrait(S_S0_M)->express();
		s.alpha = getTrait(S_ALPHA_M)->express();
		s.beta = getTrait(S_BETA_M)->express();
	}
	else {
		s.s0 = getTrait(S_S0_F)->express();
		s.alpha = getTrait(S_ALPHA_F)->express();
		s.beta = getTrait(S_BETA_F)->express();
	}

	pSettleTraits = make_unique<settleTraits>();
	pSettleTraits->s0 = (float)(s.s0);
	pSettleTraits->alpha = (float)(s.alpha);
	pSettleTraits->beta = (float)(s.beta);
	if (pSettleTraits->s0 < 0.0) pSettleTraits->s0 = 0.0;
	if (pSettleTraits->s0 > 1.0) pSettleTraits->s0 = 1.0;
	return;
}


// Inherit genome from parent(s), diploid
void Individual::inheritTraits(Species* pSpecies, Individual* mother, Individual* father, int resol)
{
	inherit(pSpecies, mother, father);
	setDispersalPhenotypes(pSpecies, resol);
}

// Inherit genome from mother, haploid
void Individual::inheritTraits(Species* pSpecies, Individual* mother, int resol)
{
	inherit(pSpecies, mother);
	setDispersalPhenotypes(pSpecies, resol);
}

//---------------------------------------------------------------------------

// Identify whether an individual is a potentially breeding female -
// if so, return her stage, otherwise return 0
int Individual::breedingFem(void) {
	if (sex == FEM) {
		if (status == 0 || status == 4 || status == 5) return stage;
		else return 0;
	}
	else return 0;
}

int Individual::getId(void) { return indId; }

sex_t Individual::getSex(void) { return sex; }

int Individual::getStatus(void) { return status; }

float Individual::getGeneticFitness(void) { return geneticFitness; }

bool Individual::isViable() const {
	float probViability = geneticFitness > 1.0 ? 1.0 : geneticFitness;
	return probViability >= pRandom->Random();
}

indStats Individual::getStats(void) {
	indStats s;
	s.stage = stage; s.sex = sex; s.age = age; s.status = status; s.fallow = fallow;
	s.isDeveloping = isDeveloping;
	return s;
}

Cell* Individual::getLocn(const short option) {
	if (option == 0) { // return previous location
		return pPrevCell;
	}
	else { // return current location
		return pCurrCell;
	}
}

Patch* Individual::getNatalPatch(void) { return pNatalPatch; }

void Individual::setYearSteps(int t) {
	if (path != 0 && t >= 0) {
		if (t >= 0) path->year = t;
		else path->year = 666;
	}
}

pathSteps Individual::getSteps(void) {
	pathSteps s;
	if (path == 0) {
		s.year = 0; s.total = 0; s.out = 0;
	}
	else {
		s.year = path->year; s.total = path->total; s.out = path->out;
	}
	return s;
}

settlePatch Individual::getSettPatch(void) {
	settlePatch s;
	if (path == 0) {
		s.pSettPatch = 0; s.settleStatus = 0;
	}
	else {
		s.pSettPatch = path->pSettPatch; s.settleStatus = path->settleStatus;
	}
	return s;
}

void Individual::setSettPatch(const settlePatch s) {
	if (path == 0) {
		path = new pathData;
		path->year = 0; path->total = 0; path->out = 0; path->settleStatus = 0;
#if RS_RCPP
		path->pathoutput = 1;
#endif
	}
	if (s.settleStatus >= 0 && s.settleStatus <= 2) path->settleStatus = s.settleStatus;
	path->pSettPatch = s.pSettPatch;
}

void Individual::setEmigTraits(Species* pSpecies, bool sexDep, bool densityDep) {
	emigTraits e; e.d0 = e.alpha = e.beta = 0.0;
	if (sexDep) {
		if (this->getSex() == MAL) {
			e.d0 = this->getTrait(E_D0_M)->express();
			if (densityDep) {
				e.alpha = getTrait(E_ALPHA_M)->express();
				e.beta = getTrait(E_BETA_M)->express();
			}
		}
		else if (this->getSex() == FEM) {
			e.d0 = this->getTrait(E_D0_F)->express();
			if (densityDep) {
				e.alpha = getTrait(E_ALPHA_F)->express();
				e.beta = getTrait(E_BETA_F)->express();
			}
		}
		else {
			throw runtime_error("Attempt to express invalid emigration trait sex.");
		}
	}	
	else {
		e.d0 = this->getTrait(E_D0)->express();
		if (densityDep) {
			e.alpha = getTrait(E_ALPHA)->express();
			e.beta = getTrait(E_BETA)->express();
		}
	}

	pEmigTraits = make_unique<emigTraits>();
	pEmigTraits->d0 = e.d0;
	pEmigTraits->alpha = e.alpha;
	pEmigTraits->beta = e.beta;

	// Below must never trigger, phenotype is bounded in express()
	if (pEmigTraits->d0 < 0.0) throw runtime_error("d0 value has become negative.");
	if (pEmigTraits->d0 > 1.0) throw runtime_error("d0 value has exceeded 1.");
	return;
}

// Get phenotypic emigration traits
emigTraits Individual::getIndEmigTraits(void) {
	emigTraits e; 
	e.d0 = e.alpha = e.beta = 0.0;
	if (pEmigTraits != 0) {
		e.d0 = pEmigTraits->d0;
		e.alpha = pEmigTraits->alpha;
		e.beta = pEmigTraits->beta;
	}
	return e;
}
// Set phenotypic transfer by kernel traits
void Individual::setIndKernelTraits(Species* pSpecies, bool sexDep, bool twinKernel, int resol) {

	trfrKernelParams k; 
	k.meanDist1 = k.meanDist2 = k.probKern1 = 0.0;

	if (sexDep) {
		if (this->sex == MAL) {
			k.meanDist1 = getTrait(KERNEL_MEANDIST_1_M)->express();

			if (twinKernel) { // twin kernel
				k.meanDist2 = getTrait(KERNEL_MEANDIST_2_M)->express();
				k.probKern1 = getTrait(KERNEL_PROBABILITY_M)->express();
			}
		}
		else if (this->sex == FEM) {
			k.meanDist1 = getTrait(KERNEL_MEANDIST_1_F)->express();

			if (twinKernel) { // twin kernel
				k.meanDist2 = getTrait(KERNEL_MEANDIST_2_F)->express();
				k.probKern1 = getTrait(KERNEL_PROBABILITY_F)->express();
			}
		}
		else {
			throw runtime_error("Attempt to express invalid kernel transfer trait sex.");
		}
	}
	else {
		k.meanDist1 = getTrait(KERNEL_MEANDIST_1)->express();

		if (twinKernel) { // twin kernel
			k.meanDist2 = getTrait(KERNEL_MEANDIST_2)->express();
			k.probKern1 = getTrait(KERNEL_PROBABILITY)->express();
		}
	}
	
	float meanDist1 = (float)(k.meanDist1);
	float meanDist2 = (float)(k.meanDist2);
	float probKern1 = (float)(k.probKern1);

	if (!pSpecies->useFullKernel()) {
		// kernel mean(s) may not be less than landscape resolution
		if (meanDist1 < resol) meanDist1 = (float)resol;
		if (meanDist2 < resol) meanDist2 = (float)resol;
	}
	if (probKern1 < 0.0) probKern1 = 0.0;
	if (probKern1 > 1.0) probKern1 = 1.0;
	auto& pKernel = dynamic_cast<kernelData&>(*pTrfrData);
	pKernel.meanDist1 = meanDist1;
	pKernel.meanDist2 = meanDist2;
	pKernel.probKern1 = probKern1;

	return;
}



// Get phenotypic emigration traits
trfrKernelParams Individual::getIndKernTraits(void) {
	trfrKernelParams k; k.meanDist1 = k.meanDist2 = k.probKern1 = 0.0;
	if (pTrfrData != 0) {

		auto& pKernel = dynamic_cast<const kernelData&>(*pTrfrData);

		k.meanDist1 = pKernel.meanDist1;
		k.meanDist2 = pKernel.meanDist2;
		k.probKern1 = pKernel.probKern1;
	}

	return k;
}

void Individual::setIndSMSTraits(Species* pSpecies) {

	trfrSMSTraits s = pSpecies->getSpSMSTraits();

	double dp, gb, alphaDB, betaDB;
	dp = gb = alphaDB = betaDB = 0.0;
	dp = getTrait(SMS_DP)->express();
	gb = getTrait(SMS_GB)->express();
	if (s.goalType == 2) {
		alphaDB = getTrait(SMS_ALPHADB)->express();
		betaDB = getTrait(SMS_BETADB)->express();
	}

	auto& pSMS = dynamic_cast<smsData&>(*pTrfrData);
	pSMS.dp = (float)(dp);
	pSMS.gb = (float)(gb);
	if (s.goalType == 2) {
		pSMS.alphaDB = (float)(alphaDB);
		pSMS.betaDB = (int)(betaDB);
	}
	else {
		pSMS.alphaDB = s.alphaDB;
		pSMS.betaDB = s.betaDB;
	}
	if (pSMS.dp < 1.0) pSMS.dp = 1.0;
	if (pSMS.gb < 1.0) pSMS.gb = 1.0;
	if (pSMS.alphaDB <= 0.0) pSMS.alphaDB = 0.000001f;
	if (pSMS.betaDB < 1) pSMS.betaDB = 1;
	return;
}

trfrData* Individual::getTrfrData(void) {
	return pTrfrData.get();
}

// Get phenotypic transfer by SMS traits
trfrSMSTraits Individual::getIndSMSTraits(void) {

	trfrSMSTraits s; s.dp = s.gb = s.alphaDB = 1.0; s.betaDB = 1;
	if (pTrfrData != 0) {

		auto& pSMS = dynamic_cast<const smsData&>(*pTrfrData);

		s.dp = pSMS.dp; s.gb = pSMS.gb;
		s.alphaDB = pSMS.alphaDB; s.betaDB = pSMS.betaDB;
	}

	return s;
}


// Set phenotypic transfer by CRW traits
void Individual::setIndCRWTraits(Species* pSpecies) {
	trfrCRWTraits c; c.stepLength = c.rho = 0.0;

	c.stepLength = getTrait(CRW_STEPLENGTH)->express();
	c.rho = getTrait(CRW_STEPCORRELATION)->express();

	auto& pCRW = dynamic_cast<crwData&>(*pTrfrData);
	pCRW.stepLength = (float)(c.stepLength);
	pCRW.rho = (float)(c.rho);
	if (pCRW.stepLength < 1.0) pCRW.stepLength = 1.0;
	if (pCRW.rho < 0.0) pCRW.rho = 0.0;
	if (pCRW.rho > 0.999) pCRW.rho = 0.999f;
	return;
}

// Get phenotypic transfer by CRW traits
trfrCRWTraits Individual::getIndCRWTraits(void) {

	trfrCRWTraits c; 
	c.stepLength = c.rho = 0.0;
	if (pTrfrData != 0) {
		auto& pCRW = dynamic_cast<const crwData&>(*pTrfrData);
		c.stepLength = pCRW.stepLength;
		c.rho = pCRW.rho;
	}
	return c;

}

// Get phenotypic settlement traits
settleTraits Individual::getIndSettTraits(void) {
	settleTraits s; s.s0 = s.alpha = s.beta = 0.0;
	if (pSettleTraits != 0) {
		s.s0 = pSettleTraits->s0;
		s.alpha = pSettleTraits->alpha;
		s.beta = pSettleTraits->beta;
	}

	return s;
}


void Individual::setStatus(short s) {
	if (s >= 0 && s <= 9) status = s;
	status = s;
}

void Individual::developing(void) {
	isDeveloping = true;
}

void Individual::develop(void) {
	stage++; isDeveloping = false;
}

void Individual::ageIncrement(short maxage) {
	if (status < 6) { // alive
		age++;
		if (age > maxage) status = 9;			// exceeds max. age - dies
		else {
			if (path != 0) path->year = 0;	// reset annual step count for movement models
			if (status == 3) // waiting to continue dispersal
				status = 1;
		}
	}
}

void Individual::incFallow(void) { fallow++; }

void Individual::resetFallow(void) { fallow = 0; }

//---------------------------------------------------------------------------
// Move to a specified neighbouring cell
void Individual::moveto(Cell* newCell) {
	// check that location is indeed a neighbour of the current cell
	locn currloc = pCurrCell->getLocn();
	locn newloc = newCell->getLocn();
	double d = sqrt(((double)currloc.x - (double)newloc.x) * ((double)currloc.x - (double)newloc.x)
		+ ((double)currloc.y - (double)newloc.y) * ((double)currloc.y - (double)newloc.y));
	if (d >= 1.0 && d < 1.5) { // ok
		pCurrCell = newCell; status = 5;
	}
}

//---------------------------------------------------------------------------
// Move to a new cell by sampling a dispersal distance from a single or double
// negative exponential kernel
// Returns 1 if still dispersing (including having found a potential patch), otherwise 0
int Individual::moveKernel(Landscape* pLandscape, Species* pSpecies, const bool absorbing)
{
	intptr patch;
	int patchNum = 0;
	int newX = 0, newY = 0;
	int dispersing = 1;
	double xrand, yrand, meandist, dist, r1, rndangle, nx, ny;
	float localK;
	trfrKernelParams kern;
	Cell* pCell;
	Patch* pPatch;
	locn loc = pCurrCell->getLocn();

	landData land = pLandscape->getLandData();

	bool usefullkernel = pSpecies->useFullKernel();
	transferRules trfr = pSpecies->getTransferRules();
	settleRules sett = pSpecies->getSettRules(stage, sex);

	pCell = NULL;
	pPatch = NULL;

	if (trfr.indVar) { // get individual's kernel parameters
		kern.meanDist1 = kern.meanDist2 = kern.probKern1 = 0.0;

		auto& pKernel = dynamic_cast<const kernelData&>(*pTrfrData);

		kern.meanDist1 = pKernel.meanDist1;
		if (trfr.twinKern)
		{
			kern.meanDist2 = pKernel.meanDist2;
			kern.probKern1 = pKernel.probKern1;
		}
	}
	else { // get kernel parameters for the species
		if (trfr.sexDep) {
			if (trfr.stgDep) {
				kern = pSpecies->getSpKernTraits(stage, sex);
			}
			else {
				kern = pSpecies->getSpKernTraits(0, sex);
			}
		}
		else {
			if (trfr.stgDep) {
				kern = pSpecies->getSpKernTraits(stage, 0);
			}
			else {
				kern = pSpecies->getSpKernTraits(0, 0);
			}
		}
	}

	// scale the appropriate kernel mean to the cell size
	if (trfr.twinKern)
	{
		if (pRandom->Bernoulli(kern.probKern1))
			meandist = kern.meanDist1 / (float)land.resol;
		else
			meandist = kern.meanDist2 / (float)land.resol;
	}
	else
		meandist = kern.meanDist1 / (float)land.resol;
	// scaled mean may not be less than 1 unless emigration derives from the kernel
	// (i.e. the 'use full kernel' option is applied)
	if (!usefullkernel && meandist < 1.0) meandist = 1.0;

	int loopsteps = 0; // new counter to prevent infinite loop added 14/8/15
	do {
		do {
			do {
				// randomise the cell within the patch, provided that the individual is still in
				// its natal cell (i.e. not waiting in the matrix)
				// this is because, if the patch is very large, the individual is near the centre
				// and the (single) kernel mean is (not much more than) the cell size, an infinite
				// loop could otherwise result, as the individual never reaches the patch edge
				// (in a cell-based model, this has no effect, other than as a processing overhead)
				if (status == 1) {
					pCell = pNatalPatch->getRandomCell();
					if (pCell != 0) {
						loc = pCell->getLocn();
					}
				}
				// randomise the position of the individual inside the cell
				// so x and y are a corner of the cell?
				xrand = (double)loc.x + pRandom->Random() * 0.999;
				yrand = (double)loc.y + pRandom->Random() * 0.999;

				// draw factor r1 0 < r1 <= 1
				r1 = 0.0000001 + pRandom->Random() * (1.0 - 0.0000001);
				dist = (-1.0 * meandist) * log(r1);

				rndangle = pRandom->Random() * 2.0 * PI;
				nx = xrand + dist * sin(rndangle);
				ny = yrand + dist * cos(rndangle);
				if (nx < 0.0) newX = -1; else newX = (int)nx;
				if (ny < 0.0) newY = -1; else newY = (int)ny;
#if RSDEBUG
				if (path != 0) (path->year)++;
#endif
				loopsteps++;
			} while (loopsteps < 1000 &&
				// keep drawing if out of bounds of landscape or same cell
				((!absorbing && (newX < land.minX || newX > land.maxX
					|| newY < land.minY || newY > land.maxY))
					|| (!usefullkernel && newX == loc.x && newY == loc.y))
				);
			if (loopsteps < 1000) {
				if (newX < land.minX || newX > land.maxX
					|| newY < land.minY || newY > land.maxY) { // beyond absorbing boundary
					// this cannot be reached if not absorbing?
					pCell = 0;
					patch = 0;
					patchNum = -1;
				}
				else {
					pCell = pLandscape->findCell(newX, newY);
					if (pCell == 0) { // no-data cell
						patch = 0;
						patchNum = -1;
					}
					else {
						patch = pCell->getPatch();
						if (patch == 0) { // matrix
							pPatch = 0;
							patchNum = 0;
						}
						else {
							pPatch = (Patch*)patch;
							patchNum = pPatch->getPatchNum();
						}
					}
				}
			}
			else { // exceeded 1000 attempts
				patch = 0;
				patchNum = -1;
			}
		} while (!absorbing && patchNum < 0 && loopsteps < 1000); 			 // in a no-data region
	} while (!usefullkernel && pPatch == pNatalPatch && loopsteps < 1000); 	// still in the original (natal) patch

	if (loopsteps < 1000) {
		if (pCell == 0) { // beyond absorbing boundary or in no-data cell
			// only if absorbing=true and out of bounddaries
			pCurrCell = 0;
			status = 6;
			dispersing = 0;
		}
		else {
			pCurrCell = pCell;
			if (pPatch == 0) localK = 0.0; // matrix
			else localK = pPatch->getK();
			if (patchNum > 0 && localK > 0.0) { // found a new patch
				status = 2; // record as potential settler
			}
			else {
				// unsuitable patch
				dispersing = 0;
				// can wait in matrix if population is stage structured ...
				if (pSpecies->stageStructured()) {
					// ... and wait option is applied ...
					if (sett.wait) { // ... it is
						status = 3; // waiting
					}
					else // ... it is not
						status = 6; // dies (unless there is a suitable neighbouring cell)
				}
				else status = 6; // dies (unless there is a suitable neighbouring cell)
			}
		}
	}
	else { // exceeded 1000 attempts
		status = 6;
		dispersing = 0;
	}

	// apply dispersal-related mortality, which may be distance-dependent
	dist *= (float)land.resol; // re-scale distance moved to landscape scale
	if (status < 7) {
		double dispmort;
		trfrMortParams mort = pSpecies->getMortParams();
		if (trfr.distMort) {
			dispmort = 1.0 / (1.0 + exp(-(dist - mort.mortBeta) * mort.mortAlpha));
		}
		else {
			dispmort = mort.fixedMort;
		}
		if (pRandom->Bernoulli(dispmort)) {
			status = 7; // dies
			dispersing = 0;
		}
	}

	return dispersing;
}

//---------------------------------------------------------------------------
// Make a single movement step according to a mechanistic movement model
// Returns 1 if still dispersing (including having found a potential patch), otherwise 0
int Individual::moveStep(Landscape* pLandscape, Species* pSpecies,
	const short landIx, const bool absorbing)
{
	if (status != 1) return 0; // not currently dispersing

	intptr patch;
	int patchNum;
	int newX, newY;
	locn loc;
	int dispersing = 1;
	double xcnew, ycnew;
	double angle;
	double mortprob, rho, steplen;
	movedata move;
	Patch* pPatch = 0;
	bool absorbed = false;
	//int popsize;

	landData land = pLandscape->getLandData();
	simParams sim = paramsSim->getSim();

	transferRules trfr = pSpecies->getTransferRules();
	trfrCRWTraits movt = pSpecies->getSpCRWTraits();
	settleSteps settsteps = pSpecies->getSteps(stage, sex);

	patch = pCurrCell->getPatch();

	if (patch == 0) { // matrix
		pPatch = 0;
		patchNum = 0;
	}
	else {
		pPatch = (Patch*)patch;
		patchNum = pPatch->getPatchNum();
	}
	// apply step-dependent mortality risk ...
	if (trfr.habMort)
	{ // habitat-dependent
		int h = pCurrCell->getHabIndex(landIx);
		if (h < 0) { // no-data cell - should not occur, but if it does, individual dies
			mortprob = 1.0;
		}
		else mortprob = pSpecies->getHabMort(h);
	}
	else mortprob = movt.stepMort;
	// ... unless individual has not yet left natal patch in emigration year
	if (pPatch == pNatalPatch && path->out == 0 && path->year == path->total) {
		mortprob = 0.0;
	}
	if (pRandom->Bernoulli(mortprob)) { // individual dies
		status = 7;
		dispersing = 0;
	}
	else { // take a step
		(path->year)++;
		(path->total)++;
		//	if (pPatch != pNatalPatch || path->out > 0) (path->out)++;
		if (patch == 0 || pPatch == 0 || patchNum == 0) { // not in a patch
			if (path != 0) path->settleStatus = 0; // reset path settlement status
			(path->out)++;
		}
		loc = pCurrCell->getLocn();
		newX = loc.x; newY = loc.y;

		switch (trfr.moveType) {

		case 1: // SMS
			move = smsMove(pLandscape, pSpecies, landIx, pPatch == pNatalPatch, trfr.indVar, absorbing);
			if (move.dist < 0.0) {
				// either INTERNAL ERROR CONDITION - INDIVIDUAL IS IN NO-DATA SQUARE
				// or individual has crossed absorbing boundary ...
				// ... individual dies
				status = 6;
				dispersing = 0;
			}
			else {

				// WOULD IT BE MORE EFFICIENT FOR smsMove TO RETURN A POINTER TO THE NEW CELL? ...

				patch = pCurrCell->getPatch();
				//int patchnum;
				if (patch == 0) {
					pPatch = 0;
					//patchnum = 0;
				}
				else {
					pPatch = (Patch*)patch;
					//patchnum = pPatch->getPatchNum();
				}
				if (sim.saveVisits && pPatch != pNatalPatch) {
					pCurrCell->incrVisits();
				}
			}
			break;

		case 2: // CRW

			auto & pCRW = dynamic_cast<crwData&>(*pTrfrData);

			if (trfr.indVar) {
				movt.stepLength = pCRW.stepLength;
				movt.rho = pCRW.rho;
			}

			steplen = movt.stepLength; 
			rho = movt.rho;
			if (pPatch == pNatalPatch) {
				rho = 0.99; // to promote leaving natal patch
				path->out = 0;
			}
			if (movt.straightenPath && path->settleStatus > 0) {
				// individual is in a patch and has already determined whether to settle
				rho = 0.99; // to promote leaving the patch
				path->out = 0;
			}
			int loopsteps = 0; // new counter to prevent infinite loop added 14/8/15
			do {
				do {
					// new direction
					if (newX < land.minX || newX > land.maxX || newY < land.minY || newY > land.maxY
						|| pCurrCell == 0) {
						// individual has tried to go out-of-bounds or into no-data area
						// allow random move to prevent repeated similar move
						angle = wrpcauchy(pCRW.prevdrn, 0.0);
					}
					else
						angle = wrpcauchy(pCRW.prevdrn, rho);
					// new continuous cell coordinates
					xcnew = pCRW.xc + sin(angle) * steplen / (float)land.resol;
					ycnew = pCRW.yc + cos(angle) * steplen / (float)land.resol;
					if (xcnew < 0.0) newX = -1; else newX = (int)xcnew;
					if (ycnew < 0.0) newY = -1; else newY = (int)ycnew;
					loopsteps++;
				} while (!absorbing && loopsteps < 1000 &&
					(newX < land.minX || newX > land.maxX || newY < land.minY || newY > land.maxY));
				if (newX < land.minX || newX > land.maxX || newY < land.minY || newY > land.maxY)
					pCurrCell = 0;
				else
					pCurrCell = pLandscape->findCell(newX, newY);
				if (pCurrCell == 0) { // no-data cell or beyond absorbing boundary
					patch = 0;
					if (absorbing) absorbed = true;
				}
				else
					patch = pCurrCell->getPatch();
			} while (!absorbing && pCurrCell == 0 && loopsteps < 1000);
			pCRW.prevdrn = (float)angle;
			pCRW.xc = (float)xcnew; pCRW.yc = (float)ycnew;
			if (absorbed) { // beyond absorbing boundary or in no-data square
				status = 6;
				dispersing = 0;
				pCurrCell = 0;
			}
			else {
				if (loopsteps >= 1000) { // unable to make a move
					// INTERNAL ERROR CONDITION - INDIVIDUAL IS IN NO-DATA SQUARE
					// NEED TO TAKE SOME FORM OF INFORMATIVE ACTION ...
					// ... individual dies as it cannot move
					status = 6;
					dispersing = 0;
					// current cell will be invalid (zero), so set back to previous cell
					pCurrCell = pPrevCell;
				}
			}
			break;

		} // end of switch (trfr.moveType)

		if (patch > 0  // not no-data area or matrix
			&& path->total >= settsteps.minSteps) {
			pPatch = (Patch*)patch;
			if (pPatch != pNatalPatch)
			{
				// determine whether the new patch is potentially suitable
				if (pPatch->getK() > 0.0)
				{ // patch is suitable
					status = 2;
				}
			}
		}
		if (status != 2 && status != 6) { // suitable patch not found, not already dead
			if (path->year >= settsteps.maxStepsYr) {
				status = 3; // waits until next year
			}
			if (path->total >= settsteps.maxSteps) {
				status = 6; // dies
				dispersing = 0;
			}
		}
	} // end of single movement step

	return dispersing;
}

//---------------------------------------------------------------------------

// Functions to implement the SMS algorithm

// Move to a neighbouring cell according to the SMS algorithm
movedata Individual::smsMove(Landscape* pLand, Species* pSpecies,
	const short landIx, const bool natalPatch, const bool indvar, const bool absorbing)
{
	array3x3d nbr; 	// to hold weights/costs/probs of moving to neighbouring cells
	array3x3d goal;	// to hold weights for moving towards a goal location
	array3x3f hab;	// to hold weights for habitat (includes percep range)
	int x2, y2; 			// x index from 0=W to 2=E, y index from 0=N to 2=S
	int newX = 0, newY = 0;
	Cell* pCell;
	Cell* pNewCell = NULL;
	double sum_nbrs = 0.0;
	movedata move;
	int cellcost, newcellcost;
	locn current;

	auto& pSMS = dynamic_cast<smsData&>(*pTrfrData);
	if (pCurrCell == 0)
	{
		// x,y is a NODATA square - this should not occur here
		// return a negative distance to indicate an error
		move.dist = -69.0; move.cost = 0.0;
		return move;
	}

	landData land = pLand->getLandData();
	trfrSMSTraits movt = pSpecies->getSpSMSTraits();
	current = pCurrCell->getLocn();

	//get weights for directional persistence....
	if ((path->out > 0 && path->out <= (movt.pr + 1))
		|| natalPatch
		|| (movt.straightenPath && path->settleStatus > 0)) {
		// inflate directional persistence to promote leaving the patch
		if (indvar) nbr = getSimDir(current.x, current.y, 10.0f * pSMS.dp);
		else nbr = getSimDir(current.x, current.y, 10.0f * movt.dp);
	}
	else {
		if (indvar) nbr = getSimDir(current.x, current.y, pSMS.dp);
		else nbr = getSimDir(current.x, current.y, movt.dp);
	}
	if (natalPatch || path->settleStatus > 0) path->out = 0;

	//get weights for goal bias....
	double gb;
	if (movt.goalType == 2) { // dispersal bias
		int nsteps = 0;
		if (path->year == path->total) { // first year of dispersal - use no. of steps outside natal patch
			nsteps = path->out;
		}
		else { // use total no. of steps
			nsteps = path->total;
		}
		if (indvar) {
			double exp_arg = -((double)nsteps - (double)pSMS.betaDB) * (-pSMS.alphaDB);
			if (exp_arg > 100.0) exp_arg = 100.0; // to prevent exp() overflow error
			gb = 1.0 + (pSMS.gb - 1.0) / (1.0 + exp(exp_arg));
		}
		else {
			double exp_arg = -((double)nsteps - (double)movt.betaDB) * (-movt.alphaDB);

			if (exp_arg > 100.0) exp_arg = 100.0; // to prevent exp() overflow error
			gb = 1.0 + (movt.gb - 1.0) / (1.0 + exp(exp_arg));
		}
	}
	else gb = movt.gb;
	goal = getGoalBias(current.x, current.y, movt.goalType, (float)gb);

	// get habitat-dependent weights (mean effective costs, given perceptual range)
	// first check if costs have already been calculated

	hab = pCurrCell->getEffCosts();
	if (hab.cell[0][0] < 0.0) { // costs have not already been calculated
		hab = getHabMatrix(pLand, pSpecies, current.x, current.y, movt.pr, movt.prMethod,
			landIx, absorbing);
		pCurrCell->setEffCosts(hab);
	}
	else { // they have already been calculated - no action required

	}

	// determine weighted effective cost for the 8 neighbours
	// multiply directional persistence, goal bias and habitat habitat-dependent weights
	for (y2 = 2; y2 > -1; y2--) {
		for (x2 = 0; x2 < 3; x2++) {
			if (x2 == 1 && y2 == 1) nbr.cell[x2][y2] = 0.0;
			else {
				if (x2 == 1 || y2 == 1) //not diagonal
					nbr.cell[x2][y2] = nbr.cell[x2][y2] * goal.cell[x2][y2] * hab.cell[x2][y2];
				else // diagonal
					nbr.cell[x2][y2] = (float)SQRT2 * nbr.cell[x2][y2] * goal.cell[x2][y2] * hab.cell[x2][y2];
			}
		}
	}

	// determine reciprocal of effective cost for the 8 neighbours
	for (y2 = 2; y2 > -1; y2--) {
		for (x2 = 0; x2 < 3; x2++) {
			if (nbr.cell[x2][y2] > 0.0) nbr.cell[x2][y2] = 1.0f / nbr.cell[x2][y2];
		}
	}

	// set any cells beyond the current landscape limits and any no-data cells
	// to have zero probability
	// increment total for re-scaling to sum to unity

	for (y2 = 2; y2 > -1; y2--) {
		for (x2 = 0; x2 < 3; x2++) {
			if (!absorbing) {
				if ((current.y + y2 - 1) < land.minY || (current.y + y2 - 1) > land.maxY
					|| (current.x + x2 - 1) < land.minX || (current.x + x2 - 1) > land.maxX)
					// cell is beyond current landscape limits
					nbr.cell[x2][y2] = 0.0;
				else { // check if no-data cell
					pCell = pLand->findCell((current.x + x2 - 1), (current.y + y2 - 1));
					if (pCell == 0) nbr.cell[x2][y2] = 0.0; // no-data cell
				}
			}
			sum_nbrs += nbr.cell[x2][y2];
		}
	}

	// scale effective costs as probabilities summing to 1
	if (sum_nbrs > 0.0) { // should always be the case, but safest to check...
		for (y2 = 2; y2 > -1; y2--) {
			for (x2 = 0; x2 < 3; x2++) {
				nbr.cell[x2][y2] = nbr.cell[x2][y2] / (float)sum_nbrs;
			}
		}
	}

	// set up cell selection probabilities
	double cumulative[9];
	int j = 0;
	cumulative[0] = nbr.cell[0][0];
	for (y2 = 0; y2 < 3; y2++) {
		for (x2 = 0; x2 < 3; x2++) {
			if (j != 0) cumulative[j] = cumulative[j - 1] + nbr.cell[x2][y2];
			j++;
		}
	}

	//to prevent very rare bug that random draw is greater than 0.999999999
	if (cumulative[8] != 1) cumulative[8] = 1;
	// select direction at random based on cell selection probabilities
	// landscape boundaries and no-data cells may be reflective or absorbing
	cellcost = pCurrCell->getCost();
	int loopsteps = 0; // new counter to prevent infinite loop added 14/8/15
	do {
		do {
			double rnd = pRandom->Random();
			j = 0;
			for (y2 = 0; y2 < 3; y2++) {
				for (x2 = 0; x2 < 3; x2++) {
					if (rnd < cumulative[j]) {
						newX = current.x + x2 - 1;
						newY = current.y + y2 - 1;
						if (x2 == 1 || y2 == 1) move.dist = (float)(land.resol);
						else move.dist = (float)(land.resol) * (float)SQRT2;
						y2 = 999; x2 = 999; //to break out of x2 and y2 loops.
					}
					j++;
				}
			}
			loopsteps++;
		} while (loopsteps < 1000
			&& (!absorbing && (newX < land.minX || newX > land.maxX
				|| newY < land.minY || newY > land.maxY)));
		if (loopsteps >= 1000) pNewCell = 0;
		else {
			if (newX < land.minX || newX > land.maxX
				|| newY < land.minY || newY > land.maxY) {
				pNewCell = 0;
			}
			pNewCell = pLand->findCell(newX, newY);
		}
	} while (!absorbing && pNewCell == 0 && loopsteps < 1000); // no-data cell
	if (loopsteps >= 1000 || pNewCell == 0) {
		// unable to make a move or crossed absorbing boundary
		// flag individual to die
		move.dist = -123.0;
		if (pNewCell == 0) pCurrCell = pNewCell;
	}
	else {
		newcellcost = pNewCell->getCost();
		move.cost = move.dist * 0.5f * ((float)cellcost + (float)newcellcost);
		// make the selected move
		if ((short)memory.size() == movt.memSize) {
			memory.pop(); // remove oldest memory element
		}
		memory.push(current); // record previous location in memory
		//if (write_out) out << "queue length is " << memory.size() << endl;
		pCurrCell = pNewCell;
	}
	return move;
}

// Weight neighbouring cells on basis of current movement direction
array3x3d Individual::getSimDir(const int x, const int y, const float dp)
{

	array3x3d d;
	locn prev;
	double theta;
	int xx, yy;

	//if (write_out) out<<"step 0"<<endl;
	if (memory.empty())
	{ // no previous movement, set matrix to unity
		for (xx = 0; xx < 3; xx++) {
			for (yy = 0; yy < 3; yy++) {
				d.cell[xx][yy] = 1;
			}
		}
	}
	else { // set up the matrix dependent on relationship of previous location to current
		//  if (write_out) out<<"step 1"<<endl;
		d.cell[1][1] = 0;
		prev = memory.front();
		//  if (write_out) out<<"step 2"<<endl;
		if ((x - prev.x) == 0 && (y - prev.y) == 0) {
			// back to 'square 1' (first memory location) - use previous step drn only
			prev = memory.back();
			if ((x - prev.x) == 0 && (y - prev.y) == 0) { // STILL HAVE A PROBLEM!
				for (xx = 0; xx < 3; xx++) {
					for (yy = 0; yy < 3; yy++) {
						d.cell[xx][yy] = 1.0;
					}
				}
				return d;
			}
		}
		else {
			//    if (write_out) out<<"step 5"<<endl;
		}
		//  if (write_out) out<<"step 6"<<endl;
		theta = atan2(((double)x - (double)prev.x), ((double)y - (double)prev.y));
		//  if (write_out) out<<"prev.x,prev.y: "<<prev.x<<","<<prev.y<<" theta: "<<theta<<endl;
		d = calcWeightings(dp, (float)theta);

	}
	return d;
}

// Weight neighbouring cells on basis of goal bias
//array3x3d Individual::getGoalBias(const int x,const int y,
//	const int goaltype,const float gb)
array3x3d Individual::getGoalBias(const int x, const int y,
	const int goaltype, const float gb)
{

	array3x3d d;
	double theta;
	int xx, yy;
	auto& pSMS = dynamic_cast<const smsData&>(*pTrfrData);

	if (goaltype == 0) { // no goal set
		for (xx = 0; xx < 3; xx++) {
			for (yy = 0; yy < 3; yy++) {
				d.cell[xx][yy] = 1.0;
			}
		}
	}
	else {
		d.cell[1][1] = 0;
		if ((x - pSMS.goal.x) == 0 && (y - pSMS.goal.y) == 0) {
			// at goal, set matrix to unity
			for (xx = 0; xx < 3; xx++) {
				for (yy = 0; yy < 3; yy++) {
					d.cell[xx][yy] = 1.0;
				}
			}
			return d;
		}
		if (goaltype == 1) {
			// TEMPORARY CODE - GOAL TYPE 1 NOT YET IMPLEMENTED, AS WE HAVE NO MEANS OF
			// CAPTURING THE GOAL LOCATION OF EACH INDIVIDUAL
			for (xx = 0; xx < 3; xx++) {
				for (yy = 0; yy < 3; yy++) {
					d.cell[xx][yy] = 1.0;
				}
			}
			return d;
		}
		else // goaltype == 2
			theta = atan2(((double)x - (double)pSMS.goal.x), ((double)y - (double)pSMS.goal.y));
		//  if (write_out) out<<"goalx,goaly: "<<goalx<<","<<goaly<<" theta: "<<theta<<endl;
		d = calcWeightings(gb, (float)theta);
	}

	return d;
}

// Calculate weightings for neighbouring cells
array3x3d Individual::calcWeightings(const double base, const double theta) {

	array3x3d d; // 3x3 array indexed from SW corner by xx and yy
	int dx, dy, xx, yy;

	double i0 = 1.0; 					// direction of theta - lowest cost bias
	double i1 = base;
	double i2 = base * base;
	double i3 = i2 * base;
	double i4 = i3 * base;		// opposite to theta - highest cost bias

	if (fabs(theta) > 7.0 * PI / 8.0) { dx = 0; dy = -1; }
	else {
		if (fabs(theta) > 5.0 * PI / 8.0) { dy = -1; if (theta > 0) dx = 1; else dx = -1; }
		else {
			if (fabs(theta) > 3.0 * PI / 8.0) { dy = 0; if (theta > 0) dx = 1; else dx = -1; }
			else {
				if (fabs(theta) > PI / 8.0) { dy = 1; if (theta > 0) dx = 1; else dx = -1; }
				else { dy = 1; dx = 0; }
			}
		}
	}
	d.cell[1][1] = 0; // central cell has zero weighting
	d.cell[dx + 1][dy + 1] = (float)i0;
	d.cell[-dx + 1][-dy + 1] = (float)i4;
	if (dx == 0 || dy == 0) { // theta points to a cardinal direction
		d.cell[dy + 1][dx + 1] = (float)i2; d.cell[-dy + 1][-dx + 1] = (float)i2;
		if (dx == 0) { // theta points N or S
			xx = dx + 1; if (xx > 1) dx -= 2; yy = dy;
			d.cell[xx + 1][yy + 1] = (float)i1; d.cell[-xx + 1][yy + 1] = (float)i1;
			d.cell[xx + 1][-yy + 1] = (float)i3; d.cell[-xx + 1][-yy + 1] = (float)i3;
		}
		else { // theta points W or E
			yy = dy + 1; if (yy > 1) dy -= 2; xx = dx;
			d.cell[xx + 1][yy + 1] = (float)i1; d.cell[xx + 1][-yy + 1] = (float)i1;
			d.cell[-xx + 1][yy + 1] = (float)i3; d.cell[-xx + 1][-yy + 1] = (float)i3;
		}
	}
	else { // theta points to an ordinal direction
		d.cell[dx + 1][-dy + 1] = (float)i2; d.cell[-dx + 1][dy + 1] = (float)i2;
		xx = dx + 1; if (xx > 1) xx -= 2; d.cell[xx + 1][dy + 1] = (float)i1;
		yy = dy + 1; if (yy > 1) yy -= 2; d.cell[dx + 1][yy + 1] = (float)i1;
		d.cell[-xx + 1][-dy + 1] = (float)i3; d.cell[-dx + 1][-yy + 1] = (float)i3;
	}

	return d;
}

// Weight neighbouring cells on basis of (habitat) costs
array3x3f Individual::getHabMatrix(Landscape* pLand, Species* pSpecies,
	const int x, const int y, const short pr, const short prmethod, const short landIx,
	const bool absorbing)
{

	array3x3f w; // array of effective costs to be returned
	int ncells, x4, y4;
	double weight, sumweights;
	// NW and SE corners of effective cost array relative to the current cell (x,y):
	int xmin = 0, ymin = 0, xmax = 0, ymax = 0;
	int cost, nodatacost, h;
	Cell* pCell;

	landData land = pLand->getLandData();
	if (absorbing) nodatacost = gAbsorbingNoDataCost;
	else nodatacost = gNoDataCost;

	for (int x2 = -1; x2 < 2; x2++) {   // index of relative move in x direction
		for (int y2 = -1; y2 < 2; y2++) { // index of relative move in x direction

			w.cell[x2 + 1][y2 + 1] = 0.0; // initialise costs array to zeroes

			// set up corners of perceptual range relative to current cell
			if (x2 == 0 && y2 == 0) { // current cell - do nothing
				xmin = 0; ymin = 0; xmax = 0; ymax = 0;
			}
			else {
				if (x2 == 0 || y2 == 0) { // not diagonal (rook move)
					if (x2 == 0) { // vertical (N-S) move
						if (pr % 2 == 0) { xmin = -pr / 2; xmax = pr / 2; ymin = y2; ymax = y2 * pr; } // PR even
						else { xmin = -(pr - 1) / 2; xmax = (pr - 1) / 2; ymin = y2; ymax = y2 * pr; } // PR odd
					}
					if (y2 == 0) { // horizontal (E-W) move
						if (pr % 2 == 0) { xmin = x2; xmax = x2 * pr; ymin = -pr / 2; ymax = pr / 2; } // PR even
						else { xmin = x2; xmax = x2 * pr; ymin = -(pr - 1) / 2; ymax = (pr - 1) / 2; } // PR odd
					}
				}
				else { // diagonal (bishop move)
					xmin = x2; xmax = x2 * pr; ymin = y2; ymax = y2 * pr;
				}
			}
			if (xmin > xmax) { int z = xmax; xmax = xmin; xmin = z; } // swap xmin and xmax
			if (ymin > ymax) { int z = ymax; ymax = ymin; ymin = z; } // swap ymin and ymax

			// calculate effective mean cost of cells in perceptual range
			ncells = 0; weight = 0.0; sumweights = 0.0;
			//		targetseen = 0;
			if (x2 != 0 || y2 != 0) { // not central cell (i.e. current cell)
				for (int x3 = xmin; x3 <= xmax; x3++) {
					for (int y3 = ymin; y3 <= ymax; y3++) {
						// if cell is out of bounds, treat landscape as a torus
						// for purpose of obtaining a cost,
						if ((x + x3) < 0) x4 = x + x3 + land.maxX + 1;
						else { if ((x + x3) > land.maxX) x4 = x + x3 - land.maxX - 1; else x4 = x + x3; }
						if ((y + y3) < 0) y4 = y + y3 + land.maxY + 1;
						else { if ((y + y3) > land.maxY) y4 = y + y3 - land.maxY - 1; else y4 = y + y3; }
						//					if (write_out && (x4 < 0 || y4 < 0)) {
						//						out<<"ERROR: x "<<x<<" y "<<y<<" x3 "<<x3<<" y3 "<<y3
						//							<<" xbound "<<xbound<<" ybound "<<ybound<<" x4 "<<x4<<" y4 "<<y4<<endl;
						//					}
						if (x4 < 0 || x4 > land.maxX || y4 < 0 || y4 > land.maxY) {
							// unexpected problem - e.g. due to ridiculously large PR
							// treat as a no-data cell
							cost = nodatacost;
						}
						else {
							// add cost of cell to total PR cost
							pCell = pLand->findCell(x4, y4);
							if (pCell == 0) { // no-data cell
								cost = nodatacost;
							}
							else {
								cost = pCell->getCost();
								if (cost < 0) cost = nodatacost;
								else {
									if (cost == 0) { // cost not yet set for the cell
										h = pCell->getHabIndex(landIx);
										cost = pSpecies->getHabCost(h);
										pCell->setCost(cost);
									}
									else {
										// nothing?
									}
								}
							}
						}
						if (prmethod == 1) { // arithmetic mean
							w.cell[x2 + 1][y2 + 1] += cost;
							ncells++;
						}
						if (prmethod == 2) { // harmonic mean
							if (cost > 0) {
								w.cell[x2 + 1][y2 + 1] += (1.0f / (float)cost);
								ncells++;
							}
						}
						if (prmethod == 3) { // arithmetic mean weighted by inverse distance
							if (cost > 0) {
								// NB distance is still given by (x3,y3)
								weight = 1.0f / (double)sqrt((pow((double)x3, 2) + pow((double)y3, 2)));
								w.cell[x2 + 1][y2 + 1] += (float)(weight * (double)cost);
								ncells++; sumweights += weight;
							}
						}
					} //end of y3 loop
				}  //end of x3 loop
				if (ncells > 0) {
					if (prmethod == 1) w.cell[x2 + 1][y2 + 1] /= ncells; // arithmetic mean
					if (prmethod == 2) w.cell[x2 + 1][y2 + 1] = ncells / w.cell[x2 + 1][y2 + 1]; // hyperbolic mean
					if (prmethod == 3 && sumweights > 0)
						w.cell[x2 + 1][y2 + 1] /= (float)sumweights; // weighted arithmetic mean
				}
			}
			else { // central cell
				// record cost if not already recorded
				// has effect of preparing for storing effective costs for the cell
				pCell = pLand->findCell(x, y);
				cost = pCell->getCost();
				if (cost < 0) cost = nodatacost;
				else {
					if (cost == 0) { // cost not yet set for the cell
						h = pCell->getHabIndex(landIx);
						cost = pSpecies->getHabCost(h);
						pCell->setCost(cost);
					}
				}
			}
			//		if (write_out2) out2<<"effective mean cost "<<w.cell[x2+1][y2+1]<<endl;

		}//end of y2 loop
	}//end of x2 loop

	return w;

}

#if RS_RCPP
//---------------------------------------------------------------------------
// Write records to movement paths file
void Individual::outMovePath(const int year)
{
	locn loc, prev_loc;

	//if (pPatch != pNatalPatch) {
	loc = pCurrCell->getLocn();
	// if still dispersing...
	if (status == 1) {
		// at first step, record start cell first
		if (path->total == 1) {
			prev_loc = pPrevCell->getLocn();
			outMovePaths << year << "\t" << indId << "\t"
				<< "0\t" << prev_loc.x << "\t" << prev_loc.y << "\t"
				<< "0\t"	// status at start cell is 0
				<< endl;
		}
		// then record current step
		outMovePaths << year << "\t" << indId << "\t"
			<< path->total << "\t" << loc.x << "\t" << loc.y << "\t"
			<< status << "\t"
			<< endl;
	}
	// if not anymore dispersing...
	if (status > 1 && status < 10) {
		prev_loc = pPrevCell->getLocn();
		// record only if this is the first step as non-disperser
		if (path->pathoutput) {
			// if this is also the first step taken at all, record the start cell first
			if (path->total == 1) {
				outMovePaths << year << "\t" << indId << "\t"
					<< "0\t" << prev_loc.x << "\t" << prev_loc.y << "\t"
					<< "0\t"	// status at start cell is 0
					<< endl;
			}
			outMovePaths << year << "\t" << indId << "\t"
				<< path->total << "\t" << loc.x << "\t" << loc.y << "\t"
				<< status << "\t"
				<< endl;
			// current cell will be invalid (zero), so set back to previous cell
			//pPrevCell = pCurrCell;
			path->pathoutput = 0;
		}
	}
}
#endif

//---------------------------------------------------------------------------

//---------------------------------------------------------------------------

double wrpcauchy(double location, double rho) {
	double result;

	if (rho < 0.0 || rho > 1.0) {
		result = location;
	}

	if (rho == 0)
		result = pRandom->Random() * M_2PI;
	else
		if (rho == 1) result = location;
		else {
			result = fmod(cauchy(location, -log(rho)), M_2PI);
		}
	return result;
}

double cauchy(double location, double scale) {
	if (scale < 0) return location;
	return location + scale * tan(PI * pRandom->Random());
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

#if RSDEBUG
// Testing utilities

Cell* Individual::getCurrCell() const {
	return pCurrCell;
}

void Individual::setInitAngle(const float angle) {
	auto pCRW = dynamic_cast<crwData*>(pTrfrData.get());
	pCRW->prevdrn = angle;
}

// Force mutations to trigger for all traits
void Individual::triggerMutations(Species* pSp) {
	for (auto const& [trType, indTrait] : spTraitTable) {
		indTrait->mutate();
		if (trType == GENETIC_LOAD1 
			|| trType == GENETIC_LOAD2
			|| trType == GENETIC_LOAD3
			|| trType == GENETIC_LOAD4
			|| trType == GENETIC_LOAD5)
			geneticFitness *= indTrait->express();
	}
	this->setDispersalPhenotypes(pSp, 1.0);
}

// Shorthand function to edit a genotype with custom values
void Individual::overrideGenotype(TraitType whichTrait, const map<int, vector<shared_ptr<Allele>>>& newGenotype) {

	GeneticFitnessTrait* pGenFitTrait;
	DispersalTrait* pDispTrait;

	switch (whichTrait)
	{
	case GENETIC_LOAD1: case GENETIC_LOAD2: case GENETIC_LOAD3: case GENETIC_LOAD4: case GENETIC_LOAD5: 
		pGenFitTrait = dynamic_cast<GeneticFitnessTrait*>(this->getTrait(whichTrait));
		pGenFitTrait->getGenes() = newGenotype;
		break;
	case E_D0: case E_ALPHA: case E_BETA:
	case S_S0: case S_ALPHA: case S_BETA:
	case E_D0_F: case E_ALPHA_F: case E_BETA_F:
	case S_S0_F: case S_ALPHA_F: case S_BETA_F: 
	case E_D0_M: case E_ALPHA_M: case E_BETA_M: 
	case S_S0_M: case S_ALPHA_M: case S_BETA_M: 
	case CRW_STEPLENGTH: case CRW_STEPCORRELATION: 
	case KERNEL_MEANDIST_1: case KERNEL_MEANDIST_2: case KERNEL_PROBABILITY: 
	case KERNEL_MEANDIST_1_F: case KERNEL_MEANDIST_2_F: case KERNEL_PROBABILITY_F: 
	case KERNEL_MEANDIST_1_M: case KERNEL_MEANDIST_2_M: case KERNEL_PROBABILITY_M: 
	case SMS_DP: case SMS_GB: case SMS_ALPHADB: case SMS_BETADB:
		pDispTrait = dynamic_cast<DispersalTrait*>(this->getTrait(whichTrait));
		pDispTrait->getGenes() = newGenotype;
		break;
	default:
		throw logic_error("Wrong trait type: please choose a valid dispersal or genetic fitness trait.");
		break;
	}
};

#endif // RSDEBUG

