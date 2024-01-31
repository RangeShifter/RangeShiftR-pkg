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
#if RSDEBUG
	//DEBUGLOG << "Individual::Individual(): indId=" << indId
	//	<< " stg=" << stg << " a=" << a << " probmale=" << probmale
	//	<< endl;
#endif
	fitness = 1.0;
	stage = stg;
	if (probmale <= 0.0) sex = FEM;
	else sex = pRandom->Bernoulli(probmale) ? FEM : MAL;
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
		//	path->leftNatalPatch = false;
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
	//pEmigTraits = make_unique<emigTraits>(nullptr);
	//pSettleTraits = make_unique<settleTraits>(nullptr);
#if RSDEBUG
//locn currloc = pCurrCell->getLocn();
//DEBUGLOG << "Individual::Individual(): indId=" << indId
//	<< " x=" << currloc.x << " y=" << currloc.y
////	<< " smsData=" << smsData << " dp=" << smsData->dp
//	<< endl;
#endif
}

Individual::~Individual(void) {
	if (path != 0) delete path;
	//if (crw != 0) delete crw;
	//if (smsData != 0) delete smsData;
	//if (emigtraits != 0) delete emigtraits;
	//if (kerntraits != 0) delete kerntraits;
	//if (setttraits != 0) delete setttraits;

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


TTrait* Individual::getTrait(TraitType trait) const {
	return this->traitTable.find(trait)->second.get();
}

//map<TraitType, std::unique_ptr<TTrait>>  Individual::getTraitTable(void) const
//{
//	return traitTable;
//}

set<TraitType> Individual::getTraitTypes() {
	auto kv = std::views::keys(this->traitTable);
	set< TraitType > keys{ kv.begin(), kv.end() };
	return keys;
}

//---------------------------------------------------------------------------

//---------------------------------------------------------------------------


void Individual::inherit(Species* pSpecies, const Individual* mother, const Individual* father) {

	int events = 0;
	const set<int> chromosomeEnds = pSpecies->getChromosomeEnds();
	const int genomeSize = pSpecies->getGenomeSize();

	int maternalStartingChromosome = pRandom->Bernoulli(0.5);
	int paternalStartingChromosome = pRandom->Bernoulli(0.5);

	set<unsigned int> maternalRecomPositions;
	set<unsigned int> paternalRecomPositions;

	// Determine at which chromosome ends does the genome break
	for (int pos : chromosomeEnds) {
		if (pRandom->Bernoulli(0.5))
			maternalRecomPositions.insert(pos);
		if (pRandom->Bernoulli(0.5))
			paternalRecomPositions.insert(pos);
	}

	// Draw recombination events for maternal genome
	if (pSpecies->getRecombinationRate() > 0.0)
		events = pRandom->Poisson(genomeSize * pSpecies->getRecombinationRate());
	int nbrCrossOvers = events + maternalRecomPositions.size();
	while (maternalRecomPositions.size() < nbrCrossOvers) {
		// Sample recombination sites
		maternalRecomPositions.insert(pRandom->IRandom(0, genomeSize));
	}

	// Draw recombination events for paternal genome
	if (pSpecies->getRecombinationRate() > 0.0)
		events = pRandom->Poisson(genomeSize * pSpecies->getRecombinationRate());
	nbrCrossOvers = events + paternalRecomPositions.size();
	while (paternalRecomPositions.size() < nbrCrossOvers) {
		paternalRecomPositions.insert(pRandom->IRandom(0, genomeSize));
	}

	// End of genome always recombines
	maternalRecomPositions.insert(genomeSize - 1);
	paternalRecomPositions.insert(genomeSize - 1);

	const auto& spTraits = pSpecies->getTraitTypes();

	for (auto const& trait : spTraits)
	{
		const auto motherTrait = mother->getTrait(trait);
		const auto fatherTrait = father->getTrait(trait);
		auto newTrait = motherTrait->clone(); // shallow copy, pointer to proto trait initialised and empty sequence

		newTrait->inherit(motherTrait, maternalRecomPositions, FEM, maternalStartingChromosome);
		if (newTrait->isInherited()) {
			newTrait->inherit(fatherTrait, paternalRecomPositions, MAL, paternalStartingChromosome);
			if (newTrait->getMutationRate() > 0 && pSpecies->areMutationsOn())
				newTrait->mutate();
		}
		if (trait == ADAPTIVE1 || trait == ADAPTIVE2 || trait == ADAPTIVE3 || trait == ADAPTIVE4 || trait == ADAPTIVE5)
			fitness *= newTrait->express();

		traitTable.insert(make_pair(trait, move(newTrait)));
	}
}

void Individual::inherit(Species* pSpecies, const Individual* mother) {
	set<unsigned int> recomPositions; //not used here cos haploid but need it for inherit function, not ideal 
	int startingChromosome = 0;
	//const auto mumTraitTable = mother->getTraitTable(); //assuming mother and father share the same genetic structure..

	const auto& mumTraits = getTraitTypes();

	for (auto const& trait : mumTraits)
	{
		const auto motherTrait = mother->getTrait(trait);

		auto newTrait = motherTrait->clone(); //shallow copy, pointer to proto trait initialised and empty sequence

		newTrait->inherit(motherTrait, recomPositions, FEM, startingChromosome);
		if (newTrait->isInherited()) {
			if (newTrait->getMutationRate() > 0 && pSpecies->areMutationsOn())
				newTrait->mutate();
		}
		if (trait == ADAPTIVE1 || trait == ADAPTIVE2 || trait == ADAPTIVE3 || trait == ADAPTIVE4 || trait == ADAPTIVE5)
			fitness *= newTrait->express();

		traitTable.insert(make_pair(trait, move(newTrait)));
	}
}

// Set genes for individual variation from species initialisation parameters
void Individual::setUpGenes(Species* pSpecies, int resol) {

	// this way to keep spp trait table immutable i.e. not able to call getTraitTable, 
	// could pass it back by value (copy) instead but could be heavy if large map
	const auto& speciesTraits = pSpecies->getTraitTypes();

	for (auto const& trait : speciesTraits)
	{
		const auto spTrait = pSpecies->getTrait(trait);
		this->traitTable.emplace(trait, traitFactory.Create(trait, spTrait));
	}
	setQTLPhenotypes(pSpecies, resol);
}

void Individual::setQTLPhenotypes(Species* pSpecies, int resol) {

	const emigRules emig = pSpecies->getEmig();
	const trfrRules trfr = pSpecies->getTrfr();
	const settleType sett = pSpecies->getSettle();

	// record phenotypic traits
	if (emig.indVar)
		this->setEmigTraits(pSpecies, emig.sexDep, emig.densDep);
	if (trfr.indVar)
		this->setTransferTraits(pSpecies, trfr, resol);
	if (sett.indVar)
		this->setSettlementTraits(pSpecies, sett.sexDep);
}

void Individual::setTransferTraits(Species* pSpecies, trfrRules trfr, int resol) {
	if (trfr.moveModel) {
		if (trfr.moveType == 1) {
			setSMSTraits(pSpecies);
		}
		else
			setCRWTraits(pSpecies, trfr.sexDep);
	}
	else
		setKernelTraits(pSpecies, trfr.sexDep, trfr.twinKern, resol);
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

#if RSDEBUG
	//DEBUGLOG << "Individual::setSettTraits(): indId=" << indId
	//	<< " s.s0=" << s.s0 << " s.alpha=" << s.alpha << " s.beta=" << s.beta
	//	<< endl;
#endif
#if RSDEBUG
		//DEBUGLOG << "Individual::setSettTraits(): indId=" << indId
		//	<< " sparams.s0Mean=" << sparams.s0Mean << " sparams.s0SD=" << sparams.s0SD 
		//	<< " sparams.s0Scale=" << sparams.s0Scale
		//	<< endl;
#endif
	pSettleTraits = make_unique<settleTraits>();
	pSettleTraits->s0 = (float)(s.s0);
	pSettleTraits->alpha = (float)(s.alpha);
	pSettleTraits->beta = (float)(s.beta);
#if RSDEBUG
	//DEBUGLOG << "Individual::setSettTraits(): indId=" << indId
	//	<< " setttraits->s0=" << setttraits->s0
	//	<< " setttraits->alpha=" << setttraits->alpha << " setttraits->beta=" << setttraits->beta
	//	<< endl;
#endif
	if (pSettleTraits->s0 < 0.0) pSettleTraits->s0 = 0.0;
	if (pSettleTraits->s0 > 1.0) pSettleTraits->s0 = 1.0;
#if RSDEBUG
	//DEBUGLOG << "Individual::setSettTraits(): indId=" << indId
	//	<< " setttraits->s0=" << setttraits->s0
	//	<< " setttraits->alpha=" << setttraits->alpha << " setttraits->beta=" << setttraits->beta
	//	<< endl;
#endif
	return;
}


// Inherit genome from parent(s)
void Individual::inheritTraits(Species* pSpecies, Individual* mother, Individual* father, int resol)
{
	inherit(pSpecies, mother, father);
	setQTLPhenotypes(pSpecies, resol);
}

// Inherit genome from mother, haploid
void Individual::inheritTraits(Species* pSpecies, Individual* mother, int resol)
{
	inherit(pSpecies, mother);
	setQTLPhenotypes(pSpecies, resol);
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

float Individual::getFitness(void) { return fitness; }

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
#if RSDEBUG
	//DEBUGLOG << "Individual::setYearSteps(): indId=" << indId
	//	<< " t=" << t << " path->year=" << path->year
	//	<< endl;
#endif
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
	if (sexDep && this->getSex() == MAL) {
		e.d0 = this->getTrait(E_D0_M)->express();
		if (densityDep) {
			e.alpha = getTrait(E_ALPHA_M)->express();
			e.beta = getTrait(E_BETA_M)->express();
		}
	}
	else {
		e.d0 = this->getTrait(E_D0_F)->express();
		if (densityDep) {
			e.alpha = getTrait(E_ALPHA_F)->express();
			e.beta = getTrait(E_BETA_F)->express();
		}
	}

#if RSDEBUG
	//DEBUGLOG << "Individual::setEmigTraits(): indId=" << indId
	//	<< " eparams.betaMean=" << eparams.betaMean << " eparams.betaSD=" << eparams.betaSD 
	//	<< " eparams.betaScale=" << eparams.betaScale
	//	<< endl;
#endif
	pEmigTraits = make_unique<emigTraits>();
	pEmigTraits->d0 = (float)(e.d0);
	pEmigTraits->alpha = (float)(e.alpha);
	pEmigTraits->beta = (float)(e.beta);
#if RSDEBUG
	//DEBUGLOG << "Individual::setEmigTraits(): indId=" << indId
	//	<< " emigtraits->d0=" << emigtraits->d0
	//	<< " emigtraits->alpha=" << emigtraits->alpha << " emigtraits->beta=" << emigtraits->beta
	//	<< endl;
#endif
	if (pEmigTraits->d0 < 0.0) pEmigTraits->d0 = 0.0;
	if (pEmigTraits->d0 > 1.0) pEmigTraits->d0 = 1.0;
#if RSDEBUG
	//DEBUGLOG << "Individual::setEmigTraits(): indId=" << indId
	//	<< " emigtraits->d0=" << emigtraits->d0
	//	<< " emigtraits->alpha=" << emigtraits->alpha << " emigtraits->beta=" << emigtraits->beta
	//	<< endl;
#endif
	return;
}

// Get phenotypic emigration traits
emigTraits Individual::getEmigTraits(void) {
#if RSDEBUG
	//DEBUGLOG << "Individual::getEmigTraits(): indId=" << indId
	//	<< endl;
#endif
	emigTraits e; e.d0 = e.alpha = e.beta = 0.0;
	if (pEmigTraits != 0) {
		e.d0 = pEmigTraits->d0;
		e.alpha = pEmigTraits->alpha;
		e.beta = pEmigTraits->beta;
	}
#if RSDEBUG
	//DEBUGLOG << "Individual::getEmigTraits(): indId=" << indId
	//	<< " e.d0=" << e.d0 << " e.alpha=" << e.alpha << " e.beta=" << e.beta
	//	<< endl;
#endif

	return e;
}
// Set phenotypic transfer by kernel traits
void Individual::setKernelTraits(Species* pSpecies, bool sexDep, bool twinKernel, int resol) {

	trfrKernTraits k; k.meanDist1 = k.meanDist2 = k.probKern1 = 0.0;
	if (sexDep && this->sex == MAL) {
		k.meanDist1 = getTrait(KERNEL_MEANDIST_1_M)->express();

		if (twinKernel) { // twin kernel
			k.meanDist2 = getTrait(KERNEL_MEANDIST_2_M)->express();
			k.probKern1 = getTrait(KERNEL_PROBABILITY_M)->express();
		}

	}
	else {
		k.meanDist1 = getTrait(KERNEL_MEANDIST_1_F)->express();

		if (twinKernel) { // twin kernel
			k.meanDist2 = getTrait(KERNEL_MEANDIST_2_F)->express();
			k.probKern1 = getTrait(KERNEL_PROBABILITY_F)->express();
		}

	}
	float meanDist1 = (float)(k.meanDist1);
	float meanDist2 = (float)(k.meanDist2);
	float probKern1 = (float)(k.probKern1);

#if RSDEBUG
	//DEBUGLOG << "Individual::setKernTraits(): indId=" << indId
	//	<< " kerntraits->meanDist1=" << kerntraits->meanDist1
	//	<< " kerntraits->meanDist2=" << kerntraits->meanDist2
	//	<< " kerntraits->probKern1=" << kerntraits->probKern1
	//	<< endl;
#endif
	if (!pSpecies->useFullKernel()) {
		// kernel mean(s) may not be less than landscape resolution
		if (meanDist1 < resol) meanDist1 = (float)resol;
		if (meanDist2 < resol) meanDist2 = (float)resol;
	}
	if (probKern1 < 0.0) probKern1 = 0.0;
	if (probKern1 > 1.0) probKern1 = 1.0;
#if RSDEBUG
	//DEBUGLOG << "Individual::setKernTraits(): indId=" << indId
	//	<< " kerntraits->meanDist1=" << kerntraits->meanDist1
	//	<< " kerntraits->meanDist2=" << kerntraits->meanDist2
	//	<< " kerntraits->probKern1=" << kerntraits->probKern1
	//	<< endl;
#endif
	auto& pKernel = dynamic_cast<kernelData&>(*pTrfrData);
	pKernel.meanDist1 = meanDist1;
	pKernel.meanDist2 = meanDist2;
	pKernel.probKern1 = probKern1;

	return;
}



// Get phenotypic emigration traits
trfrKernTraits Individual::getKernTraits(void) {
#if RSDEBUG
	//DEBUGLOG << "Individual::getKernTraits(): indId=" << indId
	//	<< endl;
#endif
	trfrKernTraits k; k.meanDist1 = k.meanDist2 = k.probKern1 = 0.0;
	if (pTrfrData != 0) {

		auto& pKernel = dynamic_cast<const kernelData&>(*pTrfrData);

		k.meanDist1 = pKernel.meanDist1;
		k.meanDist2 = pKernel.meanDist2;
		k.probKern1 = pKernel.probKern1;
	}
#if RSDEBUG
	//DEBUGLOG << "Individual::getKernTraits(): indId=" << indId
	//	<< " k.meanDist1=" << k.meanDist1 << " k.meanDist2=" << k.meanDist1
	//	<< " k.probKern1=" << k.probKern1
	//	<< endl;
#endif

	return k;
}

void Individual::setSMSTraits(Species* pSpecies) {

	trfrSMSTraits s = pSpecies->getSMSTraits();

	double dp, gb, alphaDB, betaDB;
	dp = gb = alphaDB = betaDB = 0.0;
	dp = getTrait(SMS_DP)->express();
	gb = getTrait(SMS_GB)->express();
	if (s.goalType == 2) {
		alphaDB = getTrait(SMS_ALPHADB)->express();
		betaDB = getTrait(SMS_BETADB)->express();
	}

#if RSDEBUG
	//DEBUGLOG << "Individual::setSMSTraits(): indId=" << indId
	//	<< " dp=" << dp << " gb=" << gb
	//	<< " alphaDB=" << alphaDB << " betaDB=" << betaDB
	//	<< endl;
#endif
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
#if RSDEBUG
	//DEBUGLOG << "Individual::setSMSTraits() 1111: indId=" << indId
	//	<< " smsData->dp=" << smsData->dp	<< " smsData->gb=" << smsData->gb
	//	<< " smsData->alphaDB=" << smsData->alphaDB	<< " smsData->betaDB=" << smsData->betaDB
	//	<< endl;
#endif
	if (pSMS.dp < 1.0) pSMS.dp = 1.0;
	if (pSMS.gb < 1.0) pSMS.gb = 1.0;
	if (pSMS.alphaDB <= 0.0) pSMS.alphaDB = 0.000001f;
	if (pSMS.betaDB < 1) pSMS.betaDB = 1;
#if RSDEBUG
	//DEBUGLOG << "Individual::setSMSTraits() 2222: indId=" << indId
	//	<< " smsData->dp=" << smsData->dp	<< " smsData->gb=" << smsData->gb
	//	<< " smsData->alphaDB=" << smsData->alphaDB	<< " smsData->betaDB=" << smsData->betaDB
	//	<< endl;
#endif
	return;
}

trfrData* Individual::getTrfrData(void) {
	return pTrfrData.get();
}

// Get phenotypic transfer by SMS traits
trfrSMSTraits Individual::getSMSTraits(void) {
#if RSDEBUG
	//DEBUGLOG << "Individual::getSMSTraits(): indId=" << indId << " smsData=" << smsData
	//	<< endl;
#endif
	trfrSMSTraits s; s.dp = s.gb = s.alphaDB = 1.0; s.betaDB = 1;
	if (pTrfrData != 0) {

		auto& pSMS = dynamic_cast<const smsData&>(*pTrfrData);

		s.dp = pSMS.dp; s.gb = pSMS.gb;
		s.alphaDB = pSMS.alphaDB; s.betaDB = pSMS.betaDB;
	}
#if RSDEBUG
	//DEBUGLOG << "Individual::getSMSTraits(): indId=" << indId
	//	<< " s.dp=" << s.dp << " s.gb=" << s.gb
	//	<< " s.alphaDB=" << s.alphaDB << " s.betaDB=" << s.betaDB
	//	<< endl;
#endif
	return s;
}


// Set phenotypic transfer by CRW traits
void Individual::setCRWTraits(Species* pSpecies, bool sexDep) {
	trfrCRWTraits c; c.stepLength = c.rho = 0.0;
	if (sexDep && this->sex == MAL) {
		c.stepLength = getTrait(CRW_STEPLENGTH_M)->express();
		c.rho = getTrait(CRW_STEPCORRELATION_M)->express();
	}
	else {
		c.stepLength = getTrait(CRW_STEPLENGTH_F)->express();
		c.rho = getTrait(CRW_STEPCORRELATION_F)->express();
	}

#if RSDEBUG
	//DEBUGLOG << "Individual::setCRWTraits(): indId=" << indId
	//	<< " c.stepLength=" << c.stepLength << " c.rho=" << c.rho
	//	<< endl;
#endif

	auto& pCRW = dynamic_cast<crwData&>(*pTrfrData);
	pCRW.stepLength = (float)(c.stepLength);
	pCRW.rho = (float)(c.rho);
#if RSDEBUG
	//DEBUGLOG << "Individual::setCRWTraits(): indId=" << indId
	//	<< " crw->stepL=" << crw->stepL	<< " crw->rho=" << crw->rho
	//	<< endl;
#endif
	if (pCRW.stepLength < 1.0) pCRW.stepLength = 1.0;
	if (pCRW.rho < 0.0) pCRW.rho = 0.0;
	if (pCRW.rho > 0.999) pCRW.rho = 0.999f;
#if RSDEBUG
	//DEBUGLOG << "Individual::setCRWTraits(): indId=" << indId
	//	<< " crw->stepL=" << crw->stepL	<< " crw->rho=" << crw->rho
	//	<< endl;
#endif
	return;
}

// Get phenotypic transfer by CRW traits
trfrCRWTraits Individual::getCRWTraits(void) {
#if RSDEBUG
	//DEBUGLOG << "Individual::getCRWTraits(): indId=" << indId
	//	<< endl;
#endif
	trfrCRWTraits c; c.stepLength = c.rho = 0.0;


	if (pTrfrData != 0) {

		auto& pCRW = dynamic_cast<const crwData&>(*pTrfrData);

		c.stepLength = pCRW.stepLength;
		c.rho = pCRW.rho;
	}
#if RSDEBUG
	//DEBUGLOG << "Individual::getCRWTraits(): indId=" << indId
	//	<< " c.stepLength=" << c.stepLength << " c.rho=" << c.rho
	//	<< endl;
#endif

	return c;

}

// Get phenotypic settlement traits
settleTraits Individual::getSettTraits(void) {
#if RSDEBUG
	//DEBUGLOG << "Individual::getSettTraits(): indId=" << indId
	//	<< endl;
#endif
	settleTraits s; s.s0 = s.alpha = s.beta = 0.0;
	if (pSettleTraits != 0) {
		s.s0 = pSettleTraits->s0;
		s.alpha = pSettleTraits->alpha;
		s.beta = pSettleTraits->beta;
	}
#if RSDEBUG
	//DEBUGLOG << "Individual::getSettTraits(): indId=" << indId
	//	<< " s.s0=" << s.s0 << " s.alpha=" << s.alpha << " s.beta=" << s.beta
	//	<< endl;
#endif

	return s;
}

/*
locus Individual::getAlleles(int g) {
locus l; l.allele[0] = l.allele[1] = 0.0;
if (pGenome != 0) l = pGenome->getAlleles(g);
return l;
}
*/

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
int Individual::moveKernel(Landscape* pLandscape, Species* pSpecies,
	const short repType, const bool absorbing)
{

	intptr patch;
	int patchNum = 0;
	int newX = 0, newY = 0;
	int dispersing = 1;
	double xrand, yrand, meandist, dist, r1, rndangle, nx, ny;
	float localK;
	trfrKernTraits kern;
	Cell* pCell;
	Patch* pPatch;
	locn loc = pCurrCell->getLocn();

	landData land = pLandscape->getLandData();

	bool usefullkernel = pSpecies->useFullKernel();
	trfrRules trfr = pSpecies->getTrfr();
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
				kern = pSpecies->getKernTraits(stage, sex);
			}
			else {
				kern = pSpecies->getKernTraits(0, sex);
			}
		}
		else {
			if (trfr.stgDep) {
				kern = pSpecies->getKernTraits(stage, 0);
			}
			else {
				kern = pSpecies->getKernTraits(0, 0);
			}
		}
	}
#if RSDEBUG
	//Patch *startPatch = (Patch*)startpatch;
	//DEBUGLOG << "Individual::moveKernel(): indId=" << indId << " x=" << loc.x << " y=" << loc.y
	////	<< " natalPatch = " << natalPatch
	////	<< " startpatch = " << startpatch << " patchNum = " << startPatch->getPatchNum()
	//	<< " kern.meanDist1=" << kern.meanDist1;
	//if (trfr.twinKern) {
	//	DEBUGLOG << " meanDist2=" << kern.meanDist2 << " probKern1=" << kern.probKern1;
	//}
	//DEBUGLOG << endl;
#endif

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
#if RSDEBUG
	//DEBUGLOG << "Individual::moveKernel(): indId=" << indId << " meandist=" << meandist << endl;
#endif
// scaled mean may not be less than 1 unless emigration derives from the kernel
// (i.e. the 'use full kernel' option is applied)
	if (!usefullkernel && meandist < 1.0) meandist = 1.0;
#if RSDEBUG
	//DEBUGLOG << "Individual::moveKernel(): indId=" << indId << " meandist=" << meandist << endl;
#endif

#if RSDEBUG
//Patch *startPatch = (Patch*)startpatch;
//DEBUGLOG << "Individual::moveKernel(): indId = " << indId << " x = " << x << " y = " << y
//	<< " natalPatch = " << natalPatch
////	<< " startpatch = " << startpatch << " patchNum = " << startPatch->getPatchNum()
//	<< " meanDist1 = " << kern.meanDist1;
//if (trfr.twinKern) {
//	DEBUGLOG << " probKern1 = " << kern.probKern1 << " meanDist2 = " << kern.meanDist2;
//}
//DEBUGLOG << " meandist = " << meandist << endl;
#endif

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
				xrand = (double)loc.x + pRandom->Random() * 0.999;
				yrand = (double)loc.y + pRandom->Random() * 0.999;

				r1 = 0.0000001 + pRandom->Random() * (1.0 - 0.0000001);
				//			dist = (-1.0*meandist)*std::log(r1);
				dist = (-1.0 * meandist) * log(r1);  // for LINUX_CLUSTER

				rndangle = pRandom->Random() * 2.0 * PI;
				nx = xrand + dist * sin(rndangle);
				ny = yrand + dist * cos(rndangle);
				if (nx < 0.0) newX = -1; else newX = (int)nx;
				if (ny < 0.0) newY = -1; else newY = (int)ny;
#if RSDEBUG
				if (path != 0) (path->year)++;
#endif
				loopsteps++;
#if RSDEBUG
				//DEBUGLOG << "Individual::moveKernel(): indId=" << indId << " status=" << status
				//	<< " loopsteps=" << loopsteps << " newX=" << newX << " newY=" << newY
				//	<< " loc.x=" << loc.x << " loc.y=" << loc.y
				//	<< endl;
#endif
			} while (loopsteps < 1000 &&
				((!absorbing && (newX < land.minX || newX > land.maxX
					|| newY < land.minY || newY > land.maxY))
					|| (!usefullkernel && newX == loc.x && newY == loc.y))
				);
			if (loopsteps < 1000) {
				if (newX < land.minX || newX > land.maxX
					|| newY < land.minY || newY > land.maxY) { // beyond absorbing boundary
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
			else {
				patch = 0;
				patchNum = -1;
			}
#if RSDEBUG
			//DEBUGLOG << "Individual::moveKernel(): indId=" << indId << " status=" << status
			//	<< " loopsteps=" << loopsteps << " newX=" << newX << " newY=" << newY
			//	<< " pCell=" << pCell << " patch=" << patch << " patchNum=" << patchNum
			//	<< endl;
#endif
		} while (!absorbing && patchNum < 0 && loopsteps < 1000); 			 // in a no-data region
	} while (!usefullkernel && pPatch == pNatalPatch && loopsteps < 1000); 	// still in the original (natal) patch

	if (loopsteps < 1000) {
		if (pCell == 0) { // beyond absorbing boundary or in no-data cell
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
				else
					status = 6; // dies (unless there is a suitable neighbouring cell)
			}
		}
	}
	else {
		status = 6;
		dispersing = 0;
	}
#if RSDEBUG
	//DEBUGLOG << "Individual::moveKernel(): indId=" << indId
	//	<< " newX=" << newX << " newY=" << newY
	//	<< " patch=" << patch
	//	<< " patchNum=" << patchNum << " status=" << status;
	//DEBUGLOG << endl;
#endif

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

	trfrRules trfr = pSpecies->getTrfr();
	trfrCRWTraits movt = pSpecies->getCRWTraits();
	settleSteps settsteps = pSpecies->getSteps(stage, sex);

	patch = pCurrCell->getPatch();
#if RSDEBUG
	//DEBUGLOG << "Individual::moveStep() AAAA: indId=" << indId
	//	<< " pCurrCell=" << pCurrCell << " patch=" << patch
	//	<< endl;
#endif

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
#if RSDEBUG
		//locn temploc = pCurrCell->getLocn();
		//DEBUGLOG << "Individual::moveStep(): x=" << temploc.x << " y=" << temploc.x
		//	<< " landIx=" << landIx << " h=" << h << " mortprob=" << mortprob
		//	<< endl;
#endif
	}
	else mortprob = movt.stepMort;
	// ... unless individual has not yet left natal patch in emigration year
	if (pPatch == pNatalPatch && path->out == 0 && path->year == path->total) {
		mortprob = 0.0;
	}
#if RSDEBUG
	locn loc0, loc1, loc2;
	//loc0 = pCurrCell->getLocn();
	//DEBUGLOG << "Individual::moveStep() BBBB: indId=" << indId << " status=" << status
	//	<< " path->year=" << path->year << " path->out=" << path->out
	//	<< " settleStatus=" << path->settleStatus
	//	<< " x=" << loc0.x << " y=" << loc0.y
	////	<< " patch=" << patch
	//	<< " pPatch=" << pPatch
	//	<< " patchNum=" << patchNum;
	////	<< " natalPatch=" << natalPatch;
	////if (crw != 0) {
	////	DEBUGLOG << " xc=" << crw->xc << " yc=" << crw->yc;
	////	DEBUGLOG << " rho=" << movt.rho << " stepLength=" << movt.stepLength;
	////}
	//DEBUGLOG << endl;
#endif
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
#if RSDEBUG                   
			//loc1 = pCurrCell->getLocn();      
			//DEBUGLOG << "Individual::moveStep() FFFF: indId=" << indId << " status=" << status
			////	<< " path->year=" << path->year
			//	<< " path->season=" << path->season
			//	<< " x=" << loc1.x << " y=" << loc1.y
			//	<< " smsData->goalType=" << smsData->goalType
			//	<< " goal.x=" << smsData->goal.x
			//	<< " goal.y=" << smsData->goal.y
			//	<< endl;
#endif
			move = smsMove(pLandscape, pSpecies, landIx, pPatch == pNatalPatch, trfr.indVar, absorbing);
#if RSDEBUG
			//DEBUGLOG << "Individual::moveStep() GGGG: indId=" << indId << " status=" << status
			//	<< " move.dist=" << move.dist
			//	<< endl;
#endif
			if (move.dist < 0.0) {
				// either INTERNAL ERROR CONDITION - INDIVIDUAL IS IN NO-DATA SQUARE
				// or individual has crossed absorbing boundary ...
				// ... individual dies
				status = 6;
				dispersing = 0;
			}
			else {
#if RSDEBUG
				//loc1 = pCurrCell->getLocn();
				//DEBUGLOG << "Individual::moveStep() HHHH: indId=" << indId << " status=" << status
				//	<< " path->year=" << path->year
				//	<< " x=" << loc1.x << " y=" << loc1.y
				////	<< " smsData = " << smsData
				//	<< endl;
#endif

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

			steplen = movt.stepLength; if (steplen < 0.2 * land.resol) steplen = 0.2 * land.resol;
			rho = movt.rho; if (rho > 0.99) rho = 0.99;
			if (pPatch == pNatalPatch) {
				rho = 0.99; // to promote leaving natal patch
				path->out = 0;
			}
			if (movt.straigtenPath && path->settleStatus > 0) {
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
#if RSDEBUG
					//DEBUGLOG << "Individual::moveStep(): indId=" << indId
					//	<< " xc=" << crw->xc << " yc=" << crw->yc << " pCurrCell=" << pCurrCell
					//	<< " steps=" << path->year << " loopsteps=" << loopsteps
					//	<< " steplen=" << steplen << " rho=" << rho << " angle=" << angle
					//	<< " xcnew=" << xcnew << " ycnew=" << ycnew << " newX=" << newX << " newY=" << newY << endl;
#endif
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
#if RSDEBUG
				//DEBUGLOG << "Individual::moveStep(): indId=" << indId
				//	<< " loopsteps=" << loopsteps << " absorbed=" << absorbed
				//	<< " pCurrCell=" << pCurrCell << " patch=" << patch << endl;
#endif
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
#if RSDEBUG
			//DEBUGLOG << "Individual::moveStep(): indId=" << indId
			//	<< " status=" << status
			//	<< " pCurrCell=" << pCurrCell << " patch=" << patch << endl;
#endif
			break;

		} // end of switch (trfr.moveType)

#if RSDEBUG
//locn loc2;
//if (pCurrCell > 0) {
//	loc2 = pCurrCell->getLocn();
//}
//else {
//	loc2.x = -9999; loc2.y = -9999;
//}
//DEBUGLOG << "Individual::moveStep() ZZZZ: indId=" << indId
//	<< " status=" << status
//	<< " path->total=" << path->total
//	<< " x=" << loc2.x << " y=" << loc2.y
//	<< " patch=" << patch;
//if (patch > 0) {
//	pPatch = (Patch*)patch;
//	DEBUGLOG << " patchNum=" << pPatch->getPatchNum()
//		<< " getK()=" << pPatch->getK()
//		<< " popn=" << pPatch->getPopn((int)pSpecies);
//}
//	DEBUGLOG << endl;
#endif
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
	//if (write_out) {
	//	out<<endl<<"ind = "<<get_id()<<" pr = "<<get_pr()<<" step = "<<step;
	//	if (step == 0) out<<" start at "<<x<<" "<<y;
	//	out<<endl;
	//}
	//pCell = pLand->findCell(x,y);
	if (pCurrCell == 0)
	{
		// x,y is a NODATA square - this should not occur here
		// return a negative distance to indicate an error
		move.dist = -69.0; move.cost = 0.0;
		return move;
	}

#if RSDEBUG
	//DEBUGLOG << "Individual::smsMove(): this=" << this << endl;
#endif

	landData land = pLand->getLandData();
	trfrSMSTraits movt = pSpecies->getSMSTraits();
	current = pCurrCell->getLocn();

	//get weights for directional persistence....
	//if ((path->out > 0 && path->out < 10 && path->out < 2*movt.pr)
	if ((path->out > 0 && path->out <= (movt.pr + 1))
		|| natalPatch
		|| (movt.straigtenPath && path->settleStatus > 0)) {
		// inflate directional persistence to promote leaving the patch
		if (indvar) nbr = getSimDir(current.x, current.y, 10.0f * pSMS.dp);
		else nbr = getSimDir(current.x, current.y, 10.0f * movt.dp);
	}
	else {
		if (indvar) nbr = getSimDir(current.x, current.y, pSMS.dp);
		else nbr = getSimDir(current.x, current.y, movt.dp);
	}
	if (natalPatch || path->settleStatus > 0) path->out = 0;
	//if (natalPatch) path->out = 0;
#if RSDEBUG
//DEBUGLOG << "Individual::smsMove() 0000: nbr matrix" << endl;
//for (y2 = 2; y2 > -1; y2--) {
//	for (x2 = 0; x2 < 3; x2++) DEBUGLOG << nbr.cell[x2][y2] << " ";
//	DEBUGLOG << endl;
//}
#endif
//if (write_out) {
//	out<<endl<<"directional persistence weights:"<<endl;
//	for (y2 = 2; y2 > -1; y2--) {
//		for (x2 = 0; x2 < 3; x2++) out<<nbr.cell[x2][y2]<<" ";
//		out<<endl;
//	}
//}
//if (write_out2)
//	out2<<endl<<"ind = "<<get_id()<<" pr = "<<get_pr()<<" step = "<<step<<endl<<endl;

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
#if RSDEBUG
			//DEBUGLOG << "Individual::smsMove(): exp_arg=" << exp_arg;
#endif
			if (exp_arg > 100.0) exp_arg = 100.0; // to prevent exp() overflow error
			gb = 1.0 + (pSMS.gb - 1.0) / (1.0 + exp(exp_arg));
		}
		else {
			double exp_arg = -((double)nsteps - (double)movt.betaDB) * (-movt.alphaDB);
#if RSDEBUG
			//DEBUGLOG << "Individual::smsMove(): exp_arg=" << exp_arg;
#endif
			if (exp_arg > 100.0) exp_arg = 100.0; // to prevent exp() overflow error
			gb = 1.0 + (movt.gb - 1.0) / (1.0 + exp(exp_arg));
		}
	}
	else gb = movt.gb;
	goal = getGoalBias(current.x, current.y, movt.goalType, (float)gb);
	//if (write_out) {
	//	out<<"goal bias weights:"<<endl;
	//	for (y2 = 2; y2 > -1; y2--) {
	//		for (x2 = 0; x2 < 3; x2++) out<<goal.cell[x2][y2]<<" ";
	//		out<<endl;
	//	}
	//}
	//if (write_out2)
	//  out2<<endl<<"ind = "<<get_id()<<" pr = "<<get_pr()<<" step = "<<step<<endl<<endl;

	// get habitat-dependent weights (mean effective costs, given perceptual range)
	// first check if costs have already been calculated

	hab = pCurrCell->getEffCosts();
#if RSDEBUG
	//if (hab.cell[0][0] >= 0) {
	//	DEBUGLOG << "Individual::smsMove() 1111: x=" << current.x << " y=" << current.y << endl;
	//	for (y2 = 2; y2 > -1; y2--) {
	//		for (x2 = 0; x2 < 3; x2++) DEBUGLOG << hab.cell[x2][y2] << " ";
	//		DEBUGLOG << endl;
	//	}
	//}
#endif
//if (write_out) {
//  out<<"stored effective costs:"<<endl;
//	for (y2 = 2; y2 > -1; y2--) {
//    for (x2 = 0; x2 < 3; x2++) out<<hab.cell[x2][y2]<<" ";
//	  out<<endl;
//	}
//}
	if (hab.cell[0][0] < 0.0) { // costs have not already been calculated
		hab = getHabMatrix(pLand, pSpecies, current.x, current.y, movt.pr, movt.prMethod,
			landIx, absorbing);
#if RSDEBUG
		//DEBUGLOG << "Individual::smsMove() 2222: " << endl;
		//for (y2 = 2; y2 > -1; y2--) {
		//	for (x2 = 0; x2 < 3; x2++) DEBUGLOG << hab.cell[x2][y2] << " ";
		//	DEBUGLOG << endl;
		//}
#endif
		pCurrCell->setEffCosts(hab);
	}
	else { // they have already been calculated - no action required
		//	if (write_out) {
		//		out<<"*** using previous effective costs ***"<<endl;
		//	}
	}
	//if (write_out) {
	//	out<<"mean effective costs:"<<endl;
	//	for (y2 = 2; y2 > -1; y2--) {
	//		for (x2 = 0; x2 < 3; x2++) {
	//			out<<hab.cell[x2][y2]<<" ";
	//		}
	//		out<<endl;
	//	}
	//	out<<"weighted effective costs:"<<endl;
	//}

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
			//		if (write_out) {
			//			out<<nbr.cell[x2][y2]<<" "; if (x2==1 && y2==1) out<<"         ";
			//		}
		}
		//	if (write_out) out<<endl;
	}
#if RSDEBUG
	//DEBUGLOG << "Individual::smsMove() 3333: " << endl;
	//for (y2 = 2; y2 > -1; y2--) {
	//	for (x2 = 0; x2 < 3; x2++) DEBUGLOG << nbr.cell[x2][y2] << " ";
	//	DEBUGLOG << endl;
	//}
#endif

// determine reciprocal of effective cost for the 8 neighbours
//if (write_out) out<<"reciprocal weighted effective costs:"<<endl;
	for (y2 = 2; y2 > -1; y2--) {
		for (x2 = 0; x2 < 3; x2++) {
			if (nbr.cell[x2][y2] > 0.0) nbr.cell[x2][y2] = 1.0f / nbr.cell[x2][y2];
			//		if (write_out) {
			//			out<<nbr.cell[x2][y2]<<" "; if (x2==1 && y2==1) out<<"         ";
			//		}
		}
		//	if (write_out) out<<endl;
	}

	// set any cells beyond the current landscape limits and any no-data cells
	// to have zero probability
	// increment total for re-scaling to sum to unity

#if RSDEBUG
//array3x3d temp;
//for (y2 = 2; y2 > -1; y2--) {
//	for (x2 = 0; x2 < 3; x2++) {
//		temp.cell[x2][y2] = nbr.cell[x2][y2];
//		if (current.x == 488 && current.y == 422) {
//			pCell = pLand->findCell((current.x+x2-1),(current.y+y2-1));
//			DEBUGLOG << "Individual::smsMove(): this=" << this
//				<< " IN THE PROBLEM CELL"
//				<< " y=" << current.y << " x=" << current.x
//				<< " y2=" << y2 << " x2=" << x2
//				<< " pCell=" << pCell;
//			if (pCell != 0) DEBUGLOG << " pCell->getCost=" << pCell->getCost();
//			DEBUGLOG << endl;
//		}
//	}
//}
#endif

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
#if RSDEBUG
			//DEBUGLOG << "Individual::smsMove(): this=" << this
			//	<< " y=" << current.y << " x=" << current.x
			//	<< " y2=" << y2 << " x2=" << x2
			//	<< " pCell=" << pCell
			//	<< endl;
#endif
//		if (write_out) {
//			out<<nbr.cell[x2][y2]<<" "; if (x2==1 && y2==1) out<<"      ";
//		}
			sum_nbrs += nbr.cell[x2][y2];
		}
		//	if (write_out) out<<endl;
	}

	// scale effective costs as probabilities summing to 1
	//if (write_out) out<<"probabilities:"<<endl;
	if (sum_nbrs > 0.0) { // should always be the case, but safest to check...
		for (y2 = 2; y2 > -1; y2--) {
			for (x2 = 0; x2 < 3; x2++) {
				nbr.cell[x2][y2] = nbr.cell[x2][y2] / (float)sum_nbrs;
				//		if (write_out) {
				//			out<<nbr.cell[x2][y2]<<" "; if (x2==1 && y2==1) out<<"      ";
				//		}
			}
			//	if (write_out) out<<endl;
		}
	}
#if RSDEBUG
	//DEBUGLOG << "Individual::smsMove() 4444: " << endl;
	//for (y2 = 2; y2 > -1; y2--) {
	//	for (x2 = 0; x2 < 3; x2++) DEBUGLOG << nbr.cell[x2][y2] << " ";
	//	DEBUGLOG << endl;
	//}
#endif

// set up cell selection probabilities
//if (write_out) out<<"rnd = "<<rnd<<endl;
	double cumulative[9];
	int j = 0;
	cumulative[0] = nbr.cell[0][0];
	for (y2 = 0; y2 < 3; y2++) {
		for (x2 = 0; x2 < 3; x2++) {
			if (j != 0) cumulative[j] = cumulative[j - 1] + nbr.cell[x2][y2];
			j++;
			//    if (write_out) out<<"dx = "<<x2-1<<" dy = "<<y2-1<<" sum_rnd = "<<sum_rnd<<endl;
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
#if RSDEBUG
					//DEBUGLOG << "Individual::smsMove() 7777: rnd=" << rnd
					//	<< " j=" << j	<< " cumulative[j]=" << cumulative[j]
					//	<< endl;
#endif
					if (rnd < cumulative[j]) {
						newX = current.x + x2 - 1;
						newY = current.y + y2 - 1;
						if (x2 == 1 || y2 == 1) move.dist = (float)(land.resol);
						else move.dist = (float)(land.resol) * (float)SQRT2;
						//			if (write_out) {
						//				out<<"relative x and y "<<x2-1<<" "<<y2-1<<endl;
						//				out<<"cost of move "<<move.cost<<endl;
						//				out<<"move to: x = "<<x<<" y = "<<y;
						//				if (oob) out<<" ***** OUT OF BOUNDS *****";
						//				out<<endl;
						//			}
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
#if RSDEBUG
	//DEBUGLOG << "Individual::smsMove() 8888: pNewCell=" << pNewCell
	//	<< " loopsteps=" << loopsteps
	//	<< " current.x=" << current.x << " current.y=" << current.y
	//	<< " newX=" << newX << " newY=" << newY
	//	<< " land.minX=" << land.minX << " land.minY=" << land.minY
	//	<< " land.maxX=" << land.maxX << " land.maxY=" << land.maxY
	//	<< endl;
#endif
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
			//    if (write_out) out<<"step 3"<<endl;
			//		if (write_out) out<<"*** using last step only: x,y = "<<prev.x<<","<<prev.y<<endl;
			if ((x - prev.x) == 0 && (y - prev.y) == 0) { // STILL HAVE A PROBLEM!
				for (xx = 0; xx < 3; xx++) {
					for (yy = 0; yy < 3; yy++) {
						d.cell[xx][yy] = 1.0;
					}
				}
				//      if (write_out) out<<"step 4"<<endl;
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
	//    if (write_out) out<<"*** at goal: x,y = "<<goalx<<","<<goaly<<endl;
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
	//  if (write_out) out<<"goalx,goaly: "<<goalx<<","<<goaly<<" dx,dy: "<<dx<<","<<dy
	//    <<" theta: "<<theta<<endl;
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
	if (absorbing) nodatacost = ABSNODATACOST;
	else nodatacost = NODATACOST;

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
						//out<<"ROOK N-S: x2 = "<<x2<<" y2 = "<<y2<<endl;
						if (pr % 2 == 0) { xmin = -pr / 2; xmax = pr / 2; ymin = y2; ymax = y2 * pr; } // PR even
						else { xmin = -(pr - 1) / 2; xmax = (pr - 1) / 2; ymin = y2; ymax = y2 * pr; } // PR odd
					}
					if (y2 == 0) { // horizontal (E-W) move
						//out<<"ROOK E-W: x2 = "<<x2<<" y2 = "<<y2<<endl;
						if (pr % 2 == 0) { xmin = x2; xmax = x2 * pr; ymin = -pr / 2; ymax = pr / 2; } // PR even
						else { xmin = x2; xmax = x2 * pr; ymin = -(pr - 1) / 2; ymax = (pr - 1) / 2; } // PR odd
					}
				}
				else { // diagonal (bishop move)
					//out<<"BISHOP: x2 = "<<x2<<" y2 = "<<y2<<endl;
					xmin = x2; xmax = x2 * pr; ymin = y2; ymax = y2 * pr;
				}
			}
			//out<<"pre  swap: xmin = "<<xmin<<" ymin = "<<ymin<<" xmax = "<<xmax<<" ymax = "<<ymax<<endl;
			if (xmin > xmax) { int z = xmax; xmax = xmin; xmin = z; } // swap xmin and xmax
			if (ymin > ymax) { int z = ymax; ymax = ymin; ymin = z; } // swap ymin and ymax
			//out<<"post swap: xmin = "<<xmin<<" ymin = "<<ymin<<" xmax = "<<xmax<<" ymax = "<<ymax<<endl;
	//		if (write_out2) {
	//			out2<<"current x and y  "<<x<<" "<<y<<endl;
	//			out2<<"x2 and y2  "<<x2<<" "<<y2<<endl;
	//			out2<<"xmin,ymin "<<xmin<<","<<ymin<<" xmax,ymax "<<xmax<<","<<ymax<<endl;
	//		}

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
#if RSDEBUG
										//DEBUGLOG << "Individual::getHabMatrix(): x4=" << x4 << " y4=" << y4
										//	<< " landIx=" << landIx << " h=" << h << " cost=" << cost
										//	<< endl;
#endif
										pCell->setCost(cost);
									}
									else {
#if RSDEBUG
										//DEBUGLOG << "Individual::getHabMatrix(): x4=" << x4 << " y4=" << y4
										//	<< " cost=" << cost
										//	<< endl;
#endif

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
						//          if (write_out2) out2<<x+x3<<","<<y+y3<<","<<cost<<" wt. "<<weight<<endl;
						//#if GO2TARGET
						//					if (sq != 0) {
						//						if (sq->get_target() > 50) targetseen++;
						//					}
						//#endif
					} //end of y3 loop
				}  //end of x3 loop
	//		if (write_out) out<<"ncells in PR = "<<ncells<<" tot.wt. = "<<sumweights
	//      <<" w.cell = "<<w.cell[x2+1][y2+1]<<endl;
				if (ncells > 0) {
					if (prmethod == 1) w.cell[x2 + 1][y2 + 1] /= ncells; // arithmetic mean
					if (prmethod == 2) w.cell[x2 + 1][y2 + 1] = ncells / w.cell[x2 + 1][y2 + 1]; // hyperbolic mean
					if (prmethod == 3 && sumweights > 0)
						w.cell[x2 + 1][y2 + 1] /= (float)sumweights; // weighted arithmetic mean
				}
				//#if GO2TARGET
				//      if (targetseen > 0) // target is within PR - set to a very low score
				//        w.cell[x2+1][y2+1] = (1/(1000000*(double)targetseen));
				//#endif
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
