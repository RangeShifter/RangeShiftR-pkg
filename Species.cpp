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

#include "Species.h"
//---------------------------------------------------------------------------

Species::Species(void)
{
	// initialise demographic parameters
	repType = 0; nStages = 2;
	stageStruct = false;
	propMales = 0.5; 
	harem = 1.0; 
	bc = 1.0; 
	lambda = 1.5; 
	probRep = 1.0;
	repSeasons = 1;
	repInterval = 0; 
	maxAge = 1000; 
	survival = 1;
	fecDens = false; 
	fecStageDens = false;
	devDens = false; 
	devStageDens = false;
	survDens = false; 
	survStageDens = false;
	disperseOnLoss = false;
	for (int i = 0; i < gMaxNbStages; i++) {
		for (int j = 0; j < gMaxNbSexes; j++) {
			fec[i][j] = 0.0; dev[i][j] = 0.0; surv[i][j] = 0.0;
			minAge[i][j] = 0;
		}
	}
	devCoeff = survCoeff = 1.0;
	ddwtFec = ddwtDev = ddwtSurv = 0; 
	ddwtFecDim = ddwtDevDim = ddwtSurvDim = 0;
	habK = 0; 
	habDimK = 0;
	minRK = 1.0; 
	maxRK = 2.0;

	// initialise emigration parameters
	densDepEmig = false; stgDepEmig = false; sexDepEmig = false; indVarEmig = false;
	emigStage = 0;
	for (int i = 0; i < gMaxNbStages; i++) {
		for (int j = 0; j < gMaxNbSexes; j++) {
			d0[i][j] = 0.0; alphaEmig[i][j] = 0.0; betaEmig[i][j] = 1.0;
		}
	}
	// initialise transfer parameters
	moveModel = false; stgDepTrfr = false; sexDepTrfr = false; distMort = false;
	indVarTrfr = false;
	twinKern = false;
	habMort = false;
	costMap = false;
	moveType = 1;
	for (int i = 0; i < gMaxNbStages; i++) {
		for (int j = 0; j < gMaxNbSexes; j++) {
			meanDist1[i][j] = 100.0f; meanDist2[i][j] = 1000.0f; probKern1[i][j] = 0.99f;
		}
	}
	pr = 1; prMethod = 1; memSize = 1; goalType = 0;
	dp = 1.0; gb = 1.0; alphaDB = 1.0; betaDB = 100000;
	stepMort = 0.0; stepLength = 10.0; rho = 0.9f;
	habStepMort = 0; habCost = 0;
	//costMapFile = "NULL";
	fixedMort = 0.0; mortAlpha = 0.0; mortBeta = 1.0;
	habDimTrfr = 0;
	straightenPath = false;
	fullKernel = false;
	// initialise settlement parameters
	stgDepSett = false; sexDepSett = false; indVarSett = false;
	for (int i = 0; i < gMaxNbStages; i++) {
		for (int j = 0; j < gMaxNbSexes; j++) {
			densDepSett[i][j] = false; wait[i][j] = false; go2nbrLocn[i][j] = false; findMate[i][j] = false;
			maxStepsYr[i][j] = 99999999; 	minSteps[i][j] = 0; maxSteps[i][j] = 99999999;
			s0[i][j] = 1.0; alphaS[i][j] = 0.0; betaS[i][j] = 1.0;
		}
	}
	// initialise attribute defaults
	spNum = 0;
	resetGeneticParameters();
}

Species::~Species() {
	// demographic parameters
	if (habK != NULL) deleteHabK();
	if (ddwtFec != 0) deleteDDwtFec();
	if (ddwtDev != 0) deleteDDwtDev();
	if (ddwtSurv != 0) deleteDDwtSurv();
	// transfer parameters
	if (habCost != 0 || habStepMort != 0) deleteHabCostMort();
}

short Species::getSpNum(void) { return spNum; }

//---------------------------------------------------------------------------

// Demographic functions

void Species::setDemogr(const demogrParams d) {
	if (d.repType >= 0 && d.repType <= 2) repType = d.repType;
	if (d.repType == 1 || d.repType == 2) diploid = true;
	else diploid = false;
	if (d.repSeasons >= 1) repSeasons = d.repSeasons;
	stageStruct = d.stageStruct;
	if (d.propMales > 0.0 && d.propMales < 1.0) propMales = d.propMales;
	if (d.harem > 0.0) harem = d.harem;
	if (d.bc > 0.0) bc = d.bc;
	if (d.lambda > 0.0) lambda = d.lambda;
}

demogrParams Species::getDemogrParams(void) {
	demogrParams d;
	d.repType = repType;
	d.repSeasons = repSeasons;
	d.stageStruct = stageStruct;
	d.propMales = propMales;
	d.harem = harem;
	d.bc = bc;
	d.lambda = lambda;
	return d;
}

short Species::getRepType(void) { return repType; }

bool Species::stageStructured(void) { return stageStruct; }

void Species::createHabK(short nhab) {
	if (nhab >= 0) {
		habDimK = nhab;
		if (habK != 0) deleteHabK();
		habK = new float[nhab];
		for (int i = 0; i < nhab; i++) habK[i] = 0.0;
	}
}

void Species::setHabK(short hx, float k) {
	if (hx >= 0 && hx < habDimK) {
		if (k >= 0.0) habK[hx] = k;
	}
}

float Species::getHabK(short hx) {
	float k = 0.0;
	if (hx >= 0 && hx < habDimK) k = habK[hx];
	return k;
}

float Species::getMaxK(void) {
	float k = 0.0;
	for (int i = 0; i < habDimK; i++) {
		if (habK[i] > k) k = habK[i];
	}
	return k;
}

void Species::deleteHabK(void) {
	if (habK != 0) {
		delete[] habK; habK = 0;
	}
}

void Species::setStage(const stageParams s) {
	if (s.nStages > 1) nStages = s.nStages;
	if (s.repInterval >= 0) repInterval = s.repInterval;
	if (s.maxAge >= 1) maxAge = s.maxAge;
	if (s.survival >= 0 && s.survival <= 2) survival = s.survival;
	if (s.probRep > 0.0 && s.probRep <= 1.0) probRep = s.probRep;
	fecDens = s.fecDens;		
	fecStageDens = s.fecStageDens;
	devDens = s.devDens;		
	devStageDens = s.devStageDens;
	survDens = s.survDens;	
	survStageDens = s.survStageDens;
	disperseOnLoss = s.disperseOnLoss;
}

stageParams Species::getStageParams(void) {
	stageParams s;
	s.nStages = nStages; s.repInterval = repInterval; s.maxAge = maxAge;
	s.survival = survival; s.probRep = probRep;
	s.fecDens = fecDens; s.fecStageDens = fecStageDens;
	s.devDens = devDens; s.devStageDens = devStageDens;
	s.survDens = survDens; s.survStageDens = survStageDens;
	s.disperseOnLoss = disperseOnLoss;
	return s;
}

void Species::setFec(short stg, short sex, float f) {
	// NB fecundity for stage 0 must always be zero
	if (stg > 0 && stg < gMaxNbStages && sex >= 0 && sex < gMaxNbSexes && f >= 0)
		fec[stg][sex] = f;
}

float Species::getFec(short stg, short sex) {
	if (stg >= 0 && stg < gMaxNbStages && sex >= 0 && sex < gMaxNbSexes)
		return fec[stg][sex];
	else return 0.0;
}

float Species::getMaxFec(void) {
	float maxfec = 0.0;
	if (stageStruct) {
		for (int stg = 1; stg < gMaxNbStages; stg++) {
			if (fec[stg][0] > maxfec) maxfec = fec[stg][0];
		}
	}
	else maxfec = lambda;
	return maxfec;
}

void Species::setDev(short stg, short sex, float d) {
	if (stg >= 0 && stg < gMaxNbStages && sex >= 0 && sex < gMaxNbSexes && d >= 0)
		dev[stg][sex] = d;
}

float Species::getDev(short stg, short sex) {
	if (stg >= 0 && stg < gMaxNbStages && sex >= 0 && sex < gMaxNbSexes)
		return dev[stg][sex];
	else return 0.0;
}

void Species::setSurv(short stg, short sex, float s) {
	if (stg >= 0 && stg < gMaxNbStages && sex >= 0 && sex < gMaxNbSexes && s >= 0)
		surv[stg][sex] = s;
}

float Species::getSurv(short stg, short sex) {
	if (stg >= 0 && stg < gMaxNbStages && sex >= 0 && sex < gMaxNbSexes)
		return surv[stg][sex];
	else return 0.0;
}

void Species::setMinAge(short stg, short sex, int age) {
	// NB min age for stages 0 & 1 must always be zero
	if (stg > 1 && stg < gMaxNbStages && sex >= 0 && sex < gMaxNbSexes && age >= 0)
		minAge[stg][sex] = age;
}

short Species::getMinAge(short stg, short sex) {
	if (stg >= 0 && stg < gMaxNbStages && sex >= 0 && sex < gMaxNbSexes)
		return minAge[stg][sex];
	else return 0;
}

void Species::setDensDep(float d, float s) {
	if (d > 0.0) devCoeff = d;
	if (s > 0.0) survCoeff = s;
}

densDepParams Species::getDensDep(void) {
	densDepParams d;
	d.devCoeff = devCoeff; d.survCoeff = survCoeff;
	return d;
}

void Species::createDDwtFec(short mSize) {
	if (mSize >= 0 && mSize < (gMaxNbStages * gMaxNbSexes)) {
		if (ddwtFec != 0) deleteDDwtFec();
		ddwtFecDim = mSize;
		ddwtFec = new float* [mSize];
		for (int i = 0; i < mSize; i++) {
			ddwtFec[i] = new float[mSize];
			for (int j = 0; j < mSize; j++) ddwtFec[i][j] = 1.0;
		}
	}
}

void Species::setDDwtFec(short row, short col, float f) {
	if (row >= 0 && row < ddwtFecDim && col >= 0 && col < ddwtFecDim)
		ddwtFec[row][col] = f;
}

float Species::getDDwtFec(short row, short col) {
	if (row >= 0 && row < ddwtFecDim && col >= 0 && col < ddwtFecDim)
		return ddwtFec[row][col];
	else return 0.0;
}

void Species::deleteDDwtFec(void) {
	if (ddwtFec != 0) {
		for (int i = 0; i < ddwtFecDim; i++) if (ddwtFec[i] != 0) {
			delete[] ddwtFec[i];
		}
		delete[] ddwtFec; ddwtFec = 0;
	}
}

void Species::createDDwtDev(short mSize) {
	if (mSize >= 0 && mSize < (gMaxNbStages * gMaxNbSexes)) {
		if (ddwtDev != 0) deleteDDwtDev();
		ddwtDevDim = mSize;
		ddwtDev = new float* [mSize];
		for (int i = 0; i < mSize; i++) {
			ddwtDev[i] = new float[mSize];
			for (int j = 0; j < mSize; j++) ddwtDev[i][j] = 1.0;
		}
	}
}

void Species::setDDwtDev(short row, short col, float f) {
	if (row >= 0 && row < ddwtDevDim && col >= 0 && col < ddwtDevDim)
		ddwtDev[row][col] = f;
}

float Species::getDDwtDev(short row, short col) {
	if (row >= 0 && row < ddwtDevDim && col >= 0 && col < ddwtDevDim)
		return ddwtDev[row][col];
	else return 0.0;
}

void Species::deleteDDwtDev(void) {
	if (ddwtDev != 0) {
		for (int i = 0; i < ddwtDevDim; i++) if (ddwtDev[i] != 0) {
			delete[] ddwtDev[i];
		}
		delete[] ddwtDev; ddwtDev = 0;
	}
}

void Species::createDDwtSurv(short mSize) {
	if (mSize >= 0 && mSize < (gMaxNbStages * gMaxNbSexes)) {
		if (ddwtSurv != 0) deleteDDwtSurv();
		ddwtSurvDim = mSize;
		ddwtSurv = new float* [mSize];
		for (int i = 0; i < mSize; i++) {
			ddwtSurv[i] = new float[mSize];
			for (int j = 0; j < mSize; j++) ddwtSurv[i][j] = 1.0;
		}
	}
}

void Species::setDDwtSurv(short row, short col, float f) {
	if (row >= 0 && row < ddwtSurvDim && col >= 0 && col < ddwtSurvDim)
		ddwtSurv[row][col] = f;
}

float Species::getDDwtSurv(short row, short col) {
	if (row >= 0 && row < ddwtSurvDim && col >= 0 && col < ddwtSurvDim)
		return ddwtSurv[row][col];
	else return 0.0;
}

void Species::deleteDDwtSurv(void) {
	if (ddwtSurv != 0) {
		for (int i = 0; i < ddwtSurvDim; i++) if (ddwtSurv[i] != 0) {
			delete[] ddwtSurv[i];
		}
		delete[] ddwtSurv; ddwtSurv = 0;
	}
}
//void Species::setMinMax(float min,float max) {
void Species::setMinMax(float min, float max) {
	if (min >= 0.0 && max > min) {
		minRK = min; maxRK = max;
	}
}
float Species::getMinMax(short opt) {
	if (opt == 0) return minRK;
	else return maxRK;
}

//---------------------------------------------------------------------------

bool Species::areMutationsOn(void) {
	return mutationsOn;
}

void Species::resetGeneticParameters() {
	mutationsOn = true;
	nbGeneticFitnessTraits = 0;
	genomeSize = -9999;
	recombinationRate = -9999;
	nPatchesToSample = 0;
	nIndsToSample = "";
	chromosomeEnds.clear();
	samplePatchList.clear();
}

bool Species::isDiploid() const {
	return diploid;
}

int Species::incrNbGenLoadTraits()
{
	nbGeneticFitnessTraits++;
	return nbGeneticFitnessTraits;
}

int Species::getNbGenLoadTraits() const
{
	return nbGeneticFitnessTraits;
}

void Species::addTrait(TraitType traitType, const SpeciesTrait& trait) {

	TraitType traitT = traitType;
	// hack to deal with multiple genetic load traits, could be handled better
	if (traitType == GENETIC_LOAD) {
		int n = incrNbGenLoadTraits();

		switch (n) {
		case 1:
		{
			traitT = GENETIC_LOAD1;
			break;
		}
		case 2:
		{
			traitT = GENETIC_LOAD2;
			break;
		}
		case 3:
		{
			traitT = GENETIC_LOAD3;
			break;
		}
		case 4:
		{
			traitT = GENETIC_LOAD4;
			break;
		}
		case 5:
		{
			traitT = GENETIC_LOAD5;
			break;
		}
		default:
		{
			cout << endl << ("Error:: Too many genetic load traits in Traits file, max = 5 \n");
			break;
		}
		}
	}
	spTraitTable.emplace(traitT, make_unique<SpeciesTrait>(trait));
}

SpeciesTrait* Species::getSpTrait(TraitType trait) const {
	return spTraitTable.find(trait)->second.get();
}

void Species::clearTraitTable() {
	spTraitTable.clear();
}

set<TraitType> Species::getTraitTypes() {
	auto kv = std::ranges::views::keys(spTraitTable);
	set<TraitType> keys{ kv.begin(), kv.end() };
	return keys;
}

int Species::getNTraits() const {
	return static_cast<int>(spTraitTable.size());
}

int Species::getNPositionsForTrait(const TraitType trait) const {
	return this->getSpTrait(trait)->getPositionsSize();
}

int Species::getGenomeSize() const {
	return genomeSize;
}

float Species::getRecombinationRate() const {
	return recombinationRate;
}

set<int> Species::getChromosomeEnds() const {
	return chromosomeEnds;
}

//---------------------------------------------------------------------------

// Emigration functions
void Species::setEmigRules(const emigRules e) {
	densDepEmig = e.densDep; 
	stgDepEmig = e.stgDep; 
	sexDepEmig = e.sexDep;
	indVarEmig = e.indVar;
	if (e.emigStage >= 0) emigStage = e.emigStage;
}

emigRules Species::getEmigRules(void) {
	emigRules e;
	e.densDep = densDepEmig; 
	e.stgDep = stgDepEmig; 
	e.sexDep = sexDepEmig;
	e.indVar = indVarEmig; 
	e.emigStage = emigStage;
	return e;
}

void Species::setSpEmigTraits(const short stg, const short sex, const emigTraits e) {
	if (stg >= 0 && stg < gMaxNbStages && sex >= 0 && sex < gMaxNbSexes) {
		if (e.d0 >= 0.0 && e.d0 <= 1.0) d0[stg][sex] = e.d0;
		alphaEmig[stg][sex] = e.alpha; 
		betaEmig[stg][sex] = e.beta;
	}
}

emigTraits Species::getSpEmigTraits(short stg, short sex) {
	emigTraits e;
	if (stg >= 0 && stg < gMaxNbStages && sex >= 0 && sex < gMaxNbSexes) {
		e.d0 = d0[stg][sex];
		e.alpha = alphaEmig[stg][sex];
		e.beta = betaEmig[stg][sex];
	}
	else {
		e.d0 = e.alpha = e.beta = 0.0;
	}
	return e;
}

float Species::getSpEmigD0(short stg, short sex) {
	if (stg >= 0 && stg < gMaxNbStages && sex >= 0 && sex < gMaxNbSexes) {
		return d0[stg][sex];
	}
	else {
		return 0.0;
	}
}

//---------------------------------------------------------------------------

// Transfer functions

void Species::setTrfrRules(const transferRules t) {
	moveModel = t.usesMovtProc; 
	stgDepTrfr = t.stgDep; 
	sexDepTrfr = t.sexDep;
	distMort = t.distMort;
	indVarTrfr = t.indVar;
	twinKern = t.twinKern;
	habMort = t.habMort;
	moveType = t.moveType; 
	costMap = t.costMap;
}

transferRules Species::getTransferRules(void) {
	transferRules t;
	t.usesMovtProc = moveModel; t.stgDep = stgDepTrfr; t.sexDep = sexDepTrfr;
	t.distMort = distMort; t.indVar = indVarTrfr;
	t.twinKern = twinKern;
	t.habMort = habMort;
	t.moveType = moveType; t.costMap = costMap;
	return t;
}

void Species::setFullKernel(bool k) {
	fullKernel = k;
}

bool Species::useFullKernel(void) { return fullKernel; }

void Species::setSpKernTraits(const short stg, const short sex,
	const trfrKernelParams k, const int resol)
{
	if (stg >= 0 && stg < gMaxNbStages && sex >= 0 && sex < gMaxNbSexes) {
		if (k.meanDist1 > 0.0 && k.meanDist1 >= (float)resol) meanDist1[stg][sex] = k.meanDist1;
		if (k.meanDist2 >= (float)resol) meanDist2[stg][sex] = k.meanDist2;
		if (k.probKern1 >= 0.0 && k.probKern1 <= 1.0) probKern1[stg][sex] = k.probKern1;
	}
}

trfrKernelParams Species::getSpKernTraits(short stg, short sex) {
	trfrKernelParams k;
	if (stg >= 0 && stg < gMaxNbStages && sex >= 0 && sex < gMaxNbSexes) {
		k.meanDist1 = meanDist1[stg][sex];
		k.meanDist2 = meanDist2[stg][sex];
		k.probKern1 = probKern1[stg][sex];
	}
	else {
		k.meanDist1 = 0.0; k.meanDist2 = 0.0; k.probKern1 = 1.0;
	}
	return k;
}

void Species::setMortParams(const trfrMortParams m) {
	if (m.fixedMort >= 0.0 && m.fixedMort <= 1.0) fixedMort = m.fixedMort;
	mortAlpha = m.mortAlpha;
	mortBeta = m.mortBeta;
}

trfrMortParams Species::getMortParams(void) {
	trfrMortParams m;
	m.fixedMort = fixedMort; m.mortAlpha = mortAlpha; m.mortBeta = mortBeta;
	return m;
}

void Species::setSpMovtTraits(const trfrMovtParams m) {
	if (m.pr >= 1) pr = m.pr;
	if (m.prMethod >= 1 && m.prMethod <= 3) prMethod = m.prMethod;
	if (m.memSize >= 1 && m.memSize <= 14) memSize = m.memSize;
	if (m.goalType >= 0 && m.goalType <= 2) goalType = m.goalType;
	if (m.dp >= 1.0) dp = m.dp;
	if (m.gb >= 1.0) gb = m.gb;
	if (m.alphaDB > 0.0) alphaDB = m.alphaDB;
	if (m.betaDB > 0) betaDB = m.betaDB;
	if (m.stepMort >= 0.0 && m.stepMort <= 1.0) stepMort = m.stepMort;
	if (m.stepLength > 0.0) stepLength = m.stepLength;
	if (m.rho > 0.0 && m.rho <= 1.0) rho = m.rho;
	straightenPath = m.straightenPath;
}

trfrMovtParams Species::getSpMovtTraits(void) {
	trfrMovtParams m;
	m.pr = pr; m.prMethod = prMethod; m.memSize = memSize; m.goalType = goalType;
	m.dp = dp; m.gb = gb; m.alphaDB = alphaDB;  m.betaDB = betaDB;
	m.stepMort = stepMort; m.stepLength = stepLength; m.rho = rho;
	return m;
}

trfrCRWTraits Species::getSpCRWTraits(void) {
	trfrCRWTraits m;
	m.stepMort = stepMort; m.stepLength = stepLength; m.rho = rho;
	m.straightenPath = straightenPath;
	return m;
}

trfrSMSTraits Species::getSpSMSTraits(void) {
	trfrSMSTraits m;
	m.pr = pr; m.prMethod = prMethod; m.memSize = memSize; m.goalType = goalType;
	m.dp = dp; m.gb = gb; m.alphaDB = alphaDB;  m.betaDB = betaDB; m.stepMort = stepMort;
	m.straightenPath = straightenPath;
	return m;
}

short Species::getMovtHabDim() { return habDimTrfr; }

void Species::createHabCostMort(short nhab) {
	if (nhab >= 0) {
		habDimTrfr = nhab;
		if (habCost != 0 || habStepMort != 0) deleteHabCostMort();
		habCost = new int[nhab];
		habStepMort = new double[nhab];
		for (int i = 0; i < nhab; i++) {
			habCost[i] = 1; habStepMort[i] = 0.0;
		}
	}
}

void Species::setHabCost(short hab, int cost) {
	if (hab >= 0 && hab < habDimTrfr) {
		if (cost >= 1) habCost[hab] = cost;
	}
}

void Species::setHabMort(short hab, double mort) {
	if (hab >= 0 && hab < habDimTrfr) {
		if (mort >= 0.0 && mort < 1.0) habStepMort[hab] = mort;
	}
}

int Species::getHabCost(short hab) {
	int cost = 0;
	if (hab >= 0 && hab < habDimTrfr) cost = habCost[hab];
	return cost;
}

double Species::getHabMort(short hab) {
	double pmort = 0.0;
	if (hab >= 0 && hab < habDimTrfr) pmort = habStepMort[hab];
	return pmort;
}

void Species::deleteHabCostMort(void) {
	if (habCost != 0) {
		delete[] habCost; habCost = 0;
	}
	if (habStepMort != 0) {
		delete[] habStepMort; habStepMort = 0;
	}
}

//---------------------------------------------------------------------------

// Settlement functions

void Species::setSettle(const settleType s) {
	stgDepSett = s.stgDep; sexDepSett = s.sexDep; indVarSett = s.indVar;
}

settleType Species::getSettle(void) {
	settleType s;
	s.stgDep = stgDepSett; s.sexDep = sexDepSett; s.indVar = indVarSett;
	return s;
}

void Species::setSettRules(const short stg, const short sex, const settleRules s) {
	if (stg >= 0 && stg < gMaxNbStages && sex >= 0 && sex < gMaxNbSexes) {
		densDepSett[stg][sex] = s.densDep; wait[stg][sex] = s.wait;
		go2nbrLocn[stg][sex] = s.go2nbrLocn; findMate[stg][sex] = s.findMate;
	}
}

settleRules Species::getSettRules(short stg, short sex) {
	settleRules s;
	s.densDep = false;
	s.findMate = false;
	s.go2nbrLocn = false;
	s.wait = false;
	if (stg >= 0 && stg < gMaxNbStages && sex >= 0 && sex < gMaxNbSexes) {
		s.densDep = densDepSett[stg][sex]; s.wait = wait[stg][sex];
		s.go2nbrLocn = go2nbrLocn[stg][sex]; s.findMate = findMate[stg][sex];
	}
	return s;
}

void Species::setSteps(const short stg, const short sex, const settleSteps s) {
	if (stg >= 0 && stg < gMaxNbStages && sex >= 0 && sex < gMaxNbSexes) {
		if (s.maxStepsYr >= 1) maxStepsYr[stg][sex] = s.maxStepsYr;
		else maxStepsYr[stg][sex] = 99999999;
		if (s.minSteps >= 0) minSteps[stg][sex] = s.minSteps;
		else minSteps[stg][sex] = 0;
		if (s.maxSteps >= 1) maxSteps[stg][sex] = s.maxSteps;
		else maxSteps[stg][sex] = 99999999;
	}
}

settleSteps Species::getSteps(short stg, short sex) {
	settleSteps s;
	if (stg >= 0 && stg < gMaxNbStages && sex >= 0 && sex < gMaxNbSexes) {
	    s.maxStepsYr = maxStepsYr[stg][sex];
	    s.minSteps = minSteps[stg][sex];
	    s.maxSteps = maxSteps[stg][sex];
	}
	else {
	    s.maxStepsYr = 99999999;
	    s.minSteps = 0;
	    s.maxSteps = 99999999;
	}
	return s;
}

void Species::setSpSettTraits(const short stg, const short sex, const settleTraits dd) {
	if (stg >= 0 && stg < gMaxNbStages && sex >= 0 && sex < gMaxNbSexes) {
		if (dd.s0 > 0.0 && dd.s0 <= 1.0) s0[stg][sex] = dd.s0;
		alphaS[stg][sex] = dd.alpha; betaS[stg][sex] = dd.beta;
	}
}

settleTraits Species::getSpSettTraits(short stg, short sex) {
	settleTraits dd;
	if (stg >= 0 && stg < gMaxNbStages && sex >= 0 && sex < gMaxNbSexes) {
		dd.s0 = s0[stg][sex]; dd.alpha = alphaS[stg][sex]; dd.beta = betaS[stg][sex];
	}
	else { dd.s0 = 1.0; dd.alpha = dd.beta = 0.0; }
	return dd;
}

void Species::setGeneticParameters(const std::set<int>& chromosomeEnds, const int genomeSize, const float recombinationRate,
	const std::set<int>& samplePatchList, const string nIndsToSample, const std::set<int>& stagesToSampleFrom, int nPatchesToSampleFrom)
{
	this->genomeSize = genomeSize;
	this->chromosomeEnds = chromosomeEnds;
	this->recombinationRate = recombinationRate;
	this->samplePatchList = samplePatchList;
	this->nPatchesToSample = nPatchesToSampleFrom;
	this->nIndsToSample = nIndsToSample;
	this->stagesToSampleFrom = stagesToSampleFrom;
}

// only called for cell based landscape
void Species::setSamplePatchList(const set<int>& samplePatchList) {
	this->samplePatchList = samplePatchList;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

#ifndef NDEBUG
// For testing purposes only

demogrParams createDefaultHaploidDemogrParams() {
	demogrParams d;
	d.repType = 0;
	d.repSeasons = 1;
	d.stageStruct = false;
	d.propMales = 0.0;
	d.harem = 1.0;
	d.bc = 1.0;
	d.lambda = 2.0;
	return d;
}

demogrParams createDefaultDiploidDemogrParams() {
	demogrParams d;
	d.repType = 1;
	d.repSeasons = 1;
	d.stageStruct = false;
	d.propMales = 0.5;
	d.harem = 1.0;
	d.bc = 1.0;
	d.lambda = 2.0;
	return d;
}

#endif // NDEBUG
