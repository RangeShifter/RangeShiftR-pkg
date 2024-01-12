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

//---------------------------------------------------------------------------

// Individual constructor
Individual::Individual(Cell* pCell, Patch* pPatch, short stg, short a, short repInt,
	float probmale, bool movt, short moveType)
{
	indId = indCounter;
	indCounter++; // unique identifier for each individual

	stage = stg;
	if (probmale <= 0.0) sex = 0;
	else sex = pRandom->Bernoulli(probmale);
	age = a;
	status = 0;

	if (sex == 0 && repInt > 0) { // set no. of fallow seasons for female
		fallow = pRandom->IRandom(0, repInt);
	}
	else fallow = 9999;
	isDeveloping = false;
	pPrevCell = pCurrCell = pCell;
	pNatalPatch = pPatch;
	if (movt) {
		locn loc = pCell->getLocn();
		path = new pathData;
		path->year = 0; path->total = 0; path->out = 0;
		path->pSettPatch = 0; path->settleStatus = 0;
#if RS_RCPP
		path->pathoutput = 1;
#endif
		if (moveType == 1) { // SMS
			// set up location data for SMS
			smsData = new smsdata;
			smsData->dp = smsData->gb = smsData->alphaDB = 1.0;
			smsData->betaDB = 1;
			smsData->prev.x = loc.x; 
			smsData->prev.y = loc.y; // previous location
			smsData->goal.x = loc.x; 
			smsData->goal.y = loc.y; // goal location - initialised for dispersal bias
		}
		else smsData = 0;
		if (moveType == 2) { // CRW
			// set up continuous co-ordinates etc. for CRW movement
			crw = new crwParams;
			crw->xc = ((float)pRandom->Random() * 0.999f) + (float)loc.x;
			crw->yc = (float)(pRandom->Random() * 0.999f) + (float)loc.y;
			crw->prevdrn = (float)(pRandom->Random() * 2.0 * PI);
			crw->stepL = crw->rho = 0.0;
		}
		else crw = 0;
	}
	else {
		path = 0; crw = 0; smsData = 0;
	}
	emigtraits = 0;
	kerntraits = 0;
	setttraits = 0;
	pGenome = 0;
}

Individual::~Individual(void) {
	if (path != 0) delete path;
	if (crw != 0) delete crw;
	if (smsData != 0) delete smsData;
	if (emigtraits != 0) delete emigtraits;
	if (kerntraits != 0) delete kerntraits;
	if (setttraits != 0) delete setttraits;

	if (pGenome != 0) delete pGenome;

}

//---------------------------------------------------------------------------

//---------------------------------------------------------------------------

// Set genes for individual variation from species initialisation parameters
void Individual::setGenes(Species* pSpecies, int resol) {
	demogrParams dem = pSpecies->getDemogr();
	emigRules emig = pSpecies->getEmig();
	trfrRules trfr = pSpecies->getTrfr();
	settleType sett = pSpecies->getSettle();
	genomeData gen = pSpecies->getGenomeData();
	simParams sim = paramsSim->getSim();
	int ntraits;	// first trait for all/female expression, second for male expression

	if (gen.trait1Chromosome) {
		pGenome = new Genome(pSpecies->getNChromosomes(), pSpecies->getNLoci(0),
			pSpecies->isDiploid());
	}
	else {
		pGenome = new Genome(pSpecies);
	}

	int gposn = 0;	// current position on genome
	int expr = 0;		// gene expression type - NOT CURRENTLY USED

	if (emig.indVar) { // set emigration genes
		int emigposn = gposn;
		double d0, alpha, beta;
		emigParams eparams;
		if (emig.sexDep) { // must be a sexual species
			ntraits = 2;
		}
		else {
			if (dem.repType == 0) { // asexual reproduction (haploid)
				ntraits = 1;
			}
			else { // sexual reproduction
				ntraits = 1;
			}
		}
		for (int g = 0; g < ntraits; g++) { // first trait for females/all, second for males
			eparams = pSpecies->getEmigParams(0, g);
			d0 = pRandom->Normal(0.0, eparams.d0SD) / eparams.d0Scale;
			if (emig.densDep) {
				alpha = pRandom->Normal(0.0, eparams.alphaSD) / eparams.alphaScale;
				beta = pRandom->Normal(0.0, eparams.betaSD) / eparams.betaScale;
			}
			if (gen.trait1Chromosome) {
				pGenome->setGene(gposn++, expr, d0, gen.alleleSD);
				if (emig.densDep) {
					pGenome->setGene(gposn++, expr, alpha, gen.alleleSD);
					pGenome->setGene(gposn++, expr, beta, gen.alleleSD);
				}
			}
			else {
				pGenome->setTrait(pSpecies, gposn++, d0, gen.alleleSD);
				if (emig.densDep) {
					pGenome->setTrait(pSpecies, gposn++, alpha, gen.alleleSD);
					pGenome->setTrait(pSpecies, gposn++, beta, gen.alleleSD);
				}
			}
		}
		// record phenotypic traits
		if (emig.densDep) {
			setEmigTraits(pSpecies, emigposn, 3, emig.sexDep);
		}
		else {
			setEmigTraits(pSpecies, emigposn, 1, emig.sexDep);
		}
	}

	if (trfr.indVar) { // set transfer genes
		int trfrposn = gposn;
		if (trfr.sexDep) { // must be a sexual species
			ntraits = 2;
		}
		else {
			if (dem.repType == 0) { // asexual reproduction
				ntraits = 1;
			}
			else { // sexual reproduction
				ntraits = 1;
			}
		}
		if (trfr.moveModel) {
			if (trfr.moveType == 1) { // set SMS genes
				double dp, gb, alphaDB, betaDB;
				trfrSMSParams smsparams = pSpecies->getSMSParams(0, 0); // only traits for females/all
				trfrSMSTraits smstraits = pSpecies->getSMSTraits();
				dp = pRandom->Normal(0.0, smsparams.dpSD) / smsparams.dpScale;
				gb = pRandom->Normal(0.0, smsparams.gbSD) / smsparams.gbScale;
				if (smstraits.goalType == 2) {
					alphaDB = pRandom->Normal(0.0, smsparams.alphaDBSD) / smsparams.alphaDBScale;
					betaDB = pRandom->Normal(0.0, smsparams.betaDBSD) / smsparams.betaDBScale;
				}
				if (gen.trait1Chromosome) {
					pGenome->setGene(gposn++, expr, dp, gen.alleleSD);
					pGenome->setGene(gposn++, expr, gb, gen.alleleSD);
					if (smstraits.goalType == 2) {
						pGenome->setGene(gposn++, expr, alphaDB, gen.alleleSD);
						pGenome->setGene(gposn++, expr, betaDB, gen.alleleSD);
					}
				}
				else {
					pGenome->setTrait(pSpecies, gposn++, dp, gen.alleleSD);
					pGenome->setTrait(pSpecies, gposn++, gb, gen.alleleSD);
					if (smstraits.goalType == 2) {
						pGenome->setTrait(pSpecies, gposn++, alphaDB, gen.alleleSD);
						pGenome->setTrait(pSpecies, gposn++, betaDB, gen.alleleSD);
					}
				}
				// record phenotypic traits
				if (smstraits.goalType == 2)
					setSMSTraits(pSpecies, trfrposn, 4, false);
				else
					setSMSTraits(pSpecies, trfrposn, 2, false);
			}
			if (trfr.moveType == 2) { // set CRW genes
				double stepL, rho;
				trfrCRWParams m = pSpecies->getCRWParams(0, 0); // only traits for females/all
				stepL = pRandom->Normal(0.0, m.stepLgthSD) / m.stepLScale;
				rho = pRandom->Normal(0.0, m.rhoSD) / m.rhoScale;
				if (gen.trait1Chromosome) {
					pGenome->setGene(gposn++, expr, stepL, gen.alleleSD);
					pGenome->setGene(gposn++, expr, rho, gen.alleleSD);
				}
				else {
					pGenome->setTrait(pSpecies, gposn++, stepL, gen.alleleSD);
					pGenome->setTrait(pSpecies, gposn++, rho, gen.alleleSD);
				}
				// record phenotypic traits
				setCRWTraits(pSpecies, trfrposn, 2, false);
			}
		}
		else { // set kernel genes
			double dist1, dist2, prob1;
			trfrKernParams k;
			for (int g = 0; g < ntraits; g++) { // first traits for females/all, second for males
				k = pSpecies->getKernParams(0, g);
				dist1 = pRandom->Normal(0.0, k.dist1SD) / k.dist1Scale;
				if (trfr.twinKern)
				{
					dist2 = pRandom->Normal(0.0, k.dist2SD) / k.dist2Scale;
					prob1 = pRandom->Normal(0.0, k.PKern1SD) / k.PKern1Scale;
				}
				if (gen.trait1Chromosome) {
					pGenome->setGene(gposn++, expr, dist1, gen.alleleSD);
					if (trfr.twinKern)
					{
						pGenome->setGene(gposn++, expr, dist2, gen.alleleSD);
						pGenome->setGene(gposn++, expr, prob1, gen.alleleSD);
					}
				}
				else {
					pGenome->setTrait(pSpecies, gposn++, dist1, gen.alleleSD);
					if (trfr.twinKern)
					{
						pGenome->setTrait(pSpecies, gposn++, dist2, gen.alleleSD);
						pGenome->setTrait(pSpecies, gposn++, prob1, gen.alleleSD);
					}
				}
			}
			// record phenotypic traits
			if (trfr.twinKern)
			{
				setKernTraits(pSpecies, trfrposn, 3, resol, trfr.sexDep);
			}
			else {
				setKernTraits(pSpecies, trfrposn, 1, resol, trfr.sexDep);
			}
		}
	}

	if (sett.indVar) {
		int settposn = gposn;
		double s0, alpha, beta;
		settParams sparams;
		if (sett.sexDep) { // must be a sexual species
			ntraits = 2;
		}
		else {
			if (dem.repType == 0) { // asexual reproduction
				ntraits = 1;
			}
			else { // sexual reproduction
				ntraits = 1;
			}
		}
		for (int g = 0; g < ntraits; g++) { // first trait for females/all, second for males
			if (sim.batchMode) {
				sparams = pSpecies->getSettParams(0, g);
			}
			else { // individual variability not (yet) implemented as sex-dependent in GUI
				sparams = pSpecies->getSettParams(0, 0);
			}
			s0 = pRandom->Normal(0.0, sparams.s0SD) / sparams.s0Scale;
			alpha = pRandom->Normal(0.0, sparams.alphaSSD) / sparams.alphaSScale;
			beta = pRandom->Normal(0.0, sparams.betaSSD) / sparams.betaSScale;

			if (gen.trait1Chromosome) {
				pGenome->setGene(gposn++, expr, s0, gen.alleleSD);
				pGenome->setGene(gposn++, expr, alpha, gen.alleleSD);
				pGenome->setGene(gposn++, expr, beta, gen.alleleSD);
			}
			else {
				pGenome->setTrait(pSpecies, gposn++, s0, gen.alleleSD);
				pGenome->setTrait(pSpecies, gposn++, alpha, gen.alleleSD);
				pGenome->setTrait(pSpecies, gposn++, beta, gen.alleleSD);
			}
		}
		// record phenotypic traits
		setSettTraits(pSpecies, settposn, 3, sett.sexDep);
	}

	if (!gen.trait1Chromosome) {
		if (gen.neutralMarkers || pSpecies->getNNeutralLoci() > 0) {
			pGenome->setNeutralLoci(pSpecies, gen.alleleSD);
		}
	}
}

// Inherit genome from parent(s)
void Individual::setGenes(Species* pSpecies, Individual* mother, Individual* father,
	int resol)
{
	emigRules emig = pSpecies->getEmig();
	trfrRules trfr = pSpecies->getTrfr();
	settleType sett = pSpecies->getSettle();

	Genome* pFatherGenome;
	if (father == 0) pFatherGenome = 0; else pFatherGenome = father->pGenome;

	pGenome = new Genome(pSpecies, mother->pGenome, pFatherGenome);

	if (emig.indVar) {
		// record emigration traits
		if (father == 0) { // haploid
			if (emig.densDep) {
				setEmigTraits(pSpecies, 0, 3, 0);
			}
			else {
				setEmigTraits(pSpecies, 0, 1, 0);
			}
		}
		else { // diploid
			if (emig.densDep) {
				setEmigTraits(pSpecies, 0, 3, emig.sexDep);
			}
			else {
				setEmigTraits(pSpecies, 0, 1, emig.sexDep);
			}
		}
	}

	if (trfr.indVar) {
		// record movement model traits
		if (trfr.moveModel) {
			if (trfr.moveType == 1) { // SMS
				trfrSMSTraits s = pSpecies->getSMSTraits();
				if (s.goalType == 2)
					setSMSTraits(pSpecies, trfr.movtTrait[0], 4, 0);
				else
					setSMSTraits(pSpecies, trfr.movtTrait[0], 2, 0);
			}
			if (trfr.moveType == 2) { // CRW
				setCRWTraits(pSpecies, trfr.movtTrait[0], 2, 0);
			}
		}
		else { // kernel
			if (father == 0) { // haploid
				if (trfr.twinKern)
				{
					setKernTraits(pSpecies, trfr.movtTrait[0], 3, resol, 0);
				}
				else {
					setKernTraits(pSpecies, trfr.movtTrait[0], 1, resol, 0);
				}
			}
			else { // diploid
				if (trfr.twinKern)
				{
					setKernTraits(pSpecies, trfr.movtTrait[0], 3, resol, trfr.sexDep);
				}
				else {
					setKernTraits(pSpecies, trfr.movtTrait[0], 1, resol, trfr.sexDep);
				}
			}
		}
	}

	if (sett.indVar) {
		// record settlement traits
		if (father == 0) { // haploid
			setSettTraits(pSpecies, sett.settTrait[0], 3, 0);
		}
		else { // diploid
			setSettTraits(pSpecies, sett.settTrait[0], 3, sett.sexDep);
		}
	}
}

//---------------------------------------------------------------------------

// Identify whether an individual is a potentially breeding female -
// if so, return her stage, otherwise return 0
int Individual::breedingFem(void) {
	if (sex == 0) {
		if (status == 0 || status == 4 || status == 5) return stage;
		else return 0;
	}
	else return 0;
}

int Individual::getId(void) { return indId; }

int Individual::getSex(void) { return sex; }

int Individual::getStatus(void) { return status; }

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

// Set phenotypic emigration traits
void Individual::setEmigTraits(Species* pSpecies, short emiggenelocn, short nemigtraits,
	bool sexdep) {
	emigTraits e; e.d0 = e.alpha = e.beta = 0.0;
	if (pGenome != 0) {
		if (pSpecies->has1ChromPerTrait()) {
			if (sexdep) {
				if (nemigtraits == 3) { // emigration is density-dependent
					e.d0 = (float)pGenome->express(emiggenelocn + 3 * sex, 0, 0);
					e.alpha = (float)pGenome->express(emiggenelocn + 3 * sex + 1, 0, 0);
					e.beta = (float)pGenome->express(emiggenelocn + 3 * sex + 2, 0, 0);
				}
				else {
					e.d0 = (float)pGenome->express(emiggenelocn + sex, 0, 0);
				}
			}
			else {
				e.d0 = (float)pGenome->express(emiggenelocn, 0, 0);
				if (nemigtraits == 3) { // emigration is density-dependent
					e.alpha = (float)pGenome->express(emiggenelocn + 1, 0, 0);
					e.beta = (float)pGenome->express(emiggenelocn + 2, 0, 0);
				}
			}
		}
		else {
			if (sexdep) {
				if (nemigtraits == 3) { // emigration is density-dependent
					e.d0 = (float)pGenome->express(pSpecies, emiggenelocn + 3 * sex);
					e.alpha = (float)pGenome->express(pSpecies, emiggenelocn + 3 * sex + 1);
					e.beta = (float)pGenome->express(pSpecies, emiggenelocn + 3 * sex + 2);
				}
				else {
					e.d0 = (float)pGenome->express(pSpecies, emiggenelocn + sex);
				}
			}
			else {
				e.d0 = (float)pGenome->express(pSpecies, emiggenelocn);
				if (nemigtraits == 3) { // emigration is density-dependent
					e.alpha = (float)pGenome->express(pSpecies, emiggenelocn + 1);
					e.beta = (float)pGenome->express(pSpecies, emiggenelocn + 2);
				}
			}
		}
	}

	emigParams eparams;
	if (sexdep) {
		eparams = pSpecies->getEmigParams(0, sex);
	}
	else {
		eparams = pSpecies->getEmigParams(0, 0);
	}
	emigtraits = new emigTraits;
	emigtraits->d0 = (float)(e.d0 * eparams.d0Scale + eparams.d0Mean);
	emigtraits->alpha = (float)(e.alpha * eparams.alphaScale + eparams.alphaMean);
	emigtraits->beta = (float)(e.beta * eparams.betaScale + eparams.betaMean);
	if (emigtraits->d0 < 0.0) emigtraits->d0 = 0.0;
	if (emigtraits->d0 > 1.0) emigtraits->d0 = 1.0;
	return;
}

// Get phenotypic emigration traits
emigTraits Individual::getEmigTraits(void) {
	emigTraits e; e.d0 = e.alpha = e.beta = 0.0;
	if (emigtraits != 0) {
		e.d0 = emigtraits->d0;
		e.alpha = emigtraits->alpha;
		e.beta = emigtraits->beta;
	}
	return e;
}

// Set phenotypic transfer by kernel traits
void Individual::setKernTraits(Species* pSpecies, short kerngenelocn, short nkerntraits,
	int resol, bool sexdep) {
	trfrKernTraits k; k.meanDist1 = k.meanDist2 = k.probKern1 = 0.0;
	if (pGenome != 0) {
		if (pSpecies->has1ChromPerTrait()) {
			if (sexdep) {
				if (nkerntraits == 3) { // twin kernel
					k.meanDist1 = (float)pGenome->express(kerngenelocn + 3 * sex, 0, sex);
					k.meanDist2 = (float)pGenome->express(kerngenelocn + 3 * sex + 1, 0, sex);
					k.probKern1 = (float)pGenome->express(kerngenelocn + 3 * sex + 2, 0, sex);
				}
				else {
					k.meanDist1 = (float)pGenome->express(kerngenelocn + sex, 0, sex);
				}
			}
			else {
				k.meanDist1 = (float)pGenome->express(kerngenelocn, 0, 0);
				if (nkerntraits == 3) { // twin kernel
					k.meanDist2 = (float)pGenome->express(kerngenelocn + 1, 0, 0);
					k.probKern1 = (float)pGenome->express(kerngenelocn + 2, 0, 0);
				}
			}
		}
		else {
			if (sexdep) {
				if (nkerntraits == 3) { // twin kernel
					k.meanDist1 = (float)pGenome->express(pSpecies, kerngenelocn + 3 * sex);
					k.meanDist2 = (float)pGenome->express(pSpecies, kerngenelocn + 3 * sex + 1);
					k.probKern1 = (float)pGenome->express(pSpecies, kerngenelocn + 3 * sex + 2);
				}
				else {
					k.meanDist1 = (float)pGenome->express(pSpecies, kerngenelocn + sex);
				}
			}
			else {
				k.meanDist1 = (float)pGenome->express(pSpecies, kerngenelocn);
				if (nkerntraits == 3) { // twin kernel
					k.meanDist2 = (float)pGenome->express(pSpecies, kerngenelocn + 1);
					k.probKern1 = (float)pGenome->express(pSpecies, kerngenelocn + 2);
				}
			}
		}
	}

	trfrKernParams kparams;
	if (sexdep) {
		kparams = pSpecies->getKernParams(0, sex);
	}
	else {
		kparams = pSpecies->getKernParams(0, 0);
	}
	kerntraits = new trfrKernTraits;
	kerntraits->meanDist1 = (float)(k.meanDist1 * kparams.dist1Scale + kparams.dist1Mean);
	kerntraits->meanDist2 = (float)(k.meanDist2 * kparams.dist2Scale + kparams.dist2Mean);
	kerntraits->probKern1 = (float)(k.probKern1 * kparams.PKern1Scale + kparams.PKern1Mean);
	if (!pSpecies->useFullKernel()) {
		// kernel mean(s) may not be less than landscape resolution
		if (kerntraits->meanDist1 < resol) kerntraits->meanDist1 = (float)resol;
		if (kerntraits->meanDist2 < resol) kerntraits->meanDist2 = (float)resol;
	}
	if (kerntraits->probKern1 < 0.0) kerntraits->probKern1 = 0.0;
	if (kerntraits->probKern1 > 1.0) kerntraits->probKern1 = 1.0;
	return;
}

// Get phenotypic emigration traits
trfrKernTraits Individual::getKernTraits(void) {
	trfrKernTraits k; k.meanDist1 = k.meanDist2 = k.probKern1 = 0.0;
	if (kerntraits != 0) {
		k.meanDist1 = kerntraits->meanDist1;
		k.meanDist2 = kerntraits->meanDist2;
		k.probKern1 = kerntraits->probKern1;
	}
	return k;
}

// Set phenotypic transfer by SMS traits
void Individual::setSMSTraits(Species* pSpecies, short SMSgenelocn, short nSMStraits,
	bool sexdep) {
	trfrSMSTraits s = pSpecies->getSMSTraits();
	double dp, gb, alphaDB, betaDB;
	dp = gb = alphaDB = betaDB = 0.0;
	if (pGenome != 0) {
		if (pSpecies->has1ChromPerTrait()) {
			if (sexdep) {
				dp = pGenome->express(SMSgenelocn, 0, 0);
				gb = pGenome->express(SMSgenelocn + 1, 0, 0);
				if (nSMStraits == 4) {
					alphaDB = pGenome->express(SMSgenelocn + 2, 0, 0);
					betaDB = pGenome->express(SMSgenelocn + 3, 0, 0);
				}
			}
			else {
				dp = pGenome->express(SMSgenelocn, 0, 0);
				gb = pGenome->express(SMSgenelocn + 1, 0, 0);
				if (nSMStraits == 4) {
					alphaDB = pGenome->express(SMSgenelocn + 2, 0, 0);
					betaDB = pGenome->express(SMSgenelocn + 3, 0, 0);
				}
			}
		}
		else {
			if (sexdep) {
				dp = pGenome->express(pSpecies, SMSgenelocn);
				gb = pGenome->express(pSpecies, SMSgenelocn + 1);
				if (nSMStraits == 4) {
					alphaDB = pGenome->express(pSpecies, SMSgenelocn + 2);
					betaDB = pGenome->express(pSpecies, SMSgenelocn + 3);
				}
			}
			else {
				dp = pGenome->express(pSpecies, SMSgenelocn);
				gb = pGenome->express(pSpecies, SMSgenelocn + 1);
				if (nSMStraits == 4) {
					alphaDB = pGenome->express(pSpecies, SMSgenelocn + 2);
					betaDB = pGenome->express(pSpecies, SMSgenelocn + 3);
				}
			}
		}
	}
	trfrSMSParams smsparams;
	if (sexdep) {
		smsparams = pSpecies->getSMSParams(0, 0);
	}
	else {
		smsparams = pSpecies->getSMSParams(0, 0);
	}
	smsData->dp = (float)(dp * smsparams.dpScale + smsparams.dpMean);
	smsData->gb = (float)(gb * smsparams.gbScale + smsparams.gbMean);
	if (s.goalType == 2) {
		smsData->alphaDB = (float)(alphaDB * smsparams.alphaDBScale + smsparams.alphaDBMean);
		smsData->betaDB = (int)(betaDB * smsparams.betaDBScale + smsparams.betaDBMean + 0.5);
	}
	else {
		smsData->alphaDB = s.alphaDB;
		smsData->betaDB = s.betaDB;
	}
	if (smsData->dp < 1.0) smsData->dp = 1.0;
	if (smsData->gb < 1.0) smsData->gb = 1.0;
	if (smsData->alphaDB <= 0.0) smsData->alphaDB = 0.000001f;
	if (smsData->betaDB < 1) smsData->betaDB = 1;
	return;
}

// Get phenotypic transfer by SMS traits
trfrSMSTraits Individual::getSMSTraits(void) {
	trfrSMSTraits s; s.dp = s.gb = s.alphaDB = 1.0; s.betaDB = 1;
	if (smsData != 0) {
		s.dp = smsData->dp; s.gb = smsData->gb;
		s.alphaDB = smsData->alphaDB; s.betaDB = smsData->betaDB;
	}
	return s;
}

// Set phenotypic transfer by CRW traits
void Individual::setCRWTraits(Species* pSpecies, short CRWgenelocn, short nCRWtraits,
	bool sexdep) {
	trfrCRWTraits c; c.stepLength = c.rho = 0.0;
	if (pGenome != 0) {
		if (pSpecies->has1ChromPerTrait()) {
			if (sexdep) {
				c.stepLength = (float)pGenome->express(CRWgenelocn + sex, 0, sex);
				c.rho = (float)pGenome->express(CRWgenelocn + 2 + sex, 0, sex);
			}
			else {
				c.stepLength = (float)pGenome->express(CRWgenelocn, 0, 0);
				c.rho = (float)pGenome->express(CRWgenelocn + 1, 0, 0);
			}
		}
		else {
			if (sexdep) {
				c.stepLength = (float)pGenome->express(pSpecies, CRWgenelocn + sex);
				c.rho = (float)pGenome->express(pSpecies, CRWgenelocn + 2 + sex);
			}
			else {
				c.stepLength = (float)pGenome->express(pSpecies, CRWgenelocn);
				c.rho = (float)pGenome->express(pSpecies, CRWgenelocn + 1);
			}
		}
	}

	trfrCRWParams cparams;
	if (sexdep) {
		cparams = pSpecies->getCRWParams(0, sex);
	}
	else {
		cparams = pSpecies->getCRWParams(0, 0);
	}
	crw->stepL = (float)(c.stepLength * cparams.stepLScale + cparams.stepLgthMean);
	crw->rho = (float)(c.rho * cparams.rhoScale + cparams.rhoMean);
	if (crw->stepL < 1.0) crw->stepL = 1.0;
	if (crw->rho < 0.0) crw->rho = 0.0;
	if (crw->rho > 0.999) crw->rho = 0.999f;
	return;
}

// Get phenotypic transfer by CRW traits
trfrCRWTraits Individual::getCRWTraits(void) {
	trfrCRWTraits c; c.stepLength = c.rho = 0.0;
	if (crw != 0) {
		c.stepLength = crw->stepL;
		c.rho = crw->rho;
	}
	return c;
}

// Set phenotypic settlement traits
void Individual::setSettTraits(Species* pSpecies, short settgenelocn, short nsetttraits,
	bool sexdep) {
	settleTraits s; s.s0 = s.alpha = s.beta = 0.0;
	if (pGenome != 0) {
		if (pSpecies->has1ChromPerTrait()) {
			if (sexdep) {
				s.s0 = (float)pGenome->express(settgenelocn + 3 * sex, 0, 0);
				s.alpha = (float)pGenome->express(settgenelocn + 3 * sex + 1, 0, 0);
				s.beta = (float)pGenome->express(settgenelocn + 3 * sex + 2, 0, 0);
			}
			else {
				s.s0 = (float)pGenome->express(settgenelocn, 0, 0);
				s.alpha = (float)pGenome->express(settgenelocn + 1, 0, 0);
				s.beta = (float)pGenome->express(settgenelocn + 2, 0, 0);
			}
		}
		else {
			if (sexdep) {
				s.s0 = (float)pGenome->express(pSpecies, settgenelocn + 3 * sex);
				s.alpha = (float)pGenome->express(pSpecies, settgenelocn + 3 * sex + 1);
				s.beta = (float)pGenome->express(pSpecies, settgenelocn + 3 * sex + 2);
			}
			else {
				s.s0 = (float)pGenome->express(pSpecies, settgenelocn);
				s.alpha = (float)pGenome->express(pSpecies, settgenelocn + 1);
				s.beta = (float)pGenome->express(pSpecies, settgenelocn + 2);
			}

		}
	}

	settParams sparams;
	if (sexdep) {
		sparams = pSpecies->getSettParams(0, sex);
	}
	else {
		sparams = pSpecies->getSettParams(0, 0);
	}
	setttraits = new settleTraits;
	setttraits->s0 = (float)(s.s0 * sparams.s0Scale + sparams.s0Mean);
	setttraits->alpha = (float)(s.alpha * sparams.alphaSScale + sparams.alphaSMean);
	setttraits->beta = (float)(s.beta * sparams.betaSScale + sparams.betaSMean);
	if (setttraits->s0 < 0.0) setttraits->s0 = 0.0;
	if (setttraits->s0 > 1.0) setttraits->s0 = 1.0;
	return;
}

// Get phenotypic settlement traits
settleTraits Individual::getSettTraits(void) {
	settleTraits s; s.s0 = s.alpha = s.beta = 0.0;
	if (setttraits != 0) {
		s.s0 = setttraits->s0;
		s.alpha = setttraits->alpha;
		s.beta = setttraits->beta;
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
		pCurrCell = newCell; 
		status = 5;
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
		if (pGenome != 0) {
			kern.meanDist1 = kerntraits->meanDist1;
			if (trfr.twinKern)
			{
				kern.meanDist2 = kerntraits->meanDist2;
				kern.probKern1 = kerntraits->probKern1;
			}
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
			// this should never be true??
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

	landData land = pLandscape->getLandData();
	simParams sim = paramsSim->getSim();

	trfrRules trfr = pSpecies->getTrfr();
	trfrCRWTraits movt = pSpecies->getCRWTraits();
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
				if (patch == 0) {
					pPatch = 0;
				}
				else {
					pPatch = (Patch*)patch;
				}
				if (sim.saveVisits && pPatch != pNatalPatch) {
					pCurrCell->incrVisits();
				}
			}
			break;

		case 2: // CRW
			if (trfr.indVar) {
				if (crw != 0) {
					movt.stepLength = crw->stepL;
					movt.rho = crw->rho;
				}
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
						angle = wrpcauchy(crw->prevdrn, 0.0);
					}
					else
						angle = wrpcauchy(crw->prevdrn, rho);
					// new continuous cell coordinates
					xcnew = crw->xc + sin(angle) * steplen / (float)land.resol;
					ycnew = crw->yc + cos(angle) * steplen / (float)land.resol;
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
			crw->prevdrn = (float)angle;
			crw->xc = (float)xcnew; crw->yc = (float)ycnew;
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

	if (pCurrCell == 0)
	{
		// x,y is a NODATA square - this should not occur here
		// return a negative distance to indicate an error
		move.dist = -69.0; move.cost = 0.0;
		return move;
	}

	landData land = pLand->getLandData();
	trfrSMSTraits movt = pSpecies->getSMSTraits();
	current = pCurrCell->getLocn();

	//get weights for directional persistence....
	if ((path->out > 0 && path->out <= (movt.pr + 1))
		|| natalPatch
		|| (movt.straigtenPath && path->settleStatus > 0)) {
		// inflate directional persistence to promote leaving the patch
		if (indvar) nbr = getSimDir(current.x, current.y, 10.0f * smsData->dp);
		else nbr = getSimDir(current.x, current.y, 10.0f * movt.dp);
	}
	else {
		if (indvar) nbr = getSimDir(current.x, current.y, smsData->dp);
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
			double exp_arg = -((double)nsteps - (double)smsData->betaDB) * (-smsData->alphaDB);
			if (exp_arg > 100.0) exp_arg = 100.0; // to prevent exp() overflow error
			gb = 1.0 + (smsData->gb - 1.0) / (1.0 + exp(exp_arg));
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
	else {
		// they have already been calculated - no action required
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

	if (memory.empty())
	{ // no previous movement, set matrix to unity
		for (xx = 0; xx < 3; xx++) {
			for (yy = 0; yy < 3; yy++) {
				d.cell[xx][yy] = 1;
			}
		}
	}
	else { // set up the matrix dependent on relationship of previous location to current
		d.cell[1][1] = 0;
		prev = memory.front();
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
		}
		theta = atan2(((double)x - (double)prev.x), ((double)y - (double)prev.y));
		d = calcWeightings(dp, (float)theta);

	}
	return d;
}

// Weight neighbouring cells on basis of goal bias
array3x3d Individual::getGoalBias(const int x, const int y,
	const int goaltype, const float gb)
{

	array3x3d d;
	double theta;
	int xx, yy;

	if (goaltype == 0) { // no goal set
		for (xx = 0; xx < 3; xx++) {
			for (yy = 0; yy < 3; yy++) {
				d.cell[xx][yy] = 1.0;
			}
		}
	}
	else {
		d.cell[1][1] = 0;
		if ((x - smsData->goal.x) == 0 && (y - smsData->goal.y) == 0) {
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
			theta = atan2(((double)x - (double)smsData->goal.x), ((double)y - (double)smsData->goal.y));
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
			if (x2 != 0 || y2 != 0) { // not central cell (i.e. current cell)
				for (int x3 = xmin; x3 <= xmax; x3++) {
					for (int y3 = ymin; y3 <= ymax; y3++) {
						// if cell is out of bounds, treat landscape as a torus
						// for purpose of obtaining a cost,
						if ((x + x3) < 0) x4 = x + x3 + land.maxX + 1;
						else { if ((x + x3) > land.maxX) x4 = x + x3 - land.maxX - 1; else x4 = x + x3; }
						if ((y + y3) < 0) y4 = y + y3 + land.maxY + 1;
						else { if ((y + y3) > land.maxY) y4 = y + y3 - land.maxY - 1; else y4 = y + y3; }
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
		}//end of y2 loop
	}//end of x2 loop

	return w;

}

//---------------------------------------------------------------------------
// Write records to individuals file
void Individual::outGenetics(const int rep, const int year, const int spnum,
	const int landNr, const bool xtab)
{
	if (landNr == -1) {
		if (pGenome != 0) {
			pGenome->outGenetics(rep, year, spnum, indId, xtab);
		}
	}
	else { // open/close file
		pGenome->outGenHeaders(rep, landNr, xtab);
	}

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
Cell* Individual::getCurrCell() const {
	return pCurrCell;
}

void testIndividual() {

	Patch* pPatch = new Patch(0, 0);
	int cell_x = 2;
	int cell_y = 5;
	int cell_hab = 2;
	Cell* pCell = new Cell(cell_x, cell_y, (intptr)pPatch, cell_hab);

	// Create an individual
	short stg = 0;
	short age = 0;
	short repInt = 0;
	float probmale = 0;
	bool uses_movt_process = true;
	short moveType = 1;
	Individual ind(pCell, pPatch, stg, age, repInt, probmale, uses_movt_process, moveType);

	// Gets its sex drawn from pmale
	
	// Can age or develop

	// 

	// Reproduces
	// depending on whether it is sexual or not
	// depending on the stage
	// depending on the trait inheritance


	// Disperses
	// Emigrates
	// Transfers
	// Kernel transfer
	{
		// Simple cell-based landscape layout
		// oo
		// oo
		landParams ls_params;
		ls_params.dimX = ls_params.dimY = 2;
		vector <Cell*> cells{ new Cell(0, 0, 0, 0), new Cell(1, 1, 0, 0) };
		// Set up species for habitat codes
		Species sp;
		sp.createHabK(1);
		sp.setHabK(0, 100.0); // one habitat with K = 100

		// Landscape ls = createLandscapeFromCells(cells, ls_params, sp);
			Landscape ls;
			ls.setLandParams(ls_params, true);
			// Add cells
			ls.setCellArray();
			for (auto c : cells) {
				ls.addCellToLand(c);
			}
			ls.allocatePatches(&sp);

		Patch* p = (Patch*)cells[0]->getPatch();
		Individual ind(cells[0], p, 1, 0, 0, 0.0, false, 0);
		Cell* init_cell = ind.getCurrCell();

		// ind.status
		// land.resol
		// species.trfrKernTraits.meanDist1
		// species.useFullKernel
		// pathData *path; ?
		// land.minX etc. dimX dimY
		// patch.localK
		// bool species.stageStruct
		// species trfrMortParams m; m.fixedMort = fixedMort; m.mortAlpha = mortAlpha; m.mortBeta = mortBeta;

		int isDispersing = ind.moveKernel(&ls, &sp, false);

		// After movement, individual should be...
		// in a different cell
		Cell* curr_cell = ind.getCurrCell();
		// not in a no-data cell
		assert(curr_cell != 0);
		assert(curr_cell != init_cell);
		// (non-absorbing) still within landscape boundaries
		


		// Arrival cell 

		// An individual with a small dispersal distance is unlikely to reach a distant cell

		// An individual with a large dispersal distance should be able to reach a distant cell

		// If no cell is available beyond initial cell, individual should die

		// 

	}
	// Individual::moveKernel(Landscape * pLandscape, Species * pSpecies, const short repType, const bool absorbing)
	// Settles

	// Survives

	// Develops

}
#endif // RSDEBUG

