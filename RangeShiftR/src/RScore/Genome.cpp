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

#include "Genome.h"
//---------------------------------------------------------------------------

ofstream outGenetic;

//---------------------------------------------------------------------------

Chromosome::Chromosome(int nloc)
{
	if (nloc > 0) nloci = nloc; else nloci = 1;
	pLoci = new locus[nloci];
	for (int i = 0; i < nloci; i++) {
		pLoci[i].allele[0] = pLoci[i].allele[1] = 0;
	}
}

Chromosome::~Chromosome() {
	if (pLoci != 0) {
		delete[] pLoci; pLoci = NULL;
	}
}

short Chromosome::nLoci(void) { return nloci; }

locus Chromosome::alleles(const int loc) { // return allele values at a specified locus
	locus l; l.allele[0] = l.allele[1] = 0;
	if (loc >= 0 && loc < nloci) {
		l.allele[0] = pLoci[loc].allele[0]; l.allele[1] = pLoci[loc].allele[1];
	}
	return l;
}

double Chromosome::additive(const bool diploid) {
	int sum = 0;
	for (int i = 0; i < nloci; i++) {
		sum += pLoci[i].allele[0];
		if (diploid) sum += pLoci[i].allele[1];
	}
	return (double)sum / INTBASE;
}

double Chromosome::meanvalue(const bool diploid) {
	int sum = 0;
	double mean;
	for (int i = 0; i < nloci; i++) {
		sum += pLoci[i].allele[0];
		if (diploid) sum += pLoci[i].allele[1];
	}
	mean = (double)sum / (double)nloci;
	if (diploid) mean /= 2.0;
	mean /= INTBASE;
	return  mean;
}

double Chromosome::additive(const short loc, const bool diploid) {
	int sum = 0;
	sum += pLoci[loc].allele[0];
	if (diploid) sum += pLoci[loc].allele[1];
	return (double)sum / INTBASE;
}

// Set up chromosome at simulation initialisation
void Chromosome::initialise(const double mean, const double sd,
	const bool diploid) {
	double avalue;
	double intbase = INTBASE;

	for (int i = 0; i < nloci; i++) {
		avalue = pRandom->Normal(mean, sd);
		if (avalue > 0.0)
			pLoci[i].allele[0] = (int)(avalue * intbase + 0.5);
		else
			pLoci[i].allele[0] = (int)(avalue * intbase - 0.5);
		if (diploid) {
			avalue = pRandom->Normal(mean, sd);
			if (avalue > 0.0)
				pLoci[i].allele[1] = (int)(avalue * intbase + 0.5);
			else
				pLoci[i].allele[1] = (int)(avalue * intbase - 0.5);
		}
	}

}

// Set up specified locus at simulation initialisation
void Chromosome::initialise(const short locus, const short posn, const int aval)
{
	// note that initialising value is ADDED to current value to allow for pleiotropy
	pLoci[locus].allele[posn] += aval;
}

// Inherit from specified parent
void Chromosome::inherit(const Chromosome* parentChr, const short posn, const short nloc,
	const double probmutn, const double probcross, const double mutnSD, const bool diploid)
{
	// NOTE: At present for diploid genome, presence of crossover is determined at each
	// locus (except first). However, Roslyn has shown that it is more efficient to sample
	// crossover locations from geometric distribution if number of loci is large.
	// HOW LARGE IS 'LARGE' IN THIS CASE?...

	int ix = 0; // indexes maternal and paternal strands
	if (diploid) ix = pRandom->Bernoulli(0.5); // start index at random
	for (int i = 0; i < nloc; i++) {
		if (diploid) {
			pLoci[i].allele[posn] = parentChr->pLoci[i].allele[ix];
			if (pRandom->Bernoulli(probcross)) { // crossover occurs
				if (ix == 0) ix = 1; else ix = 0;
			}
		}
		else
			pLoci[i].allele[posn] = parentChr->pLoci[i].allele[0];
		if (pRandom->Bernoulli(probmutn)) { // mutation occurs
			double intbase = INTBASE;
#if RSDEBUG
			int oldval = pLoci[i].allele[posn];
#endif
			double mutnvalue = pRandom->Normal(0, mutnSD);
			if (mutnvalue > 0.0)
				pLoci[i].allele[posn] += (int)(intbase * mutnvalue + 0.5);
			else
				pLoci[i].allele[posn] += (int)(intbase * mutnvalue - 0.5);
#if RSDEBUG
			MUTNLOG << mutnvalue << " " << oldval << " " << pLoci[i].allele[posn] << " " << endl;
#endif
		}
	}
}


//---------------------------------------------------------------------------

// NB THIS FUNCTION IS CURRENTLY NOT BEING CALLED TO CONSTRUCT AN INSTANCE OF Genome
// Genome(int) IS USED INSTEAD

Genome::Genome() {
	pChromosome = NULL;
	nChromosomes = 0;
}

// Set up new genome at initialisation for 1 chromosome per trait
Genome::Genome(int nchromosomes, int nloci, bool d) {

	diploid = d;
	if (nchromosomes > 0) nChromosomes = nchromosomes; else nChromosomes = 1;
	pChromosome = new Chromosome * [nChromosomes];
	for (int i = 0; i < nChromosomes; i++) {
		pChromosome[i] = new Chromosome(nloci);
	}

}

// Set up new genome at initialisation for trait mapping
Genome::Genome(Species* pSpecies) {
	int nloci;
	nChromosomes = pSpecies->getNChromosomes();
	diploid = pSpecies->isDiploid();
	pChromosome = new Chromosome * [nChromosomes];
	for (int i = 0; i < nChromosomes; i++) {
		nloci = pSpecies->getNLoci(i);
		pChromosome[i] = new Chromosome(nloci);
	}
}

// Inherit genome from parent(s)
Genome::Genome(Species* pSpecies, Genome* mother, Genome* father)
{
	genomeData gen = pSpecies->getGenomeData();

	nChromosomes = mother->nChromosomes;
	diploid = mother->diploid;
	pChromosome = new Chromosome * [nChromosomes];

	for (int i = 0; i < nChromosomes; i++) {
		pChromosome[i] = new Chromosome(mother->pChromosome[i]->nLoci());
		inherit(mother, 0, i, gen.probMutn, gen.probCrossover, gen.mutationSD);
		if (diploid) {
			if (father == 0) { // species is hermaphrodite - inherit again from mother
				inherit(mother, 1, i, gen.probMutn, gen.probCrossover, gen.mutationSD);
			}
			else inherit(father, 1, i, gen.probMutn, gen.probCrossover, gen.mutationSD);
		}
	}

}

Genome::~Genome() {

	if (pChromosome == NULL) return;

	for (int i = 0; i < nChromosomes; i++) {
		delete pChromosome[i];
	}
	delete[] pChromosome;

}

//---------------------------------------------------------------------------

void Genome::setDiploid(bool dip) { diploid = dip; }
bool Genome::isDiploid(void) { return diploid; }
short Genome::getNChromosomes(void) { return nChromosomes; }

//---------------------------------------------------------------------------

// Inherit from specified parent
void Genome::inherit(const Genome* parent, const short posn, const short chr,
	const double probmutn, const double probcross, const double mutnSD)
{
	pChromosome[chr]->inherit(parent->pChromosome[chr], posn, parent->pChromosome[chr]->nLoci(),
		probmutn, probcross, mutnSD, diploid);

}

void Genome::outGenHeaders(const int rep, const int landNr, const bool xtab)
{

	if (landNr == -999) { // close file
		if (outGenetic.is_open()) {
			outGenetic.close(); outGenetic.clear();
		}
		return;
	}

	string name;
	simParams sim = paramsSim->getSim();

	if (sim.batchMode) {
		name = paramsSim->getDir(2)
			+ "Batch" + Int2Str(sim.batchNum) + "_"
			+ "Sim" + Int2Str(sim.simulation)
			+ "_Land" + Int2Str(landNr) + "_Rep" + Int2Str(rep) + "_Genetics.txt";
	}
	else {
		name = paramsSim->getDir(2) + "Sim" + Int2Str(sim.simulation)
			+ "_Rep" + Int2Str(rep) + "_Genetics.txt";
	}
	outGenetic.open(name.c_str());

	outGenetic << "Rep\tYear\tSpecies\tIndID";
	if (xtab) {
		for (int i = 0; i < nChromosomes; i++) {
			int nloci = pChromosome[i]->nLoci();
			for (int j = 0; j < nloci; j++) {
				outGenetic << "\tChr" << i << "Loc" << j << "Allele0";
				if (diploid) outGenetic << "\tChr" << i << "Loc" << j << "Allele1";
			}
		}
		outGenetic << endl;
	}
	else {
		outGenetic << "\tChromosome\tLocus\tAllele0";
		if (diploid) outGenetic << "\tAllele1";
		outGenetic << endl;
	}

}

void Genome::outGenetics(const int rep, const int year, const int spnum,
	const int indID, const bool xtab)
{
	locus l;
	if (xtab) {
		outGenetic << rep << "\t" << year << "\t" << spnum << "\t" << indID;
		for (int i = 0; i < nChromosomes; i++) {
			int nloci = pChromosome[i]->nLoci();
			for (int j = 0; j < nloci; j++) {
				l = pChromosome[i]->alleles(j);
				outGenetic << "\t" << l.allele[0];
				if (diploid) outGenetic << "\t" << l.allele[1];
			}
		}
		outGenetic << endl;
	}
	else {
		for (int i = 0; i < nChromosomes; i++) {
			int nloci = pChromosome[i]->nLoci();
			for (int j = 0; j < nloci; j++) {
				outGenetic << rep << "\t" << year << "\t" << spnum << "\t"
					<< indID << "\t" << i << "\t" << j;
				l = pChromosome[i]->alleles(j);
				outGenetic << "\t" << l.allele[0];
				if (diploid) outGenetic << "\t" << l.allele[1];
				outGenetic << endl;
			}
		}
	}
}

//---------------------------------------------------------------------------

// Set up new gene at initialisation for 1 chromosome per trait
void Genome::setGene(const short chr, const short exp,
	const double traitval, const double alleleSD)
	// NB PARAMETER exp FOR EXPRESSION TYPE IS NOT CURRENTLY USED...
{
	if (chr >= 0 && chr < nChromosomes) {
		pChromosome[chr]->initialise(traitval, alleleSD, diploid);
	}
}

// Set up trait at initialisation for trait mapping
void Genome::setTrait(Species* pSpecies, const int trait,
	const double traitval, const double alleleSD)
{
	traitAllele allele;
	int nalleles = pSpecies->getNTraitAlleles(trait);
	int ntraitmaps = pSpecies->getNTraitMaps();

	int avalue;
	double intbase = INTBASE;
	if (trait < ntraitmaps) {
		for (int i = 0; i < nalleles; i++) {
			allele = pSpecies->getTraitAllele(trait, i);
			avalue = (int)(pRandom->Normal(traitval, alleleSD) * intbase);
			pChromosome[allele.chromo]->initialise(allele.locus, 0, avalue);
			if (diploid) {
				avalue = (int)(pRandom->Normal(traitval, alleleSD) * intbase);
				pChromosome[allele.chromo]->initialise(allele.locus, 1, avalue);
			}
		}
	}
	else { // insufficient traits were defined
		// alleles cannot be initialised - all individuals have mean phenotype
	}

}

// Set up trait at initialisation for trait mapping
void Genome::setNeutralLoci(Species* pSpecies, const double alleleSD)
{
	traitAllele allele;
	int nneutral = pSpecies->getNNeutralLoci();

	double avalue;
	double intbase = INTBASE;
	for (int i = 0; i < nneutral; i++) {
		allele = pSpecies->getNeutralAllele(i);
		avalue = pRandom->Normal(0.0, alleleSD);
		if (avalue > 0.0)
			pChromosome[allele.chromo]->initialise(allele.locus, 0, (int)(avalue * intbase + 0.5));
		else
			pChromosome[allele.chromo]->initialise(allele.locus, 0, (int)(avalue * intbase - 0.5));
		if (diploid) {
			avalue = pRandom->Normal(0.0, alleleSD);
			if (avalue > 0.0)
				pChromosome[allele.chromo]->initialise(allele.locus, 1, (int)(avalue * intbase + 0.5));
			else
				pChromosome[allele.chromo]->initialise(allele.locus, 1, (int)(avalue * intbase - 0.5));
		}
	}
}

// Return the expressed value of a gene when species has one chromosome per trait
double Genome::express(short chr, short expr, short indsex)
{
	double genevalue = 0.0;
	genevalue = pChromosome[chr]->meanvalue(diploid);
	return genevalue;
}

// Return the expressed value of a trait when genetic architecture is defined
double Genome::express(Species* pSpecies, short traitnum)
{
	double genevalue = 0.0;

	traitAllele allele;
	int nalleles = pSpecies->getNTraitAlleles(traitnum);
	if (nalleles > 0) {
		for (int i = 0; i < nalleles; i++) {
			allele = pSpecies->getTraitAllele(traitnum, i);
			genevalue += pChromosome[allele.chromo]->additive(allele.locus, diploid);
		}
		genevalue /= (double)nalleles;
		if (diploid) genevalue /= 2.0;
	}
	return genevalue;
}


locusOK Genome::getAlleles(short chr, short loc) {
	locusOK l;
	l.allele[0] = l.allele[1] = 0; l.ok = false;
	if (chr >= 0 && chr < nChromosomes) {
		if (pChromosome[chr] != 0) {
			if (loc >= 0 && loc < pChromosome[chr]->nLoci()) {
				locus a = pChromosome[chr]->alleles(loc);
				l.allele[0] = a.allele[0]; l.allele[1] = a.allele[1]; l.ok = true;
			}
		}
	}

	return l;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

