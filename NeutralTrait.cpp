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
 * File Created by Roslyn Henry March 2023. Code adapted from NEMO (https://nemo2.sourceforge.io/)
 --------------------------------------------------------------------------*/

#include "NeutralTrait.h"

// ----------------------------------------------------------------------------------------
//  Initialisation constructor
//  Called when initialising community
//  Sets up initial values, and immutable attributes (distributions and parameters)
//	that are defined at the species-level
// ----------------------------------------------------------------------------------------
NeutralTrait::NeutralTrait(SpeciesTrait* P)
{
	pSpeciesTrait = P;
	ploidy = pSpeciesTrait->getPloidy();

	DistributionType mutationDistribution = pSpeciesTrait->getMutationDistribution();
	const map<GenParamType, float>& mutationParameters = pSpeciesTrait->getMutationParameters();

	if (mutationParameters.count(MAX) != 1)
		throw logic_error("KAM or SSM mutation distribution parameter must contain max value, max cannot exceed 256  \n");

	// Set default value to user-specified max
	wildType = (int)mutationParameters.find(MAX)->second;
	if (wildType > NeutralValUpperBound)
		throw logic_error("max number of alleles cannot exceed " + to_string(NeutralValUpperBound) + ".\n");

	_inherit_func_ptr = (pSpeciesTrait->getPloidy() == 1) ? &NeutralTrait::inheritHaploid : &NeutralTrait::inheritDiploid; //this could be changed if we wanted some alternative form of inheritance

	if (mutationDistribution == SSM)
		_mutate_func_ptr = &NeutralTrait::mutate_SSM;

	if (mutationDistribution == KAM)
		_mutate_func_ptr = &NeutralTrait::mutate_KAM;

	if (mutationDistribution != SSM && mutationDistribution != KAM)
		throw logic_error("wrong mutation distribution for neutral markers, must be KAM or SSM \n");

	if (mutationParameters.count(MAX) != 1)
		throw logic_error("KAM or SSM mutation distribution parameter must contain max value, max cannot exceed 256  \n");

	DistributionType initialDistribution = pSpeciesTrait->getInitialDistribution();
	const map<GenParamType, float>& initialParameters = pSpeciesTrait->getInitialParameters();

	if (mutationDistribution == SSM && initialDistribution != UNIFORM)
		throw logic_error("If using SSM mutation model for neutral trait, must use uniform initial distribution.\n");

	switch (initialDistribution) {
	case UNIFORM:
	{
		if (initialParameters.count(MAX) != 1)
			throw logic_error("initial distribution parameter must contain one max value if set to UNIFORM (e.g. max= ), max cannot exceed " + to_string(NeutralValUpperBound) + "\n");

		float maxNeutralVal = initialParameters.find(MAX)->second;
		if (maxNeutralVal > NeutralValUpperBound) {
			throw logic_error("initial distribution parameter max cannot exceed " + to_string(NeutralValUpperBound) + ", resetting to " + to_string(NeutralValUpperBound) + "\n");
			maxNeutralVal = NeutralValUpperBound; //reserve 255 for wildtype
		}
		initialiseUniform(maxNeutralVal);
		break;
	}
	default:
	{
		throw logic_error("wrong parameter value for parameter \"initialisation of neutral trait\", must be left uniform \n");
		break; //should return false
	}
	}
}

// ----------------------------------------------------------------------------------------
// Inheritance constructor
// Copies immutable features from a parent trait
// Only called via clone()
// ----------------------------------------------------------------------------------------
NeutralTrait::NeutralTrait(const NeutralTrait& T) :
	ploidy(T.ploidy), pSpeciesTrait(T.pSpeciesTrait), _mutate_func_ptr(T._mutate_func_ptr), _inherit_func_ptr(T._inherit_func_ptr) {
}

// ----------------------------------------------------------------------------------------
// Index of a gene position in the flat genes array (throws if not a locus of this trait)
// ----------------------------------------------------------------------------------------
int NeutralTrait::locusIndexOrThrow(int position, const char* msg) const {
	const int l = pSpeciesTrait->getLocusIndex(position);
	if (l < 0 || static_cast<size_t>(l) * ploidy >= genes.size())
		throw runtime_error(msg);
	return l;
}

// ----------------------------------------------------------------------------------------
// mutate options
// ----------------------------------------------------------------------------------------

// ----------------------------------------------------------------------------------------
// Draw and apply mutations according to a KAM process
// 
// Mutations drawn only for existing positions, 
// that is no new genes are created during simulation
// KAM = randomly drawn value in 0-MAX, differs from previous value
// ----------------------------------------------------------------------------------------
void NeutralTrait::mutate_KAM()
{
	const int positionsSize = pSpeciesTrait->getPositionsSize();
	const auto& genePositions = pSpeciesTrait->getGenePositions();
	const short ploidy = pSpeciesTrait->getPloidy();
	const float mutationRate = pSpeciesTrait->getMutationRate();
	auto& rng = pRandom->getRNG();// bugfix: auto& - plain auto would copy the generator again
	unsigned char mut;

	const map<GenParamType, float>& mutationParameters = pSpeciesTrait->getMutationParameters();

	int maxNeutralVal = (int)mutationParameters.find(MAX)->second;
	if (maxNeutralVal > NeutralValUpperBound) maxNeutralVal = NeutralValUpperBound; //reserve max value for wildtype

	for (int whichChromosome = 0; whichChromosome < ploidy; whichChromosome++) {

		unsigned int NbMut = pRandom->Binomial(positionsSize, mutationRate);
		if (NbMut > 0) {
			vector<int> mutationPositions;
			sample(genePositions.begin(), genePositions.end(), std::back_inserter(mutationPositions),
				NbMut, rng); // without replacement

			for (int m : mutationPositions) {
				mut = (unsigned char)pRandom->IRandom(0, maxNeutralVal); // draw new mutation, could draw wildtype

				const int l = locusIndexOrThrow(m, "Locus selected for mutation doesn't exist.");
				unsigned char& allele = genes[static_cast<size_t>(l) * ploidy + whichChromosome];
				auto currentChar = allele; // current allele
				if (maxNeutralVal > 0) { // dodge the infinite loop
					do {
						mut = (unsigned char)pRandom->IRandom(0, maxNeutralVal);
					} while (mut == currentChar); // new allele value is different
				}
				else mut = 0; 
				allele = mut; //overwrite with new value
			}
		}
	}
}


// ----------------------------------------------------------------------------------------
// Draw and apply single-step mutations (SSM)
// 
// Mutations drawn only for existing positions, 
// that is no new genes are created during simulation
// Increment previous value by 1 or -1,
// unless already 0 (then always +1) or MAX (then always -1)
// ----------------------------------------------------------------------------------------
void NeutralTrait::mutate_SSM()
{
	const int positionsSize = pSpeciesTrait->getPositionsSize();
	const auto& genePositions = pSpeciesTrait->getGenePositions();
	const short ploidy = pSpeciesTrait->getPloidy();
	const float mutationRate = pSpeciesTrait->getMutationRate();
	auto& rng = pRandom->getRNG();// bugfix: auto& - plain auto would copy the generator again

	const map<GenParamType, float>& mutationParameters = pSpeciesTrait->getMutationParameters();

	int maxNeutralVal = (int)mutationParameters.find(MAX)->second;
	if (maxNeutralVal > NeutralValUpperBound) maxNeutralVal = NeutralValUpperBound; //reserved max value for wildtype

	for (int whichChromosome = 0; whichChromosome < ploidy; whichChromosome++) {

		unsigned int NbMut = pRandom->Binomial(positionsSize, mutationRate);
		if (NbMut > 0) {
			vector<int> mutationPositions;
			sample(genePositions.begin(), genePositions.end(), std::back_inserter(mutationPositions),
				NbMut, rng);

			for (int m : mutationPositions) {
				int mutateUp = pRandom->Bernoulli(0.5);
				const int l = locusIndexOrThrow(m, "Locus selected for mutation doesn't exist.");
				unsigned char& allele = genes[static_cast<size_t>(l) * ploidy + whichChromosome];
				auto currentAllele = allele;
				if (mutateUp == 1 && currentAllele < maxNeutralVal)
					allele += 1; // one step up
				else if (currentAllele > 0) // step down or already max
					allele -= 1; // one step down
				else // current allele is 0, step up
					allele += 1;
			}
		}
	}
}

// ----------------------------------------------------------------------------------------
//  Wrapper to inheritance function
// ----------------------------------------------------------------------------------------
void NeutralTrait::inheritGenes(const bool& fromMother, QuantitativeTrait* parent, set<unsigned int> const& recomPositions, int startingChromosome)
{
	auto parentCast = dynamic_cast<NeutralTrait*> (parent); // must convert QuantitativeTrait to NeutralTrait
	const auto& parent_seq = parentCast->getGenes();
	(this->*_inherit_func_ptr) (fromMother, parent_seq, recomPositions, startingChromosome);
}

// ----------------------------------------------------------------------------------------
// Inheritance for diploid, sexual species
// Called once for each parent. Given a list of recombinant sites, 
// populates offspring genes with appropriate parent alleles
// Assumes mother genes are inherited first
// ----------------------------------------------------------------------------------------
void NeutralTrait::inheritDiploid(const bool& fromMother, vector<unsigned char> const& parentGenes, set<unsigned int> const& recomPositions, int parentChromosome) {

	// all individuals carry all loci of the trait, in the order of the gene positions
	const vector<int>& positions = pSpeciesTrait->getGenePositionsVec();
	const int nLoci = static_cast<int>(positions.size());
	if (parentGenes.size() != static_cast<size_t>(nLoci) * 2)
		throw runtime_error("Parent neutral genotype has the wrong number of loci.");

	const int lastPosition = positions.back();
	auto recomIt = recomPositions.lower_bound(positions.front());
	// If no recombination sites, only breakpoint is last position
	// i.e., no recombination occurs
	int nextBreakpoint = recomIt == recomPositions.end() ? lastPosition : *recomIt;
	
	// Is the first parent gene position already recombinant?
	auto distance = std::distance(recomPositions.begin(), recomIt);
	if (distance % 2 != 0)
		parentChromosome = 1 - parentChromosome; //switch chromosome

	if (fromMother) {
		// loci do not exist yet: create and initialise them (paternal allele = wild type)
		if (!genes.empty()) throw runtime_error("Mother-inherited locus already exists.");
		genes.assign(static_cast<size_t>(nLoci) * 2, static_cast<unsigned char>(wildType));
	}
	else if (genes.empty()) throw runtime_error("Father-inherited locus does not exist.");
	const int offspringChromosome = fromMother ? sex_t::FEM : sex_t::MAL;

	for (int l = 0; l < nLoci; l++) {
		const int locus = positions[l];

		// Switch chromosome if locus is past recombination site
		while (locus > nextBreakpoint) {
			parentChromosome = 1 - parentChromosome;
			std::advance(recomIt, 1); // go to next recombination site
			nextBreakpoint = recomIt == recomPositions.end() ? lastPosition : *recomIt;
		}

		if (locus <= nextBreakpoint) {
			genes[static_cast<size_t>(l) * 2 + offspringChromosome] = parentGenes[static_cast<size_t>(l) * 2 + parentChromosome];
		}
	}
}

// ----------------------------------------------------------------------------------------
// Inheritance for haploid, asexual species
// Simply pass down parent genes
// Arguments are still needed to match overloaded function in base class
// ----------------------------------------------------------------------------------------
void NeutralTrait::inheritHaploid(const bool& fromMother, vector<unsigned char> const& parentGenes, set<unsigned int> const& recomPositions, int parentChromosome)
{
	genes = parentGenes;
}

// ----------------------------------------------------------------------------------------
// Initialise neutral loci
// ----------------------------------------------------------------------------------------
void NeutralTrait::initialiseUniform(int maxAlleleVal)
{
	const auto& genePositions = pSpeciesTrait->getGenePositions();
	const auto& initPositions = pSpeciesTrait->getInitPositions();
	short ploidy = pSpeciesTrait->getPloidy();

	genes.clear();
	genes.reserve(genePositions.size() * ploidy);
	for (auto position : genePositions) { // ascending order = locus index order
		for (int i = 0; i < ploidy; i++) {
			unsigned char alleleVal = char(0);
			if (initPositions.contains(position)) {
				//  allele values span 0 - max inclusive, max is wildtype
				alleleVal = (unsigned char)pRandom->IRandom(0, maxAlleleVal);
			}
			genes.push_back(alleleVal);
		}
	}
}


// ----------------------------------------------------------------------------------------
// Check if particular loci is heterozygote
// ----------------------------------------------------------------------------------------
bool NeutralTrait::isHeterozygoteAtLocus(int locus) const {
	// assumes diploidy
	const int l = locusIndexOrThrow(locus, "Neutral gene queried for heterozygosity does not exist.");
	return isHeterozygoteAtIndex(l);
}

// ----------------------------------------------------------------------------------------
// Count heterozygote loci in genome 
// ----------------------------------------------------------------------------------------
int NeutralTrait::countHeterozygoteLoci() const {
	// assumes diploidy
	int count = 0;
	for (size_t i = 0; i + 1 < genes.size(); i += 2) {
		count += (genes[i] != genes[i + 1]);
	}
	return count;
}

// ----------------------------------------------------------------------------------------
// Get allele value at loci 
// ----------------------------------------------------------------------------------------
float NeutralTrait::getAlleleValueAtLocus(short whichChromosome, int position) const {

	const int l = locusIndexOrThrow(position, "The neutral locus queried for its allele value does not exist.");
	return genes[static_cast<size_t>(l) * ploidy + whichChromosome];
}

#ifdef UNIT_TESTS // Testing only

// Set the genotype from a map position -> alleles (as used by the unit tests)
void NeutralTrait::setGenesFromMap(const map<int, vector<unsigned char>>& newGenes) {
	genes.clear();
	if (newGenes.empty()) return;
	ploidy = static_cast<short>(newGenes.begin()->second.size());
	for (int position : pSpeciesTrait->getGenePositionsVec()) {
		auto it = newGenes.find(position);
		if (it == newGenes.end()) throw logic_error("Test genotype misses a neutral locus.");
		genes.insert(genes.end(), it->second.begin(), it->second.end());
	}
}

// Create a default set of neutral alleles for testing
//
// Shorthand function to manually set genotypes for neutral
// traits, instead of having to manipulate mutations.
map<int, vector<unsigned char>> createTestNeutralGenotype(
	const int genomeSz, const bool isDiploid,
	const unsigned char valAlleleA,
	const unsigned char valAlleleB
) {
	vector<unsigned char> gene(isDiploid ? 2 : 1);
	gene[0] = valAlleleA;
	if (isDiploid) gene[1] = valAlleleB;

	map<int, vector<unsigned char>> genotype;
	for (int i = 0; i < genomeSz; i++) {
		genotype.emplace(i, gene);
	}
	return genotype;
}

#endif // UNIT_TESTS
