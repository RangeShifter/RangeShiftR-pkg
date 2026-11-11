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

#include "DispersalTrait.h"

// ----------------------------------------------------------------------------------------
//  Initialisation constructor
//  Called when initialising community
//  Sets up initial values, and immutable attributes (distributions and parameters)
//	that are defined at the species-level
// ----------------------------------------------------------------------------------------
DispersalTrait::DispersalTrait(SpeciesTrait* P)
{
	pSpeciesTrait = P;
	ploidy = pSpeciesTrait->getPloidy();
	ExpressionType expressionType = pSpeciesTrait->getExpressionType();

	if (!pSpeciesTrait->isInherited()) // there is a trait for individual variation but this isn't inherited variation it's sampled from initial distribution 
		_inherit_func_ptr = &DispersalTrait::reInitialiseGenes;
	else {
		_inherit_func_ptr = (pSpeciesTrait->getPloidy() == 1) ? &DispersalTrait::inheritHaploid : &DispersalTrait::inheritDiploid; //this could be changed if we wanted some alternative form of inheritance

		DistributionType mutationDistribution = pSpeciesTrait->getMutationDistribution();
		const map<GenParamType, float>& mutationParameters = pSpeciesTrait->getMutationParameters();

		// Set mutation parameters
		switch (mutationDistribution) {
		case UNIFORM:
		{
			if (mutationParameters.count(MAX) != 1)
				throw logic_error("mutation uniform distribution parameter must contain max value (e.g. max= ) \n");
			if (mutationParameters.count(MIN) != 1)
				throw logic_error("mutation uniform distribution parameter must contain min value (e.g. min= ) \n");
			_mutate_func_ptr = &DispersalTrait::mutateUniform;
			break;
		}
		case NORMAL:
		{
			if (mutationParameters.count(MEAN) != 1)
				throw logic_error("mutation distribution set to normal so parameters must contain mean value (e.g. mean= ) \n");
			if (mutationParameters.count(SD) != 1)
				throw logic_error("mutation distribution set to normal so parameters must contain sdev value (e.g. sdev= ) \n");
			_mutate_func_ptr = &DispersalTrait::mutateNormal;
			break;
		}
		default:
		{
			throw logic_error("wrong parameter value for mutation model, must be uniform/normal \n"); //unless want to add gamma or negative exp 
			break;
		}
		}
	}

	// Set initialisation parameters
	DistributionType initialDistribution = pSpeciesTrait->getInitialDistribution();
	const map<GenParamType, float>& initialParameters = pSpeciesTrait->getInitialParameters();
	switch (initialDistribution) {
	case UNIFORM:
	{
		if (initialParameters.count(MAX) != 1)
			throw logic_error("initial uniform distribution parameter must contain max value (e.g. max= ) \n");
		if (initialParameters.count(MIN) != 1)
			throw logic_error("initial uniform distribution parameter must contain min value (e.g. min= ) \n");
		float maxD = initialParameters.find(MAX)->second;
		float minD = initialParameters.find(MIN)->second;
		initialiseUniform(minD, maxD);
		break;
	}
	case NORMAL:
	{
		if (initialParameters.count(MEAN) != 1)
			throw logic_error("initial normal distribution parameter must contain mean value (e.g. mean= ) \n");
		if (initialParameters.count(SD) != 1)
			throw logic_error("initial normal distribution parameter must contain sdev value (e.g. sdev= ) \n");
		float mean = initialParameters.find(MEAN)->second;
		float sd = initialParameters.find(SD)->second;
		initialiseNormal(mean, sd);
		break;
	}
	default:
	{
		throw logic_error("wrong parameter value for parameter \"initialisation of dispersal traits\", must be uniform/normal \n");
		break;
	}
	}

	// Set expression mode parameters
	switch (expressionType) {
	case AVERAGE:
	{
		_express_func_ptr = &DispersalTrait::expressAverage;
		break;
	}
	case ADDITIVE:
	{
		_express_func_ptr = &DispersalTrait::expressAdditive;
		break;
	}
	default:
	{
		throw logic_error("wrong parameter value for parameter \"expression of dispersal trait\", must be average/additive \n");
		break;
	}
	}
}

// ----------------------------------------------------------------------------------------
// Inheritance constructor
// Copies immutable features from a parent trait
// Only called via clone()
// ----------------------------------------------------------------------------------------
DispersalTrait::DispersalTrait(const DispersalTrait& T) : 
	ploidy(T.ploidy),
	pSpeciesTrait(T.pSpeciesTrait), 
	_mutate_func_ptr(T._mutate_func_ptr), 
	_inherit_func_ptr(T._inherit_func_ptr), 
	_express_func_ptr(T._express_func_ptr) {}

// ----------------------------------------------------------------------------------------
// Index of a gene position in the flat genes array (throws if not a locus of this trait)
// ----------------------------------------------------------------------------------------
int DispersalTrait::locusIndexOrThrow(int position, const char* msg) const {
	const int l = pSpeciesTrait->getLocusIndex(position);
	if (l < 0 || static_cast<size_t>(l) * ploidy >= genes.size())
		throw runtime_error(msg);
	return l;
}

// ----------------------------------------------------------------------------------------
// Sample and apply mutations from a uniform distribution
// 
// Mutations drawn only for existing positions, 
// that is no new genes are created during simulation
// ----------------------------------------------------------------------------------------
void DispersalTrait::mutateUniform()
{
	const int positionsSize = pSpeciesTrait->getPositionsSize();
	const auto& genePositions = pSpeciesTrait->getGenePositions();
	const short ploidy = pSpeciesTrait->getPloidy();
	const float mutationRate = pSpeciesTrait->getMutationRate();
	float newAlleleVal;

	auto& rng = pRandom->getRNG();// bugfix: auto& - plain auto would copy the generator again

	const map<GenParamType, float>& mutationParameters = pSpeciesTrait->getMutationParameters();
	float maxD = mutationParameters.find(MAX)->second;
	float minD = mutationParameters.find(MIN)->second;

	for (int p = 0; p < ploidy; p++) {

		unsigned int NbMut = pRandom->Binomial(positionsSize, mutationRate);

		if (NbMut > 0) {
			vector<int> mutationPositions;
			sample(genePositions.begin(), genePositions.end(), std::back_inserter(mutationPositions),
				NbMut, rng);

			for (int m : mutationPositions) {
				const int l = locusIndexOrThrow(m, "Locus sampled for mutation doesn't exist.");
				float& allele = genes[static_cast<size_t>(l) * ploidy + p];
				float currentAlleleVal = allele;//current
				newAlleleVal = pRandom->FRandom(minD, maxD) + currentAlleleVal;
				allele = newAlleleVal;
			}
		}
	}
}

// ----------------------------------------------------------------------------------------
// Sample and apply mutations from a normal distribution
// Mutations drawn only for existing positions, 
// that is no new genes are created during simulation
// ----------------------------------------------------------------------------------------
void DispersalTrait::mutateNormal()
{
	const int positionsSize = pSpeciesTrait->getPositionsSize();
	const auto& genePositions = pSpeciesTrait->getGenePositions();
	const short ploidy = pSpeciesTrait->getPloidy();
	const float mutationRate = pSpeciesTrait->getMutationRate();

	auto& rng = pRandom->getRNG();

	const map<GenParamType, float>& mutationParameters = pSpeciesTrait->getMutationParameters();
	const float mean = mutationParameters.find(MEAN)->second;
	const float sd = mutationParameters.find(SD)->second;
	float newAlleleVal;

	for (int p = 0; p < ploidy; p++) {

		unsigned int NbMut = pRandom->Binomial(positionsSize, mutationRate);

		if (NbMut > 0) {
			vector<int> mutationPositions;
			sample(genePositions.begin(), genePositions.end(), std::back_inserter(mutationPositions),
				NbMut, rng);

			for (int m : mutationPositions) {
				const int l = locusIndexOrThrow(m, "Locus sampled for mutation doesn't exist.");
				float& allele = genes[static_cast<size_t>(l) * ploidy + p];
				float currentAlleleVal = allele; //current
				newAlleleVal = pRandom->Normal(mean, sd) + currentAlleleVal;
				allele = newAlleleVal;
			}
		}
	}
}

// ----------------------------------------------------------------------------------------
//  Wrapper to inheritance function
// ----------------------------------------------------------------------------------------
void DispersalTrait::inheritGenes(const bool& fromMother, QuantitativeTrait* parentTrait, set<unsigned int> const& recomPositions, int startingChromosome)
{
	auto parentCast = dynamic_cast<DispersalTrait*>(parentTrait); // must convert QuantitativeTrait to DispersalTrait
	const auto& parent_seq = parentCast->getGenes();
	(this->*_inherit_func_ptr)(fromMother, parent_seq, recomPositions, startingChromosome);
}

// ----------------------------------------------------------------------------------------
// Inheritance for diploid, sexual species
// Called once for each parent. 
// Pass the correct parental strand, resolving crossing-overs
// after each recombinant site e.g. if parent genotype is
// 0000
// 1111
// and position 2 is selected to recombine, then offspring inherits
// 0001
// Assumes mother genes are inherited first.
// ----------------------------------------------------------------------------------------
void DispersalTrait::inheritDiploid(const bool& fromMother, vector<float> const& parentGenes, set<unsigned int> const& recomPositions, int parentChromosome) {

	// all individuals carry all loci of the trait, in the order of the gene positions
	const vector<int>& positions = pSpeciesTrait->getGenePositionsVec();
	const int nLoci = static_cast<int>(positions.size());
	if (parentGenes.size() != static_cast<size_t>(nLoci) * 2)
		throw runtime_error("Parent dispersal genotype has the wrong number of loci.");

	const int lastPosition = positions.back();
	auto recomIt = recomPositions.lower_bound(positions.front());
	// If no recombination sites, only breakpoint is last position
	// i.e., no recombination occurs
	int nextBreakpoint = recomIt == recomPositions.end() ? lastPosition : *recomIt;

	// Is the first parent gene position already recombinant?
	auto distance = std::distance(recomPositions.begin(), recomIt);
	if (distance % 2 != 0) // odd positions = switch, even = switch back
		parentChromosome = 1 - parentChromosome; //switch chromosome

	if (fromMother) {
		// loci do not exist yet: create them (paternal alleles are set by the father's call)
		if (!genes.empty()) throw runtime_error("Mother-inherited locus already exists.");
		genes.assign(static_cast<size_t>(nLoci) * 2, 0.0f);
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
void DispersalTrait::inheritHaploid(const bool& fromMother, vector<float> const& parentGenes, set<unsigned int> const& recomPositions, int parentChromosome)
{
	genes = parentGenes;
}

// ----------------------------------------------------------------------------------------
// Non-inheritance 
// For cases where isInherited option is turned off
// In this case, offspring alleles are populated using the initialise functions
// Arguments are still needed to match overloaded function in base class
// ----------------------------------------------------------------------------------------
void DispersalTrait::reInitialiseGenes(const bool& fromMother, vector<float> const& parentGenes, set<unsigned int> const& recomPositions, int parentChromosome)
{
	DistributionType initialDistribution = pSpeciesTrait->getInitialDistribution();
	const map<GenParamType, float>& initialParameters = pSpeciesTrait->getInitialParameters();

	switch (initialDistribution) {
	case UNIFORM:
	{
		if (initialParameters.count(MAX) != 1)
			throw logic_error("initial uniform distribution parameter must contain max value (e.g. max= ) \n");
		if (initialParameters.count(MIN) != 1)
			throw logic_error("initial uniform distribution parameter must contain min value (e.g. min= ) \n");
		float maxD = initialParameters.find(MAX)->second;
		float minD = initialParameters.find(MIN)->second;
		initialiseUniform(minD, maxD);
		break;
	}
	case NORMAL:
	{
		if (initialParameters.count(MEAN) != 1)
			throw logic_error("initial normal distribution parameter must contain mean value (e.g. mean= ) \n");
		if (initialParameters.count(SD) != 1)
			throw logic_error("initial normal distribution parameter must contain sdev value (e.g. sdev= ) \n");
		float mean = initialParameters.find(MEAN)->second;
		float sd = initialParameters.find(SD)->second;
		initialiseNormal(mean, sd);
		break;
	}
	default:
	{
		throw logic_error("wrong parameter value for parameter \"initialisation of dispersal trait\", must be uniform/normal \n");
		break; //should return false
	}
	}
}

// ----------------------------------------------------------------------------------------
// Dispersal initialisation options
// ----------------------------------------------------------------------------------------
void DispersalTrait::initialiseNormal(float mean, float sd) {

	const int nLoci = pSpeciesTrait->getPositionsSize();

	genes.clear();
	genes.reserve(static_cast<size_t>(nLoci) * ploidy);
	for (int l = 0; l < nLoci; l++) { // loci in ascending position order
		for (int i = 0; i < ploidy; i++) {
			float alleleVal = pRandom->Normal(mean, sd);
			genes.push_back(alleleVal);
		}
	}
}

void DispersalTrait::initialiseUniform(float min, float max) {

	const int nLoci = pSpeciesTrait->getPositionsSize();

	genes.clear();
	genes.reserve(static_cast<size_t>(nLoci) * ploidy);
	for (int l = 0; l < nLoci; l++) { // loci in ascending position order
		for (int i = 0; i < ploidy; i++) {
			float alleleVal = pRandom->FRandom(min, max);
			genes.push_back(alleleVal);
		}
	}
}

// ----------------------------------------------------------------------------------------
// Dispersal gene expression options
// ----------------------------------------------------------------------------------------
float DispersalTrait::expressAdditive() {

	float phenotype = 0.0;

	for (float allele : genes) // same summation order as before: locus by locus, chromosome A then B
		phenotype += allele;
	trimPhenotype(phenotype);
	return phenotype;
}

float DispersalTrait::expressAverage() {

	int positionsSize = pSpeciesTrait->getPositionsSize();
	short ploidy = pSpeciesTrait->getPloidy();
	float phenotype = 0.0;

	for (float allele : genes) // same summation order as before
		phenotype += allele;
	phenotype /= positionsSize * ploidy;
	trimPhenotype(phenotype);
	return phenotype;
}

void DispersalTrait::trimPhenotype(float& val) {
	const float minPositiveVal = 1e-06;
	switch (pSpeciesTrait->getTraitType())
	{
	// Values bound between 0 and 1
	case E_D0_F: case E_D0_M: case E_D0:
	case S_S0_F: case S_S0_M: case S_S0:
	case KERNEL_PROBABILITY_F: case KERNEL_PROBABILITY_M: case KERNEL_PROBABILITY:
	case CRW_STEPCORRELATION:
	{
		if (val < 0.0) val = 0;
		else if (val > 1.0) val = 1.0;
		break;
	}
	// Positive values
	case KERNEL_MEANDIST_1_F: case KERNEL_MEANDIST_1_M: case KERNEL_MEANDIST_1:
	case KERNEL_MEANDIST_2_F: case KERNEL_MEANDIST_2_M: case KERNEL_MEANDIST_2:
	case CRW_STEPLENGTH:
	{
		if (val < 0.0) val = 0;
		break;
	}
	// Strictly positive values
	case SMS_ALPHADB:
	{
		if (val <= 0.0) val = minPositiveVal;
		break;
	}
	// Minimum 1
	case SMS_DP:
	case SMS_GB:
	{
		if (val <= 1.0) val = 1.0;
		break;
	}
	// Not bound
	// Slopes of the density reaction norms may be negative: a negative settlement
	// slope gives settlement that decreases with density (eqn. 16), a negative
	// emigration slope negative density-dependent emigration (cf. isValidTraitVal())
	case E_ALPHA_F: case E_ALPHA_M: case E_ALPHA:
	case S_ALPHA_F: case S_ALPHA_M: case S_ALPHA:
	case E_BETA_F: case E_BETA_M: case E_BETA:
	case S_BETA_F: case S_BETA_M: case S_BETA:
	case SMS_BETADB:
	{
		break;
	}
	default:
		break;
	}
}

// ----------------------------------------------------------------------------------------
// Get allele value at locus
// ----------------------------------------------------------------------------------------
float DispersalTrait::getAlleleValueAtLocus(short whichChromosome, int position) const {

	const int l = locusIndexOrThrow(position, "The Dispersal locus queried for its allele value does not exist.");
	return genes[static_cast<size_t>(l) * ploidy + whichChromosome];
}

float DispersalTrait::getDomCoefAtLocus(short whichChromosome, int position) const {
	locusIndexOrThrow(position, "The genetic load locus queried for its dominance coefficient does not exist.");
	return static_cast<float>(dispDominanceFactor); // all dispersal alleles have the same dominance
}

#ifdef UNIT_TESTS

// Convert a test genotype (position -> alleles) to the flat representation
vector<float> DispersalTrait::flattenTestGenotype(const map<int, vector<shared_ptr<Allele>>>& genSeq) const {
	vector<float> flat;
	for (int position : pSpeciesTrait->getGenePositionsVec()) {
		auto it = genSeq.find(position);
		if (it == genSeq.end()) throw logic_error("Test genotype misses a dispersal locus.");
		for (const auto& pAllele : it->second)
			flat.push_back(pAllele ? pAllele->getAlleleValue() : 0.0f);
	}
	return flat;
}

// Create a default set of alleles for testing
//
// Shorthand function to manually set genotypes for dispersal and 
// genetic fitness traits, instead of having to manipulate mutations.
map<int, vector<shared_ptr<Allele>>> createTestGenotype(
	const int genomeSz, const bool isDiploid,
	const float valAlleleA,
	const float valAlleleB,
	const float domCoeffA,
	const float domCoeffB
) {
	vector<shared_ptr<Allele>> gene(isDiploid ? 2 : 1);
	if (isDiploid) {
		gene[0] = make_shared<Allele>(valAlleleA, domCoeffA);
		gene[1] = make_shared<Allele>(valAlleleB, domCoeffB);
	}
	else {
		gene[0] = make_shared<Allele>(valAlleleA, domCoeffA);
	}
	map<int, vector<shared_ptr<Allele>>> genotype;
	for (int i = 0; i < genomeSz; i++) {
		genotype.emplace(i, gene);
	}
	return genotype;
}
#endif // UNIT_TESTS
