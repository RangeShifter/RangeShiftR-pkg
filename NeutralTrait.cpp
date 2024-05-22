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

	DistributionType mutationDistribution = pSpeciesTrait->getMutationDistribution();
	map<GenParamType, float> mutationParameters = pSpeciesTrait->getMutationParameters();

	// Set default value to user-specified max
	wildType = (int)mutationParameters.find(MAX)->second - 1;
	if (wildType > NeutralValUpperBound)
		throw logic_error("Error:: max number of alleles cannot exceed " + to_string(NeutralValUpperBound) + ".\n");

	_inherit_func_ptr = (pSpeciesTrait->getPloidy() == 1) ? &NeutralTrait::inheritHaploid : &NeutralTrait::inheritDiploid; //this could be changed if we wanted some alternative form of inheritance

	if (mutationDistribution == SSM)
		_mutate_func_ptr = &NeutralTrait::mutate_SSM;

	if (mutationDistribution == KAM)
		_mutate_func_ptr = &NeutralTrait::mutate_KAM;

	if (mutationDistribution != SSM && mutationDistribution != KAM)
		throw logic_error("Error:: wrong mutation distribution for neutral markers, must be KAM or SSM \n");

	if (mutationParameters.count(MAX) != 1)
		throw logic_error("Error:: KAM or SSM mutation distribution parameter must contain max value (e.g. max= ), max cannot exceed 256  \n");

	DistributionType initialDistribution = pSpeciesTrait->getInitialDistribution();
	map<GenParamType, float> initialParameters = pSpeciesTrait->getInitialParameters();

	if (mutationDistribution == SSM && initialDistribution != UNIFORM)
		throw logic_error("Error:: If using SSM mutation model must initialise genome with alleles (microsats) \n");

	switch (initialDistribution) {
	case UNIFORM:
	{
		if (initialParameters.count(MAX) != 1)
			throw logic_error("Error:: initial distribution parameter must contain one max value if set to UNIFORM (e.g. max= ), max cannot exceed " + to_string(NeutralValUpperBound) + "\n");

		float maxNeutralVal = initialParameters.find(MAX)->second;
		if (maxNeutralVal > NeutralValUpperBound) {
			throw logic_error("Warning:: initial distribution parameter max cannot exceed " + to_string(NeutralValUpperBound) + ", resetting to " + to_string(NeutralValUpperBound) + "\n");
			maxNeutralVal = NeutralValUpperBound; //reserve 255 for wildtype
		}
		initialiseUniform(maxNeutralVal);
		break;
	}
	case NONE: 
		break;
	default:
	{
		throw logic_error("wrong parameter value for parameter \"initialisation of neutral trait\", must be left as default (#) or uniform \n");
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
	pSpeciesTrait(T.pSpeciesTrait), _mutate_func_ptr(T._mutate_func_ptr), _inherit_func_ptr(T._inherit_func_ptr) {
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
	auto rng = pRandom->getRNG();
	unsigned char mut;

	map<GenParamType, float> mutationParameters = pSpeciesTrait->getMutationParameters();

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

				auto it = genes.find(m);
				if (it == genes.end())
					throw runtime_error("Locus selected for mutation doesn't exist.");
				auto currentChar = it->second[whichChromosome]; // current allele
				do {
					mut = (unsigned char)pRandom->IRandom(0, maxNeutralVal);
				} while (mut == currentChar); // new allele value is different

				it->second[whichChromosome] = mut; //overwrite with new value
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
	auto rng = pRandom->getRNG();

	map<GenParamType, float> mutationParameters = pSpeciesTrait->getMutationParameters();

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
				auto it = genes.find(m);
				if (it == genes.end())
					throw runtime_error("Locus selected for mutation doesn't exist.");
				auto currentAllele = it->second[whichChromosome];
				if (mutateUp == 1 && currentAllele < maxNeutralVal)
					it->second[whichChromosome] += 1; // one step up
				else if (currentAllele > 0) // step down or already max
					it->second[whichChromosome] -= 1; // one step down
				else // current allele is 0, step up
					it->second[whichChromosome] += 1;
			}
		}
	}
}

// ----------------------------------------------------------------------------------------
//  Wrapper to inheritance function
// ----------------------------------------------------------------------------------------
void NeutralTrait::inheritGenes(const bool& fromMother, TTrait* parent, set<unsigned int> const& recomPositions, int startingChromosome)
{
	auto parentCast = dynamic_cast<NeutralTrait*> (parent); // must convert TTrait to NeutralTrait
	const auto& parent_seq = parentCast->getGenes();
	(this->*_inherit_func_ptr) (fromMother, parent_seq, recomPositions, startingChromosome);
}

// ----------------------------------------------------------------------------------------
// Inheritance for diploid, sexual species
// Called once for each parent. Given a list of recombinant sites, 
// populates offspring genes with appropriate parent alleles
// Assumes mother genes are inherited first
// ----------------------------------------------------------------------------------------
void NeutralTrait::inheritDiploid(const bool& fromMother, map<int, vector<unsigned char>> const& parentGenes, set<unsigned int> const& recomPositions, int parentChromosome) {

	auto it = recomPositions.lower_bound(parentGenes.begin()->first);
	unsigned int nextBreakpoint = *it;
	auto distance = std::distance(recomPositions.begin(), it);
	if (distance - 1 % 2 != 0)
		parentChromosome = 1 - parentChromosome; //switch chromosome

	for (auto const& [locus, allelePair] : parentGenes) {

		// Switch chromosome if locus is past recombination site
		while (locus > nextBreakpoint) {
			std::advance(it, 1);
			nextBreakpoint = *it;
			parentChromosome = 1 - parentChromosome; //switch chromosome
		}

		if (locus <= nextBreakpoint) {
			unsigned char sp = allelePair[parentChromosome];
			auto it = genes.find(locus);
			if (it == genes.end()) {
				// locus does not exist yet, create and initialise it
				if (!fromMother) throw runtime_error("Father-inherited locus does not exist.");
				vector<unsigned char> newAllelePair(2, wildType);
				newAllelePair[sex_t::FEM] = sp;
				genes.insert(make_pair(locus, newAllelePair));
			}
			else { // father, locus already exists
				if (fromMother) throw runtime_error("Mother-inherited locus already exists.");
				it->second[sex_t::MAL] = sp;
			}
		}
	}
}

// ----------------------------------------------------------------------------------------
// Inheritance for haploid, asexual species
// Simply pass down parent genes
// Arguments are still needed to match overloaded function in base class
// ----------------------------------------------------------------------------------------
void NeutralTrait::inheritHaploid(const bool& fromMother, map<int, vector<unsigned char>> const& parentGenes, set<unsigned int> const& recomPositions, int parentChromosome)
{
	genes = parentGenes;
}

// ----------------------------------------------------------------------------------------
// Initialise neutral loci
// ----------------------------------------------------------------------------------------
void NeutralTrait::initialiseUniform(int maxAlleleVal)
{
	const auto& genePositions = pSpeciesTrait->getGenePositions();
	short ploidy = pSpeciesTrait->getPloidy();

	for (auto position : genePositions) {
		vector<unsigned char> allelePair;

		for (int i = 0; i < ploidy; i++) {
			//  allele values span 0 - max inclusive, max is wildtype
			auto alleleVal = (unsigned char)pRandom->IRandom(0, maxAlleleVal); 
			allelePair.emplace_back(alleleVal);
		}
		genes.insert(make_pair(position, allelePair));
	}
}


// ----------------------------------------------------------------------------------------
// Check if particular loci is heterozygote
// ----------------------------------------------------------------------------------------
bool NeutralTrait::isHeterozygoteAtLocus(int locus) const {
	// assumes diploidy
	auto it = genes.find(locus);
	if (it == genes.end()) 
		throw runtime_error("Neutral gene queried for heterozygosity does not exist.");
	else
		return(it->second[0] != it->second[1]);
}

// ----------------------------------------------------------------------------------------
// Count heterozygote loci in genome 
// ----------------------------------------------------------------------------------------
int NeutralTrait::countHeterozygoteLoci() const {
	// assumes diploidy
	int count = 0;
	for (auto const& [locus, allelePair] : genes) {
		count += (allelePair[0] != allelePair[1]);
	}
	return count;
}

// ----------------------------------------------------------------------------------------
// Get allele value at loci 
// ----------------------------------------------------------------------------------------
float NeutralTrait::getAlleleValueAtLocus(short whichChromosome, int position) const {

	auto it = genes.find(position);
	if (it == genes.end()) //no mutations there
		throw runtime_error("The neutral locus queried for its allele value does not exist.");
	return it->second[whichChromosome];
}