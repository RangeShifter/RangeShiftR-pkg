#include "SNPTrait.h"

// ----------------------------------------------------------------------------------------
// for initialising population
// ----------------------------------------------------------------------------------------
SNPTrait::SNPTrait(SpeciesTrait* P)
{
	pSpeciesTrait = P;

	DistributionType mutationDistribution = pSpeciesTrait->getMutationDistribution();
	map<GenParamType, float> mutationParameters = pSpeciesTrait->getMutationParameters();

	_inherit_func_ptr = (pSpeciesTrait->getPloidy() == 1) ? &SNPTrait::inheritHaploid : &SNPTrait::inheritDiploid; //this could be changed if we wanted some alternative form of inheritance

	if (mutationDistribution == SSM)
		_mutate_func_ptr = &SNPTrait::mutate_SSM;

	if (mutationDistribution == KAM)
		_mutate_func_ptr = &SNPTrait::mutate_KAM;

	if (mutationDistribution != SSM && mutationDistribution != KAM)
		cout << endl << ("Error:: wrong mutation distribution for neutral markers, must be KAM or SSM \n");

	if (mutationParameters.count(MAX) != 1)
		cout << endl << ("Error:: KAM or SSM mutation distribution parameter must contain max value (e.g. max= ), max cannot exceed 256  \n");

	if (wildType == -999)
		wildType = (int)mutationParameters.find(MAX)->second - 1;

	if (wildType > SNPvalUpperBound)
		cout << endl << "Error:: max number of alleles cannot exceed " << SNPvalUpperBound << ".\n";

	DistributionType initialDistribution = pSpeciesTrait->getInitialDistribution();
	map<GenParamType, float> initialParameters = pSpeciesTrait->getInitialParameters();

	if (mutationDistribution == SSM && initialDistribution != UNIFORM)
		cout << endl << ("Error:: If using SSM mutation model must initialise genome with alleles (microsats) \n");

	switch (initialDistribution) {
	case UNIFORM:
	{
		if (initialParameters.count(MAX) != 1)
			cout << endl << "Error:: initial SNP/Microsat distribution parameter must contain one max value if set to UNIFORM (e.g. max= ), max cannot exceed " << SNPvalUpperBound << "\n";

		float maxSNPval = initialParameters.find(MAX)->second;
		if (maxSNPval > SNPvalUpperBound) {
			cout << endl << "Warning:: initial SNP/Microsat distribution parameter max cannot exceed " << SNPvalUpperBound << ", resetting to " << SNPvalUpperBound << "\n";
			maxSNPval = SNPvalUpperBound; //reserve 255 for wildtype
		}
		initialiseUniform(maxSNPval);

		break;
	}
	case NONE:
	{ 
		break; 
	}
	default:
	{
		cout << endl << ("wrong parameter value for parameter \"initialisation of snp/microsat\", must be left as default (#) or uniform \n");
		break; //should return false
	}
	}
}

// ----------------------------------------------------------------------------------------
// for creating new individuals
// ----------------------------------------------------------------------------------------

SNPTrait::SNPTrait(const SNPTrait& T) :
	pSpeciesTrait(T.pSpeciesTrait), _mutate_func_ptr(T._mutate_func_ptr), _inherit_func_ptr(T._inherit_func_ptr) {

}


// ----------------------------------------------------------------------------------------
// mutate options
// ----------------------------------------------------------------------------------------

// ----------------------------------------------------------------------------------------
// mutate_KAM
// ----------------------------------------------------------------------------------------
void SNPTrait::mutate_KAM()
{
	const int positionsSize = pSpeciesTrait->getPositionsSize();
	const auto& genePositions = pSpeciesTrait->getGenePositions();
	const short ploidy = pSpeciesTrait->getPloidy();
	const float mutationRate = pSpeciesTrait->getMutationRate();
	auto rng = pRandom->getRNG();
	unsigned char mut;

	map<GenParamType, float> mutationParameters = pSpeciesTrait->getMutationParameters();

	int maxSNPval = (int)mutationParameters.find(MAX)->second;

	if (maxSNPval > SNPvalUpperBound) maxSNPval = SNPvalUpperBound; //reserve max value for wildtype

	for (int whichChromosome = 0; whichChromosome < ploidy; whichChromosome++) {

		unsigned int NbMut = pRandom->Poisson(positionsSize * mutationRate);
		if (NbMut > positionsSize) NbMut = positionsSize;

		if (NbMut > 0) {
			vector<int> mutationPositions;
			sample(genePositions.begin(), genePositions.end(), std::back_inserter(mutationPositions),
				NbMut, rng); // without replacement

			for (int m : mutationPositions) {

				mut = (unsigned char)pRandom->IRandom(0, maxSNPval); // draw new mutation, could draw wildtype

				auto it = genes.find(m);
				if (it == genes.end())
					throw runtime_error("Locus selected for mutation doesn't exist.");

				auto currentChar = it->second[whichChromosome]; //current mutation
				do {
					mut = (unsigned char)pRandom->IRandom(0, maxSNPval); //make sure new value differs from old , could be a problem here with infinite loops
				} while (mut == currentChar);

				it->second[whichChromosome] = mut; //overwrite with new value
			}
		}
	}
}


// ----------------------------------------------------------------------------------------
// mutate_SSM
// ----------------------------------------------------------------------------------------
void SNPTrait::mutate_SSM()
{
	const int positionsSize = pSpeciesTrait->getPositionsSize();
	const auto& genePositions = pSpeciesTrait->getGenePositions();
	const short ploidy = pSpeciesTrait->getPloidy();
	const float mutationRate = pSpeciesTrait->getMutationRate();
	auto rng = pRandom->getRNG();

	map<GenParamType, float> mutationParameters = pSpeciesTrait->getMutationParameters();

	int maxSNPval = (int)mutationParameters.find(MAX)->second;
	if (maxSNPval > SNPvalUpperBound) maxSNPval = SNPvalUpperBound; //reserved max value for wildtype

	for (int whichChromosome = 0; whichChromosome < ploidy; whichChromosome++) {

		unsigned int NbMut = pRandom->Poisson(positionsSize * mutationRate);
		if (NbMut > positionsSize) NbMut = positionsSize;

		if (NbMut > 0) {
			vector<int> mutationPositions;
			sample(genePositions.begin(), genePositions.end(), std::back_inserter(mutationPositions),
				NbMut, rng);

			for (int m : mutationPositions) {
				int mutateForward = pRandom->Bernoulli(0.5);
				auto it = genes.find(m);
				if (it == genes.end())
					throw runtime_error("Locus selected for mutation doesn't exist.");

				auto currentAllele = it->second[whichChromosome];
				if (mutateForward == 1 && currentAllele < maxSNPval)
					it->second[whichChromosome] += 1; // shift one char to the right
				else if (currentAllele > 0)
					it->second[whichChromosome] -= 1; // shift one char to the left
				else // current allele is 0 and mutate backwards
					it->second[whichChromosome] += 1;
			}
		}
	}
}

// ----------------------------------------------------------------------------------------
// inheritance options
// ----------------------------------------------------------------------------------------

void SNPTrait::inherit(const bool& fromMother, TTrait* parent, set<unsigned int> const& recomPositions, int startingChromosome)
{
	auto parentCast = dynamic_cast<SNPTrait*> (parent); //horrible

	const auto& parent_seq = parentCast->getGenes();
	if (parent_seq.size() > 0) //else nothing to inherit
		(this->*_inherit_func_ptr) (fromMother, parent_seq, recomPositions, startingChromosome);
}



void SNPTrait::inheritDiploid(const bool& fromMother, map<int, vector<unsigned char>> const& parentGenes, set<unsigned int> const& recomPositions, int parentChromosome) {

	if (parentGenes.size() > 0) {
		auto it = recomPositions.lower_bound(parentGenes.begin()->first);

		unsigned int nextBreakpoint = *it;

		auto distance = std::distance(recomPositions.begin(), it);
		if (distance - 1 % 2 != 0)
			parentChromosome = 1 - parentChromosome; //switch chromosome


		for (auto const& [locus, allelePair] : parentGenes) {

			while (locus > nextBreakpoint) {
				std::advance(it, 1);
				nextBreakpoint = *it;
				parentChromosome = 1 - parentChromosome; //switch chromosome
			}

			if (locus <= nextBreakpoint) {
				unsigned char sp = allelePair[parentChromosome];

				auto it = genes.find(locus);
				if (it == genes.end()) {
					if (!fromMother) throw runtime_error("Father-inherited locus does not exist.");
					// not found
					vector<unsigned char> newAllelePair(2, wildType);
					newAllelePair[sex_t::FEM] = sp;
					genes.insert(make_pair(locus, newAllelePair));
				}
				else {
					if (fromMother) throw runtime_error("Mother-inherited locus already exists.");
					it->second[sex_t::MAL] = sp;
				}
			}
		}
	}
}

void SNPTrait::inheritHaploid(const bool& fromMother, map<int, vector<unsigned char>> const& parentGenes, set<unsigned int> const& recomPositions, int parentChromosome)
{
	genes = parentGenes;
}

// ----------------------------------------------------------------------------------------
// Initialise neutral loci
// ----------------------------------------------------------------------------------------

void SNPTrait::initialiseUniform(int maxAlleleVal)
{
	const auto& genePositions = pSpeciesTrait->getGenePositions();
	short ploidy = pSpeciesTrait->getPloidy();

	for (auto position : genePositions) {
		vector<unsigned char> allelePair;

		for (int i = 0; i < ploidy; i++) {
			auto alleleVal = (unsigned char)pRandom->IRandom(0, maxAlleleVal); //  allele values span 0 - max (255 ceiling) inclusive, max == wildtype
			allelePair.emplace_back(alleleVal);
		}
		genes.insert(make_pair(position, allelePair));
	}
}


// ----------------------------------------------------------------------------------------
// check if particular loci is heterozygote
// ----------------------------------------------------------------------------------------


bool SNPTrait::isHeterozygoteAtLocus(int locus) const {

	auto it = genes.find(locus);

	if (it == genes.end()) //not found
		return false;
	else
		return(it->second[0] != it->second[1]);
}

// ----------------------------------------------------------------------------------------
// count heterozygote loci in genome 
// ----------------------------------------------------------------------------------------

int SNPTrait::countHeterozygoteLoci() const {

	int count = 0;
	for (auto const& [locus, allelePair] : genes) {
			count += (allelePair[0] != allelePair[1]);
	}
	return count;
}

// ----------------------------------------------------------------------------------------
// get allele value at loci 
// ----------------------------------------------------------------------------------------


float SNPTrait::getAlleleValueAtLocus(short whichChromosome, int position) const {

	auto it = genes.find(position);

	if (it == genes.end()) //no mutations there
		return wildType; //must still be wildtype at loci
	else
		return it->second[whichChromosome];
}