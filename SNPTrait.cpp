#include "SNPTrait.h"

// ----------------------------------------------------------------------------------------
// for initialising population
// ----------------------------------------------------------------------------------------
SNPTrait::SNPTrait(SpeciesTrait* P)
{
	pProtoTrait = P;

	DistributionType mutationDistribution = pProtoTrait->getMutationDistribution();
	map<parameter_t, float> mutationParameters = pProtoTrait->getMutationParameters();

	_inherit_func_ptr = (pProtoTrait->getPloidy() == 1) ? &SNPTrait::inheritHaploid : &SNPTrait::inheritDiploid; //this could be changed if we wanted some alternative form of inheritance

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

	if (wildType > maxSNPAlleles - 1)
		cout << endl << "Error:: max number of alleles cannot exceed " << maxSNPAlleles << ".\n";

	DistributionType initialDistribution = pProtoTrait->getInitialDistribution();
	map<parameter_t, float> initialParameters = pProtoTrait->getInitialParameters();

	if (mutationDistribution == SSM && initialDistribution != UNIFORM)
		cout << endl << ("Error:: If using SSM mutation model must initialise genome with alleles (microsats) \n");

	switch (initialDistribution) {
	case UNIFORM:
	{
		if (initialParameters.count(MAX) != 1)
			cout << endl << "Error:: initial SNP/Microsat distribution parameter must contain one max value if set to UNIFORM (e.g. max= ), max cannot exceed " << maxSNPAlleles << "\n";

		float maxD = initialParameters.find(MAX)->second;
		if (maxD > maxSNPAlleles) {
			cout << endl << "Warning:: initial SNP/Microsat distribution parameter max cannot exceed " << maxSNPAlleles << ", resetting to " << maxSNPAlleles << "\n";

			maxD = maxSNPAlleles; //reserve 255 for wildtype
		}
		initialiseFull(maxD);

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
	pProtoTrait(T.pProtoTrait), _mutate_func_ptr(T._mutate_func_ptr), _inherit_func_ptr(T._inherit_func_ptr) {

}


// ----------------------------------------------------------------------------------------
// mutate options
// ----------------------------------------------------------------------------------------

// ----------------------------------------------------------------------------------------
// mutate_KAM
// ----------------------------------------------------------------------------------------
void SNPTrait::mutate_KAM()
{
	const int positionsSize = pProtoTrait->getPositionsSize();
	const auto& positions = pProtoTrait->getPositions();
	const short ploidy = pProtoTrait->getPloidy();
	const float mutationRate = pProtoTrait->getMutationRate();
	auto rng = pRandom->getRNG();
	unsigned char mut;

	map<parameter_t, float> mutationParameters = pProtoTrait->getMutationParameters();


	int maxD = (int)mutationParameters.find(MAX)->second;

	if (maxD > maxSNPAlleles) maxD = maxSNPAlleles; //reserve max value -1 for wildtype

	for (int p = 0; p < ploidy; p++) {

		unsigned int NbMut = pRandom->Poisson(positionsSize * mutationRate);

		if (NbMut > 0) {
			vector<int> mutationPositions;
			sample(positions.begin(), positions.end(), std::back_inserter(mutationPositions),
				NbMut, rng);

			for (int m : mutationPositions) {

				mut = (unsigned char)pRandom->IRandom(0, maxD - 1); //draw new mutation, could draw wildtype
				auto it = genes.find(m); //find if position in map already has mutations there

				if (it == genes.end()) {   // not found so create new entry in map with wildtype as char default

					vector<unsigned char> allelePair(2, wildType);
					allelePair[p] = mut; //put new mutation value in 

					genes.insert(make_pair(m, allelePair));
				}
				else { //position found, already mutations there

					auto currentChar = it->second[p]; //current mutation
					do {
						mut = (unsigned char)pRandom->IRandom(0, maxD - 1); //make sure new value differs from old , could be a problem here with infinite loops
					} while (mut == currentChar);

					it->second[p] = mut; //overwrite with new value
				}
			}
		}
	}
}


// ----------------------------------------------------------------------------------------
// mutate_SSM
// ----------------------------------------------------------------------------------------
void SNPTrait::mutate_SSM()
{
	const int positionsSize = pProtoTrait->getPositionsSize();
	const auto& positions = pProtoTrait->getPositions();
	const short ploidy = pProtoTrait->getPloidy();
	const float mutationRate = pProtoTrait->getMutationRate();
	auto rng = pRandom->getRNG();

	map<parameter_t, float> mutationParameters = pProtoTrait->getMutationParameters();

	int maxD = (int)mutationParameters.find(MAX)->second;
	if (maxD > maxSNPAlleles) maxD = maxSNPAlleles; //reserved max value for wildtype

	for (int p = 0; p < ploidy; p++) {

		unsigned int NbMut = pRandom->Poisson(positionsSize * mutationRate);

		if (NbMut > 0) {
			vector<int> mutationPositions;
			sample(positions.begin(), positions.end(), std::back_inserter(mutationPositions),
				NbMut, rng);

			for (int m : mutationPositions) {
				int mutateForward = pRandom->Bernoulli(0.5);
				auto it = genes.find(m);
				auto currentAllele = it->second[p];//current
				 //alleles values are from 0 to maxD - 1
				if (mutateForward && currentAllele < maxD - 1)
					it->second[p] += 1; //one step to the right
				else if (currentAllele > 0) // !direction || all==_allele_num
					it->second[p] -= 1; //one step to the left
				else //!direction && all == 0
					it->second[p] += 1;
			}
		}
	}
}

// ----------------------------------------------------------------------------------------
// inheritance options
// ----------------------------------------------------------------------------------------

void SNPTrait::inherit(TTrait* parent, set<unsigned int> const& recomPositions, sex_t chromosome, int startingChromosome)
{
	auto parentCast = dynamic_cast<SNPTrait*> (parent); //horrible

	const auto& parent_seq = parentCast->getGenes();
	if (parent_seq.size() > 0) //else nothing to inherit
		(this->*_inherit_func_ptr) (chromosome, parent_seq, recomPositions, startingChromosome);
}



void SNPTrait::inheritDiploid(sex_t whichChromosome, map<int, vector<unsigned char>> const& parentGenes, set<unsigned int> const& recomPositions, int parentChromosome) {

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
					// not found
					vector<unsigned char> newAllelePair(2, wildType);
					newAllelePair[whichChromosome] = sp;
					genes.insert(make_pair(locus, newAllelePair));
				}
				else {
					it->second[whichChromosome] = sp;
				}
			}
		}
	}
}

void SNPTrait::inheritHaploid(sex_t chromosome, map<int, vector<unsigned char>> const& parentGenes, set<unsigned int> const& recomPositions, int parentChromosome)
{
	genes = parentGenes;
}

// ----------------------------------------------------------------------------------------
// Initialise neutral loci
// ----------------------------------------------------------------------------------------

void SNPTrait::initialiseFull(int maxAlleleVal)
{
	const auto& positions = pProtoTrait->getPositions();
	short ploidy = pProtoTrait->getPloidy();

	for (auto position : positions) {
		vector<unsigned char> allelePair;

		for (int i = 0; i < ploidy; i++) {
			auto alleleVal = (unsigned char)pRandom->IRandom(0, maxAlleleVal - 1); //  allele values span 0 - max (255 ceiling) inclusive, max == wildtype
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