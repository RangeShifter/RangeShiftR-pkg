#include "QTLTrait.h"

// ----------------------------------------------------------------------------------------
//  Initialisation constructor
//  Called when initialising community
//  Sets up initial values, and immutable attributes (distributions and parameters)
//	that are defined at the species-level
// ----------------------------------------------------------------------------------------
QTLTrait::QTLTrait(SpeciesTrait* P)
{
	pSpeciesTrait = P;
	ExpressionType expressionType = pSpeciesTrait->getExpressionType();

	if (!pSpeciesTrait->isInherited()) // there is a trait for individual variation but this isn't inherited variation it's sampled from initial distribution 
		_inherit_func_ptr = &QTLTrait::reInitialiseGenes;
	else {
		_inherit_func_ptr = (pSpeciesTrait->getPloidy() == 1) ? &QTLTrait::inheritHaploid : &QTLTrait::inheritDiploid; //this could be changed if we wanted some alternative form of inheritance

		DistributionType mutationDistribution = pSpeciesTrait->getMutationDistribution();
		map<GenParamType, float> mutationParameters = pSpeciesTrait->getMutationParameters();

		// Set mutation parameters
		switch (mutationDistribution) {
		case UNIFORM:
		{
			if (mutationParameters.count(MAX) != 1)
				throw logic_error("Error:: mutation uniform qtl distribution parameter must contain max value (e.g. max= ) \n");
			if (mutationParameters.count(MIN) != 1)
				throw logic_error("Error:: mutation uniform qtl distribution parameter must contain min value (e.g. min= ) \n");
			_mutate_func_ptr = &QTLTrait::mutateUniform;
			break;
		}
		case NORMAL:
		{
			if (mutationParameters.count(MEAN) != 1)
				throw logic_error("Error:: qtl mutation distribution set to normal so parameters must contain mean value (e.g. mean= ) \n");
			if (mutationParameters.count(SD) != 1)
				throw logic_error("Error::qtl mutation distribution set to normal so parameters must contain sdev value (e.g. sdev= ) \n");
			_mutate_func_ptr = &QTLTrait::mutateNormal;
			break;
		}
		default:
		{
			throw logic_error("Error:: wrong parameter value for qtl mutation model, must be uniform/normal \n"); //unless want to add gamma or negative exp 
			break;
		}
		}
	}

	// Set initialisation parameters
	DistributionType initialDistribution = pSpeciesTrait->getInitialDistribution();
	map<GenParamType, float> initialParameters = pSpeciesTrait->getInitialParameters();
	switch (initialDistribution) {
	case UNIFORM:
	{
		if (initialParameters.count(MAX) != 1)
			throw logic_error("Error:: initial uniform qtl distribution parameter must contain max value (e.g. max= ) \n");
		if (initialParameters.count(MIN) != 1)
			throw logic_error("Error:: initial uniform qtl distribution parameter must contain min value (e.g. min= ) \n");
		float maxD = initialParameters.find(MAX)->second;
		float minD = initialParameters.find(MIN)->second;
		initialiseUniform(minD, maxD);
		break;
	}
	case NORMAL:
	{
		if (initialParameters.count(MEAN) != 1)
			throw logic_error("Error:: initial normal qtl distribution parameter must contain mean value (e.g. mean= ) \n");
		if (initialParameters.count(SD) != 1)
			throw logic_error("Error:: initial normal qtl distribution parameter must contain sdev value (e.g. sdev= ) \n");
		float mean = initialParameters.find(MEAN)->second;
		float sd = initialParameters.find(SD)->second;
		initialiseNormal(mean, sd);
		break;
	}
	default:
	{
		throw logic_error("wrong parameter value for parameter \"initialisation of qtl\", must be uniform/normal \n");
		break;
	}
	}

	// Set expression mode parameters
	switch (expressionType) {
	case AVERAGE:
	{
		_express_func_ptr = &QTLTrait::expressAverage;
		break;
	}
	case ADDITIVE:
	{
		_express_func_ptr = &QTLTrait::expressAdditive;
		break;
	}
	default:
	{
		throw logic_error("wrong parameter value for parameter \"expression of qtl\", must be average/additive \n");
		break;
	}
	}
}

// ----------------------------------------------------------------------------------------
// Inheritance constructor
// Copies immutable features from a parent trait
// Only called via clone()
// ----------------------------------------------------------------------------------------
QTLTrait::QTLTrait(const QTLTrait& T) : pSpeciesTrait(T.pSpeciesTrait), _mutate_func_ptr(T._mutate_func_ptr), _inherit_func_ptr(T._inherit_func_ptr), _express_func_ptr(T._express_func_ptr)
{}

// ----------------------------------------------------------------------------------------
// Sample and apply mutations from a uniform distribution
// 
// Mutations drawn only for existing positions, 
// that is no new genes are created during simulation
// ----------------------------------------------------------------------------------------
void QTLTrait::mutateUniform()
{
	const int positionsSize = pSpeciesTrait->getPositionsSize();
	const auto& genePositions = pSpeciesTrait->getGenePositions();
	const short ploidy = pSpeciesTrait->getPloidy();
	const float mutationRate = pSpeciesTrait->getMutationRate();
	float newAlleleVal;

	auto rng = pRandom->getRNG();

	map<GenParamType, float> mutationParameters = pSpeciesTrait->getMutationParameters();
	float maxD = mutationParameters.find(MAX)->second;
	float minD = mutationParameters.find(MIN)->second;

	for (int p = 0; p < ploidy; p++) {

		unsigned int NbMut = pRandom->Binomial(positionsSize, mutationRate);

		if (NbMut > 0) {
			vector<int> mutationPositions;
			sample(genePositions.begin(), genePositions.end(), std::back_inserter(mutationPositions),
				NbMut, rng);

			for (int m : mutationPositions) {
				auto it = genes.find(m);
				if (it == genes.end())
					throw runtime_error("Locus sampled for mutation doesn't exist.");
				float currentAlleleVal = it->second[p].get()->getAlleleValue();//current
				newAlleleVal = pRandom->FRandom(minD, maxD) + currentAlleleVal;
				it->second[p] = make_shared<Allele>(newAlleleVal, QTLDominanceFactor);
			}
		}
	}
}

// ----------------------------------------------------------------------------------------
// Sample and apply mutations from a normal distribution
// Mutations drawn only for existing positions, 
// that is no new genes are created during simulation
// ----------------------------------------------------------------------------------------
void QTLTrait::mutateNormal()
{
	const int positionsSize = pSpeciesTrait->getPositionsSize();
	const auto& genePositions = pSpeciesTrait->getGenePositions();
	const short ploidy = pSpeciesTrait->getPloidy();
	const float mutationRate = pSpeciesTrait->getMutationRate();

	auto rng = pRandom->getRNG();

	const map<GenParamType, float> mutationParameters = pSpeciesTrait->getMutationParameters();
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
				auto it = genes.find(m);
				if (it == genes.end())
					throw runtime_error("Locus sampled for mutation doesn't exist.");
				float currentAlleleVal = it->second[p].get()->getAlleleValue(); //current
				newAlleleVal = pRandom->Normal(mean, sd) + currentAlleleVal;
				it->second[p] = make_shared<Allele>(newAlleleVal, QTLDominanceFactor);
			}
		}
	}
}

// ----------------------------------------------------------------------------------------
//  Wrapper to inheritance function
// ----------------------------------------------------------------------------------------
void QTLTrait::inheritGenes(const bool& fromMother, TTrait* parentTrait, set<unsigned int> const& recomPositions, int startingChromosome)
{
	auto parentCast = dynamic_cast<QTLTrait*>(parentTrait); // must convert TTrait to QTLTrait
	const auto& parent_seq = parentCast->getGenes();
	(this->*_inherit_func_ptr)(fromMother, parent_seq, recomPositions, startingChromosome);
}

// ----------------------------------------------------------------------------------------
// Inheritance for diploid, sexual species
// Called once for each parent. Given a list of recombinant sites, 
// populates offspring genes with appropriate parent alleles
// Assumes mother genes are inherited first
// ----------------------------------------------------------------------------------------
void QTLTrait::inheritDiploid(const bool& fromMother, map<int, vector<shared_ptr<Allele>>> const& parentGenes, set<unsigned int> const& recomPositions, int parentChromosome) {

	auto it = recomPositions.lower_bound(parentGenes.begin()->first);
	int nextBreakpoint = *it;
	// Is the first parent gene position already recombinant?
	auto distance = std::distance(recomPositions.begin(), it);
	if (distance % 2 != 0) // odd positions = switch, even = switch back
		parentChromosome = 1 - parentChromosome; //switch chromosome

	for (auto const& [locus, allelePair] : parentGenes) {

		// Switch chromosome if locus is past recombination site
		while (locus > nextBreakpoint) {
			parentChromosome = 1 - parentChromosome;
			std::advance(it, 1); // go to next recombination site
			nextBreakpoint = *it;
		}

		if (locus <= nextBreakpoint) {
			auto& parentAllele = allelePair[parentChromosome];
			auto itGene = genes.find(locus);
			if (itGene == genes.end()) {
				// locus does not exist yet, create and initialise it
				if (!fromMother) throw runtime_error("Father-inherited locus does not exist.");
				vector<shared_ptr<Allele>> newAllelePair(2);
				newAllelePair[sex_t::FEM] = parentAllele;
				genes.insert(make_pair(locus, newAllelePair));
			}
			else { // father, locus already exists
				if (fromMother) throw runtime_error("Mother-inherited locus already exists.");
				itGene->second[sex_t::MAL] = parentAllele;
			}
		}
	}
}

// ----------------------------------------------------------------------------------------
// Inheritance for haploid, asexual species
// Simply pass down parent genes
// Arguments are still needed to match overloaded function in base class
// ----------------------------------------------------------------------------------------
void QTLTrait::inheritHaploid(const bool& fromMother, map<int, vector<shared_ptr<Allele>>> const& parentGenes, set<unsigned int> const& recomPositions, int parentChromosome)
{
	genes = parentGenes;
}

// ----------------------------------------------------------------------------------------
// Non-inheritance 
// For cases where isInherited option is turned off
// In this case, offspring alleles are populated using the initialise functions
// Arguments are still needed to match overloaded function in base class
// ----------------------------------------------------------------------------------------
void QTLTrait::reInitialiseGenes(const bool& fromMother, map<int, vector<shared_ptr<Allele>>> const& parentGenes, set<unsigned int> const& recomPositions, int parentChromosome)
{
	DistributionType initialDistribution = pSpeciesTrait->getInitialDistribution();
	map<GenParamType, float> initialParameters = pSpeciesTrait->getInitialParameters();

	switch (initialDistribution) {
	case UNIFORM:
	{
		if (initialParameters.count(MAX) != 1)
			throw logic_error("Error:: initial uniform qtl distribution parameter must contain max value (e.g. max= ) \n");
		if (initialParameters.count(MIN) != 1)
			throw logic_error("Error:: initial uniform qtl distribution parameter must contain min value (e.g. min= ) \n");
		float maxD = initialParameters.find(MAX)->second;
		float minD = initialParameters.find(MIN)->second;
		initialiseUniform(minD, maxD);
		break;
	}
	case NORMAL:
	{
		if (initialParameters.count(MEAN) != 1)
			throw logic_error("Error:: initial normal qtl distribution parameter must contain mean value (e.g. mean= ) \n");
		if (initialParameters.count(SD) != 1)
			throw logic_error("Error:: initial normal qtl distribution parameter must contain sdev value (e.g. sdev= ) \n");
		float mean = initialParameters.find(MEAN)->second;
		float sd = initialParameters.find(SD)->second;
		initialiseNormal(mean, sd);
		break;
	}
	default:
	{
		throw logic_error("wrong parameter value for parameter \"initialisation of qtl\", must be uniform/normal \n");
		break; //should return false
	}
	}
}

// ----------------------------------------------------------------------------------------
// QTL initialisation options
// ----------------------------------------------------------------------------------------
void QTLTrait::initialiseNormal(float mean, float sd) {

	const set<int> genePositions = pSpeciesTrait->getGenePositions();
	short ploidy = pSpeciesTrait->getPloidy();

	for (auto position : genePositions) {
		vector<shared_ptr<Allele>> newAllelePair;
		for (int i = 0; i < ploidy; i++) {
			float alleleVal = pRandom->Normal(mean, sd);
			newAllelePair.emplace_back(make_shared<Allele>(alleleVal, QTLDominanceFactor));
		}
		genes.insert(make_pair(position, newAllelePair));
	}
}

void QTLTrait::initialiseUniform(float min, float max) {

	const set<int> genePositions = pSpeciesTrait->getGenePositions();
	short ploidy = pSpeciesTrait->getPloidy();

	for (auto position : genePositions) {
		vector<shared_ptr<Allele>> newAllelePair;
		for (int i = 0; i < ploidy; i++) {
			float alleleVal = pRandom->FRandom(min, max);
			newAllelePair.emplace_back(make_shared<Allele>(alleleVal, QTLDominanceFactor));
		}
		genes.insert(make_pair(position, newAllelePair));
	}
}

// ----------------------------------------------------------------------------------------
// QTL gene expression options
// ----------------------------------------------------------------------------------------
float QTLTrait::expressAdditive() {

	float phenotype = 0.0;

	for (auto const& [locus, allelePair] : genes)
	{
		for (const std::shared_ptr<Allele> m : allelePair)
			phenotype += m->getAlleleValue();
	}
	trimQTLPhenotype(phenotype);
	return phenotype;
}

float QTLTrait::expressAverage() {

	int positionsSize = pSpeciesTrait->getPositionsSize();
	short ploidy = pSpeciesTrait->getPloidy();
	float phenotype = 0.0;

	for (auto const& [locus, allelePair] : genes)
	{
		for (auto& m : allelePair)
			phenotype += m->getAlleleValue();
	}
	phenotype /= positionsSize * ploidy;
	trimQTLPhenotype(phenotype);
	return phenotype;
}

void QTLTrait::trimQTLPhenotype(float& val) {
	const float minPositiveVal = 1e-06;
	switch (pSpeciesTrait->getTraitType())
	{
	// Values bound between 0 and 1
	case E_D0_F: case E_D0_M: 
	case S_S0_F: case S_S0_M:
	case KERNEL_PROBABILITY_F: case KERNEL_PROBABILITY_M:
	case CRW_STEPCORRELATION:
	{
		if (val < 0.0) val = 0;
		else if (val > 1.0) val = 1.0;
		break;
	}
	// Positive values
	case KERNEL_MEANDIST_1_F: case KERNEL_MEANDIST_1_M:
	case KERNEL_MEANDIST_2_F: case KERNEL_MEANDIST_2_M:
	case CRW_STEPLENGTH:
	{
		if (val < 0.0) val = 0;
		break;
	}
	// Strictly positive values
	case E_ALPHA_F: case E_ALPHA_M: 
	case S_ALPHA_F: case S_ALPHA_M:
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
	case E_BETA_F: case E_BETA_M:
	case S_BETA_F: case S_BETA_M:
	case SMS_BETADB:
	{
		break;
	}
	default:
		break;
	}
}

// ----------------------------------------------------------------------------------------
// Check if specific locus is heterozygote
// ----------------------------------------------------------------------------------------
bool QTLTrait::isHeterozygoteAtLocus(int locus) const {
	// assumes diploidy
	auto it = genes.find(locus);
	if (it == genes.end()) //not found
		throw runtime_error("QTL gene queried for heterozygosity does not exist.");
	else
		return(it->second[0].get()->getId() != it->second[1].get()->getId());
}

// ----------------------------------------------------------------------------------------
// Count heterozygote loci in genome 
// ----------------------------------------------------------------------------------------
int QTLTrait::countHeterozygoteLoci() const {
	// assumes diploidy
	int count = 0;
	for (auto const& [locus, allelePair] : genes) {
		if (allelePair.size() == 2)
			count += (allelePair[0].get()->getId() != allelePair[1].get()->getId());
	}
	return count;
}

// ----------------------------------------------------------------------------------------
// Get allele value at locus
// ----------------------------------------------------------------------------------------
float QTLTrait::getAlleleValueAtLocus(short whichChromosome, int position) const {

	auto it = genes.find(position);
	if (it == genes.end())
		throw runtime_error("The QTL locus queried for its allele value does not exist.");
	return it->second[whichChromosome].get()->getAlleleValue();
}