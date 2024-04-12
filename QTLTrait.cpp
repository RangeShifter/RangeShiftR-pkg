#include "QTLTrait.h"

// ----------------------------------------------------------------------------------------
// for initialising population
// ----------------------------------------------------------------------------------------

QTLTrait::QTLTrait(SpeciesTrait* P)
{
	pSpeciesTrait = P;
	ExpressionType expressionType = pSpeciesTrait->getExpressionType();

	if (!pSpeciesTrait->isInherited()) //there is a trait for individual variation but this isn't inherited variation it's sampled from initial distribution 
		_inherit_func_ptr = &QTLTrait::inheritInitialParameters;
	else {
		_inherit_func_ptr = (pSpeciesTrait->getPloidy() == 1) ? &QTLTrait::inheritHaploid : &QTLTrait::inheritDiploid; //this could be changed if we wanted some alternative form of inheritance

		DistributionType mutationDistribution = pSpeciesTrait->getMutationDistribution();
		map<GenParamType, float> mutationParameters = pSpeciesTrait->getMutationParameters();

		switch (mutationDistribution) {
		case UNIFORM:
		{
			if (mutationParameters.count(MAX) != 1)
				cout << endl << ("Error:: mutation uniform qtl distribution parameter must contain max value (e.g. max= ) \n");

			if (mutationParameters.count(MIN) != 1)
				cout << endl << ("Error:: mutation uniform qtl distribution parameter must contain min value (e.g. min= ) \n");

			_mutate_func_ptr = &QTLTrait::mutateUniform;
			break;
		}
		case NORMAL:
		{
			if (mutationParameters.count(MEAN) != 1)
				cout << endl << ("Error:: qtl mutation distribution set to normal so parameters must contain mean value (e.g. mean= ) \n");

			if (mutationParameters.count(SD) != 1)
				cout << endl << ("Error::qtl mutation distribution set to normal so parameters must contain sdev value (e.g. sdev= ) \n");

			_mutate_func_ptr = &QTLTrait::mutateNormal;
			break;
		}

		default:
		{
			cout << endl << ("Error:: wrong parameter value for qtl mutation model, must be uniform/normal \n"); //unless want to add gamma or negative exp 
			break; //should return false
		}
		}
	}

	DistributionType initialDistribution = pSpeciesTrait->getInitialDistribution();
	map<GenParamType, float> initialParameters = pSpeciesTrait->getInitialParameters();

	switch (initialDistribution) {
	case UNIFORM:
	{
		if (initialParameters.count(MAX) != 1)
			cout << endl << ("Error:: initial uniform qtl distribution parameter must contain max value (e.g. max= ) \n");

		if (initialParameters.count(MIN) != 1)
			cout << endl << ("Error:: initial uniform qtl distribution parameter must contain min value (e.g. min= ) \n");

		float maxD = initialParameters.find(MAX)->second;
		float minD = initialParameters.find(MIN)->second;

		initialiseUniform(minD, maxD);
		break;
	}
	case NORMAL:
	{
		if (initialParameters.count(MEAN) != 1)
			cout << endl << ("Error:: initial normal qtl distribution parameter must contain mean value (e.g. mean= ) \n");

		if (initialParameters.count(SD) != 1)
			cout << endl << ("Error:: initial normal qtl distribution parameter must contain sdev value (e.g. sdev= ) \n");

		float mean = initialParameters.find(MEAN)->second;
		float sd = initialParameters.find(SD)->second;

		initialiseNormal(mean, sd);
		break;
	}

	default:
	{
		cout << endl << ("wrong parameter value for parameter \"initialisation of qtl\", must be uniform/normal \n");
		break; //should return false
	}
	}

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
		cout << endl << ("wrong parameter value for parameter \"expression of qtl\", must be average/additive \n");
		break; //should return false
	}
	}
}


// ----------------------------------------------------------------------------------------
// for creating new individuals
// ----------------------------------------------------------------------------------------

QTLTrait::QTLTrait(const QTLTrait& T) : pSpeciesTrait(T.pSpeciesTrait), _mutate_func_ptr(T._mutate_func_ptr), _inherit_func_ptr(T._inherit_func_ptr), _express_func_ptr(T._express_func_ptr)
{}



// ----------------------------------------------------------------------------------------
// mutate uniform
// ----------------------------------------------------------------------------------------
void QTLTrait::mutateUniform()
{
	const int positionsSize = pSpeciesTrait->getPositionsSize();
	const auto& genePositions = pSpeciesTrait->getGenePositions();
	const short ploidy = pSpeciesTrait->getPloidy();
	const float mutationRate = pSpeciesTrait->getMutationRate();

	auto rng = pRandom->getRNG();

	map<GenParamType, float> mutationParameters = pSpeciesTrait->getMutationParameters();
	float maxD = mutationParameters.find(MAX)->second;
	float minD = mutationParameters.find(MIN)->second;

	for (int p = 0; p < ploidy; p++) {

		unsigned int NbMut = pRandom->Poisson(positionsSize * mutationRate);

		if (NbMut > 0) {
			vector<int> mutationPositions;
			sample(genePositions.begin(), genePositions.end(), std::back_inserter(mutationPositions),
				NbMut, rng);

			for (int m : mutationPositions) {
				auto it = genes.find(m);
				float currentAlleleVal = it->second[p].get()->getAlleleValue();//current
				float newAlleleVal = pRandom->FRandom(minD, maxD) + currentAlleleVal;
				it->second[p] = make_shared<Allele>(newAlleleVal, 1.0);
			}
		}
	}
}

// ----------------------------------------------------------------------------------------
// mutate normal
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

	for (int p = 0; p < ploidy; p++) {

		unsigned int NbMut = pRandom->Poisson(positionsSize * mutationRate);

		if (NbMut > 0) {
			vector<int> mutationPositions;
			sample(genePositions.begin(), genePositions.end(), std::back_inserter(mutationPositions),
				NbMut, rng);

			for (int m : mutationPositions) {
				auto it = genes.find(m);
				float currentAlleleVal = it->second[p].get()->getAlleleValue();//current
				float newAlleleVal = pRandom->Normal(mean, sd) + currentAlleleVal;
				it->second[p] = make_shared<Allele>(newAlleleVal, 1.0);
			}
		}
	}
}

// ----------------------------------------------------------------------------------------
// inheritance options
// ----------------------------------------------------------------------------------------


void QTLTrait::inherit(const bool& fromMother, TTrait* parentTrait, set<unsigned int> const& recomPositions, int startingChromosome)
{
	auto parentCast = dynamic_cast<QTLTrait*> (parentTrait); //horrible

	const auto& parent_seq = parentCast->getGenes();
	if (parent_seq.size() > 0) //else nothing to inherit, should always be something to inherit with QTL
		(this->*_inherit_func_ptr) (fromMother, parent_seq, recomPositions, startingChromosome);
}

void QTLTrait::inheritDiploid(const bool& fromMother, map<int, vector<shared_ptr<Allele>>> const& parentGenes, set<unsigned int> const& recomPositions, int parentChromosome) {

	auto it = recomPositions.lower_bound(parentGenes.begin()->first);
	int nextBreakpoint = *it;
	auto distance = std::distance(recomPositions.begin(), it);
	if (distance % 2 != 0)
		parentChromosome = 1 - parentChromosome; //switch chromosome

	for (auto const& [locus, allelePair] : parentGenes) {

		while (locus > nextBreakpoint) {
			std::advance(it, 1);
			nextBreakpoint = *it;
			parentChromosome = 1 - parentChromosome; //switch chromosome
		}

		if (locus <= nextBreakpoint) {
			auto& sp = allelePair[parentChromosome];
			auto it = genes.find(locus);
			if (it == genes.end()) {
				// locus does not exist yet, initiate it
				if (!fromMother) throw runtime_error("Father-inherited locus does not exist.");
				vector<shared_ptr<Allele>> newAllelePair(2);
				newAllelePair[sex_t::FEM] = allele;
				genes.insert(make_pair(locus, newAllelePair));
			}
			else { // father, locus already exists
				if (fromMother) throw runtime_error("Mother-inherited locus already exists.");
				it->second[sex_t::MAL] = allele;
			}
		}
	}
}

void QTLTrait::inheritHaploid(const bool& fromMother, map<int, vector<shared_ptr<Allele>>> const& parentGenes, set<unsigned int> const& recomPositions, int parentChromosome)
{
	genes = parentGenes;
}

// ----------------------------------------------------------------------------------------
// 'Inherit' from initialisation parameters, for simulations with individual variation but no inheritance
// ----------------------------------------------------------------------------------------

void QTLTrait::inheritInitialParameters(sex_t whichChromosome, map<int, vector<shared_ptr<Allele>>> const& parentGenes, set<unsigned int> const& recomPositions, int parentChromosome)
{
	DistributionType initialDistribution = pSpeciesTrait->getInitialDistribution();
	map<GenParamType, float> initialParameters = pSpeciesTrait->getInitialParameters();

	switch (initialDistribution) {
	case UNIFORM:
	{
		if (initialParameters.count(MAX) != 1)
			cout << endl << ("Error:: initial uniform qtl distribution parameter must contain max value (e.g. max= ) \n");

		if (initialParameters.count(MIN) != 1)
			cout << endl << ("Error:: initial uniform qtl distribution parameter must contain min value (e.g. min= ) \n");

		float maxD = initialParameters.find(MAX)->second;
		float minD = initialParameters.find(MIN)->second;

		initialiseUniform(minD, maxD);

		break;
	}
	case NORMAL:
	{
		if (initialParameters.count(MEAN) != 1)
			cout << endl << ("Error:: initial normal qtl distribution parameter must contain mean value (e.g. mean= ) \n");

		if (initialParameters.count(SD) != 1)
			cout << endl << ("Error:: initial normal qtl distribution parameter must contain sdev value (e.g. sdev= ) \n");

		float mean = initialParameters.find(MEAN)->second;
		float sd = initialParameters.find(SD)->second;

		initialiseNormal(mean, sd);

		break;
	}

	default:
	{
		cout << endl << ("wrong parameter value for parameter \"initialisation of qtl\", must be uniform/normal \n");
		break; //should return false
	}
	}
}


// ----------------------------------------------------------------------------------------
// Initialisation options
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
// expression options
// ----------------------------------------------------------------------------------------

float QTLTrait::expressAdditive() {

	float phenotype = 0.0;

	for (auto const& [locus, allelePair] : genes)
	{
		for (const std::shared_ptr<Allele> m : allelePair)
			phenotype += m->getAlleleValue();
	}
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
	return phenotype;
}

// ----------------------------------------------------------------------------------------
// check if particular loci is heterozygote
// ----------------------------------------------------------------------------------------

bool QTLTrait::isHeterozygoteAtLocus(int locus) const {

	auto it = genes.find(locus);

	if (it == genes.end()) //not found
		return false;
	else {
		if (it->second.size() == 1) //only one ptr there so must be allele and wildtype at loci (heterozygote)
			return true;
		else
			return(it->second[0].get()->getId() != it->second[1].get()->getId());
	}
}

// ----------------------------------------------------------------------------------------
// count heterozygote loci in genome 
// ----------------------------------------------------------------------------------------

int QTLTrait::countHeterozygoteLoci() const {

	int count = 0;
	for (auto const& [locus, allelePair] : genes) {

		if (allelePair.size() == 1) //only one ptr there so must be allele and wildtype at loci (heterozygote)
			count++;

		if (allelePair.size() == 2)
			count += (allelePair[0].get()->getId() != allelePair[1].get()->getId());
	}
	return count;
}

// ----------------------------------------------------------------------------------------
// get allele value at loci 
// ----------------------------------------------------------------------------------------

float QTLTrait::getAlleleValueAtLocus(short whichChromosome, int position) const {

	auto it = genes.find(position);

	if (it == genes.end()) { //no mutations there, should never happen at QTLs should always hold a value 
		return 0; //must still be wildtype at loci
		cout << endl << ("Error:: trying to find QTL at ", position, " but doesn't exist \n");
	}
	else
		return it->second[whichChromosome].get()->getAlleleValue();

}