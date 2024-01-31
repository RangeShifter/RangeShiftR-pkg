#include "QTLTrait.h"

// ----------------------------------------------------------------------------------------
// for initialising population
// ----------------------------------------------------------------------------------------

QTLTrait::QTLTrait(ProtoTrait* P)
{
	pProtoTrait = P;
	ExpressionType expressionType = pProtoTrait->getExpressionType();

	if (!pProtoTrait->isInherited()) //there is a trait for individual variation but this isn't inherited variation it's sampled from initial distribution 
		_inherit_func_ptr = &QTLTrait::inheritInitialParameters;
	else {
		_inherit_func_ptr = (pProtoTrait->getPloidy() == 1) ? &QTLTrait::inheritHaploid : &QTLTrait::inheritDiploid; //this could be changed if we wanted some alternative form of inheritance

		DistributionType mutationDistribution = pProtoTrait->getMutationDistribution();
		map<parameter_t, float> mutationParameters = pProtoTrait->getMutationParameters();

		switch (mutationDistribution) {
		case UNIFORM:
		{
			if (!mutationParameters.count(MAX))
				cout << endl << ("Error:: mutation uniform qtl distribution parameter must contain max value (e.g. max= ) \n");

			if (!mutationParameters.count(MIN))
				cout << endl << ("Error:: mutation uniform qtl distribution parameter must contain min value (e.g. min= ) \n");

			_mutate_func_ptr = &QTLTrait::mutateUniform;
			break;
		}
		case NORMAL:
		{
			if (!mutationParameters.count(MEAN))
				cout << endl << ("Error:: qtl mutation distribution set to normal so parameters must contain mean value (e.g. mean= ) \n");

			if (!mutationParameters.count(SDEV))
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

	DistributionType initialDistribution = pProtoTrait->getInitialDistribution();
	map<parameter_t, float> initialParameters = pProtoTrait->getInitialParameters();

	switch (initialDistribution) {
	case UNIFORM:
	{
		if (!initialParameters.count(MAX))
			cout << endl << ("Error:: initial uniform qtl distribution parameter must contain max value (e.g. max= ) \n");

		if (!initialParameters.count(MIN))
			cout << endl << ("Error:: initial uniform qtl distribution parameter must contain min value (e.g. min= ) \n");

		float maxD = initialParameters.find(MAX)->second;
		float minD = initialParameters.find(MIN)->second;

		initialiseUniform(minD, maxD);
		break;
	}
	case NORMAL:
	{
		if (!initialParameters.count(MEAN))
			cout << endl << ("Error:: initial normal qtl distribution parameter must contain mean value (e.g. mean= ) \n");

		if (!initialParameters.count(SDEV))
			cout << endl << ("Error:: initial normal qtl distribution parameter must contain sdev value (e.g. sdev= ) \n");

		float mean = initialParameters.find(MEAN)->second;
		float sd = initialParameters.find(SDEV)->second;

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

QTLTrait::QTLTrait(const QTLTrait& T) : pProtoTrait(T.pProtoTrait), _mutate_func_ptr(T._mutate_func_ptr), _inherit_func_ptr(T._inherit_func_ptr), _express_func_ptr(T._express_func_ptr)
{}



// ----------------------------------------------------------------------------------------
// mutate uniform
// ----------------------------------------------------------------------------------------
void QTLTrait::mutateUniform()
{
	const int positionsSize = pProtoTrait->getPositionsSize();
	const auto& positions = pProtoTrait->getPositions();
	const short ploidy = pProtoTrait->getPloidy();
	const float mutationRate = pProtoTrait->getMutationRate();

	auto rng = pRandom->getRNG();

	map<parameter_t, float> mutationParameters = pProtoTrait->getMutationParameters();
	float maxD = mutationParameters.find(MAX)->second;
	float minD = mutationParameters.find(MIN)->second;

	for (int p = 0; p < ploidy; p++) {

		unsigned int NbMut = pRandom->Poisson(positionsSize * mutationRate);

		if (NbMut > 0) {
			vector<int> mutationPositions;
			sample(positions.begin(), positions.end(), std::back_inserter(mutationPositions),
				NbMut, rng);

			for (int m : mutationPositions) {

				auto it = mutations.find(m);
				float currentSelCoef = it->second[p].get()->getSelectionCoef();//current
				float newSelectionCoef = pRandom->FRandom(minD, maxD) + currentSelCoef;
				it->second[p] = make_shared<Mutation>(newSelectionCoef, 1.0);

			}
		}

	}

}

// ----------------------------------------------------------------------------------------
// mutate normal
// ----------------------------------------------------------------------------------------

void QTLTrait::mutateNormal()
{

	const int positionsSize = pProtoTrait->getPositionsSize();
	const auto& positions = pProtoTrait->getPositions();
	const short ploidy = pProtoTrait->getPloidy();
	const float mutationRate = pProtoTrait->getMutationRate();

	auto rng = pRandom->getRNG();


	const map<parameter_t, float> mutationParameters = pProtoTrait->getMutationParameters();
	const float mean = mutationParameters.find(MEAN)->second;
	const float sd = mutationParameters.find(SDEV)->second;


	for (int p = 0; p < ploidy; p++) {

		unsigned int NbMut = pRandom->Poisson(positionsSize * mutationRate);

		if (NbMut > 0) {
			vector<int> mutationPositions;
			sample(positions.begin(), positions.end(), std::back_inserter(mutationPositions),
				NbMut, rng);

			for (int m : mutationPositions) {

				auto it = mutations.find(m);
				float currentSelCoef = it->second[p].get()->getSelectionCoef();//current
				float newSelectionCoef = pRandom->Normal(mean, sd) + currentSelCoef;
				it->second[p] = make_shared<Mutation>(newSelectionCoef, 1.0);
			}
		}


	}

}

// ----------------------------------------------------------------------------------------
// inheritance options
// ----------------------------------------------------------------------------------------


void QTLTrait::inherit(TTrait* parent, set<unsigned int> const& recomPositions, sex_t chromosome, int startingChromosome)
{

	auto parentCast = dynamic_cast<QTLTrait*> (parent); //horrible

	const auto& parent_seq = parentCast->get_mutations();
	if (parent_seq.size() > 0) //else nothing to inherit, should always be something to inherit with QTL
		(this->*_inherit_func_ptr) (chromosome, parent_seq, recomPositions, startingChromosome);


}

void QTLTrait::inheritDiploid(sex_t chromosome, map<int, vector<shared_ptr<Mutation>>>  const& parentMutations, set<unsigned int> const& recomPositions, int parentChromosome) {

	auto it = recomPositions.lower_bound(parentMutations.begin()->first);

	unsigned int nextBreakpoint = *it;

	auto distance = std::distance(recomPositions.begin(), it);
	if (distance % 2 != 0)
		parentChromosome = !parentChromosome; //switch chromosome

	for (auto const& [key, val] : parentMutations) {

		while (key > nextBreakpoint) {
			std::advance(it, 1);
			nextBreakpoint = *it;
			parentChromosome = !parentChromosome; //switch chromosome
		}

		if (key <= nextBreakpoint) {
			auto& sp = val[parentChromosome];


			auto it = mutations.find(key);
			if (it == mutations.end()) {
				// not found
				vector<shared_ptr<Mutation>> vect(2);
				vect[chromosome] = sp;

				mutations.insert(make_pair(key, vect));
			}
			else {
				it->second[chromosome] = sp;


			}



		}



	}





}

void QTLTrait::inheritHaploid(sex_t chromosome, map<int, vector<shared_ptr<Mutation>>> const& parentMutations, set<unsigned int> const& recomPositions, int parentChromosome)
{
	mutations = parentMutations;
}

// ----------------------------------------------------------------------------------------
// 'Inherit' from initialisation parameters, for simulations with individual variation but no inheritance
// ----------------------------------------------------------------------------------------

void QTLTrait::inheritInitialParameters(sex_t chromosome, map<int, vector<shared_ptr<Mutation>>> const& parentMutations, set<unsigned int> const& recomPositions, int parentChromosome)
{
	DistributionType initialDistribution = pProtoTrait->getInitialDistribution();
	map<parameter_t, float> initialParameters = pProtoTrait->getInitialParameters();


	switch (initialDistribution) {
	case UNIFORM:
	{
		if (!initialParameters.count(MAX))
			cout << endl << ("Error:: initial uniform qtl distribution parameter must contain max value (e.g. max= ) \n");

		if (!initialParameters.count(MIN))
			cout << endl << ("Error:: initial uniform qtl distribution parameter must contain min value (e.g. min= ) \n");

		float maxD = initialParameters.find(MAX)->second;
		float minD = initialParameters.find(MIN)->second;

		initialiseUniform(minD, maxD);

		break;
	}
	case NORMAL:
	{
		if (!initialParameters.count(MEAN))
			cout << endl << ("Error:: initial normal qtl distribution parameter must contain mean value (e.g. mean= ) \n");

		if (!initialParameters.count(SDEV))
			cout << endl << ("Error:: initial normal qtl distribution parameter must contain sdev value (e.g. sdev= ) \n");

		float mean = initialParameters.find(MEAN)->second;
		float sd = initialParameters.find(SDEV)->second;

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

	const auto positions = pProtoTrait->getPositions();
	short ploidy = pProtoTrait->getPloidy();

	for (auto position : positions) {
		vector<shared_ptr<Mutation>> vect;

		for (int i = 0; i < ploidy; i++) {
			float selectionCoef = pRandom->Normal(mean, sd);
			vect.emplace_back(make_shared<Mutation>(selectionCoef, 1.0));//crude here to set dominance to 1 for qtl but can be changed later
		}
		mutations.insert(make_pair(position, vect));
	}
}

void QTLTrait::initialiseUniform(float min, float max) {

	const auto positions = pProtoTrait->getPositions();
	short ploidy = pProtoTrait->getPloidy();

	for (auto position : positions) {
		vector<shared_ptr<Mutation>> vect;

		for (int i = 0; i < ploidy; i++) {
			float selectionCoef = pRandom->FRandom(min, max);
			vect.emplace_back(make_shared<Mutation>(selectionCoef, 1.0)); //crude here to set dominance to 1 for qtl but can be changed later
		}
		mutations.insert(make_pair(position, vect));
	}
}

// ----------------------------------------------------------------------------------------
// expression options
// ----------------------------------------------------------------------------------------

float QTLTrait::expressAdditive() {

	float phenotype = 0.0;

	for (auto const& [key, val] : mutations)
	{
		for (auto m : val)
			phenotype += m->getSelectionCoef();
	}

	return phenotype;
}

float QTLTrait::expressAverage() {

	int positionsSize = pProtoTrait->getPositionsSize();
	short ploidy = pProtoTrait->getPloidy();
	float phenotype = 0.0;

	for (auto const& [key, val] : mutations)
	{
		for (auto& m : val)
			phenotype += m->getSelectionCoef();
	}
	phenotype = phenotype / (positionsSize * ploidy);
	return phenotype;
}

// ----------------------------------------------------------------------------------------
// check if particular loci is heterozygote
// ----------------------------------------------------------------------------------------


bool QTLTrait::isHeterozygoteAtLoci(int loci) const {

	auto it = mutations.find(loci);

	if (it == mutations.end()) //not found
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

	for (auto const& [key, val] : mutations) {

		if (val.size() == 1) //only one ptr there so must be allele and wildtype at loci (heterozygote)
			count++;

		if (val.size() == 2)
			count += (val[0].get()->getId() != val[1].get()->getId());
	}

	return count;
}

// ----------------------------------------------------------------------------------------
// get allele value at loci 
// ----------------------------------------------------------------------------------------


float QTLTrait::getSelectionCoefAtLoci(short chromosome, int position) const {

	auto it = mutations.find(position);

	if (it == mutations.end()) { //no mutations there, should never happen at QTLs should always hold a value 
		return 0; //must still be wildtype at loci
		cout << endl << ("Error:: trying to find QTL at ", position, " but doesn't exist \n");
	}
	else
		return it->second[chromosome].get()->getSelectionCoef();

}