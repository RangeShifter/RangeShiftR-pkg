#include "GeneticLoad.h"

// ----------------------------------------------------------------------------------------
// for initialising population
// ----------------------------------------------------------------------------------------

GeneticLoad::GeneticLoad(ProtoTrait* P)
{

	pProtoTrait = P;

	if (wildType.get() == nullptr)
		wildType = make_shared<Mutation>(0.0, 0.0);

	ExpressionType expressionType = pProtoTrait->getExpressionType();

	_inherit_func_ptr = (pProtoTrait->getPloidy() == 1) ? &GeneticLoad::inheritHaploid : &GeneticLoad::inheritDiploid; //this could be changed if we wanted some alternative form of inheritance

	DistributionType mutationDistribution = pProtoTrait->getMutationDistribution();
	map<parameter_t, float> mutationParameters = pProtoTrait->getMutationParameters();

	switch (mutationDistribution) {
	case UNIFORM:
	{
		if (!mutationParameters.count(MAX))
			cout << endl << ("Error:: adaptive mutation uniform distribution parameter must contain max value (e.g. max= ) \n");

		if (!mutationParameters.count(MIN))
			cout << endl << ("Error:: adaptive mutation uniform distribution parameter must contain min value (e.g. min= ) \n");

		break;
	}
	case NORMAL:
	{

		if (!mutationParameters.count(MEAN))
			cout << endl << ("Error:: adaptive mutation distribution set to normal so parameters must contain mean value (e.g. mean= ) \n");

		if (!mutationParameters.count(SDEV))
			cout << endl << ("Error:: adaptive mutation distribution set to normal so parameters must contain sdev value (e.g. sdev= ) \n");

		break;
	}
	case GAMMA:
	{
		if (!mutationParameters.count(SHAPE))
			cout << endl << ("Error:: adaptive mutation distribution set to gamma so parameters must contain shape value (e.g. shape= ) \n");

		if (!mutationParameters.count(SCALE))
			cout << endl << ("Error:: adaptive mutation distribution set to gamma so parameters must contain scale value (e.g. scale= ) \n");

		break;
	}
	case NEGEXP:
	{
		if (!mutationParameters.count(MEAN))
			cout << endl << ("Error:: adaptive mutation distribution set to negative exponential (negative decay) so parameters must contain mean value (e.g. mean= ) \n");


		break;
	}

	default:
	{
		cout << endl << ("Error:: wrong parameter value for adaptive mutation model, must be uniform/normal/gamma/negExp \n");
	}
	}

	DistributionType dominanceDistribution = pProtoTrait->getDominanceDistribution();
	map<parameter_t, float> dominanceParameters = pProtoTrait->getDominanceParameters();

	switch (dominanceDistribution) {
	case UNIFORM:
	{
		if (!dominanceParameters.count(MAX))
			cout << endl << ("Error:: adaptive dominance uniform distribution parameter must contain max value (e.g. max= ) \n");

		if (!dominanceParameters.count(MIN))
			cout << endl << ("Error:: adaptive dominance uniform distribution parameter must contain min value (e.g. min= ) \n");

		break;
	}
	case NORMAL:
	{

		if (!dominanceParameters.count(MEAN))
			cout << endl << ("Error:: adaptive dominance distribution set to normal so parameters must contain mean value (e.g. mean= ) \n");

		if (!dominanceParameters.count(SDEV))
			cout << endl << ("Error:: adaptive dominance distribution set to normal so parameters must contain sdev value (e.g. sdev= ) \n");

		break;
	}
	case GAMMA:
	{
		if (!dominanceParameters.count(SHAPE))
			cout << endl << ("Error:: adaptive dominance distribution set to gamma so parameters must contain shape value (e.g. shape= ) \n");

		if (!dominanceParameters.count(SCALE))
			cout << endl << ("Error:: adaptive dominance distribution set to gamma so parameters must contain scale value (e.g. scale= ) \n");

		break;
	}
	case NEGEXP:
	{
		if (!dominanceParameters.count(MEAN))
			cout << endl << ("Error:: adaptive dominance distribution set to negative exponential (negative decay) so parameters must contain mean value (e.g. mean= ) \n");

		break;
	}
	case SCALED:
	{
		break;
	}

	default:
	{
		cout << endl << ("Error:: wrong parameter value for adaptive dominance model, must be uniform/normal/gamma/negExp/scaled \n");
		break; //should return false
	}
	}

	DistributionType initialDistribution = pProtoTrait->getInitialDistribution();
	map<parameter_t, float> initialParameters = pProtoTrait->getInitialParameters();

	//switch (expressionType) {
	//case MULTIPLICATIVE:
	//{
	//    _express_func_ptr = &AdaptiveTrait::expressMulti;
	//    break;
	//}
	//default:
	//{
	//    cout << endl << ("wrong parameter value for parameter \"expression of adaptive mutations\", must be multiplicative  \n");
	//    break; //should return false
	//}

	//}

}


// ----------------------------------------------------------------------------------------
// for creating new individuals
// ----------------------------------------------------------------------------------------

GeneticLoad::GeneticLoad(const GeneticLoad& T) : pProtoTrait(T.pProtoTrait), _inherit_func_ptr(T._inherit_func_ptr)
{}

// ----------------------------------------------------------------------------------------
// mutate uniform
// ----------------------------------------------------------------------------------------
void GeneticLoad::mutate()
{
	const int positionsSize = pProtoTrait->getPositionsSize();
	const auto& positions = pProtoTrait->getPositions();
	const short ploidy = pProtoTrait->getPloidy();
	const float mutationRate = pProtoTrait->getMutationRate();

	auto rng = pRandom->getRNG();

	for (int p = 0; p < ploidy; p++) {

		// Determine nb of mutations
		unsigned int NbMut = pRandom->Poisson(positionsSize * mutationRate);

		if (NbMut > 0) {
			vector<int> mutationPositions;
			// Draw which positions mutate
			sample(positions.begin(), positions.end(), std::back_inserter(mutationPositions),
				NbMut, rng);

			for (int m : mutationPositions) {

				float newSelectionCoef = drawSelectionCoef();
				float newDominanceCoef = drawDominance(newSelectionCoef);

				auto it = mutations.find(m); //find if position in map already has mutations there

				if (it == mutations.end()) {   // not found so create new entry in map with wildtype as char default
					vector<shared_ptr<Mutation>> newAllelePair(2);
					newAllelePair[p] = make_shared<Mutation>(newSelectionCoef, newDominanceCoef); //put new mutation value in 
					mutations.insert(make_pair(m, newAllelePair));
				}
				else { //position found, already mutations there
					it->second[p] = make_shared<Mutation>(newSelectionCoef, newDominanceCoef);
				}
			}
		}
	}
}
// ----------------------------------------------------------------------------------------
// get dominance value for new mutation
// ----------------------------------------------------------------------------------------


float GeneticLoad::drawDominance(float selCoef) {

	DistributionType dominanceDistribution = pProtoTrait->getDominanceDistribution();
	map<parameter_t, float> dominanceParameters = pProtoTrait->getDominanceParameters();

	float h = 1.0; //default dominance is  1

	switch (dominanceDistribution) {
	case UNIFORM:
	{
		float maxD = dominanceParameters.find(MAX)->second;
		float minD = dominanceParameters.find(MIN)->second;
		h = pRandom->FRandom(minD, maxD);
		break;
	}
	case NORMAL:
	{
		const float mean = dominanceParameters.find(MEAN)->second;
		const float sd = dominanceParameters.find(SDEV)->second;
		h = pRandom->Normal(mean, sd);
		break;
	}
	case GAMMA:
	{
		const float shape = dominanceParameters.find(SHAPE)->second;
		const float scale = dominanceParameters.find(SCALE)->second;
		h = pRandom->Gamma(shape, scale);
		break;
	}
	case NEGEXP:
	{
		const float mean = dominanceParameters.find(MEAN)->second;
		h = pRandom->NegExp(mean);
		break;
	}
	case SCALED:
	{
		const float min = 0;
		const float max = exp((-log(2 * 0.36) / 0.05) * selCoef);
		h = pRandom->FRandom(min, max);
		break;
	}

	default:
	{
		cout << endl << ("Error:: wrong parameter value for adaptive dominance model, must be uniform/normal/gamma/negExp/scaled \n");
		break; //should return false
	}
	}

	return h;
}

// ----------------------------------------------------------------------------------------
// get selection coefficient for new mutation
// ----------------------------------------------------------------------------------------


float GeneticLoad::drawSelectionCoef() {

	DistributionType mutationDistribution = pProtoTrait->getMutationDistribution();
	map<parameter_t, float> mutationParameters = pProtoTrait->getMutationParameters();

	float s = 0.0; //default selection coefficient is 0

	switch (mutationDistribution) {
	case UNIFORM:
	{
		float maxD = mutationParameters.find(MAX)->second;
		float minD = mutationParameters.find(MIN)->second;
		s = pRandom->FRandom(minD, maxD);

		break;
	}
	case NORMAL:
	{
		const float mean = mutationParameters.find(MEAN)->second;
		const float sd = mutationParameters.find(SDEV)->second;
		s = pRandom->Normal(mean, sd);

		break;
	}
	case GAMMA:
	{
		const float shape = mutationParameters.find(SHAPE)->second;
		const float scale = mutationParameters.find(SCALE)->second;
		s = pRandom->Gamma(shape, scale);
		break;
	}
	case NEGEXP:
	{
		const float mean = mutationParameters.find(MEAN)->second;
		s = pRandom->NegExp(mean);
		break;
	}
	default:
	{
		cout << endl << ("Error:: wrong parameter value for adaptive mutation model, must be uniform/normal/gamma/negExp/scaled \n");
		break; //should return false
	}
	}
	return s;
}


// ----------------------------------------------------------------------------------------
// inheritance options
// ----------------------------------------------------------------------------------------


void GeneticLoad::inherit(TTrait* parent, set<unsigned int> const& recomPositions, sex_t whichChromosome, int startingChromosome)
{

	auto parentCast = dynamic_cast<GeneticLoad*> (parent); //horrible

	const auto& parent_seq = parentCast->get_mutations();
	if (parent_seq.size() > 0) //else nothing to inherit
		(this->*_inherit_func_ptr) (whichChromosome, parent_seq, recomPositions, startingChromosome);


}

void GeneticLoad::inheritDiploid(sex_t whichChromosome, map<int, vector<shared_ptr<Mutation>>> const& parentMutations, set<unsigned int> const& recomPositions, int parentChromosome) {

	auto it = recomPositions.lower_bound(parentMutations.begin()->first);

	unsigned int nextBreakpoint = *it;

	auto distance = std::distance(recomPositions.begin(), it);
	if (distance % 2 != 0)
		parentChromosome = !parentChromosome; // switch to the other one
		// use 1-parentChromosome, or switch to a sex_t ?

	for (auto const& [locus, allelePair] : parentMutations) {

		while (locus > nextBreakpoint) {
			std::advance(it, 1);
			nextBreakpoint = *it;
			parentChromosome = !parentChromosome; // switch to the other one
		}

		if (locus <= nextBreakpoint) {
			auto& allele = allelePair[parentChromosome];
			auto it = mutations.find(locus);
			if (it == mutations.end()) {
				// locus does not exist yet, initiate it
				vector<shared_ptr<Mutation>> newAllelePair(2);
				newAllelePair[whichChromosome] = allele;
				mutations.insert(make_pair(locus, newAllelePair));
			}
			else {
				// locus already exists
				// presumably bc already inherited from other parent
				// set corresponding allele
				it->second[whichChromosome] = allele;
			}
		}
	}
}

void GeneticLoad::inheritHaploid(sex_t chromosome, map<int, vector<shared_ptr<Mutation>>> const& parentMutations, set<unsigned int> const& recomPositions, int parentChromosome)
{
	mutations = parentMutations;
}

// ----------------------------------------------------------------------------------------
// expression options
// ----------------------------------------------------------------------------------------

float GeneticLoad::express() {

	float phenotype = 1.0;

	for (auto const& [locus, pAllelePair] : mutations)
	{
		auto pAlleleLeft  = (!pAllelePair[0]) ? wildType : pAllelePair[0];
		auto pAlleleRight = (!pAllelePair[1]) ? wildType : pAllelePair[1];

		if (pAlleleLeft.get()->getId() != pAlleleRight.get()->getId()) // heterozygote
		{
			phenotype *= 1 + pAlleleLeft->getSelectionCoef()  * pAlleleLeft->getDominanceCoef();
			phenotype *= 1 + pAlleleRight->getSelectionCoef() * pAlleleRight->getDominanceCoef();
		}
		else { // homozygote
			phenotype *= 1 + pAlleleLeft->getSelectionCoef();
			phenotype *= 1 + pAlleleRight->getSelectionCoef();
		}
	}
	return phenotype;
}

// ----------------------------------------------------------------------------------------
// check if particular loci is heterozygote
// ----------------------------------------------------------------------------------------


bool GeneticLoad::isHeterozygoteAtLoci(int loci) const {

	auto it = mutations.find(loci);

	if (it == mutations.end()) //not found so must be wildtype homozygous
		return false;
	else {
		auto a = (!it->second[0]) ? wildType : it->second[0];
		auto b = (!it->second[1]) ? wildType : it->second[1];
		return a != b;

	}
}

// ----------------------------------------------------------------------------------------
// count heterozygote loci in genome 
// ----------------------------------------------------------------------------------------


int GeneticLoad::countHeterozygoteLoci() const {

	int count = 0;

	for (auto const& [key, val] : mutations) {

		auto a = (!val[0]) ? wildType : val[0];
		auto b = (!val[1]) ? wildType : val[1];
		count += a != b;
	}

	return count;
}

// ----------------------------------------------------------------------------------------
// get allele value at loci 
// ----------------------------------------------------------------------------------------


float GeneticLoad::getSelectionCoefAtLoci(short chromosome, int position) const {

	auto it = mutations.find(position);

	if (it == mutations.end()) {
		return wildType->getSelectionCoef(); //must still be wildtype at loci
	}
	else
		return (!it->second[chromosome]) ? wildType->getSelectionCoef() : it->second[chromosome]->getSelectionCoef();
}