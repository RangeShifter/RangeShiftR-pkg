
#include "SpeciesTrait.h"

// Species trait constructor
SpeciesTrait::SpeciesTrait(
	const TraitType& traitType, const sex_t& sx, 
	const set<int>& pos, const ExpressionType& expr,
	const DistributionType& initDist, const map<GenParamType, float> initParams,
	const DistributionType& dominanceDist, const map<GenParamType, float> dominanceParams,
	bool isInherited, const float& mutRate,
	const DistributionType& mutationDist, const map<GenParamType, float> mutationParams,
	Species* pSpecies) :
	sex{sx},
	genePositions{pos},
	expressionType{expr},
	initialDistribution{initDist},
	initialParameters{initParams},
	dominanceDistribution{dominanceDist},
	dominanceParameters{dominanceParams},
	inherited{isInherited},
	mutationDistribution{mutationDist},
	mutationParameters{mutationParams},
	mutationRate{mutRate}
{
	// Initialise ploidy only once per species
	if (ploidy == NULL) this->ploidy = pSpecies->isDiploid() ? 2 : 1;
}
