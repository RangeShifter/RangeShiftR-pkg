
#include "SpeciesTrait.h"

//could be handled in header file but here for now for flexibility
SpeciesTrait::SpeciesTrait(
	const TraitType& traitType, const sex_t& sx, 
	const set<int>& pos, const ExpressionType& expr,
	const DistributionType& initDist, const map<parameter_t, float> initParams,
	const DistributionType& dominanceDist, const map<parameter_t, float> dominanceParams,
	bool isInherited, const float& mutRate,
	const DistributionType& mutationDist, const map<parameter_t, float> mutationParams,
	Species* pSpecies) :
	sex{sx},
	positions{pos},
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
	if (ploidy == 0) this->ploidy = pSpecies->isDiploid() ? 2 : 1;
}
