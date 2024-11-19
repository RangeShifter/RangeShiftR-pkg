#ifndef TRAITFACTORYH
#define TRAITFACTORYH

#include <memory>

#include "SpeciesTrait.h"
#include "NeutralTrait.h"
#include "DispersalTrait.h"
#include "GeneticFitnessTrait.h"

// Create handled pointers to a new trait
class TraitFactory
{
public:
	TraitFactory() {};

	unique_ptr<QuantitativeTrait> Create(const TraitType traitType, SpeciesTrait* pSpTrait)
	{
		if (traitType == NEUTRAL) {
			return make_unique<NeutralTrait>(pSpTrait);
		}
		else if (traitType == GENETIC_LOAD1 
			|| traitType == GENETIC_LOAD2 
			|| traitType == GENETIC_LOAD3 
			|| traitType == GENETIC_LOAD4 
			|| traitType == GENETIC_LOAD5) {
			return make_unique<GeneticFitnessTrait>(pSpTrait);
		}
		else {
			return make_unique<DispersalTrait>(pSpTrait);
		}
	}
};
#endif