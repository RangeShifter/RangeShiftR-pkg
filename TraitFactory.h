#ifndef TRAITFACTORYH
#define TRAITFACTORYH

#include <memory>

#include "SpeciesTrait.h"
#include "SNPTrait.h"
#include "QTLTrait.h"
#include "GeneticLoad.h"

class TraitFactory
{
public:
	TraitFactory() {};

	unique_ptr<TTrait> Create(const TraitType traitType, SpeciesTrait* protoTrait)
	{
		if (traitType == SNP) {
			return make_unique<SNPTrait>(protoTrait);
		}
		else if (traitType == GENETIC_LOAD1 
			|| traitType == GENETIC_LOAD2 
			|| traitType == GENETIC_LOAD3 
			|| traitType == GENETIC_LOAD4 
			|| traitType == GENETIC_LOAD5) {
			return make_unique<GeneticLoad>(protoTrait);
		}
		else {
			return make_unique<QTLTrait>(protoTrait);
		}
	}
};
#endif