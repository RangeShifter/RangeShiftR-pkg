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
		else if (traitType == ADAPTIVE1 || traitType == ADAPTIVE2 || traitType == ADAPTIVE3 || traitType == ADAPTIVE4 || traitType == ADAPTIVE5) {
			return make_unique<GeneticLoad>(protoTrait);
		}
		else {
			return make_unique<QTLTrait>(protoTrait);
		}
	}
};
#endif