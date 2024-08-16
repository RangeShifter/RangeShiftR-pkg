#ifndef TRAITTH
#define TRAITTH

#include <vector>
#include <set>
#include <memory>

#include "Parameters.h"
#include "Allele.h"
#include "SpeciesTrait.h"

using namespace std;

// Base interface for all genetic traits
class QuantitativeTrait {
public:
    virtual void mutate() = 0;
    virtual unique_ptr<QuantitativeTrait> clone() const = 0; //copies parameters (if not static) not gene seqeunces
    virtual void inheritGenes(const bool&, QuantitativeTrait*, set<unsigned int> const& , int ) = 0;
    virtual int getNLoci() const = 0;
    virtual float getMutationRate() const = 0;
    virtual bool isInherited() const = 0;
    virtual float getAlleleValueAtLocus(short chromosome, int i) const = 0;
    virtual int countHeterozygoteLoci() const = 0;
    virtual bool isHeterozygoteAtLocus(int loci) const = 0;
    virtual float express() = 0;
    virtual ~QuantitativeTrait() { }
#if RSDEBUG // for testing only
	virtual int getAlleleIDAtLocus(short whichChromosome, int position) const = 0;
#endif
};

extern RSrandom* pRandom;

#endif