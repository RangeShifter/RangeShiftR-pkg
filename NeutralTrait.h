#ifndef NeutralTRAITH
#define NeutralTRAITH

#include <vector>
#include <string>
#include <memory>
#include <algorithm>

#include "QuantitativeTrait.h"

using namespace std;

// Neutral traits
// 
// Not expressed and are only used to compute neutral statistics 
// e.g. the Fst.
// To save on mem usage, allele values are represented by character types,
// taking a value between 0 and a user-specified max >= 255
class NeutralTrait : public QuantitativeTrait {

public:

	// Initialisation constructor, set initial values and immutable features 
	NeutralTrait(SpeciesTrait* P);

	// Inheritance constructor, copies pointers to immutable features when cloning from parent
	NeutralTrait(const NeutralTrait& T);

	// Make a shallow copy to pass to offspring trait
	// Return new pointer to new trait created by inheritance c'tor 
	// This avoids copying shared attributes: distributions and parameters
	virtual unique_ptr<QuantitativeTrait> clone() const override { return std::make_unique<NeutralTrait>(*this); }

	virtual ~NeutralTrait() { }

	// Getters
	virtual int getNLoci()  const override { return pSpeciesTrait->getPositionsSize(); }
	float getMutationRate() const override { return pSpeciesTrait->getMutationRate(); }
	bool isInherited() const override { return pSpeciesTrait->isInherited(); }
	map<int, vector<unsigned char>>& getGenes() { return genes; } //returning reference, reciever must be const

	virtual void mutate() override { (this->*_mutate_func_ptr) (); }
	virtual void inheritGenes(const bool& fromMother, QuantitativeTrait* parent, set<unsigned int> const& recomPositions, int startingChromosome) override;
	virtual float express() {
		throw runtime_error("Neutral trait shouldn't be expressed.");
		return -9999;
	}

	virtual float getAlleleValueAtLocus(short chromosome, int position) const override;
	virtual float getDomCoefAtLocus(short chromosome, int position) const override {
		return 0.0;
	}

	virtual int countHeterozygoteLoci() const;
	virtual bool isHeterozygoteAtLocus(int locus) const override;
#if RSDEBUG // for testing only
	int getAlleleIDAtLocus(short whichChromosome, int position) const;
#endif

private:

	inline static int wildType; // default allele value, value set at construction
	const int NeutralValUpperBound = UCHAR_MAX; // alleles are char, can take value 0-255

	// <Locus position, <Allele A, Allele B>>
	map<int, vector<unsigned char>> genes;

	// Initialisation
	void initialiseUniform(int max); //other option is that mutations map is empty until a mutation happens, default when empty is to return a 0 value for wildtype

	// Immutable features, set at initialisation
	// and passed down to every subsequent trait copy
	//// Species-level trait attributes, invariant across individuals
	SpeciesTrait* pSpeciesTrait;
	//// Species-level trait functions
	void (NeutralTrait::* _mutate_func_ptr) (void);
	void (NeutralTrait::* _inherit_func_ptr) (const bool& fromMother, map<int, vector<unsigned char>> const& parent, set<unsigned int> const& recomPositions, int parentChromosome);

	// Possible values for immutable functions
	//// Inheritance
	void inheritDiploid(const bool& fromMother, map<int, vector<unsigned char>> const&, set<unsigned int> const& recomPositions, int parentChromosome);
	void inheritHaploid(const bool& fromMother, map<int, vector<unsigned char>> const& parentMutations, set<unsigned int> const& recomPositions, int parentChromosome);
	//// Mutation
	void mutate_KAM();
	void mutate_SSM(); // single-step mutations

};
#endif