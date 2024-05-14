#ifndef QTLTRAITH
#define QTLTRAITH

#include <vector>
#include <string>
#include <memory>
#include <algorithm>

#include "TTrait.h"

using namespace std;

// Quantitative-trait-loci (QTL) trait
// 
// That is, all dispersal traits in RangeShifter,
// which individual phenotype is determined by small contributions 
// from (potentially) many quantitative loci
class QTLTrait : public TTrait {

public:
	
	// Initialisation constructor, set initial values and immutable features 
	QTLTrait(SpeciesTrait* P);

	// Inheritance constructor, copies pointers to immutable features when cloning from parent
	QTLTrait(const QTLTrait& T);
	
	// Make a shallow copy to pass to offspring trait
	// Return new pointer to new trait created by inheritance c'tor 
	// This avoids copying shared attributes: distributions and parameters
	virtual unique_ptr<TTrait> clone() const override { return std::make_unique<QTLTrait>(*this); }

	virtual ~QTLTrait() { }
	
	// Getters
	int getNLoci()  const override { return pSpeciesTrait->getPositionsSize(); }
	float getMutationRate() const override { return pSpeciesTrait->getMutationRate(); }
	bool isInherited() const override { return pSpeciesTrait->isInherited(); }
	
	map<int, vector<shared_ptr<Allele>>>& getGenes() { return genes; } // returning reference, receiver must be const
	
	void mutate() override { (this->*_mutate_func_ptr) (); }
	float express() override { return (this->*_express_func_ptr) (); }
	void inheritGenes(const bool& fromMother, TTrait* parent, set<unsigned int> const& recomPositions, int startingChromosome) override;

	float getAlleleValueAtLocus(short chromosome, int i) const override;
	int countHeterozygoteLoci() const;
	bool isHeterozygoteAtLocus(int locus) const override;

private:

	const double QTLDominanceFactor = 1.0; // no dominance for QTL traits (yet?)

	// <Locus position, <Allele A, Allele B> >
	map<int, vector<shared_ptr<Allele>>> genes;

	// Initialisation
	void initialiseUniform(float min, float max);
	void initialiseNormal(float mean, float sd);

	// Immutable features, set at initialisation
	// and passed down to every subsequent trait copy
	//// Species-level trait attributes, invariant across individuals
	SpeciesTrait* pSpeciesTrait;
	//// Species-level trait functions
	void (QTLTrait::* _mutate_func_ptr) (void);
	void (QTLTrait::* _inherit_func_ptr) (const bool& fromMother, map<int, vector<shared_ptr<Allele>>> const& parent, set<unsigned int> const& recomPositions, int parentChromosome);
	float (QTLTrait::* _express_func_ptr) (void);

	// Possible values for immutable functions
	//// Inheritance
	void inheritDiploid(const bool& fromMother, map<int, vector<shared_ptr<Allele>>> const& parent, set<unsigned int> const& recomPositions, int parentChromosome);
	void inheritHaploid(const bool& fromMother, map<int, vector<shared_ptr<Allele>>> const& parent, set<unsigned int> const& recomPositions, int parentChromosome);
	void reInitialiseGenes(const bool& fromMother, map<int, vector<shared_ptr<Allele>>> const& parentMutations, set<unsigned int> const& recomPositions, int parentChromosome);
	//// Mutation
	void mutateUniform();
	void mutateNormal();
	void trimQTLPhenotype(float& phenotype);
	//// Gene expression
	float expressAverage();
	float expressAdditive();
};

#endif