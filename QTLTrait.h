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
	
	// Species-level constructor, initialise
	QTLTrait(SpeciesTrait* P);

	// Individual-level constructor, called when new individual is created
	QTLTrait(const QTLTrait& T);
	
	virtual ~QTLTrait() { }

	virtual unique_ptr<TTrait> clone() const override { return std::make_unique<QTLTrait>(*this); }
	
	// Getters
	int getNLoci()  const override { return pSpeciesTrait->getPositionsSize(); }
	float getMutationRate() const override { return pSpeciesTrait->getMutationRate(); }
	bool isInherited() const override { return pSpeciesTrait->isInherited(); }
	
	map<int, vector<shared_ptr<Allele>>>& getGenes() { return genes; } // returning reference, receiver must be const
	
	void mutate() override { (this->*_mutate_func_ptr) (); }
	float express() override { return (this->*_express_func_ptr) (); }
	void inherit(const bool& fromMother, TTrait* parent, set<unsigned int> const& recomPositions, int startingChromosome) override;

	float getAlleleValueAtLocus(short chromosome, int i) const override;
	int countHeterozygoteLoci() const;
	bool isHeterozygoteAtLocus(int locus) const override;

private:

	const double QTLDominanceFactor = 1.0; // no dominance for QTL traits (yet?)

	// Species-level trait attributes, invariant across individuals
	SpeciesTrait* pSpeciesTrait; 

	map<int, vector<shared_ptr<Allele>>> genes; // position <strand A , strand B>>

	void (QTLTrait::* _mutate_func_ptr) (void);
	void (QTLTrait::* _inherit_func_ptr) (const bool& fromMother, map<int, vector<shared_ptr<Allele>>> const& parent, set<unsigned int> const& recomPositions, int parentChromosome);
	float (QTLTrait::* _express_func_ptr) (void);

	void initialiseUniform(float min, float max);
	void initialiseNormal(float mean, float sd);
	void inheritDiploid(const bool& fromMother, map<int, vector<shared_ptr<Allele>>> const& parent, set<unsigned int> const& recomPositions, int parentChromosome);
	void inheritHaploid(const bool& fromMother, map<int, vector<shared_ptr<Allele>>> const& parent, set<unsigned int> const& recomPositions, int parentChromosome);
	void reInitialiseGenes(const bool& fromMother, map<int, vector<shared_ptr<Allele>>> const& parentMutations, set<unsigned int> const& recomPositions, int parentChromosome);
	void mutateUniform();
	void mutateNormal();
	float expressAverage();
	float expressAdditive();
};

#endif