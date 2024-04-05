#ifndef QTLTRAITH
#define QTLTRAITH

#include <vector>
#include <string>
#include <memory>
#include <algorithm>

#include "TTrait.h"

using namespace std;


class QTLTrait : public TTrait {

private:

	const double QTLDominanceFactor = 1.0; // that is, no dominance

	SpeciesTrait* pSpeciesTrait; // would be better as const so immutable, but means passing positions list is heavy and can't be passed by reference

	map<int, vector<shared_ptr<Allele>>> genes; //position <strand A , strand B>>

	void (QTLTrait::* _mutate_func_ptr) (void);
	void (QTLTrait::* _inherit_func_ptr) (sex_t chromosome, map<int, vector<shared_ptr<Allele>>> const& parent, set<unsigned int> const& recomPositions, int parentChromosome);
	float (QTLTrait::* _express_func_ptr) (void);

	void initialiseUniform(float min, float max);
	void initialiseNormal(float mean, float sd);
	void inheritDiploid(sex_t chromosome, map<int, vector<shared_ptr<Allele>>> const& parent, set<unsigned int> const& recomPositions, int parentChromosome);
	void inheritHaploid(sex_t chromosome, map<int, vector<shared_ptr<Allele>>> const& parent, set<unsigned int> const& recomPositions, int parentChromosome);
	void inheritInitialParameters(sex_t chromosome, map<int, vector<shared_ptr<Allele>>> const& parentMutations, set<unsigned int> const& recomPositions, int parentChromosome);
	void mutateUniform();
	void mutateNormal();
	float expressAverage();
	float expressAdditive();

public:

	//this one for species held trait table, e.g. prototype table, sets static members
	QTLTrait(SpeciesTrait* P);
	QTLTrait(const QTLTrait& T);
	virtual unique_ptr<TTrait> clone() const override { return std::make_unique<QTLTrait>(*this); }
	int getNLoci()  const override { return pSpeciesTrait->getPositionsSize(); }
	float getMutationRate() const override { return pSpeciesTrait->getMutationRate(); }
	bool isInherited() const override { return pSpeciesTrait->isInherited(); }
	void mutate() override { (this->*_mutate_func_ptr) (); }
	float express() override { return (this->*_express_func_ptr) (); }
	void inherit(TTrait* parent, set<unsigned int> const& recomPositions, sex_t chromosome, int startingChromosome) override;
	map<int, vector<shared_ptr<Allele>>>& getGenes() { return genes; } // returning reference, receiver must be const

	// virtual float getAlleleValueAtLocus(short chromosome, int i) const override { return (double)mutations[chromosome][i]->getSelectionCoef(); }

	float getAlleleValueAtLocus(short chromosome, int i) const override;
	int countHeterozygoteLoci() const;
	bool isHeterozygoteAtLocus(int locus) const override;
	virtual ~QTLTrait() { }
};

#endif