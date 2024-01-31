#ifndef QTLTRAITH
#define QTLTRAITH

#include <vector>
#include <string>

#include "TTrait.h"

using namespace std;


class QTLTrait : public TTrait {

private:

	ProtoTrait* pProtoTrait; // would be better as const so immutable, but means passing positions list is heavy and can't be passed by reference

	map<int, vector<shared_ptr<Mutation>>> mutations; //position <strand A , strand B>>

	void (QTLTrait::* _mutate_func_ptr) (void);
	void (QTLTrait::* _inherit_func_ptr) (sex_t chromosome, map<int, vector<shared_ptr<Mutation>>> const& parent, set<unsigned int> const& recomPositions, int parentChromosome);
	float (QTLTrait::* _express_func_ptr) (void);

	void initialiseUniform(float min, float max);
	void initialiseNormal(float mean, float sd);
	void inheritDiploid(sex_t chromosome, map<int, vector<shared_ptr<Mutation>>> const& parent, set<unsigned int> const& recomPositions, int parentChromosome);
	void inheritHaploid(sex_t chromosome, map<int, vector<shared_ptr<Mutation>>> const& parent, set<unsigned int> const& recomPositions, int parentChromosome);
	void inheritInitialParameters(sex_t chromosome, map<int, vector<shared_ptr<Mutation>>> const& parentMutations, set<unsigned int> const& recomPositions, int parentChromosome);
	void mutateUniform();
	void mutateNormal();
	float expressAverage();
	float expressAdditive();

public:

	//this one for species held trait table, e.g. prototype table, sets static members
	QTLTrait(ProtoTrait* P);
	QTLTrait(const QTLTrait& T);
	virtual unique_ptr<TTrait> clone() const override { return std::make_unique<QTLTrait>(*this); }
	int getNLoci()  const override { return pProtoTrait->getPositionsSize(); }
	float getMutationRate() const override { return pProtoTrait->getMutationRate(); }
	bool isInherited() const override { return pProtoTrait->isInherited(); }
	void mutate() override { (this->*_mutate_func_ptr) (); }
	float express() override { return (this->*_express_func_ptr) (); }
	void inherit(TTrait* parent, set<unsigned int> const& recomPositions, sex_t chromosome, int startingChromosome) override;
	map<int, vector<shared_ptr<Mutation>>>& get_mutations() { return mutations; } // returning reference, receiver must be const

	// virtual float getSelectionCoefAtLoci(short chromosome, int i) const override { return (double)mutations[chromosome][i]->getSelectionCoef(); }

	float getSelectionCoefAtLoci(short chromosome, int i) const override;
	int countHeterozygoteLoci() const;
	bool isHeterozygoteAtLoci(int loci) const override;
	virtual ~QTLTrait() { }
};

#endif