#ifndef GENETICLOADH
#define GENETICLOADH

#include <vector>
#include <string>
#include <memory>
#include <map>
#include <algorithm>

#include "TTrait.h"

using namespace std;

class GeneticLoad : public TTrait {

private:

    SpeciesTrait* pSpeciesTrait;

    map<int, vector<shared_ptr<Allele>>> genes; // position <strand A , strand B>>

    inline static shared_ptr<Allele> wildType = make_shared<Allele>(0.0, 0.0);

    void (GeneticLoad::* _inherit_func_ptr) (const bool& fromMother, map<int, vector<shared_ptr<Allele>>> const& parent, set<unsigned int> const& recomPositions, int parentChromosome);

    void initialise();

    void inheritDiploid(const bool& fromMother, map<int, vector<shared_ptr<Allele>>> const& parent, set<unsigned int> const& recomPositions, int parentChromosome);
    void inheritHaploid(const bool& fromMother, map<int, vector<shared_ptr<Allele>>> const& parent, set<unsigned int> const& recomPositions, int parentChromosome);

    float drawDominance(float);
    float drawSelectionCoef();

public:

    //this one for species held trait table, e.g. prototype table, sets static members

    GeneticLoad(SpeciesTrait* P);

    //this one for individuals, static members are not reset
    GeneticLoad(const GeneticLoad& T);

    virtual unique_ptr<TTrait> clone() const override { return std::make_unique<GeneticLoad>(*this); }

    virtual void inherit(const bool& fromMother, TTrait* parent, set<unsigned int> const& recomPositions, int startingChromosome) override;

    virtual void  mutate() override;

    virtual int getNLoci()  const override { return pSpeciesTrait->getPositionsSize(); }

    float getMutationRate() const override { return pSpeciesTrait->getMutationRate(); }

    bool isInherited() const override { return pSpeciesTrait->isInherited(); }

    map<int, vector<shared_ptr<Allele>>>& getGenes() { return genes; } //returning reference, reciever must be const

    virtual float getAlleleValueAtLocus(short chromosome, int position) const override;

    virtual int countHeterozygoteLoci() const;

    virtual bool isHeterozygoteAtLocus(int locus) const override;

    virtual float express();
    virtual ~GeneticLoad() { }
};
#endif // GENETICLOADH
