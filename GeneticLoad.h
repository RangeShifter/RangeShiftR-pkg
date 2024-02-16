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

    inline static shared_ptr<Allele> wildType;

   // void (AdaptiveTrait::* _mutate_func_ptr) (void);
    void (GeneticLoad::* _inherit_func_ptr) (sex_t chromosome, map<int, vector<shared_ptr<Allele>>> const& parent, set<unsigned int> const& recomPositions, int parentChromosome);
    //float (AdaptiveTrait::* _express_func_ptr) (void);

    void inheritDiploid(sex_t chromosome, map<int, vector<shared_ptr<Allele>>> const& parent, set<unsigned int> const& recomPositions, int parentChromosome);
    void inheritHaploid(sex_t chromosome, map<int, vector<shared_ptr<Allele>>> const& parent, set<unsigned int> const& recomPositions, int parentChromosome);

    float drawDominance(float);
    float drawSelectionCoef();

public:

    //this one for species held trait table, e.g. prototype table, sets static members

    GeneticLoad(SpeciesTrait* P);

    //this one for individuals, static members are not reset
    GeneticLoad(const GeneticLoad& T);

    virtual unique_ptr<TTrait> clone() const override { return std::make_unique<GeneticLoad>(*this); }

    virtual void  inherit(TTrait* parent, set<unsigned int> const& recomPositions, sex_t chromosome, int startingChromosome) override;

    virtual void  mutate() override;

    virtual int getNLoci()  const override { return pSpeciesTrait->getPositionsSize(); }

    float getMutationRate() const override { return pSpeciesTrait->getMutationRate(); }

    bool isInherited() const override { return pSpeciesTrait->isInherited(); }

    map<int, vector<shared_ptr<Allele>>>& get_mutations() { return genes; } //returning reference, reciever must be const

    virtual float getSelectionCoefAtLoci(short chromosome, int position) const override;

    virtual int countHeterozygoteLoci() const;

    virtual bool isHeterozygoteAtLocus(int locus) const override;

    virtual float express();

    virtual ~GeneticLoad() { }
};
#endif // GENETICLOADH
