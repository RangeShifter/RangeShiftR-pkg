#ifndef GENETICFITNESSH
#define GENETICFITNESSH

#include <vector>
#include <string>
#include <memory>
#include <map>
#include <algorithm>

#include "QuantitativeTrait.h"

using namespace std;

// Genetic Load trait
// 
// Alleles contribute to their value to reducing
// offspring viability. Alleles start with value
// zero but increase through simulation via mutations.
// There can be up to five genetic load traits.
class GeneticFitnessTrait : public QuantitativeTrait {

public:

    // Initialisation constructor, set initial values and immutable features 
    GeneticFitnessTrait(SpeciesTrait* P);

    // Inheritance constructor, copies pointers to immutable features when cloning from parent
    GeneticFitnessTrait(const GeneticFitnessTrait& T);

   // Make a shallow copy to pass to offspring trait
    // Return new pointer to new trait created by inheritance c'tor 
    // This avoids copying shared attributes: distributions and parameters
    virtual unique_ptr<QuantitativeTrait> clone() const override { return std::make_unique<GeneticFitnessTrait>(*this); }

    virtual ~GeneticFitnessTrait() { }

    // Getters
    virtual int getNLoci()  const override { return pSpeciesTrait->getPositionsSize(); }
    float getMutationRate() const override { return pSpeciesTrait->getMutationRate(); }
    bool isInherited() const override { return pSpeciesTrait->isInherited(); }

    map<int, vector<shared_ptr<Allele>>>& getGenes() { return genes; } // returning reference, reciever must be const

    virtual void mutate() override;
    virtual float express();
    virtual void inheritGenes(const bool& fromMother, QuantitativeTrait* parent, set<unsigned int> const& recomPositions, int startingChromosome) override;

    virtual float getAlleleValueAtLocus(short chromosome, int position) const override;
    virtual int countHeterozygoteLoci() const;
    virtual bool isHeterozygoteAtLocus(int locus) const override;

private:

    // Default allele has value 0 and dominance 0
    inline static shared_ptr<Allele> wildType = make_shared<Allele>(0.0, 0.0);

    // <Locus position, <Allele A, Allele B> >
    map<int, vector<shared_ptr<Allele>>> genes; // position <strand A , strand B>>

    // Initialisation
    void initialise();

    // Mutation
    float scaledDomMeanSelCoeff; // s_d, only for scaled dominance distribution
    float drawDominance(float);
    float drawSelectionCoef();

    // Immutable features, set at initialisation
    // and passed down to every subsequent trait copy
    //// Species-level trait attributes, invariant across individuals
    SpeciesTrait* pSpeciesTrait;
    //// Species-level trait functions
    void (GeneticFitnessTrait::* _inherit_func_ptr) (const bool& fromMother, map<int, vector<shared_ptr<Allele>>> const& parent, set<unsigned int> const& recomPositions, int parentChromosome);

    // Possible values for immutable functions
     //// Inheritance
    void inheritDiploid(const bool& fromMother, map<int, vector<shared_ptr<Allele>>> const& parent, set<unsigned int> const& recomPositions, int parentChromosome);
    void inheritHaploid(const bool& fromMother, map<int, vector<shared_ptr<Allele>>> const& parent, set<unsigned int> const& recomPositions, int parentChromosome);

};
#endif // GENETICFITNESSH
