#ifndef SPECIESTRAITH
#define SPECIESTRAITH

#include "Parameters.h"
#include "Species.h"

#include <map>
#include <vector>
#include <set>
#include <memory>

class Species; // forward declaration to overcome circularity issue

// Species-level traits
// Features of traits that are shared across all individuals in the same species
class SpeciesTrait {
public:
    SpeciesTrait(
        const TraitType& traitType,
        const sex_t& sex,
        const set<int>& pos,
        const ExpressionType& expr,
        const DistributionType& initDist,
        const map<GenParamType, float> initParams,
        const DistributionType& initDomDist,
        const map<GenParamType, float> initDomParams,
        bool isInherited,
        const float& mutationRate,
        const DistributionType& mutationDist,
        const map<GenParamType, float> mutationParams,
        const DistributionType& dominanceDist,
        const map<GenParamType, float> dominanceParams,
        const int ploidy,
        const bool isOutput
    );

    bool isValidTraitVal(const float& val) const;
    TraitType getTraitType() const { return traitType; }
    bool isOutput() const { return traitIsOutput; }

    // Getters
    sex_t getSex() const { return sex; }
    float getMutationRate() const { return mutationRate; }
    short getPloidy() const { return ploidy; }
    set<int>& getGenePositions() { return genePositions; } // returning by reference, make sure receiver is const
    int getPositionsSize() const { return static_cast<int>(genePositions.size()); }
    bool isInherited() const { return inherited; }

    DistributionType getInitialDistribution() const { return initialDistribution; };
    map<GenParamType, float> getInitialParameters() const { return initialParameters; };
    DistributionType getInitDomDistribution() const { return initialDomDistribution; };
    map<GenParamType, float> getInitDomParameters() const { return initialDomParameters; };
    DistributionType getMutationDistribution() const { return mutationDistribution; };
    map<GenParamType, float> getMutationParameters() const { return mutationParameters; };
    DistributionType getDominanceDistribution() const { return dominanceDistribution; };
    map<GenParamType, float> getDominanceParameters() const { return dominanceParameters; };

    ExpressionType getExpressionType() const { return expressionType; };

    int getNbNeutralAlleles() const {
        if (!traitType == NEUTRAL) throw logic_error("getNbNeutralAlleles() should only be called for neutral traits.");
        else {
            int maxAlleleVal = max(
                getMutationParameters().find(MAX)->second + 1,
                getInitialParameters().find(MAX)->second + 1
            ); // possible values range from 0 to MAX
            return maxAlleleVal;
        }
    }

private:

    int ploidy;
    float mutationRate;
    TraitType traitType;
    bool traitIsOutput;
    sex_t sex;

    // Positions in the genome of all genes (loci) pertaining to this trait
    // The genome itself is not modelled explicitly
    set<int> genePositions;

    ExpressionType expressionType;
    DistributionType initialDistribution;
    map<GenParamType, float> initialParameters;
    DistributionType initialDomDistribution;
    map<GenParamType, float> initialDomParameters;
    DistributionType dominanceDistribution;
    map<GenParamType, float> dominanceParameters;
    bool inherited;
    DistributionType mutationDistribution;
    map<GenParamType, float> mutationParameters;
};

#ifndef NDEBUG // Testing only
// Create a default set of gene positions ranging from zero to genome size
set<int> createTestGenePositions(const int genomeSz);
SpeciesTrait* createTestEmigSpTrait(const set<int>& genePositions, const bool& isDiploid);
SpeciesTrait* createTestGenLoadTrait(const set<int>& genePositions, const bool& isDiploid);
SpeciesTrait* createTestNeutralSpTrait(const float& maxAlleleVal, const set<int>& genePositions, const bool& isDiploid);
#endif // NDEBUG

#endif // SPECIESTRAITH