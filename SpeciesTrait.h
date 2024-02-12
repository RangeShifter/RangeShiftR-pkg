#ifndef SPECIESTRAITH
#define SPECIESTRAITH

#include "Parameters.h"
#include "Species.h"

#include <map>
#include <vector>
#include <set>
#include <memory>

class Species;

class SpeciesTrait {

private:
    inline static int ploidy = 0;
    float mutationRate;
    sex_t sex;
    set<int> positions;
    ExpressionType expressionType;
    DistributionType initialDistribution;
    map<parameter_t, float> initialParameters;
    DistributionType dominanceDistribution;
    map<parameter_t, float> dominanceParameters;
    bool inherited;
    DistributionType mutationDistribution;
    map<parameter_t, float> mutationParameters;

public:

    SpeciesTrait(vector<string> parameters, Species* pSpecies);

    // Getters
    sex_t getSex() const { return sex; }
    float getMutationRate() const { return mutationRate; }
    short getPloidy() const { return ploidy; }
    set<int>& getPositions() { return positions; } // returning by reference, make sure receiver is const
    int getPositionsSize() const { return positions.size(); }
    bool isInherited() const { return inherited;  }
    DistributionType getMutationDistribution() const { return mutationDistribution; };
    map<parameter_t, float> getMutationParameters() const { return mutationParameters; };
    DistributionType getDominanceDistribution() const { return dominanceDistribution; };
    map<parameter_t, float> getDominanceParameters() const { return dominanceParameters; };
    DistributionType getInitialDistribution() const { return initialDistribution; };
    map<parameter_t, float> getInitialParameters() const { return initialParameters; };
    ExpressionType getExpressionType() const { return expressionType; };

    DistributionType stringToDistributionType(const std::string& str) const;
    ExpressionType stringToExpressionType(const std::string& str) const;
    map<parameter_t, float> stringToParameterMap(string parameters) const;

    set<int> selectRandomLociPositions(int noLoci, Species* pSpecies) const;

    set<int> stringToLoci(string pos, string nLoci, Species* pSpecies) const;
    TraitType stringToTraitType(const std::string& str, sex_t sex) const;
};
#endif