/*----------------------------------------------------------------------------
 *
 *	Copyright (C) 2026 Greta Bocedi, Stephen C.F. Palmer, Justin M.J. Travis, Anne-Kathleen Malchow, Roslyn Henry, Théo Pannetier, Jette Wolff, Damaris Zurell
 *
 *	This file is part of RangeShifter.
 *
 *	RangeShifter is free software: you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *	the Free Software Foundation, either version 3 of the License, or
 *	(at your option) any later version.
 *
 *	RangeShifter is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *	GNU General Public License for more details.
 *
 *	You should have received a copy of the GNU General Public License
 *	along with RangeShifter. If not, see <https://www.gnu.org/licenses/>.
 *
 * File Created by Roslyn Henry March 2023. Code adapted from NEMO (https://nemo2.sourceforge.io/)
 --------------------------------------------------------------------------*/

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
    virtual float getDomCoefAtLocus(short chromosome, int position) const override;

private:

    // Default allele has value 0 and dominance 0
    inline static shared_ptr<Allele> wildType = make_shared<Allele>(0.0, 0.0);

    // <Locus position, <Allele A, Allele B> >
    map<int, vector<shared_ptr<Allele>>> genes; // position <strand A , strand B>>

    // Initialisation
    void initialise();

    void setScaledCoeff(const DistributionType& selCoeffDist, const map<GenParamType, float>& selCoeffParams);

    // Mutation
    float scaledDomMeanSelCoeff = 0; // s_d, only for scaled dominance distribution
    float drawDominance(
        float selCoef, 
        const DistributionType& domDist, 
        const map<GenParamType, float>& domParams
    );
    float drawSelectionCoef(
        const DistributionType& mutationDistribution, 
        const map<GenParamType, float>& mutationParameters
    );

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
