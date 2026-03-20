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

#ifndef DISPTRAITH
#define DISPTRAITH

#include <vector>
#include <string>
#include <memory>
#include <algorithm>

#include "QuantitativeTrait.h"

using namespace std;

// Dispersal trait
// 
// That is, all evolvable that control emigration, transfer and settlement
class DispersalTrait : public QuantitativeTrait {

public:
	
	// Initialisation constructor, set initial values and immutable features 
	DispersalTrait(SpeciesTrait* P);

	// Inheritance constructor, copies pointers to immutable features when cloning from parent
	DispersalTrait(const DispersalTrait& T);
	
	// Make a shallow copy to pass to offspring trait
	// Return new pointer to new trait created by inheritance c'tor 
	// This avoids copying shared attributes: distributions and parameters
	virtual unique_ptr<QuantitativeTrait> clone() const override { return std::make_unique<DispersalTrait>(*this); }

	virtual ~DispersalTrait() { }
	
	// Getters
	int getNLoci()  const override { return pSpeciesTrait->getPositionsSize(); }
	float getMutationRate() const override { return pSpeciesTrait->getMutationRate(); }
	bool isInherited() const override { return pSpeciesTrait->isInherited(); }
	
	map<int, vector<shared_ptr<Allele>>>& getGenes() { return genes; } // returning reference, receiver must be const
	
	void mutate() override { (this->*_mutate_func_ptr) (); }
	float express() override { return (this->*_express_func_ptr) (); }
	void inheritGenes(const bool& fromMother, QuantitativeTrait* parent, set<unsigned int> const& recomPositions, int startingChromosome) override;

	float getAlleleValueAtLocus(short chromosome, int i) const override;
	float getDomCoefAtLocus(short chromosome, int position) const override;

#ifdef UNIT_TESTS // for testing only
	void overwriteGenes(map<int, vector<shared_ptr<Allele>>> genSeq) {
		genes = genSeq;
	}
	void triggerInherit(
		// inheritGenes requires passing a QuantitativeTrait, unfeasible in tests
		const bool& fromMother, 
		map<int, vector<shared_ptr<Allele>>> const& parentGenes, 
		set<unsigned int> const& recomPositions, 
		int startChr) {
		(this->*_inherit_func_ptr)(fromMother, parentGenes, recomPositions, startChr);
	}
#endif

private:

	const double dispDominanceFactor = 1.0; // no dominance for Dispersal traits (yet?)

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
	void (DispersalTrait::* _mutate_func_ptr) (void);
	void (DispersalTrait::* _inherit_func_ptr) (const bool& fromMother, map<int, vector<shared_ptr<Allele>>> const& parent, set<unsigned int> const& recomPositions, int parentChromosome);
	float (DispersalTrait::* _express_func_ptr) (void);

	// Possible values for immutable functions
	//// Inheritance
	void inheritDiploid(const bool& fromMother, map<int, vector<shared_ptr<Allele>>> const& parent, set<unsigned int> const& recomPositions, int parentChromosome);
	void inheritHaploid(const bool& fromMother, map<int, vector<shared_ptr<Allele>>> const& parent, set<unsigned int> const& recomPositions, int parentChromosome);
	void reInitialiseGenes(const bool& fromMother, map<int, vector<shared_ptr<Allele>>> const& parentMutations, set<unsigned int> const& recomPositions, int parentChromosome);
	//// Mutation
	void mutateUniform();
	void mutateNormal();
	void trimPhenotype(float& phenotype);
	//// Gene expression
	float expressAverage();
	float expressAdditive();
};

#ifdef UNIT_TESTS
// Test utilities

map<int, vector<shared_ptr<Allele>>> createTestGenotype(
	const int genomeSz, const bool isDiploid,
	const float valAlleleA,
	const float valAlleleB = -99.9, // allow no value for haploids
	const float domCoeffA = 1.0, // default for dispersal traits
	const float domCoeffB = 1.0
);
#endif // UNIT_TESTS

#endif // DISPTRAITH
