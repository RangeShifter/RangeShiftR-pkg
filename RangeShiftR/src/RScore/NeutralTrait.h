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

#ifndef NeutralTRAITH
#define NeutralTRAITH

#include <vector>
#include <string>
#include <memory>
#include <algorithm>

#include "QuantitativeTrait.h"

using namespace std;

// Neutral traits
// 
// Not expressed and are only used to compute neutral statistics 
// e.g. the Fst.
// To save on mem usage, allele values are represented by character types,
// taking a value between 0 and a user-specified max >= 255
class NeutralTrait : public QuantitativeTrait {

public:

	// Initialisation constructor, set initial values and immutable features 
	NeutralTrait(SpeciesTrait* P);

	// Inheritance constructor, copies pointers to immutable features when cloning from parent
	NeutralTrait(const NeutralTrait& T);

	// Make a shallow copy to pass to offspring trait
	// Return new pointer to new trait created by inheritance c'tor 
	// This avoids copying shared attributes: distributions and parameters
	virtual unique_ptr<QuantitativeTrait> clone() const override { return std::make_unique<NeutralTrait>(*this); }

	virtual ~NeutralTrait() { }

	// Getters
	virtual int getNLoci()  const override { return pSpeciesTrait->getPositionsSize(); }
	float getMutationRate() const override { return pSpeciesTrait->getMutationRate(); }
	bool isInherited() const override { return pSpeciesTrait->isInherited(); }
	// Genes are stored in a flat array: allele of locus index l on chromosome c is
	// genes[l * ploidy + c], where locus index l is the rank of the gene position
	// in the (sorted) gene positions of the species trait.
	const vector<unsigned char>& getGenes() const { return genes; }
	unsigned char getAlleleAtIndex(short whichChromosome, int locusIndex) const {
		return genes[static_cast<size_t>(locusIndex) * ploidy + whichChromosome];
	}
#ifdef UNIT_TESTS
	void setGenesFromMap(const map<int, vector<unsigned char>>& newGenes); // testing only
#endif

	virtual void mutate() override { (this->*_mutate_func_ptr) (); }
	virtual void inheritGenes(const bool& fromMother, QuantitativeTrait* parent, set<unsigned int> const& recomPositions, int startingChromosome) override;
	virtual float express() {
		throw runtime_error("Neutral trait shouldn't be expressed.");
		return -9999;
	}

	virtual float getAlleleValueAtLocus(short chromosome, int position) const override;
	virtual float getDomCoefAtLocus(short chromosome, int position) const override {
		return 0.0;
	}

	int countHeterozygoteLoci() const;
	bool isHeterozygoteAtLocus(int locus) const;
	bool isHeterozygoteAtIndex(int locusIndex) const { // assumes diploidy
		return genes[static_cast<size_t>(locusIndex) * 2] != genes[static_cast<size_t>(locusIndex) * 2 + 1];
	}

private:

	inline static int wildType; // default allele value, value set at construction
	const int NeutralValUpperBound = UCHAR_MAX; // alleles are char, can take value 0-255

	// <Locus position, <Allele A, Allele B>>
	// Formerly map<int, vector<unsigned char>>: one map node plus one heap-allocated
	// vector per locus (~100 bytes per locus). Now ploidy bytes per locus.
	vector<unsigned char> genes;
	short ploidy = 2;
	int locusIndexOrThrow(int position, const char* msg) const;

	// Initialisation
	void initialiseUniform(int max);

	// Immutable features, set at initialisation
	// and passed down to every subsequent trait copy
	//// Species-level trait attributes, invariant across individuals
	SpeciesTrait* pSpeciesTrait;
	//// Species-level trait functions
	void (NeutralTrait::* _mutate_func_ptr) (void);
	void (NeutralTrait::* _inherit_func_ptr) (const bool& fromMother, vector<unsigned char> const& parent, set<unsigned int> const& recomPositions, int parentChromosome);

	// Possible values for immutable functions
	//// Inheritance
	void inheritDiploid(const bool& fromMother, vector<unsigned char> const&, set<unsigned int> const& recomPositions, int parentChromosome);
	void inheritHaploid(const bool& fromMother, vector<unsigned char> const& parentMutations, set<unsigned int> const& recomPositions, int parentChromosome);
	//// Mutation
	void mutate_KAM();
	void mutate_SSM(); // single-step mutations

};

#ifdef UNIT_TESTS // for testing purposes only
map<int, vector<unsigned char>> createTestNeutralGenotype(
	const int genomeSz, const bool isDiploid,
	const unsigned char valAlleleA,
	const unsigned char valAlleleB = char(0) // if haploid
);
#endif

#endif // NeutralTraitH
