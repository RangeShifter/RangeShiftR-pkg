#if RSDEBUG

#include "../Community.h"

void testNeutralStats() {

	// In haploid settings, Ho is always zero
	{
	const int genomeSz = 4;
	Species* pSpecies = new Species();
	pSpecies->setGeneticParameters(
		set<int>{genomeSz - 1}, // single chromosome
		genomeSz,
		0.0, // no recombination
		set<int>{}, "none", set<int>{}, 0 // no output so no sampling
	);
	const set<int> genePositions = { 0, 1, 3 }; // arbitrary
	const bool isDiploid{ false };
	SpeciesTrait* spTr = createTestNeutralSpTrait(genePositions, isDiploid);
	pSpecies->addTrait(TraitType::NEUTRAL, *spTr);

	const int nbInds = 2;
	for (int i = 0; i < nbInds; i++) {
		Individual ind = Individual(pCell, pPatch, 0, 0, 0, 0.0, false, 0);
		ind.setUpGenes(pSpecies, 1.0);
		pPop->recruit(&ind);
	}

		// create a patch with populations with individuals with traits

		auto pNeutralStatistics = make_unique<NeutralStatsManager>(patchList.size(), nLoci);
		pNeutralStatistics->calculateHo(
			//set<int> const& patchList, 
			//const int nbInds, 
			//const int nbrLoci, 
			//Species * pSpecies, 
			//Landscape * pLandscape
		);
		assert(pNeutralStatistics->getHo() == 0.0);
	}
}

#endif //RSDEBUG