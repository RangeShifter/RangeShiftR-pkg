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
	const float maxAlleleVal = 0.0; // sample from zero to zero, i.e. all same allele
	SpeciesTrait* spTr = createTestNeutralSpTrait(maxAlleleVal, genePositions, isDiploid);
	pSpecies->addTrait(TraitType::NEUTRAL, *spTr);

	// not sure if patch zero (matrix) is ok
	Patch* pPatch = new Patch(0, 0);
	Cell* pCell = new Cell(0, 0, (intptr)pPatch, 0);
	pPatch->addCell(pCell, 0, 0);

	const int nbInds = 2;
	Population* pPop = new Population(pSpecies, pPatch, nbInds, 1);
	/*
	for (int i = 0; i < nbInds; i++) {
		Individual ind = Individual(pCell, pPatch, 0, 0, 0, 0.0, false, 0);
		ind.setUpGenes(pSpecies, 1.0);
		pPop->recruit(&ind);
	}

		auto pNeutralStatistics = make_unique<NeutralStatsManager>(patchList.size(), nLoci);
		pNeutralStatistics->calculateHo(
			//set<int> const& patchList, 
			//const int nbInds, 
			//const int nbrLoci, 
			//Species * pSpecies, 
			//Landscape * pLandscape
		);
		assert(pNeutralStatistics->getHo() == 0.0);
		*/
	}
}

#endif //RSDEBUG