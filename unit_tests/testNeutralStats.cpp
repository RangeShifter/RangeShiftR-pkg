#if RSDEBUG

#include "../Community.h"

void testNeutralStats() {

	// In haploid settings, Ho is always zero
	{
		// Create genetic structure
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
		const float maxAlleleVal = 255; // arbitrary
		SpeciesTrait* spTr = createTestNeutralSpTrait(maxAlleleVal, genePositions, isDiploid);
		pSpecies->addTrait(TraitType::NEUTRAL, *spTr);

		// Create empty 1-patch, 1-cell landscape
		Landscape* pLandscape = new Landscape;
		Patch* pPatch = pLandscape->newPatch(0);
		Cell* pCell = new Cell(0, 0, (intptr)pPatch, 0);
		pPatch->addCell(pCell, 0, 0);

		// Initialise population
		const int nbInds = 2;
		Population* pPop = new Population(pSpecies, pPatch, nbInds, 1);

		// Calculate heterozygosity
		set<int> patchList{ pPatch->getPatchNum()};
		auto pNeutralStatistics = make_unique<NeutralStatsManager>(patchList.size(), genePositions.size());
		pNeutralStatistics->calculateHo(
			patchList,
			nbInds, 
			genePositions.size(),
			pSpecies, 
			pLandscape
		);
		assert(pNeutralStatistics->getHo() == 0.0);
	}

	// If every individual in a sample is homozygote, Ho is zero
	{
		// Create genetic structure
		const int genomeSz = 4;
		Species* pSpecies = new Species();
		pSpecies->setGeneticParameters(
			set<int>{genomeSz - 1}, // single chromosome
			genomeSz,
			0.0, // no recombination
			set<int>{}, "none", set<int>{}, 0 // no output so no sampling
		);
		const set<int> genePositions = { 0, 1, 3 }; // arbitrary
		const bool isDiploid{ true };
		const float maxAlleleVal = 0; // sample in uniform(0,0) so every individual is homozygote
		SpeciesTrait* spTr = createTestNeutralSpTrait(maxAlleleVal, genePositions, isDiploid);
		pSpecies->addTrait(TraitType::NEUTRAL, *spTr);

		// Create empty 1-patch, 1-cell landscape
		Landscape* pLandscape = new Landscape;
		Patch* pPatch = pLandscape->newPatch(0);
		Cell* pCell = new Cell(0, 0, (intptr)pPatch, 0);
		pPatch->addCell(pCell, 0, 0);

		// Initialise population
		const int nbInds = 2;
		Population* pPop = new Population(pSpecies, pPatch, nbInds, 1);

		// Calculate heterozygosity
		set<int> patchList{ pPatch->getPatchNum() };
		auto pNeutralStatistics = make_unique<NeutralStatsManager>(patchList.size(), genePositions.size());
		pNeutralStatistics->calculateHo(
			patchList,
			nbInds,
			genePositions.size(),
			pSpecies,
			pLandscape
		);
		assert(pNeutralStatistics->getHo() == 0.0);
	}
}

#endif //RSDEBUG