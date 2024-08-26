#if RSDEBUG

#include "../Community.h"

demogrParams createDefltHaploidDemogrParams() {
	demogrParams d;
	d.repType = 0;
	d.repSeasons = 1;
	d.stageStruct = false;
	d.propMales = 0.5;
	d.harem = 1.0;
	d.bc = 1.0;
	d.lambda = 1.0;
	return d;
}

demogrParams createDefaultDiploidDemogrParams() {
	demogrParams d;
	d.repType = 1;
	d.repSeasons = 1;
	d.stageStruct = false;
	d.propMales = 0.5;
	d.harem = 1.0;
	d.bc = 1.0;
	d.lambda = 1.0;
	return d;
}

void testNeutralStats() {

	// In haploid settings, Ho is always zero
	{
		// Create empty 1-patch, 1-cell landscape
		Landscape* pLandscape = new Landscape;
		Patch* pPatch = pLandscape->newPatch(0);
		const set<int> patchList{ pPatch->getPatchNum() };
		const string indSampling = "all";
		const set<int> stgToSample = { 1 };
		Cell* pCell = new Cell(0, 0, (intptr)pPatch, 0);
		pPatch->addCell(pCell, 0, 0);

		// Create genetic structure
		const int genomeSz = 4;
		Species* pSpecies = new Species();
		pSpecies->setDemogr(createDefltHaploidDemogrParams());
		pSpecies->setGeneticParameters(
			set<int>{genomeSz - 1}, // single chromosome
			genomeSz,
			0.0, // no recombination
			patchList,
			indSampling,
			stgToSample,
			1
		);
		const set<int> genePositions = { 0, 1, 3 }; // arbitrary
		const bool isDiploid{ false };
		const float maxAlleleVal = 0; // sample in uniform(0,0) so every individual is homozygote
		SpeciesTrait* spTr = createTestNeutralSpTrait(maxAlleleVal, genePositions, isDiploid);
		pSpecies->addTrait(TraitType::NEUTRAL, *spTr);

		// Initialise population
		const int nbInds = 2;
		Population* pPop = new Population(pSpecies, pPatch, nbInds, 1);
		pPop->sampleIndsWithoutReplacement(indSampling, stgToSample);

		// Calculate heterozygosity
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
		// Create empty 1-patch, 1-cell landscape
		Landscape* pLandscape = new Landscape;
		Patch* pPatch = pLandscape->newPatch(0);
		const set<int> patchList{ pPatch->getPatchNum() };
		const string indSampling = "all";
		const set<int> stgToSample = { 1 };
		Cell* pCell = new Cell(0, 0, (intptr)pPatch, 0);
		pPatch->addCell(pCell, 0, 0);

		// Create genetic structure
		const int genomeSz = 4;
		Species* pSpecies = new Species();
		pSpecies->setDemogr(createDefaultDiploidDemogrParams());
		pSpecies->setGeneticParameters(
			set<int>{genomeSz - 1}, // single chromosome
			genomeSz,
			0.0, // no recombination
			patchList, 
			indSampling, 
			stgToSample,
			1
		);
		const set<int> genePositions = { 0, 1, 3 }; // arbitrary
		const bool isDiploid{ true };
		const float maxAlleleVal = 0; // sample in uniform(0,0) so every individual is homozygote
		SpeciesTrait* spTr = createTestNeutralSpTrait(maxAlleleVal, genePositions, isDiploid);
		pSpecies->addTrait(TraitType::NEUTRAL, *spTr);

		// Initialise population
		const int nbInds = 2;
		Population* pPop = new Population(pSpecies, pPatch, nbInds, 1);
		pPop->sampleIndsWithoutReplacement(indSampling, stgToSample);

		// Calculate heterozygosity
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

	// If every individual in a sample is heterozygote, Ho is one.
	{
		// Create empty 1-patch, 1-cell landscape
		Landscape* pLandscape = new Landscape;
		Patch* pPatch = pLandscape->newPatch(0);
		const set<int> patchList{ pPatch->getPatchNum() };
		const string indSampling = "all";
		const set<int> stgToSample = { 1 };
		Cell* pCell = new Cell(0, 0, (intptr)pPatch, 0);
		pPatch->addCell(pCell, 0, 0);

		// Create genetic structure
		const int genomeSz = 4;
		Species* pSpecies = new Species();
		pSpecies->setDemogr(createDefaultDiploidDemogrParams());
		pSpecies->setGeneticParameters(
			set<int>{genomeSz - 1}, // single chromosome
			genomeSz,
			0.0, // no recombination
			patchList,
			indSampling,
			stgToSample,
			1
		);
		const set<int> genePositions = { 0, 1, 3 }; // arbitrary
		const bool isDiploid{ true };
		const float maxAlleleVal = 255; // highly unlikely that same value sampled twice
		SpeciesTrait* spTr = createTestNeutralSpTrait(maxAlleleVal, genePositions, isDiploid);
		pSpecies->addTrait(TraitType::NEUTRAL, *spTr);

		// Initialise population
		const int nbInds = 2;
		Population* pPop = new Population(pSpecies, pPatch, nbInds, 1);
		pPop->sampleIndsWithoutReplacement(indSampling, stgToSample);

		// Calculate heterozygosity
		auto pNeutralStatistics = make_unique<NeutralStatsManager>(patchList.size(), genePositions.size());
		pNeutralStatistics->calculateHo(
			patchList,
			nbInds,
			genePositions.size(),
			pSpecies,
			pLandscape
		);
		assert(pNeutralStatistics->getHo() == 1.0);
	}

	// If two sampled pops have zero allele in common, Fst is 1
	// - global Fst
	// - per-locus Fst
	// - pairwise per-patch Fst
	{
		// Create two-patches landscape
		Landscape* pLandscape = new Landscape;
		const int nbPatches = 2;
		vector<Patch*> patches(nbPatches);
		vector<Cell*> cells(nbPatches);
		set<int> patchList;
		for (int i = 0; i < nbPatches; i++) {
			patches[i] = pLandscape->newPatch(i);
			cells[i] = new Cell(0, 0, (intptr)patches[i], 0);
			patches[i]->addCell(cells[i], 0, 0);
			patchList.insert(patches[i]->getPatchNum());
		}
		const string indSampling = "all";
		const set<int> stgToSample = { 1 };

		// Create genetic structure
		const int genomeSz = 2;
		Species* pSpecies = new Species();
		pSpecies->setDemogr(createDefaultDiploidDemogrParams());
		pSpecies->setGeneticParameters(
			set<int>{genomeSz - 1}, // single chromosome
			genomeSz,
			0.0, // no recombination
			patchList,
			indSampling,
			stgToSample,
			1
		);
		const set<int> genePositions = { 0, 1 };
		const int nbLoci = genePositions.size();
		const bool isDiploid{ true };
		const float maxAlleleVal = 255; // highly unlikely that same value sampled twice
		SpeciesTrait* spTr = createTestNeutralSpTrait(maxAlleleVal, genePositions, isDiploid);
		pSpecies->addTrait(TraitType::NEUTRAL, *spTr);

		const int nbIndsPerPop = 2;
		unsigned char alleleValPopA = char(0);
		unsigned char alleleValPopB = char(1);
		auto popAGenotype = createTestNeutralGenotype(genomeSz, true, alleleValPopA, alleleValPopA);
		auto popBGenotype = createTestNeutralGenotype(genomeSz, true, alleleValPopB, alleleValPopB);
		vector<
			map<int,std::vector<unsigned char>>
		> genotypes = { popAGenotype, popBGenotype };
		
		// Initialise populations
		for (int p = 0; p < patches.size(); p++) {
			Population* pPop = new Population(pSpecies, patches[p], 0, 1);
			// create individuals and add to pop 
			for (int i = 0; i < nbIndsPerPop; i++) {
				Individual ind = Individual(cells[p], patches[p], 0, 0, 0, 0.0, false, 1);
				ind.setUpGenes(pSpecies, 1.0);
				ind.overrideGenotype(NEUTRAL, genotypes[p]);
				pPop->recruit(&ind);
			}
			pPop->sampleIndsWithoutReplacement(indSampling, stgToSample);
		} // do individuals still exist?
		
		// Calculate heterozygosity
		auto pNeutralStatistics = make_unique<NeutralStatsManager>(nbPatches, nbLoci);
		pNeutralStatistics->updateAllNeutralTables(
			pSpecies, 
			pLandscape, 
			patchList
		);
		const int maxNbNeutralAlleles = static_cast<int>(maxAlleleVal);
		pNeutralStatistics->calculateFstatWC(
			patchList,
			nbIndsPerPop * patchList.size(),
			nbLoci,
			maxNbNeutralAlleles,
			pSpecies,
			pLandscape
		);
		// assert(pNeutralStatistics->getFstWC() == 0.4);
	}
}

#endif //RSDEBUG