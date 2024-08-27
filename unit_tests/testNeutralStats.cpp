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

	// If the sample is too small, Fst evaluates to zero
	// More importantly, the neutral stats module still runs without error 
	{
		// 1 - if there is zero population in the sample
		{
			// Patch setup
			const int nbPatches = 2;
			
			// Create two-patches landscape
			Landscape* pLandscape = new Landscape;
			vector<Patch*> patches(nbPatches);
			vector<Cell*> cells(nbPatches);
			set<int> patchList;
			for (int i = 0; i < nbPatches; i++) {
				patches[i] = pLandscape->newPatch(i);
				cells[i] = new Cell(0, 0, (intptr)patches[i], 0);
				patches[i]->addCell(cells[i], 0, 0);
				patchList.insert(patches[i]->getPatchNum());
			}
			// Create species trait structure
			Species* pSpecies = new Species();
			const bool isDiploid = true;
			const set<int> genePositions = { 0 };
			const float maxAlleleVal = 0;
			SpeciesTrait* spTr = createTestNeutralSpTrait(maxAlleleVal, genePositions, isDiploid);
			pSpecies->addTrait(TraitType::NEUTRAL, *spTr);

			// Oops, no populations in patches

			// Compute F-stats
			const int nbLoci = 0;
			auto pNeutralStatistics = make_unique<NeutralStatsManager>(nbPatches, nbLoci);
			pNeutralStatistics->updateAllNeutralTables(
				pSpecies,
				pLandscape,
				patchList
			);
			const int maxNbNeutralAlleles = 0;
			pNeutralStatistics->calculateFstatWC(
				patchList,
				0,
				nbLoci,
				maxNbNeutralAlleles,
				pSpecies,
				pLandscape
			);
			assert(pNeutralStatistics->getFstWC() == 0.0);
		} // end case 1 - zero pop in sample

		// 2 - If there is only one population in the sample
		{
			// Patch setup
			const int nbPatches = 2;
			const int nbPops = nbPatches - 1; // oops, 1 population is missing
			const int nbIndsPerPop = 5;
			// Genetic setup
			const int genomeSz = 2;
			const bool isDiploid{ true };
			const set<int> genePositions = { 0, 1 };
			const float maxAlleleVal = 0;

			// Create two-patches landscape
			Landscape* pLandscape = new Landscape;
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

			// Create species trait structure
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
			const int nbLoci = genePositions.size();
			SpeciesTrait* spTr = createTestNeutralSpTrait(maxAlleleVal, genePositions, isDiploid);
			pSpecies->addTrait(TraitType::NEUTRAL, *spTr);

			// Initialise populations
			const int indStg = 1;
			for (int p = 0; p < nbPops; p++) {
				Population* pPop = new Population(pSpecies, patches[p], 0, 1);
				// create individuals and add to pop 
				for (int i = 0; i < nbIndsPerPop; i++) {
					Individual* pInd = new Individual(cells[p], patches[p], indStg, 0, 0, 0.0, false, 1);
					pInd->setUpGenes(pSpecies, 1.0);
					pPop->recruit(pInd);
				}
				pPop->sampleIndsWithoutReplacement(indSampling, { indStg });
			}

			// Compute F-stats
			auto pNeutralStatistics = make_unique<NeutralStatsManager>(nbPatches, nbLoci);
			pNeutralStatistics->updateAllNeutralTables(
				pSpecies,
				pLandscape,
				patchList
			);
			const int maxNbNeutralAlleles = static_cast<int>(maxAlleleVal) + 1;
			pNeutralStatistics->calculateFstatWC(
				patchList,
				nbIndsPerPop * patchList.size(),
				nbLoci,
				maxNbNeutralAlleles,
				pSpecies,
				pLandscape
			);
			assert(pNeutralStatistics->getFstWC() == 0.0);
		} // end case 2, only one population in sample
		
		// 3 - Two empty populations
		{
			// Patch setup
			const int nbPatches = 2;
			const int nbIndsPerPop = 0; // oops

			// Genetic setup
			const int genomeSz = 2;
			const bool isDiploid{ true };
			const set<int> genePositions = { 0, 1 };
			const float maxAlleleVal = 0;

			// Create two-patches landscape
			Landscape* pLandscape = new Landscape;
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

			// Create species trait structure
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
			const int nbLoci = genePositions.size();
			SpeciesTrait* spTr = createTestNeutralSpTrait(maxAlleleVal, genePositions, isDiploid);
			pSpecies->addTrait(TraitType::NEUTRAL, *spTr);

			// Initialise populations
			const int indStg = 1;
			for (int p = 0; p < nbPatches; p++) {
				Population* pPop = new Population(pSpecies, patches[p], nbIndsPerPop, 1);
				// create individuals and add to pop 
				pPop->sampleIndsWithoutReplacement(indSampling, { indStg });
			}

			// Compute F-stats
			auto pNeutralStatistics = make_unique<NeutralStatsManager>(nbPatches, nbLoci);
			pNeutralStatistics->updateAllNeutralTables(
				pSpecies,
				pLandscape,
				patchList
			);
			const int maxNbNeutralAlleles = static_cast<int>(maxAlleleVal) + 1;
			pNeutralStatistics->calculateFstatWC(
				patchList,
				nbIndsPerPop * patchList.size(),
				nbLoci,
				maxNbNeutralAlleles,
				pSpecies,
				pLandscape
			);
			assert(pNeutralStatistics->getFstWC() == 0.0);
		}
	}

	// If two sampled pops are monomorphic for different alleles,
	// then Fst is 1
	{
		// Patch setup
		const int nbPatches = 2;
		const int nbIndsPerPop = 2;
		// Genetic setup
		const int genomeSz = 2;
		const bool isDiploid{ true };
		const set<int> genePositions = { 0, 1 };
		const float maxAlleleVal = 1;
		unsigned char alleleValPopA = char(0);
		unsigned char alleleValPopB = char(1);
		auto popAGenotype = createTestNeutralGenotype(genomeSz, true, alleleValPopA, alleleValPopA);
		auto popBGenotype = createTestNeutralGenotype(genomeSz, true, alleleValPopB, alleleValPopB);

		// Create two-patches landscape
		Landscape* pLandscape = new Landscape;
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

		// Create species trait structure
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
		const int nbLoci = genePositions.size();
		SpeciesTrait* spTr = createTestNeutralSpTrait(maxAlleleVal, genePositions, isDiploid);
		pSpecies->addTrait(TraitType::NEUTRAL, *spTr);

		// Initialise populations
		vector<map<int, std::vector<unsigned char>>> genotypes = { popAGenotype, popBGenotype };
		const int indStg = 1;
		for (int p = 0; p < patches.size(); p++) {
			Population* pPop = new Population(pSpecies, patches[p], 0, 1);
			// create individuals and add to pop 
			for (int i = 0; i < nbIndsPerPop; i++) {
				Individual* pInd = new Individual(cells[p], patches[p], indStg, 0, 0, 0.0, false, 1);
				pInd->setUpGenes(pSpecies, 1.0);
				pInd->overrideGenotype(NEUTRAL, genotypes[p]);
				pPop->recruit(pInd);
			}
			pPop->sampleIndsWithoutReplacement(indSampling, { indStg });
		}
		
		// Compute F-stats
		auto pNeutralStatistics = make_unique<NeutralStatsManager>(nbPatches, nbLoci);
		pNeutralStatistics->updateAllNeutralTables(
			pSpecies, 
			pLandscape, 
			patchList
		);
		const int maxNbNeutralAlleles = static_cast<int>(maxAlleleVal) + 1;
		pNeutralStatistics->calculateFstatWC(
			patchList,
			nbIndsPerPop * patchList.size(),
			nbLoci,
			maxNbNeutralAlleles,
			pSpecies,
			pLandscape
		);
		assert(pNeutralStatistics->getFstWC() == 1.0);
		
		pNeutralStatistics->calcPairwiseWeightedFst(
			patchList,
			nbIndsPerPop* patchList.size(),
			nbLoci,
			pSpecies,
			pLandscape
		);
		assert(pNeutralStatistics->getPairwiseFst(0, 1) == 0.0);
	}

	// In strictly homozygote samples:
	// 1. Fis = 1 (maximum inbreeding)
	// 2. The sign of the Fst depends on the ratio of variation within vs between populations
	{
		// Case 1/2: variation within > between
		// Within: frequencies of allele A and B are quite different
		// Between: no variation
		{
			// Patch setup
			const int nbPatches = 2;
			const int nbIndsPerPop = 10;
			// Genetic setup
			const int genomeSz = 1;
			const bool isDiploid{ true };
			const set<int> genePositions = { 0 };
			const float maxAlleleVal = 1;
			unsigned char alleleValPopA = char(0);
			unsigned char alleleValPopB = char(1);
			const auto genotypeAA = createTestNeutralGenotype(genomeSz, true, alleleValPopA, alleleValPopA);
			const auto genotypeBB = createTestNeutralGenotype(genomeSz, true, alleleValPopB, alleleValPopB);
			vector<map<int, vector<unsigned char>>> genotypes = {
				// 8 AA, 2 BB (no heterozygotes)
				genotypeAA, genotypeAA, genotypeAA, genotypeAA, genotypeAA,
				genotypeAA, genotypeAA, genotypeAA, genotypeBB, genotypeBB
			};

			// Create two-patches landscape
			Landscape* pLandscape = new Landscape;
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

			// Create species trait structure
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
			const int nbLoci = genePositions.size();
			SpeciesTrait* spTr = createTestNeutralSpTrait(maxAlleleVal, genePositions, isDiploid);
			pSpecies->addTrait(TraitType::NEUTRAL, *spTr);

			// Initialise populations
			const int indStg = 1;
			for (int p = 0; p < patches.size(); p++) {
				Population* pPop = new Population(pSpecies, patches[p], 0, 1);
				// create individuals and add to pop 
				for (int i = 0; i < nbIndsPerPop; i++) {
					Individual* pInd = new Individual(cells[p], patches[p], indStg, 0, 0, 0.0, false, 1);
					pInd->setUpGenes(pSpecies, 1.0);
					pInd->overrideGenotype(NEUTRAL, genotypes[i]);
					pPop->recruit(pInd);
				}
				pPop->sampleIndsWithoutReplacement(indSampling, { indStg });
			}

			// Compute F-stats
			auto pNeutralStatistics = make_unique<NeutralStatsManager>(nbPatches, nbLoci);
			pNeutralStatistics->updateAllNeutralTables(
				pSpecies,
				pLandscape,
				patchList
			);
			const int maxNbNeutralAlleles = static_cast<int>(maxAlleleVal) + 1;
			pNeutralStatistics->calculateFstatWC(
				patchList,
				nbIndsPerPop * patchList.size(),
				nbLoci,
				maxNbNeutralAlleles,
				pSpecies,
				pLandscape
			);
			assert(pNeutralStatistics->getFstWC() < 0.0);
			assert(pNeutralStatistics->getFisWC() == 1.0);
		}

		// Case 2/2: variation within < between
		// Within: no variation for pop1, same allelic frequencies
		// Between: larger variation
		{
			// Patch setup
			const int nbPatches = 2;
			const int nbIndsPerPop = 10;
			// Genetic setup
			const int genomeSz = 1;
			const bool isDiploid{ true };
			const set<int> genePositions = { 0 };
			const float maxAlleleVal = 1;
			unsigned char alleleValPopA = char(0);
			unsigned char alleleValPopB = char(1);
			const auto genotypeAA = createTestNeutralGenotype(genomeSz, true, alleleValPopA, alleleValPopA);
			const auto genotypeBB = createTestNeutralGenotype(genomeSz, true, alleleValPopB, alleleValPopB);
			vector<vector<map<int, vector<unsigned char>>>> genotypeList = {
				{ // Pop 1: 5 AA, 5 BB (no heterozygotes)
					genotypeAA, genotypeAA, genotypeAA, genotypeAA, genotypeAA,
					genotypeBB, genotypeBB, genotypeBB, genotypeBB, genotypeBB
				}, 
				{ // Pop 2: 8 AA, 2 BB (no heterozygotes)
					genotypeAA, genotypeAA, genotypeAA, genotypeAA, genotypeAA,
					genotypeAA, genotypeAA, genotypeAA, genotypeBB, genotypeBB
				}
			};

			// Create two-patches landscape
			Landscape* pLandscape = new Landscape;
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

			// Create species trait structure
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
			const int nbLoci = genePositions.size();
			SpeciesTrait* spTr = createTestNeutralSpTrait(maxAlleleVal, genePositions, isDiploid);
			pSpecies->addTrait(TraitType::NEUTRAL, *spTr);

			// Initialise populations
			const int indStg = 1;
			for (int p = 0; p < patches.size(); p++) {
				Population* pPop = new Population(pSpecies, patches[p], 0, 1);
				// create individuals and add to pop 
				for (int i = 0; i < nbIndsPerPop; i++) {
					Individual* pInd = new Individual(cells[p], patches[p], indStg, 0, 0, 0.0, false, 1);
					pInd->setUpGenes(pSpecies, 1.0);
					pInd->overrideGenotype(NEUTRAL, genotypeList[p][i]);
					pPop->recruit(pInd);
				}
				pPop->sampleIndsWithoutReplacement(indSampling, { indStg });
			}

			// Compute F-stats
			auto pNeutralStatistics = make_unique<NeutralStatsManager>(nbPatches, nbLoci);
			pNeutralStatistics->updateAllNeutralTables(
				pSpecies,
				pLandscape,
				patchList
			);
			const int maxNbNeutralAlleles = static_cast<int>(maxAlleleVal) + 1;
			pNeutralStatistics->calculateFstatWC(
				patchList,
				nbIndsPerPop * patchList.size(),
				nbLoci,
				maxNbNeutralAlleles,
				pSpecies,
				pLandscape
			);
			assert(pNeutralStatistics->getFstWC() > 0.0);
			assert(pNeutralStatistics->getFisWC() == 1.0);
		}
	}

	// Fst calculation is correct for an ordinary sample
	{
		// Two 10-individual pops,
		// 1 tri-allelic locus

		// Patch setup
		const int nbPatches = 2;
		const int nbIndsPerPop = 10;

		// Genetic setup
		const int genomeSz = 1;
		const bool isDiploid{ true };
		const set<int> genePositions = { 0 };
		const float maxAlleleVal = 2;
		unsigned char alleleA = char(0);
		unsigned char alleleB = char(1);
		unsigned char alleleC = char(2);
		// Create genotypes
		auto genotypeAA = createTestNeutralGenotype(genomeSz, true, alleleA, alleleA);
		auto genotypeAB = createTestNeutralGenotype(genomeSz, true, alleleA, alleleB);
		auto genotypeAC = createTestNeutralGenotype(genomeSz, true, alleleA, alleleC);
		auto genotypeBB = createTestNeutralGenotype(genomeSz, true, alleleB, alleleB);
		auto genotypeBC = createTestNeutralGenotype(genomeSz, true, alleleB, alleleC);
		auto genotypeCC = createTestNeutralGenotype(genomeSz, true, alleleC, alleleC);
		vector<vector<map<int, std::vector<unsigned char>>>> genotypeList = {
			{ // Genotypes of population 1
				genotypeAA, genotypeAA, genotypeAB, genotypeAB, genotypeAB,
				genotypeAB, genotypeAC, genotypeAC, genotypeBB, genotypeCC
			},
			{ // Genotypes of population 2
				genotypeAA, genotypeAB, genotypeAC, genotypeBB, genotypeBC,
				genotypeBC, genotypeBC, genotypeBC, genotypeCC, genotypeCC
			}
		};

		// Create two-patches landscape
		Landscape* pLandscape = new Landscape;
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

		// Create species trait structure
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
		const int nbLoci = genePositions.size();
		SpeciesTrait* spTr = createTestNeutralSpTrait(maxAlleleVal, genePositions, isDiploid);
		pSpecies->addTrait(TraitType::NEUTRAL, *spTr);

		// Initialise populations
		const int indStg = 1;
		for (int p = 0; p < patches.size(); p++) {
			Population* pPop = new Population(pSpecies, patches[p], 0, 1);
			// create individuals and add to pop 
			for (int i = 0; i < nbIndsPerPop; i++) {
				Individual* pInd = new Individual(cells[p], patches[p], indStg, 0, 0, 0.0, false, 1);
				pInd->setUpGenes(pSpecies, 1.0);
				pInd->overrideGenotype(NEUTRAL, genotypeList[p][i]);
				pPop->recruit(pInd);
			}
			pPop->sampleIndsWithoutReplacement(indSampling, { indStg });
		}

		// Compute F-stats
		auto pNeutralStatistics = make_unique<NeutralStatsManager>(nbPatches, nbLoci);
		pNeutralStatistics->updateAllNeutralTables(
			pSpecies,
			pLandscape,
			patchList
		);
		const int maxNbNeutralAlleles = static_cast<int>(maxAlleleVal) + 1;
		pNeutralStatistics->calculateFstatWC(
			patchList,
			nbIndsPerPop * patchList.size(),
			nbLoci,
			maxNbNeutralAlleles,
			pSpecies,
			pLandscape
		);
		const double expectedFst = 0.0583; // calculated by hand from Weir and Cockerham 1984
		double calcError = abs(pNeutralStatistics->getFstWC() - expectedFst);
		assert(calcError < 0.0001);
	}
}

#endif //RSDEBUG