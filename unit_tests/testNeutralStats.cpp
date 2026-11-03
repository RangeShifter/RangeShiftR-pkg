#ifdef UNIT_TESTS

#include "../Community.h"

// ----------------------------------------------------------------------------------------
// Per-patch Fis
// ----------------------------------------------------------------------------------------

// Genotype of one individual: one (alleleA, alleleB) pair per locus (alleleB ignored if haploid)
using TestGenotype = vector<pair<int, int>>;

struct PerPatchFisSetup {
	Landscape* pLandscape;
	Species* pSpecies;
	vector<Patch*> patches;
	set<int> patchList;
	int nbInds = 0;
	int nbLoci = 0;
};

// One patch per entry of genotypesPerPatch, loci at genome positions 0 .. nbLoci-1
PerPatchFisSetup makePerPatchFisSetup(const vector<vector<TestGenotype>>& genotypesPerPatch, const int nbLoci, const bool isDiploid) {
	PerPatchFisSetup s;
	s.nbLoci = nbLoci;
	s.pLandscape = new Landscape;
	for (int i = 0; i < (int)genotypesPerPatch.size(); i++) {
		Patch* p = s.pLandscape->newPatch(i + 1);
		Cell* c = new Cell(i, 0, p, 0);
		p->addCell(c, i, 0);
		s.patches.push_back(p);
		s.patchList.insert(p->getPatchNum());
	}
	s.pSpecies = createDefaultSpecies();
	s.pSpecies->setDemogr(isDiploid ? createDefaultDiploidDemogrParams() : createDefaultHaploidDemogrParams());
	s.pSpecies->setGeneticParameters(set<int>{nbLoci - 1}, nbLoci, 0.0, s.patchList, "all", { 1 }, 1);
	set<int> positions;
	for (int l = 0; l < nbLoci; l++) positions.insert(l);
	const float maxAlleleVal = 3; // alleles 0..3
	SpeciesTrait* spTr = createTestNeutralSpTrait(maxAlleleVal, positions, isDiploid);
	s.pSpecies->addTrait(TraitType::NEUTRAL, *spTr);

	for (int p = 0; p < (int)genotypesPerPatch.size(); p++) {
		Population* pPop = new Population(s.pSpecies, s.patches[p], 0, 1);
		for (const TestGenotype& gt : genotypesPerPatch[p]) {
			Individual* pInd = new Individual(s.pSpecies, s.patches[p]->getCell(0), s.patches[p], 1, 0, 0, 0.0, false, 1);
			pInd->setUpGenes(s.pSpecies, 1);
			map<int, vector<unsigned char>> genotype;
			for (int l = 0; l < nbLoci; l++) {
				vector<unsigned char> gene(isDiploid ? 2 : 1);
				gene[0] = (unsigned char)gt[l].first;
				if (isDiploid) gene[1] = (unsigned char)gt[l].second;
				genotype.emplace(l, gene);
			}
			pInd->overrideGenotype(NEUTRAL, genotype);
			pPop->recruit(pInd);
			s.nbInds++;
		}
		pPop->sampleIndsWithoutReplacement("all", { 1 });
	}
	return s;
}

vector<PatchFisStats> computePerPatchFis(PerPatchFisSetup& s) {
	NeutralStatsManager nsm(s.patchList.size(), s.nbLoci);
	nsm.updateAllNeutralTables(s.pSpecies, s.pLandscape, s.patchList);
	nsm.calculatePerPatchFis(s.patchList, s.pSpecies, s.pLandscape);
	return nsm.getPerPatchFis();
}

void testPerPatchFis() {

	const double tol = 1e-9;

	// Hand-computed example, one locus, n = 4:
	// genotypes 0/1, 0/0, 0/1, 1/1 -> p = 0.5, Ho = 0.5,
	// Hs = 4/3 * (1 - 0.5 - 0.5/8) = 7/12, Fis = 1 - 0.5 / (7/12) = 1/7
	{
		auto s = makePerPatchFisSetup({ { {{0,1}}, {{0,0}}, {{0,1}}, {{1,1}} } }, 1, true);
		auto res = computePerPatchFis(s);
		assert(res.size() == 1);
		assert(res[0].patchId == s.patches[0]->getPatchNum());
		assert(res[0].nbSampledInds == 4);
		assert(fabs(res[0].ho - 0.5) < tol);
		assert(fabs(res[0].hs - 7.0 / 12.0) < tol);
		assert(fabs(res[0].fis - 1.0 / 7.0) < tol);
	}

	// Per-patch Fis is the Weir & Cockerham f-hat for a single sample:
	// with two identical patches (s^2 = 0), the global WC Fis must equal the per-patch Fis
	{
		const vector<TestGenotype> pop = { {{0,1},{2,2}}, {{0,0},{1,2}}, {{1,1},{1,1}}, {{0,1},{2,1}}, {{0,0},{2,2}} };
		auto s = makePerPatchFisSetup({ pop, pop }, 2, true);
		NeutralStatsManager nsm(s.patchList.size(), s.nbLoci);
		nsm.updateAllNeutralTables(s.pSpecies, s.pLandscape, s.patchList);
		nsm.calculateFstatWC(s.patchList, s.nbInds, s.nbLoci, 4, s.pSpecies, s.pLandscape, false);
		nsm.calculatePerPatchFis(s.patchList, s.pSpecies, s.pLandscape);
		const auto& res = nsm.getPerPatchFis();
		assert(res.size() == 2);
		assert(fabs(res[0].fis - res[1].fis) < tol);
		assert(fabs(res[0].fis - nsm.getFisWC()) < tol);
	}

	// Multi-locus Fis is the ratio of sums over loci, not the mean of per-locus ratios.
	// Locus 0: 0/1, 0/1, 0/1 -> p = 0.5, Ho = 1,   Hs = 3/2 * (0.5 - 1/6)  = 0.5
	// Locus 1: 0/0, 1/1, 0/1 -> p = 0.5, Ho = 1/3, Hs = 3/2 * (0.5 - 1/18) = 2/3
	// Fis = 1 - (4/3) / (7/6) = -1/7
	{
		auto s = makePerPatchFisSetup({ { {{0,1},{0,0}}, {{0,1},{1,1}}, {{0,1},{0,1}} } }, 2, true);
		auto res = computePerPatchFis(s);
		assert(fabs(res[0].ho - 2.0 / 3.0) < tol);
		assert(fabs(res[0].hs - 7.0 / 12.0) < tol);
		assert(fabs(res[0].fis - (-1.0 / 7.0)) < tol);
	}

	// Edge cases: monomorphic patch (Fis undefined), single individual (Hs, Fis undefined),
	// patch without sampled individuals (skipped)
	{
		auto s = makePerPatchFisSetup({ { {{0,0}}, {{0,0}} }, { {{0,1}} }, {} }, 1, true);
		auto res = computePerPatchFis(s);
		assert(res.size() == 2); // empty patch skipped
		assert(res[0].ho == 0.0 && res[0].hs == 0.0);
		assert(std::isnan(res[0].fis));
		assert(res[1].nbSampledInds == 1);
		assert(res[1].ho == 1.0);
		assert(std::isnan(res[1].hs));
		assert(std::isnan(res[1].fis));
	}

	// Haploid: Ho and Fis undefined, Hs is the unbiased gene diversity n/(n-1) * (1 - sum p^2)
	{
		auto s = makePerPatchFisSetup({ { {{0,0}}, {{1,0}}, {{1,0}}, {{2,0}} } }, 1, false);
		auto res = computePerPatchFis(s);
		assert(std::isnan(res[0].ho));
		assert(std::isnan(res[0].fis));
		// p = (0.25, 0.5, 0.25) -> 1 - sum p^2 = 0.625 -> Hs = 4/3 * 0.625
		assert(fabs(res[0].hs - 4.0 / 3.0 * 0.625) < tol);
	}
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
		Cell* pCell = new Cell(0, 0, pPatch, 0);
		pPatch->addCell(pCell, 0, 0);

		// Create genetic structure
		const int genomeSz = 4;
		Species* pSpecies = createDefaultSpecies();
		pSpecies->setDemogr(createDefaultHaploidDemogrParams());
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
		Cell* pCell = new Cell(0, 0, pPatch, 0);
		pPatch->addCell(pCell, 0, 0);

		// Create genetic structure
		const int genomeSz = 4;
		Species* pSpecies = createDefaultSpecies();
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
		Cell* pCell = new Cell(0, 0, pPatch, 0);
		pPatch->addCell(pCell, 0, 0);

		// Create genetic structure
		const int genomeSz = 4;
		Species* pSpecies = createDefaultSpecies();
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
		const set<int> genePositions = { 0, 2 }; // arbitrary
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

	// The correct number of individuals is sampled
	{
		const int nbIndsPopA = 15;
		const int nbIndsPopB = 11;
		const string indSampling = "13";
		const set<int> stgToSample = { 1 };

		// Patch setup
		const int nbPatches = 2;
		// Genetic setup
		const int genomeSz = 2;
		const bool isDiploid{ true };
		const set<int> genePositions = { 0, 1 };
		const float maxAlleleVal = 1;
		
		// Create two-patches landscape
		Landscape* pLandscape = new Landscape;
		vector<Patch*> patches(nbPatches);
		vector<Cell*> cells(nbPatches);
		set<int> patchList;
		for (int i = 0; i < nbPatches; i++) {
			patches[i] = pLandscape->newPatch(i);
			cells[i] = new Cell(0, 0, patches[i], 0);
			patches[i]->addCell(cells[i], 0, 0);
			patchList.insert(patches[i]->getPatchNum());
		}
		
		// Create species trait structure
		Species* pSpecies = createDefaultSpecies();
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
		Population* pPopA = new Population(pSpecies, patches[0], nbIndsPopA, 1);
		Population* pPopB = new Population(pSpecies, patches[1], nbIndsPopB, 1);
		pPopA->sampleIndsWithoutReplacement(indSampling, { indStg });
		pPopB->sampleIndsWithoutReplacement(indSampling, { indStg });
		
		assert(pPopA->sampleSize() == stoi(indSampling)); // the correct nb of individuals has been sampled
		assert(pPopB->sampleSize() == nbIndsPopB); // too small; all individuals have been sampled
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
				cells[i] = new Cell(0, 0, patches[i], 0);
				patches[i]->addCell(cells[i], 0, 0);
				patchList.insert(patches[i]->getPatchNum());
			}
			// Create species trait structure
			Species* pSpecies = createDefaultSpecies();
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
				pLandscape,
				false
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
				cells[i] = new Cell(0, 0, patches[i], 0);
				patches[i]->addCell(cells[i], 0, 0);
				patchList.insert(patches[i]->getPatchNum());
			}
			const string indSampling = "all";
			const set<int> stgToSample = { 1 };

			// Create species trait structure
			Species* pSpecies = createDefaultSpecies();
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
					Individual* pInd = new Individual(pSpecies, cells[p], patches[p], indStg, 0, 0, 0.0, false, 1);
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
				pLandscape,
				false
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
				cells[i] = new Cell(0, 0, patches[i], 0);
				patches[i]->addCell(cells[i], 0, 0);
				patchList.insert(patches[i]->getPatchNum());
			}
			const string indSampling = "all";
			const set<int> stgToSample = { 1 };

			// Create species trait structure
			Species* pSpecies = createDefaultSpecies();
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
				pLandscape,
				false
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
			cells[i] = new Cell(0, 0, patches[i], 0);
			patches[i]->addCell(cells[i], 0, 0);
			patchList.insert(patches[i]->getPatchNum());
		}
		const string indSampling = "all";
		const set<int> stgToSample = { 1 };

		// Create species trait structure
		Species* pSpecies = createDefaultSpecies();
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
				Individual* pInd = new Individual(pSpecies, cells[p], patches[p], indStg, 0, 0, 0.0, false, 1);
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
			pLandscape,
			false
		);
		assert(pNeutralStatistics->getFstWC() == 1.0);
		
//		pNeutralStatistics->calcPairwiseWeightedFst(
//			patchList,
//			nbIndsPerPop* patchList.size(),
//			nbLoci,
//			pSpecies,
//			pLandscape
//		);
//		assert(pNeutralStatistics->getPairwiseFst(0, 1) == 0.0);
	}

	double refWeirCockerhamDiploidFst; // for use in further tests below

	// In strictly homozygote samples:
	// 1. Fis = 1 (maximum inbreeding)
	// 2. The sign of the Fst depends on the ratio of variation within vs between populations
	// 3. (if sample sizes are equal) The Weir&Cockerham 1984, Weir&Hill 2002 estimators yield the same values
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
				genotypeAA, genotypeAA, genotypeAA, genotypeAA, genotypeBB,
				genotypeAA, genotypeAA, genotypeAA, genotypeAA, genotypeBB
			};

			// Create two-patches landscape
			Landscape* pLandscape = new Landscape;
			vector<Patch*> patches(nbPatches);
			vector<Cell*> cells(nbPatches);
			set<int> patchList;
			for (int i = 0; i < nbPatches; i++) {
				patches[i] = pLandscape->newPatch(i);
				cells[i] = new Cell(0, 0, patches[i], 0);
				patches[i]->addCell(cells[i], 0, 0);
				patchList.insert(patches[i]->getPatchNum());
			}
			const string indSampling = "all";
			const set<int> stgToSample = { 1 };

			// Create species trait structure
			Species* pSpecies = createDefaultSpecies();
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
					Individual* pInd = new Individual(pSpecies, cells[p], patches[p], indStg, 0, 0, 0.0, false, 1);
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
				pLandscape,
				false
					);
			assert(pNeutralStatistics->getFstWC() < 0.0);
			assert(pNeutralStatistics->getFisWC() == 1.0);

	//		pNeutralStatistics->calcPairwiseWeightedFst(
	//			patchList,
	//			nbIndsPerPop * patchList.size(),
	//			nbLoci,
	//			pSpecies,
	//			pLandscape
	//		);
	//		const double tol = 0.000001;
	//		assert(abs(pNeutralStatistics->getWeightedFst() - pNeutralStatistics->getFstWC()) < tol);

			refWeirCockerhamDiploidFst = pNeutralStatistics->getFstWC(); // for use in further tests below
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
				cells[i] = new Cell(0, 0, patches[i], 0);
				patches[i]->addCell(cells[i], 0, 0);
				patchList.insert(patches[i]->getPatchNum());
			}
			const string indSampling = "all";
			const set<int> stgToSample = { 1 };

			// Create species trait structure
			Species* pSpecies = createDefaultSpecies();
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
					Individual* pInd = new Individual(pSpecies, cells[p], patches[p], indStg, 0, 0, 0.0, false, 1);
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
				nbIndsPerPop* patchList.size(),
				nbLoci,
				maxNbNeutralAlleles,
				pSpecies,
				pLandscape,
				false
			);
			assert(pNeutralStatistics->getFstWC() > 0.0);
			assert(pNeutralStatistics->getFisWC() == 1.0);

	//		// Weir & Hill population-specific estimates average to the (Weir & Hill) global estimator
	//		pNeutralStatistics->calcPairwiseWeightedFst(
	//			patchList,
	//			nbIndsPerPop* patchList.size(),
	//			nbLoci,
	//			pSpecies,
	//			pLandscape
	//		);
	//		const double pop1Fst = pNeutralStatistics->getPairwiseFst(0, 0);
	//		const double pop2Fst = pNeutralStatistics->getPairwiseFst(1, 1);
	//		assert((pop1Fst + pop2Fst) / 2.0 == pNeutralStatistics->getWeightedFst());

		}
	}

	// In a strictly heterozygote sample, 
	// 1. Fis = -1
	// 2. If there is no variation between individuals, Fst = 0
	// 3. Weir & Cockerham 1984 estimator > Weir & Hill 2002
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
		const auto genotypeAB = createTestNeutralGenotype(genomeSz, true, alleleValPopA, alleleValPopB);
		// all individuals have genotype AB

		// Create two-patches landscape
		Landscape* pLandscape = new Landscape;
		vector<Patch*> patches(nbPatches);
		vector<Cell*> cells(nbPatches);
		set<int> patchList;
		for (int i = 0; i < nbPatches; i++) {
			patches[i] = pLandscape->newPatch(i);
			cells[i] = new Cell(0, 0, patches[i], 0);
			patches[i]->addCell(cells[i], 0, 0);
			patchList.insert(patches[i]->getPatchNum());
		}
		const string indSampling = "all";
		const set<int> stgToSample = { 1 };

		// Create species trait structure
		Species* pSpecies = createDefaultSpecies();
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
				Individual* pInd = new Individual(pSpecies, cells[p], patches[p], indStg, 0, 0, 0.0, false, 1);
				pInd->setUpGenes(pSpecies, 1.0);
				pInd->overrideGenotype(NEUTRAL, genotypeAB);
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
			nbIndsPerPop* patchList.size(),
			nbLoci,
			maxNbNeutralAlleles,
			pSpecies,
			pLandscape,
			false
		);
		assert(pNeutralStatistics->getFstWC() == 0.0);
		assert(pNeutralStatistics->getFisWC() == -1.0);

	//	pNeutralStatistics->calcPairwiseWeightedFst(
	//		patchList,
	//		nbIndsPerPop* patchList.size(),
	//		nbLoci,
	//		pSpecies,
	//		pLandscape
	//	);
	//	assert(pNeutralStatistics->getWeightedFst() < pNeutralStatistics->getFstWC());
		// Weir and Hill is still equal to Weir and Cockerham full homozygote case
	//	const double tol = 0.000001;
	//	assert(abs(pNeutralStatistics->getWeightedFst() - refWeirCockerhamDiploidFst) < tol);
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
			cells[i] = new Cell(0, 0, patches[i], 0);
			patches[i]->addCell(cells[i], 0, 0);
			patchList.insert(patches[i]->getPatchNum());
		}
		const string indSampling = "all";
		const set<int> stgToSample = { 1 };

		// Create species trait structure
		Species* pSpecies = createDefaultSpecies();
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
				Individual* pInd = new Individual(pSpecies, cells[p], patches[p], indStg, 0, 0, 0.0, false, 1);
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
			pLandscape,
			false
		);
		const double expectedFst = 0.0583; // calculated by hand from Weir and Cockerham 1984
		double calcError = abs(pNeutralStatistics->getFstWC() - expectedFst);
		assert(calcError < 0.0001);
	}

	// The Fst is a haploid sample equals that of a diploid sample if:
	// 1 - allelic frequencies are equal
	// 2 - the diploid sample has no heterozygotes
	// (reference diploid Fst has already been calculated above)
	{
		const bool isDiploid{ false };
		// Patch setup
		const int nbPatches = 2;
		const int nbIndsPerPop = 10;
		// Genetic setup
		const int genomeSz = 1;
		const set<int> genePositions = { 0 };
		const float maxAlleleVal = 1;
		unsigned char alleleValPopA = char(0);
		unsigned char alleleValPopB = char(1);
		const auto genotypeA = createTestNeutralGenotype(genomeSz, true, alleleValPopA);
		const auto genotypeB = createTestNeutralGenotype(genomeSz, true, alleleValPopB);
		vector<map<int, vector<unsigned char>>> genotypes = {
			// 8 A, 2B (no heterozygotes)
			genotypeA, genotypeA, genotypeA, genotypeA, genotypeB,
			genotypeA, genotypeA, genotypeA, genotypeA, genotypeB
		};

		// Create two-patches landscape
		Landscape* pLandscape = new Landscape;
		vector<Patch*> patches(nbPatches);
		vector<Cell*> cells(nbPatches);
		set<int> patchList;
		for (int i = 0; i < nbPatches; i++) {
			patches[i] = pLandscape->newPatch(i);
			cells[i] = new Cell(0, 0, patches[i], 0);
			patches[i]->addCell(cells[i], 0, 0);
			patchList.insert(patches[i]->getPatchNum());
		}
		const string indSampling = "all";
		const set<int> stgToSample = { 1 };

		// Create species trait structure
		Species* pSpecies = createDefaultSpecies();
		pSpecies->setDemogr(createDefaultHaploidDemogrParams());
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
				Individual* pInd = new Individual(pSpecies, cells[p], patches[p], indStg, 0, 0, 0.0, false, 1);
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
			pLandscape,
			false
		);
		assert(pNeutralStatistics->getFstWC() == refWeirCockerhamDiploidFst);
	}

	// Multi-locus case
	// First locus has no variation between individuals, all heterozygotes (Fst = 0)
	// Second locus has different fixed alleles between both populations (Fst = 1)
	{
		// Patch setup
		const int nbPatches = 2;
		const int nbIndsPerPop = 10;
		// Genetic setup
		const int genomeSz = 3;
		const bool isDiploid{ true };
		const set<int> genePositions = { 0, 2 }; // position 1 is arbitrarily empty
		const float maxAlleleVal = 1;
		unsigned char alleleValPopA = char(0);
		unsigned char alleleValPopB = char(1);
		map<int, vector<unsigned char>> genotypePop1;
		map<int, vector<unsigned char>> genotypePop2;
		// First locus: all heterozygotes, same genotype
		genotypePop1.emplace(0, vector<unsigned char>{alleleValPopA, alleleValPopB});
		genotypePop2.emplace(0, vector<unsigned char>{alleleValPopA, alleleValPopB});
		// Second locus: different fixed alleles
		genotypePop1.emplace(2, vector<unsigned char>{alleleValPopA, alleleValPopA});
		genotypePop2.emplace(2, vector<unsigned char>{alleleValPopB, alleleValPopB});

		// Create two-patches landscape
		Landscape* pLandscape = new Landscape;
		vector<Patch*> patches(nbPatches);
		vector<Cell*> cells(nbPatches);
		set<int> patchList;
		for (int i = 0; i < nbPatches; i++) {
			patches[i] = pLandscape->newPatch(i);
			cells[i] = new Cell(0, 0, patches[i], 0);
			patches[i]->addCell(cells[i], 0, 0);
			patchList.insert(patches[i]->getPatchNum());
		}
		const string indSampling = "all";
		const set<int> stgToSample = { 1 };

		// Create species trait structure
		Species* pSpecies = createDefaultSpecies();
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
				Individual* pInd = new Individual(pSpecies, cells[p], patches[p], indStg, 0, 0, 0.0, false, 1);
				pInd->setUpGenes(pSpecies, 1.0);
				pInd->overrideGenotype(NEUTRAL, p == 0 ? genotypePop1 : genotypePop2);
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
			pLandscape,
			false
		);
		assert(pNeutralStatistics->getPerLocusFst(0) == 0.0);
		assert(pNeutralStatistics->getPerLocusFst(1) == 1.0);
	}

	testPerPatchFis();
}

#endif // UNIT_TESTS
