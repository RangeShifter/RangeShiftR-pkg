#if RSDEBUG

#include "../Individual.h"
#include "../Population.h"

void testPopulation()
{
	// Given a genetic load trait, offspring 
	// Survival is (inversely) proportional to the mutation rate
	{
		vector<float> mutationRates = { 0.0, 0.05, 0.1 };
		vector<int> survivingInds;
		const int initialNbInds = 1000;
		const float localK = 10000; // not limiting

		// Simple genetic layout
		const bool isDiploid{ false }; // haploid suffices
		const int genomeSz = 1;
		const set<int> genePositions = { 0 };

		// All mutations are dominant lethal
		const map<GenParamType, float> mutParams{
				pair<GenParamType, float>{GenParamType::MIN, 1},
				pair<GenParamType, float>{GenParamType::MAX, 1}
		};
		const map<GenParamType, float> domParams{
			pair<GenParamType, float>{GenParamType::MIN, 1},
			pair<GenParamType, float>{GenParamType::MAX, 1}
		};

		for (float mutationRate : mutationRates) {
			Landscape* pLandscape = new Landscape;
			Patch* pPatch = pLandscape->newPatch(1);
			Cell* pCell = new Cell(0, 0, (intptr)pPatch, 0);
			pPatch->addCell(pCell, 0, 0);

			Species* pSpecies = new Species();
			pSpecies->setDemogr(createDefaultHaploidDemogrParams());
			pSpecies->setGeneticParameters(
				set<int>{genomeSz - 1}, // single chromosome
				genomeSz,
				0.0, // no recombination
				{ }, "none", { }, 0 // no sampling
			);
			SpeciesTrait* spTr = new SpeciesTrait(
				TraitType::GENETIC_LOAD,
				sex_t::NA,
				genePositions,
				ExpressionType::MULTIPLICATIVE,
				DistributionType::NONE, map<GenParamType, float>{},
				DistributionType::UNIFORM, domParams,
				true, // isInherited
				mutationRate, // mutation rate
				DistributionType::UNIFORM, mutParams,
				isDiploid ? 2 : 1,
				false
			);
			pSpecies->addTrait(TraitType::GENETIC_LOAD, *spTr);

			Population pop = Population(pSpecies, pPatch, initialNbInds, 1);
			pop.reproduction(localK, 1, 1); // juveniles are checked for viability at birth
			pop.fledge(); // non-overlapping: adults are replaced with juveniles
			survivingInds.push_back(pop.getNInds());
		}
		assert(survivingInds[0] > survivingInds[1] && survivingInds[1] > survivingInds[2]);
	}

	// Dispersal is proportional to the mutation rate
	{
		vector<float> mutationRates = { 0.0, 0.05, 0.1 };
		vector<int> emigratingInds;
		const int initialNbInds = 1000;
		const float localK = 10000; // not limiting

		// Simple genetic layout
		const bool isDiploid{ false }; // haploid suffices
		const int genomeSz = 1;
		const set<int> genePositions = { 0 };

		// Wild-types nver emigrate, mutants always do
		const map<GenParamType, float> initParams{
			pair<GenParamType, float>{GenParamType::MIN, 0},
			pair<GenParamType, float>{GenParamType::MAX, 0}
		};
		const map<GenParamType, float> mutParams{
			pair<GenParamType, float>{GenParamType::MIN, 1},
			pair<GenParamType, float>{GenParamType::MAX, 1}
		};

		for (float mutationRate : mutationRates) {
			Landscape* pLandscape = new Landscape;
			Patch* pPatch = pLandscape->newPatch(1);
			Cell* pCell = new Cell(0, 0, (intptr)pPatch, 0);
			pPatch->addCell(pCell, 0, 0);

			Species* pSpecies = new Species();
			emigRules emig;
			emig.indVar = true;
			emig.sexDep = false;
			emig.densDep = false;
			emig.stgDep = false;
			pSpecies->setEmigRules(emig);
			stageParams stg;
			stg.nStages = 2;
			pSpecies->setStage(stg);
			pSpecies->setDemogr(createDefaultHaploidDemogrParams());
			pSpecies->setGeneticParameters(
				set<int>{genomeSz - 1}, // single chromosome
				genomeSz,
				0.0, // no recombination
				{ }, "none", { }, 0 // no sampling
			);
			SpeciesTrait* spTr = new SpeciesTrait(
				TraitType::E_D0,
				sex_t::NA,
				genePositions,
				ExpressionType::ADDITIVE,
				DistributionType::UNIFORM, initParams,
				DistributionType::NONE, map<GenParamType, float>{}, // no dominance
				true, // isInherited
				mutationRate, // mutation rate
				DistributionType::UNIFORM, mutParams,
				isDiploid ? 2 : 1,
				false
			);
			pSpecies->addTrait(TraitType::E_D0, *spTr);

			Population pop = Population(pSpecies, pPatch, initialNbInds, 1);
			pop.reproduction(localK, 1, 1);
			pop.fledge(); // replace initial pop with juveniles
			pop.emigration(localK); // select and flag emigrants
			int popSize = pop.totalPop();
			for (int i = 0; i < popSize; i++) {
				pop.extractDisperser(i); // rm emigrants from pop
			}
			int nbEmigrating = popSize - pop.totalPop(); // diff is nb of emigrants
			if (mutationRate == 0.0)
				assert(nbEmigrating == 0);
			emigratingInds.push_back(nbEmigrating);
		}
		assert(emigratingInds[0] < emigratingInds[1] && emigratingInds[1] < emigratingInds[2]);
	}

	// In the absence of evolutionary forces, neutral gene 
	// frequencies conform to Hardy-Weinberg principle, i.e.:
	// 1 - Allele frequencies p and q remain constant through generations
	// 2 - Genotype frequencies conform to fAA = p^2, fAB = 2pq, fBB = q^2
	{
		float mutationRate = 0.0;
		const float localK = 10000.0;
		const int initialNbInds = localK;
		const float initFreqA = 0.7;
		const float exptdFreqA = initFreqA; // Allelic freqs are constant under HW
		const float exptdFreqB = 1 - exptdFreqA;
		const float exptdFreqHeteroZ = 2 * exptdFreqA * exptdFreqB; // according to HW
		const int nbGens = 10;
		float obsFreqA = 0.0;
		float obsFreqB = 0.0;
		float obsFreqHeteroZ = 0.0;
		const float tolerance = 0.02; // fairly high tolerance, I expect a bit of drift to act.

		// Simple genetic layout
		const bool isDiploid{ true }; // HW only applies to diploids
		const int genomeSz = 1;
		const set<int> genePositions = { 0 };
		const float maxAlleleVal = 1;
		unsigned char alleleA = char(0);
		unsigned char alleleB = char(1);
		auto genotypeAA = createTestNeutralGenotype(genomeSz, true, alleleA, alleleA);
		auto genotypeBB = createTestNeutralGenotype(genomeSz, true, alleleB, alleleB);

		Landscape* pLandscape = new Landscape;
		Patch* pPatch = pLandscape->newPatch(1);
		Cell* pCell = new Cell(0, 0, (intptr)pPatch, 0);
		pPatch->addCell(pCell, 0, 0);

		Species* pSpecies = new Species();
		pSpecies->setDemogr(createDefaultDiploidDemogrParams());
		pSpecies->setGeneticParameters(
			set<int>{genomeSz - 1}, // single chromosome
			genomeSz,
			0.0, // no recombination
			{ }, "none", { }, 0 // no sampling
		);
		SpeciesTrait* spTr = createTestNeutralSpTrait(maxAlleleVal, genePositions, isDiploid);
		pSpecies->addTrait(TraitType::NEUTRAL, *spTr);

		// Initialise population with 
		Population pop = Population(pSpecies, pPatch, 0, 1);
		for (int i = 0; i < initialNbInds; i++) {
			Individual* pInd = new Individual(pCell, pPatch, 1, 0, 0, 0.5, false, 1);
			pInd->setUpGenes(pSpecies, 1.0);
			if (i < initialNbInds * initFreqA)
				pInd->overrideGenotype(NEUTRAL, genotypeAA);
			else 
				pInd->overrideGenotype(NEUTRAL, genotypeBB);
			pop.recruit(pInd);
		}

		// Check allele frequencies conform to HW through generations
		for (int yr = 0; yr < nbGens; yr++) {
			pop.reproduction(localK, 1, 1);
			pop.fledge(); // replace initial pop with juveniles
			pop.survival0(localK, 0, 0); // flag juveniles for development
			pop.survival1(); // develop to stage 1 (breeders)

			// Count allele and heterozygote frequencies
			pop.sampleIndsWithoutReplacement("all", { 1 });
			pop.updatePopNeutralTables();
			obsFreqA = pop.getAlleleFrequency(0, alleleA);
			obsFreqB = pop.getAlleleFrequency(0, alleleB);
			float nbHeteroZ = pop.getHeteroTally(0, alleleA);
			int nbInds = pop.getNInds();
			obsFreqHeteroZ = nbHeteroZ / nbInds;
			assert(abs(obsFreqA - exptdFreqA) < tolerance);
			assert(abs(obsFreqB - exptdFreqB) < tolerance);
			assert(abs(obsFreqHeteroZ - exptdFreqHeteroZ) < tolerance);
		}
	}

	// Genetic load meets Hardy-Weinberg expectation
	// If a lethal (s = 1) recessive (h = 0) allele starts at freq 0.6,
	// then (if no mutations) the prop. of unviable homozygote offspring should be 0.36
	{
		const float initFreqA = 0.6;
		const float sA = 1.0; // lethal
		const float hA = 0.0; // fully recessive
		const float sB = 0.0; // benign
		const float hB = 1.0; // fully dominant
		float mutationRate = 0.0;
		const float localK = 10000.0;
		const int initialNbInds = localK;
		const float tolerance = 0.015; // high tolerance, still a lot of stochasticity
		const float expectedFreqAA = initFreqA * initFreqA;

		// Simple genetic layout
		const bool isDiploid{ true }; // HW only applies to diploids
		const int genomeSz = 1;
		const set<int> genePositions = { 0 };
		const float maxAlleleVal = 1;
		unsigned char alleleA = char(0);
		unsigned char alleleB = char(1);
		auto genotypeAA = createTestGenotype(genomeSz, true, sA, sA, hA, hA);
		auto genotypeBB = createTestGenotype(genomeSz, true, sB, sB, hB, hB);

		Landscape* pLandscape = new Landscape;
		Patch* pPatch = pLandscape->newPatch(1);
		Cell* pCell = new Cell(0, 0, (intptr)pPatch, 0);
		pPatch->addCell(pCell, 0, 0);

		Species* pSpecies = new Species();
		pSpecies->setDemogr(createDefaultDiploidDemogrParams());
		pSpecies->setGeneticParameters(
			set<int>{genomeSz - 1}, // single chromosome
			genomeSz,
			0.0, // no recombination
			{ }, "none", { }, 0 // no sampling
		);
		SpeciesTrait* spTr = createTestGenLoadTrait(genePositions, isDiploid);
		pSpecies->addTrait(TraitType::GENETIC_LOAD, *spTr);

		// Initialise population with 
		Population pop = Population(pSpecies, pPatch, 0, 1);
		for (int i = 0; i < initialNbInds; i++) {
			Individual* pInd = new Individual(pCell, pPatch, 1, 0, 0, 0.5, false, 1);
			pInd->setUpGenes(pSpecies, 1.0);
			if (i < initialNbInds * initFreqA)
				pInd->overrideGenotype(GENETIC_LOAD1, genotypeAA);
			else
				pInd->overrideGenotype(GENETIC_LOAD1, genotypeBB);
			pop.recruit(pInd);
		}

		// Check allele frequencies conform to HW
		pop.reproduction(localK, 1, 1);
		pop.fledge(); // replace initial pop with juveniles
		double obsFreqUnviable = 1 - pop.getNInds() / localK;
		assert(abs(obsFreqUnviable - expectedFreqAA) < tolerance);
	}
}

#endif //RSDEBUG