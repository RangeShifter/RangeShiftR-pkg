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
}

#endif //RSDEBUG