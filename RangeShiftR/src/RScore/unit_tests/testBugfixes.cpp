#ifdef UNIT_TESTS

// Regression tests for issues found in the code review of October 2026.
// Each test documents the expected (correct) behaviour. On the unfixed code
// base the tests either fail an assertion or trigger an AddressSanitizer /
// UndefinedBehaviourSanitizer report (build with -fsanitize=address,undefined).
//
// Tests are registered in a map so that they can be run all at once
// (testBugfixes(), called from Main.cpp) or individually
// (runBugfixTest("name"), useful to demonstrate one failure at a time).

#include <numeric>
#include <functional>
#include <map>
#include <string>
#include <cmath>

#include "../Individual.h"
#include "../Population.h"
#include "../Community.h"
#include "../SubCommunity.h"
#include "../FractalGenerator.h"

using namespace std;

namespace {

	// ---------------------------------------------------------------------
	// helpers
	// ---------------------------------------------------------------------

	// Kernel-dispersing, non-structured asexual species with one suitable
	// habitat (index 0, K = 100) and one unsuitable habitat (index 1, K = 0)
	void setUpKernelSpecies(Species& sp, float meanDist, bool wait, bool go2nbr) {
		sp.createHabK(2);
		sp.setHabK(0, 100.0);
		sp.setHabK(1, 0.0);
		demogrParams d = createDefaultHaploidDemogrParams();
		sp.setDemogr(d);
		transferRules trfr;
		trfr.usesMovtProc = false;
		trfr.indVar = trfr.sexDep = trfr.stgDep = false;
		trfr.twinKern = trfr.distMort = trfr.habMort = trfr.costMap = false;
		trfr.moveType = 1;
		sp.setTrfrRules(trfr);
		sp.setFullKernel(false);
		trfrMortParams mort;
		mort.fixedMort = 0.0; mort.mortAlpha = 0.0; mort.mortBeta = 1.0;
		sp.setMortParams(mort);
		settleType st; st.stgDep = st.sexDep = st.indVar = false;
		sp.setSettle(st);
		settleRules sr;
		sr.densDep = false; sr.findMate = false;
		sr.wait = wait; sr.go2nbrLocn = go2nbr;
		sp.setSettRules(0, 0, sr);
		trfrKernelParams kern;
		kern.meanDist1 = meanDist; kern.meanDist2 = meanDist; kern.probKern1 = 1.0;
		sp.overrideKernels(0, 0, kern);
	}

	landParams makeLandParams(int dimX, int dimY) {
		landParams lp = createDefaultLandParams(dimX);
		lp.dimX = dimX; lp.dimY = dimY;
		lp.maxX = dimX - 1; lp.maxY = dimY - 1;
		return lp;
	}

	// ---------------------------------------------------------------------
	// RSrandom
	// ---------------------------------------------------------------------

	// getRNG() must hand out the generator itself, not a copy: otherwise
	// two consecutive shuffles (or samples) reuse the same random stream.
	void test_getRNG_advances_state() {
		RSrandom rsr;
		vector<int> a(50), b(50);
		iota(a.begin(), a.end(), 0);
		iota(b.begin(), b.end(), 0);
		shuffle(a.begin(), a.end(), rsr.getRNG());
		shuffle(b.begin(), b.end(), rsr.getRNG());
		assert(a != b); // P(equal | correct RNG) = 1/50!
	}

	// ---------------------------------------------------------------------
	// Species
	// ---------------------------------------------------------------------

	// The ploidy flag must be set by the constructor, not only by setDemogr()
	void test_species_diploid_initialised() {
		Species asexual(0);
		assert(!asexual.isDiploid());
		Species sexual(1);
		assert(sexual.isDiploid());
	}

	// 10 stages x 2 sexes = 20 is a valid (maximal) weights matrix dimension
	void test_ddweights_maximal_dimension() {
		Species sp;
		const short dim = gMaxNbStages * gMaxNbSexes;
		sp.createDDwtFec(dim);
		sp.setDDwtFec(dim - 1, dim - 1, 0.5f);
		assert(sp.getDDwtFec(dim - 1, dim - 1) == 0.5f);
	}

	// Requesting a trait that was never added must not dereference end()
	void test_getSpTrait_missing_trait_throws() {
		Species sp;
		bool thrown = false;
		try { sp.getSpTrait(NEUTRAL); }
		catch (const exception&) { thrown = true; }
		assert(thrown);
	}

	// ---------------------------------------------------------------------
	// Parameters
	// ---------------------------------------------------------------------

	// saveVisits (SMS heat map) must survive a setSim()/getSim() round trip
	// also outside of the R package
	void test_paramSim_roundtrip_saveVisits() {
		paramSim p;
		simParams s = p.getSim();
		s.saveVisits = true;
		s.batchNum = 0; s.simulation = 0; s.reps = 1; s.years = 1;
		p.setSim(s);
		assert(p.getSim().saveVisits == true);
	}

	// ---------------------------------------------------------------------
	// Cell / Patch
	// ---------------------------------------------------------------------

	// An invalid landscape-change index must not write out of bounds
	void test_cell_changeHabIndex_out_of_range() {
		Cell c(0, 0, (Patch*)nullptr, 1);
		c.changeHabIndex(5, 2); // only index 0 exists
		assert(c.getHabIndex(0) == 1);
	}

	// Cells without a demographic scaling vector for the current landscape
	// must not be read out of bounds
	void test_patch_demoscaling_missing_layer() {
		Patch p(1, 1);
		Cell* c = new Cell(0, 0, &p, 0);
		p.addCell(c, 0, 0);
		patchLimits lim; lim.xMin = lim.yMin = 0; lim.xMax = lim.yMax = 10;
		p.setPatchDemoScaling(0, lim); // the cell has no scaling layers at all
		vector<float> ds = p.getDemoScaling();
		assert((int)ds.size() == nDSlayer);
		for (float v : ds) assert(v >= 0.0f && v <= 1.0f);
		delete c;
	}

	// ---------------------------------------------------------------------
	// Population
	// ---------------------------------------------------------------------

	// Translocation filter: max_age must remove individuals OLDER than max_age
	void test_translocation_max_age_filter() {
		Landscape* pLandscape = new Landscape;
		Patch* pPatch = pLandscape->newPatch(1);
		Cell* pCell = new Cell(0, 0, pPatch, 0);
		pPatch->addCell(pCell, 0, 0);
		Species* pSpecies = createDefaultSpecies();
		pSpecies->setDemogr(createDefaultHaploidDemogrParams());
		Population pop(pSpecies, pPatch, 0, 1);
		for (short age = 1; age <= 5; age++)
			pop.recruit(new Individual(pSpecies, pCell, pPatch, 1, age, 0, 0.0, false, 1));

		vector<Individual*> selected = pop.getIndsWithCharacteristics(-9, 3, -9, -9);
		assert(selected.size() == 3);
		for (Individual* ind : selected) assert(ind->getStats().age <= 3);
	}

	// Offspring of an asexual (female-only) species must all be female,
	// even if a proportion of males has been set
	void test_asexual_offspring_are_female() {
		Landscape* pLandscape = new Landscape;
		Patch* pPatch = pLandscape->newPatch(1);
		Cell* pCell = new Cell(0, 0, pPatch, 0);
		pPatch->addCell(pCell, 0, 0);
		Species* pSpecies = createDefaultSpecies();
		demogrParams d = createDefaultHaploidDemogrParams();
		d.propMales = 0.5; // e.g. left over from a sexual parameterisation
		pSpecies->setDemogr(d);
		Population pop(pSpecies, pPatch, 200, 1);
		vector<float> scaling(nDSlayer, 1.0f);
		pop.reproduction(100000.0f, 1.0f, 1, scaling);
		pop.fledge();
		assert(pop.getNbInds() > 0);
		for (Individual* ind : pop.getIndividualsInStage(0))
			assert(ind->getSex() == FEM);
		assert(pop.getNbInds(0, 1) == 0); // no males counted
	}

	// ---------------------------------------------------------------------
	// Neutral statistics
	// ---------------------------------------------------------------------

	struct TwoPatchNeutralSetup {
		Landscape* pLandscape;
		Species* pSpecies;
		vector<Patch*> patches;
		set<int> patchList;
		int nbInds = 0;
	};

	// one locus, two patches, genotypes given per individual
	TwoPatchNeutralSetup makeNeutralSetup(const vector<vector<pair<int, int>>>& genotypes) {
		TwoPatchNeutralSetup s;
		const int genomeSz = 1;
		s.pLandscape = new Landscape;
		for (int i = 0; i < (int)genotypes.size(); i++) {
			Patch* p = s.pLandscape->newPatch(i + 1);
			Cell* c = new Cell(i, 0, p, 0);
			p->addCell(c, i, 0);
			s.patches.push_back(p);
			s.patchList.insert(p->getPatchNum());
		}
		s.pSpecies = createDefaultSpecies();
		s.pSpecies->setDemogr(createDefaultDiploidDemogrParams());
		s.pSpecies->setGeneticParameters(set<int>{genomeSz - 1}, genomeSz, 0.0,
			s.patchList, "all", { 1 }, 1);
		SpeciesTrait* spTr = createTestNeutralSpTrait(1, { 0 }, true);
		s.pSpecies->addTrait(TraitType::NEUTRAL, *spTr);

		for (int p = 0; p < (int)genotypes.size(); p++) {
			Population* pPop = new Population(s.pSpecies, s.patches[p], 0, 1);
			for (auto [a, b] : genotypes[p]) {
				Individual* pInd = new Individual(s.pSpecies, s.patches[p]->getCell(0), s.patches[p], 1, 0, 0, 0.0, false, 1);
				pInd->setUpGenes(s.pSpecies, 1);
				pInd->overrideGenotype(NEUTRAL, createTestNeutralGenotype(genomeSz, true, (unsigned char)a, (unsigned char)b));
				pPop->recruit(pInd);
				s.nbInds++;
			}
			pPop->sampleIndsWithoutReplacement("all", { 1 });
		}
		return s;
	}

	// With a single locus, the per-locus F-statistics must equal the global ones
	void test_perLocusFit_equals_global_for_one_locus() {
		auto s = makeNeutralSetup({ { {0,1}, {0,0}, {0,1} }, { {1,1}, {0,1}, {1,1} } });
		NeutralStatsManager nsm(s.patchList.size(), 1);
		nsm.updateAllNeutralTables(s.pSpecies, s.pLandscape, s.patchList);
		nsm.calculateFstatWC(s.patchList, s.nbInds, 1, 2, s.pSpecies, s.pLandscape, false);
		assert(fabs(nsm.getPerLocusFst(0) - nsm.getFstWC()) < 1e-9);
		assert(fabs(nsm.getPerLocusFis(0) - nsm.getFisWC()) < 1e-9);
		assert(fabs(nsm.getPerLocusFit(0) - nsm.getFitWC()) < 1e-9);
	}

	// A patch in which every sampled individual is heterozygous has Ho = 1
	void test_patch_heterozygosity_not_doubled() {
		auto s = makeNeutralSetup({ { {0,1}, {0,1}, {1,0} } });
		for (Patch* p : s.patches) p->getPopn(s.pSpecies)->updatePopNeutralTables();
		Community comm(s.pLandscape);
		const float het = comm.getPatchHet(s.pSpecies, s.patches[0]->getPatchNum(), 0);
		assert(fabs(het - 1.0f) < 1e-6);
	}

	// All loci monomorphic: F-statistics are undefined and must not be NaN
	void test_fstat_monomorphic_no_nan() {
		auto s = makeNeutralSetup({ { {0,0}, {0,0} }, { {0,0}, {0,0} } });
		NeutralStatsManager nsm(s.patchList.size(), 1);
		nsm.updateAllNeutralTables(s.pSpecies, s.pLandscape, s.patchList);
		nsm.calculateFstatWC(s.patchList, s.nbInds, 1, 2, s.pSpecies, s.pLandscape, false);
		assert(std::isnan(nsm.getFstWC()));
		assert(std::isnan(nsm.getFitWC()));
		assert(std::isnan(nsm.getPerLocusFst(0)));
		assert(std::isnan(nsm.getPerLocusFit(0)));
	}

	// ---------------------------------------------------------------------
	// Genetic traits
	// ---------------------------------------------------------------------

	// Ecological plausibility: a heterozygote carrying one deleterious allele
	// (s, h) and one wild-type allele must have fitness 1 - h*s.
	// With the relative-dominance formula h_A / (h_A + h_B) and h_wildtype = 0,
	// every mutation becomes fully dominant (fitness 1 - s) irrespective of h.
	void test_genetic_load_heterozygote_with_wildtype() {
		const int genomeSz = 1;
		Landscape ls;
		Patch* pPatch = ls.newPatch(1);
		Cell* pCell = new Cell(0, 0, pPatch, 0);
		pPatch->addCell(pCell, 0, 0);
		Species* pSpecies = createDefaultSpecies();
		pSpecies->setDemogr(createDefaultDiploidDemogrParams());
		pSpecies->setGeneticParameters(set<int>{genomeSz - 1}, genomeSz, 0.0, set<int>{}, "none", set<int>{}, 0);
		const map<GenParamType, float> distParams{ {MIN, 0.0}, {MAX, 0.0} };
		SpeciesTrait* spTr = new SpeciesTrait(TraitType::GENETIC_LOAD1, sex_t::NA, { 0 },
			ExpressionType::MULTIPLICATIVE, { 0 },
			DistributionType::NONE, distParams, DistributionType::UNIFORM, distParams,
			true, 0.0, DistributionType::UNIFORM, distParams, DistributionType::UNIFORM, distParams,
			2, false);
		pSpecies->addTrait(TraitType::GENETIC_LOAD1, *spTr);
		Individual ind(pSpecies, pCell, pPatch, 0, 0, 0, 0.0, false, 0);
		ind.setUpGenes(pSpecies, 1);

		const float s = 0.5f, h = 0.2f; // partially recessive deleterious allele
		ind.overrideGenotype(GENETIC_LOAD1, createTestGenotype(genomeSz, true, s, 0.0f, h, 0.0f));
		const float w = ind.getTrait(GENETIC_LOAD1)->express();
		assert(fabs(w - (1.0f - h * s)) < 1e-6); // 0.9, not 0.5
	}

	// Ecological plausibility: evolving settlement slopes must be allowed to
	// become negative. With p_s = S0 / (1 + exp(-(N/K - beta) * alpha)),
	// settlement only decreases with density if alpha < 0 (eqn. 16 of the manual).
	void test_settlement_alpha_can_be_negative() {
		const map<GenParamType, float> init{ {MIN, -2.0}, {MAX, -2.0} };
		const map<GenParamType, float> mut{ {MIN, 0.0}, {MAX, 0.0} };
		SpeciesTrait spTr(TraitType::S_ALPHA, sex_t::NA, { 0 }, ExpressionType::ADDITIVE, { 0 },
			DistributionType::UNIFORM, init, DistributionType::NONE, {},
			true, 0.0, DistributionType::UNIFORM, mut, DistributionType::NONE, {}, 1, false);
		DispersalTrait trait(&spTr);
		assert(trait.express() < 0.0f);
	}

	// A neutral trait without the mandatory MAX parameter must raise a clear
	// error instead of dereferencing map::end()
	void test_neutral_trait_missing_max_throws() {
		SpeciesTrait spTr(TraitType::NEUTRAL, sex_t::NA, { 0 }, ExpressionType::NOTEXPR, { 0 },
			DistributionType::UNIFORM, { {MAX, 3.0} }, DistributionType::NONE, {},
			true, 0.0, DistributionType::KAM, {} /* no MAX */, DistributionType::NONE, {}, 2, false);
		bool thrown = false;
		try { NeutralTrait t(&spTr); }
		catch (const logic_error&) { thrown = true; }
		assert(thrown);
	}

	// ---------------------------------------------------------------------
	// Dispersal
	// ---------------------------------------------------------------------

	// Ecological plausibility: an individual that waited in the matrix (status 3)
	// must resume kernel dispersal from its CURRENT location (manual 2.5.3),
	// not be put back into its natal patch.
	void test_kernel_waiting_individual_restarts_from_current_cell() {
		const int dimX = 31;
		landParams lp = makeLandParams(dimX, 1);
		Landscape ls;
		ls.setLandParams(lp, true);
		ls.setCellArray();
		Cell* natal = new Cell(0, 0, (Patch*)nullptr, 0); // suitable
		ls.addCellToLand(natal);
		vector<Cell*> row{ natal };
		for (int x = 1; x < dimX; x++) {
			Cell* c = new Cell(x, 0, (Patch*)nullptr, 1); // matrix
			ls.addCellToLand(c);
			row.push_back(c);
		}
		Species sp;
		setUpKernelSpecies(sp, 1.0f, true, false);
		demogrParams d = sp.getDemogrParams(); d.stageStruct = true; sp.setDemogr(d);
		stageParams stg = sp.getStageParams(); stg.nStages = 2; sp.setStage(stg);
		ls.allocatePatches(&sp);
		ls.updateCarryingCapacity(&sp, 0, 0);

		int nearNatal = 0;
		const int nReps = 50;
		for (int r = 0; r < nReps; r++) {
			Individual ind(&sp, natal, natal->getPatch(), 1, 1, 0, 0.0, false, 0);
			for (int x = 1; x <= 20; x++) ind.moveto(row[x]); // walk to x = 20
			ind.setStatus(3);         // waiting in the matrix ...
			ind.ageIncrement(1000);   // ... and released for the next dispersal event
			assert(ind.getStatus() == 1);
			ind.moveKernel(&ls, &sp, false);
			if (ind.getCurrCell() != nullptr && ind.getCurrCell()->getLocn().x < 10) nearNatal++;
		}
		// mean dispersal distance is one cell: from x = 20, reaching x < 10 is (nearly) impossible
		assert(nearNatal == 0);
	}

	// An individual that crossed an absorbing boundary has no current cell;
	// the "go to neighbouring cell" settlement rule must not dereference it
	void test_settlement_neighbour_rule_absorbed_individual() {
		landParams lp = makeLandParams(3, 3);
		Landscape ls;
		ls.setLandParams(lp, true);
		ls.setCellArray();
		for (int x = 0; x < 3; x++) for (int y = 0; y < 3; y++)
			ls.addCellToLand(new Cell(x, y, (Patch*)nullptr, (x == 0 && y == 0) ? 0 : 1));
		Species sp;
		setUpKernelSpecies(sp, 1000.0f, false, true); // huge kernel -> leaves landscape
		ls.allocatePatches(&sp);
		ls.updateCarryingCapacity(&sp, 0, 0);
		Cell* natal = ls.findCell(0, 0);

		Individual* pInd = nullptr;
		for (int attempt = 0; attempt < 100; attempt++) {
			delete pInd;
			pInd = new Individual(&sp, natal, natal->getPatch(), 1, 0, 0, 0.0, false, 0);
			pInd->setStatus(1);
			pInd->moveKernel(&ls, &sp, true); // absorbing boundaries
			if (pInd->getCurrCell() == nullptr) break;
		}
		assert(pInd->getCurrCell() == nullptr);
		assert(pInd->getStatus() == 6);

		std::map<Species*, vector<Individual*>> pool{ { &sp, { pInd } } };
		SubCommunity::resolveSettlement(pool, &ls); // must not crash
		assert(pInd->getStatus() == 6);
		delete pInd;
	}

	// ---------------------------------------------------------------------
	// Landscape generation
	// ---------------------------------------------------------------------

	// Fractal landscapes with X != Y (allowed by the manual, e.g. 33 x 65)
	// must not read outside the arena (diamond step used ky2 for the x index)
	void test_fractal_non_square() {
		vector<land> a = fractal_landscape(9, 17, 0.5, 0.5, 100.0, 1.0);
		assert((int)a.size() == 9 * 17);
		vector<land> b = fractal_landscape(17, 9, 0.5, 0.5, 100.0, 1.0);
		assert((int)b.size() == 9 * 17);
	}

	// new[] must be paired with delete[] (ASan: alloc-dealloc-mismatch)
	void test_connect_matrix_delete() {
		Landscape ls;
		ls.newPatch(0); ls.newPatch(1); ls.newPatch(2);
		ls.createConnectMatrix();
		ls.incrConnectMatrix(1, 2);
		ls.deleteConnectMatrix();
	}

	const std::map<string, std::function<void()>>& bugfixTests() {
		static const std::map<string, std::function<void()>> tests{
			{ "getRNG_advances_state", test_getRNG_advances_state },
			{ "species_diploid_initialised", test_species_diploid_initialised },
			{ "ddweights_maximal_dimension", test_ddweights_maximal_dimension },
			{ "getSpTrait_missing_trait_throws", test_getSpTrait_missing_trait_throws },
			{ "paramSim_roundtrip_saveVisits", test_paramSim_roundtrip_saveVisits },
			{ "cell_changeHabIndex_out_of_range", test_cell_changeHabIndex_out_of_range },
			{ "patch_demoscaling_missing_layer", test_patch_demoscaling_missing_layer },
			{ "translocation_max_age_filter", test_translocation_max_age_filter },
			{ "asexual_offspring_are_female", test_asexual_offspring_are_female },
			{ "perLocusFit_equals_global_for_one_locus", test_perLocusFit_equals_global_for_one_locus },
			{ "patch_heterozygosity_not_doubled", test_patch_heterozygosity_not_doubled },
			{ "fstat_monomorphic_no_nan", test_fstat_monomorphic_no_nan },
			{ "genetic_load_heterozygote_with_wildtype", test_genetic_load_heterozygote_with_wildtype },
			{ "settlement_alpha_can_be_negative", test_settlement_alpha_can_be_negative },
			{ "neutral_trait_missing_max_throws", test_neutral_trait_missing_max_throws },
			{ "kernel_waiting_individual_restarts_from_current_cell", test_kernel_waiting_individual_restarts_from_current_cell },
			{ "settlement_neighbour_rule_absorbed_individual", test_settlement_neighbour_rule_absorbed_individual },
			{ "fractal_non_square", test_fractal_non_square },
			{ "connect_matrix_delete", test_connect_matrix_delete },
		};
		return tests;
	}
} // namespace

// Run a single regression test by name; returns false if the name is unknown
bool runBugfixTest(const string& name) {
	auto it = bugfixTests().find(name);
	if (it == bugfixTests().end()) return false;
	it->second();
	return true;
}

vector<string> listBugfixTests() {
	vector<string> names;
	for (auto& [name, f] : bugfixTests()) names.push_back(name);
	return names;
}

void testBugfixes() {
	for (auto& [name, f] : bugfixTests()) f();
}

#endif // UNIT_TESTS
