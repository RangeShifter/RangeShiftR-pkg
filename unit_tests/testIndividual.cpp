#ifndef NDEBUG

#include "../Individual.h"
#include "../Population.h"

void testTransferKernels() {
	// Simple 5*5 cell-based landscape layout
	int lsDim = 5;
	landParams ls_params = createDefaultLandParams(lsDim);

	Landscape ls;
	ls.setLandParams(ls_params, true);

	// Two suitable cells in opposite corners
	Cell* init_cell = new Cell(0, 0, 0, 0);
	Cell* final_cell = new Cell(ls_params.dimX - 1, ls_params.dimY - 1, 0, 0);
	ls.setCellArray();
	ls.addCellToLand(init_cell);
	ls.addCellToLand(final_cell);

	// Set up species
	Species sp;
	// Habitat codes
	sp.createHabK(1);
	sp.setHabK(0, 100.0); // one habitat with K = 100
	// Demography
	demogrParams d;
	d.stageStruct = false;
	sp.setDemogr(d);
	// Transfer rules
	transferRules trfr;
	trfr.indVar = trfr.sexDep = trfr.stgDep = false;
	trfr.twinKern = trfr.distMort = false;
	sp.setTrfrRules(trfr);
	sp.setFullKernel(false);
	// Transfer traits
	trfrKernelParams kern;
	kern.meanDist1 = static_cast<float>(ls_params.dimX); // can reach destination cell reasonably often
	sp.setSpKernTraits(0, 0, kern, ls_params.resol);
	// Transfer mortality params
	trfrMortParams mort;
	mort.fixedMort = 0.0;
	sp.setMortParams(mort);
	// Settlement
	settleRules sett;
	sett.wait = false;
	sp.setSettRules(0, 0, sett);

	// Set up patches
	ls.allocatePatches(&sp);
	ls.updateCarryingCapacity(&sp, 0, 0);
	Patch* init_patch = (Patch*)init_cell->getPatch();

	// Create and set up individual
	Individual ind1(&sp,init_cell, init_patch, 1, 0, 0, 0.0, false, 0);
	int isDispersing = ind1.moveKernel(&ls, &sp, false);

	// After moving, individual should be in the only available cell
	Cell* curr_cell = ind1.getCurrCell();
	assert(curr_cell != init_cell);
	assert(curr_cell == final_cell);
	assert(ind1.getStatus() == 2); // potential settler

	// If no cell within reasonable dispersal reach, individual does not move and dies
	kern.meanDist1 = 1.0;
	sp.setSpKernTraits(0, 0, kern, ls_params.resol);
	Individual ind2(&sp, init_cell, init_patch, 1, 0, 0, 0.0, false, 0); // reset individual
	isDispersing = ind2.moveKernel(&ls, &sp, false);
	curr_cell = ind2.getCurrCell();
	assert(ind2.getStatus() == 6); // RIP in peace
	assert(curr_cell == init_cell);

	// Twin kernels
	trfr.twinKern = true;
	sp.setTrfrRules(trfr);
	kern.meanDist1 = 1.0; // very unlikely to reach suitable cell
	kern.meanDist2 = 5.0; // easily reaches suitable cell...
	kern.probKern1 = 1.0; // ... but never used
	sp.setSpKernTraits(0, 0, kern, ls_params.resol);
	Individual ind3(&sp, init_cell, init_patch, 1, 0, 0, 0.0, false, 0);
	isDispersing = ind3.moveKernel(&ls, &sp, false);
	assert(ind3.getStatus() == 6); // dead, could not reach destination cell
	kern.probKern1 = 0.0; // always use second kernel
	sp.setSpKernTraits(0, 0, kern, ls_params.resol);
	Individual ind4(&sp, init_cell, init_patch, 1, 0, 0, 0.0, false, 0);
	isDispersing = ind4.moveKernel(&ls, &sp, false);
	assert(ind4.getStatus() == 2);
	// reset
	trfr.twinKern = false;
	sp.setTrfrRules(trfr);
	kern.probKern1 = 1.0;
	sp.setSpKernTraits(0, 0, kern, ls_params.resol);

	// Sex-dependent dispersal distances
	trfr.sexDep = true;
	sp.setTrfrRules(trfr);
	trfrKernelParams kern_f = kern;
	kern_f.meanDist1 = 1.0; // female very unlikely to reach suitable cell
	sp.setSpKernTraits(0, 0, kern_f, ls_params.resol);
	trfrKernelParams kern_m = kern;
	kern_m.meanDist1 = 5.0; // male easily reaches suitable cell
	sp.setSpKernTraits(0, 1, kern_m, ls_params.resol);

	Individual ind5(&sp, init_cell, init_patch, 1, 0, 0, 0.0, false, 0); // female as default
	isDispersing = ind5.moveKernel(&ls, &sp, false);
	assert(ind5.getStatus() == 6); // dead, could not reach destination

	Individual ind6(&sp, init_cell, init_patch, 1, 0, 0, 1.0, false, 0); // male
	assert(ind6.getSex() == 1);
	isDispersing = ind6.moveKernel(&ls, &sp, false);
	assert(ind6.getStatus() == 2);
	// reset
	trfr.sexDep = false;
	sp.setTrfrRules(trfr);

	// Stage-dependent
	trfr.stgDep = true;
	sp.setTrfrRules(trfr);
	trfrKernelParams kern_juv = kern;
	kern_juv.meanDist1 = 1.0; // juveniles very unlikely to reach suitable cell
	sp.setSpKernTraits(0, 0, kern_juv, ls_params.resol);
	trfrKernelParams kern_adult = kern;
	kern_adult.meanDist1 = 5.0; // adults easily reach suitable cell
	sp.setSpKernTraits(1, 0, kern_adult, ls_params.resol);

	Individual ind7(&sp, init_cell, init_patch, 0, 0, 0, 0.0, false, 0); // juvenile
	isDispersing = ind7.moveKernel(&ls, &sp, false);
	assert(ind7.getStatus() == 6);

	Individual ind8(&sp, init_cell, init_patch, 1, 0, 0, 0.0, false, 0); // adult by default
	isDispersing = ind8.moveKernel(&ls, &sp, false);
	assert(ind8.getStatus() == 2);
	// reset
	trfr.stgDep = false;
	sp.setTrfrRules(trfr);

	/* Boundaries: dispersal distance overshoots
	Only adjacent cells are available
	-----
	-ooo-
	-oio-
	-ooo-
	-----
	*/
	ls.setCellArray(); // reset cells
	vector <Cell*> cells;
	// Set central cell and all adjacent
	for (int x = ls_params.minX + 1; x < ls_params.maxX; ++x) {
		for (int y = ls_params.minY + 1; y < ls_params.maxY; ++y) {
			cells.push_back(new Cell(x, y, 0, 0));
		}
	}

	for (auto c : cells) ls.addCellToLand(c);
	ls.allocatePatches(&sp);
	ls.updateCarryingCapacity(&sp, 0, 0);
	init_cell = cells[4]; // that is, the center
	init_patch = (Patch*)init_cell->getPatch();

	kern.meanDist1 = 10; // overshoots *most* of the time...
	sp.setSpKernTraits(0, 0, kern, ls_params.resol);
	Individual ind9(&sp, init_cell, init_patch, 1, 0, 0, 0.0, false, 0); // reset individual

	// Non-absorbing boundaries
	bool absorbing_boundaries{ false };
	isDispersing = ind9.moveKernel(&ls, &sp, absorbing_boundaries);
	curr_cell = ind9.getCurrCell();
	assert(curr_cell != init_cell); // ...should be able to move eventually
	assert(ind9.getStatus() == 2);

	// Absorbing boundaries
	absorbing_boundaries = true;
	Individual ind10(&sp, init_cell, init_patch, 1, 0, 0, 0.0, false, 0); // reset individual
	isDispersing = ind10.moveKernel(&ls, &sp, absorbing_boundaries);
	curr_cell = ind10.getCurrCell();
	assert(ind10.getStatus() == 6);
	assert(curr_cell == 0); // out of the landscape

	// Dispersal-related mortality
	// Fixed mortality
	mort.fixedMort = 1.0; // Individual *will* die after any step
	sp.setMortParams(mort);
	trfr.distMort = false;
	sp.setTrfrRules(trfr);
	Individual ind11(&sp, init_cell, init_patch, 1, 0, 0, 0.0, false, 0);
	isDispersing = ind11.moveKernel(&ls, &sp, false);
	assert(ind11.getStatus() == 7);
	// Distance-dependent mortality
	trfr.distMort = true;
	sp.setTrfrRules(trfr);
	mort.mortAlpha = 1000.0; // very steep threshold
	mort.mortBeta = 0.5; // very small distance
	sp.setMortParams(mort);
	kern.meanDist1 = 5; // very likely to go over threshold
	sp.setSpKernTraits(0, 0, kern, ls_params.resol);
	Individual ind12(&sp, init_cell, init_patch, 1, 0, 0, 0.0, false, 0);
	isDispersing = ind12.moveKernel(&ls, &sp, false);
	assert(ind12.getStatus() == 7);
	mort.mortBeta = 30; // very large distance, unlikely to draw
	sp.setMortParams(mort);
	Individual ind13(&sp, init_cell, init_patch, 1, 0, 0, 0.0, false, 0);
	isDispersing = ind13.moveKernel(&ls, &sp, false);
	assert(ind13.getStatus() != 7);

	// Reset mortality params
	trfr.distMort = false;
	mort.fixedMort = 0.0;
	sp.setTrfrRules(trfr);
	sp.setMortParams(mort);
}

void testTransferCRW() {
	// Simple 5*5 cell-based landscape layout
	int lsDim = 5;
	landParams ls_params = createDefaultLandParams(lsDim);
	double cellDiagLength = ls_params.resol * SQRT2;

	Landscape ls;
	ls.setLandParams(ls_params, true);

	// All cells are suitable
	const int hab_index = 0;
	vector<Cell*> cell_vec;
	for (int x = ls_params.minX; x < ls_params.dimX; ++x) {
		for (int y = ls_params.minY; y < ls_params.dimY; ++y) {
			cell_vec.push_back(new Cell(x, y, 0, hab_index));
		}
	}
	Cell* init_cell = cell_vec[12]; // central
	ls.setCellArray();
	for (auto c : cell_vec) ls.addCellToLand(c);

	// Set up species
	Species sp;

	// Habitat codes
	sp.createHabK(1);
	sp.setHabK(hab_index, 100.0); // one habitat with K = 100

	// Habitat-dependent mortality
	sp.createHabCostMort(1);

	// Transfer rules
	transferRules trfr;
	trfr.indVar = false;
	trfr.habMort = false;
	trfr.moveType = 2; // CRW
	sp.setTrfrRules(trfr);

	// Transfer CRW traits
	trfrMovtParams m;
	m.stepMort = 0.0;
	m.stepLength = cellDiagLength; // guaranteed to move out
	m.rho = 1.0;
	m.straightenPath = false;
	sp.setSpMovtTraits(m);

	// Settlement rules
	settleRules sett;
	sett.wait = false;
	sp.setSettRules(0, 0, sett);
	settleSteps steps;
	steps.maxSteps = 1;
	steps.minSteps = 1;
	steps.maxStepsYr = 1;
	sp.setSteps(0, 0, steps);

	// Set up patches
	ls.allocatePatches(&sp);
	ls.updateCarryingCapacity(&sp, 0, 0);
	Patch* init_patch = (Patch*)init_cell->getPatch();

	// Create and set up individual
	Individual ind0(&sp, init_cell, init_patch, 1, 0, 0, 0.0, true, 2);

	// Set status
	assert(ind0.getStatus() == 0); // default status, not emigrating
	int isDispersing = ind0.moveStep(&ls, &sp, hab_index, false);
	assert(ind0.getStatus() == 0); // status didn't change
	assert(ind0.getCurrCell() == init_cell); // not emigrating so didn't move

	// Per-step mortality
	m.stepMort = 1.0; // should die
	sp.setSpMovtTraits(m);
	Individual ind1(&sp, init_cell, init_patch, 0, 0, 0, 0.0, true, 2);
	// force set path bc for some reason path gets deallocated upon exiting constructor??
	ind1.setStatus(1);
	isDispersing = ind1.moveStep(&ls, &sp, hab_index, false);
	// Individual begins in natal patch so mortality is disabled
	assert(ind1.getStatus() != 7);
	// Individual should be in a different patch
	Cell* first_step_cell = ind1.getCurrCell();
	assert(first_step_cell != init_cell);

	assert((Patch*)first_step_cell->getPatch() != init_patch);
	ind1.setStatus(1); // emigrating again

	// Individual should die on second step
	isDispersing = ind1.moveStep(&ls, &sp, hab_index, false);
	assert(ind1.getCurrCell() == first_step_cell); // shouldn't have moved
	assert(ind1.getStatus() == 7); // died by transfer
	m.stepMort = 0.0; // not dying
	sp.setSpMovtTraits(m);

	// Habitat-dep mortality
	// ...

	// Step size

	ls = Landscape();
	ls.setLandParams(ls_params, true);
	sp.createHabK(2);
	const int hab_suitable = 0;
	const int hab_unsuitable = 1;
	sp.setHabK(hab_suitable, 100.0);
	sp.setHabK(hab_unsuitable, 0.0);

	// Initial cell unsuitable, suitable cell in opposite corner
	init_cell = new Cell(0, 0, 0, hab_unsuitable);
	Cell* final_cell = new Cell(ls_params.dimX - 1, ls_params.dimY - 1, 0, hab_suitable);
	ls.setCellArray();
	ls.addCellToLand(init_cell);
	ls.addCellToLand(final_cell);
	ls.allocatePatches(&sp);
	ls.updateCarryingCapacity(&sp, 0, 0);
	// Init cell is NOT in natal patch
	Patch* natalPatch = new Patch(0, 0);
	init_patch = (Patch*)init_cell->getPatch();

	// Step length too short
	m.stepLength = 0.1; // will not reach final cell
	m.rho = 0.0; // random angle
	sp.setSpMovtTraits(m);
	steps.minSteps = 1;
	steps.maxStepsYr = 2;
	steps.maxSteps = 3;
	sp.setSteps(0, 0, steps);
	Individual ind2(&sp, init_cell, natalPatch, 0, 0, 0, 0.0, true, 2);
	ind2.setStatus(1); // dispersing
	// First step - still in unsuitable cell so still dispersing
	isDispersing = ind2.moveStep(&ls, &sp, hab_index, false);
	assert(ind2.getCurrCell() == init_cell);
	assert(ind2.getStatus() == 1);
	// Second step - reaching max steps this year, wait next year
	isDispersing = ind2.moveStep(&ls, &sp, hab_index, false);
	assert(ind2.getCurrCell() == init_cell);
	assert(ind2.getStatus() == 3);
	ind2.setStatus(1); // dispersing again
	// Third step - reaching max steps, dies in unsuitable cell
	isDispersing = ind2.moveStep(&ls, &sp, hab_index, false);
	assert(ind2.getCurrCell() == init_cell);
	assert(ind2.getStatus() == 6);

	// Step length too long
	m.stepLength = ls_params.dimX * SQRT2 * 1.5; // overshoots
	sp.setSpMovtTraits(m);
	Individual ind3(&sp, init_cell, init_patch, 0, 0, 0, 0.0, true, 2);
	ind3.setStatus(1); // dispersing
	steps.minSteps = 1;
	steps.maxStepsYr = 1;
	steps.maxSteps = 1; // no need to test more than one step this time
	sp.setSteps(0, 0, steps);
	isDispersing = ind3.moveStep(&ls, &sp, hab_index, false);
	assert(ind3.getCurrCell() == init_cell);
	assert(ind3.getStatus() == 6);

	// Adequate step length
	m.stepLength = (ls_params.dimX - 1) * SQRT2;
	sp.setSpMovtTraits(m);
	Individual ind4(&sp, init_cell, natalPatch, 0, 0, 0, 0.0, true, 2);
	ind4.setStatus(1); // dispersing
	// Initial angle still random but should eventually reach the suitable cell
	isDispersing = ind4.moveStep(&ls, &sp, hab_index, false);
	assert(ind4.getStatus() == 2);
	assert(ind4.getCurrCell() == final_cell);

	// If boundaries are absorbing however, most likely to die
	Individual ind5(&sp, init_cell, natalPatch, 0, 0, 0, 0.0, true, 2);
	ind5.setStatus(1); // dispersing
	bool absorbing_boundaries = true;
	isDispersing = ind5.moveStep(&ls, &sp, hab_index, absorbing_boundaries);
	assert(ind5.getStatus() == 6);
	assert(ind5.getCurrCell() == 0); // deref apparently

	// Correlation parameter
	// If rho = 1 should move in a straight line
	// Many cells, all suitable
	lsDim = 10;
	ls_params = createDefaultLandParams(lsDim);
	ls = Landscape();
	ls.setLandParams(ls_params, true);
	cell_vec.clear();
	for (int x = ls_params.minX; x < ls_params.dimX; ++x) {
		for (int y = ls_params.minY; y < ls_params.dimY; ++y) {
			cell_vec.push_back(new Cell(x, y, 0, hab_suitable));
		}
	}
	ls.setCellArray();
	for (auto c : cell_vec) ls.addCellToLand(c);
	ls.allocatePatches(&sp);
	ls.updateCarryingCapacity(&sp, 0, 0);
	// Individual moves by 1 along the diagonal
	m.stepLength = cellDiagLength; // 1 diagonal cell
	m.rho = 1; // angle = previous angle
	sp.setSpMovtTraits(m);
	steps.maxStepsYr = steps.maxSteps = ls_params.dimX;
	sp.setSteps(0, 0, steps);
	Individual ind6(&sp, cell_vec[0], natalPatch, 0, 0, 0, 0.0, true, 2);
	const float diag_angle = PI / 4.0; // 45 degrees
	ind6.setInitAngle(diag_angle);
	// Individual moves only along diagonal cells
	for (int i = 1; i < ls_params.dimX; ++i) {
		ind6.setStatus(1); // dispersing
		isDispersing = ind6.moveStep(&ls, &sp, hab_index, false);
		assert(ind6.getStatus() == 2);
		assert(ind6.getCurrCell() == cell_vec[i * (ls_params.dimX + 1)]);
	}
}

void testGenetics() {

	// Individuals inherit alleles from their parents
	{
		// 1 - Diploid case, 2 parents with diff. alleles
		{
			const int genomeSz = 5;
			const bool isDiploid{ true };
			SpeciesTrait* spTr = createTestEmigSpTrait(createTestGenePositions(genomeSz), isDiploid);

			// Heterozygote parent genotypes
			const float valAlleleMotherA(1.0), valAlleleMotherB(2.0);
			const float valAlleleFatherA(3.0), valAlleleFatherB(4.0);
			auto motherGenotype = createTestGenotype(genomeSz, isDiploid, valAlleleMotherA, valAlleleMotherB);
			auto fatherGenotype = createTestGenotype(genomeSz, isDiploid, valAlleleFatherA, valAlleleFatherB);
			const int startMotherChr = 0; // Strand A
			const int startFatherChr = 1; // Strand B

			// Create individual trait objects
			DispersalTrait dispTrParent(spTr); // initialisation constructor
			DispersalTrait dispTrChild(dispTrParent); // inheritance constructor

			set<unsigned int> recomSites{}; // no recombination

			dispTrChild.triggerInherit(true, motherGenotype, recomSites, startMotherChr);
			dispTrChild.triggerInherit(false, fatherGenotype, recomSites, startFatherChr);

			// Child should have inherited strand A from mother and strand B from father
			float valChildAlleleA, valChildAlleleB;
			for (int i = 0; i < genomeSz; i++) {
				valChildAlleleA = dispTrChild.getAlleleValueAtLocus(0, i);
				assert(valChildAlleleA == valAlleleMotherA);
				valChildAlleleB = dispTrChild.getAlleleValueAtLocus(1, i);
				assert(valChildAlleleB == valAlleleFatherB);
			}
		}

		// 2 - Haploid case, simply copy parent
		{
			const int genomeSz = 5;
			const bool isDiploid{ false };
			SpeciesTrait* spTr = createTestEmigSpTrait(createTestGenePositions(genomeSz), isDiploid);

			// Heterozygote parent genotypes
			const float valAlleleMother(5.0);
			auto motherGenotype = createTestGenotype(genomeSz, isDiploid, valAlleleMother);
			const int startMotherChr = 999; // doesn't matter, not used

			// Create individual trait objects
			DispersalTrait dispTrParent(spTr); // initialisation constructor
			DispersalTrait dispTrChild(dispTrParent); // inheritance constructor

			set<unsigned int> recomSites;
			for (unsigned int i = 0; i < genomeSz; i++)
				recomSites.insert(i); // recombination should be ignored

			dispTrChild.triggerInherit(true, motherGenotype, recomSites, startMotherChr);

			// Child should have inherited strand B from mother
			float valChildAllele;
			for (int i = 0; i < genomeSz; i++) {
				valChildAllele = dispTrChild.getAlleleValueAtLocus(0, i);
				assert(valChildAllele == valAlleleMother);
			}
		}
	}

	// Recombination occurs just after selected recombination site
	{
		const int genomeSz = 3;
		const bool isDiploid{ true };
		SpeciesTrait* spTr = createTestEmigSpTrait(createTestGenePositions(genomeSz), isDiploid);

		// Create individual trait objects
		DispersalTrait dispTrParent(spTr); // initialisation constructor

		// Fill mother gene sequence
		const float valAlleleA(0.0), valAlleleB(1.0);
		// Mother genotype:
		// 000
		// 111
		auto motherGenotype = createTestGenotype(genomeSz, isDiploid, valAlleleA, valAlleleB);

		// Trigger inheritance from mother
		const vector<unsigned int> recombinationSites{0, 1, genomeSz - 1}; // should work for any position
		const int startingChr{0}; // Chromosome A

		for (unsigned int site : recombinationSites) {
			DispersalTrait dispTrChild(dispTrParent); // inheritance constructor
			dispTrChild.triggerInherit(true, motherGenotype, set<unsigned int>{site}, startingChr);
			// for this test we ignore inheritance from father

			// Mother-inherited alleles should have value of chr A before recombination site, chr B after
			float valMotherAllele;
			for (int i = 0; i < genomeSz; i++) {
				valMotherAllele = dispTrChild.getAlleleValueAtLocus(0, i);
				assert(valMotherAllele == (i <= site ? valAlleleA : valAlleleB));
				// don't check other chromosome, empty bc we did not resolve father inheritance
			}
		}
	}

	// Traits that are not inherited sample their mutations in initial distribution
	{
		const float initialAlleleVal = 0.5;
		const float parentAlleleVal = 0.8;
		const bool isInherited = false;

		const int genomeSz = 5;
		const bool isDiploid{ false };

		// Create species trait
		const map<GenParamType, float> distParams{
			// all alleles initialised at a single value
			pair<GenParamType, float>{GenParamType::MIN, initialAlleleVal},
			pair<GenParamType, float>{GenParamType::MAX, initialAlleleVal}
		};
		SpeciesTrait* spTr = new SpeciesTrait(
			TraitType::E_D0,
			sex_t::NA,
			createTestGenePositions(genomeSz),
			ExpressionType::ADDITIVE,
			DistributionType::UNIFORM, distParams,
			DistributionType::NONE, distParams, // no dominance, not used
			isInherited,
			0.0, // no mutations
			DistributionType::UNIFORM, distParams, // ignored
			DistributionType::NONE, distParams, // no dominance, not used
			isDiploid ? 2 : 1,
			false
		);

		// Create individual trait objects
		DispersalTrait dispTrParent(spTr); // initialisation constructor
		auto parentGenotype = createTestGenotype(genomeSz, isDiploid, parentAlleleVal);
		const int startMotherChr = 999; // doesn't matter, not used
		dispTrParent.overwriteGenes(parentGenotype);
		assert(dispTrParent.getAlleleValueAtLocus(0, 0) == parentAlleleVal);

		DispersalTrait dispTrChild(dispTrParent); // inheritance constructor

		set<unsigned int> recomSites{genomeSz - 1};
		dispTrChild.triggerInherit(true, parentGenotype, recomSites, startMotherChr);

		// Child does not inherit from mother,
		// instead allele value is initial value
		float valChildAllele;
		for (int i = 0; i < genomeSz; i++) {
			valChildAllele = dispTrChild.getAlleleValueAtLocus(0, i);
			assert(valChildAllele == initialAlleleVal);
		}
	}

	// Genetic fitness mutations are constrained between -1 or 1
	// + sampled dominance coefficients are always positive
	{
		{
			// Most likely (~96%) to sample a mutation > 1
			const float gammaMutShapeParam = 5.0;
			const float gammaMutScaleParam = 1.0;
			//  Normal centered on 0 : ~50% of sampling negative dominance coefficient
			const float dominanceMeanParam = 0.0;
			const float dominanceSdParam = 1.0;

			const int genomeSz = 5;
			const bool isDiploid{ false };

			// Create species trait
			const map<GenParamType, float> mutationParams{
				pair<GenParamType, float>{GenParamType::SHAPE, gammaMutShapeParam},
				pair<GenParamType, float>{GenParamType::SCALE, gammaMutScaleParam}
			};
			const map<GenParamType, float> dominanceParams{
				pair<GenParamType, float>{GenParamType::MEAN, dominanceMeanParam},
				pair<GenParamType, float>{GenParamType::SD, dominanceSdParam}
			};
			const map<GenParamType, float> placeholderParams = mutationParams;

			SpeciesTrait* spTr = new SpeciesTrait(
				TraitType::GENETIC_LOAD1,
				sex_t::NA,
				createTestGenePositions(genomeSz),
				ExpressionType::MULTIPLICATIVE,
				DistributionType::NONE, placeholderParams, // not used for genetic load
				DistributionType::NONE, dominanceParams,
				true,
				1.0, // every site mutates
				DistributionType::GAMMA, mutationParams,
				DistributionType::NORMAL, dominanceParams,
				isDiploid ? 2 : 1,
				false
			);

			// Create individual trait object
			GeneticFitnessTrait traitInd(spTr); // initialisation constructor
			traitInd.mutate();
			for (int i = 0; i < genomeSz; i++) {
				float valAllele = traitInd.getAlleleValueAtLocus(0, i);
				assert(valAllele <= 1.0);
				float domCoef = traitInd.getDomCoefAtLocus(0, i);
				assert(domCoef >= 0.0);
			}
		}

		{
			// 1/2 chance to sample a mutation < -1
			const float normalMutMeanParam = 0;
			const float normalMutSdParam = 1.0;

			const int genomeSz = 10;
			const bool isDiploid{ false };

			// Create species trait
			const map<GenParamType, float> mutationParams{
				// all alleles initialised at a single value
				pair<GenParamType, float>{GenParamType::MEAN, normalMutMeanParam},
				pair<GenParamType, float>{GenParamType::SD, normalMutSdParam}
			};
			const map<GenParamType, float> placeholderParams = mutationParams;

			SpeciesTrait* spTr = new SpeciesTrait(
				TraitType::GENETIC_LOAD1,
				sex_t::NA,
				createTestGenePositions(genomeSz),
				ExpressionType::MULTIPLICATIVE,
				DistributionType::NONE, placeholderParams, // not used for genetic load
				DistributionType::NONE, placeholderParams, // doesn't matter for this test
				true,
				1.0, // every site mutates
				DistributionType::NORMAL, mutationParams,
				DistributionType::NORMAL, placeholderParams, // doesn't matter for this test
				isDiploid ? 2 : 1,
				false
			);

			// Create individual trait object
			GeneticFitnessTrait traitInd(spTr); // initialisation constructor
			traitInd.mutate();
			for (int i = 0; i < genomeSz; i++) {
				float valAllele = traitInd.getAlleleValueAtLocus(0, i);
				assert(valAllele > -1.0);
			}
		}
	}
}

bool haveSameEmigD0Allele(const Individual& indA, const Individual& indB, const int& position, short whichHaplo = 0) {
	return indA.getTrait(E_D0)->getAlleleValueAtLocus(whichHaplo, position)
		== indB.getTrait(E_D0)->getAlleleValueAtLocus(whichHaplo, position);
}

void testIndividual() {

	// Kernel-based transfer
	testTransferKernels();

	// Correlated random walk (CRW)
	testTransferCRW();

	testGenetics();

	// Genetic linkage + Chromosome breaks
		// Considering diallelic genes A, B, C, D with:
		// A, B, C are on chr.1, D is on chr.2
		// A, B are adjacent, C sits on the other end of chr.1
		// C, D have adjacent positions in genome (but are on separate chr.)
		// AB------------C//D
	// We simulate 100 inheritance + recombination processes and expect that:
		// 1. freq(A,B have same alleles) >> freq(A,C have same alleles)
		// 2. 0.65 > freq(C,D have same alleles) > 0.35 despite being adjacent because of chrom. break
		// (both freq. have p < 0.001 from a binomial with p 0.5 and 100 trials)
	{
		Patch* pPatch = new Patch(0, 0);
		Cell* pCell = new Cell(0, 0, (intptr)pPatch, 0);

		const float recombinationRate = 0.01;
		const int genomeSz = 10;

		Species* pSpecies = createDefaultSpecies();
		pSpecies->setGeneticParameters(
			set<int>{genomeSz-2, genomeSz - 1}, // two chromosomes
			genomeSz,
			recombinationRate,
			set<int>{}, "none", set<int>{}, 0 // no output so no sampling
		);

		int posA = 0;
		int posB = 1;
		int posC = genomeSz - 2;
		int posD = genomeSz - 1;
		const set<int> genePositions = {posA, posB, posC, posD};
		const bool isDiploid{ true };
		SpeciesTrait* spTr = createTestEmigSpTrait(genePositions, isDiploid);
		pSpecies->addTrait(TraitType::E_D0, *spTr);

		Individual indMother = Individual(pCell, pPatch, 0, 0, 0, 0.0, false, 0);
		Individual indFather = Individual(pCell, pPatch, 0, 0, 0, 1.0, false, 0);
		indMother.setUpGenes(pSpecies, 1.0);
		indFather.setUpGenes(pSpecies, 1.0);

		int countRecombineTogetherAB = 0;
		int countRecombineTogetherAC = 0;
		int countRecombineTogetherCD = 0;

		const int nbTrials = 100;
		for (int i = 0; i < nbTrials; ++i)
		{
			Individual indChild = Individual(pCell, pPatch, 0, 0, 0, 0.0, false, 0);
			indChild.inheritTraits(pSpecies, &indMother, &indFather, 1.0);

			bool hasInheritedA0 = haveSameEmigD0Allele(indChild, indMother, posA);
			bool hasInheritedB0 = haveSameEmigD0Allele(indChild, indMother, posB);
			bool hasInheritedC0 = haveSameEmigD0Allele(indChild, indMother, posC);
			bool hasInheritedD0 = haveSameEmigD0Allele(indChild, indMother, posD);

			countRecombineTogetherAB += (hasInheritedA0 && hasInheritedB0)
				|| (!hasInheritedA0 && !hasInheritedB0);
			countRecombineTogetherAC += (hasInheritedA0 && hasInheritedC0)
				|| (!hasInheritedA0 && !hasInheritedC0);
			countRecombineTogetherCD += (hasInheritedC0 && hasInheritedD0)
				|| (!hasInheritedC0 && !hasInheritedD0);
		}
		assert(countRecombineTogetherAB > countRecombineTogetherAC);
		assert(35 < countRecombineTogetherCD && countRecombineTogetherCD < 65);
	}

	// Sex-specific traits use the correct genes
	/// Set up a sex-dependent emigration probability with male and female loci
	/// Emigration probability is 1 initially, but female trait mutates.
	{
		Patch* pPatch = new Patch(0, 0);
		Cell* pCell = new Cell(0, 0, (intptr)pPatch, 0);

		// Species-level paramters
		const int genomeSz = 6;
		Species* pSpecies = createDefaultSpecies();
		pSpecies->setGeneticParameters(
			set<int>{genomeSz - 1}, // one chromosome
			genomeSz,
			0.0, // no recombination
			set<int>{}, "none", set<int>{}, 0 // no output so no sampling
		);
		emigRules emig;
		emig.indVar = true;
		emig.sexDep = true;
		emig.densDep = false;
		pSpecies->setEmigRules(emig);

		// Create species trait
		//// Shared params
		const map<GenParamType, float> initParams{
			// all values are 1
			pair<GenParamType, float>{GenParamType::MIN, 1.0},
			pair<GenParamType, float>{GenParamType::MAX, 1.0}
		};
		const map<GenParamType, float> mutationParams{
			// reduction of emigration probability
			pair<GenParamType, float>{GenParamType::MIN, -0.5},
			pair<GenParamType, float>{GenParamType::MAX, -0.25}
		};
		const bool isDiploid{ true };

		//// Sex-specific params
		const set<int> maleGenePositions = { 0, 2, 4 };
		const set<int> femaleGenePositions = { 1, 3, 5 };
		const float maleMutationRate = 0.0;
		const float femaleMutationRate = 1.0;

		SpeciesTrait* spTrM = new SpeciesTrait(
			TraitType::E_D0_M,
			sex_t::MAL,
			maleGenePositions,
			ExpressionType::AVERAGE,
			// Set all initial alleles values to 1
			DistributionType::UNIFORM, initParams,
			DistributionType::NONE, initParams, // no dominance, params are ignored
			true, // isInherited
			maleMutationRate, // does not mutate
			DistributionType::UNIFORM, mutationParams, // not used
			DistributionType::NONE, initParams, // no dominance, params are ignored
			isDiploid ? 2 : 1,
			false
		);
		pSpecies->addTrait(TraitType::E_D0_M, *spTrM);

		SpeciesTrait* spTrF = new SpeciesTrait(
			TraitType::E_D0_F,
			sex_t::FEM,
			femaleGenePositions,
			ExpressionType::AVERAGE,
			// Set all initial alleles values to 1
			DistributionType::UNIFORM, initParams,
			DistributionType::NONE, initParams, // no dominance, params are ignored
			true, // isInherited
			femaleMutationRate, // does mutate
			DistributionType::UNIFORM, mutationParams, // not used
			DistributionType::NONE, initParams,
			isDiploid ? 2 : 1,
			false
		);
		pSpecies->addTrait(TraitType::E_D0_F, *spTrF);

		// Set up male and female individuals, trigger mutations
		Individual indFemale = Individual(pCell, pPatch, 0, 0, 0, 0.0, false, 0);
		Individual indMale = Individual(pCell, pPatch, 0, 0, 0, 1.0, false, 0);
		indFemale.setUpGenes(pSpecies, 1.0);
		indMale.setUpGenes(pSpecies, 1.0);
		indFemale.triggerMutations(pSpecies);
		indMale.triggerMutations(pSpecies);

		// Male should use male trait, still 1
		// Female should use female trait, has mutated
		emigTraits femaleEmig = indFemale.getIndEmigTraits();
		emigTraits maleEmig = indMale.getIndEmigTraits();
		assert(femaleEmig.d0 != 1.0);
		assert(maleEmig.d0 == 1.0);
	}

	// Individuals use species-level trait when not individual variable,
	// and individual-level trait when trait is individual variable
	{
		float spEmigProb = 1.0;
		float indEmigProb = 0.0;

		Patch* pPatch = new Patch(0, 0);
		Cell* pCell = new Cell(0, 0, (intptr)pPatch, 0);

		// Species-level paramters
		const int genomeSz = 1;
		Species* pSpecies = createDefaultSpecies();
		pSpecies->setGeneticParameters(
			set<int>{genomeSz - 1}, // one chromosome
			genomeSz,
			0.0, // no recombination
			set<int>{}, "none", set<int>{}, 0 // no output so no sampling
		);
		emigRules emig;
		emig.indVar = false;
		emig.stgDep = false; emig.sexDep = false; emig.densDep = false;
		pSpecies->setEmigRules(emig);

		emigTraits spEmigTraits; spEmigTraits.d0 = spEmigProb;
		pSpecies->setSpEmigTraits(0, 0, spEmigTraits);

		// Create species trait
		const map<GenParamType, float> initParams{
			// Emigration probability is always 0
			pair<GenParamType, float>{GenParamType::MIN, indEmigProb},
			pair<GenParamType, float>{GenParamType::MAX, indEmigProb}
		};
		const map<GenParamType, float> mutationParams = initParams; // doesn't matter, not used
		SpeciesTrait* spTr = new SpeciesTrait(
			TraitType::E_D0,
			sex_t::NA,
			set<int>{ 0 }, // only one locus
			ExpressionType::AVERAGE,
			DistributionType::UNIFORM, initParams,
			DistributionType::NONE, initParams, // no dominance, params are ignored
			true, // isInherited
			0.0, // no mutation
			DistributionType::UNIFORM, mutationParams, // not used
			DistributionType::NONE, initParams, // no dominance, params are ignored
			2, // diploid
			false
		);
		pSpecies->addTrait(TraitType::E_D0, *spTr);

		Individual ind = Individual(pCell, pPatch, 0, 0, 0, 0.0, false, 0);
		ind.setUpGenes(pSpecies, 1.0);

		// Create population to trigger emigration selection
		Population pop(pSpecies, pPatch, 0, 1.0);
		pop.recruit(&ind);
		assert(ind.getStatus() == 0);
		pop.emigration(100.0);

		// Individual is using the species-wide emigration prob,
		// so should be selected to emigrate (status 1)
		assert(ind.getStatus() == 1);

		// Change rules to use individual-variable trait
		ind.setStatus(0);
		emig.indVar = true;
		pSpecies->setEmigRules(emig);
		pop.emigration(100.0);

		// Individual-level emig prob is zero, must not be emigrating;
		assert(ind.getStatus() == 0);

		pop.clearInds(); // empty inds vector to not deallocate individual twice
	}

	// Individuals with genetic fitness = 1 are always viable
	// Individuals with genetic fitness = 0 are never viable
	{
		Patch* pPatch = new Patch(0, 0);
		Cell* pCell = new Cell(0, 0, (intptr)pPatch, 0);

		// Species-level paramters
		const int genomeSz = 1;
		Species* pSpecies = createDefaultSpecies();
		pSpecies->setGeneticParameters(
			set<int>{genomeSz - 1}, // one chromosome
			genomeSz,
			0.0, // no recombination
			set<int>{}, "none", set<int>{}, 0 // no output so no sampling
		);

		// Create species trait
		const map<GenParamType, float> mutationParams{
			// Always sample 1 = a lethal mutation
			pair<GenParamType, float>{GenParamType::MIN, 1.0},
			pair<GenParamType, float>{GenParamType::MAX, 1.0}
		};
		const map<GenParamType, float> domParams = mutationParams; // all dominance parameters are equal
		const map<GenParamType, float> initParams = mutationParams; // doesn't matter, not used

		SpeciesTrait* spTr = new SpeciesTrait(
			TraitType::GENETIC_LOAD1,
			sex_t::NA,
			set<int>{ 0 }, // only one locus
			ExpressionType::MULTIPLICATIVE,
			DistributionType::NONE, initParams,
			DistributionType::UNIFORM, domParams, // no dominance, params are ignored
			true, // isInherited
			1.0, // will mutate
			DistributionType::UNIFORM, mutationParams, // lethal mutation
			DistributionType::NONE, initParams,
			2, // diploid
			false
		);
		pSpecies->addTrait(TraitType::GENETIC_LOAD1, *spTr);

		Individual ind = Individual(pCell, pPatch, 0, 0, 0, 0.0, false, 0);
		ind.setUpGenes(pSpecies, 1.0);

		// By default, all loci are initialised at 0 so individual is viable
		assert(ind.getGeneticFitness() == 1.0);
		assert(ind.isViable());

		ind.triggerMutations(pSpecies);

		// Individual now bears a lethal allele
		assert(ind.getGeneticFitness() == 0.0);
		assert(!ind.isViable());
	}

	// A largely dominant alleles overrides the expression of its homologue
	{
		Patch* pPatch = new Patch(0, 0);
		Cell* pCell = new Cell(0, 0, (intptr)pPatch, 0);

		// Species-level paramters
		const int genomeSz = 1;
		Species* pSpecies = createDefaultSpecies();
		pSpecies->setGeneticParameters(
			set<int>{genomeSz - 1}, // one chromosome
			genomeSz,
			0.0, // no recombination
			set<int>{}, "none", set<int>{}, 0 // no output so no sampling
		);

		// Create species trait
		const map<GenParamType, float> distParams{
			pair<GenParamType, float>{GenParamType::MIN, 0.0},
			pair<GenParamType, float>{GenParamType::MAX, 0.0}
		};

		SpeciesTrait* spTr = new SpeciesTrait(
			TraitType::GENETIC_LOAD1,
			sex_t::NA,
			set<int>{ 0 }, // only one locus
			ExpressionType::MULTIPLICATIVE,
			DistributionType::NONE, distParams,
			DistributionType::UNIFORM, distParams, // no dominance, params are ignored
			true, // isInherited
			0.0, // no mutation
			DistributionType::UNIFORM, distParams, // lethal mutation
			DistributionType::NONE, distParams,
			2, // diploid
			false
		);
		pSpecies->addTrait(TraitType::GENETIC_LOAD1, *spTr);

		Individual ind = Individual(pCell, pPatch, 0, 0, 0, 0.0, false, 0);
		ind.setUpGenes(pSpecies, 1.0);

		const float valAlleleA = 1.0; // lethal
		const float valAlleleB = 0.0; // mild
		float domCoefA = 0.00001; // highly recessive
		float domCoefB = 1.0;
		auto viableGenotype = createTestGenotype(genomeSz, true, valAlleleA, valAlleleB, domCoefA, domCoefB);
		ind.overrideGenotype(GENETIC_LOAD1, viableGenotype);
		ind.triggerMutations(pSpecies); // no mutations (rate = 0) but updates genetic fitness
		assert(ind.isViable());

		domCoefA = 10000.0; // oh gosh now it's lethal AND dominant
		auto lethalGenotype = createTestGenotype(genomeSz, true, valAlleleA, valAlleleB, domCoefA, domCoefB);
		ind.overrideGenotype(GENETIC_LOAD1, lethalGenotype);
		ind.triggerMutations(pSpecies);
		assert(!ind.isViable());
	}

	// Dispersal trait alleles can take any value, but
	// phenotypes are constrained to valid values
	{
		// Case 1 - Settlement + Kernel-based transfer
		{
			const int genomeSz = 4;

			const bool isDiploid{ true }; // haploid, simpler check
			const float mutationRate = 0.0; // no mutations

			Patch* pPatch = new Patch(0, 0);
			Cell* pCell = new Cell(0, 0, (intptr)pPatch, 0);

			// Genome-level settings
			Species* pSpecies = createDefaultSpecies();
			pSpecies->setGeneticParameters(
				set<int>{genomeSz - 1}, // one chromosome
				genomeSz,
				0.0, // no recombination
				set<int>{}, "none", set<int>{}, 0 // no output so no sampling
			);

			settleType sett;
			sett.indVar = true;
			sett.sexDep = false;
			sett.stgDep = false;
			pSpecies->setSettle(sett);
			settleRules settRules;
			settRules.densDep = true;
			pSpecies->setSettRules(0, 0, settRules);

			transferRules trfr;
			trfr.indVar = true;
			trfr.usesMovtProc = false;
			trfr.moveType = 0; // kernels
			trfr.twinKern = false;
			trfr.sexDep = false;
			trfr.stgDep = false;
			pSpecies->setTrfrRules(trfr);

			// Species-level traits
			const map<GenParamType, float> distParams{
				// prameters don't matter, allele values are overwritten below
				pair<GenParamType, float>{GenParamType::MIN, 1.0},
				pair<GenParamType, float>{GenParamType::MAX, 1.0}
			};
			SpeciesTrait* trSettProb = new SpeciesTrait(
				TraitType::S_S0,
				sex_t::NA,
				set<int>{ 0 },
				ExpressionType::AVERAGE,
				DistributionType::UNIFORM, distParams,
				DistributionType::NONE, distParams, // no dominance, params are ignored
				true, // isInherited
				mutationRate, // does not mutate
				DistributionType::UNIFORM, distParams, // not used
				DistributionType::NONE, distParams, // no dominance, params are ignored
				isDiploid ? 2 : 1,
				false
			);
			SpeciesTrait* trSettAlpha = new SpeciesTrait(
				TraitType::S_ALPHA,
				sex_t::NA,
				set<int>{ 1 },
				ExpressionType::AVERAGE,
				DistributionType::UNIFORM, distParams,
				DistributionType::NONE, distParams, // no dominance, params are ignored
				true, // isInherited
				mutationRate, // does not mutate
				DistributionType::UNIFORM, distParams, // not used
				DistributionType::NONE, distParams,
				isDiploid ? 2 : 1,
				false
			);
			SpeciesTrait* trSettBeta = new SpeciesTrait(
				TraitType::S_BETA,
				sex_t::NA,
				set<int>{ 2 },
				ExpressionType::AVERAGE,
				DistributionType::UNIFORM, distParams,
				DistributionType::NONE, distParams, // no dominance, params are ignored
				true, // isInherited
				mutationRate, // does not mutate
				DistributionType::UNIFORM, distParams, // not used
				DistributionType::NONE, distParams, // no dominance, params are ignored
				isDiploid ? 2 : 1,
				false
			);
			SpeciesTrait* trMeanKern = new SpeciesTrait(
				TraitType::KERNEL_MEANDIST_1,
				sex_t::NA,
				set<int>{ 3 },
				ExpressionType::ADDITIVE,
				DistributionType::UNIFORM, distParams,
				DistributionType::NONE, distParams, // no dominance, params are ignored
				true, // isInherited
				mutationRate, // does not mutate
				DistributionType::UNIFORM, distParams, // not used
				DistributionType::NONE, distParams, // no dominance, params are ignored
				isDiploid ? 2 : 1,
				false
			);

			pSpecies->addTrait(TraitType::S_S0, *trSettProb);
			pSpecies->addTrait(TraitType::S_ALPHA, *trSettAlpha);
			pSpecies->addTrait(TraitType::S_BETA, *trSettBeta);
			pSpecies->addTrait(TraitType::KERNEL_MEANDIST_1, *trMeanKern);

			Individual ind = Individual(pCell, pPatch, 0, 0, 0, 0.0, false, 0);
			ind.setUpGenes(pSpecies, 1.0);

			// Overwrite genotypes with alleles resulting in invalid phenotypes
			auto sProbGenotype = createTestGenotype(genomeSz, true, 1.1, 1.2);
			auto sAlphaGenotype = createTestGenotype(genomeSz, true, -1.5, -0.5);
			auto kernelDistGenotype = createTestGenotype(genomeSz, true, -0.2, -0.4);
			ind.overrideGenotype(S_S0, sProbGenotype);
			ind.overrideGenotype(S_ALPHA, sAlphaGenotype);
			ind.overrideGenotype(KERNEL_MEANDIST_1, kernelDistGenotype);
			// any value is valid for settlement beta

			ind.triggerMutations(pSpecies); // no mutations, but trigger expression

			settleTraits settTr = ind.getIndSettTraits();
			assert(settTr.s0 <= 1.0);
			assert(settTr.alpha > 0.0);
			kernelData trfrTr = *(static_cast<kernelData*>(ind.getTrfrData()));
			assert(trfrTr.meanDist1 >= 0.0);
		}

		// Case 2 - Correlated random walk transfer
		{
			const int genomeSz = 2;

			const bool isDiploid{ true }; // haploid, simpler check
			const float mutationRate = 0.0; // no mutations

			Patch* pPatch = new Patch(0, 0);
			Cell* pCell = new Cell(0, 0, (intptr)pPatch, 0);

			// Genome-level settings
			Species* pSpecies = createDefaultSpecies();
			pSpecies->setGeneticParameters(
				set<int>{genomeSz - 1}, // one chromosome
				genomeSz,
				0.0, // no recombination
				set<int>{}, "none", set<int>{}, 0 // no output so no sampling
			);

			settleType sett;
			sett.indVar = false;
			pSpecies->setSettle(sett);

			transferRules trfr;
			trfr.indVar = true;
			trfr.usesMovtProc = true;
			trfr.moveType = 2; // CRW
			pSpecies->setTrfrRules(trfr);

			// Species-level traits
			const map<GenParamType, float> distParams{
				// prameters don't matter, allele values are overwritten below
				pair<GenParamType, float>{GenParamType::MIN, 1.0},
				pair<GenParamType, float>{GenParamType::MAX, 1.0}
			};
			SpeciesTrait* trCRWLen = new SpeciesTrait(
				TraitType::CRW_STEPLENGTH,
				sex_t::NA,
				set<int>{ 0 },
				ExpressionType::ADDITIVE,
				DistributionType::UNIFORM, distParams,
				DistributionType::NONE, distParams, // no dominance, params are ignored
				true, // isInherited
				mutationRate, // does not mutate
				DistributionType::UNIFORM, distParams, // not used
				DistributionType::NONE, distParams, // no dominance, params are ignored
				isDiploid ? 2 : 1,
				false
			);
			SpeciesTrait* trCRWCorr = new SpeciesTrait(
				TraitType::CRW_STEPCORRELATION,
				sex_t::NA,
				set<int>{ 1 },
				ExpressionType::AVERAGE,
				DistributionType::UNIFORM, distParams,
				DistributionType::NONE, distParams, // no dominance, params are ignored
				true, // isInherited
				mutationRate, // does not mutate
				DistributionType::UNIFORM, distParams, // not used
				DistributionType::NONE, distParams, // no dominance, params are ignored
				isDiploid ? 2 : 1,
				false
			);

			pSpecies->addTrait(TraitType::CRW_STEPLENGTH, *trCRWLen);
			pSpecies->addTrait(TraitType::CRW_STEPCORRELATION, *trCRWCorr);

			bool usesMovtProcess = true;
			short whichMovtProcess = 2; // CRW
			Individual ind = Individual(pCell, pPatch, 0, 0, 0, 0.0, usesMovtProcess, whichMovtProcess);
			ind.setUpGenes(pSpecies, 1.0);

			// Overwrite genotypes with alleles resulting in invalid phenotypes
			auto crwLenGenotype = createTestGenotype(genomeSz, true, -1.1, -1.2);
			auto crwCorrGenoType = createTestGenotype(genomeSz, true, 1.5, 2.0);
			ind.overrideGenotype(CRW_STEPLENGTH, crwLenGenotype);
			ind.overrideGenotype(CRW_STEPCORRELATION, crwCorrGenoType);

			ind.triggerMutations(pSpecies); // no mutations, but trigger expression

			crwData trfrTr = *(static_cast<crwData*>(ind.getTrfrData()));
			assert(trfrTr.stepLength >= 0.0);
			assert(trfrTr.rho <= 1.0);
		}

		// Case 3 - Transfer with Stochastic Movement Simulator
		{
			const int genomeSz = 2;
			const bool isDiploid{ true }; // haploid, simpler check
			const float mutationRate = 0.0; // no mutations

			Patch* pPatch = new Patch(0, 0);
			Cell* pCell = new Cell(0, 0, (intptr)pPatch, 0);

			// Genome-level settings
			Species* pSpecies = createDefaultSpecies();
			pSpecies->setGeneticParameters(
				set<int>{genomeSz - 1}, // one chromosome
				genomeSz,
				0.0, // no recombination
				set<int>{}, "none", set<int>{}, 0 // no output so no sampling
			);

			settleType sett;
			sett.indVar = false;
			pSpecies->setSettle(sett);

			transferRules trfr;
			trfr.indVar = true;
			trfr.usesMovtProc = true;
			trfr.moveType = 1; // SMS
			pSpecies->setTrfrRules(trfr);

			// Species-level traits
			const map<GenParamType, float> distParams{
				// prameters don't matter, allele values are overwritten below
				pair<GenParamType, float>{GenParamType::MIN, 1.0},
				pair<GenParamType, float>{GenParamType::MAX, 1.0}
			};
			SpeciesTrait* trSmsDirPers = new SpeciesTrait(
				TraitType::SMS_DP,
				sex_t::NA,
				set<int>{ 0 },
				ExpressionType::ADDITIVE,
				DistributionType::UNIFORM, distParams,
				DistributionType::NONE, distParams, // no dominance, params are ignored
				true, // isInherited
				mutationRate, // does not mutate
				DistributionType::UNIFORM, distParams, // not used
				DistributionType::NONE, distParams, // no dominance, params are ignored
				isDiploid ? 2 : 1,
				false
			);
			SpeciesTrait* trSmsGoalBias = new SpeciesTrait(
				TraitType::SMS_GB,
				sex_t::NA,
				set<int>{ 1 },
				ExpressionType::AVERAGE,
				DistributionType::UNIFORM, distParams,
				DistributionType::NONE, distParams, // no dominance, params are ignored
				true, // isInherited
				mutationRate, // does not mutate
				DistributionType::UNIFORM, distParams, // not used
				DistributionType::NONE, distParams, // no dominance, params are ignored
				isDiploid ? 2 : 1,
				false
			);
			SpeciesTrait* trSmsAlpha = new SpeciesTrait(
				TraitType::SMS_ALPHADB,
				sex_t::NA,
				set<int>{ 0 },
				ExpressionType::ADDITIVE,
				DistributionType::UNIFORM, distParams,
				DistributionType::NONE, distParams, // no dominance, params are ignored
				true, // isInherited
				mutationRate, // does not mutate
				DistributionType::UNIFORM, distParams, // not used
				DistributionType::NONE, distParams, // no dominance, params are ignored
				isDiploid ? 2 : 1,
				false
			);
			SpeciesTrait* trSmsBeta = new SpeciesTrait(
				TraitType::SMS_BETADB,
				sex_t::NA,
				set<int>{ 1 },
				ExpressionType::AVERAGE,
				DistributionType::UNIFORM, distParams,
				DistributionType::NONE, distParams, // no dominance, params are ignored
				true, // isInherited
				mutationRate, // does not mutate
				DistributionType::UNIFORM, distParams, // not used
				DistributionType::NONE, distParams, // no dominance, params are ignored
				isDiploid ? 2 : 1,
				false
			);

			pSpecies->addTrait(TraitType::SMS_DP, *trSmsDirPers);
			pSpecies->addTrait(TraitType::SMS_GB, *trSmsGoalBias);
			pSpecies->addTrait(TraitType::SMS_ALPHADB, *trSmsAlpha);
			pSpecies->addTrait(TraitType::SMS_BETADB, *trSmsBeta);

			bool usesMovtProcess = true;
			short whichMovtProcess = 1; // SMS
			Individual ind = Individual(pCell, pPatch, 0, 0, 0, 0.0, usesMovtProcess, whichMovtProcess);
			ind.setUpGenes(pSpecies, 1.0);

			// Overwrite genotypes with alleles resulting in invalid phenotypes
			auto smsDPGenotype = createTestGenotype(genomeSz, true, 0.1, 0.2);
			auto smsGBGenotype = createTestGenotype(genomeSz, true, 0.5, 0.0);
			auto smsAlphaGenotype = createTestGenotype(genomeSz, true, 0.0, 0.0);
			ind.overrideGenotype(SMS_DP, smsDPGenotype);
			ind.overrideGenotype(SMS_GB, smsGBGenotype);
			ind.overrideGenotype(SMS_ALPHADB, smsAlphaGenotype);
			// any value is valid for SMS beta

			ind.triggerMutations(pSpecies); // no mutations, but trigger expression

			smsData trfrTr = *(static_cast<smsData*>(ind.getTrfrData()));
			assert(trfrTr.dp >= 1.0);
			assert(trfrTr.gb >= 1.0);
			assert(trfrTr.alphaDB > 0.0);
		}
	}

}

#endif //NDEBUG
