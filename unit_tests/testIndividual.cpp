#if RSDEBUG

#include "../Individual.h"

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
	Individual ind1(init_cell, init_patch, 1, 0, 0, 0.0, false, 0);
	int isDispersing = ind1.moveKernel(&ls, &sp, false);

	// After moving, individual should be in the only available cell
	Cell* curr_cell = ind1.getCurrCell();
	assert(curr_cell != init_cell);
	assert(curr_cell == final_cell);
	assert(ind1.getStatus() == 2); // potential settler

	// If no cell within reasonable dispersal reach, individual does not move and dies
	kern.meanDist1 = 1.0;
	sp.setSpKernTraits(0, 0, kern, ls_params.resol);
	Individual ind2(init_cell, init_patch, 1, 0, 0, 0.0, false, 0); // reset individual
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
	Individual ind3(init_cell, init_patch, 1, 0, 0, 0.0, false, 0);
	isDispersing = ind3.moveKernel(&ls, &sp, false);
	assert(ind3.getStatus() == 6); // dead, could not reach destination cell
	kern.probKern1 = 0.0; // always use second kernel
	sp.setSpKernTraits(0, 0, kern, ls_params.resol);
	Individual ind4(init_cell, init_patch, 1, 0, 0, 0.0, false, 0);
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

	Individual ind5(init_cell, init_patch, 1, 0, 0, 0.0, false, 0); // female as default
	isDispersing = ind5.moveKernel(&ls, &sp, false);
	assert(ind5.getStatus() == 6); // dead, could not reach destination

	Individual ind6(init_cell, init_patch, 1, 0, 0, 1.0, false, 0); // male
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

	Individual ind7(init_cell, init_patch, 0, 0, 0, 0.0, false, 0); // juvenile
	isDispersing = ind7.moveKernel(&ls, &sp, false);
	assert(ind7.getStatus() == 6);

	Individual ind8(init_cell, init_patch, 1, 0, 0, 0.0, false, 0); // adult by default
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
	Individual ind9(init_cell, init_patch, 1, 0, 0, 0.0, false, 0); // reset individual

	// Non-absorbing boundaries
	bool absorbing_boundaries{ false };
	isDispersing = ind9.moveKernel(&ls, &sp, absorbing_boundaries);
	curr_cell = ind9.getCurrCell();
	assert(curr_cell != init_cell); // ...should be able to move eventually
	assert(ind9.getStatus() == 2);

	// Absorbing boundaries
	absorbing_boundaries = true;
	Individual ind10(init_cell, init_patch, 1, 0, 0, 0.0, false, 0); // reset individual
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
	Individual ind11(init_cell, init_patch, 1, 0, 0, 0.0, false, 0);
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
	Individual ind12(init_cell, init_patch, 1, 0, 0, 0.0, false, 0);
	isDispersing = ind12.moveKernel(&ls, &sp, false);
	assert(ind12.getStatus() == 7);
	mort.mortBeta = 30; // very large distance, unlikely to draw
	sp.setMortParams(mort);
	Individual ind13(init_cell, init_patch, 1, 0, 0, 0.0, false, 0);
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
	Individual ind0(init_cell, init_patch, 1, 0, 0, 0.0, true, 2);

	// Set status
	assert(ind0.getStatus() == 0); // default status, not emigrating
	int isDispersing = ind0.moveStep(&ls, &sp, hab_index, false);
	assert(ind0.getStatus() == 0); // status didn't change
	assert(ind0.getCurrCell() == init_cell); // not emigrating so didn't move

	// Per-step mortality
	m.stepMort = 1.0; // should die 
	sp.setSpMovtTraits(m);
	Individual ind1(init_cell, init_patch, 0, 0, 0, 0.0, true, 2);
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
	Individual ind2(init_cell, natalPatch, 0, 0, 0, 0.0, true, 2);
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
	Individual ind3(init_cell, init_patch, 0, 0, 0, 0.0, true, 2);
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
	Individual ind4(init_cell, natalPatch, 0, 0, 0, 0.0, true, 2);
	ind4.setStatus(1); // dispersing
	// Initial angle still random but should eventually reach the suitable cell
	isDispersing = ind4.moveStep(&ls, &sp, hab_index, false);
	assert(ind4.getStatus() == 2);
	assert(ind4.getCurrCell() == final_cell);

	// If boundaries are absorbing however, most likely to die
	Individual ind5(init_cell, natalPatch, 0, 0, 0, 0.0, true, 2);
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
	Individual ind6(cell_vec[0], natalPatch, 0, 0, 0, 0.0, true, 2);
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
			SpeciesTrait* spTr = createTestSpTrait(createTestGenePositions(genomeSz), isDiploid);

			// Heterozygote parent genotypes
			const float valAlleleMotherA(1.0), valAlleleMotherB(2.0);
			const float valAlleleFatherA(3.0), valAlleleFatherB(4.0);
			auto motherGenotype = createTestEmigTrGenotype(genomeSz, isDiploid, valAlleleMotherA, valAlleleMotherB);
			auto fatherGenotype = createTestEmigTrGenotype(genomeSz, isDiploid, valAlleleFatherA, valAlleleFatherB);
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
			SpeciesTrait* spTr = createTestSpTrait(createTestGenePositions(genomeSz), isDiploid);

			// Heterozygote parent genotypes
			const float valAlleleMother(5.0);
			auto motherGenotype = createTestEmigTrGenotype(genomeSz, isDiploid, valAlleleMother);
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

	// Recombination occurs just after recombination site
	{
		const int genomeSz = 5;
		const bool isDiploid{ true };
		SpeciesTrait* spTr = createTestSpTrait(createTestGenePositions(genomeSz), isDiploid);

		// Create individual trait objects
		DispersalTrait dispTrParent(spTr); // initialisation constructor
		DispersalTrait dispTrChild(dispTrParent); // inheritance constructor

		// Fill mother gene sequence
		const float valAlleleA(0.0), valAlleleB(1.0);
		// Mother genotype:
		// 000000
		// 111111
		auto motherGenotype = createTestEmigTrGenotype(genomeSz, isDiploid, valAlleleA, valAlleleB);

		// Trigger inheritance from mother
		const unsigned int recombinationSite = pRandom->IRandom(0, genomeSz); // should work for any position
		const int startingChr{0}; // Chromosome A
		dispTrChild.triggerInherit(true, motherGenotype, set<unsigned int>{recombinationSite}, startingChr);
		// for this test we ignore inheritance from father

		// Mother-inherited alleles should have value of chr A before recombination site, chr B after
		float valMotherAllele;
		for (int i = 0; i < genomeSz; i++) {
			valMotherAllele = dispTrChild.getAlleleValueAtLocus(0, i);
			assert(valMotherAllele == (i <= recombinationSite ? valAlleleA : valAlleleB));
			// don't check other chromosome, empty bc we did not resolve father inheritance 
		}
	}

	// Genetic linkage does occur
	{
		// Should also test for interactions with chromosome breaks?
		// 
		// Sample X vectors of recombination sites
		// Create traits and resolve inheritance for each vector
		// Count times c1 = A & B; c2 = B & C have same alleles
		// Assert c1 > c2
		/*
		const int genomeSz = 100;
		// Genes A, B, C with dist(A,B) << dist(B,C);
		const set<int> genePositions = {0, 1, genomeSz - 1 };
		const bool isDiploid{ true };
		SpeciesTrait* spTr = createTestSpTrait(genePositions, isDiploid);
		// All genes heterozygote with same alleles for each gene
		const float valAlleleA(0.0), valAlleleB(1.0);
		const int startingChr{ 0 }; // Chromosome A

		set<unsigned int> recombinationSites;
		const int nbDraws = 10;
		
		for (int i = 0; i < nbDraws; i++) {

			if (i > 0) recombinationSites.clear();
			drawRecombinantSites(recombinationSites);

			// Create individual trait objects
			DispersalTrait dispTrParent(spTr); // initialisation constructor
			DispersalTrait dispTrChild(dispTrParent); // inheritance constructor
			auto motherGenotype = createTestEmigTrGenotype(genomeSz, isDiploid, valAlleleA, valAlleleB);

			// Trigger inheritance from mother
			dispTrChild.triggerInherit(true, motherGenotype, recombinationSites, startingChr);
			// inheritance from father not needed

			// Check linked genes have same alleles more often
			// ...

		}
		*/
		
	}
}

void testIndividual() {

	// Kernel-based transfer
	testTransferKernels();

	// Correlated random walk (CRW)
	testTransferCRW();

	testGenetics();

	// Genetic linkage does occur
		// Should also test for interactions with chromosome breaks?
		// 
		// Sample X vectors of recombination sites
		// Create traits and resolve inheritance for each vector
		// Count times c1 = A & B; c2 = B & C have same alleles
		// Assert c1 > c2
	{

		Patch* pPatch = new Patch(0, 0);
		Cell* pCell = new Cell(0, 0, (intptr)pPatch, 0);

		const float recombinationRate = 0.1;
		const int genomeSz = 10;
		Species* pSpecies = new Species();
		
		pSpecies->setGeneticParameters(
			set<int>{genomeSz - 1}, // a single chromosome
			genomeSz,
			recombinationRate,
			set<int>{}, "none", set<int>{}, 0 // no output so no sampling
		);
		const set<int> genePositions = {0, 1, genomeSz - 1};
		const bool isDiploid{ true };
		SpeciesTrait* spTr = createTestSpTrait(genePositions, isDiploid);
		pSpecies->addTrait(TraitType::E_D0, *spTr);
		
		Individual indMother = Individual(pCell, pPatch, 0, 0, 0, 0.0, false, 0);
		Individual indFather = Individual(pCell, pPatch, 0, 0, 0, 1.0, false, 0);

		DispersalTrait dispTrParent(spTr); // initialisation constructor
		// All genes heterozygote with same alleles for each gene
		const float valAlleleA(0.0), valAlleleB(1.0);
		const float domCoef(0.0); // no dominance for dispersal traits
		vector<shared_ptr<Allele>> gene(2);
		gene[0] = make_shared<Allele>(valAlleleA, domCoef);
		gene[1] = make_shared<Allele>(valAlleleB, domCoef);

		map<int, vector<shared_ptr<Allele>>> genotype;
		for (int i : genePositions) {
			genotype.emplace(i, gene);
		}
		dispTrParent.overwriteGenes(genotype);
		indMother.insertIndDispTrait(TraitType::E_D0, dispTrParent);
		indFather.insertIndDispTrait(TraitType::E_D0, dispTrParent);

		Individual indChild = Individual(pCell, pPatch, 0, 0, 0, 0.0, false, 0);
		indChild.inheritTraits(pSpecies, &indMother, &indFather, 1.0);
	}
}

#endif //RSDEBUG