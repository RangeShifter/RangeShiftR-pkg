
#include "NeutralStatsManager.h"
#include "Population.h"

/*----------------------------------------------------------------------------
 *
 *	Copyright (C) 2020 Greta Bocedi, Stephen C.F. Palmer, Justin M.J. Travis, Anne-Kathleen Malchow, Damaris Zurell
 *
 *	This file is part of RangeShifter.
 *
 *	RangeShifter is free software: you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *	the Free Software Foundation, either version 3 of the License, or
 *	(at your option) any later version.
 *
 *	RangeShifter is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *	GNU General Public License for more details.
 *
 *	You should have received a copy of the GNU General Public License
 *	along with RangeShifter. If not, see <https://www.gnu.org/licenses/>.
 *
 * File Created by Roslyn Henry March 2023. Code adapted from NEMO (https://nemo2.sourceforge.io/)
 --------------------------------------------------------------------------*/

 // ----------------------------------------------------------------------------------------
 // Constructor
 // ----------------------------------------------------------------------------------------
NeutralStatsManager::NeutralStatsManager(const int& nbSampledPatches, const int nLoci) {
	this->pairwiseFstMatrix = PatchMatrix(nbSampledPatches, nbSampledPatches);
	commNeutralCountTables.reserve(nLoci); //don't have to be pointers, not shared or moved
}

// ----------------------------------------------------------------------------------------
// Populate population and community-level NEUTRAL count tables
// Update allele occurrence and heterozygosity counts, and allele frequencies
// ----------------------------------------------------------------------------------------
void NeutralStatsManager::updateAllNeutralTables(Species* pSpecies, Landscape* pLandscape, set<int> const& patchList) {

	const int nLoci = pSpecies->getNPositionsForTrait(NEUTRAL);
	const int nAlleles = pSpecies->getSpTrait(NEUTRAL)->getNbNeutralAlleles();
	const int ploidy = pSpecies->isDiploid() ? 2 : 1;

	// Create / Update community-level NEUTRAL counts table
	if (!commNeutralCountTables.empty()) {
		resetCommNeutralTables();
	}
	else { // populate the tables with default values
		for (int thisLocus = 0; thisLocus < nLoci; thisLocus++) {
			NeutralCountsTable newNeutralTbl = NeutralCountsTable(nAlleles);
			commNeutralCountTables.push_back(newNeutralTbl);
		}
	}

	int nbSampledInds = 0;
	int patchAlleleCount;

	// Update counts for each population
	for (int patchId : patchList) {
		const auto patch = pLandscape->findPatch(patchId);
		const auto pPop = (Population*)patch->getPopn((intptr)pSpecies);
		if (pPop != 0) {
			// Update this population's NEUTRAL counts tables
			pPop->updatePopNeutralTables();
			nbSampledInds += pPop->sampleSize();
		}
		// Add population-level counts to community-level counts 
		for (int thisLocus = 0; thisLocus < nLoci; thisLocus++) {
			for (int allele = 0; allele < nAlleles; allele++) {

				if (pPop != 0) {
					patchAlleleCount = pPop->getAlleleTally(thisLocus, allele);
				}
				else {
					patchAlleleCount = 0;
				}
				commNeutralCountTables[thisLocus].incrementTallyBy(patchAlleleCount, allele);
			}
		}
	}

	// Update community-level frequencies
	std::for_each(commNeutralCountTables.begin(),
		commNeutralCountTables.end(),
		[&](NeutralCountsTable& v) -> void {
			v.setFrequencies(nbSampledInds * ploidy);
		});
}

// ----------------------------------------------------------------------------------------
// Reset allele tables in NeutralTable structs
// ----------------------------------------------------------------------------------------
void NeutralStatsManager::resetCommNeutralTables() {
	for (auto& entry : commNeutralCountTables) {
		entry.reset();
	}
}

// ----------------------------------------------------------------------------------------
//	Calculate allelic diversity metrics
// ----------------------------------------------------------------------------------------
void NeutralStatsManager::calcAllelicDiversityMetrics(set<int> const& patchList, const int nInds, Species* pSpecies, Landscape* pLandscape)
{
	int i, j;
	const int nLoci = pSpecies->getNPositionsForTrait(NEUTRAL);
	const int nAlleles = pSpecies->getSpTrait(NEUTRAL)->getNbNeutralAlleles();
	const int ploidy = pSpecies->isDiploid() ? 2 : 1;
	unsigned int nbPopulatedPatches = 0;
	int nbAllelesInPatch = 0;
	double meanAllelicDivInPatch = 0;
	bool alleleExistsInPop = 0;

	bool** alleleExistsInCommTable;
	alleleExistsInCommTable = new bool* [nLoci];
	for (i = 0; i < nLoci; ++i) {
		alleleExistsInCommTable[i] = new bool[nAlleles];
		for (j = 0; j < nAlleles; ++j)
			alleleExistsInCommTable[i][j] = 0;
	}

	// Compute mean nb alleles per locus per patch
	for (int patchId : patchList) {
		const auto patch = pLandscape->findPatch(patchId);
		const auto pPop = (Population*)patch->getPopn((intptr)pSpecies);
		if (pPop != 0) {
			if (pPop->sampleSize() > 0) {
				nbPopulatedPatches++;
				nbAllelesInPatch = 0;
				for (i = 0; i < nLoci; ++i)
					for (j = 0; j < nAlleles; ++j) {
						alleleExistsInPop = pPop->getAlleleTally(i, j) != 0;
						nbAllelesInPatch += alleleExistsInPop;
						alleleExistsInCommTable[i][j] |= alleleExistsInPop; // OR operator
					}
				// add mean nb of alleles per locus for Patch k to the pop mean
				meanAllelicDivInPatch += static_cast<double>(nbAllelesInPatch) / nLoci;
			}
		}
	}
	meanNbAllelesPerLocusPerPatch = nbPopulatedPatches > 0 ? meanAllelicDivInPatch / nbPopulatedPatches : 0;

	// Compute mean nb alleles per locus
	meanNbAllelesPerLocus = 0;
	for (i = 0; i < nLoci; ++i)
		for (j = 0; j < nAlleles; ++j)
			meanNbAllelesPerLocus += alleleExistsInCommTable[i][j];
	meanNbAllelesPerLocus /= nLoci;
	// Clear table 
	for (i = 0; i < nLoci; ++i)
		delete[] alleleExistsInCommTable[i];
	delete[] alleleExistsInCommTable;

	// Compute number of fixed loci per patch
	// mean number of loci that are fixed at pop level per pop
	meanNbFixedAllelesPerPatch = 0;
	if (nbPopulatedPatches > 0) {
		for (int patchId : patchList) {
			const auto patch = pLandscape->findPatch(patchId);
			const auto pPop = (Population*)patch->getPopn((intptr)pSpecies);
			if (pPop != 0) {
				for (i = 0; i < nLoci; ++i)
					for (j = 0; j < nAlleles; ++j)
						meanNbFixedAllelesPerPatch += pPop->getAlleleFrequency(i, j) == 1;
			}
		}
		meanNbFixedAllelesPerPatch /= nbPopulatedPatches;
	}

	// Compute number of fixed loci
	nbGloballyFixedAlleles = 0;
	for (i = 0; i < nLoci; ++i)
		for (j = 0; j < nAlleles; ++j)
			nbGloballyFixedAlleles += commNeutralCountTables[i].getFrequency(j) == 1;
}

// ----------------------------------------------------------------------------------------
// Calculate Ho per Nei and Chesser
// Average (observed) heterozygosity per individual
// Sum (nb of heterozygote loci) across individuals / nb individuals / nb loci
// ----------------------------------------------------------------------------------------
void NeutralStatsManager::calculateHo(set<int> const& patchList, const int nbInds, const int nbrLoci, Species* pSpecies, Landscape* pLandscape) {

	int nbHetero = 0;

	if (nbInds != 0 && pSpecies->isDiploid()) {
		for (int patchId : patchList) {
			const auto patch = pLandscape->findPatch(patchId);
			const auto pPop = (Population*)patch->getPopn((intptr)pSpecies);
			if (pPop != 0) nbHetero += pPop->countHeterozygoteLoci();
		}
		ho = static_cast<double>(nbHetero) / (nbInds * nbrLoci);
	}
	else ho = 0.0;
}

// ----------------------------------------------------------------------------------------
// Calculate Hs per Nei and Chesser
// Average expected population-level heterozygosity per locus per population
// currently not used but may be useful
// ----------------------------------------------------------------------------------------
void NeutralStatsManager::calculateHs(set<int> const& patchList, const int nbrLoci, Species* pSpecies, Landscape* pLandscape) {

	double hs = 0;
	int nPatches = 0;

	for (int patchId : patchList) {
		const auto patch = pLandscape->findPatch(patchId);
		const auto pPop = (Population*)patch->getPopn((intptr)pSpecies);
		if (pPop->sampleSize() > 0) {
			nPatches++;
			hs += pPop->computeHs();
		}
	}
	hs = (nPatches != 0 ? hs / (nbrLoci * nPatches) : 0.0);
}

// ----------------------------------------------------------------------------------------
// Calculate Ht per Nei and Chesser
// Average expected community-level heterozygosity per locus
// Currently not used but may be useful
// ----------------------------------------------------------------------------------------
void NeutralStatsManager::calculateHt(Species* pSpecies, Landscape* pLandscape, const int nLoci, const int nAlleles) {

	double ht = 0;
	int nPatches = 0;
	vector<double>locihet(nLoci, 1);
	double freq;

	for (int thisLocus = 0; thisLocus < nLoci; ++thisLocus) {
		for (int allele = 0; allele < nAlleles; ++allele) {
			freq = commNeutralCountTables[thisLocus].getFrequency(allele);
			freq *= freq; //squared frequencies
			locihet[thisLocus] -= freq;  //1 - sum of p2 = expected heterozygosity
		}
		ht += locihet[thisLocus];
	}
	ht = ht / nLoci;
}

// ----------------------------------------------------------------------------------------
// Calculate Ho per locus as per Nei and Chesser
// Observed proportion of heterozygote individuals for each locus
// Sum (nb of heterozygote individuals) / nb individuals for each locus
// ----------------------------------------------------------------------------------------
void NeutralStatsManager::calculatePerLocusHo(set<int> const& patchList, const int nbInds, const int nbrLoci, Species* pSpecies, Landscape* pLandscape) {

	vector<int> nbHeterosInComm(nbrLoci, 0);
	vector<int> nbHeterosInPop(nbrLoci);

	if (pSpecies->isDiploid()) {
		for (int patchId : patchList) {
			const auto patch = pLandscape->findPatch(patchId);
			const auto pPop = (Population*)patch->getPopn((intptr)pSpecies);
			if (pPop != 0) {
				if (pPop->sampleSize() > 0) {
					nbHeterosInPop = pPop->countNbHeterozygotesEachLocus();
					// Add counts to community total
					transform(nbHeterosInComm.begin(), nbHeterosInComm.end(), nbHeterosInPop.begin(),
						nbHeterosInComm.begin(), plus<int>());
				}
			}
		}
	}

	perLocusHo = vector<double>(nbrLoci, 0);
	if (nbInds != 0) {
		for (int i = 0; i < nbHeterosInComm.size(); i++) {
			perLocusHo[i] = static_cast<double>(nbHeterosInComm[i]) / nbInds;
		}
	}
}

// ----------------------------------------------------------------------------------------
// Fstat Weir & Cockerham
// ----------------------------------------------------------------------------------------
void NeutralStatsManager::calculateFstatWC(set<int> const& patchList, const int nbSampledIndsInComm, const int nLoci, const int nAlleles, Species* pSpecies, Landscape* pLandscape) {

	double inverseNtotal;
	double sumWeights = 0;
	double nBar, nC, inverseNbar;
	unsigned int nbPops = 0;
	const int ploidy = pSpecies->isDiploid() ? 2 : 1;
	const int globalSampleSize = nbSampledIndsInComm * ploidy; // total nb of alleles

	for (int patchId : patchList) {
		const auto patch = pLandscape->findPatch(patchId);
		const auto pPop = (Population*)patch->getPopn((intptr)pSpecies);
		if (pPop != 0) {
			int popSampleSize = pPop->sampleSize() * ploidy;
			if (popSampleSize > 0) {
				nbPops++;
				sumWeights += static_cast<double>(popSampleSize * popSampleSize) / globalSampleSize;
			}
		}
	}

	nbExtantPops = nbPops;
	totalNbSampledInds = nbSampledIndsInComm;

	if (nbPops > 1) {

		// Calculate F stats
		nBar = static_cast<double>(globalSampleSize) / nbPops; // average sample size, cannot be less than 1
		nC = (globalSampleSize - sumWeights) / nbPops - 1;
		double nBarMinusOne = (nBar == 1.0) ? 1.0 : nBar - 1.0; // avoid / 0 if exactly 1 ind per pop
		inverseNbar = 1.0 / nBarMinusOne;
		inverseNtotal = 1.0 / globalSampleSize;

		double var;
		double s2, pBar, hBar;
		double s2Denom = 1.0 / ((nbPops - 1) * nBar);
		double rTerm = static_cast<double> (nbPops - 1) / nbPops;
		double hBarFactor = (2 * nBar - 1) / (4 * nBar);

		double a = 0, b = 0, c = 0, intermediateTerm;
		for (int thisLocus = 0; thisLocus < nLoci; ++thisLocus) {
			for (int allele = 0; allele < nAlleles; ++allele) {

				s2 = hBar = 0;
				pBar = commNeutralCountTables[thisLocus].getFrequency(allele);

				for (int patchId : patchList) {
					const auto patch = pLandscape->findPatch(patchId);
					const auto pPop = (Population*)patch->getPopn((intptr)pSpecies);
					if (pPop != 0) {
						var = pPop->getAlleleFrequency(thisLocus, allele) - pBar;
						var *= var;
						s2 += var * pPop->sampleSize() * ploidy;
						hBar += pPop->getHeteroTally(thisLocus, allele); // n_i * h_i
					}
				} //end for pop

				s2 *= s2Denom;
				hBar *= inverseNtotal; // / (r * n_bar)

				intermediateTerm = pBar * (1 - pBar) - rTerm * s2;
				a += s2 - inverseNbar * (intermediateTerm - 0.25 * hBar);
				b += intermediateTerm - hBarFactor * hBar;
				c += hBar;
			} // end for allele 
		} // end for locus

		a *= nBar / nC;
		b *= nBar / nBarMinusOne;
		c *= 0.5;

		fst = a / (a + b + c); // theta hat in eq. 1 in WC 1984
		fis = b / (b + c); // f hat
		fit = (a + b) / (a + b + c); // F hat
	}
	else { // zero or one sampled pops, cannot compute F stats
		fst = 0.0;
		fis = 0.0;
		fit = 0.0;
	}
}

// ----------------------------------------------------------------------------------------
// Fstat Weir & Cockerham using Mean square approach. Similar to implementation in Hierfstat
// ----------------------------------------------------------------------------------------
void NeutralStatsManager::calcPerLocusMeanSquaresFst(set<int> const& patchList, const int nInds, const int nLoci, const int maxNbAllelesPerLocus, Species* pSpecies, Landscape* pLandscape) {

	double sumWeights = 0;
	unsigned int nbExtantPops = 0;
	const int ploidy = pSpecies->isDiploid() ? 2 : 1;
	const int totSampleSize = nInds * ploidy;

	for (int patchId : patchList) {
		const auto patch = pLandscape->findPatch(patchId);
		const auto pPop = (Population*)patch->getPopn((intptr)pSpecies);
		if (pPop != 0) {
			int ni = pPop->sampleSize() * ploidy;
			if (ni > 0) {
				nbExtantPops++;
				sumWeights += static_cast<double>(ni * ni) / totSampleSize;
			}
		}
	}

	// per locus stats, resize should only happen in first timestep of calculation:
	if (perLocusFst.size() == 0)
		perLocusFst.resize(nLoci);
	if (perLocusFis.size() == 0)
		perLocusFis.resize(nLoci);
	if (perLocusFit.size() == 0)
		perLocusFit.resize(nLoci);

	if (nbExtantPops > 1) {
		vector<int> nbAllelesEachLocus(nLoci);
		bool** alleleExistsMatrix = new bool* [nLoci];
		for (int i = 0; i < nLoci; ++i)
			alleleExistsMatrix[i] = new bool[maxNbAllelesPerLocus];

		int nbAllelesInComm = 0;
		for (int locus = 0; locus < nLoci; ++locus) {
			nbAllelesEachLocus[locus] = 0;
			for (int allele = 0; allele < maxNbAllelesPerLocus; ++allele) {
				int count = 0;
				for (int patchId : patchList) {
					const auto patch = pLandscape->findPatch(patchId);
					const auto pPop = (Population*)patch->getPopn((intptr)pSpecies);
					if (pPop != 0) count += pPop->getAlleleTally(locus, allele);
				}
				alleleExistsMatrix[locus][allele] = count != 0;
				nbAllelesEachLocus[locus] += count != 0;
			}
			nbAllelesInComm += nbAllelesEachLocus[locus];
		}

		// n, and nal are given by pop_sizes, same num ind typed at all loci in each patch
		// nc is the same for each locus
		// nt is given by tot_size, same tot num of ind typed for all loci

		//SSG: het/2 for each allele
		vector<double> SSG(nbAllelesInComm);
		vector<double> SSP(nbAllelesInComm);
		vector<double> SSi(nbAllelesInComm);

		int totalAlleleCounter = 0;
		double het, pi, var, pBar;
		int popSize;

		for (int locus = 0; locus < nLoci; ++locus) {
			for (int allele = 0; allele < maxNbAllelesPerLocus && totalAlleleCounter < nbAllelesInComm; ++allele) {

				if (alleleExistsMatrix[locus][allele] == false) continue; //do not consider alleles not present in the pop
				SSG[totalAlleleCounter] = 0;
				SSi[totalAlleleCounter] = 0;
				SSP[totalAlleleCounter] = 0;

				for (int patchId : patchList) {
					const auto patch = pLandscape->findPatch(patchId);
					const auto pPop = (Population*)patch->getPopn((intptr)pSpecies);
					if (pPop == 0) popSize = 0;
					else popSize = pPop->sampleSize();
					if (popSize == 0) continue; // skip empty patches

					het = pPop->getHeteroTally(locus, allele); // ni * h_i
					pi = pPop->getAlleleFrequency(locus, allele);
					pBar = commNeutralCountTables[locus].getFrequency(allele);
					var = pi - pBar; //(p_liu - pbar_u)^2
					var *= var;

					SSG[totalAlleleCounter] += het; // numerator MSG
					SSi[totalAlleleCounter] += 2 * popSize * pi * (1 - pi) - het / 2; // numerator
					SSP[totalAlleleCounter] += 2 * popSize * var; // numerator MSP
				}
				totalAlleleCounter++;
			}
		}

		if (totalAlleleCounter != nbAllelesInComm)
			throw runtime_error("Error:: allele counter and total number of alleles differ in WC mean squared Fstat calculation \n");

		vector<double> MSG(nbAllelesInComm);
		vector<double> MSP(nbAllelesInComm);
		vector<double> MSI(nbAllelesInComm);
		vector<double> sigw(nbAllelesInComm);
		vector<double> siga(nbAllelesInComm);
		vector<double> sigb(nbAllelesInComm);

		double SIGA = 0, SIGB = 0, SIGW = 0;

		if (nbAllelesInComm != nLoci) { // more than one allele per locus
			double nc = (totSampleSize - sumWeights) / (nbExtantPops - 1);
			int MSiDenom = nInds == nbExtantPops ? 1 : nInds - nbExtantPops; // avoid /0 if exactly 1 ind per pop

			for (int i = 0; i < nbAllelesInComm; ++i) {

				MSG[i] = SSG[i] / (2 * nInds);
				MSP[i] = SSP[i] / (nbExtantPops - 1);
				MSI[i] = SSi[i] / MSiDenom;

				siga[i] = (MSP[i] - MSI[i]) / (2 * nc);
				sigb[i] = 0.5 * (MSI[i] - MSG[i]);
				sigw[i] = MSG[i];

				SIGA += siga[i];
				SIGB += sigb[i];
				SIGW += sigw[i];
			}

			double locusSIGA, locusSIGB, locusSIGW;
			int alleleCounter;
			for (int locus = 0; locus < nLoci; ++locus) {
				alleleCounter = 0;
				locusSIGA = locusSIGB = locusSIGW = 0;

				for (int allele = 0; allele < nbAllelesEachLocus[locus]; ++allele) {
					locusSIGA += siga[alleleCounter];
					locusSIGB += sigb[alleleCounter];
					locusSIGW += sigw[alleleCounter];
					alleleCounter++;
				}
				perLocusFst[locus] = locusSIGA / (locusSIGA + locusSIGB + locusSIGW);
				perLocusFis[locus] = locusSIGB / (locusSIGB + locusSIGW);
				perLocusFit[locus] = (locusSIGA + locusSIGB) / (locusSIGA + locusSIGB + locusSIGW);
			}

			// Total F-stats
			fst = SIGA / (SIGA + SIGB + SIGW);
			fis = (SIGA + SIGB) / (SIGA + SIGB + SIGW);
			fit = SIGB / (SIGB + SIGW);
		}
		else { // no variation: only 1 allele (wildtype) at each locus 
			// so don't calculate to avoid division by zero
			fst = 0;
			fis = 0;
			fit = 0;
		}

		// Deallocate matrix
		for (int i = 0; i < nLoci; ++i)
			delete[]alleleExistsMatrix[i];
		delete[]alleleExistsMatrix;
	}
	else { // zero or one sampled pops, cannot calculate Fst
		for (int locus = 0; locus < nLoci; ++locus) {
			perLocusFst[locus] = 0.0;
			perLocusFis[locus] = 0.0;
			perLocusFit[locus] = 0.0;
		}
		fst = 0;
		fis = 0;
		fit = 0;
	}
}

// ----------------------------------------------------------------------------------------
// Patch pairwise Fst 
// Computes the weighted within and between patch Fst's as well as the overall Fst (Theta).
// The method used here is that of Weir& Hill 2002, Ann.Rev.Genet. 36:721 - 750.
// The weighting is done for samples(patches) of unequal sizes.
// ----------------------------------------------------------------------------------------
void NeutralStatsManager::calcPairwiseWeightedFst(set<int> const& patchList, const int nInds, const int nLoci, Species* pSpecies, Landscape* pLandscape) {

	const int nAlleles = (int)pSpecies->getSpTrait(NEUTRAL)->getNbNeutralAlleles();
	const int ploidy = pSpecies->isDiploid() ? 2 : 1;

	// Needs to be in vector to iterate over, copy preserves order
	vector<int> patchVect;
	copy(patchList.begin(), patchList.end(), std::back_inserter(patchVect));

	int nPatches = static_cast<int>(patchList.size());
	int nbPopulatedPatches = 0;

	// Initialise 
	if (pairwiseFstMatrix.getNbCells() != nPatches * nPatches)
		pairwiseFstMatrix = PatchMatrix(nPatches, nPatches);

	// Reset table
	pairwiseFstMatrix.setAll(0.0); // or nanf("NULL")?

	//init
	vector<double> popWeights(nPatches);
	vector<double> popSizes(nPatches);
	double** numeratorPairwiseFst = new double* [nPatches];
	for (int i = 0; i < nPatches; i++) numeratorPairwiseFst[i] = new double[nPatches];
	double totSize;
	double numeratorWeightedFst = 0;
	double denominator = 0;
	double sumWeights = 0;

	totSize = nInds * ploidy;

	// Calculate weight (n_ic) terms
	for (int i = 0; i < nPatches; ++i) {
		const auto patch = pLandscape->findPatch(patchVect[i]);
		const auto pPop = (Population*)patch->getPopn((intptr)pSpecies);
		if (pPop != 0) {
			popSizes[i] = pPop->sampleSize() * ploidy;
		} // else popSizes[i] remain default init value 0, safe
		popWeights[i] = popSizes[i] - (popSizes[i] * popSizes[i] / totSize); // n_ic in Weir & Hill 2002
		sumWeights += popWeights[i];
		if (popSizes[i] > 0) nbPopulatedPatches++;

		// Fill the pairwise Fst matrix with default value 0
		for (int j = 0; j < nPatches; j++)
			numeratorPairwiseFst[i][j] = 0;
	}

	if (nbPopulatedPatches > 1) {
		// Calculate Fst numerators and denominators
		double p, pq, pBar, sqDist, num;
		for (int i = 0; i < nPatches; ++i) {
			if (popSizes[i] == 0) continue;
			const auto patch = pLandscape->findPatch(patchVect[i]);
			const auto pPop = (Population*)patch->getPopn((intptr)pSpecies);

			for (int l = 0; l < nLoci; ++l) {
				for (int u = 0; u < nAlleles; ++u) {
					p = pPop->getAlleleFrequency(l, u); //p_liu
					pq = p * (1 - p);
					pBar = commNeutralCountTables[l].getFrequency(u);
					sqDist = p - pBar; //(p_liu - pbar_u)^2 
					sqDist *= sqDist;

					num = pq * popSizes[i] / (popSizes[i] - 1); // eq. 8 Weir & Hill 2002
					numeratorPairwiseFst[i][i] += num;
					numeratorWeightedFst += num * popSizes[i]; // see equ. 9, Weir & Hill 2002
					denominator += popSizes[i] * sqDist + popWeights[i] * pq; //common denominator

				} // end for allele
			} // end for locus
		} // end for pop

		// Diagonals
		double pairwiseFst;
		for (int i = 0; i < nPatches; ++i) {
			if (popSizes[i] == 0) continue;
			else if (denominator != 0)
			{
				pairwiseFst = 1 - (numeratorPairwiseFst[i][i] * sumWeights / denominator);
				pairwiseFstMatrix.set(i, i, pairwiseFst);
			}
			// else remain 0
		}

		// Add allele frequencies to numerators
		double pi, pj;
		for (int l = 0; l < nLoci; ++l)
			for (int u = 0; u < nAlleles; ++u)
				for (int i = 0; i < nPatches - 1; ++i) { // nPatches-1 bc bottom row not filled
					if (popSizes[i] == 0) continue;
					const auto patch = pLandscape->findPatch(patchVect[i]);
					const auto pPopI = (Population*)patch->getPopn((intptr)pSpecies);

					for (int j = i + 1; j < nPatches; ++j) { // fill only upper half of matrix
						if (popSizes[j] == 0) continue;
						const auto patch = pLandscape->findPatch(patchVect[j]);
						const auto pPopJ = (Population*)patch->getPopn((intptr)pSpecies);

						pi = pPopI->getAlleleFrequency(l, u);
						pj = pPopJ->getAlleleFrequency(l, u);
						numeratorPairwiseFst[i][j] += pi * (1 - pj) + pj * (1 - pi); // equ. 7 of Weir & Hill 2002
					}
				}

		// Final estimates of pairwise Fst (beta_ii' in eq. 7 in WC 2002)
		for (int i = 0; i < nPatches - 1; ++i) {
			if (popSizes[i] == 0) continue; // Fst for this pair remains NULL
			for (int j = i + 1; j < nPatches; ++j) {
				if (popSizes[j] == 0) continue;
				else if (denominator != 0) {
					pairwiseFst = 1 - (numeratorPairwiseFst[i][j] * sumWeights) / (2 * denominator);
					pairwiseFstMatrix.set(i, j, pairwiseFst);
				}
				// else remain 0
			}
		}

		// Estimator of global Fst weighted by sample sizes (beta_W in eq. 9 in WH 2002)
		if (denominator != 0) {
			weightedFst = 1 - (numeratorWeightedFst * sumWeights) / (denominator * totSize); // beta_w in Eq. 9 in WH 2002
		}
		else {
			weightedFst = 0.0;
		}

		// Deallocate pairwise Fst matrix
		for (int i = 0; i < nPatches; i++) delete[] numeratorPairwiseFst[i];
		delete[] numeratorPairwiseFst;
	}
	else { // zero or one pop, cannot calculate Fst
		// pairwiseFstMatrix keeps default values (0)
		weightedFst = 0.0;
	}
}

