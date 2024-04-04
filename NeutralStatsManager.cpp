
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
 // Cstor
 // ----------------------------------------------------------------------------------------

NeutralStatsManager::NeutralStatsManager(const int& nbSampledPatches, const int nLoci) {
	this->pairwiseFstMatrix = PatchMatrix(nbSampledPatches, nbSampledPatches);
	commSNPTables.reserve(nLoci); //don't have to be pointers, not shared or moved
}

// ----------------------------------------------------------------------------------------
// Set allele tables in SNPtable structs
// ----------------------------------------------------------------------------------------

void NeutralStatsManager::updateAllSNPTables(Species* pSpecies, Landscape* pLandscape, set<int> const& patchList) {

	const int nLoci = pSpecies->getNPositionsForTrait(SNP);
	const int nAlleles = (int)pSpecies->getSpTrait(SNP)->getMutationParameters().find(MAX)->second;
	const int ploidy = (pSpecies->isDiploid() ? 2 : 1);

	if (!commSNPTables.empty()) {
		resetCommSNPtables();
	}
	else { // populate the tables with default values
		for (int thisLocus = 0; thisLocus < nLoci; thisLocus++) {
			SNPtable newSNPtbl = SNPtable(nAlleles);
			commSNPTables.push_back(newSNPtbl);
		}
	}

	int populationSize = 0;
	int patchAlleleCount;

	for (int patchId : patchList) {
		const auto patch = pLandscape->findPatch(patchId);
		const auto pPop = (Population*)patch->getPopn((intptr)pSpecies);
		if (pPop != 0) {
			// Update this population's SNP counts tables
			pPop->updatePopSNPtables();
			populationSize += pPop->sampleSize();
		}
		// Update global SNP counts tables
		for (int thisLocus = 0; thisLocus < nLoci; thisLocus++) {
			for (int allele = 0; allele < nAlleles; allele++) {

				if (pPop != 0) {
					patchAlleleCount = pPop->getAlleleTally(thisLocus, allele);
				}
				else {
					patchAlleleCount = 0;
				}
				commSNPTables[thisLocus].incrementTallyBy(patchAlleleCount, allele);
			}
		}
	}

	// Update global frequency
	populationSize *= ploidy;
	std::for_each(commSNPTables.begin(),
		commSNPTables.end(),
		[&](SNPtable& v) -> void {
			v.setFrequencies(populationSize);
		});
}

// ----------------------------------------------------------------------------------------
// Reset allele tables in SNPtable structs
// ----------------------------------------------------------------------------------------

void NeutralStatsManager::resetCommSNPtables() {
	for (auto& entry : commSNPTables) {
		entry.reset();
	}
}


// ----------------------------------------------------------------------------------------
// set loci diversity
// ----------------------------------------------------------------------------------------

void NeutralStatsManager::setLociDiversityCounter(set<int> const& patchList, const int nInds, Species* pSpecies, Landscape* pLandscape)
{
	int i, j;
	const int nLoci = pSpecies->getNPositionsForTrait(SNP);
	const int nAlleles = (int)pSpecies->getSpTrait(SNP)->getMutationParameters().find(MAX)->second;
	const int ploidy = (pSpecies->isDiploid() ? 2 : 1);
	unsigned int nbPopulatedPatches = 0;
	int nbAllelesInPatch = 0;
	double meanAllelicDivInPatch = 0;
	bool alleleExistsInPop = 0;

	bool** alleleExistsInCommTable;

	// number of alleles per locus, Patch and pop counters:
	alleleExistsInCommTable = new bool* [nLoci];

	for (i = 0; i < nLoci; ++i) {
		alleleExistsInCommTable[i] = new bool[nAlleles];
		for (j = 0; j < nAlleles; ++j)
			alleleExistsInCommTable[i][j] = 0;
	}

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

	// Compute mean allelic diversity per locus
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
			nbGloballyFixedAlleles += commSNPTables[i].getFrequency(j) == 1;
}

// ----------------------------------------------------------------------------------------
// calculate Ho per Nei and Chesser
// ----------------------------------------------------------------------------------------
void NeutralStatsManager::calculateHo(set<int> const& patchList, const int nbInds, const int nbrLoci, Species* pSpecies, Landscape* pLandscape) {

	int nbHetero = 0;
	int nLoci = nbInds * nbrLoci;

	if (nLoci != 0) {
		for (int patchId : patchList) {
			const auto patch = pLandscape->findPatch(patchId);
			const auto pPop = (Population*)patch->getPopn((intptr)pSpecies);
			if (pPop != 0) {
				nbHetero += pPop->countHeterozygoteLoci();
			}
		}
		_ho = static_cast<double>(nbHetero) / static_cast<double>(nLoci);
	}
	else _ho = 0.0;
}

// ----------------------------------------------------------------------------------------
// calculate Hs per Nei and Chesser, currently not used but may be useful
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

	_hs = (nPatches != 0 ? hs / (nbrLoci * nPatches) : 0.0);

}

// ----------------------------------------------------------------------------------------
// calculate Ht per Nei and Chesser, currently not used but may be useful
// ----------------------------------------------------------------------------------------


void NeutralStatsManager::calculateHt(Species* pSpecies, Landscape* pLandscape, const int nLoci, const int nAlleles) {

	double ht = 0;
	int nPatches = 0;
	vector<double>locihet(nLoci, 1);
	double freq;

	for (int thisLocus = 0; thisLocus < nLoci; ++thisLocus) {

		for (int allele = 0; allele < nAlleles; ++allele) {

			freq = commSNPTables[thisLocus].getFrequency(allele);

			freq *= freq; //squared frequencies

			locihet[thisLocus] -= freq;  //1 - sum of p2 = expected heterozygosity
		}

		ht += locihet[thisLocus];
	}

	_ht = ht / nLoci;
}

// ----------------------------------------------------------------------------------------
// calculate Ho per locus as per Nei and Chesser
// ----------------------------------------------------------------------------------------

void NeutralStatsManager::calculateHo2(set<int> const& patchList, const int nbInds, const int nbrLoci, Species* pSpecies, Landscape* pLandscape) {

	vector<double> hetero(nbrLoci, 0);
	double nLoci = nbInds * nbrLoci;

	for (int patchId : patchList) {
		const auto patch = pLandscape->findPatch(patchId);
		const auto pPop = (Population*)patch->getPopn((intptr)pSpecies);
		if (pPop != 0) {
			if (pPop->sampleSize() > 0) {
				const vector<double> heteroPatch = pPop->countLociHeterozyotes();
				transform(hetero.begin(), hetero.end(), heteroPatch.begin(),
					hetero.begin(), plus<double>());
			}
		}
	}

	if (nbInds != 0)
		for (double h : hetero) {
			h /= nbInds;
		}
	perLocusHo = hetero;
}


// ----------------------------------------------------------------------------------------
// Fstat Weir & Cockerham
// ----------------------------------------------------------------------------------------


void NeutralStatsManager::calculateFstatWC(set<int> const& patchList, const int nbSampledIndsInComm, const int nLoci, const int nAlleles, Species* pSpecies, Landscape* pLandscape) {

	double inverseNtotal;
	double sumWeights = 0;
	double nBar, nC, inverseNbar;
	unsigned int nbPopulatedPatches = 0;

	for (int patchId : patchList) {
		const auto patch = pLandscape->findPatch(patchId);
		const auto pPop = (Population*)patch->getPopn((intptr)pSpecies);
		if (pPop != 0) {
			int nbSampledIndsinPop = pPop->sampleSize();
			if (nbSampledIndsinPop > 0) {
				nbPopulatedPatches++;
				sumWeights += static_cast<double>(nbSampledIndsinPop * nbSampledIndsinPop) / nbSampledIndsInComm;
			}
		}
	}

	nbExtantPops = nbPopulatedPatches;
	totalNbSampledInds = nbSampledIndsInComm; // r * nBar

	if (nbPopulatedPatches > 1) {

		// Calculate F stats
		nBar = static_cast<double>(nbSampledIndsInComm) / nbPopulatedPatches; // average sample size, cannot be less than 1
		nC = (nbSampledIndsInComm - sumWeights) / nbPopulatedPatches - 1;
		double nBarMinusOne = (nBar == 1.0) ? 1.0 : nBar - 1.0; // avoid / 0 if exactly 1 ind per pop
		inverseNbar = 1.0 / nBarMinusOne;
		inverseNtotal = 1.0 / nbSampledIndsInComm;

		double var;
		double s2, pBar, hBar;
		double s2Denom = 1.0 / ((nbPopulatedPatches - 1) * nBar);
		double rTerm = static_cast<double>(nbPopulatedPatches - 1) / nbPopulatedPatches;
		double hBarFactor = (2 * nBarMinusOne) / (4 * nBar);
		double a = 0, b = 0, c = 0, intermediateTerm;

		for (int thisLocus = 0; thisLocus < nLoci; ++thisLocus) {
			for (int allele = 0; allele < nAlleles; ++allele) {

				s2 = hBar = 0;
				pBar = commSNPTables[thisLocus].getFrequency(allele);

				for (int patchId : patchList) {
					const auto patch = pLandscape->findPatch(patchId);
					const auto pPop = (Population*)patch->getPopn((intptr)pSpecies);
					if (pPop != 0) {
						var = pPop->getAlleleFrequency(thisLocus, allele) - pBar;
						var *= var;
						s2 += var * pPop->sampleSize();
						hBar += pPop->getHeteroTally(thisLocus, allele);
					}
				} //end for pop

				s2 *= s2Denom;
				hBar *= inverseNtotal;

				intermediateTerm = pBar * (1 - pBar) - rTerm * s2;
				a += s2 - inverseNbar * (intermediateTerm - 0.25 * hBar);
				b += intermediateTerm - hBarFactor * hBar;
				c += hBar;
			} // end for allele 
		} // end for locus

		a *= nBar / nC;
		b *= nBar / nBarMinusOne;
		c *= 0.5;

		FstWeirCockerham = a / (a + b + c); // theta hat in eq. 1 in WC 1984
		FitWeirCockerham = (a + b) / (a + b + c); // F hat
		FisWeirCockerham = b / (b + c); // f hat
	}
	else { // zero or one sampled pops, cannot compute F stats
		FstWeirCockerham = 0.0;
		FitWeirCockerham = 0.0;
		FisWeirCockerham = 0.0;
	}
}


// ----------------------------------------------------------------------------------------
// Fstat Weir & Cockerham using Mean square approach. Similar to implementation in Hierfstat
// ----------------------------------------------------------------------------------------

void NeutralStatsManager::calculateFstatWC_MS(set<int> const& patchList, const int nInds, const int nLoci, const int maxNbAllelesPerLocus, Species* pSpecies, Landscape* pLandscape) {

	double sumWeights = 0;
	unsigned int nbExtantPops = 0;

	for (int patchId : patchList) {
		const auto patch = pLandscape->findPatch(patchId);
		const auto pPop = (Population*)patch->getPopn((intptr)pSpecies);
		if (pPop != 0) {
			int patchSize = pPop->sampleSize();
			if (patchSize > 0) {
				nbExtantPops++;
				sumWeights += static_cast<double>(patchSize * patchSize) / nInds;
			}
		}
	}

	// per locus stats, resize should only happen in first timestep of calculation:
	if (perLocusFstWeirCockerham.size() == 0)
		perLocusFstWeirCockerham.resize(nLoci);
	if (perLocusFisWeirCockerham.size() == 0)
		perLocusFisWeirCockerham.resize(nLoci);
	if (perLocusFitWeirCockerham.size() == 0)
		perLocusFitWeirCockerham.resize(nLoci);

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
					pBar = commSNPTables[locus].getFrequency(allele);
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
			double nc = (nInds - sumWeights) / (nbExtantPops - 1);
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
				perLocusFstWeirCockerham[locus] = locusSIGA / (locusSIGA + locusSIGB + locusSIGW);
				perLocusFisWeirCockerham[locus] = locusSIGB / (locusSIGB + locusSIGW);
				perLocusFitWeirCockerham[locus] = (locusSIGA + locusSIGB) / (locusSIGA + locusSIGB + locusSIGW);
			}

			// Total F-stats
			FstWeirCockerham = SIGA / (SIGA + SIGB + SIGW);
			FitWeirCockerham = (SIGA + SIGB) / (SIGA + SIGB + SIGW);
			FisWeirCockerham = SIGB / (SIGB + SIGW);
		}
		else { // no variation: only 1 allele (wildtype) at each locus 
			// so don't calculate to avoid division by zero
			FstWeirCockerham = 0;
			FitWeirCockerham = 0;
			FisWeirCockerham = 0;
		}

		// Deallocate matrix
		for (int i = 0; i < nLoci; ++i)
			delete[]alleleExistsMatrix[i];
		delete[]alleleExistsMatrix;
	}
	else { // zero or one sampled pops, cannot calculate Fst
		for (int locus = 0; locus < nLoci; ++locus) {
			perLocusFstWeirCockerham[locus] = 0.0;
			perLocusFisWeirCockerham[locus] = 0.0;
			perLocusFitWeirCockerham[locus] = 0.0;
		}
		FstWeirCockerham = 0;
		FitWeirCockerham = 0;
		FisWeirCockerham = 0;
	}
}


// ----------------------------------------------------------------------------------------
// Patch pairwise Fst 
// Computes the weighted within and between patch Fst's as well as the overall Fst (Theta).
//	The method used here is that of Weir& Hill 2002, Ann.Rev.Genet. 36:721 - 750.
// The weighting is done for samples(patches) of unequal sizes.
// ----------------------------------------------------------------------------------------

void NeutralStatsManager::setFstMatrix(set<int> const& patchList, const int nInds, const int nLoci, Species* pSpecies, Landscape* pLandscape) {

	const int nAlleles = (int)pSpecies->getSpTrait(SNP)->getMutationParameters().find(MAX)->second;

	// Needs to be in vector to iterate over, copy preserves order
	vector<int> patchVect;
	copy(patchList.begin(), patchList.end(), std::back_inserter(patchVect));

	int nPatches = static_cast<int>(patchList.size());
	int nbPopulatedPatches = 0;

	// Initialise 
	if (pairwiseFstMatrix.getNbCells() != nPatches * nPatches)
		pairwiseFstMatrix = PatchMatrix(nPatches, nPatches);

	// Reset table
	pairwiseFstMatrix.assign(0.0); // or nanf("NULL")?

	//init
	vector<double> popWeights(nPatches);
	vector<double> popSizes(nPatches);
	double** numeratorPairwiseFst = new double* [nPatches];
	for (int i = 0; i < nPatches; i++) numeratorPairwiseFst[i] = new double[nPatches];
	double totSize;
	double numeratorWeightedFst = 0;
	double denominator = 0;
	double sumWeights = 0;

	totSize = nInds * 2; // diploid

	// Calculate weight (n_ic) terms
	for (int i = 0; i < nPatches; ++i) {
		const auto patch = pLandscape->findPatch(patchVect[i]);
		const auto pPop = (Population*)patch->getPopn((intptr)pSpecies);
		if (pPop != 0) {
			popSizes[i] = pPop->sampleSize() * 2.0;
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
					pBar = commSNPTables[l].getFrequency(u);
					sqDist = p - pBar; //(p_liu - pbar_u)^2 
					sqDist *= sqDist;

					num = pq * popSizes[i] / (popSizes[i] - 1);
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
				for (int i = 0; i < nPatches - 1; ++i) {
					if (popSizes[i] == 0) continue;
					const auto patch = pLandscape->findPatch(patchVect[i]);
					const auto pPopI = (Population*)patch->getPopn((intptr)pSpecies);

					for (int j = i + 1; j < nPatches; ++j) {
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
			weightedFstWeirHill = 1 - (numeratorWeightedFst * sumWeights) / (denominator * totSize); // beta_w in Eq. 9 in WH 2002
		}
		else {
			weightedFstWeirHill = 0.0;
		}

		// Deallocate pairwise Fst matrix
		for (int i = 0; i < nPatches; i++) delete[] numeratorPairwiseFst[i];
		delete[] numeratorPairwiseFst;
	}
	else { // zero or one pop, cannot calculate Fst
		// pairwiseFstMatrix keeps default values (0)
		weightedFstWeirHill = 0.0;
	}
}





