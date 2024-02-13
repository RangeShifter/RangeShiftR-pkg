
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

NeutralStatsManager::NeutralStatsManager(set<int> const& patchList, const int nLoci) {
	this->_fst_matrix = PatchMatrix(static_cast<int>(patchList.size()), static_cast<int>(patchList.size()));
	globalAlleleTable.reserve(nLoci); //don't have to be pointers, not shared or moved
}

// ----------------------------------------------------------------------------------------
// Set allele tables in NeutralData structs
// ----------------------------------------------------------------------------------------


void NeutralStatsManager::updateAlleleTables(Species* pSpecies, Landscape* pLandscape, set<int> const& patchList) {

	const int nLoci = pSpecies->getNPositionsForTrait(SNP);
	const int nAlleles = (int)pSpecies->getSpTrait(SNP)->getMutationParameters().find(MAX)->second;
	const int chromosomes = (pSpecies->isDiploid() ? 2 : 1);

	if (!globalAlleleTable.empty())
		resetGlobalAlleleTable();

	int populationSize = 0;

	for (int patchId : patchList) {
		const auto patch = pLandscape->findPatch(patchId);
		const auto pPop = (Population*)patch->getPopn((intptr)pSpecies);
		pPop->updateAlleleTable();
		populationSize += pPop->sampleSize();

		for (int thisLocus = 0; thisLocus < nLoci; thisLocus++) {
			for (int allele = 0; allele < nAlleles; allele++) {

				int patchAlleleCount = pPop->getAlleleCount(thisLocus, allele);

				if (globalAlleleTable.size() <= thisLocus) { //if first allele of new loci (should only happen in first calculation step)
					NeutralData n = NeutralData(nAlleles, allele, patchAlleleCount);
					globalAlleleTable.push_back(n);
				}
				else globalAlleleTable[thisLocus].incrementCountBy(patchAlleleCount, allele);
			}
		}
	}

	populationSize *= chromosomes;

	std::for_each(globalAlleleTable.begin(),
		globalAlleleTable.end(),
		[&](NeutralData &v) -> void {
			v.setFrequencies(populationSize);
		});
}

// ----------------------------------------------------------------------------------------
// Reset allele tables in NeutralData structs
// ----------------------------------------------------------------------------------------

void NeutralStatsManager::resetGlobalAlleleTable() {
	for (auto& entry : globalAlleleTable) {
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
	const int chromosomes = (pSpecies->isDiploid() ? 2 : 1);
	unsigned int nbpatch = 0;
	double patch_mean, pop_mean = 0;

	bool** pop_div;

	// number of alleles per locus, Patch and pop counters:
	pop_div = new bool* [nLoci];

	for (i = 0; i < nLoci; ++i) {
		pop_div[i] = new bool[nAlleles];
		for (j = 0; j < nAlleles; ++j)
			pop_div[i][j] = 0;
	}

	for (int patchId : patchList) {
		const auto patch = pLandscape->findPatch(patchId);
		const auto pPop = (Population*)patch->getPopn((intptr)pSpecies);

		nbpatch += (pPop->sampleSize() != 0);
		if (pPop->sampleSize() > 0) {
			patch_mean = 0;
			for (i = 0; i < nLoci; ++i)
				for (j = 0; j < nAlleles; ++j) {
					patch_mean += (pPop->getAlleleCount(i, j) != 0);
					pop_div[i][j] |= (pPop->getAlleleCount(i, j) != 0); // OR
				}
			// add mean nb of alleles per locus for Patch k to the pop mean
			pop_mean += patch_mean / nLoci;
		}
	}

	_nb_alleles_local = (nbpatch ? pop_mean / nbpatch : nanf("NULL"));
	_nb_alleles_global = 0;

	for (i = 0; i < nLoci; ++i)
		for (j = 0; j < nAlleles; ++j)
			_nb_alleles_global += pop_div[i][j];

	_nb_alleles_global /= nLoci;

	for (i = 0; i < nLoci; ++i)
		delete[] pop_div[i];
	delete[] pop_div;

	//number of fixed loci, local and global counters:
	_fix_loc_local = 0;

	for (int patchId : patchList) {
		const auto patch = pLandscape->findPatch(patchId);
		const auto pPop = (Population*)patch->getPopn((intptr)pSpecies);
		for (i = 0; i < nLoci; ++i)
			for (j = 0; j < nAlleles; ++j)
				_fix_loc_local += (pPop->getAlleleFrequency(i, j) == 1);
	}

	_fix_loc_local /= nbpatch;
	_fix_loc_global = 0;

	//globally:  
	for (i = 0; i < nLoci; ++i)
		for (j = 0; j < nAlleles; ++j)
			_fix_loc_global += (globalAlleleTable[i].getFrequency(j) == 1);
}

// ----------------------------------------------------------------------------------------
// calculate Ho per Nei and Chesser
// ----------------------------------------------------------------------------------------
void NeutralStatsManager::calculateHo(set<int> const& patchList, const int nbInds, const int nbrLoci, Species* pSpecies, Landscape* pLandscape) {

	double hetero = 0;
	double nLoci = nbInds * nbrLoci;

	if (nLoci != 0) {
		for (int patchId : patchList) {
			const auto patch = pLandscape->findPatch(patchId);
			const auto pPop = (Population*)patch->getPopn((intptr)pSpecies);
			hetero += pPop->countHeterozygoteLoci();
		}
		_ho = hetero / nLoci;
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

			freq = globalAlleleTable[thisLocus].getFrequency(allele);

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
		if (pPop->sampleSize() > 0) {
			const vector<double> heteroPatch = pPop->countLociHeterozyotes();

			transform(hetero.begin(), hetero.end(), heteroPatch.begin(),
				hetero.begin(), plus<double>());
		}
	}

	if (nbInds != 0)
		for (double h : hetero) {
			h /= nbInds;
		}
	ho_loc = hetero;
}


// ----------------------------------------------------------------------------------------
// Fstat Weir & Cockerham
// ----------------------------------------------------------------------------------------


void NeutralStatsManager::calculateFstatWC(set<int> const& patchList, const int nInds, const int nLoci, const int nAlleles, Species* pSpecies, Landscape* pLandscape) {

	double inverse_n_total;
	double sum_weights = 0;
	double n_bar, n_c, inverse_n_bar;
	unsigned int extantPs = 0;

	for (int patchId : patchList) {
		const auto patch = pLandscape->findPatch(patchId);
		const auto pPop = (Population*)patch->getPopn((intptr)pSpecies);
		int patchSize = pPop->sampleSize();
		if (patchSize) {
			extantPs++;
			sum_weights += (patchSize * patchSize / static_cast<double>(nInds));
		}
	}

	_n_extantPopulations = extantPs;
	_n_individuals = nInds;

	n_bar = nInds / static_cast<double>(extantPs);
	n_c = (nInds - sum_weights) / (extantPs - 1);
	inverse_n_bar = 1.0 / (n_bar - 1);
	inverse_n_total = 1.0 / nInds;

	double var;
	double s2, p_bar, h_bar;
	double s2_denom = 1.0 / ((extantPs - 1) * n_bar),
		r = (double)(extantPs - 1) / extantPs,
		h_bar_factor = (2 * n_bar - 1) / (4 * n_bar);

	double a = 0, b = 0, c = 0, x;

	for (int thisLocus = 0; thisLocus < nLoci; ++thisLocus) {

		for (int allele = 0; allele < nAlleles; ++allele) {

			s2 = p_bar = h_bar = 0;

			for (int patchId : patchList) {
				const auto patch = pLandscape->findPatch(patchId);
				const auto pPop = (Population*)patch->getPopn((intptr)pSpecies);

				var = pPop->getAlleleFrequency(thisLocus, allele) - globalAlleleTable[thisLocus].getFrequency(allele); //(p_liu - pbar_u)^2 

				var *= var;

				s2 += var * pPop->sampleSize();

				h_bar += pPop->getHetero(thisLocus, allele);

			}//end for pop

			s2 *= s2_denom;
			p_bar = globalAlleleTable[thisLocus].getFrequency(allele);
			h_bar *= inverse_n_total;

			x = p_bar * (1 - p_bar) - r * s2;
			a += s2 - inverse_n_bar * (x - 0.25 * h_bar);
			b += x - h_bar_factor * h_bar;
			c += h_bar; 
		} // end for allele 
	} // end for locus

	a *= n_bar / n_c;
	b *= n_bar / (n_bar - 1);
	c *= 0.5;

	_fst_WC = a / (a + b + c);
	_fit_WC = (a + b) / (a + b + c);
	_fis_WC = b / (b + c);
}


// ----------------------------------------------------------------------------------------
// Fstat Weir & Cockerham using Mean square approach. Similar to implementation in Hierfstat
// ----------------------------------------------------------------------------------------

void NeutralStatsManager::calculateFstatWC_MS(set<int> const& patchList, const int nInds, const int nLoci, const int nAlleles, Species* pSpecies, Landscape* pLandscape) {

	double sum_weights = 0;
	double nc;
	unsigned int extantPs = 0;

	for (int patchId : patchList) {
		const auto patch = pLandscape->findPatch(patchId);
		const auto pPop = (Population*)patch->getPopn((intptr)pSpecies);
		int patchSize = pPop->sampleSize();
		if (patchSize) {
			extantPs++;
			sum_weights += (patchSize * patchSize / static_cast<double>(nInds));
		}

	}

	nc = (nInds - sum_weights) / (extantPs - 1);

	unsigned int npl = extantPs; //all loci typed in all patches

	//p = _alleleFreqTable
	//pb = _globalAlleleFreq

	vector<int> alploc(nLoci);

	unsigned int** alploc_table = new unsigned int* [nLoci];

	for (int i = 0; i < nLoci; ++i)
		alploc_table[i] = new unsigned int[nAlleles];

	int tot_num_allele = 0;

	for (int l = 0; l < nLoci; ++l) {

		alploc[l] = 0;

		for (int cnt, a = 0; a < nAlleles; ++a) {

			cnt = 0;

			for (int patchId : patchList) {
				const auto patch = pLandscape->findPatch(patchId);
				const auto pPop = (Population*)patch->getPopn((intptr)pSpecies);

				cnt += pPop->getAlleleCount(l, a);

			}
			alploc_table[l][a] = (cnt != 0);
			alploc[l] += (cnt != 0);
		}

		tot_num_allele += alploc[l];
	}

	//n, and nal are given by pop_sizes, same num ind typed at all loci in each patch
//nc is the same for each locus
//nt is given by tot_size, same tot num of ind typed for all loci

//SSG: het/2 for each allele
	vector<double> SSG(tot_num_allele);
	vector<double> SSP(tot_num_allele);
	vector<double> SSi(tot_num_allele);

	int all_cntr = 0;

	double het, freq, var;

	for (int l = 0; l < nLoci; ++l) {

		for (int a = 0; a < nAlleles && all_cntr < tot_num_allele; ++a) {

			if (alploc_table[l][a] == 0) continue; //do not consider alleles not present in the pop

			SSG[all_cntr] = 0;
			SSi[all_cntr] = 0;
			SSP[all_cntr] = 0;

			for (int patchId : patchList) {
				const auto patch = pLandscape->findPatch(patchId);
				const auto pPop = (Population*)patch->getPopn((intptr)pSpecies);
				int popSize = pPop->sampleSize();
				if (!popSize) continue; //skip empty patches

				het = pPop->getHetero(l, a);

				freq = pPop->getAlleleFrequency(l, a);

				var = freq - globalAlleleTable[l].getFrequency(a); //(p_liu - pbar_u)^2

				var *= var;

				SSG[all_cntr] += het;

				SSi[all_cntr] += 2 * popSize * freq * (1 - freq) - het / 2;

				SSP[all_cntr] += 2 * popSize * var;
			}
			all_cntr++;
		}
	}

	if (all_cntr != tot_num_allele)
		cout << endl << ("Error:: allele counter and total number of alleles differ in WC mean squared Fstat calculation \n");

	//these shouldn't have to be dynamically allocated manually, if using a vector from stl then would do allocation for you 

	vector<double> MSG(tot_num_allele);
	vector<double> MSP(tot_num_allele);
	vector<double> MSI(tot_num_allele);
	vector<double> sigw(tot_num_allele);
	vector<double> siga(tot_num_allele);
	vector<double> sigb(tot_num_allele);

	//	double *FST_pal = new double[tot_num_allele];
	//	double *FIS_pal = new double[tot_num_allele];

	double SIGA = 0, SIGB = 0, SIGW = 0;

	//per locus stats, resize should only happen in first timestep of calculation:
	if (_fst_WC_loc.size() == 0)
		_fst_WC_loc.resize(nLoci);
	if (_fis_WC_loc.size() == 0)
		_fis_WC_loc.resize(nLoci);
	if (_fit_WC_loc.size() == 0)
		_fit_WC_loc.resize(nLoci);
	
	if (tot_num_allele != nLoci) { 
		for (int i = 0; i < tot_num_allele; ++i) {

			MSG[i] = SSG[i] / (2 * nInds);
			sigw[i] = MSG[i]; //wasted!

			MSP[i] = SSP[i] / (npl - 1);

			MSI[i] = SSi[i] / (nInds - npl);

			sigb[i] = 0.5 * (MSI[i] - MSG[i]);

			siga[i] = (MSP[i] - MSI[i]) / (2 * nc);

			//		FST_pal[i] = siga[i]/(siga[i]+sigb[i]+sigw[i]);
			//		FIS_pal[i] = sigb[i]/(sigb[i]+sigw[i]);

			SIGA += siga[i];
			SIGB += sigb[i];
			SIGW += sigw[i];
		}
		double lsiga, lsigb, lsigw;

		//	cout<<"  computing sigma per locus\n";

		for (int allcntr = 0, i = 0; i < nLoci; ++i) {

			lsiga = lsigb = lsigw = 0;

			for (int l = 0; l < alploc[i]; ++l) {
				lsiga += siga[allcntr];
				lsigb += sigb[allcntr];
				lsigw += sigw[allcntr];
				allcntr++;
			}

			_fst_WC_loc[i] = lsiga / (lsiga + lsigb + lsigw);
			_fis_WC_loc[i] = lsigb / (lsigb + lsigw);
			_fit_WC_loc[i] = (lsiga + lsigb) / (lsiga + lsigb + lsigw);

		}

		// Total F-stats
		_fst_WC = SIGA / (SIGA + SIGB + SIGW);
		_fit_WC = (SIGA + SIGB) / (SIGA + SIGB + SIGW);
		_fis_WC = SIGB / (SIGB + SIGW);

	}
	else { //then there is no variation at any locus, only 1 allele (wildtype) at each locus so don't calculate to avoid division by zero issues
		// Total F-stats
		_fst_WC = 0;
		_fit_WC = 0;
		_fis_WC = 0;
	} 

	for (int i = 0; i < nLoci; ++i)
		delete[]alploc_table[i];
	delete[]alploc_table;
}


// ----------------------------------------------------------------------------------------
// Patch pairwise Fst 
// Computes the weighted within and between patch Fst's as well as the overall Fst (Theta).
//	The method used here is that of Weir& Hill 2002, Ann.Rev.Genet. 36:721 - 750.
// The weighting is done for samples(patches) of unequal sizes.
// ----------------------------------------------------------------------------------------

void NeutralStatsManager::setFstMatrix(set<int> const& patchList, const int nInds, const int nLoci, Species* pSpecies, Landscape* pLandscape) {

	const int nAlleles = (int)pSpecies->getSpTrait(SNP)->getMutationParameters().find(MAX)->second;

	vector<int> patchVect;

	copy(patchList.begin(), patchList.end(), std::back_inserter(patchVect)); //needs to be in vector to iterate over, copy preserves order

	int nPatches = static_cast<int>(patchList.size());

	//initialise 

	 if (_fst_matrix.length() != nPatches * nPatches)
		_fst_matrix = PatchMatrix(nPatches, nPatches);

	//reset table
	_fst_matrix.assign(nanf("NULL"));

	//init
	double* pop_weights = new double[nPatches];
	double* pop_sizes = new double[nPatches];
	double** numerator = new double* [nPatches];
	for (int i = 0; i < nPatches; i++) numerator[i] = new double[nPatches];
	double tot_size;
	double numerator_W = 0;
	double denominator = 0;
	double sum_weights = 0;

	tot_size = nInds * 2; //diploid

	for (int i = 0; i < nPatches; ++i) {

		const auto patch = pLandscape->findPatch(patchVect[i]);
		const auto pPop = (Population*)patch->getPopn((intptr)pSpecies);
		pop_sizes[i] = pPop->sampleSize() * 2;
		pop_weights[i] = pop_sizes[i] - (pop_sizes[i] * pop_sizes[i] / tot_size); //n_ic in Weir & Hill 2002
		sum_weights += pop_weights[i];
		for (int j = 0; j < nPatches; j++)
			numerator[i][j] = 0;
	}

	double p, pq, var, num;

	for (int i = 0; i < nPatches; ++i) {

		if (!pop_sizes[i]) continue;

		const auto patch = pLandscape->findPatch(patchVect[i]);
		const auto pPop = (Population*)patch->getPopn((intptr)pSpecies);

//#if RSDEBUG
//		DEBUGLOG << "OutputPairwiseStats: 644 "
//			<< " population " << patch->getPatchNum() << " popSize[i] " << pop_sizes[i] << endl;
//		//	<< endl;
//#endif

		for (int l = 0; l < nLoci; ++l) {

			for (int u = 0; u < nAlleles; ++u) {

				p = pPop->getAlleleFrequency(l, u); //p_liu

				pq = p * (1 - p);

				var = p - globalAlleleTable[l].getFrequency(u); //(p_liu - pbar_u)^2 

				var *= var;

				num = pq * pop_sizes[i] / (pop_sizes[i] - 1);

				numerator[i][i] += num;

				numerator_W += num * pop_sizes[i]; //see equ. 9, Weir & Hill 2002

				denominator += pop_sizes[i] * var + pop_weights[i] * pq; //common denominator

//#if RSDEBUG
//				DEBUGLOG << "OutputPairwiseStats: 671 "
//					<< " loci " << l << " allele " << u << " frequency " << p << " pq " << pq << " var " << var << " num " << num << 
//					" numerator[i][i] " << numerator[i][i] << " numerator_W " << numerator_W << " denominator " << denominator
//					<< endl;
//#endif


			} // end for allele
		}// end for locus
	}//end for pop

	for (int i = 0; i < nPatches; ++i) {
		if (!pop_sizes[i]) continue;
		if(denominator != 0)
			_fst_matrix.set(i, i, 1 - (numerator[i][i] * sum_weights / denominator));
		else
			_fst_matrix.set(i, i, 0.0);
//#if RSDEBUG
//		DEBUGLOG << "OutputPairwiseStats: 717 "
//			<< " result " << 1 - (numerator[i][i] * sum_weights / denominator) << " in matrix " << getPairwiseFst(i, i)
//			<< endl;
//#endif
	}
	_fst_WH = 1 - ((numerator_W * sum_weights) / (denominator * tot_size)); //equ. 9 Weir & Hill 2002

	//pairwise Fst:
	double pi, pj;
	for (int l = 0; l < nLoci; ++l)
		for (int u = 0; u < nAlleles; ++u)
			for (int i = 0; i < nPatches - 1; ++i) {
				if (!pop_sizes[i]) continue;

				const auto patch = pLandscape->findPatch(patchVect[i]);
				const auto pPopI = (Population*)patch->getPopn((intptr)pSpecies);

				for (int j = i + 1; j < nPatches; ++j) {
					if (!pop_sizes[j]) continue;
					const auto patch = pLandscape->findPatch(patchVect[j]);
					const auto pPopJ = (Population*)patch->getPopn((intptr)pSpecies);

					pi = pPopI->getAlleleFrequency(l, u);
					pj = pPopJ->getAlleleFrequency(l, u);
					numerator[i][j] += pi * (1 - pj) + pj * (1 - pi); //equ. 7 of Weir & Hill 2002
				}
			}

	for (int i = 0; i < nPatches - 1; ++i) {
		if (!pop_sizes[i]) continue;
		for (int j = i + 1; j < nPatches; ++j) {
			if (!pop_sizes[j]) continue;
			if (denominator != 0)
				_fst_matrix.set(i, j, 1 - ((numerator[i][j] * sum_weights) / (2 * denominator)));
			else
				_fst_matrix.set(i, j, 0.0);

//#if RSDEBUG
//			DEBUGLOG << "OutputPairwiseStats: 717 "
//				<< " result " << 1 - ((numerator[i][j] * sum_weights) / (2 * denominator)) << " in matrix " << getPairwiseFst(i, j)
//				<< endl;
//#endif

		}
	} 
	delete[] pop_weights;
	delete[] pop_sizes;
	for (int i = 0; i < nPatches; i++) delete[] numerator[i];
	delete[] numerator; 
}





