#ifndef NEUTRALSTATSH
#define NEUTRALSTATSH

#include "Species.h"
#include "Landscape.h"

using namespace std;

// Patch * patch matrix to store pairwise Fst calculations
/** Creates an array of doubles of size = rows*cols, taken from NEMO **/
struct PatchMatrix
{
public:
	PatchMatrix(int rows = 0, int cols = 0) : rows(0), cols(0), nbCells(0), value(0) {
		nbCells = rows * cols;
		value.resize(nbCells);
		rows = rows; cols = cols;
	};

	// Get value at specified position
	double get(unsigned int i, unsigned int j) {
		if (!((i + 1) * (j + 1) > nbCells))
			return value[i * cols + j];
		else throw runtime_error("Error: PatchMatrix::get() out of range!\n");
		return 0;
	}

	int getNbCells() { return nbCells; };

	/** Sets element at row i and column j to value val **/
	void set(unsigned int i, unsigned int j, double val) {
		if (i * j < nbCells)
			value[i * cols + j] = val;
		else throw runtime_error("Error: PatchMatrix::set() out of range!\n");
	}

	/** Assigns a value to all elements of the matrix. */
	void setAll(double val)
	{
		for (unsigned int i = 0; i < nbCells; ++i) value[i] = val;
	}

private:
	unsigned int rows, cols, nbCells;
	vector<double> value;
};

// Counts of SNP allele occurrences in populations
// for neutral statistics calculations
struct SNPCountsTable {

public:
	SNPCountsTable(int nAllele) : alleleTallies(nAllele), alleleFrequencies(nAllele), alleleHeterozygoteTallies(nAllele) {};
	
	void reset() {
		fill(alleleTallies.begin(), alleleTallies.end(), 0); fill(alleleFrequencies.begin(), alleleFrequencies.end(), 0);
		fill(alleleHeterozygoteTallies.begin(), alleleHeterozygoteTallies.end(), 0);
	}

	// Getters
	int getTally(int whichAllele) { return alleleTallies[whichAllele]; };
	double getFrequency(int whichAllele) { return alleleFrequencies[whichAllele]; };
	int getHeteroTally(int whichAllele) { return alleleHeterozygoteTallies[whichAllele]; };

	// Setters / increments
	void incrementTally(int whichAllele) { alleleTallies[whichAllele]++; };
	void incrementTallyBy(int count, int whichAllele) { this->alleleTallies[whichAllele] += count; }
	void incrementHeteroTally(int whichAllele) { this->alleleHeterozygoteTallies[whichAllele]++; }
	void setFrequencies(int sampleSize) {
		for (int i = 0; i < alleleFrequencies.size(); i++) {
			alleleFrequencies[i] = sampleSize > 0 ? static_cast<double>(alleleTallies[i]) / sampleSize : 0.0;
		}
	};

private:
	// Tallies, one for each possible allele (so absolute max size is 255)
	vector<int> alleleTallies; // number of occurrences of each allele in pop
	vector<double> alleleFrequencies; // frequency of each allele in pop
	vector<int> alleleHeterozygoteTallies; // nb of times each allele is found in a heterozygous pair
};

class NeutralStatsManager {

private:
	int  nbExtantPops, totalNbSampledInds;
	/** F-statistics */
	double _ho, _hs, _ht, _hsnei, _htnei, meanNbAllelesPerLocusPerPatch, meanNbAllelesPerLocus,
		_fst, _fis, _fit, meanNbFixedAllelesPerPatch, nbGloballyFixedAlleles;
	/** Weir & Hill (2002) F-stat estimates. */
	double weightedFstWeirHill;
	/** Weir & Cockerham (1984) F-stat estimates. */
	double FstWeirCockerham, FisWeirCockerham, FitWeirCockerham;
	/** Per-locus F-stats (Weir&Cockerham). */
	vector<double> perLocusFstWeirCockerham, perLocusFisWeirCockerham, perLocusFitWeirCockerham, perLocusHo; //no need for pointers because shouldn't be copied or moved, resized 

	/** Pairwise Fst matrix. */
	PatchMatrix pairwiseFstMatrix;
	vector<SNPCountsTable> commSNPTables; //don't have to be pointers, not shared or moved

public: 

	NeutralStatsManager(const int& nbSampledPatches, const int nLoci);

	void updateAllSNPTables(Species* pSpecies, Landscape* pLandscape, set<int> const& patchList);
	void resetCommSNPtables();
	void setLociDiversityCounter(set<int> const& patchList, const int nInds, Species* pSpecies, Landscape* pLandscape);
	void setFstMatrix(set<int> const& patchList, const int nInds, const int nLoci, Species* pSpecies, Landscape* pLandscape);
	void calculateHo(set<int> const& patchList, const int totalNbSampledInds, const int nbrLoci, Species* pSpecies, Landscape* pLandscape);
	void calculateHs(set<int> const& patchList, const int nbrLoci, Species* pSpecies, Landscape* pLandscape);
	void calculateHt(Species* pSpecies, Landscape* pLandscape, const int nLoci, const int nAlleles);
	void calculateHo2(set<int> const& patchList, const int totalNbSampledInds, const int nbrLoci, Species* pSpecies, Landscape* pLandscape);
	void calculateFstatWC(set<int> const& patchList, const int nInds, const int nLoci, const int nAlleles, Species* pSpecies, Landscape* pLandscape);
	void calculateFstatWC_MS(set<int> const& patchList, const int nInds, const int nLoci, const int nAlleles, Species* pSpecies, Landscape* pLandscape);

	double getHsnei() const { return _hsnei; }
	double getHtnei() const { return _htnei; }
	double getHo() const { return _ho; }
	double getHs() const { return _hs; }
	double getHt() const { return _ht; }
	double getFst() const { return _fst; }
	double getFis() const { return _fis; }
	double getFit() const { return _fit; }
	double getFstWC() const { return FstWeirCockerham; }
	double getFisWC() const { return FisWeirCockerham; }
	double getFitWC() const { return FitWeirCockerham; }
	int getNbPopulatedSampledPatches() const { return nbExtantPops;  }
	int getTotalNbSampledInds() const { return totalNbSampledInds;  }
	double getWeightedFst() { return weightedFstWeirHill; }
	double getMeanNbAllPerLocusPerPatch() const { return meanNbAllelesPerLocusPerPatch; }
	double getMeanNbAllPerLocus() const { return meanNbAllelesPerLocus; }
	double getMeanFixdAllelesPerPatch() const { return meanNbFixedAllelesPerPatch; }
	double getTotalFixdAlleles() const { return nbGloballyFixedAlleles; }
	double getPairwiseFst(int i, int j) { return pairwiseFstMatrix.get(i, j); }
	double get_fst_WC_loc(int i) const { return perLocusFstWeirCockerham[i]; }
	double get_fis_WC_loc(int i) const { return perLocusFisWeirCockerham[i]; }
	double get_fit_WC_loc(int i) const { return perLocusFitWeirCockerham[i]; }
	double get_ho_loc(int i) const { return perLocusHo[i]; }
};

#endif




