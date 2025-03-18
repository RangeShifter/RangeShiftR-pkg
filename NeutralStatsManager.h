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
	PatchMatrix(unsigned int nRows = 0, unsigned int nCols = 0) : 
		rows{ nRows }, 
		cols{ nCols }, 
		nbCells{nCols * nRows} {
		value.resize(nbCells);
	};

	// Get value at specified position
	double get(unsigned int i, unsigned int j) {
		if (i >= cols || j >= rows)
			throw runtime_error("PatchMatrix::get() out of range!\n");
		else return value[i * cols + j];
	}

	int getNbCells() { return nbCells; };

	/** Sets element at row i and column j to value val **/
	void set(unsigned int i, unsigned int j, double val) {
		if (i >= cols || j >= rows)
			throw runtime_error("PatchMatrix::set() out of range!\n");
		else value[i * cols + j] = val;
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

// Counts of NEUTRAL allele occurrences in populations
// for neutral statistics calculations
struct NeutralCountsTable {

public:
	NeutralCountsTable(int maxNbNeutralAlleles) : alleleTallies(maxNbNeutralAlleles), alleleFrequencies(maxNbNeutralAlleles), alleleHeterozygoteTallies(maxNbNeutralAlleles) {};
	
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


// Handles calculations of neutral statistics
class NeutralStatsManager {

public:

	NeutralStatsManager(const int& nbSampledPatches, const int nLoci);

	// Count alleles and their frequencies in all pops and community
	void updateAllNeutralTables(Species* pSpecies, Landscape* pLandscape, set<int> const& patchList);
	void resetCommNeutralTables();

	void calcAllelicDiversityMetrics(set<int> const& patchList, const int nInds, Species* pSpecies, Landscape* pLandscape);
	
	// Heterozygosity calculations
	void calculateHo(set<int> const& patchList, const int totalNbSampledInds, const int nbrLoci, Species* pSpecies, Landscape* pLandscape);
	void calculateHs(set<int> const& patchList, const int nbrLoci, Species* pSpecies, Landscape* pLandscape);
	void calculateHt(Species* pSpecies, Landscape* pLandscape, const int nLoci, const int nAlleles);
	void calculatePerLocusHo(set<int> const& patchList, const int totalNbSampledInds, const int nbrLoci, Species* pSpecies, Landscape* pLandscape);
	
	// F-stats calculations
	void calculateFstatWC(set<int> const& patchList, const int nInds, const int nLoci, const int nAlleles, Species* pSpecies, Landscape* pLandscape);
	void calcPairwiseWeightedFst(set<int> const& patchList, const int nInds, const int nLoci, Species* pSpecies, Landscape* pLandscape);

	// Getters
	int getNbPopulatedSampledPatches() const { return nbExtantPops; }
	int getTotalNbSampledInds() const { return totalNbSampledInds; }
	
	double getMeanNbAllPerLocusPerPatch() const { return meanNbAllelesPerLocusPerPatch; }
	double getMeanNbAllPerLocus() const { return meanNbAllelesPerLocus; }
	double getMeanFixdAllelesPerPatch() const { return meanNbFixedLociPerPatch; }
	double getTotalFixdAlleles() const { return meanFixedLoci; }
	
	double getHo() const { return ho; }
	double getHs() const { return hs; }
	double getHt() const { return ht; }

	double getPerLocusHo(int i) const { return perLocusHo[i]; }

	double getFstWC() const { return fst; }
	double getFisWC() const { return fis; }
	double getFitWC() const { return fit; }

	double getWeightedFst() { return weightedFst; }

	double getPerLocusFst(int i) const { return perLocusFst[i]; }
	double getPerLocusFis(int i) const { return perLocusFis[i]; }
	double getPerLocusFit(int i) const { return perLocusFit[i]; }

	double getPairwiseFst(int i, int j) { return pairwiseFstMatrix.get(i, j); }

private:

	int nbExtantPops, totalNbSampledInds;
	vector<NeutralCountsTable> commNeutralCountTables; // community-level tallies of allele counts and freqs

	double meanNbAllelesPerLocusPerPatch, meanNbAllelesPerLocus;
	double meanNbFixedLociPerPatch, meanFixedLoci;

	double ho; // observed heterozygosity 
	double hs; // expected population-level heterozygosity
	double ht; // expected community-level heterozygosity
	
	vector<double> perLocusHo; // Per-locus observed heterozygosity

	// F-statistics
	// Weir & Cockerham (1984) F-stat estimates.
	double fst, fis, fit;

	// Weir & Hill (2002) F-stat estimates 
	double weightedFst;

	// Per-locus F-stats (Weir & Cockerham).
	vector<double> perLocusFst, perLocusFis, perLocusFit;

	// Pairwise Fst matrix
	PatchMatrix pairwiseFstMatrix;
};

#endif




