#ifndef NEUTRALSTATSH
#define NEUTRALSTATSH

#include "Species.h"
#include "Landscape.h"

using namespace std;

/**Creates an array of doubles of size = rows*cols, taken from NEMO**/
struct PatchMatrix
{
private:
	unsigned int _rows, _cols, _length;
	vector<double> _val;

public:

	PatchMatrix(int rows = 0, int cols = 0) : _rows(0), _cols(0), _length(0), _val(0) {
		_length = rows * cols;
		_val.resize(_length);
		_rows = rows; _cols = cols;
	};


	/**Assigns a value to all element of the matrix.*/
	void assign(double val)
	{
		for (unsigned int i = 0; i < _length; ++i) _val[i] = val;
	}

	int length() { return _length; };

	/**Sets element at row i and column j to value val**/
	void set(unsigned int i, unsigned int j, double val) {
		if (i * j < _length)
			_val[i * _cols + j] = val;
		else
			cout << endl << ("Error: PatchMatrix::set() out of range!\n");
	}

	double get(unsigned int i, unsigned int j) {
		if (!((i + 1) * (j + 1) > _length))
			return _val[i * _cols + j];
		else
			cout << endl << ("Error: PatchMatrix::get() out of range!\n");
		return 0;
	}
};


struct SNPtable {

private:
	vector<int> alleleTallies;
	vector<double> alleleFrequencies;
	vector<int> alleleHeterozygoteTallies;

public:
	//for community allele table, heteros not technically needed so don't reserve 
	SNPtable(int nAllele, int allele, int alleleCount) : alleleTallies(nAllele), alleleFrequencies(nAllele) {
		this->incrementTallyBy(alleleCount, allele);
	};

	//for population allele tables 
	SNPtable(int nAllele) : alleleTallies(nAllele), alleleFrequencies(nAllele), alleleHeterozygoteTallies(nAllele) {};

	void setFrequencies(int populationSize) {
		for (int i = 0; i < alleleFrequencies.size(); i++) {
			alleleFrequencies[i] = alleleTallies[i] / static_cast<double>(populationSize);
		}
	};

	void reset() {
		fill(alleleTallies.begin(), alleleTallies.end(), 0); fill(alleleFrequencies.begin(), alleleFrequencies.end(), 0);
		fill(alleleHeterozygoteTallies.begin(), alleleHeterozygoteTallies.end(), 0);
	}

	int getTally(int whichAllele) { return alleleTallies[whichAllele]; };
	double getFrequency(int whichAllele) { return alleleFrequencies[whichAllele]; };
	int getHeteroTally(int whichAllele) { return alleleHeterozygoteTallies[whichAllele]; };

	void incrementTally(int whichAllele) { alleleTallies[whichAllele]++; };
	void incrementTallyBy(int count, int whichAllele) { this->alleleTallies[whichAllele] += count; }
	void incrementHeteroTally(int whichAllele) { this->alleleHeterozygoteTallies[whichAllele]++; }
};



class NeutralStatsManager {

private:
	int  _n_extantPopulations, _n_individuals;
	/**F-statistics*/
	double _ho, _hs, _ht, _hsnei, _htnei, meanNbAllelesPerLocusPerPatch, meanNbAllelesPerLocus,
		_fst, _fis, _fit, meanNbFixedAllelesPerPatch, nbGloballyFixedAlleles;
	/**Weir & Hill (2002) F-stat estimates.*/
	double _fst_WH;
	/**Weir & Cockerham (1984) F-stat estimates.*/
	double _fst_WC, _fis_WC, _fit_WC;
	/**Per-locus F-stats (Weir&Cockerham).*/
	vector<double> _fst_WC_loc, _fis_WC_loc, _fit_WC_loc, ho_loc; //no need for pointers because shouldn't be copied or moved, resized 

	/**Pairwise Fst matrix.*/
	PatchMatrix _fst_matrix;
	vector<SNPtable> globalSNPtables; //don't have to be pointers, not shared or moved

public: 

	NeutralStatsManager(const int& nbSampledPatches, const int nLoci);

	void updateAllSNPTables(Species* pSpecies, Landscape* pLandscape, set<int> const& patchList);

	void resetGlobalSNPtables();

	void setLociDiversityCounter(set<int> const& patchList, const int nInds, Species* pSpecies, Landscape* pLandscape);

	void calculateHo(set<int> const& patchList, const int nbInds, const int nbrLoci, Species* pSpecies, Landscape* pLandscape);
	void calculateHs(set<int> const& patchList, const int nbrLoci, Species* pSpecies, Landscape* pLandscape);
	void calculateHt(Species* pSpecies, Landscape* pLandscape, const int nLoci, const int nAlleles);
	void calculateHo2(set<int> const& patchList, const int nbInds, const int nbrLoci, Species* pSpecies, Landscape* pLandscape);

	void calculateFstatWC(set<int> const& patchList, const int nInds, const int nLoci, const int nAlleles, Species* pSpecies, Landscape* pLandscape);
	void calculateFstatWC_MS(set<int> const& patchList, const int nInds, const int nLoci, const int nAlleles, Species* pSpecies, Landscape* pLandscape);

	void setFstMatrix(set<int> const& patchList, const int nInds, const int nLoci, Species* pSpecies, Landscape* pLandscape);

	double getHsnei() const { return _hsnei; }
	double getHtnei() const { return _htnei; }
	double getHo() const { return _ho; }
	double getHs() const { return _hs; }
	double getHt() const { return _ht; }
	double getFst() const { return _fst; }
	double getFis() const { return _fis; }
	double getFit() const { return _fit; }
	double getFstWC() const { return _fst_WC; }
	double getFisWC() const { return _fis_WC; }
	double getFitWC() const { return _fit_WC; }
	int getNExtantPatchs() const { return _n_extantPopulations;  }
	int getNIndividuals() const { return _n_individuals;  }
	double getWeightedFst() { return _fst_WH; }
	double getMeanNbAllPerLocusPerPatch() const { return meanNbAllelesPerLocusPerPatch; }
	double getMeanNbAllPerLocus() const { return meanNbAllelesPerLocus; }
	double getMeanFixdAllelesPerPatch() const { return meanNbFixedAllelesPerPatch; }
	double getTotalFixdAlleles() const { return nbGloballyFixedAlleles; }
	double getPairwiseFst(int i, int j) { return _fst_matrix.get(i, j); }
	double get_fst_WC_loc(int i) const { return _fst_WC_loc[i]; }
	double get_fis_WC_loc(int i) const { return _fis_WC_loc[i]; }
	double get_fit_WC_loc(int i) const { return _fit_WC_loc[i]; }
	double get_ho_loc(int i) const { return ho_loc[i]; }
};

#endif




