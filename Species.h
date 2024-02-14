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
 --------------------------------------------------------------------------*/


 /*------------------------------------------------------------------------------

 RangeShifter v2.0 Species

 Implements the Species class

 There is ONE instance of a Species for each species within the Community
 AND THIS IS CURRENTLY LIMITED TO A SINGLE SPECIES.
 The class holds all the demographic and dispersal parameters of the species.

 For full details of RangeShifter, please see:
 Bocedi G., Palmer S.C.F., Pe’er G., Heikkinen R.K., Matsinos Y.G., Watts K.
 and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
 eco-evolutionary dynamics and species’ responses to environmental changes.
 Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

 Authors: Greta Bocedi & Steve Palmer, University of Aberdeen

 Last updated: 28 July 2021 by Greta Bocedi

 ------------------------------------------------------------------------------*/

#ifndef SpeciesH
#define SpeciesH

#include "Parameters.h"
#include "SpeciesTrait.h"
#include "TTrait.h"
#include <ranges>
#include <map>
#include <set>

class SpeciesTrait;

 // structures for demographic parameters

struct demogrParams {
	short repType;
	short repSeasons;
	float propMales; float harem; float bc; float lambda;
	bool stageStruct;
};
struct stageParams {
	short nStages; short repInterval; short maxAge; short survival;
	float probRep;
	bool fecDens;  bool fecStageDens; bool devDens; bool devStageDens;
	bool survDens; bool survStageDens; bool disperseOnLoss;
};
struct densDepParams {
	float devCoeff; float survCoeff;
};

// structures for emigration parameters

struct emigRules {
	bool densDep; bool stgDep; bool sexDep; bool indVar;
	short emigStage;
	short emigTrait[2];
};
struct emigTraits {
	float d0; float alpha; float beta;

	emigTraits() : d0(0.0), alpha(0.0), beta(0.0) {};

	emigTraits(const emigTraits& e) : d0(e.d0), alpha(e.alpha), beta(e.beta) {};

	emigTraits* clone() { return new emigTraits(*this); }

	void divideTraitsBy(int i) {

		d0 /= i;
		alpha /= i;
		beta /= i;
	}
};

// structures for transfer parameters

struct trfrRules {
	bool moveModel; bool stgDep; bool sexDep;
	bool distMort; bool indVar;
	bool twinKern;
	bool habMort;
	short moveType; bool costMap;
	short movtTrait[2];
};
struct trfrKernTraits {
	float	meanDist1;
	float	meanDist2;
	float	probKern1;
};
struct trfrMortParams {
	float fixedMort; float mortAlpha; float mortBeta;
};
struct trfrMovtTraits {
	short	pr; short	prMethod; short	memSize; short goalType;
	float	dp; float	gb; float alphaDB; int betaDB;
	float	stepMort; float	stepLength; float	rho;
	bool straigtenPath;
};
struct trfrCRWTraits {
	float	stepMort;
	float	stepLength;
	float	rho;
	bool	straigtenPath;
};
struct trfrSMSTraits {
	short	pr; short	prMethod; short	memSize; short goalType;
	float	dp; float	gb; float alphaDB; int betaDB; float	stepMort;
	bool straigtenPath;
};

// structures for settlement parameters

struct settleType {
	bool stgDep; bool sexDep; bool indVar;
	short settTrait[2];
};
struct settleRules {
	bool densDep; bool wait; bool go2nbrLocn; bool findMate;
};
struct settleSteps {
	int minSteps; int maxSteps; int maxStepsYr;
};
struct settleTraits {
	float s0; float alpha; float beta;

	settleTraits() : s0(0.0), alpha(0.0), beta(0.0) {};

	settleTraits(const settleTraits& e) : s0(e.s0), alpha(e.alpha), beta(e.beta) {};

	void divideTraitsBy(int i) {
		s0 /= i;
		alpha /= i;
		beta /= i;
	}
};


//---------------------------------------------------------------------------

class Species {

public:
	Species(void);
	~Species(void);
	short getSpNum(void);

	// demographic parameter functions

	void createHabK( // Create habitat carrying capacity table
		short	// no. of habitats
	);
	void setHabK(
		short,	// habitat index no. (NB may differ from habitat no. supplied by user)
		float		// carrying capacity (inds/cell)
	);
	float getHabK(
		short		// habitat index no. (NB may differ from habitat no. supplied by user)
	);
	float getMaxK(void); // return highest carrying capacity over all habitats
	void deleteHabK(void); // Delete habitat carrying capacity table
	void setStage( // Set stage structure parameters
		const stageParams	// structure holding stage structure parameters
	);
	stageParams getStage(void); // Get stage structure parameters
	void setDemogr( // Set general demographic parameters
		const demogrParams	// structure holding general demographic parameters
	);
	demogrParams getDemogr(void); // Get general demographic parameters
	short getRepType(void);
	bool stageStructured(void);
	void setDensDep( // Set demographic density dependence coefficients
		float,	// development coefficient
		float		// survival coefficient
	);
	densDepParams getDensDep(void); // Get development and survival coefficients

	void setFec( // Set fecundity
		short,	// stage (must be > 0)
		short,	// sex
		float		// fecundity
	);
	float getFec( // Get fecundity
		short,	// stage
		short		// sex
	);
	void setDev( // Set development probability
		short,	// stage
		short,	// sex
		float		// development probability
	);
	float getDev( // Get development probability
		short,	// stage
		short		// sex
	);
	void setSurv( // Set survival probability
		short,	// stage
		short,	// sex
		float		// survival probability
	);
	float getSurv( // Get survival probability
		short,	// stage
		short		// sex
	);

	float getMaxFec(void); // Get highest fecundity of any stage
	void setMinAge( // Set minimum age
		short,	// stage
		short,	// sex
		int			// minimum age (years) (must be zero for stages 0 and 1)
	);
	short getMinAge( // Get minimum age
		short,	// stage
		short		// sex
	);
	void createDDwtFec( // Create fecundity weights matrix
		short		// matrix dimension - no. of stages * no. of sexes
	);
	void setDDwtFec( // Set fecundity weights matrix element
		short,	// row
		short,	// column
		float		// weight
	);
	float getDDwtFec( // Get fecundity weights matrix element
		short,	// row
		short		// column
	);
	void deleteDDwtFec(void); // Delete fecundity weights matrix
	void createDDwtDev( // Create development weights matrix
		short		// matrix dimension - no. of stages * no. of sexes
	);
	void setDDwtDev( // Set development weights matrix element
		short,	// row
		short,	// column
		float		// weight
	);
	float getDDwtDev( // Get development weights matrix element
		short,	// row
		short		// column
	);
	void deleteDDwtDev(void); // Delete development weights matrix
	void createDDwtSurv( // Create survival weights matrix
		short		// matrix dimension - no. of stages * no. of sexes
	);
	void setDDwtSurv( // Set survival weights matrix element
		short,	// row
		short,	// column
		float		// weight
	);
	float getDDwtSurv( // Get survival weights matrix element
		short,	// row
		short		// column
	);
	void deleteDDwtSurv(void); // Delete survival weights matrix
	// Functions to handle min/max R or K (under environmental stochasticity)
	void setMinMax( // Set min and max values
		float,	// min
		float		// max
	);
	float getMinMax( // Get min/max value
		short		// option: 0 = return minimum, otherwise = return maximum
	);

	std::set<int>& getSamplePatches() {
		return samplePatchList;
	};

	string getNIndsToSample() {
		return nIndsToSample;
	};

	std::set<int>& getStagesToSample() {
		return stagesToSampleFrom;
	}

	string getNSampleCellsFst() {
		return nSampleCellsFst;
	}

	// Genetic functions
	void resetGeneticParameters();
	bool areMutationsOn(void);
	bool isDiploid() const;
	void setNumberOfNeutralLoci(int);
	int getNumberOfNeutralLoci() const;
	int incrementAdaptiveTraits();
	int getNumberOfAdaptiveTraits() const;

	// emigration parameter functions

	void setEmig( // Set emigration rules
		const emigRules	// structure holding emigration rules
	);
	emigRules getEmig(void); // Get emigration rules
	void setEmigTraits( // Set emigration trait parameters
		const short,			// stage
		const short,			// sex
		const emigTraits	// structure holding emigration trait parameters
	);
	emigTraits getEmigTraits( // Get emigration trait parameters
		short,	// stage
		short		// sex
	);
	float getEmigD0( // Get (maximum) emigration probability
		short,	// stage
		short		// sex
	);

	// transfer parameter functions

	void setTrfr( // Set transfer rules
		const trfrRules	// structure holding transfer rules
	);
	trfrRules getTrfr(void); // Get transfer rules
	void setFullKernel( // Set fullKernel condition
		bool						// fullKernel value
	);
	bool useFullKernel(void);
	void setKernTraits( // Set transfer by kernel parameters
		const short,					// stage
		const short,					// sex
		const trfrKernTraits,	// structure holding transfer by kernel parameters
		const int							// Landscape resolution
	);
	trfrKernTraits getKernTraits( // Get transfer by kernel parameters
		short,	// stage
		short		// sex
	);
	void setMortParams( // Set transfer mortality parameters
		const trfrMortParams	// structure holding transfer mortality parameters
	);
	trfrMortParams getMortParams(void); // Get transfer mortality parameters
	void setMovtTraits( // Set transfer movement model parameters
		const trfrMovtTraits	// structure holding transfer movement model parameters
	);
	trfrMovtTraits getMovtTraits(void); // Get transfer movement model traits
	trfrCRWTraits getCRWTraits(void);		// Get CRW traits
	trfrSMSTraits getSMSTraits(void);		// Get SMS traits
	// Return dimension of habitat-dependent step mortality and costs matrices
	short getMovtHabDim(void);
	void createHabCostMort( // Create habitat-dependent costs and mortality matrices
		short	// no. of habitats
	);
	void setHabCost( // Set habitat-dependent cost
		short,	// habitat index no.
		int			// cost value
	);
	void setHabMort( // Set habitat-dependent per-step mortality
		short,	// habitat index no.
		double	// mortality rate
	);
	int getHabCost( // Get habitat-dependent cost
		short		// habitat index no.
	);
	double getHabMort( // Get habitat-dependent per-step mortality
		short		// habitat index no.
	);
	void deleteHabCostMort(void); // Delete habitat-dependent costs and mortality matrices

	// settlement parameter functions

	void setSettle( // Set settlement type
		const settleType	// structure holding settlement type (stage- and/or sex-dependent)
	);
	settleType getSettle(void); // Get settlement type
	void setSettRules( // Set settlement rules
		const short,			// stage
		const short,			// sex
		const settleRules	// structure holding settlement rules
	);
	settleRules getSettRules( // Get settlement rules
		short,	// stage
		short		// sex
	);
	void setSteps( // Set path step limit parameters
		const short,			// stage
		const short,			// sex
		const settleSteps	// structure holding path step limit parameters
	);
	settleSteps getSteps( // Set path step limit parameters
		short,	// stage
		short		// sex
	);
	void setSettTraits( // Set settlement density dependence traits
		const short,					// stage
		const short,					// sex
		const settleTraits	// structure holding density dependence traits
	);
	settleTraits getSettTraits( // Get settlement density dependence traits
		short,	// stage
		short		// sex
	);

	void addTrait(TraitType traitType, const SpeciesTrait& trait);

	void clearTraitTable();

	SpeciesTrait* getSpTrait(TraitType trait) const;

	//map<TraitType, std::unique_ptr<ProtoTrait>>& getTraitTable(void); //return by reference so ensure variable recieving is const

	std::set<TraitType> getTraitTypes();

	int getNTraits() const;
	int getNPositionsForTrait(const TraitType trait) const;
	int getGenomeSize() const;
	float getRecombinationRate() const;
	std::set<int> getChromosomeEnds() const;
	void setGeneticParameters(const std::set<int>& chromosomeEnds, const int genomeSize, const float recombinationRate,
		const std::set<int>& samplePatchList, const string nIndsToSample, const std::set<int>& stagesToSampleFrom, string nSampleCellsFst);
	void setSamplePatchList(const std::set<int>& samplePatchList);

private:

	// NOTE: SEQUENCE OF PARAMETER VARIABLES MAY NEED TO BE ALTERED FOR EFFICIENTCY ...
	// ... but that is of low importance, as there will only be one (or a few) instance(s)

	// demographic parameters

	short repType;			// 0 = asexual, 1 = simple two sex, 2 = complex two sex
	short nStages;      // no. of stages (incl. juveniles) in structured population
	float propMales;    // proportion of males at birth in sexual model
	float harem;        // max harem size in complex sexual model
	float bc;						// competition coefficient for non-structured population
	float lambda;       // max intrinsic growth rate for non-structured population
	float probRep; 			// probability of reproducing in subsequent seasons
	short repSeasons;		// no. of reproductive seasons per year
	short repInterval;	// no. of reproductive seasons between subsequent reproductions
	short maxAge;       // max age in structured population
	short survival;			// survival timing: 0 = at reprodn, 1 = between reprodns, 2 = anually
	bool stageStruct;
	bool fecDens;
	bool fecStageDens;
	bool devDens;
	bool devStageDens;
	bool survDens;
	bool survStageDens;
	bool disperseOnLoss;	// individuals disperse on complete loss of patch
	// (otherwise they die)
	short habDimK;			// dimension of carrying capacities matrix
	float* habK;				// habitat-specific carrying capacities (inds/cell)
	float devCoeff; 		// density-dependent development coefficient
	float survCoeff; 		// density-dependent survival coefficient
	float** ddwtFec;    // density-dependent weights matrix for fecundity
	float** ddwtDev;    // density-dependent weights matrix for development
	float** ddwtSurv;   // density-dependent weights matrix for survival
	// NB for the following arrays, sex 0 is females, sex 1 is males
	float fec[NSTAGES][NSEXES];			// fecundities
	float dev[NSTAGES][NSEXES];			// development probabilities
	float surv[NSTAGES][NSEXES];		// survival probabilities
	short minAge[NSTAGES][NSEXES];	// minimum age to enter stage
	// NOTE - IN THEORY, NEXT 3 VARIABLES COULD BE COMMON, BUT WE WOULD NEED TO ENSURE THAT
	// ALL MATRICES ARE DELETED IF THERE IS A CHANGE IN NO. OF STAGES OR REPRODUCTION TYPE
	// ***** TO BE RECONSIDERED LATER *****
	short ddwtFecDim;		// dimension of density-dependent weights matrix for fecundity
	short ddwtDevDim;		// dimension of density-dependent weights matrix for fecundity
	short ddwtSurvDim;	// dimension of density-dependent weights matrix for fecundity
	float minRK; 				// minimum ) growth rate OR carrying capacity
	float maxRK; 				// maximum ) (under environmental stochasticity)

	// genome parameters

	/**The traits table.*/
	std::map<TraitType, std::unique_ptr<SpeciesTrait>> spTraitTable;
	std::set<int> chromosomeEnds;
	int genomeSize;
	bool diploid;
	bool mutationsOn;
	int numberOfNeutralLoci;
	int numberOfAdaptiveTraits;
	float recombinationRate;
	std::set<int> samplePatchList;
	string nSampleCellsFst; //for cell based landscape
	std::set<int> stagesToSampleFrom;
	string nIndsToSample; //could be integer or 'all', all means in in selected patches not necessarily all in population

	// emigration parameters

	bool	densDepEmig;	// density-dependent emigration
	bool	stgDepEmig;   // stage-dependent emigration
	bool	sexDepEmig;   // sex-dependent emigration
	bool	indVarEmig;   // individual variation in emigration
	short emigStage;		// stage which emigrates (used for stage-strucutred population
	// having individual variability in emigration probability)
// NB for the following arrays, sex 0 is females, sex 1 is males
	float	d0[NSTAGES][NSEXES];				 // maximum emigration probability
	float	alphaEmig[NSTAGES][NSEXES];	 // slope of density-dependent reaction norm
	float	betaEmig[NSTAGES][NSEXES];	 // inflection point of reaction norm (in terms of N/K)
	// NB Initialisation parameters are made double to avoid conversion errors (reason unclear)
	// on traits maps using FloatToStr()

	// transfer parameters

	bool moveModel;
	bool stgDepTrfr;
	bool sexDepTrfr;
	bool distMort;
	bool indVarTrfr;
	bool twinKern;
	bool habMort;		// habitat-dependent mortality
	float	meanDist1[NSTAGES][NSEXES];	// mean of 1st dispersal kernel (m)
	float	meanDist2[NSTAGES][NSEXES]; // mean of 2nd dispersal kernel (m)
	float	probKern1[NSTAGES][NSEXES]; // probability of dispersing with the 1st kernel
	// NB INITIAL limits are made double to avoid conversion errors (reason unclear)
	// on traits maps using FloatToStr()
	// As evolving traits are are not stage-dependent, no. of rows can be 1
	float fixedMort;		// constant mortality probability
	float mortAlpha;		// slope for mortality distance dependence function
	float mortBeta;			// inflection point for mortality distance dependence function
	short moveType; 		// 1 = SMS, 2 = CRW
	short pr;						// SMS perceptual range (cells)
	short prMethod;			// SMS perceptual range evaluation method:
	// 1 = arith. mean, 2 = harmonic mean, 3 = inverse weighted arith. mean
	short memSize;			// SMS memory size (1-14 steps)
	short goalType;			// SMS goal bias type: 0 = none, 1 = towards goal, 2 = dispersal bias
	float dp;						// SMS directional persistence
	float gb;						// SMS goal bias
	float alphaDB; 			// SMS dispersal bias decay rate
	int betaDB; 				// SMS dispersal bias decay inflection point (no. of steps)
	float stepMort;			// constant per-step mortality probability for movement models
	double* habStepMort;	// habitat-dependent per-step mortality probability
	float stepLength;		// CRW step length (m)
	float rho;					// CRW correlation coefficient
	short habDimTrfr;		// dimension of habitat-dependent step mortality and costs matrices
	int* habCost;				// habitat costs
	bool costMap;				// import cost map from file?
	bool straigtenPath;	// straighten path after decision not to settle
	bool fullKernel;		// used to indicate special case when density-independent emigration
	// is 1.0, and kernel-based movement within the natal cell is used
	// to determine philopatry

// settlement parameters

	bool stgDepSett;
	bool sexDepSett;
	bool indVarSett;   								// individual variation in settlement
	bool densDepSett[NSTAGES][NSEXES];
	bool wait[NSTAGES][NSEXES];				// wait to continue moving next season (stage-structured model only)
	bool go2nbrLocn[NSTAGES][NSEXES];	// settle in neighbouring cell/patch if available (ditto)
	bool findMate[NSTAGES][NSEXES];
	int minSteps;     								// minimum no. of steps
	int maxSteps;											// maximum total no. of steps
	int maxStepsYr[NSTAGES][NSEXES]; 	// maximum no. of steps in any one dispersal period
	float	s0[NSTAGES][NSEXES];				// maximum settlement probability
	float alphaS[NSTAGES][NSEXES];		// slope of the settlement reaction norm to density
	float betaS[NSTAGES][NSEXES];			// inflection point of the settlement reaction norm to density

	// other attributes
	int spNum;

};


//---------------------------------------------------------------------------

#if RSDEBUG
//extern ofstream DEBUGLOG;
extern void DebugGUI(string);
#endif

//---------------------------------------------------------------------------
#endif
