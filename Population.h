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

 RangeShifter v2.0 Population

 Implements the Population class

 There is ONE instance of a Population for each Species within each SubCommunity
 (including the matrix). The Population holds a list of all the Individuals in
 the Population.

 The matrix Population(s) hold(s) Individuals which are currently in the process
 of transfer through the matrix.

 For full details of RangeShifter, please see:
 Bocedi G., Palmer S.C.F., Pe’er G., Heikkinen R.K., Matsinos Y.G., Watts K.
 and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
 eco-evolutionary dynamics and species’ responses to environmental changes.
 Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

 Authors: Greta Bocedi & Steve Palmer, University of Aberdeen

 Last updated: 25 June 2021 by Steve Palmer

 ------------------------------------------------------------------------------*/

#ifndef PopulationH
#define PopulationH

#include <vector>
#include <algorithm>
using namespace std;

#include "Parameters.h"
#include "Individual.h"
#include "Species.h"
#include "Landscape.h"
#include "Patch.h"
#include "Cell.h"
#include "NeutralStatsManager.h"

//---------------------------------------------------------------------------

struct popStats {
	Species* pSpecies; Patch* pPatch; int spNum, nInds, nNonJuvs, nAdults; bool breeding;
};
struct disperser {
	Individual* pInd; Cell* pCell; bool yes;
};
struct zombie {
	Individual* pInd;
};
struct traitsums { // sums of trait genes for dispersal
	int ninds[maxNbSexes];				// no. of individuals
	double sumD0[maxNbSexes];			// sum of maximum emigration probability
	double ssqD0[maxNbSexes];			// sum of squares of maximum emigration probability
	double sumAlpha[maxNbSexes];	// sum of slope of emigration dens-dep reaction norm
	double ssqAlpha[maxNbSexes];	// sum of squares of slope of emigration den-dep reaction norm
	double sumBeta[maxNbSexes]; 	// sum of inflection point of emigration reaction norm
	double ssqBeta[maxNbSexes]; 	// sum of squares of inflection point of emigration reaction norm
	double sumDist1[maxNbSexes]; 	// sum of kernel I mean
	double ssqDist1[maxNbSexes]; 	// sum of squares of kernel I mean
	double sumDist2[maxNbSexes]; 	// sum of kernel II mean
	double ssqDist2[maxNbSexes]; 	// sum of squares of kernel II mean
	double sumProp1[maxNbSexes]; 	// sum of propn using kernel I
	double ssqProp1[maxNbSexes]; 	// sum of squares of propn using kernel I
	double sumDP[maxNbSexes]; 		// sum of SMS directional persistence
	double ssqDP[maxNbSexes]; 		// sum of squares of SMS directional persistence
	double sumGB[maxNbSexes]; 		// sum of SMS goal bias
	double ssqGB[maxNbSexes]; 		// sum of squares of SMS goal bias
	double sumAlphaDB[maxNbSexes];	// sum of SMS dispersal bias decay rate
	double ssqAlphaDB[maxNbSexes]; 	// sum of squares of SMS dispersal bias decay rate
	double sumBetaDB[maxNbSexes];		// sum of SMS dispersal bias decay infl. pt.
	double ssqBetaDB[maxNbSexes]; 	// sum of squares of SMS dispersal bias decay infl. pt.
	double sumStepL[maxNbSexes]; 	// sum of CRW step length
	double ssqStepL[maxNbSexes]; 	// sum of squares of CRW step length
	double sumRho[maxNbSexes]; 		// sum of CRW correlation coefficient
	double ssqRho[maxNbSexes]; 		// sum of squares of CRW correlation coefficient
	double sumS0[maxNbSexes];			// sum of maximum settlement probability
	double ssqS0[maxNbSexes];			// sum of squares of maximum settlement probability
	double sumAlphaS[maxNbSexes];	// sum of slope of settlement den-dep reaction norm
	double ssqAlphaS[maxNbSexes];	// sum of squares of slope of settlement den-dep reaction norm
	double sumBetaS[maxNbSexes]; 	// sum of inflection point of settlement reaction norm
	double ssqBetaS[maxNbSexes]; 	// sum of squares of inflection point of settlement reaction norm
	double sumFitness[maxNbSexes];
	double ssqFitness[maxNbSexes];
};

class Population {

public:
	Population(void); // default constructor
	Population( // constructor for a Population of a specified size
		Species*,	// pointer to Species
		Patch*,		// pointer to Patch
		int,			// no. of Individuals
		int				// Landscape resolution
	);
	~Population(void);
	traitsums getIndTraitsSums(Species*);
	popStats getStats(void);
	Species* getSpecies(void);
	int getNInds(void);
	int totalPop(void);
	int stagePop( // return no. of Individuals in a specified stage
		int	// stage
	);
	void extirpate(void); // Remove all individuals
	void reproduction(
		const float,	// local carrying capacity
		const float,	// effect of environmental gradient and/or stochasticty
		const int			// Landscape resolution
	);
	// Following reproduction of ALL species, add juveniles to the population
	void fledge(void);
	void emigration( // Determine which individuals will disperse
		float   // local carrying capacity
	);
	void allEmigrate(void); // All individuals emigrate after patch destruction
	// If an individual has been identified as an emigrant, remove it from the Population
	disperser extractDisperser(
		int		// index no. to the Individual in the inds vector
	);
	// For an individual identified as being in the matrix population:
	// if it is a settler, return its new location and remove it from the current population
	// otherwise, leave it in the matrix population for possible reporting before deletion
	disperser extractSettler(
		int   // index no. to the Individual in the inds vector
	);
	void recruit( // Add a specified individual to the population
		Individual*	// pointer to Individual
	);
	Individual* sampleInd() const;
	void sampleIndsWithoutReplacement(string n, const set<int>& sampleStages);
	int sampleSize() const;
	set<Individual*> getIndividualsInStage(int stage);
#if RS_RCPP
	int transfer( // Executed for the Population(s) in the matrix only
		Landscape*,	// pointer to Landscape
		short,				// landscape change index
		short				// year
	);
	// Determine whether there is a potential mate present in a patch which a potential
	// settler has reached
	bool matePresent(
		Cell*,	// pointer to the Cell which the potential settler has reached
		short		// sex of the required mate (0 = female, 1 = male)
	);
#else
	int transfer( // Executed for the Population(s) in the matrix only
		Landscape*,	// pointer to Landscape
		short				// landscape change index
	);
	// Determine whether there is a potential mate present in a patch which a potential
	// settler has reached
	bool matePresent(
		Cell*,	// pointer to the Cell which the potential settler has reached
		short		// sex of the required mate (0 = female, 1 = male)
	);
#endif // RS_RCPP
	// Determine survival and development and record in individual's status code
	// Changes are NOT applied to the Population at this stage
	void survival0(
		float,	// local carrying capacity
		short,	// option0:	0 - stage 0 (juveniles) only
		//	  			1 - all stages
		//					2 - stage 1 and above (all non-juveniles)
		short 	// option1:	0 - development only (when survival is annual)
						//	  	 		1 - development and survival
						//	  	 		2 - survival only (when survival is annual)
	);
	void survival1(void); // Apply survival changes to the population
	void ageIncrement(void);
	bool outPopHeaders( // Open population file and write header record
		int,	// Landscape number (-999 to close the file)
		bool	// TRUE for a patch-based model, FALSE for a cell-based model
	);
	void outPopulation( // Write record to population file
		int,		// replicate
		int,		// year
		int,		// generation
		float,	// epsilon - global stochasticity value
		bool,		// TRUE for a patch-based model, FALSE for a cell-based model
		bool,		// TRUE to write environmental data
		bool		// TRUE if there is a gradient in carrying capacity
	);

	void outIndsHeaders( // Open individuals file and write header record
		int,	// replicate
		int,	// Landscape number (-999 to close the file)
		bool	// TRUE for a patch-based model, FALSE for a cell-based model
	);
	void outIndividual( // Write records to individuals file
		Landscape*,	// pointer to Landscape
		int,				// replicate
		int,				// year
		int,				// generation
		int					// Patch number
	);

	void clean(void); // Remove zero pointers to dead or dispersed individuals

	void updatePopSNPtables();
	double getAlleleFrequency(int locus, int allele);
	int getAlleleTally(int locus, int allele);
	int getHeteroTally(int locus, int allele);
	int countHeterozygoteLoci();
	vector<double> countLociHeterozyotes();
	double computeHs();

private:
	short nStages;
	short nSexes;
	Species* pSpecies;	// pointer to the species
	Patch* pPatch;			// pointer to the patch
	int nInds[maxNbStages][maxNbSexes];		// no. of individuals in each stage/sex

	std::vector <Individual*> inds; // all individuals in population except ...
	std::vector <Individual*> juvs; // ... juveniles until reproduction of ALL species
	// has been completed

	std::set <Individual*> sampledInds;
	vector<SNPtable> popSNPtables;
	void resetPopSNPtables();
};

//---------------------------------------------------------------------------

extern paramGrad* paramsGrad;
extern paramStoch* paramsStoch;
extern paramInit* paramsInit;
extern paramSim* paramsSim;
extern RSrandom* pRandom;

#if RSDEBUG
extern ofstream DEBUGLOG;
#endif

//---------------------------------------------------------------------------
#endif

