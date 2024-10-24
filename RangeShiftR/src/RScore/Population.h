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
Bocedi G., Palmer S.C.F., Pe?er G., Heikkinen R.K., Matsinos Y.G., Watts K.
and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
eco-evolutionary dynamics and species? responses to environmental changes.
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
	Species *pSpecies; Patch *pPatch; int spNum,nInds,nNonJuvs,nAdults; bool breeding;
};
struct disperser {
	Individual *pInd; Cell *pCell; bool yes;
};
struct zombie {
	Individual* pInd;
};
struct traitsums { // sums of trait genes for dispersal
	int ninds[gMaxNbSexes];				// no. of individuals
	double sumD0[gMaxNbSexes];			// sum of maximum emigration probability
	double ssqD0[gMaxNbSexes];			// sum of squares of maximum emigration probability
	double sumAlpha[gMaxNbSexes];	// sum of slope of emigration dens-dep reaction norm
	double ssqAlpha[gMaxNbSexes];	// sum of squares of slope of emigration den-dep reaction norm
	double sumBeta[gMaxNbSexes]; 	// sum of inflection point of emigration reaction norm
	double ssqBeta[gMaxNbSexes]; 	// sum of squares of inflection point of emigration reaction norm
	double sumDist1[gMaxNbSexes]; 	// sum of kernel I mean
	double ssqDist1[gMaxNbSexes]; 	// sum of squares of kernel I mean
	double sumDist2[gMaxNbSexes]; 	// sum of kernel II mean
	double ssqDist2[gMaxNbSexes]; 	// sum of squares of kernel II mean
	double sumProp1[gMaxNbSexes]; 	// sum of propn using kernel I
	double ssqProp1[gMaxNbSexes]; 	// sum of squares of propn using kernel I
	double sumDP[gMaxNbSexes]; 		// sum of SMS directional persistence
	double ssqDP[gMaxNbSexes]; 		// sum of squares of SMS directional persistence
	double sumGB[gMaxNbSexes]; 		// sum of SMS goal bias
	double ssqGB[gMaxNbSexes]; 		// sum of squares of SMS goal bias
	double sumAlphaDB[gMaxNbSexes];	// sum of SMS dispersal bias decay rate
	double ssqAlphaDB[gMaxNbSexes]; 	// sum of squares of SMS dispersal bias decay rate
	double sumBetaDB[gMaxNbSexes];		// sum of SMS dispersal bias decay infl. pt.
	double ssqBetaDB[gMaxNbSexes]; 	// sum of squares of SMS dispersal bias decay infl. pt.
	double sumStepL[gMaxNbSexes]; 	// sum of CRW step length
	double ssqStepL[gMaxNbSexes]; 	// sum of squares of CRW step length
	double sumRho[gMaxNbSexes]; 		// sum of CRW correlation coefficient
	double ssqRho[gMaxNbSexes]; 		// sum of squares of CRW correlation coefficient
	double sumS0[gMaxNbSexes];			// sum of maximum settlement probability
	double ssqS0[gMaxNbSexes];			// sum of squares of maximum settlement probability
	double sumAlphaS[gMaxNbSexes];	// sum of slope of settlement den-dep reaction norm
	double ssqAlphaS[gMaxNbSexes];	// sum of squares of slope of settlement den-dep reaction norm
	double sumBetaS[gMaxNbSexes]; 	// sum of inflection point of settlement reaction norm
	double ssqBetaS[gMaxNbSexes]; 	// sum of squares of inflection point of settlement reaction norm
	double sumGeneticFitness[gMaxNbSexes];
	double ssqGeneticFitness[gMaxNbSexes];
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
	vector<Individual*> getIndividualsInStage(int stage);
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
	void outputGeneValues(ofstream& ofsGenes, const int& yr, const int& gen) const;
	void clean(void); // Remove zero pointers to dead or dispersed individuals

	void updatePopNeutralTables();
	double getAlleleFrequency(int locus, int allele);
	int getAlleleTally(int locus, int allele);
	int getHeteroTally(int locus, int allele);
	int countHeterozygoteLoci();
	vector<int> countNbHeterozygotesEachLocus();
	double computeHs();
	std::vector<Individual*> getIndsWithCharacteristics( // Return a set of individuals with specified characteristics
		int,	// min age
		int,    // max age
		int,    // stage
		int     //sex
	);
	void cleanSampledInds(
	    Individual* // individual to remove from sampled individuals vector
	); // clean sampled individuals vector

	int sampleIndividuals( // Select a set of individuals with specified characteristics; return the number of individuals with those characteristics
	// void sampleIndividuals( // Select a set of individuals with specified characteristics; return the number of individuals with those characteristics
	        int, //number of individuals to sample
        	int,	// min age (0 if not set)
        	int,    // max age (max age if not set)
        	int,    // stage
        	int     //sex
	);

	Individual* catchIndividual(
	    double, // catching rate
	    int
	);

	// void completeTranslocation(
	//         std::vector <Individual*> // catched individuals
	// );

	// void recruitTranslocated(
	//         Individual*
	// );

	bool getSizeSampledInds(
	);

#ifndef NDEBUG
	// Testing only
	void clearInds() { inds.clear(); } // empty inds vector to avoid deallocating individual is used separately in test
#endif // NDEBUG

private:
	short nStages;
	short nSexes;
	Species *pSpecies;	// pointer to the species
	Patch *pPatch;			// pointer to the patch
	int nInds[gMaxNbStages][gMaxNbSexes];		// no. of individuals in each stage/sex

	vector <Individual*> inds; // all individuals in population except ...
	vector <Individual*> juvs; // ... juveniles until reproduction of ALL species
																	// has been completed

	vector<Individual*> sampledInds;
	//std::vector <Individual*> sampledInds; // individuals with specified characteristics from translocation!!! 
	vector<NeutralCountsTable> popNeutralCountTables;
	void resetPopNeutralTables();
};

//---------------------------------------------------------------------------

extern paramGrad *paramsGrad;
extern paramStoch *paramsStoch;
extern paramInit *paramsInit;
extern paramSim *paramsSim;
extern RSrandom *pRandom;

//---------------------------------------------------------------------------
#endif

