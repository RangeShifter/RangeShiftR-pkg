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

 RangeShifter v2.0 Individual

 Implements the Individual class

 Various optional attributes (genes for traits, movement parameters, etc.) are
 allocated dynamically and accessed by pointers if required.

 For full details of RangeShifter, please see:
 Bocedi G., Palmer S.C.F., Pe’er G., Heikkinen R.K., Matsinos Y.G., Watts K.
 and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
 eco-evolutionary dynamics and species’ responses to environmental changes.
 Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

 Authors: Greta Bocedi & Steve Palmer, University of Aberdeen

 Last updated: 28 July 2021 by Greta Bocedi

 ------------------------------------------------------------------------------*/

#ifndef IndividualH
#define IndividualH


#include <queue>
#include <algorithm>
#include <ranges>
#include <memory>
using namespace std;

#include "Parameters.h"
#include "Species.h"
#include "Landscape.h"
#include "Patch.h"
#include "Cell.h"
#include "TraitFactory.h"

//---------------------------------------------------------------------------

struct indStats {
	short stage; sex_t sex; short age; short status; short fallow;
	bool isDeveloping;
};
struct pathData { // to hold path data common to SMS and CRW models
	int year, total, out; // nos. of steps
	Patch* pSettPatch;		// pointer to most recent patch tested for settlement
	short settleStatus; 	// whether ind may settle in current patch
	// 0 = not set, 1 = debarred through density dependence rule
	// 2 = OK to settle subject to finding a mate
//	bool leftNatalPatch;	// individual has moved out of its natal patch
#if RS_RCPP
	short pathoutput;
#endif
};
struct pathSteps { // nos. of steps for movement model
	int year, total, out;
};
struct settlePatch {
	Patch* pSettPatch; short settleStatus;
};

struct trfrData {

	virtual void addMyself(trfrData& toAdd) = 0;

	virtual void clone(const trfrData& copyFrom) = 0;

	virtual void divideTraitsBy(int) = 0;

	virtual movement_t getType() = 0;

	virtual ~trfrData() {}

};

struct crwData : trfrData { // to hold data for CRW movement model

	float prevdrn;	// direction of previous step (UNITS)
	float xc, yc;		// continuous cell co-ordinates	
	float rho;			// phenotypic step correlation coefficient
	float stepLength; // phenotypic step length (m)
	//static bool straigtenPath; //does not vary between individuals, shared
	//static float stepMort; //does not vary between individuals, shared

	crwData(float prevdrnA, float xcA, float ycA) : prevdrn(prevdrnA), xc(xcA), yc(ycA), rho(0.0), stepLength(0.0) {}
	~crwData() {}

	void addMyself(trfrData& toAdd) {

		auto& CRW = dynamic_cast<crwData&>(toAdd);

		CRW.stepLength += stepLength;
		CRW.rho += rho;

		//stepLength += pCRW.stepLength;
	//	rho += pCRW.rho;
	}

	movement_t getType() { return CRW; }

	void clone(const trfrData& copyFrom) {


		const crwData& pCopy = dynamic_cast<const crwData&>(copyFrom);

		stepLength = pCopy.stepLength;
		rho = pCopy.rho;
	}

	void divideTraitsBy(int i) {

		stepLength /= i;
		rho /= i;
	}

};
struct array3x3d { double cell[3][3]; };
struct movedata { float dist; float cost; };
struct smsData : trfrData {
	locn prev;			// location of previous cell
	locn goal;			// location of goal
	float dp;				// directional persistence
	float gb;				// goal bias
	float alphaDB;	// dispersal bias decay rate
	int betaDB;			// dispersal bias decay inflection point (no. of steps)

	//below are shared
	//static short pr;
	//static short prMethod;
	//static short memSize;
	//static short goalType;
	//static float stepMort;
	//static bool straigtenPath;

	smsData(locn prevA, locn goalA) : prev(prevA), goal(goalA), dp(0.0), gb(0.0), alphaDB(0.0), betaDB(0) {}
	~smsData() {}


	void addMyself(trfrData& toAdd) {
		auto& SMS = dynamic_cast<smsData&>(toAdd);
		SMS.dp += dp;
		SMS.gb += gb;
		SMS.alphaDB += alphaDB;
		SMS.betaDB += betaDB;
	}

	movement_t getType() { return SMS; }

	void clone(const trfrData& copyFrom) {
		auto& pCopy = dynamic_cast<const smsData&>(copyFrom);
		dp = pCopy.dp;
		gb = pCopy.gb;
		alphaDB = pCopy.alphaDB;
		betaDB = pCopy.betaDB;
	}

	void divideTraitsBy(int i) {

		dp /= i;
		gb /= i;
		alphaDB /= i;
		betaDB /= i;
	}

};

struct kernelData : trfrData {
	float	meanDist1;
	float	meanDist2;
	float	probKern1;

	kernelData(float meanDist1A, float meanDist2A, float probKern1A) : meanDist1(meanDist1A), meanDist2(meanDist2A), probKern1(probKern1A) {}
	~kernelData() {}

	void addMyself(trfrData& toAdd) {

		auto& Kernel = dynamic_cast<kernelData&>(toAdd);

		Kernel.meanDist1 += meanDist1;
		Kernel.meanDist2 += meanDist2;
		Kernel.probKern1 += probKern1;
	}

	movement_t getType() { return KERNEL; }

	void clone(const trfrData& copyFrom) {
		const kernelData& pCopy = dynamic_cast<const kernelData&>(copyFrom);
		meanDist1 = pCopy.meanDist1;
		meanDist2 = pCopy.meanDist2;
		probKern1 = pCopy.probKern1;
	}

	void divideTraitsBy(int i) {
		meanDist1 /= i;
		meanDist2 /= i;
		probKern1 /= i;
	}
};

class Individual {

public:
	static int indCounter; // used to create ID, held by class, not members of class
	static TraitFactory traitFactory;
	Individual( // Individual constructor
		Cell*,	// pointer to Cell
		Patch*,	// pointer to patch
		short,	// stage
		short,	// age
		short,	// reproduction interval (no. of years/seasons between breeding attempts)
		float,	// probability that sex is male
		bool,		// TRUE for a movement model, FALSE for kernel-based transfer
		short		// movement type: 1 = SMS, 2 = CRW
	);
	~Individual(void);
	void setUpGenes( // Set genes for individual variation from species initialisation parameters
		Species*,			// pointer to Species
		int						// Landscape resolution
	);
	void inheritTraits( // Inherit genome from parents
		Species*,			// pointer to Species
		Individual*,	// pointer to mother
		Individual*,	// pointer to father (must be 0 for an asexual Species)
		int						// Landscape resolution
	);

	void inheritTraits(Species* pSpecies, Individual* mother, int resol); //haploid

	void setDispersalPhenotypes(Species* pSpecies, int resol);

	TTrait* getTrait(TraitType trait) const;

	set<TraitType> getTraitTypes();

	void inherit(Species* pSpecies, const Individual* mother, const Individual* father);

	void inherit(Species* pSpecies, const Individual* mother); // haploid

	Individual* traitClone(Cell*, Patch*, float, bool, short);

	void setEmigTraits(Species* pSpecies, bool sexDep, bool densityDep);
	void setTransferTraits(Species* pSpecies, transferRules trfr, int resol);

	emigTraits getEmigTraits(void); // Get phenotypic emigration traits

	void setKernelTraits(Species* pSpecies, bool sexDep, bool twinKernel, int resol);

	trfrKernelParams getKernTraits(void); // Get phenotypic transfer by kernel traits

	void setSMSTraits(Species* pSpecies);

	trfrSMSTraits getSMSTraits(void); // Get phenotypic transfer by SMS traits

	void setCRWTraits(Species* pSpecies);

	trfrCRWTraits getCRWTraits(void); // Get phenotypic transfer by CRW traits

	void setSettlementTraits(Species* pSpecies, bool sexDep);

	settleTraits getSettTraits(void); // Get phenotypic settlement traits

	trfrData* getTrfrData(void);
	void setEmigTraits(const emigTraits& emig);
	void setSettleTraits(const settleTraits& settle);


	// Identify whether an individual is a potentially breeding female -
	// if so, return her stage, otherwise return 0
	int breedingFem(void);
	int getId(void);
	sex_t getSex(void);
	int getStatus(void);
	float getGeneticFitness(void);
	indStats getStats(void);
	Cell* getLocn( // Return location (as pointer to Cell)
		const short	// option: 0 = get natal locn, 1 = get current locn
	); //
	Patch* getNatalPatch(void);
	void setYearSteps(int);
	pathSteps getSteps(void);
	settlePatch getSettPatch(void);
	void setSettPatch(const settlePatch);
	void setStatus(short);
	void developing(void);
	void develop(void);
	void ageIncrement( // Age by one year
		short	// maximum age - if exceeded, the Individual dies
	);
	void incFallow(void); // Inrement no. of reproductive seasons since last reproduction
	void resetFallow(void);
	void moveto( // Move to a specified neighbouring cell
		Cell*	// pointer to the new cell
	);
	// Move to a new cell by sampling a dispersal distance from a single or double
	// negative exponential kernel
	// Returns 1 if still dispersing (including having found a potential patch), otherwise 0
	int moveKernel(
		Landscape*,		// pointer to Landscape
		Species*,			// pointer to Species
		const short,	// reproduction type (see Species)
		const bool    // absorbing boundaries?
	);
	// Make a single movement step according to a mechanistic movement model
	// Returns 1 if still dispersing (including having found a potential patch), otherwise 0
	int moveStep(
		Landscape*,		// pointer to Landscape
		Species*,			// pointer to Species
		const short,	// landscape change index
		const bool    // absorbing boundaries?
	);
	movedata smsMove( // Move to a neighbouring cell according to the SMS algorithm
		Landscape*,		// pointer to Landscape
		Species*,			// pointer to Species
		const short,	// landscape change index
		const bool,		// TRUE if still in (or returned to) natal patch
		const bool,   // individual variability?
		const bool    // absorbing boundaries?
	);
	array3x3d getSimDir( // Weight neighbouring cells on basis of current movement direction
		const int,	// current x co-ordinate
		const int,	// current y co-ordinate
		const float	// directional persistence value
	);
	array3x3d getGoalBias( // Weight neighbouring cells on basis of goal bias
		const int,	// current x co-ordinate
		const int,	// current y co-ordinate
		const int,	// goal type: 0 = none, 1 = towards goal (NOT IMPLEMENTED), 2 = dispersal bias
		const float	// GOAL BIAS VALUE
	);
	array3x3d calcWeightings( // Calculate weightings for neighbouring cells
		const double,	// base for power-law (directional persistence or goal bias value)
		const double	// direction in which lowest (unit) weighting is to be applied
	);
	array3x3f getHabMatrix( // Weight neighbouring cells on basis of (habitat) costs
		Landscape*,		// pointer to Landscape
		Species*,			// pointer to Species
		const int,		// current x co-ordinate
		const int,		// current y co-ordinate
		const short,	// perceptual range (cells)
		const short,	// perceptual range evaluation method (see Species)
		const short,	// landscape change index
		const bool    // absorbing boundaries?
	);
#if RS_RCPP
	void outMovePath( // Write records to movement paths file
		const int		 	// year
	);
#endif

private:
	int indId;
	float geneticFitness;
	short stage;
	sex_t sex;
	short age;
	short status;	// 0 = initial status in natal patch / philopatric recruit
	// 1 = disperser
	// 2 = disperser awaiting settlement in possible suitable patch
	// 3 = waiting between dispersal events
	// 4 = completed settlement
	// 5 = completed settlement in a suitable neighbouring cell
	// 6 = died during transfer by failing to find a suitable patch
	//     (includes exceeding maximum number of steps or crossing
	//			absorbing boundary)
	// 7 = died during transfer by constant, step-dependent,
	//     habitat-dependent or distance-dependent mortality
	// 8 = failed to survive annual (demographic) mortality
	// 9 = exceeded maximum age
	short fallow; // reproductive seasons since last reproduction
	bool isDeveloping;
	Cell* pPrevCell;						// pointer to previous Cell
	Cell* pCurrCell;						// pointer to current Cell
	Patch* pNatalPatch;					// pointer to natal Patch
	pathData* path; 						// pointer to path data for movement model
	std::unique_ptr <emigTraits> pEmigTraits;			// pointer to emigration traits
	std::unique_ptr <settleTraits> pSettleTraits;		// pointer to settlement traits
	std::unique_ptr <trfrData> pTrfrData; //can be sms, kernel, crw
	std::queue <locn> memory;		// memory of last N squares visited for SMS
	map<TraitType, unique_ptr<TTrait>> spTraitTable;
};


//---------------------------------------------------------------------------

double cauchy(double location, double scale);
double wrpcauchy(double location, double rho = exp(double(-1)));

extern RSrandom* pRandom;

#if RSDEBUG
extern ofstream DEBUGLOG;
#endif

#if RS_RCPP
extern ofstream outMovePaths;
#endif

#if RSDEBUG
void testIndividual();
#endif

//---------------------------------------------------------------------------
#endif // IndividualH
