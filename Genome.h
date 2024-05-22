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

#ifndef GenomeH
#define GenomeH

#include <vector>
#include <algorithm>

#include "Parameters.h"
#include "Species.h"

#define INTBASE 100.0; // to convert integer alleles into continuous traits

struct locus { short allele[2]; };
struct locusOK { short allele[2]; bool ok; };

//---------------------------------------------------------------------------

class Chromosome {

public:
	Chromosome(int);
	~Chromosome();
	short nLoci(void);
	double additive( // Return trait value on normalised genetic scale
		const bool	// diploid
	);
	double meanvalue( // Return trait value on normalised genetic scale
		const bool	// diploid
	);
	double additive( // Return trait value on normalised genetic scale
		const short,	// locus
		const bool		// diploid
	);
	locus alleles(	// Return allele values at a specified locus
		const int			// position of locus on chromosome
	);
	void initialise( // Set up chromosome at simulation initialisation
		const double,	// normalised phenotypic trait value			
		const double,	// s.d. of allelic variance (genetic scale)
		const bool		// diploid
	);
	void initialise( // Set up specified locus at simulation initialisation
		const short,	// locus
		const short,	// position: 0 from mother, 1 from father
		const int			// allele value
	);
	void inherit( // Inherit chromosome from specified parent
		const Chromosome*,	// pointer to parent's chromosome
		const short,				// position: 0 from mother, 1 from father
		const short,				// no. of loci
		const double,				// mutation probability
		const double,				// crossover probability
		const double,				// s.d. of mutation magnitude (genetic scale)
		const bool					// diploid
	);

protected:

private:
	short nloci;
	locus* pLoci;

};

//---------------------------------------------------------------------------

class Genome {

public:
	Genome();
	Genome(int, int, bool);
	Genome(Species*);
	Genome(Species*, Genome*, Genome*);
	~Genome();
	void setGene( // Set up new gene at initialisation for 1 chromosome per trait
		const short,	// chromosome number
		const short,	// expression type (NOT CURRENTLY USED)
		const double,	// normalised trait value		
		const double		// s.d. of allelic variance
	);
	void setTrait( // Set up trait at initialisation for trait mapping
		Species*,			// pointer to Species
		const int,		// trait number			
		const double,	// normalised trait value		
		const double		// s.d. of allelic variance
	);
	void setNeutralLoci( // Set up neutral loci at initialisation
		Species*,			// pointer to Species
		const double		// s.d. of allelic variance
	);
	double express(
		// Return the expressed value of a gene when species has one chromosome per trait
		short,	// chromosome number
		short,	// expression type (NOT CURRENTLY USED)
		short		// individual's sex (NOT CURRENTLY USED)
	);
	double express(
		// Return the expressed value of a trait when genetic architecture is defined
		Species*,	// pointer to Species
		short			// trait number
	);
	locusOK getAlleles( // Get allele values at a specified locus
		short,	// chromosome number
		short		// locus position on chromosome
	);
	// SCFP NEW DECLARATIONS
	void setDiploid(bool);
	bool isDiploid(void);
	void inherit( // Inherit from specified parent
		const Genome*,	// pointer to parent's genome
		const short,		// position: 0 from mother, 1 from father
		const short,		// chromasome number
		const double,		// mutation probability
		const double,		// crossover probability
		const double			// s.d. of mutation magnitude (genetic scale)
	);
	short getNChromosomes(void);
	void outGenHeaders(
		const int,	// replicate
		const int,	// landscape number
		const bool	// output as cross table?
	);
	void outGenetics(
		const int,	// replicate
		const int,	// year
		const int,	// species number
		const int, 	// individual ID
		const bool 	// output as cross table?
	);


private:
	short nChromosomes;						// no. of chromosomes
	bool diploid;
	Chromosome** pChromosome;

};

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

extern paramSim* paramsSim;
extern RSrandom* pRandom;

#if RSDEBUG
extern ofstream DEBUGLOG;
extern ofstream MUTNLOG;
#endif

//---------------------------------------------------------------------------

#endif
