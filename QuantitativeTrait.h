/*----------------------------------------------------------------------------
 *
 *	Copyright (C) 2026 Greta Bocedi, Stephen C.F. Palmer, Justin M.J. Travis, Anne-Kathleen Malchow, Roslyn Henry, Théo Pannetier, Jette Wolff, Damaris Zurell
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

#ifndef TRAITTH
#define TRAITTH

#include <vector>
#include <set>
#include <memory>

#include "Parameters.h"
#include "Allele.h"
#include "SpeciesTrait.h"

using namespace std;

// Base interface for all genetic traits
class QuantitativeTrait {
public:
    virtual void mutate() = 0;
    virtual unique_ptr<QuantitativeTrait> clone() const = 0; //copies parameters (if not static) not gene seqeunces
    virtual void inheritGenes(const bool&, QuantitativeTrait*, set<unsigned int> const& , int ) = 0;
    virtual int getNLoci() const = 0;
    virtual float getMutationRate() const = 0;
    virtual bool isInherited() const = 0;
    virtual float getAlleleValueAtLocus(short chromosome, int i) const = 0;
    virtual float getDomCoefAtLocus(short whichChromosome, int position) const = 0;
    virtual float express() = 0;
    virtual ~QuantitativeTrait() { }
};

extern RSrandom* pRandom;

#endif
