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
#ifndef TRAITFACTORYH
#define TRAITFACTORYH

#include <memory>

#include "SpeciesTrait.h"
#include "NeutralTrait.h"
#include "DispersalTrait.h"
#include "GeneticFitnessTrait.h"

// Create handled pointers to a new trait
class TraitFactory
{
public:
	TraitFactory() {};

	unique_ptr<QuantitativeTrait> Create(const TraitType traitType, SpeciesTrait* pSpTrait)
	{
		if (traitType == NEUTRAL) {
			return make_unique<NeutralTrait>(pSpTrait);
		}
		else if (traitType == GENETIC_LOAD1 
			|| traitType == GENETIC_LOAD2 
			|| traitType == GENETIC_LOAD3 
			|| traitType == GENETIC_LOAD4 
			|| traitType == GENETIC_LOAD5) {
			return make_unique<GeneticFitnessTrait>(pSpTrait);
		}
		else {
			return make_unique<DispersalTrait>(pSpTrait);
		}
	}
};
#endif
