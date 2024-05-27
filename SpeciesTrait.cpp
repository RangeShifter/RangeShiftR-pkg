
#include "SpeciesTrait.h"

// Species trait constructor
SpeciesTrait::SpeciesTrait(
	const TraitType& trType, const sex_t& sx,
	const set<int>& pos, const ExpressionType& expr,
	const DistributionType& initDist, const map<GenParamType, float> initParams,
	const DistributionType& dominanceDist, const map<GenParamType, float> dominanceParams,
	bool isInherited, const float& mutRate,
	const DistributionType& mutationDist, const map<GenParamType, float> mutationParams,
	Species* pSpecies) :
	traitType{ trType },
	sex{ sx },
	genePositions{ pos },
	expressionType{ expr },
	initialDistribution{ initDist },
	initialParameters{ initParams },
	dominanceDistribution{ dominanceDist },
	dominanceParameters{ dominanceParams },
	inherited{ isInherited },
	mutationDistribution{ mutationDist },
	mutationParameters{ mutationParams },
	mutationRate{ mutRate }
{
	// Initialise ploidy only once per species
	if (ploidy == NULL) this->ploidy = pSpecies->isDiploid() ? 2 : 1;

	// Check distribution parameters
	// Initial distribution
	for (auto [paramType, paramVal] : initParams) {
		switch (paramType)
		{
		case MIN: case MAX: case MEAN:
			if (!isValidTraitVal(paramVal))
				throw logic_error("Invalid parameter value: initial parameter " + to_string(paramType) + " must have a valid value for trait" + to_string(traitType) + ".");
			break;
		case SD:
			if (paramVal <= 0.0)
				throw logic_error("Invalid parameter value: initial parameter " + to_string(paramType) + " must be strictly positive");
			break;
		default:
			break;
		}
	}

	// Mutation distribution
	for (auto [paramType, paramVal] : mutationParams) {
		switch (paramType)
		{
		case MIN: case MAX: case MEAN:
			if (!isValidTraitVal(paramVal))
				throw logic_error("Invalid parameter value: mutation parameter " + to_string(paramType) + " must have a valid value for trait" + to_string(traitType) + ".");
			break;
		case SD: case SHAPE: case SCALE:
			if (paramVal <= 0.0)
				throw logic_error("Invalid parameter value: mutation parameter " + to_string(paramType) + " must be strictly positive");
			break;
		default:
			break;
		}
	}

	// Dominance distribution
	for (auto [paramType, paramVal] : dominanceParams) {
		switch (paramType)
		{
		case MIN: case MAX: case MEAN:
			if (paramVal < 0.0)
				throw logic_error("Invalid parameter value: dominance parameter " + to_string(paramType) + " must not be negative.");
			break;
		case SD: case SHAPE: case SCALE:
			if (paramVal <= 0.0)
				throw logic_error("Invalid parameter value: dominance parameter " + to_string(paramType) + " must be strictly positive");
			break;
		default:
			break;
		}
	}
}

bool SpeciesTrait::isValidTraitVal(const float& val) const {
	switch (traitType)
	{
	// Neutral trait
	case NEUTRAL: // only need to check for input parameters
	{
		return val >= 0.0 && val <= 255.0;
	}
	// Genetic Load
	case GENETIC_LOAD: case GENETIC_LOAD1: case GENETIC_LOAD2: case GENETIC_LOAD3: case GENETIC_LOAD4: case GENETIC_LOAD5:
	{
		return val >= -1.0 // genetic fitness traits can be beneficial
			&& val <= 1.0;
		break;
	}
	// Dispersal traits
	/// Emigration
	case E_D0_F: case E_D0_M: case E_D0: {
		return val >= 0.0 && val <= 1.0; // is a probability
		break;
	}
	case E_ALPHA_F: case E_ALPHA_M: case E_ALPHA:
	{
		return val > 0.0;
		break;
	}
	case E_BETA_F: case E_BETA_M: case E_BETA:
	{
		return true; // inflexion point can be any value
		break;
	}
	/// Settlement
	case S_S0_F: case S_S0_M: case S_S0:
	{
		return val >= 0.0 && val <= 1.0;
		break;
	}
	case S_ALPHA_F: case S_ALPHA_M: case S_ALPHA:
	{
		return val > 0.0;
		break;
	}
	case S_BETA_F: case S_BETA_M: case S_BETA:
	{
		return true;
		break;
	}
	/// Transfer - Kernels
	case KERNEL_MEANDIST_1_F: case KERNEL_MEANDIST_1_M: case KERNEL_MEANDIST_1:
	case KERNEL_MEANDIST_2_F: case KERNEL_MEANDIST_2_M: case KERNEL_MEANDIST_2:
	{
		return val >= 0.0; // is a distance
		break;
	}
	case KERNEL_PROBABILITY_F: case KERNEL_PROBABILITY_M: case KERNEL_PROBABILITY:
	{
		return val >= 0.0 && val <= 1.0;
		break;
	}
	/// Transfer - Correlated random walk
	case CRW_STEPLENGTH:
	{
		return val >= 0.0;
		break;
	}
	case CRW_STEPCORRELATION:
	{
		return val >= 0.0 && val <= 1.0;
		break;
	}
	/// Transfer - Stochastic Movement Simulator
	case SMS_DP:
	{
		return val >= 1.0; // according to parameter doc
		break;
	}
	case SMS_GB:
	{
		return val >= 1.0; // according to parameter doc
		break;
	}
	case SMS_ALPHADB:
	{
		return val > 0.0;
		break;
	}
	case SMS_BETADB:
	{
		return true;
		break;
	}
	default:
		throw logic_error("Invalid trait type " + to_string(traitType) + " passed to isValidTraitVal().");
		break;
	}
}

