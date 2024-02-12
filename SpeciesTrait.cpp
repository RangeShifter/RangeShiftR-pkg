
#include "SpeciesTrait.h"

//could be handled in header file but here for now for flexibility
SpeciesTrait::SpeciesTrait(vector<string> parameters, Species* pSpecies) {

	bool neutralPresent = false;
	if (ploidy == 0) this->ploidy = pSpecies->isDiploid() ? 2 : 1;

	this->sex = stringToSex(parameters[2]);
	TraitType traitType = stringToTraitType(parameters[1], this->sex);
	this->positions = stringToLoci(parameters[3], parameters[4], pSpecies);
	this->expressionType = stringToExpressionType(parameters[5]);
	this->initialDistribution = stringToDistributionType(parameters[6]);
	this->initialParameters = stringToParameterMap(parameters[7]);
	this->dominanceDistribution = stringToDistributionType(parameters[8]);
	this->dominanceParameters = stringToParameterMap(parameters[9]);

	if (traitType == SNP || traitType == ADAPTIVE)
		this->inherited = true;
	else
		this->inherited = iequals(parameters[10], "true") ? true : false;

	if (this->isInherited()) {
		this->mutationDistribution = stringToDistributionType(parameters[11]);
		this->mutationParameters = stringToParameterMap(parameters[12]);
		this->mutationRate = stof(parameters[13]);
	}

	// error outputting for different traits 
	if (traitType == SNP) {
		if (mutationDistribution != KAM && mutationDistribution != SSM)
			cout << endl << "Traits file: ERROR - Neutral marker mutation distribution must be KAM or SSM (max = 256))" << endl;

		if (pSpecies->getNumberOfNeutralLoci() > 0)
			cout << endl << "Traits file: WARNING - can only have one set of neutral markers, overwriting previous" << endl;
		else pSpecies->setNumberOfNeutralLoci(positions.size());
	}
}

TraitType SpeciesTrait::stringToTraitType(const std::string& str, sex_t sex) const {

	if (sex == MAL) {
		if (str == "emigration_d0") return E_D0_M;
		else if (str == "emigration_alpha") return E_ALPHA_M;
		else if (str == "emigration_beta") return E_BETA_M;
		else if (str == "settlement_s0") return S_S0_M;
		else if (str == "settlement_alpha") return S_ALPHA_M;
		else if (str == "settlement_beta") return S_BETA_M;
		else if (str == "kernel_meanDistance1") return KERNEL_MEANDIST_1_M;
		else if (str == "kernel_meanDistance2") return KERNEL_MEANDIST_2_M;
		else if (str == "kernel_probability") return KERNEL_PROBABILITY_M;
		else if (str == "crw_stepLength") return CRW_STEPLENGTH_M;
		else if (str == "crw_stepCorrelation") return CRW_STEPCORRELATION_M;
		else throw logic_error(str + " is not a valid trait type.");
	} else {
		if (str == "emigration_d0") return E_D0_F;
		else if (str == "emigration_alpha") return E_ALPHA_F;
		else if (str == "emigration_beta") return E_BETA_F;
		else if (str == "settlement_s0") return S_S0_F;
		else if (str == "settlement_alpha") return S_ALPHA_F;
		else if (str == "settlement_beta") return S_BETA_F;
		else if (str == "kernel_meanDistance1") return KERNEL_MEANDIST_1_F;
		else if (str == "kernel_meanDistance2") return KERNEL_MEANDIST_2_F;
		else if (str == "kernel_probability") return KERNEL_PROBABILITY_F;
		else if (str == "crw_stepLength") return CRW_STEPLENGTH_F;
		else if (str == "crw_stepCorrelation") return CRW_STEPCORRELATION_F;
		else if (str == "sms_directionalPersistence") return SMS_DP;
		else if (str == "sms_goalBias") return SMS_GB;
		else if (str == "sms_alphaDB") return SMS_ALPHADB;
		else if (str == "sms_betaDB") return SMS_BETADB;
		else if (str == "neutral") return SNP;
		else if (str == "adaptive") return ADAPTIVE;
		else throw logic_error(str + " is not a valid trait type.");
	}
}

ExpressionType SpeciesTrait::stringToExpressionType(const std::string& str) const {
	if (str == "average") return AVERAGE;
	else if (str == "additive") return ADDITIVE;
	else if (str == "multiplicative") return MULTIPLICATIVE;
	else if (str == "#") return NEUTRAL;
	else throw logic_error(str + " is not a valid gene expression type.");
}

DistributionType SpeciesTrait::stringToDistributionType(const std::string& str) const
{
	if (str == "#") return NONE;
	else if (str == "uniform") return UNIFORM;
	else if (str == "normal") return NORMAL;
	else if (str == "gamma") return GAMMA;
	else if (str == "scaled") return SCALED;
	else if (str == "negExp") return NEGEXP;
	else if (str == "KAM") return KAM;
	else if (str == "SSM") return SSM;
	else throw logic_error(str + " is not a valid distribution type.");
}

map<parameter_t, float> SpeciesTrait::stringToParameterMap(string parameters) const {

	map<parameter_t, float> paramMap;
	if (parameters != "#") {
		parameters.erase(remove(parameters.begin(), parameters.end(), '\"'), parameters.end());
		stringstream ss(parameters);

		string value, valueWithin;
		while (std::getline(ss, value, ',')) {
			stringstream sss(value);
			vector<string> paramValue;
			while (std::getline(sss, valueWithin, '=')) {
				paramValue.push_back(valueWithin);
			}

			if (paramValue.size() == 2) {
				parameter_t parameterT = paramValue[0];
				float value = stof(paramValue[1]);
				paramMap.emplace(parameterT, value);
			}
			else
				cout << endl << "Traits file: ERROR - parameter values for a distribution missing, should be e.g. 'mean=0,standard_deviation=0.5' or if not applicable put #" << endl;
		}
	}
	return paramMap;
}

set<int> SpeciesTrait::selectRandomLociPositions(int nbLoci, Species* pSpecies) const {

	int genomeSize = pSpecies->getGenomeSize();
	set<int> positions;
	for (int i = 0; i < nbLoci; ++i)
		positions.insert(pRandom->IRandom(0, genomeSize));
	return positions;
}


set<int> SpeciesTrait::stringToLoci(string pos, string nLoci, Species* pSpecies) const {

	set<int> positions;

	if (pos != "random") {
		//stringstream ss(pos);

		//string value, valueWithin;
		//while (std::getline(ss, value, ',')) {
		//	stringstream sss(value);
		//	vector<int> minMax;
		//	while (std::getline(sss, valueWithin, '-')) {
		//		minMax.push_back(stoi(valueWithin));
		//	}
		//	if (minMax[0] > pSpecies->getGenomeSize() || minMax[1] > pSpecies->getGenomeSize()) {
		//		cout << endl << "Traits file: ERROR - trait positions must not exceed genome size" << endl;
		//	}
		//	else {
		//		if (minMax.size() > 1) {
		//			for (int i = minMax[0]; i < minMax[1] + 1; ++i) {
		//				positions.insert(i);
		//			}
		//		}
		//		else
		//			positions.insert(minMax[0]);
		//	}
		//}
		positions = convertStringToSet(pos);

		for (auto position : positions) {
			if (position > pSpecies->getGenomeSize())
				cout << endl << "Traits file: ERROR - trait positions " << position << " must not exceed genome size" << endl;
		}
	}
	else {
		positions = selectRandomLociPositions(stoi(nLoci), pSpecies);
	}
	return positions;
}
