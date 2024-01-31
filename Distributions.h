#ifndef DISTRIBUTIONSH
#define DISTRIBUTIONSH


#include <stdlib.h>
#include <fstream>
//#include <iostream>
using namespace std;

#include <cmath>
#include <random>
#if !LINUX_CLUSTER
#include <ctime>
#endif

#include "RSRandom.h"

struct distributionType {

	virtual ~distributionType() {}

};

struct NormalDistribution : distributionType {

	// Set up standard normal distribution
	std::normal_distribution<>* pNormal;
	
	NormalDistribution(float mean, float sd) {
		pNormal = new normal_distribution<double>(mean, sd);
	}

	float sample(void) {


	mt19937 gen = pRandom->getRNG();
	return pNormal->operator()(gen);
	}

	~NormalDistribution() {
		if (pNormal != 0) delete pNormal;
	}
};

struct NormalDistribution : distributionType {

	// Set up standard normal distribution
	std::normal_distribution<>* pNormal;

	NormalDistribution(float mean, float sd) {
		pNormal = new normal_distribution<double>(mean, sd);
	}

	float sample(void) {


		mt19937 gen = pRandom->getRNG();
		return pNormal->operator()(gen);
	}

	~NormalDistribution() {
		if (pNormal != 0) delete pNormal;
	}
};



extern RSrandom* pRandom;

#endif