//---------------------------------------------------------------------------

#pragma hdrstop

#include "RandomCheck.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)

ifstream inRandom;
ofstream outRandom;
ofstream outBernoulli;
ofstream outNormal;
ofstream outPoisson;
ofstream outIRandom;

void randomCheck(void) {

int samplesize,irandMin,irandMax;
double bernMean,normMean,normSD,poisMean;
string name,header;
simParams sim = paramsSim->getSim();

name = paramsSim->getDir(1) + "RandomCheck.txt";
inRandom.open(name.c_str());
if (!inRandom.is_open()) {
	cout << endl << "***** Error opening input file RandomCheck.txt" << endl;
	inRandom.clear();
	return;
}
for (int i = 0; i < 7; i++) {
	inRandom >> header;
}
inRandom >> samplesize >> bernMean >> normMean >> normSD >> poisMean >> irandMin >> irandMax;

name = paramsSim->getDir(2) + "Random.txt";
outRandom.open(name.c_str());
name = paramsSim->getDir(2) + "Bernoulli.txt";
outBernoulli.open(name.c_str());
name = paramsSim->getDir(2) + "Normal.txt";
outNormal.open(name.c_str());
name = paramsSim->getDir(2) + "Poisson.txt";
outPoisson.open(name.c_str());
name = paramsSim->getDir(2) + "IRandom.txt";
outIRandom.open(name.c_str());

for (int i = 0; i < samplesize; i++) {
	outRandom << pRandom->Random() << endl;
	outBernoulli << pRandom->Bernoulli(bernMean) << endl;
	outNormal << pRandom->Normal(normMean,normSD) << endl;
	outPoisson << pRandom->Poisson(poisMean) << endl;
	outIRandom << pRandom->IRandom(irandMin,irandMax) << endl;
}

inRandom.close(); inRandom.clear();
outRandom.close(); outRandom.clear();
outBernoulli.close(); outBernoulli.clear();
outNormal.close(); outNormal.clear();
outPoisson.close(); outPoisson.clear();
outIRandom.close(); outIRandom.clear();

}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


