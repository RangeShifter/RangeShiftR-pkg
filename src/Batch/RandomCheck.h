//---------------------------------------------------------------------------

#ifndef RandomCheckH
#define RandomCheckH

#include <fstream>
using namespace std;

#include "Version.h"
#include "Parameters.h"
#include "RSrandom.h"

void randomCheck(void);

extern paramSim *paramsSim;
#if RSRANDOM
extern RSrandom *pRandom;
#else
extern StochasticLib1 *pRandom;
#endif

//---------------------------------------------------------------------------
#endif
