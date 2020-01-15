/*------------------------------------------------------------------------------

RangeShifter v2.0 Main

Entry level function for BATCH MODE version

For compilation in Embarcadero

For full details of RangeShifter, please see:
Bocedi G., Palmer S.C.F., Peer G., Heikkinen R.K., Matsinos Y.G., Watts K.
and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
eco-evolutionary dynamics and species responses to environmental changes.
Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

Author: Steve Palmer, University of Aberdeen

Last updated: 6 January 2020 by Steve Palmer

------------------------------------------------------------------------------*/

#pragma hdrstop
#pragma argsused

//#include <tchar.h>
#include <string>
#include <stdio.h>
//#include <iostream.h>
#include <sstream>
#include <iostream>
//#include <io.h>
#include <iomanip>
#include <stdlib.h>
using namespace std;

#include "Parameters.h"
#include "Landscape.h"
#include "Species.h"
#include "SubCommunity.h"
#include "BatchMode.h"
#include "RandomCheck.h"

const string Int2Str(const int x)
{
	ostringstream o;
	if (!(o << x)) return "ERROR";
	return o.str();
}
const string Float2Str(const float x) {
	ostringstream o;
	if (!(o << x)) return "ERROR";
	return o.str();
}
const string Double2Str(const double x) {
	ostringstream o;
	if (!(o << x)) return "ERROR";
	return o.str();
}

void MemoLine(string msg) {
// dummy function for batch version
}

#if RSDEBUG
void DebugGUI(string msg) {
// dummy function for batch version
}
#endif

string habmapname,patchmapname,distnmapname;	// req'd for compilation, but not used
string costmapname,genfilename;					 			// ditto
vector <string> hfnames;											// ditto

paramGrad *paramsGrad;			// pointer to environmental gradient parameters
paramStoch *paramsStoch;		// pointer to environmental stochasticity parameters
paramInit *paramsInit;			// pointer to initialisation parameters
paramSim *paramsSim;				// pointer to simulation parameters

Species *pSpecies;  				// pointer to species
Community *pComm;						// pointer to community
RSrandom *pRandom;          // pointer to random number routines

#if RSDEBUG
ofstream DEBUGLOG;
ofstream MUTNLOG;
#endif

//---------------------------------------------------------------------------
#if CLUSTER || RCPP
int main(int argc, char* argv[])
#else
int _tmain(int argc, _TCHAR* argv[])
#endif
{

//int i,t0,t1,Nruns;
int i,t0,t1;
int nSimuls, nLandscapes; // no. of simulations and landscapes in batch

t0 = time(0);

// set up parameter objects
paramsGrad = new paramGrad;
paramsStoch = new paramStoch;
paramsInit = new paramInit;
paramsSim = new paramSim;                

// set up working directory and control file name
string cname;
#if CLUSTER || RCPP
if (argc > 1) {
	// full path name of directory passed as a parameter
	paramsSim->setDir(argv[1]);
	if (argc > 2) {
		// control file number also passed as a parameter
		int i = atoi(argv[2]);
		cname  = paramsSim->getDir(0) + "Inputs/CONTROL" + Int2Str(i) + ".txt";
	}
	else {
		// default name is CONTROL.txt
		cname  = paramsSim->getDir(0) + "Inputs/CONTROL.txt";
	}
}
else {
	// use current directory - get name from first (automatic) parameter
	string nameS = argv[0];
	string path  = argv[0];
	unsigned int loc = nameS.find("/",0);
	while(loc < 999999) {
		nameS = nameS.substr(loc+1);
		loc = nameS.find("/",0);
	}
	path = path.substr(0,path.length()-nameS.length());
	paramsSim->setDir(path);
	// control file name is forced to be CONTROL.txt
	cname  = paramsSim->getDir(0) + "Inputs/CONTROL.txt";
}
#else
if (_argc > 1) {
	// full path name of directory passed as a parameter
	paramsSim->setDir(_argv[1]);
	if (_argc > 2) {
		// control file name also passed as a parameter
		cname  = paramsSim->getDir(0) + "Inputs\\" + _argv[2];
	}
	else {
		// default name is CONTROL.txt
		cname  = paramsSim->getDir(0) + "Inputs\\CONTROL.txt";
	}
}
else {
	// use current directory - get name from first (automatic) parameter
	string nameS = _argv[0];
	string path  = _argv[0];
	unsigned int loc = nameS.find("\\",0);
	while(loc < 999999) {
		nameS = nameS.substr(loc+1);
		loc = nameS.find("\\",0);
	}
	path = path.substr(0,path.length()-nameS.length());
	paramsSim->setDir(path);
	// control file name is forced to be CONTROL.txt
	cname  = paramsSim->getDir(0) + "Inputs\\CONTROL.txt";
}
#endif
#if RSDEBUG
cout << endl << "Working directory: " << paramsSim->getDir(0) << endl;
cout << endl << "Inputs folder:     " << paramsSim->getDir(1) << endl;
cout << endl << "Control file:      " << cname << endl << endl;
#endif

bool errorfolder = CheckDirectory();
if (errorfolder) {
	cout << endl << "***** Invalid working directory: " << paramsSim->getDir(0)
		<< endl << endl;
	cout << "***** Working directory must contain Inputs, Outputs and Output_Maps folders"
		<< endl << endl;
	cout << "*****" << endl;
	cout << "***** Simulation ABORTED - enter any number to terminate program" << endl;
	cout << "*****" << endl;
	cin >> i;
	return 666;
}

#if RSDEBUG
// set up debugging log file
string name = paramsSim->getDir(2) + "DebugLog.txt";
DEBUGLOG.open(name.c_str());
name = paramsSim->getDir(2) + "MutnLog.txt";
MUTNLOG.open(name.c_str());
//DEBUGLOG << "Main(): random integers:";
//for (int i = 0; i < 5; i++) {
//	int rrrr = pRandom->IRandom(1000,2000); DEBUGLOG << " " << rrrr;
//}
//DEBUGLOG << endl;
//DEBUGLOG << "Main(): paramsSim = " << paramsSim << endl;
if (DEBUGLOG.is_open())
	cout << endl << "Main(): DEBUGLOG is open" << endl << endl;
else
	cout << endl << "Main(): DEBUGLOG is NOT open" << endl << endl;
#endif

/*
for (int i = 0; i < 10; i++) {
//	DEBUGLOG << pRandom->Random() << endl;
//	DEBUGLOG << pRandom->IRandom(5,55) << endl;
//	DEBUGLOG << pRandom->Poisson(4.2) << endl;
//	DEBUGLOG << pRandom->Bernoulli(0.6045) << endl;
	DEBUGLOG << pRandom->Normal(-564.7,123.4) << endl;
}
*/

/*

DEBUGLOG << endl << "Random():" << endl;
for (int i = 0; i < 5; i++) {
	for (int j = 0; j < 10; j++) {
		DEBUGLOG << pRandom->Random() << " ";
	}
	DEBUGLOG << endl;
}
DEBUGLOG << endl << "IRandom(5,55):" << endl;
for (int i = 0; i < 5; i++) {
	for (int j = 0; j < 50; j++) {
		DEBUGLOG << pRandom->IRandom(5,55) << " ";
	}
	DEBUGLOG << endl;
}
DEBUGLOG << endl << "Poisson(4.2):" << endl;
for (int i = 0; i < 5; i++) {
	for (int j = 0; j < 10; j++) {
		DEBUGLOG << pRandom->Poisson(4.2) << " ";
	}
	DEBUGLOG << endl;
}
DEBUGLOG << endl << "Bernoulli(0.6):" << endl;
for (int i = 0; i < 5; i++) {
	for (int j = 0; j < 20; j++) {
		DEBUGLOG << pRandom->Bernoulli(0.6) << " ";
	}
	DEBUGLOG << endl;
}
DEBUGLOG << endl << "Normal(0.0,1.0):" << endl;
for (int i = 0; i < 5; i++) {
	for (int j = 0; j < 10; j++) {
		DEBUGLOG << pRandom->Normal(0.0,1.0) << " ";
	}
	DEBUGLOG << endl;
}
DEBUGLOG << endl << "Normal(2.5,0.35):" << endl;
for (int i = 0; i < 5; i++) {
	for (int j = 0; j < 10; j++) {
		DEBUGLOG << pRandom->Normal(2.5,0.35) << " ";
	}
	DEBUGLOG << endl;
}
DEBUGLOG << endl << "Normal(-564.7,123.4):" << endl;
for (int i = 0; i < 5; i++) {
	for (int j = 0; j < 10; j++) {
		DEBUGLOG << pRandom->Normal(-564.7,123.4) << " ";
	}
	DEBUGLOG << endl;
}

*/

/*
DEBUGLOG.close();
DEBUGLOG.clear();

cout << "*****" << endl;
cout << "***** Simulation completed - enter any number to terminate program" << endl;
cout << "*****" << endl;
cin >> i;

return 0;
*/


// set up species
// FOR MULTI-SPECIES MODEL, THERE WILL BE AN ARRAY OF SPECIES POINTERS
// OR A COMMUNITY CLASS TO HOLD THE SPECIES
pSpecies = new Species;
demogrParams dem = pSpecies->getDemogr();
stageParams sstruct = pSpecies->getStage();
trfrRules trfr = pSpecies->getTrfr();

batchfiles b;
string indir  = paramsSim->getDir(1);
string outdir = paramsSim->getDir(2);
b = ParseControlFile(cname,indir,outdir);       
if (b.ok) { 
	nSimuls = b.nSimuls;
	nLandscapes = b.nLandscapes;
	dem.repType = b.reproductn;
	dem.repSeasons = b.repseasons;
	if (b.stagestruct == 0) dem.stageStruct = false; else dem.stageStruct = true;
	sstruct.nStages = b.stages;
	if (b.transfer == 0) trfr.moveModel = false;
	else {
		trfr.moveModel = true;
		trfr.moveType = b.transfer;
	}
	cout << endl << "Batch input files OK" << endl;
	pSpecies->setDemogr(dem);
	pSpecies->setStage(sstruct);
	pSpecies->setTrfr(trfr);
	simParams sim = paramsSim->getSim();
	sim.batchMode = true;
	sim.batchNum = b.batchNum;  
	paramsSim->setSim(sim);
}
else {
	cout << endl << "Error in parsing batch input files - see BatchLog file for details" << endl;
}
#if RSDEBUG
DEBUGLOG << "Main(): dem.repType = " << dem.repType << endl;
#endif

// set up random number class
#if RCPP
#if RSDEBUG
pRandom = new RSrandom(666);
#else
pRandom = new RSrandom(-1);  // need to be replaced with parameter from control file
#endif
#else
pRandom = new RSrandom();
#endif


#if RANDOMCHECK
randomCheck();
#else
if (b.ok) {
	RunBatch(nSimuls,nLandscapes);
}
#endif

#if RSDEBUG
if (DEBUGLOG.is_open()) {
	DEBUGLOG.close(); DEBUGLOG.clear();
}
if (MUTNLOG.is_open()) {
	MUTNLOG.close(); MUTNLOG.clear();
}
#endif

t1 = time(0);
cout << endl << "***** Elapsed time " << t1-t0 << " seconds" << endl << endl;

cout << "*****" << endl;
cout << "***** Simulation completed - enter any number to terminate program" << endl;
cout << "*****" << endl;
cin >> i;

return 0;
}

//---------------------------------------------------------------------------

// Dummy functions corresponding to those used in GUI version

/* Batch mode of v2.0 currently has no facility to save maps (unless initiated from GUI).
To do so, we would need a form of bit map which is portable across platforms
and operating systems, rather than the Embarcadero VCL classes.
Does such exist?
*/

traitCanvas SetupTraitCanvas(void) {
traitCanvas tcanv;
for (int i = 0; i < NTRAITS; i++) { tcanv.pcanvas[i] = 0; }
return tcanv;
}

void Landscape::setLandMap(void) { }
void Landscape::drawLandscape(int rep,int yr,int landnum) { }
void Community::viewOccSuit(int year,double mn,double se) { }
void Community::draw(int rep,int yr,int gen,int landNum) { }

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

