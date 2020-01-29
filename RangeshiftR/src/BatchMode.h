/*------------------------------------------------------------------------------

RangeShifter v2.0 BatchMode

Functions for running in BATCH MODE

For full details of RangeShifter, please see:
Bocedi G., Palmer S.C.F., Peer G., Heikkinen R.K., Matsinos Y.G., Watts K.
and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
eco-evolutionary dynamics and species responses to environmental changes.
Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

Authors: Greta Bocedi & Steve Palmer, University of Aberdeen

Last updated: 7 January 2020 by Steve Palmer

------------------------------------------------------------------------------*/

#ifndef BatchModeH
#define BatchModeH

#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
using namespace std;

#include "Parameters.h"
#include "Landscape.h"
#include "Species.h"
#include "Model.h"
#if RS_ABC
#include "ABC.h"
#endif
#if RS_CONTAIN
#include "Control.h"
#endif // RS_CONTAIN 

struct batchfiles {
	bool ok;
	int batchNum;
	int nSimuls;
	int nLandscapes;
	int patchmodel,resolution,landtype,maxNhab,speciesdist,distresolution;
	int reproductn;
#if SEASONAL
	int nseasons;
#else
	int repseasons;
#endif // SEASONAL  
	int stagestruct,stages,transfer;
	int sexesDem;		// no. of explicit sexes for demographic model
	int sexesDisp;	// no. of explicit sexes for dispersal model
#if SEASONAL
//	string seasonFile;
#endif // SEASONAL 
	string parameterFile;
	string landFile;
	string stageStructFile;
	string emigrationFile;
	string transferFile;
	string settleFile;
	string geneticsFile;
#if RS_CONTAIN
	string manageFile;
#endif // RS_CONTAIN 
	string initFile;
#if VIRTUALECOLOGIST
	string virtEcolFile;
#endif
#if RS_ABC
	string abcParamsFile,abcObsFile;
#endif
};

struct simCheck {
	bool newsimul;
	int simul,simlines,reqdsimlines,errors;
};

batchfiles ParseControlFile(string,string,string);
#if BUTTERFLYDISP
int ParseParameterFile(string);
#else
int ParseParameterFile(void);
#endif
int ParseLandFile(int,string);
int ParseDynamicFile(string);
int ParseStageFile(string);
int ParseTransitionFile(short,short);
int ParseWeightsFile(string);
#if RS_CONTAIN
int ParseHabDemFile(short,short,string);
int ParseManageFile(string);
#endif // RS_CONTAIN 
#if SEASONAL
int ParseSeasonFile(string);
//#if PARTMIGRN
int ParseExtremeFile(string);
//#endif // PARTMIGRN
#endif // SEASONAL
int ParseEmigFile(void);
int ParseTransferFile(string);
int ParseSettleFile(void);
int ParseGeneticsFile(string);
int ParseArchFile(void);
int ParseInitFile(string);
int ParseInitIndsFile(void);
#if VIRTUALECOLOGIST
int ParseVirtEcolFile(string);
int ParseSampleFile(void);
int ParsePatchFile(void);
#endif // VIRTUALECOLOGIST
#if RS_ABC
int ParseABCParamsFile(void);
int ParseABCObsFile(void);
#endif // RS_ABC
#if TEMPMORT
int ParseMortFile(void);
int ReadMortalities(string);
#endif // TEMPMORT 
int CheckCostRaster(string,string);
simCheck CheckStageSex(string,int,int,simCheck,int,int,int,int,int,bool,bool);

void BatchError(
    string,	// file name
    int,		// line number
    int,		// option
    string	// fieldname
);
/* Write error message to batch log file. Options are as follows:
0 - general message only, no reference to field name
1 - fieldname must be 0 or 1
2 - fieldname must be 0, 1 or 2
3 - fieldname must be 0, 1, 2 or 3
4 - fieldname must be from 0 to 4
5 - fieldname must be from 0 to 5
6 - fieldname must be from 0 to 6
7 - fieldname must be from 0 to 7
10 - fieldname must be >0
11 - fieldname must be >=1
12 - fieldname must be >=2
13 - fieldname must be >=3
18 - fieldname must be >1
19 - fieldname must be >=0
20 - fieldname must be between 0 and 1
21 - fieldname must be >1
33 - fieldname must be 1, 2 or 3
44 - fieldname must be from 1 to 4
55 - fieldname must be from 1 to 5
66 - fieldname must be from 1 to 6
100 - fieldname must be between 0 and 100
111 - simulation must match first simulation in ParameterFile
222 - simulations must be sequential integers
223 - seasons must be sequential integers
333 - columns must match no. of habitats
444 - columns must be one fewer than no. of stages
555 - columns must match no. of stages
666 - fieldname must be a unique positive integer
*/
void BatchError(
    string,	// file name
    int,		// line number
    int,		// option
    string,	// fieldname
    string	// fieldname2
);
/* Write error message to batch log file. Options are as follows:
1 - fieldname must be greater than fieldname2
2 - fieldname must be greater than or equal to fieldname2
3 - fieldname must be less than or equal to fieldname2
4 - fieldname must be less than fieldname2
*/

int power2check(int x);

void CtrlFormatError(void);
void ArchFormatError(void);
#if VIRTUALECOLOGIST
void SampleFormatError(void);
#endif // VIRTUALECOLOGIST
void FormatError(string,int);
void OpenError(string,string);
void EOFerror(string);
void FileOK(string,int,int);
void FileHeadersOK(string);
void SimulnCountError(string);

void RunBatch(int,int);
int ReadParameters(int,Landscape*);
int ReadLandFile(int);
int ReadLandFile(int,Landscape*);
int ReadDynLandFile(Landscape*);
#if SEASONAL
int ReadStageStructure(int,Landscape*);
#else
int ReadStageStructure(int);
#endif // SEASONAL   
#if RS_CONTAIN
int ReadHabDemFile(const short,const short);
int ReadManageFile(int,Landscape*);
#endif // RS_CONTAIN 
int ReadTransitionMatrix(
    short,	// no. of stages
    short,	// no. of sexes represented for demography
    short,	// habitat index
    short		// season
);
int ReadStageWeights(int);
int ReadEmigration(int);
int ReadTransfer(int,Landscape*);
int ReadSettlement(int);
int ReadGenetics(int);
int ReadArchFile(string);
int ReadInitialisation(int,Landscape*);
int ReadInitIndsFile(int,Landscape*,string);
#if SEASONAL
int ReadSeasonFile(const short,const short);
//#if PARTMIGRN
int ReadExtremeFile(Landscape*,const short);
//#endif // PARTMIGRN
#endif // SEASONAL
#if VIRTUALECOLOGIST
int ReadVirtEcol(int);
int ReadSampleFile(string);
int ReadPatchFile(string);
#endif // VIRTUALECOLOGIST

#if RSDEBUG
extern ofstream DEBUGLOG;
#endif

// external pointers to parameter sets
extern paramGrad *paramsGrad;
extern paramStoch *paramsStoch;
extern paramInit *paramsInit;
extern paramSim *paramsSim;

extern Species *pSpecies;
extern string costmapname;	// see FormMove.cpp (VCL) OR Main.cpp (batch)
extern string genfilename;	// see FormGenetics.cpp (VCL) OR Main.cpp (batch)
#if RS_CONTAIN
extern Cull *pCull;
#endif // RS_CONTAIN 
#if VIRTUALECOLOGIST
extern string locfilename;		// see FormVirtEcol.cpp (VCL) OR Main.cpp (batch)
extern string patchfilename;	// see [NOT YET CODED FOR GUI] (VCL) OR Main.cpp (batch)
#endif // VIRTUALECOLOGIST 
#if TEMPMORT
extern string mortfilename;	// see [NOT YET CODED FOR GUI] (VCL) OR Main.cpp (batch)
#endif // TEMPMORT
#if !CLUSTER || RS_RCPP
extern int RS_random_seed;			// see RSrandom.cpp
#if RS_RCPP
void EOFerrorR(string);
void StreamErrorR(string);
#endif // RS_RCPP
#endif // !CLUSTER || RS_RCPP

//---------------------------------------------------------------------------
#endif
