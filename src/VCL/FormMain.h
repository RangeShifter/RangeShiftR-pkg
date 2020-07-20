/*------------------------------------------------------------------------------

RangeShifter v2.0 FormMain

Handles the entry-level form for the GUI version of the program.

Also contains methods for the Landscape, Patch, Community, Sub-community and
Individual classes to display output on the screen.

A random number stream is set up, which is available to all other units of the
program through the pointer pRandom.

For full details of RangeShifter, please see:
Bocedi G., Palmer S.C.F., Pe’er G., Heikkinen R.K., Matsinos Y.G., Watts K.
and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
eco-evolutionary dynamics and species’ responses to environmental changes.
Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

Authors: Greta Bocedi & Steve Palmer, University of Aberdeen

Last updated: 1 July 2020 by Steve Palmer

------------------------------------------------------------------------------*/

#ifndef FormMainH
#define FormMainH

#include <System.Classes.hpp>
#include <Vcl.Controls.hpp>
#include <Vcl.StdCtrls.hpp>
#include <Vcl.Forms.hpp>
#include <Vcl.ExtCtrls.hpp>
#include <Buttons.hpp>
#include <Dialogs.hpp>
#include <ExtDlgs.hpp>
#include <ComCtrls.hpp>
#include <Menus.hpp>
#include <ActnList.hpp>
#include <VCLTee.Chart.hpp>
#include <VCLTee.Series.hpp>
#include <VCLTee.TeEngine.hpp>
#include <VCLTee.TeeProcs.hpp>
#include <VclTee.TeeGDIPlus.hpp>

#include <algorithm>
#include <fstream>
#include <io.h>
#include <stdio.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>
#include <utility>
#include <tchar.h>

//#include <sys/types.h>
//#include <sys/stat.h>

#include <float.h>
//#include <math.h>
//#include <Math.hpp>
#include <time.h>
#include <numeric>
using namespace std;

//#include "mathlib.h"

#include "Parameters.h"
#include "Model.h"
#include "Landscape.h"
#include "Species.h"
#include "Population.h"
#include "BatchMode.h"
#include "FormLand.h"
#include "FormArtificialLand.h"
#include "FormEnvGradient.h"
#include "FormSpecies.h"
#include "FormGenetics.h"
#include "FormSim.h"
#include "FormSeeding.h"
#include "FormVisualPatch.h"
#include "FormVisualCost.h"
#include "FormVisualTraits0.h"
#include "FormVisualTraits1.h"
#include "FormVisualTraits2.h"

//---------------------------------------------------------------------------
class TfrmMain : public TForm
{
__published:	// IDE-managed Components
	TMainMenu *MainMenu1;
	TMenuItem *File1;
	TMenuItem *Parameterset;
	TMenuItem *SetDirectory;
	TMenuItem *Simulations;
	TMenuItem *Pause;
	TMenuItem *Stop;
	TMenuItem *Run1;
	TMenuItem *Refresh;
	TMenuItem *LandMenu;
	TMenuItem *RasterLand;
	TMenuItem *Artificial;
	TMenuItem *SpeciesMenu;
	TMenuItem *GeneticsMenu;
	TMenuItem *BatchMode;
	TMenuItem *EnvGradient;
	TOpenDialog *OpenDialog1;
	TOpenDialog *OpenDialog2;
	TImage *LandImage;		// landscape map
	TImage *CommImage;		// community map
	TImage *SpDistImage;  // loaded species distribution map
	TImage *PopLegend;
	TPaintBox *MovtPaintBox;
	TMemo *Memo1;
	TChart *ChartOccSuit;
	TLineSeries *OccSuitmean;
	TLineSeries *OccSuitplusSE;
	TLineSeries *OccSuitminusSE;
	TChart *ChartPop;
	TLineSeries *Population;
	TLineSeries *Cells;
	TChart *ChartNoise;
	TLineSeries *EnvNoise;
	TLabel *LabelPop;
	TLabel *LabelMinPop;
	TLabel *LabelMaxPop;
	TScrollBox *LandScrollBox;
	TButton *BtnZoomIn;
	TButton *BtnZoomOut;
	void __fastcall SetDirectoryClick(TObject *Sender);
	void __fastcall BatchModeClick(TObject *Sender);
	void __fastcall LandMenuClick(TObject *Sender);
	void __fastcall RasterLandClick(TObject *Sender);
	void __fastcall ArtificialClick(TObject *Sender);
	void __fastcall EnvGradientClick(TObject *Sender);
	void __fastcall SpeciesMenuClick(TObject *Sender);
	void __fastcall GeneticsMenuClick(TObject *Sender);
	void __fastcall SimulationsClick(TObject *Sender);
	void __fastcall Run1Click(TObject *Sender);
	void __fastcall StopClick(TObject *Sender);
	void __fastcall PauseClick(TObject *Sender);
	void __fastcall RefreshClick(TObject *Sender);
	void __fastcall BtnZoomInClick(TObject *Sender);
	void __fastcall BtnZoomOutClick(TObject *Sender);

private:	// User declarations
public:		// User declarations
	__fastcall TfrmMain(TComponent* Owner);
	void __fastcall SetMapDimensions(void);
};
//---------------------------------------------------------------------------
extern PACKAGE TfrmMain *frmMain;
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
bool stopRun = false;
bool pause = false;
//bool drawZoomedOut = false; // POTENTIAL OPTION TO DISPLAY IN V1.0 FORMAT
// ***** NB WOULD REQUIRE Landscape attribute pix to be float, not int *****
// default dimensions of main landscape image and of gradient image
// ENSURE THESE MATCH THE SIZES OF THE IMAGES IN THE FORMS
const int def_land_width  = 600;
const int def_land_height = 1000;
const int def_grad_width  = 200;
const int def_grad_height = 3000;

bool directoryChosen = false;
bool readyToRun = false;
bool cancelled = false;
string map_name;

// species' distribution
int Sp_resol_ratio;

// artificial landscape
int artLandStatus;  // 0 = not created, 1 = fractal series created, 2 = random series created
										// 9 = generation parameters successfully entered

// time series parameters
bool chartNoise = false;
float *epsGlobal;

// stage-structure
bool forceInit = false; // condition to force initialisation if landscape or
												// stage structure has been changed
bool speciesOK; 				// to check whether species data have been validated
bool geneticsOK; 				// to check whether genetics data have been validated

//Simulation parameters-------------------------------------
bool batchMode = false;
string controlFile_name;
int nSimuls, nLandscapes;
bool landscapeLoaded = false; // indicates whether any landscape has been loaded
bool newLandscape = false;		// indicates whether a new landscape has been loaded
bool costsSet = false;				// indicates whether valid costs have been set
int simseq = 0;								// sequential simulation count on current landscape

ofstream RSlog; // performance log for recording simulation times, etc.

//---------------------------------------------------------------------------

// Graphics objects

Graphics::TBitmap *bmpLand; 		TCanvas *canvasLand;			// landscape map
Graphics::TBitmap *bmpComm; 		TCanvas *canvasComm;			// community map
Graphics::TBitmap *bmpSpDist; 	TCanvas *canvasSpDist;		// species' distribution
Graphics::TBitmap *bmpPatches;	TCanvas *canvasPatches; 	// for visualising patches
Graphics::TBitmap *bmpCosts;		TCanvas *canvasCosts; 		// for visualising costs
Graphics::TBitmap *bmpPopLgd;		TCanvas *canvasPopLgd; 		// population size legend
// for visualising gradient
Graphics::TBitmap *bmpGrad; TCanvas *canvasGrad;
Graphics::TBitmap *bmpGradL; TCanvas *canvasGradL;
// for visualising variable traits
Graphics::TBitmap *bmpT[18]; TCanvas *canvasT[18];
Graphics::TBitmap *bmpL[18]; TCanvas *canvasL[18];

//---------------------------------------------------------------------------

void SetGradTraitsForms(void);
void SetTraitsImage(
	const int,						// image no. on form
	const int,						// trait image no.
	const UnicodeString,	// caption
	const int,						// image height
	const int,						// image width
	const int,						// legend height
	const float,					// legend pixel scaling factor
	const float,					// trait mean
	const float,					// trait scaling factor
	const float,					// minimum value
	const float						// maximum value
);
void RedrawTraitImage(
	const int,					// trait image no
	const bool 					// environmental gradient
);
void SaveTraitMap(int,int,int,int,string,Graphics::TBitmap *map,string);

void DrawLandscape(void);
void DrawInitial(bool loaded);
void DrawPatches(void);
const rgb draw_wheel(int);

// these functions to have different version for GUI and batch applications ...
void MemoLine(UnicodeString);
#if RSDEBUG
void DebugGUI(string);
#endif

//---------------------------------------------------------------------------

#endif

