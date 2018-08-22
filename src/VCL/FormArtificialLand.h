/*------------------------------------------------------------------------------

RangeShifter v2.0 FormArtificialLand

Input from the Artificial Landscapes form

For full details of RangeShifter, please see:
Bocedi G., Palmer S.C.F., Pe’er G., Heikkinen R.K., Matsinos Y.G., Watts K.
and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
eco-evolutionary dynamics and species’ responses to environmental changes.
Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

Authors: Greta Bocedi & Steve Palmer, University of Aberdeen

Last updated: 9 November 2016 by Steve Palmer

------------------------------------------------------------------------------*/

#ifndef FormArtificialLandH
#define FormArtificialLandH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <ExtCtrls.hpp>
#include <Grids.hpp>
#include <ComCtrls.hpp>

#include <vector>
#include <io.h>
#include <iomanip>
//#include <iostream>
#include <fstream>
#include <sstream>
//#include <stdlib.h>
//#include <string>
using namespace std;

//#include <math.h>
//#include "mathlib.h"
//#include <Math.hpp>
#include "Parameters.h"
#include "BatchMode.h"
#include "Landscape.h"
#include "FractalGenerator.h"

//---------------------------------------------------------------------------
class TfrmGenerateLand : public TForm
{
__published:	// IDE-managed Components
	TLabeledEdit *edtXdim;
	TLabeledEdit *edtYdim;
	TLabeledEdit *edtMaxProp;
	TLabeledEdit *edtMinProp;
	TButton *BtnOK;
	TLabeledEdit *edtResolution;
	TLabel *Label7;
	TLabeledEdit *edtSeriesNr;
	TLabeledEdit *edtNLand;
	TStringGrid *StringGridLand;
	TLabel *GridLandLabel;
	TButton *BtnCreateSeries;
	TRadioGroup *RGLandSeries;
	TStatusBar *StatusBar1;
	TRadioGroup *RGLandType;
	TButton *BtnCancel;
	TRadioGroup *RGFract;
	TMemo *Memo2;
	TCheckBox *CBpatch;
	TLabel *Label10;
	TLabeledEdit *edtMaxCells;
	TStringGrid *StringGridMatrix;
	TLabel *GridMatrixLabel;
	TCheckBox *CBmatrix;
	void __fastcall BtnOKClick(TObject *Sender);
	void __fastcall BtnCreateSeriesClick(TObject *Sender);
	void __fastcall RGLandSeriesClick(TObject *Sender);
	void __fastcall RGFractClick(TObject *Sender);
	void __fastcall RGLandTypeClick(TObject *Sender);
	void __fastcall BtnCancelClick(TObject *Sender);
	void __fastcall CBpatchClick(TObject *Sender);
	void __fastcall CBmatrixClick(TObject *Sender);

private:	// User declarations
	bool __fastcall CheckDimensions(void);
	int LandSeries(int,bool,bool);
	void CreatePatchLand(void);
	void DeletePatchLand(void);
	void AddNeighbours(int,int,int);
	void LandSeriesPatches(const string,float**);
	void __fastcall SetMemo(void);
public:		// User declarations
	__fastcall TfrmGenerateLand(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TfrmGenerateLand *frmGenerateLand;
extern Landscape *pLandscape;
extern paramSim *paramsSim;
extern string Inputs;
extern int artLandStatus;
	// 0 = not created, 1 = fractal series created, 2 = random series created

#if RSDEBUG
extern void DebugGUI(string);
#endif

//---------------------------------------------------------------------------
#endif
