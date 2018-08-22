/*------------------------------------------------------------------------------

RangeShifter v2.0 FormDensity

Input from the Stage's Weights form

For full details of RangeShifter, please see:
Bocedi G., Palmer S.C.F., Pe’er G., Heikkinen R.K., Matsinos Y.G., Watts K.
and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
eco-evolutionary dynamics and species’ responses to environmental changes.
Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

Authors: Greta Bocedi & Steve Palmer, University of Aberdeen

Last updated: 12 February 2016 by Steve Palmer

------------------------------------------------------------------------------*/

#ifndef FormDensityH
#define FormDensityH
//---------------------------------------------------------------------------
#include <System.Classes.hpp>
#include <Vcl.Controls.hpp>
#include <Vcl.StdCtrls.hpp>
#include <Vcl.Forms.hpp>
#include <Vcl.Grids.hpp>
//---------------------------------------------------------------------------
class TfrmDensity : public TForm
{
__published:	// IDE-managed Components
	TStringGrid *FecundityDens;
	TLabel *Label4;
	TStringGrid *DevelopDens;
	TLabel *Label1;
	TStringGrid *SurvivalDens;
	TLabel *Label2;
	TButton *BtnOK;
	TButton *BtnCancel;
	void __fastcall BtnOKClick(TObject *Sender);
	void __fastcall BtnCancelClick(TObject *Sender);
private:	// User declarations
public:		// User declarations
	__fastcall TfrmDensity(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TfrmDensity *frmDensity;
//---------------------------------------------------------------------------

// declarations of functions defined in FormSpecies.cpp
void resetFecundityDens(TfrmDensity*,int);
void resetDevelopDens(TfrmDensity*,int);
void resetSurvivalDens(TfrmDensity*,int);

#endif
