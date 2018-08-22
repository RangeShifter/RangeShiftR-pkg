/*------------------------------------------------------------------------------

RangeShifter v2.0 FormEnvGradientH

Input from the Environmental Gradient form

For full details of RangeShifter, please see:
Bocedi G., Palmer S.C.F., Pe’er G., Heikkinen R.K., Matsinos Y.G., Watts K.
and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
eco-evolutionary dynamics and species’ responses to environmental changes.
Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

Authors: Greta Bocedi & Steve Palmer, University of Aberdeen

Last updated: 12 February 2016 by Steve Palmer

------------------------------------------------------------------------------*/

#ifndef FormEnvGradientH
#define FormEnvGradientH
//---------------------------------------------------------------------------
#include <System.Classes.hpp>
#include <Vcl.Controls.hpp>
#include <Vcl.StdCtrls.hpp>
#include <Vcl.Forms.hpp>
#include <Vcl.ExtCtrls.hpp>
//#include <sstream> // SCFP 2/11/13 TEMPORARY
//#include <stdlib.h> // SCFP 2/11/13 TEMPORARY
//#include <string> // SCFP 2/11/13 TEMPORARY
//using namespace std; // SCFP 2/11/13 TEMPORARY

#include "Parameters.h"
#include "Landscape.h"

//---------------------------------------------------------------------------
class TfrmEnvGradient : public TForm
{
__published:	// IDE-managed Components
	TRadioGroup *RGgradientType;
	TLabeledEdit *edtGradInc;
	TLabeledEdit *edtOptRow;
	TLabeledEdit *edtDev;
	TLabeledEdit *edtOptEXT;
	TCheckBox *CBshift;
	TLabeledEdit *edtShiftRate;
	TLabeledEdit *edtStartShift;
	TLabeledEdit *edtStopShift;
	TBevel *Bevel3;
	TButton *BtnOK;
	TButton *BtnCancel;
	void __fastcall BtnOKClick(TObject *Sender);
	void __fastcall CBshiftClick(TObject *Sender);
	void __fastcall RGgradientTypeClick(TObject *Sender);
	void __fastcall BtnCancelClick(TObject *Sender);
	void __fastcall Refresh(void);
//	void __fastcall CBenvGradClick(TObject *Sender);
private:	// User declarations
public:		// User declarations
	__fastcall TfrmEnvGradient(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TfrmEnvGradient *frmEnvGradient;
//extern bool gradient;
extern bool patch_model;
extern int SimType;
extern bool fractal;
extern bool continuous;
extern paramGrad *paramsGrad;
extern Landscape *pLandscape;

//---------------------------------------------------------------------------
#endif
