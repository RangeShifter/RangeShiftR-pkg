/*------------------------------------------------------------------------------

RangeShifter v2.0 FormDynLand

Input from the Dynamic Landscape form

For full details of RangeShifter, please see:
Bocedi G., Palmer S.C.F., Pe’er G., Heikkinen R.K., Matsinos Y.G., Watts K.
and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
eco-evolutionary dynamics and species’ responses to environmental changes.
Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

Authors: Steve Palmer, University of Aberdeen

Last updated: 14 November 2016 by Steve Palmer

------------------------------------------------------------------------------*/

#ifndef FormDynLandH
#define FormDynLandH
//---------------------------------------------------------------------------
#include <System.Classes.hpp>
#include <Vcl.Controls.hpp>
#include <Vcl.StdCtrls.hpp>
#include <Vcl.Forms.hpp>
#include <Vcl.ExtCtrls.hpp>
#include <Vcl.Dialogs.hpp>

#include "Version.h"
#include "Landscape.h"

//---------------------------------------------------------------------------
class TfrmDynLand : public TForm
{
__published:	// IDE-managed Components
	TButton *BtnFinished;
	TButton *BtnAdd;
	TLabeledEdit *edtYear;
	TLabeledEdit *edtChange;
	TOpenDialog *OpenDialog1;
	TButton *BtnReset;
	void __fastcall BtnAddClick(TObject *Sender);
	void __fastcall BtnFinishedClick(TObject *Sender);
	void __fastcall BtnResetClick(TObject *Sender);
private:	// User declarations
public:		// User declarations
	__fastcall TfrmDynLand(TComponent* Owner);
	void __fastcall refresh(void);
};
//---------------------------------------------------------------------------
extern PACKAGE TfrmDynLand *frmDynLand;

extern Landscape *pLandscape;
// error message defined in FormLand.cpp ...
extern string msgResol;
extern string msgHeaders;
extern string msgPatchHdr;
extern string msgLandHdr;
extern string msgReadError;

//---------------------------------------------------------------------------
#endif
