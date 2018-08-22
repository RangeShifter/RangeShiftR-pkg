/*------------------------------------------------------------------------------

RangeShifter v2.0 FormVisualPatch

Displays the Patch Landscape form

For full details of RangeShifter, please see:
Bocedi G., Palmer S.C.F., Pe’er G., Heikkinen R.K., Matsinos Y.G., Watts K.
and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
eco-evolutionary dynamics and species’ responses to environmental changes.
Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

Authors: Greta Bocedi & Steve Palmer, University of Aberdeen

Last updated: 12 February 2016 by Steve Palmer

------------------------------------------------------------------------------*/

#ifndef FormVisualPatchH
#define FormVisualPatchH

#include "FormLand.h"

//---------------------------------------------------------------------------
#include <System.Classes.hpp>
#include <Vcl.Controls.hpp>
#include <Vcl.StdCtrls.hpp>
#include <Vcl.Forms.hpp>
#include <Vcl.ExtCtrls.hpp>
//---------------------------------------------------------------------------
class TfrmVisualPatch : public TForm
{
__published:	// IDE-managed Components
	TScrollBox *PatchScrollBox;
	TImage *PatchImage;
	void __fastcall FormClose(TObject *Sender, TCloseAction &Action);
private:	// User declarations
public:		// User declarations
	__fastcall TfrmVisualPatch(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TfrmVisualPatch *frmVisualPatch;
extern PACKAGE TfrmLand *frmLand;

//---------------------------------------------------------------------------
#endif
