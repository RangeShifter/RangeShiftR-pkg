/*------------------------------------------------------------------------------

RangeShifter v2.0 FormLand

Input from the Landscape form

For full details of RangeShifter, please see:
Bocedi G., Palmer S.C.F., Pe’er G., Heikkinen R.K., Matsinos Y.G., Watts K.
and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
eco-evolutionary dynamics and species’ responses to environmental changes.
Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

Authors: Greta Bocedi & Steve Palmer, University of Aberdeen

Last updated: 14 July 2020 by Steve Palmer

------------------------------------------------------------------------------*/

#ifndef FormLandH
#define FormLandH
#include <System.Classes.hpp>
#include <Vcl.Controls.hpp>
#include <Vcl.Dialogs.hpp>
#include <Vcl.ExtCtrls.hpp>
#include <Vcl.ExtDlgs.hpp>
#include <Vcl.Grids.hpp>
#include <Vcl.StdCtrls.hpp>
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <ExtCtrls.hpp>
#include <Grids.hpp>
#include <Buttons.hpp>
#include <Vcl.Dialogs.hpp>
#include <Vcl.ExtDlgs.hpp>

#include <vector>
#include <string>
using namespace std;

#include "FormDynLand.h"
#include "FormSeeding.h"
#include "Landscape.h"
#include "BatchMode.h"

//---------------------------------------------------------------------------
class TfrmLand : public TForm
{
__published:	// IDE-managed Components
	TEdit *edtRes;
	TButton *BtnImport;
	TLabel *LabelResolution;
	TLabeledEdit *edtNhabitat;
	TRadioGroup *RGhab;
	TButton *BtnCancel;
	TStringGrid *StringGridHab;
	TOpenTextFileDialog *OpenTextFileDialog1;
	TRadioGroup *RGCellPatch;
	TCheckBox *CBVisualPatch;
	TLabel *LabelNhab;
	TLabel *LabelSpDist;
	TBevel *Bevel1;
	TButton *BtnImportSp;
	TEdit *edtSpResol;
	TLabel *LabelSpResol;
	TButton *BtnChangeColours;
	TButton *BtnOK;
	TButton *BtnDynamic;
	void __fastcall FormCreate(TObject *Sender);
	void __fastcall BtnImportClick(TObject *Sender);
	void __fastcall BtnCancelClick(TObject *Sender);
	void __fastcall edtNhabitatExit(TObject *Sender);
	void __fastcall RGhabClick(TObject *Sender);
	void __fastcall RGCellPatchClick(TObject *Sender);
	void __fastcall BtnImportSpClick(TObject *Sender);
//	void __fastcall CBVisualPatchClick(TObject *Sender);
	void __fastcall BtnChangeColoursClick(TObject *Sender);
	void __fastcall BtnOKClick(TObject *Sender);
	void __fastcall BtnDynamicClick(TObject *Sender);

private:	// User declarations
public:		// User declarations
	__fastcall TfrmLand(TComponent* Owner);
	void __fastcall refresh(bool);
	void __fastcall displayColourTable(int);
};
//---------------------------------------------------------------------------
extern PACKAGE TfrmLand *frmLand;

extern Landscape *pLandscape;
extern bool landscapeLoaded;
extern bool newLandscape;
extern bool costsSet;

//---------------------------------------------------------------------------
#endif
