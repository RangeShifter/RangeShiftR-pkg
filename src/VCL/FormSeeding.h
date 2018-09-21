/*------------------------------------------------------------------------------

RangeShifter v2.0 FormSeeding

Input from the Initialisation Rules form

For full details of RangeShifter, please see:
Bocedi G., Palmer S.C.F., Pe’er G., Heikkinen R.K., Matsinos Y.G., Watts K.
and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
eco-evolutionary dynamics and species’ responses to environmental changes.
Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

Authors: Greta Bocedi & Steve Palmer, University of Aberdeen

Last updated: 21 September 2018 by Steve Palmer

------------------------------------------------------------------------------*/

#ifndef FormSeedingH
#define FormSeedingH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <Buttons.hpp>
#include <ExtCtrls.hpp>
#include <Dialogs.hpp>
#include <ExtDlgs.hpp>
#include <Grids.hpp>
#include <ComCtrls.hpp>
#include <Vcl.Menus.hpp>

#include "Landscape.h"
#include "BatchMode.h"

//---------------------------------------------------------------------------
class TfrmSeeding : public TForm
{
__published:	// IDE-managed Components
	TBitBtn *BtnOK;
	TRadioGroup *RG1;
	TRadioGroup *RGRules;
	TRadioGroup *RGCells;
	TPanel *PanelSpDist;
	TLabel *LabelListBox;
	TListBox *LBSpPres;
	TBitBtn *BtnCancel;
	TLabeledEdit *edtNSpCells;
	TOpenTextFileDialog *OpenInitialisationFile;
	TBitBtn *BtnSaveInitial;
	TStatusBar *StatusBar1;
	TPanel *PanelNInds;
	TLabel *LabelNinds;
	TLabel *LabelPropn;
	TBevel *Bevel6;
	TComboBox *cboNinds;
	TEdit *edtNinds;
	TStringGrid *SGinitialStage;
	TLabel *LabelSpDistn;
	TPanel *PanelPlus;
	TLabel *LabelLandRes;
	TBevel *Bevel3;
	TCheckBox *CBAdditional;
	TLabeledEdit *edtAddX;
	TLabeledEdit *edtAddY;
	TButton *BtnAdd;
	TMemo *InitialList;
	TLabeledEdit *edtPatchID;
	TSaveTextFileDialog *SaveTextFileDialog1;
	TLabel *init_unit;
	TPanel *PanelFree;
	TLabel *LabInitPos;
	TLabel *LabelFreeInit;
	TBevel *Bevel5;
	TLabeledEdit *edtMaxY;
	TLabeledEdit *edtMaxX;
	TLabeledEdit *edtMinY;
	TLabeledEdit *edtMinX;
	TLabeledEdit *edtTotRandomCells;
	TRadioGroup *RGinitAges;
	TMainMenu *MainMenu_Seeding;
	TMenuItem *Refresh1;
	TLabeledEdit *edtInitFreezeYear;
	TLabeledEdit *edtRestrictRows;
	TLabeledEdit *edtRestrictFreq;
	TCheckBox *CBrestrictRange;
	TLabeledEdit *edtFinalFreezeYear;
	void __fastcall NindsChange(TObject *Sender);
	void __fastcall BtnOKClick(TObject *Sender);
	void __fastcall RG1Click(TObject *Sender);
	void __fastcall RGCellsClick(TObject *Sender);
	void __fastcall RGRulesClick(TObject *Sender);
	void __fastcall BtnCancelClick(TObject *Sender);
	void __fastcall BtnAddClick(TObject *Sender);
	void __fastcall CBAdditionalClick(TObject *Sender);
	void __fastcall BtnSaveInitialClick(TObject *Sender);
	void __fastcall Refresh1Click(TObject *Sender);
	void __fastcall CBrestrictRangeClick(TObject *Sender);

private:	// User declarations
public:		// User declarations
	__fastcall TfrmSeeding(TComponent* Owner);
	void __fastcall initialise(void);
	void __fastcall setDefaults(void);
	int  __fastcall CheckInitIndsFile(string);

};

//---------------------------------------------------------------------------
extern PACKAGE TfrmSeeding *frmSeeding;

extern Landscape *pLandscape;
extern Species *pSpecies;
extern bool landscapeLoaded;
extern bool forceInit;

//---------------------------------------------------------------------------
#endif
