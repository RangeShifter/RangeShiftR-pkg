/*------------------------------------------------------------------------------

RangeShifter v2.0 FormMove

Input from the Movement Processes form

For full details of RangeShifter, please see:
Bocedi G., Palmer S.C.F., Pe’er G., Heikkinen R.K., Matsinos Y.G., Watts K.
and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
eco-evolutionary dynamics and species’ responses to environmental changes.
Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

Authors: Greta Bocedi & Steve Palmer, University of Aberdeen

Last updated: 22 August 2018 by Steve Palmer

------------------------------------------------------------------------------*/

#ifndef FormMoveH
#define FormMoveH

//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <ComCtrls.hpp>
#include <ExtCtrls.hpp>
#include <Grids.hpp>
#include <Vcl.Dialogs.hpp>

#include "Parameters.h"
#include "Landscape.h"
#include "Species.h"

//---------------------------------------------------------------------------
class TfrmMove : public TForm
{
__published:	// IDE-managed Components
	TButton *BtnOK;
	TButton *BtnCancel;
	TPanel *PanelRW;
	TPanel *PanelHab;
	TStringGrid *SGhabCosts;
	TCheckBox *CBCosts;
	TLabeledEdit *edtStepLength;
	TRadioGroup *RGMoveModel;
	TRadioGroup *RGstepM;
	TPanel *PanelSMS;
	TLabeledEdit *edtPRMethod;
	TLabeledEdit *edtPR;
	TLabeledEdit *edtDP;
	TLabeledEdit *edtRho;
	TCheckBox *CBIndVarCRW;
	TLabeledEdit *edtStepLMean;
	TLabeledEdit *edtStepLSD;
	TLabeledEdit *edtStepMort;
	TLabeledEdit *edtSM;
	TEdit *edtRhoMean;
	TEdit *edtRhoSD;
	TLabel *PRLabel;
	TLabel *PRmethodLabel;
	TLabel *DPLabel;
	TMemo *Memo2;
	TOpenDialog *OpenDialog1;
	TCheckBox *CBVisualCosts;
	TLabeledEdit *edtStepLScale;
	TEdit *edtRhoScale;
	TLabeledEdit *edtMemorySize;
	TLabeledEdit *edtGoalBias;
	TLabeledEdit *edtGoalType;
	TLabel *GoalBiasLabel;
	TLabel *MemorySizeLabel;
	TLabel *GoalTypeLabel;
	TCheckBox *CBstraightPath;
	TLabeledEdit *edtAlphaDB;
	TLabeledEdit *edtBetaDB;
	TLabel *AlphaDBlabel;
	TLabel *BetaDBlabel;
	TPanel *PanelDB;
	void __fastcall BtnOKClick(TObject *Sender);
	void __fastcall RGstepMClick(TObject *Sender);
	void __fastcall CBCostsClick(TObject *Sender);
	void __fastcall RGMoveModelClick(TObject *Sender);
	void __fastcall CBIndVarCRWClick(TObject *Sender);
	void __fastcall BtnCancelClick(TObject *Sender);
	void __fastcall edtGoalTypeChange(TObject *Sender);
//	void __fastcall CBEvolRWClick(TObject *Sender);
private:	// User declarations
public:		// User declarations
	__fastcall TfrmMove(TComponent* Owner);
	void __fastcall refresh(void);
	void __fastcall setHabTable(void);
};
//---------------------------------------------------------------------------

extern Graphics::TBitmap *bmpCosts;
extern TCanvas *canvasCosts;

extern PACKAGE TfrmMove *frmMove;

extern Landscape *pLandscape;	// see FormMain.cpp
extern Species   *pSpecies;		// see FormMain.cpp
extern paramSim  *paramsSim;	// see FormMain.cpp

//extern int import_CostsLand(string);
extern bool movtModelOK;	// see FormSpecies.cpp
extern bool newLandscape;	// see FormMain.h
extern bool costsSet;			// see FormMain.h

//---------------------------------------------------------------------------
#endif
