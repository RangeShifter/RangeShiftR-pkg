/*------------------------------------------------------------------------------

RangeShifter v2.0 FormSim

Input from the Simulation Parameters form

For full details of RangeShifter, please see:
Bocedi G., Palmer S.C.F., Pe’er G., Heikkinen R.K., Matsinos Y.G., Watts K.
and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
eco-evolutionary dynamics and species’ responses to environmental changes.
Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

Authors: Greta Bocedi & Steve Palmer, University of Aberdeen

Last updated: 3 March 2017 by Steve Palmer

------------------------------------------------------------------------------*/

#ifndef FormSimH
#define FormSimH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <ExtCtrls.hpp>
#include <Buttons.hpp>
#include <Grids.hpp>
#include <Vcl.ActnColorMaps.hpp>
#include <Vcl.ActnMan.hpp>

#include "Parameters.h"
#include "Landscape.h"
#include "FormSeeding.h"

//---------------------------------------------------------------------------
class TfrmSim : public TForm
{
__published:	// IDE-managed Components
	TButton *BtnOK;
	TPanel *PanelOutputs;
	TRadioGroup *RGMap;
	TPanel *PanelMaps;
	TLabeledEdit *edtMap_interval;
	TCheckBox *CBDrawInit;
	TPanel *PanelEnvStoch;
	TLabeledEdit *edtAC;
	TLabeledEdit *edtStd;
	TLabel *LabelAC;
	TLabel *LabelStd;
	TCheckBox *CBEnvStoch;
	TLabeledEdit *edtMaxR;
	TLabeledEdit *edtMinR;
	TButton *BtnCancel;
	TGroupBox *GroupBoxOutputs;
	TCheckBox *CBoutPop;
	TCheckBox *CBoutTraitCell;
	TCheckBox *CBoutOcc;
	TCheckBox *CBoutTraitRow;
	TCheckBox *CBoutInd;
	TGroupBox *GroupBoxVisual;
	TCheckBox *CBVisualLand;
	TCheckBox *CBVisualTraits;
	TCheckBox *CBVisualGraph;
	TCheckBox *CBVisualPop;
	TCheckBox *CBVisualGrad;
	TCheckBox *CBVisualMove;
	TBevel *Bevel3;
	TCheckBox *CBLocExt;
	TEdit *edtLocExt;
	TCheckBox *CBoutRange;
	TRadioGroup *RGEnvStoch;
	TRadioGroup *RGEnvStochType;
	TBevel *Bevel2;
	TLabeledEdit *edtYears;
	TLabeledEdit *edtRep;
	TLabeledEdit *edtSimNr;
	TLabel *SimulationDesc;
	TBitBtn *BtnSeeding;
	TCheckBox *CBoutConnect;
	TRadioGroup *RGTraitsMap;
	TLabeledEdit *edtTraitsMap_Int;
	TLabeledEdit *edtMovtSlow;
	TEdit *edtFreqRange;
	TEdit *edtFreqOcc;
	TEdit *edtFreqPop;
	TEdit *edtFreqInd;
	TEdit *edtFreqTraitCell;
	TEdit *edtFreqTraitRow;
	TEdit *edtFreqConn;
	TLabel *LabelFreq;
	TEdit *edtStartRange;
	TEdit *edtStartOcc;
	TEdit *edtStartPop;
	TEdit *edtStartInd;
	TEdit *edtStartTraitCell;
	TEdit *edtStartTraitRow;
	TEdit *edtStartConn;
	TLabel *LabelStart;
	TEdit *edtStartGen;
	TEdit *edtFreqGen;
	TCheckBox *CBoutGen;
	TRadioGroup *RGGen;
	TCheckBox *CBGenCrosstab;
	TCheckBox *CBabsorbing;
	TCheckBox *CBSMSheatmap;
	void __fastcall BtnOKClick(TObject *Sender);
	void __fastcall BtnSeedingClick(TObject *Sender);
	void __fastcall RGMapClick(TObject *Sender);
	void __fastcall CBEnvStochClick(TObject *Sender);
	void __fastcall BtnCancelClick(TObject *Sender);
	void __fastcall CBLocExtClick(TObject *Sender);
	void __fastcall CBoutOccClick(TObject *Sender);
	void __fastcall RGEnvStochTypeClick(TObject *Sender);
	void __fastcall RGEnvStochClick(TObject *Sender);
	void __fastcall CBVisualLandClick(TObject *Sender);
	void __fastcall CBoutConnectClick(TObject *Sender);
	void __fastcall CBVisualTraitsClick(TObject *Sender);
	void __fastcall CBVisualMoveClick(TObject *Sender);
	void __fastcall edtRepChange(TObject *Sender);
	void __fastcall CBoutTraitCellClick(TObject *Sender);
	void __fastcall CBoutTraitRowClick(TObject *Sender);
	void __fastcall CBoutGenClick(TObject *Sender);
	void __fastcall RGTraitsMapClick(TObject *Sender);
	void __fastcall CBoutPopClick(TObject *Sender);
	void __fastcall CBoutIndClick(TObject *Sender);

private:	// User declarations
public:		// User declarations
	__fastcall TfrmSim(TComponent* Owner);
	void __fastcall refresh(bool,bool);
};
//---------------------------------------------------------------------------
extern PACKAGE TfrmSim *frmSim;

extern bool chartNoise;
extern bool patch_model;
extern bool readyToRun;
extern bool cancelled;
extern paramGrad *paramsGrad;
extern Landscape *pLandscape;
extern Species *pSpecies;

//---------------------------------------------------------------------------
#endif
