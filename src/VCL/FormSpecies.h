/*------------------------------------------------------------------------------

RangeShifter v2.0 FormSpecies

Input from the Species Parameters form

For full details of RangeShifter, please see:
Bocedi G., Palmer S.C.F., Pe’er G., Heikkinen R.K., Matsinos Y.G., Watts K.
and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
eco-evolutionary dynamics and species’ responses to environmental changes.
Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

Authors: Greta Bocedi & Steve Palmer, University of Aberdeen

Last updated: 27 June 2017 by Steve Palmer

------------------------------------------------------------------------------*/

#ifndef FormSpeciesH
#define FormSpeciesH

//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <ComCtrls.hpp>
#include <ExtCtrls.hpp>
#include <Buttons.hpp>
#include <Grids.hpp>
#include <Dialogs.hpp>
#include <ExtDlgs.hpp>

#include <VCLTee.Chart.hpp>
#include <VCLTee.Series.hpp>
#include <VCLTee.TeEngine.hpp>
#include <VCLTee.TeeProcs.hpp>

#include "FormDensity.h"
#include "FormMove.h"
#include "FormSeeding.h"
#include "FormSim.h"
#include "Parameters.h"
#include "Species.h"
#include <VclTee.TeeGDIPlus.hpp>

//---------------------------------------------------------------------------
class TfrmSpecies : public TForm
{
__published:	// IDE-managed Components
	TPageControl *PageControl1;
	TTabSheet *TSPop;
	TTabSheet *TSDispersal;
	TButton *BtnOK;
	TPanel *PanelK;
	TStringGrid *SGhab;
	TButton *BtnCancel;
	TTabSheet *TSDispersal2;
	TRadioGroup *RGReproduction;
	TLabeledEdit *edtR;
	TCheckBox *CBStageModel;
	TPanel *PanelStage;
	TLabeledEdit *edtNstages;
	TStringGrid *transMatrix;
	TLabel *transMatrixLabel;
	TLabeledEdit *edtMaxAge;
	TBevel *Bevel1;
	TButton *BtnSetCoeff;
	TStringGrid *MinAges;
	TLabel *MinAgesLabel;
	TLabeledEdit *edtSexRatio;
	TGroupBox *GroupBox1;
	TCheckBox *CBFecundity;
	TCheckBox *CBDevelopment;
	TCheckBox *CBSurvival;
	TLabeledEdit *edtHarem;
	TPanel *PanelEmigration;
	TRadioGroup *RGEmigProb;
	TPanel *PanelTransfer;
	TRadioGroup *RGMovements;
	TBitBtn *BtnMovements;
	TRadioGroup *RGKernel;
	TPanel *PanelKernels;
	TButton *BtnUpdKernel;
	TCheckBox *CBIndVarKernel;
	TPanel *PanelMortality;
	TLabeledEdit *edtMortProb;
	TLabeledEdit *edtMortSlope;
	TLabeledEdit *edtMortInfl;
	TPanel *PanelSettProcess;
	TLabeledEdit *edtMinSteps;
	TRadioGroup *RGSteps;
	TLabeledEdit *edtSettS0;
	TLabeledEdit *edtSettAlpha;
	TLabeledEdit *edtSettBeta;
	TLabeledEdit *edtNsteps;
	TLabeledEdit *edtMaxStepYear;
	TMemo *Memo2;
	TCheckBox *CBSexKernels;
	TCheckBox *CBSexSettMovt;
	TCheckBox *CBStageKernels;
	TCheckBox *CBStageSettMovt;
	TPanel *PanelSexEmig;
	TPanel *PanelSexTransfer;
	TPanel *PanelSexSettle;
	TStringGrid *SexEmigPar;
	TStringGrid *SexKernPar;
	TStringGrid *SexSettlePar;
	TStringGrid *SexEmigDensPar;
	TLabel *LabelSexEmigPar;
	TLabel *LabelSexEmigDensPar;
	TPanel *PanelSettKernels;
	TRadioGroup *RGSettKern;
	TLabel *EmigrationLabel;
	TLabel *TransferLabel;
	TLabel *SettProcessLabel;
	TLabel *SettKernelsLabel;
	TLabel *SexEmigLabel;
	TLabel *SexTransferLabel;
	TLabel *SexSettleParLabel;
	TPanel *PanelEP;
	TLabel *LabelDensIndept;
	TLabeledEdit *edtEP;
	TLabeledEdit *edtEPmean;
	TLabeledEdit *edtEPsd;
	TPanel *PanelDensEmig;
	TLabeledEdit *edtD0Mean;
	TLabeledEdit *edtD0SD;
	TLabeledEdit *edtD0;
	TLabeledEdit *edtAlpha;
	TLabeledEdit *edtBeta;
	TEdit *edtAlphaMean;
	TEdit *edtBetaMean;
	TEdit *edtAlphaSD;
	TEdit *edtBetaSD;
	TPanel *PanelDist;
	TLabeledEdit *edtDist1;
	TLabeledEdit *edtDist2;
	TLabeledEdit *edtPKern1;
	TEdit *edtDist2Mean;
	TLabeledEdit *edtDist1Mean;
	TLabeledEdit *edtDist1SD;
	TEdit *edtDist2SD;
	TEdit *edtPKern1SD;
	TEdit *edtPKern1Mean;
	TCheckBox *CBSexSettKern;
	TCheckBox *CBStageSettKern;
	TCheckBox *CBSettKernMate;
	TLabeledEdit *edtC;
	TLabeledEdit *edtDist1Scale;
	TEdit *edtDist2Scale;
	TEdit *edtPKern1Scale;
	TRadioGroup *RGMortality;
	TLabeledEdit *edtEPscale;
	TLabeledEdit *edtD0Scale;
	TEdit *edtAlphaScale;
	TEdit *edtBetaScale;
	TLabel *SGhabLabel;
	TLabeledEdit *edtK;
	TLabel *edtKLabel;
	TLabeledEdit *edtPRep;
	TLabeledEdit *edtRepInterval;
	TGroupBox *GroupBox2;
	TCheckBox *CBweightFec;
	TCheckBox *CBweightDev;
	TCheckBox *CBweightSurv;
	TLabeledEdit *edtDevCoeff;
	TLabeledEdit *edtSurvCoeff;
	TBevel *Bevel2;
	TRadioGroup *RGSurvival;
	TLabel *SurvivalLabel;
	TLabeledEdit *edtGen;
	TPanel *VarTextPanel;
	TLabel *VarLabel1;
	TLabel *VarLabel2;
	TLabel *VarLabel3;
	TPanel *EmigCheckboxPanel;
	TCheckBox *CBFullKernel;
	TCheckBox *CBSexEmig;
	TCheckBox *CBStageEmig;
	TCheckBox *CBIndVarEmig;
	TChart *KernelGraph;
	TAreaSeries *Kernel1;
	TAreaSeries *Kernel2;
	TCheckBox *CBIndVarSettle;
	TPanel *PanelSettleDD;
	TLabeledEdit *edtSettS0Mean;
	TLabeledEdit *edtSettS0SD;
	TEdit *edtSettAlphaMean;
	TEdit *edtSettBetaMean;
	TEdit *edtSettAlphaSD;
	TEdit *edtSettBetaSD;
	TLabeledEdit *edtSettS0Scale;
	TEdit *edtSettAlphaScale;
	TEdit *edtSettBetaScale;
	TPanel *PanelMateDD;
	TCheckBox *CBDensDepSettle;
	TCheckBox *CBFindMate;
	TLabeledEdit *edtEmigStage;
	TRadioGroup *RGaction;
	void __fastcall BtnOKClick(TObject *SendeSexTranferParr);
	void __fastcall RGMovementsClick(TObject *Sender);
	void __fastcall BtnMovementsClick(TObject *Sender);
	void __fastcall RGEmigProbClick(TObject *Sender);
	void __fastcall RGKernelClick(TObject *Sender);
	void __fastcall BtnUpdKernelClick(TObject *Sender);
	void __fastcall RGMortalityClick(TObject *Sender);
	void __fastcall RGStepsClick(TObject *Sender);
	void __fastcall CBStageModelClick(TObject *Sender);
	void __fastcall RGReproductionClick(TObject *Sender);
	void __fastcall edtNstagesExit(TObject *Sender);
	void __fastcall MinAgesExit(TObject *Sender);
	void __fastcall BtnSetCoeffClick(TObject *Sender);
	void __fastcall CBFecundityClick(TObject *Sender);
	void __fastcall CBDevelopmentClick(TObject *Sender);
	void __fastcall CBSurvivalClick(TObject *Sender);
	void __fastcall CBSexEmigClick(TObject *Sender);
	void __fastcall CBStageEmigClick(TObject *Sender);
	void __fastcall CBIndVarEmigClick(TObject *Sender);
	void __fastcall CBSexKernelsClick(TObject *Sender);
	void __fastcall CBStageKernelsClick(TObject *Sender);
	void __fastcall CBIndVarKernelClick(TObject *Sender);
	void __fastcall CBSexSettMovtClick(TObject *Sender);
	void __fastcall CBStageSettMovtClick(TObject *Sender);
	void __fastcall TSDispersalExit(TObject *Sender);
	void __fastcall TSDispersal2Enter(TObject *Sender);
	void __fastcall CBSexSettKernClick(TObject *Sender);
	void __fastcall CBStageSettKernClick(TObject *Sender);
	void __fastcall RGSettKernClick(TObject *Sender);
	void __fastcall CBweightFecClick(TObject *Sender);
	void __fastcall CBweightDevClick(TObject *Sender);
	void __fastcall CBweightSurvClick(TObject *Sender);
	void __fastcall BtnCancelClick(TObject *Sender);
	void __fastcall emigCheckboxPanelExit(TObject *Sender);
	void __fastcall CBIndVarSettleClick(TObject *Sender);
	void __fastcall CBDensDepSettleClick(TObject *Sender);

private:	// User declarations
public:		// User declarations
	__fastcall TfrmSpecies(TComponent* Owner);
	void __fastcall refresh(bool);
	void __fastcall setEmigFields(void);
	void __fastcall setKernelFields(void);
	void __fastcall setSettleFields(void);
//	void __fastcall setSettleDD(void);
	void __fastcall resetTransMatrix(int);
	void __fastcall resetSexEmigPar(int);
	void __fastcall resetSexEmigDensPar(int);
	void __fastcall resetSexTransferPar(int);
	void __fastcall resetSexSettlePar(int);
	void __fastcall updateSpecies(void);
};
//---------------------------------------------------------------------------
extern PACKAGE TfrmSpecies *frmSpecies;
//---------------------------------------------------------------------------

extern bool forceInit;
extern bool speciesOK;
extern bool newLandscape;
//extern Species *pSpecies;

#endif
