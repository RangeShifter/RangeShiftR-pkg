//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop
#include <stdlib.h>
#include <string>
#include <string.h>
#include <sstream>
using namespace std;

#include "FormSim.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TfrmSim *frmSim;

//---------------------------------------------------------------------------
__fastcall TfrmSim::TfrmSim(TComponent* Owner)
	: TForm(Owner)
{
Left = 100;
Top = 60;
}

//---------------------------------------------------------------------------
void __fastcall TfrmSim::refresh(bool patchModel,bool generated)
{
envGradParams grad = paramsGrad->getGradient();
demogrParams dem = pSpecies->getDemogr();
emigRules emig = pSpecies->getEmig();
trfrRules trfr = pSpecies->getTrfr();
settleType sett = pSpecies->getSettle();
genomeData gen = pSpecies->getGenomeData();

if (patchModel) {
	RGEnvStoch->ItemIndex = 0;
	CBoutTraitCell->Caption = "Mean Traits by patch";
	CBLocExt->Enabled = false; CBLocExt->Visible = false;
	CBLocExt->Checked = false;
	edtLocExt->Enabled = false; edtLocExt->Visible = false;
	CBoutTraitRow->Checked = false; CBoutTraitRow->Enabled = false;
	edtStartTraitRow->Enabled = false;
	edtFreqTraitRow->Enabled = false;
	RGEnvStoch->Enabled = false;
	CBoutConnect->Enabled = true;   CBoutConnect->Visible = true;
}
else { // cell-based model
	CBoutTraitCell->Caption = "Mean Traits by cell";
	if (grad.gradient && grad.gradType == 3) { // gradient in local extinction probability
		// additional local extinction is not permitted
		CBLocExt->Checked = false;
		CBLocExt->Enabled = false; CBLocExt->Visible = true;
		edtLocExt->Enabled = false; edtLocExt->Visible = false;
	}
	else {
		CBLocExt->Enabled = true; CBLocExt->Visible = true;
		if (CBLocExt->Checked) {
			edtLocExt->Enabled = true; edtLocExt->Visible = true;
		}
		else {
			edtLocExt->Enabled = false; edtLocExt->Visible = false;
		}
	}
	if (emig.indVar || trfr.indVar || sett.indVar) {
		CBoutTraitRow->Enabled = true;
		if (CBoutTraitRow->Checked) {
			edtStartTraitRow->Enabled = true;
			edtFreqTraitRow->Enabled = true;
		}
		else {
			edtStartTraitRow->Enabled = false;
			edtFreqTraitRow->Enabled = false;
		}
	}
	else {
		CBoutTraitRow->Checked = false; CBoutTraitRow->Enabled = false;
		edtStartTraitRow->Enabled = false;
		edtFreqTraitRow->Enabled = false;
	}
	CBoutConnect->Enabled = false; CBoutConnect->Visible = false;
}

if (generated) {
	CBoutOcc->Enabled = false;
	edtStartOcc->Enabled = false;
	edtFreqOcc->Enabled = false;
}
else {
	if (StrToInt(edtRep->Text) == 1) { // single replicate
		CBoutOcc->Checked = false; CBoutOcc->Enabled = false;
		edtStartOcc->Enabled = false;
		edtFreqOcc->Enabled = false;
	}
	else {
		CBoutOcc->Enabled = true;
		if (CBoutOcc->Checked) {
			edtStartOcc->Enabled = true;
			edtFreqOcc->Enabled = true;
		}
	}
}

if (grad.gradient)
	CBVisualGrad->Enabled = true;
else {
	CBVisualGrad->Checked = false; CBVisualGrad->Enabled = false;
}

if (dem.stageStruct) {
	RGEnvStochType->Items->Strings[0] = "in fecundities";
	RGEnvStochType->Items->Strings[1] = "in dens. dependence";
	edtMinR->EditLabel->Caption = "Min. fecundity ";
	edtMaxR->EditLabel->Caption = "Max. fecundity ";
}
else {
	RGEnvStochType->Items->Strings[0] = "in growth rate";
	RGEnvStochType->Items->Strings[1] = "in carrying capacity";
	edtMinR->EditLabel->Caption = "Min. growth rate ";
	edtMaxR->EditLabel->Caption = "Max. growth rate ";
}

if (CBoutPop->Checked) {
	edtStartPop->Enabled = true; edtFreqPop->Enabled = true;
}
else {
	edtStartPop->Enabled = false; edtFreqPop->Enabled = false;
}

if (CBoutInd->Checked) {
	edtStartInd->Enabled = true; edtFreqInd->Enabled = true;
}
else {
	edtStartInd->Enabled = false; edtFreqInd->Enabled = false;
}

if (emig.indVar || trfr.indVar || sett.indVar || gen.neutralMarkers) {
	CBoutGen->Enabled = true;
	if (CBoutGen->Checked) {
		edtStartGen->Enabled = true; edtFreqGen->Enabled = true;
		CBGenCrosstab->Enabled = true; CBGenCrosstab->Visible = true;
	}
	else {
		edtStartGen->Enabled = false; edtFreqGen->Enabled = false;
		CBGenCrosstab->Enabled = false; CBGenCrosstab->Visible = false;
	}
}
else {
	CBoutGen->Checked = false; CBoutGen->Enabled = false;
	edtStartGen->Enabled = false; edtFreqGen->Enabled = false;
}

if (emig.indVar || trfr.indVar || sett.indVar) {
	CBoutTraitCell->Enabled = true;
	if (CBoutTraitCell->Checked) {
		edtStartTraitCell->Enabled = true; edtFreqTraitCell->Enabled = true;
	}
	else {
		edtStartTraitCell->Enabled = false; edtFreqTraitCell->Enabled = false;
	}
	CBVisualTraits->Enabled = true;
}
else {
	CBoutTraitCell->Checked = false; CBoutTraitCell->Enabled = false;
	edtStartTraitCell->Enabled = false; edtFreqTraitCell->Enabled = false;
	CBVisualTraits->Checked = false; CBVisualTraits->Enabled = false;
}

if (trfr.moveModel) {
	CBVisualMove->Enabled = true;
	if (CBVisualMove->Checked) edtMovtSlow->Enabled = true;
	else edtMovtSlow->Enabled = false;
}
else {
	CBVisualMove->Enabled = false;
	edtMovtSlow->Enabled = false;
}

}

//---------------------------------------------------------------------------
void __fastcall TfrmSim::edtRepChange(TObject *Sender)
{
landParams ppLand = pLandscape->getLandParams();
if (ppLand.generated || StrToInt(edtRep->Text) < 2) {
	CBoutOcc->Checked = false; CBoutOcc->Enabled = false;
}
else
	CBoutOcc->Enabled = true;
}

//---------------------------------------------------------------------------
void __fastcall TfrmSim::CBEnvStochClick(TObject *Sender)
{
landParams ppLand = pLandscape->getLandParams();

if (CBEnvStoch->Checked) {
	edtAC->Enabled = true;		LabelAC->Enabled = true;
	edtStd->Enabled = true;		LabelStd->Enabled = true;
	edtMinR->Enabled = true;	edtMaxR->Enabled = true;
	if (ppLand.patchModel) {
		RGEnvStoch->ItemIndex = 0; RGEnvStoch->Enabled = false;
	}
	else RGEnvStoch->Enabled = true;
	RGEnvStochType->Enabled = true;
}
else {
	edtAC->Enabled = false;		LabelAC->Enabled = false;
	edtStd->Enabled = false;	LabelStd->Enabled = false;
	edtMinR->Enabled = false;	edtMaxR->Enabled = false;
	RGEnvStoch->Enabled = false;
	RGEnvStochType->Enabled = false;
}
}

//---------------------------------------------------------------------------
void __fastcall TfrmSim::BtnCancelClick(TObject *Sender)
{
cancelled = true;
Close();
}

//---------------------------------------------------------------------------
void __fastcall TfrmSim::CBLocExtClick(TObject *Sender)
{
envGradParams g;

if (CBLocExt->Checked) {
	g = paramsGrad->getGradient();
	if (g.gradient && g.gradType == 3) { // gradient in local extinction probability
		// NOTE - this combination should not occur, as local extinction check-box is unchecked
		// and disabled for gradient in local extinction probability
		MessageDlg("Gradient in local extinction probability exists.\nAdditional extinction is not permitted.",
				mtError, TMsgDlgButtons() << mbRetry,0);
		CBLocExt->Checked = false;
	}
	else {
		edtLocExt->Enabled = true; edtLocExt->Visible = true;
	}
}
else { edtLocExt->Enabled = false; edtLocExt->Visible = false; }
}
//---------------------------------------------------------------------------
void __fastcall TfrmSim::CBoutOccClick(TObject *Sender)
{
if (CBoutOcc->Checked && StrToInt(edtRep->Text) == 1) {
	MessageDlg("This output can be computed only in case of multiple replicates",
		mtWarning, TMsgDlgButtons() << mbRetry,0);
	CBoutOcc->Checked = false;
}
if (CBoutOcc->Checked) {
//	edtStartOcc->Enabled = true;
	edtFreqOcc->Enabled = true;
}
else {
//	edtStartOcc->Enabled = false;
	edtFreqOcc->Enabled = false;
}
}
//---------------------------------------------------------------------------
void __fastcall TfrmSim::CBoutConnectClick(TObject *Sender)
{
if (CBoutConnect->Checked) {
	edtStartConn->Enabled = true; edtStartConn->Visible = true;
	edtFreqConn->Enabled  = true; edtFreqConn->Visible  = true;
}
else {
	edtStartConn->Enabled = false; edtStartConn->Visible = false;
	edtFreqConn->Enabled  = false; edtFreqConn->Visible  = false;
}
}

//---------------------------------------------------------------------------
void __fastcall TfrmSim::BtnSeedingClick(TObject *Sender)
{
frmSeeding->initialise();
frmSeeding->Show();
}

//---------------------------------------------------------------------------
void __fastcall TfrmSim::RGEnvStochClick(TObject *Sender)
{
landParams ppLand = pLandscape->getLandParams();

if (ppLand.patchModel) {
	if (RGEnvStoch->ItemIndex == 1) {
		MessageDlg("This option is not active for a patch-based model",
			mtError, TMsgDlgButtons() << mbOK,0);
		RGEnvStoch->ItemIndex = 0;
	}
}
}

//---------------------------------------------------------------------------
void __fastcall TfrmSim::RGEnvStochTypeClick(TObject *Sender)
{
demogrParams dem = pSpecies->getDemogr();

if (RGEnvStochType->ItemIndex == 0) { // in growth rate / fecundity
	if (dem.stageStruct) {
		edtMaxR->EditLabel->Caption = "Max. fecundity ";
		edtMinR->EditLabel->Caption = "Min. fecundity ";
	}
	else{
		edtMaxR->EditLabel->Caption = "Max. growth rate ";
		edtMinR->EditLabel->Caption = "Min. growth rate ";
	}
}
else { // in carrying capacity / density dependence
	if (dem.stageStruct) {
		edtMaxR->EditLabel->Caption = "Maximum 1/b ";
		edtMinR->EditLabel->Caption = "Minimum 1/b ";
	}
	else {
		edtMaxR->EditLabel->Caption = "Maximum K ";
		edtMinR->EditLabel->Caption = "Minimum K ";
	}
}
}

//---------------------------------------------------------------------------
void __fastcall TfrmSim::CBVisualLandClick(TObject *Sender)
{
trfrRules trfr = pSpecies->getTrfr();
simView v = paramsSim->getViews();

if (CBVisualLand->Checked) {
	CBVisualPop->Enabled = true;
	if (trfr.moveModel) CBVisualMove->Enabled = true;
	else CBVisualMove->Enabled = false;
	RGMap->Enabled = true;   RGMap->ItemIndex = 1;
}
else {
	CBVisualPop->Enabled = false; CBVisualPop->Checked = false;
	if (trfr.moveModel) {
		if (v.viewCosts)
			CBVisualMove->Enabled = true;
		else CBVisualMove->Enabled = false;
	}
	else CBVisualMove->Enabled = false;
	if (CBVisualMove->Enabled) edtMovtSlow->Enabled = true;
	else edtMovtSlow->Enabled = false;

	RGMap->Enabled = false;  RGMap->ItemIndex = 0;
}
if (CBVisualMove->Enabled && CBVisualMove->Checked) edtMovtSlow->Enabled = true;
else edtMovtSlow->Enabled = false;
}

//---------------------------------------------------------------------------
void __fastcall TfrmSim::RGMapClick(TObject *Sender)
{
landParams ppLand = pLandscape->getLandParams();
//if (ppLand.spDist) {
//	MessageDlg("ppLand.spDist is TRUE",mtWarning, TMsgDlgButtons() << mbOK,0);
//}
//else {
//	MessageDlg("ppLand.spDist is FALSE",mtWarning, TMsgDlgButtons() << mbOK,0);
//}
if (RGMap->ItemIndex == 1) {
	PanelMaps->Visible = true; PanelMaps->Enabled = true;
	if (ppLand.spDist)
		CBDrawInit->Enabled = true;
	else CBDrawInit->Enabled = false;
}
else {
	PanelMaps->Visible = false; PanelMaps->Enabled = false;
	CBDrawInit->Enabled = false;
}
}

//---------------------------------------------------------------------------

void __fastcall TfrmSim::RGTraitsMapClick(TObject *Sender)
{
if (RGTraitsMap->ItemIndex == 0) edtTraitsMap_Int->Enabled = false;
else edtTraitsMap_Int->Enabled = true;
}

//---------------------------------------------------------------------------
void __fastcall TfrmSim::CBVisualTraitsClick(TObject *Sender)
{
if (CBVisualTraits->Checked) {
	RGTraitsMap->Enabled = true;  RGTraitsMap->Visible = true;
	edtTraitsMap_Int->Visible = true;
	if (RGTraitsMap->ItemIndex == 0) edtTraitsMap_Int->Enabled = false;
	else edtTraitsMap_Int->Enabled = true;}
else {
	RGTraitsMap->Enabled = false;  RGTraitsMap->Visible = false;
	edtTraitsMap_Int->Enabled = false; edtTraitsMap_Int->Visible = false;
}
}


//---------------------------------------------------------------------------
void __fastcall TfrmSim::CBVisualMoveClick(TObject *Sender)
{
if (CBVisualMove->Checked) edtMovtSlow->Enabled = true;
else edtMovtSlow->Enabled = false;
}

//---------------------------------------------------------------------------
void __fastcall TfrmSim::BtnOKClick(TObject *Sender)
{
landParams ppLand = pLandscape->getLandParams();
demogrParams dem = pSpecies->getDemogr();
simParams sim = paramsSim->getSim();
simView v = paramsSim->getViews();

if (frmSeeding->RG1->ItemIndex == -1
|| (ppLand.generated && frmSeeding->RGRules->ItemIndex == -1) ) {
	MessageDlg("You have NOT chosen any initialisation rule.\nPlease select one.",
		mtWarning, TMsgDlgButtons() << mbOK,0);
	return;
}
if (forceInit) {
	MessageDlg("Please specify initialisation rule",
		mtWarning, TMsgDlgButtons() << mbOK,0);
	return;
}

sim.simulation = StrToInt(edtSimNr->Text);
sim.reps = StrToInt(edtRep->Text);
sim.years = StrToInt(edtYears->Text);
sim.absorbing = CBabsorbing->Checked;

//Parameters for Environmental stochastcity
envStochParams env;
if (CBEnvStoch->Checked) {
	float ffff,minLimit,maxLimit;
	env.stoch = true;
	if (RGEnvStochType->ItemIndex == 0) env.inK = false; else env.inK = true;
	if (RGEnvStoch->ItemIndex == 0) env.local = false; else env.local = true;
	ffff = StrToFloat(edtAC->Text);
	if (ffff < 0.0 || ffff >= 1.0) {
		MessageDlg("Temporal autocorrelation must be >= 0.0 and < 1.0",
			mtError, TMsgDlgButtons() << mbOK,0);
		return;
	}
	else env.ac = ffff;
	ffff = StrToFloat(edtStd->Text);
	if (ffff <= 0.0 || ffff > 1.0) {
		MessageDlg("Amplitude must be > 0 and <= 1.0",
			mtError, TMsgDlgButtons() << mbOK,0);
		return;
	}
	else env.std = ffff;
	minLimit = StrToFloat(edtMinR->Text);
	maxLimit = StrToFloat(edtMaxR->Text);
	string msg,limittype;
	if (env.inK)
	{
		if (dem.stageStruct) limittype = " 1/b ";
		else limittype = " K ";
	}
	else {
		if (dem.stageStruct) limittype = " fecundity ";
		else limittype = " growth rate ";
	}
	// min value of env stoch in r or K must be zero if there is a gradient
	// in the same parameter
	envGradParams grad;
	grad = paramsGrad->getGradient();
	if (minLimit < 0.0) {
		msg = "Minimum" + limittype + "must be >= 0.0";
		MessageDlg(msg.c_str(), mtError, TMsgDlgButtons() << mbOK,0);
		return;
	}
	else {
		if (grad.gradient) {
			if ((grad.gradType == 1 && env.inK && minLimit > 0.0)
			||  (grad.gradType == 2 && !env.inK && minLimit > 0.0))
			{
				msg = "Minimum" + limittype + "must be zero as there is a gradient in" + limittype;
				MessageDlg(msg.c_str(), mtError, TMsgDlgButtons() << mbOK,0);
				return;
			}
		}
  }
	if (maxLimit <= minLimit) {
		msg = "Maximum" + limittype + "must be > minimum" + limittype;
		MessageDlg(msg.c_str(), mtError, TMsgDlgButtons() << mbOK,0);
		return;
	}
	if (env.inK) {
		float minK,maxK;
		minK = minLimit * (((float)(ppLand.resol*ppLand.resol))/10000.0);
		maxK = maxLimit * (((float)(ppLand.resol*ppLand.resol))/10000.0);
		pSpecies->setMinMax(minK,maxK);
	}
	else {
		pSpecies->setMinMax(minLimit,maxLimit);
	}
	if (env.local) chartNoise = false;
	else chartNoise = true;
}
else {
	env.stoch = false;
	env.ac = 0.0; env.std = 1.0;
	chartNoise = false;
}

if (CBLocExt->Checked) {
	env.localExt = true;
	float locExtProb = -1.0;
	locExtProb = StrToFloat(edtLocExt->Text);
	if (locExtProb <= 0.0 || locExtProb >= 1.0) {
		MessageDlg("Local extinction probability must be > 0.0 and < 1.0",
			mtError, TMsgDlgButtons() << mbOK,0);
		return;
	}
	else
		env.locExtProb = locExtProb;
}
else {
	env.localExt = false; env.locExtProb = 0.0;
}
paramsStoch->setStoch(env);

if (CBVisualMove->Checked && edtMovtSlow->Enabled) {
	if (StrToInt(edtMovtSlow->Text) < 1) {
		MessageDlg("Movement slow factor must be > 0",
			mtError, TMsgDlgButtons() << mbOK,0);
		return;
	}
}

//Outputs ------------------------------------------------------------------
sim.outRange = CBoutRange->Checked;
sim.outOccup = CBoutOcc->Checked;
sim.outPop = CBoutPop->Checked;
sim.outInds = CBoutInd->Checked;
sim.outGenetics = CBoutGen->Checked;
sim.outGenType = RGGen->ItemIndex;
sim.outGenXtab = CBGenCrosstab->Checked;
sim.outTraitsCells = CBoutTraitCell->Checked;
if (ppLand.patchModel) sim.outTraitsRows = false;
else sim.outTraitsRows = CBoutTraitRow->Checked;
sim.outConnect = CBoutConnect->Checked;
//sim.outStartRange = StrToInt(edtStartRange->Text);
//sim.outStartOcc   = StrToInt(edtStartOcc->Text);
sim.outStartPop   = StrToInt(edtStartPop->Text);
sim.outStartInd   = StrToInt(edtStartInd->Text);
sim.outStartGenetic = StrToInt(edtStartGen->Text);
sim.outStartTraitCell = StrToInt(edtStartTraitCell->Text);
sim.outStartTraitRow  = StrToInt(edtStartTraitRow->Text);
sim.outStartConn  = StrToInt(edtStartConn->Text);
sim.outIntRange = StrToInt(edtFreqRange->Text);
sim.outIntOcc   = StrToInt(edtFreqOcc->Text);
sim.outIntPop   = StrToInt(edtFreqPop->Text);
sim.outIntInd   = StrToInt(edtFreqInd->Text);
sim.outIntGenetic = StrToInt(edtFreqGen->Text);
sim.outIntTraitCell = StrToInt(edtFreqTraitCell->Text);
sim.outIntTraitRow  = StrToInt(edtFreqTraitRow->Text);
sim.outIntConn  = StrToInt(edtFreqConn->Text);

sim.saveMaps = RGMap->ItemIndex;
//Map outputs parameters
if (sim.saveMaps) {
	sim.mapInt = StrToInt(edtMap_interval->Text);
//	sim.saveInitMap = CBInitMap->Checked;
	if (ppLand.spDist) // species initial distribution exists
		sim.drawLoaded = CBDrawInit->Checked;
	else sim.drawLoaded = false;
}

sim.saveTraitMaps = RGTraitsMap->ItemIndex;
if (sim.saveTraitMaps) sim.traitInt = StrToInt(edtTraitsMap_Int->Text);

paramsSim->setSim(sim);

//Visualisations -----------------------------------------------------------
v.viewLand = CBVisualLand->Checked;
v.viewPop = CBVisualPop->Checked;
v.viewGraph = CBVisualGraph->Checked;
v.viewGrad = CBVisualGrad->Checked;
v.viewTraits = CBVisualTraits->Checked;
v.viewPaths = CBVisualMove->Checked;
v.slowFactor = StrToInt(edtMovtSlow->Text);
paramsSim->setViews(v);

MemoLine("Simulation parameters saved.");
MemoLine("Simulation ready to run.");

readyToRun = true;
Close();

}

//---------------------------------------------------------------------------

void __fastcall TfrmSim::CBoutPopClick(TObject *Sender)
{
if (CBoutPop->Checked) {
	edtStartPop->Enabled = true; edtFreqPop->Enabled = true;
}
else {
	edtStartPop->Enabled = false; edtFreqPop->Enabled = false;
}
}
//---------------------------------------------------------------------------

void __fastcall TfrmSim::CBoutIndClick(TObject *Sender)
{
if (CBoutInd->Checked) {
	edtStartInd->Enabled = true; edtFreqInd->Enabled = true;
}
else {
	edtStartInd->Enabled = false; edtFreqInd->Enabled = false;
}
}
//---------------------------------------------------------------------------

void __fastcall TfrmSim::CBoutGenClick(TObject *Sender)
{
if (CBoutGen->Checked) {
	edtStartGen->Enabled = true; edtFreqGen->Enabled = true;
	if (pSpecies->stageStructured()) {
		RGGen->Enabled = true; RGGen->Visible = true;
	}
	else {
		RGGen->Enabled = false; RGGen->Visible = false;
	}
	CBGenCrosstab->Enabled = true; CBGenCrosstab->Visible = true;
}
else {
	edtStartGen->Enabled = false; edtFreqGen->Enabled = false;
	RGGen->Enabled = false; RGGen->Visible = false;
	CBGenCrosstab->Enabled = false; CBGenCrosstab->Visible = false;
}
}
//---------------------------------------------------------------------------

void __fastcall TfrmSim::CBoutTraitCellClick(TObject *Sender)
{
if (CBoutTraitCell->Checked) {
	edtStartTraitCell->Enabled = true; edtFreqTraitCell->Enabled = true;
}
else {
	edtStartTraitCell->Enabled = false; edtFreqTraitCell->Enabled = false;
}
}
//---------------------------------------------------------------------------

void __fastcall TfrmSim::CBoutTraitRowClick(TObject *Sender)
{
if (CBoutTraitRow->Checked) {
	edtStartTraitRow->Enabled = true; edtFreqTraitRow->Enabled = true;
}
else {
	edtStartTraitRow->Enabled = false; edtFreqTraitRow->Enabled = false;
}
}
//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


