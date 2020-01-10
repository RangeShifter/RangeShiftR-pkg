//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "FormMove.h"
#include "FormSpecies.h"
#include "FormVisualCost.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TfrmMove *frmMove;

string costmapname;

//---------------------------------------------------------------------------
__fastcall TfrmMove::TfrmMove(TComponent* Owner)
	: TForm(Owner)
{
Left = 100;
Top = 60;
}

//---------------------------------------------------------------------------

void __fastcall TfrmMove::refresh(void)
{
landParams ppLand = pLandscape->getLandParams();
trfrRules trfr = pSpecies->getTrfr();

#if RSDEBUG
//if (costsSet) MemoLine("costsSet = true");
//else MemoLine("costsSet = false");
#endif

if (ppLand.generated) { // importing cost map is not permitted
	CBCosts->Enabled = false; CBCosts->Checked = false;
}
else {
	CBCosts->Enabled = true;
	if (!costsSet) { // any previous cost map may be invalid
		trfr.costMap = false;
		pSpecies->setTrfr(trfr);
		CBCosts->Checked = false;
		CBVisualCosts->Checked = false; CBVisualCosts->Visible = false;
	}
}

if (trfr.moveType == 1) {
	RGMoveModel->ItemIndex = 0; // SMS
	trfrSMSTraits sms = pSpecies->getSMSTraits();
	if (sms.goalType == 2) { // dispersal bias
		edtAlphaDB->Enabled = true; edtBetaDB->Enabled = true;
	}
	else {
		edtAlphaDB->Enabled = false; edtBetaDB->Enabled = false;
	}
}
else RGMoveModel->ItemIndex = 1; // CRW

RGstepM->Enabled = true; 
#if TEMPMORT              
if (trfr.smType == 1) RGstepM->ItemIndex = 1; // habitat-dependent mortality
#else
if (trfr.habMort) RGstepM->ItemIndex = 1; // habitat-dependent mortality
#endif // TEMPMORT
else RGstepM->ItemIndex = 0; // constant mortality
#if TEMPMORT              
if (trfr.moveType == 2 && !trfr.smType == 1)
#else
if (trfr.moveType == 2 && !trfr.habMort) 
#endif // TEMPMORT
{
	// CRW and constant per step mortality - do not show costs/mortality matrix
	PanelHab->Enabled = false;   PanelHab->Visible = false;
}
else {
	PanelHab->Enabled = true;   PanelHab->Visible = true;
}
if (trfr.costMap) CBCosts->Checked = true; else CBCosts->Checked = false;

if (!ppLand.generated && ppLand.rasterType != 0) {
	// for imported habitat cover or quality, user is required to upload costs map
	RGstepM->ItemIndex = 0;
	RGstepM->Enabled = false;
	CBCosts->Enabled = false; CBCosts->Visible = true; CBCosts->Checked = true;
	CBVisualCosts->Enabled = true; CBVisualCosts->Visible = true;
}

//if (pSpecies->stageStructed()) {
//	// individual variability in CRW is not allowed
//	CBIndVarCRW->Checked = false; CBIndVarCRW->Enabled = false;
//}
//else {
	CBIndVarCRW->Enabled = true;
//}

setHabTable();

if (!costsSet || !trfr.costMap) frmVisualCost->Close();

}


//---------------------------------------------------------------------------

void __fastcall TfrmMove::setHabTable(void)
{

//MessageDlg("Starting setHabTable()",mtWarning,TMsgDlgButtons() << mbOK,0);

SGhabCosts->Enabled = true; SGhabCosts->Visible = true;

int habcode,habcost;
int mortcol = 0;
landParams ppLand = pLandscape->getLandParams();
trfrRules trfr = pSpecies->getTrfr();

if (RGMoveModel->ItemIndex == 0 && !CBCosts->Checked) { // SMS and NOT cost map
	// costs column is required
	if (RGstepM->ItemIndex == 0) { // constant mortality
		SGhabCosts->ColCount = 2;
	}
	else { // habitat-dependent mortality
		SGhabCosts->ColCount = 3;
		SGhabCosts->Cells[2][0] = "Step mort.";
		mortcol = 2;
	}
	SGhabCosts->Cells[1][0] = "Costs";
}
else { // CRW or SMS with loaded cost map
	// costs column is not required
	if (RGstepM->ItemIndex == 0) { // constant mortality
		// neither column is required
		SGhabCosts->Enabled = false; SGhabCosts->Visible = false;
		return;
	}
	else { // habitat-dependent mortality
		SGhabCosts->ColCount = 2;
		SGhabCosts->Cells[1][0] = "Step mort.";
		mortcol = 1;
	}
}

if (ppLand.generated) {
	SGhabCosts->RowCount = 3;
	SGhabCosts->Cells[0][1] = "Matrix";
	SGhabCosts->Cells[0][2] = "Habitat";
}
else {
	SGhabCosts->RowCount = ppLand.nHab + 1;
	for (int i = 0; i < ppLand.nHab; i++) {
		habcode = pLandscape->getHabCode(i);
		SGhabCosts->Cells[0][i+1] = habcode;
	}
}

if (RGMoveModel->ItemIndex == 0 && !CBCosts->Checked) { // SMS and NOT cost map
	// populate habitat cost values
	for (int i = 0; i < ppLand.nHab; i++) {
		if (costsSet) {
			habcost = pSpecies->getHabCost(i);
		}
		else {
			habcost = pLandscape->getHabCode(i);
		}
		if (habcost < 1) habcost = 1;
		SGhabCosts->Cells[1][i+1] = habcost;
	}
}
if (mortcol > 0) { // habitat-dependent mortality
	// populate habitat-dependent mortality values
	for (int i = 0; i < ppLand.nHab; i++) {
		SGhabCosts->Cells[mortcol][i+1] = pSpecies->getHabMort(i);
  }
}

}

//---------------------------------------------------------------------------

void __fastcall TfrmMove::RGMoveModelClick(TObject *Sender)
{
landParams ppLand = pLandscape->getLandParams();
demogrParams dem = pSpecies->getDemogr();
trfrRules trfr = pSpecies->getTrfr();

//if (dem.stageStruct)
//	MessageDlg("Stage structure is TRUE",mtWarning, TMsgDlgButtons() << mbRetry,0);
//else
//	MessageDlg("Stage structure is FALSE",mtWarning, TMsgDlgButtons() << mbRetry,0);

switch (RGMoveModel->ItemIndex) {

case 0: // SMS
	PanelHab->Enabled = true; PanelHab->Visible = true;
	PanelSMS->Enabled = true; PanelSMS->Visible = true;
	PanelRW->Enabled = false; PanelRW->Visible = false;
	CBCosts->Enabled = true; CBCosts->Visible = true;
//	CBCosts->Checked = false;
	if (StrToInt(edtGoalType->Text) == 2) { // dispersal bias
		edtAlphaDB->Enabled = true; edtBetaDB->Enabled = true;
	}
	else {
		edtAlphaDB->Enabled = false; edtBetaDB->Enabled = false;
	}

	if (RGstepM->ItemIndex == 0) { // constant per-step mortality
		edtSM->Enabled = true; edtSM->Visible = true;
		if (!ppLand.generated && ppLand.rasterType != 0) {
			// real landscape and not habitat codes - cost map must be supplied
			CBCosts->Enabled = false; CBCosts->Visible = true; CBCosts->Checked = true;
			CBVisualCosts->Enabled = true; CBVisualCosts->Visible = true;
		}
	}
	else { // habitat-dependent per-step mortality
		edtSM->Enabled = false; edtSM->Visible = false;
	}
	// individual variability is not implemented for SMS
	trfr.indVar = false;
	pSpecies->setTrfr(trfr);
	break;

case 1: // CRW
	PanelHab->Enabled = false; PanelHab->Visible = false;
	PanelRW->Enabled = true; PanelRW->Visible = true;
	PanelSMS->Enabled = false; PanelSMS->Visible = false;
	CBCosts->Enabled = false; CBCosts->Visible = false; CBCosts->Checked = false;
//	if (dem.stageStruct) {
//		// individual variability not allowed
//		CBIndVarCRW->Checked = false; CBIndVarCRW->Enabled = false;
//	}
//	else {
		CBIndVarCRW->Enabled = true;
//	}
	if (CBIndVarCRW->Checked == false) {
		edtStepLength->Enabled = true; edtRho->Enabled = true;
		edtStepLMean->Enabled = false; edtStepLSD->Enabled = false;
		edtRhoMean->Enabled =false; edtRhoSD->Enabled = false;
	}
	else{
		edtStepLength->Enabled = false; edtRho->Enabled = false;
		edtStepLMean->Enabled = true; edtStepLSD->Enabled = true;
		edtRhoMean->Enabled =true; edtRhoSD->Enabled = true;
	}
	if (RGstepM->ItemIndex == 1) {
		PanelHab->Enabled = true; PanelHab->Visible = true;
		edtStepMort->Enabled = false; edtSM->Visible = false;
	}
	else PanelHab->Enabled = false;
	break;

}

setHabTable();

}

//---------------------------------------------------------------------------
void __fastcall TfrmMove::RGstepMClick(TObject *Sender)
{

if (RGstepM->ItemIndex == 0) { // constant step mortality
	if (RGMoveModel->ItemIndex == 0) { // SMS
		edtSM->Visible = true; edtSM->Enabled = true;
		if (CBCosts->Checked) { // import costs map - table is not required
			PanelHab->Enabled = false; PanelHab->Visible = false;
		}
		else { // costs input in table
			PanelHab->Enabled = true; PanelHab->Visible = true;
		}
	}
	else { // CRW
		PanelHab->Enabled = false; PanelHab->Visible = false;
		edtStepMort->Enabled = true; edtStepMort->Visible = true;
	}
}
else { // habitat-dependent step mortality
	PanelHab->Enabled = true; PanelHab->Visible = true;
	if (RGMoveModel->ItemIndex == 0) { // SMS
		edtSM->Visible = false; edtSM->Enabled = false;
	}
	else { // CRW
	 edtStepMort->Enabled = false; edtStepMort->Visible = false;
	}
}

setHabTable();

}
//---------------------------------------------------------------------------

void __fastcall TfrmMove::edtGoalTypeChange(TObject *Sender)
{
if (StrToInt(edtGoalType->Text) == 2) { // dispersal bias
	if (CBIndVarSMS->Checked) {
		PanelDBIndVar->Enabled = true; PanelDBIndVar->Visible = true;
		edtAlphaDB->Enabled = false; edtBetaDB->Enabled = false;
	}
	else {
		PanelDBIndVar->Enabled = false; PanelDBIndVar->Visible = false;
		edtAlphaDB->Enabled = true; edtBetaDB->Enabled = true;
	}
}
else {
	PanelDBIndVar->Enabled = false; PanelDBIndVar->Visible = false;
	edtAlphaDB->Enabled = false; edtBetaDB->Enabled = false;
}
}

//---------------------------------------------------------------------------
void __fastcall TfrmMove::CBIndVarSMSClick(TObject *Sender)
{
if (CBIndVarSMS->Checked) {
	edtDP->Enabled = false; edtGoalBias->Enabled = false;
	edtAlphaDB->Enabled = false; edtBetaDB->Enabled = false;
	PanelDPGBIndVar->Enabled = true; PanelDPGBIndVar->Visible = true;
	if (StrToInt(edtGoalType->Text) == 2) { // dispersal bias
		PanelDBIndVar->Enabled = true; PanelDBIndVar->Visible = true;
	}
	else {
		PanelDBIndVar->Enabled = false; PanelDBIndVar->Visible = false;
	}
}
else {
	edtDP->Enabled = true; edtGoalBias->Enabled = true;
	edtAlphaDB->Enabled = true; edtBetaDB->Enabled = true;
	PanelDPGBIndVar->Enabled = false; PanelDPGBIndVar->Visible = false;
	PanelDBIndVar->Enabled = false;   PanelDBIndVar->Visible = false;
	if (StrToInt(edtGoalType->Text) == 2) { // dispersal bias
		edtAlphaDB->Enabled = true; edtBetaDB->Enabled = true;
	}
	else {
		edtAlphaDB->Enabled = false; edtBetaDB->Enabled = false;
	}
}
}

//---------------------------------------------------------------------------
void __fastcall TfrmMove::CBIndVarCRWClick(TObject *Sender)
{
if (CBIndVarCRW->Checked) {
	edtStepLength->Enabled = false; edtRho->Enabled = false;
	edtStepLMean->Enabled = true;   edtStepLMean->Visible = true;
	edtStepLSD->Enabled = true;     edtStepLSD->Visible = true;
	edtStepLScale->Enabled = true;  edtStepLScale->Visible = true;
	edtRhoMean->Enabled = true;     edtRhoMean->Visible = true;
	edtRhoSD->Enabled = true;       edtRhoSD->Visible = true;
	edtRhoScale->Enabled = true;    edtRhoScale->Visible = true;
}
else {
	edtStepLength->Enabled = true; edtRho->Enabled = true;
	edtStepLMean->Enabled = false;   edtStepLMean->Visible = false;
	edtStepLSD->Enabled = false;     edtStepLSD->Visible = false;
	edtStepLScale->Enabled = false;  edtStepLScale->Visible = false;
	edtRhoMean->Enabled = false;     edtRhoMean->Visible = false;
	edtRhoSD->Enabled = false;       edtRhoSD->Visible = false;
	edtRhoScale->Enabled = false;    edtRhoScale->Visible = false;
}
}
/*
//---------------------------------------------------------------------------
void __fastcall TfrmMove::CBEvolRWClick(TObject *Sender)
{
 if (CBEvolRW->Checked){
	edtMutRW->Enabled = true;
	edtMutRW->Visible = true;
	edtMutSizeSL->Enabled = true;
	edtMutSizeSL->Visible = true;
	edtMutSizeRho->Enabled = true;
	edtMutSizeRho->Visible = true;
 }
 else{
	edtMutRW->Enabled = false;
	edtMutRW->Visible = false;
	edtMutSizeSL->Enabled = false;
	edtMutSizeSL->Visible = false;
	edtMutSizeRho->Enabled = false;
	edtMutSizeRho->Visible = false;
 }
}
*/
//---------------------------------------------------------------------------
void __fastcall TfrmMove::CBCostsClick(TObject *Sender)
{
simView v = paramsSim->getViews();
if (CBCosts->Checked) {
	CBVisualCosts->Enabled = true; CBVisualCosts->Visible = true;
}
else { // not imported cost map
	CBVisualCosts->Enabled = false; CBVisualCosts->Visible = false;
	v.viewCosts = false;
	paramsSim->setViews(v);
}

setHabTable();

// as source of costs has been changed, any existing costs in the landscape must be cleared
pLandscape->resetCosts();
}

//---------------------------------------------------------------------------
void __fastcall TfrmMove::BtnOKClick(TObject *Sender)
{
int rtncode;
Cell *pCell;
landParams ppLand = pLandscape->getLandParams();
trfrRules trfr = pSpecies->getTrfr();
trfrMovtTraits move = pSpecies->getMovtTraits();
trfrSMSParams smsparams;
trfrCRWParams crwparams;
movtModelOK = false;
bool costsEntered = false;
float tttt0,tttt1;

string msg;
string msgdp = "Directional persistence ";
string msggb = "Goal bias ";
string msgalpha = "Alpha DB ";
string msgbeta = "Beta DB ";
string msgsteplgh = "Step length ";
string msgrho = "Step correlation ";
string msgmean = "mean ";
string msgsd = "s.d. ";
string msgscale = "scaling factor ";
string msggt0 = "must be greater than zero ";
string msgnotlt = "may not be less than ";
string msg01 = "must be in the range 0 to 1 ";

trfr.moveType = RGMoveModel->ItemIndex + 1; // so that ... 1 = SMS, 2 = CRW
if (trfr.moveType == 1) { // SMS
	if (RGstepM->ItemIndex == 0) { // constant mortality
		if (StrToFloat(edtSM->Text) >= 0.0 && StrToFloat(edtSM->Text) < 1.0 ) {
			move.stepMort = StrToFloat(edtSM->Text);
		}
		else {
			MessageDlg("Step mortality must be >= 0 and < 1",mtError,TMsgDlgButtons() << mbOK,0);
			return;
		}
	}
	if (StrToInt(edtPR->Text) >= 1) {
		move.pr = StrToInt(edtPR->Text);
	}
	else {
		MessageDlg("Error in Perceptual range",mtError,TMsgDlgButtons() << mbOK,0);
		return;
	}
	if (StrToInt(edtPRMethod->Text) >= 1 && StrToInt(edtPRMethod->Text) <= 3) {
		move.prMethod = StrToInt(edtPRMethod->Text);
	}
	else {
		MessageDlg("Error in PR method",mtError,TMsgDlgButtons() << mbOK,0);
		return;
	}
	if (StrToFloat(edtDP->Text) >= 1.0) {
		move.dp = StrToFloat(edtDP->Text);
	}
	else {
		MessageDlg("Error in Directional persistence",mtError,TMsgDlgButtons() << mbOK,0);
		return;
	}
	if (StrToInt(edtMemorySize->Text) >= 1 && StrToInt(edtMemorySize->Text) <= 14) {
		move.memSize = StrToInt(edtMemorySize->Text);
	}
	else {
		MessageDlg("Error in Memory size",mtError,TMsgDlgButtons() << mbOK,0);
		return;
	}
	if (StrToFloat(edtGoalBias->Text) >= 1.0) {
		move.gb = StrToFloat(edtGoalBias->Text);
	}
	else {
		MessageDlg("Error in Goal bias",mtError,TMsgDlgButtons() << mbOK,0);
		return;
	}
	if (StrToInt(edtGoalType->Text) == 0 || StrToInt(edtGoalType->Text) == 2) {
		move.goalType = StrToInt(edtGoalType->Text);
	}
	else {
		if (StrToInt(edtGoalType->Text) == 1) {
			MessageDlg("Goal type 1 is not implemented in this version",
				mtError,TMsgDlgButtons() << mbOK,0);
		}
		else {
			MessageDlg("Error in Goal type",mtError,TMsgDlgButtons() << mbOK,0);
		}
		return;
	}
	if (StrToInt(edtGoalType->Text) == 2) { // dispersal bias
		if (StrToFloat(edtAlphaDB->Text) > 0.0) {
			move.alphaDB = StrToFloat(edtAlphaDB->Text);
		}
		else {
			msg = msgalpha + msggt0;
			MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
			return;
		}
		if (StrToInt(edtBetaDB->Text) > 0) {
			move.betaDB = StrToInt(edtBetaDB->Text);
		}
		else {
			msg = msgbeta + msggt0;
			MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
			return;
		}
	}
	if (CBIndVarSMS->Checked) {
		if (StrToFloat(edtDPMean->Text) < 1.0) {
			msg = msgdp + msgmean + msgnotlt + "1.0";
			MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
			return;
		}
		tttt0 = StrToFloat(edtDPSD->Text);
		if (tttt0 <= 0.0) {
			msg = msgdp + msgsd + msggt0;
			MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
			return;
		}
		tttt1 = StrToFloat(edtDPScale->Text);
		if (tttt1 <= 0.0) {
			msg = msgdp + msgscale + msggt0;
			MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
			return;
		}
		if (tttt1 < tttt0) {
			msg = msgdp + msgscale + msgnotlt + msgsd;
			MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
			return;
		}
		if (StrToFloat(edtGBMean->Text) < 1.0) {
			msg = msggb + msgmean + msgnotlt + "1.0";
			MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
			return;
		}
		tttt0 = StrToFloat(edtGBSD->Text);
		if (tttt0 <= 0.0) {
			msg = msggb + msgsd + msggt0;
			MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
			return;
		}
		tttt1 = StrToFloat(edtGBScale->Text);
		if (tttt1 <= 0.0) {
			msg = msggb + msgscale + msggt0;
			MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
			return;
		}
		if (tttt1 < tttt0) {
			msg = msggb + msgscale + msgnotlt + msgsd;
			MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
			return;
		}
		if (StrToInt(edtGoalType->Text) == 2) { // dispersal bias
			if (StrToFloat(edtAlphaDBMean->Text) <= 0.0) {
				msg = msgalpha + msgmean + msggt0;
				MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
				return;
			}
			tttt0 = StrToFloat(edtAlphaDBSD->Text);
			if (tttt0 <= 0.0) {
				msg = msgalpha + msgsd + msggt0;
				MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
				return;
			}
			tttt1 = StrToFloat(edtAlphaDBScale->Text);
			if (tttt1 <= 0.0) {
				msg = msgalpha + msgscale + msggt0;
				MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
				return;
			}
			if (tttt1 < tttt0) {
				msg = msgalpha + msgscale + msgnotlt + msgsd;
				MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
				return;
			}
			if (StrToFloat(edtBetaDBMean->Text) < 1.0) {
				msg = msgbeta + msgmean + msgnotlt + "1.0";
				MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
				return;
			}
			tttt0 = StrToFloat(edtBetaDBSD->Text);
			if (tttt0 <= 0.0) {
				msg = msgbeta + msgsd + msggt0;
				MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
				return;
			}
			tttt1 = StrToFloat(edtBetaDBScale->Text);
			if (tttt1 <= 0.0) {
				msg = msgbeta + msgscale + msggt0;
				MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
				return;
			}
			if (tttt1 < tttt0) {
				msg = msgalpha + msgscale + msgnotlt + msgsd;
				MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
				return;
			}
		}
	} // end of (CBIndVarSMS->Checked)
	move.straigtenPath = CBstraightPath->Checked;
//	maxcost = 0;
	if (!ppLand.generated && ppLand.rasterType != 0) // real landscape and not habitat codes
		trfr.costMap = true; // cost map is compulsory
	else
		trfr.costMap = CBCosts->Checked;
	if (CBIndVarSMS->Checked) {
		trfrScales scales = pSpecies->getTrfrScales();
		trfr.indVar = true;
		smsparams.dpMean  = StrToFloat(edtDPMean->Text);
		smsparams.dpSD    = StrToFloat(edtDPSD->Text);
		scales.dpScale 		= StrToFloat(edtDPScale->Text);
		smsparams.gbMean  = StrToFloat(edtGBMean->Text);
		smsparams.gbSD    = StrToFloat(edtGBSD->Text);
		scales.gbScale 	  = StrToFloat(edtGBScale->Text);
		if (StrToInt(edtGoalType->Text) == 2) { // dispersal bias
			smsparams.alphaDBMean  = StrToFloat(edtAlphaDBMean->Text);
			smsparams.alphaDBSD    = StrToFloat(edtAlphaDBSD->Text);
			scales.alphaDBScale    = StrToFloat(edtAlphaDBScale->Text);
			smsparams.betaDBMean   = StrToFloat(edtBetaDBMean->Text);
			smsparams.betaDBSD     = StrToFloat(edtBetaDBSD->Text);
			scales.betaDBScale     = StrToFloat(edtBetaDBScale->Text);
		}
		pSpecies->setSMSParams(0,0,smsparams);
		pSpecies->setTrfrScales(scales);
	}
	else {
		trfr.indVar = false;
	}
} // end of SMS
else { // CRW
	trfr.costMap = false; // cost map is not applicable
	if (CBIndVarCRW->Checked) {
		if (StrToFloat(edtStepLMean->Text) <= 0.0) {
			msg = msgsteplgh + msgmean + msggt0;
			MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
			return;
		}
		tttt0 = StrToFloat(edtStepLSD->Text);
		if (tttt0 <= 0.0) {
			msg = msgsteplgh + msgsd + msggt0;
			MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
			return;
		}
		tttt1 = StrToFloat(edtStepLScale->Text);
		if (tttt1 <= 0.0) {
			msg = msgsteplgh + msgscale + msggt0;
			MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
			return;
		}
		if (tttt1 < tttt0) {
			msg = msgsteplgh + msgscale + msgnotlt + msgsd;
			MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
			return;
		}
		tttt0 = StrToFloat(edtRhoMean->Text);
		if (tttt0 <= 0.0 || tttt0 >= 1.0) {
			msg = msgrho + msgmean + msg01;
			MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
			return;
		}
		tttt0 = StrToFloat(edtRhoSD->Text);
		if (tttt0 <= 0.0 || tttt0 >= 1.0) {
			msg = msgrho + msgsd + msg01;
			MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
			return;
		}
		tttt1 = StrToFloat(edtRhoScale->Text);
		if (tttt1 <= 0.0 || tttt1 >= 1.0) {
			msg = msgrho + msgscale + msg01;
			MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
			return;
		}
		if (tttt1 < tttt0) {
			msg = msgrho + msgscale + msgnotlt + msgsd;
			MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
			return;
		}
	}
	else {
		if (StrToFloat(edtStepLength->Text) <= 0.0) {
			msg = msgsteplgh + msggt0;
			MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
			return;
		}
		if (StrToFloat(edtRho->Text) <= 0.0 || StrToFloat(edtRho->Text) >= 1.0) {
			MessageDlg("Step correlation must be > 0 and < 1",mtError,TMsgDlgButtons() << mbOK,0);
			return;
		}
	}
	if (RGstepM->ItemIndex == 0) {
		if (StrToFloat(edtStepMort->Text) >= 0.0 && StrToFloat(edtStepMort->Text) < 1.0 ) {
			move.stepMort = StrToFloat(frmMove->edtStepMort->Text);
		}
		else {
			MessageDlg("Step mortality must be >= 0 and < 1",mtError,TMsgDlgButtons() << mbOK,0);
			return;
		}
	}
	if (CBIndVarCRW->Checked) {
		trfr.indVar = true;
		crwparams.stepLgthMean = StrToFloat(edtStepLMean->Text);
		crwparams.stepLgthSD = StrToFloat(edtStepLSD->Text);
		crwparams.rhoMean = StrToFloat(edtRhoMean->Text);
		crwparams.rhoSD = StrToFloat(edtRhoSD->Text);
		pSpecies->setCRWParams(0,0,crwparams);
		trfrScales scales = pSpecies->getTrfrScales();
		scales.stepLScale = StrToFloat(edtStepLScale->Text);
		scales.rhoScale   = StrToFloat(edtRhoScale->Text);
		pSpecies->setTrfrScales(scales);
	}
	else {
		trfr.indVar = false;
		move.stepLength = StrToFloat(edtStepLength->Text);
		move.rho = StrToFloat(edtRho->Text);
	}
}

// record habitat-dependent step mortality indicator
#if TEMPMORT              
trfr.smType = RGstepM->ItemIndex;
#else
if (RGstepM->ItemIndex == 0) trfr.habMort = false; else trfr.habMort = true;
#endif // TEMPMORT

// Check the habitat-dependent costs / step mortality table, if activated
bool checktable = false;
int nCols = 0;

if (trfr.moveType == 1) { // SMS
	if (trfr.costMap) {
#if TEMPMORT              
		if (trfr.smType == 1)
#else
		if (trfr.habMort) 
#endif // TEMPMORT
		{
			checktable = true; nCols = 1;
		}
	}
	else { // no cost map
		checktable = true; nCols = 1; costsEntered = true;
#if TEMPMORT              
		if (trfr.smType == 1) nCols = 2;
#else
		if (trfr.habMort) nCols = 2;
#endif // TEMPMORT
	}
}
else { // CRW
#if TEMPMORT              
	if (trfr.smType == 1) 
#else
	if (trfr.habMort)
#endif // TEMPMORT 
	{
		checktable = true; nCols = 1;
	}
}
bool cost_error = false;
bool mort_error = false;
if (checktable) {
	int cost; double mort;
	for (int i = 0; i < ppLand.nHab; i++) {
		cost = 1; mort = 0.0;
		if (nCols == 2) {
			cost = (SGhabCosts->Cells[1][i+1]).ToInt();
			mort = (SGhabCosts->Cells[2][i+1]).ToDouble();
		}
		else { // one column
			if ((trfr.moveType == 1 && trfr.costMap) || trfr.moveType == 2)
				mort = (SGhabCosts->Cells[1][i+1]).ToDouble();
			else
				cost = (SGhabCosts->Cells[1][i+1]).ToInt();
		}
		if (cost < 1) cost_error = true;
		if (mort < 0.0 || mort >= 1.0) mort_error = true;
	}
}
if (cost_error) {
	MessageDlg("Error in table: costs may not be less than 1"
				,mtError,TMsgDlgButtons() << mbOK,0);
	return;
}
if (mort_error) {
	MessageDlg("Error in table: step mortality must be >= 0.0 and < 1.0"
				,mtError,TMsgDlgButtons() << mbOK,0);
	return;
}

if (!ppLand.generated) {
	pLandscape->resetCosts(); // in case costs, PR or PR method have been changed
}
pSpecies->setMovtTraits(move);
if (costsEntered) costsSet = true;
pSpecies->setTrfr(trfr);
simView v = paramsSim->getViews();

if (trfr.costMap) { // read the costs map
	bool fileSelected = false;
	OpenDialog1->InitialDir = paramsSim->getDir(1).c_str();
	OpenDialog1->Title = "Select costs map";
	if (OpenDialog1->Execute()) {
		costmapname = AnsiString(OpenDialog1->FileName).c_str();
		fileSelected = true;
	}
	else {
		MessageDlg("The costs file selection has been cancelled",
			mtWarning, TMsgDlgButtons() << mbOK,0);
		return;
	}
	rasterdata costsraster;
	costsraster.ok = false;
	if (fileSelected) {
		// check selected file
		costsraster = CheckRasterFile(costmapname);
		if (!costsraster.ok) {
			MessageDlg("Format error in header lines of the selected costs file",
				mtError, TMsgDlgButtons() << mbOK,0);
			return;
		}
		if (costsraster.ncols != ppLand.maxX+1 || costsraster.nrows != ppLand.maxY+1) {
			MessageDlg("The dimensions of the costs map differ from those of the landscape"
				,mtError, TMsgDlgButtons() << mbOK,0);
			return;
		}
		landOrigin origin = pLandscape->getOrigin();
		if (costsraster.xllcorner != origin.minEast || costsraster.yllcorner != origin.minNorth) {
			MessageDlg("The origin of the costs map differs from the landscape origin"
				,mtError, TMsgDlgButtons() << mbOK,0);
			return;
		}
		// read costs file
		rtncode = pLandscape->readCosts(costmapname);
		if (rtncode <= 0) {
			string msgReadError = "Error returned from reading costs file: "
				+ Int2Str(rtncode);
			MessageDlg(msgReadError.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
			return;
		}
		// costs file read successfully
		landParams ppLand = pLandscape->getLandParams();
		landPix p = pLandscape->getLandPix();
		if (CBVisualCosts->Checked) {
			v.viewCosts = true;
			frmVisualCost->Show();
			// setting image
			frmVisualCost->CostImage->Height = (ppLand.maxY + 1) * p.pix;
			frmVisualCost->CostImage->Width  = (ppLand.maxX + 1) * p.pix;
			frmVisualCost->PaintBox1->Height = frmVisualCost->CostImage->Height;
			frmVisualCost->PaintBox1->Width  = frmVisualCost->CostImage->Width;
			frmVisualCost->CostScrollBox->VertScrollBar->Range = frmVisualCost->CostImage->Height;
			frmVisualCost->CostScrollBox->HorzScrollBar->Range = frmVisualCost->PaintBox1->Width;
			bmpCosts = new Graphics::TBitmap();
			canvasCosts = bmpCosts->Canvas;
			bmpCosts->Height = frmVisualCost->CostImage->Height;
			bmpCosts->Width  = frmVisualCost->CostImage->Width;
		}
		else
			v.viewCosts = false;
	}
}
else {
	v.viewCosts = false;
}
paramsSim->setViews(v);

if (v.viewCosts) { // draw the cost map
	double maxcost = std::log((double)rtncode);
	int cost,col;
	landPix p = pLandscape->getLandPix();
	for (int y = ppLand.maxY; y>-1; y--){
		for (int x = 0; x <= ppLand.maxX; x++){
			pCell = pLandscape->findCell(x,y);
//			if (pCell->getHabIndex() < 0)
			if (pCell == 0)
				cost = -99; // nodata cell
			else cost = pCell->getCost();
#if RSDEBUG
//string msg = "x=" + Int2Str(x) + " y=" + Int2Str(y)  + " cost=" + Int2Str(cost);
//MemoLine(msg.c_str());
#endif
			if (cost < 1) { // nodata cell (or invalid)
				canvasCosts->Brush->Color = (TColor)RGB(176,196,222);
			}
			else {
				col = 255 - (int)((std::log((double)(cost))*255)/maxcost);
				canvasCosts->Brush->Color = (TColor)RGB(col,col,col);
			}
			canvasCosts->FillRect(Rect((int)(x*p.pix),(int)((ppLand.maxY-y)*p.pix),
				(int)((x+1)*p.pix),(int)((ppLand.maxY-y+1)*p.pix)));
		}
	}
	frmVisualCost->CostImage->Picture->Bitmap = bmpCosts;
	frmVisualCost->Repaint();
}
else {
	frmVisualCost->Close();
}

// read the table again and store costs / mortalities
if (checktable) {
	int cost; double mort;
	for (int i = 0; i < ppLand.nHab; i++) {
		cost = 1; mort = 0.0;
		if (nCols == 2) {
			cost = (SGhabCosts->Cells[1][i+1]).ToInt();
			pSpecies->setHabCost(i,cost);
			mort = (SGhabCosts->Cells[2][i+1]).ToDouble();
			pSpecies->setHabMort(i,mort);
		}
		else { // one column
			if ((trfr.moveType == 1 && trfr.costMap) || trfr.moveType == 2) {
				mort = (SGhabCosts->Cells[1][i+1]).ToDouble();
				pSpecies->setHabMort(i,mort);
			}
			else {
				cost = (SGhabCosts->Cells[1][i+1]).ToInt();
				pSpecies->setHabCost(i,cost);
			}
		}
	}
}

movtModelOK = true;
Close();

}

//---------------------------------------------------------------------------
void __fastcall TfrmMove::BtnCancelClick(TObject *Sender)
{
frmMove->Close();
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


