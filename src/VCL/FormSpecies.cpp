//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "Version.h"
#include "FormSpecies.h"

#if RSWIN64
#include<iostream>
#include<iomanip>
#else
#include<iostream.h>
#include<iomanip.h>
#endif
#include<string>

//---------------------------------------------------------------------------

#pragma package(smart_init)
#pragma resource "*.dfm"
TfrmSpecies *frmSpecies;

//---------------------------------------------------------------------------

int Ns; //number of stages
bool movtModelOK = false;
bool formCancelled = false;

//---------------------------------------------------------------------------
__fastcall TfrmSpecies::TfrmSpecies(TComponent* Owner)
	: TForm(Owner)
{
PageControl1->ActivePageIndex = 0;
Left = 100; Top = 60;
}

//---------------------------------------------------------------------------
// Refresh the form from retained parameter values
void __fastcall TfrmSpecies::refresh(bool patchModel)
{
int habcode;
landParams ppLand = pLandscape->getLandParams();
demogrParams dem = pSpecies->getDemogr();
TObject *dummy = NULL;

formCancelled = false;

// Population panel
RGReproduction->ItemIndex = dem.repType;
if (dem.stageStruct) {
	stageParams sstruct = pSpecies->getStage();
	if (ppLand.dynamic) {
		RGaction->Enabled = true; RGaction->Visible = true;
	}
	else {
		RGaction->Enabled = false; RGaction->Visible = false;
	}
	RGSurvival->ItemIndex = sstruct.survival;
	CBFecundity->Checked = sstruct.fecDens;
	CBFecundityClick(dummy);
	CBDevelopment->Checked = sstruct.devDens;
	CBDevelopmentClick(dummy);
	CBSurvival->Checked = sstruct.survDens;
	CBSurvivalClick(dummy);
	CBweightFec->Checked = sstruct.fecStageDens;
	CBweightFecClick(dummy);
	CBweightDev->Checked = sstruct.devStageDens;
	CBweightDevClick(dummy);
	CBweightSurv->Checked = sstruct.survStageDens;
	CBweightSurvClick(dummy);
}
else {
	RGaction->Enabled = false; RGaction->Visible = false;
}
if (ppLand.generated) {
	SGhabLabel->Visible = false;
	SGhab->Visible = false;
	edtKLabel->Caption = "inds / ha "; edtKLabel->Visible = true;
	edtK->Visible = true; edtK->Enabled = true;
}
else { // imported landscape
	if (ppLand.rasterType == 2) { // landscape quality
		SGhabLabel->Visible = false;
		SGhab->Visible = false;
		edtKLabel->Caption = "inds / ha (assuming 100% quality) "; edtKLabel->Visible = true;
		edtK->Visible = true; edtK->Enabled = true;
	}
	else { // habitats
		SGhabLabel->Visible = true;
		SGhab->Visible = true;
		edtKLabel->Visible = false;
		edtK->Visible = false; edtK->Enabled = false;
		if (newLandscape) {
			SGhab->RowCount = ppLand.nHab + 1;
			SGhab->Cells[1][0] = "K (inds/ha)";
			if (ppLand.rasterType == 0) { // habitat codes
				SGhab->Cells[0][0] = "Hab. code";
				for (int i = 0; i < ppLand.nHab; i++) {
					habcode = pLandscape->getHabCode(i);
					SGhab->Cells[0][i+1] = habcode;
					SGhab->Cells[1][i+1] = "0.0";
				}
			}
			else { // % cover (rasterType == 1)
				SGhab->Cells[0][0] = "Hab. nr.";
				for (int i = 0; i < ppLand.nHab; i++){
					SGhab->Cells[0][i+1] = IntToStr(i+1).c_str();
					SGhab->Cells[1][i+1] = "0.0";
				}
			}
		}
	}
}
// Dispersal panel - check boxes
emigRules emig = pSpecies->getEmig();
if (emig.densDep) RGEmigProb->ItemIndex = 1; else RGEmigProb->ItemIndex = 0;
trfrRules trfr = pSpecies->getTrfr();
if (trfr.moveModel) RGMovements->ItemIndex = 1;
else {
	RGMovements->ItemIndex = 0;
	if (trfr.twinKern) RGKernel->ItemIndex = 1; else RGKernel->ItemIndex = 0;
}
if (trfr.distMort) RGMortality->ItemIndex = 1; else RGMortality->ItemIndex = 0;
/*
CBSexEmig->Checked = emig.sexDep;
//CBSexEmigClick(dummy);
CBStageEmig->Checked = emig.stgDep;
//CBStageEmigClick(dummy);
CBIndVarEP->Checked = emig.indVar;
//CBIndVarEPClick(dummy);
*/

// Dispersal panel - initial kernel means
if (newLandscape) {
	edtDist1->Text = ppLand.resol;
	edtDist2->Text = ppLand.resol;
	edtDist1Mean->Text = ppLand.resol; edtDist2Mean->Text = ppLand.resol;
	edtDist1SD->Text = ppLand.resol / 10; edtDist2SD->Text = ppLand.resol / 10;
	edtDist1Scale->Text = ppLand.resol / 10;
	edtDist2Scale->Text = ppLand.resol / 10;
}

// Dispersal panel - settlement options
settleRules srules = pSpecies->getSettRules(0,0);
if (srules.wait) {
	if (srules.go2nbrLocn) RGSettKern->ItemIndex = 3;
	else RGSettKern->ItemIndex = 1;
}
else {
	if (srules.go2nbrLocn) RGSettKern->ItemIndex = 2;
	else RGSettKern->ItemIndex = 0;
}
if (dem.repType == 0) { // asexual model
	CBFindMate->Enabled = false; CBFindMate->Checked = false;
}
else {
	if (CBSexSettMovt->Checked) {
		CBFindMate->Enabled = false; CBFindMate->Checked = false;
	}
	else {
		CBFindMate->Enabled = true;  CBFindMate->Checked = srules.findMate;
	}
}

//UnicodeString find0 = "Find a suitable";
UnicodeString find1;
if (patchModel) {
	find1 = "Recruit to the nearest patch (1 cell) / ";
}
else{
	find1 = "Randomly choose a suitable neigh. cell / ";
}
RGSettKern->Items->Strings[2] = find1 + "die ";
RGSettKern->Items->Strings[3] = find1 + "wait ";

settleSteps ssteps = pSpecies->getSteps(0,0);
if (ssteps.maxSteps == 99999999) RGSteps->ItemIndex = 1; else RGSteps->ItemIndex = 0;

newLandscape = false;
}

//---------------------------------------------------------------------------

void __fastcall TfrmSpecies::resetTransMatrix(int option)
{
if (option == 1) { // asexual or simple sexual model
	transMatrix->ColCount = Ns+1; transMatrix->RowCount = Ns+1;
	transMatrix->Cells[1][0] = "juv"; transMatrix->Cells[0][1] = "juv";
	for (int j = 1; j < Ns+1; j++) transMatrix->Cells[1][j] = "0.0";
	for (int i = 2; i < Ns+1; i++) {
		transMatrix->Cells[i][0] = i-1;
		transMatrix->Cells[0][i] = i-1;
		for (int j = 1; j < Ns+1; j++) transMatrix->Cells[i][j] = "0.0";
	}
	MinAges->RowCount = Ns;
	for (int i = 1; i < Ns; i++) {
		MinAges->Cells[0][i] = i;
		MinAges->Cells[1][i] = 0;
	}
}
if (option == 2) { // complex sexual model
	transMatrix->ColCount = Ns*2+1; transMatrix->RowCount = Ns*2;
	transMatrix->Cells[0][1] = "juv";
	transMatrix->Cells[1][0] = "juv m";  transMatrix->Cells[2][0] = "juv f";
	transMatrix->Cells[1][1] = "0.0";  transMatrix->Cells[2][1] = "0.0";
	int ii = 1;
	for (int i = 3; i < Ns*2+1; i++) {
		if (i%2 != 0) transMatrix->Cells[i][0] = IntToStr(ii)+" m";
		else{
			transMatrix->Cells[i][0] = IntToStr(ii)+" f";
			ii++;
		}
	}
	ii = 1;
	for (int i = 2; i < Ns*2+1; i++) {
		if (i%2 == 0) transMatrix->Cells[0][i] = IntToStr(ii)+" m";
		else{
			transMatrix->Cells[0][i] = IntToStr(ii)+" f";
			ii++;
		}
	}
	for (int i = 1; i < Ns*2+1; i++){
		for (int j = 1; j < Ns*2; j++) transMatrix->Cells[i][j] = "0.0";
	}
	MinAges->RowCount = (Ns-1)*2+1;
	ii = 1;
	for (int i = 1; i < (Ns-1)*2+1; i++) {
		if (i%2 != 0) MinAges->Cells[0][i] = IntToStr(ii)+" m";
		else{
			MinAges->Cells[0][i] = IntToStr(ii)+" f";
			ii++;
		}
		MinAges->Cells[1][i] = 0;
	}
}
}

void __fastcall TfrmSpecies::resetSexEmigPar(int option)
{
if (option == 1) { // asexual or simple sexual model
//	if (SexEmigPar->RowCount == (Ns+1) && SexEmigPar->Cells[0][1] == "juv") {
//		//  format does not need to be changed
//		MessageDlg("RETURNING...",mtWarning,TMsgDlgButtons() << mbOK,0);
//		return;
//	}
//	else {
//		MessageDlg("NOT RETURNING...",mtWarning,TMsgDlgButtons() << mbOK,0);
//	}
	SexEmigPar->RowCount = Ns+1;
	SexEmigPar->Cells[0][0] = "Stage"; SexEmigPar->Cells[0][1] = "juv";
	SexEmigPar->Cells[1][1] = "0.0";
	for (int i = 2; i < Ns+1; i++){
		SexEmigPar->Cells[0][i] = i-1;
		for (int j = 1; j < SexEmigPar->ColCount; j++)
			SexEmigPar->Cells[j][i] = "0.0";
	}
}
if (option == 2) { // complex sexual model
//	if (SexEmigPar->RowCount == (Ns*2+1) && SexEmigPar->Cells[0][1] == "juv f") {
//		//  format does not need to be changed
//		return;
//	}
	SexEmigPar->RowCount = Ns*2+1;
	SexEmigPar->Cells[0][0] = "Stage/Sex";
	SexEmigPar->Cells[0][1] = "juv f";  SexEmigPar->Cells[0][2] = "juv m";
	for (int j = 1; j < SexEmigPar->ColCount; j++){
		SexEmigPar->Cells[j][1] = "0.0";    SexEmigPar->Cells[j][2] = "0.0";
	}
	for (int i = 1; i < Ns; i++) {
		for (int j = 0; j < NSEXES; j++) {
			if (j == 0) SexEmigPar->Cells[0][2*i+1+j] = IntToStr(i) + " f";
			else SexEmigPar->Cells[0][2*i+1+j] = IntToStr(i) + " m";
			for (int k = 1; k < SexEmigPar->ColCount; k++)
				SexEmigPar->Cells[k][2*i+1+j] = "0.0";
		}
	}
}
}

void __fastcall TfrmSpecies::resetSexEmigDensPar(int option)
{
if (option == 1) { // asexual or simple sexual model
	SexEmigDensPar->RowCount = Ns+1;
	SexEmigDensPar->Cells[0][0] = "Stage";
	for (int i = 1; i < Ns+1; i++) {
		if (i == 1) SexEmigDensPar->Cells[0][i] = "juv";
		else SexEmigDensPar->Cells[0][i] = i-1;
		for (int j = 1; j < SexEmigDensPar->ColCount; j++)
			SexEmigDensPar->Cells[j][i] = "0.0";
	}
}
if (option == 2) { // complex sexual model
	SexEmigDensPar->RowCount = Ns*2+1;
	SexEmigDensPar->Cells[0][0] = "Stage/Sex";
	SexEmigDensPar->Cells[0][1] = "juv f";  SexEmigDensPar->Cells[0][2] = "juv m";
	for (int j = 1; j < SexEmigDensPar->ColCount; j++){
		SexEmigDensPar->Cells[j][1] = "0.0";    SexEmigDensPar->Cells[j][2] = "0.0";
	}
	for (int i = 1; i < Ns; i++) {
		for (int j = 0; j < NSEXES; j++) {
			if (j == 0) SexEmigDensPar->Cells[0][2*i+1+j] = IntToStr(i) + " f";
			else SexEmigDensPar->Cells[0][2*i+1+j] = IntToStr(i) + " m";
			for (int k = 1; k < SexEmigDensPar->ColCount; k++)
				SexEmigDensPar->Cells[k][2*i+1+j] = "0.0";
		}
	}
}

}

void __fastcall TfrmSpecies::resetSexTransferPar(int option)
{
if (option == 1) { // asexual or simple sexual model
	SexKernPar->RowCount = Ns + 1;
	SexKernPar->Cells[0][0] = "Stage";
	for (int i = 1; i < Ns+1; i++){
		if (i == 1) SexKernPar->Cells[0][1] = "juv";
		else SexKernPar->Cells[0][i] = i-1;
		for (int j = 1; j < SexKernPar->ColCount; j++)
			SexKernPar->Cells[j][i] = "0.0";
	}
}
if (option == 2) { // complex sexual model
	SexKernPar->RowCount = Ns*2+1;
	SexKernPar->Cells[0][0] = "Stage/Sex";
	SexKernPar->Cells[0][1] = "juv f";  SexKernPar->Cells[0][2] = "juv m";
	for (int j = 1; j < SexKernPar->ColCount; j++){
		SexKernPar->Cells[j][1] = "0.0";  SexKernPar->Cells[j][2] = "0.0";
	}
	for (int i = 1; i < Ns; i++) {
		for (int j = 0; j < NSEXES; j++) {
			if (j == 0) SexKernPar->Cells[0][2*i+1+j] = IntToStr(i) + " f";
			else SexKernPar->Cells[0][2*i+1+j] = IntToStr(i) + " m";
			for (int k = 1; k < SexKernPar->ColCount; k++)
				SexKernPar->Cells[k][2*i+1+j] = "0.0";
		}
	}
}
}

void __fastcall TfrmSpecies::resetSexSettlePar(int option)
{
if (option == 1) { // dispersal kernel - asexual or simple sexual model
	SexSettlePar->RowCount = Ns + 1;
	SexSettlePar->Cells[0][0] = "Stage"; SexSettlePar->Cells[0][1] = "juv";
	for (int i = 1; i < Ns+1; i++){
		if (i > 1) SexSettlePar->Cells[0][i] = i-1;
		SexSettlePar->Cells[1][i] = "0"; SexSettlePar->Cells[2][i] = "0";
	}
}

if (option == 2) { // dispersal kernel - complex sexual model
	SexSettlePar->RowCount = Ns*2 + 1;
	SexSettlePar->Cells[0][0] = "Stage/Sex";
	SexSettlePar->Cells[0][1] = "juv f"; SexSettlePar->Cells[0][2] = "juv m";
	for (int i = 1; i < Ns; i++) {
		for (int j = 0; j < NSEXES; j++) {
			if (j == 0) SexSettlePar->Cells[0][2*i+1+j] = IntToStr(i) + " f";
			else SexSettlePar->Cells[0][2*i+1+j] = IntToStr(i) + " m";
			for (int k = 1; k < SexSettlePar->ColCount; k++)
				SexSettlePar->Cells[k][2*i+1+j] = "0";
		}
	}
}

if (option == 3) { // movement process - asexual or simple sexual model
	SexSettlePar->RowCount = Ns + 1;
	SexSettlePar->Cells[0][0] = "Stage"; SexSettlePar->Cells[0][1] = "juv";
	for (int i = 1; i < Ns+1; i++){
		if (i > 1) SexSettlePar->Cells[0][i] = i-1;
		SexSettlePar->Cells[1][i] = "0";
		SexSettlePar->Cells[2][i] = "0";
		SexSettlePar->Cells[3][i] = "1.0";
		SexSettlePar->Cells[4][i] = "0.0";
		SexSettlePar->Cells[5][i] = "1.0";
	}
}

if (option == 4) { // // movement process - complex sexual model
	SexSettlePar->RowCount = Ns*2 + 1;
	SexSettlePar->Cells[0][0] = "Stage/Sex";
	SexSettlePar->Cells[0][1] = "juv f"; SexSettlePar->Cells[0][2] = "juv m";
	int ii = 1;
	for (int i = 1; i < Ns*2+1; i++) {
		if (i > 2) {
			if (i%2 != 0) SexSettlePar->Cells[0][i] = IntToStr(ii)+" m";
			else{
				SexSettlePar->Cells[0][i] = IntToStr(ii)+" f";
				ii++;
			}
		}
		SexSettlePar->Cells[1][i] = "0";
		SexSettlePar->Cells[2][i] = "0";
		SexSettlePar->Cells[3][i] = "1.0";
		SexSettlePar->Cells[4][i] = "0.0";
		SexSettlePar->Cells[5][i] = "1.0";
	}
	for (int i = 1; i < Ns; i++) {
		for (int j = 0; j < NSEXES; j++) {
			if (j == 0) SexSettlePar->Cells[0][2*i+1+j] = IntToStr(i) + " f";
			else SexSettlePar->Cells[0][2*i+1+j] = IntToStr(i) + " m";
			SexSettlePar->Cells[1][2*i+1+j] = "0";
			SexSettlePar->Cells[2][2*i+1+j] = "0";
			SexSettlePar->Cells[3][2*i+1+j] = "1.0";
			SexSettlePar->Cells[4][2*i+1+j] = "0.0";
			SexSettlePar->Cells[5][2*i+1+j] = "1.0";
		}
	}
}
}

void resetFecundityDens(TfrmDensity* frmDensity, int option)
{
if (option < 2) { // asexual or simple sexual model
	if (frmDensity->FecundityDens->ColCount == (Ns+1)
	&&  frmDensity->FecundityDens->Cells[1][0] == "juv") {
		// matrix format does not need to be changed
		return;
	}
	frmDensity->FecundityDens->ColCount = Ns+1;
	frmDensity->FecundityDens->RowCount = Ns+1;
	frmDensity->FecundityDens->Cells[1][0] = "juv";
	frmDensity->FecundityDens->Cells[0][1] = "juv";
	for (int j = 1; j < Ns+1; j++) frmDensity->FecundityDens->Cells[1][j] = "0.0";
	for (int i = 2; i < Ns+1; i++){
		frmDensity->FecundityDens->Cells[i][0] = i-1;
		frmDensity->FecundityDens->Cells[0][i] = i-1;
		frmDensity->FecundityDens->Cells[i][1] = "0.0";
		for (int j = 2; j < Ns+1; j++) frmDensity->FecundityDens->Cells[i][j] = "1.0";
	}
}
if (option == 2) { // complex sexual model
	if (frmDensity->FecundityDens->ColCount == (Ns*2+1)
	&&  frmDensity->FecundityDens->Cells[1][0] == "juv m") {
		// matrix format does not need to be changed
		return;
	}
	frmDensity->FecundityDens->ColCount = Ns*2+1;
	frmDensity->FecundityDens->RowCount = Ns*2+1;
	frmDensity->FecundityDens->Cells[1][0] = "juv m";
	frmDensity->FecundityDens->Cells[0][1] = "juv m";
	frmDensity->FecundityDens->Cells[2][0] = "juv f";
	frmDensity->FecundityDens->Cells[0][2] = "juv f";
	int ii = 1;
	for (int i = 3; i < Ns*2+1; i++) {
		if (i%2 != 0) {
			frmDensity->FecundityDens->Cells[i][0] = IntToStr(ii)+" m";
			frmDensity->FecundityDens->Cells[0][i] = IntToStr(ii)+" m";
		}
		else{
			frmDensity->FecundityDens->Cells[i][0] = IntToStr(ii)+" f";
			frmDensity->FecundityDens->Cells[0][i] = IntToStr(ii)+" f";
			ii++;
		}
	}
	for (int i = 1; i < Ns*2+1; i++) {
		for (int j = 1; j < Ns*2+1; j++) {
			if (i <= 2 || j <= 2) frmDensity->FecundityDens->Cells[i][j] = "0.0";
			else frmDensity->FecundityDens->Cells[i][j] = "1.0";
		}
	}
}
}

void resetDevelopDens(TfrmDensity* frmDensity, int option)
{
if (option < 2) { // asexual or simple sexual model
	if (frmDensity->DevelopDens->ColCount == (Ns+1)
	&&  frmDensity->DevelopDens->Cells[1][0] == "juv") {
		// matrix format does not need to be changed
		return;
	}
	frmDensity->DevelopDens->ColCount = Ns+1;
	frmDensity->DevelopDens->RowCount = Ns+1;
	frmDensity->DevelopDens->Cells[1][0] = "juv";
	frmDensity->DevelopDens->Cells[0][1] = "juv";
	for (int j = 1; j < Ns+1; j++) frmDensity->DevelopDens->Cells[1][j] = "1.0";
	for (int i = 2; i < Ns+1; i++){
		frmDensity->DevelopDens->Cells[i][0] = i-1;
		frmDensity->DevelopDens->Cells[0][i] = i-1;
		for (int j = 1; j < Ns+1; j++) {
			if (i == Ns) frmDensity->DevelopDens->Cells[i][j] = "0.0";
			else frmDensity->DevelopDens->Cells[i][j] = "1.0";
		}
	}
}
if (option == 2) { // complex sexual model
	if (frmDensity->DevelopDens->ColCount == (Ns*2+1)
	&&  frmDensity->DevelopDens->Cells[1][0] == "juv m") {
		// matrix format does not need to be changed
		return;
	}
	frmDensity->DevelopDens->ColCount = Ns*2+1;
	frmDensity->DevelopDens->RowCount = Ns*2+1;
	frmDensity->DevelopDens->Cells[1][0] = "juv m";
	frmDensity->DevelopDens->Cells[0][1] = "juv m";
	frmDensity->DevelopDens->Cells[2][0] = "juv f";
	frmDensity->DevelopDens->Cells[0][2] = "juv f";
	int ii = 1;
	for (int i = 3; i < Ns*2+1; i++) {
		if (i%2 != 0) {
			frmDensity->DevelopDens->Cells[i][0] = IntToStr(ii)+" m";
			frmDensity->DevelopDens->Cells[0][i] = IntToStr(ii)+" m";
		}
		else{
			frmDensity->DevelopDens->Cells[i][0] = IntToStr(ii)+" f";
			frmDensity->DevelopDens->Cells[0][i] = IntToStr(ii)+" f";
			ii++;
		}
	}
	for (int i = 1; i < Ns*2+1; i++){
		for (int j = 1; j < Ns*2+1; j++) {
			if (i >= (Ns*2-1)) frmDensity->DevelopDens->Cells[i][j] = "0.0";
			else frmDensity->DevelopDens->Cells[i][j] = "1.0";
		}
	}

}
}

void resetSurvivalDens(TfrmDensity* frmDensity, int option)
{
//UnicodeString cellval;
//cellval = frmDensity->SurvivalDens->Cells[1][0];
//MessageDlg(cellval.c_str(),mtWarning,TMsgDlgButtons() << mbOK,0);

if (option < 2) { // asexual or simple sexual model
	if (frmDensity->SurvivalDens->ColCount == (Ns+1)
	&&  frmDensity->SurvivalDens->Cells[1][0] == "juv") {
		// matrix format does not need to be changed
		return;
	}
	frmDensity->SurvivalDens->ColCount = Ns+1;
	frmDensity->SurvivalDens->RowCount = Ns+1;
	frmDensity->SurvivalDens->Cells[1][0] = "juv";
	frmDensity->SurvivalDens->Cells[0][1] = "juv";
	for (int j = 1; j < Ns+1; j++) frmDensity->SurvivalDens->Cells[1][j] = "1.0";
	for (int i = 2; i < Ns+1; i++){
		frmDensity->SurvivalDens->Cells[i][0] = i-1;
		frmDensity->SurvivalDens->Cells[0][i] = i-1;
		for (int j = 1; j < Ns+1; j++) frmDensity->SurvivalDens->Cells[i][j] = "1.0";
	}
}
if (option == 2) { // complex sexual model
	if (frmDensity->SurvivalDens->ColCount == (Ns*2+1)
	&&  frmDensity->SurvivalDens->Cells[1][0] == "juv m") {
		// matrix format does not need to be changed
		return;
	}
	frmDensity->SurvivalDens->ColCount = Ns*2+1;
	frmDensity->SurvivalDens->RowCount = Ns*2+1;
	frmDensity->SurvivalDens->Cells[1][0] = "juv m";
	frmDensity->SurvivalDens->Cells[0][1] = "juv m";
	frmDensity->SurvivalDens->Cells[2][0] = "juv f";
	frmDensity->SurvivalDens->Cells[0][2] = "juv f";
	int ii = 1;
	for (int i = 3; i < Ns*2+1; i++) {
		if (i%2 != 0) {
			frmDensity->SurvivalDens->Cells[i][0] = IntToStr(ii)+" m";
			frmDensity->SurvivalDens->Cells[0][i] = IntToStr(ii)+" m";
		}
		else{
			frmDensity->SurvivalDens->Cells[i][0] = IntToStr(ii)+" f";
			frmDensity->SurvivalDens->Cells[0][i] = IntToStr(ii)+" f";
			ii++;
		}
	}
	for (int i = 1; i < Ns*2+1; i++){
		for (int j = 1; j < Ns*2+1; j++) frmDensity->SurvivalDens->Cells[i][j] = "1.0";
	}
}
}

//---------------------------------------------------------------------------
void __fastcall TfrmSpecies::TSDispersal2Enter(TObject *Sender)
{
//MessageDlg("EXECUTING TSDispersal2Enter()",mtWarning,TMsgDlgButtons() << mbOK,0);
if (RGMovements->ItemIndex == 1) { // movement process
	PanelSexTransfer->Visible = false;
}
}

//---------------------------------------------------------------------------
void __fastcall TfrmSpecies::TSDispersalExit(TObject *Sender)
{
//if (formCancelled) {
//	MessageDlg("EXECUTING TSDispersalexit() after cancel",mtWarning,TMsgDlgButtons() << mbOK,0);
//}
//else {
//	MessageDlg("EXECUTING TSDispersalexit() normally",mtWarning,TMsgDlgButtons() << mbOK,0);
//}
}

//---------------------------------------------------------------------------

// Demographic functions

void __fastcall TfrmSpecies::CBStageModelClick(TObject *Sender)
{
landParams ppLand = pLandscape->getLandParams();
demogrParams dem = pSpecies->getDemogr();

Ns = StrToInt(edtNstages->Text);
if (Ns < 2) Ns = 2;
forceInit = true;

//CBSexEmig->Checked = false;
//CBSexKernels->Checked = false;
//CBSexSettMovt->Checked = false;
//CBSexSettKern->Checked = false;
CBStageEmig->Checked = false;
CBStageKernels->Checked = false;
CBStageSettMovt->Checked = false;
CBStageSettKern->Checked = false;

if (CBStageModel->Checked) {
  PanelStage->Enabled = true; PanelStage->Visible = true;
  edtR->Enabled = false; edtR->Visible = false;
  edtC->Enabled = false; edtC->Visible = false;
  edtPRep->Enabled = true; edtPRep->Visible = true;
	edtRepInterval->Enabled = true; edtRepInterval->Visible = true;
	if (ppLand.dynamic) {
		RGaction->Enabled = true; RGaction->Visible = true;
	}
	else {
		RGaction->Enabled = false; RGaction->Visible = false;
	}
	if (CBIndVarEmig->Checked) {
		CBStageEmig->Enabled = false; CBStageEmig->Visible = true;
	}
	else {
		CBStageEmig->Enabled = true; CBStageEmig->Visible = true;
  }
	CBStageKernels->Enabled = true; CBStageKernels->Visible = true;
	CBStageSettMovt->Enabled = true; CBStageSettMovt->Visible = true;
	CBStageSettKern->Enabled = true; CBStageSettKern->Visible = true;

	if (CBStageEmig->Checked) {
		CBIndVarEmig->Enabled = false; CBIndVarEmig->Checked = false;
	}
	else {
		CBIndVarEmig->Enabled = true;
	}
	if (CBStageKernels->Checked) {
		CBIndVarKernel->Enabled = false; CBIndVarKernel->Checked = false;
	}
	else {
		CBIndVarKernel->Enabled = true;
	}
	if (CBStageSettMovt->Checked) {
		CBIndVarSettle->Enabled = false; CBIndVarSettle->Checked = false;
	}
	else {
		CBIndVarSettle->Enabled = true;
  }

  MinAges->Cells[0][0] = "Stage";
  MinAges->Cells[1][0] = "Age";

	if (RGReproduction->ItemIndex != 2) {
		resetTransMatrix(1);
	}
  // explicit 2 sexes model
  else{
		resetTransMatrix(2);
	} // end of explicit 2 sexes model

	edtMaxStepYear->Visible = true; edtMaxStepYear->Enabled = true;
  Memo2->Visible = true;

	edtK->EditLabel->Caption = "1/b: ";
	edtKLabel->Caption = "Strength of density dependence ";
	SGhabLabel->Caption = "Habitat-specific strength\nof density dependence (1/b): ";
	SGhabLabel->Top = 30; SGhabLabel->Width = 260; SGhabLabel->WordWrap = false;
	SGhab->Cells[1][0] = "1/b";

	dem.stageStruct = true;
	pSpecies->setDemogr(dem);
} // end of CBStageModel->Checked
else { // !CBStageModel->Checked
  PanelStage->Enabled = false; PanelStage->Visible = false;
  edtR->Enabled = true; edtR->Visible = true;
  edtC->Enabled = true; edtC->Visible = true;
  edtPRep->Enabled = false; edtPRep->Visible = false;
  edtRepInterval->Enabled = false; edtRepInterval->Visible = false;
	RGaction->Enabled = false; RGaction->Visible = false;
	edtMaxStepYear->Visible = false; edtMaxStepYear->Enabled = false;
  Memo2->Visible = false;
	CBStageEmig->Enabled = false;     CBStageEmig->Visible = false;
	edtEmigStage->Enabled = false;    edtEmigStage->Visible = false;
	CBStageKernels->Enabled = false;  CBStageKernels->Visible = false;
	CBStageSettMovt->Enabled = false; CBStageSettMovt->Visible = false;
	CBStageSettKern->Enabled = false; CBStageSettKern->Visible = false;
	CBIndVarEmig->Enabled = true;
	CBIndVarKernel->Enabled = true;
	CBIndVarSettle->Enabled = true;

	edtK->EditLabel->Caption = "K: ";
	edtKLabel->Caption = "inds / ha (assuming 100% quality) ";
	SGhabLabel->Caption = "Habitat-dependent carrying capacities: ";
	SGhabLabel->Top = 50; SGhabLabel->Width = 260; SGhabLabel->WordWrap = true;
	SGhab->Cells[1][0] = "K (inds/ha)";

	dem.stageStruct = false;
	pSpecies->setDemogr(dem);
}
}
//---------------------------------------------------------------------------
void __fastcall TfrmSpecies::RGReproductionClick(TObject *Sender)
{
initParams init = paramsInit->getInit();

Ns = StrToInt(edtNstages->Text);
if (Ns < 2) Ns = 2;

CBSexEmig->Checked = false;
//CBStageEmig->Checked = false;
CBSexKernels->Checked = false;
//CBStageKernels->Checked = false;
CBSexSettKern->Checked = false;
//CBStageSettKern->Checked = false;
CBSettKernMate->Checked = false;
CBSexSettMovt->Checked = false;
//CBStageSettMovt->Checked = false;

if (frmDensity->Showing) frmDensity->Close();

switch (RGReproduction->ItemIndex) {
case 0:
	edtHarem->Enabled = false;       edtHarem->Visible = false;
	edtSexRatio->Enabled = false;    edtSexRatio->Visible = false;
	CBSexEmig->Enabled = false;      CBSexEmig->Visible = false;
	CBSexKernels->Enabled = false;   CBSexKernels->Visible = false;
	CBSexSettKern->Enabled = false;  CBSexSettKern->Visible = false;
	CBSettKernMate->Enabled = false; CBSettKernMate->Visible = false;
	CBSettKernMate->Checked = false;
	CBSexSettMovt->Enabled = false;  CBSexSettMovt->Visible = false;
	CBFindMate->Enabled = false;     CBFindMate->Checked = false;
	if (CBStageModel->Checked) {
		edtR->Enabled = false; edtR->Visible = false;
		resetTransMatrix(1);
	}
	else {
		edtR->Enabled = true; edtR->Visible = true;
	}
	break;
case 1:
	edtHarem->Enabled = false;      edtHarem->Visible = false;
	edtSexRatio->Enabled = true;    edtSexRatio->Visible = true;
	CBSexEmig->Enabled = true;      CBSexEmig->Visible = true;
	CBSexKernels->Enabled = true;   CBSexKernels->Visible = true;
	CBSexSettKern->Enabled = true;  CBSexSettKern->Visible = true;
	CBSettKernMate->Enabled = true; CBSettKernMate->Visible = true;
	CBSexSettMovt->Enabled = true;  CBSexSettMovt->Visible = true;
	CBFindMate->Enabled = true;
	if (CBStageModel->Checked) {
		edtR->Enabled = false; edtR->Visible = false;
		resetTransMatrix(1);
		if (CBweightFec->Checked)  resetFecundityDens(frmDensity,1);
		if (CBweightDev->Checked)  resetDevelopDens(frmDensity,1);
		if (CBweightSurv->Checked) resetSurvivalDens(frmDensity,1);
	}
	else {
		edtR->Enabled = true; edtR->Visible = true;
	}
	break;
case 2:
	edtHarem->Enabled = true;       edtHarem->Visible = true;
	CBSettKernMate->Enabled = true; CBSettKernMate->Visible = true;
	edtSexRatio->Enabled = true;    edtSexRatio->Visible = true;
	CBSexEmig->Enabled = true;      CBSexEmig->Visible = true;
	CBSexKernels->Enabled = true;   CBSexKernels->Visible = true;
	CBSexSettKern->Enabled = true;  CBSexSettKern->Visible = true;
	CBSettKernMate->Enabled = true; CBSettKernMate->Visible = true;
	CBSexSettMovt->Enabled = true;  CBSexSettMovt->Visible = true;
	CBFindMate->Enabled = true;
	if (CBStageModel->Checked) {
		edtR->Enabled = false; edtR->Visible = false;
		resetTransMatrix(2);
		if (CBweightFec->Checked)  resetFecundityDens(frmDensity,2);
		if (CBweightDev->Checked)  resetDevelopDens(frmDensity,2);
		if (CBweightSurv->Checked) resetSurvivalDens(frmDensity,2);
	}
	else {
		edtR->Enabled = true; edtR->Visible = true;
	}
	break;
}

if (init.seedType == 2) { // from initial individuals file
	forceInit = true;
}

}
//---------------------------------------------------------------------------
void __fastcall TfrmSpecies::edtNstagesExit(TObject *Sender)
{
if (Ns != StrToInt(edtNstages->Text)) // no. of stages has been changed
	forceInit = true;
Ns = StrToInt(edtNstages->Text);
if (Ns < 2) {
	Ns = 2; edtNstages->Text = "2";
	MessageDlg("The minimum number of stages is 2",mtError,
				TMsgDlgButtons() << mbRetry,0);
}
if (Ns > 10) {
	Ns = 10; edtNstages->Text = "10";
	MessageDlg("The maximum number of stages is 10",mtError,
				TMsgDlgButtons() << mbRetry,0);
}

if (frmDensity->Showing) frmDensity->Close();

if (RGReproduction->ItemIndex == 2) { //explicit 2 sexes model
	resetTransMatrix(2);
	if (CBweightFec->Checked)  resetFecundityDens(frmDensity,2);
	if (CBweightDev->Checked)  resetDevelopDens(frmDensity,2);
	if (CBweightSurv->Checked) resetSurvivalDens(frmDensity,2);
}
else {
	resetTransMatrix(1);
	if (CBweightFec->Checked)  resetFecundityDens(frmDensity,1);
	if (CBweightDev->Checked)  resetDevelopDens(frmDensity,1);
	if (CBweightSurv->Checked) resetSurvivalDens(frmDensity,1);
}
CBStageEmig->Checked = false;
CBStageKernels->Checked = false;
CBStageSettKern->Checked = false;
CBStageSettMovt->Checked = false;
}
//---------------------------------------------------------------------------
void __fastcall TfrmSpecies::MinAgesExit(TObject *Sender)
{
int Nrows,i;
bool error = false;
if (RGReproduction->ItemIndex == 2) Nrows = 2; else Nrows = 1;

for (i = 0; i < Nrows; i++) {
	if (StrToInt(MinAges->Cells[1][i+1]) != 0) error = true;
}
if (error) {
	MessageDlg("The minimum age for the first stage must be zero.",mtError,
				TMsgDlgButtons() << mbOK,0);
	for (i = 0; i < Nrows; i++) MinAges->Cells[1][i+1] = 0;
}

}
//---------------------------------------------------------------------------
void __fastcall TfrmSpecies::CBweightFecClick(TObject *Sender)
{
Ns = StrToInt(edtNstages->Text);
if (Ns < 2) Ns = 2;

if (CBweightFec->Checked) {
	if (CBFecundity->Checked) {
		BtnSetCoeff->Enabled = true;
		//set coefficient matrices------------------------------------
		//Fecundity
		frmDensity->FecundityDens->Enabled = true;
		frmDensity->FecundityDens->Visible = true;
		if (RGReproduction->ItemIndex != 2) {
			resetFecundityDens(frmDensity,1);
		}
		else {
			resetFecundityDens(frmDensity,2);
		}
	}
	else {
		MessageDlg("Density dependence in fecundity has not been selected.",mtError,
				TMsgDlgButtons() << mbOK,0);
		CBweightFec->Checked = false;
	}
}
else {
	if (CBweightSurv == false && CBweightDev == false) BtnSetCoeff->Enabled = false;
	frmDensity->FecundityDens->Enabled = false;
	frmDensity->FecundityDens->Visible = false;
}
}
//---------------------------------------------------------------------------
void __fastcall TfrmSpecies::CBweightDevClick(TObject *Sender)
{
Ns = StrToInt(edtNstages->Text);
if (Ns < 2) Ns = 2;

if (CBweightDev->Checked) {
	if (CBDevelopment->Checked) {
		BtnSetCoeff->Enabled = true;
		//set coefficient matrices------------------------------------
		//Development
		frmDensity->DevelopDens->Enabled = true;
		frmDensity->DevelopDens->Visible = true;
		if (RGReproduction->ItemIndex != 2) {
			resetDevelopDens(frmDensity,1);
		}
		else {
			resetDevelopDens(frmDensity,2);
		}
	}
	else{
		MessageDlg("Density dependence in development has not been selected.",mtError,
				TMsgDlgButtons() << mbOK,0);
		CBweightDev->Checked = false;
	}
}
else {
	if (CBweightSurv == false && CBweightFec == false) BtnSetCoeff->Enabled = false;
	frmDensity->DevelopDens->Enabled = false;
	frmDensity->DevelopDens->Visible = false;
}
}
//---------------------------------------------------------------------------
void __fastcall TfrmSpecies::CBweightSurvClick(TObject *Sender)
{
Ns = StrToInt(edtNstages->Text);
if (Ns < 2) Ns = 2;

if (CBweightSurv->Checked) {
	if (CBSurvival->Checked) {
		BtnSetCoeff->Enabled = true;
		//set coefficient matrices------------------------------------
		//Survival
		if (CBSurvival->Checked) {
			frmDensity->SurvivalDens->Enabled = true;
			frmDensity->SurvivalDens->Visible = true;
			if (RGReproduction->ItemIndex != 2) {
				resetSurvivalDens(frmDensity,1);
			}
			else {
				resetSurvivalDens(frmDensity,2);
			}
		}
	}
	else {
		MessageDlg("Density dependence in survival has not been selected.",mtError,
				TMsgDlgButtons() << mbOK,0);
		CBweightSurv->Checked = false;
	}
}
else {
	if (CBweightDev == false && CBweightFec == false) BtnSetCoeff->Enabled = false;
	frmDensity->SurvivalDens->Enabled = false;
	frmDensity->SurvivalDens->Visible = false;
}
}
//---------------------------------------------------------------------------
void __fastcall TfrmSpecies::CBFecundityClick(TObject *Sender)
{
Ns = StrToInt(edtNstages->Text);
if (Ns < 2) Ns = 2;

if (CBFecundity->Checked) {
	frmDensity->FecundityDens->Enabled = true;
	frmDensity->FecundityDens->Visible = true;
	resetFecundityDens(frmDensity,RGReproduction->ItemIndex);
}
else {
	frmDensity->FecundityDens->Enabled = false;
	frmDensity->FecundityDens->Visible = false;
	if (CBweightFec->Checked) CBweightFec->Checked = false;
	if (CBweightDev->Checked == false && CBweightSurv->Checked == false) {
		if (frmDensity->Showing) frmDensity->Close();
		BtnSetCoeff->Enabled = false;
	}
}
}
//---------------------------------------------------------------------------
void __fastcall TfrmSpecies::CBDevelopmentClick(TObject *Sender)
{
Ns = StrToInt(edtNstages->Text);
if (Ns < 2) Ns = 2;

if (CBDevelopment->Checked) {
	edtDevCoeff->Enabled = true; edtDevCoeff->Visible = true;
	frmDensity->DevelopDens->Enabled = true;
	frmDensity->DevelopDens->Visible = true;
	resetDevelopDens(frmDensity,RGReproduction->ItemIndex);
}
else {
	edtDevCoeff->Enabled = false; edtDevCoeff->Visible = false;
	frmDensity->DevelopDens->Enabled = false;
	frmDensity->DevelopDens->Visible = false;
	if (CBweightDev->Checked) CBweightDev->Checked = false;
	if (CBweightFec->Checked == false && CBweightSurv->Checked == false) {
		if (frmDensity->Showing) frmDensity->Close();
		BtnSetCoeff->Enabled = false;
	}
}
}
//---------------------------------------------------------------------------
void __fastcall TfrmSpecies::CBSurvivalClick(TObject *Sender)
{
Ns = StrToInt(edtNstages->Text);
if (Ns < 2) Ns = 2;

if (CBSurvival->Checked) {
	edtSurvCoeff->Enabled = true; edtSurvCoeff->Visible = true;
	frmDensity->SurvivalDens->Enabled = true;
	frmDensity->SurvivalDens->Visible = true;
	resetSurvivalDens(frmDensity,RGReproduction->ItemIndex);
}
else {
	edtSurvCoeff->Enabled = false; edtSurvCoeff->Visible = false;
	frmDensity->SurvivalDens->Enabled = false;
	frmDensity->SurvivalDens->Visible = false;
	if (CBweightSurv->Checked) CBweightSurv->Checked = false;
	if (CBweightDev->Checked == false && CBweightFec->Checked == false) {
		if (frmDensity->Showing) frmDensity->Close();
		BtnSetCoeff->Enabled = false;
	}
}
}
//---------------------------------------------------------------------------
void __fastcall TfrmSpecies::BtnSetCoeffClick(TObject *Sender)
{
Ns = StrToInt(edtNstages->Text);
if (Ns < 2) Ns = 2;

if (CBweightFec->Checked) {
	frmDensity->FecundityDens->Enabled = true;
	frmDensity->FecundityDens->Visible = true;
	if (RGReproduction->ItemIndex != 2) {
		resetFecundityDens(frmDensity,1);
	}
	else{
		resetFecundityDens(frmDensity,2);
	}
}
else {
  frmDensity->FecundityDens->Enabled = false;
  frmDensity->FecundityDens->Visible = false;
}

if (CBweightDev->Checked) {
	frmDensity->Height = 302;
	if (RGReproduction->ItemIndex != 2) {
		resetDevelopDens(frmDensity,1);
	}
	else{
		resetDevelopDens(frmDensity,2);
	}
}
else{
	frmDensity->Height = 302;
	frmDensity->DevelopDens->Enabled = false;
	frmDensity->DevelopDens->Visible = false;
}

if (CBweightSurv->Checked) {
	frmDensity->SurvivalDens->Enabled = true;
	frmDensity->SurvivalDens->Visible = true;
	if (RGReproduction->ItemIndex != 2) {
		resetSurvivalDens(frmDensity,1);
	}
	else{
		resetSurvivalDens(frmDensity,2);
	}
}
else{
	frmDensity->SurvivalDens->Enabled = false;
	frmDensity->SurvivalDens->Visible = false;
}

 frmDensity->Show();
}
//---------------------------------------------------------------------------

// Emigration functions

void __fastcall TfrmSpecies::RGEmigProbClick(TObject *Sender)
{
setEmigFields();
}

//---------------------------------------------------------------------------
void __fastcall TfrmSpecies::CBSexEmigClick(TObject *Sender)
{

setEmigFields();

}
//---------------------------------------------------------------------------
void __fastcall TfrmSpecies::CBStageEmigClick(TObject *Sender)
{

//if (CBStageEmig->Checked && CBIndVarEP->Checked) {
//	// stage-dependent parameters and individual variability are incompatible
//	CBIndVarEP->Checked = false;
//	MessageDlg("Stage dependent dispersal parameters may not vary between individuals",
//		mtWarning, TMsgDlgButtons() << mbRetry,0);
//	return;
//}
if (CBStageEmig->Checked) {
	CBIndVarEmig->Enabled = false; CBIndVarEmig->Checked = false;
	edtEmigStage->Enabled = false; edtEmigStage->Visible = false;
}
else {
	CBIndVarEmig->Enabled = true;
}

setEmigFields();

}

void __fastcall TfrmSpecies::setEmigFields(void) {

//MessageDlg("EXECUTING setEmigFields()",mtWarning,TMsgDlgButtons() << mbOK,0);

if (RGEmigProb->ItemIndex == 0) { // density-independent
	if (RGMovements->ItemIndex == 0) { // dispersel kernels
		CBFullKernel->Enabled = true;
	}
	PanelEP->Visible = true; PanelEP->Enabled = true;
	if (CBIndVarEmig->Checked) {
		edtEP->Visible = true; edtEP->Enabled = false;
		edtEPmean->Visible = true; edtEPmean->Enabled = true;
		edtEPsd->Visible = true; edtEPsd->Enabled = true;
		edtEPscale->Visible = true; edtEPscale->Enabled = true;
	}
	else {
		edtEP->Visible = true; edtEP->Enabled = true;
		edtEPmean->Visible = false; edtEPmean->Enabled = false;
		edtEPsd->Visible = false; edtEPsd->Enabled = false;
		edtEPscale->Visible = false; edtEPscale->Enabled = false;
	}
	PanelDensEmig->Visible = false; PanelDensEmig->Enabled = false;
}
else { // density-dependent
	CBFullKernel->Enabled = false; CBFullKernel->Checked = false;
	PanelEP->Visible = false; PanelEP->Enabled = false;
	if (CBIndVarEmig->Checked) {
		edtD0->Visible = true; edtD0->Enabled = false;
		edtD0Mean->Visible = true; edtD0Mean->Enabled = true;
		edtD0SD->Visible = true; edtD0SD->Enabled = true;
		edtD0Scale->Visible = true; edtD0Scale->Enabled = true;
		edtAlpha->Visible = true; edtAlpha->Enabled = false;
		edtAlphaMean->Visible = true; edtAlphaMean->Enabled = true;
		edtAlphaSD->Visible = true; edtAlphaSD->Enabled = true;
		edtAlphaScale->Visible = true; edtAlphaScale->Enabled = true;
		edtBeta->Visible = true; edtBeta->Enabled = false;
		edtBetaMean->Visible = true; edtBetaMean->Enabled = true;
		edtBetaSD->Visible = true; edtBetaSD->Enabled = true;
		edtBetaScale->Visible = true; edtBetaScale->Enabled = true;
	}
	else {
		edtD0->Visible = true; edtD0->Enabled = true;
		edtD0Mean->Visible = false; edtD0Mean->Enabled = false;
		edtD0SD->Visible = false; edtD0SD->Enabled = false;
		edtD0Scale->Visible = false; edtD0Scale->Enabled = false;
		edtAlpha->Visible = true; edtAlpha->Enabled = true;
		edtAlphaMean->Visible = false; edtAlphaMean->Enabled = false;
		edtAlphaSD->Visible = false; edtAlphaSD->Enabled = false;
		edtAlphaScale->Visible = false; edtAlphaScale->Enabled = false;
		edtBeta->Visible = true; edtBeta->Enabled = true;
		edtBetaMean->Visible = false; edtBetaMean->Enabled = false;
		edtBetaSD->Visible = false; edtBetaSD->Enabled = false;
		edtBetaScale->Visible = false; edtBetaScale->Enabled = false;
	}
	PanelDensEmig->Visible = true; PanelDensEmig->Enabled = true;
}

if (CBSexEmig->Checked || CBStageEmig->Checked) {
	PanelSexEmig->Visible = true; PanelSexEmig->Enabled = true;

	PanelEP->Visible = false; PanelDensEmig->Visible = false;
	if (RGEmigProb->ItemIndex == 0) { // density-independent
		SexEmigPar->Enabled = true; SexEmigPar->Visible = true;
		LabelSexEmigPar->Visible = true;
		SexEmigDensPar->Enabled = false; SexEmigDensPar->Visible = false;
		LabelSexEmigDensPar->Visible = false;
		if (CBIndVarEmig->Checked) { // individual variability
			SexEmigPar->ColCount = 3;
			SexEmigPar->Cells[1][0] = "d mean"; SexEmigPar->Cells[2][0] = "d s.d.";
				PanelEP->Visible = true; PanelEP->Enabled = true;
				edtEP->Enabled = false;
				edtEPmean->Enabled = false; edtEPsd->Enabled = false;
				edtEPscale->Enabled = true;
		}
		else { // no individual variability
			PanelEP->Visible = false;
			SexEmigPar->ColCount = 2; SexEmigPar->Cells[1][0] = "d";
		}
	}
	else { // density-dependent
		SexEmigPar->Enabled = false; SexEmigPar->Visible = false;
		LabelSexEmigPar->Visible = false;
		SexEmigDensPar->Enabled = true; SexEmigDensPar->Visible = true;
		LabelSexEmigDensPar->Visible = true;
		if (CBIndVarEmig->Checked)	{ // Individual variability
			SexEmigDensPar->ColCount = 7;
			SexEmigDensPar->Cells[1][0] = "D0 mean";
			SexEmigDensPar->Cells[2][0] = "D0 s.d.";
			SexEmigDensPar->Cells[3][0] = "Alpha mean";
			SexEmigDensPar->Cells[4][0] = "Alpha s.d.";
			SexEmigDensPar->Cells[5][0] = "Beta mean";
			SexEmigDensPar->Cells[6][0] = "Beta s.d.";
				PanelDensEmig->Visible = true;
				edtD0->Enabled = false; edtD0Mean->Enabled = false; edtD0SD->Enabled = false;
				edtAlpha->Enabled = false; edtAlphaMean->Enabled = false; edtAlphaSD->Enabled = false;
				edtBeta->Enabled = false; edtBetaMean->Enabled = false; edtBetaSD->Enabled = false;
				edtD0Scale->Enabled = true;
				edtAlphaScale->Enabled = true; edtBetaScale->Enabled = true;
		}
		else { // no individual variability
			PanelDensEmig->Visible = false;
			SexEmigDensPar->ColCount = 4;
			SexEmigDensPar->Cells[1][0] = "D0";
			SexEmigDensPar->Cells[2][0] = "Alpha";
			SexEmigDensPar->Cells[3][0] = "Beta";
		}
	}

	if (CBSexEmig->Checked && !CBStageEmig->Checked) { // only sex-dependent
		if (RGEmigProb->ItemIndex == 0) { // density-independent
			SexEmigPar->RowCount = 3;
			SexEmigPar->Cells[0][0] = "Sex";
			SexEmigPar->Cells[0][1] = "f";
			SexEmigPar->Cells[0][2] = "m";
			for (int i = 1; i < SexEmigPar->ColCount; i++) {
				for (int j = 1; j < SexEmigPar->RowCount; j++) SexEmigPar->Cells[i][j] = "0.0";
			}
		}
		else { // density-dependent
			SexEmigDensPar->RowCount = 3;
			SexEmigDensPar->Cells[0][0] = "Sex";
			SexEmigDensPar->Cells[0][1] = "f";
			SexEmigDensPar->Cells[0][2] = "m";
			for (int i = 1; i < SexEmigDensPar->ColCount; i++) {
				for (int j = 1; j < SexEmigDensPar->RowCount; j++)
											SexEmigDensPar->Cells[i][j] = "0.0";
			}
		}
	}
	else {
		if (!CBSexEmig->Checked && CBStageEmig->Checked) { // only stage-dependent
			if (RGEmigProb->ItemIndex == 0) { // density-independent
				resetSexEmigPar(1);
			}
			else { // density-dependent
				resetSexEmigDensPar(1);
			}
		}
		else {
			if (CBSexEmig->Checked && CBStageEmig->Checked) { // both dependencies
				if (RGEmigProb->ItemIndex == 0) { // density-independent
					resetSexEmigPar(2);
				}
				else { // density-dependent
					resetSexEmigDensPar(2);
				}
			}
		}
	}
}
else { // neither sex- nor stage-dependent
	PanelSexEmig->Visible = false; PanelSexEmig->Enabled = false;
	PanelSexEmig->Enabled = false; PanelSexEmig->Visible = false;
	if (RGEmigProb->ItemIndex == 0) { // density-independent
		PanelEP->Visible = true; PanelEP->Enabled = true;
		PanelDensEmig->Visible = false; PanelDensEmig->Enabled = false;
		if (CBIndVarEmig->Checked)	{ // individual variability
			edtEPmean->Enabled = true; edtEPsd->Enabled = true;
		}
		else {
			edtEP->Enabled = true;
		}
	}
	else { // density-dependent
		PanelEP->Visible = false; PanelEP->Enabled = false;
		PanelDensEmig->Visible = true; PanelDensEmig->Enabled = true;
		if (CBIndVarEmig->Checked)	{ // individual variability
			edtD0Mean->Enabled = true; edtD0SD->Enabled = true;
			edtAlphaMean->Enabled = true; edtAlphaSD->Enabled = true;
			edtBetaMean->Enabled = true; edtBetaSD->Enabled = true;
		}
		else {
			edtD0->Enabled = true; edtAlpha->Enabled = true; edtBeta->Enabled = true;
		}
	}
}

}

// Reset table(s) on third panel if any of the checkboxes have been changed
void __fastcall TfrmSpecies::emigCheckboxPanelExit(TObject *Sender)
{
/*
bool changed = false;
emigRules emig = pSpecies->getEmig();
if ((CBSexEmig->Checked && !emig.sexDep) || (!CBSexEmig->Checked && emig.sexDep))
	changed = true;
if ((CBStageEmig->Checked && !emig.stgDep) || (!CBStageEmig->Checked && emig.stgDep))
	changed = true;
if ((CBIndVarEP->Checked && !emig.indVar) || (!CBIndVarEP->Checked && emig.indVar))
	changed = true;
if (changed) {
	MessageDlg("EXECUTING emigCheckboxPanelExit() CHANGED",mtWarning,TMsgDlgButtons() << mbOK,0);
}
else {
	MessageDlg("EXECUTING emigCheckboxPanelExit() no change",mtWarning,TMsgDlgButtons() << mbOK,0);
}
*/
}
//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
void __fastcall TfrmSpecies::CBIndVarEmigClick(TObject *Sender)
{

//if (CBStageEmig->Checked && CBIndVarEP->Checked) {
//	// stage-dependent parameters and individual variability are incompatible
//	CBStageEmig->Checked = false;
//	MessageDlg("Stage-dependent dispersal parameters may not vary between individuals",
//		mtWarning, TMsgDlgButtons() << mbRetry,0);
//	return;
//}
if (CBIndVarEmig->Checked) {
	CBStageEmig->Enabled = false; CBStageEmig->Checked = false;
	if (CBStageModel->Checked) {
		edtEmigStage->Enabled = true; edtEmigStage->Visible = true;
	}
	else {
		edtEmigStage->Enabled = false; edtEmigStage->Visible = false;
	}
}
else {
	edtEmigStage->Enabled = false; edtEmigStage->Visible = false;
	if (CBStageModel->Checked) {
		CBStageEmig->Enabled = true;
	}
	else {
		CBStageEmig->Enabled = false; CBStageEmig->Checked = false;
	}
}

setEmigFields();

}

//---------------------------------------------------------------------------

// Transfer functions

void __fastcall TfrmSpecies::RGMovementsClick(TObject *Sender)
{
switch (RGMovements->ItemIndex){
case 0: // dispersal kernel(s)
	if (RGEmigProb->ItemIndex == 0) { // density-independent emigration
		CBFullKernel->Enabled = true;
	}
	PanelKernels->Enabled = true;      PanelKernels->Visible = true;
	RGKernel->Enabled = true;          RGKernel->Visible = true;
	PanelMortality->Enabled = true;    PanelMortality->Visible = true;
	RGMortality->Enabled = true;       RGMortality->Visible = true;
	BtnMovements->Enabled = false;     BtnMovements->Visible = false;
	CBIndVarKernel->Enabled = true;    CBIndVarKernel->Visible = true;
	PanelSettKernels->Enabled = true;  PanelSettKernels->Visible = true;
	PanelSettProcess->Enabled = false; PanelSettProcess->Visible = false;
	edtMaxStepYear->Enabled = false;   edtMaxStepYear->Visible = false;
	break;
case 1: // movement process
	CBFullKernel->Enabled = false;     CBFullKernel->Checked = false;
	PanelKernels->Enabled = false;     PanelKernels->Visible = false;
	RGKernel->Enabled = false;         RGKernel->Visible = false;
	PanelMortality->Enabled = false;   PanelMortality->Visible = false;
	RGMortality->Enabled = false;      RGMortality->Visible = false;
	BtnMovements->Enabled = true;      BtnMovements->Visible = true;
	CBIndVarKernel->Enabled = false;   CBIndVarKernel->Visible = false;
	PanelSettKernels->Enabled = false; PanelSettKernels->Visible = false;
	PanelSettProcess->Enabled = true;  PanelSettProcess->Visible = true;
	if (CBStageModel->Checked) {
		edtMaxStepYear->Enabled = true; edtMaxStepYear->Visible = true;
	}
	else{
		edtMaxStepYear->Enabled = false; edtMaxStepYear->Visible = false;
	}
	if (CBDensDepSettle->Checked) {
//		if (CBStageSettMovt->Checked) {
//			CBIndVarSettle->Enabled = false; CBIndVarSettle->Checked = false;
//		}
//		else {
//			CBIndVarSettle->Enabled = true;
//		}
		CBIndVarSettle->Enabled = true;
	}
	else {
		CBIndVarSettle->Enabled = false; CBIndVarSettle->Checked = false;
	}
	break;
}
// whenever chage is made between kernel and movement model, reset sex- and stage-
// dependencies in settlement to ensure that table of panel 3 is re-formatted
CBSexSettKern->Checked = false;   CBSexSettMovt->Checked = false;
CBStageSettKern->Checked = false; CBStageSettMovt->Checked = false;
}

//---------------------------------------------------------------------------
void __fastcall TfrmSpecies::BtnMovementsClick(TObject *Sender)
{
landParams ppLand = pLandscape->getLandParams();

if (pSpecies->getMovtHabDim() == 0) // set up habitat-dependent cost and mortality matrices
	pSpecies->createHabCostMort(ppLand.nHab);

frmMove->refresh();
frmMove->Show();
}

//---------------------------------------------------------------------------
void __fastcall TfrmSpecies::RGKernelClick(TObject *Sender)
{
PanelDist->Visible = true;

setKernelFields();

//CBSexKernels->Checked = false;
//CBStageKernels->Checked = false;
}

//---------------------------------------------------------------------------
void __fastcall TfrmSpecies::CBSexKernelsClick(TObject *Sender)
{
setKernelFields();
}
//---------------------------------------------------------------------------
void __fastcall TfrmSpecies::CBStageKernelsClick(TObject *Sender)
{
if (CBStageKernels->Checked) {
	CBIndVarKernel->Enabled = false; CBIndVarKernel->Checked = false;
}
else {
	CBIndVarKernel->Enabled = true;
}

setKernelFields();

}

//---------------------------------------------------------------------------
void __fastcall TfrmSpecies::CBIndVarKernelClick(TObject *Sender)
{

if (CBIndVarKernel->Checked) {
	CBStageKernels->Enabled = false; CBStageKernels->Checked = false;
}
else {
  CBStageKernels->Enabled = true;
}

setKernelFields();

}

//---------------------------------------------------------------------------
void __fastcall TfrmSpecies::setKernelFields(void) {

if (CBIndVarKernel->Checked) { // individual variability
	KernelGraph->Visible = false;
	BtnUpdKernel->Visible = false; BtnUpdKernel->Enabled = false;
	edtDist1->Enabled = false;
	edtDist1Mean->Enabled = true;  edtDist1Mean->Visible = true;
	edtDist1SD->Enabled = true;    edtDist1SD->Visible = true;
	edtDist1Scale->Enabled = true; edtDist1Scale->Visible = true;
	if (RGKernel->ItemIndex == 1) { // double kernel
		edtDist2->Enabled = false;	    edtDist2->Visible = true;
		edtPKern1->Enabled = false;	    edtPKern1->Visible = true;
		edtDist2Mean->Enabled = true;		edtDist2Mean->Visible = true;
		edtDist2SD->Enabled = true;		  edtDist2SD->Visible = true;
		edtDist2Scale->Enabled = true;	edtDist2Scale->Visible = true;
		edtPKern1Mean->Enabled = true;  edtPKern1Mean->Visible = true;
		edtPKern1SD->Enabled = true;    edtPKern1SD->Visible = true;
		edtPKern1Scale->Enabled = true; edtPKern1Scale->Visible = true;
	}
	else { // single kernel
		edtDist2->Enabled = false;	    edtDist2->Visible = false;
		edtPKern1->Enabled = false;	    edtPKern1->Visible = false;
		edtDist2Mean->Enabled = false;	 edtDist2Mean->Visible = false;
		edtDist2SD->Enabled = false;		 edtDist2SD->Visible = false;
		edtDist2Scale->Enabled = false;	 edtDist2Scale->Visible = false;
		edtPKern1Mean->Enabled = false;  edtPKern1Mean->Visible = false;
		edtPKern1SD->Enabled = false;    edtPKern1SD->Visible = false;
		edtPKern1Scale->Enabled = false; edtPKern1Scale->Visible = false;
	}
}
else { // no individual variability
	KernelGraph->Visible = true;
	BtnUpdKernel->Visible = true; BtnUpdKernel->Enabled = true;
	edtDist1->Enabled = true;
	edtDist1Mean->Enabled = false;	 edtDist1Mean->Visible = false;
	edtDist1SD->Enabled = false;		 edtDist1SD->Visible = false;
	edtDist1Scale->Enabled = false;	 edtDist1Scale->Visible = false;
	edtDist2Mean->Enabled = false;	 edtDist2Mean->Visible = false;
	edtDist2SD->Enabled = false;	   edtDist2SD->Visible = false;
	edtDist2Scale->Enabled = false;	 edtDist2Scale->Visible = false;
	edtPKern1Mean->Enabled = false;  edtPKern1Mean->Visible = false;
	edtPKern1SD->Enabled = false;    edtPKern1SD->Visible = false;
	edtPKern1Scale->Enabled = false; edtPKern1Scale->Visible = false;
	if (RGKernel->ItemIndex == 1) { // double kernel
		edtDist2->Enabled = true;	 edtDist2->Visible = true;
		edtPKern1->Enabled = true; edtPKern1->Visible = true;
	}
	else { // single kernel
		edtDist2->Enabled = false;	edtDist2->Visible = false;
		edtPKern1->Enabled = false; edtPKern1->Visible = false;
	}
}

if (CBSexKernels->Checked || CBStageKernels->Checked) {
	PanelSexTransfer->Enabled = true; PanelSexTransfer->Visible = true;
	PanelDist->Visible = false;
	if (CBIndVarKernel->Checked) {
		VarTextPanel->Visible = true;
		if (RGKernel->ItemIndex == 0) { // single neg. exp.
			VarLabel2->Visible = false; VarLabel3->Visible = false;
		}
		else { // double neg. exp.
			VarLabel2->Visible = true;  VarLabel3->Visible = true;
		}
	}
	else {
		VarTextPanel->Visible = false;
  }
	KernelGraph->Visible = false; BtnUpdKernel->Visible = false;
	if (RGKernel->ItemIndex == 0) { // single neg. exp.
		if (CBIndVarKernel->Checked) { // individual variability
			SexKernPar->ColCount = 3;
			SexKernPar->Cells[1][0] = "dist mean (m)";
			SexKernPar->Cells[2][0] = "dist s.d. (m)";
		}
		else { // no individual variability
			SexKernPar->ColCount = 2;
			SexKernPar->Cells[1][0] = "mean dist (m)";
		}
	}
	else { // double neg. exp.
		if (CBIndVarKernel->Checked) { // individual variability
			SexKernPar->ColCount = 7;
			SexKernPar->Cells[1][0] = "dist I mean (m)";
			SexKernPar->Cells[2][0] = "dist I s.d. (m)";
			SexKernPar->Cells[3][0] = "dist II mean (m)";
			SexKernPar->Cells[4][0] = "dist II s.d. (m)";
			SexKernPar->Cells[5][0] = "P kernel I mean";
			SexKernPar->Cells[6][0] = "P kernel I s.d.";
		}
		else { // no individual variability
			SexKernPar->ColCount = 4;
			SexKernPar->Cells[1][0] = "mean Dist I(m)";
			SexKernPar->Cells[2][0] = "mean Dist II(m)";
			SexKernPar->Cells[3][0] = "P kernel I";
		}
	}

	if (CBSexKernels->Checked && !CBStageKernels->Checked) { // only sex-dependent
		SexKernPar->RowCount = 3;
		SexKernPar->Cells[0][0] = "Sex";
		SexKernPar->Cells[0][1] = "f";
		SexKernPar->Cells[0][2] = "m";
		for (int i = 1; i < SexKernPar->ColCount; i++) {
			for (int j = 1; j < SexKernPar->RowCount; j++)
				SexKernPar->Cells[i][j] = "0.0";
		}
	}
	if (!CBSexKernels->Checked && CBStageKernels->Checked) { // only stage-dependent
		resetSexTransferPar(1);
	}
	if (CBSexKernels->Checked && CBStageKernels->Checked) { // sex- and stage-dependent
		resetSexTransferPar(2);
	}
}
else { // not sex-dependent and not stage-dependent
	PanelSexTransfer->Enabled = false; PanelSexTransfer->Visible = false;
	PanelDist->Visible = true;
	VarTextPanel->Visible = false;
	if (CBIndVarKernel->Checked) {
		KernelGraph->Visible = false; BtnUpdKernel->Visible = false;
	}
	else {
		KernelGraph->Visible = true; BtnUpdKernel->Visible = true;
	}
}

}

//---------------------------------------------------------------------------
void __fastcall TfrmSpecies::BtnUpdKernelClick(TObject *Sender)
{
double k1,k2;
double max0,max1,max2,pK1;
double meanD1 = (double)StrToFloat(edtDist1->Text);
double meanD2 = (double)StrToFloat(edtDist2->Text);
if (RGKernel->ItemIndex == 0) pK1 = 1.0; // single negative exponential
else pK1 = (double)StrToFloat(edtPKern1->Text);

if (meanD1 < 1.0 || meanD2 < 1.0) { // invalid values - do not update chart
	return;
}

Kernel1->Clear();
Kernel2->Clear();

// NEW METHOD to show weighted curves
// take maximum to be the 99th percentile
max0 = -(std::log(0.01)*meanD1);
if (RGKernel->ItemIndex == 1) { // double negative exponential
	max2 = -(std::log(0.00001*meanD2)*meanD2);
	if (max2 > max0) max0 = max2;
}
#if RSDEBUG
//MemoLine(("TfrmSpecies::btnUpdKernelClick(): meanD1=" + Float2Str(meanD1)
//	+ " meanD2="	+ Float2Str(meanD2)
//	+ " max0="	+ Float2Str(max0)
//	).c_str());
#endif
for (int i = 0; i <= 100; i++) // draw 101 points on graph regardless of maximum
{
	k1 = (1.0/meanD1)*exp(-((double)i*max0/100.0)/meanD1);
	if (RGKernel->ItemIndex == 1) { // double negative exponential
		k1 *= pK1;
		k2 = (1.0-pK1)*(1.0/meanD2)*exp(-((double)i*max0/100.0)/meanD2);
	}
	Kernel1->AddXY(((double)i*max0/100.0),k1,"",clBlue);
	if (RGKernel->ItemIndex == 1) // double negative exponential
		Kernel2->AddXY(((double)i*max0/100.0),k2);
}

}
//---------------------------------------------------------------------------
void __fastcall TfrmSpecies::RGMortalityClick(TObject *Sender)
{
switch(RGMortality->ItemIndex){
case 0:
	edtMortProb->Enabled = true;
	edtMortSlope->Enabled = false; edtMortSlope->Visible = false;
	edtMortInfl->Enabled = false;  edtMortInfl->Visible = false;
	break;
case 1:
	edtMortProb->Enabled = false;
	edtMortSlope->Enabled = true; edtMortSlope->Visible = true;
	edtMortInfl->Enabled = true;  edtMortInfl->Visible = true;
	break;
}
}
//---------------------------------------------------------------------------

// Settlement functions

void __fastcall TfrmSpecies::RGSettKernClick(TObject *Sender)
{
if (CBStageModel->Checked == false) {
	if (RGSettKern->ItemIndex == 1 || RGSettKern->ItemIndex == 3) {
		MessageDlg("This option is possible only with overlapping generations.\n"
			"Please choose between the first and third options.",
			mtError,TMsgDlgButtons() << mbOK,0);
		RGSettKern->ItemIndex = 0;
	}
}
}
//---------------------------------------------------------------------------
void __fastcall TfrmSpecies::RGStepsClick(TObject *Sender)
{
if (RGSteps->ItemIndex == 0) {
	edtNsteps->Enabled = true; edtNsteps->Visible = true;
	edtMinSteps->Enabled = true;
}
else {
	edtNsteps->Enabled = false; edtNsteps->Visible = false;
	edtMinSteps->Enabled = true;
 }
}

//---------------------------------------------------------------------------
void __fastcall TfrmSpecies::CBIndVarSettleClick(TObject *Sender)
{
setSettleFields();
}

//---------------------------------------------------------------------------
void __fastcall TfrmSpecies::CBSexSettMovtClick(TObject *Sender)
{
if (CBSexSettMovt->Checked || CBStageSettMovt->Checked) {
	CBFindMate->Enabled = false; CBFindMate->Checked = false;
}
else {
	CBFindMate->Enabled = true;
}
setSettleFields();
}

//---------------------------------------------------------------------------
void __fastcall TfrmSpecies::CBStageSettMovtClick(TObject *Sender)
{
if (CBStageSettMovt->Checked || CBSexSettMovt->Checked) {
	CBFindMate->Enabled = false; CBFindMate->Checked = false;
//	if (CBDensDepSettle->Checked) {
//		CBIndVarSettle->Enabled = false; CBIndVarSettle->Checked = false;
//	}
}
else {
	if (RGReproduction->ItemIndex == 0) CBFindMate->Enabled = false;
	else CBFindMate->Enabled = true;
//	if (CBDensDepSettle->Checked) {
//		CBIndVarSettle->Enabled = true;
//	}
}
setSettleFields();
}

//---------------------------------------------------------------------------
void __fastcall TfrmSpecies::CBDensDepSettleClick(TObject *Sender)
{
if (CBDensDepSettle->Checked) {
//	if (CBStageSettMovt->Checked) {
//		CBIndVarSettle->Enabled = false; CBIndVarSettle->Checked = false;
//	}
//	else {
//		CBIndVarSettle->Enabled = true;
//	}
	CBIndVarSettle->Enabled = true;
}
else {
	CBIndVarSettle->Enabled = false; CBIndVarSettle->Checked = false;
}
setSettleFields();
}

//---------------------------------------------------------------------------
void __fastcall TfrmSpecies::setSettleFields(void) {

int ii;

// If there is individual variability in density dependence for a movement model,
// then it may NOT be sex- or stage-dependent. However, sex- and/or stage-dependent
// mate-finding requirement and maximum number of steps are permitted.

//MessageDlg("EXECUTING setSettleFields()",mtWarning,TMsgDlgButtons() << mbOK,0);

// Tab sheet Dispersal
if (CBSexSettMovt->Checked || CBStageSettMovt->Checked) {
	edtMaxStepYear->Visible = false; Memo2->Visible = false;
	if (CBIndVarSettle->Checked) {
		PanelSettleDD->Enabled = true; PanelSettleDD->Visible = true;
		edtSettS0->Enabled = false;
		edtSettAlpha->Enabled = false;      edtSettBeta->Enabled = false;
		edtSettS0Mean->Enabled = true;      edtSettS0Mean->Visible = true;
		edtSettAlphaMean->Enabled = true;   edtSettAlphaMean->Visible = true;
		edtSettBetaMean->Enabled = true;    edtSettBetaMean->Visible = true;
		edtSettS0SD->Enabled = true;        edtSettS0SD->Visible = true;
		edtSettAlphaSD->Enabled = true;     edtSettAlphaSD->Visible = true;
		edtSettBetaSD->Enabled = true;      edtSettBetaSD->Visible = true;
		edtSettS0Scale->Enabled = true;     edtSettS0Scale->Visible = true;
		edtSettAlphaScale->Enabled = true;  edtSettAlphaScale->Visible = true;
		edtSettBetaScale->Enabled = true;   edtSettBetaScale->Visible = true;
	}
	else {
		PanelSettleDD->Enabled = false; PanelSettleDD->Visible = false;
	}
}
else { // !(CBSexSettle->Checked || CBStageSettle->Checked)
	if (CBStageModel->Checked) {
		edtMaxStepYear->Visible = true; Memo2->Visible = true;
	}
	else {
		edtMaxStepYear->Visible = false; Memo2->Visible = false;
	}
	if (CBDensDepSettle->Checked) {
		PanelSettleDD->Enabled = true; PanelSettleDD->Visible = true;
		if (CBIndVarSettle->Checked) {
			edtSettS0->Enabled = false;
			edtSettAlpha->Enabled = false;     edtSettBeta->Enabled = false;
			edtSettS0Mean->Enabled = true;     edtSettS0Mean->Visible = true;
			edtSettAlphaMean->Enabled = true;  edtSettAlphaMean->Visible = true;
			edtSettBetaMean->Enabled = true;   edtSettBetaMean->Visible = true;
			edtSettS0SD->Enabled = true;       edtSettS0SD->Visible = true;
			edtSettAlphaSD->Enabled = true;    edtSettAlphaSD->Visible = true;
			edtSettBetaSD->Enabled = true;     edtSettBetaSD->Visible = true;
			edtSettS0Scale->Enabled = true;    edtSettS0Scale->Visible = true;
			edtSettAlphaScale->Enabled = true; edtSettAlphaScale->Visible = true;
			edtSettBetaScale->Enabled = true;  edtSettBetaScale->Visible = true;
		}
		else {
			edtSettS0->Enabled = true;
			edtSettAlpha->Enabled = true;       edtSettBeta->Enabled = true;
			edtSettS0Mean->Enabled = false;     edtSettS0Mean->Visible = false;
			edtSettAlphaMean->Enabled = false;  edtSettAlphaMean->Visible = false;
			edtSettBetaMean->Enabled = false;   edtSettBetaMean->Visible = false;
			edtSettS0SD->Enabled = false;       edtSettS0SD->Visible = false;
			edtSettAlphaSD->Enabled = false;    edtSettAlphaSD->Visible = false;
			edtSettBetaSD->Enabled = false;     edtSettBetaSD->Visible = false;
			edtSettS0Scale->Enabled = false;    edtSettS0Scale->Visible = false;
			edtSettAlphaScale->Enabled = false; edtSettAlphaScale->Visible = false;
			edtSettBetaScale->Enabled = false;  edtSettBetaScale->Visible = false;
		}
	}
	else {
		PanelSettleDD->Enabled = false; PanelSettleDD->Visible = false;
	}
}

// Tab sheet Sex-/stage-dependent dispersal
SexSettleParLabel->Caption = "SETTLEMENT - Movement Processes  ";

PanelSexSettle->Enabled = true; PanelSexSettle->Visible = true;
if (CBIndVarSettle->Checked) {
	SexSettlePar->ColCount = 3;
	SexSettlePar->Cells[1][0] = "find mate";
	SexSettlePar->Cells[2][0] = "max. step/year";
}
else {
	if (CBDensDepSettle->Checked) {
		SexSettlePar->ColCount = 6;
		SexSettlePar->Cells[1][0] = "find mate";
		SexSettlePar->Cells[2][0] = "max. step/year";
		SexSettlePar->Cells[3][0] = "S0";
		SexSettlePar->Cells[4][0] = "alphaS";
		SexSettlePar->Cells[5][0] = "BetaS";
	}
	else {
		SexSettlePar->ColCount = 3;
		SexSettlePar->Cells[1][0] = "find mate";
		SexSettlePar->Cells[2][0] = "max. step/year";
	}
}

if (CBSexSettMovt->Checked) {
	if (CBStageSettMovt->Checked) {
		resetSexSettlePar(4);
	}
	else {
		SexSettlePar->RowCount = 3;
		SexSettlePar->Cells[0][0] = "Sex";
		SexSettlePar->Cells[0][1] = "f";
		SexSettlePar->Cells[0][2] = "m";
		for (int j = 1; j < SexSettlePar->RowCount; j++){
			SexSettlePar->Cells[1][j] = "0";
			SexSettlePar->Cells[2][j] = "0";
			if (!CBIndVarSettle->Checked) {
				if (CBDensDepSettle->Checked) {
					SexSettlePar->Cells[3][j] = "1.0";
					SexSettlePar->Cells[4][j] = "0.0";
					SexSettlePar->Cells[5][j] = "1.0";
				}
			}
		}
	}
}
else {
	if (CBStageSettMovt->Checked) {
		resetSexSettlePar(3);
	}
	else {
		PanelSexSettle->Enabled = false; PanelSexSettle->Visible = false;
	}
}

}

//---------------------------------------------------------------------------
void __fastcall TfrmSpecies::CBSexSettKernClick(TObject *Sender)
{
int ii;
SexSettleParLabel->Caption = "SETTLEMENT - Dispersal Kernels  ";
CBSettKernMate->Checked = false;

if (CBSexSettKern->Checked == false && CBStageSettKern->Checked == false) {
	PanelSexSettle->Enabled = false;
	PanelSexSettle->Visible = false;
	RGSettKern->Visible = true; RGSettKern->ItemIndex = 0;
	CBSettKernMate->Visible = true;
	if (RGReproduction->ItemIndex != 0) CBSettKernMate->Visible = true;
}
else {
	PanelSexSettle->Enabled = true;
  PanelSexSettle->Visible = true;
	SexSettlePar->ColCount = 3;
  SexSettlePar->Cells[1][0] = "If unsuitable...";
  SexSettlePar->Cells[2][0] = "Mating requir.";
  RGSettKern->Visible = false; RGSettKern->ItemIndex = 0;
  CBSettKernMate->Visible = false;
  //Only sex dependent ---------------------------------------------------
  if (CBSexSettKern->Checked && CBStageSettKern->Checked == false) {
		SexSettlePar->RowCount = 3;
		SexSettlePar->Cells[0][0] = "Sex";
		SexSettlePar->Cells[0][1] = "f";
		SexSettlePar->Cells[0][2] = "m";
		for (int j = 1; j < SexSettlePar->RowCount; j++){
			SexSettlePar->Cells[1][j] = "0"; SexSettlePar->Cells[2][j] = "0";
		}
	}
	else{
		//Only stage dependent -----------------------------------------------
		if (CBSexSettKern->Checked == false && CBStageSettKern->Checked) {
			resetSexSettlePar(1);
		}
		else{
			//Both dependencies -------------------------------------------------
			if (CBSexSettKern->Checked && CBStageSettKern->Checked) {
				resetSexSettlePar(2);
			}
		}
	}
}
}

//---------------------------------------------------------------------------
void __fastcall TfrmSpecies::CBStageSettKernClick(TObject *Sender)
{
int ii;
SexSettleParLabel->Caption = "SETTLEMENT - Dispersal Kernels  ";
CBSettKernMate->Checked = false;

if (CBSexSettKern->Checked == false && CBStageSettKern->Checked == false) {
	PanelSexSettle->Enabled = false;
	PanelSexSettle->Visible = false;
	RGSettKern->Visible = true; RGSettKern->ItemIndex = 0;
	CBSettKernMate->Visible = true;
}
else{
	PanelSexSettle->Enabled = true;
	PanelSexSettle->Visible = true;
	RGSettKern->Visible = false; RGSettKern->ItemIndex = 0;
	CBSettKernMate->Visible = false;
	if (CBSexSettKern->Checked && !CBStageSettKern->Checked) { // only sex-dependent
		SexSettlePar->ColCount = 3;
		SexSettlePar->Cells[1][0] = "If unsuitable...";
		SexSettlePar->Cells[2][0] = "Mating requir.";
		SexSettlePar->RowCount = 3;
		SexSettlePar->Cells[0][0] = "Sex";
		SexSettlePar->Cells[0][1] = "f";
		SexSettlePar->Cells[0][2] = "m";
		for (int j = 1; j < SexSettlePar->RowCount; j++){
			SexSettlePar->Cells[1][j] = "0"; SexSettlePar->Cells[2][j] = "0";
		}
	}
	else{
		if (!CBSexSettKern->Checked && CBStageSettKern->Checked) { // only stage-dependent
			SexSettlePar->ColCount = 2;
			SexSettlePar->Cells[1][0] = "If unsuitable...";
			resetSexSettlePar(1);
		}
		else { // sex- and stage-dependent
			SexSettlePar->ColCount = 3;
			SexSettlePar->Cells[1][0] = "If unsuitable...";
			SexSettlePar->Cells[2][0] = "Mating requir.";
			if (CBSexSettKern->Checked && CBStageSettKern->Checked) {
				resetSexSettlePar(2);
			}
		}
	}
}
}

//---------------------------------------------------------------------------
void __fastcall TfrmSpecies::BtnOKClick(TObject *Sender)
{

double tttt0, tttt1;
string msg;
bool changed;
landParams ppLand = pLandscape->getLandParams();

float resol = (float)ppLand.resol;

formCancelled = false;

// PANEL 1 - Population dynamics
// =======   ===================

// Check Rmax > 1 for non-stage-structured model
if (!CBStageModel->Checked && StrToFloat(edtR->Text) <= 0.0) {
	MessageDlg("Rmax must be greater than zero",mtError,TMsgDlgButtons() << mbOK,0);
	return;
}

// Check carrying capacities (inds/ha)
float k;
float sum_K = 0.0;
bool any_error = false;
if (ppLand.generated || ppLand.rasterType == 2) { // artificial landscape or habitat quality
	sum_K = StrToFloat(edtK->Text);
	if (sum_K <= 0.0) {
		PageControl1->ActivePageIndex = 0;
		MessageDlg("Carrying capacity must be greater than zero"
			 ,mtError,TMsgDlgButtons() << mbOK,0);
		return;
	}
}
else {
	for (int i = 0; i < ppLand.nHab; i++) {
		k = (float)(SGhab->Cells[1][i+1]).ToDouble();
		if (k < 0.0) any_error = true;
		sum_K += k;
	}
	if (any_error || sum_K <= 0.0) {
		PageControl1->ActivePageIndex = 0;
		if (any_error)
			MessageDlg("Error in carrying capacities",mtError,TMsgDlgButtons() << mbOK,0);
		else
			MessageDlg("Total carrying capacity must be greater than zero"
				,mtError,TMsgDlgButtons() << mbOK,0);
		return;
	}
}

// Read carrying capacities (inds/ha) & convert them into inds/cell assuming 100% cover
if (ppLand.generated) {
	k = StrToFloat(edtK->Text);
	k *= ((float)(ppLand.resol*ppLand.resol))/10000.0;
	genLandParams gen = pLandscape->getGenLandParams();
	if (gen.continuous) {
		pSpecies->createHabK(1);
		pSpecies->setHabK(0,k);
	}
	else { // discrete
		pSpecies->createHabK(2);
		pSpecies->setHabK(0,0.0); // matrix (0) has K = 0;
		pSpecies->setHabK(1,k);		// habitat (1) has K = k
	}
}
else {
	if (ppLand.rasterType == 2) { // habitat quality
		pSpecies->createHabK(1);
		k = StrToFloat(edtK->Text);
		k *= ((float)(ppLand.resol*ppLand.resol))/10000.0;
		pSpecies->setHabK(0,k);
	}
	else {
		pSpecies->createHabK(ppLand.nHab);
		for (int i = 0; i < ppLand.nHab; i++) {
			k = (float)(SGhab->Cells[1][i+1]).ToDouble();
			k *= ((float)(ppLand.resol*ppLand.resol))/10000.0;
			pSpecies->setHabK(i,k);
		}
	}
}

// check that there are at least 2 stages including juveniles
int Ns = StrToInt(edtNstages->Text);
if (CBStageModel->Checked && Ns < 2) {
	PageControl1->ActivePageIndex = 0;
	Ns = 2; edtNstages->Text = "2";
	MessageDlg("There must be at least two stages including juveniles. "
		,mtError,TMsgDlgButtons() << mbOK,0);
	return;
}

if (CBStageModel->Checked) {

	// check maximum age
	int max = StrToInt(edtMaxAge->Text);
	string maxmsg0 = "The maximum age must be ";
	string maxmsg1 = "higher than the minimum age of the last stage-class";
//	string maxmsg2 = "not lower than the no. of stages for annual reproduction";
	if (RGReproduction->ItemIndex  < 2) {
		if (max <= StrToInt(MinAges->Cells[1][Ns-1])) {
			PageControl1->ActivePageIndex = 0;
			MessageDlg((maxmsg0+maxmsg1).c_str(),mtError,TMsgDlgButtons() << mbOK,0);
			return;
		}
	}
	else {
		if (max <= StrToInt(MinAges->Cells[1][Ns-1])|| max <= StrToInt(MinAges->Cells[1][Ns])){
			PageControl1->ActivePageIndex = 0;
			MessageDlg((maxmsg0+maxmsg1).c_str(),mtError,TMsgDlgButtons() << mbOK,0);
			return;
		}
	}
	// check fecundities
	tttt0 = 0.0;
	if (RGReproduction->ItemIndex  < 2) {
		for (int i = 2; i < Ns+1; i++)
			tttt0 += StrToFloat(transMatrix->Cells[i][1]);
		if (tttt0 <= 0.0) {
			PageControl1->ActivePageIndex = 0;
			MessageDlg("Adult fecundity must be > 0 to avoid extinction",mtError,
				TMsgDlgButtons() << mbOK,0);
			return;
		}
	}
	else {
		for (int i = 4; i < Ns*2+1; i++)
			if (i%2 == 0) tttt0 += StrToFloat(transMatrix->Cells[i][1]);
		if (tttt0 <= 0.0) {
			PageControl1->ActivePageIndex = 0;
			MessageDlg("Adult female fecundity must be > 0 to avoid extinction",mtError,
				TMsgDlgButtons() << mbOK,0);
			return;
		}
		tttt0 = 0.0;
		for (int i = 3; i < Ns*2+1; i++)
			if (i%2 != 0) tttt0 += StrToFloat(transMatrix->Cells[i][1]);
		if (tttt0 <= 0.0) {
			PageControl1->ActivePageIndex = 0;
			MessageDlg("At least one male adult stage must be reproductive",mtError,
				TMsgDlgButtons() << mbOK,0);
			return;
		}
	}
	// check Survival
	string juvmsg =
		"Juvenile development must be greater than 0 and less than or equal to 1";
	if (RGReproduction->ItemIndex < 2) {
		tttt0 =  StrToFloat(transMatrix->Cells[1][1]);
		if (tttt0 > 0.0) {
			MessageDlg("Juveniles are only allowed to develop to stage 1",mtError,
				TMsgDlgButtons() << mbOK,0);
			transMatrix->Cells[1][1] = "0.0";
			return;
		}
		tttt0 =  StrToFloat(transMatrix->Cells[1][2]);
		if (tttt0 <= 0.0 || tttt0 > 1.0) {
			MessageDlg(juvmsg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
			transMatrix->Cells[1][1] = "0.0";
			return;
		}
	}
	else { // complex sexual model
		tttt0 =  transMatrix->Cells[1][1].ToDouble();
		tttt1 =  transMatrix->Cells[2][1].ToDouble();
		if (tttt0 > 0.0 || tttt1 > 0.0) {
			MessageDlg("Juveniles are only allowed to develop to stage 1",mtError,
				TMsgDlgButtons() << mbOK,0);
			transMatrix->Cells[1][1] = "0.0";
			transMatrix->Cells[2][1] = "0.0";
			return;
		}
		tttt0 =  StrToFloat(transMatrix->Cells[1][2]);
		tttt1 =  StrToFloat(transMatrix->Cells[2][3]);
		if (tttt0 <= 0.0 || tttt0 > 1.0 || tttt1 <= 0.0 || tttt1 > 1.0) {
			MessageDlg(juvmsg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
			transMatrix->Cells[1][1] = "0.0";
			return;
		}
	}

	double sum;
	changed =  false;
	if (RGReproduction->ItemIndex < 2) {
		// check juveniles
		sum = 0.0;
		for (int j = 1; j < Ns+1; j++) sum += transMatrix->Cells[1][j].ToDouble();
		if (sum > 1.0) changed = true;
		// check adults
		for (int i = 2; i < Ns+1; i++){
			sum = 0.0;
			for (int j = 2; j < Ns+1; j++) sum += transMatrix->Cells[i][j].ToDouble();
			if (sum > 1.0) changed = true;
		}
	}
	else {
		// check juveniles
		// males
		sum = transMatrix->Cells[1][1].ToDouble();
		for (int j = 2; j < Ns+1; j++) if(j%2 == 0) sum += transMatrix->Cells[1][j].ToDouble();
		if (sum > 1.0) changed = true;
		// females
		sum = transMatrix->Cells[2][1].ToDouble();
		for (int j = 2; j < Ns+1; j++) if(j%2 != 0) sum += transMatrix->Cells[2][j].ToDouble();
		if (sum > 1.0) changed = true;
		// check adults
		for (int i = 3; i < Ns*2+1; i++){
			sum = 0.0;
			for (int j = 2; j < Ns*2; j++) {
				// males
				if(i%2 != 0 && j%2 == 0) sum += transMatrix->Cells[i][j].ToDouble();
				// females
				if(i%2 == 0 && j%2 != 0) sum += transMatrix->Cells[i][j].ToDouble();
			}
			if (sum > 1.0) changed = true;
		}
	}
	if (changed) {
		PageControl1->ActivePageIndex = 0;
		MessageDlg("Error: The sum of the transition probabilities of each column (stage), "
				 "excluding fecundity, must be <= 1.0.", mtError,TMsgDlgButtons() << mbOK,0);
		return;
	}
	// check density dependence coefficients (if appropriate)
	if (CBDevelopment->Checked) {
		if (StrToFloat(edtDevCoeff->Text) <= 0.0) {
			PageControl1->ActivePageIndex = 0;
			MessageDlg("Development coefficient must by > 0.0",
				mtError,TMsgDlgButtons() << mbOK,0);
			return;
		}
	}
	if (CBSurvival->Checked) {
		if (StrToFloat(edtSurvCoeff->Text) <= 0.0) {
			PageControl1->ActivePageIndex = 0;
			MessageDlg("Survival coefficient must by > 0.0",
				mtError,TMsgDlgButtons() << mbOK,0);
			return;
		}
	}
}

// PANELS 2 and 3 - Dispersal
// ==============   =========

string msgemig = "Emigration ";
string msgdd = "Dispersal distance ";
string msgmeandd = "mean dispersal distance ";
string msg4k2 = "for kernel II ";
string msgmean = "mean ";
string msgsd = "s.d. ";
string msgscale = "scaling factor ";
string msgprob = "probability ";
string msgmaxsett = "S0 ";
string msgsettletype = "Settlement type ";
string msgmaxsteps = "Maximum number of steps ";
string msgrelsd = "may not be less than the s.d. ";
string msgrelmaxsd = "may not be less than the maximum corresponding s.d. ";
string msgpkernI = "Proportion using kernel I ";
string range01incl = "must be between 0 and 1 ";
string range01excl = "must be greater than 0 and less than 1 ";
string range01asym = "must be greater than 0 and less than or equal to 1 ";
string msggt0 = "must be greater than 0 ";
string msgge0 = "may not be less than 0 ";
string msg03 = "must be an integer between 0 and 3";
string msg0or1 = "must be either 0 or 1";

// Emigration

emigRules emig = pSpecies->getEmig();

bool probError = false;
bool sdError = false;

if (CBSexEmig->Checked || CBStageEmig->Checked) {
	// check table on panel 3
	int nsexes = 1; 	// no. of sexes if not sex-dependent
	if (CBSexEmig->Checked) nsexes = 2; // sex-dependent
	int nstages = 1; 	// no. of stages if not stage-dependent
	if (CBStageEmig->Checked) nstages = Ns; // stage-dependent
	if (RGEmigProb->ItemIndex == 0) { // density-independent
		msg = msgemig + msgprob + range01incl;
		if (CBIndVarEmig->Checked) { // individual variability
			double maxSD = 0.0;
			for (int i = 1; i < nsexes*nstages+1; i++) {
				tttt0 = SexEmigPar->Cells[1][i].ToDouble(); // mean
				if (tttt0 < 0.0 || tttt0 > 1.0) {
					PageControl1->ActivePageIndex = 2;
					msg = "d " + msgmean + range01incl;
					MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
					return;
				}
				tttt1 = SexEmigPar->Cells[2][i].ToDouble(); // s.d.
				if (tttt1 <= 0.0 || tttt1 >= 1.0) {
					PageControl1->ActivePageIndex = 2;
					msg = "d " + msgsd + range01excl;
					MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
					return;
				}
				if (tttt1 > maxSD) maxSD = tttt1;
			}
			if (edtEPscale->Text.ToDouble() < maxSD) {
				PageControl1->ActivePageIndex = 1;
				msg = msgemig + msgscale + msgrelmaxsd;
				MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
				return;
			}
		}
		else { // no individual variability
			for (int i = 1; i < nsexes*nstages+1; i++) {
				tttt0 = SexEmigPar->Cells[1][i].ToDouble();
				if (tttt0 < 0.0 || tttt0 > 1.0) probError = true;
			}
		}
	}
	else { // density-dependent
		if (CBIndVarEmig->Checked) { // individual variability
			double maxD0SD,maxAlphaSD,maxBetaSD;
			maxD0SD = maxAlphaSD = maxBetaSD = 0.0;
			for (int i = 1; i < nsexes*nstages+1; i++) {
				tttt0 = SexEmigDensPar->Cells[1][i].ToDouble(); // D0 mean
				if (tttt0 < 0.0 || tttt0 > 1.0) {
					PageControl1->ActivePageIndex = 2;
					msg = "D0 " + msgmean + range01incl;
					MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
					return;
				}
				tttt1 = SexEmigDensPar->Cells[2][i].ToDouble(); // D0 s.d.
				if (tttt1 <= 0.0 || tttt1 >= 1.0) {
					PageControl1->ActivePageIndex = 2;
					msg = "D0 " + msgsd + range01excl;
					MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
					return;
				}
				if (tttt1 > maxD0SD) maxD0SD = tttt1;
				tttt1 = SexEmigDensPar->Cells[4][i].ToDouble(); // alpha s.d.
				if (tttt1 <= 0.0) sdError = true;
				if (tttt1 > maxAlphaSD) maxAlphaSD = tttt1;
				tttt1 = SexEmigDensPar->Cells[6][i].ToDouble(); // beta s.d.
				if (tttt1 <= 0.0) sdError = true;
				if (tttt1 > maxBetaSD) maxBetaSD = tttt1;
			}
			if (sdError) {
				PageControl1->ActivePageIndex = 2;
				msg = "Any " + msgsd + msggt0;
				MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
				return;
			}
			if (edtD0Scale->Text.ToDouble() < maxD0SD) {
				PageControl1->ActivePageIndex = 1;
				msg = "D0 " + msgscale + msgrelmaxsd;
				MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
				return;
			}
			if (edtAlphaScale->Text.ToDouble() < maxAlphaSD) {
				PageControl1->ActivePageIndex = 1;
				msg = "Alpha " + msgscale + msgrelsd;
				MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
				return;
			}
			if (edtBetaScale->Text.ToDouble() < maxBetaSD) {
				PageControl1->ActivePageIndex = 1;
				msg = "Beta " + msgscale + msgrelsd;
				MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
				return;
			}
		}
		else { // no individual variability
			for (int i = 1; i < nsexes*nstages+1; i++) {
				for (int j = 1; j < 2; j++) {
					if (RGEmigProb->ItemIndex == 0) { // density-independent
						tttt0 = SexEmigDensPar->Cells[j][i].ToDouble();
						msg = msgemig + msgprob + range01incl;
					}
					else { // density-dependent
						tttt0 = SexEmigDensPar->Cells[j][i].ToDouble();
						msg = "D0 " + range01incl;
					}
					if (tttt0 < 0.0 || tttt0 > 1.0) probError = true;
				}
			}
		}
	}
	if (probError) {
		PageControl1->ActivePageIndex = 2;
		MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
		return;
	}
}
else { // not sex or stage dependence
	// check entries on panel 2
	if (CBIndVarEmig->Checked) { // individual variability
		if (RGEmigProb->ItemIndex == 0) { // density-independent
			tttt0 = edtEPmean->Text.ToDouble();
			if (tttt0 < 0.0 || tttt0 > 1.0) {
				PageControl1->ActivePageIndex = 1;
				msg = msgemig + msgmean + range01incl;
				MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
				return;
			}
			tttt0 = edtEPsd->Text.ToDouble();
			if (tttt0 <= 0.0 || tttt0 >= 1.0) {
				PageControl1->ActivePageIndex = 1;
				msg = msgemig + msgsd + range01excl;
				MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
				return;
			}
			tttt1 = edtEPscale->Text.ToDouble();
			if (tttt1 <= 0.0 || tttt1 >= 1.0) {
				PageControl1->ActivePageIndex = 1;
				msg = msgemig + msgscale + range01excl;
				MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
				return;
			}
			if (tttt1 < tttt0) {
				PageControl1->ActivePageIndex = 1;
				msg = msgemig + msgscale + msgrelsd;
				MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
				return;
			}
		}
		else { // density-dependent
			tttt0 = edtD0Mean->Text.ToDouble();
			if (tttt0 < 0.0 || tttt0 > 1.0) {
				PageControl1->ActivePageIndex = 1;
				msg = "D0 " + msgmean + range01incl;
				MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
				return;
			}
			tttt0 = edtD0SD->Text.ToDouble();
			if (tttt0 <= 0.0 || tttt0 >= 1.0) {
				PageControl1->ActivePageIndex = 1;
				msg = "D0 " + msgsd + range01excl;
				MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
				return;
			}
			tttt1 = edtD0Scale->Text.ToDouble();
			if (tttt1 <= 0.0 || tttt1 >= 1.0) {
				PageControl1->ActivePageIndex = 1;
				msg = "D0 " + msgscale + range01excl;
				MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
				return;
			}
			if (tttt1 < tttt0) {
				PageControl1->ActivePageIndex = 1;
				msg = "D0 " + msgscale + msgrelsd;
				MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
				return;
			}
			tttt0 = edtAlphaSD->Text.ToDouble();
			if (tttt0 <= 0.0) {
				PageControl1->ActivePageIndex = 1;
				msg = "Alpha " + msgsd + msggt0;
				MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
				return;
			}
			tttt1 = edtAlphaScale->Text.ToDouble();
			if (tttt1 < tttt0) {
				PageControl1->ActivePageIndex = 1;
				msg = "Alpha " + msgscale + msgrelsd;
				MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
				return;
			}
			tttt0 = edtBetaSD->Text.ToDouble();
			if (tttt0 <= 0.0) {
				PageControl1->ActivePageIndex = 1;
				msg = "Beta " + msgsd + msggt0;
				MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
				return;
			}
			tttt1 = edtBetaScale->Text.ToDouble();
			if (tttt1 < tttt0) {
				PageControl1->ActivePageIndex = 1;
				msg = "Beta " + msgscale + msgrelsd;
				MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
				return;
			}
		}
	}
	else { // no individual variability
		if (RGEmigProb->ItemIndex == 0) { // density-independent
			tttt0 = edtEP->Text.ToDouble();
			if (tttt0 < 0.0 || tttt0 > 1.0) {
				PageControl1->ActivePageIndex = 1;
				msg = msgemig + msgprob + range01incl;
				MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
				return;
			}
		}
		else { // density-independent
			tttt0 = edtD0->Text.ToDouble();
			if (tttt0 < 0.0 || tttt0 > 1.0) {
				PageControl1->ActivePageIndex = 1;
				msg = "D0 " + range01incl;
				MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
				return;
			}
		}
	}
}
if (CBStageModel->Checked && CBIndVarEmig->Checked) {
	int emigstage = edtEmigStage->Text.ToInt();
	int nstages   = edtNstages->Text.ToInt();
	if (emigstage < 0 || emigstage >= nstages) {
		PageControl1->ActivePageIndex = 1;
		msg = "Emigration stage must be between 0 and " + Int2Str(nstages-1);
		MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
		return;
	}
}

// Transfer

trfrRules trfr = pSpecies->getTrfr();

//bool rangeError = false;
string gtrres;

if (RGMovements->ItemIndex == 0) { // dispersal kernel(s)
	double minMean; // minimum value allowed for mean of negative exponential kernel
	if (CBFullKernel->Checked) { // full kernel used - no constraint on mean
		minMean = 1;
		gtrres = "must be >= 1m";
	}
	else {
		minMean = (double)resol;
		gtrres = "must be >= the cell resolution";
	}
		if (CBSexKernels->Checked || CBStageKernels->Checked) {
			// check table on panel 3
			int nsexes = 1; 	// no. of sexes if not sex-dependent
			if (CBSexKernels->Checked) nsexes = 2; // sex-dependent
			int nstages = 1; 	// no. of stages if not stage-dependent
			if (CBStageKernels->Checked) nstages = Ns; // stage-dependent
			double maxDist1SD,maxDist2SD,maxPKern1SD;
			maxDist1SD = maxDist2SD = maxPKern1SD = 0.0;
			for (int i = 1; i < nsexes*nstages+1; i++) {
				if (CBIndVarKernel->Checked) { // individual variability
					// NB cannot occur if model is stage-structured - hence there must be 2 lines
					if (SexKernPar->Cells[1][i].ToDouble() < minMean) { // kernel I mean
						PageControl1->ActivePageIndex = 2;
						if (i == 1) msg = "The fe"; else msg = "The ";
						msg += "males' " + msgdd + "mean " + gtrres;
						MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
						return;
					}
					tttt0 = SexKernPar->Cells[2][i].ToDouble();
					if (tttt0 <= 0.0) { // kernel I s.d.
						PageControl1->ActivePageIndex = 2;
						if (i == 1) msg = "The fe"; else msg = "The ";
						PageControl1->ActivePageIndex = 2;
						msg += "males' " + msgdd + "s.d. " + msggt0;
						MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
						return;
					}
					if (tttt0 > maxDist1SD) maxDist1SD = tttt0;
					if (RGKernel->ItemIndex == 1) { // twin kernels
						if (SexKernPar->Cells[3][i].ToDouble() < minMean) { // kernel II mean
							PageControl1->ActivePageIndex = 2;
							if (i == 1) msg = "The fe"; else msg = "The ";
							msg += "males' " + msgdd + msg4k2 + "mean " + gtrres;
							MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
							return;
						}
						tttt0 = SexKernPar->Cells[4][i].ToDouble();
						if (tttt0 <= 0.0) { // kernel II s.d.
							PageControl1->ActivePageIndex = 2;
							if (i == 1) msg = "The fe"; else msg = "The ";
							msg += "males' " + msgdd + msg4k2 + "s.d." + msggt0;
							MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
							return;
						}
						if (tttt0 > maxDist2SD) maxDist2SD = tttt0;
						tttt0 = (float)SexKernPar->Cells[5][i].ToDouble();
						if (tttt0 < 0.0 || tttt0 > 1.0) {
							PageControl1->ActivePageIndex = 2;
							if (i == 1) msg = "The fe"; else msg = "The ";
							msg += "males' " + msgpkernI + "mean " + range01incl;
							MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
							return;
						}
						tttt0 = (float)SexKernPar->Cells[6][i].ToDouble();
						if (tttt0 <= 0.0 || tttt0 >= 1.0) {
							PageControl1->ActivePageIndex = 2;
							if (i == 1) msg = "The fe"; else msg = "The ";
							msg += "males' " + msgpkernI + "s.d. " + range01excl;
							MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
							return;
						}
						if (tttt0 > maxPKern1SD) maxPKern1SD = tttt0;
					}
				}
				else { // no individual variability - can be sex- AND/OR stage-structured
					string msg0;
					if (nstages == 1) msg0 = "The ";
					else {
						msg0 = "Stage ";
						if (nsexes == 1) msg0 += Int2Str((int)(i-1));
						else msg0 += Int2Str((int)((i-1)/2));
					}
					if (nsexes > 1) {
						if (i%2 == 0) msg0 += " "; else msg0 += " fe";
						msg0 += "males' ";
					}
					// kernel I mean
					if (SexKernPar->Cells[1][i].ToDouble() < minMean) {
						if (true)
						// a conditional check on the equivalent emigration proportion is required
						// here (as in case 2), but it would be very complex, as emigration could
						// be constant, stage-dependent, sex-dependent or stage-sex-dependent -
						// a better approach would be to always have a NSxP  matrix (N = no. of stages,
						// S = no. of sexes, P = set of all transfer parameters)) which is populated
						// during edit checking (rather than after the form is fully checked,
						// thereby enabling direct comparison between corresponding transfer parameters
						// SCFP 28/6/13
						{
							PageControl1->ActivePageIndex = 2;
							msg = msg0 + msgmeandd + gtrres;
							MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
							return;
						}
					}
					if (RGKernel->ItemIndex == 1) { // kernel II
						if (SexKernPar->Cells[2][i].ToDouble() < minMean) {
							PageControl1->ActivePageIndex = 2;
							msg = msg0 + msgmeandd + msg4k2 + gtrres;
							MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
							return;
						}
						tttt0 = (float)SexKernPar->Cells[3][i].ToDouble();
						if (tttt0 < 0.0 || tttt0 > 1.0) {
							PageControl1->ActivePageIndex = 2;
							msg = msg0 + msgpkernI + range01incl;
							MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
							return;
						}
					}
				}
			} // end of for loop
			// check standard deviations against corresponding scaling factors
			if (edtDist1Scale->Text.ToDouble() < maxDist1SD) {
				PageControl1->ActivePageIndex = 1;
				msg = msgdd + msgscale + msgrelmaxsd;
//				msg = Double2Str(edtEPscale->Text.ToDouble()) + " v " + Double2Str(maxSD);
				MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
				return;
			}
			if (edtDist2Scale->Text.ToDouble() < maxDist2SD) {
				PageControl1->ActivePageIndex = 1;
				msg = msgdd + msgscale + msg4k2 + msgrelmaxsd;
				MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
				return;
			}
			if (edtPKern1Scale->Text.ToDouble() < maxPKern1SD) {
				PageControl1->ActivePageIndex = 1;
				msg = msgpkernI + msgscale + msgrelmaxsd;
				MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
				return;
			}
		} // end of CBSexKernels->Checked || CBStageKernels->Checked
		else { // !CBSexKernels->Checked && !CBStageKernels->Checked
			// check panel 2
			if (CBIndVarKernel->Checked) { // individual variability
				if(StrToFloat(edtDist1Mean->Text) < minMean) {
					PageControl1->ActivePageIndex = 1;
					msg = msgdd + " mean " + gtrres;
					MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
					return;
				}
				tttt0 = StrToFloat(edtDist1SD->Text);
				if(tttt0 <= 0.0) {
					PageControl1->ActivePageIndex = 1;
					msg = msgdd + " s.d. " + msggt0;
					MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
					return;
				}
				tttt1 = StrToFloat(edtDist1Scale->Text);
				if(tttt1 <= 0.0) {
					PageControl1->ActivePageIndex = 1;
					msg = msgdd + msgscale + msggt0;
					MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
					return;
				}
				if(tttt1 < tttt0) {
					PageControl1->ActivePageIndex = 1;
					msg = msgdd + msgscale + msgrelsd;
					MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
					return;
				}
				if (RGKernel->ItemIndex == 1){ // double exponential kernel
					if(StrToFloat(edtDist2Mean->Text) < minMean){
						PageControl1->ActivePageIndex = 1;
						msg = msgdd + msg4k2 + "mean " + gtrres;
						MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
						return;
					}
					tttt0 = StrToFloat(edtDist2SD->Text);
					if(tttt0 <= 0.0){
						PageControl1->ActivePageIndex = 1;
						msg = msgdd + msg4k2 + "s.d. " + gtrres;
						MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
						return;
					}
					tttt1 = StrToFloat(edtDist2Scale->Text);
					if(tttt1 <= 0.0){
						PageControl1->ActivePageIndex = 1;
						msg = msgdd + msg4k2 + msgscale + gtrres;
						MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
						return;
					}
					if(tttt1 < tttt0){
						PageControl1->ActivePageIndex = 1;
						msg = msgdd + msg4k2 + msgscale + msgrelsd;
						MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
						return;
					}
					tttt0 = StrToFloat(edtPKern1Mean->Text);
					if(tttt0 < 0.0 || tttt0 > 1.0){
						PageControl1->ActivePageIndex = 1;
						msg = msgpkernI + "mean " + range01incl;
						MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
						return;
					}
					tttt0 = StrToFloat(edtPKern1SD->Text);
					if(tttt0 <= 0.0 || tttt0 >= 1.0){
						PageControl1->ActivePageIndex = 1;
						msg = msgpkernI + "s.d. " + range01excl;
						MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
						return;
					}
					tttt1 = StrToFloat(edtPKern1Scale->Text);
					if(tttt1 <= 0.0 || tttt1 >= 1.0){
						PageControl1->ActivePageIndex = 1;
						msg = msgpkernI + msgscale + range01excl;
						MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
						return;
					}
					if(tttt1 < tttt0){
						PageControl1->ActivePageIndex = 1;
						msg = msgpkernI + msgscale + msgrelsd;
						MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
						return;
					}
				}
			}
			else { // no individual variability
				if(StrToFloat(edtDist1->Text) < minMean) {
					PageControl1->ActivePageIndex = 1;
					msg = msgdd + gtrres;
					MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
					return;
				}
				if (RGKernel->ItemIndex == 1){
					if(StrToFloat(edtDist2->Text) < minMean) {
						PageControl1->ActivePageIndex = 1;
						msg = msgdd + msg4k2 + gtrres,
						MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
						return;
					}
					tttt0 = StrToFloat(edtPKern1->Text);
					if (tttt0 <= 0.0 || tttt0 >= 1.0) {
						PageControl1->ActivePageIndex = 1;
						msg = msgpkernI + range01excl;
						MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
						return;
					}
				}
			}
		}

		// SCFP 28/7/13 - the following edits are not correct and do not allow for sex-
		// dependent transfer, but for similar reasons as above, correct checks conditional
		// on emigration proportions will be easier to implement following code revision
		if (CBStageKernels->Checked == true && RGKernel->ItemIndex == 1) {
			// multiple stage double exponential dispersal kernel
			tttt0 = StrToFloat(SexKernPar->Cells[3][1]);
			if (tttt0 <= 0.0 || tttt0 >= 1.0)
			{
				if ((CBStageEmig->Checked && SexEmigPar->Cells[1][1].ToDouble() > 0.0)
				||  (!CBStageEmig->Checked && StrToFloat(edtEP->Text) > 0.0) )
				{
					PageControl1->ActivePageIndex = 2;
					msg = msgpkernI + "for stage 0 " + range01excl;
					MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
					return;
				}
			}
			tttt0 = StrToFloat(SexKernPar->Cells[3][2]);
			if (tttt0 <= 0.0 || tttt0 >= 1.0)
			{
				if ((CBStageEmig->Checked && SexEmigPar->Cells[1][2].ToDouble() > 0.0)
				||  (!CBStageEmig->Checked && StrToFloat(edtEP->Text) > 0.0) ){
				  	PageControl1->ActivePageIndex = 2;
					msg = msgpkernI + "for stage 1 " + range01excl;
					MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
					return;
				}
			}
		}
} // end of kernel-based transfer condition
else { // movement processes
	// all data are on separate form FormMove
	if (!movtModelOK) {
	 	PageControl1->ActivePageIndex = 1;
		msg = "Movement process parameters have not been correctly specified";
		MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
		return;
	}
}

// Settlement

settleType sett;
if (RGMovements->ItemIndex == 0) { // transfer is by dispersal kernel(s)
	// checks on settlement rules for kernel-based dispersal
	int nsexes = 1; 	// no. of sexes if not sex-dependent
	if (CBSexSettKern->Checked) {
		nsexes = 2; // sex-dependent
		if (CBStageSettKern->Checked) {
			sett.sexDep = true; sett.stgDep = true;
		}
		else {
			sett.sexDep = true; sett.stgDep = false;
		}
	}
	else {
		if (CBStageSettKern->Checked) {
			sett.sexDep = false; sett.stgDep = true;
		}
		else {
			sett.sexDep = false; sett.stgDep = false;
		}
	}
	int nstages = 1; 	// no. of stages if not stage-dependent
	if (CBStageSettKern->Checked) nstages = Ns; // stage-dependent
	if (CBSexSettKern->Checked || CBStageSettKern->Checked) {
		// check panel 3 (2 columns)
		for (int i = 1; i < nsexes*nstages+1; i++) {
			if (SexSettlePar->Cells[1][i].ToInt() < 0
			||  SexSettlePar->Cells[1][i].ToInt() > 3) {
				PageControl1->ActivePageIndex = 2;
				msg = msgsettletype;
				if (nsexes > 1) {
					if (i%2 == 0) msg += "for males "; else msg += "for females ";
				}
				msg += msg03;
				MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
				return;
			}
			if (SexSettlePar->Cells[2][i].ToInt() < 0
			||  SexSettlePar->Cells[2][i].ToInt() > 1) {
				PageControl1->ActivePageIndex = 2;
				msg = "The mate-finding condition ";
				if (nsexes > 1) {
					if (i%2 == 0) msg += "for males "; else msg += "for females ";
				}
				msg += msg0or1;
				MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
				return;
			}
		}
	}
	else {
		// check panel 2
		// no checks required - selection by radio-button only

		// WHAT ABOUT MAX NO OF STEPS?....

	}
}
else { // transfer is by movement rule
	// checks on sex- and stage-dependent settlement rules for movement rules
	sett.indVar = CBIndVarSettle->Checked;
	int nsexes = 1; 	// no. of sexes if not sex-dependent
	if (CBSexSettMovt->Checked) {
		nsexes = 2; // sex-dependent
		if (CBStageSettMovt->Checked) {
			sett.sexDep = true; sett.stgDep = true;
		}
		else {
			sett.sexDep = true; sett.stgDep = false;
		}
	}
	else {
		if (CBStageSettMovt->Checked) {
			sett.sexDep = false; sett.stgDep = true;
		}
		else {
			sett.sexDep = false; sett.stgDep = false;
		}
	}
	int nstages = 1; 	// no. of stages if not stage-dependent
	if (CBStageSettMovt->Checked) nstages = Ns; // stage-dependent
	if (CBSexSettMovt->Checked || CBStageSettMovt->Checked) {
		int settletype;
//		double maxS0SD,maxAlphaSD,maxBetaSD;
//		maxS0SD = maxAlphaSD = maxBetaSD = 0.0;
		// check panel 3
		for (int i = 1; i < nsexes*nstages+1; i++) {
			// find mate
			int findmate = SexSettlePar->Cells[1][i].ToInt();
			if (findmate < 0 ||  findmate > 1) {
				PageControl1->ActivePageIndex = 2;
				msg = "find mate " + msg0or1;
				MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
				return;
			}
			// max steps per year
			if (SexSettlePar->Cells[2][i].ToInt() < 0) {
				PageControl1->ActivePageIndex = 2;
				msg = msgmaxsteps + "per year ";
				if (nsexes > 1) {
					if (i%2 == 0) msg += "for males"; else msg += "for females";
				}
				msg += msgge0;
				MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
				return;
			}
			// density dependence parameters
			if (CBIndVarSettle->Checked) {
				/*  these checks are not required as individual variability in density
						dependence precludes sex dependence
						HOWEVER, RETAIN CODE IN CASE SEX DEPENDENCE IS REINTRODUCED
				if (SexSettlePar->Cells[3][i].ToDouble() <= 0.0
				||  SexSettlePar->Cells[3][i].ToDouble() > 1.0) {
					PageControl1->ActivePageIndex = 2;
					msg = "S0 " + msgmean;
					if (nsexes > 1) {
						if (i%2 == 0) msg += "for males "; else msg += "for females ";
					}
					msg += range01asym;
					MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
					return;
				}
				tttt1 = SexSettlePar->Cells[4][i].ToDouble(); // S0 s.d.
				if (tttt1 <= 0.0 || tttt1 > 1.0) {
					PageControl1->ActivePageIndex = 2;
					msg = "S0 " + msgsd;
					if (nsexes > 1) {
						if (i%2 == 0) msg += "for males "; else msg += "for females ";
					}
					msg += range01asym;
					MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
					return;
				}
				if (tttt1 > maxS0SD) maxS0SD = tttt1;
				tttt1 = SexSettlePar->Cells[6][i].ToDouble(); // alpha s.d.
				if (tttt1 <= 0.0) {
					PageControl1->ActivePageIndex = 2;
					msg = "alphaS " + msgsd;
					if (nsexes > 1) {
						if (i%2 == 0) msg += "for males "; else msg += "for females ";
					}
					msg += msggt0;
					MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
					return;
				}
				if (tttt1 > maxAlphaSD) maxAlphaSD = tttt1;
				tttt1 = SexSettlePar->Cells[8][i].ToDouble(); // beta s.d.
				if (tttt1 <= 0.0) {
					PageControl1->ActivePageIndex = 2;
					msg = "betaS " + msgsd;
					if (nsexes > 1) {
						if (i%2 == 0) msg += "for males "; else msg += "for females ";
					}
					msg += msggt0;
					MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
					return;
				}
				if (tttt1 > maxBetaSD) maxBetaSD = tttt1;
				*/
			}
			else { // no individual variation
				if (CBDensDepSettle->Checked) {
					// max. settlement prob. (s0) is the only constrained parameter
					if (SexSettlePar->Cells[3][i].ToDouble() <= 0.0
					||  SexSettlePar->Cells[3][i].ToDouble() > 1.0) {
						PageControl1->ActivePageIndex = 2;
						msg = "S0 ";
						if (nsexes > 1) {
							if (i%2 == 0) msg += "for males "; else msg += "for females ";
						}
						msg += range01asym;
						MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
						return;
					}
				}
			}
		}
		/* AS ABOVE - RETAIN CODE...
		if (edtSettS0Scale->Text.ToDouble() < maxS0SD) {
			PageControl1->ActivePageIndex = 1;
			msg = "S0 " + msgscale + msgrelmaxsd;
			MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
			return;
		}
		if (edtSettAlphaScale->Text.ToDouble() < maxAlphaSD) {
			PageControl1->ActivePageIndex = 1;
			msg = "AlphaS " + msgscale + msgrelsd;
			MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
			return;
		}
		if (edtSettBetaScale->Text.ToDouble() < maxBetaSD) {
			PageControl1->ActivePageIndex = 1;
			msg = "BetaS " + msgscale + msgrelsd;
			MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
			return;
		}
		*/
	}
	else { // not sex- or stage-dependent
		// check panel 2
		if (CBDensDepSettle->Checked) {
			// density dependent settlement
			if (!CBIndVarSettle->Checked) { // no individual variability
				if (StrToFloat(edtSettS0->Text) <= 0.0 || StrToFloat(edtSettS0->Text) > 1.0) {
//					msg = msgmaxsett + msgprob + range01asym;
					msg = msgmaxsett + range01asym;
					MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
					return;
				}
			}
		}
		msg = msgmaxsteps;
		if (RGSteps->ItemIndex == 0) {
			if (StrToInt(edtNsteps->Text) < 1) {
				msg += msggt0;
				MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
				return;
			}
		}
		if (CBStageModel->Checked) {
			if (StrToInt(edtMaxStepYear->Text) < 0) {
				msg += "per year must be zero or greater ";
				MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
				return;
			}
		}
	}
	if (CBIndVarSettle->Checked) { // individual variability
		tttt0 = edtSettS0Mean->Text.ToDouble();
		if (tttt0 < 0.0 || tttt0 > 1.0) {
			PageControl1->ActivePageIndex = 1;
			msg = msgmaxsett + msgmean + range01incl;
			MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
			return;
		}
		tttt0 = edtSettS0SD->Text.ToDouble();
		if (tttt0 <= 0.0 || tttt0 >= 1.0) {
			PageControl1->ActivePageIndex = 1;
			msg = msgmaxsett + msgsd + range01excl;
			MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
			return;
		}
		tttt1 = edtSettS0Scale->Text.ToDouble();
		if (tttt1 <= 0.0 || tttt1 >= 1.0) {
			PageControl1->ActivePageIndex = 1;
			msg = msgmaxsett + msgscale + range01excl;
			MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
			return;
		}
		if (tttt1 < tttt0) {
			PageControl1->ActivePageIndex = 1;
			msg = msgmaxsett + msgscale + msgrelsd;
			MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
			return;
		}
		tttt0 = edtSettAlphaSD->Text.ToDouble();
		if (tttt0 <= 0.0) {
			PageControl1->ActivePageIndex = 1;
			msg = "AlphaS " + msgsd + range01asym;
			MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
			return;
		}
		tttt1 = edtSettAlphaScale->Text.ToDouble();
		if (tttt1 < tttt0) {
			PageControl1->ActivePageIndex = 1;
			msg = "AlphaS " + msgscale + msgrelsd;
			MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
			return;
		}
		tttt0 = edtSettBetaSD->Text.ToDouble();
		if (tttt0 <= 0.0) {
			PageControl1->ActivePageIndex = 1;
			msg = "BetaS" + msgsd + range01asym;
			MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
			return;
		}
		tttt1 = edtSettBetaScale->Text.ToDouble();
		if (tttt1 < tttt0) {
			PageControl1->ActivePageIndex = 1;
			msg = "BetaS " + msgscale + msgrelsd;
			MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
			return;
		}
	}
}

pSpecies->setTrfr(trfr);
pSpecies->setSettle(sett);

if (CBStageModel->Checked) { // stage-structured
	stageParams sstruct = pSpecies->getStage();
	sstruct.nStages = StrToInt(edtNstages->Text);
	pSpecies->setStage(sstruct);
	stageParams sss;
	sss = pSpecies->getStage();
}

// All species checks have been passed, and values can now be saved...
updateSpecies();
speciesOK = true;

Close();

}

//---------------------------------------------------------------------------

void __fastcall TfrmSpecies::updateSpecies(void) {

int z,h,stg;
float k;
float fff,surv;
double ss,dd;
float devCoeff,survCoeff;

landParams ppLand = pLandscape->getLandParams();
demogrParams dem = pSpecies->getDemogr();
stageParams sstruct = pSpecies->getStage();

// Reproduction parameters

dem.repSeasons = StrToFloat(edtGen->Text);
dem.repType = RGReproduction->ItemIndex;

if (dem.repType != 0) dem.propMales = StrToFloat(edtSexRatio->Text);
if (dem.repType == 2) dem.harem = StrToInt(edtHarem->Text);

if (CBStageModel->Checked) { // stage-structured
	dem.stageStruct = true;
	if (ppLand.dynamic) {
		if (RGaction->ItemIndex == 0) sstruct.disperseOnLoss = false;
		else sstruct.disperseOnLoss = true;
	}
	else sstruct.disperseOnLoss = false;

	sstruct.probRep = StrToFloat(edtPRep->Text);
	sstruct.repInterval = StrToInt(edtRepInterval->Text);

	sstruct.maxAge = StrToInt(edtMaxAge->Text);
	sstruct.survival = RGSurvival->ItemIndex;

	sstruct.fecDens = CBFecundity->Checked;
	sstruct.devDens = CBDevelopment->Checked;
	sstruct.survDens = CBSurvival->Checked;
	sstruct.fecStageDens = CBweightFec->Checked;
	sstruct.devStageDens = CBweightDev->Checked;
	sstruct.survStageDens = CBweightSurv->Checked;

	if (dem.repType != 2) { // ASEXUAL or IMPLICIT 2 sex model
		for (int i = 2; i < sstruct.nStages+1; i++) {
			fff = (float)(transMatrix->Cells[i][1]).ToDouble();
			pSpecies->setFec(i-1,0,fff);
			if (fff > dem.lambda) dem.lambda = fff;
		}
		for (int i = 1; i < sstruct.nStages; i++)
			pSpecies->setMinAge(i,0,(MinAges->Cells[1][i]).ToInt());

		for (int i = 1; i < sstruct.nStages+1; i++) {
			ss = 0.0; dd = 0.0;
			for (int j = 1; j < sstruct.nStages+1; j++) {
				if (j == i)     ss = (transMatrix->Cells[i][j]).ToDouble();
				if (j == (i+1)) dd = (transMatrix->Cells[i][j]).ToDouble();
			}
			pSpecies->setSurv(i-1,0,ss+dd);
			if ((ss+dd) > 0)
				pSpecies->setDev(i-1,0,dd/(ss+dd));
			else
				pSpecies->setDev(i-1,0,0.0);
		}

		// stage-specific density dependence coefficients & stages' weights
		if (sstruct.fecDens && sstruct.fecStageDens) { // fecundity
			// create stage weights matrix
			pSpecies->createDDwtFec(sstruct.nStages);
			// read coefficients
			for (int i = 0; i < sstruct.nStages; i++){
				for (int j = 0; j < sstruct.nStages; j++)
					pSpecies->setDDwtFec(i,j,(float)frmDensity->FecundityDens->Cells[i+1][j+1].ToDouble());
			}
		}
		if (sstruct.devDens) { // development
			devCoeff = StrToFloat(edtDevCoeff->Text);
			if (sstruct.devStageDens) {
				// create stage weights matrix
				pSpecies->createDDwtDev(sstruct.nStages);
				// read coefficients
				for (int i = 0; i < sstruct.nStages; i++){
					for (int j = 0; j < sstruct.nStages; j++)
						pSpecies->setDDwtDev(i,j,(float)frmDensity->DevelopDens->Cells[i+1][j+1].ToDouble());
				}
			}
		}
		else devCoeff = 1.0;
		if (sstruct.survDens) { // survival
			survCoeff = StrToFloat(edtSurvCoeff->Text);
			if (sstruct.survStageDens) {
				// create stage weights matrix
				pSpecies->createDDwtSurv(sstruct.nStages);
				// read coefficients
				for (int i = 0; i < sstruct.nStages; i++){
					for (int j = 0; j < sstruct.nStages; j++)
						pSpecies->setDDwtSurv(i,j,(float)frmDensity->SurvivalDens->Cells[i+1][j+1].ToDouble());
				}
			}
		}
		else survCoeff = 1.0;
		if (sstruct.devDens || sstruct.survDens) {
			pSpecies->setDensDep(devCoeff,survCoeff);
		}
	}  // end of dem.repType != 2
	else { // complex 2 sex model
		stg = 1;
		string msg;
		for (int i = 1; i < sstruct.nStages*2-1; i++) {
			if (i%2 != 0) { // male
				pSpecies->setMinAge(stg,1,(MinAges->Cells[1][i]).ToInt());
			}
			else { // female
				pSpecies->setMinAge(stg,0,(MinAges->Cells[1][i]).ToInt());
				stg++;
			}
		}

		stg = 1;
		for (int i = 3; i < sstruct.nStages*2+1; i++) {
			fff = (float)(transMatrix->Cells[i][1]).ToDouble();
			if (i%2 != 0) { // males
				pSpecies->setFec(stg,1,fff);
			}
			else { // females
				pSpecies->setFec(stg,0,fff);
				if (fff > dem.lambda) dem.lambda = fff;
				stg++;
			}
		}

		ss = (transMatrix->Cells[1][1]).ToDouble();
		dd = (transMatrix->Cells[1][2]).ToDouble();
		pSpecies->setSurv(0,1,ss+dd);
		if ((ss+dd) > 0)
			pSpecies->setDev(0,1,dd/(ss+dd));
		else
			pSpecies->setDev(0,1,0.0);
		ss = (transMatrix->Cells[2][1]).ToDouble();
		dd = (transMatrix->Cells[2][3]).ToDouble();
		pSpecies->setSurv(0,0,ss+dd);
		if ((ss+dd) > 0)
			pSpecies->setDev(0,0,dd/(ss+dd));
		else
			pSpecies->setDev(0,0,0.0);
		stg = 1;
		for (int i = 3; i < sstruct.nStages*2+1; i++) {
			ss = 0.0; dd = 0.0;

			if (i%2 != 0) { // males
				for (int j = 2; j < sstruct.nStages*2; j++) {
					if (i == j+1) ss = (transMatrix->Cells[i][j]).ToDouble();
					if (i == j-1) dd = (transMatrix->Cells[i][j]).ToDouble();
				}
				pSpecies->setSurv(stg,1,ss+dd);
				if ((ss+dd) > 0)
					pSpecies->setDev(stg,1,dd/(ss+dd));
				else
					pSpecies->setDev(stg,1,0.0);
			}
			else { // females
				for (int j = 1; j < sstruct.nStages*2; j++) {
					if (i == j+1) ss = (transMatrix->Cells[i][j]).ToDouble();
					if (i == j-1) dd = (transMatrix->Cells[i][j]).ToDouble();
				}
				pSpecies->setSurv(stg,0,ss+dd);
				if ((ss+dd) > 0)
					pSpecies->setDev(stg,0,dd/(ss+dd));
				else
					pSpecies->setDev(stg,0,0.0);
				stg++;
			}
		}

		// stage-specific density dependence coefficients
		if (sstruct.fecDens && sstruct.fecStageDens) { // fecundity
			// create stage weights matrix
			pSpecies->createDDwtFec(sstruct.nStages*NSEXES);
			// read coefficients
			for (int i = 0; i < sstruct.nStages*2; i++){
				for (int j = 0; j < sstruct.nStages*2; j++)
					pSpecies->setDDwtFec(i,j,(float)frmDensity->FecundityDens->Cells[i+1][j+1].ToDouble());
			}
		}
		if (sstruct.devDens) { // development
			devCoeff = StrToFloat(edtDevCoeff->Text);
			if (sstruct.devStageDens) {
				// create stage weights matrix
				pSpecies->createDDwtDev(sstruct.nStages*NSEXES);
				// read coefficients
				for (int i = 0; i < sstruct.nStages*2; i++){
					for (int j = 0; j < sstruct.nStages*2; j++)
						pSpecies->setDDwtDev(i,j,(float)frmDensity->DevelopDens->Cells[i+1][j+1].ToDouble());
				}
			}
		}
		else devCoeff = 1.0;
		if (sstruct.survDens) { // survival
			survCoeff = StrToFloat(edtSurvCoeff->Text);
			if (sstruct.survStageDens) {
				// create stage weights matrix
				pSpecies->createDDwtSurv(sstruct.nStages*NSEXES);
				// read coefficients
				for (int i = 0; i < sstruct.nStages*2; i++){
					for (int j = 0; j < sstruct.nStages*2; j++)
						pSpecies->setDDwtSurv(i,j,(float)frmDensity->SurvivalDens->Cells[i+1][j+1].ToDouble());
				}
			}
		}
		else survCoeff = 1.0;
		if (sstruct.devDens || sstruct.survDens) {
			pSpecies->setDensDep(devCoeff,survCoeff);
		}
	} // end of explicit 2 sex option
	pSpecies->setStage(sstruct);
} // end of stage-structured
else { // NO stage structure
	dem.stageStruct = false;
	dem.lambda = StrToFloat(edtR->Text);
	dem.bc =  StrToFloat(edtC->Text);
}
pSpecies->setDemogr(dem);

// Dispersal parameters

// Emigration

emigRules emig = pSpecies->getEmig();
emigTraits etraits;

if (RGEmigProb->ItemIndex == 1) emig.densDep = true;
else emig.densDep = false;
if (CBIndVarEmig->Checked == false) emig.indVar = false;
else emig.indVar = true;
if (emig.indVar) {
	emigScales s;
	if (emig.densDep) {
		s.d0Scale    = StrToFloat(edtD0Scale->Text);
		s.alphaScale = StrToFloat(edtAlphaScale->Text);
		s.betaScale  = StrToFloat(edtBetaScale->Text);
	}
	else {
		s.d0Scale = StrToFloat(edtEPscale->Text);
	}
	pSpecies->setEmigScales(s);
}

if (CBSexEmig->Checked) emig.sexDep = true; else emig.sexDep = false;
if (CBStageEmig->Checked) emig.stgDep = true;
else { // not stage-dependent
	emig.stgDep = false;
	if (emig.indVar) emig.emigStage = StrToInt(edtEmigStage->Text);
}

pSpecies->setFullKernel(CBFullKernel->Checked);
pSpecies->setEmig(emig);

emigParams eparams;

if (emig.sexDep) {
	if (emig.stgDep) {
		if (emig.densDep) {
			for (int i = 0; i < sstruct.nStages; i++) {
				for (int sex = 0; sex < NSEXES; sex++) {
					etraits.d0    = (float)SexEmigDensPar->Cells[1][2*i+sex+1].ToDouble();
					etraits.alpha = (float)SexEmigDensPar->Cells[2][2*i+sex+1].ToDouble();
					etraits.beta  = (float)SexEmigDensPar->Cells[3][2*i+sex+1].ToDouble();
					pSpecies->setEmigTraits(i,sex,etraits);
				}
			}
		}
		else {
			for (int i = 0; i < sstruct.nStages; i++) {
				for (int sex = 0; sex < NSEXES; sex++) {
					etraits.d0 = (float)SexEmigPar->Cells[1][2*i+sex+1].ToDouble();
					etraits.alpha = etraits.beta = 0.0;
					pSpecies->setEmigTraits(i,sex,etraits);
				}
			}
		}
	}
	else { // !emig.stgDep
		if (emig.indVar) {
			if (emig.densDep) {
				for (int sex = 0; sex < NSEXES; sex++) {
					eparams.d0Mean    = SexEmigDensPar->Cells[1][sex+1].ToDouble();
					eparams.d0SD      = SexEmigDensPar->Cells[2][sex+1].ToDouble();
					eparams.alphaMean = SexEmigDensPar->Cells[3][sex+1].ToDouble();
					eparams.alphaSD   = SexEmigDensPar->Cells[4][sex+1].ToDouble();
					eparams.betaMean  = SexEmigDensPar->Cells[5][sex+1].ToDouble();
					eparams.betaSD    = SexEmigDensPar->Cells[6][sex+1].ToDouble();
					pSpecies->setEmigParams(0,sex,eparams);
				}
			}
			else {
				for (int sex = 0; sex < NSEXES; sex++) {
					eparams.d0Mean = SexEmigPar->Cells[1][sex+1].ToDouble();
					eparams.d0SD   = SexEmigPar->Cells[2][sex+1].ToDouble();
					eparams.alphaMean = eparams.betaMean = 0.0;
					eparams.alphaSD = eparams.betaSD = 0.000001;
					pSpecies->setEmigParams(0,sex,eparams);
				}
			}
		}
		else { // !emig.indVar
			if (emig.densDep) {
				for (int sex = 0; sex < NSEXES; sex++) {
					etraits.d0    = (float)SexEmigDensPar->Cells[1][sex+1].ToDouble();
					etraits.alpha = (float)SexEmigDensPar->Cells[2][sex+1].ToDouble();
					etraits.beta  = (float)SexEmigDensPar->Cells[3][sex+1].ToDouble();
					pSpecies->setEmigTraits(0,sex,etraits);
				}
			}
			else {
				for (int sex = 0; sex < NSEXES; sex++) {
					etraits.d0 = (float)SexEmigPar->Cells[1][sex+1].ToDouble();
					etraits.alpha = etraits.beta = 0.0;
					pSpecies->setEmigTraits(0,sex,etraits);
				}
			}
		}
	}
}
else { // !emig.sexDep
	if (emig.stgDep) {
		if (emig.densDep) {
			for (int i = 0; i < sstruct.nStages; i++) {
				etraits.d0    = (float)SexEmigDensPar->Cells[1][i+1].ToDouble();
				etraits.alpha = (float)SexEmigDensPar->Cells[2][i+1].ToDouble();
				etraits.beta  = (float)SexEmigDensPar->Cells[3][i+1].ToDouble();
				pSpecies->setEmigTraits(i,0,etraits);
			}
		}
		else {
			for (int i = 0; i < sstruct.nStages; i++) {
				etraits.d0 = (float)SexEmigPar->Cells[1][i+1].ToDouble();
				etraits.alpha = etraits.beta = 0.0;
				pSpecies->setEmigTraits(i,0,etraits);
			}
		}
	}
	else { // !emig.stgDep
		if (emig.indVar) {
			if (emig.densDep) {
				eparams.d0Mean    = StrToFloat(edtD0Mean->Text);
				eparams.d0SD      = StrToFloat(edtD0SD->Text);
				eparams.alphaMean = StrToFloat(edtAlphaMean->Text);
				eparams.alphaSD   = StrToFloat(edtAlphaSD->Text);
				eparams.betaMean  = StrToFloat(edtBetaMean->Text);
				eparams.betaSD    = StrToFloat(edtBetaSD->Text);
			}
			else {
				eparams.d0Mean = StrToFloat(edtEPmean->Text);
				eparams.d0SD   = StrToFloat(edtEPsd->Text);
				eparams.alphaMean = eparams.betaMean = 0.0;
				eparams.alphaSD = eparams.betaSD = 0.000001;
			}
			pSpecies->setEmigParams(0,0,eparams);
		}
		else { // !emig.indVar
			if (emig.densDep) {
				etraits.d0 = StrToFloat(edtD0->Text);
				etraits.alpha = StrToFloat(edtAlpha->Text);
				etraits.beta  = StrToFloat(edtBeta->Text);
			}
			else {
				etraits.d0 = StrToFloat(frmSpecies->edtEP->Text);
				etraits.alpha = etraits.beta = 0.0;
			}

			// TEMPORARY FIX TO ENSURE EMIGRATION PROBABILITY LIES WITHIN BOUNDS
			// NEEDS PROPER CHECK ON VALUE ENTERED INTO THIS FIELD
			// SIMILAR CHECKS ARE REQUIRED FOR SEX- AND STAGE-DEPENDENT PARAMETERS
			if (etraits.d0 < 0.0) etraits.d0 = 0.0;
			if (etraits.d0 > 1.0) etraits.d0 = 1.0;

			pSpecies->setEmigTraits(0,0,etraits);
		}
	}
}

// Transfer

trfrRules trfr = pSpecies->getTrfr();
trfrScales scale = pSpecies->getTrfrScales();

if (RGMovements->ItemIndex == 0) trfr.moveModel = false;
else trfr.moveModel = true;

// NOTE: movement model parameters are handled by FormMove
if (!trfr.moveModel) { // dispersal kernels
	if (RGKernel->ItemIndex == 0) trfr.twinKern = false;
	else trfr.twinKern = true;
	if (CBIndVarKernel->Checked) trfr.indVar = true;
	else trfr.indVar = false;
	if (trfr.indVar) {
		scale.dist1Scale  = StrToFloat(edtDist1Scale->Text);
		scale.dist2Scale  = StrToFloat(edtDist2Scale->Text);
		scale.PKern1Scale = StrToFloat(edtPKern1Scale->Text);
		pSpecies->setTrfrScales(scale);
	}

	if (CBSexKernels->Checked) trfr.sexDep = true; else trfr.sexDep = false;
	if (CBStageKernels->Checked) trfr.stgDep = true; else trfr.stgDep = false;

	trfrKernTraits k;
	trfrKernParams kparams;
	float minkernelmean; // minimum value permitted for kernel means
	if (pSpecies->useFullKernel()) minkernelmean = 1.0;
	else minkernelmean = (float)ppLand.resol;

	if (trfr.sexDep) {
		if (trfr.stgDep) {
			for (int i = 0; i < sstruct.nStages; i++) {
				for (int sex = 0; sex < NSEXES; sex++) {
					k.meanDist1 = (float)SexKernPar->Cells[1][2*i+sex+1].ToDouble();
					if (trfr.twinKern) {
						k.meanDist2 = (float)SexKernPar->Cells[2][2*i+sex+1].ToDouble();
						k.probKern1 = (float)SexKernPar->Cells[3][2*i+sex+1].ToDouble();
					}
					else {
						k.meanDist2 = k.meanDist1; k.probKern1 = 1.0;
					}
					pSpecies->setKernTraits(i,sex,k,minkernelmean);
				}
			}
		}
		else { // !trfr.stgDep
			if (trfr.indVar) {
				for (int sex = 0; sex < NSEXES; sex++) {
					kparams.dist1Mean = (float)SexKernPar->Cells[1][sex+1].ToDouble();
					kparams.dist1SD   = (float)SexKernPar->Cells[2][sex+1].ToDouble();
					if (trfr.twinKern) {
						kparams.dist2Mean  = (float)SexKernPar->Cells[3][sex+1].ToDouble();
						kparams.dist2SD    = (float)SexKernPar->Cells[4][sex+1].ToDouble();
						kparams.PKern1Mean = (float)SexKernPar->Cells[5][sex+1].ToDouble();
						kparams.PKern1SD   = (float)SexKernPar->Cells[6][sex+1].ToDouble();
					}
					else {
						kparams.dist2Mean = kparams.dist1Mean; kparams.dist2SD = kparams.dist1SD;
						kparams.PKern1Mean = 0.999; kparams.PKern1SD = 0.001;
					}
					pSpecies->setKernParams(0,sex,kparams,minkernelmean);
				}
			}
			else {
				for (int sex = 0; sex < NSEXES; sex++) {
					k.meanDist1 = (float)SexKernPar->Cells[1][sex+1].ToDouble();
					if (trfr.twinKern) {
						k.meanDist2 = (float)SexKernPar->Cells[2][sex+1].ToDouble();
						k.probKern1 = (float)SexKernPar->Cells[3][sex+1].ToDouble();
					}
					else {
						k.meanDist2 = k.meanDist1; k.probKern1 = 1.0;
					}
					pSpecies->setKernTraits(0,sex,k,minkernelmean);
				}
			}
		}
	}
	else { // !trfr.sexDep
		if (trfr.stgDep) {
			for (int i = 0; i < sstruct.nStages; i++) {
				k.meanDist1 = (float)SexKernPar->Cells[1][i+1].ToDouble();
				if (trfr.twinKern) {
					k.meanDist2 = (float)SexKernPar->Cells[2][i+1].ToDouble();
					k.probKern1 = (float)SexKernPar->Cells[3][i+1].ToDouble();
				}
				else {
					k.meanDist2 = k.meanDist1; k.probKern1 = 1.0;
				}
				pSpecies->setKernTraits(i,0,k,minkernelmean);
			}
		}
		else { // !trfr.stgDep
			if (trfr.indVar) {
				kparams.dist1Mean  = StrToFloat(edtDist1Mean->Text);
				kparams.dist1SD    = StrToFloat(edtDist1SD->Text);
				kparams.dist1Scale = StrToFloat(edtDist1Scale->Text);
				if (trfr.twinKern) {
					kparams.dist2Mean   = StrToFloat(edtDist2Mean->Text);
					kparams.dist2SD     = StrToFloat(edtDist2SD->Text);
					kparams.dist2Scale  = StrToFloat(edtDist2Scale->Text);
					kparams.PKern1Mean  = StrToFloat(edtPKern1Mean->Text);
					kparams.PKern1SD    = StrToFloat(edtPKern1SD->Text);
					kparams.PKern1Scale = StrToFloat(edtPKern1Scale->Text);
				}
				else {
					kparams.dist2Mean  = kparams.dist1Mean; kparams.dist2SD = kparams.dist1SD;
					kparams.dist2Scale =  kparams.dist1Scale;
					kparams.PKern1Mean = 0.999; kparams.PKern1SD = 0.001;
					kparams.PKern1Scale = 0.001;
				}
//				else MAXDist = kparams.maxDist1;
				pSpecies->setKernParams(0,0,kparams,minkernelmean);
			}
			else {
				k.meanDist1 = StrToFloat(edtDist1->Text);
				if (trfr.twinKern) {
					k.meanDist2 = StrToFloat(edtDist2->Text);
					k.probKern1 = StrToFloat(edtPKern1->Text);
				}
				else {
					k.meanDist2 = k.meanDist1; k.probKern1 = 1.0;
				}
				pSpecies->setKernTraits(0,0,k,minkernelmean);
			}
		}
	}

	// Mortality
	if (RGMortality->ItemIndex == 0) trfr.distMort = false;
	else trfr.distMort = true;
	trfrMortParams mort;
	mort.fixedMort = StrToFloat(edtMortProb->Text);
	mort.mortAlpha = StrToFloat(edtMortSlope->Text);
	mort.mortBeta  = StrToFloat(edtMortInfl->Text);
	pSpecies->setMortParams(mort);
} // end of dispersal kernels

pSpecies->setTrfr(trfr);

// Settlement

settleType sett;
settleRules srules;
settleSteps ssteps;
settleTraits settleDD;
settParams sparams;

if (trfr.moveModel) {
	sett.indVar = CBIndVarSettle->Checked;
	if (CBSexSettMovt->Checked) {
		if (CBStageSettMovt->Checked) {
			sett.sexDep = true; sett.stgDep = true;
		}
		else {
			sett.sexDep = true; sett.stgDep = false;
		}
	}
	else {
		if (CBSexSettMovt->Checked) {
			sett.sexDep = false; sett.stgDep = true;
		}
		else {
			sett.sexDep = false; sett.stgDep = false;
		}
	}
	pSpecies->setSettle(sett);

	ssteps.minSteps = StrToInt(edtMinSteps->Text);
	// max steps must be set to very large number for per-step mortality only
	if (RGSteps->ItemIndex == 0) ssteps.maxSteps = StrToInt(edtNsteps->Text);
	else ssteps.maxSteps = 99999999;
	// max steps per year should be set for all stage-structured populations
//	if (dem.stageStruct){
//		ssteps.maxStepsYr = StrToInt(edtMaxStepYear->Text);
//	}
//	else ssteps.maxStepsYr = 99999999;
	pSpecies->setSteps(0,0,ssteps);

//	int settleType;
	int nstages,nsexes;
	if (dem.stageStruct) nstages = sstruct.nStages; else nstages = 1;
	if (dem.repType == 0) nsexes = 1; else nsexes = 2;
	for (int stg = 0; stg < nstages; stg++) {
		for (int sex = 0; sex < nsexes; sex++) {
			srules = pSpecies->getSettRules(stg,sex);
			srules.wait = srules.go2nbrLocn = false;
			srules.densDep = CBDensDepSettle->Checked;
			if (CBSexSettMovt->Checked) { // sex-dependent
				if (CBStageSettMovt->Checked) { // stage-dependent
					srules.findMate   = SexSettlePar->Cells[1][2*stg+sex+1].ToInt();
					ssteps.maxStepsYr = SexSettlePar->Cells[2][2*stg+sex+1].ToInt();
					if (srules.densDep && !sett.indVar) {
						settleDD.s0       = SexSettlePar->Cells[3][2*stg+sex+1].ToDouble();
						settleDD.alpha    = SexSettlePar->Cells[4][2*stg+sex+1].ToDouble();
						settleDD.beta     = SexSettlePar->Cells[5][2*stg+sex+1].ToDouble();
					}
				}
				else { // not stage-dependent
					srules.findMate   = SexSettlePar->Cells[1][sex+1].ToInt();
					ssteps.maxStepsYr = SexSettlePar->Cells[2][sex+1].ToInt();
					if (srules.densDep && !sett.indVar) {
						settleDD.s0       = SexSettlePar->Cells[3][sex+1].ToDouble();
						settleDD.alpha    = SexSettlePar->Cells[4][sex+1].ToDouble();
						settleDD.beta     = SexSettlePar->Cells[5][sex+1].ToDouble();
					}
				}
			}
			else { // not sex-dependent
				if (CBStageSettMovt->Checked) { // stage-dependent
					srules.findMate   = SexSettlePar->Cells[1][stg+1].ToInt();
					ssteps.maxStepsYr = SexSettlePar->Cells[2][stg+1].ToInt();
					if (srules.densDep && !sett.indVar) {
						settleDD.s0       = SexSettlePar->Cells[3][stg+1].ToDouble();
						settleDD.alpha    = SexSettlePar->Cells[4][stg+1].ToDouble();
						settleDD.beta     = SexSettlePar->Cells[5][stg+1].ToDouble();
					}
				}
				else { // not stage-dependent - parameters from panel 2
					srules.findMate   = CBFindMate->Checked;
					ssteps.maxStepsYr = StrToInt(edtMaxStepYear->Text);
					if (srules.densDep && !sett.indVar) {
						settleDD.s0       = StrToFloat(edtSettS0->Text);
						settleDD.alpha    = StrToFloat(edtSettAlpha->Text);
						settleDD.beta     = StrToFloat(edtSettBeta->Text);
					}
				}
			}
			pSpecies->setSettRules(stg,sex,srules);
			pSpecies->setSteps(stg,sex,ssteps);
			if (sett.indVar) {
				if (stg == 0) {
					// initial DD parameters are always from panel 2
					// and applied to both sexes in the GUI version
					sparams.s0Mean      = StrToFloat(edtSettS0Mean->Text);
					sparams.s0SD        = StrToFloat(edtSettS0SD->Text);
					sparams.s0Scale     = StrToFloat(edtSettS0Scale->Text);
					sparams.alphaSMean  = StrToFloat(edtSettAlphaMean->Text);
					sparams.alphaSSD    = StrToFloat(edtSettAlphaSD->Text);
					sparams.alphaSScale = StrToFloat(edtSettAlphaScale->Text);
					sparams.betaSMean   = StrToFloat(edtSettBetaMean->Text);
					sparams.betaSSD     = StrToFloat(edtSettBetaSD->Text);
					sparams.betaSScale  = StrToFloat(edtSettBetaScale->Text);
					pSpecies->setSettParams(stg,sex,sparams);
				}
			}
			else {
				pSpecies->setSettTraits(stg,sex,settleDD);
			}
		} // end of for (int sex = 0; sex < nsexes; sex++)
	} // end of for (int stg = 0; stg < nstages; stg++)
}
else { // dispersal kernels
	sett.indVar = false;
	if (CBSexSettKern->Checked) {
		if (CBStageSettKern->Checked) {
			sett.sexDep = true; sett.stgDep = true;
		}
		else {
			sett.sexDep = true; sett.stgDep = false;
		}
	}
	else {
		if (CBStageSettKern->Checked) {
			sett.sexDep = false; sett.stgDep = true;
		}
		else {
			sett.sexDep = false; sett.stgDep = false;
		}
	}
	pSpecies->setSettle(sett);

	int settKernType,findMate,nstages,nsexes;
	if (dem.stageStruct) nstages = sstruct.nStages; else nstages = 1;
	if (dem.repType == 0) nsexes = 1; else nsexes = 2;
	for (int stg = 0; stg < nstages; stg++) {
		for (int sex = 0; sex < nsexes; sex++) {
			srules = pSpecies->getSettRules(stg,sex);
			srules.densDep = false;
			if (CBSexSettKern->Checked) { // sex-structured
				if (CBStageSettKern->Checked) { // stage-structured
					settKernType = SexSettlePar->Cells[1][2*stg+sex+1].ToInt();
					findMate = SexSettlePar->Cells[2][2*stg+sex+1].ToInt();
				}
				else { // not stage-structured
					settKernType = SexSettlePar->Cells[1][sex+1].ToInt();
					findMate = SexSettlePar->Cells[2][sex+1].ToInt();
				}
			}
			else { // not sex-structured
				if (CBStageSettKern->Checked) { // stage-structured
					settKernType = SexSettlePar->Cells[1][stg+1].ToInt();
					findMate = SexSettlePar->Cells[2][stg+1].ToInt();
				}
				else { // not stage-structured - parameters from panel 2
					settKernType = RGSettKern->ItemIndex;
					findMate = CBSettKernMate->Checked;
				}
			}
			switch (settKernType) {
			case 0:
				srules.wait = false; srules.go2nbrLocn = false;
				break;
			case 1:
				srules.wait = true; srules.go2nbrLocn = false;
				break;
			case 2:
				srules.wait = false; srules.go2nbrLocn = true;
				break;
			case 3:
				srules.wait = true; srules.go2nbrLocn = true;
				break;
			}
			if (findMate == 0) srules.findMate = false; else srules.findMate = true;
			pSpecies->setSettRules(stg,sex,srules);
		}
	}
}

pSpecies->setTraits();

}

//---------------------------------------------------------------------------
void __fastcall TfrmSpecies::BtnCancelClick(TObject *Sender)
{
formCancelled = true;
Close();
}


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

