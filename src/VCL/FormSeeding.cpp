//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include <stdlib.h>
#include <fstream>
#include <io.h>
#include <iostream>
#include <sstream>
#include <string>
#include <string.h>
using namespace std;

#include "FormSeeding.h"
#include "FormSpecies.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TfrmSeeding *frmSeeding;
//---------------------------------------------------------------------------
__fastcall TfrmSeeding::TfrmSeeding(TComponent* Owner)
	: TForm(Owner)
{
Left = 100; Top = 60;
}
//---------------------------------------------------------------------------
void __fastcall TfrmSeeding::NindsChange(TObject *Sender)
{
if (cboNinds->ItemIndex == 2) {
	edtNinds->Enabled = true; edtNinds->Visible = true;
}
else {
	edtNinds->Enabled = false; edtNinds->Visible = false;
}
}
//---------------------------------------------------------------------------
void __fastcall TfrmSeeding::RGCellsClick(TObject *Sender)
{
locn loc;

switch (RGCells->ItemIndex) {
case 0:
	edtNSpCells->Enabled = false;
	LabelListBox->Enabled = false;   LabelListBox->Visible = false;
	LBSpPres->Enabled = false; LBSpPres->Visible = false;
	PanelPlus->Enabled = true; PanelPlus->Visible = true;
	CBAdditional->Checked = false;
	BtnOK->Enabled = true;
	break;
case 1:
	edtNSpCells->Enabled = true;
	LabelListBox->Enabled = false;   LabelListBox->Visible = false;
	LBSpPres->Enabled = false; LBSpPres->Visible = false;
	PanelPlus->Enabled = false; PanelPlus->Visible = false;
	CBAdditional->Checked = false;
	BtnOK->Enabled = true;
	break;
case 2:
	edtNSpCells->Enabled = false;
	LabelListBox->Enabled = true;   LabelListBox->Visible = true;
	LBSpPres->Enabled = true; LBSpPres->Visible = true;
	PanelPlus->Enabled = true; PanelPlus->Visible = true;
	CBAdditional->Checked = false;
	BtnOK->Enabled = true;
	if (LBSpPres->Items->Count == 0) {
		string list;
		if (landscapeLoaded) {
			int nDistns = pLandscape->distnCount();
			// NOTE: CHANGE FOR MULTIPLE SPECIES ...
			if (nDistns > 0) {
				int nCells = pLandscape->distCellCount(0);
				for (int i = 0; i < nCells; i++) {
					loc = pLandscape->getDistnCell(0,i);
					list = Int2Str(loc.x) + " " + Int2Str(loc.y);
					LBSpPres->Items->Add(list.c_str());
				}
			}
		}
	}
	break;
}
}
//---------------------------------------------------------------------------
void __fastcall TfrmSeeding::RGRulesClick(TObject *Sender)
{
landParams ppLand = pLandscape->getLandParams();

switch (RGRules->ItemIndex) {
case 0:
	PanelFree->Enabled = true;  PanelFree->Visible = true;
	PanelPlus->Enabled = false; PanelPlus->Visible = false;
	edtTotRandomCells->Enabled = true; edtTotRandomCells->Visible = true;
	break;
case 1:
	PanelFree->Enabled = true;  PanelFree->Visible = true;
	PanelPlus->Enabled = false; PanelPlus->Visible = false;
	edtTotRandomCells->Enabled = false; edtTotRandomCells->Visible = false;
	break;
case 2:
	if (ppLand.generated) {
		MessageDlg("Manual selection of cells is not permitted for a generated landscape",
			mtError, TMsgDlgButtons() << mbOK,0);
		RGRules->ItemIndex = 1;
		PanelFree->Enabled = true;	PanelFree->Visible = true;
		PanelPlus->Enabled = false;	PanelPlus->Visible = false;
	}
	else {
		PanelFree->Enabled = false;	PanelFree->Visible = false;
		PanelPlus->Enabled = true; 	PanelPlus->Visible = true;
	}
	break;
}
}



//------------------------------------------------------------------------------
void __fastcall TfrmSeeding::BtnAddClick(TObject *Sender)
{
landParams ppLand = pLandscape->getLandParams();

if (ppLand.patchModel) {
	int id = StrToInt(frmSeeding->edtPatchID->Text);
	if (id == 0) {
		MessageDlg("PatchID = 0 corresponds to no patch",
				mtError, TMsgDlgButtons() << mbOK,0);
		return;
	}
	if (id < 0) {
		MessageDlg("PatchID must be a positive integer",
				mtError, TMsgDlgButtons() << mbOK,0);
		return;
	}
	InitialList->Lines->Add((frmSeeding->edtPatchID->Text).c_str());
	pLandscape->addInitCell(id,-1);
	edtPatchID->Clear();
}
else {
	int x = StrToInt(frmSeeding->edtAddX->Text);
	int y = StrToInt(frmSeeding->edtAddY->Text);
	if (x < 0 || x > ppLand.maxX || y < 0 || y > ppLand.maxY) {
		MessageDlg("Invalid co-ordinates",
				mtError, TMsgDlgButtons() << mbOK,0);
		return;
	}
	InitialList->Lines->Add((frmSeeding->edtAddX->Text
		+ "  " + frmSeeding->edtAddY->Text).c_str());
	pLandscape->addInitCell(x,y);
	edtAddX->Clear(); frmSeeding->edtAddY->Clear();
}
}


//---------------------------------------------------------------------------
void __fastcall TfrmSeeding::BtnCancelClick(TObject *Sender)
{
frmSeeding->Close();
}

//------------------------------------------------------------------------------
void __fastcall TfrmSeeding::Refresh1Click(TObject *Sender)
{
setDefaults();
}

//---------------------------------------------------------------------------
void __fastcall TfrmSeeding::initialise(void)
{
UnicodeString units,spprescells,manually;
landParams ppLand = pLandscape->getLandParams();
demogrParams dem = pSpecies->getDemogr();

if (forceInit) setDefaults();

if (ppLand.generated) {
	RG1->ItemIndex = 0;
	RG1->Enabled = false;
	units = "cells";
}
else { // Imported landscape
	RG1->Enabled = true;
	spprescells = "species presence cells";
	manually = "Manually include ";
	if (ppLand.patchModel) {
		units = "patches";
		init_unit->Caption = "ha";
		CBAdditional->Caption = manually + units;
		edtAddX->Visible = false; edtAddY->Visible = false;
		edtPatchID->Visible = true;
	}
	else {
		units = "cells";
		init_unit->Caption = "cell";
		CBAdditional->Caption = manually + "landscape " + units;
		edtAddX->Visible = true; edtAddY->Visible = true;
		edtPatchID->Visible = false;
	}
	LabelLandRes->Visible = false;
	edtAddX->Enabled = false; edtAddY->Enabled = false;
	edtPatchID->Enabled = false;
	CBAdditional->Checked = false;
	BtnAdd->Enabled = false;
}

if (dem.stageStruct) {
	LabelPropn->Visible = true;
	SGinitialStage->Visible = true;
	SGinitialStage->Enabled = true;
	RGinitAges->Visible = true;
	RGinitAges->Enabled = true;
}
else {
	LabelPropn->Visible = false;
	SGinitialStage->Visible = false;
	SGinitialStage->Enabled = false;
	RGinitAges->Visible = false;
	RGinitAges->Enabled = false;
}

RGRules->Items->Strings[0] = "Random (given no. of " + units + ")";
RGRules->Items->Strings[1] = "All suitable " + units;
RGRules->Items->Strings[2] = "Manually select " + units;
edtTotRandomCells->EditLabel->Caption = "No. of randomly selected " + units + ": ";
RGCells->Items->Strings[0] = "All " + units + " within all " + spprescells;
RGCells->Items->Strings[1] = "All " + units + " within some " + spprescells + " (randomly chosen) ";
RGCells->Items->Strings[2] = "All " + units + " within selected " + spprescells;

}

//---------------------------------------------------------------------------
void __fastcall TfrmSeeding::setDefaults(void)
{
landParams ppLand = pLandscape->getLandParams();
demogrParams dem = pSpecies->getDemogr();

RG1->ItemIndex = 0;
RGRules->ItemIndex = 1;
RGRules->Enabled = true; RGRules->Visible = true;
PanelFree->Enabled = true; PanelFree->Visible = true;
RGCells->ItemIndex = 0;
RGCells->Enabled = false; RGCells->Visible = false;
PanelSpDist->Enabled = false; PanelSpDist->Visible = false;
BtnSaveInitial->Enabled = false; BtnSaveInitial->Visible = false;
PanelPlus->Enabled = false; PanelPlus->Visible = false;

edtMinX->Text = "0"; edtMaxX->Text = IntToStr(ppLand.maxX);
edtMinY->Text = "0"; edtMaxY->Text = IntToStr(ppLand.maxY);
edtNinds->Text = "0";
edtNSpCells->Text = "0";
edtTotRandomCells->Text = "0";
edtInitFreezeYear->Text = "0";
CBrestrictRange->Checked = false;
edtRestrictRows->Text = "100";
edtRestrictRows->Enabled = false; edtRestrictRows->Visible = false;
edtRestrictFreq->Text = "10";
edtRestrictFreq->Enabled = false; edtRestrictFreq->Visible = false;
edtFinalFreezeYear->Text = "0";
edtFinalFreezeYear->Enabled = false; edtFinalFreezeYear->Visible = false;

LBSpPres->ClearSelection();
CBAdditional->Checked = false;
InitialList->Clear();
pLandscape->clearInitCells();

if (ppLand.generated){
	BtnSaveInitial->Enabled = false; BtnSaveInitial->Visible = false;
}
else {
	BtnSaveInitial->Enabled = true;  BtnSaveInitial->Visible = true;
}

if (dem.stageStruct) {
	stageParams sstruct = pSpecies->getStage();
	frmSeeding->SGinitialStage->ColCount = sstruct.nStages;
	frmSeeding->SGinitialStage->Cells[0][0] = "Stage";
	frmSeeding->SGinitialStage->Cells[0][1] = "Propn.";
	for (int i = 1; i < sstruct.nStages+1; i++) {
		frmSeeding->SGinitialStage->Cells[i][0] = i;
		if (sstruct.nStages == 2) frmSeeding->SGinitialStage->Cells[i][1] = "1.0";
		else frmSeeding->SGinitialStage->Cells[i][1] = "0.0";
	}
}

BtnOK->Enabled = true;

}

//---------------------------------------------------------------------------
void __fastcall TfrmSeeding::CBrestrictRangeClick(TObject *Sender)
{
if (CBrestrictRange->Checked) {
	edtRestrictRows->Enabled = true;     edtRestrictRows->Visible = true;
	edtRestrictFreq->Enabled = true;     edtRestrictFreq->Visible = true;
	edtFinalFreezeYear->Enabled = true;  edtFinalFreezeYear->Visible = true;
}
else {
	edtRestrictRows->Enabled = false;    edtRestrictRows->Visible = false;
	edtRestrictFreq->Enabled = false;    edtRestrictFreq->Visible = false;
	edtFinalFreezeYear->Enabled = false; edtFinalFreezeYear->Visible = false;

}
}

//---------------------------------------------------------------------------
void __fastcall TfrmSeeding::CBAdditionalClick(TObject *Sender)
{
landParams ppLand = pLandscape->getLandParams();

if (CBAdditional->Checked) {
	if (ppLand.patchModel) {
		edtAddX->Enabled = false; edtAddY->Enabled = false;
		edtPatchID->Enabled = true;
		LabelLandRes->Enabled = false;
	}
	else{
		edtAddX->Enabled = true; edtAddY->Enabled = true;
		edtPatchID->Enabled = false;
		LabelLandRes->Enabled = true;
	}
	BtnAdd->Enabled = true;
}
else {
	edtAddX->Enabled = false; edtAddY->Enabled = false;
	edtPatchID->Enabled = false;
	LabelLandRes->Enabled = false;
	BtnAdd->Enabled = false;
}
}

//------------------------------------------------------------------------------
void __fastcall TfrmSeeding::RG1Click(TObject *Sender)
{
//simParams sim = paramsSim->getSim();
landParams ppLand = pLandscape->getLandParams();
switch (RG1->ItemIndex) {
case 0: // free initialisaton
	RGRules->Enabled = true; RGRules->Visible = true;
	RGCells->Enabled = false; RGCells->Visible = false;
	PanelNInds->Enabled = true; PanelNInds->Visible = true;
	PanelSpDist->Enabled = false; PanelSpDist->Visible = false;
	BtnSaveInitial->Enabled = true; BtnSaveInitial->Visible = true;
	if (RGRules->ItemIndex == 2) { // manually select cells
		PanelFree->Enabled = false; PanelFree->Visible = false;
		PanelPlus->Enabled = true; PanelPlus->Visible = true;
		CBAdditional->Checked = true;
	}
	else {
		PanelFree->Enabled = true; PanelFree->Visible = true;
		PanelPlus->Enabled = false; PanelPlus->Visible = false;
		CBAdditional->Checked = false;
	}
	break;
case 1: // from species distribution
//	if (sim.initDistLoaded)
	if (ppLand.spDist)
	{
		RGRules->Enabled = false; RGRules->Visible = false;
		RGCells->Enabled = true; RGCells->Visible = true;
		PanelNInds->Enabled = true; PanelNInds->Visible = true;
		PanelFree->Enabled = false; PanelFree->Visible = false;
		PanelSpDist->Enabled = true; PanelSpDist->Visible = true;
		BtnSaveInitial->Enabled = true; BtnSaveInitial->Visible = true;
		PanelPlus->Enabled = true; PanelPlus->Visible = true;
		CBAdditional->Checked = false;
	}
	else {
		MessageDlg("The species distribution has not been loaded.",
			mtError, TMsgDlgButtons() << mbRetry,0);
		RG1->ItemIndex = 0;
	}
	break;
case 2: // from initial individuals file
//	break;
case 3: // from initialisation file
	RGRules->Enabled = false; RGRules->Visible = false;
	RGCells->Enabled = false; RGCells->Visible = false;
	PanelNInds->Enabled = false; PanelNInds->Visible = false;
	PanelFree->Enabled = false; PanelFree->Visible = false;
	PanelSpDist->Enabled = false; PanelSpDist->Visible = false;
	BtnSaveInitial->Enabled = false; BtnSaveInitial->Visible = false;
	PanelPlus->Enabled = false; PanelPlus->Visible = false;
	CBAdditional->Checked = false;
	break;
}
}


//------------------------------------------------------------------------------
void __fastcall TfrmSeeding::BtnSaveInitialClick(TObject *Sender)
{
string initname,initcellsname,initspdistcellsname;
int xx, yy;
landParams ppLand = pLandscape->getLandParams();
demogrParams dem = pSpecies->getDemogr();
stageParams sstruct = pSpecies->getStage();

ofstream initfile;

SaveTextFileDialog1->InitialDir = paramsSim->getDir(1).c_str();
if (SaveTextFileDialog1->Execute()) {
	StatusBar1->Panels->Items[0]->Text = "Saving Initialisation file. Please wait...";
	StatusBar1->Refresh();
	initname = AnsiString(SaveTextFileDialog1->FileName).c_str();
	// extract separate file and path names
	string fname = initname;
	string path = fname;
	unsigned int loc = fname.find("\\",0);
#if RSWIN64
	while(loc < 999999)
#else
	while(loc != string::npos)
#endif
	{
		fname = fname.substr(loc+1);
		loc =  fname.find("\\",0);
	}
	path = path.substr(0,path.length()-fname.length());

//string msg = path;
//MessageDlg(msg.c_str(),mtWarning, TMsgDlgButtons() << mbOK,0);

	initfile.open(initname.c_str());

	initfile << "SeedType\tFreeType\tSpType\tInitDens";
	if (ppLand.patchModel) initfile << "\tIndsHa";
	else  initfile << "\tIndsCell";
	initfile << "\tminX\tmaxX\tminY\tmaxY\tNCells\tNSpCells";
	initfile << "\tInitFreezeYear\tRestrictRows\tRestrictFreq\tFinalFreezeYear";
	initfile << "\tInitIndsFile\tInitAge";
	if (dem.stageStruct) {
		for (int i = 1; i < sstruct.nStages; i++) initfile << "\t" << "PropStage_" << i;
	}
	initfile << "\tInitCells_File\tInitSpDistCells_File" << endl;

	initfile << RG1->ItemIndex;
	if (RG1->ItemIndex == 0) initfile << "\t" << RGRules->ItemIndex << "\t-9";
	else initfile << "\t-9\t" << RGCells->ItemIndex;
	initfile << "\t" << cboNinds->ItemIndex;
	if (cboNinds->ItemIndex == 2) initfile << "\t" << StrToInt(edtNinds->Text);
	else initfile <<"\t-9";
	if (RG1->ItemIndex == 0) {
		initfile << "\t" << StrToInt(edtMinX->Text) << "\t" << StrToInt(edtMaxX->Text);
		initfile << "\t" << StrToInt(edtMinY->Text) << "\t" << StrToInt(edtMaxY->Text);
	}
	else
		initfile << "\t-9\t-9\t-9\t-9";
	if (RG1->ItemIndex == 0 && RGRules->ItemIndex == 0)
		initfile << "\t" << StrToInt(edtTotRandomCells->Text);
	else initfile << "\t-9";
	if (RG1->ItemIndex == 1) initfile << "\t" << StrToInt(edtNSpCells->Text);
	else initfile << "\t-9";
	if (RG1->ItemIndex == 0) initfile << "\t" << StrToInt(edtInitFreezeYear->Text);
	else initfile << "\t-9";

	if (CBrestrictRange->Checked) {
		initfile << "\t" << StrToInt(edtRestrictRows->Text)
			<< "\t" << StrToInt(edtRestrictFreq->Text) << "\t" << StrToInt(edtFinalFreezeYear->Text);
	}
	else initfile << "\t-9\t-9\t-9";

  // initial individuals file name is never saved, as option is not allowed
	initfile << "\t-9";    
	
	if (dem.stageStruct) {
		initfile <<"\t" << RGinitAges->ItemIndex;
		for (int i = 1; i < sstruct.nStages; i++)
			initfile <<"\t" << SGinitialStage->Cells[i][1].ToDouble();
	}
	else  initfile << "\t-9";

	if (pLandscape->initCellCount() > 0) { // additional cells selected
		initcellsname = fname.substr(0,fname.length()-4) + "_InitCells.txt";
		initfile <<"\t" << initcellsname;
	}
	else initfile <<"\t-9";
	if (RG1->ItemIndex == 1 && RGCells->ItemIndex == 2)
	{
		initspdistcellsname = fname.substr(0,fname.length()-4) + "_InitSpDistCells.txt";
		initfile <<"\t" << initspdistcellsname;
	}
	else initfile <<"\t-9";

	initfile.close(); initfile.clear();

	if (RG1->ItemIndex == 1) { // from species distribution
		if (RGCells->ItemIndex == 2) { // selected cells where the species is present
			initfile.open((path + initspdistcellsname).c_str());
			initfile << "x\ty\tpres" << endl;
			// TO BE REVISED FOR MULTIPLE SPECIES ...
			int ndistns = pLandscape->distnCount();
			if (ndistns > 0) {
				locn distloc;
				int ndistcells = (int)pLandscape->distCellCount(0);
				for (int i = 0; i < ndistcells; i++){
					distloc = pLandscape->getSelectedDistnCell(0,i);
					if (distloc.x >= 0) { // location is valid and selected
						initfile << distloc.x << "\t" << distloc.y << "\t1" << endl;
					}
				}
			}
			initfile.close(); initfile.clear();
		}
	}

	if (pLandscape->initCellCount() > 0) { // additional cells selected
		locn loc;
		initfile.open((path + initcellsname).c_str());
		if (ppLand.patchModel) initfile << "PatchID" << endl;
		else initfile << "x\ty" << endl;
		int ncells = pLandscape->initCellCount();
		for (int i = 0; i < ncells; i++) {
			loc = pLandscape->getInitCell(i);
			initfile << loc.x;
			if (!ppLand.patchModel) initfile << "\t" << loc.y;
			initfile << endl;
		}
		initfile.close(); initfile.clear();
	}

	StatusBar1->Panels->Items[0]->Text = "Initialisation file saved.";
	StatusBar1->Refresh();
}
InitialList->Clear();

}

//------------------------------------------------------------------------------
void __fastcall TfrmSeeding::BtnOKClick(TObject *Sender)
{
//int maxcells;
int nspcells,h;
float check,k;
Cell *pCell;
landParams ppLand = pLandscape->getLandParams();
demogrParams dem = pSpecies->getDemogr();
stageParams sstruct = pSpecies->getStage();
initParams init = paramsInit->getInit();
ifstream initfile,initcellsfile,initspdistcellsfile;
string initname,initcellsname,initspdistcellsname;

if (RG1->ItemIndex == -1) {
	MessageDlg("Please select an initialisation option to continue",
		mtError,TMsgDlgButtons() << mbOK,0);
	return;
}

init.seedType = RG1->ItemIndex;
if (init.seedType == 0) init.freeType = RGRules->ItemIndex;
if (init.seedType == 1) init.spDistType = RGCells->ItemIndex;

init.initDens = cboNinds->ItemIndex;
if (init.initDens == 2 && init.seedType < 2) {
	if (ppLand.patchModel) {
		init.indsHa = StrToFloat(edtNinds->Text);
		if (init.indsHa <= 0.0) {
			MessageDlg("Initial density must be greater than zero",
				mtError,TMsgDlgButtons() << mbOK,0);
			return;
		}
	}
	else {
		init.indsCell = StrToInt(edtNinds->Text);
		if (init.indsCell < 1) {
			MessageDlg("Initial individuals per cell must be 1 or more",
				mtError,TMsgDlgButtons() << mbOK,0);
			return;
		}
	}
}

if (dem.stageStruct && init.seedType < 2) {
	// check initialisation proportions
	float p;
	init.initAge = RGinitAges->ItemIndex;
	check = 0.0;
	for (int i = 1; i < sstruct.nStages; i++){
		p = frmSeeding->SGinitialStage->Cells[i][1].ToDouble();
		check += p;
	}
	if (check < 0.99999 || check > 1.00001) {
		MessageDlg("Please set the proportions of individuals per stage class"
				" to sum up to 1.0", mtError, TMsgDlgButtons() << mbOK,0);
		return;
	}
	else {
		for (int i = 1; i < sstruct.nStages; i++){
			p = frmSeeding->SGinitialStage->Cells[i][1].ToDouble();
			paramsInit->setProp(i,p);
		}
	}
}

init.initFrzYr = 0;
init.restrictRange = false;

switch (init.seedType) {

case 0: // free initialisation
	init.minSeedX = StrToInt(edtMinX->Text); init.maxSeedX = StrToInt(edtMaxX->Text);
	init.minSeedY = StrToInt(edtMinY->Text); init.maxSeedY = StrToInt(edtMaxY->Text);
	if (init.minSeedX < 0 || init.maxSeedX >= ppLand.dimX
	||  init.minSeedY < 0 || init.maxSeedY >= ppLand.dimY )  {
		MessageDlg("The initialisation limits exceed the landscape boundaries",
			mtError, TMsgDlgButtons() << mbOK,0);
		return;
	}
	if (init.minSeedX > init.maxSeedX || init.minSeedY > init.maxSeedY) {
		MessageDlg("Minimum limit exceeds maximum limit",
			mtError, TMsgDlgButtons() << mbOK,0);
		return;

	}
	if (init.freeType == 0) { // random cells
		init.nSeedPatches = StrToInt(edtTotRandomCells->Text);
		if (init.nSeedPatches < 1) {
			string randmsg = "Number of ";
			if (ppLand.patchModel) randmsg += "patches ";
			else randmsg += "cells ";
			randmsg += "to initialise must be 1 or more";
			MessageDlg(randmsg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
			return;
		}
		// determine no. of patches / cells within the specified limits
		int npatches = pLandscape->patchCount();
		patchData pch;
		patchLimits limits;
		limits.xMin = init.minSeedX; limits.xMax = init.maxSeedX;
		limits.yMin = init.minSeedY; limits.yMax = init.maxSeedY;
		if (ppLand.patchModel) {
			int maxpatches = 0;
//string aaa1 = ("npatches=" + Int2Str(npatches));
//MessageDlg(aaa1.c_str(), mtWarning, TMsgDlgButtons() << mbOK,0);
			for (int i = 0; i < npatches; i++) {
				pch = pLandscape->getPatchData(i);
//if (i < 3) {
//string aaa1 = ("i=" + Int2Str(i) + " patchNum=" + Int2Str(pch.patchNum)
//	+ " nCells=" + Int2Str(pch.nCells) + " x=" + Int2Str(pch.x) + " y=" + Int2Str(pch.y));
//MessageDlg(aaa1.c_str(), mtWarning, TMsgDlgButtons() << mbOK,0);
//}
				if (pch.pPatch->withinLimits(limits)) {
					if (pch.pPatch->getPatchNum() != 0) maxpatches++;
				}
			}
			if (init.nSeedPatches > maxpatches) {
				string initmsg = "You cannot initialise more than " + Int2Str(maxpatches)
					+ " patch(es) within the specified limits";
				MessageDlg(initmsg.c_str(), mtError, TMsgDlgButtons() << mbOK,0);
				return;
			}
		}
		else { // cell-based model
			// it is not possible to check max. no. of suitable cells, as they are not
			// identified until model is run (because they depend on K values for habitats)
			// however, if more are specified that are available, all are initialised
		}
	} // end of FreeType == 0
	init.initFrzYr = StrToInt(edtInitFreezeYear->Text);
	if (init.initFrzYr < 0) {
		MessageDlg("Freeze initial range year must be zero or greater",
			mtError, TMsgDlgButtons() << mbOK,0);
		return;
	}
	if (CBrestrictRange->Checked) {
		init.restrictRange = true;
		init.restrictRows = StrToInt(edtRestrictRows->Text);
		if (init.restrictRows < 1 || init.restrictRows >= ppLand.dimY) {
			MessageDlg("No. of rows must be > 0 and < Y dimension",
				mtError, TMsgDlgButtons() << mbOK,0);
			return;
		}
		init.restrictFreq = StrToInt(edtRestrictFreq->Text);
		if (init.restrictFreq < 1) {
			MessageDlg("Frequency must be > 0",
				mtError, TMsgDlgButtons() << mbOK,0);
			return;
		}
		init.finalFrzYr = StrToInt(edtFinalFreezeYear->Text);
		if (init.finalFrzYr <= init.initFrzYr) {
			MessageDlg("Freeze range year must be > initial freeze year",
				mtError, TMsgDlgButtons() << mbOK,0);
			return;
		}
	}
	else {
		init.restrictRange = false;
	}
	if (init.freeType == 2) { // manually selected cells/patches
		// at least one patch/cell must have been selected
		if (pLandscape->initCellCount() < 1) {
			string unittype;
			if (ppLand.patchModel) unittype = "patch"; else unittype = "cell";
			string label = ("At least one " + unittype + " must be selected manually");
			MessageDlg(label.c_str(), mtError, TMsgDlgButtons() << mbOK,0);
			return;
		}
	} // end of FreeType == 2
	break;

case 1: // from species distribution
	nspcells = pLandscape->distCellCount(0);
	if (init.spDistType == 1) { // random initialisation
		int nchosen = StrToInt(edtNSpCells->Text);
		if (nchosen < 1 || nchosen >= nspcells) {
			string msg = "The chosen nr. of cells must be in the range 1 to "
				+ Int2Str(nspcells-1);
			MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
			return;
		}
		init.nSpDistPatches = nchosen;
	}
	if (init.spDistType == 2) { // selected dsitribution cells
		int totselected = 0;
		// identify cells to be initialised
		for (int i = 0; i < nspcells; i++) {
			if (frmSeeding->LBSpPres->Selected[i]) { // mark cell to be initialised
				pLandscape->setDistnCell(0,i,true);
				totselected++;
			}
			else { // entry was de-selected - de-activate the cell
				pLandscape->setDistnCell(0,i,false);
			}
		}
		if (totselected == 0) {
			MessageDlg("At least one species distribution cell must be selected",
				mtError, TMsgDlgButtons() << mbOK,0);
			return;
		}
	}
	break;

case 2: // from initial individuals file
	OpenInitialisationFile->InitialDir = paramsSim->getDir(1).c_str();
	OpenInitialisationFile->Title = "Select the Initial Individuals File";
	if (OpenInitialisationFile->Execute()){
		initname = (AnsiString(OpenInitialisationFile->FileName).c_str());
		// check that selected file is the correct format
		int errcode = CheckInitIndsFile(initname);
		if (errcode < 0) {
			return;
		}
		else {
			int ninds = ReadInitIndsFile(0,pLandscape,initname);
#if RSDEBUG
DebugGUI(("TfrmSeeding::BtnOKClick): initname=" + initname
	+ " ninds=" + Int2Str(ninds)
	).c_str());
#endif
			if (ninds < 1) {
				MessageDlg("Problem encountered reading selected file",
					mtError, TMsgDlgButtons() << mbOK,0);
				return;
			}
			// extract file name from path+file
			string fname = initname;
			unsigned int loc = fname.find("\\",0);
#if RSWIN64
			while(loc < 999999)
#else
			while(loc != string::npos)
#endif
			{
				fname = fname.substr(loc+1);
				loc =  fname.find("\\",0);
			}
			init.indsFile = fname;
		}
	}
	else return;
	break;

case 3: // from initialisation file
	OpenInitialisationFile->InitialDir = paramsSim->getDir(1).c_str();
	OpenInitialisationFile->Title = "Select the main Initialisation File";
	if (OpenInitialisationFile->Execute()){
		bool hdrError = false;
		initname = (AnsiString(OpenInitialisationFile->FileName).c_str());
		// check that selected file is the correct format
		initfile.open(initname.c_str());
		string h0,h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12,h13,h14,h15,h16,h17,h18,h19;
		h0 = h1 = h2 = h3 = h4 = h5 = h6 = h7 = h8 = h9 = "x";
		h10 = h11 = h12 = h13 = h14 = h15 = h16 = h17 = h18 = h19 = "x";
		initfile >> h0 >> h1 >> h2 >> h3 >> h4 >> h5 >> h6 >> h7 >> h8 >> h9 >> h10;                
		initfile >> h11 >> h12 >> h13 >> h14 >> h15 >> h16;
		if (dem.stageStruct) {
			for (int i = 1; i < sstruct.nStages; i++) initfile >> h17;
		}
		initfile >> h18 >> h19;
		if (h0 != "SeedType" || h1 != "FreeType" || h2 != "SpType" || h3 != "InitDens")
			hdrError = true;
		if (ppLand.patchModel) {
			if (h4 != "IndsHa") hdrError = true;
		}
		else {
			if (h4 != "IndsCell") hdrError = true;
		}                      
		if (h5 != "minX" || h6 != "maxX" || h7 != "minY" || h8 != "maxY"
		||  h9 != "NCells" ||  h10 != "NSpCells" || h11 != "InitFreezeYear" 
		|| h12 != "RestrictRows" || h13 != "RestrictFreq" || h14 != "FinalFreezeYear" 
		|| h15 != "InitIndsFile" || h16 != "InitAge" 
		|| h18 != "InitCells_File" || h19 != "InitSpDistCells_File" ) hdrError = true;      
		if (hdrError) {
			initfile.close(); initfile.clear();
			MessageDlg("Format error in header line of selected file",
				mtError, TMsgDlgButtons() << mbRetry,0);
			return;
		}
		// read the initialisation data from the single data line
		int iiii;
		float ffff;
		string tttt;
		initfile >> init.seedType >> init.freeType >> init.spDistType >> init.initDens;
		if (ppLand.patchModel) initfile >> init.indsHa;
		else initfile >> init.indsCell;
		initfile >> init.minSeedX >> init.maxSeedX >> init.minSeedY >> init.maxSeedY;
		initfile >> init.nSeedPatches >> init.nSpDistPatches >> init.initFrzYr >> init.restrictRows;		
		initfile >> init.restrictFreq >> init.finalFrzYr >> tttt >> init.initAge; 
		if (init.restrictRows > 0) init.restrictRange = true; else init.restrictRange = false;
		if (dem.stageStruct) {
			for (int i = 1; i < sstruct.nStages; i++) {
				initfile >> ffff;
				paramsInit->setProp(i,ffff);
			}
		}
		initfile >> initcellsname >> initspdistcellsname;
//string msg = "seedtype=" + Int2Str(init.seedType) + " =" + Int2Str(init.freeType) + " =" + Int2Str() + " =" + Int2Str() +
//				+ Int2Str(nspcells-1);
//MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
//		paramsInit->setInit(init);
		initfile.close(); initfile.clear();
		// extract path name from full file name
		string fname = initname;
		string path = fname;
		unsigned int loc = fname.find("\\",0);
#if RSWIN64
		while(loc < 999999)
#else
		while(loc != string::npos)
#endif
		{
			fname = fname.substr(loc+1);
			loc =  fname.find("\\",0);
		}
		path = path.substr(0,path.length()-fname.length());
//string msg = path;
//MessageDlg(msg.c_str(),mtWarning, TMsgDlgButtons() << mbOK,0);

		string spdist = " species distribution";
		string cellsfile = " cells file";
		string pchfile = " patches file";
		string formatError = "Format error in header line of";
		if (init.seedType == 1) { // from species distribution
			// check that a distribution has been loaded
			// (NB there is no way of checking that it is the correct one)
			if (pLandscape->distnCount() < 1) {
				string msg = "Initial" + spdist + " has not been loaded";
				MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
				return;
			}
		}

		if (initcellsname != "-9") { // read additional patches/cells
			initcellsfile.open((path + initcellsname).c_str());
			if (!initcellsfile.is_open()) {
				string msg = "Error opening additional";
				if (ppLand.patchModel) msg += pchfile;
				else msg += cellsfile;
				MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
				initcellsfile.clear();
				return;
			}
			h0 = h1 = "q";
			if (ppLand.patchModel) {
				initcellsfile >> h0;
				if (h0 != "PatchID") {
					string msg = formatError + " additional" + pchfile;
					MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
					initcellsfile.close(); initcellsfile.clear();
					return;
				}
			}
			else {
				initcellsfile >> h0 >> h1;
				if (h0 != "x" || h1 != "y") {
					string msg = formatError + " additional" + cellsfile;
					MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
					initcellsfile.close(); initcellsfile.clear();
					return;
				}
			}
			pLandscape->clearInitCells();
			int x,y;
			if (ppLand.patchModel) {
				while (initcellsfile >> x) {
					pLandscape->addInitCell(x,-1);
				}
			}
			else {
				while (initcellsfile >> x >> y) {
					pLandscape->addInitCell(x,y);
				}
      }
			initcellsfile.close(); initcellsfile.clear();
		}

		if (initspdistcellsname != "-9") { // read species distribution cells
			// NOTE: it is assumed that the correct landscape and species distribution map
			// have been loaded
			initspdistcellsfile.open((path + initspdistcellsname).c_str());
			if (!initspdistcellsfile.is_open()) {
				string msg = "Error opening" + spdist + cellsfile;
				MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
				initspdistcellsfile.clear();
				return;
			}
			h0 = h1 = h2 = "q";
			initspdistcellsfile >> h0 >> h1 >> h2;
			if (h0 != "x" || h1 != "y" || h2 != "pres") {
				string msg = formatError + spdist + cellsfile;
				MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
				initspdistcellsfile.close(); initspdistcellsfile.clear();
				return;
			}
			int pres;
			locn loc;
			pLandscape->resetDistribution(pSpecies);
			while (initspdistcellsfile >> loc.x >> loc.y >> pres){
//string msg = "x=" + Int2Str(loc.x) + " y=" + Int2Str(loc.y) + " pres=" + Int2Str(pres) ;
//MessageDlg(msg.c_str(),mtWarning, TMsgDlgButtons() << mbOK,0);
				if (pres > 0) {
					pLandscape->setDistnCell(0,loc,true);
				}
			}
			initspdistcellsfile.close(); initspdistcellsfile.clear();
		}

//			MessageDlg("File OK",
//				mtWarning, TMsgDlgButtons() << mbOK,0);
	}
	break;

}

// no errors detected - close form
paramsInit->setInit(init);
forceInit = false;
Close();

}

//---------------------------------------------------------------------------
int __fastcall TfrmSeeding::CheckInitIndsFile(string fname) {
string header,msg;
string msg0 = "Error in selected file at line ";
int year,species,patchID,x,y,ninds,sex,age,stage,prevyear;
landParams ppLand = pLandscape->getLandParams();
demogrParams dem = pSpecies->getDemogr();

int errors = 0;

ifstream initindsfile;

// Open specified file
initindsfile.open(fname.c_str());
if (!initindsfile.is_open()) {
	initindsfile.clear();
	return -999;
}

// Check header line;
initindsfile >> header; if (header != "Year" ) errors++;
initindsfile >> header; if (header != "Species" ) errors++;
if (ppLand.patchModel) {
	initindsfile >> header; if (header != "PatchID" ) errors++;
}
else {
	initindsfile >> header; if (header != "X" ) errors++;
	initindsfile >> header; if (header != "Y" ) errors++;
}
initindsfile >> header; if (header != "Ninds" ) errors++;
if (dem.repType > 0) {
	initindsfile >> header; if (header != "Sex" ) errors++;
}
if (dem.stageStruct) {
	initindsfile >> header; if (header != "Age" ) errors++;
	initindsfile >> header; if (header != "Stage" ) errors++;
}
if (errors > 0) {
	msg =
		"Error in header line of selected file\nPlease ensure format matches model settings";
	MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
	initindsfile.close(); initindsfile.clear();
	return -1;
}

// Check data lines
int line = 1;
year = prevyear = -98765;
initindsfile >> year;
while (year != -98765) {
	if (year < 0) {
		msg = msg0 + Int2Str(line) + "\nYear must be >=0";
		MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
		initindsfile.close(); initindsfile.clear();
		return -1;
	}
	else {
		if (year < prevyear) {
			msg = msg0 + Int2Str(line) + "\nYear must be >= previous Year";
			MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
			initindsfile.close(); initindsfile.clear();
			return -1;
		}
	}
	prevyear = year;
	initindsfile >> species;
	if (species != 0) {
		msg = msg0 + Int2Str(line) + "\nSpecies must be 0";
		MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
		initindsfile.close(); initindsfile.clear();
		return -1;
	}
	if (ppLand.patchModel) {
		initindsfile >> patchID;
		if (patchID < 1) {
			msg = msg0 + Int2Str(line) + "\nPatchID must be >= 1";
			MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
			initindsfile.close(); initindsfile.clear();
			return -1;
		}
	}
	else {
		initindsfile >> x >> y ;
		if (x < 0 || y < 0) {
			msg = msg0 + Int2Str(line) + "\nX and Y must be >= 0";
			MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
			initindsfile.close(); initindsfile.clear();
			return -1;
		}
	}
	initindsfile >> ninds;
	if (ninds < 1) {
		msg = msg0 + Int2Str(line) + "\nNinds must be > 0";
		MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
		initindsfile.close(); initindsfile.clear();
		return -1;
	}
	if (dem.repType > 0) {
		initindsfile >> sex;
		if (sex < 0 || sex > 1) {
			msg = msg0 + Int2Str(line) + "\nSex must be 0 or 1";
			MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
			initindsfile.close(); initindsfile.clear();
			return -1;
		}
	}
	if (dem.stageStruct) {
		stageParams sstruct = pSpecies->getStage();
		initindsfile >> age >> stage;
		if (age < 1) {
			msg = msg0 + Int2Str(line) + "\nAge must be >= 1";
			MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
			initindsfile.close(); initindsfile.clear();
			return -1;
		}
		if (stage < 1) {
			msg = msg0 + Int2Str(line) + "\nStage must be >= 1";
			MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
			initindsfile.close(); initindsfile.clear();
			return -1;
		}
		if (stage >= sstruct.nStages) {
			msg = msg0 + Int2Str(line) + "\nStage must be < no. of stages";
			MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
			initindsfile.close(); initindsfile.clear();
			return -1;
		}
	}

	// read next line
	line++;
	year = -98765;
	initindsfile >> year;
	if (initindsfile.eof()) year = -98765;
} // end of while loop

if (initindsfile.is_open()) initindsfile.close();
initindsfile.clear();

return 0;
}



//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


