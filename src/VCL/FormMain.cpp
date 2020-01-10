#include <vcl.h>
#pragma hdrstop

#include "FormMain.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TfrmMain *frmMain;

paramGrad *paramsGrad;			// pointer to environmental gradient parameters
paramStoch *paramsStoch;		// pointer to environmental stochasticity parameters
paramInit *paramsInit;			// pointer to initialisation parameters
paramSim *paramsSim;				// pointer to simulation parameters
Landscape *pLandscape;  		// pointer to landscape
Species *pSpecies;					// pointer to species
Community *pComm;						// pointer to community

RSrandom *pRandom;

#if RSDEBUG
ofstream DEBUGLOG;
ofstream MUTNLOG;
ofstream DEBUGLOGGUI;
#endif

//---------------------------------------------------------------------------
__fastcall TfrmMain::TfrmMain(TComponent* Owner)
	: TForm(Owner)
{
pLandscape = new Landscape;
paramsGrad = new paramGrad;
paramsStoch = new paramStoch;
pSpecies = new Species;
paramsInit = new paramInit;
paramsSim = new paramSim;

#if RSWIN64
Caption = "RangeShifter v2.0  -  64 bit implementation";
#else
Caption = "RangeShifter v2.0  -  32 bit implementation";
#endif

// set up random number stream
pRandom = new RSrandom();

//MemoLine(("TfrmMain::TfrmMain(): frmGenerateLand = "
//	+ Int2Str((int)frmGenerateLand)).c_str());
//MemoLine(("TfrmMain::TfrmMain(): frmSpecies = "
//	+ Int2Str((int)frmSpecies)).c_str());

}
//---------------------------------------------------------------------------
const string Int2Str(const int x) {
	ostringstream o;
	if (!(o << x)) return "ERROR";
	return o.str();
}
const string Float2Str(const float x) {
	ostringstream o;
	if (!(o << x)) return "ERROR";
	return o.str();
}
const string Double2Str(const double x) {
	ostringstream o;
	if (!(o << x)) return "ERROR";
	return o.str();
}
//---------------------------------------------------------------------------

void MemoLine(UnicodeString msg) {
frmMain->Memo1->Lines->Add(msg);
}

#if RSDEBUG
void DebugGUI(string msg) {
DEBUGLOGGUI << msg << endl;
}
#endif

//---------------------------------------------------------------------------
const rgb draw_wheel(int nFactor)
{
rgb colours;

if (nFactor < 256) {
	if (nFactor < 0) nFactor = 0;
	// red 0, green increasing, blue 255 (BLUE to LIGHT BLUE)
	colours.r = 0; colours.g = nFactor; colours.b = 255;
}
else {
	if (nFactor < 512) {
		// red 0, green 255, blue decreasing.  (LIGHT BLUE to GREEN)
		colours.r = 0; colours.g = 255; colours.b = 511 - nFactor;
	}
	else {
		if (nFactor > 767) nFactor = 767;
			// red increasing, green 255, blue 0 (GREEN to YELLOW)
			colours.r = nFactor - 512; colours.g = 255; colours.b = 0;
	}
}

return colours;
}

//---------------------------------------------------------------------------
void __fastcall TfrmMain::SetDirectoryClick(TObject *Sender)
{

//MemoLine(("TfrmMain::SetDirectoryClick(): 0000 frmGenerateLand = "
//	+ Int2Str((int)frmGenerateLand)).c_str());
//MemoLine(("TfrmMain::SetDirectoryClick(): 0000 frmSpecies = "
//	+ Int2Str((int)frmSpecies)).c_str());

OpenDialog1->Title = "Select working directory";
if (OpenDialog1->Execute()) {
	String filename = OpenDialog1->FileName;
	string nameS = AnsiString(filename).c_str();
	string path = nameS;

	unsigned int loc = nameS.find( "\\", 0 );
#if RSWIN64
	while(loc < 999999)
#else
	while(loc != string::npos)
#endif
	{
		nameS = nameS.substr(loc+1);
		loc =  nameS.find( "\\", 0 );
	}
	path = path.substr(0,path.length()-nameS.length());
//	MemoLine(("SetDirectoryClick(): Path = " + path).c_str());
	paramsSim->setDir(path);

	// check that Inputs, Outputs and Output_Maps folders exist
	bool errorfolder = CheckDirectory();
	if (errorfolder) {
		string msg = "Working directory must contain Inputs, Outputs and Output_Maps folders";
		MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
		return;
	}

	directoryChosen = true;
//		frmMain->Memo1->Lines->Add("CHOSEN ... ");
	frmMain->BatchMode->Enabled = true;
	frmMain->LandMenu->Enabled= true;
	if (RSlog.is_open()) {
		RSlog.close(); RSlog.clear();
	}
}

if (directoryChosen) {
	string name = paramsSim->getDir(2) + "RS_log.csv";
	if (!RSlog.is_open()) {
		RSlog.open(name.c_str());
	}
	RSlog << "Event,Number,Reps,Years,Time" << endl;
#if RSDEBUG
	RSlog << "WARNING,***** RSDEBUG mode is active *****,,," << endl;
	name = paramsSim->getDir(2) + "DebugGUI.txt";
	if (!DEBUGLOGGUI.is_open()) {
		DEBUGLOGGUI.open(name.c_str());
	}
	DEBUGLOGGUI << "TfrmMain::SetDirectoryClick(): ";
	for (int i = 0; i < 5; i++) {
		int rrrr = pRandom->IRandom(1000,2000);
		DEBUGLOGGUI << " " << rrrr;
	}
	DEBUGLOGGUI << endl;
//	if (DEBUGLOGGUI.is_open())
//		MessageDlg("TfrmMain: DEBUGLOGGUI is open",mtWarning,TMsgDlgButtons() << mbOK,0);
//	else
//		MessageDlg("TfrmMain: DEBUGLOGGUI is NOT open",mtWarning,TMsgDlgButtons() << mbOK,0);
//	DEBUGLOGGUI.close(); DEBUGLOGGUI.clear();
//	if (DEBUGLOGGUI.is_open())
//		MessageDlg("TfrmMain: DEBUGLOGGUI is open",mtWarning,TMsgDlgButtons() << mbOK,0);
//	else
//		MessageDlg("TfrmMain: DEBUGLOGGUI is NOT open",mtWarning,TMsgDlgButtons() << mbOK,0);
#endif

}
else {
//	frmMain->Memo1->Lines->Add("NOT CHOSEN ... ");
	frmMain->BatchMode->Enabled = false;
}
//MemoLine(("SetDirectoryClick(): Path = " + paramsSim->getDir(0)).c_str());

//MemoLine(("TfrmMain::SetDirectoryClick(): 9999 frmGenerateLand = "
//	+ Int2Str((int)frmGenerateLand)).c_str());
//MemoLine(("TfrmMain::SetDirectoryClick(): 9999 frmSpecies = "
//	+ Int2Str((int)frmSpecies)).c_str());

}

//---------------------------------------------------------------------------
void __fastcall TfrmMain::BatchModeClick(TObject *Sender)
{
batchfiles b;
landParams ppLand;
demogrParams dem;
stageParams sstruct;
trfrRules trfr;
simParams sim = paramsSim->getSim();

// create a new landscape
if (pLandscape != NULL) delete pLandscape;
pLandscape = new Landscape;
//ppLand = pLandscape->getLandParams();

if (BatchMode->Checked) {
	batchMode = false;
	BatchMode->Checked = false;
	LandMenu->Enabled = true;
	Parameterset->Enabled = false;
	Run1->Enabled = false;
}
else {
	OpenDialog2->InitialDir = paramsSim->getDir(1).c_str();
	OpenDialog2->Title = "Select the Control file.";
	if (OpenDialog2->Execute()){
		controlFile_name = AnsiString(OpenDialog2->FileName).c_str();
		MemoLine(("Checking the Control file " + controlFile_name).c_str());
		b = ParseControlFile(controlFile_name,paramsSim->getDir(1),paramsSim->getDir(2));
		if (b.ok) {
			ppLand = pLandscape->getLandParams();
			nSimuls = b.nSimuls;
			nLandscapes = b.nLandscapes;
			ppLand.patchModel = b.patchmodel;
			ppLand.resol = b.resolution;
			ppLand.nHabMax = b.maxNhab;
			if (b.landtype == 9) {
				ppLand.generated = true;
				ppLand.nHab = 2;
			}
			else {
				ppLand.generated = false;
				ppLand.rasterType = b.landtype;
				if (b.landtype == 2) ppLand.nHab = 1; // habitat quality landscape
			}
			ppLand.spDist = b.speciesdist;
			ppLand.spResol = b.distresolution;
			dem.repType = b.reproductn;
			dem.repSeasons = b.repseasons;
			dem.stageStruct = b.stagestruct;
			sstruct.nStages = b.stages;
			if (b.transfer == 0) trfr.moveModel = false;
			else {
				trfr.moveModel = true;
				trfr.moveType = b.transfer;
			}
			batchMode = true;
			BatchMode->Checked = true;
			LandMenu->Enabled = false;
			Parameterset->Enabled = false;
			Run1->Enabled = true;
			pLandscape->setLandParams(ppLand,true);
			pSpecies->setDemogr(dem);
			pSpecies->setStage(sstruct);
			pSpecies->setTrfr(trfr);
			MemoLine("Control file valid");
		}
		else { // ERROR
			MessageDlg("Error in batch input.\nSee batch log file in Outputs folder for details."
				,mtError, TMsgDlgButtons() << mbOK,0);
			batchMode = false;
			BatchMode->Checked = false;
			LandMenu->Enabled = true;
			Parameterset->Enabled = false;
			Run1->Enabled = false;
		}
	}
}
sim.batchNum = b.batchNum;
sim.batchMode = batchMode;
paramsSim->setSim(sim);
}

//---------------------------------------------------------------------------
void __fastcall TfrmMain::LandMenuClick(TObject *Sender)
{
if (!directoryChosen)
{
MessageDlg("Please select the work Directory to proceed.",
	mtWarning, TMsgDlgButtons() << mbOK,0);
}

//MemoLine(("TfrmMain::LandMenuClick(): frmGenerateLand = "
//	+ Int2Str((int)frmGenerateLand)).c_str());
//MemoLine(("TfrmMain::LandMenuClick(): frmSpecies = "
//	+ Int2Str((int)frmSpecies)).c_str());

}

//---------------------------------------------------------------------------
void __fastcall TfrmMain::RasterLandClick(TObject *Sender)
{
#if RSDEBUG
//	if (DEBUGLOGGUI.is_open())
//		MessageDlg("RasterLandClick: DEBUGLOGGUI is open",mtWarning,TMsgDlgButtons() << mbOK,0);
//	else
//		MessageDlg("RasterLandClick: DEBUGLOGGUI is NOT open",mtWarning,TMsgDlgButtons() << mbOK,0);
#endif
landParams ppLand;

Run1->Enabled = false;
frmLand->refresh(landscapeLoaded);
frmLand->ShowModal();

if (landscapeLoaded) {
	ppLand = pLandscape->getLandParams();
	if (ppLand.rasterType == 0) { // habitat codes
		if (pLandscape->habitatsIndexed()) {
			LandImage->Visible = true;
		}
		else {
			Parameterset->Enabled = false;
			LandImage->Visible = false;
			MemoLine("Landscape habitat indices not set, please re-load");
			return;
		}
	}
	else
		LandImage->Visible = true;
}
else {
	Parameterset->Enabled = false;
	LandImage->Visible = false;
	MemoLine("Landscape erased, please re-load");
	return;
}
//forceInit = true;
simseq = 0; // ensure patches are reallocated for new cell-based landscape

UnicodeString units;
if (ppLand.patchModel) {
	units = "patches";
	// environmental gradient is not permitted for a patch-based model
	EnvGradient->Enabled = false;
	paramsGrad->noGradient();
}
else {
	units = "cells";
	EnvGradient->Enabled = true;
}
ChartPop->RightAxis->Title->Caption = "Occupied " + units;
ChartOccSuit->Title->Caption = "Proportion of suitable " + units + " occupied";

File1->Enabled = false;
Artificial->Enabled = false;
Parameterset->Enabled = true;
SpeciesMenu->Enabled = true;
GeneticsMenu->Enabled = false;
Simulations->Enabled = false;
//Refresh->Enabled = true;
DrawLandscape();
if (pLandscape->distnCount() > 0) { // initial distribution has been loaded
	// set the scale for drawing the initial distribution
	Sp_resol_ratio = ppLand.spResol/ppLand.resol;
#if RSDEBUG
//UnicodeString pmsg2 = "TfrmMain::RasterLandClick() 2: Sp_resol_ratio = " + IntToStr(Sp_resol_ratio);
//MemoLine(pmsg2);
#endif
	DrawInitial(true);
	MemoLine("Initial distribution loaded");
}
else {
	DrawInitial(false);
}

if (frmLand->CBVisualPatch->Checked) {
	if (pLandscape->patchCount() > 0) {
		//setting image
		frmVisualPatch->PatchImage->Height = frmMain->LandImage->Height;
		frmVisualPatch->PatchImage->Width  = frmMain->LandImage->Width;
		frmVisualPatch->PatchScrollBox->VertScrollBar->Range = frmMain->LandImage->Height;
		frmVisualPatch->PatchScrollBox->HorzScrollBar->Range = frmMain->LandImage->Width;
		if (bmpPatches == NULL) bmpPatches = new Graphics::TBitmap();
		canvasPatches = bmpPatches->Canvas;
		bmpPatches->Height = bmpLand->Height;
		bmpPatches->Width  = bmpLand->Width;
		frmVisualPatch->Show();
		DrawPatches();
	}
}
else
	if (frmVisualPatch->Showing) frmVisualPatch->Close();

}

//---------------------------------------------------------------------------
void __fastcall TfrmMain::ArtificialClick(TObject *Sender)
{
landParams ppLand;

//MemoLine(("TfrmMain::ArtificialClick(): frmGenerateLand = "
//	+ Int2Str((int)frmGenerateLand)).c_str());
//MemoLine(("TfrmMain::ArtificialClick(): frmSpecies = "
//	+ Int2Str((int)frmSpecies)).c_str());
//MemoLine(("TfrmMain::ArtificialClick(): pLandscape = "
//	+ Int2Str((int)pLandscape)).c_str());

if (pLandscape != 0) {
	delete pLandscape;
	pLandscape = new Landscape;
	frmMain->Parameterset->Enabled = false;
	frmMain->SpeciesMenu->Enabled = false;
}
artLandStatus = 0;
while (frmGenerateLand->ShowModal() == 1) {} // show form until OK or cancelled

//string msg;
//msg = "ArtificialClick: artform = " + Int2Str(artform);
//MessageDlg(msg.c_str(),mtWarning, TMsgDlgButtons() << mbOK,0);

ppLand = pLandscape->getLandParams();
forceInit = true;

//string msg;
//msg = "ArtificialClick: artLandStatus = " + Int2Str(artLandStatus);
//MessageDlg(msg.c_str(),mtWarning, TMsgDlgButtons() << mbOK,0);
if (artLandStatus) {
	if (artLandStatus == 1) {
		MemoLine("Fractal maps series has been created");
	}
	if (artLandStatus == 2) {
		MemoLine("Random maps series has been created");
	}
	if (artLandStatus == 9) {
		MemoLine("Artificial landscape parameters OK");
	}
}
else {
	return;
}

//string msg;
//msg = "ArtificialClick: maxX " + Int2Str(ppLand.maxX) + " maxY " + Int2Str(ppLand.maxY);
//MessageDlg(msg.c_str(),mtWarning, TMsgDlgButtons() << mbOK,0);

frmMain->File1->Enabled = false;

if (artLandStatus < 9) {
	// no landscape is loaded, therefore user may proceed no further,
	// but may generate further series or generate a single landscape for immediate use
	return;
}
frmMain->RasterLand->Enabled = false;
frmMain->Parameterset->Enabled = true;
frmMain->SpeciesMenu->Enabled = true;
frmMain->GeneticsMenu->Enabled = false;
frmMain->Simulations->Enabled = false;
if (!ppLand.patchModel) frmMain->EnvGradient->Enabled = true;

newLandscape = true;
costsSet = false; // force costs to be respecified (for SMS only)
frmVisualCost->Close();

}

//---------------------------------------------------------------------------
void __fastcall TfrmMain::EnvGradientClick(TObject *Sender)
{
frmMain->Run1->Enabled = false;
frmEnvGradient->Refresh();
frmEnvGradient->Show();
}

//------------------------------------------------------------------------------
void __fastcall TfrmMain::SpeciesMenuClick(TObject *Sender)
{
string msg1 = "Species parameters saved. Please set ";
string msg2 = " parameters.";
landParams ppLand = pLandscape->getLandParams();

Simulations->Enabled = false;
Run1->Enabled = false;

frmSpecies->refresh(ppLand.patchModel);

speciesOK = false;
while (frmSpecies->ShowModal() == 1) {} // show form until OK or cancelled
if (!speciesOK)
	// user has cancelled or closed Species form - prevent potentially invalid parameters
	// being set and Genetics / Simulation form being activated
	return;

emigRules emig = pSpecies->getEmig();
trfrRules trfr = pSpecies->getTrfr();
settleType sett = pSpecies->getSettle();
if (emig.indVar || trfr.indVar || sett.indVar) {
	// individual variation - genetics must be specified
	MemoLine((msg1 + "genetics" + msg2).c_str());
	GeneticsMenu->Enabled = true;
//	GeneticsMenu->Caption = "Genetics";
	Simulations->Enabled = false;
}
else { // no individual variation
	// optional neutral genetics OR go straight to simulation form
	MemoLine((msg1 + "neutral genetics (optional) and simulation" + msg2).c_str());
	genomeData gen = pSpecies->getGenomeData();
	GeneticsMenu->Enabled = true;
	if (gen.neutralMarkers) Simulations->Enabled = false;
	else Simulations->Enabled = true;
}

Application->ProcessMessages();
}

//------------------------------------------------------------------------------
void __fastcall TfrmMain::GeneticsMenuClick(TObject *Sender)
{
Run1->Enabled = false;

demogrParams dem = pSpecies->getDemogr();
bool sexual;
if (dem.repType == 0) sexual = false; else sexual = true;

frmGenetics->refresh(sexual);

geneticsOK = false;
while (frmGenetics->ShowModal() == 1) {} // show form until OK or cancelled
//MemoLine("RETURNED FROM GENETICS FORM");
if (!geneticsOK) {
	// user has cancelled or closed Genetics form - prevent potentially invalid parameters
	// being set and Simulation form being activated
//	MemoLine("MAIN FORM - NOT geneticsOK");
	return;
}

MemoLine("Genetics parameters saved. Please set simulation parameters.");
frmMain->Simulations->Enabled = true;

Application->ProcessMessages();
}

//---------------------------------------------------------------------------

void __fastcall TfrmMain::SimulationsClick(TObject *Sender)
{
landParams ppLand = pLandscape->getLandParams();

frmSim->refresh(ppLand.patchModel,ppLand.generated);
readyToRun = false;
cancelled = false;
do {
	frmSim->ShowModal();
} while (!readyToRun && !cancelled);

if (chartNoise) frmMain->ChartNoise->Visible = true;
else frmMain->ChartNoise->Visible = false;

if (readyToRun) frmMain->Run1->Enabled = true;

}

//---------------------------------------------------------------------------
// Main function - Run
void __fastcall TfrmMain::Run1Click(TObject *Sender)
{
simParams sim = paramsSim->getSim();
#if RSDEBUG
string name = paramsSim->getDir(2) + "Sim" + Int2Str(sim.simulation)
//	+ "_" + Int2Str(batch_line)
	+ "_DebugLog.txt";
DEBUGLOG.open(name.c_str());
if (!DEBUGLOG.is_open()) {
	MessageDlg("Unable to open DEBUGLOG file",
		mtError, TMsgDlgButtons() << mbOK,0);
	return;
}
DEBUGLOG << "Run1Click(): random integers:";
for (int i = 0; i < 5; i++) {
	int rrrr = pRandom->IRandom(1000,2000);
DEBUGLOG << " " << rrrr;
}
DEBUGLOG << endl;
//if (DEBUGLOG.is_open())
//	MessageDlg("Run1Click 222: DEBUGLOG is open",mtWarning,TMsgDlgButtons() << mbOK,0);
//else
//	MessageDlg("Run1Click 222: DEBUGLOG is NOT open",mtWarning,TMsgDlgButtons() << mbOK,0);
//DEBUGLOG.close(); DEBUGLOG.clear();
//return;
#endif
landParams ppLand = pLandscape->getLandParams();
if (ppLand.spResol == ppLand.resol) {
	// prevent initial distribution being shown further during simulation
	frmMain->SpDistImage->Visible = false;
}

Run1->Enabled = false; Stop->Enabled = true; Pause->Enabled = true;
// other menu options are disabled whilst running until Refresh
File1->Enabled = false;
LandMenu->Enabled = false;
Parameterset->Enabled = false;
Refresh->Enabled = false;

if (pause) {
	pause = false;
	frmMain->Run1->Enabled = false;
	MemoLine("Simulation running...");
}
else {
	if (batchMode) {
		RunBatch(nSimuls,nLandscapes);
	}
	else {
		MemoLine("----------------------");
		MemoLine(("Running Simulation nr. " + Int2Str(sim.simulation)).c_str());

		OutParameters(pLandscape);

		int t0,t1;
		t0 = time(0);

#if RSDEBUG
//if (DEBUGLOG.is_open())
//	MessageDlg("Run1Click 111: DEBUGLOG is open",mtWarning,TMsgDlgButtons() << mbOK,0);
//else
//	MessageDlg("Run1Click 111: DEBUGLOG is NOT open",mtWarning,TMsgDlgButtons() << mbOK,0);
//emigLimits eee = pSpecies->getEmigLimits(0,1);
//DEBUGLOG << "RunReal(): eee.minD0 = " << eee.minD0 << " eee.maxD0 = " << eee.maxD0 << endl;
#endif

		int retcode = RunModel(pLandscape,simseq++);
		if (retcode == 666) {
			MessageDlg("Unable to open one or more output files",
				mtError, TMsgDlgButtons() << mbOK,0);
			Stop->Enabled = false;
			Pause->Enabled = false;
			Refresh->Enabled = true;
			return;
		}

#if RSDEBUG
//if (DEBUGLOG.is_open())
//	MessageDlg("RunReal 999: DEBUGLOG is open",mtWarning,TMsgDlgButtons() << mbOK,0);
//else
//	MessageDlg("RunReal 999: DEBUGLOG is NOT open",mtWarning,TMsgDlgButtons() << mbOK,0);
t1 = time(0);
DEBUGLOG << endl << endl << "Run1Click(): TOTAL TIME = " << t1-t0 << " seconds"
	<< endl << endl;
DEBUGLOG << "Run1Click(): end" << endl;
#endif

		// Write performance data
		t1 = time(0);
		RSlog << "Simulation," << sim.simulation << "," << sim.reps << "," << sim.years
			<< "," << t1-t0 << endl;

	} // end of batchMode == false

	if (stopRun) { // simulation stopped
		stopRun = false;
		Stop->Enabled = false;
	}
	else { // normal termination
		Stop->Enabled = false;
		Pause->Enabled = false;
		MemoLine("SIMULATION COMPLETED.");
		if (batchMode) {
			File1->Enabled = false;
			Pause->Enabled = false;
			MemoLine("Please Close the program.");
		}
		else {
			MemoLine("Please Refresh or Close the program.");
		}
	}
	if (!batchMode) Refresh->Enabled = true;
}

#if RSDEBUG
//	 DEBUGLOG.close();
#endif

}

//------------------------------------------------------------------------------
void __fastcall TfrmMain::BtnZoomInClick(TObject *Sender)
{
landPix p = pLandscape->getLandPix();
landData land = pLandscape->getLandData();
#if RSDEBUG
//DebugGUI(("TfrmMain::BtnZoomInClick(): 00000 p.pix=" + Int2Str(p.pix)
//	).c_str());
#endif
//MemoLine(("BtnZoomInClick(): p.pix = " + Int2Str(p.pix)).c_str());
if (p.pix < 5) p.pix++;
else {
	p.pix = (int)((float)p.pix * 1.25);
}
if (p.pix > 600) {
	p.pix = 600;
	BtnZoomIn->Enabled = false;
}
pLandscape->setLandPix(p);
#if RSDEBUG
//DebugGUI(("TfrmMain::BtnZoomInClick(): 99999 p.pix=" + Int2Str(p.pix)
//	).c_str());
#endif

bmpLand->Height = land.dimY * p.pix;
bmpLand->Width  = land.dimX * p.pix;

SetMapDimensions();
pLandscape->drawLandscape(0,0,0);
if (pLandscape->distnCount() > 0) DrawInitial(true);
BtnZoomOut->Enabled = true;
}

//---------------------------------------------------------------------------
void __fastcall TfrmMain::BtnZoomOutClick(TObject *Sender)
{
landPix p = pLandscape->getLandPix();
landData land = pLandscape->getLandData();
#if RSDEBUG
//DebugGUI(("TfrmMain::BtnZoomOutClick(): 00000 p.pix=" + Int2Str(p.pix)
//	).c_str());
#endif
p.pix = (int)((float)p.pix * 0.8);
if (p.pix < 2) {
	p.pix = 1;
	BtnZoomOut->Enabled = false;
}
pLandscape->setLandPix(p);
#if RSDEBUG
//DebugGUI(("TfrmMain::BtnZoomOutClick(): 99999 p.pix=" + Int2Str(p.pix)
//	).c_str());
#endif

bmpLand->Height = land.dimY * p.pix;
bmpLand->Width  = land.dimX * p.pix;

SetMapDimensions();
pLandscape->drawLandscape(0,0,0);
if (pLandscape->distnCount() > 0) DrawInitial(true);
BtnZoomIn->Enabled = true;
}

void __fastcall TfrmMain::SetMapDimensions(void) {
simView v = paramsSim->getViews();
landParams ppLand = pLandscape->getLandParams();
landPix p = pLandscape->getLandPix();

#if RSDEBUG
//DebugGUI(("TfrmMain::SetMapDimensions(): bmpLand->Height=" + Int2Str(bmpLand->Height)
//	).c_str());
#endif

LandImage->Height = bmpLand->Height;
LandImage->Width  = bmpLand->Width;
CommImage->Height = bmpLand->Height;
CommImage->Width  = bmpLand->Width;
SpDistImage->Height = bmpLand->Height;
SpDistImage->Width  = bmpLand->Width;
MovtPaintBox->Height = bmpLand->Height;
MovtPaintBox->Width  = bmpLand->Width;
LandScrollBox->VertScrollBar->Range = bmpLand->Height;
LandScrollBox->HorzScrollBar->Range = bmpLand->Width;
if (!ppLand.generated) {
	if (p.pix > 1) BtnZoomOut->Enabled = true;
	else BtnZoomOut->Enabled = false;
	if (p.pix < 600) BtnZoomIn->Enabled = true;
	else BtnZoomIn->Enabled = false;
}

if (v.viewLand) {
	LandImage->Picture->Bitmap = bmpLand;
	LandImage->Repaint();
}
if (v.viewPop) {
	CommImage->Refresh();
	CommImage->Picture->Bitmap = NULL;
	CommImage->Transparent = true;
	CommImage->Repaint();
}
}

//---------------------------------------------------------------------------
// Called from RunModel() prior to looping through replicates
void SetupVisualOutput(void)
{
landParams ppLand = pLandscape->getLandParams();
simParams sim = paramsSim->getSim();
simView v = paramsSim->getViews();
emigRules emig = pSpecies->getEmig();
trfrRules trfr = pSpecies->getTrfr();
envStochParams env = paramsStoch->getStoch();

if (v.viewGrad || v.viewTraits) {
	SetGradTraitsForms();
}

if (v.viewPop || sim.saveMaps) {
	int rs, gs;
	float lpix = (float)(frmMain->PopLegend->Height) / 280.0;
	if (bmpPopLgd == NULL) bmpPopLgd = new Graphics::TBitmap(); canvasPopLgd = bmpPopLgd->Canvas;
	bmpPopLgd->Width = frmMain->PopLegend->Width; bmpPopLgd->Height = frmMain->PopLegend->Height;
	if (v.viewPop) {
		float maxK = pSpecies->getMaxK();
		frmMain->PopLegend->Visible = true; frmMain->LabelPop->Visible = true;
		frmMain->LabelMaxPop->Visible = true; frmMain->LabelMinPop->Visible = true;
		if (ppLand.patchModel) frmMain->LabelMaxPop->Caption =
			(">" + Float2Str(maxK*10000/(float)(ppLand.resol*ppLand.resol)) + " ind/ha").c_str();
		else
			frmMain->LabelMaxPop->Caption = (">" + Int2Str(maxK+1.0) + " ind/cell").c_str();
	}
	for (int i = 100; i < 381; i++) {
		gs = 0;
		if(i <= 256) rs = (int)((double)255*(sin(((double)i/256)*PI/2)));
		else{
			rs = 255;
			gs = (int)((double)255*(sin(((double)(i-256)/256)*PI/2)));
		}
		canvasPopLgd->Brush->Color = (TColor)RGB(rs,gs,0);
		canvasPopLgd->FillRect(Rect(0,(int)((i-100)*lpix),10,(int)((i-100)*lpix+lpix)));
	}
	if (v.viewPop) {
		frmMain->PopLegend->Picture->Bitmap = bmpPopLgd;
		frmMain->PopLegend->Repaint();
	}
}
else {
	frmMain->PopLegend->Visible = false; frmMain->LabelPop->Visible = false;
	frmMain->LabelMaxPop->Visible = false; frmMain->LabelMinPop->Visible = false;
}

if (sim.outOccup) {
	frmMain->ChartOccSuit->Visible = true;
}
else {
	frmMain->ChartOccSuit->Visible = false;
}

}

void ResetVisualOutput(void)
{
envGradParams grad = paramsGrad->getGradient();
envStochParams env = paramsStoch->getStoch();
simView v = paramsSim->getViews();

if (v.viewGraph) { frmMain->Population->Clear(); frmMain->Cells->Clear(); }
if (env.stoch && !env.local) {
	if (v.viewGraph) frmMain->EnvNoise->Clear();
}
if (grad.gradient && v.viewGrad) {
	frmVisualTraits0->Image0->Picture->Bitmap = bmpGrad;
	frmVisualTraits0->Repaint();
}

}

// OUTPUTS ------------------------------------------------------------------

//---------------------------------------------------------------------------
// Update the occupancy graph on the screen
void Community::viewOccSuit(int year,double mn,double se)
{
simParams sim = paramsSim->getSim();
frmMain->OccSuitmean->AddXY(year,mn,"",clBlack);
frmMain->OccSuitplusSE->AddXY(year,(mn+se),"",clBlack);
frmMain->OccSuitminusSE->AddXY(year,(mn-se),"",clBlack);
}

//---------------------------------------------------------------------------
void SaveTraitMap(int rr, int yy, int gg, int landNr, string folder,
	Graphics::TBitmap *map, string trait)
{
#if RSDEBUG
//DEBUGLOG << "SaveTraitMap(): rr=" << rr << " yy=" << yy << " gg= " << gg
//	<< " landNr=" << landNr << " folder=" << folder << " map=" << map << " trait=" << trait
//	<< endl;
#endif
string mapName;
int simnum = paramsSim->getSimNum();

mapName = (folder + "Sim" + Int2Str(simnum) +
	+ "_land" + Int2Str(landNr) + "_rep" + Int2Str(rr) + "_" + trait);
if (yy == -1) {
	mapName = (mapName + "_initial.bmp");
}
else{
	mapName = (mapName + "_year" + Int2Str(yy)+ "_rs" + Int2Str(gg) + ".bmp");
}
map->SaveToFile(mapName.c_str());

#if RSDEBUG
//DEBUGLOG << "SaveTraitMap(): finished"
//	<< endl;
#endif
}

//---------------------------------------------------------------------------
void __fastcall TfrmMain::StopClick(TObject *Sender)
{
stopRun = true;
if (pause) pause = false;
frmMain->Pause->Enabled = false;
frmMain->Run1->Enabled = false;
MemoLine("Simulation stopped.");
if (batchMode) {
	frmMain->File1->Enabled = false;
	MemoLine("Please Close the program.");
}
else {
	frmMain->Refresh->Enabled = true;
	MemoLine("Please Refresh or Close the program.");
}
frmMain->Stop->Enabled = false;
}

//---------------------------------------------------------------------------
void __fastcall TfrmMain::PauseClick(TObject *Sender)
{
pause = true;
frmMain->Run1->Enabled = true;
MemoLine("Simulation paused. Press Run to continue...");
int p = 1;
while (p > 0){
	Application->ProcessMessages();
	if (stopRun) break;
	if(pause == false) break;
}
}

//---------------------------------------------------------------------------
void __fastcall TfrmMain::RefreshClick(TObject *Sender)
{
landParams ppLand = pLandscape->getLandParams();
envGradParams grad = paramsGrad->getGradient();
envStochParams env = paramsStoch->getStoch();
trfrRules trfr = pSpecies->getTrfr();
simParams sim = paramsSim->getSim();
simView v = paramsSim->getViews();

#if RSDEBUG
//DEBUGLOG << "RefreshClick() starting: Inds.size() " << Inds.size()
//	<< " Temp_Inds.size() " << Temp_Inds.size() << " Offs.size() " << Offs.size()
//	<< endl;
#endif

#if RSDEBUG
//DEBUGLOG << "RefreshClick(): 00000 " << endl;
#endif

if (env.stoch && !env.local) delete epsGlobal;

if (grad.gradient) paramsGrad->resetOptY();

if (sim.outConnect && ppLand.patchModel)
	pLandscape->deleteConnectMatrix();
#if RSDEBUG
//DEBUGLOG << "RefreshClick(): 11111 " << endl;
#endif

pLandscape->resetPatchPopns();
pLandscape->resetLandLimits();
if (pComm != 0) {
	delete pComm;
}

#if RSDEBUG
//DEBUGLOG << "RefreshClick(): 22222 " << endl;
#endif

// Visualisations -----------------------------------------------------------
if (v.viewLand) {
	if (ppLand.generated) LandImage->Picture->Bitmap = NULL;
	else
//		DrawLandscape();
		pLandscape->drawLandscape(0,0,0);
	LandImage->Repaint();
}
if (v.viewGraph) {
	Population->Clear(); Cells->Clear();
	if (sim.outOccup) {
		OccSuitmean->Clear(); OccSuitplusSE->Clear(); OccSuitminusSE->Clear();
//		ChartOccSuit->Visible = false;
	}
	if (env.stoch && !env.local)
		EnvNoise->Clear();
}


//if (v.viewLand) {
//	LandImage->Picture->Bitmap = bmpLand;
//	LandImage->Repaint();
//}
if (v.viewPop) {
	CommImage->Refresh();
	CommImage->Picture->Bitmap = NULL;
	CommImage->Transparent = true;
	CommImage->Repaint();
}

if (v.viewPop) {
	if (ppLand.generated) {
		CommImage->Refresh();
		CommImage->Picture->Bitmap = NULL;
		CommImage->Transparent = true;
		CommImage->Repaint();
	}
	if (bmpPopLgd != NULL){ delete bmpPopLgd; bmpPopLgd = NULL;}
	PopLegend->Picture->Bitmap = NULL;
	PopLegend->Visible = false; frmMain->LabelPop->Visible = false;
	LabelMaxPop->Visible = false; frmMain->LabelMinPop->Visible = false;
}

if(v.viewGrad || v.viewTraits) {
	for (int i = 0; i < NTRAITS; i++) {
		if (bmpT[i] != NULL){ delete bmpT[i]; bmpT[i] = NULL;}
		if (bmpL[i] != NULL){ delete bmpL[i]; bmpL[i] = NULL;}
	}
	if (frmVisualTraits0->Showing) frmVisualTraits0->Close();
	if (frmVisualTraits1->Showing) frmVisualTraits1->Close();
	if (frmVisualTraits2->Showing) frmVisualTraits2->Close();
}
#if RSDEBUG
//DEBUGLOG << "RefreshClick(): 88888 " << endl;
#endif

//if (sim.initDistLoaded && ppLand.spResol == ppLand.resol)
if (ppLand.spDist && ppLand.spResol == ppLand.resol)
{
	// prevent initial distribution being shown further during simulation
	SpDistImage->Visible = true;
}

MemoLine("Refreshed.");

Refresh->Enabled = false;
// Other menu options disabled whilst running can now be enabled
LandMenu->Enabled = true;
Parameterset->Enabled = true;
stopRun = false;

#if RSDEBUG
//DEBUGLOG << "RefreshClick(): 99999 " << endl;
#endif

#if RSDEBUG
if (DEBUGLOG.is_open()) { DEBUGLOG.close(); DEBUGLOG.clear(); }
#endif

}

//---------------------------------------------------------------------------
void DrawLandscape(void)
{
simView v = paramsSim->getViews();

pLandscape->setLandMap();
pLandscape->drawLandscape(0,0,0);
frmMain->SetMapDimensions();

//if (v.viewLand) {
//	frmMain->LandImage->Picture->Bitmap = bmpLand;
//	frmMain->LandImage->Repaint();
//}
//
//if (v.viewPop) {
//	frmMain->CommImage->Refresh();
//	frmMain->CommImage->Picture->Bitmap = NULL;
//	frmMain->CommImage->Transparent = true;
//	frmMain->CommImage->Repaint();
//}

}

void DrawPopnGraph(Community *pComm,int yr)
{
frmMain->Population->AddXY(yr,pComm->totalInds(),"",clRed);
commStats s = pComm->getStats();
frmMain->Cells->AddXY(yr,s.occupied,"",clBlue);
}

//---------------------------------------------------------------------------
void DrawInitial(bool loaded)
{
locn distCell;
landData ppLand = pLandscape->getLandData();
landPix p = pLandscape->getLandPix();

if (bmpSpDist != NULL) delete bmpSpDist;
bmpSpDist = new Graphics::TBitmap();
canvasSpDist = bmpSpDist->Canvas;
bmpSpDist->Height = ppLand.dimY * p.pix;
bmpSpDist->Width  = ppLand.dimX * p.pix;
bmpSpDist->Transparent = true;

if (loaded) {
	// CHANGE NEXT LINE FOR MULTIPLE SPECIES DISTRIBUTIONS ...
	locn distDim = pLandscape->getDistnDimensions(0);
	int diffY = (ppLand.dimY) - (distDim.y+1) * Sp_resol_ratio;
#if RSDEBUG
//MemoLine(("diffY = " + Int2Str(diffY)).c_str());
#endif
	canvasSpDist->Brush->Color=(TColor) RGB(255,255,0);

#if RSDEBUG
//UnicodeString pmsg0 = "DrawInitial() 0: p.pix = " + FloatToStr(p.pix);
//MemoLine(pmsg0);
//UnicodeString pmsg1 = "DrawInitial() 1: Sp_resol_ratio = " + IntToStr(Sp_resol_ratio);
//MemoLine(pmsg1);
#endif

	// NOTE: CURRENTLY DRAWING THE FIRST SPECIES DISTRIBUTION ONLY
	int ndistcells = pLandscape->distCellCount(0);
	for (int i = 0; i < ndistcells; i++) {
		distCell = pLandscape->getDistnCell(0,i);
#if RSDEBUG
//if (i == 0) {
//UnicodeString loc0 = "DrawInitial(): x0 = "
//	+ FloatToStr(loc.x * p.sppix);
//MemoLine(loc0);
//UnicodeString loc1 = "DrawInitial(): y0 = "
//	+ FloatToStr((ppLand.maxY - (loc.y+1)*Sp_resol_ratio) * p.pix);
//MemoLine(loc1);
//UnicodeString loc2 = "DrawInitial(): x1 = "
//	+ FloatToStr((loc.x+1) * p.sppix);
//MemoLine(loc2);
//UnicodeString loc3 = "DrawInitial(): y1 = "
//	+ FloatToStr((ppLand.maxY - (loc.y)*Sp_resol_ratio) * p.pix);
//MemoLine(loc3);
//}
#endif
		if (distCell.x >= 0) { // location is valid
				canvasSpDist->FrameRect(
					Rect((int)(distCell.x * Sp_resol_ratio * p.pix),
							 (int)(((distDim.y-distCell.y)*Sp_resol_ratio + diffY) * p.pix),
							 (int)((distCell.x+1) * Sp_resol_ratio * p.pix),
							 (int)(((distDim.y-distCell.y+1)*Sp_resol_ratio + diffY) * p.pix))
				);
		}
	}
}

frmMain->SpDistImage->Visible = true;
frmMain->SpDistImage->Picture->Bitmap = bmpSpDist;
//frmMain->Repaint();
frmMain->SpDistImage->Repaint();
}

//---------------------------------------------------------------------------
void DrawPatches(void)
{
intptr patchnum;
int *rs, *gs, *bs;
Patch *pPatch;
Cell *pCell;
landParams ppLand = pLandscape->getLandParams();
landPix p = pLandscape->getLandPix();
int npatches = pLandscape->patchCount();

#if RSDEBUG
//UnicodeString pmsg2 = "DrawPatches(): p.pix = " + FloatToStr(p.pix);
//MemoLine(pmsg2);
//UnicodeString pmsg3 = "DrawPatches(): npatches = " + FloatToStr((float)npatches);
//MemoLine(pmsg3);
#endif

rs = new int[npatches]; gs = new int[npatches]; bs = new int[npatches];

#if RSDEBUG
//UnicodeString pmsg0 = "DrawPatches(): npatches = " + IntToStr(npatches);
//MemoLine(pmsg0);
#endif
for (int i = 0; i < npatches; i++) {
	rs[i] = gs[i] = bs[i] = 0;
	do { // loop to avoid dark colours
		rs[i] = pRandom->IRandom(0,255);
		gs[i] = pRandom->IRandom(0,255);
		bs[i] = pRandom->IRandom(0,255);
#if RSDEBUG
//UnicodeString pmsg2 = "DrawPatches(): i = " + IntToStr(i)
//	+ " r = " + IntToStr(rs[i]) + " g = " + IntToStr(gs[i]) + " b = " + IntToStr(bs[i]);
//MemoLine(pmsg2);
#endif
	}
	while (rs[i] < 50 && gs[i] < 50 && bs[i] < 50);
//	while ((rs[i] + gs[i] + bs[i]) < 100);

}

for (int y = ppLand.maxY; y >= 0; y--) {
	for (int x = 0; x <= ppLand.maxX; x++) {
		pCell = pLandscape->findCell(x,y);
		if (pCell == 0) {
			canvasPatches->Brush->Color = (TColor)RGB(40,40,40); // grey
		}
		else {
			patchnum = pCell->getPatch();
			if (patchnum == 0) {
				canvasPatches->Brush->Color = (TColor)RGB(40,40,40); // grey
			}
			else {
				pPatch = (Patch*)patchnum;
				patchnum = pPatch->getSeqNum();
				if (patchnum == 0) { // matrix
					canvasPatches->Brush->Color = (TColor)RGB(0,0,0); // grey
				}
				else {
					// allow for patches not being consecutively numbered
					// SHOULD NOT NOW OCCUR, AS PATCH SEQUENTIAL NUMBER IS USED
					// (may result in replicated colours)
					if (patchnum > npatches) patchnum -= npatches;
					canvasPatches->Brush->Color = (TColor)RGB(rs[patchnum],gs[patchnum],bs[patchnum]);
				}
			}
		}
		canvasPatches->FillRect(Rect((int)(x*p.pix),(int)((ppLand.maxY-y)*p.pix),
			(int)((x+1)*p.pix),(int)((ppLand.maxY-y+1)*p.pix)));
	}
}

frmVisualPatch->PatchImage->Picture->Bitmap = bmpPatches;
frmVisualPatch->Repaint();
delete[] rs; delete[] bs; delete[] gs;
}

//---------------------------------------------------------------------------
traitCanvas SetupTraitCanvas(void) {
traitCanvas tcanv;
emigRules emig = pSpecies->getEmig();
trfrRules trfr = pSpecies->getTrfr();
settleType sett = pSpecies->getSettle();

int ncanvas = 0;

for (int i = 0; i < NTRAITS; i++) {
	tcanv.pcanvas[i] = 0;
}
if (emig.indVar) {
	if (emig.sexDep) {
		if (emig.densDep) ncanvas += 6; else ncanvas += 2;;
	}
	else {
		if (emig.densDep) ncanvas += 3; else ncanvas += 1;
	}
}
if (trfr.indVar) {
	if (trfr.moveModel) {
		if (trfr.moveType == 2) ncanvas += 2; // CRW
	}
	else { // kernel
		if (trfr.sexDep) {
			if (trfr.twinKern) ncanvas += 6; else ncanvas += 2;
		}
		else {
			if (trfr.twinKern) ncanvas += 3; else ncanvas += 1;
		}
	}
}
if (sett.indVar) ncanvas += 3;

for (int i = 0; i < ncanvas; i++) { tcanv.pcanvas[i] = canvasT[i]; }
#if RSDEBUG
//DEBUGLOG << "SetupTraitCanvas(): ncanvas=" << ncanvas
//	<< " canvasT[0]=" << (int)canvasT[0]
//	<< " canvasT[1]=" << (int)canvasT[1]
//	<< " tcanv.pcanvas[0]=" << (int)tcanv.pcanvas[0]
//	<< endl;
#endif

return tcanv;
}

//---------------------------------------------------------------------------
void Outputs_Visuals_B(int rep,int yr,int gen,int landNr)
{
emigRules emig = pSpecies->getEmig();
trfrRules trfr = pSpecies->getTrfr();
settleType sett = pSpecies->getSettle();
simParams sim = paramsSim->getSim();
simView v = paramsSim->getViews();
string DirMap = paramsSim->getDir(3);
string mapname;
int ixt = 0;
bool savemap;
if (sim.saveTraitMaps && yr%sim.traitInt == 0) savemap = true; else savemap = false;

#if RSDEBUG
//DEBUGLOG << "Outputs_Visuals_B(): rep=" << rep << " yr=" << yr
//	<< " gen=" << gen << " landNr=" << landNr << " ixt=" << ixt << endl;
#endif

if (emig.indVar) {
	RedrawTraitImage(ixt,v.viewGrad);
	if (savemap) {
		if (emig.sexDep) {
			if (emig.densDep) mapname = "D0_F"; else mapname = "EP_F";
		}
		else {
			if (emig.densDep) mapname = "D0"; else mapname = "EP";
		}
		SaveTraitMap(rep,yr,gen,landNr,DirMap,bmpT[ixt],mapname);
	}
	ixt++;
	if (emig.sexDep) {
		RedrawTraitImage(ixt,v.viewGrad);
		if (savemap) {
			if (emig.densDep) mapname = "D0_M"; else mapname = "EP_M";
			SaveTraitMap(rep,yr,gen,landNr,DirMap,bmpT[ixt],mapname);
		}
		ixt++;
		if (emig.densDep) {
			RedrawTraitImage(ixt,v.viewGrad);
			if (savemap) {
				mapname = "alpha_F";
				SaveTraitMap(rep,yr,gen,landNr,DirMap,bmpT[ixt],mapname);
			}
			ixt++;
			RedrawTraitImage(ixt,v.viewGrad);
			if (savemap) {
				mapname = "alpha_M";
				SaveTraitMap(rep,yr,gen,landNr,DirMap,bmpT[ixt],mapname);
			}
			ixt++;
			RedrawTraitImage(ixt,v.viewGrad);
			if (savemap) {
				mapname = "beta_F";
				SaveTraitMap(rep,yr,gen,landNr,DirMap,bmpT[ixt],mapname);
			}
			ixt++;
			RedrawTraitImage(ixt,v.viewGrad);
			if (savemap) {
				mapname = "beta_M";
				SaveTraitMap(rep,yr,gen,landNr,DirMap,bmpT[ixt],mapname);
			}
			ixt++;
		}
	}
	else { // !emig.sexDep
		if (emig.densDep) {
			RedrawTraitImage(ixt,v.viewGrad);
			if (savemap) {
				mapname = "alpha";
				SaveTraitMap(rep,yr,gen,landNr,DirMap,bmpT[ixt],mapname);
			}
			ixt++;
			RedrawTraitImage(ixt,v.viewGrad);
			if (savemap) {
				mapname = "beta";
				SaveTraitMap(rep,yr,gen,landNr,DirMap,bmpT[ixt],mapname);
			}
			ixt++;
		}
	}
}

if (trfr.indVar) {
	if (trfr.moveModel) {
		if (trfr.moveType == 2) { // CRW
			RedrawTraitImage(ixt,v.viewGrad);
			if (savemap) {
				mapname = "StepLength";
				SaveTraitMap(rep,yr,gen,landNr,DirMap,bmpT[ixt],mapname);
			}
			ixt++;
			RedrawTraitImage(ixt,v.viewGrad);
			if (savemap) {
				mapname = "Rho";
				SaveTraitMap(rep,yr,gen,landNr,DirMap,bmpT[ixt],mapname);
			}
			ixt++;
		}
	}
	else { // kernel
		RedrawTraitImage(ixt,v.viewGrad);
		if (savemap) {
			if (trfr.sexDep) {
				if (trfr.twinKern) mapname = "meanDistI_F"; else mapname = "meanDist_F";
			}
			else {
				if (trfr.twinKern) mapname = "meanDistI"; else mapname = "meanDist";
			}
			SaveTraitMap(rep,yr,gen,landNr,DirMap,bmpT[ixt],mapname);
		}
		ixt++;
		if (trfr.sexDep) {
			RedrawTraitImage(ixt,v.viewGrad);
			if (savemap) {
				if (trfr.twinKern) mapname = "meanDistI_M"; else mapname = "meanDist_M";
				SaveTraitMap(rep,yr,gen,landNr,DirMap,bmpT[ixt],mapname);
			}
			ixt++;
			if (trfr.twinKern) {
				RedrawTraitImage(ixt,v.viewGrad);
				if (savemap) {
					mapname = "meanDistII_F";
					SaveTraitMap(rep,yr,gen,landNr,DirMap,bmpT[ixt],mapname);
				}
				ixt++;
				RedrawTraitImage(ixt,v.viewGrad);
				if (savemap) {
					mapname = "meanDistII_M";
					SaveTraitMap(rep,yr,gen,landNr,DirMap,bmpT[ixt],mapname);
				}
				ixt++;
				RedrawTraitImage(ixt,v.viewGrad);
				if (savemap) {
					mapname = "meanProbKernelI_F";
					SaveTraitMap(rep,yr,gen,landNr,DirMap,bmpT[ixt],mapname);
				}
				ixt++;
				RedrawTraitImage(ixt,v.viewGrad);
				if (savemap) {
					mapname = "meanProbKernelI_M";
					SaveTraitMap(rep,yr,gen,landNr,DirMap,bmpT[ixt],mapname);
				}
				ixt++;
			}
		}
		else { // !trfr.sexDep
			if (trfr.twinKern) {
				RedrawTraitImage(ixt,v.viewGrad);
				if (savemap) {
					mapname = "meanDistII";
					SaveTraitMap(rep,yr,gen,landNr,DirMap,bmpT[ixt],mapname);
				}
				ixt++;
				RedrawTraitImage(ixt,v.viewGrad);
				if (savemap) {
					mapname = "meanProbKernelI";
					SaveTraitMap(rep,yr,gen,landNr,DirMap,bmpT[ixt],mapname);
				}
				ixt++;
			}
    }
	}
}
if (sett.indVar) {
	RedrawTraitImage(ixt,v.viewGrad);
	if (savemap) {
		mapname = "S0";
		SaveTraitMap(rep,yr,gen,landNr,DirMap,bmpT[ixt],mapname);
	}
	ixt++;
	RedrawTraitImage(ixt,v.viewGrad);
	if (savemap) {
		mapname = "alphaS";
		SaveTraitMap(rep,yr,gen,landNr,DirMap,bmpT[ixt],mapname);
	}
	ixt++;
	RedrawTraitImage(ixt,v.viewGrad);
	if (savemap) {
		mapname = "betaS";
		SaveTraitMap(rep,yr,gen,landNr,DirMap,bmpT[ixt],mapname);
	}
	ixt++;
}

frmVisualTraits0->Repaint();
if (ixt > 5) frmVisualTraits1->Repaint();
if (ixt > 11) frmVisualTraits2->Repaint();

#if RSDEBUG
//DEBUGLOG << "Outputs_Visuals_B(): finished " << endl;
#endif
}

//---------------------------------------------------------------------------

// Refresh visual cost screen
void RefreshVisualCost(void) {
#if RSDEBUG
//DEBUGLOG << "RefreshVisualCost(): running " << endl;
#endif
frmVisualCost->PaintBox1->Refresh();
frmVisualCost->Refresh();
}

//---------------------------------------------------------------------------

void SetGradTraitsForms(void)
{
int rs,gs,bs;
float lpix; // scaling factor for image labels (colour keys)
double minLabel,maxLabel;
int height,width,lgndht;
rgb col;
UnicodeString caption;
landParams ppLand = pLandscape->getLandParams();
landPix p = pLandscape->getLandPix();
envGradParams grad = paramsGrad->getGradient();
demogrParams dem = pSpecies->getDemogr();
emigRules emig = pSpecies->getEmig();
trfrRules trfr = pSpecies->getTrfr();
settleType sett = pSpecies->getSettle();
simView v = paramsSim->getViews();
int ixf = 0; // indexes form images
int ixt = 0; // indexes trait canvases
UnicodeString envgradtitle = "Environmental gradient ";
UnicodeString traits1title = "Mean dispersal traits (1) ";

#if RSDEBUG
//DEBUGLOG << "SetGradTraitsForms(): p.pix = " << p.pix << endl;
#endif

p.gpix = (float)def_grad_width / (float)(ppLand.maxX+1);
height = (int)((ppLand.maxY+1)*p.gpix);
width = frmVisualTraits0->Image0->Width;
lgndht = frmVisualTraits0->Image0L->Height;
frmVisualTraits0->refresh(height);
frmVisualTraits1->refresh(height);
frmVisualTraits2->refresh(height);

pLandscape->setLandPix(p);
#if RSDEBUG
MemoLine(("SetGradTraitsForms(): p.pix=" + Int2Str(p.pix)
	+ " p.gpix=" + Int2Str(p.gpix)
	).c_str());
#endif

lpix = (float)(width) / 768.0;

if (grad.gradient && v.viewGrad) {
	frmVisualTraits0->Caption = envgradtitle;
	if (emig.indVar || trfr.indVar || sett.indVar)
		frmVisualTraits0->Caption = envgradtitle + "& " + traits1title;

	if (bmpGrad == NULL) bmpGrad = new Graphics::TBitmap();
	canvasGrad = bmpGrad->Canvas;
	bmpGrad->Width = width; bmpGrad->Height = height;

	// SCFP 17/9/13
	// Legend for environmental gradient map
	// THIS DOES NOT WORK, AS draw_wheel() IS NOT SET UP FOR BLUE-RED RAMP...
	// NEED TO DISCUSS WITH GRETA...

	if (bmpGradL == NULL) bmpGradL = new Graphics::TBitmap(); canvasGradL = bmpGradL->Canvas;
	bmpGradL->Width = width; bmpGradL->Height = lgndht;
	UnicodeString caption;
	switch (grad.gradType) {
	case 1:
		caption = "Carrying capacity  ";
		maxLabel =
			((double)((int)(10*pSpecies->getMaxK()*10000.0)/(ppLand.resol*ppLand.resol)))/10.0;
		break;
	case 2:
		caption = "Intrinsic growth rate  ";
		maxLabel = ((double)((int)(10*pSpecies->getMaxFec())))/10.0;
		break;
	case 3:
		caption = "Local extinction probability  ";
		maxLabel = 1.0;
		break;
	}
	for (int i = 0; i < 768; i++) {
		bs = 255; rs = (int)(i*510/768);
		if (rs>255){ rs=255; bs = 510 -(int)(i*510/768); }
		canvasGradL->Brush->Color = (TColor)RGB(rs,0,bs);
		canvasGradL->FillRect(Rect((int)(i*lpix),0,(int)(i*lpix+lpix),10));
	}
	frmVisualTraits0->setImage(0,caption,maxLabel/2.0,maxLabel/6.0,NSD,0.0,maxLabel,bmpGradL);
	ixf++;
}
else {  // no gradient
	frmVisualTraits0->Caption = traits1title;
}

if (emig.indVar && v.viewTraits) { // set up images for emigration traits
	emigParams elim0,elim1;
	elim0 = pSpecies->getEmigParams(0,0); // females / only sex
	UnicodeString meanEP = "Mean Emigration Probability";
	if (emig.sexDep) {
		if (emig.densDep) caption = "Mean D0 (females) ";
		else caption = meanEP + " (females) ";
	}
	else {
		if (emig.densDep) caption = "Mean D0  ";
		else caption = meanEP;
	}
	SetTraitsImage(ixf++,ixt++,caption,height,width,lgndht,lpix,
		elim0.d0Mean,elim0.d0Scale,0.0,1.0);
	if (emig.sexDep) {
		elim1 = pSpecies->getEmigParams(0,1); // males
		if (emig.densDep) caption = "Mean D0 (males)  ";
		else caption = meanEP + " (males) ";
		SetTraitsImage(ixf++,ixt++,caption,height,width,lgndht,lpix,
			elim1.d0Mean,elim1.d0Scale,0.0,1.0);
		if (emig.densDep) {
			caption = "Mean alpha (females)  ";
			SetTraitsImage(ixf++,ixt++,caption,height,width,lgndht,lpix,
				elim0.alphaMean,elim0.alphaScale,1.0,0.0);
			caption = "Mean alpha (males) ";
			SetTraitsImage(ixf++,ixt++,caption,height,width,lgndht,lpix,
				elim1.alphaMean,elim1.alphaScale,1.0,0.0);
			caption = "Mean beta (females) ";
			SetTraitsImage(ixf++,ixt++,caption,height,width,lgndht,lpix,
				elim0.betaMean,elim0.betaScale,1.0,0.0);
			caption = "Mean beta (males) ";
			SetTraitsImage(ixf++,ixt++,caption,height,width,lgndht,lpix,
				elim1.betaMean,elim1.betaScale,1.0,0.0);
		}
	}
	else { // not sex-dependent
		if (emig.densDep) {
			caption = "Mean alpha ";
			SetTraitsImage(ixf++,ixt++,caption,height,width,lgndht,lpix,
				elim0.alphaMean,elim0.alphaScale,1.0,0.0);
			caption = "Mean beta ";
			SetTraitsImage(ixf++,ixt++,caption,height,width,lgndht,lpix,
				elim0.betaMean,elim0.betaScale,1.0,0.0);
		}
	}
}

if (trfr.indVar && v.viewTraits) { // set up images for transfer traits
	if (trfr.moveModel) {
#if RSDEBUG
//DEBUGLOG << "SetGradientForm(): trfr.moveType = " << trfr.moveType
//	<< " trfr.indVar = " << trfr.indVar << endl;
#endif
		if (trfr.moveType == 2) {
			trfrCRWParams mparams = pSpecies->getCRWParams(0,0);
			caption = "Mean step length (m) ";
			SetTraitsImage(ixf++,ixt++,caption,height,width,lgndht,lpix,
				mparams.stepLgthMean,mparams.stepLScale,1.0,1000000.0);
			caption = "Mean step correlation ";
			SetTraitsImage(ixf++,ixt++,caption,height,width,lgndht,lpix,
				mparams.rhoMean,mparams.rhoScale,0.0,1.0);
		}
	}
	else { // transfer by kernel
		float minDist,maxDist; // min and max distance limits to trait maps
		minDist = (float)ppLand.resol;
		maxDist = minDist * 100000.0; // effectively infinite
		if (pSpecies->useFullKernel()) minDist = 0.0;
		trfrKernParams k0,k1;
		k0 = pSpecies->getKernParams(0,0);
		UnicodeString meandd = "Mean dispersal distance ";
		if (trfr.sexDep) {
			if (trfr.twinKern) caption = meandd + "(1st kernel - females) (m) ";
			else caption = meandd + "(m) (females) ";
		}
		else {
			if (trfr.twinKern) caption = meandd + "(1st kernel) (m) ";
			else caption = meandd + "(m) ";
		}
		SetTraitsImage(ixf++,ixt++,caption,height,width,lgndht,lpix,
			k0.dist1Mean,k0.dist1Scale,minDist,maxDist);
		if (trfr.sexDep) {
			k1 = pSpecies->getKernParams(0,1);
			if (trfr.twinKern) caption = meandd + "(1st kernel - males (m) ";
			else caption = meandd + "(m) (males)";
			SetTraitsImage(ixf++,ixt++,caption,height,width,lgndht,lpix,
				k1.dist1Mean,k1.dist1Scale,minDist,maxDist);
			if (trfr.twinKern) {
				caption = meandd + "(2nd kernel - females) (m) ";
				SetTraitsImage(ixf++,ixt++,caption,height,width,lgndht,lpix,
					k0.dist2Mean,k0.dist2Scale,minDist,maxDist);
				caption = meandd + "(2nd kernel - males) (m) ";
				SetTraitsImage(ixf++,ixt++,caption,height,width,lgndht,lpix,
					k1.dist2Mean,k1.dist2Scale,minDist,maxDist);
				caption = "1st kernel probabilty (females) ";
				SetTraitsImage(ixf++,ixt++,caption,height,width,lgndht,lpix,
					k0.PKern1Mean,k0.PKern1Scale,0.0,1.0);
				caption = "1st kernel probabilty (males) ";
				SetTraitsImage(ixf++,ixt++,caption,height,width,lgndht,lpix,
					k1.PKern1Mean,k1.PKern1Scale,0.0,1.0);
			}
		}
		else { // sex-independent
			if (trfr.twinKern) {
				caption = meandd + "(2nd kernel) (m) ";
				SetTraitsImage(ixf++,ixt++,caption,height,width,lgndht,lpix,
					k0.dist2Mean,k0.dist2Scale,minDist,maxDist);
				caption = "1st kernel probabilty ";
				SetTraitsImage(ixf++,ixt++,caption,height,width,lgndht,lpix,
					k0.PKern1Mean,k0.PKern1Scale,0.0,1.0);
			}
		}
	}
}

if (sett.indVar && v.viewTraits) {
	settParams sparams = pSpecies->getSettParams(0,0);
	caption = "Mean S0 ";
	SetTraitsImage(ixf++,ixt++,caption,height,width,lgndht,lpix,
		sparams.s0Mean,sparams.s0Scale,0.0,1.0);
	caption = "Mean alphaS ";
	SetTraitsImage(ixf++,ixt++,caption,height,width,lgndht,lpix,
		sparams.alphaSMean,sparams.alphaSScale,1.0,0.0);
	caption = "Mean betaS ";
	SetTraitsImage(ixf++,ixt++,caption,height,width,lgndht,lpix,
		sparams.betaSMean,sparams.betaSScale,1.0,0.0);
}

if (ixf > 0) frmVisualTraits0->Show();
if (ixf > 6) frmVisualTraits1->Show();
if (ixf > 12) frmVisualTraits2->Show();

}

void RedrawTraitImage(const int ixt,bool gradient)
{
int ixf = ixt;
if (gradient) ixf++;
if (ixf < 6) {
	frmVisualTraits0->drawImage(ixf,bmpT[ixt]);
}
else {
	if (ixf < 12) {
		frmVisualTraits1->drawImage(ixf-6,bmpT[ixt]);
	}
	else {
		frmVisualTraits2->drawImage(ixf-12,bmpT[ixt]);
	}
}

}

void SetTraitsImage(const int ixf,const int ixt,
	const UnicodeString caption,const int h,const int w,const int lgdh,const float lpix,
	const float m,const float s,const float min,const float max)
{
int rs,gs,bs;
rgb col;

if (bmpT[ixt] == NULL) bmpT[ixt] = new Graphics::TBitmap(); canvasT[ixt] = bmpT[ixt]->Canvas;
bmpT[ixt]->Width = w; bmpT[ixt]->Height = h;

if (bmpL[ixt] == NULL) bmpL[ixt] = new Graphics::TBitmap(); canvasL[ixt] = bmpL[ixt]->Canvas;
bmpL[ixt]->Width = w; bmpL[ixt]->Height = lgdh;
for (int i = 0; i < 768; i++) {
	col = draw_wheel(i); rs = col.r; gs = col.g; bs = col.b;
	canvasL[ixt]->Brush->Color=(TColor)RGB(rs,gs,bs);
	canvasL[ixt]->FillRect(Rect((int)(i*lpix),0,(int)(i*lpix+lpix),10));
}
if (ixf < 6) {
	frmVisualTraits0->setImage(ixf,caption,m,s,NSD,min,max,bmpL[ixt]);
}
else {
	if (ixf < 12) {
		frmVisualTraits1->setImage(ixf-6,caption,m,s,NSD,min,max,bmpL[ixt]);
	}
	else {
		frmVisualTraits2->setImage(ixf-12,caption,m,s,NSD,min,max,bmpL[ixt]);
	}
}

}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

// FUNCTION DEFINITIONS FOR GUI VERSION ONLY

void Landscape::setLandMap(void) {
float hpix,wpix;

if (bmpLand == NULL) bmpLand = new Graphics::TBitmap();
canvasLand = bmpLand->Canvas;

// set the scale for drawing the landscape
#if RSDEBUG
//UnicodeString pmsg1 = "Landscape::drawLandscape() 1: bmpLand=" + IntToStr((int)bmpLand);
//MemoLine(pmsg1);
#endif
hpix = (float)def_land_height / (float)(dimY);
wpix = (float)def_land_width / (float)(dimX);
#if RSDEBUG
DebugGUI(("Landscape::setLandMap(): def_land_height=" + Int2Str(def_land_height)
	+ " dimY=" + Int2Str(dimY) + " hpix=" + Float2Str(hpix)
	+ " dimX=" + Int2Str(dimX) + " wpix=" + Float2Str(wpix)
	).c_str());
#endif
//if (wpix < hpix) hpix = wpix;
if (wpix > 1.0) pix = (int)wpix; else pix = 1;

bmpLand->Height = dimY * pix;
bmpLand->Width  = dimX * pix;
// try to keep map height within some reasonable maximum (arbitrarily 10000 pixels)
// otherwise saved maps can become too large
while (bmpLand->Height > 10000 && pix > 1) {
	pix /= 2;
	bmpLand->Height = dimY * pix;
	bmpLand->Width  = dimX * pix;
}
//if (bmpLand->Height < 1) bmpLand->Height = 1;
//if (bmpLand->Width < 1)  bmpLand->Width = 1;

//frmMain->LandImage->Height = bmpLand->Height;
//frmMain->LandImage->Width  = bmpLand->Width;
//frmMain->CommImage->Height = bmpLand->Height;
//frmMain->CommImage->Width  = bmpLand->Width;
//frmMain->SpDistImage->Height = bmpLand->Height;
//frmMain->SpDistImage->Width  = bmpLand->Width;
//frmMain->MovtPaintBox->Height = bmpLand->Height;
//frmMain->MovtPaintBox->Width  = bmpLand->Width;
//frmMain->LandScrollBox->VertScrollBar->Range = bmpLand->Height;
//frmMain->LandScrollBox->HorzScrollBar->Range = bmpLand->Width;
frmMain->SetMapDimensions();

#if RSDEBUG
DebugGUI(("Landscape::setLandMap(): nHab=" + Int2Str(nHab)
	+ " rasterType=" + Int2Str(rasterType)
//	+ " maxX=" + Int2Str(maxX) + " maxY=" + Int2Str(maxY)
	+ " dimX=" + Int2Str(dimX) + " dimY=" + Int2Str(dimY)
	+ " pix=" + Int2Str(pix)
	+ " bmpLand->Height=" + Int2Str(bmpLand->Height)
	+ " bmpLand->Width=" + Int2Str(bmpLand->Width)
	).c_str());
DebugGUI(("Landscape::setLandMap(): nHab=" + Int2Str(nHab)
	+ " LandImage->Height=" + Int2Str(frmMain->LandImage->Height)
	+ " LandImage->Width" + Int2Str(frmMain->LandImage->Width)
	).c_str());
DebugGUI(("Landscape::setLandMap(): nHab=" + Int2Str(nHab)
	+ " LandScrollBox->VertScrollBar->Range=" + Int2Str(frmMain->LandScrollBox->VertScrollBar->Range)
	+ " LandScrollBox->HorzScrollBar->Range=" + Int2Str(frmMain->LandScrollBox->HorzScrollBar->Range)
	).c_str());
#endif

}

// NOTE: do not use DEBUGLOG << here, as this function is run for an imported
// landscape before the stream is opened

void Landscape::drawLandscape(int rep,int landIx,int landnum)
{
#if RSDEBUG
DebugGUI(("Landscape::drawLandscape(): rep=" + Int2Str(rep)
	+ " nHab=" + Int2Str(nHab)
	+ " landIx=" + Int2Str(landIx)).c_str());
#endif
int hx; // habitat index
//float q,hpix,wpix;
float q;
int gscale;
locn loc;
MemoLine("Drawing landscape...");

rgb *colours = 0;
envGradParams grad = paramsGrad->getGradient();
simView v = paramsSim->getViews();

if (!generated && rasterType != 2) {
	// load colours into local table to increase speed of drawing cells
	colours = new rgb [nHab];
	for (int i = 0; i < nHab; i++) {
		colours[i] = getColour(i);
//	MemoLine(("Red " + Int2Str(i) + " = " + Int2Str(colours[i].r)).c_str());
	}
}

// NOTE: MAY BE MORE EFFICIENT TO EMBED rasterType CONDITION WITHIN FOR LOOPS ...

switch (rasterType) {

case 0: // habitat codes
	for (int y = dimY-1; y >= 0; y--)
	{
#if RSDEBUG
//DebugGUI(("Landscape::drawLandscape(): y=" + Int2Str(y)
//	+ " cells[y]=" + Int2Str((int)cells[y])).c_str());
#endif
		for (int x = 0; x < dimX; x++)
		{
			if (cells[y][x] == 0) { // no-data cell
				canvasLand->Brush->Color = (TColor)RGB(176,196,222);
			}
			else {
//				hx = cells[y][x]->getHabIndex();
				hx = cells[y][x]->getHabIndex(landIx);
#if RSDEBUG
//DebugGUI(("Landscape::drawLandscape(): y=" + Int2Str(y) + " x=" + Int2Str(x)
//	+ " hx=" + Int2Str(hx)).c_str());
#endif
				if (generated) {
					if (hx == 0) // matrix
						canvasLand->Brush->Color = (TColor)RGB(1,1,1); // black
					else // habitat
						canvasLand->Brush->Color = (TColor)RGB(255,255,255); // white
				}
				else { // imported landscape
						canvasLand->Brush->Color = (TColor)RGB(colours[hx].r,colours[hx].g,colours[hx].b);
//					}
				}
			}
#if RSDEBUG
//DebugGUI(("Landscape::drawLandscape(): y=" + Int2Str(y) + " x=" + Int2Str(x)
//	+ " cells[y][x]=" + Int2Str((int)cells[y][x]) + " hx=" + Int2Str(hx)).c_str());

//DebugGUI(("Landscape::drawLandscape(): y=" + Int2Str(y) + " x=" + Int2Str(x)
//	+ " hx=" + Int2Str(hx)
//	+ " x*pix=" + Int2Str((int)x*pix)
//	+ " (x+1)*pix=" + Int2Str((int)(x+1)*pix)
//	+ " (maxY-y)*pix=" + Int2Str((int)(maxY-y)*pix)
//	+ " (maxY-y+1)*pix=" + Int2Str((int)(maxY-y+1)*pix)
//	).c_str());
#endif
			canvasLand->FillRect(Rect((int)(x*pix),(int)((dimY-1-y)*pix),
				(int)((x+1)*pix),(int)((dimY-y)*pix)));
		}
	}
	break;

case 1: // % cover
	for (int y = dimY-1; y >= 0; y--) {
		for (int x = 0; x < dimX; x++) {
#if RSDEBUG
DebugGUI(("Landscape::drawLandscape(): y=" + Int2Str(y) + " x=" + Int2Str(x)
	+ " cells[y][x]=" + Int2Str((int)cells[y][x])
	).c_str());
#endif
			if (cells[y][x] == 0) { // no-data cell
				canvasLand->Brush->Color = (TColor)RGB(176,196,222);
			}
			else {
				hx = 0;
				float maxcover = 0.0;
				for (int i = 0; i < nHab; i++) {
					q = cells[y][x]->getHabitat(i);
					if (q > maxcover) {
						hx = i; maxcover = q;
					}
#if RSDEBUG
DebugGUI(("Landscape::drawLandscape(): y=" + Int2Str(y) + " x=" + Int2Str(x)
	+ " i=" + Int2Str(i)
	+ " q=" + Float2Str(q)
	+ " maxcover=" + Float2Str(maxcover)
	+ " hx=" + Int2Str(hx)).c_str());
#endif
					if (maxcover > 0.0) {
						canvasLand->Brush->Color = (TColor)RGB(colours[hx].r,colours[hx].g,colours[hx].b);
					}
					else { // no cover in any layer, i.e. matrix
						canvasLand->Brush->Color = (TColor)RGB(0,0,0); // black
					}
				}
			}
			canvasLand->FillRect(Rect((int)(x*pix),(int)((dimY-1-y)*pix),
				(int)((x+1)*pix),(int)((dimY-y)*pix)));
		}
	}
	break;

case 2: // habitat quality
	for (int y = dimY-1; y >= 0; y--)
	{
		for (int x = 0; x < dimX; x++)
		{
			if (cells[y][x] == 0) { // no-data cell
				canvasLand->Brush->Color = (TColor)RGB(176,196,222);
			}
			else {
//				q = cells[y][x]->getQuality();
				q = cells[y][x]->getHabitat(landIx);
				if (q > 0.0) {
					gscale = (int)(55 + q*2.0);
					if (gscale > 255) gscale = 255;
				}
				else // matrix in artificial continuous landscape
					gscale = 0; // black
				canvasLand->Brush->Color = (TColor)RGB(gscale,gscale,gscale);
			}
			canvasLand->FillRect(Rect((int)(x*pix),(int)((dimY-1-y)*pix),
				(int)((x+1)*pix),(int)((dimY-y)*pix)));
		}
	}
	break;
}

if (colours != 0) delete colours;

if (v.viewLand) {
	frmMain->LandImage->Picture->Bitmap = bmpLand;
	frmMain->LandImage->Repaint();
}
MemoLine("...finished");

}

// Draw environmental gradient map
void Landscape::drawGradient(void)
{
float envval,val;
int bs,rs;
int gtype;
landParams ppLand = pLandscape->getLandParams();
landPix p = pLandscape->getLandPix();
envGradParams grad = paramsGrad->getGradient();
envStochParams env = paramsStoch->getStoch();
demogrParams dem = pSpecies->getDemogr();

for (int y = dimY-1; y >= 0; y--) {
	for (int x = 0; x < dimX; x++) {
		if (cells[y][x] == 0) { // no-data area - force black
			gtype = 99;
		}
		else {
			envval = cells[y][x]->getEnvVal();
			if (grad.gradType == 3) { // extinction probability
				// need to distinguish between a matrix cell and a cell having zero envval
				// which implies extinction probability of 1
//				int hab = cells[y][x]->getHabIndex();
				int hab = cells[y][x]->getHabIndex(0);
				if (hab < 0) gtype = 99; // no-data area - force black
				else {
					if (pSpecies->getHabK(hab) > 0.0) gtype = grad.gradType;
					else gtype = 99; // matrix - force black
				}
			}
			else {
				if (envval <= 0.0) gtype = 99; // force black
				else {
					gtype = grad.gradType;
				}
			}
		}
		switch (gtype) {
		case 1: // K
			val = envval;
			if (env.stoch && env.local && env.inK) {
				val += cells[y][x]->getEps();
			}
			bs = 255; rs = (int)(510*val); if (rs > 255){ bs = 510-rs; rs = 255; }
			if (val > 0.0) canvasGrad->Brush->Color=(TColor)RGB(rs,0,bs);
			else canvasGrad->Brush->Color = (TColor)RGB(0,0,0);
//			if (val > 0.0) canvasT[0]->Brush->Color=(TColor)RGB(rs,0,bs);
//			else canvasT[0]->Brush->Color = (TColor)RGB(0,0,0);
			break;
		case 2: // r
			val = envval;
			if (env.stoch && env.local && !env.inK) {
				val += cells[y][x]->getEps();
			}
			bs = 255; rs = (int)(510*val); if (rs > 255){ bs = 510-rs; rs = 255; }
			if (val > 0.0) canvasGrad->Brush->Color=(TColor)RGB(rs,0,bs);
			else canvasGrad->Brush->Color = (TColor)RGB(0,0,0);
//			if (val > 0.0) canvasT[0]->Brush->Color=(TColor)RGB(rs,0,bs);
//			else canvasT[0]->Brush->Color = (TColor)RGB(0,0,0);
			break;
		case 3: // extinction probability
			val = 1.0 - envval + grad.extProbOpt;
			if (val > 1.0) val = 1.0;
			bs = 255; rs = (int)(510*val); if (rs > 255){ bs = 510-rs; rs = 255; }
			if (val < 1.0) canvasGrad->Brush->Color=(TColor)RGB(rs,0,bs);
			else canvasGrad->Brush->Color = (TColor)RGB(0,0,0);
//			if (val < 1.0) canvasT[0]->Brush->Color=(TColor)RGB(rs,0,bs);
//			else canvasT[0]->Brush->Color = (TColor)RGB(0,0,0);
			break;
		default:
			canvasGrad->Brush->Color = (TColor)RGB(0,0,0);
//			canvasT[0]->Brush->Color = (TColor)RGB(0,0,0);
			;
		}
//		canvasGrad->FillRect(Rect((int)(x*p.gpix),(int)((ppLand.maxY-y)*p.gpix),
//			(int)((x+1)*p.gpix),(int)((ppLand.maxY-y+1)*p.gpix)));
		canvasGrad->FillRect(Rect((int)(x*p.gpix),(int)((ppLand.dimY-1-y)*p.gpix),
			(int)((x+1)*p.gpix),(int)((ppLand.dimY-y)*p.gpix)));
//		canvasT[0]->FillRect(Rect((int)(x*p.gpix),(int)((ppLand.maxY-y)*p.gpix),
//			(int)((x+1)*p.gpix),(int)((ppLand.maxY-y+1)*p.gpix)));
	}
}

frmVisualTraits0->Image0->Picture->Bitmap = bmpGrad;
frmVisualTraits0->Repaint();

}

// Draw environmental stochasticity time-series
void Landscape::drawGlobalStoch(int nyears){
for (int i = 0; i < nyears; i++){
	frmMain->EnvNoise->AddXY(i,epsGlobal[i],"");
}
}

//---------------------------------------------------------------------------

// Save SMS path visits map to .bmp file
void Landscape::saveVisits(int rep, int landNr)
{

int rs,gs;
int nFactor;
unsigned long int visits,maxvisits;
string mapName;
simParams sim = paramsSim->getSim();

// find max no of visits to scale colour ramp;
maxvisits = 0;
for (int y = dimY-1; y >= 0; y--) {
	for (int x = 0; x < dimX; x++) {
		if (cells[y][x] != 0) { // not a no-data cell
			visits = cells[y][x]->getVisits();
			if (visits > maxvisits) maxvisits = visits;
		}
	}
}
if (maxvisits < 1) return;

// set up new bitmap for the map
Graphics::TBitmap *mapVisits = new Graphics::TBitmap();
TCanvas *canvasVisits = mapVisits->Canvas;
mapVisits->Height = bmpLand->Height;
mapVisits->Width  = bmpLand->Width;
mapVisits->Canvas->CopyRect(Rect(0,0,mapVisits->Width,mapVisits->Height),
	bmpLand->Canvas,Rect(0,0,bmpLand->Width,bmpLand->Height));

// draw visits on red-blue colour ramp
for (int y = dimY-1; y >= 0; y--) {
	for (int x = 0; x < dimX; x++) {
		if (cells[y][x] != 0) { // not a no-data cell
			visits = cells[y][x]->getVisits();
			nFactor = 100 + (int) (visits * 280.0 / maxvisits);
			if (visits > 0.0) {
				if (nFactor <= 256) { rs = nFactor; gs = 0; } // red increasing, green 0
				else {
					if (nFactor <= 380) {
						rs = 255; gs = nFactor - 255; // red 255, green increasing until 380
					}
					else { rs = 255; gs = 175; }
				}
				canvasVisits->Brush->Color=(TColor)RGB(rs,gs,0);
				canvasVisits->FillRect(Rect((int)(x*pix),(int)((dimY-1-y)*pix),
					(int)((x+1)*pix),(int)((dimY-y)*pix)));
			}
		}
	}
}

if (sim.batchMode) {
	mapName = paramsSim->getDir(3)
		+ "Batch" + Int2Str(sim.batchNum) + "_"
		+ "Sim" + Int2Str(sim.simulation)
		+ "_land" + Int2Str(landNr) + "_rep" + Int2Str(rep)
		+ "_Visits.bmp";
}
else {
	mapName = paramsSim->getDir(3)
		+ "Sim" + Int2Str(sim.simulation)
		+ "_land" + Int2Str(landNr) + "_rep" + Int2Str(rep)
		+ "_Visits.bmp";
}
mapVisits->SaveToFile(mapName.c_str());
delete mapVisits;

}

//---------------------------------------------------------------------------

void Patch::drawCells(TCanvas *pCanvas,float pix,int dimY,rgb col) {
locn loc;

pCanvas->Brush->Color=(TColor)RGB(col.r,col.g,col.b);
int ncells = (int)cells.size();
for (int i = 0; i < ncells; i++) {
	loc = cells[i]->getLocn();
#if RSDEBUG
//DEBUGLOG << "Patch::drawCells(): pCanvas=" << pCanvas
//	<< " dimY=" << dimY << " loc.x=" << loc.x << " loc.y=" << loc.y
//	<< endl;
#endif
	pCanvas->FillRect(Rect((int)(loc.x*pix),(int)((dimY-1-loc.y)*pix),
		(int)((loc.x+1)*pix),(int)((dimY-loc.y)*pix)));
}
}

//---------------------------------------------------------------------------
// Draw the community on the landscape map
void Community::draw(int rep, int yr, int gen, int landNr)
{
simParams sim = paramsSim->getSim();
simView v = paramsSim->getViews();
#if RSDEBUG
//DEBUGLOG << "Community::draw(): rep = " <<  rep << "  yr = " << yr << endl;
#endif

if (bmpComm != NULL) delete bmpComm;
bmpComm = NULL;
bmpComm = new Graphics::TBitmap();
bmpComm->Height = bmpLand->Height;
bmpComm->Width  = bmpLand->Width;
//canvasComm = bmpComm->Canvas;
// first paint entire canvas black to be used as transparent colour
//canvasComm->Brush->Color=(TColor)RGB(0,0,0);
//canvasComm->FillRect(Rect(0,0,bmpComm->Width,bmpComm->Height));
bmpComm->Canvas->Brush->Color = (TColor)RGB(0,0,0);
bmpComm->Canvas->FillRect(Rect(0,0,bmpComm->Width,bmpComm->Height));
bmpComm->Transparent = true;
bmpComm->TransparentColor = (TColor)RGB(0,0,0);

// generate output for each sub-community (patch) in the community
#if RSDEBUG
//DEBUGLOG << "Community::draw(): bmpComm->Height = " << bmpComm->Height << endl;
#endif
int ncomms = (int)subComms.size();
for (int i = 0; i < ncomms; i++) { // all sub-communities
//	subComms[i]->draw(canvasComm,pLandscape);
	subComms[i]->draw(bmpComm->Canvas,pLandscape);
}

//if (v.viewLand) {
//	frmMain->LandImage->Picture->Bitmap = bmpLand;
//	frmMain->Repaint();
//}
if (v.viewPop) {
	frmMain->CommImage->Picture->Bitmap->Height = frmMain->LandImage->Picture->Bitmap->Height;
	frmMain->CommImage->Picture->Bitmap->Width  = frmMain->LandImage->Picture->Bitmap->Width;
	frmMain->CommImage->Picture->Bitmap = bmpComm;
	frmMain->CommImage->Transparent = true;
	frmMain->CommImage->Repaint();
}

if (sim.saveMaps && yr%sim.mapInt == 0) {
	// save to file compound map of community superimposed on landscape
	string mapName;
	locn distCell;
	landParams ppLand = pLandscape->getLandParams();

	// set up new bitmap for the compound map
	Graphics::TBitmap *mapI = new Graphics::TBitmap();
	TCanvas *canvasI = mapI->Canvas;
	mapI->Height = bmpLand->Height;
	mapI->Width  = bmpLand->Width;
	mapI->Canvas->CopyRect(Rect(0,0,mapI->Width,mapI->Height),
		bmpLand->Canvas,Rect(0,0,bmpLand->Width,bmpLand->Height));

	// add the community
#if RSDEBUG
//DEBUGLOG << "Community::draw(): mapI->Height = " << mapI->Height << endl;
#endif
	for (int i = 0; i < ncomms; i++) { // all sub-communities
		subComms[i]->draw(canvasI,pLandscape);
//		subComms[i]->draw(mapI->Canvas,pLandscape);
	}

	if (ppLand.spDist && sim.drawLoaded) { // add initial distribution
		canvasI->Brush->Color = (TColor)RGB(255,255,0); // yellow
		locn distDim = pLandscape->getDistnDimensions(0);
		int diffY = (ppLand.dimY) - (distDim.y+1) * Sp_resol_ratio;
		landPix p = pLandscape->getLandPix();
		int nDistns = pLandscape->distnCount();
		if (nDistns > 0) {
			int nCells = pLandscape->distCellCount(0);
			for (int i = 0; i < nCells; i++) {
				distCell = pLandscape->getDistnCell(0,i);
				if (distCell.x >= 0) { // cell is valid
					canvasI->FrameRect(
						Rect((int)(distCell.x * Sp_resol_ratio * p.pix),
								 (int)(((distDim.y-distCell.y)*Sp_resol_ratio + diffY) * p.pix),
								 (int)((distCell.x+1) * Sp_resol_ratio * p.pix),
							 (int)(((distDim.y-distCell.y+1)*Sp_resol_ratio + diffY) * p.pix))
					);
				}
			}
		}
	}

#if RSDEBUG
//DEBUGLOG << "Community::draw(): saving map..." << endl;
#endif
	mapName = paramsSim->getDir(3) + "Sim" + Int2Str(sim.simulation)
		+ "_land" + Int2Str(landNr) + "_rep" + Int2Str(rep)
		+ "_yr" + Int2Str(yr) + "_rs" + Int2Str(gen) + ".bmp";
	mapI->SaveToFile(mapName.c_str());
	delete mapI;
}

#if RSDEBUG
//DEBUGLOG << "Community::draw(): finished" << endl;
#endif
}

//---------------------------------------------------------------------------

void SubCommunity::draw(TCanvas *pCanvas,Landscape *pLandscape) {
#if RSDEBUG
//DEBUGLOG << "SubCommunity::draw(): pLandscape = " <<  pLandscape
//	<< "  subCommNum = " << subCommNum << endl;
#endif
popStats pop;
int ninds = 0;
landPix p = pLandscape->getLandPix();
landData land = pLandscape->getLandData();
rgb colour;
Species *pSpecies;
double density,maxDensity;

if (subCommNum == 0) // matrix sub-community
	return;

// NOTE: MAX DENSITY SHOULD BE SET OUTWITH THIS FUNCTION, AND REMAIN CONSTANT FOR THE
// WHOLE SIMULATION ...
maxDensity = 0.0;

int npops = (int)popns.size();
for (int i = 0; i < npops; i++) { // all populations
	pop = popns[i]->getStats();
	ninds += pop.nInds;
	pSpecies = popns[i]->getSpecies();
	density = pSpecies->getMaxK() * 10000 / (double)(land.resol*land.resol);
	if (density > maxDensity) maxDensity = density;
}
int ncells = pPatch->getNCells();
// NB no. of cells could be zero in a dynamic landscape if patch has been destroyed
#if RSDEBUG
//DEBUGLOG << "SubCommunity::draw(): ncells=" <<  ncells
//	<< " npops=" << npops << " ninds=" << ninds << endl;
#endif

if (ninds > 0 && ncells > 0) {
	double density = ((double)ninds * 10000.0 /
			 (double)(ncells * land.resol * land.resol));
	int nFactor = 100 + (int) (density * 280.0 / maxDensity);

	if (nFactor <= 256) { // red increasing, green 0
		colour.r = nFactor; colour.g = 0;
	}
	else {
		if (nFactor <= 380)
		{ // red 255, green increasing until 380
			colour.r = 255; colour.g = nFactor - 255;
		}
		else {
			colour.r = 255; colour.g = 175;
		}
	}
	colour.b = 0;
#if RSDEBUG
//DEBUGLOG << "SubCommunity::draw(): density = " << density
//	<< " maxDensity = " << maxDensity
//	<< " nFactor = " << nFactor
//	<< " r = " << colour.r
//	<< " g = " << colour.g
//	<< " b = " << colour.b
//	<< endl;
#endif
	pPatch->drawCells(pCanvas,p.pix,land.dimY,colour);
}
//else {
//	colour.r = colour.g = colour.b = 0;
//	pPatch->drawCells(canvas,p.pix,land.maxY,colour);
//}

}

//---------------------------------------------------------------------------

// Visualise paths resulting from movement simulation model
void Individual::drawMove(const float x0,const float y0,const float x1,const float y1)
{
landParams ppLand = pLandscape->getLandParams();
demogrParams dem = pSpecies->getDemogr();
simView v = paramsSim->getViews();
landPix p = pLandscape->getLandPix();

if (v.viewCosts) {
	if (dem.repType == 0) // asexual model
		frmVisualCost->PaintBox1->Canvas->Pen->Color = (TColor)RGB(200,0,200);
	else{
		if (sex == 1) frmVisualCost->PaintBox1->Canvas->Pen->Color = (TColor)RGB(0,0,200);
		else frmVisualCost->PaintBox1->Canvas->Pen->Color = (TColor)RGB(200,0,200);
	}
//	frmVisualCost->PaintBox1->Canvas->MoveTo((int)(x0*p.pix),(int)((ppLand.maxY-y0+1)*p.pix));
//	frmVisualCost->PaintBox1->Canvas->LineTo((int)(x1*p.pix),(int)((ppLand.maxY-y1+1)*p.pix));
	frmVisualCost->PaintBox1->Canvas->MoveTo((int)(x0*p.pix),(int)((ppLand.dimY-y0)*p.pix));
	frmVisualCost->PaintBox1->Canvas->LineTo((int)(x1*p.pix),(int)((ppLand.dimY-y1)*p.pix));
}
else {
	if (dem.repType == 0) // asexual model
		frmMain->MovtPaintBox->Canvas->Pen->Color = (TColor)RGB(200,0,200);
	else{
		if (sex == 1) frmMain->MovtPaintBox->Canvas->Pen->Color = (TColor)RGB(0,0,200);
		else frmMain->MovtPaintBox->Canvas->Pen->Color = (TColor)RGB(200,0,200);
	}
//	frmMain->MovtPaintBox->Canvas->MoveTo((int)(x0*p.pix),(int)((ppLand.maxY-y0+1)*p.pix));
//	frmMain->MovtPaintBox->Canvas->LineTo((int)(x1*p.pix),(int)((ppLand.maxY-y1+1)*p.pix));
	frmMain->MovtPaintBox->Canvas->MoveTo((int)(x0*p.pix),(int)((ppLand.dimY-y0)*p.pix));
	frmMain->MovtPaintBox->Canvas->LineTo((int)(x1*p.pix),(int)((ppLand.dimY-y1)*p.pix));
}

}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


