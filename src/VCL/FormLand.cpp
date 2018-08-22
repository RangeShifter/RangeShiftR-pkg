//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "FormLand.h"
#include "FormVisualCost.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TfrmLand *frmLand;

vector <string> hfnames; // list of additional habitat file names for the % cover option
string habmapname;
string patchmapname;
string distnmapname;
int hab;
string msgResol = "The raster resolution is different from the set landscape resolution";
string msgHeaders = "The origin and/or dimensions of the current habitat\n"
	"file differ from those of the first habitat file";
string msgPatchHdr = "Patch raster dimensions differ from habitat raster dimensions";
string msgLandHdr = "Format error in header lines of the selected landscape file";
string msgReadError = "Error returned from reading ";

int RED[21]   = {  0,250,200,100,200,150,153,155,128,230,  0,  0,  0,  0,200,60,0,204,255,128,  0};
int GREEN[21] = {200,200,200,250,150,150,128,100, 26,140,100,128,  0,180,200,60,0,179,255,102,  0};
int BLUE[21]  = { 50,150,100,100,250,150,  0, 60,128,166,  0,115,255,190,200,60,0,  0,128,255,128};

//---------------------------------------------------------------------------
__fastcall TfrmLand::TfrmLand(TComponent* Owner)
	: TForm(Owner)
{
frmLand->Left = 100;
frmLand->Top = 60;
}
//---------------------------------------------------------------------------
void __fastcall TfrmLand::FormCreate(TObject *Sender)
{
// initialise form
RGhab->ItemIndex = 0;
RGCellPatch->ItemIndex = 0;
StringGridHab->Cells[1][0] = "R";
StringGridHab->Cells[2][0] = "G";
StringGridHab->Cells[3][0] = "B";
BtnOK->Enabled = false;
BtnCancel->Enabled = true;
}

//---------------------------------------------------------------------------
void __fastcall TfrmLand::BtnCancelClick(TObject *Sender)
{
frmLand->Close();
}

//---------------------------------------------------------------------------
void __fastcall TfrmLand::refresh(bool landscapeLoaded)
{
landParams ppLand = pLandscape->getLandParams();
if (landscapeLoaded && ppLand.rasterType < 2) {
	displayColourTable(ppLand.nHab);
	StringGridHab->Visible = true;
	StringGridHab->Enabled = true;
}
else {
	StringGridHab->Visible = false;
	StringGridHab->Enabled = false;
}
}

//------------------------------------------------------------------------------
void __fastcall TfrmLand::edtNhabitatExit(TObject *Sender)
{
hab = StrToInt(edtNhabitat->Text);
//setHabitatTable(hab);
}
//---------------------------------------------------------------------------
void __fastcall TfrmLand::RGhabClick(TObject *Sender)
{
switch (RGhab->ItemIndex) {
case 0:
	edtNhabitat->Enabled = false; edtNhabitat->Visible = false;
	LabelNhab->Visible = false;
	hab = 1;
	break;

case 1:
	edtNhabitat->Enabled = true; edtNhabitat->Visible = true;
	LabelNhab->Visible = true;
	hab = StrToInt(edtNhabitat->Text);
	break;

case 2:
	hab = 1;
	edtNhabitat->Enabled = false; edtNhabitat->Visible = false;
	LabelNhab->Visible = false;
	StringGridHab->Enabled = false; StringGridHab->Visible = false;
	BtnChangeColours->Enabled = false; BtnChangeColours->Visible = false;
	break;
}
}
//---------------------------------------------------------------------------
void __fastcall TfrmLand::RGCellPatchClick(TObject *Sender)
{
if (RGCellPatch->ItemIndex == 0) { // cell-based model
	CBVisualPatch->Enabled = false;
	CBVisualPatch->Visible = false;
  CBVisualPatch->Checked = false;
}
else { // patch-based model
	CBVisualPatch->Enabled = true;
	CBVisualPatch->Visible = true;
}
}

//------------------------------------------------------------------------------
void __fastcall TfrmLand::BtnImportClick(TObject *Sender)
{
string name,msg;
int m;
int imported = 1;
rasterdata landraster,patchraster;
simParams sim = paramsSim->getSim();
hfnames.clear(); // clear list of habitat file names for the % cover option

StringGridHab->Enabled = false; StringGridHab->Visible = false;
BtnChangeColours->Enabled = false; BtnChangeColours->Visible = false;
BtnDynamic->Enabled = false;
BtnOK->Enabled = false;
BtnCancel->Enabled = true;

int *setRED;
int *setGREEN;
int *setBLUE;

if (pLandscape != NULL) {
	delete pLandscape;
	landscapeLoaded = false;
	//	sim.initDistLoaded = false;
//	paramsSim->setSim(sim);
}
pLandscape = new Landscape();
//string msgResol = "The raster resolution is different from the set landscape resolution.";

landParams ppLand = pLandscape->getLandParams();
simView v = paramsSim->getViews();

ppLand.generated = 0; //imported landscape
ppLand.resol = -1;
ppLand.resol = StrToInt(frmLand->edtRes->Text);
if (ppLand.resol < 1) {
	MessageDlg("Resolution must be greater than zero",
		mtError, TMsgDlgButtons() << mbOK,0);
	return;
}
ppLand.rasterType = frmLand->RGhab->ItemIndex;
ppLand.nHab = hab;

v.viewLand = true;
if (RGCellPatch->ItemIndex == 0) { // cell-based model
	ppLand.patchModel = false;
}
else { // patch-based model
	ppLand.patchModel = true;
	v.viewPatch = CBVisualPatch->Checked;
}
paramsSim->setViews(v);

//string msg;
//msg = "BEFORE WRITE: generated " + Int2Str(ppLand.generated)
//	+ " RasterType " + Int2Str(ppLand.RasterType);
//MessageDlg(msg.c_str(),mtWarning, TMsgDlgButtons() << mbOK,0);

pLandscape->setLandParams(ppLand,false);

//string msgRead = "Before reading landscape nHab = " + Int2Str(ppLand.nHab);
//MessageDlg(msgRead.c_str(),mtWarning, TMsgDlgButtons() << mbOK,0);

landraster.ok = false;
int nloaded = 0;

// Read the only or first habitat raster

// WHY DOES THE DIALOG SCREEN SHOW THE WRONG FOLDER WHEN A NEW DIRECTORY
// HAS BEEN SET SECOND TIME ROUND?...
//	MemoLine(("btnImportClick(): Path = " + paramsSim->getDir(0)).c_str());

OpenTextFileDialog1->InitialDir = paramsSim->getDir(1).c_str();
//	MemoLine(("btnImportClick(): Path = " + paramsSim->getDir(1)).c_str());
msg = "Select raster map";
if (ppLand.rasterType == 1 && ppLand.nHab > 1) // % cover and multiple layers
	msg += " for habitat nr. 1";
OpenTextFileDialog1->Title = msg.c_str();
if (OpenTextFileDialog1->Execute()) {
	habmapname = AnsiString(OpenTextFileDialog1->FileName).c_str();
	landraster = CheckRasterFile(habmapname);
	if (landraster.ok) {
		if (landraster.cellsize != ppLand.resol) {
			landraster.ok = false;
			MessageDlg(msgResol.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
		}
	}
	else {
		MessageDlg(msgLandHdr.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
	}
}
else {
	landraster.ok = false;
	MessageDlg("Selection of habitat raster file cancelled",
		mtWarning, TMsgDlgButtons() << mbOK,0);
}
if (!landraster.ok) return;

//	dialogmsg;
//	dialogmsg = "Return code = " + Int2Str(getfile) + " habmapname = " + habmapname;
//	MessageDlg(dialogmsg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);

bool otherRastersOK = true;
if (ppLand.rasterType == 1 && ppLand.nHab > 1) { // % cover and multiple layers
	// read remaining habitat layers
	string fname;
	for (int i = 1; i < ppLand.nHab; i++) {
//		string msgRead1 = "% cover for loop, i = " + Int2Str(i);
//		MessageDlg(msgRead1.c_str(),mtWarning, TMsgDlgButtons() << mbOK,0);
		OpenTextFileDialog1->InitialDir = paramsSim->getDir(1).c_str();
		msg = "Select raster map for habitat nr. " + Int2Str(i+1);
		OpenTextFileDialog1->Title = msg.c_str();
		if (OpenTextFileDialog1->Execute()){
			fname = AnsiString(OpenTextFileDialog1->FileName).c_str();
			rasterdata raster2;
			raster2 = CheckRasterFile(fname);
			if (raster2.ok) {
				if (raster2.cellsize == ppLand.resol) {
					if (raster2.ncols == landraster.ncols
					&&  raster2.nrows == landraster.nrows
					&&  raster2.xllcorner == landraster.xllcorner
					&&  raster2.yllcorner == landraster.yllcorner) {
						otherRastersOK = true;
						hfnames.push_back(fname);
					}
					else {
						otherRastersOK = false;
						MessageDlg(msgHeaders.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
					}
				}
				else {
					otherRastersOK = false;
					MessageDlg(msgResol.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
				}
			}
			else {
				otherRastersOK = false;
				MessageDlg(msgLandHdr.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
			}
			if (!otherRastersOK) {
				imported = 1;
				i = ppLand.nHab+1; break; // to break out of habitats loop
			}
		}
		else {
			// file selection dialog has been cancelled
			MessageDlg("The landscape file selection has been cancelled",
						mtWarning, TMsgDlgButtons() << mbOK,0);
			otherRastersOK = false;
			i = ppLand.nHab+1; break; // to break out of habitats loop
		}
  }
}
if (!otherRastersOK) { hfnames.clear(); return; }

if (landraster.ok && otherRastersOK && ppLand.patchModel) {
	// get the name of the patch raster
	msg = "Habitat raster";
	if (ppLand.rasterType == 1 && ppLand.nHab > 1) // % cover and multiple layers
		msg += "s";
	msg += " OK - now specify patch IDs file ...";
	MessageDlg(msg.c_str(),mtWarning, TMsgDlgButtons() << mbOK,0);
	patchraster.ok = false;
	OpenTextFileDialog1->InitialDir = paramsSim->getDir(1).c_str();
	OpenTextFileDialog1->Title = "Select patch raster map";
	if (OpenTextFileDialog1->Execute()){
		patchmapname = AnsiString(OpenTextFileDialog1->FileName).c_str();
		patchraster = CheckRasterFile(patchmapname);
		if (patchraster.ok) {
			if (patchraster.cellsize != ppLand.resol) {
				patchraster.ok = false;
				MessageDlg(msgResol.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
			}
			if (patchraster.nrows != landraster.nrows || patchraster.ncols != landraster.ncols) {
				patchraster.ok = false;
				MessageDlg(msgPatchHdr.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
			}
		}
		else {
			MessageDlg(msgLandHdr.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
		}
	}
	else {
		MessageDlg("Selection of patch raster file cancelled",
			mtWarning, TMsgDlgButtons() << mbOK,0);
	}
}

if (!landraster.ok || !otherRastersOK) return;

// Read the landscape file(s) if checks have been passed
if (landraster.ok) {
	if (ppLand.patchModel) {
		if (patchraster.ok) {
			MemoLine("Reading landscape files...");
			imported = pLandscape->readLandscape(0,habmapname,patchmapname);
		}
		else return;
	}
	else {
		MemoLine("Reading landscape file...");
		imported = pLandscape->readLandscape(0,habmapname," ");
	}
	if (imported != 0) {
		msg = msgReadError + "landscape file: " + Int2Str(imported);
		MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
		return;
	}
	else nloaded++;
}
else
	return;

if (ppLand.rasterType == 1 && ppLand.nHab > 1) { // % cover and multiple layers
	// read additional habitat layers
	MemoLine("Reading additional landscape files...");
	for (int i = 1; i < ppLand.nHab; i++) {
		MemoLine(("Loading map for habitat nr. " + Int2Str(i+1) + ". Please wait...").c_str());
		string fname = hfnames[i-1];
//		string hfsize = Int2Str(hfnames.size());
//		MemoLine(("i=" + Int2Str(i) + " hfsize=" + hfsize + " fname=" + fname).c_str());
		imported = pLandscape->readLandscape(i,fname," ");
		if (imported != 0) {
			msg = msgReadError + "habitat " + Int2Str(i+1) + " file: " + Int2Str(imported);
			MessageDlg(msgReadError.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
			i = ppLand.nHab+1; break; // to break out of habitats loop
		}
		else nloaded++;
	}
	// check that no cell exceeds total cover of 100%
	int nCells = pLandscape->checkTotalCover();
	if (nCells > 0) { // error condition
		string msgCoverError = Int2Str(nCells) + " cell";
		if (nCells == 1) msgCoverError += " has";
		else  msgCoverError += "s have";
		msgCoverError += " total cover exceeding 100%";
		MessageDlg(msgCoverError.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
		imported = -1;
	}
}

if (imported == 0) { // landscape data loaded successfully

	ppLand = pLandscape->getLandParams();

	// display the no. of habitat types
	LabelNhab->Visible = true;
	edtNhabitat->Text = ppLand.nHab;
	edtNhabitat->Visible = true;

	frmDynLand->refresh();
	if (RGhab->ItemIndex == 1) { // habitat % cover
		// dynamic landscape is not available - colours may be changed
		BtnDynamic->Enabled = false; BtnDynamic->Visible = false;
		displayColourTable(ppLand.nHab);
//		edtNhabitat->Text = ppLand.nHab;
//		edtNhabitat->Visible = true;
		BtnChangeColours->Enabled = true; BtnChangeColours->Visible = true;
	}
	else {
		// dynamic landscape is available - colours may not (yet) be changed
		BtnDynamic->Enabled = true;  BtnDynamic->Visible = true;
		BtnChangeColours->Enabled = false; BtnChangeColours->Visible = false;
	}
//	if (RGhab->ItemIndex == 2) { // habitat quality
//		BtnChangeColours->Enabled = false; BtnChangeColours->Visible = false;
//	}
//	else {
//		// display the default habitat colours
//		// or amended colours if already stored
//		displayColourTable(ppLand.nHab);
//		BtnChangeColours->Enabled = true; BtnChangeColours->Visible = true;
//	}
	LabelSpResol->Enabled = true;
	edtSpResol->Enabled = true;
	BtnImportSp->Enabled = true;

//string msgRead2 = "After reading landscape nHab = " + Int2Str(ppLand.nHab);
//MessageDlg(msgRead2.c_str(),mtWarning, TMsgDlgButtons() << mbOK,0);

//	if (!ppLand.patchModel) {
//		frmMain->EnvGradient->Enabled = true;
//	}
	MemoLine("...Landscape loaded");

	landscapeLoaded = true; newLandscape = true;
	forceInit = true;
	costsSet = false; // force costs to be respecified (for SMS only)
	frmVisualCost->Close();

	BtnOK->Enabled = true;
	BtnCancel->Enabled = false;

}
#if RSDEBUG
DebugGUI(("TfrmLand::BtnImportClick(): ppLand.nHab=" + Int2Str(ppLand.nHab)
//	+ " h=" + Int2Str(h) + " p=" + Int2Str(p)
).c_str());
#endif

}

//---------------------------------------------------------------------------
void __fastcall TfrmLand::BtnDynamicClick(TObject *Sender)
{
//frmDynLand->refresh();
while (frmDynLand->ShowModal() == 1) {} // show form until OK or cancelled
}

//------------------------------------------------------------------------------
void __fastcall TfrmLand::BtnImportSpClick(TObject *Sender)
{
int loaded;
rasterdata spdist;
landParams ppLand = pLandscape->getLandParams();
landOrigin origin = pLandscape->getOrigin();
//simParams sim = paramsSim->getSim();
UnicodeString msg;
UnicodeString txtspdist = "species distribution ";

ppLand.spResol = StrToInt(frmLand->edtSpResol->Text);

if (ppLand.spResol < ppLand.resol || ppLand.spResol%ppLand.resol != 0) {
	msg = "The resolution of the "
		+ txtspdist + "must be\n an integer multiple of the landscape resolution";
	MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
	return;
}

// species distribution resolution is OK - get file name
frmLand->OpenTextFileDialog1->InitialDir = paramsSim->getDir(0).c_str();
frmLand->OpenTextFileDialog1->Title = "Select the species distribution file";
if (frmLand->OpenTextFileDialog1->Execute()){
	distnmapname = AnsiString(frmLand->OpenTextFileDialog1->FileName).c_str();
}
else {
	msg = "Selection of " + txtspdist + "raster file cancelled";
	MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
	return;
}

// check that species distribution raster headers are ok
spdist = CheckRasterFile(distnmapname);
if (!spdist.ok) {
	msg = "Format error in header lines of the " + txtspdist + "file";
	MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
	return;
}

// check that headers are compatible with landscape
if (spdist.cellsize != ppLand.spResol) {
	msg = "The raster resolution " + IntToStr(spdist.cellsize)
		+ "m differs from the specified resolution for the " + txtspdist
		+ IntToStr(ppLand.spResol) +"m";
	MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
	return;
}
if ((int)spdist.xllcorner != (int)origin.minEast
||  (int)spdist.yllcorner != (int)origin.minNorth) {
	msg = "Origin co-ordinates of the " + txtspdist + "file\n"
		"do not match those of the landscape file";
	MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
	return;
}
if (spdist.ncols > ppLand.maxX+1 || spdist.nrows > ppLand.maxY+1) {
	msg = "The dimensions of the " + txtspdist + "exceed those of the landscape";
	MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
	return;
}

MemoLine("Reading species distribution file...");
// WILL NEED TO READ DISTRIBUTION FOR CORRECT SPECIES ...
// ... CURRENTLY IT IS THE ONLY ONE
if (pLandscape->distnCount() > 0) pLandscape->deleteDistribution(pSpecies);
loaded = pLandscape->newDistribution(pSpecies,distnmapname);
if (loaded == 0) {
//	sim.initDistLoaded = true;
	ppLand.spDist = true;
	MemoLine("Species distribution loaded");
}
else {
//	sim.initDistLoaded = false;
	ppLand.spDist = false;
	string msg = "Error in distribution file: ";
	if (loaded == 22) { // invalid cell in the distribution file
		msg += "values must be 0 or 1";
		MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbRetry,0);
	}
	else {
		msg += "return code = " + Int2Str(loaded);
		MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbRetry,0);
	}
}

//paramsSim->setSim(sim);
pLandscape->setLandParams(ppLand,false);
//string msg = "ppLand.spResol = " + Int2Str(ppLand.spResol);
//MessageDlg(msg.c_str(),mtWarning, TMsgDlgButtons() << mbOK,0);

}

//------------------------------------------------------------------------------

// Display the default habitat colours or amended colours if already stored
void __fastcall TfrmLand::displayColourTable(int nhab) {
int ncolours = (int)pLandscape->colourCount();
//MemoLine(("No. of stored colours is " + Int2Str(ncolours)).c_str());
//MemoLine(("nhab = " + Int2Str(nhab)).c_str());
StringGridHab->Visible = true;
StringGridHab->Enabled = true;
rgb colour;
StringGridHab->RowCount = nhab + 1;
for (int i = 0; i < nhab; i++) {
	if (RGhab->ItemIndex == 0) { // habitat codes
		StringGridHab->Cells[0][i+1] = pLandscape->getHabCode(i);
	}
	else { // habitat % cover - layers take consecutive integers
		StringGridHab->Cells[0][i+1] = IntToStr(i+1);
  }
	if (i < ncolours) { // colour already exists
		colour = pLandscape->getColour(i);
		StringGridHab->Cells[1][i+1] = colour.r;
		StringGridHab->Cells[2][i+1] = colour.g;
		StringGridHab->Cells[3][i+1] = colour.b;
	}
	else {
		if (i < 21) {
			StringGridHab->Cells[1][i+1] = colour.r = RED[i];
			StringGridHab->Cells[2][i+1] = colour.g = GREEN[i];
			StringGridHab->Cells[3][i+1] = colour.b = BLUE[i];
		}
		else {
			StringGridHab->Cells[1][i+1] = colour.r = 0;
			StringGridHab->Cells[2][i+1] = colour.g = 0;
			StringGridHab->Cells[3][i+1] = colour.b = 0;
    }
		pLandscape->addColour(colour);
	}
}
}

//---------------------------------------------------------------------------

void __fastcall TfrmLand::BtnChangeColoursClick(TObject *Sender)
{
int col;
rgb colours;
bool coloursOK = true;

for (int i = 0; i < StringGridHab->RowCount-1; i++) {
	for (int j = 1; j < 4; j++) {
		col = -999;
		col = StringGridHab->Cells[j][i+1].ToInt();
		if (col < 0 || col > 255) coloursOK = false;
	}
}
if (coloursOK) {
	for (int i = 0; i < StringGridHab->RowCount-1; i++) {
		colours.r = StringGridHab->Cells[1][i+1].ToInt();
		colours.g = StringGridHab->Cells[2][i+1].ToInt();
		colours.b = StringGridHab->Cells[3][i+1].ToInt();
		pLandscape->changeColour(i,colours);
	}
	MessageDlg("Colour values updated",
		mtWarning, TMsgDlgButtons() << mbOK,0);
}
else {
	MessageDlg("Colour values must be from 0 to 255 inclusive",
				mtError, TMsgDlgButtons() << mbOK,0);
}
}
//---------------------------------------------------------------------------

void __fastcall TfrmLand::BtnOKClick(TObject *Sender)
{
landParams ppLand = pLandscape->getLandParams();
BtnDynamic->Enabled = false;

switch (RGhab->ItemIndex) {
case 0: // habitat codes
	if (pLandscape->habitatsIndexed()) {
		Close(); // return control to main form
	}
	else {
		pLandscape->updateHabitatIndices();
		displayColourTable(ppLand.nHab);
		edtNhabitat->Text = ppLand.nHab;
		edtNhabitat->Visible = true;
		BtnChangeColours->Enabled = true; BtnChangeColours->Visible = true;
	}
	break;
case 1: // habitat % cover
	Close(); // return control to main form
	break;
case 2: // habitat quality
	Close(); // return control to main form
	break;
default:
	;
}

}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


