//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "FormArtificialLand.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"

ofstream outLandFile;
ofstream Fmaps;
vector<land> artLand; 		// used in generating fractal landscapes
TfrmGenerateLand *frmGenerateLand;

int **patch_land = NULL;	// used in generating patches
int cellCount;

//---------------------------------------------------------------------------
__fastcall TfrmGenerateLand::TfrmGenerateLand(TComponent* Owner)
	: TForm(Owner)
{
Top = 60; Left = 100;
StringGridLand->Cells[0][1] = "p";
StringGridLand->Cells[1][1] = "0.2";
SetMemo();
}

//---------------------------------------------------------------------------
// Common edit checks for both generation options
bool __fastcall TfrmGenerateLand::CheckDimensions(void) {
string msg,msgtxt1,msgtxt2;
int xDim = StrToInt(edtXdim->Text);
int yDim = StrToInt(edtYdim->Text);
float maxPct = StrToFloat(edtMaxProp->Text);
float minPct = StrToFloat(edtMinProp->Text);
int fractal =  RGFract->ItemIndex;
int continuous = RGLandType->ItemIndex;

msgtxt1 = "Invalid dimensions: ";
if (fractal) { // conditions for fractal landcapes
	if (xDim < 3 || yDim < 3)
	{
		msg = msgtxt1 + "X and Y must be 3 or more";
		MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
		return false;
	}
	msgtxt2 = "please ensure ";
	if (xDim > yDim)
	{
		msg = msgtxt1 + msgtxt2 + "Y >= X";
		MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
		return false;
	}
	if ((xDim > 2 && power2check(xDim-1) != 1)
	||  (yDim > 2 && power2check(yDim-1) != 1)) {
		msg = msgtxt1 + msgtxt2 + "that X & Y are powers of 2 + 1";
		MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
		return false;
	}
}
else { // conditions for random landscapes
	if (xDim < 1 || yDim < 1) {
		msg = msgtxt1  + "X and Y must be 1 or more";
		MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
		return false;
	}
}

if (continuous) {
	if (minPct <= 0.0) {
		MessageDlg("The minimum % of habitat must be greater than 0.0",
				mtError, TMsgDlgButtons() << mbOK,0);
		return false;
	}
	if (maxPct > 100.0) {
		MessageDlg("The maximum % of habitat may not be greater than 100.0",
				mtError, TMsgDlgButtons() << mbOK,0);
		return false;
	}
	if (maxPct < minPct) {
		MessageDlg("The maximum % of habitat may not be greater than the minimum",
				mtError, TMsgDlgButtons() << mbOK,0);
		return false;
	}
}
return true;
}

//---------------------------------------------------------------------------
void __fastcall TfrmGenerateLand::CBpatchClick(TObject *Sender)
{
if (CBpatch->Checked) {
	edtMaxCells->Enabled = true; edtMaxCells->Visible = true;
}
else {
	edtMaxCells->Enabled = false; edtMaxCells->Visible = false;
}
}
//---------------------------------------------------------------------------
void __fastcall TfrmGenerateLand::RGLandSeriesClick(TObject *Sender)
{
switch (RGLandSeries->ItemIndex) {

case 0: // new landscape for each replicate
	StatusBar1->Panels->Items[0]->Text = ""; StatusBar1->Refresh();
	edtSeriesNr->Enabled = false; edtSeriesNr->Visible = false;
	edtNLand->Enabled = false; edtNLand->Visible = false;
	CBpatch->Enabled = false; CBpatch->Visible = false;
	BtnCreateSeries->Enabled = false; BtnCreateSeries->Visible = false;
	edtMaxCells->Enabled = false; edtMaxCells->Visible = false;
	CBmatrix->Enabled = false; CBmatrix->Visible = false;
	GridMatrixLabel->Visible = false;
	StringGridMatrix->Enabled = false; StringGridMatrix->Visible = false;
	SetMemo();
	StringGridLand->ColCount = 2;
	if (RGFract->ItemIndex == 0) { // random
		StringGridLand->RowCount = 2;
		StringGridLand->Cells[0][1] = "p";
		StringGridLand->Cols[1]->Clear();
		StringGridLand->Cells[1][1] = "0.2";
	}
	else { // fractal
		StringGridLand->RowCount = 3;
		StringGridLand->Cols[1]->Clear();
		StringGridLand->Cells[1][1] = "0.2";
		StringGridLand->Cells[1][2] = "0.1";
		StringGridLand->Cells[0][1] = "p";
		StringGridLand->Cells[0][2] = "H";
	}
	BtnOK->Enabled = true;
	break;

case 1: // generate and save a series
	edtSeriesNr->Enabled = true; edtSeriesNr->Visible = true;
	edtNLand->Enabled = true; edtNLand->Visible = true;
	CBpatch->Enabled = true; CBpatch->Visible = true;
	BtnCreateSeries->Enabled = true; BtnCreateSeries->Visible = true;
	if (CBpatch->Checked) {
		edtMaxCells->Enabled = true; edtMaxCells->Visible = true;
	}
	else {
		edtMaxCells->Enabled = false; edtMaxCells->Visible = false;
	}
	if (RGLandType->ItemIndex == 0) { // discrete
		CBmatrix->Enabled = true; CBmatrix->Visible = true;
		if (CBmatrix->Checked) {
			GridMatrixLabel->Visible = true;
			StringGridMatrix->Enabled = true; StringGridMatrix->Visible = true;
		}
	}
	else { // continuous
		CBmatrix->Checked = false;
		CBmatrix->Enabled = false; CBmatrix->Visible = false;
		GridMatrixLabel->Visible = false;
		StringGridMatrix->Enabled = false; StringGridMatrix->Visible = false;
	}
	SetMemo();
	if (RGFract->ItemIndex == 0) { // random
		StringGridLand->ColCount = 11;
		StringGridLand->RowCount = 2;
		StringGridLand->Cells[0][1] = "p";
		for (int i = 1; i <= StringGridLand->ColCount; i++){
		 StringGridLand->Cols[i]->Clear();
		 StringGridLand->Cells[i][0] = IntToStr(i);
		}
	}
	else { // fractal
		StringGridLand->ColCount = 11;
		StringGridLand->RowCount = 3;
		StringGridLand->Cells[0][1] = "p";
		StringGridLand->Cells[0][2] = "H";
		for (int i = 1; i <= StringGridLand->ColCount; i++){
		 StringGridLand->Cols[i]->Clear();
		 StringGridLand->Cells[i][0] = IntToStr(i);
		}
	}
	StringGridMatrix->ColCount = 11;
	StringGridMatrix->RowCount = 2;
	StringGridMatrix->Cells[0][1] = "m";
	for (int i = 1; i <= StringGridMatrix->ColCount; i++){
		StringGridMatrix->Cols[i]->Clear();
		StringGridMatrix->Cells[i][0] = IntToStr(i);
	}
	BtnOK->Enabled = false;
	break;

}
if (RGLandType->ItemIndex == 0) { // discrete
	edtMaxProp->Enabled = false; edtMaxProp->Visible = false;
	edtMinProp->Enabled = false; edtMinProp->Visible = false;
}
else { // continuous
	edtMaxProp->Enabled = true; edtMaxProp->Visible = true;
	edtMinProp->Enabled = true; edtMinProp->Visible = true;
}
}
//---------------------------------------------------------------------------
void __fastcall TfrmGenerateLand::RGFractClick(TObject *Sender)
{
if (RGFract->ItemIndex == 0) { // random landcape
	if (RGLandSeries->ItemIndex == 0) { // new landscape for each replicate
		StringGridLand->ColCount = 2;
		StringGridLand->RowCount = 2;
		StringGridLand->Cells[0][1] = "p";
		StringGridLand->Cols[1]->Clear();
		StringGridLand->Cells[1][1] = "0.2";
	}
	else { // generate and save a series
		StringGridLand->ColCount = 11;
		StringGridLand->RowCount = 2;
		StringGridLand->Cells[0][1] = "p";
		for (int i = 1; i <= StringGridLand->ColCount; i++){
			StringGridLand->Cols[i]->Clear();
			StringGridLand->Cells[i][0] = IntToStr(i);
		}
	}
	SetMemo();
}
else { // fractal landscape
	if (RGLandSeries->ItemIndex == 0) { // new landscape for each replicate
		StringGridLand->ColCount = 2;
		StringGridLand->RowCount = 3;
		StringGridLand->Cols[1]->Clear();
		StringGridLand->Cells[1][1] = "0.2";
		StringGridLand->Cells[1][2] = "0.1";
		StringGridLand->Cells[0][1] = "p";
		StringGridLand->Cells[0][2] = "H";
	}
	else { // generate and save a series
	  StringGridLand->ColCount = 11;
	  StringGridLand->RowCount = 3;
	  StringGridLand->Cells[0][1] = "p";
	  StringGridLand->Cells[0][2] = "H";
		for (int i = 1; i <= StringGridLand->ColCount; i++){
	   StringGridLand->Cols[i]->Clear();
	   StringGridLand->Cells[i][0] = IntToStr(i);
	  }
	}
	SetMemo();
}
}
//---------------------------------------------------------------------------
void __fastcall TfrmGenerateLand::RGLandTypeClick(TObject *Sender)
{
if (RGLandType->ItemIndex == 0) { // discrete
	edtMaxProp->Enabled = false; edtMaxProp->Visible = false;
	edtMinProp->Enabled = false; edtMinProp->Visible = false;
	if (RGLandSeries->ItemIndex == 1) {
		CBmatrix->Enabled = true; CBmatrix->Visible = true;
		if (CBmatrix->Checked) {
			GridMatrixLabel->Visible = true;
			StringGridMatrix->Enabled = true; StringGridMatrix->Visible = true;
		}
	}
}
else { // continuous
	edtMaxProp->Enabled = true; edtMaxProp->Visible = true;
	edtMinProp->Enabled = true; edtMinProp->Visible = true;
	CBmatrix->Checked = false;
	CBmatrix->Enabled = false; CBmatrix->Visible = false;
	GridMatrixLabel->Visible = false;
	StringGridMatrix->Enabled = false; StringGridMatrix->Visible = false;
}
}

//---------------------------------------------------------------------------
void __fastcall TfrmGenerateLand::CBmatrixClick(TObject *Sender)
{
if (CBmatrix->Checked) {
	GridMatrixLabel->Visible = true;
	StringGridMatrix->Visible = true; StringGridMatrix->Enabled = true;
}
else {
	GridMatrixLabel->Visible = false;
	StringGridMatrix->Visible = false; StringGridMatrix->Enabled = false;
}
SetMemo();
}

//---------------------------------------------------------------------------
void __fastcall TfrmGenerateLand::SetMemo(void) {
Memo2->Visible = true;
Memo2->Clear();
Memo2->Lines->Add("NOTES:");
Memo2->Lines->Add("p = proportion of suitable habitat in the landscape");
if (RGFract->ItemIndex == 1) { // fractal landcape
	Memo2->Lines->Add("H = Hurst exponent");
}
if (RGLandSeries->ItemIndex == 1 && CBmatrix->Checked) {
	Memo2->Lines->Add("m = proportion of second randomly distributed matrix habitat (max. 0.5)");
}
}

//---------------------------------------------------------------------------
void __fastcall TfrmGenerateLand::BtnOKClick(TObject *Sender)
{
string msg,msgtxt1,msgtxt2;
landParams ppLand = pLandscape->getLandParams();
genLandParams ppGenLand = pLandscape->getGenLandParams();

//msg = "ppGenLand.maxY " + Int2Str(ppLand.maxY)
//	+ " ppGenLand.maxX " + Int2Str(ppLand.maxX);
//MessageDlg(msg.c_str(),mtWarning, TMsgDlgButtons() << mbOK,0);

if (!CheckDimensions()) return;

ppLand.generated = true;
ppLand.nHab = 2;

ppLand.dimX = StrToInt(edtXdim->Text);
ppLand.dimY = StrToInt(edtYdim->Text);
ppGenLand.maxPct = StrToFloat(edtMaxProp->Text);
ppGenLand.minPct = StrToFloat(edtMinProp->Text);
ppLand.resol = StrToInt(edtResolution->Text);
ppGenLand.fractal =  RGFract->ItemIndex;
ppGenLand.continuous = RGLandType->ItemIndex;
ppGenLand.propSuit = (StringGridLand->Cells[1][1]).ToDouble();

msgtxt1 = " must lie in the range ";
if (ppGenLand.propSuit <= 0.0 || ppGenLand.propSuit > 1.0) {
	msg = "p" + msgtxt1 + ">0 to 1";
	MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
	return;
}

// no. of habitat cells must be > 0 for random landscape
float habCells = ppGenLand.propSuit * (float)(ppLand.dimX) * (float)(ppLand.dimY);
msgtxt2 = "No. of habitat cells (i.e. p * no. cells) must be ";
#if RSDEBUG
//DebugGUI(("TfrmGenerateLand::ButtonOKClick(): dimX=" + Int2Str(ppLand.dimX)
//	+ " dimY=" + Float2Str(ppLand.dimY)
//	+ " propSuit=" + Float2Str(ppGenLand.propSuit)
//	+ " habCells=" + Int2Str(habCells)).c_str());
#endif
if (!ppGenLand.fractal && habCells < 1.0) {
	msg = msgtxt2 + "1 or more";
	MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
	return;
}

if (ppGenLand.fractal) {
	// no. of habitat cells must be > 1 for fractal algorithm
	if (habCells < 2.0) {
		msg = msgtxt2 + "2 or more";
		MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
		return;
	}
	ppGenLand.hurst = (StringGridLand->Cells[1][2]).ToDouble();
	if (ppGenLand.hurst <= 0.0 || ppGenLand.hurst >= 1.0) {
		msg = "H" + msgtxt1 + ">0 to <1";
		MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
		return;
	}
}

ppLand.minX = ppLand.minY = 0;
ppLand.maxX = ppLand.dimX - 1; ppLand.maxY = ppLand.dimY - 1;
pLandscape->setLandParams(ppLand,false);
pLandscape->setGenLandParams(ppGenLand);

artLandStatus = 9;
frmGenerateLand->Close();

}

//---------------------------------------------------------------------------
void __fastcall TfrmGenerateLand::BtnCreateSeriesClick(TObject *Sender)
{
string LandFileName;
string msg;
string msg00 = "Please ensure that valid ";
string msg01 = " values have been entered";
string msg10 = " values must be between ";
string msg11 = "0 and 1 (exclusive)";
int SeriesNr;
bool p_error,cells_error,H_error,m_error;
double ppp,habCells,p_sum,H_sum,m_sum;

landParams ppLand = pLandscape->getLandParams();
genLandParams ppGenLand = pLandscape->getGenLandParams();

SeriesNr = StrToInt(edtSeriesNr->Text);
ppLand.landNum = StrToInt(edtNLand->Text);
if (ppLand.landNum < 1) {
	MessageDlg("No. of replicates must be 1 or more"
		,mtError, TMsgDlgButtons() << mbOK,0);
	return;
}

if (!CheckDimensions()) return;

ppLand.generated = true;
ppLand.dimX = StrToInt(edtXdim->Text);
ppLand.dimY = StrToInt(edtYdim->Text);
ppGenLand.maxPct = StrToFloat(edtMaxProp->Text);
ppGenLand.minPct = StrToFloat(edtMinProp->Text);
ppLand.resol = StrToInt(edtResolution->Text);
ppGenLand.fractal =  RGFract->ItemIndex;
ppGenLand.continuous = RGLandType->ItemIndex;

if (CBpatch->Checked) {
	ppLand.patchModel = true;
	ppGenLand.maxCells = StrToInt(edtMaxCells->Text);
	if (ppGenLand.maxCells < 1) {
		MessageDlg("Max. cells per patch must be 1 or more",
			mtError, TMsgDlgButtons() << mbOK,0);
		return;
	}
}
else ppLand.patchModel = false;

p_error = cells_error = H_error = false;
p_sum = H_sum = 0.0;
for (int i = 1; i < 11; i++) {
	if (!(StringGridLand->Cells[i][1]).IsEmpty()) {
		p_sum += ppp = (StringGridLand->Cells[i][1]).ToDouble();
		if (ppp <= 0.0 || ppp >= 1.0) p_error = true;
		habCells = ppp * (double)(ppLand.dimX) * (double)(ppLand.dimY);
		// no. of suitable habitat cells must be > 0 for random landscape
		// or >1 for fractal landscape
		if (ppGenLand.fractal && habCells < 2.0 ) cells_error = true;
		if (!ppGenLand.fractal && habCells < 1.0 ) cells_error = true;
	}
#if RSDEBUG
//DebugGUI(("TfrmGenerateLand::BtnCreateSeriesClick(): dimX=" + Int2Str(ppLand.dimX)
//	+ " dimY=" + Float2Str(ppLand.dimY)
//	+ " ppp=" + Float2Str(ppp)
//	+ " habCells=" + Int2Str(habCells)).c_str());
#endif
	if (ppGenLand.fractal) {
		if (!(StringGridLand->Cells[i][2]).IsEmpty()) {
			H_sum += ppp = (StringGridLand->Cells[i][2]).ToDouble();
			if (ppp <= 0.0 || ppp >= 1.0) H_error = true;
		}
	}
}
if (p_sum <= 0.0) {
	msg = msg00 + "p" + msg01;
	MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
	return;
}
if (p_error) {
	msg = "p" + msg10 + msg11;
	MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
	return;
}
if (cells_error) {
	msg = "No. of habitat cells (i.e. p * no. cells) must be ";
	if (ppGenLand.fractal) msg += "2 or more"; else msg += "1 or more";
	MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
	return;
}
if (ppGenLand.fractal){
	if (H_sum <= 0.0) {
		msg = msg00 + "H" + msg01;
		MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
		return;
	}
	if (H_error) {
		msg = "H" + msg10 + msg11;
		MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
		return;
	}
}

if (CBmatrix->Checked) {
	// check table of proportions of second matrix habitat
	m_error = false;
	m_sum = 0.0;
	for (int i = 1; i < 11; i++) {
		if (!(StringGridMatrix->Cells[i][1]).IsEmpty()) {
			m_sum += ppp = (StringGridMatrix->Cells[i][1]).ToDouble();
			if (ppp <= 0.0 || ppp > 0.5) m_error = true;
		}
	}
	if (m_sum <= 0.0) {
		msg = msg00 + "m" + msg01;
		MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
		return;
	}
	if (m_error) {
		msg = "m" + msg10 + "0 and 0.5";
		MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
		return;
	}
}

// Check that LandFile.txt does not already exist in Inputs folder
ifstream dummyLandFile;
//LandFileName = paramsSim->getDir(1) + "LandFile" + artInt2Str(SeriesNr) + ".txt";
LandFileName = paramsSim->getDir(1) + "LandFile" + Int2Str(SeriesNr) + ".txt";
dummyLandFile.open(LandFileName.c_str());
if (dummyLandFile.is_open()) { // file exists
//	msg = "LandFile" + artInt2Str(SeriesNr) + ".txt" + " already exists in Inputs folder.\n"
	msg = "LandFile" + Int2Str(SeriesNr) + ".txt" + " already exists in Inputs folder.\n"
		"Please either move, rename or delete it,\n or use a different Series number.";
	MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
	dummyLandFile.close(); dummyLandFile.clear();
	return;
}
else { // open LandFile for output
	outLandFile.open(LandFileName.c_str());
	outLandFile << "LandNum\tNhabitats\tLandscapeFile\tPatchFile\tCostMapFile\tDynLandFile\tSpDistFile"
		<< endl;
}

ppLand.minX = ppLand.minY = 0;
ppLand.maxX = ppLand.dimX - 1; ppLand.maxY = ppLand.dimY - 1;
pLandscape->setLandParams(ppLand,false);
pLandscape->setGenLandParams(ppGenLand);

// Generate the landscapes
artLandStatus = 0;
if (!LandSeries(SeriesNr,ppGenLand.fractal,CBmatrix->Checked)) {
	if (ppGenLand.fractal) artLandStatus = 1; else artLandStatus = 2;
}
/*
if (ppGenLand.fractal) {
	if (!LandSeries(SeriesNr,true)) {
		artLandStatus = 1;
	}
}
else {
	if (!LandSeries(SeriesNr,false,CBmatrix->Checked)) {
		artLandStatus = 2;
	}
}
*/
outLandFile.close(); outLandFile.clear();

frmGenerateLand->Close();

}

//---------------------------------------------------------------------------
void __fastcall TfrmGenerateLand::BtnCancelClick(TObject *Sender)
{
frmGenerateLand->Close();
}

//---------------------------------------------------------------------------
// Generate a series of fractal landscapes
int TfrmGenerateLand::LandSeries(int SeriesNr,bool fractal,bool randmatrix)
{
#if RSDEBUG
//DebugGUI(("LandSeries(): SeriesNr=" + Int2Str(SeriesNr)
//	+ " fractal=" + Int2Str(fractal) + " randmatrix=" + Int2Str(randmatrix)
//	).c_str());
#endif
int x,y,jmax,kmax,ncells,nmatrixcells;
double p,H,pmatrix;
bool firstcolumn;
float **art_land;
float **random_matrix;
// string FragMapName;
vector<land>::iterator iter;
string MapName,patchMapName,name,contdisc;
landParams ppLand = pLandscape->getLandParams();
genLandParams ppGenLand = pLandscape->getGenLandParams();
string Inputs = paramsSim->getDir(1);

StatusBar1->Panels->Items[0]->Text =
	"Generating maps series. Please wait...";
StatusBar1->Refresh();

SeriesNr = StrToInt(edtSeriesNr->Text);
ppLand.landNum = StrToInt(edtNLand->Text);

// create a 2D array
art_land = new float *[ppLand.dimX];
for (int xx = 0; xx < ppLand.dimX; xx++) {
	art_land[xx] = new float[ppLand.dimY];
}
int nhabitats = 2;
if (randmatrix) { // second matrix habitat
	random_matrix = new float *[ppLand.dimX];
	for (int xx = 0; xx < ppLand.dimX; xx++) {
		random_matrix[xx] = new float[ppLand.dimY];
	}
	nhabitats = 3;
}
if (ppLand.patchModel) CreatePatchLand();

int landnum = 1;
for (int i = 1; i <= 10; i++) { // each propn. of suitable habitat
	if (!(StringGridLand->Cells[i][1]).IsEmpty()) {
		p = (StringGridLand->Cells[i][1]).ToDouble(); // propn. suitable habitat
		if (fractal) jmax = 10; else jmax = 1; // no. of possible Hurst exponents
		for (int j = 1; j <= jmax; j++) { // each Hurst exponent
			if (!fractal || !(StringGridLand->Cells[j][2]).IsEmpty()) {
				if (fractal) H = (StringGridLand->Cells[j][2]).ToDouble(); // Hurst exponent
				else H = 0.0; // random landscape
				if (randmatrix) kmax = 10; else kmax = 1; // no. of possible random matrix propns.
				for (int k = 1; k <= kmax; k++) { // each matrix propn.
					if (!randmatrix || !(StringGridMatrix->Cells[k][1].IsEmpty())) {
						if (randmatrix) pmatrix = (StringGridMatrix->Cells[k][1]).ToDouble(); // matrix propn.
						else pmatrix = 0.0;

						for (int rr = 0; rr < ppLand.landNum; rr++) { // each replicate
							// reset the 2D array(s)
							for (int xx = 0; xx < ppLand.dimX; xx++) {
								for (int yy = 0; yy < ppLand.dimY; yy++) {
									art_land[xx][yy] = 0.0;
									if (randmatrix) random_matrix[xx][yy] = 0.0;
								}
							}
							// calculate no. of suitable habitat cells
							ncells = (int)((ppLand.dimX) * (ppLand.dimY) * (p + 0.000001));
							if (fractal) {
								// NOTE: x & y are swapped in the function fractal_landscape()
								artLand = fractal_landscape(ppLand.dimY,ppLand.dimX,H,(1.0-p),
									ppGenLand.maxPct,ppGenLand.minPct);
								iter = artLand.begin();
								while (iter != artLand.end()){
									x = iter->y_coord; y = iter->x_coord;
									art_land[x][y] = iter->value;
#if RSDEBUG
//DebugGUI(("LandSeries(): x=" + Int2Str(x)
//	+ " y=" + Int2Str(y) 
//	+ " value=" + Float2Str(art_land[x][y])).c_str());
#endif
									iter++;
								}
								if (!ppGenLand.continuous) {
									if (iter->avail == 0) art_land[x][y] = 0.0;
								}
								artLand.clear();
							} // end of fractal
							else { // random
#if RSDEBUG
//DebugGUI(("LandSeries(): rr=" + Int2Str(rr)
//	+ " dimX=" + Int2Str(ppLand.dimX) + " dimY=" + Int2Str(ppLand.dimY)
//	+ " p=" + Float2Str(p) + " ncells=" + Int2Str(ncells)).c_str());
#endif
								int i = 0;
								do {
									do {
										x = pRandom->IRandom(0,ppLand.dimX-1); y = pRandom->IRandom(0,ppLand.dimY-1);
									} while (art_land[x][y] > 0.0);
									if (ppGenLand.continuous) {
										art_land[x][y] =
											ppGenLand.minPct + pRandom->Random() * (ppGenLand.maxPct - ppGenLand.minPct);
									}
									else { // discrete
										art_land[x][y] = 1.0;
									}
#if RSDEBUG
//DebugGUI(("LandSeries(): i=" + Int2Str(i)
//	+ " x=" + Int2Str(x) + " y=" + Int2Str(y)
//	+ " art_land=" + Float2Str(art_land[x][y])).c_str());
#endif
									i++;
								} while(i < ncells);
							} // end of random

							if (randmatrix) {
								// change specified propn. of matrix cells to habitat code 3
								nmatrixcells = (int)((ppLand.dimX * ppLand.dimY - ncells) * (pmatrix + 0.000001));
#if RSDEBUG
//DebugGUI(("LandSeries(): rr=" + Int2Str(rr)
//	+ " dimX=" + Int2Str(ppLand.dimX) + " dimY=" + Int2Str(ppLand.dimY)
//	+ " pmatrix=" + Float2Str(pmatrix) + " nmatrixcells=" + Int2Str(nmatrixcells)).c_str());
#endif
								int i = 0;
								do {
									do {
										x = pRandom->IRandom(0,ppLand.dimX-1); y = pRandom->IRandom(0,ppLand.dimY-1);
									} while (art_land[x][y] > 0.0);
									art_land[x][y] = 999999.0;
#if RSDEBUG
//DebugGUI(("LandSeries(): i=" + Int2Str(i)
//	+ " x=" + Int2Str(x) + " y=" + Int2Str(y)
//	+ " art_land=" + Float2Str(art_land[x][y])).c_str());
#endif
									i++;
								} while(i < nmatrixcells);
							}

							if (ppGenLand.continuous) contdisc = "cont"; else contdisc = "disc";
							MapName = "Series" + Int2Str(SeriesNr) + contdisc;
							if (!fractal) MapName += "Random";
							MapName += "_X" + Int2Str(ppLand.dimX) + "Y" + Int2Str(ppLand.dimY) + "_p" + Float2Str(p);
							if (fractal) MapName += "H" + Float2Str(H);
							if (randmatrix) MapName += "_matrix" + Float2Str(pmatrix);
							MapName += "_nr" + Int2Str(rr) + ".txt";
							if (ppLand.patchModel) {
								patchMapName = "Series" + Int2Str(SeriesNr) + contdisc + "Patch";
								patchMapName += "_X" + Int2Str(ppLand.dimX) + "Y" + Int2Str(ppLand.dimY) + "_p" + Float2Str(p);
								if (fractal) patchMapName += "H" + Float2Str(H);
								if (randmatrix) patchMapName += "_matrix" + Float2Str(pmatrix);
								patchMapName += "_nr" + Int2Str(rr) + ".txt";
							}
							name = Inputs + MapName;
							Fmaps.open(name.c_str());
//							streamsize prec = Fmaps.precision();
//						 if (CBFragstats->Checked) {
//							FragMapName = FractMap + "Series" + Int2Str(SeriesNr) + "_Frag/s"
//								+ Int2Str(SeriesNr) + "_discrete_p" + Float2Str(1-p)
//								+ H" + Float2Str(H)+ "_nr"+ Int2Str(jj)+"_Frag.txt";
//							FmapsFrag.open(FragMapName.c_str());
//						 }
							// headers for raster format
							Fmaps << "ncols\t" << ppLand.dimX << endl;
							Fmaps << "nrows\t" << ppLand.dimY << endl;
							Fmaps << "xllcorner\t" << 0.0 << endl;
							Fmaps << "yllcorner\t" << 0.0 << endl;
							Fmaps << "cellsize\t" << ppLand.resol << endl;
							Fmaps << "NODATA_value\t" << -9999 << endl;
							for (y = ppLand.dimY-1; y >= 0; y--) {
								firstcolumn = true;
								for (x = 0; x < ppLand.dimX; x++) {
									if (firstcolumn) firstcolumn = false;
									else Fmaps << "\t";
									if (ppGenLand.continuous) {
										Fmaps << setprecision(5) << art_land[x][y];
									}
									else { // discrete
										// habitat codes are 1 and 2 (and 3) for consistency with batch mode input
										if (art_land[x][y] > 9999.0) { // second matrix class
											Fmaps << "3";
										}
										else {
											if (fractal) {
												if (art_land[x][y] <= 0.0) Fmaps << "1"; // matrix
												else Fmaps << "2"; // suitable habitat
											}
											else { // random
												if (art_land[x][y] < 0.5) Fmaps << "1"; // matrix
												else Fmaps << "2"; // suitable habitat
											}
										}
									}
//									if (CBFragstats->Checked) {
//										if (art_land[x][y] > 0) FmapsFrag << 1 <<"\t";
//										else FmapsFrag << 999 <<"\t";
//									}
								}
								Fmaps << endl;
							}
//							if (CBFragstats->Checked) FmapsFrag << endl;
							Fmaps.close(); Fmaps.clear();

							if (ppLand.patchModel) LandSeriesPatches(patchMapName,art_land);

//							if (CBFragstats->Checked) FmapsFrag.close();
							// write entry to LandFile
							if (ppLand.patchModel)
								outLandFile << landnum++ << "\t" << nhabitats
									<< "\t" << MapName << "\t" << patchMapName;
							else
								outLandFile << landnum++ << "\t" << nhabitats
									<< "\t" << MapName << "\tNULL";
							outLandFile << "\tNULL\tNULL\tNULL" << endl;

						} //end of replicate loop
					}
				} // end of matrix propn. loop
			}
		} // end of Hurst exponent loop
	}
} // end of suitable habitat loop


//RGLandSeries->ItemIndex = 0;

if (ppLand.patchModel) DeletePatchLand();

return 0;
}

//---------------------------------------------------------------------------
// Utility functions for artifical landscape patches

void TfrmGenerateLand::CreatePatchLand(void)
{
landParams ppLand = pLandscape->getLandParams();
#if RSDEBUG
//DebugGUI("CreatePatchLand(): patch_land=" + Int2Str((int)patch_land)
//	+ " dimX=" + Int2Str(ppLand.dimX) + " dimY=" + Int2Str(ppLand.dimY)
//);
#endif

patch_land = new int *[ppLand.dimX];
for (int xx = 0; xx < ppLand.dimX; xx++) {
	patch_land[xx] = new int[ppLand.dimY];
	for (int yy = 0; yy < ppLand.dimY; yy++) {
		patch_land[xx][yy] = 0;
	}
}
}

void TfrmGenerateLand::DeletePatchLand(void)
{
landParams ppLand = pLandscape->getLandParams();
#if RSDEBUG
//DebugGUI("DeletePatchLand(): patch_land=" + Int2Str((int)patch_land)
//	+ " dimX=" + Int2Str(ppLand.dimX) + " dimY=" + Int2Str(ppLand.dimY)
//);
//DebugGUI("DeletePatchLand(): patch_land[0]=" + Int2Str((int)patch_land[0])
//	+ " patch_land[dimX-1]=" + Int2Str((int)patch_land[ppLand.dimX-1])
//);
#endif

if (patch_land != NULL) {
	for (int xx = 0; xx < ppLand.dimX; xx++) {
		delete[] patch_land[xx];
	}
#if RSDEBUG
//	DebugGUI("DeletePatchLand(): patch_land=" + Int2Str((int)patch_land)
//		);
//	for (int xx = 0; xx < ppLand.dimX; xx++) {
//		DebugGUI("DeletePatchLand(): xx=" + Int2Str(xx)
//			+ " patch_land[xx]=" + Int2Str((int)patch_land[xx])
//		);
//	}
#endif
	delete[] patch_land;
	patch_land = NULL;
}
}

//---------------------------------------------------------------------------
// Recursive function to add neighbouring cells to a patch
void TfrmGenerateLand::AddNeighbours(int xx,int yy,int patchnum) {

int xxx,yyy,dx,dy,ddx,ddy,dx0,dx1,dy0,dy1;
landParams ppLand = pLandscape->getLandParams();
genLandParams ppGenLand = pLandscape->getGenLandParams();

// randomise directions searched to prevent quasi-linear patches when
// max patch size is set
if (pRandom->Bernoulli(0.5)) ddx = 1; else ddx = -1;
if (pRandom->Bernoulli(0.5)) ddy = 1; else ddy = -1;
if (ddx > 0) { dx0 = -1; dx1 = 1; } else { dx0 = 1; dx1 = -1; }
if (ddy > 0) { dy0 = -1; dy1 = 1; } else { dy0 = 1; dy1 = -1; }

#if RSDEBUG
//int ID = pRandom->IRandom(10000,99999);
//DebugGUI((Int2Str(0)).c_str());
//DebugGUI(("AddNeighbours(): NEW ID=" + Int2Str(ID)
//	+ " xx=" + Int2Str(xx) + " yy=" + Int2Str(yy)
//	+ " patchnum=" + Int2Str(patchnum)
//	+ " cellCount=" + Int2Str(cellCount)
//	).c_str());
#endif

// mark immediately neighbouring cells for addition to the patch (subject to maximum)
dx = dx0; dy = dy0;
do {
	do {
		if (cellCount >= ppGenLand.maxCells) { // add no more cells
			dx = 2*dx1; dy = 2*dy1;
		}
		else {
			if (dx != 0 || dy != 0) {
				xxx = xx+dx; yyy = yy+dy;
				if (xxx >= 0 && xxx < ppLand.dimX && yyy >= 0 && yyy < ppLand.dimY) {
					if (patch_land[xxx][yyy] == -1) { // suitable and not yet allocated to a patch
						// add to the current patch
//						patch_land[xxx][yyy] = patchnum;
						patch_land[xxx][yyy] = -2;
						cellCount++;
						if (cellCount >= ppGenLand.maxCells) { // add no more cells
							dx = 2*dx1; dy = 2*dy1;
						}
					}
#if RSDEBUG
//DebugGUI(("AddNeighbours(): ID=" + Int2Str(ID)
//	+ " xx=" + Int2Str(xx) + " yy=" + Int2Str(yy)
//	+ " ddx=" + Int2Str(ddx) + " ddy=" + Int2Str(ddy)
//	+ " xxx=" + Int2Str(xxx) + " yyy=" + Int2Str(yyy)
//	+ " patch_land[" + Int2Str(xxx) + "][" + Int2Str(yyy) + "]=" + Int2Str(patch_land[xxx][yyy])
//	+ " patchnum=" + Int2Str(patchnum)
//	+ " cellCount=" + Int2Str(cellCount)
//	).c_str());
#endif
				}
			}
			dy += ddy;
		}
	} while ((ddy > 0 && dy <= dy1) || (ddy < 0 && dy >= dy1));
	dx += ddx; dy = dy0;
} while ((ddx > 0 && dx <= dx1) || (ddx < 0 && dx >= dx1));

// add any of their neighbours
dx = dx0; dy = dy0;
do {
	do {
		if (dx != 0 || dy != 0) {
			xxx = xx+dx; yyy = yy+dy;
			if (xxx >= 0 && xxx < ppLand.dimX && yyy >= 0 && yyy < ppLand.dimY) {
				if (patch_land[xxx][yyy] == -2) { // allocated to current patch
					// add to the current patch
					patch_land[xxx][yyy] = patchnum;
					AddNeighbours(xxx,yyy,patchnum); // ... of the new cell ...
				}
			}
		}
		dy += ddy;
	} while ((ddy > 0 && dy <= dy1) || (ddy < 0 && dy >= dy1));
	dx += ddx; dy = dy0;
} while ((ddx > 0 && dx <= dx1) || (ddx < 0 && dx >= dx1));

}

//---------------------------------------------------------------------------
// Generate patches for an artifical landscape
void TfrmGenerateLand::LandSeriesPatches(const string patchMapName, float **art_land)
{
bool firstcolumn;
int patchnum,x,y;
string name;
landParams ppLand = pLandscape->getLandParams();
genLandParams ppGenLand = pLandscape->getGenLandParams();
string Inputs = paramsSim->getDir(1);

// reset the 2D array - suitable habitat is temporaily identified as -1
for (int xx = 0; xx < ppLand.dimX; xx++) {
	for (int yy = 0; yy < ppLand.dimY; yy++) {
		patch_land[xx][yy] = 0;
		if (ppGenLand.continuous) {
			if (art_land[xx][yy] > 0.0) patch_land[xx][yy] = -1;
		}
		else {
			if (ppGenLand.fractal) {
				// discrete habitat lies in range 1 to 100
				if (art_land[xx][yy] > 0.0 && art_land[xx][yy] < 9999.0)
					patch_land[xx][yy] = -1;
			}
			else { // random
				// discrete habitat codes are still set to 0 and 1
				if (art_land[xx][yy] >= 0.5 && art_land[xx][yy] <= 1.5)
					patch_land[xx][yy] = -1;
			}
		}
	}
}

patchnum = 1;
x = y = 0;
do {
	do {
		if (patch_land[x][y] == -1) { // create new patch
			patch_land[x][y] = patchnum;
			cellCount = 1;
#if RSDEBUG
//DebugGUI((Int2Str(0)).c_str());
//DebugGUI(("LandSeriesPatches(): x=" + Int2Str(x) + " y=" + Int2Str(y)
//	+ " patchnum=" + Int2Str(patchnum)
//	+ " cellCount=" + Int2Str(cellCount)
//	).c_str());
#endif
			// expand patch to adjacent suitable cells
			AddNeighbours(x,y,patchnum);
			patchnum++;
		}
		x++;
	} while (x < ppLand.dimX);
	x = 0;
	y++;
} while (y < ppLand.dimY);

name = Inputs + patchMapName;
Fmaps.open(name.c_str());

// headers for raster format
Fmaps << "ncols\t" << ppLand.dimX << endl;
Fmaps << "nrows\t" << ppLand.dimY << endl;
Fmaps << "xllcorner\t" << 0.0 << endl;
Fmaps << "yllcorner\t" << 0.0 << endl;
Fmaps << "cellsize\t" << ppLand.resol << endl;
Fmaps << "NODATA_value\t" << -9999 << endl;
for (int y = ppLand.dimY-1; y > -1; y--){
	firstcolumn = true;
	for (int x = 0; x < ppLand.dimX; x++) {
		if (firstcolumn) firstcolumn = false;
		else Fmaps << "\t";;
		Fmaps << patch_land[x][y];
	}
	Fmaps << endl;
}
Fmaps.close(); Fmaps.clear();

}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

