//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "FormDynLand.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TfrmDynLand *frmDynLand;

int chgnum,chgyear,maxchgyr;

//---------------------------------------------------------------------------
__fastcall TfrmDynLand::TfrmDynLand(TComponent* Owner)
	: TForm(Owner)
{
}

//---------------------------------------------------------------------------
void __fastcall TfrmDynLand::refresh(void)
{
chgnum = chgyear = 1;
maxchgyr = 0;
edtChange->Text = chgnum;
edtYear->Text = chgyear;
}

//---------------------------------------------------------------------------
void __fastcall TfrmDynLand::BtnAddClick(TObject *Sender)
{
landChange chg;
rasterdata habraster,pchraster;
string habchgmapname,pchchgmapname;
UnicodeString msg;
//string msgResol = "The raster resolution is different from the set landscape resolution.";
landParams ppLand = pLandscape->getLandParams();
landOrigin origin = pLandscape->getOrigin();

chgyear = StrToInt(edtYear->Text);
if (chgyear <= maxchgyr) {
	MessageDlg("Year must be greater than year of previous change",
		mtError, TMsgDlgButtons() << mbOK,0);
	return;
}

OpenDialog1->InitialDir = paramsSim->getDir(1).c_str();
msg = "Select habitat map for change no. " + edtChange->Text;
OpenDialog1->Title = msg.c_str();
if (OpenDialog1->Execute()) {
	habchgmapname = AnsiString(OpenDialog1->FileName).c_str();
	habraster = CheckRasterFile(habchgmapname);
	if (habraster.ok) {
		if (habraster.cellsize != ppLand.resol) {
			habraster.ok = false;
			MessageDlg(msgResol.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
		}
		else {
			if (habraster.ncols == ppLand.dimX
			&&  habraster.nrows == ppLand.dimY
			&&  habraster.xllcorner == origin.minEast
			&&  habraster.yllcorner == origin.minNorth) {
				// habitat raster matches original landscape - no action yet
			}
			else {
				habraster.ok = false;
				MessageDlg(msgHeaders.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
			}
		}
	}
	else {
		habraster.ok = false;
		MessageDlg(msgLandHdr.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
	}
}
else {
	habraster.ok = false;
	MessageDlg("Selection of habitat file cancelled",
		mtWarning, TMsgDlgButtons() << mbOK,0);
}
if (!habraster.ok) return;

if (ppLand.patchModel) {
	// get the name of the patch raster
	msg = "Habitat raster OK - now specify corresponding patch IDs file ...";
	MessageDlg(msg.c_str(),mtWarning, TMsgDlgButtons() << mbOK,0);
	pchraster.ok = false;
	OpenDialog1->InitialDir = paramsSim->getDir(1).c_str();
	OpenDialog1->Title = "Select patch raster map";
	if (OpenDialog1->Execute()){
		pchchgmapname = AnsiString(OpenDialog1->FileName).c_str();
		pchraster = CheckRasterFile(pchchgmapname);
		if (pchraster.ok) {
			if (pchraster.cellsize != ppLand.resol) {
				pchraster.ok = false;
				MessageDlg(msgResol.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
			}
			if (pchraster.nrows != habraster.nrows || pchraster.ncols != habraster.ncols) {
				pchraster.ok = false;
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
	if (!pchraster.ok) return;
}

if (chgnum == 1) {
	pLandscape->setDynamicLand(true);
}
maxchgyr = chgyear;
chg.chgnum = chgnum; chg.chgyear = chgyear;
chg.habfile = habchgmapname;
chg.pchfile = pchchgmapname;
chg.costfile = "NULL";
pLandscape->addLandChange(chg);
chgnum++; chgyear++;
edtChange->Text = chgnum;
edtYear->Text = chgyear;

}

//---------------------------------------------------------------------------
void __fastcall TfrmDynLand::BtnResetClick(TObject *Sender)
{
pLandscape->deleteLandChanges();
//pLandscape->setDynamicLand(false,false);
pLandscape->setDynamicLand(false);
refresh();
}

//---------------------------------------------------------------------------
// Read the selected landscape changes
void __fastcall TfrmDynLand::BtnFinishedClick(TObject *Sender)
{
landParams ppLand = pLandscape->getLandParams();
int imported;
string msg;
if (ppLand.patchModel) {
	MemoLine("Creating patch change matrix...");
	pLandscape->createPatchChgMatrix();
	MemoLine("...finished");
}
MemoLine("Recording landscape changes...");
int nchanges = pLandscape->numLandChanges();
for (int i = 0; i < nchanges; i++) {
	MemoLine(("..." + Int2Str(i+1) + "...").c_str());
	imported = pLandscape->readLandChange(i,false);
	if (imported != 0) {
		msg = msgReadError + "change " + Int2Str(i+1) + " file: " + Int2Str(imported);
		MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
		return;
	}
	if (ppLand.patchModel) {
		pLandscape->recordPatchChanges(i+1);
	}
}
if (ppLand.patchModel) {
	// record changes back to original landscape for multiple replicates
	pLandscape->recordPatchChanges(0);
	pLandscape->deletePatchChgMatrix();
}
MemoLine("...finished");
Close();
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


