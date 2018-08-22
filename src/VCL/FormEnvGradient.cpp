//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "FormEnvGradient.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TfrmEnvGradient *frmEnvGradient;

//---------------------------------------------------------------------------

__fastcall TfrmEnvGradient::TfrmEnvGradient(TComponent* Owner)
	: TForm(Owner)
{
frmEnvGradient->Top = 60;
frmEnvGradient->Left = 100;
}

//-----------------------------------------------------------------------------
// Refresh the form from retained parameter values
void __fastcall TfrmEnvGradient::Refresh(void)
{
envGradParams g;
landParams ppLand = pLandscape->getLandParams();

g = paramsGrad->getGradient();
if (g.gradient) {
	RGgradientType->ItemIndex = g.gradType;
	edtGradInc->Enabled = true;
	edtOptRow->Enabled = true;
	edtDev->Enabled = true;
	CBshift->Enabled = true;
	if (ppLand.patchModel){
		edtOptEXT->Enabled = false;
	}
	else {
		edtOptEXT->Enabled = true;
	}
	if (g.shifting) {
		CBshift->Checked = true;
		edtShiftRate->Enabled = true;
		edtStartShift->Enabled = true;
		edtStopShift->Enabled = true;
	}
	else {
		CBshift->Checked = false;
		edtShiftRate->Enabled = false;
		edtStartShift->Enabled = false;
		edtStopShift->Enabled = false;
	}
}
else {
	RGgradientType->ItemIndex = 0;
	edtGradInc->Enabled = false;
	edtOptRow->Enabled = false;
	edtDev->Enabled = false;
	edtOptEXT->Enabled = false;
	CBshift->Enabled = false;
	CBshift->Checked = false;
	edtShiftRate->Enabled = false;
	edtStartShift->Enabled = false;
	edtStopShift->Enabled = false;
}

}

//-----------------------------------------------------------------------------

void __fastcall TfrmEnvGradient::CBshiftClick(TObject *Sender)
{
if (CBshift->Checked){
	edtShiftRate->Enabled = true;
	edtStartShift->Enabled = true;
	edtStopShift->Enabled = true;
}
else{
	edtShiftRate->Enabled = false;
	edtStartShift->Enabled = false;
	edtStopShift->Enabled = false;
}
}

//---------------------------------------------------------------------------
/*
// SCFP 2/11/13 TEMPORARY
const string Int22Str(const int x)
{
	ostringstream o;
	if (!(o << x)) return "ERROR";
	return o.str();
}
*/

//---------------------------------------------------------------------------

void __fastcall TfrmEnvGradient::RGgradientTypeClick(TObject *Sender)
{
landParams ppLand = pLandscape->getLandParams();
genLandParams ppGenLand = pLandscape->getGenLandParams();

/*
if (RGgradientType->ItemIndex == 1) { // habitat suitability
//	string msg;
//	msg = Int22Str(SimType) + " " + Int22Str(fractal) + " " + Int22Str(continuous);
//		MessageDlg(msg.c_str(),mtWarning, TMsgDlgButtons() << mbRetry,0);
	if (!ppLand.generated || ppGenLand.fractal || ppGenLand.continuous) {
		MessageDlg("A gradient in habitat suitability is only possible with"
			" discrete random landscapes generated at each replicate.",
			mtError, TMsgDlgButtons() << mbRetry,0);
		RGgradientType->ItemIndex = 0;
	}
}
*/

if (RGgradientType->ItemIndex == 0) { // no gradient
	edtGradInc->Enabled = false;
	edtOptRow->Enabled = false;
	edtDev->Enabled = false;
	edtOptEXT->Enabled = false;
	CBshift->Enabled = false;
	CBshift->Checked = false;
	edtShiftRate->Enabled = false;
	edtStartShift->Enabled = false;
	edtStopShift->Enabled = false;
}
else { // gradient
	edtGradInc->Enabled = true;
	edtOptRow->Enabled = true;
	edtDev->Enabled = true;
	frmEnvGradient->CBshift->Enabled = true;
	if (ppLand.patchModel){
		edtOptEXT->Enabled = false;
	}
	else {
		edtOptEXT->Enabled = true;
	}
}

//if (RGgradientType->ItemIndex == 4)
if (RGgradientType->ItemIndex == 3)
{ // gradient in local extinction probability
	edtOptEXT->Visible = true;  edtOptEXT->Enabled = true;
}
else {
	edtOptEXT->Visible = false; edtOptEXT->Enabled = false;
}
}

//----------------------------------------------------------------------------
//Read parameters for environmental gradient

void __fastcall TfrmEnvGradient::BtnOKClick(TObject *Sender)
{
int gradType,shift_begin,shift_stop;
bool gradient;
float gradSteep,optY,f,optEXT,shift_rate;

if (RGgradientType->ItemIndex == 0) gradient = false; else gradient = true;
if (gradient) {
	gradType  = RGgradientType->ItemIndex;
	gradSteep = StrToFloat(edtGradInc->Text);
	optY = StrToFloat(edtOptRow->Text); // opt_y0 = opt_y;
	f = StrToFloat(edtDev->Text);
	optEXT = StrToFloat(edtOptEXT->Text);
//	if (grad_inc < 0.0)
	if (gradSteep < 0.0 || gradSteep > 1.0)
	{
//		MessageDlg("Gradient steepness may not be less than zero",
		MessageDlg("Gradient steepness must lie between 0 and 1",
			mtError, TMsgDlgButtons() << mbRetry,0);
		return;
	}
	if (optY < 0) {
		MessageDlg("Optimum Y may not be less than zero",
			mtError, TMsgDlgButtons() << mbRetry,0);
		return;
	}
	if (f < 0.0) {
		MessageDlg("Scaling factor may not be less than zero",
			mtError, TMsgDlgButtons() << mbRetry,0);
		return;
	}
//	if (RGgradientType->ItemIndex == 4)
	if (RGgradientType->ItemIndex == 3)
	{
		if (optEXT < 0.0 || optEXT >= 1.0) {
			MessageDlg("Extinction probability must be in the range  0.0 <= P < 1.0",
				mtError, TMsgDlgButtons() << mbRetry,0);
			return;
		}
	}
	paramsGrad->setGradient(gradType,gradSteep,optY,f,optEXT);
	if (CBshift->Checked == true) {
		shift_rate = StrToFloat(edtShiftRate->Text);
		shift_begin = StrToInt(edtStartShift->Text);
		shift_stop = StrToInt(edtStopShift->Text);
		if (shift_rate <= 0.0) {
			MessageDlg("Shifting rate must be greater than zero",
				mtError, TMsgDlgButtons() << mbRetry,0);
			return;
		}
		if (shift_begin < 0) {
			MessageDlg("Starting year may not be less than zero",
				mtError, TMsgDlgButtons() << mbRetry,0);
			return;
		}
		if (shift_stop <= shift_begin) {
			MessageDlg("Stopping year must be greater than starting year",
				mtError, TMsgDlgButtons() << mbRetry,0);
			return;
		}
		paramsGrad->setShifting(shift_rate,shift_begin,shift_stop);
	}
	else paramsGrad->noShifting();
}
else {
	paramsGrad->noGradient();
	paramsGrad->noShifting();
}

Close();
}

//---------------------------------------------------------------------------
void __fastcall TfrmEnvGradient::BtnCancelClick(TObject *Sender)
{
Close();
}

//---------------------------------------------------------------------------

