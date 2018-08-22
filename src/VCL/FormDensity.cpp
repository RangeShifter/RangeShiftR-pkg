//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "FormDensity.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TfrmDensity *frmDensity;
//---------------------------------------------------------------------------
__fastcall TfrmDensity::TfrmDensity(TComponent* Owner)
	: TForm(Owner)
{
 frmDensity->Left = 100;
 frmDensity->Top = 60;
}
//---------------------------------------------------------------------------
void __fastcall TfrmDensity::BtnOKClick(TObject *Sender)
{
 frmDensity->Close();
}
//---------------------------------------------------------------------------

void __fastcall TfrmDensity::BtnCancelClick(TObject *Sender)
{
 frmDensity->Close();
}
//---------------------------------------------------------------------------
