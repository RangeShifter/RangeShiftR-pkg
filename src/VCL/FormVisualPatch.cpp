//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "FormVisualPatch.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TfrmVisualPatch *frmVisualPatch;
//---------------------------------------------------------------------------
__fastcall TfrmVisualPatch::TfrmVisualPatch(TComponent* Owner)
	: TForm(Owner)
{
 frmVisualPatch->Left = 200;
 frmVisualPatch->Top = 200;
}

//---------------------------------------------------------------------------
void __fastcall TfrmVisualPatch::FormClose(TObject *Sender, TCloseAction &Action)
{
 frmLand->CBVisualPatch->Checked = false;
}

//---------------------------------------------------------------------------


