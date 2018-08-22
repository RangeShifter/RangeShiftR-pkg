//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "FormVisualCost.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TfrmVisualCost *frmVisualCost;
//---------------------------------------------------------------------------
__fastcall TfrmVisualCost::TfrmVisualCost(TComponent* Owner)
	: TForm(Owner)
{
 frmVisualCost->Left = 200;
 frmVisualCost->Top = 200;
}
//---------------------------------------------------------------------------
