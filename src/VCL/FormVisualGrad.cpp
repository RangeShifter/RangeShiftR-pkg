//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "FormVisualGrad.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TfrmVisualGrad *frmVisualGrad;
//---------------------------------------------------------------------------
__fastcall TfrmVisualGrad::TfrmVisualGrad(TComponent* Owner)
	: TForm(Owner)
{
 frmVisualGrad->Left = 200;
 frmVisualGrad->Top = 200;
}

void __fastcall TfrmVisualGrad::refresh(int height) {
// set image heights
Image1->Height = height;
Image2->Height = height;
Image3->Height = height;
Image4->Height = height;
Image5->Height = height;
Image6->Height = height;
Image7->Height = height;
VertScrollBar->Range = height + 200;

// set captions of charts to blank
Label1->Caption = " ";
Label2->Caption = " ";
Label3->Caption = " ";
Label4->Caption = " ";
Label5->Caption = " ";
Label6->Caption = " ";
Label7->Caption = " ";
}

// Set labels and legend for a specified trait image
void __fastcall TfrmVisualGrad::setImage(int panel,UnicodeString titletext,
	double mean,double sd,float nsd,float min,float max,Graphics::TBitmap *map)
{
TImage *image;
TLabel *title;
TLabel *left;
TLabel *right;
double minLabel,maxLabel;

switch (panel) {
case 1:
	title = Label1; left = Label1L; right = Label1R; image = Image1L;
	break;
case 2:
	title = Label2; left = Label2L; right = Label2R; image = Image2L;
	break;
case 3:
	title = Label3; left = Label3L; right = Label3R; image = Image3L;
	break;
case 4:
	title = Label4; left = Label4L; right = Label4R; image = Image4L;
	break;
case 5:
	title = Label5; left = Label5L; right = Label5R; image = Image5L;
	break;
case 6:
	title = Label6; left = Label6L; right = Label6R; image = Image6L;
	break;
case 7:
	title = Label7; left = Label7L; right = Label7R; image = Image7L;
	break;
default:
	return;
}

title->Caption = titletext;
title->Visible = true; title->Width = 150;

minLabel = (double)((int)(1000.0 * (mean - nsd*sd) + 0.5)) / 1000.0;
if (max > min) { if (minLabel < min) minLabel = min; }
maxLabel = (double)((int)(1000.0 * (mean + nsd*sd) + 0.5)) / 1000.0;
if (max > min) { if (maxLabel > max) maxLabel = max; }
left->Visible = true; right->Visible = true;
left->Caption = FloatToStr(minLabel).c_str();
right->Caption = FloatToStr(maxLabel).c_str();
image->Picture->Bitmap = map;
image->Repaint();

}
//---------------------------------------------------------------------------
