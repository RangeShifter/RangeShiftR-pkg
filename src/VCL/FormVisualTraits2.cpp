//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "FormVisualTraits2.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TfrmVisualTraits2 *frmVisualTraits2;
//---------------------------------------------------------------------------
__fastcall TfrmVisualTraits2::TfrmVisualTraits2(TComponent* Owner)
	: TForm(Owner)
{
 frmVisualTraits2->Left = 200;
 frmVisualTraits2->Top = 200;
}

void __fastcall TfrmVisualTraits2::refresh(int height) {

// reset images
Image0->Picture->Bitmap = NULL;
Image1->Picture->Bitmap = NULL;
Image2->Picture->Bitmap = NULL;
Image3->Picture->Bitmap = NULL;
Image4->Picture->Bitmap = NULL;
Image5->Picture->Bitmap = NULL;

// reset legends
Image0L->Picture->Bitmap = NULL;
Image1L->Picture->Bitmap = NULL;
Image2L->Picture->Bitmap = NULL;
Image3L->Picture->Bitmap = NULL;
Image4L->Picture->Bitmap = NULL;
Image5L->Picture->Bitmap = NULL;

// set image heights
Image0->Height = height;
Image1->Height = height;
Image2->Height = height;
Image3->Height = height;
Image4->Height = height;
Image5->Height = height;
VertScrollBar->Range = height + 200;

// set captions of charts to blank
Label0->Caption = " ";
Label1->Caption = " ";
Label2->Caption = " ";
Label3->Caption = " ";
Label4->Caption = " ";
Label5->Caption = " ";

// remove legend labels
Label0L->Visible = false; Label0R->Visible = false;
Label1L->Visible = false; Label1R->Visible = false;
Label2L->Visible = false; Label2R->Visible = false;
Label3L->Visible = false; Label3R->Visible = false;
Label4L->Visible = false; Label4R->Visible = false;
Label5L->Visible = false; Label5R->Visible = false;

}

// Set labels and legend for a specified trait image
void __fastcall TfrmVisualTraits2::setImage(int panel,UnicodeString titletext,
	double mean,double sd,float nsd,float min,float max,Graphics::TBitmap *map)
{
TImage *image;
TLabel *title;
TLabel *left;
TLabel *right;
double minLabel,maxLabel;

switch (panel) {
case 0:
	title = Label0; left = Label0L; right = Label0R; image = Image0L;
	break;
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

// Redraw a specified trait image
void __fastcall TfrmVisualTraits2::drawImage(int panel,Graphics::TBitmap *map) {
TImage *image;

switch (panel) {
case 0:
	image = Image0;
	break;
case 1:
	image = Image1;
	break;
case 2:
	image = Image2;
	break;
case 3:
	image = Image3;
	break;
case 4:
	image = Image4;
	break;
case 5:
	image = Image5;
	break;
default:
	return;
}

image->Picture->Bitmap = map;
image->Repaint();
}

//---------------------------------------------------------------------------
