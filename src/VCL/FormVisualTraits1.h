/*------------------------------------------------------------------------------

RangeShifter v2.0 FormVisualTraits1H

Displays Traits form 1

For full details of RangeShifter, please see:
Bocedi G., Palmer S.C.F., Pe’er G., Heikkinen R.K., Matsinos Y.G., Watts K.
and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
eco-evolutionary dynamics and species’ responses to environmental changes.
Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

Authors: Greta Bocedi & Steve Palmer, University of Aberdeen

Last updated: 12 February 2016 by Steve Palmer

------------------------------------------------------------------------------*/

#ifndef FormVisualTraits1H
#define FormVisualTraits1H
//---------------------------------------------------------------------------
#include <System.Classes.hpp>
#include <Vcl.Controls.hpp>
#include <Vcl.StdCtrls.hpp>
#include <Vcl.Forms.hpp>
#include <Vcl.ExtCtrls.hpp>
//---------------------------------------------------------------------------
class TfrmVisualTraits1 : public TForm
{
__published:	// IDE-managed Components
	TImage *Image0;
	TLabel *Label0;
	TImage *Image1;
	TLabel *Label1;
	TImage *Image2;
	TLabel *Label2;
	TImage *Image3;
	TLabel *Label3;
	TImage *Image4;
	TLabel *Label4;
	TImage *Image5;
	TLabel *Label5;
	TLabel *Label0L;
	TLabel *Label0R;
	TImage *Image0L;
	TLabel *Label1L;
	TImage *Image1L;
	TLabel *Label1R;
	TImage *Image2L;
	TLabel *Label2R;
	TLabel *Label2L;
	TImage *Image3L;
	TLabel *Label3R;
	TLabel *Label3L;
	TLabel *Label4L;
	TLabel *Label4R;
	TImage *Image4L;
	TLabel *Label5L;
	TLabel *Label5R;
	TImage *Image5L;
private:	// User declarations
public:		// User declarations
	__fastcall TfrmVisualTraits1(TComponent* Owner);
	void __fastcall refresh(int);
	void __fastcall setImage( // Set labels and legend for a specified trait image
		int,								// image number (0...5)
		UnicodeString,			// title text
		double,							// mean initial trait value
		double,							// s.d. initial trait value
		float,							// no. of s.d. used to determine min and max legend values
		float,							// min value ) only applied if
		float,							// max value )  max > min
		Graphics::TBitmap*	// pointer to legend image
	);
	void __fastcall drawImage( // Redraw a specified trait image
		int,								// image number (0...5)
		Graphics::TBitmap*	// pointer to trait image
	);
};
//---------------------------------------------------------------------------
extern PACKAGE TfrmVisualTraits1 *frmVisualTraits1;
//---------------------------------------------------------------------------
#endif
