//---------------------------------------------------------------------------

#ifndef FormVisualGradH
#define FormVisualGradH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <ExtCtrls.hpp>
//---------------------------------------------------------------------------
class TfrmVisualGrad : public TForm
{
__published:	// IDE-managed Components
	TImage *Image1;
	TImage *Image2;
	TImage *Image3;
	TImage *Image4;
	TLabel *Label1;
	TLabel *Label2;
	TLabel *Label3;
	TLabel *Label4;
	TImage *Image5;
	TLabel *Label5;
	TImage *Image6;
	TLabel *Label6;
	TImage *Image7;
	TLabel *Label7;
	TLabel *Label1R;
	TLabel *Label1L;
	TImage *Image1L;
	TLabel *Label2R;
	TLabel *Label2L;
	TImage *Image2L;
	TLabel *Label3R;
	TLabel *Label3L;
	TImage *Image3L;
	TLabel *Label4L;
	TLabel *Label4R;
	TImage *Image4L;
	TLabel *Label5L;
	TLabel *Label5R;
	TImage *Image5L;
	TLabel *Label6L;
	TLabel *Label6R;
	TImage *Image6L;
	TLabel *Label7L;
	TLabel *Label7R;
	TImage *Image7L;
private:	// User declarations
public:		// User declarations
	__fastcall TfrmVisualGrad(TComponent* Owner);
	void __fastcall refresh(int);
	void __fastcall setImage( // Set labels and legend for a specified trait image
		int,								// image number (1...7)
		UnicodeString,			// title text
		double,							// mean initial trait value
		double,							// s.d. initial trait value
		float,							// no. of s.d. used to determine min and max legend values
		float,							// min value ) only applied if
		float,							// max value )  max > min
		Graphics::TBitmap*	// pointer to legend image
	);
};
//---------------------------------------------------------------------------
extern PACKAGE TfrmVisualGrad *frmVisualGrad;
//---------------------------------------------------------------------------
#endif
