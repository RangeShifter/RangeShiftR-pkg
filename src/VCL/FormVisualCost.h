/*------------------------------------------------------------------------------

RangeShifter v2.0 FormVisualCost

Displays the Cost Landscape form

For full details of RangeShifter, please see:
Bocedi G., Palmer S.C.F., Pe’er G., Heikkinen R.K., Matsinos Y.G., Watts K.
and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
eco-evolutionary dynamics and species’ responses to environmental changes.
Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

Authors: Greta Bocedi & Steve Palmer, University of Aberdeen

Last updated: 12 February 2016 by Steve Palmer

------------------------------------------------------------------------------*/

#ifndef FormVisualCostH
#define FormVisualCostH
//---------------------------------------------------------------------------
#include <System.Classes.hpp>
#include <Vcl.Controls.hpp>
#include <Vcl.StdCtrls.hpp>
#include <Vcl.Forms.hpp>
#include <Vcl.ExtCtrls.hpp>
//---------------------------------------------------------------------------
class TfrmVisualCost : public TForm
{
__published:	// IDE-managed Components
	TScrollBox *CostScrollBox;
	TImage *CostImage;
	TPaintBox *PaintBox1;
private:	// User declarations
public:		// User declarations
	__fastcall TfrmVisualCost(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TfrmVisualCost *frmVisualCost;
//---------------------------------------------------------------------------
#endif
