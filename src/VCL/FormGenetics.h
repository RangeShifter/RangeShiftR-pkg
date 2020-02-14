/*------------------------------------------------------------------------------

RangeShifter v2.0 FormGenetics

Input from the Genetics Parameters form

For full details of RangeShifter, please see:
Bocedi G., Palmer S.C.F., Pe’er G., Heikkinen R.K., Matsinos Y.G., Watts K.
and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
eco-evolutionary dynamics and species’ responses to environmental changes.
Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

Authors: Greta Bocedi & Steve Palmer, University of Aberdeen

Last updated: 13 February 2020 by Steve Palmer

------------------------------------------------------------------------------*/

#ifndef FormGeneticsH
#define FormGeneticsH
//---------------------------------------------------------------------------

#include <System.Classes.hpp>
#include <Vcl.Controls.hpp>
#include <Vcl.StdCtrls.hpp>
#include <Vcl.Forms.hpp>
#include <Vcl.ExtCtrls.hpp>

#include "Genome.h"
#include "Parameters.h"
#include <Vcl.Dialogs.hpp>

//---------------------------------------------------------------------------
class TfrmGenetics : public TForm
{
__published:	// IDE-managed Components
	TLabeledEdit *edtSize;
	TLabeledEdit *edtMutnProb;
	TLabeledEdit *edtXoverProb;
	TLabeledEdit *edtAlleleSD;
	TLabeledEdit *edtMutnSD;
	TLabel *LabelChromosome;
	TButton *BtnOK;
	TRadioGroup *RGploidy;
	TLabeledEdit *nTraits;
	TRadioGroup *RGarchitecture;
	TButton *BtnReadFile;
	TOpenDialog *OpenDialog1;
	TCheckBox *CBneutral;
	void __fastcall BtnOKClick(TObject *Sender);
	void __fastcall RGarchitectureClick(TObject *Sender);
	void __fastcall BtnReadFileClick(TObject *Sender);
	void __fastcall CBneutralClick(TObject *Sender);
private:	// User declarations
public:		// User declarations
	__fastcall TfrmGenetics(TComponent* Owner);
	void __fastcall neutralGenetics(void);
	void __fastcall refresh(bool);
	int __fastcall checkArchFile(string);
	void __fastcall archFormatError(void);
};
//---------------------------------------------------------------------------
extern PACKAGE TfrmGenetics *frmGenetics;
//---------------------------------------------------------------------------

extern Species *pSpecies;			// see FormMain.cpp
extern int fileNtraits;				// see BatchMode.cpp

extern bool forceInit;
extern bool geneticsOK;
extern void MemoLine(UnicodeString);

#endif
