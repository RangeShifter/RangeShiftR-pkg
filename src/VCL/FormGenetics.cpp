//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "FormGenetics.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TfrmGenetics *frmGenetics;

string genfilename;
bool geneticfileOK = false;
string msgfile = "genetic architecture file ";

//---------------------------------------------------------------------------
__fastcall TfrmGenetics::TfrmGenetics(TComponent* Owner)
	: TForm(Owner)
{
Left = 100; Top = 60;
}

//---------------------------------------------------------------------------
// Refresh the form before display
void __fastcall TfrmGenetics::refresh(bool sexualSpecies)
{
if (sexualSpecies) { // species must be diploid
	RGploidy->ItemIndex = 0; RGploidy->Enabled = false;
}
else { // species may be haploid or self-fertilising
//	RGploidy->Enabled = true;
	// currently there is no mating system for hermaphrodites
	RGploidy->ItemIndex = 1; RGploidy->Enabled = false;
}
int ntraits = pSpecies->getNTraits();
nTraits->Text = ntraits;
if (ntraits == 0) {
	edtSize->Enabled = false;
	RGarchitecture->ItemIndex = 1; RGarchitecture->Enabled = false;
	CBneutral->Enabled = true;  CBneutral->Visible = true;
	neutralGenetics();
}
else {
	CBneutral->Enabled = false; CBneutral->Visible = false;
	CBneutral->Checked = false;
	RGarchitecture->Enabled = true;
	BtnReadFile->Enabled = true;
	edtMutnProb->Enabled = true;  edtXoverProb->Enabled = true;
	edtAlleleSD->Enabled = true;  edtMutnSD->Enabled = true;
	if (RGarchitecture->ItemIndex == 0) { // one chromosome per trait
		edtSize->Enabled = true;  BtnReadFile->Enabled = false;
	}
	else { // read from file
		edtSize->Enabled = false; BtnReadFile->Enabled = true;
	}
}
}

//---------------------------------------------------------------------------
void __fastcall TfrmGenetics::CBneutralClick(TObject *Sender)
{
neutralGenetics();
}

//---------------------------------------------------------------------------
void __fastcall TfrmGenetics::neutralGenetics(void)
{
if (CBneutral->Checked) {
	BtnReadFile->Enabled = true;
	edtMutnProb->Enabled = true;  edtXoverProb->Enabled = true;
	edtAlleleSD->Enabled = true;  edtMutnSD->Enabled = true;
}
else {
	BtnReadFile->Enabled = false;
	edtMutnProb->Enabled = false;  edtXoverProb->Enabled = false;
	edtAlleleSD->Enabled = false;  edtMutnSD->Enabled = false;
}
}

//---------------------------------------------------------------------------
void __fastcall TfrmGenetics::RGarchitectureClick(TObject *Sender)
{
if (RGarchitecture->ItemIndex == 0) { // one chromosome per trait
	edtSize->Enabled = true;  BtnReadFile->Enabled = false;
}
else { // read from file
	edtSize->Enabled = false; BtnReadFile->Enabled = true;
}
}

//---------------------------------------------------------------------------
void __fastcall TfrmGenetics::BtnReadFileClick(TObject *Sender)
{
bool fileSelected = false;
string msg;

OpenDialog1->InitialDir = paramsSim->getDir(1).c_str();
OpenDialog1->Title = "Select genetic architecture file";
if (OpenDialog1->Execute()) {
	genfilename = AnsiString(OpenDialog1->FileName).c_str();
	fileSelected = true;
}
else {
	msg = "The " + msgfile + "selection has been cancelled";
	MessageDlg(msg.c_str(),mtWarning, TMsgDlgButtons() << mbOK,0);
	return;
}
fileNtraits = 0;
if (fileSelected) {
	// check selected file
	geneticfileOK = false;
	fileNtraits = checkArchFile(genfilename);
	if (fileNtraits >= 0) {
		geneticfileOK = true;
	}
	if (fileNtraits == 0 && CBneutral->Checked) {
		geneticfileOK = true;
	}
}

if (geneticfileOK) { // read file and set up trait mapping to genome
	ifstream archfile;
	archfile.open(genfilename.c_str());
	string paramname;
	int nchromosomes,nloci,traitnum,chrom,locus;
	// set no. of chromosomes
	archfile >> paramname >> nchromosomes;
	pSpecies->setNChromosomes(nchromosomes);
	int nchromset = pSpecies->getNChromosomes();
#if RSDEBUG
	msg = "No. of chromosomes set is " + Int2Str(nchromset);
	MessageDlg(msg.c_str(),mtWarning, TMsgDlgButtons() << mbOK,0);
#endif
	if (nchromset <= 0) geneticfileOK = false;
	pSpecies->setTraitData(fileNtraits);
	// set no. of loci for each chromosome
	archfile >> paramname;
	for (int i = 0; i < nchromosomes; i++) {
		archfile >> nloci;
		pSpecies->setNLoci(i,nloci);
	}
	if (!CBneutral->Checked) {
		// set trait maps
		paramname = "XXXyyyZZZ";
		archfile >> paramname;
		do {
			archfile >> traitnum >> paramname >> nloci;
			pSpecies->setTraitMap(traitnum,nloci);
			for (int allele = 0; allele < nloci; allele++) {
				chrom = locus = -999999;
				archfile >> chrom >> locus;
				pSpecies->setTraitAllele(traitnum,allele,chrom,locus);
			}
			paramname = "XXXyyyZZZ";
			archfile >> paramname;
		} while (paramname != "XXXyyyZZZ");
	}
	archfile.close(); archfile.clear();
	// NB neural markers are established when OK button is clicked, in case no. of
	// traits has been changed but new file has not been selected
}

if (!geneticfileOK) {
	msg = "Problem encountered in setting up trait mapping ";
	MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
	return;
}

}

//---------------------------------------------------------------------------
int __fastcall TfrmGenetics::checkArchFile(string fname)
{
string paramname,msg;
int nchromosomes,nloci;
bool formatError = false;
int *chromsize;

ifstream archfile;

// open specified file
archfile.open(fname.c_str());
if (!archfile.is_open()) return 11;

// check no. of chromosomes, and terminate if in error
archfile >> paramname >> nchromosomes;
if (paramname == "NChromosomes") {
	if (nchromosomes < 1) {
		msg = "No. of chromosomes must be >=1 ";
		MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
		archfile.close(); archfile.clear();
		return -1;
	}
}
else {
	archFormatError();
	archfile.close(); archfile.clear();
	return -1;
}
chromsize = new int[nchromosomes];
for (int i = 0; i < nchromosomes; i++) chromsize[i] = 0;

// check no. of loci on each chromosome, and terminate if in error
archfile >> paramname;
if (paramname != "NLoci") {
	archFormatError();
	archfile.close(); archfile.clear();
	return -1;
}
int locerrors = 0;
for (int i = 0; i < nchromosomes; i++) {
	nloci = -999;
	archfile >> nloci;
	if (nloci < 1) locerrors++; else chromsize[i] = nloci;
}
if (locerrors) {
	msg = "No. of loci must be >=1 ";
	MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
	archfile.close(); archfile.clear();
	return -1;
}

if (CBneutral->Checked) { // neutral genetics only
	// no further checking required
	return 0;
}

// check unspecified no. of traits
int traitnum,prevtrait,chrom,locus;
traitnum = prevtrait = -1;
int ntraits = 0;
bool traitError = false;
bool lociError = false;
bool chromError = false;
bool locusError = false;
paramname = "XXXyyyZZZ";
archfile >> paramname;
if (paramname != "Trait") formatError = true;
do {
	archfile >> traitnum;
	if (paramname != "Trait") formatError = true;
	if (traitnum == (prevtrait+1)) prevtrait = traitnum;
	else traitError = true;
	archfile >> paramname >> nloci;
	if (paramname != "NLoci") formatError = true;
	if (nloci < 1) lociError = true;
	for (int i = 0; i < nloci; i++) {
		chrom = locus = -999999;
		archfile >> chrom >> locus;
		if (chrom == -999999 || locus == -999999) {
			msg = "Too few loci listed for trait " + Int2Str(traitnum);
			MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
			archfile.close(); archfile.clear();
			return -1;
		}
		else {
			if (chrom >= 0 && chrom < nchromosomes) {
				if (locus < 0 || locus >= chromsize[chrom]) locusError = true;
			}
			else chromError = true;
		}
	}
	ntraits++;
	paramname = "XXXyyyZZZ";
	archfile >> paramname;
} while (paramname != "XXXyyyZZZ");

if (formatError) {
	archFormatError();
	archfile.close(); archfile.clear();
	return -1;
}
if (traitError) {
	msg = "Traits must be sequentially numbered starting at 0 ";
	MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
	archfile.close(); archfile.clear();
	return -1;
}
if (lociError) {
	msg = "No. of trait loci must be >=1 ";
	MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
	archfile.close(); archfile.clear();
	return -1;
}
if (chromError) {
	msg = "Chromosome no. must be from 0 to " + Int2Str(nchromosomes-1);
	MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
	archfile.close(); archfile.clear();
	return -1;
}
if (locusError) {
	msg = "Locus no. must not exceed no. of loci on specified chromosome ";
	MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
	archfile.close(); archfile.clear();
	return -1;
}

// final read should hit EOF

if (!archfile.eof()) {
	msg = "Failed to read to end of file ";
	MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
	archfile.close(); archfile.clear();
	return -1;
}

archfile.close(); archfile.clear();

return ntraits;

}

void __fastcall TfrmGenetics::archFormatError(void)
{
string msg = "Format error in architecture file: case-sensitive parameter names "
	"must match the specification exactly";
MessageDlg(msg.c_str(),mtError, TMsgDlgButtons() << mbOK,0);
}

//---------------------------------------------------------------------------
void __fastcall TfrmGenetics::BtnOKClick(TObject *Sender)
{
//MemoLine("GENETICS FORM - start of btnOKClick()");

string msg;
string msg01 = "probability must be between 0 and 1 ";
string msggt0 = "must be greater than 0 ";

if (RGarchitecture->ItemIndex == 0) { // one chromosome per trait
	if (StrToInt(edtSize->Text) < 1) {
		msg = "No. of loci per chromosome " + msggt0;
		MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
		return;
	}
}
else { // read from file
	if (!geneticfileOK) {
		msg = "Please select valid " + msgfile;
		MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
		return;
	}
	int ntraits = pSpecies->getNTraits();
	if (fileNtraits != ntraits) {
		msg = "No. of traits defined in " + msgfile
			+ "does not match no. set for the species";
		MessageDlg(msg.c_str(),mtWarning, TMsgDlgButtons() << mbOK,0);
	}
}
if (StrToFloat(edtMutnProb->Text) < 0.0 || StrToFloat(edtMutnProb->Text) > 1.0) {
	msg = "Mutation " + msg01;
	MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
	return;
}
if (StrToFloat(edtXoverProb->Text) < 0.0 || StrToFloat(edtXoverProb->Text) > 1.0) {
	msg = "Crossover " + msg01;
	MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
	return;
}
if (StrToFloat(edtAlleleSD->Text) <= 0.0) {
	msg = "Allele s.d. " + msggt0;
	MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
	return;
}
if (StrToFloat(edtMutnSD->Text) <= 0.0) {
	msg = "Mutation s.d. " + msggt0;
	MessageDlg(msg.c_str(),mtError,TMsgDlgButtons() << mbOK,0);
	return;
}

// Update genome data
genomeData d;
if (RGploidy->ItemIndex == 0) d.diploid = true; else d.diploid = false;
d.neutralMarkers = CBneutral->Checked;
if (RGarchitecture->ItemIndex == 0) {
	d.trait1Chromosome = true;
	d.nLoci = StrToInt(edtSize->Text);
	pSpecies->set1ChromPerTrait(d.nLoci);
}
else {
	d.trait1Chromosome = false;
	d.nLoci = 0;
}
d.probMutn = StrToFloat(edtMutnProb->Text);
d.probCrossover = StrToFloat(edtXoverProb->Text);
d.alleleSD = StrToFloat(edtAlleleSD->Text);
d.mutationSD = StrToFloat(edtMutnSD->Text);
pSpecies->setGenomeData(d);
if (!d.trait1Chromosome) { // architechture read from file
	// any loci not contributing to a trait are recorded as neutral
	pSpecies->setNeutralLoci(CBneutral->Checked);
}

geneticsOK = true;
//if (geneticsOK) MemoLine("COMPLETED EDIT CHECKS IN GENETICS FORM - geneticsOK");
//else MemoLine("COMPLETED EDIT CHECKS IN GENETICS FORM - NOT geneticsOK");

Close();

}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


