/*----------------------------------------------------------------------------
 *	
 *	Copyright (C) 2020 Greta Bocedi, Stephen C.F. Palmer, Justin M.J. Travis, Anne-Kathleen Malchow, Damaris Zurell 
 *	
 *	This file is part of RangeShifter.
 *	
 *	RangeShifter is free software: you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *	the Free Software Foundation, either version 3 of the License, or
 *	(at your option) any later version.
 *	
 *	RangeShifter is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *	GNU General Public License for more details.
 *	
 *	You should have received a copy of the GNU General Public License
 *	along with RangeShifter. If not, see <https://www.gnu.org/licenses/>.
 *	
 --------------------------------------------------------------------------*/
 
 
//---------------------------------------------------------------------------
#if RS_EMBARCADERO
#pragma hdrstop
#endif

#include "BatchMode.h"
//---------------------------------------------------------------------------
#if RS_EMBARCADERO
#pragma package(smart_init) 
#endif

ifstream controlfile;
// Note - all batch files are prefixed 'b' here for reasons concerned with RS v1.0
ifstream bParamFile,bLandFile,bDynLandFile;
ifstream bSpDistFile,bStageStructFile,bTransMatrix;
ifstream bStageWeightsFile;
ifstream bEmigrationFile,bTransferFile,bSettlementFile;
ifstream bGeneticsFile,bArchFile,bInitFile,bInitIndsFile;
#if VIRTUALECOLOGIST
ifstream bVirtEcolFile,bSampleFile,bPatchFile;
#endif
#if RS_ABC
ifstream bABCpFile;
ifstream bABCoFile;
#endif
#if TEMPMORT
ifstream bMortFile;
#endif
#if SEASONAL
ifstream bSeasonFile;
//#if PARTMIGRN
ifstream bExtremeFile;
//#endif // PARTMIGRN 
#endif // SEASONAL
#if RS_CONTAIN
ifstream bHabDemFile,bManageFile;
#endif // RS_CONTAIN 

ofstream batchlog;

ofstream rsLog; // performance log for recording simulation times, etc.

// NOTE: THE STREAMS USED TO READ THE DATA AT RUN TIME COULD TAKE THE SAME NAMES AS
// USED DURING PARSING (ABOVE)
ifstream parameters;
ifstream ssfile,tmfile,fdfile,ddfile,sdfile;
#if RS_CONTAIN
ifstream hdfile,amfile;
#endif // RS_CONTAIN 
ifstream emigFile,transFile,settFile,genFile,archFile,initFile,initIndsFile;
ifstream landfile,dynlandfile;
#if TEMPMORT
ifstream mortFile;
#endif
#if SEASONAL
ifstream seasonfile;
//#if PARTMIGRN
ifstream extremefile;
//#endif // PARTMIGRN 
#endif // SEASONAL 

// global variables passed between parsing functions...
int batchnum;
int patchmodel,resolution,landtype,maxNhab,speciesdist,distresolution;
int reproductn;
#if SEASONAL
int nseasons;
#else
int repseasons;
#endif // SEASONAL 
int stagestruct,stages,transfer;
int sexesDem;		// no. of explicit sexes for demographic model
int sexesDisp;	// no. of explicit sexes for dispersal model
int firstsimul = 0;
int fileNtraits; // no. of traits defined in genetic architecture file
//rasterdata landraster,patchraster,spdistraster;
#if !RS_RCPP
rasterdata landraster;
#else
rasterdata landraster,patchraster,spdistraster,costsraster;
#endif
// ...including names of the input files
string parameterFile;
string landFile;
string name_landscape,name_patch,name_dynland,name_sp_dist,name_costfile;
#if RS_CONTAIN
string name_damagefile;
string name_managefile;
#endif // RS_CONTAIN 
#if SPATIALMORT
string name_mortfile[2];
#endif // SPATIALMORT 
#if SEASONAL
string seasonFile;
#endif // SEASONAL 
string stageStructFile,transMatrix;
string emigrationFile,transferFile,settleFile,geneticsFile,initialFile;
#if VIRTUALECOLOGIST
string virtEcolFile;
#endif // VIRTUALECOLOGIST 
#if RS_ABC
string abcParamsFile,abcObsFile;
int nABCsamples;
#endif // RS_ABC 
string prevInitialIndsFile = " ";

string msgnlines = "No. of lines for final Simulation ";
string msgshldbe = " should be ";
string msgresol0 = "*** Resolution of ";
string msgresol1 = " does not match Resolution in Control file ";
string msghdrs0  = "*** Headers of ";
string msghdrs1  = " do not match headers of LandscapeFile";
string msgpatch  = " is required for patch-based model";
string msgmatch  = " must match the specification exactly";
string msgcase   = " case-sensitive parameter names";

float **matrix = NULL;	// temporary matrix used in batch mode
int matrixsize = 0; 		// size of temporary matrix

//---------------------------------------------------------------------------
// Returns input value less next highest power of 2 (for x > 2)
int power2check(int x) {
if (x < 2) return 0;
int r = x%2;
while (r == 0) {
	x /= 2; r = x%2;
}
return x;
}

//---------------------------------------------------------------------------
batchfiles ParseControlFile(string ctrlfile, string indir, string outdir)
{
batchfiles b;
int lines,nSimuls;
int errors = 0;
string paramname,filename,fname,logname,header;
string filetype = "Control file";
bool controlFormatError = false;
b.ok = true; b.nSimuls = 0; b.nLandscapes = 0;

// open batch log file
logname = outdir + "BatchLog.txt";
batchlog.open(logname.c_str());
if (!batchlog.is_open()) {
//	MessageDlg("Error opening batch output log file",mtError, TMsgDlgButtons() << mbOK,0);
	cout << "Error opening batch output log file " << logname << endl;
	b.ok = false;
	return b;
}

controlfile.open(ctrlfile.c_str());

if (!controlfile.is_open()) {
//	MessageDlg("Error opening Control file",mtError, TMsgDlgButtons() << mbOK,0);
	cout << "Error opening Control file: " << ctrlfile << endl;
	batchlog << "Error opening Control file: " << ctrlfile << endl;
	b.ok = false;
	if (batchlog.is_open()) { batchlog.close(); batchlog.clear(); }
	return b;
}
else {
	batchlog << "Checking Control file " << ctrlfile << endl;
}

// Check fixed model parameters

controlfile >> paramname >> batchnum;
if (paramname == "BatchNum") {
	if (batchnum < 0) {
		BatchError(filetype,-999,19,"BatchNum"); errors++;
	}
	else b.batchNum = batchnum;
}
else controlFormatError = true; // wrong control file format

controlfile >> paramname >> patchmodel;
if (paramname == "PatchModel") {
	if (patchmodel < 0 || patchmodel > 1) {
		BatchError(filetype,-999,1,"PatchModel"); errors++;
	}
	else b.patchmodel = patchmodel;
}
else controlFormatError = true; // wrong control file format

controlfile >> paramname >> resolution;
if (paramname == "Resolution") {
	if (resolution < 1) {
		BatchError(filetype,-999,11,"Resolution"); errors++;
	}
	else b.resolution = resolution;
}
else controlFormatError = true; // wrong control file format

controlfile >> paramname >> landtype;
if (paramname == "LandType") {
#if SEASONAL
	if (landtype != 0) {
		BatchError(filetype,-999,0,"LandType");
		batchlog << "LandType must be 0 for a seasonal model" << endl;
		errors++;
	}
	else b.landtype = landtype;
#else
	if (landtype != 0 && landtype != 2 && landtype != 9) {
		BatchError(filetype,-999,0,"LandType");
		batchlog << "LandType must be 0, 2 or 9" << endl;
		errors++;
	}
	else {
		if (landtype == 9 && patchmodel) {
			BatchError(filetype,-999,0,"LandType");
			batchlog << "LandType may not be 9 for a patch-based model" << endl;
			errors++;
		}
		else b.landtype = landtype;
	}
#endif // SEASONAL 
}
else controlFormatError = true; // wrong control file format

controlfile >> paramname >> maxNhab;
if (paramname == "MaxHabitats") {
	if (landtype == 0) { // raster with unique habitat codes
		if (maxNhab < 2) {
			BatchError(filetype,-999,12,"MaxHabitats"); errors++;
		}
		else b.maxNhab = maxNhab;
	}
	else { // raster with habitat quality OR artificial landscape
		if (maxNhab != 1) {
			BatchError(filetype,-999,0," "); errors++;
			batchlog << "MaxHabitats must be 1 for LandType = " << landtype << endl;
		}
		else {
			if (landtype == 9) // artificial landscape
				// although the user enters 1, the actual number of habitats is 2
				b.maxNhab = 2;
			else
				b.maxNhab = maxNhab;
		}
	}
}
else controlFormatError = true; // wrong control file format

controlfile >> paramname >> speciesdist;
if (paramname == "SpeciesDist") {
	if (speciesdist < 0 || speciesdist > 1) {
		BatchError(filetype,-999,1,"SpeciesDist"); errors++;
	}
	else {
		if (speciesdist != 0 && landtype == 9) {
			BatchError(filetype,-999,0,"SpeciesDist");
			batchlog << "SpeciesDist must be 0 for an artificial landscape" << endl;
			errors++;
			
		}
		else b.speciesdist = speciesdist;
	}
}
else controlFormatError = true; // wrong control file format

controlfile >> paramname >> distresolution;
if (paramname == "DistResolution") {
	if (speciesdist == 1) { // distribution resolution is required
		if (distresolution < resolution) {
			BatchError(filetype,-999,0,"DistResolution");
			batchlog << "DistResolution may not be less than Resolution" << endl;
			errors++;
		}
		else {
			if (distresolution%resolution) {
				BatchError(filetype,-999,0,"DistResolution");
				batchlog << "DistResolution must be an integer multiple of Resolution" << endl;
				errors++;
			}
			else b.distresolution = distresolution;
		}
	}
}
else controlFormatError = true; // wrong control file format

controlfile >> paramname >> reproductn;
sexesDem = sexesDisp = 0;
if (paramname == "Reproduction") {
#if GROUPDISP
	if (reproductn < 0 || reproductn > 3) {
		BatchError(filetype,-999,3,"Reproduction"); errors++;
	}
#else
	if (reproductn < 0 || reproductn > 2) {
		BatchError(filetype,-999,2,"Reproduction"); errors++;
	}
#endif
	else {
		switch (reproductn) {
		case 0: { sexesDem = 1; sexesDisp = 1; break; }
		case 1: { sexesDem = 1; sexesDisp = 2; break; }
		case 2: { sexesDem = 2; sexesDisp = 2; break; }
#if GROUPDISP
		case 3: { sexesDem = 1; sexesDisp = 1; break; }
#endif
		}
		b.reproductn = reproductn; b.sexesDem = sexesDem; b.sexesDisp = sexesDisp;
	}
}
else controlFormatError = true; // wrong control file format

#if SEASONAL
controlfile >> paramname >> nseasons;
if (paramname == "NSeasons") {
	if (nseasons < 2) {
		BatchError(filetype,-999,12,"NSeasons"); errors++;
	}
	else b.nseasons = nseasons;
}
#else
controlfile >> paramname >> repseasons;
if (paramname == "RepSeasons") {
	if (repseasons < 1) {
		BatchError(filetype,-999,11,"RepSeasons"); errors++;
	}
	else b.repseasons = repseasons;
}
#endif
else controlFormatError = true; // wrong control file format

controlfile >> paramname >> stagestruct;
if (paramname == "StageStruct") {
#if SEASONAL
	if (stagestruct != 1) {
		BatchError(filetype,-999,0," "); errors++;
		batchlog << "StageStruct must be 1 for a partial migration model" << endl;
	}
#else
	if (stagestruct < 0 || stagestruct > 1) {
		BatchError(filetype,-999,1,"StageStruct"); errors++;
	}
#endif
	else b.stagestruct = stagestruct;
}
else controlFormatError = true; // wrong control file format

controlfile >> paramname >> stages;
if (paramname == "Stages") {
	if (stagestruct) {
		if (stages < 2 || stages > 10) {
			BatchError(filetype,-999,0," "); errors++;
			batchlog << "Stages must be between 2 and 10" << endl;
		}
		b.stages = stages;
	}
	else { // non-stage-structured model must have 2 stages
		b.stages = stages = 2;
	}
}
else controlFormatError = true; // wrong control file format

controlfile >> paramname >> transfer;
if (paramname == "Transfer") {
#if RS_CONTAIN
	if (transfer < 0 || transfer > 4) {
		BatchError(filetype,-999,4,"Transfer"); errors++;
	}
	else {
		if (stagestruct || transfer < 3) b.transfer = transfer;
		else {
			BatchError(filetype,-999,0," "); errors++;
			batchlog << "Transfer option " << transfer 
				<< " is not allowed for a non-structured species" << endl;
		} 
	}
#else
	if (transfer < 0 || transfer > 2) {
		BatchError(filetype,-999,2,"Transfer"); errors++;
	}
	else b.transfer = transfer;
#endif // RS_CONTAIN 
}
else controlFormatError = true; // wrong control file format

#if RS_ABC

controlfile >> paramname >> nABCsamples;
if (paramname == "ABCsamplesize") {
	if (nABCsamples < 1) {
		BatchError(filetype,-999,10,"ABCsamplesize"); errors++;
	}
}
else controlFormatError = true; // wrong control file format

#endif

if (controlFormatError || errors > 0) { // terminate batch error checking
	if (controlFormatError) {
		CtrlFormatError();
	}
	batchlog << endl
		<< "*** Model parameters in Control file must be corrected before further input file checks are conducted"
		<< endl;
	batchlog.close(); batchlog.clear();
	b.ok = false;
	controlfile.close(); controlfile.clear();
	return b;
}

// Check parameter file
controlfile >> paramname >> filename;
if (paramname == "ParameterFile" && !controlFormatError) {
	fname = indir + filename;
	batchlog << endl << "Checking " << paramname << " " << fname << endl;
	bParamFile.open(fname.c_str());
	if (bParamFile.is_open()) {
#if BUTTERFLYDISP
		b.nSimuls = ParseParameterFile(indir);
#else
		b.nSimuls = ParseParameterFile();
#endif
		if (b.nSimuls < 0) {
			b.ok = false;
		}
		else {
			FileOK(paramname,b.nSimuls,0);
			b.parameterFile = fname;
		}
		bParamFile.close();
	}
	else {
		OpenError(paramname,fname); b.ok = false;
		cout << "Unable to open ParameterFile" << endl;
	}
	bParamFile.clear();
	if (!b.ok) {
		batchlog << endl
			<< "*** ParameterFile must be corrected before further input file checks are conducted"
			<< endl;
		batchlog.close(); batchlog.clear();
		b.ok = false;
		controlfile.close(); controlfile.clear();
		return b;
	}
}
else controlFormatError = true; // wrong control file format
if (bParamFile.is_open()) bParamFile.close();
bParamFile.clear();

// Check land file
controlfile >> paramname >> filename;
if (paramname == "LandFile" && !controlFormatError) {
	fname = indir + filename;
	batchlog << endl << "Checking " << paramname << " " << fname << endl;
	bLandFile.open(fname.c_str());
	if (bLandFile.is_open()) {
		lines = ParseLandFile(landtype,indir);
		if (lines < 0) {
			b.ok = false;
			if (lines < -111)
				batchlog << "*** Format error in " << paramname << endl;
		}
		else {
			FileOK(paramname,lines,1);
			b.landFile = fname; b.nLandscapes = lines;
		}
		bLandFile.close();
	}
	else {
		OpenError(paramname,fname); b.ok = false;
	}
	bLandFile.clear();
}
else controlFormatError = true; // wrong control file format

/*
#if SEASONAL
// Check seasonal file
controlfile >> paramname >> filename;
if (paramname == "SeasonFile" && !controlFormatError) {
	fname = indir + filename;
	batchlog << endl << "Checking " << paramname << " " << fname << endl;
	bSeasonFile.open(fname.c_str());
	if (bSeasonFile.is_open()) {
		lines = ParseSeasonFile(indir);
		if (lines < 0) {
			b.ok = false;
			if (lines < -111)
				batchlog << "*** Format error in " << paramname << endl;
		}
		else {
			if (lines == b.nseasons) {
				FileOK(paramname,lines,3);    
				b.seasonFile = fname;
			}
			else {
				b.ok = false;
				batchlog << "*** No. of seasons in " << filename
					<< " does not match no. in Control file" << endl;
			}
		}
		bSeasonFile.close();
	}
	else {
		OpenError(paramname,fname); b.ok = false;
	}
	bSeasonFile.clear();
}
else controlFormatError = true; // wrong control file format
#endif
*/

// Check stage structure file if required file
controlfile >> paramname >> filename;
batchlog << endl;
if (paramname == "StageStructFile" && !controlFormatError) {
	if (filename == "NULL") {
		if (stagestruct) {
			batchlog << "*** File name is required for " << paramname << endl;
			b.ok = false;
		}
		else b.stageStructFile = filename;
	}
	else { // filename is not NULL
		if (stagestruct) { // check file only if it is required
			fname = indir + filename;
			batchlog << "Checking " << paramname << " " << fname << endl;
			bStageStructFile.open(fname.c_str());
			if (bStageStructFile.is_open()) {
				nSimuls = ParseStageFile(indir);
				if (nSimuls < 0) {
					b.ok = false;
				}
				else {
					FileOK(paramname,nSimuls,0);
					if (nSimuls != b.nSimuls) {
						SimulnCountError(filename); b.ok = false;
					}
					else b.stageStructFile = fname;
				}
				bStageStructFile.close();
			}
			else {
				OpenError(paramname,fname); b.ok = false;
			}
			bStageStructFile.clear();
		} // end of required
		else { // file is not required, and filename should be NULL
			if (filename != "NULL") {
				batchlog << "*** File name for stageStructFile should be NULL as StageStruct = "
					<< stagestruct << endl;
				b.ok = false;
			}
		}
	}
}
else controlFormatError = true; // wrong control file format

// Check emigration file
controlfile >> paramname >> filename;
if (paramname == "EmigrationFile" && !controlFormatError) {
	fname = indir + filename;
	batchlog << endl << "Checking " << paramname << " " << fname << endl;
	bEmigrationFile.open(fname.c_str());
	if (bEmigrationFile.is_open()) {
		nSimuls = ParseEmigFile();
		if (nSimuls < 0) {
			b.ok = false;
		}
		else {
			FileOK(paramname,nSimuls,0);
			if (nSimuls != b.nSimuls) {
				SimulnCountError(filename); b.ok = false;
			}
			else b.emigrationFile = fname;
		}
		bEmigrationFile.close();
	}
	else {
		OpenError(paramname,fname); b.ok = false;
	}
	bEmigrationFile.clear();
}
else controlFormatError = true; // wrong control file format

// Check transfer file
controlfile >> paramname >> filename;
if (paramname == "TransferFile" && !controlFormatError) {
	fname = indir + filename;
	batchlog << endl << "Checking " << paramname << " " << fname << endl;
	bTransferFile.open(fname.c_str());
	if (bTransferFile.is_open()) {
		nSimuls = ParseTransferFile(indir);
		if (nSimuls < 0) {
			b.ok = false;
		}
		else {
			FileOK(paramname,nSimuls,0);
			if (nSimuls != b.nSimuls) {
				SimulnCountError(filename); b.ok = false;
			}
			else b.transferFile = fname;
		}
		bTransferFile.close(); bTransferFile.clear();
	}
	else {
		OpenError(paramname,fname); b.ok = false;
	}
  bTransferFile.clear();
}
else controlFormatError = true; // wrong control file format

// Check settlement file
controlfile >> paramname >> filename;
if (paramname == "SettlementFile" && !controlFormatError) {
	fname = indir + filename;
	batchlog << endl << "Checking " << paramname << " " << fname << endl;
	bSettlementFile.open(fname.c_str());
	if (bSettlementFile.is_open()) {
		nSimuls = ParseSettleFile();
		if (nSimuls < 0) {
			b.ok = false;
		}
		else {
			FileOK(paramname,nSimuls,0);
			if (nSimuls != b.nSimuls) {
				SimulnCountError(filename); b.ok = false;
			}
			else b.settleFile = fname;
		}
		bSettlementFile.close();
	}
	else {
		OpenError(paramname,fname); b.ok = false;
	}
	bSettlementFile.clear();
}
else controlFormatError = true; // wrong control file format

// Check genetics file (optional)
controlfile >> paramname >> filename;
batchlog << endl;
if (paramname == "GeneticsFile" && !controlFormatError) {
	if (filename == "NULL") {
		// this is allowed, because at this stage we do not know whether any simulation
		// includes individual variability - if so, default genetics settings are applied
		b.geneticsFile = filename;
	}
	else { // filename is not NULL
		fname = indir + filename;
		batchlog << "Checking " << paramname << " " << fname << endl;
		bGeneticsFile.open(fname.c_str());
		if (bGeneticsFile.is_open()) {
			nSimuls = ParseGeneticsFile(indir);
			if (nSimuls < 0) {
				b.ok = false;
			}
			else {
				FileOK(paramname,nSimuls,0);
				if (nSimuls != b.nSimuls) {
					SimulnCountError(filename); b.ok = false;
				}
				else b.geneticsFile = fname;
			}
			bGeneticsFile.close();
		}
		else {
			OpenError(paramname,fname); b.ok = false;
		}
		bGeneticsFile.clear();
	}
}
else controlFormatError = true; // wrong control file format

#if RS_CONTAIN

// Check adaptive management file (optional)
controlfile >> paramname >> filename;
batchlog << endl;
if (paramname == "ManagementFile" && !controlFormatError) {
	if (filename == "NULL") {
		// management cull is not applied
		b.manageFile = filename;
	}
	else { // filename is not NULL
		fname = indir + filename;
		batchlog << "Checking " << paramname << " " << fname << endl;
		bManageFile.open(fname.c_str());
		if (bManageFile.is_open()) {
			nSimuls = ParseManageFile(indir);
			if (nSimuls < 0) {
				b.ok = false;
			}
			else {
				FileOK(paramname,nSimuls,0);
				if (nSimuls != b.nSimuls) {
					SimulnCountError(filename); b.ok = false;
				}
				else b.manageFile = fname;
			}
			bManageFile.close();
		}
		else {
			OpenError(paramname,fname); b.ok = false;
		}
		bManageFile.clear();
	}
}
else controlFormatError = true; // wrong control file format

#endif // RS_CONTAIN 

// Check initialisation file
controlfile >> paramname >> filename;
if (paramname == "InitialisationFile" && !controlFormatError) {
	fname = indir + filename;
	batchlog << endl << "Checking " << paramname << " " << fname << endl;
	bInitFile.open(fname.c_str());
	if (bInitFile.is_open()) {
		nSimuls = ParseInitFile(indir);
		if (nSimuls < 0) {
			b.ok = false;
		}
		else {
			FileOK(paramname,nSimuls,0);
			if (nSimuls != b.nSimuls) {
				SimulnCountError(filename); b.ok = false;
			}
			else b.initFile = fname;
		}
		bInitFile.close();
	}
	else {
		OpenError(paramname,fname); b.ok = false;
	}
	bInitFile.clear();
}
else controlFormatError = true; // wrong control file format

#if VIRTUALECOLOGIST
// Check virtual ecologist file (optional)
// CURRENTLY IMPLEMENTED FOR LANDSCAPE GENETICS ONLY...
controlfile >> paramname >> filename;
batchlog << endl;
if (paramname == "VirtualEcolFile" && !controlFormatError) {
	if (filename == "NULL") {
		// this is allowed, because at this stage we do not know whether any simulation
		// includes individual variability - if so, default genetics settings are applied
		b.virtEcolFile = filename;
	}
	else { // filename is not NULL
		fname = indir + filename;
		batchlog << "Checking " << paramname << " " << fname << endl;
		bVirtEcolFile.open(fname.c_str());
		if (bVirtEcolFile.is_open()) {
			nSimuls = ParseVirtEcolFile(indir);
			if (nSimuls < 0) {
				b.ok = false;
			}
			else {
				FileOK(paramname,nSimuls,0);
				if (nSimuls != b.nSimuls) {
					SimulnCountError(filename); b.ok = false;
				}
				else b.virtEcolFile = fname;
			}
			bVirtEcolFile.close();
		}
		else {
			OpenError(paramname,fname); b.ok = false;
		}
		bVirtEcolFile.clear();
	}
}
else controlFormatError = true; // wrong control file format
#endif

#if RS_ABC

// Check ABC parameters file
controlfile >> paramname >> filename;
batchlog << endl;
if (paramname == "ABCParamsFile" && !controlFormatError) {
	fname = indir + filename;
	batchlog << "Checking " << paramname << " " << fname << endl;
	bABCpFile.open(fname.c_str());
	if (bABCpFile.is_open()) {
		int nparams = ParseABCParamsFile();
		if (nparams < 1) {
			b.ok = false;
		}
		else {
			FileOK(paramname,nparams,2);
			b.abcParamsFile = fname;
		}
		bABCpFile.close();
	}
	else {
		OpenError(paramname,fname); b.ok = false;
	}
	bABCpFile.clear();
}
else controlFormatError = true; // wrong control file format

// Check ABC observations file
controlfile >> paramname >> filename;
batchlog << endl;
if (paramname == "ABCObsFile" && !controlFormatError) {
	fname = indir + filename;
	batchlog << "Checking " << paramname << " " << fname << endl;
	bABCoFile.open(fname.c_str());
	if (bABCoFile.is_open()) {
		int nparams = ParseABCObsFile();
		if (nparams < 1) {
			b.ok = false;
		}
		else {
			FileOK(paramname,nparams,2);
			b.abcObsFile = fname;
		}
		bABCoFile.close();
	}
	else {
		OpenError(paramname,fname); b.ok = false;
	}
	bABCoFile.clear();
}
else controlFormatError = true; // wrong control file format

#endif // RS_ABC

if (controlFormatError) {
	CtrlFormatError();
	b.ok = false;
}

if (controlfile.is_open()) 	{ controlfile.close(); controlfile.clear(); }
if (batchlog.is_open()) 		{ batchlog.close(); batchlog.clear(); }

// NOTE: THE FOLLOWING ELEMENTS COULD BE REMOVED FROM b ...
parameterFile = b.parameterFile;
landFile = b.landFile;
#if SEASONAL
//seasonFile = b.seasonFile;
#endif
stageStructFile = b.stageStructFile;
emigrationFile = b.emigrationFile;
transferFile = b.transferFile;
settleFile = b.settleFile;
geneticsFile = b.geneticsFile;
#if RS_CONTAIN
name_managefile = b.manageFile;
#endif // RS_CONTAIN 
initialFile = b.initFile;
#if VIRTUALECOLOGIST
virtEcolFile = b.virtEcolFile;
#endif
#if RS_ABC
abcParamsFile = b.abcParamsFile;
abcObsFile = b.abcObsFile;
#endif

return b;

}

//---------------------------------------------------------------------------
#if BUTTERFLYDISP
int ParseParameterFile(string indir)
#else
int ParseParameterFile(void)
#endif
{
string header,Kheader,intext;
int i,inint,replicates,years;
int absorb,gradient,shifting,shiftstart,shiftend,envstoch,stochtype;
int localext,savemaps;
int prevsimul = 0;
float infloat,minR,maxR,minK,maxK,sum_K,min_K,max_K;
int errors = 0;
int Kerrors = 0;
string filetype = "ParameterFile";

//batchlog << "ParseParametersFile(): starting " << endl;
// Parse header line;
bParamFile >> header; if (header != "Simulation" ) errors++;
bParamFile >> header; if (header != "Replicates" ) errors++;
bParamFile >> header; if (header != "Years" ) errors++;
#if SPATIALMORT
bParamFile >> header; if (header != "MortChgYr" ) errors++;
#endif
bParamFile >> header; if (header != "Absorbing" ) errors++;
bParamFile >> header; if (header != "Gradient" ) errors++;
bParamFile >> header; if (header != "GradSteep" ) errors++;
bParamFile >> header; if (header != "Optimum" ) errors++;
bParamFile >> header; if (header != "f" ) errors++;
bParamFile >> header; if (header != "LocalExtOpt" ) errors++;
bParamFile >> header; if (header != "Shifting" ) errors++;
bParamFile >> header; if (header != "ShiftRate" ) errors++;
bParamFile >> header; if (header != "ShiftStart" ) errors++;
bParamFile >> header; if (header != "ShiftEnd" ) errors++;
bParamFile >> header; if (header != "EnvStoch" ) errors++;
bParamFile >> header; if (header != "EnvStochType" ) errors++;
#if BUTTERFLYDISP
bParamFile >> header; if (header != "EnvStochFromFile" ) errors++;
bParamFile >> header; if (header != "EnvStochFileName" ) errors++;
#endif
bParamFile >> header; if (header != "ac" ) errors++;
bParamFile >> header; if (header != "std" ) errors++;
bParamFile >> header; if (header != "minR" ) errors++;
bParamFile >> header; if (header != "maxR" ) errors++;
bParamFile >> header; if (header != "minK" ) errors++;
bParamFile >> header; if (header != "maxK" ) errors++;
bParamFile >> header; if (header != "LocalExt" ) errors++;
bParamFile >> header; if (header != "LocalExtProb" ) errors++;
#if BUTTERFLYDISP
bParamFile >> header; if (header != "Dispersal" ) errors++;
#endif
bParamFile >> header; if (header != "PropMales" ) errors++;
bParamFile >> header; if (header != "Harem" ) errors++;
bParamFile >> header; if (header != "bc" ) errors++;
bParamFile >> header; if (header != "Rmax" ) errors++;
#if SEASONAL
for (int j = 0; j < nseasons; j++) {
	for (i = 0; i < maxNhab; i++) {
		Kheader = "K" + Int2Str(i+1) + "_" + Int2Str(j);
		bParamFile >> header; if (header != Kheader ) Kerrors++;
	}
}
#else
for (i = 0; i < maxNhab; i++) {
	Kheader = "K" + Int2Str(i+1);
	bParamFile >> header; if (header != Kheader ) Kerrors++;
}
#endif // SEASONAL 
#if GROUPDISP
// ADDITIONAL PARAMETERS FOR GROUP DISPERSAL MODEL
bParamFile >> header; if (header != "Selfing" ) errors++;
bParamFile >> header; if (header != "Paternity" ) errors++;
bParamFile >> header; if (header != "PropLocal" ) errors++;
bParamFile >> header; if (header != "PropNghbr" ) errors++;
#if PEDIGREE
bParamFile >> header; if (header != "PedMatSize" ) errors++;
#endif // PEDIGREE
#endif // GROUPDISP 
#if SOCIALMODEL
// ADDITIONAL PARAMETERS FOR PROBIS SOCIAL POLYMORPHISM MODEL
bParamFile >> header; if (header != "phenMean" ) errors++;
bParamFile >> header; if (header != "phenSD" ) errors++;
bParamFile >> header; if (header != "phenScale" ) errors++;
bParamFile >> header; if (header != "ratioK" ) errors++;
bParamFile >> header; if (header != "ratioRmax" ) errors++;
bParamFile >> header; if (header != "ratioBc" ) errors++;
//bParamFile >> header; if (header != "Ra" ) errors++;
//bParamFile >> header; if (header != "Rs" ) errors++;
bParamFile >> header; if (header != "Ta" ) errors++;
bParamFile >> header; if (header != "Ts" ) errors++;
bParamFile >> header; if (header != "ca" ) errors++;
bParamFile >> header; if (header != "cs" ) errors++;
bParamFile >> header; if (header != "ba" ) errors++;
bParamFile >> header; if (header != "bs" ) errors++;
bParamFile >> header; if (header != "dK" ) errors++;
bParamFile >> header; if (header != "alpha" ) errors++;
#endif // SOCIALMODEL
bParamFile >> header; if (header != "OutStartPop" ) errors++;
bParamFile >> header; if (header != "OutStartInd" ) errors++;
bParamFile >> header; if (header != "OutStartGenetic" ) errors++;
bParamFile >> header; if (header != "OutStartTraitCell" ) errors++;
bParamFile >> header; if (header != "OutStartTraitRow" ) errors++;
bParamFile >> header; if (header != "OutStartConn" ) errors++;
#if RS_CONTAIN
bParamFile >> header; if (header != "OutStartDamage" ) errors++;
#endif // RS_CONTAIN 
bParamFile >> header; if (header != "OutIntRange" ) errors++;
bParamFile >> header; if (header != "OutIntOcc" ) errors++;
bParamFile >> header; if (header != "OutIntPop" ) errors++;
bParamFile >> header; if (header != "OutIntInd" ) errors++;
bParamFile >> header; if (header != "OutIntGenetic" ) errors++;
bParamFile >> header; if (header != "OutGenType" ) errors++;
bParamFile >> header; if (header != "OutGenCrossTab" ) errors++;
bParamFile >> header; if (header != "OutIntTraitCell" ) errors++;
bParamFile >> header; if (header != "OutIntTraitRow" ) errors++;
bParamFile >> header; if (header != "OutIntConn" ) errors++;
#if RS_CONTAIN
bParamFile >> header; if (header != "OutIntDamage" ) errors++;
#endif // RS_CONTAIN 
bParamFile >> header; if (header != "SaveMaps" ) errors++;
bParamFile >> header; if (header != "MapsInterval" ) errors++;
bParamFile >> header; if (header != "SMSHeatMap" ) errors++;
bParamFile >> header; if (header != "DrawLoadedSp" ) errors++;
if (errors > 0 || Kerrors > 0) {
	FormatError(filetype,errors);
	batchlog << "*** Ensure column headers are correct to continue checking data" << endl;
	if (Kerrors > 0) {
		BatchError(filetype,-999,333,"K");
	}
	return -111;
}

// Parse data lines
int line = 1;
int nSimuls = 0;
inint = -98765;
bParamFile >> inint; // first simulation number
if (inint < 0) {
	batchlog << "*** Error in ParameterFile - first simulation number must be >= 0" << endl;
	errors++;
}
else {
	prevsimul = firstsimul = inint; nSimuls++;
}
while (inint != -98765) {
	bParamFile >> replicates; if (replicates <= 0) { BatchError(filetype,line,11,"Replicates"); errors++; }
	bParamFile >> years; if (years <= 0) { BatchError(filetype,line,11,"Years"); errors++; }
#if SPATIALMORT
	bParamFile >> inint;
	if (landtype != 9 && inint <= 0) { BatchError(filetype,line,11,"MortChgYr"); errors++; }
#endif
	bParamFile >> absorb;
	if (absorb < 0 || absorb > 1) { BatchError(filetype,line,1,"Absorbing"); errors++; }
	bParamFile >> gradient;
	if (patchmodel) {
		if (gradient != 0) {
			BatchError(filetype,line,0," ");
			batchlog << "Gradient must be 0 for patch-based model" << endl;
			errors++;
			gradient = 0; // to prevent checking of subsequent fields
		}
	gradient = 0; // to prevent unnecessary checking of subsequent fields
	}
	else { // cell-based model
		if (gradient < 0 || gradient > 3) {
			BatchError(filetype,line,0," ");
			batchlog << "Gradient must be between 0 and 3 for cell-based model" << endl;
			errors++;
		}
	}
	bParamFile >> infloat;
	if (gradient && infloat < 0.0) { BatchError(filetype,line,19,"GradSteep"); errors++; }
	bParamFile >> inint;
	if (gradient && inint < 0) { BatchError(filetype,line,19,"Optimum"); errors++; }
	bParamFile >> infloat;
	if (gradient && infloat < 0.0) { BatchError(filetype,line,19,"f"); errors++; }
	bParamFile >> infloat;
	if (gradient == 4 && (infloat < 0.0 || infloat >= 1.0))
	{ BatchError(filetype,line,20,"LocalExtOpt"); errors++; }
	bParamFile >> shifting;
	if (gradient && (shifting < 0 || shifting > 1)) { BatchError(filetype,line,1,"Shifting"); errors++; }
	bParamFile >> infloat;
	if (gradient && shifting && infloat <= 0.0) { BatchError(filetype,line,10,"ShiftRate"); errors++; }
	bParamFile >> shiftstart;
	if (gradient && shifting && shiftstart <= 0) { BatchError(filetype,line,10,"ShiftStart"); errors++; }
	bParamFile >> shiftend;
	if (gradient && shifting && shiftend <= shiftstart) {
		BatchError(filetype,line,0," ");
		batchlog << "ShiftEnd must be greater than ShiftStart" << endl;
		errors++;
	}
	bParamFile >> envstoch;
	if (patchmodel == 0) { // cell-based model
		if (envstoch < 0 || envstoch > 2) {
			BatchError(filetype,line,0," ");
			batchlog << "EnvStoch must be 0, 1 or 2 for cell-based model" << endl;
			errors++;
//			envstoch = 0; // to prevent checking of subsequent fields
		}
	}
	else { // patch-based model
		if (envstoch < 0 || envstoch > 1) {
			BatchError(filetype,line,0," ");
			batchlog << "EnvStoch must be 0 or 1 for patch-based model" << endl;
			errors++;
//			envstoch = 0; // to prevent checking of subsequent fields
		}
	}
	bParamFile >> stochtype;
	if (envstoch && (stochtype < 0 || stochtype > 1)) {
		BatchError(filetype,line,1,"EnvStochType"); errors++;
	}
#if BUTTERFLYDISP
	ifstream bStochFile;
	int envstochfromfile;
	string envstochfilename;
	string ftype = "EnvStochFileName";
	bParamFile >> envstochfromfile >> envstochfilename;
	if (envstoch == 1) { // global stochasticity
		if (envstochfromfile < 0 || envstochfromfile > 1) {
			BatchError(filetype,line,1,"EnvStochFromFile"); errors++;
		}
	}
	else {
		if (envstochfromfile != 0) {
			BatchError(filetype,line,0," ");
			batchlog << "EnvStochFromFile must be 0 if EnvStoch is not 1" << endl;
			errors++;
		}
	}
	if (envstochfromfile == 1) {
		if (envstochfilename == "NULL") {
			BatchError(filetype,line,0," ");
			batchlog << ftype << " must be supplied if EnvStochFromFile is 1" << endl;
			errors++;
		}
		else {
			string fname = indir + envstochfilename;
			batchlog << "Checking " << ftype << " " << fname << endl;
			bStochFile.open(fname.c_str());
			if (bStochFile.is_open()) {
				// check selected file
				string hdr0,hdr1;
				int year,prevyear;
				int nyears = 0;
				float epsilon;
				bStochFile >> hdr0 >> hdr1;
				if (hdr0 != "Year" || hdr1 != "Epsilon") {
					batchlog << "*** Format error in EnvStochFileName:" << msgcase << msgmatch << endl;
					errors++;
				}
				year = prevyear = -98765;
				bStochFile >> year >> epsilon;
				while (year != -98765) {
					if (nyears == 0) {
						if (year != 0) {
							BatchError(filetype,line,0," ");
							batchlog << "First year in file must be 0" << endl;
							errors++;
						}
					}
					else {
						if (year != prevyear+1) {
							BatchError(filetype,line,0," ");
							batchlog << "Years must be sequential" << endl;
							errors++;
						}
					}
					nyears++;
					prevyear = year;
					year = -98765;
					bStochFile >> year >> epsilon;
				}
				if (nyears < years+1 ) {
					BatchError(filetype,line,0," ");
					batchlog << "Too few years specified in selected file" << endl;
					errors++;
				}
				bStochFile.close();
			}
			else {
				OpenError(ftype,fname); errors++;
			}
			if (bStochFile.is_open()) bStochFile.close();
			bStochFile.clear();
		}
	}
	else {
		if (envstochfilename != "NULL") {
			BatchError(filetype,line,0," ");
			batchlog << ftype << " must be NULL if EnvStochFromFile is 0" << endl;
			errors++;
		}
	}
#endif
	bParamFile >> infloat;
#if BUTTERFLYDISP
	if (envstochfromfile != 1) {
#endif
	if (envstoch && (infloat < 0.0 || infloat >= 1.0)) { BatchError(filetype,line,20,"ac"); errors++; }
#if BUTTERFLYDISP
	}
#endif
	bParamFile >> infloat;
#if BUTTERFLYDISP
	if (envstochfromfile != 1) {
#endif
	if (envstoch && (infloat <= 0.0 || infloat > 1.0)) { BatchError(filetype,line,20,"std"); errors++; }
#if BUTTERFLYDISP
	}
#endif
	bParamFile >> minR;
	if (envstoch && stochtype == 0 && minR <= 0.0) { BatchError(filetype,line,10,"minR"); errors++; }
	bParamFile >> maxR;
	if (envstoch && stochtype == 0 && maxR <= minR) {
		BatchError(filetype,line,0," ");
		batchlog << "maxR must be greater than minR" << endl;
		errors++;
	}
	bParamFile >> minK >> maxK;
	if (envstoch && stochtype == 1) {
		if (minK <= 0.0) { BatchError(filetype,line,10,"minK"); errors++; }
		if (maxK <= minK) {
			BatchError(filetype,line,0," ");
			batchlog << "maxK must be greater than minK" << endl;
			errors++;
		}
	}
	bParamFile >> localext;
	if (patchmodel == 0) { // cell-based model
		if (localext < 0 || localext > 1) {
			BatchError(filetype,line,1,"LocalExt");
			errors++;
		}
		else {
			if (gradient == 4) { // gradient in local extinction probability
				if (localext != 0) {
					BatchError(filetype,line,0," ");
					batchlog << "LocalExt must be zero if Gradient is 4" << endl;
					errors++;
				}
			}
		}
	}
	else { // patch-based model
		if (localext != 0) {
			BatchError(filetype,line,0,"null");
			batchlog << "LocalExt must be 0 for patch-based model" << endl;
			errors++;
		}
	}
	bParamFile >> infloat;
	if (patchmodel == 0 && localext == 1 && (infloat <= 0.0 || infloat >= 1.0))
	{ BatchError(filetype,line,20,"LocalExtProb"); errors++; }
#if BUTTERFLYDISP
	bParamFile >> inint;
	if (stagestruct == 0 && (inint < 0 || inint > 1)) {
		BatchError(filetype,line,1,"Dispersal"); errors++;
	}
#endif
	bParamFile >> infloat;
	if (reproductn && (infloat <= 0.0 || infloat >= 1.0)) {
		BatchError(filetype,line,20,"PropMales"); errors++;
	}
	bParamFile >> infloat;
	if (reproductn == 2 && infloat <= 0.0) { BatchError(filetype,line,10,"Harem"); errors++; }
	bParamFile >> infloat;
	if (stagestruct == 0 && infloat <= 0.0) { BatchError(filetype,line,10,"bc"); errors++; }
	bParamFile >> infloat;
	if (stagestruct == 0 && infloat <= 0.0) { BatchError(filetype,line,10,"Rmax"); errors++; }
#if SEASONAL
	for (int s = 0; s < nseasons; s++) {
		sum_K = 0.0; min_K = 9999999.0; max_K = 0.0;
		for (i = 0; i < maxNhab; i++) {
			bParamFile >> infloat;
			if (infloat < 0.0) {
				Kheader = "K" + Int2Str(i+1) + "_" + Int2Str(s);
				BatchError(filetype,line,19,Kheader); errors++;
			}
			else {
				sum_K += infloat;
				if (infloat > 0.0) {
					if (infloat < min_K) min_K = infloat;
					if (infloat > max_K) max_K = infloat;
				}
			}
		}
		if (sum_K <= 0.0) {
			BatchError(filetype,line,0," "); errors++;
			batchlog << "At least one K column per season must be non-zero" << endl;
		}
		else {
			if (envstoch && stochtype == 1) { // environmental stochasticity in K
				if (min_K < minK || max_K > maxK) {
					BatchError(filetype,line,0," "); errors++;
					batchlog << "Non-zero K values must lie between minK and maxK" << endl;
				}
			}
		}
	}
#else
	sum_K = 0.0; min_K = 9999999.0; max_K = 0.0;
	for (i = 0; i < maxNhab; i++) {
		bParamFile >> infloat;
		if (infloat < 0.0) {
			Kheader = "K" + Int2Str(i+1);
			BatchError(filetype,line,19,Kheader); errors++;
		}
		else {
			sum_K += infloat;
			if (infloat > 0.0) {
				if (infloat < min_K) min_K = infloat;
				if (infloat > max_K) max_K = infloat;
			}
		}
	}
	if (sum_K <= 0.0) {
		BatchError(filetype,line,0," "); errors++;
		batchlog << "At least one K column must be non-zero" << endl;
	}
	else {
		if (envstoch && stochtype == 1) { // environmental stochasticity in K
			if (min_K < minK || max_K > maxK) {
				BatchError(filetype,line,0," "); errors++;
				batchlog << "Non-zero K values must lie between minK and maxK" << endl;
			}
		}
	}
#endif // SEASONAL 

#if GROUPDISP
// ADDITIONAL PARAMETERS FOR GROUP DISPERSAL MODEL
	bParamFile >> inint;
	if (inint < 0 || inint > 1) {
		BatchError(filetype,line,1,"Selfing"); errors++;
	}
	int paternity;
	bParamFile >> paternity;
	if (paternity < 0 || paternity > 2) {
		BatchError(filetype,line,2,"Paternity"); errors++;
	}
	if (paternity == 2 && patchmodel) {
		BatchError(filetype,line,0," "); errors++;
		batchlog << "Paternity may not be 2 for a patch-based model" << endl;
	}
	float proplocal,propnghbr;
	bParamFile >> proplocal >> propnghbr;
	if (paternity == 2) {
		if (proplocal < 0.0 || proplocal > 1.0) {
			BatchError(filetype,line,20,"PropLocal"); errors++;
		}
		if (propnghbr < 0.0 || propnghbr > 1.0) {
			BatchError(filetype,line,20,"PropNghbr"); errors++;
		}
		if (proplocal+propnghbr > 1.0) {
			BatchError(filetype,line,0," "); errors++;
			batchlog << "Local and neighbourhood proportions may not sum to more than 1.0" << endl;
		}
	}
#if PEDIGREE
	bParamFile >> inint; 
	if (inint < 1 ) {
		BatchError(filetype,line,11,"PedMatSize"); errors++;
	}
#endif // PEDIGREE
#endif // GROUPDISP 
#if SOCIALMODEL
	// ADDITIONAL PARAMETERS FOR PROBIS SOCIAL POLYMORPHISM MODEL
	socialParams soc;
//	bParamFile >> soc.asocK >> soc.asocRmax >> soc.asocBc >> soc.ra >> soc.rs
//		>> soc.Ta >> soc.Ts >> soc.dK >> soc.alpha;
//	bParamFile >> soc.asocK >> soc.asocRmax >> soc.asocBc
//		>> soc.Ta >> soc.Ts >> soc.Ca >> soc.Cs >> soc.dK >> soc.alpha;
	bParamFile >> soc.socMean >> soc.socSD >> soc.socScale
		>> soc.asocK >> soc.asocRmax >> soc.asocBc
		>> soc.Ta >> soc.Ts >> soc.ca >> soc.cs >> soc.ba >> soc.bs >> soc.dK >> soc.alpha;
	if (soc.socSD <= 0.0) {
		BatchError(filetype,line,10,"phenSD"); errors++;
	}
	if (soc.socScale <= 0.0) {
		BatchError(filetype,line,10,"phenScale"); errors++;
	}
	if (soc.socSD > soc.socScale) {
		BatchError(filetype,line,3,"phenSD","phenScale"); errors++;
	}
	if (soc.asocK <= 0.0) {
		BatchError(filetype,line,10,"ratioK"); errors++;
	}
	if (soc.asocRmax <= 0.0) {
		BatchError(filetype,line,10,"ratioRmax"); errors++;
	}
	if (soc.asocBc <= 0.0) {
		BatchError(filetype,line,10,"ratioBc"); errors++;
	}
//	if (soc.ra <= 0.0) {
//		BatchError(filetype,line,10,"Ra"); errors++;
//	}
//	if (soc.rs <= 0.0) {
//		BatchError(filetype,line,10,"Rs"); errors++;
//	}
	if (soc.Ta <= 0.0 || soc.Ta >= 1.0) {
		BatchError(filetype,line,20,"Ta"); errors++;
	}
	if (soc.Ts <= 0.0 || soc.Ts >= 1.0) {
		BatchError(filetype,line,20,"Ts"); errors++;
	}
//	if (soc.Ca <= 0.0) {
//		BatchError(filetype,line,10,"Ca"); errors++;
//	}
//	if (soc.Cs <= 0.0) {
//		BatchError(filetype,line,10,"Cs"); errors++;
//	}
	if (soc.dK <= 0.0 || soc.dK >= 1.0) {
		BatchError(filetype,line,20,"dK"); errors++;
	}
	if (soc.alpha <= 0.0) {
		BatchError(filetype,line,10,"alpha"); errors++;
	}
#endif // SOCIALMODEL
	bParamFile >> inint;
	if (inint < 0) { BatchError(filetype,line,19,"OutStartPop"); errors++; }
	bParamFile >> inint;
	if (inint < 0) { BatchError(filetype,line,19,"OutStartInd"); errors++; }
	bParamFile >> inint;
	if (inint < 0) { BatchError(filetype,line,19,"OutStartGenetic"); errors++; }
	bParamFile >> inint;
	if (inint < 0) { BatchError(filetype,line,19,"OutStartTraitCell"); errors++; }
	bParamFile >> inint;
	if (inint < 0) { BatchError(filetype,line,19,"OutStartTraitRow"); errors++; }
	bParamFile >> inint;
	if (inint < 0) { BatchError(filetype,line,19,"OutStartConn"); errors++; }
#if RS_CONTAIN
	bParamFile >> inint;
	if (inint < 0) { BatchError(filetype,line,19,"OutStartDamage"); errors++; }
#endif // RS_CONTAIN 
	bParamFile >> inint;
	if (inint < 0) { BatchError(filetype,line,19,"OutIntRange"); errors++; }
	bParamFile >> inint;
	if (inint < 0) { BatchError(filetype,line,19,"OutIntOcc"); errors++; }
	else {
		if (landtype == 9) {
			if (inint > 0) {
				BatchError(filetype,line,0," "); errors++;
				batchlog << "OutIntOcc must be zero for a generated landscape" << endl;
			}
		}
		else {
			if (replicates < 2 && inint > 0) {
				BatchError(filetype,line,0," "); errors++;
				batchlog << "OutIntOcc may be non-zero only if Replicates >= 2" << endl;
			}
		}
	}
	bParamFile >> inint;
	if (inint < 0) { BatchError(filetype,line,19,"OutIntPop"); errors++; }
	bParamFile >> inint;
	if (inint < 0) { BatchError(filetype,line,19,"OutIntInd"); errors++; }
	bParamFile >> inint;
	if (inint < 0) { BatchError(filetype,line,19,"OutIntGenetic"); errors++; }
	bParamFile >> inint;
	if (inint < 0 || inint > 2) { BatchError(filetype,line,2,"OutGenType"); errors++; }
	bParamFile >> inint;
	if (inint < 0 || inint > 1) { BatchError(filetype,line,1,"OutGenCrossTab"); errors++; }
	bParamFile >> inint;
	if (inint < 0) { BatchError(filetype,line,19,"OutIntTraitCell"); errors++; }
	bParamFile >> inint;
	if (inint < 0) { BatchError(filetype,line,19,"OutIntTraitRow"); errors++; }
	bParamFile >> inint;
	if (inint < 0) { BatchError(filetype,line,19,"OutIntConn"); errors++; }
	else {
		if (patchmodel != 1 && inint > 0) {
			BatchError(filetype,line,0," ");
			batchlog << "OutIntConn may be >0 only if PatchModel is 1" << endl;
			errors++;
		}
	}
#if RS_CONTAIN
	bParamFile >> inint;
	if (landtype == 9) {
		if (inint != 0) {
			BatchError(filetype,line,0," ");
			batchlog << "OutIntDamage must be 0 for an artificial landscape" << endl;
			errors++;					
		}
	}
	else {
		if (inint < 0) { BatchError(filetype,line,19,"OutIntDamage"); errors++; }
	}
#endif // RS_CONTAIN 
	bParamFile >> savemaps; if (savemaps < 0 || savemaps > 1)
	{ BatchError(filetype,line,1,"SaveMaps"); errors++; }
	bParamFile >> inint; if (savemaps == 1 && inint < 1) {
		BatchError(filetype,line,11,"MapsInterval");
		errors++;
	}
	bParamFile >> inint; if (inint < 0 || inint > 1) {
		BatchError(filetype,line,1,"SMSHeatMap");
		errors++;
	}
	bParamFile >> inint; if (savemaps == 1 && (inint < 0 || inint > 1)) {
		BatchError(filetype,line,1,"DrawLoadedSp");
		errors++;
	}
	line++;
	// read next simulation number
	inint = -98765;
	bParamFile >> inint;
	if (bParamFile.eof()) {
		inint = -98765;
	}
	else { // check for valid simulation number
		if (inint != prevsimul+1) {
			BatchError(filetype,line,222," ");
			errors++;
		}
		prevsimul = inint; nSimuls++;
	}
//	batchlog << "ParseParametersFile(): First item of next line = " << inint << endl;
} // end of while loop
if (!bParamFile.eof()) {
	EOFerror(filetype);
	errors++;
}

if (errors > 0) return -111;
else return nSimuls;
}

int ParseLandFile(int landtype, string indir)
{
string fname,header,intext,ftype,costfile;
int j,inint,line;
float infloat;
#if !RS_RCPP
rasterdata patchraster,spdistraster,costraster;
#endif
#if RS_CONTAIN
rasterdata damageraster;
#endif // RS_CONTAIN 
#if SPATIALMORT
rasterdata mortraster;
#endif // SPATIALMORT 
int errors = 0;
//int Kerrors = 0;
int totlines = 0;
//bool errorshown = false;
vector <int> landlist;
string filetype = "LandFile";

//batchlog << "ParseLandFile(): starting " << endl;
if (landtype == 0 || landtype == 2) { // real landscape
	// Parse header line;
	bLandFile >> header; if (header != "LandNum" ) errors++;
//	batchlog << "ParseLandFile(): header = " << header << endl;
	bLandFile >> header; if (header != "Nhabitats" ) errors++;
	bLandFile >> header; if (header != "LandscapeFile" ) errors++;
	bLandFile >> header; if (header != "PatchFile" ) errors++;
	bLandFile >> header; if (header != "CostMapFile" ) errors++;
	bLandFile >> header; if (header != "DynLandFile" ) errors++;
	bLandFile >> header; if (header != "SpDistFile" ) errors++;
#if RS_CONTAIN
	bLandFile >> header; if (header != "DamageFile" ) errors++;
#endif // RS_CONTAIN 
#if SPATIALMORT
	bLandFile >> header; if (header != "MortFile1" ) errors++;
	bLandFile >> header; if (header != "MortFile2" ) errors++;
#endif
	if (errors > 0) {
		FormatError(filetype,0);
		batchlog << "*** Ensure format is correct for real landscape" << endl;
		return -111;
	}
	// Parse data lines
	line = 1;
	inint = -98765;
	bLandFile >> inint;
	while (inint != -98765) {
//		 batchlog << "ParseLandFile(): Landscape no. = " << inint << endl;
		if (inint < 1) {
			BatchError(filetype,line,11,"LandNum"); errors++;
		}
		else {
			// landscape number must be unique - retain in list to check
			for (j = 0; j < (int)landlist.size(); j++) {
				if (inint == landlist[j]) {
					BatchError(filetype,line,666,"LandNum"); j = (int)landlist.size() + 1; errors++;
				}
			}
//		batchlog << "ParseLandFile(): Adding landscape no. " << inint
//			<< " to landscape list" << endl;
			landlist.push_back(inint);
		}
		bLandFile >> inint;
		if (landtype == 0) { // raster map with unique habitat codes
			if (inint < 0) {
				BatchError(filetype,line,10,"Nhabitats"); errors++;
			}
			if (inint > maxNhab) {
				BatchError(filetype,line,0," ");
				batchlog << "Nhabitats may not exceed MaxHabitats in Control file" << endl;
				errors++;
			}
		}

		// check landscape filename
		ftype = "LandscapeFile";
		bLandFile >> intext;
		fname = indir + intext;
		landraster = CheckRasterFile(fname);
		if (landraster.ok) {
			if (landraster.cellsize == resolution)
				batchlog << ftype << " headers OK: " << fname << endl;
			else {
				errors++;
				batchlog << msgresol0 << ftype << " " << fname
						<< msgresol1 << endl;
			}
		}
		else {
			errors++;
			if (landraster.errors == -111)
				OpenError(ftype,fname);
			else
				FormatError(fname,landraster.errors);
		}

		// check patch map filename
		ftype = "PatchFile";
		bLandFile >> intext;
		if (intext == "NULL") {
			if (patchmodel) {
				BatchError(filetype,line,0," "); errors++;
				batchlog << ftype << msgpatch << endl;
			}
		}
		else {
			if (patchmodel) {
				fname = indir + intext;
				patchraster = CheckRasterFile(fname);
				if (patchraster.ok) {
					if (patchraster.cellsize == resolution) {
						if (patchraster.ncols == landraster.ncols
						&&  patchraster.nrows == landraster.nrows
						&&  patchraster.cellsize == landraster.cellsize
						&&  (int)patchraster.xllcorner == (int)landraster.xllcorner
						&&  (int)patchraster.yllcorner == (int)landraster.yllcorner) {
							batchlog << ftype << " headers OK: " << fname << endl;
						}
						else {
							batchlog << msghdrs0 << ftype << " " << fname
								<< msghdrs1 << endl;
							errors++;
						}
					}
					else {
						batchlog << msgresol0 << ftype << " " << fname
							<< msgresol1 << endl;
						errors++;
					}
				}
				else {
					errors++;
					if (patchraster.errors == -111)
						OpenError(ftype,fname);
					else
						FormatError(fname,patchraster.errors);
				}
			}
		}

		// check initial distribution map filename
		ftype = "CostMapFile";
		bLandFile >> costfile;
		if (costfile == "NULL") {
			if (transfer == 1) { // SMS
				if (landtype == 2) {
					BatchError(filetype,line,0," "); errors++;
					batchlog << ftype << " is required for a habitat quality landscape" << endl;
				}
			}
		}
		else {
			if (transfer == 1) { // SMS
				fname = indir + costfile;
				costraster = CheckRasterFile(fname);
				if (costraster.ok) {
					if (costraster.cellsize == resolution) {
						if (costraster.ncols == landraster.ncols
						&&  costraster.nrows == landraster.nrows
						&&  costraster.cellsize == landraster.cellsize
						&&  (int)costraster.xllcorner == (int)landraster.xllcorner
						&&  (int)costraster.yllcorner == (int)landraster.yllcorner) {
							batchlog << ftype << " headers OK: " << fname << endl;
						}
						else {
							batchlog << msghdrs0 << ftype << " " << fname
								<< msghdrs1 << endl;
							errors++;
						}
					}
					else {
						batchlog << msgresol0 << ftype << " " << fname
							<< msgresol1 << endl;
						errors++;
					}
				}
				else {
					errors++;
					if (costraster.errors == -111)
						OpenError(ftype,fname);
					else
						FormatError(fname,costraster.errors);
				}
			}
			else {
				BatchError(filetype,line,0," "); errors++;
				batchlog << ftype << " must be NULL if transfer model is not SMS" << endl;
			}
		}

		// check dynamic landscape filename
		ftype = "DynLandFile";
		bLandFile >> intext;
		if (intext != "NULL") { // landscape is dynamic
			fname = indir + intext;
			batchlog << "Checking " << ftype << " " << fname << endl;
			bDynLandFile.open(fname.c_str());
			if (bDynLandFile.is_open()) {
				int something = ParseDynamicFile(indir,costfile);
				if (something < 0) {
						errors++;
				}
				bDynLandFile.close(); bDynLandFile.clear();
			}
			else {
				bDynLandFile.clear();
				errors++;
				OpenError(ftype,fname);
			}
		}

		// check initial distribution map filename
		ftype = "SpDistFile";
		bLandFile >> intext;
		if (intext == "NULL") {
			if (speciesdist) {
				BatchError(filetype,line,0," "); errors++;
				batchlog << ftype << " is required as SpeciesDist is 1 in Control file" << endl;
			}
		}
		else {
			if (speciesdist) {
				fname = indir + intext;
				spdistraster = CheckRasterFile(fname);
				if (spdistraster.ok) {
					if (spdistraster.cellsize == distresolution) {
						if (spdistraster.cellsize == landraster.cellsize) {
							// check that extent matches landscape extent
							if (spdistraster.ncols != landraster.ncols
							||  spdistraster.nrows != landraster.nrows) {
								batchlog << "*** Extent of " << ftype
									<< " does not match extent of LandscapeFile" << endl;
								errors++;
							}
							else {
								// check origins match
								if ((int)spdistraster.xllcorner == (int)landraster.xllcorner
								&&  (int)spdistraster.yllcorner == (int)landraster.yllcorner) {
									batchlog << ftype << " headers OK: " << fname << endl;
								}
								else {
									batchlog << "*** Origin co-ordinates of " << ftype
										<< " do not match those of LandscapeFile" << endl;
									errors++;
								}
							}
						}
						else { // not able to check extents match
							// check origins match
							if ((int)spdistraster.xllcorner == (int)landraster.xllcorner
							&&  (int)spdistraster.yllcorner == (int)landraster.yllcorner) {
								batchlog << ftype << " headers OK: " << fname << endl;
							}
							else {
								batchlog << "*** Origin co-ordinates of " << ftype
									<< " do not match those of LandscapeFile" << endl;
								errors++;
							}
						}
					}
					else {
						batchlog << "*** Resolution of " << ftype << " " << fname
							<< " does not match DistResolution in Control file" << endl;
						errors++;
					}
				}
				else {
					errors++;
					if (spdistraster.errors == -111)
						OpenError(ftype,fname);
					else
						FormatError(fname,spdistraster.errors);
				}
			}
		}
		
#if RS_CONTAIN
		// check economic / environmental damage map filename
		ftype = "DamageFile";
		bLandFile >> intext;
		if (intext != "NULL") { 
			if (true) {
				fname = indir + intext;
				damageraster = CheckRasterFile(fname);
				if (damageraster.ok) {
					if (damageraster.cellsize == resolution) {
						if (damageraster.ncols == landraster.ncols
						&&  damageraster.nrows == landraster.nrows
						&&  damageraster.cellsize == landraster.cellsize
						&&  (int)damageraster.xllcorner == (int)landraster.xllcorner
						&&  (int)damageraster.yllcorner == (int)landraster.yllcorner) {
							batchlog << ftype << " headers OK: " << fname << endl;
						}
						else {
							batchlog << "*** Headers of " << ftype << " " << fname
								<< " do not match headers of LandscapeFile" << endl;
							errors++;
						}
					}
					else {
						batchlog << msgresol0 << ftype << " " << fname
							<< msgresol1 << endl;
						errors++;
					}
				}
				else {
					errors++;
					if (damageraster.errors == -111)
						OpenError(ftype,fname);
					else
						FormatError(fname,damageraster.errors);
				}
			}
		}
#endif // RS_CONTAIN 
		
#if SPATIALMORT
		// check mortality map filenames
		bool filenull[2];
		for (int mm = 0; mm < 2; mm++) {
			int fnum = mm+1;
			ftype = "MortFile" + Int2Str(fnum);
			bLandFile >> intext;
			if (intext == "NULL") {
				filenull[mm] = true;
			}
			else {
				filenull[mm] = false;
				fname = indir + intext;
				mortraster = CheckRasterFile(fname);
				if (mortraster.ok) {
					if (mortraster.cellsize == resolution) {
						if (mortraster.ncols == landraster.ncols
						&&  mortraster.nrows == landraster.nrows
						&&  mortraster.cellsize == landraster.cellsize
						&&  (int)mortraster.xllcorner == (int)landraster.xllcorner
						&&  (int)mortraster.yllcorner == (int)landraster.yllcorner) {
							batchlog << ftype << " headers OK: " << fname << endl;
						}
						else {
							batchlog << "*** Headers of " << ftype << " " << fname
								<< " do not match headers of LandscapeFile" << endl;
							errors++;
						}
					}
					else {
						batchlog << "*** Resolution of " << ftype << " " << fname
							<< " does not match Resolution in Control file" << endl;
						errors++;
					}
				}
				else {
					errors++;
					if (mortraster.errors == -111)
						OpenError(ftype,fname);
					else
						FormatError(fname,mortraster.errors);
				}
			}
		}
		if ((filenull[0] && !filenull[1]) || (!filenull[0] && filenull[1])) {
			batchlog << "*** There must be either two spatial mortality files,"
				<< " or both must be NULL to omit spatial mortality" << endl;
			errors++;
		}
#endif // SPATIALMORT

		totlines++; line++;
		// read first field on next line
		inint = -98765;
		bLandFile >> inint;
//		batchlog << "ParseLandFile(): first item of next line = " << inint << endl;
	} // end of while loop
	landlist.clear();
} // end of real landscape
else {
	if (landtype == 9) { // artificial landscape
		int fractal,type,Xdim,Ydim;
		float minhab,maxhab;
		// Parse header line;
		bLandFile >> header; if (header != "LandNum" ) errors++;
		bLandFile >> header; if (header != "Fractal" ) errors++;
		bLandFile >> header; if (header != "Type" ) errors++;
		bLandFile >> header; if (header != "Xdim" ) errors++;
		bLandFile >> header; if (header != "Ydim" ) errors++;
		bLandFile >> header; if (header != "MinHab" ) errors++;
		bLandFile >> header; if (header != "MaxHab" ) errors++;
		bLandFile >> header; if (header != "Psuit" ) errors++;
		bLandFile >> header; if (header != "H" ) errors++;
		if (errors > 0) {
			FormatError(filetype,0);
			batchlog << "*** Ensure format is correct for artificial landscape" << endl;
			return -111;
		}
		// Parse data lines
		line = 1;
		inint = -98765;
		bLandFile >> inint;
		while (inint != -98765) {
			for (j = 0; j < (int)landlist.size(); j++) {
				if (inint < 1 || inint == landlist[j]) {
					BatchError(filetype,line,666,"LandNum"); j = (int)landlist.size() + 1; errors++;
				}
			}
//		batchlog << "ParseLandFile(): Adding landscape no. " << inint
//			<< " to landscape list" << endl;
			landlist.push_back(inint);
			bLandFile >> fractal;
			if (fractal < 0 || fractal > 1) {
				BatchError(filetype,line,1,"Fractal"); errors++;
			}
			bLandFile >> type;
			if (type < 0 || type > 1) {
				BatchError(filetype,line,1,"Type"); errors++;
			}
			bLandFile >> Xdim >> Ydim;
			if (fractal == 1) {
				if (Xdim < 3) {
					BatchError(filetype,line,13,"Xdim"); errors++;
				}
				if (Ydim < 3) {
					BatchError(filetype,line,13,"Ydim"); errors++;
				}
			}
			else {
				if (Xdim < 1) {
					BatchError(filetype,line,11,"Xdim"); errors++;
				}
				if (Ydim < 1) {
					BatchError(filetype,line,11,"Ydim"); errors++;
				}
      }
			if (fractal == 1) {
				if (Ydim < Xdim) {
					BatchError(filetype,line,0," ");
					batchlog << "Y dimension may not be less than X dimension" << endl; errors++;
				}
				if ((Xdim > 2 && power2check(Xdim-1) != 1)
				||  (Ydim > 2 && power2check(Ydim-1) != 1)) {
					BatchError(filetype,line,0," ");
					batchlog << "X and Y dimensions must be a power of 2 plus 1" << endl; errors++;
				}
			}
			bLandFile >> minhab >> maxhab;
			if (type == 1) { // continuous landscape
				if (minhab <= 0.0 || minhab >= 100.0) {
					BatchError(filetype,line,100,"MinHab"); errors++;
				}
				if (maxhab <= 0.0 || maxhab > 100.0) {
					BatchError(filetype,line,100,"MaxHab"); errors++;
				}
				if (maxhab <= minhab) {
					BatchError(filetype,line,0," ");
					batchlog << "MaxHab must exceed MinHab" << endl; errors++;
				}
			}
			bLandFile >> infloat;
			if (infloat < 0.0 || infloat > 1.0) {
				BatchError(filetype,line,20,"Psuit"); errors++;
			}
			bLandFile >> infloat;
			if (fractal == 1) {
				if (infloat <= 0.0 || infloat >= 1.0) {
					BatchError(filetype,line,20,"H"); errors++;
				}
			}
			totlines++; line++;
			// read first field on next line
			inint = -98765;
			bLandFile >> inint;
		} // end of while loop
	} // end of artificial landscape
	else { // ERROR condition which should not occur
		batchlog << "*** Critical error in land file. "
			<< "Invalid value of landscape type passed to function ParseLandFile()" << endl;
		errors++;
	}
}
if (!bLandFile.eof()) {
	EOFerror(filetype);
	errors++;
}

if (errors > 0) return -111;
else return totlines;

}

int ParseDynamicFile(string indir,string costfile) {
#if RSDEBUG
DEBUGLOG << "ParseDynamicFile(): costfile=" << costfile << endl;
#endif
string header,filename,fname,ftype,intext;
int change,prevchange,year,prevyear=0;
rasterdata landchgraster,patchchgraster,costchgraster;
int errors = 0;
string filetype = "DynLandFile";
//int totlines = 0;

bDynLandFile >> header; if (header != "Change" ) errors++;
bDynLandFile >> header; if (header != "Year" ) errors++;
bDynLandFile >> header; if (header != "LandChangeFile" ) errors++;
bDynLandFile >> header; if (header != "PatchChangeFile" ) errors++;
bDynLandFile >> header; if (header != "CostChangeFile" ) errors++;

if (errors > 0) {
	FormatError(filetype,errors);
	return -111;
}

// Parse data lines
int line = 1;
change = -98765;
bDynLandFile >> change; // first change number
if (change != 1) {
	batchlog << "*** Error in DynLandFile - first change number must be 1" << endl;
	errors++;
}
else {
	prevchange = change;
}
while (change != -98765) {

	bDynLandFile >> year; if (year <= 0) { BatchError(filetype,line,10,"Year"); errors++; }
	if (line > 1) {
		if (year <= prevyear) {
			BatchError(filetype,line,1,"Year","previous Year"); errors++;
		}
	}
	prevyear = year;

	// check landscape filename
	ftype = "LandChangeFile";
	bDynLandFile >> intext;
//batchlog << "***** indir=" << indir << " intext=" << intext << endl;
	fname = indir + intext;
	landchgraster = CheckRasterFile(fname);
	if (landchgraster.ok) {
		if (landchgraster.cellsize == resolution)
			if (landchgraster.ncols == landraster.ncols
			&&  landchgraster.nrows == landraster.nrows
			&&  landchgraster.cellsize == landraster.cellsize
			&&  (int)landchgraster.xllcorner == (int)landraster.xllcorner
			&&  (int)landchgraster.yllcorner == (int)landraster.yllcorner) {
				batchlog << ftype << " headers OK: " << fname << endl;
			}
			else {
				batchlog << msghdrs0 << ftype << " " << fname
					<< msghdrs1 << endl;
				errors++;
			}
		else {
			errors++;
			batchlog << msgresol0 << ftype << " " << fname << msgresol1 << endl;
		}
	}
	else {
		errors++;
		if (landchgraster.errors == -111)
			OpenError(ftype,fname);
		else
			FormatError(fname,landchgraster.errors);
	}

	// check patch filename
	ftype = "PatchChangeFile";
	bDynLandFile >> intext;
	if (intext == "NULL") {
		if (patchmodel) {
			BatchError(filetype,line,0," "); errors++;
			batchlog << ftype << msgpatch << endl;
		}
	}
	else {
		if (patchmodel) {
			fname = indir + intext;
			patchchgraster = CheckRasterFile(fname);
			if (patchchgraster.ok) {
				if (patchchgraster.cellsize == resolution) {
					if (patchchgraster.ncols == landraster.ncols
					&&  patchchgraster.nrows == landraster.nrows
					&&  patchchgraster.cellsize == landraster.cellsize
					&&  (int)patchchgraster.xllcorner == (int)landraster.xllcorner
					&&  (int)patchchgraster.yllcorner == (int)landraster.yllcorner) {
						batchlog << ftype << " headers OK: " << fname << endl;
					}
					else {
						batchlog << msghdrs0 << ftype << " " << fname
							<< msghdrs1 << endl;
						errors++;
					}
				}
				else {
					batchlog << msgresol0 << ftype << " " << fname
						<< msgresol1 << endl;
					errors++;
				}
			}
			else {
				errors++;
				if (patchchgraster.errors == -111)
					OpenError(ftype,fname);
				else
					FormatError(fname,patchchgraster.errors);
			}
		}
	}

	// check costs filename
	ftype = "CostChangeFile";
	bDynLandFile >> intext;
	if (intext == "NULL") {
		if (costfile != "NULL") {
			BatchError(filetype,line,0," "); errors++;
			batchlog << ftype << " must be supplied " << endl;
		}
	}
	else {
		if (costfile == "NULL") {
			BatchError(filetype,line,0," "); errors++;
			batchlog << ftype << " must be NULL to match LandFile " << endl;
		}
		else {
			fname = indir + intext;
			costchgraster = CheckRasterFile(fname);
			if (costchgraster.ok) {
				if (costchgraster.cellsize == resolution) {
					if (costchgraster.ncols == landraster.ncols
					&&  costchgraster.nrows == landraster.nrows
					&&  costchgraster.cellsize == landraster.cellsize
					&&  (int)costchgraster.xllcorner == (int)landraster.xllcorner
					&&  (int)costchgraster.yllcorner == (int)landraster.yllcorner) {
						batchlog << ftype << " headers OK: " << fname << endl;
					}
					else {
						batchlog << msghdrs0 << ftype << " " << fname
							<< msghdrs1 << endl;
						errors++;
					}
				}
				else {
					batchlog << msgresol0 << ftype << " " << fname
						<< msgresol1 << endl;
					errors++;
				}
			}
			else {
				errors++;
				if (costchgraster.errors == -111)
					OpenError(ftype,fname);
				else
					FormatError(fname,costchgraster.errors);
			}
		}
	}

	line++;
	// read first field on next line
	change = -98765;
	bDynLandFile >> change;
	if (bDynLandFile.eof()) {
		change = -98765;
	}
	else { // check for valid change number
		if (change != prevchange+1) {
			BatchError(filetype,line,0," ");
			batchlog << "Change numbers must be sequential integers" << endl;
			errors++;
		}
		prevchange = change;
	}
}

if (errors > 0) return -111;
else return 0;

}

//---------------------------------------------------------------------------
int ParseStageFile(string indir)
{
string header,filename,fname,ftype2;
int inint,i,err,fecdensdep,fecstagewts,devdensdep,devstagewts,survdensdep,survstagewts;
#if GOBYMODEL
float initPhenSD,initPhenScale,fasocial;
#endif
float infloat;
int errors = 0;
int simuls = 0;
int prevsimul;
bool checkfile;
vector <string> transfiles,wtsfiles;
#if RS_CONTAIN
vector <string> habdemfiles;
bool habdemfileNULL;
#endif // RS_CONTAIN 
#if SEASONAL
vector <string> seasonfiles;
//#if PARTMIGRN
vector <string> extremefiles;
//#endif // PARTMIGRN 
#endif // SEASONAL
string filetype = "StageStructFile";

// Parse header line;
bStageStructFile >> header; if (header != "Simulation" ) errors++;
#if GOBYMODEL
bStageStructFile >> header; if (header != "InitPhenMean" ) errors++;
bStageStructFile >> header; if (header != "InitPhenSD" ) errors++;
bStageStructFile >> header; if (header != "InitPhenScale" ) errors++;
bStageStructFile >> header; if (header != "Fasocial" ) errors++;
#endif
bStageStructFile >> header; if (header != "PostDestructn" ) errors++;
bStageStructFile >> header; if (header != "PRep" ) errors++;
bStageStructFile >> header; if (header != "RepInterval" ) errors++;
bStageStructFile >> header; if (header != "MaxAge" ) errors++;
#if RS_CONTAIN
bStageStructFile >> header; if (header != "HabDemFile" ) errors++;
bStageStructFile >> header; if (header != "TransMatrixFile" ) errors++;
#else
#if SEASONAL
bStageStructFile >> header; if (header != "SeasonFile" ) errors++;
#else
bStageStructFile >> header; if (header != "TransMatrixFile" ) errors++;
#endif // SEASONAL 
#endif // RS_CONTAIN 
bStageStructFile >> header; if (header != "SurvSched" ) errors++;
bStageStructFile >> header; if (header != "FecDensDep" ) errors++;
bStageStructFile >> header; if (header != "FecStageWts" ) errors++;
bStageStructFile >> header; if (header != "FecStageWtsFile" ) errors++;
bStageStructFile >> header; if (header != "DevDensDep" ) errors++;
bStageStructFile >> header; if (header != "DevDensCoeff" ) errors++;
bStageStructFile >> header; if (header != "DevStageWts" ) errors++;
bStageStructFile >> header; if (header != "DevStageWtsFile" ) errors++;
bStageStructFile >> header; if (header != "SurvDensDep" ) errors++;
bStageStructFile >> header; if (header != "SurvDensCoeff" ) errors++;
bStageStructFile >> header; if (header != "SurvStageWts" ) errors++;
bStageStructFile >> header; if (header != "SurvStageWtsFile" ) errors++;
#if PARTMIGRN
bStageStructFile >> header; if (header != "PropPhilRes" ) errors++;
bStageStructFile >> header; if (header != "PropPhilMigFxd" ) errors++;
bStageStructFile >> header; if (header != "PropPhilMigVar" ) errors++;
bStageStructFile >> header; if (header != "PropDispRes" ) errors++;
bStageStructFile >> header; if (header != "PropDispMigFxd" ) errors++;
bStageStructFile >> header; if (header != "PropDispMigVar" ) errors++;
bStageStructFile >> header; if (header != "ResetMigrn" ) errors++;
#endif // PARTMIGRN 
#if SEASONAL
bStageStructFile >> header; if (header != "ExtremeFile" ) errors++;
#endif // SEASONAL 
if (errors > 0) {
	FormatError(filetype,errors);
	return -111;
}

// Parse data lines
int line = 1;
inint = -98765;
bStageStructFile >> inint;
// first simulation number must match first one in parameterFile
if (inint != firstsimul) {
	BatchError(filetype,line,111,"Simulation"); errors++;
}
prevsimul = inint;
while (inint != -98765) {
	simuls++;
#if GOBYMODEL
	bStageStructFile >> infloat; // no check required on initial phenotype mean
	bStageStructFile >> initPhenSD;
	if (initPhenSD <= 0.0) { BatchError(filetype,line,10,"InitPhenSD"); errors++; }
	bStageStructFile >> initPhenScale;
	if (initPhenScale <= 0.0) { BatchError(filetype,line,10,"InitPhenScale"); errors++; }
	if (initPhenSD > initPhenScale) { BatchError(filetype,line,3,"InitPhenSD","InitPhenScale"); errors++; }
	bStageStructFile >> fasocial;
	if (fasocial <= 0.0) { BatchError(filetype,line,10,"Fasocial"); errors++; }
#endif
	bStageStructFile >> inint;
	if (inint < 0 || inint > 1) { BatchError(filetype,line,1,"PostDestructn"); errors++; }
	bStageStructFile >> infloat;
	if (infloat <= 0 || infloat > 1.0) { BatchError(filetype,line,20,"PRep"); errors++; }
	bStageStructFile >> inint;
	if (inint < 0) { BatchError(filetype,line,19,"RepInterval"); errors++; }
	bStageStructFile >> inint;
	if (inint < 2) { BatchError(filetype,line,12,"MaxAge"); errors++; }

	bStageStructFile >> filename;
#if RS_CONTAIN

	// habitat demography file - compulsory if TransMatrixFile is NULL 
	ftype2 = "HabDemFile";    
	if (filename == "NULL") {
		habdemfileNULL = true;
	}
	else {
		habdemfileNULL = false;
		if (landtype != 0) {
			batchlog << "*** HabDemFile may be specified for LandType 0 only" << endl;
			errors++;			
		}
		else {
			checkfile = true;
			for (i = 0; i < (int)habdemfiles.size(); i++) {
				if (filename == habdemfiles[i]) { // file has already been checked
//					batchlog << "*** line = " << line << " i = " << i << " filename = " << filename
//						<< " habdemfiles[i] = " << habdemfiles[i] << endl;
					checkfile = false;
				}
			}
			if (checkfile) {
				fname = indir + filename;
				batchlog << "Checking " << ftype2 << " " << fname << endl;
				bHabDemFile.open(fname.c_str());
				if (bHabDemFile.is_open()) {
					err = ParseHabDemFile(stages,sexesDem,indir);
					if (err == 0) FileHeadersOK(ftype2); else errors++;
					bHabDemFile.close();
				}
				else {
					OpenError(ftype2,fname); errors++;
				}
				if (bHabDemFile.is_open()) bHabDemFile.close();
				bHabDemFile.clear();
			}
			habdemfiles.push_back(filename);
		}
	}

	// transition matrix file - compulsory if HabDemFile is NULL 
	bStageStructFile >> filename;
	ftype2 = "TransMatrixFile";    
	if (filename == "NULL") {
		if (habdemfileNULL) {
			batchlog << "*** " << ftype2 << " is compulsory if HabDemFile is NULL" << endl;
			errors++;
		}
	}
	else {
		if (!habdemfileNULL) {
			batchlog << "*** Only one of HabDemFile and " << ftype2 << " may be specified" << endl;
			errors++;			
		}
		else {
			checkfile = true;
			for (i = 0; i < (int)transfiles.size(); i++) {
				if (filename == transfiles[i]) { // file has already been checked
//					batchlog << "*** line = " << line << " i = " << i << " filename = " << filename
//						<< " transfiles[i] = " << transfiles[i] << endl;
					checkfile = false;
				}
			}
			if (checkfile) {
				fname = indir + filename;
				batchlog << "Checking " << ftype2 << " " << fname << endl;
				bTransMatrix.open(fname.c_str());
				if (bTransMatrix.is_open()) {
					err = ParseTransitionFile(stages,sexesDem);
					if (err == 0) FileHeadersOK(ftype2); else errors++;
					bTransMatrix.close();
				}
				else {
					OpenError(ftype2,fname); errors++;
				}
				if (bTransMatrix.is_open()) bTransMatrix.close();
				bTransMatrix.clear();
			}
			transfiles.push_back(filename);
		}
	}

#else
#if SEASONAL
	// seasons file - compulsory
	ftype2 = "SeasonFile";
	checkfile = true;
	for (i = 0; i < (int)seasonfiles.size(); i++) {
		if (filename == seasonfiles[i]) { // file has already been checked
//			batchlog << "*** line = " << line << " i = " << i << " filename = " << filename
//				<< " seasonfiles[i] = " << seasonfiles[i] << endl;
			checkfile = false;
		}
	}
	if (checkfile) {
		if (filename == "NULL") {
			batchlog << "*** " << ftype2 << " is compulsory for seasonal model" << endl;
			errors++;
		}
		else {
			fname = indir + filename;
			batchlog << "Checking " << ftype2 << " " << fname << endl;
			bSeasonFile.open(fname.c_str());
			if (bSeasonFile.is_open()) {
				int lines = ParseSeasonFile(indir);
				if (lines < 0) {
					errors++;
					if (lines < -111)
						batchlog << "*** Format error in " << ftype2 << endl;
				}
				else {
					if (lines == nseasons) {
						FileOK(ftype2,lines,3);    
//						b.seasonFile = fname;
					}
					else {
						errors++;
						batchlog << "*** No. of seasons in " << filename
							<< " does not match no. in Control file" << endl;
					}
				}
				bSeasonFile.close();
			}
			else {
				OpenError(ftype2,fname); errors++;
			}
			bSeasonFile.clear();
		}
	}
	seasonfiles.push_back(filename);
#else
	// transition matrix file - compulsory
	ftype2 = "TransMatrixFile";
	checkfile = true;
	for (i = 0; i < (int)transfiles.size(); i++) {
		if (filename == transfiles[i]) { // file has already been checked
//			batchlog << "*** line = " << line << " i = " << i << " filename = " << filename
//				<< " transfiles[i] = " << transfiles[i] << endl;
			checkfile = false;
		}
	}
	if (checkfile) {
		if (filename == "NULL") {
			batchlog << "*** " << ftype2 << " is compulsory for stage-structured model" << endl;
			errors++;
		}
		else {
			fname = indir + filename;
			batchlog << "Checking " << ftype2 << " " << fname << endl;
			bTransMatrix.open(fname.c_str());
			if (bTransMatrix.is_open()) {
				err = ParseTransitionFile(stages,sexesDem);
				if (err == 0) FileHeadersOK(ftype2); else errors++;
				bTransMatrix.close();
			}
			else {
				OpenError(ftype2,fname); errors++;
			}
			if (bTransMatrix.is_open()) bTransMatrix.close();
			bTransMatrix.clear();
		}
	}
	transfiles.push_back(filename);
#endif // SEASONAL 
#endif // RS_CONTAIN 

	bStageStructFile >> inint;
#if SEASONAL
	if (inint < 0 || inint > 1) { BatchError(filetype,line,1,"SurvSched"); errors++; }
#else
	if (inint < 0 || inint > 2) { BatchError(filetype,line,2,"SurvSched"); errors++; }
#endif
	bStageStructFile >> fecdensdep;
	if (fecdensdep < 0 || fecdensdep > 1)
	{ BatchError(filetype,line,1,"FecDensDep"); errors++; fecdensdep = 1; }
	bStageStructFile >> fecstagewts;
	if (fecdensdep) {
		if (fecstagewts < 0 || fecstagewts > 1)
		{ BatchError(filetype,line,1,"FecStageWts"); errors++; fecstagewts = 1; }
	}
	else {
		if (fecstagewts != 0) {
			BatchError(filetype,line,0," ");
			batchlog << "FecStageWts must be 0 if FecDensDep is 0" << endl; errors++;
			errors++; fecstagewts = 1;
		}
	}

	// fecundity stage weights file - optional
	ftype2 = "FecStageWtsFile";
	bStageStructFile >> filename;
	if (filename == "NULL") {
		if (fecstagewts) {
			BatchError(filetype,line,0," ");
			batchlog << ftype2 << " is compulsory unless FecStageWts is 0" << endl;
			errors++;
		}
	}
	else {
		checkfile = true;
		for (i = 0; i < (int)wtsfiles.size(); i++) {
			if (filename == wtsfiles[i]) checkfile = false; // file has already been checked
		}
		if (checkfile) {
			fname = indir + filename;
			batchlog << "Checking " << ftype2 << " " << fname << endl;
			bStageWeightsFile.open(fname.c_str());
			if (bStageWeightsFile.is_open()) {
				err = ParseWeightsFile(ftype2);
				if (err == 0) FileHeadersOK(ftype2); else errors++;
				bStageWeightsFile.close();
			}
			else {
				OpenError(ftype2,fname); errors++;
			}
			if (bStageWeightsFile.is_open()) bStageWeightsFile.close();
			bStageWeightsFile.clear();
		}
		wtsfiles.push_back(filename);
	}

	bStageStructFile >> devdensdep;
	if (devdensdep < 0 || devdensdep > 1)
	{ BatchError(filetype,line,1,"DevDensDep"); errors++; devdensdep = 1; }
	bStageStructFile >> infloat >> devstagewts;
	if (devdensdep) {
		if (infloat <= 0.0) {
			BatchError(filetype,line,10,"DevDensCoeff"); errors++;
		}
		if (devstagewts < 0 || devstagewts > 1) {
			BatchError(filetype,line,1,"DevStageWts"); errors++; devstagewts = 1;
		}
	}
	else {
		if (devstagewts != 0) {
			BatchError(filetype,line,0," ");
			batchlog << "DevStageWts must be 0 if DevDensDep is 0" << endl; errors++;
			errors++; devstagewts = 1;
		}
	}

	// development stage weights file - optional
	ftype2 = "DevStageWtsFile";
	bStageStructFile >> filename;
	if (filename == "NULL") {
		if (devstagewts) {
			BatchError(filetype,line,0," ");
			batchlog << ftype2 << " is compulsory unless DevStageWts is 0" << endl;
			errors++;
		}
	}
	else {
		checkfile = true;
		for (i = 0; i < (int)wtsfiles.size(); i++) {
			if (filename == wtsfiles[i]) checkfile = false; // file has already been checked
		}
		if (checkfile) {
			fname = indir + filename;
			batchlog << "Checking " << ftype2 << " " << fname << endl;
			bStageWeightsFile.open(fname.c_str());
			if (bStageWeightsFile.is_open()) {
				err = ParseWeightsFile(ftype2);
				if (err == 0) FileHeadersOK(ftype2); else errors++;
				bStageWeightsFile.close();
			}
			else {
				OpenError(ftype2,fname); errors++;
			}
			if (bStageWeightsFile.is_open()) bStageWeightsFile.close();
			bStageWeightsFile.clear();
		}
		wtsfiles.push_back(filename);
	}

	bStageStructFile >> survdensdep;
	if (survdensdep < 0 || survdensdep > 1)
	{ BatchError(filetype,line,1,"SurvDensDep"); errors++; survdensdep = 1; }
	bStageStructFile >> infloat >> survstagewts;
	if (survdensdep) {
		if (infloat <= 0.0) {
			BatchError(filetype,line,10,"SurvDensCoeff"); errors++;
		}
		if (survstagewts < 0 || survstagewts > 1) {
			BatchError(filetype,line,1,"SurvStageWts"); errors++; survstagewts = 1;
		}
	}
	else {
		if (survstagewts != 0) {
			BatchError(filetype,line,0," ");
			batchlog << "SurvStageWts must be 0 if SurvDensDep is 0" << endl; errors++;
			errors++; survstagewts = 1;
		}
	}

	// survival stage weights file - optional
	ftype2 = "SurvStageWtsFile";
	bStageStructFile >> filename;
	if (filename == "NULL") {
		if (survstagewts) {
			BatchError(filetype,line,0," ");
			batchlog << ftype2 << " is compulsory unless SurvStageWts is 0" << endl;
			errors++;
		}
	}
	else {
		checkfile = true;
		for (i = 0; i < (int)wtsfiles.size(); i++) {
			if (filename == wtsfiles[i]) checkfile = false; // file has already been checked
		}
		if (checkfile) {
			fname = indir + filename;
			batchlog << "Checking " << ftype2 << " " << fname << endl;
			bStageWeightsFile.open(fname.c_str());
			if (bStageWeightsFile.is_open()) {
				err = ParseWeightsFile(ftype2);
				if (err == 0) FileHeadersOK(ftype2); else errors++;
				bStageWeightsFile.close();
			}
			else {
				OpenError(ftype2,fname); errors++;
			}
			if (bStageWeightsFile.is_open()) bStageWeightsFile.close();
			bStageWeightsFile.clear();
		}
		wtsfiles.push_back(filename);
	}

#if PARTMIGRN
	float prop;
	float cumprop = 0.0;
	for (i = 0; i < 6; i++) {
		bStageStructFile >> prop;
		cumprop += prop;
		if (prop < 0.0 || prop > 1.0) {
			BatchError(filetype,line,20,"Proportion Phil/Disp Res/Mig"); errors++;
		}
	}
	if (cumprop < 1.0 || cumprop > 1.0) {
		BatchError(filetype,line,0," "); errors++;
		batchlog << "Proportions must sum to 1.0" << endl;
	}
	int resetmigrn;
	bStageStructFile >> resetmigrn;
	if (resetmigrn < 0 || resetmigrn > 1) {
		BatchError(filetype,line,1,"ResetMigrn"); errors++;
	}
#endif // PARTMIGRN 

#if SEASONAL
	// extreme events file - optional
	ftype2 = "ExtremeFile";
	bStageStructFile >> filename;
	if (filename != "NULL") {
		checkfile = true;
		for (i = 0; i < (int)extremefiles.size(); i++) {
			if (filename == extremefiles[i]) checkfile = false; // file has already been checked
		}
		if (checkfile) {
			fname = indir + filename;
			batchlog << "Checking " << ftype2 << " " << fname << endl;
			bExtremeFile.open(fname.c_str());
			if (bExtremeFile.is_open()) {
				err = ParseExtremeFile(ftype2);
				if (err >= 0) FileHeadersOK(ftype2); else errors++;
				bExtremeFile.close();
			}
			else {
				OpenError(ftype2,fname); errors++;
			}
			if (bExtremeFile.is_open()) bExtremeFile.close();
			bExtremeFile.clear();
		}
		extremefiles.push_back(filename);
	}
#endif // SEASONAL 
	
	// read next simulation
	line++;
	inint = -98765;
	bStageStructFile >> inint;
	if (bStageStructFile.eof()) {
		inint = -98765;
	}
	else { // check for valid simulation number
		if (inint != prevsimul+1) {
			BatchError(filetype,line,222," ");
			errors++;
		}
		prevsimul = inint;
	}
}
if (!bStageStructFile.eof()) {
	EOFerror(filetype);
	errors++;
}

transfiles.clear();
wtsfiles.clear();
#if RS_CONTAIN
habdemfiles.clear();
#endif // RS_CONTAIN 
#if SEASONAL
seasonfiles.clear();
//#if PARTMIGRN
extremefiles.clear();
//#endif // PARTMIGRN 
#endif // SEASONAL

if (errors > 0) return -111;
else return simuls;

}

//---------------------------------------------------------------------------
// Check transition matrix file
int ParseTransitionFile(short nstages,short nsexesDem)
{
string header,hhh;
int i,j,stage,sex,line,minage;
//int prevminage;
float infloat;
int errors = 0;
string filetype = "TransMatrixFile";

// check header records
bTransMatrix >> header; if (header != "Transition" ) errors++;
for (i = 0; i < nstages; i++) {
	for (j = 0; j < nsexesDem; j++) {
		bTransMatrix >> header;
		if (nsexesDem == 1) hhh = Int2Str(i);
		else {
			if (j == 0) hhh = Int2Str(i) + "m"; else hhh = Int2Str(i) + "f";
		}
		if (header != hhh) errors++;
//		batchlog << "i = " << i << " j = " << j << " hhh = " << hhh << " header = " << header
//			<< " errors = " << errors << endl;
	}
}
bTransMatrix >> header; if (header != "MinAge" ) errors++;

if (errors > 0) {
	FormatError(filetype,errors);
	return -111;
}

// check matrix, including row headers

// single row for juveniles
line = 1;
bTransMatrix >> header;
if (header != "0" ) {
	BatchError(filetype,line,0," ");
	batchlog << "Invalid row header" << endl; errors++;
}
float totfecundity = 0.0;
for (i = 0; i < nstages; i++) {
	for (j = 0; j < nsexesDem; j++) {
		bTransMatrix >> infloat;
		if (i > 0) {
			if (infloat < 0.0) {
				BatchError(filetype,line,19,"Fecundity"); errors++;
			}
			totfecundity += infloat;
		}
	}
}
#if !SEASONAL
if (totfecundity <= 0.0) {
	BatchError(filetype,line,10,"Total fecundity"); errors++;
}
#endif
bTransMatrix >> minage;
//prevminage = minage;
//				batchlog << "MINAGE = " << minage << endl;
if (minage != 0) {
	BatchError(filetype,line,0," ");
	batchlog << "MinAge must be zero for juvenile stage" << endl; errors++;
}

// one row for each stage/sex combination
//				batchlog << "HEADER = " << header << endl;
for (stage = 1; stage < nstages; stage++) {
	for (sex = 0; sex < nsexesDem; sex++) {
		line++;
		// row header
		bTransMatrix >> header;
		if (nsexesDem == 1) hhh = Int2Str(stage);
		else {
			if (sex == 0) hhh = Int2Str(stage) + "m"; else hhh = Int2Str(stage) + "f";
		}
		if (header != hhh) {
			BatchError(filetype,line,0," ");
			batchlog << "Invalid row header" << endl; errors++;
		}
		for (i = 0; i < nstages; i++) {
			for (j = 0; j < nsexesDem; j++) {
				bTransMatrix >> infloat;
//				batchlog << "TRANS PROB = " << infloat << endl;
				if (infloat < 0.0 || infloat > 1) {
					BatchError(filetype,line,20,"Transition probability"); errors++;
				}
			}
		}
//		 prevminage = minage;
		 bTransMatrix >> minage;
//				batchlog << "MINAGE = " << minage << endl;
		if (stage == 1 && minage != 0) {
			BatchError(filetype,line,0," ");
			batchlog << "MinAge must be zero for stage 1" << endl; errors++;
		}
		if (stage > 1) {
			if (minage < 0) {
				BatchError(filetype,line,19,"MinAge"); errors++;
			}
			// SCFP 30/9/13 - IDEALLY OUGHT TO TEST THAT MINAGE IS NO LESS THAN PREVIOUS MINAGE
			// BUT WOULD NEED TO BE PREV MINAGE FOR THE SAME SEX
			// HOWEVER, IT IS NOT CRITICAL, AS A MINAGE OF LESS THAN PREVIOUS CANNOT CAUSE ANY
			// PROBLEM, AS PREVIOUS MINAGE GETS APPLIED EARLIER
//			if (minage < prevminage) {
//				BatchError(filetype,line,0," ");
//				batchlog << "MinAge may not be less than MinAge of previous stage" << endl; errors++;
//			}
		}
	}
}
// final read should hit EOF
bTransMatrix >> header;

if (!bTransMatrix.eof()) {
	EOFerror(filetype);
	errors++;
}

return errors;

}

//---------------------------------------------------------------------------

#if RS_CONTAIN

int ParseHabDemFile(short nstages,short nsexesDem,string indir)
{
string header,filename,fname,ftype2;
int inint,i,err;
int errors = 0;
bool checkfile;
vector <string> transfiles;
string filetype = "HabDemFile";

// Parse header line;
bHabDemFile >> header; if (header != "HabCode" ) errors++;
#if SEASONAL
bHabDemFile >> header; if (header != "Season" ) errors++;
#endif // SEASONAL 
bHabDemFile >> header; if (header != "TransMatrixFile" ) errors++;
if (errors > 0) {
	FormatError(filetype,errors);
	return -111;
}

// Parse data lines
int maxcode = NHABITATS;
if (maxcode > maxNhab) maxcode = maxNhab;
int line = 1;
inint = -98765;
bHabDemFile >> inint;
//prevseason = inint;
while (inint != -98765) {

	// check for valid habitat code
	if (inint < 1 || inint > maxcode) {
		BatchError(filetype,line,0," ");
		batchlog << "Habitat code must be from 1 to " << maxcode << endl; errors++;
	}
#if SEASONAL
	bHabDemFile >> inint;
//	if (inint < 0) {
//		BatchError(filetype,line,19,"Season");
//	}
	if (inint < 0 || inint >= nseasons) {
		BatchError(filetype,line,0," ");
		batchlog << "Season must be from 0 to " << nseasons-1 << endl; errors++;
	}
#endif // SEASONAL 
	// transition matrix file - compulsory
	ftype2 = "TransMatrixFile";
	bHabDemFile >> filename;
	checkfile = true;
	for (i = 0; i < (int)transfiles.size(); i++) {
		if (filename == transfiles[i]) { // file has already been checked
//			batchlog << "*** line = " << line << " i = " << i << " filename = " << filename
//				<< " transfiles[i] = " << transfiles[i] << endl;
			checkfile = false;
		}
	}
	if (checkfile) {
		if (filename == "NULL") {
			batchlog << "*** " << ftype2 << " is compulsory for partial migration model" << endl;
			errors++;
		}
		else {
			fname = indir + filename;
			batchlog << "Checking " << ftype2 << " " << fname << endl;
			bTransMatrix.open(fname.c_str());
			if (bTransMatrix.is_open()) {
				err = ParseTransitionFile(stages,sexesDem);
				if (err == 0) FileHeadersOK(ftype2); else errors++;
				bTransMatrix.close();
			}
			else {
				OpenError(ftype2,fname); errors++;
			}
			if (bTransMatrix.is_open()) bTransMatrix.close();
			bTransMatrix.clear();
		}
	}
	transfiles.push_back(filename);
	
	// read next habitat
	line++;
	inint = -98765;
	bHabDemFile >> inint;
	if (bHabDemFile.eof()) {
		inint = -98765;
	}
//	else { // check for valid season number
//		if (inint != prevseason+1) {
//			BatchError(filetype,line,223," ");   
//			errors++;
//		}
//		prevseason = inint;
//	}
}
if (!bHabDemFile.eof()) {
	EOFerror(filetype);
	errors++;
}

transfiles.clear();

if (errors > 0) return -111;
else return 0;

}

#endif // RS_CONTAIN 

//---------------------------------------------------------------------------
// Check stage weights matrix file
int ParseWeightsFile(string filetype)
{
string header,hhh;
int i,j,stage,sex,line;
float infloat;
int errors = 0;

// check header records
bStageWeightsFile >> header; if (header != "StageWts" ) errors++;
for (i = 0; i < stages; i++) {
	for (j = 0; j < sexesDem; j++) {
		bStageWeightsFile >> header;
		if (sexesDem == 1) hhh = Int2Str(i);
		else {
			if (j == 0) hhh = Int2Str(i) + "m"; else hhh = Int2Str(i) + "f";
		}
		if (header != hhh) errors++;
	}
}

if (errors > 0) {
	FormatError(filetype,errors);
	return -111;
}

// check matrix, including row headers
// one row for each stage/sex combination
line = 0;
for (stage = 0; stage < stages; stage++) {
	for (sex = 0; sex < sexesDem; sex++) {
		line++;
		// row header
		bStageWeightsFile >> header;
		if (sexesDem == 1) hhh = Int2Str(stage);
		else {
			if (sex == 0) hhh = Int2Str(stage) + "m"; else hhh = Int2Str(stage) + "f";
		}
		if (header != hhh) {
			BatchError(filetype,line,0," ");
			batchlog << "Invalid row header" << endl; errors++;
		}
		for (i = 0; i < stages; i++) {
			for (j = 0; j < sexesDem; j++) {
				bStageWeightsFile >> infloat;
				// NOTE - any real number is acceptable - no check required
			}
		}
	}
}
// final read should hit EOF
bStageWeightsFile >> header;

if (!bStageWeightsFile.eof()) {
	EOFerror(filetype);
	errors++;
}

return errors;

}

#if SEASONAL

//---------------------------------------------------------------------------
int ParseSeasonFile(string indir)
{
string header,filename,fname,ftype2;
//int inint,i,err,fecdensdep,fecstagewts,devdensdep,devstagewts,survdensdep,survstagewts;
int inint,i,err;
float infloat;
int errors = 0;
int seasons = 0;
int prevseason;
bool checkfile;
vector <string> transfiles;
string filetype = "SeasonFile";

// Parse header line;
bSeasonFile >> header; if (header != "Season" ) errors++;
bSeasonFile >> header; if (header != "TransMatrixFile" ) errors++;
bSeasonFile >> header; if (header != "ProbExtreme" ) errors++;
bSeasonFile >> header; if (header != "MortExtreme" ) errors++;
if (errors > 0) {
	FormatError(filetype,errors);
	return -111;
}

// Parse data lines
int line = 1;
inint = -98765;
bSeasonFile >> inint;
// first season number must be zero
if (inint != 0) {
	BatchError(filetype,line,0," ");
	batchlog << "First season must be 0" << endl; errors++;
}
prevseason = inint;
while (inint != -98765) {
	seasons++;

	// transition matrix file - compulsory
	ftype2 = "TransMatrixFile";
	bSeasonFile >> filename;
	checkfile = true;
	for (i = 0; i < (int)transfiles.size(); i++) {
		if (filename == transfiles[i]) { // file has already been checked
//			batchlog << "*** line = " << line << " i = " << i << " filename = " << filename
//				<< " transfiles[i] = " << transfiles[i] << endl;
			checkfile = false;
		}
	}
	if (checkfile) {
		if (filename == "NULL") {
			batchlog << "*** " << ftype2 << " is compulsory for partial migration model" << endl;
			errors++;
		}
		else {
			fname = indir + filename;
			batchlog << "Checking " << ftype2 << " " << fname << endl;
			bTransMatrix.open(fname.c_str());
			if (bTransMatrix.is_open()) {
				err = ParseTransitionFile(stages,sexesDem);
				if (err == 0) FileHeadersOK(ftype2); else errors++;
				bTransMatrix.close();
			}
			else {
				OpenError(ftype2,fname); errors++;
			}
			if (bTransMatrix.is_open()) bTransMatrix.close();
			bTransMatrix.clear();
		}
	}
	transfiles.push_back(filename);

	bSeasonFile >> infloat;
	if (infloat < 0.0 || infloat > 1.0) { BatchError(filetype,line,20,"ProbExtreme"); errors++; }
	bSeasonFile >> infloat;
	if (infloat < 0.0 || infloat > 1.0) { BatchError(filetype,line,20,"MortExtreme"); errors++; }

	// ADD SURVIVAL SCHEDULE AND/OR STAGE WEIGHTS HERE?????
	
	// read next season
	line++;
	inint = -98765;
	bSeasonFile >> inint;
	if (bSeasonFile.eof()) {
		inint = -98765;
	}
	else { // check for valid season number
		if (inint != prevseason+1) {
			BatchError(filetype,line,223," ");   
			errors++;
		}
		prevseason = inint;
	}
}
if (!bSeasonFile.eof()) {
	EOFerror(filetype);
	errors++;
}

transfiles.clear();

if (errors > 0) return -111;
else return seasons;

}

//---------------------------------------------------------------------------
int ParseExtremeFile(string indir)
{
string header,filename,fname,ftype2;
//int inint,i,err;
int inint,year,prevyear,season,prevseason,x,y;
float infloat;
int errors = 0;
int events = 0;
string filetype = "ExtremeFile";

// Parse header line;
bExtremeFile >> header; if (header != "Year" ) errors++;
bExtremeFile >> header; if (header != "Season" ) errors++;
if (patchmodel) {
	bExtremeFile >> header; if (header != "PatchID" ) errors++;
}
else {
	bExtremeFile >> header; if (header != "X" ) errors++;
	bExtremeFile >> header; if (header != "Y" ) errors++;
}
bExtremeFile >> header; if (header != "MortExtreme" ) errors++;
if (errors > 0) {
	FormatError(filetype,errors);
	return -111;
}

// Parse data lines
int line = 1;
year = prevyear = prevseason = -98765;
bExtremeFile >> year;
while (year != -98765) {
	events++;
	if (year < 0) {
		BatchError(filetype,line,19,"Year"); errors++;
	}
	else {
		if (year < prevyear) {
			BatchError(filetype,line,2,"Year","previous Year"); errors++;
		}
	}
	bExtremeFile >> season;
	if (season < 0) { BatchError(filetype,line,19,"Season"); errors++; }
	if (season >= nseasons) { BatchError(filetype,line,4,"Season","NSeasons in Control file"); errors++; }
	if (year == prevyear) {
		if (season < prevseason) {
			BatchError(filetype,line,2,"Season","previous Season"); errors++;
		}		
	}
	prevyear = year; prevseason = season;
	if (patchmodel) {
		bExtremeFile >> inint;
		if (inint < 1) { BatchError(filetype,line,11,"PatchID"); errors++; }
	}
	else {  
		bExtremeFile >> x >> y;
		if (x < 0 || y < 0) {
			BatchError(filetype,line,19,"X and Y"); errors++;
		}
	}
	bExtremeFile >> infloat;
	if (infloat < 0.0 || infloat > 1.0) { BatchError(filetype,line,20,"MortExtreme"); errors++; }

	// read next event
	line++;
	year = -98765;
	bExtremeFile >> year;
	if (bExtremeFile.eof()) {
		year = -98765;
	}
}
if (!bExtremeFile.eof()) {
	EOFerror(filetype);
	errors++;
}

if (errors > 0) return -111;
else return events;

}

#endif // SEASONAL 

//---------------------------------------------------------------------------
int ParseEmigFile(void)
{
string header;
int simul;
int densdep,usefullkern,stagedep,sexdep,indvar,emigstage,stage,sex;
bool densdepset,indvarset;
#if GOBYMODEL
float asocD;
#endif
float	ep,d0,alpha,beta,epMean,epSD,d0Mean,d0SD,alphaMean,alphaSD,betaMean,betaSD;
float epScale,d0Scale,alphaScale,betaScale;
float epScale0 = 0,d0Scale0 = 0,alphaScale0 = 0,betaScale0 = 0;
int errors = 0;
int simuls = 0;
string filetype = "EmigrationFile";

densdepset = false; 
indvarset = false;
ep = 0.0;

// Parse header line;
bEmigrationFile >> header; if (header != "Simulation" ) errors++;
bEmigrationFile >> header; if (header != "DensDep" ) errors++;
bEmigrationFile >> header; if (header != "UseFullKern" ) errors++;
bEmigrationFile >> header; if (header != "StageDep" ) errors++;
bEmigrationFile >> header; if (header != "SexDep" ) errors++;
bEmigrationFile >> header; if (header != "IndVar" ) errors++;
bEmigrationFile >> header; if (header != "EmigStage" ) errors++;
bEmigrationFile >> header; if (header != "Stage" ) errors++;
bEmigrationFile >> header; if (header != "Sex" ) errors++;
#if GROUPDISP
// ADDITIONAL PARAMETERS FOR GROUP DISPERSAL MODEL
bEmigrationFile >> header; if (header != "GroupDisp" ) errors++;
bEmigrationFile >> header; if (header != "GroupMean" ) errors++;
#endif
#if GOBYMODEL
bEmigrationFile >> header; if (header != "Dasocial" ) errors++;
#endif
bEmigrationFile >> header; if (header != "EP" ) errors++;
bEmigrationFile >> header; if (header != "D0" ) errors++;
bEmigrationFile >> header; if (header != "alpha" ) errors++;
bEmigrationFile >> header; if (header != "beta" ) errors++;
bEmigrationFile >> header; if (header != "EPMean" ) errors++;
bEmigrationFile >> header; if (header != "EPSD" ) errors++;
bEmigrationFile >> header; if (header != "D0Mean" ) errors++;
bEmigrationFile >> header; if (header != "D0SD" ) errors++;
bEmigrationFile >> header; if (header != "alphaMean" ) errors++;
bEmigrationFile >> header; if (header != "alphaSD" ) errors++;
bEmigrationFile >> header; if (header != "betaMean" ) errors++;
bEmigrationFile >> header; if (header != "betaSD" ) errors++;
bEmigrationFile >> header; if (header != "EPScale" ) errors++;
bEmigrationFile >> header; if (header != "D0Scale" ) errors++;
bEmigrationFile >> header; if (header != "alphaScale" ) errors++;
bEmigrationFile >> header; if (header != "betaScale" ) errors++;
if (errors > 0) {
	FormatError(filetype,errors);
	return -111;
}

// Parse data lines
int line = 1;
simCheck current{}, prev{};
simul = -98765;
prev.simul = -999;
prev.simlines = prev.reqdsimlines = 0;
bEmigrationFile >> simul;
// first simulation number must match first one in parameterFile
if (simul != firstsimul) {
	BatchError(filetype,line,111,"Simulation"); errors++;
}
while (simul != -98765) {
	// read and validate columns relating to stage and sex-dependency and to IIV
	bEmigrationFile >> densdep >> usefullkern >> stagedep >> sexdep;
	bEmigrationFile >> indvar >> emigstage >> stage >> sex;
	current = CheckStageSex(filetype,line,simul,prev,stagedep,sexdep,stage,sex,indvar,true,false);
	if (current.newsimul) simuls++;
	errors += current.errors;
	prev = current;
	// validate density dependency
	if (densdep < 0 || densdep > 1) {
		BatchError(filetype,line,1,"DensDep"); errors++;
	}
	// validate use full kernel
	if (usefullkern < 0 || usefullkern > 1) {
		BatchError(filetype,line,1,"UseFullKern"); errors++;
	}
	if (densdep != 0) {
		if (usefullkern != 0) {
			BatchError(filetype,line,0,"UseFullKern"); errors++;
			batchlog << "UseFullKern must be 0 if there is density-dependent emigration" << endl;
		}
	}
	// validate emigration stage
	if (stagestruct && !stagedep && indvar) {
		if (stage == 0 && sex == 0) {
			if (emigstage < 0 || emigstage >= stages) {
				BatchError(filetype,line,0,"EmigStage"); errors++;
				batchlog << "EmigStage must be from 0 to " << Int2Str(stages-1) << endl;
			}
		}
	}
	if (stage == 0 && sex == 0) { // first line of a simulation
		// record whether density dependence and individual variability are applied
		if (densdep == 1) densdepset = true; else densdepset = false;
		if (indvar == 1)  indvarset = true;  else indvarset = false;
	}

	// read remaining columns of the current record
#if GROUPDISP
// ADDITIONAL PARAMETERS FOR GROUP DISPERSAL MODEL
	int groupdisp;
	float groupmean;
	bEmigrationFile >> groupdisp >> groupmean;
	if (stage == 0 && sex == 0) { // first line of a simulation
		if (groupdisp < 0 || groupdisp > 1) {
			BatchError(filetype,line,1,"GroupDisp"); errors++;
		}
		if (groupdisp == 1) {
			if (groupmean <= 1.0) {
				BatchError(filetype,line,0," "); errors++;
				batchlog << "Mean group size must be > 1.0" << endl;
			}
		}
	}
#endif // GROUPDISP
#if GOBYMODEL
	bEmigrationFile >> asocD;
	if (densdep == 1 && stage == 0 && sex == 0) { // first line of a DD simulation
		if (asocD <= 0.0) {
			BatchError(filetype,line,10,"Dasocial"); errors++;
		}
	}
#endif // GOBYMODEL
	bEmigrationFile >> ep >> d0 >> alpha >> beta >> epMean >> epSD >> d0Mean >> d0SD;
	bEmigrationFile >> alphaMean >> alphaSD >> betaMean >> betaSD;
	bEmigrationFile >> epScale >> d0Scale >> alphaScale >> betaScale;
#if RSDEBUG
//DEBUGLOG << "ParseEmigFile(): simul=" << simul
//	<< " reqdsimlines=" << current.reqdsimlines
//	<< " line=" << line << " stage=" << stage << " sex=" << sex
//	<< " densdep=" << densdep << " indvar=" << indvar
//	<< " ep=" << ep << " d0=" << d0 << " alpha=" << alpha << " beta=" << beta
//	<< " epMean=" << epMean << " epSD=" << epSD
//	<< " d0Mean=" << d0Mean << " d0SD=" << d0SD
//	<< " alphaMean=" << alphaMean << " alphaSD=" << alphaSD
//	<< " betaMean=" << betaMean << " betaSD=" << betaSD
//	<< " epScale=" << epScale << " d0Scale=" << d0Scale
//	<< " alphaScale=" << alphaScale << " betaScale=" << betaScale
//	<< endl;
#endif
	if (current.newsimul) {
		// record scaling factors from first line of the simulation
		epScale0 = epScale; d0Scale0 = d0Scale;
		alphaScale0 = alphaScale; betaScale0 = betaScale;
	}

	if (densdepset) {
		if (indvarset) {
			if (d0Mean <= 0.0 || d0Mean > 1.0) {
				BatchError(filetype,line,20,"D0Mean"); errors++;
			}
			if (d0SD <= 0.0 || d0SD > 1.0) {
				BatchError(filetype,line,20,"D0SD"); errors++;
			}
			if (stage == 0 && sex == 0) {
				if (d0Scale <= 0.0 || d0Scale > 1.0) {
					BatchError(filetype,line,20,"D0Scale"); errors++;
				}
			}
			if (d0SD > d0Scale0) {
				BatchError(filetype,line,3,"D0SD","D0Scale (first line)"); errors++;
			}
			if (alphaSD <= 0.0) {
				BatchError(filetype,line,10,"alphaSD"); errors++;
			}
			if (stage == 0 && sex == 0) {
				if (alphaScale0 <= 0.0) {
					BatchError(filetype,line,10,"alphaScale"); errors++;
				}
			}
			if (alphaSD > alphaScale0) {
				BatchError(filetype,line,3,"alphaSD","alphaScale (first line)"); errors++;
			}
			if (betaSD <= 0.0) {
				BatchError(filetype,line,10,"betaSD"); errors++;
			}
			if (stage == 0 && sex == 0) {
				if (betaScale0 <= 0.0) {
					BatchError(filetype,line,10,"betaScale"); errors++;
				}
			}
			if (betaSD > betaScale0) {
				BatchError(filetype,line,3,"betaSD","betaScale (first line)"); errors++;
			}
		}
		else { // !indvarset
			if (d0 < 0.0 || d0 > 1.0) {
				BatchError(filetype,line,20,"D0"); errors++;
			}
			// NB alpha and beta may take any value
		}
	}
	else { // !densdepset
		if (indvarset) {
			if (epMean <= 0.0 || epMean > 1.0) {
				BatchError(filetype,line,20,"EPMean"); errors++;
			}
			if (epSD <= 0.0 || epSD > 1.0) {
				BatchError(filetype,line,20,"EPSD"); errors++;
			}
			if (stage == 0 && sex == 0) {
				if (epScale <= 0.0 || epScale > 1.0) {
					BatchError(filetype,line,20,"EPScale"); errors++;
				}
			}
			if (epSD > epScale0) {
				BatchError(filetype,line,3,"EPSD","EPScale (first line)"); errors++;
			}
		}
		else { // !indvarset
			if (ep < 0.0 || ep > 1.0) {
				BatchError(filetype,line,20,"EP"); errors++;
			}
		}
	}

	// read next simulation
	line++;
	simul = -98765;
	bEmigrationFile >> simul;
	if (bEmigrationFile.eof()) simul = -98765;
} // end of while loop
// check for correct number of lines for previous simulation
if (current.simlines != current.reqdsimlines) {
	BatchError(filetype,line,0," "); errors++;
	batchlog << msgnlines << current.simul
		<< msgshldbe << current.reqdsimlines << endl;
}
if (!bEmigrationFile.eof()) {
	EOFerror(filetype);
	errors++;
}

if (errors > 0) return -111;
else return simuls;

}

//---------------------------------------------------------------------------
int ParseTransferFile(string indir)
{
string header,colheader,intext,fname,ftype;
int i,simul,stagedep,sexdep,kerneltype,distmort,indvar,stage,sex;
int	prMethod,smtype,straightenPath;
float pr,dp,smconst;
int goaltype,memsize,betaDB; float gb,alphaDB;
float dpMean,dpSD,gbMean,gbSD,alphaDBMean,alphaDBSD,betaDBMean,betaDBSD;
float dpScale,gbScale,alphaDBScale,betaDBScale;
float meanDistI,meanDistII,ProbKernelI;
float DistIMean,DistISD,DistIIMean,DistIISD,ProbKernelIMean,ProbKernelISD;
float DistIScale,DistIIScale,ProbKernelIScale;
float DistIScale0,DistIIScale0,ProbKernelIScale0;
#if RS_CONTAIN
float u0Kernel1,p0Kernel1,u0Kernel2,p0Kernel2,propKernel1;
float meanU,sigma_w,hc,hr,vt,kappa,dirnmean,dirnsd;
#endif // RS_CONTAIN 
float mortProb,slope,inflPoint;
float morthab,mortmatrix;
int costhab,costmatrix;
float SL,rho;
float StepLMean,StepLSD,RhoMean,RhoSD,StepLScale,RhoScale;

vector <string> costsfiles;
#if TEMPMORT
bool checkfile; 
vector <string> mortfiles;
#endif

int errors = 0; int morthaberrors = 0; int costerrors = 0; int hrerrors = 0;
int simuls = 0;
string filetype = "TransferFile";

DistIScale0 = 0.0; 
DistIIScale0 = 0.0;
ProbKernelIScale0 = 0.0;

// Parse header line;
bTransferFile >> header; if (header != "Simulation" ) errors++;
switch (transfer) {

case 0: { // negative exponential dispersal kernel
	batchlog << "Checking dispersal kernel format file" << endl;
	bTransferFile >> header; if (header != "StageDep" ) errors++;
	bTransferFile >> header; if (header != "SexDep" ) errors++;
	bTransferFile >> header; if (header != "KernelType" ) errors++;
	bTransferFile >> header; if (header != "DistMort" ) errors++;
	bTransferFile >> header; if (header != "IndVar" ) errors++;
	bTransferFile >> header; if (header != "Stage" ) errors++;
	bTransferFile >> header; if (header != "Sex" ) errors++;
	bTransferFile >> header; if (header != "meanDistI" ) errors++;
	bTransferFile >> header; if (header != "meanDistII" ) errors++;
	bTransferFile >> header; if (header != "ProbKernelI" ) errors++;
	bTransferFile >> header; if (header != "DistIMean" ) errors++;
	bTransferFile >> header; if (header != "DistISD" ) errors++;
	bTransferFile >> header; if (header != "DistIIMean" ) errors++;
	bTransferFile >> header; if (header != "DistIISD" ) errors++;
	bTransferFile >> header; if (header != "ProbKernelIMean" ) errors++;
	bTransferFile >> header; if (header != "ProbKernelISD" ) errors++;
	bTransferFile >> header; if (header != "DistIScale" ) errors++;
	bTransferFile >> header; if (header != "DistIIScale" ) errors++;
	bTransferFile >> header; if (header != "ProbKernelIScale" ) errors++;
	bTransferFile >> header; if (header != "MortProb" ) errors++;
	bTransferFile >> header; if (header != "Slope" ) errors++;
	bTransferFile >> header; if (header != "InflPoint" ) errors++;
	break;
} // end of negative exponential dispersal kernel

case 1: { // SMS
	batchlog << "Checking SMS format file ";
	bTransferFile >> header; if (header != "IndVar" ) errors++;
	bTransferFile >> header; if (header != "PR" ) errors++;
	bTransferFile >> header; if (header != "PRMethod" ) errors++;
	bTransferFile >> header; if (header != "DP" ) errors++;
	bTransferFile >> header; if (header != "MemSize" ) errors++;
	bTransferFile >> header; if (header != "GB" ) errors++;
	bTransferFile >> header; if (header != "GoalType" ) errors++;
	bTransferFile >> header; if (header != "AlphaDB" ) errors++;
	bTransferFile >> header; if (header != "BetaDB" ) errors++;
	bTransferFile >> header; if (header != "DPMean" ) errors++;
	bTransferFile >> header; if (header != "DPSD" ) errors++;
	bTransferFile >> header; if (header != "GBMean" ) errors++;
	bTransferFile >> header; if (header != "GBSD" ) errors++;
	bTransferFile >> header; if (header != "AlphaDBMean" ) errors++;
	bTransferFile >> header; if (header != "AlphaDBSD" ) errors++;
	bTransferFile >> header; if (header != "BetaDBMean" ) errors++;
	bTransferFile >> header; if (header != "BetaDBSD" ) errors++;
	bTransferFile >> header; if (header != "DPScale" ) errors++;
	bTransferFile >> header; if (header != "GBScale" ) errors++;
	bTransferFile >> header; if (header != "AlphaDBScale" ) errors++;
	bTransferFile >> header; if (header != "BetaDBScale" ) errors++;
#if TEMPMORT
	bTransferFile >> header; if (header != "StraightenPath" ) errors++;
	bTransferFile >> header; if (header != "SMtype" ) errors++;
	bTransferFile >> header; if (header != "SMconst" ) errors++;
	bTransferFile >> header; if (header != "MortFile" ) errors++;
#else
	bTransferFile >> header; if (header != "StraightenPath" ) errors++;
	bTransferFile >> header; if (header != "SMtype" ) errors++;
	bTransferFile >> header; if (header != "SMconst" ) errors++;
#endif
	switch (landtype) {
	case 0: { // raster map with unique habitat codes
		batchlog << "for LandType = 0" << endl;
		for (i = 0; i < maxNhab; i++) {
			colheader = "MortHab" + Int2Str(i+1);
			bTransferFile >> header; if (header != colheader ) morthaberrors++;
		}
		for (i = 0; i < maxNhab; i++) {
			colheader = "CostHab" + Int2Str(i+1);
			bTransferFile >> header; if (header != colheader ) costerrors++;
		}
		break;
	} // end of raster map with unique habitat codes
	case 2: { // raster map with habitat quality
		batchlog << "for LandType = 2" << endl;
		break;
	} // end of raster map with habitat quality
	case 9: { // artificial landscape
		batchlog << "for LandType = 9" << endl;
		bTransferFile >> header; if (header != "MortHabitat" ) errors++;
		bTransferFile >> header; if (header != "MortMatrix" ) errors++;
		bTransferFile >> header; if (header != "CostHabitat" ) errors++;
		bTransferFile >> header; if (header != "CostMatrix" ) errors++;
		break;
	} // end of artificial landscape
	} // end of switch (landtype)
	break;
} // end of SMS

case 2: { // CRW
	batchlog << "Checking CRW format file" << endl;
	bTransferFile >> header; if (header != "IndVar" ) errors++;
	bTransferFile >> header; if (header != "SL" ) errors++;
	bTransferFile >> header; if (header != "Rho" ) errors++;
	bTransferFile >> header; if (header != "StepLMean" ) errors++;
	bTransferFile >> header; if (header != "StepLSD" ) errors++;
	bTransferFile >> header; if (header != "RhoMean" ) errors++;
	bTransferFile >> header; if (header != "RhoSD" ) errors++;
	bTransferFile >> header; if (header != "StepLScale" ) errors++;
	bTransferFile >> header; if (header != "RhoScale" ) errors++;
	bTransferFile >> header; if (header != "StraightenPath" ) errors++;
	bTransferFile >> header; if (header != "SMtype" ) errors++;
	bTransferFile >> header; if (header != "SMconst" ) errors++;
	if (landtype == 0) {
		for (i = 0; i < maxNhab; i++) {
			colheader = "MortHab" + Int2Str(i+1);
			bTransferFile >> header; if (header != colheader ) morthaberrors++;
		}
	}
	break;
} // end of CRW

#if RS_CONTAIN

case 3: { // 2Dt dispersal kernel
	batchlog << "Checking 2DT dispersal kernel format file" << endl;
	bTransferFile >> header; if (header != "DistMort" ) errors++;
	bTransferFile >> header; if (header != "U0Kernel1" ) errors++;
	bTransferFile >> header; if (header != "P0Kernel1" ) errors++;
	bTransferFile >> header; if (header != "U0Kernel2" ) errors++;
	bTransferFile >> header; if (header != "P0Kernel2" ) errors++;
	bTransferFile >> header; if (header != "PropKernel1" ) errors++;
	bTransferFile >> header; if (header != "MortProb" ) errors++;
	bTransferFile >> header; if (header != "Slope" ) errors++;
	bTransferFile >> header; if (header != "InflPoint" ) errors++;
	break;
} // end of 2Dt dispersal kernel

case 4: { // WALD dispersal kernel
	batchlog << "Checking WALD dispersal kernel format file" << endl;
	bTransferFile >> header; if (header != "DistMort" ) errors++;
	bTransferFile >> header; if (header != "MeanU" ) errors++;
	bTransferFile >> header; if (header != "SigmaW" ) errors++;
	bTransferFile >> header; if (header != "Hc" ) errors++;
	for (i = 1; i < stages; i++) {
		colheader = "Hr" + Int2Str(i);
		bTransferFile >> header; if (header != colheader ) hrerrors++;
	}
	bTransferFile >> header; if (header != "Vt" ) errors++;
	bTransferFile >> header; if (header != "Kappa" ) errors++;
	bTransferFile >> header; if (header != "DirnMean" ) errors++;
	bTransferFile >> header; if (header != "DirnSD" ) errors++;
	bTransferFile >> header; if (header != "MortProb" ) errors++;
	bTransferFile >> header; if (header != "Slope" ) errors++;
	bTransferFile >> header; if (header != "InflPoint" ) errors++;
	break;
} // end of WALD dispersal kernel

#endif // RS_CONTAIN 

} // end of switch (transfer)
// report any errors in headers, and if so, terminate validation
if (errors > 0 || morthaberrors > 0 || costerrors > 0 || hrerrors > 0) {
	FormatError(filetype,errors+morthaberrors+costerrors);
	if (morthaberrors > 0) BatchError(filetype,-999,333,"MortHab");
	if (costerrors > 0) BatchError(filetype,-999,333,"CostHab");
	if (hrerrors > 0) BatchError(filetype,-999,444,"Hr");
	return -111;
}

// Parse data lines
int line = 1;
simCheck current{}, prev{};
simul = -98765;
prev.simul = -999;
prev.simlines = prev.reqdsimlines = 0;
bTransferFile >> simul;
// first simulation number must match first one in parameterFile
if (simul != firstsimul) {
	BatchError(filetype,line,111,"Simulation"); errors++;
}
while (simul != -98765) {

	switch (transfer) {
	
	case 0: { // negative exponential dispersal kernel
		// read and validate columns relating to stage and sex-dependency and to IIV
		bTransferFile >> stagedep >> sexdep >> kerneltype >> distmort;
		bTransferFile >> indvar >> stage >> sex;
		current = CheckStageSex(filetype,line,simul,prev,stagedep,sexdep,stage,sex,indvar,true,false);
		if (current.newsimul) simuls++;
		errors += current.errors;
		prev = current;
		// validate kernel type
		if (kerneltype < 0 || kerneltype > 1) {
			BatchError(filetype,line,1,"KernelType"); errors++;
		}
		// validate mortality
		if (distmort < 0 || distmort > 1) {
			BatchError(filetype,line,1,"DistMort"); errors++;
		}
		// read remaining columns of the current record
		bTransferFile >> meanDistI >> meanDistII >> ProbKernelI >> DistIMean >> DistISD;
		bTransferFile	>> DistIIMean >> DistIISD >> ProbKernelIMean >> ProbKernelISD;
		bTransferFile	>> DistIScale >> DistIIScale >> ProbKernelIScale;
		bTransferFile >> mortProb	>> slope >> inflPoint;

		if (current.newsimul) {
			DistIScale0 = DistIScale; DistIIScale0 = DistIIScale;
			ProbKernelIScale0 = ProbKernelIScale;
		}

		if (!indvar) {
			if (meanDistI < resolution) {
				// NOTE - should also check whether emigration prob is constant and equal to 1
				//but checks across diffferent input files are not yet implemented
				BatchError(filetype,line,2,"meanDistI","Resolution"); errors++;
			}
			if (kerneltype != 0) {
				if (meanDistII < resolution) {
					// NOTE - DITTO
					BatchError(filetype,line,2,"meanDistII","Resolution"); errors++;
				}
				if (ProbKernelI <= 0.0 || ProbKernelI >= 1.0) {
					BatchError(filetype,line,20,"ProbKernelI"); errors++;
				}
			}
		}

		if (!stagestruct && indvar) {
			if (DistIMean < resolution) {
				// NOTE - DITTO
				BatchError(filetype,line,2,"DistIMean","Resolution"); errors++;
			}
			if (DistISD <= 0.0) {
				BatchError(filetype,line,10,"DistISD"); errors++;
			}
			if (DistISD > DistIScale0) {
				BatchError(filetype,line,3,"DistISD","DistIScale (first line)"); errors++;
			}
			if (current.newsimul) {
				if (DistIScale <= 0.0) {
					BatchError(filetype,line,10,"DistIScale"); errors++;
				}
			}
			if (kerneltype != 0) {
				if (DistIIMean < resolution) {
					// NOTE - DITTO
					BatchError(filetype,line,2,"DistIIMean","Resolution"); errors++;
				}
				if (DistIISD <= 0.0) {
					BatchError(filetype,line,10,"DistIISD"); errors++;
				}
				if (DistIISD > DistIIScale0) {
					BatchError(filetype,line,3,"DistIISD","DistIIScale (first line)"); errors++;
				}
				if (current.newsimul) {
					if (DistIIScale <= 0.0) {
						BatchError(filetype,line,10,"DistIIScale"); errors++;
					}
				}
				if (ProbKernelIMean <= 0.0 || ProbKernelIMean > 1.0) {
					BatchError(filetype,line,20,"ProbKernelIMean"); errors++;
				}
				if (ProbKernelISD <= 0.0 || ProbKernelISD > 1.0) {
					BatchError(filetype,line,20,"ProbKernelISD"); errors++;
				}
				if (ProbKernelISD > ProbKernelIScale0) {
					BatchError(filetype,line,3,"ProbKernelISD","ProbKernelIScale (first line)"); errors++;
				}
				if (current.newsimul) {
					if (ProbKernelIScale <= 0.0 || ProbKernelIScale > 1.0) {
						BatchError(filetype,line,20,"ProbKernelIScale"); errors++;
					}
				}
			}
		}

		if (stage == 0 && sex == 0) {
			if (distmort) { // distance-dependent mortality
				// WHAT CONDITIONS APPLY TO MORTALITY SLOPE AND INFLECTION POINT?
			}
			else { // constant mortality
				if (mortProb < 0.0 || mortProb >= 1.0) {
					BatchError(filetype,line,20,"MortProb"); errors++;
				}
			}
		}

		break;
	} // end of negative exponential dispersal kernel

	case 1: { // SMS
		bTransferFile >> indvar;
		bTransferFile	>> pr >> prMethod	>> dp;
		bTransferFile	>> memsize >> gb	>> goaltype >> alphaDB >> betaDB;
		current = CheckStageSex(filetype,line,simul,prev,0,0,0,0,0,true,false);
		if (current.newsimul) simuls++;
		errors += current.errors;
		prev = current;
		// validate SMS movement parameters
		if (pr < 1) {
			BatchError(filetype,line,11,"PR"); errors++;
		}
		if (prMethod < 1 || prMethod > 3) {
			BatchError(filetype,line,33,"PRmethod"); errors++;
		}
		if (!indvar && dp < 1.0) {
			BatchError(filetype,line,11,"DP"); errors++;
		}
		if (memsize < 1 || memsize > 14) {
			BatchError(filetype,line,0,"MemSize"); errors++;
			batchlog << "MemSize must be from 1 to 14" << endl;
		}
		if (!indvar && gb < 1.0) {
			BatchError(filetype,line,11,"GB"); errors++;
		}
		if (goaltype < 0 || goaltype > 2) {
			BatchError(filetype,line,2,"GoalType"); errors++;
		}
		if (!indvar && goaltype == 2) { // dispersal bias
			if (alphaDB <= 0.0) {
				BatchError(filetype,line,10,"AlphaDB"); errors++;
			}
			if (betaDB <= 0.0) {
				BatchError(filetype,line,10,"BetaDB"); errors++;
			}
		}
		bTransferFile >> dpMean >> dpSD >> gbMean >> gbSD >> alphaDBMean >> alphaDBSD
			>> betaDBMean >> betaDBSD >> dpScale >> gbScale >> alphaDBScale >> betaDBScale;
		if (indvar) {
			if (dpMean < 1.0) {
				BatchError(filetype,line,11,"DPMean"); errors++;
			}
			if (dpSD <= 0.0) {
				BatchError(filetype,line,10,"DPSD"); errors++;
			}
			if (dpScale <= 0.0) {
				BatchError(filetype,line,10,"DPScale"); errors++;
			}
			if (dpSD > dpScale) {
				BatchError(filetype,line,3,"DPSD","DPScale"); errors++;
			}
			if (gbMean < 1.0) {
				BatchError(filetype,line,11,"GBMean"); errors++;
			}
			if (gbSD <= 0.0) {
				BatchError(filetype,line,10,"GBSD"); errors++;
			}
			if (gbScale <= 0.0) {
				BatchError(filetype,line,10,"GBScale"); errors++;
			}
			if (gbSD > gbScale) {
				BatchError(filetype,line,3,"GBSD","GBScale"); errors++;
			}
			if (goaltype == 2) { // dispersal bias
				if (alphaDBMean <= 0.0) {
					BatchError(filetype,line,10,"AlphaDBMean"); errors++;
				}
				if (alphaDBSD <= 0.0) {
					BatchError(filetype,line,10,"AlphaDBSD"); errors++;
				}
				if (alphaDBScale <= 0.0) {
					BatchError(filetype,line,10,"AlphaDBScale"); errors++;
				}
				if (alphaDBSD > alphaDBScale) {
					BatchError(filetype,line,3,"AlphaDBSD","AlphaDBScale"); errors++;
				}
				if (betaDBMean < 1.0) {
					BatchError(filetype,line,11,"BetaDBMean"); errors++;
				}
				if (betaDBSD <= 0.0) {
					BatchError(filetype,line,10,"BetaDBSD"); errors++;
				}
				if (betaDBScale <= 0.0) {
					BatchError(filetype,line,10,"BetaDBScale"); errors++;
				}
				if (betaDBSD > betaDBScale) {
					BatchError(filetype,line,3,"BetaDBSD","BetaDBScale"); errors++;
				}
			}
		}
#if TEMPMORT
		bTransferFile >> straightenPath	>> smtype >> smconst >> intext;
#else
		bTransferFile >> straightenPath	>> smtype >> smconst;
#endif // TEMPMORT 
		if (straightenPath < 0 || straightenPath > 1) {
			BatchError(filetype,line,1,"StraightenPath"); errors++;
		}
		if (landtype == 2) // habitat quality landscape 
#if TEMPMORT
		{ // must have constant or temporally variable mortality 
			if (smtype != 0 && smtype != 2) {
				BatchError(filetype,line,0," "); errors++;
				batchlog << "SMtype must be 0 or 2 for LandType 2" << endl;
			}
		}
#else
		{ // must have constant mortality
			if (smtype != 0) {
				BatchError(filetype,line,0," "); errors++;
				batchlog << "SMtype must be 0 for LandType 2" << endl;
			}
		}
#endif // TEMPMORT 
		else {
#if TEMPMORT
			if (smtype < 0 || smtype > 2) {
				BatchError(filetype,line,2,"SMtype"); errors++;
			}
#else
			if (smtype < 0 || smtype > 1) {
				BatchError(filetype,line,1,"SMtype"); errors++;
			}
#endif // TEMPMORT
		}
#if TEMPMORT
		if (smtype == 0 || smtype == 2) 
#else
		if (smtype == 0) 
#endif // TEMPMORT 
		{
			if (smconst < 0.0 || smconst >= 1.0) {
				BatchError(filetype,line,20,"SMconst"); errors++;
			}
		}
#if TEMPMORT
		ftype = "MortFile";
		if (smtype == 2) {
			if (intext == "NULL") {
				BatchError(filetype,line,0," "); errors++;
				batchlog << ftype << " is required if SMtype = 2" << endl;
			}
			else {
				checkfile = true;
				for (i = 0; i < (int)mortfiles.size(); i++) {
					if (intext == mortfiles[i]) { // file has already been checked
					checkfile = false;
					}
				}
				if (checkfile) {
					fname = indir + intext;
					batchlog << "Checking " << ftype << " " << fname << endl;
					bMortFile.open(fname.c_str());
					if (bMortFile.is_open()) {
						int err = ParseMortFile();
						if (err == 0) FileHeadersOK(ftype); else errors++;
						bMortFile.close();
					}	
					else {
						OpenError(ftype,fname); errors++;
					}
				if (bMortFile.is_open()) bMortFile.close();
				bMortFile.clear();
				} // end of checkfile
				mortfiles.push_back(intext);
			} // end of intext != NULL
		}
		else { 
			if (intext != "NULL") {
				BatchError(filetype,line,0," "); errors++;
				batchlog << ftype << " should be NULL if SMtype is not 2" << endl;
			}
		}
#endif // TEMPMORT 
		switch (landtype) {
		
		case 0: { // raster map with unique habitat codes
//			batchlog << "for LandType = 0" << endl;
			for (i = 0; i < maxNhab; i++) {
				bTransferFile >> morthab;
				if (smtype == 1) 
				{
					if (morthab < 0.0 || morthab >= 1.0) {
						colheader = "MortHab" + Int2Str(i+1);
						BatchError(filetype,line,20,colheader); errors++;
					}
				}
			}
			for (i = 0; i < maxNhab; i++) {
				bTransferFile >> costhab;
				if (name_costfile == "NULL") {
					if (costhab < 1) {
						colheader = "CostHab" + Int2Str(i+1);
						BatchError(filetype,line,11,colheader); errors++;
					}
				}
			}
			break;
		} // end of raster map with unique habitat codes
		
		case 2: { // raster map with habitat quality
//			batchlog << "for LandType = 2" << endl;
			break;
		} // end of raster map with habitat quality
		
		case 9: { // artificial landscape
//			batchlog << "for LandType = 9" << endl;
			bTransferFile >> morthab >> mortmatrix;
			bTransferFile >> costhab >> costmatrix;
			if (smtype) { // validate habitat-dependent mortality
				if (morthab < 0.0 || morthab >= 1.0) {
					BatchError(filetype,line,20,"MortHabitat"); errors++;
				}
				if (mortmatrix < 0.0 || mortmatrix >= 1.0) {
					BatchError(filetype,line,20,"MortMatrix"); errors++;
				}
			}
			if (costhab < 1) {
				BatchError(filetype,line,11,"CostHabitat"); errors++;
			}
			if (costmatrix < 1) {
				BatchError(filetype,line,11,"CostMatrix"); errors++;
			}
			break;
		} // end of artificial landscape
		
		} // end of switch (landtype)
		
		break;
		
	} // end of SMS

	case 2: { // CRW
		bTransferFile	>> indvar >> SL	>> rho >> StepLMean	>> StepLSD >> RhoMean >> RhoSD;
		bTransferFile >> StepLScale >> RhoScale >> straightenPath	>> smtype >> smconst;
		current = CheckStageSex(filetype,line,simul,prev,0,0,0,0,indvar,true,false);
		if (current.newsimul) simuls++;
		errors += current.errors;
		prev = current;

		if (indvar) { // individual variability
			if (StepLMean <= 0.0) {
				BatchError(filetype,line,10,"StepLMean"); errors++;
			}
			if (StepLSD <= 0.0) {
				BatchError(filetype,line,10,"StepLSD"); errors++;
			}
			if (StepLScale <= 0.0) {
				BatchError(filetype,line,10,"StepLScale"); errors++;
			}
			if (StepLSD > StepLScale) {
				BatchError(filetype,line,3,"StepLSD","StepLScale"); errors++;
			}
			if (RhoMean <= 0.0 || RhoMean >= 1.0) {
				BatchError(filetype,line,20,"RhoMean"); errors++;
			}
			if (RhoSD <= 0.0 || RhoSD >= 1.0) {
				BatchError(filetype,line,20,"RhoSD"); errors++;
			}
			if (RhoScale <= 0.0 || RhoScale >= 1.0) {
				BatchError(filetype,line,20,"RhoScale"); errors++;
			}
			if (RhoSD > RhoScale) {
				BatchError(filetype,line,3,"RhoSD","RhoScale"); errors++;
			}
		}
		else { // no individual variability
			if (SL <= 0.0) {
				BatchError(filetype,line,10,"SL"); errors++;
			}
			if (rho <= 0.0 || rho >= 1.0) {
				BatchError(filetype,line,20,"Rho"); errors++;
			}
		}
		if (straightenPath < 0 || straightenPath > 1) {
			BatchError(filetype,line,1,"StraightenPath"); errors++;
		}
		if (landtype == 0) { // real landscape with habitat types
			if (smtype < 0 || smtype > 1) {
				BatchError(filetype,line,1,"SMtype"); errors++;
			}
			if (!smtype) {
				if (smconst < 0.0 || smconst >= 1.0) {
					BatchError(filetype,line,20,"SMconst"); errors++;
				}
			}
			for (int i = 0; i < maxNhab; i++) {
				bTransferFile >> morthab;
				if (smtype) {
					if (morthab < 0.0 || morthab >= 1.0) {
						colheader = "MortHab" + Int2Str(i+1);
						BatchError(filetype,line,20,colheader); errors++;
					}
				}
			}
		}
		else { // real landscape with quality OR artificial landscape
			if (smtype != 0) {
				BatchError(filetype,line,0," "); errors++;
				batchlog << "SMtype must be 0 for LandType 2 or 9" << endl;
			}
			if (smconst < 0.0 || smtype >= 1.0) {
				BatchError(filetype,line,20,"SMconst"); errors++;
			}
		}
		break;
	} // end of CRW

#if RS_CONTAIN

	case 3: { // 2Dt dispersal kernel
		// read and validate columns relating to stage and sex-dependency and to IIV
		current = CheckStageSex(filetype,line,simul,prev,0,0,0,0,0,true,false);
		if (current.newsimul) simuls++;
		errors += current.errors;
		prev = current;
		// validate mortality
		bTransferFile >> distmort;          
		if (distmort < 0 || distmort > 1) {
			BatchError(filetype,line,1,"DistMort"); errors++;
		}
		// read remaining columns of the current record
		bTransferFile >> u0Kernel1 >> p0Kernel1 >> u0Kernel2 >> p0Kernel2 >> propKernel1;
		bTransferFile >> mortProb	>> slope >> inflPoint;

		if (u0Kernel1 < 0.0) {
			BatchError(filetype,line,19,"U0Kernel1"); errors++;
		}
		if (p0Kernel1 < 0.0) {
			BatchError(filetype,line,19,"P0Kernel1"); errors++;
		}
		if (u0Kernel2 < 0.0) {
			BatchError(filetype,line,19,"U0Kernel2"); errors++;
		}
		if (p0Kernel2 < 0.0) {
			BatchError(filetype,line,19,"P0Kernel2"); errors++;
		}
		if (propKernel1 < 0.0 || propKernel1 > 1.0) {
			BatchError(filetype,line,20,"PropKernel1"); errors++;
		}

		if (distmort) { // distance-dependent mortality
			// WHAT CONDITIONS APPLY TO MORTALITY SLOPE AND INFLECTION POINT?
		}
		else { // constant mortality
			if (mortProb < 0.0 || mortProb >= 1.0) {
				BatchError(filetype,line,20,"MortProb"); errors++;
			}
		}

		break;
	} // end of 2Dt dispersal kernel

	case 4: { // WALD dispersal kernel
		// read and validate columns relating to stage and sex-dependency and to IIV
		current = CheckStageSex(filetype,line,simul,prev,0,0,0,0,0,true,false);
		if (current.newsimul) simuls++;
		errors += current.errors;
		prev = current;
		// validate mortality
		bTransferFile >> distmort;          
		if (distmort < 0 || distmort > 1) {
			BatchError(filetype,line,1,"DistMort"); errors++;
		}
		// read remaining columns of the current record
		bTransferFile >> meanU >> sigma_w >> hc;
		if (meanU <= 0.0) {
			BatchError(filetype,line,10,"MeanU"); errors++;
		}
		if (sigma_w <= 0.0) {
			BatchError(filetype,line,10,"SigmaW"); errors++;
		}
		if (hc <= 0.0) {
			BatchError(filetype,line,10,"Hc"); errors++;
		}
		for (int i = 1; i < stages; i++) {
			bTransferFile >> hr; 
			colheader = "Hr" + Int2Str(i);
			if (hr <= 0.0) {
				BatchError(filetype,line,10,colheader); errors++;
			}
			if (hr > hc) {
				BatchError(filetype,line,3,colheader,"Hc"); errors++;
			}
		}
		bTransferFile >> vt	>> kappa >> dirnmean >> dirnsd;
		if (vt <= 0.0) {
			BatchError(filetype,line,10,"Vt"); errors++;
		}
		if (kappa <= 0.0) {
			BatchError(filetype,line,10,"Kappa"); errors++;
		}
		if (dirnmean < 0.0 || dirnmean >= 360.0) {
			BatchError(filetype,line,0," "); errors++;
			batchlog << "DirnMean must be >= 0.0 and < 360.0" << endl;
		}
		if (dirnsd <= 0.0) {
			BatchError(filetype,line,10,"DirnSD"); errors++;
		}
		bTransferFile >> mortProb	>> slope >> inflPoint;
		if (distmort) { // distance-dependent mortality
			// WHAT CONDITIONS APPLY TO MORTALITY SLOPE AND INFLECTION POINT?
		}
		else { // constant mortality
			if (mortProb < 0.0 || mortProb >= 1.0) {
				BatchError(filetype,line,20,"MortProb"); errors++;
			}
		}

		break;
	} // end of WALD dispersal kernel

#endif // RS_CONTAIN 
	
	} // end of switch (transfer)

	// read next simulation
	line++;
	simul = -98765;
	bTransferFile >> simul;
	if (bTransferFile.eof()) simul = -98765;
} // end of while loop
// check for correct number of lines for previous simulation
if (transfer == 0 // no. of lines checked for dispersal kernel transfer method only
&&  current.simlines != current.reqdsimlines) {
	BatchError(filetype,line,0," "); errors++;
	batchlog << msgnlines << current.simul
		<< msgshldbe << current.reqdsimlines << endl;
}
if (!bTransferFile.eof()) {
	EOFerror(filetype);
	errors++;
}
costsfiles.clear();

if (errors > 0) return -111;
else return simuls;

}


#if TEMPMORT
//---------------------------------------------------------------------------

int ParseMortFile(void)
{
int year,prevyear;
float gradient;
int errors = 0;
string filetype = "MortFile";

bool yearError = false;
bool gradError = false;
prevyear = year = -98765;
bMortFile >> year;
if (year < 1) {
	BatchError(filetype,-999,0," "); errors++;
	batchlog << "First year must be 1 or more" << endl;
}
while (year != -98765) {
	bMortFile >> gradient;
//	batchlog << "year = " << year << " gradient = " << gradient << endl;
	if (year <= prevyear) yearError = true;
	if (gradient <= -1.0 || gradient >= 1.0) gradError = true;
	prevyear = year;
	year = -98765;
	bMortFile >> year;          
//	batchlog << "paramname=" << paramname << " (end of loop)" << endl;
}
//batchlog << "paramname=" << paramname << " (after loop)" << endl;

if (yearError) {
	BatchError(filetype,-999,0," "); errors++;
	batchlog << "Year must be greater than previous year" << endl;
}
if (gradError) {
	BatchError(filetype,-999,0," "); errors++;
	batchlog << "mortality gradient must be between -1 and +1 exclusive " << endl;
}

// final read should hit EOF
if (!bMortFile.eof()) {
	EOFerror(filetype);
	errors++;
}

return errors;

}
#endif // TEMPMORT

//---------------------------------------------------------------------------
int ParseSettleFile(void)
{
string header;
int simul,stagedep,sexdep,stage,sex,settletype;
int densdep,indvar,findmate,minSteps,maxSteps,maxStepsYear;
float s0,alphaS,betaS;
#if GOBYMODEL
float alphaSasoc,betaSasoc;
#endif
float s0mean,alphaSmean,betaSmean,s0sd,alphaSsd,betaSsd,s0scale,alphaSscale,betaSscale;
float s0scale0 = 0,alphaSscale0 = 0,betaSscale0 = 0;
int errors = 0;
int simuls = 0;
string filetype = "SettlementFile";

// Parse header line;
bSettlementFile >> header; if (header != "Simulation" ) errors++;
bSettlementFile >> header; if (header != "StageDep" ) errors++;
bSettlementFile >> header; if (header != "SexDep" ) errors++;
bSettlementFile >> header; if (header != "Stage" ) errors++;
bSettlementFile >> header; if (header != "Sex" ) errors++;
#if RS_CONTAIN
if (transfer == 0 || transfer == 3 || transfer == 4) 
#else
if (transfer == 0) 
#endif // RS_CONTAIN 
{ // dispersal kernel
	bSettlementFile >> header; if (header != "SettleType" ) errors++;
	bSettlementFile >> header; if (header != "FindMate" ) errors++;
}
else { // movement method
	bSettlementFile >> header; if (header != "DensDep" ) errors++;
	bSettlementFile >> header; if (header != "IndVar" ) errors++;
	bSettlementFile >> header; if (header != "FindMate" ) errors++;
	bSettlementFile >> header; if (header != "MinSteps" ) errors++;
	bSettlementFile >> header; if (header != "MaxSteps" ) errors++;
	bSettlementFile >> header; if (header != "MaxStepsYear" ) errors++;
#if GOBYMODEL
	bSettlementFile >> header; if (header != "AlphaSasocial" ) errors++;
	bSettlementFile >> header; if (header != "BetaSasocial" ) errors++;
#endif
	bSettlementFile >> header; if (header != "S0" ) errors++;
	bSettlementFile >> header; if (header != "AlphaS" ) errors++;
	bSettlementFile >> header; if (header != "BetaS" ) errors++;
	bSettlementFile >> header; if (header != "S0Mean" ) errors++;
	bSettlementFile >> header; if (header != "S0SD" ) errors++;
	bSettlementFile >> header; if (header != "AlphaSMean" ) errors++;
	bSettlementFile >> header; if (header != "AlphaSSD" ) errors++;
	bSettlementFile >> header; if (header != "BetaSMean" ) errors++;
	bSettlementFile >> header; if (header != "BetaSSD" ) errors++;
	bSettlementFile >> header; if (header != "S0Scale" ) errors++;
	bSettlementFile >> header; if (header != "AlphaSScale" ) errors++;
	bSettlementFile >> header; if (header != "BetaSScale" ) errors++;
}
if (errors > 0) {
	FormatError(filetype,errors);
	return -111;
}

// Parse data lines
int line = 1;
simCheck current{}, prev{};
simul = -98765;
prev.simul = -999;
prev.simlines = prev.reqdsimlines = 0;
bSettlementFile >> simul;
// first simulation number must match first one in parameterFile
if (simul != firstsimul) {
	BatchError(filetype,line,111,"Simulation"); errors++;
}
while (simul != -98765) {
#if RS_CONTAIN
	if (transfer == 0 || transfer == 3 || transfer == 4) 
#else
	if (transfer == 0) 
#endif // RS_CONTAIN 
	{ // dispersal kernel
		// read and validate columns relating to stage and sex-dependency (NB no IIV here)
		bSettlementFile >> stagedep >> sexdep >> stage >> sex >> settletype >> findmate;
		current = CheckStageSex(filetype,line,simul,prev,stagedep,sexdep,stage,sex,0,true,false);
		if (current.newsimul) simuls++;
		errors += current.errors;
		prev = current;
		if (settletype < 0 || settletype > 3) {
			BatchError(filetype,line,3,"SettleType"); errors++;
		}
		if (!stagestruct && (settletype == 1 || settletype == 3)) {
			BatchError(filetype,line,0," "); errors++;
			batchlog << "Invalid SettleType for a non-stage-structured population" << endl;
		}
		if (sexesDisp > 1) {
			if (findmate < 0 || findmate > 1) {
				BatchError(filetype,line,1,"FindMate"); errors++;
			}
		}
	}
	else { // movement method
		// read and validate columns relating to stage and sex-dependency (IIV psossible)
		bSettlementFile >> stagedep >> sexdep >> stage >> sex >> densdep >> indvar >> findmate;
		current = CheckStageSex(filetype,line,simul,prev,stagedep,sexdep,stage,sex,indvar,true,false);
		if (current.newsimul) simuls++;
		errors += current.errors;
		prev = current;
		if (densdep < 0 || densdep > 1) {
			BatchError(filetype,line,1,"DensDep"); errors++;
		}
		if (densdep == 0) {
			if (indvar != 0) {
				BatchError(filetype,line,0," "); errors++;
				batchlog << "IndVar must be 0 if DensDep is 0" << endl;
			}
		}
		if (reproductn != 0 && sexesDisp > 1) {
			if (findmate < 0 || findmate > 1) {
				BatchError(filetype,line,1,"FindMate"); errors++;
			}
		}
		bSettlementFile >> minSteps >> maxSteps >> maxStepsYear;
		if (stage == 0 && sex == 0) {
			if (minSteps < 0) {
				BatchError(filetype,line,19,"MinSteps"); errors++;
			}
			if (maxSteps < 0) {
				BatchError(filetype,line,19,"MaxSteps"); errors++;
			}
		}
		if (maxStepsYear < 0) {
			BatchError(filetype,line,19,"MaxStepsYear"); errors++;
		}
#if GOBYMODEL
		bSettlementFile >> alphaSasoc >> betaSasoc;
		if (densdep == 1 && stage == 0 && sex == 0) {
			if (betaSasoc <= 0.0) {
				BatchError(filetype,line,10,"BetaSasocial"); errors++;
			}
		}
#endif
		bSettlementFile >> s0 >> alphaS >> betaS;
		bSettlementFile >> s0mean >> s0sd >> alphaSmean >> alphaSsd
			>> betaSmean >> betaSsd >> s0scale >> alphaSscale >> betaSscale;
#if RSDEBUG
//DEBUGLOG << "ParseSettleFile(): simul=" << simul
//	<< " reqdsimlines=" << current.reqdsimlines
//	<< " line=" << line << " stage=" << stage << " sex=" << sex
//	<< " densdep=" << densdep << " indvar=" << indvar << " findmate=" << findmate
//	<< " s0=" << s0 << " alphaS=" << alphaS << " betaS=" << betaS
//	<< " s0mean=" << s0mean << " s0sd=" << s0sd
//	<< " alphaSmean=" << alphaSmean << " alphaSsd=" << alphaSsd
//	<< " betaSmean=" << betaSmean << " betaSsd=" << betaSsd
//	<< " s0scale=" << s0scale << " alphaSscale=" << alphaSscale << " betaSscale=" << betaSscale
//	<< endl;
#endif
		if (current.newsimul) {
			// record scaling factors from first line of the simulation
			s0scale0 = s0scale;
			alphaSscale0 = alphaSscale; betaSscale0 = betaSscale;
		}
	
		if (densdep == 1) {
			if (indvar == 1) {
//				if (stage == 0 && sex == 0) 
				if (stage == 0)
				{
#if RSDEBUG
//DEBUGLOG << "ParseSettleFile(): validating initial trait parameters" << endl;
#endif
					if (s0mean <= 0.0 || s0mean > 1.0) {
						BatchError(filetype,line,20,"S0Mean"); errors++;
					}
					if (s0sd <= 0.0 || s0sd > 1.0) {
						BatchError(filetype,line,20,"S0SD"); errors++;
					}
					if (alphaSsd <= 0.0) {
						BatchError(filetype,line,10,"AlphaSsd"); errors++;
					}
					if (betaSsd <= 0.0) {
						BatchError(filetype,line,10,"BetaSsd"); errors++;
					}
					if (sex == 0) {
						if (s0scale <= 0.0 || s0scale > 1.0) {
							BatchError(filetype,line,20,"S0Scale"); errors++;
						}
						if (alphaSscale <= 0.0) {
							BatchError(filetype,line,10,"AlphaSscale"); errors++;
						}
						if (betaSscale <= 0.0) {
							BatchError(filetype,line,10,"BetaSscale"); errors++;
						}
					}
					if (s0sd > s0scale0) {     
						BatchError(filetype,line,3,"S0SD","S0Scale (first line)"); errors++;
					}
					if (alphaSsd > alphaSscale0) {
						BatchError(filetype,line,3,"AlphaSsd","AlphaSscale (first line)"); errors++;
					}
					if (betaSsd > betaSscale0) {
						BatchError(filetype,line,3,"BetaSsd","BetaSscale (first line)"); errors++;
					}
				}
			}
			else { // no individual variation
#if RSDEBUG
//DEBUGLOG << "ParseSettleFile(): validating s0 only" << endl;
#endif
				if (s0 <= 0.0 || s0 > 1.0) {
					BatchError(filetype,line,20,"S0"); errors++;
				}
				// NOTE: alphaS and betaS can take any value
			}
		}
	}
	// read next simulation
	line++;
	simul = -98765;
	bSettlementFile >> simul;
	if (bSettlementFile.eof()) simul = -98765;
} // end of while loop
// check for correct number of lines for previous simulation
if (current.simlines != current.reqdsimlines) {
	BatchError(filetype,line,0," "); errors++;
	batchlog << msgnlines << current.simul
		<< msgshldbe << current.reqdsimlines << endl;
}
if (!bSettlementFile.eof()) {
	EOFerror(filetype);
	errors++;
}

if (errors > 0) return -111;
else return simuls;

}

//---------------------------------------------------------------------------
int ParseGeneticsFile(string indir)
{
string header,colheader;
int i,simul,err;
int arch,nLoci;
string filename,ftype,fname;
float probMutn,probCross,alleleSD,mutationSD;
bool checkfile;
int errors = 0;
int simuls = 0;
vector <string> archfiles;
string filetype = "GeneticsFile";

// Parse header line;
bGeneticsFile >> header; if (header != "Simulation" ) errors++;
bGeneticsFile >> header; if (header != "Architecture" ) errors++;
bGeneticsFile >> header; if (header != "NLoci" ) errors++;
bGeneticsFile >> header; if (header != "ArchFile" ) errors++;
bGeneticsFile >> header; if (header != "ProbMutn" ) errors++;
bGeneticsFile >> header; if (header != "ProbCross" ) errors++;
bGeneticsFile >> header; if (header != "AlleleSD" ) errors++;
bGeneticsFile >> header; if (header != "MutationSD" ) errors++;
if (errors > 0) {
	FormatError(filetype,errors);
	return -111;
}

// Parse data lines
int line = 1;
simCheck current{},prev{};
simul = -98765;
prev.simul = -999;
prev.simlines = prev.reqdsimlines = 0;
bGeneticsFile >> simul;
// first simulation number must match first one in parameterFile
if (simul != firstsimul) {
	BatchError(filetype,line,111,"Simulation"); errors++;
}
while (simul != -98765) {
	// read and validate columns relating to stage and sex-dependency (NB no IIV here)
	bGeneticsFile >> arch >> nLoci >> filename
		>> probMutn >> probCross >> alleleSD >> mutationSD;

	current = CheckStageSex(filetype,line,simul,prev,0,0,0,0,0,true,false);
	if (current.newsimul) simuls++;
	errors += current.errors;
	prev = current;

	// validate parameters

	if (arch < 0 || arch > 1) {
		BatchError(filetype,line,1,"Architecture"); errors++;
	}

	// genetic architecture file - optional
	ftype = "ArchFile";
	if (filename == "NULL") {
		if (arch != 0) {
			BatchError(filetype,line,0," ");
			batchlog << ftype << " is compulsory unless Architecture is 0" << endl;
			errors++;
		}
	}
	else {
		if (arch != 1) {
			BatchError(filetype,line,0," ");
			batchlog << ftype << " must be NULL if Architecture is 0" << endl;
			errors++;
		}
		else { // check architecture file
			checkfile = true;
			for (i = 0; i < (int)archfiles.size(); i++) {
				if (filename == archfiles[i]) checkfile = false; // file has already been checked
			}
			if (checkfile) {
				fname = indir + filename;
				batchlog << "Checking " << ftype << " " << fname << endl;
				bArchFile.open(fname.c_str());
				if (bArchFile.is_open()) {
					err = ParseArchFile();
					if (err == 0) FileHeadersOK(ftype); else errors++;
					bArchFile.close();
				}
				else {
					OpenError(ftype,fname); errors++;
				}
				if (bArchFile.is_open()) bArchFile.close();
				bArchFile.clear();
			}
			archfiles.push_back(filename);
		}
	}

	if (arch == 0) {
		if (nLoci < 1) {
			BatchError(filetype,line,11,"NLoci"); errors++;
		}
	}
	if (probMutn < 0.0 || probMutn > 1.0) {
		BatchError(filetype,line,20,"ProbMutn"); errors++;
	}
	if (probCross < 0.0 || probCross > 1.0) {
		BatchError(filetype,line,20,"ProbCross"); errors++;
	}
	if (alleleSD <= 0.0) {
		BatchError(filetype,line,10,"AlleleSD"); errors++;
	}
	if (mutationSD <= 0.0) {
		BatchError(filetype,line,10,"MutationSD"); errors++;
	}

	// read next simulation
	line++;
	simul = -98765;
	bGeneticsFile >> simul;
	if (bGeneticsFile.eof()) simul = -98765;
} // end of while loop
// check for correct number of lines for previous simulation
if (current.simlines != current.reqdsimlines) {
	BatchError(filetype,line,0," "); errors++;
	batchlog << msgnlines << current.simul
		<< msgshldbe << current.reqdsimlines << endl;
}
if (!bGeneticsFile.eof()) {
	EOFerror(filetype);
	errors++;
}

if (errors > 0) return -111;
else return simuls;

}

//---------------------------------------------------------------------------
int ParseArchFile(void)
{
string paramname;
int nchromosomes,nloci;
int errors = 0;
bool formatError = false;
string filetype = "ArchFile";
int *chromsize = 0;

// check no. of chromosomes, and terminate if in error
bArchFile >> paramname >> nchromosomes;
if (paramname == "NChromosomes") {
	if (nchromosomes < 1) {
		BatchError(filetype,-999,11,"NChromosomes"); errors++;
		return -111;
	}
}
else {
	ArchFormatError();
	return -111;
}
chromsize = new int[nchromosomes];
for (int i = 0; i < nchromosomes; i++) chromsize[i] = 0;

// check no. of loci on each chromosome, and terminate if in error
bArchFile >> paramname;
if (paramname != "NLoci") formatError = true;
int locerrors = 0;
for (int i = 0; i < nchromosomes; i++) {
	nloci = -999;
	bArchFile >> nloci;
	if (nloci < 1) locerrors++; else chromsize[i] = nloci;
}
if (locerrors) {
	BatchError(filetype,-999,11,"NLoci");
	return -111;
}
//if (formatError) batchlog << "formatError is TRUE" << endl;
//else batchlog << "formatError is FALSE" << endl;

// check unspecified no. of traits
fileNtraits = 0;
int traitnum,prevtrait,chrom,locus;
traitnum = prevtrait = -1;
bool traitError = false;
bool lociError = false;
bool chromError = false;
bool locusError = false;
paramname = "XXXyyyZZZ";
//batchlog << "paramname=" << paramname << endl;
bArchFile >> paramname;
//batchlog << "paramname=" << paramname << endl;
while (paramname != "XXXyyyZZZ") {
	bArchFile >> traitnum;
//	batchlog << "traitnum=" << traitnum << endl;
	if (paramname != "Trait") formatError = true;
	if (traitnum == (prevtrait+1)) prevtrait = traitnum;
	else traitError = true;
	bArchFile >> paramname >> nloci;
//	batchlog << "paramname=" << paramname << " nloci=" << nloci << endl;
	if (paramname != "NLoci") formatError = true;
	if (nloci < 1) lociError = true;
	for (int i = 0; i < nloci; i++) {
		chrom = locus = -999999;
		bArchFile >> chrom >> locus;
//		batchlog << "chrom=" << chrom << " locus=" << locus << endl;
		if (chrom == -999999 || locus == -999999) {
			BatchError(filetype,-999,0," "); errors++;
			batchlog << "Too few loci listed for trait " << traitnum << endl;
		}
		else {
			if (chrom >= 0 && chrom < nchromosomes) {
//				batchlog << "chromsize[" << chrom << "]=" << chromsize[chrom] << endl;
				if (locus < 0 || locus >= chromsize[chrom]) locusError = true;
			}
			else chromError = true;
		}
	}
	fileNtraits++;
	paramname = "XXXyyyZZZ";
	bArchFile >> paramname;
//	batchlog << "paramname=" << paramname << " (end of loop)" << endl;
}
//batchlog << "paramname=" << paramname << " (after loop)" << endl;

if (traitError) {
	BatchError(filetype,-999,0," "); errors++;
	batchlog << "Traits must be sequentially numbered starting at 0 " << endl;
}
if (lociError) {
	BatchError(filetype,-999,11,"Trait NLoci"); errors++;
}
if (chromError) {
	BatchError(filetype,-999,0," "); errors++;
	batchlog << "Chromosome no. must be from 0 to " << (nchromosomes-1) << endl;
}
if (locusError) {
	BatchError(filetype,-999,0," "); errors++;
	batchlog << "Locus no. must not exceed no. of loci on specified chromosome " << endl;
}
//if (formatError) batchlog << "formatError is TRUE" << endl;
//else batchlog << "formatError is FALSE" << endl;

if (formatError || errors > 0) { // terminate batch error checking
	if (formatError) ArchFormatError();
	return -111;
}

// final read should hit EOF

if (!bArchFile.eof()) {
	EOFerror(filetype);
	errors++;
}

if (chromsize != 0) delete[] chromsize;

return errors;

}

#if RS_CONTAIN

//---------------------------------------------------------------------------
int ParseManageFile(string indir)
{
string header,colheader;
int simul,err;
int inint,cullNstages,prevNstages,cullstage,prevstage,cullrate,firstcullrate;
int damagetype,occoption;
string filename,ftype,fname;
float infloat;
bool checkfile;
int errors = 0; 
//int cullerrors = 0;
int simuls = 0;
//vector <string> archfiles;
string filetype = "ManagementFile";

// Parse header line;
bManageFile >> header; if (header != "Simulation" ) errors++;
bManageFile >> header; if (header != "Method" ) errors++;
bManageFile >> header; if (header != "Timing" ) errors++;
bManageFile >> header; if (header != "DistDecay" ) errors++;
bManageFile >> header; if (header != "MaxPatches" ) errors++;
//bManageFile >> header; if (header != "ThresholdPop" ) errors++;
bManageFile >> header; if (header != "ThresholdDens" ) errors++;
bManageFile >> header; if (header != "CountCV" ) errors++;
if (stagestruct) {
	bManageFile >> header; if (header != "CullNstages" ) errors++;	
	bManageFile >> header; if (header != "CullStage" ) errors++;	
}
bManageFile >> header; if (header != "CullRate" ) errors++;
bManageFile >> header; if (header != "CullMaxRate" ) errors++;
bManageFile >> header; if (header != "CullAlpha" ) errors++;
bManageFile >> header; if (header != "CullBeta" ) errors++;
//bManageFile >> header; if (header != "EdgeBias" ) errors++;
//if (stagestruct) {
//	for (i = 0; i < stages; i++) {
//		colheader = "CullStage" + Int2Str(i);
//		bManageFile >> header; if (header != colheader ) cullerrors++;
//	}
//}
bManageFile >> header; if (header != "DamageTiming" ) errors++;
bManageFile >> header; if (header != "DamageType" ) errors++;
bManageFile >> header; if (header != "OccOption" ) errors++;
bManageFile >> header; if (header != "DamageStage" ) errors++;
bManageFile >> header; if (header != "AlphaOccupancy" ) errors++;
bManageFile >> header; if (header != "BetaOccupancy" ) errors++;
bManageFile >> header; if (header != "AlphaTraversal" ) errors++;
bManageFile >> header; if (header != "BetaTraversal" ) errors++;
// report any errors in headers, and if so, terminate validation
//if (errors > 0 || cullerrors > 0) {
//	FormatError(filetype,errors+cullerrors);
//	if (cullerrors > 0) BatchError(filetype,-999,555,"CullStage");
//	return -111;
//}
if (errors > 0) {
	FormatError(filetype,errors);
	return -111;
}

// Parse data lines
int line = 1;
int nstages = 0;
simCheck current,prev;
simul = -98765;
prev.simul = -999;
prev.simlines = prev.reqdsimlines = 0;
prevstage = prevNstages = -1;
string msgsame = " must be the same for all lines of the Simulation";

bManageFile >> simul;
// first simulation number must match first one in parameterFile
if (simul != firstsimul) {
	BatchError(filetype,line,111,"Simulation"); errors++;
}
while (simul != -98765) {

	if (stagestruct) 
		current = CheckStageSex(filetype,line,simul,prev,0,0,0,0,0,false,false);
	else 
		current = CheckStageSex(filetype,line,simul,prev,0,0,0,0,0,true,false);
	if (current.newsimul) {
		if (stagestruct && line > 1) {
			// check that no. of stages for previous simulation was correct
			if (nstages != cullNstages) {
				BatchError(filetype,(line-1),0," "); errors++;  
				batchlog << "No. of stages must match CullNstages" << endl;				
			}
		} 
		simuls++; prevstage = prevNstages = -1; nstages = 0; 
	}
	errors += current.errors;
	prev = current;

	// validate parameters

	bManageFile >> inint;
	if (current.newsimul && (inint < 0 || inint > 7)) {
		BatchError(filetype,line,7,"Method"); errors++;
	}
	bManageFile >> inint;
	if (current.newsimul && (inint < 0 || inint > 1)) {
		BatchError(filetype,line,1,"Timing"); errors++;
	}
	bManageFile >> infloat;
	if (current.newsimul && infloat < 0.0) {
		BatchError(filetype,line,10,"DistDecay"); errors++;
	}
	bManageFile >> inint;
	if (current.newsimul && inint < 2) {
		BatchError(filetype,line,12,"MaxPatches"); errors++;
	}
//	bManageFile >> inint;
//	if (inint < 1) {
//		BatchError(filetype,line,11,"ThresholdPop"); errors++;
//	}
	bManageFile >> infloat;
	if (current.newsimul && infloat <= 0.0) {
		BatchError(filetype,line,10,"ThresholdDens"); errors++;
	}
	bManageFile >> infloat;
	if (current.newsimul && infloat < 0.0) {
		BatchError(filetype,line,19,"CountCV"); errors++;
	}
	if (stagestruct) {
		bManageFile >> cullNstages;
		if (cullNstages < 1 || cullNstages > stages) {
			BatchError(filetype,line,0," "); errors++;  
			batchlog << "CullNstages must be from 1 to " << stages << endl;
		}
		if (!current.newsimul && cullNstages != prevNstages) {
			BatchError(filetype,line,0," "); errors++;  
			batchlog << "CullNstages" << msgsame << endl;
		}
		bManageFile >> cullstage;
		if (cullstage < 0 || cullstage >= stages) {
			BatchError(filetype,line,0," "); errors++;  
			batchlog << "CullStage must be from 0 to " << (stages-1) << endl;
		}
		if (cullstage <= prevstage) {
			BatchError(filetype,line,0," "); errors++;  
			batchlog << "CullStage must be greater than previous CullStage" << endl;
		}
		prevstage = cullstage; prevNstages = cullNstages; nstages++;
	}
	bManageFile >> cullrate;
	if (stagestruct) {
		if (cullrate < 0 || cullrate > 2) {
			BatchError(filetype,line,2,"CullRate"); errors++;
		}		
		if (current.newsimul) {
			firstcullrate = cullrate;
		}
		else {
			if (cullrate != firstcullrate) {
				BatchError(filetype,line,0," "); errors++;  
				batchlog << "CullRate" << msgsame << endl;
			}
		}		
	}
	else {
		if (cullrate < 0 || cullrate > 1) {
			BatchError(filetype,line,1,"CullRate"); errors++;
		}
	}
	bManageFile >> infloat;
	if (current.newsimul || cullrate == 2) {
		if (infloat < 0.0 || infloat > 1.0) {
			BatchError(filetype,line,20,"CullMaxRate"); errors++;
		}
	}
	bManageFile >> infloat >> infloat; // CullAlpha and CullBeta - no conditions
//	bManageFile >> inint;
//	if (inint < 0 || inint > 1) {
//		BatchError(filetype,line,1,"EdgeBias"); errors++;
//	}
	
//	if (stagestruct) {
////		int cullstage;
//		for (i = 0; i < stages; i++) {
//			bManageFile >> inint;
//			if (inint < 0 || inint > 1) {
//				colheader = "CullStage" + Int2Str(i);
//				BatchError(filetype,line,1,colheader); errors++;
//			}
//		}
//	}
	bManageFile >> inint;
	if (current.newsimul && (inint < 0 || inint > 1)) {
		BatchError(filetype,line,1,"DamageTiming"); errors++;
	}
	bManageFile >> damagetype;
	if (current.newsimul && (damagetype < 0 || damagetype > 2)) {
		BatchError(filetype,line,2,"DamageType"); errors++;
	}
	bManageFile >> occoption;
	if (stagestruct) {
		if (current.newsimul && (occoption < 0 || occoption > 3)) {
			BatchError(filetype,line,3,"OccOption"); errors++;
		}
	}
	else {
		if (occoption < 0 || occoption > 1) {
			BatchError(filetype,line,1,"OccOption"); errors++;
		}
	}
	bManageFile >> inint;
	if (stagestruct && occoption == 3) {
		if (inint < 0 || inint >= stages) {
			BatchError(filetype,line,0," "); errors++;  
			batchlog << "DamageStage must be from 0 to " << (stages-1) << endl;
		}
	}
	bManageFile >> infloat;
	if (damagetype == 2 && infloat < 0.0) {
		BatchError(filetype,line,10,"AlphaOccupancy"); errors++;
	}
	bManageFile >> infloat;
	if (damagetype == 1 && infloat < 0.0) {
		BatchError(filetype,line,10,"BetaOccupancy"); errors++;
	}
	bManageFile >> infloat;
	if (transfer == 1 && damagetype == 2 && infloat < 0.0) {
		BatchError(filetype,line,10,"AlphaTraversal"); errors++;
	}
	bManageFile >> infloat;
	if (transfer == 1 && damagetype == 1 && infloat < 0.0) {
		BatchError(filetype,line,10,"BetaTraversal"); errors++;
	}

	// read next simulation
	line++;
	simul = -98765;
	bManageFile >> simul;
	if (bManageFile.eof()) simul = -98765;
} // end of while loop
// check for correct number of lines for previous simulation
if (!stagestruct && current.simlines != current.reqdsimlines) {
	BatchError(filetype,line,0," "); errors++;
	batchlog << msgnlines << current.simul
		<< msgshldbe << current.reqdsimlines << endl;
}
if (stagestruct) {
	// check that no. of stages for final simulation was correct
	if (nstages != cullNstages) {
		BatchError(filetype,(line-1),0," "); errors++;  
		batchlog << "No. of stages must match CullNstages" << endl;				
	}
} 
if (!bManageFile.eof()) {
	EOFerror(filetype);
	errors++;
}

if (errors > 0) return -111;
else return simuls;

}

#endif // RS_CONTAIN 

//---------------------------------------------------------------------------
int ParseInitFile(string indir)
{
string header,colheader;
int i,simul;
int seedtype,freetype,sptype,initdens,indscell = 0,minX,maxX,minY,maxY;
int nCells,nSpCells,initAge;
int initFreezeYear,restrictRows,restrictFreq,finalFreezeYear;
float inds_per_ha;

int errors = 0; int propnerrors = 0;
int simuls = 0;
string filetype = "InitialisationFile";

// Parse header line;
bInitFile >> header; if (header != "Simulation" ) errors++;
bInitFile >> header; if (header != "SeedType" ) errors++;
bInitFile >> header; if (header != "FreeType" ) errors++;
bInitFile >> header; if (header != "SpType" ) errors++;
bInitFile >> header; if (header != "InitDens" ) errors++;
bInitFile >> header;
if (patchmodel) { if (header != "IndsHa" ) errors++; }
else 						{ if (header != "IndsCell" ) errors++; }
bInitFile >> header; if (header != "minX" ) errors++;
bInitFile >> header; if (header != "maxX" ) errors++;
bInitFile >> header; if (header != "minY" ) errors++;
bInitFile >> header; if (header != "maxY" ) errors++;
bInitFile >> header; if (header != "NCells" ) errors++;
bInitFile >> header; if (header != "NSpCells" ) errors++;
bInitFile >> header; if (header != "InitFreezeYear" ) errors++;
bInitFile >> header; if (header != "RestrictRows" ) errors++;
bInitFile >> header; if (header != "RestrictFreq" ) errors++;
bInitFile >> header; if (header != "FinalFreezeYear" ) errors++;
bInitFile >> header; if (header != "InitIndsFile" ) errors++;
if (stagestruct) {
	bInitFile >> header; if (header != "InitAge" ) errors++;
	for (i = 1; i < stages; i++) {
		colheader = "PropStage" + Int2Str(i);
		bInitFile >> header; if (header != colheader ) propnerrors++;
	}
}
// report any errors in headers, and if so, terminate validation
if (errors > 0 || propnerrors > 0) {
	FormatError(filetype,errors+propnerrors);
	if (propnerrors > 0) BatchError(filetype,-999,444,"PropStage");
	return -111;
}

// Parse data lines
int line = 1;
int err;
simCheck current{}, prev{};
bool checkfile;
string filename,ftype2,fname;
vector <string> indsfiles;
ftype2 = "InitIndsFile";
simul = -98765;
prev.simul = -999;
prev.simlines = prev.reqdsimlines = 0;
bInitFile >> simul;
// first simulation number must match first one in parameterFile
if (simul != firstsimul) {
	BatchError(filetype,line,111,"Simulation"); errors++;
}
while (simul != -98765) {
	current = CheckStageSex(filetype,line,simul,prev,0,0,0,0,0,true,false);
	if (current.newsimul) simuls++;
	errors += current.errors;
	prev = current;

	bInitFile >> seedtype >> freetype >> sptype >> initdens >> inds_per_ha;
	if (!patchmodel) indscell = (int)inds_per_ha;
	if (seedtype < 0 || seedtype > 2) {
		BatchError(filetype,line,2,"SeedType"); errors++;
	}
	if (landtype == 9 && seedtype != 0) {
		BatchError(filetype,line,0," "); errors++;
		batchlog << "SeedType must be 0 for an artificial landscape"
			<< endl;
	}
	if (!speciesdist && seedtype == 1) {
		BatchError(filetype,line,0," "); errors++;
		batchlog << "SeedType may not be 1 if there is no initial species distribution map"
			<< endl;
	}
	if (seedtype == 0) {
		if (freetype < 0 || freetype > 1) {
			BatchError(filetype,line,1,"FreeType"); errors++;
		}
	}
	if (seedtype == 1) {
		if (sptype < 0 || sptype > 1) {
			BatchError(filetype,line,1,"SpType"); errors++;
		}
	}
	if (initdens < 0 || initdens > 2) {
		BatchError(filetype,line,2,"initDens"); errors++;
	}
	if (seedtype < 2) {
		if (initdens == 2) { // specified density
			if (patchmodel) {
				 if (inds_per_ha <= 0.0) {
					BatchError(filetype,line,10,"IndsHa"); errors++;
				 }
			}
			else {
				if (indscell < 1) {
					BatchError(filetype,line,11,"IndsCell"); errors++;
				}
			}
		}
	}

	bInitFile >> minX >> maxX >> minY >> maxY >> nCells >> nSpCells;
	if (seedtype == 0) {
//		if (minX < 0) {
//			BatchError(filetype,line,19,"minX"); errors++;
//		}
		if (maxX < minX) {
			BatchError(filetype,line,2,"maxX","minX"); errors++;
		}
//		if (minY < 0) {
//			BatchError(filetype,line,19,"minY"); errors++;
//		}
		if (maxY < minY) {
			BatchError(filetype,line,2,"maxY","minY"); errors++;
		}
	}
	if (seedtype == 0 && freetype == 0) {
		if (nCells < 1) {
			BatchError(filetype,line,11,"NCells"); errors++;
		}
		int range_cells;
		range_cells = (maxX-minX)*(maxY-minY);
		if (nCells > range_cells) {
			BatchError(filetype,line,0," "); errors++;
			batchlog << "NCells may not be greater than the area specified (i.e. "
				<< range_cells << " cells)" << endl;
		}
	}
	if (seedtype == 1 && sptype == 1 && nSpCells < 1) {
		BatchError(filetype,line,11,"NSpCells"); errors++;
	}

	bInitFile >> initFreezeYear >> restrictRows >> restrictFreq >> finalFreezeYear;
	if (seedtype == 0) {
		if (initFreezeYear < 0) {
			BatchError(filetype,line,19,"InitFreezeYear"); errors++;
		}
		if (restrictRows < 0) {
			BatchError(filetype,line,19,"RestrictRows"); errors++;
		}
		if (restrictRows > 0 && restrictFreq <= 0) {
			BatchError(filetype,line,10,"RestrictFreq"); errors++;
		}
		if (finalFreezeYear < 0) {
			BatchError(filetype,line,19,"FinalFreezeYear"); errors++;
		}
		else {
			if (finalFreezeYear > 0 && finalFreezeYear <= initFreezeYear) {
				BatchError(filetype,line,1,"FinalFreezeYear","InitFreezeYear"); errors++;
			}
		}
	}

	bInitFile >> filename;
	if (filename == "NULL") {
		if (seedtype == 2) {
			BatchError(filetype,line,0," "); errors++;
			batchlog << ftype2 << " is compulsory for SeedType 2" << endl;
		}
	}
	else {
		if (seedtype == 2) {
			checkfile = true;
			for (i = 0; i < (int)indsfiles.size(); i++) {
				if (filename == indsfiles[i]) { // file has already been checked
					checkfile = false;
				}
			}
			if (checkfile) {
				fname = indir + filename;
				batchlog << "Checking " << ftype2 << " " << fname << endl;
				bInitIndsFile.open(fname.c_str());
				if (bInitIndsFile.is_open()) {
					err = ParseInitIndsFile();
					if (err == 0) FileHeadersOK(ftype2); else errors++;
					bInitIndsFile.close();
				}
				else {
					OpenError(ftype2,fname); errors++;
				}
				if (bInitIndsFile.is_open()) bInitIndsFile.close();
				bInitIndsFile.clear();
				indsfiles.push_back(filename);
			}
		}
		else {
			BatchError(filetype,line,0," "); errors++;
			batchlog << ftype2 << " must be NULL for SeedType "
				<< seedtype << endl;
		}
	}

	if (stagestruct) {
		bInitFile >> initAge;
		if (seedtype != 2 && (initAge < 0 || initAge > 2)) {
			BatchError(filetype,line,2,"initAge"); errors++;
		}
		float propstage;
		float cumprop = 0.0;
		for (i = 1; i < stages; i++) {
			bInitFile >> propstage;
			cumprop += propstage;
			if (seedtype != 2 && (propstage < 0.0 || propstage > 1.0)) {
				colheader = "PropStage" + Int2Str(i);
				BatchError(filetype,line,20,colheader); errors++;
			}
		}
		if (seedtype != 2 && (cumprop < 0.99999 || cumprop > 1.00001)) {
			BatchError(filetype,line,0," "); errors++;
			batchlog << "Initial proportions must sum to 1.0" << endl;
		}
	}

	// read next simulation
	line++;
	simul = -98765;
	bInitFile >> simul;
	if (bInitFile.eof()) simul = -98765;
} // end of while loop
// check for correct number of lines for previous simulation
if (current.simlines != current.reqdsimlines) {
	BatchError(filetype,line,0," "); errors++;
	batchlog << msgnlines << current.simul
		<< msgshldbe << current.reqdsimlines << endl;
}
if (!bInitFile.eof()) {
	EOFerror(filetype);
	errors++;
}

if (errors > 0) return -111;
else return simuls;

}

//---------------------------------------------------------------------------
int ParseInitIndsFile(void) {
string header;
int year,species,patchID,x,y,ninds,sex,age,stage,prevyear;

int errors = 0;
string filetype = "InitIndsFile";

// Parse header line
bInitIndsFile >> header; if (header != "Year" ) errors++;
bInitIndsFile >> header; if (header != "Species" ) errors++;
if (patchmodel) {
	bInitIndsFile >> header; if (header != "PatchID" ) errors++;
}
else {
	bInitIndsFile >> header; if (header != "X" ) errors++;
	bInitIndsFile >> header; if (header != "Y" ) errors++;
}
bInitIndsFile >> header; if (header != "Ninds" ) errors++;
if (reproductn > 0) {
	bInitIndsFile >> header; if (header != "Sex" ) errors++;
}
if (stagestruct) {
	bInitIndsFile >> header; if (header != "Age" ) errors++;
	bInitIndsFile >> header; if (header != "Stage" ) errors++;
}

// Report any errors in headers, and if so, terminate validation
if (errors > 0) {
	FormatError(filetype,errors);
	return -111;
}

// Parse data lines
int line = 1;
string filename,ftype2,fname;
year = prevyear = -98765;
bInitIndsFile >> year;
while (year != -98765) {
	if (year < 0) {
		BatchError(filetype,line,19,"Year"); errors++;
	}
	else {
		if (year < prevyear) {
			BatchError(filetype,line,2,"Year","previous Year"); errors++;
		}
	}
	prevyear = year;
	bInitIndsFile >> species;
	if (species != 0) {
		BatchError(filetype,line,0," "); errors++;
		batchlog << "Species must be 0" << endl;
	}
	if (patchmodel) {
		bInitIndsFile >> patchID;
		if (patchID < 1) {
			BatchError(filetype,line,11,"PatchID"); errors++;
		}
	}
	else {
		bInitIndsFile >> x >> y ;
		if (x < 0 || y < 0) {
			BatchError(filetype,line,19,"X and Y"); errors++;
		}
	}
	bInitIndsFile >> ninds;
	if (ninds < 1) {
		BatchError(filetype,line,11,"Ninds"); errors++;
	}
	if (reproductn > 0) {
		bInitIndsFile >> sex;
		if (sex < 0 || sex > 1) {
			BatchError(filetype,line,1,"Sex"); errors++;
		}
	}
	if (stagestruct) {
		bInitIndsFile >> age >> stage;
		if (age < 1) {
			BatchError(filetype,line,11,"Age"); errors++;
		}
		if (stage < 1) {
			BatchError(filetype,line,11,"Stage"); errors++;
		}
		if (stage >= stages) {
			BatchError(filetype,line,4,"Stage","no. of stages"); errors++;
		}
	}
	line++;
	year = -98765;
	bInitIndsFile >> year;
	if (bInitIndsFile.eof()) year = -98765;
} // end of while loop
if (!bInitIndsFile.eof()) {
	EOFerror(filetype);
	errors++;
}

return errors;
}

#if VIRTUALECOLOGIST

int ParseVirtEcolFile(string indir)
{
string header,colheader;
int i,simul,err;
int patchMethod,rowsFront,maxNPatches,maxPatchNum,minIndsPatch,maxIndsPatch,stageMethod;
int minX,maxX,minY,maxY,landGen,landGenStart,landGenFreq,writeSamples,sampleLoci;
string filename,ftype,fname;
bool checkfile;
int errors = 0;
int simuls = 0;
vector <string> samplefiles;
vector <string> patchfiles;
string filetype = "VirtualEcolFile";

// Parse header line;
bVirtEcolFile >> header; if (header != "Simulation" ) errors++;
bVirtEcolFile >> header; if (header != "PatchMethod" ) errors++;
bVirtEcolFile >> header; if (header != "RowsFront" ) errors++;
bVirtEcolFile >> header; if (header != "MaxNPatches" ) errors++;
bVirtEcolFile >> header; if (header != "MaxPatchNum" ) errors++;
bVirtEcolFile >> header; if (header != "PatchFile" ) errors++;
bVirtEcolFile >> header; if (header != "MinIndsPatch" ) errors++;
bVirtEcolFile >> header; if (header != "MaxIndsPatch" ) errors++;
bVirtEcolFile >> header; if (header != "StageMethod" ) errors++;
bVirtEcolFile >> header; if (header != "MinX" ) errors++;
bVirtEcolFile >> header; if (header != "MaxX" ) errors++;
bVirtEcolFile >> header; if (header != "MinY" ) errors++;
bVirtEcolFile >> header; if (header != "MaxY" ) errors++;
bVirtEcolFile >> header; if (header != "LandGen" ) errors++;
bVirtEcolFile >> header; if (header != "LandGenStart" ) errors++;
bVirtEcolFile >> header; if (header != "LandGenFreq" ) errors++;
bVirtEcolFile >> header; if (header != "WriteSamples" ) errors++;
bVirtEcolFile >> header; if (header != "SampleLoci" ) errors++;
bVirtEcolFile >> header; if (header != "SampleFile" ) errors++;
if (errors > 0) {
	FormatError(filetype,errors);
	return -111;
}

// Parse data lines
int line = 1;
simCheck current,prev;
simul = -98765;
prev.simul = -999;
prev.simlines = prev.reqdsimlines = 0;
bVirtEcolFile >> simul;
// first simulation number must match first one in parameterFile
if (simul != firstsimul) {
	BatchError(filetype,line,111,"Simulation"); errors++;
}
while (simul != -98765) {
	// read and validate columns relating to stage and sex-dependency (NB none here)
	current = CheckStageSex(filetype,line,simul,prev,0,0,0,0,0,false,false);
	if (current.newsimul) simuls++;
	errors += current.errors;
	prev = current;

	// validate parameters
	bVirtEcolFile >> patchMethod >> rowsFront >> maxNPatches >> maxPatchNum; 

	if (patchMethod < 0 || patchMethod > 3) {
		BatchError(filetype,line,3,"PatchMethod"); errors++;
	}
	if (patchMethod == 2) {
		if (rowsFront < 1) {
			BatchError(filetype,line,11,"RowsFront"); errors++;
		}
	}
	if (maxNPatches < 2) {
		BatchError(filetype,line,12,"MaxNPatches"); errors++;
	}
	if (maxPatchNum < 0) {
		BatchError(filetype,line,19,"MaxPatchNum"); errors++;
	}

	// patch sample file - optional
	bVirtEcolFile >> filename;
	ftype = "PatchFile";
	if (filename == "NULL") {
		if (patchMethod == 3) {
			BatchError(filetype,line,0," ");
			batchlog << ftype << " is compulsory for PatchMethod 3" << endl;
			errors++;
		}
	}
	else {
		if (patchMethod != 3) {
			BatchError(filetype,line,0," ");
			batchlog << ftype << " must be NULL unless PatchMethod is 3" << endl;
			errors++;
		}
		else { // check patch sample file
			checkfile = true;
			for (i = 0; i < (int)patchfiles.size(); i++) {
				if (filename == patchfiles[i]) checkfile = false; // file has already been checked
			}
			if (checkfile) {
				fname = indir + filename;
				batchlog << "Checking " << ftype << " " << fname << endl;
				bPatchFile.open(fname.c_str());
				if (bPatchFile.is_open()) {
					err = ParsePatchFile();
					if (err == 0) FileHeadersOK(ftype); else errors++;
					bPatchFile.close();
				}
				else {
					OpenError(ftype,fname); errors++;
				}
				if (bPatchFile.is_open()) bPatchFile.close();
				bPatchFile.clear();
			}
			patchfiles.push_back(filename);
		}
	}

	bVirtEcolFile >> minIndsPatch >> maxIndsPatch >> stageMethod >> minX >> maxX >> minY >> maxY
		>> landGen >> landGenStart >> landGenFreq >> writeSamples >> sampleLoci;
	if (patchMethod == 2) {
		if (minIndsPatch < 2) {
			BatchError(filetype,line,12,"MinIndsPatch"); errors++;
		}
	}
	if (maxIndsPatch < 2) {
		BatchError(filetype,line,12,"MaxIndsPatch"); errors++;
	}
	if (stagestruct) {
		if (stageMethod < 0 || stageMethod > 2) {
			BatchError(filetype,line,2,"StageMethod"); errors++;
		}
	}
	if (patchMethod == 1) {
		if (minX < 0) {
			BatchError(filetype,line,19,"MinX"); errors++;
		}
		if (maxX < minX) {
			BatchError(filetype,line,2,"MaxX","MinX"); errors++;
		}
		if (minY < 0) {
			BatchError(filetype,line,19,"MinY"); errors++;
		}
		if (maxY < minY) {
			BatchError(filetype,line,2,"MaxY","MinY"); errors++;
		}
	}
	if (landGen < 0 || landGen > 1) {
		BatchError(filetype,line,1,"LandGen"); errors++;
	}
	if (landGen == 1) {
		if (landGenStart < 0) {
			BatchError(filetype,line,19,"LandGenStart"); errors++;
		}
		if (landGenFreq < 0) {
			BatchError(filetype,line,19,"LandGenFreq"); errors++;
		}
		if (writeSamples < 0 || writeSamples > 1) {
			BatchError(filetype,line,1,"WriteSamples"); errors++;
		}
		if (sampleLoci < 0 || sampleLoci > 1) {
			BatchError(filetype,line,1,"SampleLoci"); errors++;
		}
	}

	// sample loci file - optional
	bVirtEcolFile >> filename;
	ftype = "SampleFile";
	if (landGen == 1) {
		if (filename == "NULL") {
			if (sampleLoci != 0) {
				BatchError(filetype,line,0," ");
				batchlog << ftype << " is compulsory unless SampleLoci is 0" << endl;
				errors++;
			}
		}
		else {
			if (sampleLoci != 1) {
				BatchError(filetype,line,0," ");
				batchlog << ftype << " must be NULL if SampleLoci is 0" << endl;
				errors++;
			}
			else { // check sample loci file
				checkfile = true;
				for (i = 0; i < (int)samplefiles.size(); i++) {
					if (filename == samplefiles[i]) checkfile = false; // file has already been checked
				}
				if (checkfile) {
					fname = indir + filename;
					batchlog << "Checking " << ftype << " " << fname << endl;
					bSampleFile.open(fname.c_str());
					if (bSampleFile.is_open()) {
						err = ParseSampleFile();
						if (err == 0) FileHeadersOK(ftype); else errors++;
						bSampleFile.close();
					}
					else {
						OpenError(ftype,fname); errors++;
					}
					if (bSampleFile.is_open()) bSampleFile.close();
					bSampleFile.clear();
				}
				samplefiles.push_back(filename);
			}
		}
	}

	// read next simulation
	line++;
	simul = -98765;
	bVirtEcolFile >> simul;
	if (bVirtEcolFile.eof()) simul = -98765;
} // end of while loop
// check for correct number of lines for previous simulation
if (current.simlines != current.reqdsimlines) {
	BatchError(filetype,line,0," "); errors++;
	batchlog << msgnlines << current.simul
		<< msgshldbe << current.reqdsimlines << endl;
}
if (!bVirtEcolFile.eof()) {
	EOFerror(filetype);
	errors++;
}

if (errors > 0) return -111;
else return simuls;

}

//---------------------------------------------------------------------------
int ParseSampleFile(void)
{
string paramname;
int chromosome,locus;
int nsamples;
//float infloat;
int errors = 0;
bool formatError = false;
string filetype = "SampleFile";

// check no. of samples, and terminate if in error
bSampleFile >> paramname >> nsamples;
if (paramname == "NSamples") {
	if (nsamples < 2) {
		BatchError(filetype,-999,11,"NSamples"); errors++;
		return -111;
	}
}
else {
	SampleFormatError();
	return -111;
}

// check no. of loci on each chromosome, and terminate if in error
bSampleFile >> paramname;
if (paramname != "Loci") formatError = true;
int nchromosomes = pSpecies->getNChromosomes();
//#if RSDEBUG
//DEBUGLOG << "ParseSampleFile(): nchromosomes=" << nchromosomes
//	<< endl;
//#endif
int nloci;
int locerrors = 0;
for (int i = 0; i < nsamples; i++) {
	chromosome = locus = -999;
	bSampleFile >> chromosome >> locus;
	if (chromosome < 0 || locus < 0) locerrors++;
//	if (chromosome < 0 || chromosome >= nchromosomes || locus < 0) locerrors++;
//	else {
//		// check specified locus exists
//		nloci = pSpecies->getNLoci(chromosome);
//		if (locus >= nloci) locerrors++;
//	}
}
if (locerrors) {
	BatchError(filetype,-999,19,"Chromosome and locus");
//	BatchError(filetype,2,0," "); errors++;
//	batchlog << "Invalid chromosome and locus combination - not found on genome " << endl;
	return -111;
}
//if (formatError) batchlog << "formatError is TRUE" << endl;
//else batchlog << "formatError is FALSE" << endl;

if (formatError || errors > 0) { // terminate batch error checking
	if (formatError) SampleFormatError();
	return -111;
}

// final read should hit EOF
bSampleFile >> paramname;

if (!bSampleFile.eof()) {
	EOFerror(filetype);
	errors++;
}

return errors;

}

//---------------------------------------------------------------------------
int ParsePatchFile(void)
{
int patchnum;
int errors = 0;
bool formatError = false;
string filetype = "PatchFile";

patchnum = -98765;
bPatchFile >> patchnum;
while (patchnum != -98765) {
	if (patchnum < 1) errors++;
	patchnum = -98765;
	bPatchFile >> patchnum;
}

if (errors) {    
	BatchError(filetype,-999,10,"Patch number");
	return -111;
}

return errors;

}

#endif // VIRTUALECOLOGIST

#if RS_ABC

//---------------------------------------------------------------------------
int ParseABCParamsFile(void)
{
string header;
string nm;
int i,inint,type,priortype,stage,sex;
int prevparam;
float pr0,pr1;
int errors = 0;
string filetype = "ABCParamsFile";

//batchlog << "ParseABCFile(): starting " << endl;
// Parse header line;
bABCpFile >> header; if (header != "ParamNum" ) errors++;
bABCpFile >> header; if (header != "ParamType" ) errors++;
bABCpFile >> header; if (header != "ParamName" ) errors++;
bABCpFile >> header; if (header != "Stage" ) errors++;
bABCpFile >> header; if (header != "Sex" ) errors++;
bABCpFile >> header; if (header != "PriorType" ) errors++;
bABCpFile >> header; if (header != "Pr0" ) errors++;
bABCpFile >> header; if (header != "Pr1" ) errors++;
if (errors > 0) {
	FormatError(filetype,errors);
	batchlog << "*** Ensure column headers are correct to continue checking data" << endl;
	return -111;
}

// Parse data lines
int line = 1;
int nparams = 0;
inint = -98765;
bABCpFile >> inint; // first parameter number
if (inint != 0) {
	batchlog << "*** Error in ABCParamsFile - first parameter number must be 0" << endl;
	errors++;
}
else {
	prevparam = inint; nparams++;
}
while (inint != -98765) {
	bABCpFile >> type;
	if (type < 1 || type > 5) { BatchError(filetype,line,55,"ParamType"); errors++; }
	nm = "ZZZxxxZZZxxx";
	bABCpFile >> nm;
	if (nm == "ZZZxxxZZZxxx") {
		BatchError(filetype,line,0," ");
		batchlog << "Invalid parameter name" << endl;
		errors++;
	}
	else {
		// permitted parameter names depend on parameter type
		bool valid = false;
		switch (type) {
		case 1: // demographic
			if (nm == "K1" || nm == "K2" || nm == "K3" || nm == "K4" || nm == "K5") valid = true;
			if (nm == "K6" || nm == "K7" || nm == "K8" || nm == "K9") valid = true;
			if (nm == "Rmax") valid = true;
			if (nm == "Fec" || nm == "Dev" || nm == "Surv") valid = true;
			break;
		case 2: // genome
			break;
		case 3: // emigration
			if (nm == "D0" || nm == "alpha" || nm == "beta") valid = true;
			break;
		case 4: // transfer
			// kernel
			if (nm == "meanDistI" || nm == "meanDistII" || nm == "ProbKernelI") valid = true;
			// SMS
			if (nm == "DP" || nm == "MemSize" || nm == "GB" || nm == "AlphaDB" || nm == "BetaDB") valid = true;
			if (nm == "SMconst") valid = true;
			if (nm == "MortHab1" || nm == "MortHab2" || nm == "MortHab3") valid = true;
			if (nm == "MortHab4" || nm == "MortHab5" || nm == "MortHab6") valid = true;
			if (nm == "MortHab7" || nm == "MortHab8" || nm == "MortHab9") valid = true;
			if (nm == "CostHab1" || nm == "CostHab2" || nm == "CostHab3") valid = true;
			if (nm == "CostHab4" || nm == "CostHab5" || nm == "CostHab6") valid = true;
			if (nm == "CostHab7" || nm == "CostHab8" || nm == "CostHab9") valid = true;
			// CRW
			if (nm == "SL" || nm == "Rho") valid = true;
			break;
		case 5: // settlement
			if (nm == "S0" || nm == "AlphaS" || nm == "BetaS") valid = true;
			break;
		default:
			;
		}
		if (!valid) {
			BatchError(filetype,line,0," ");
			batchlog << "Invalid parameter name for parameter type " << type << endl;
			errors++;
		}
	}
	bABCpFile >> stage >> sex;
	bABCpFile >> priortype;
	if (priortype < 1 || priortype > 6) { BatchError(filetype,line,66,"PriorType"); errors++; }
	bABCpFile >> pr0 >> pr1;

	switch (priortype) {
	case 1: // uniform
		if (pr0 >= pr1) {
			BatchError(filetype,line,1,"Pr1","Pr0"); errors++;
		}
		break;
	case 2: // normal
		if (pr1 <= 0.0) { // s.d.
			BatchError(filetype,line,10,"Pr1"); errors++;
		}
		break;
	case 3: // log-normal
		if (pr0 <= 0.0) { // mean
			BatchError(filetype,line,10,"Pr0"); errors++;
		}
		if (pr1 <= 1.0) { // s.d.
			BatchError(filetype,line,21,"Pr1"); errors++;
		}
		break;
	case 4: // beta
	case 5: // gamma
		if (pr0 <= 0.0) {
			BatchError(filetype,line,10,"Pr0"); errors++;
		}
		if (pr1 <= 0.0) {
			BatchError(filetype,line,10,"Pr1"); errors++;
		}
		break;
	case 6: // uniform integer
		if (pr0 >= pr1) {
			BatchError(filetype,line,1,"Pr1","Pr0"); errors++;
		}
		break;
	}
	line++;
	// read next simulation number
	inint = -98765;
	bABCpFile >> inint;
	if (bABCpFile.eof()) {
		inint = -98765;
	}
	else { // check for valid parameter number
		if (inint != prevparam+1) {
			BatchError(filetype,line,0," ");
			batchlog << "Parameter numbers must be sequential integers" << endl;
			errors++;
		}
		prevparam = inint; nparams++;
	}
//	batchlog << "ParseABCFile(): First item of next line = " << inint << endl;
} // end of while loop
if (!bABCpFile.eof()) {
	EOFerror(filetype);
	errors++;
}

if (errors > 0) return -111;
else return nparams;
}

//---------------------------------------------------------------------------

int ParseABCObsFile(void)
{
string header;
string nm;
int i,inint,type;
int prevobs;
int year,loc0,loc1,stage,sex;
float value,weight;
int errors = 0;
string filetype = "ABCObsFile";

//batchlog << "ParseABCFile(): starting " << endl;
// Parse header line;
bABCoFile >> header; if (header != "ObsNum" ) errors++;
bABCoFile >> header; if (header != "ObsType" ) errors++;
bABCoFile >> header; if (header != "ObsName" ) errors++;
bABCoFile >> header; if (header != "Year" ) errors++;
bABCoFile >> header; if (header != "Loc0" ) errors++;
bABCoFile >> header; if (header != "Loc1" ) errors++;
bABCoFile >> header; if (header != "Stage" ) errors++;
bABCoFile >> header; if (header != "Sex" ) errors++;
bABCoFile >> header; if (header != "Value" ) errors++;
bABCoFile >> header; if (header != "Weight" ) errors++;
if (errors > 0) {
	FormatError(filetype,errors);
	batchlog << "*** Ensure column headers are correct to continue checking data" << endl;
	return -111;
}

// Parse data lines
int line = 1;
int nparams = 0;
inint = -98765;
bABCoFile >> inint; // first observation number
if (inint != 0) {
	batchlog << "*** Error in ABCObsFile - first observation number must be 0" << endl;
	errors++;
}
else {
	prevobs = inint; nparams++;
}
while (inint != -98765) {
	bABCoFile >> type;
	if (type < 1 || type > 5) { BatchError(filetype,line,55,"ObsType"); errors++; }
	nm = "ZZZxxxZZZxxx";
	bABCoFile >> nm;
	if (nm == "ZZZxxxZZZxxx") {
		BatchError(filetype,line,0," ");
		batchlog << "Invalid observation name" << endl;
		errors++;
	}
	else {
		// permitted observation names depend on observation type
		bool valid = false;
		switch (type) {
		case 1: // range
			if (nm == "NInds" || nm == "NOccupPatches") valid = true;
			break;
		case 2: // population
			if (nm == "NInds" || nm == "Occupied") valid = true;
			break;
		case 3: // individual
			if (nm == "EmigRate" || nm == "DispSuccess" || nm == "DistMean" || nm == "DistSD")
				valid = true;
			break;
		case 4: // connectivity
			if (nm == "NInds") valid = true;
			break;
		case 5: // genome
			if (nm == "FST" || nm == "D") valid = true;
			break;
		default:
			;
		}
		if (!valid) {
			BatchError(filetype,line,0," ");
			batchlog << "Invalid observation name for observation type " << type << endl;
			errors++;
		}
	}
	bABCoFile >> year >> loc0 >> loc1 >> stage >> sex >> value >> weight;
	if (year < 1) {
		BatchError(filetype,line,11,"Year"); errors++;
	}
	switch (type) {
	case 1: // range
		if (nm == "NInds") {
			if (stagestruct) {    
				if (stage < 0 || stage >= stages) {
					BatchError(filetype,line,0," ");
					batchlog << "Stage must be 0 for all stages or equal to a valid stage number in the model" << endl;
					errors++;
				}
			}
			if (sex != 0) {
				BatchError(filetype,line,0," ");
				batchlog << "Sex must be 0 for a range-level observation" << endl;
				errors++;
			}				
		}
		break;
	case 2: // population
		if (patchmodel) {
			if (loc0 < 1) {
				BatchError(filetype,line,0," ");
				batchlog << "Loc0 (PatchID) must be 1 or more for a patch-based model" << endl;
				errors++;
			}
		}
		else {
			if (loc0 < 0 || loc1 < 0) {
				BatchError(filetype,line,0," ");
				batchlog << "Loc0 (X) and Loc1 (Y) must be 0 or more for a cell-based model" << endl;
				errors++;
			}
		}
		if (nm == "NInds") {
			if (stagestruct) {    
				if (stage < 0 || stage >= stages) {
					BatchError(filetype,line,0," ");
					batchlog << "Stage must be 0 for all stages or equal to a valid stage number in the model" << endl;
					errors++;
				}
			}
			if (reproductn == 0) {
				if (sex != 0) {
					BatchError(filetype,line,0," ");
					batchlog << "Sex must be 0 for a female-only model" << endl;
					errors++;
				}				
			}
			else {
				if (sex < 0 || sex > 2) {
					BatchError(filetype,line,2,"Sex"); errors++;					
				}
			}
		}
		if (nm == "Occupied") {
			if (value != 0.0 && value != 1.0) {
				BatchError(filetype,line,1,"Value"); errors++;
			}
		}
		break;
	case 3: // individual
		if (nm == "EmigRate" || nm == "DispSuccess") {
			if (value < 0.0 || value > 1.0) {
				BatchError(filetype,line,20,"Value"); errors++;
			}
		}
		break;
	case 4: // connectivity
		if (patchmodel) {
			if (loc0 < 1 && loc0 != -999) {
				BatchError(filetype,line,0," ");
				batchlog << "Loc0 (PatchID) must be 1 or more, or -999 for all immigrants" << endl;
				errors++;
			}
			if (loc1 < 1 && loc1 != -999) {
				BatchError(filetype,line,0," ");
				batchlog << "Loc1 (PatchID) must be 1 or more, or -999 for all emigrants" << endl;
				errors++;
			}
			if (loc0 == -999 && loc1 == -999) {
				BatchError(filetype,line,0," ");
				batchlog << "Either Loc0 or Loc1 must be a PatchID of 1 or more " << endl;
				errors++;
			}
		}
		break;
	case 5: // genome
		break;
	default:
		;
	}
	if (nm == "NInds" || nm == "NOccupPatches" || nm == "DistMean" || nm == "DistSD") {
		if (value <= 0.0) {
			BatchError(filetype,line,10,"Value"); errors++;
		}
	}
	if (weight <= 0) {
		BatchError(filetype,line,10,"Weight"); errors++;
	}
	line++;
	// read next simulation number
	inint = -98765;
	bABCoFile >> inint;
	if (bABCoFile.eof()) {
		inint = -98765;
	}
	else { // check for valid parameter number
		if (inint != prevobs+1) {
			BatchError(filetype,line,0," ");
			batchlog << "Observation numbers must be sequential integers" << endl;
			errors++;
		}
		prevobs = inint; nparams++;
	}
//	batchlog << "ParseABCFile(): First item of next line = " << inint << endl;
} // end of while loop
if (!bABCoFile.eof()) {
	EOFerror(filetype);
	errors++;
}

if (errors > 0) return -111;
else return nparams;
}

#endif // RS_ABC

//---------------------------------------------------------------------------
/*
Check stage- and sex-dependency fields for any of the dispersal files.
Check that the number of records for a simulation matches the stage-
and sex-dependency settings (unless checklines is false).
Validate the IIV field (if present).
*/
simCheck CheckStageSex(string filetype, int line, int simul, simCheck prev,
	int stagedep, int sexdep, int stage, int sex, int indvar,
	bool checklines, bool stgdepindvarok)
{
simCheck current{};
current.errors = 0;
int iii;

// has there been a change of simulation number?;
if (simul == prev.simul) { // no
	current.newsimul = false; current.simlines = prev.simlines+1;
}
else { // yes
	// check for valid simulation number
	current.newsimul = true; current.simlines = 1;
	if (line > 1 && simul != prev.simul+1) {
		BatchError(filetype,line,222," "); current.errors++;
	}
	// check for correct number of lines for previous simulation
	if (checklines && prev.simlines != prev.reqdsimlines) {
		BatchError(filetype,line,0," "); current.errors++;
		batchlog << "No. of lines for previous Simulation " << prev.simul
			<< msgshldbe << prev.reqdsimlines << endl;
	}
}
current.simul = simul;

// validate stagedep
if (stagestruct) {
	if (stagedep < 0 || stagedep > 1) {
		BatchError(filetype,line,1,"StageDep"); current.errors++;
		stagedep = 1; // to calculate required number of lines
	}
}
else {
	if (stagedep != 0) {
		BatchError(filetype,line,0," "); current.errors++;
		batchlog << "StageDep must be 0 for non-stage-structured model" << endl;
		stagedep = 0; // to calculate required number of lines
	}
}
// validate sexdep
if (sexesDisp == 2) {
	if (sexdep < 0 || sexdep > 1) {
		BatchError(filetype,line,1,"SexDep"); current.errors++;
		sexdep = 1; // to calculate required number of lines
	}
}
else {
	if (sexdep != 0) {
		BatchError(filetype,line,0," "); current.errors++;
		batchlog << "SexDep must be 0 for asexual model" << endl;
		sexdep = 0; // to calculate required number of lines
	}
}
if (current.newsimul) { // set required number of lines
	if (stagedep) {
		if (sexdep) current.reqdsimlines = stages*sexesDisp;
		else current.reqdsimlines = stages;
	}
	else {
		if (sexdep) current.reqdsimlines = sexesDisp;
		else current.reqdsimlines = 1;
	}
}
else current.reqdsimlines= prev.reqdsimlines;

// validate stage
if (stagedep) { // there must be 1 or 2 lines for each stage
	if (sexdep) { // there must be 2 lines for each stage
		if (current.simlines%2) iii = (current.simlines+1)/2; else  iii = current.simlines/2;
		if (stage != iii-1) {
			BatchError(filetype,line,0," "); current.errors++;
			batchlog << "Stages must be sequentially numbered from 0" << endl;
		}
	}
	else { // there must be 1 line for each stage
		if (stage != current.simlines-1) {
			BatchError(filetype,line,0," "); current.errors++;
			batchlog << "Stages must be sequentially numbered from 0" << endl;
		}
	}
}
else { // no stage-dependent emigration
	if (stage != 0) {
		BatchError(filetype,line,0," "); current.errors++;
		batchlog << "Stage must be 0 for non-stage-structured model" << endl;
	}
}
// validate sex
if (sexdep) {
	if (sex != (current.simlines+1)%2) {
		BatchError(filetype,line,0," "); current.errors++;
		batchlog << "Sex must be alternately 0 and 1 if SexDep is 1" << endl;
	}
}
else {
	if (sex != 0) {
		BatchError(filetype,line,0," "); current.errors++;
		batchlog << "Sex must be 0 if SexDep is 0" << endl;
	}
}

// validate indvar
if (stagedep && !stgdepindvarok) {
	if (indvar != 0) {
		BatchError(filetype,line,0," "); current.errors++;
		batchlog << "IndVar must be 0 if stage-dependent" << endl;
	}
}
else {
	if (indvar < 0 || indvar > 1) {
		BatchError(filetype,line,1,"IndVar"); current.errors++;
	}
}

return current;

}

//---------------------------------------------------------------------------

/* TEMPLATE PARSING FUNCTION
int ParseXXXXXXXXFile(void)
{
string header;
int errors = 0;
int simuls = 0;

 >> header; if (header != "" ) errors++;

if (errors > 0) return -111;
else return simuls;

}
*/

//---------------------------------------------------------------------------

// Functions to handle and report error conditions

void BatchError(string filename, int line, int option, string fieldname)
{
if (line == -999) { // message does not cite line number
	batchlog << "*** Error in " << filename << ": ";
}
else {
	batchlog << "*** Error in " << filename << " at line " << line <<": ";
}
switch (option) {
case 0:
break;
case 1:
	batchlog << fieldname << " must be 0 or 1" ;
break;
case 2:
	batchlog << fieldname << " must be 0, 1 or 2";
break;
case 3:
	batchlog << fieldname << " must be 0, 1, 2 or 3";
break;
case 4:
	batchlog << fieldname << " must be from 0 to 4";
break;
case 5:
	batchlog << fieldname << " must be from 0 to 5";
break;
case 6:
	batchlog << fieldname << " must be from 0 to 6";
break;
case 7:
	batchlog << fieldname << " must be from 0 to 7";
break;
case 10:
	batchlog << fieldname << " must be greater than zero";
break;
case 11:
	batchlog << fieldname << " must be 1 or more";
break;
case 12:
	batchlog << fieldname << " must be 2 or more";
break;
case 13:
	batchlog << fieldname << " must be 3 or more";
break;
case 18:
	batchlog << fieldname << " must be greater than 1.0";
break;
case 19:
	batchlog << fieldname << " must be 0 or more";
break;
case 20:
	batchlog << fieldname << " must be between 0 and 1";
break;
case 21:
	batchlog << fieldname << " must be greater than 1";
break;
case 33:
	batchlog << fieldname << " must be 1, 2 or 3";
break;
case 44:
	batchlog << fieldname << " must be from 1 to 4";
break;
case 55:
	batchlog << fieldname << " must be from 1 to 5";
break;
case 66:
	batchlog << fieldname << " must be from 1 to 6";
break;
case 100:
	batchlog << fieldname << " must be between 0 and 100";
break;
case 111:
	batchlog << fieldname << " must match the first Simulation in ParameterFile";
break;
case 222:
	batchlog << "Simulation numbers must be sequential integers";
break;
#if SEASONAL
case 223:
	batchlog << "Season numbers must be sequential integers";
break;
#endif
case 333:
	batchlog << "No. of " << fieldname << " columns must equal max. no. of habitats ("
		<< maxNhab << ") and be sequentially numbered starting from 1";
break;
case 444:
	batchlog << "No. of " << fieldname << " columns must be one fewer than no. of stages, i.e. "
		<< stages-1 << ", and be sequentially numbered starting from 1";
break;
case 555:
	batchlog << "No. of " << fieldname << " columns must equal no. of stages, i.e. "
		<< stages << ", and be sequentially numbered starting from 0";
break;
case 666:
	batchlog << fieldname << " must be a unique positive integer";
break;
default:
	batchlog << "*** Unspecified error regarding parameter " << fieldname;
}
if (option != 0) batchlog << endl;
}

void BatchError(string filename,int line,int option,string fieldname,string fieldname2)
{
if (line == -999) { // message does not cite line number
	batchlog << "*** Error in " << filename << ": ";
}
else {
	batchlog << "*** Error in " << filename << " at line " << line <<": ";
}
switch (option) {
case 0:
break;
case 1:
	batchlog << fieldname << " must be greater than " << fieldname2;
break;
case 2:
	batchlog << fieldname << " must be greater than or equal to " << fieldname2;
break;
case 3:
	batchlog << fieldname << " must be less than or equal to " << fieldname2;
break;
case 4:
	batchlog << fieldname << " must be less than " << fieldname2;
break;
default:
	batchlog << "*** Unspecified error regarding parameters " << fieldname
		<< " and " << fieldname2;
	;
}
if (option != 0) batchlog << endl;
}

void CtrlFormatError(void)
{
cout << "Format error in Control file" << endl;
batchlog << endl << "***" << endl << "*** Format error in Control file:"
//	<< " case-sensitive parameter and file names must match the specification exactly"
	<< msgcase << " and file names" << msgmatch
	<< endl
	<< "***" << endl;
}

void ArchFormatError(void)
{
//batchlog << "*** Format error in ArchFile: case-sensitive parameter names "
//	<< "must match the specification exactly" << endl;
batchlog << "*** Format error in ArchFile:" << msgcase << msgmatch << endl;
}

#if VIRTUALECOLOGIST
void SampleFormatError(void)
{
batchlog << "*** Format error in SampleFile:" << msgcase << msgmatch << endl;
}
#endif

void FormatError(string filename, int errors)
{
batchlog << "*** Format error in header line of ";
if (errors == 0) {
	batchlog << filename << endl;
}
else {
	batchlog << filename << ": " << errors << " error";
	if (errors > 1) batchlog << "s";
	batchlog << " detected" << endl;
}
}

void OpenError(string ftype, string fname)
{
batchlog << "*** Unable to open " << ftype << " " << fname << endl;
}

void EOFerror(string filename)
{
batchlog << "*** Failed to read to EOF in " << filename << endl;
}

void FileOK(string ftype, int n, int option)
{
batchlog << ftype << " OK: total no. of ";
switch (option) {
case 0:
	batchlog << "simulations = ";
	break;
case 1:
	batchlog << "landscapes = ";
	break;
case 2:
	batchlog << "parameters = ";
	break;
#if SEASONAL
case 3:
	batchlog << "seasons = ";
	break;
#endif
default:
	batchlog << "PROBLEMS = ";
}
//if (option == 0) batchlog << "simulations = ";
//else batchlog << "landscapes = ";
batchlog << n << endl;
}

void FileHeadersOK(string filename)
{
batchlog << filename << " OK" << endl;
}

void SimulnCountError(string filename)
{
batchlog << "*** No. of simulations in " << filename
	<< " does not match no. in ParameterFile" << endl;
}

//---------------------------------------------------------------------------
int ReadLandFile(int option)
{
#if RSDEBUG
DEBUGLOG << "ReadLandFile(): option=" << option
	<< " landFile=" << landFile << endl;
#endif

if (option == 0) { // open file and read header line
	landfile.open(landFile.c_str());
	if (landfile.is_open()) {
		string header;
		int nheaders;
		if (landtype == 9) nheaders = 9; // artificial landscape
		else { // imported raster map
//			nheaders = 6;
			nheaders = 7;
#if RS_CONTAIN
			nheaders += 1; 
#endif // RS_CONTAIN 
#if SPATIALMORT
			nheaders += 2; 
#endif // SPATIALMORT 
		}
		for (int i = 0; i < nheaders; i++) landfile >> header;
	}
	else return 1;
}

if (option == 9) { // close file
	if (landfile.is_open()) {
		landfile.close();  landfile.clear();
	}
}
return 0;
}

//---------------------------------------------------------------------------
int ReadLandFile(int option, Landscape *pLandscape)
{
#if RSDEBUG
DEBUGLOG << "ReadLandFile(): option=" << option << endl;
#endif
landParams ppLand = pLandscape->getLandParams();
genLandParams ppGenLand = pLandscape->getGenLandParams();
simParams sim = paramsSim->getSim();

if (landtype == 9) { //artificial landscape
#if RSDEBUG
DEBUGLOG << "ReadLandFile(): artificial: "  << endl;
#endif
	ppLand.nHab = 2;
  ppLand.rasterType = 9;
	landfile >> ppLand.landNum >> ppGenLand.fractal >> ppGenLand.continuous
		>> ppLand.dimX >> ppLand.dimY >> ppGenLand.minPct >> ppGenLand.maxPct
		>> ppGenLand.propSuit >> ppGenLand.hurst;
	ppLand.maxX = ppLand.dimX-1; ppLand.maxY = ppLand.dimY-1;
	if (ppGenLand.fractal && ppLand.maxX > ppLand.maxY) {
		return -901;
	}
	if (ppGenLand.fractal) {
		if ((ppLand.dimX < 3 || ppLand.dimX%2 != 1)
		||  (ppLand.dimY < 3 || ppLand.dimY%2 != 1)) {
			return -902;
		}
	}
	// SCFP 26/9/13 - min and max habitat percentages need to be set for all types of
	// fractal landscape (including discrete), as they are passed to the fractal generator
	// NOTE that will not have been checked for a discrete landscape
	if (ppGenLand.fractal && !ppGenLand.continuous) { ppGenLand.minPct = 1; ppGenLand.maxPct = 100; }
#if RSDEBUG
DEBUGLOG << "ReadLandFile(): ppLand.landNum=" << ppLand.landNum
	<< " ppLand.dimX=" << ppLand.dimX << " ppLand.dimY=" << ppLand.dimY
	<< " ppLand.maxX=" << ppLand.maxX << " ppLand.maxY=" << ppLand.maxY
	<< " propSuit=" << ppGenLand.propSuit
	<< endl;
#endif
}
else { // imported raster map
	string dummy; // no longer necessary to read no. of habitats from landFile
//	landfile >> ppGenLand.landNum >> ppLand.nHab >> name_landscape >> name_patch >> name_sp_dist;
	landfile >> ppLand.landNum >> dummy >> name_landscape >> name_patch;
	landfile >> name_costfile >> name_dynland >> name_sp_dist;
#if RS_CONTAIN
	landfile >> name_damagefile;
#endif // RS_CONTAIN 
#if SPATIALMORT
	landfile >> name_mortfile[0] >> name_mortfile[1];
#endif // SPATIALMORT 
	if (landtype == 2) ppLand.nHab = 1; // habitat quality landscape has one habitat class
#if RSDEBUG
DEBUGLOG << "ReadLandFile(): ppLand.landNum=" << ppLand.landNum
	<< " name_landscape=" << name_landscape
	<< " name_patch=" << name_patch
	<< " name_costfile=" << name_costfile
	<< " name_dynland=" << name_dynland
	<< " name_sp_dist=" << name_sp_dist
#if RS_CONTAIN
	<< " name_damagefile=" << name_damagefile
#endif // RS_CONTAIN 
	<< endl;
#endif
}

pLandscape->setLandParams(ppLand,sim.batchMode);
pLandscape->setGenLandParams(ppGenLand);

//simParams sim = paramsSim->getSim();

#if RSDEBUG
//DEBUGLOG << "ReadLandFile(): NHab=" << ppLand.nHab << endl;
DEBUGLOG << "ReadLandFile(): ppLand.landNum=" << ppLand.landNum << endl;
#endif

return ppLand.landNum;
}

//---------------------------------------------------------------------------
int ReadDynLandFile(Landscape *pLandscape) {
#if RSDEBUG
DEBUGLOG << "ReadDynLandFile(): pLandscape=" << pLandscape
	<< " name_dynland=" << name_dynland
	<< endl;
#endif
//int change,year;
string landchangefile,patchchangefile,costchangefile;
int change,imported;
int nchanges = 0;
landChange chg;
landParams ppLand = pLandscape->getLandParams();
string fname = paramsSim->getDir(1) + name_dynland;

dynlandfile.open(fname.c_str());
if (dynlandfile.is_open()) {
	string header;
	int nheaders = 5;
	for (int i = 0; i < nheaders; i++) dynlandfile >> header;
}
else {
	dynlandfile.clear();
	return 72727;
}

// read data lines
change = -98765;
dynlandfile >> change; // first change number
while (change != -98765) {
	chg.chgnum = change;                   
	dynlandfile >> chg.chgyear >> landchangefile >> patchchangefile >> costchangefile;
//	dynlandfile >> chg.chgyear >> chg.habfile >> chg.pchfile;
	chg.habfile = paramsSim->getDir(1) + landchangefile;
	chg.pchfile = paramsSim->getDir(1) + patchchangefile;
	if (costchangefile == "NULL") chg.costfile = "none";
	else chg.costfile = paramsSim->getDir(1) + costchangefile;
	nchanges++;
	pLandscape->addLandChange(chg);
	// read first field on next line
	change = -98765;
	dynlandfile >> change;
	if (dynlandfile.eof()) {
		change = -98765;
	}
}

dynlandfile.close(); dynlandfile.clear();

// read landscape change maps
if (ppLand.patchModel) {
	pLandscape->createPatchChgMatrix();
}
if (costchangefile != "NULL") {
	pLandscape->createCostsChgMatrix();
}
for (int i = 0; i < nchanges; i++) {
	if (costchangefile == "NULL") imported = pLandscape->readLandChange(i,false);
	else imported = pLandscape->readLandChange(i,true);
	if (imported != 0) {
		return imported;
	}
	if (ppLand.patchModel) {
		pLandscape->recordPatchChanges(i+1);
	}
	if (costchangefile != "NULL") {
		pLandscape->recordCostChanges(i+1);
	}
}
if (ppLand.patchModel) {
	// record changes back to original landscape for multiple replicates
	pLandscape->recordPatchChanges(0);
	pLandscape->deletePatchChgMatrix();
}
if (costchangefile != "NULL") {
	pLandscape->recordCostChanges(0);
	pLandscape->deleteCostsChgMatrix();
}

#if RSDEBUG
DEBUGLOG << "ReadDynLandFile(): finished" << endl;
#endif
return 0;
}

//---------------------------------------------------------------------------
int ReadParameters(int option, Landscape *pLandscape)
{
#if RSDEBUG
DEBUGLOG << endl << "ReadParameters(): option=" << option
//	<< " parameters=" << parameters
	<< endl;
#endif
int iiii,jjjj;
int error = 0;
landParams paramsLand = pLandscape->getLandParams();

if (option == 0) { // open file and read header line
	parameters.open(parameterFile.c_str());
	if (parameters.is_open()) {
		string header;
		int nheaders = 47 + paramsLand.nHabMax;
#if RS_CONTAIN
		nheaders += 2;	// ADDITIONAL HEADERS FOR ADAPTIVE MANAGEMENT MODEL
#endif // RS_CONTAIN 
#if SEASONAL
		nheaders +=	paramsLand.nHabMax * (nseasons-1); // ADDITIONAL HEADER FOR SEASONAL MODEL 
#endif // SEASONAL 
#if GROUPDISP
		nheaders += 4;	// ADDITIONAL HEADERS FOR GROUP DISPERSAL MODEL
#if PEDIGREE
		nheaders += 1;	// ADDITIONAL HEADER FOR PEDIGREE MATRIX 
#endif // PEDIGREE
#endif // GROUPDISP 
#if SOCIALMODEL
		nheaders += 14; // ADDITIONAL HEADERS FOR PROBIS SOCIAL POLYMORPHISM MODEL
#endif
#if SPATIALMORT
		nheaders += 1;	// ADDITIONAL HEADER FOR SPATIALLY VARYING MORTALITY
#endif
#if BUTTERFLYDISP
		nheaders += 3;	// ADDITIONAL HEADER FOR BUTTERFLY DISPERSAL MODEL
#endif
		for (int i = 0; i < nheaders; i++) parameters >> header;
		return 0;
	}
	else return 1;
}

if (option == 9) { // close file
	if (parameters.is_open()) {
		parameters.close();  parameters.clear();
	}
	return 0;
}

envStochParams env = paramsStoch->getStoch();
demogrParams dem = pSpecies->getDemogr();
simParams sim = paramsSim->getSim();
simView v = paramsSim->getViews();

if (!parameters.is_open()) {
	cout << endl << "ReadParameters(): ERROR - ParameterFile is not open" << endl;
	return 4086534;
}

int gradType,shift_begin,shift_stop;
float grad_inc,opt_y,f,optEXT,shift_rate;
bool shifting;

parameters >> sim.simulation >> sim.reps >> sim.years;
#if SPATIALMORT
parameters >> sim.mortChgYear;
#endif
parameters >> iiii;
if (iiii == 1) sim.absorbing = true; else sim.absorbing = false;
#if RSDEBUG
//DEBUGLOG << "ReadParameters(): paramsSim=" << paramsSim << endl;
DEBUGLOG << "ReadParameters(): simulation=" << sim.simulation
	<< " reps=" << sim.reps << " years=" << sim.years << endl;
#endif
parameters >> gradType;
parameters >> grad_inc >> opt_y >> f >> optEXT >> iiii >> shift_rate;
if (iiii == 1 && gradType != 0) shifting = true; else shifting = false;
parameters >> shift_begin >> shift_stop;
paramsGrad->setGradient(gradType,grad_inc,opt_y,f,optEXT);
if (shifting) paramsGrad->setShifting(shift_rate,shift_begin,shift_stop);
else paramsGrad->noShifting();

parameters >> iiii;
if (iiii == 0) env.stoch = false;
else {
	env.stoch = true;
	if (iiii == 2) env.local = true; else env.local = false;
}
if (paramsLand.patchModel && env.local) error = 101;
parameters >> iiii;
if (iiii == 1) env.inK = true; else env.inK = false;
#if BUTTERFLYDISP
parameters >> iiii;
if (iiii == 1) env.fromFile = true; else env.fromFile = false;
string envfile;
parameters >> envfile;
if (env.fromFile) {
	envstochfilename = paramsSim->getDir(1) + envfile;
}
#endif
// as from v1.1, there is just one pair of min & max values,
// which are attributes of the species
// ULTIMATELY, THE PARAMETER FILE SHOULD HAVE ONLY TWO COLUMNS ...
//parameters >> env.ac >> env.std >> env.minR >> env.maxR >> env.minK >> env.maxK;
float minR,maxR,minK,maxK; 
parameters >> env.ac >> env.std >> minR >> maxR >> minK >> maxK;
if (env.inK) {
	float minKK,maxKK;
	minKK = minK * ((((float)paramsLand.resol * (float)paramsLand.resol)) / 10000.0f);
	maxKK = maxK * ((((float)paramsLand.resol * (float)paramsLand.resol)) / 10000.0f);
	pSpecies->setMinMax(minKK,maxKK);
}
else pSpecies->setMinMax(minR,maxR);
#if RSDEBUG
//DEBUGLOG << "ReadParameters(): minR=" << env.minR << " maxR=" << env.maxR << endl;
#endif
parameters >> iiii;
if (iiii == 1) env.localExt = true; else env.localExt = false;
if (paramsLand.patchModel && env.localExt) error = 102;
parameters >> env.locExtProb;
paramsStoch->setStoch(env);

#if BUTTERFLYDISP
parameters >> iiii;
if (dem.stageStruct) dem.dispersal = 1; else dem.dispersal = iiii;
#endif
parameters >> dem.propMales >> dem.harem >> dem.bc >> dem.lambda;
#if !GROUPDISP
pSpecies->setDemogr(dem);
#endif

float k;

if (landtype == 9) { // artificial landscape
	// only one value of K is read, but the first 'habitat' is the matrix where K = 0
#if SEASONAL
	pSpecies->createHabK(2,1);
#else
	pSpecies->createHabK(2);
#endif // SEASONAL 
	parameters >> k;
	k *= (((float)paramsLand.resol*(float)paramsLand.resol))/10000.0f;
#if SEASONAL
	pSpecies->setHabK(0,1,0);
	pSpecies->setHabK(1,1,k);
#else
	pSpecies->setHabK(0,0);
	pSpecies->setHabK(1,k);
#endif // SEASONAL 
}
else {
#if SEASONAL
	pSpecies->createHabK(paramsLand.nHabMax,dem.nSeasons);
	for (int j = 0; j < dem.nSeasons; j++) {
		for (int i = 0; i < paramsLand.nHabMax; i++) {
			parameters >> k;
			k *= ((double)(paramsLand.resol*paramsLand.resol))/10000.0;
				pSpecies->setHabK(i,j,k);			
		}
	}
#else
	pSpecies->createHabK(paramsLand.nHabMax);
	for (int i = 0; i < paramsLand.nHabMax; i++) {
		parameters >> k;
		k *= (((float)paramsLand.resol*(float)paramsLand.resol))/10000.0f;
		pSpecies->setHabK(i,k);
	}
#endif // SEASONAL 
}

#if RSDEBUG
DEBUGLOG << "ReadParameters(): dem.lambda=" << dem.lambda
#if SEASONAL
	<< " habK[0] = " << pSpecies->getHabK(0,0)
#else
	<< " habK[0] = " << pSpecies->getHabK(0)
#endif // SEASONAL 
	<< " nHabMax = " << paramsLand.nHabMax
	<< endl;
#endif

#if GROUPDISP
// ADDITIONAL PARAMETERS FOR GROUP DISPERSAL MODEL
parameters >> iiii;
if (iiii == 1) dem.selfing = true; else dem.selfing = false;
parameters >> dem.paternity >> dem.propLocal >> dem.propNghbr;
pSpecies->setDemogr(dem);
#if PEDIGREE
parameters >> sim.relMatSize;
#endif // PEDIGREE
#endif
#if SOCIALMODEL
// ADDITIONAL PARAMETERS FOR PROBIS SOCIAL POLYMORPHISM MODEL
socialParams soc;
//parameters >> soc.asocK >> soc.asocRmax >> soc.asocBc >> soc.ra >> soc.rs
//	>> soc.Ta >> soc.Ts >> soc.dK >> soc.alpha;
//parameters >> soc.asocK >> soc.asocRmax >> soc.asocBc
//	>> soc.Ta >> soc.Ts >> soc.Ca >> soc.Cs >> soc.dK >> soc.alpha;
parameters >> soc.socMean >> soc.socSD >> soc.socScale
	>> soc.asocK >> soc.asocRmax >> soc.asocBc
	>> soc.Ta >> soc.Ts >> soc.ca >> soc.cs >> soc.ba >> soc.bs >> soc.dK >> soc.alpha;
pSpecies->setSocialParams(soc);
#endif

parameters >> sim.outStartPop >> sim.outStartInd >> sim.outStartGenetic
	>> sim.outStartTraitCell >> sim.outStartTraitRow >> sim.outStartConn;
#if RS_CONTAIN
parameters >> sim.outStartDamage;
#endif // RS_CONTAIN 
parameters >> sim.outIntRange >> sim.outIntOcc >> sim.outIntPop >> sim.outIntInd
	>> sim.outIntGenetic >> sim.outGenType >> jjjj
	>> sim.outIntTraitCell >> sim.outIntTraitRow >> sim.outIntConn;
#if RS_CONTAIN
parameters >> sim.outIntDamage;       
#endif // RS_CONTAIN 
if (sim.outIntRange > 0)     sim.outRange = true; else sim.outRange = false;
if (sim.outIntOcc > 0)       sim.outOccup = true; else sim.outOccup = false;
if (sim.outIntPop > 0)       sim.outPop = true; else sim.outPop = false;
if (sim.outIntInd > 0)       sim.outInds = true; else sim.outInds = false;
if (sim.outIntGenetic > 0)   sim.outGenetics = true; else sim.outGenetics = false;
if (jjjj == 1)   						 sim.outGenXtab = true; else sim.outGenXtab = false;
if (sim.outIntRange > 0)     sim.outRange = true; else sim.outRange = false;
if (sim.outIntTraitCell > 0) sim.outTraitsCells = true; else sim.outTraitsCells = false;
if (sim.outIntTraitRow > 0)  sim.outTraitsRows = true; else sim.outTraitsRows = false;
if (sim.outIntConn > 0)      sim.outConnect = true; else sim.outConnect = false;
#if RS_CONTAIN
if (sim.outIntDamage > 0)    sim.outDamage = true; else sim.outDamage = false;
#endif // RS_CONTAIN 
if (sim.outOccup && sim.reps < 2) error = 103;
if (paramsLand.patchModel) {
	if (sim.outTraitsRows) error = 104;
}
else{
	if (sim.outConnect) error = 105;
}
#if RSDEBUG
DEBUGLOG << "ReadParameters(): outRange=" << sim.outRange << " outInt=" << sim.outIntRange
	<< endl;
#endif
parameters >> iiii >> sim.mapInt;
if (iiii == 0) sim.saveMaps = false;
else sim.saveMaps = true;
parameters >> iiii;
if (iiii == 0) sim.saveVisits = false;
else sim.saveVisits = true;
parameters >> iiii;
//sim.saveInitMap = false;
if (iiii == 0) sim.drawLoaded = false; else sim.drawLoaded = true;

paramsSim->setSim(sim);
paramsSim->setViews(v);

return error;
}

//---------------------------------------------------------------------------
#if SEASONAL
int ReadStageStructure(int option,Landscape *pLandscape)
#else
int ReadStageStructure(int option)
#endif  
{
string name;
int simulation,postDestructn;
#if RS_CONTAIN || SEASONAL
demogrParams dem = pSpecies->getDemogr();
#endif
stageParams sstruct = pSpecies->getStage();
string Inputs = paramsSim->getDir(1);

if (option == 0) { // open file and read header line
	ssfile.open(stageStructFile.c_str());
	string header;
	int nheaders = 18;
#if RS_CONTAIN
	nheaders += 1;
#endif // RS_CONTAIN 
#if GOBYMODEL
	nheaders += 4;
#endif // GOBYMODEL 
#if SEASONAL
	nheaders += 1;
#if PARTMIGRN
	nheaders += 7;
#endif // PARTMIGRN 
#endif // SEASONAL
#if RSDEBUG
DEBUGLOG << "ReadStageStructure{}: nheaders=" << nheaders << endl;
#endif
	for (int i = 0; i < nheaders; i++) ssfile >> header;
	return 0;
}

if (option == 9) { // close file
	if (ssfile.is_open()) {
		ssfile.close(); ssfile.clear();
	}
	return 0;
}

#if RSDEBUG
DEBUGLOG << "ReadStageStructure{}: sstruct.nStages=" << sstruct.nStages << endl;
#endif

ssfile >> simulation;
#if GOBYMODEL
socialParams soc;
ssfile >> soc.socMean >> soc.socSD >> soc.socScale >> sstruct.asocF;
pSpecies->setSocialParams(soc);
#endif
ssfile >> postDestructn >> sstruct.probRep >> sstruct.repInterval >> sstruct.maxAge;
if (postDestructn == 1) sstruct.disperseOnLoss = true;
else sstruct.disperseOnLoss = false;

ssfile >> name; 
#if RS_CONTAIN
pSpecies->resetDem(-1); // reset demography for all habitats
// first 'name' is HabDemFile
if (name == "NULL") { 
	dem.habDepDem = false;
}
else {
	dem.habDepDem = true;
	hdfile.open((Inputs+name).c_str());
	ReadHabDemFile(sstruct.nStages,sexesDem);
	hdfile.close(); hdfile.clear();	
}
pSpecies->setDemogr(dem);
ssfile >> name; 
// second 'name' is TransMatrixFile
if (name != "NULL") {
#if SEASONAL
	for (int i = 0; i < dem.nSeasons; i++) {
		tmfile.open((Inputs+name).c_str());
		ReadTransitionMatrix(sstruct.nStages,sexesDem,0,i);		
		tmfile.close(); tmfile.clear();	
	}
#else
	tmfile.open((Inputs+name).c_str());
	ReadTransitionMatrix(sstruct.nStages,sexesDem,0,0);
	tmfile.close(); tmfile.clear();	
#endif // SEASONAL 
}
#else
#if SEASONAL
// 'name' is SeasonFile
seasonfile.open((Inputs+name).c_str());
ReadSeasonFile(1,dem.nSeasons);
seasonfile.close(); seasonfile.clear();
#else
// 'name' is TransMatrixFile
tmfile.open((Inputs+name).c_str());
ReadTransitionMatrix(sstruct.nStages,sexesDem,0,0);
tmfile.close(); tmfile.clear();
#endif // SEASONAL 
#endif // RS_CONTAIN 
ssfile >> sstruct.survival; 

float devCoeff,survCoeff;
ssfile >> sstruct.fecDens >> sstruct.fecStageDens >> name; // 'name' is FecStageWtsFile
if (name != "NULL") {
	fdfile.open((Inputs+name).c_str());
	ReadStageWeights(1);
	fdfile.close(); fdfile.clear();
}
ssfile >> sstruct.devDens >> devCoeff >> sstruct.devStageDens >> name; // 'name' is DevStageWtsFile
if (name != "NULL") {
	ddfile.open((Inputs+name).c_str());
	ReadStageWeights(2);
	ddfile.close(); ddfile.clear();
}
ssfile >> sstruct.survDens >> survCoeff >> sstruct.survStageDens >> name; // 'name' is SurvStageWtsFile
if (name != "NULL") {
	sdfile.open((Inputs+name).c_str());
	ReadStageWeights(3);
	sdfile.close(); sdfile.clear();
}
#if PARTMIGRN
float prop;
int resetmigrn;
for (short i = 1; i < 7; i++) { ssfile >> prop; pSpecies->setPropDispMigrn(i,prop); }
ssfile >> resetmigrn;
if (resetmigrn == 1) pSpecies->setResetMigrn(true);
else pSpecies->setResetMigrn(true);
#endif // PARTMIGRN 
#if SEASONAL
ssfile >> name; // 'name' is ExtremeFile
if (name != "NULL") {
	pLandscape->resetExtEvents();
	extremefile.open((Inputs+name).c_str());
	ReadExtremeFile(pLandscape,nseasons);
	extremefile.close(); extremefile.clear();
}
#endif // SEASONAL 

pSpecies->setStage(sstruct);

if (sstruct.devDens || sstruct.survDens) {
	pSpecies->setDensDep(devCoeff,survCoeff);
}

return 0;
}

#if RS_CONTAIN

int ReadHabDemFile(const short nstages,const short nsexesDem) {
#if RSDEBUG
DEBUGLOG << "ReadHabDemFile(): nstages=" << nstages << " nsexesDem=" << nsexesDem << endl;
#endif
short habcode;        
string name;
string Inputs = paramsSim->getDir(1);

// read header line
string header;
int nheaders = 2;
#if SEASONAL
nheaders += 1;
short season;
#endif
for (int i = 0; i < nheaders; i++) hdfile >> header;

habcode = -1;
hdfile >> habcode;
while (habcode != -1) 
//for (int i = 0; i < maxNhab; i++) 
{	
//	habcode = -1;
//	hdfile >> habcode;
	if (habcode < 0) { // EOF has been reached
		return 1;
	}
#if SEASONAL
	hdfile >> season; 
#endif
	hdfile >> name; // 'name' is TransMatrixFile
#if RSDEBUG
//DEBUGLOG << "ReadHabDemFile(): i=" << i << " habcode=" << habcode 
DEBUGLOG << "ReadHabDemFile(): habcode=" << habcode 
#if SEASONAL
	<< " season=" << season
#endif
	<< " name=" << name << endl;
#endif
	tmfile.open((Inputs+name).c_str());
	// user supplies habitat code, but as habitats must be sequential starting from 1
	// then habitat index is habcode-1
#if SEASONAL
	ReadTransitionMatrix(nstages,nsexesDem,habcode-1,season);    
#else
	ReadTransitionMatrix(nstages,nsexesDem,habcode-1,0);    
#endif
	tmfile.close(); tmfile.clear();
	habcode = -1;
	hdfile >> habcode;
}

return 0;
}

#endif // RS_CONTAIN 

//---------------------------------------------------------------------------
int ReadTransitionMatrix(short nstages,short nsexesDem,short hab,short season)  
{
//#if RS_CONTAIN
//int hab = 0; // TEMPORARY set suitable habitat to 0
//#endif // RS_CONTAIN 
int ii;
int minAge;
float ss, dd; 
string header;
demogrParams dem = pSpecies->getDemogr();
//stageParams sstruct = pSpecies->getStage();
#if SEASONAL
bool breeding = false;
#if RSDEBUG
//DEBUGLOG << "Read_TransitionMatrix(): season=" << season << endl;
#endif
#endif // SEASONAL

// read header line
//for (int i = 0; i < (sstruct.nStages*nsexesDem)+2; i++)
for (int i = 0; i < (nstages*nsexesDem)+2; i++)
{
	tmfile >> header;
#if RSDEBUG
//DEBUGLOG << "Read_TransitionMatrix(): i= << i << " header=" << header << endl;
#endif
}

if (matrix != NULL) {
	for (int j = 0; j < matrixsize; j++) delete[] matrix[j];
	delete[] matrix;
	matrix = NULL; matrixsize = 0;
}

if (dem.repType != 2) { // asexual or implicit sexual model
#if RSDEBUG
//DEBUGLOG << "Read_TransitionMatrix(): asexual model, sstruct.nStages = " << sstruct.nStages << endl;
#endif
	// create a temporary matrix
	matrix = new float *[nstages];
	matrixsize = nstages;
	for (int i = 0; i < nstages; i++)
		matrix[i] = new float [nstages];

//	for (int i = 0; i < sstruct.nStages; i++) 
	for (int i = 0; i < nstages; i++)
	{ // i = row; j = coloumn
		tmfile >> header;
#if RSDEBUG
//DEBUGLOG << "Read_TransitionMatrix(): i=" << i << " header=" << header << endl;
#endif
//		for (int j = 0; j < sstruct.nStages; j++) 
		for (int j = 0; j < nstages; j++)
		{
			tmfile >> matrix[j][i];
#if RSDEBUG
//DEBUGLOG << "Read_TransitionMatrix(): j=" << j << " matrix[j][i]=" << matrix[j][i] << endl;
#endif
		}
		tmfile >> minAge; pSpecies->setMinAge(i,0,minAge);
	}

//	for (int j = 1; j < sstruct.nStages; j++)
	for (int j = 1; j < nstages; j++)
#if RS_CONTAIN
#if SEASONAL
		{
		pSpecies->setFec(hab,season,j,0,matrix[j][0]);
		if (matrix[j][0] > 0.0) breeding = true;
		}
#else
		pSpecies->setFec(hab,j,0,matrix[j][0]);
#endif // SEASONAL 
#else
#if SEASONAL
		{
		pSpecies->setFec(season,j,0,matrix[j][0]);
		if (matrix[j][0] > 0.0) breeding = true;
		}
#else
		pSpecies->setFec(j,0,matrix[j][0]);
#endif // SEASONAL 
#endif // RS_CONTAIN 
//	for (int j = 0; j < sstruct.nStages; j++) 
	for (int j = 0; j < nstages; j++)
	{
		ss = 0.0; dd = 0.0;
//		for (int i = 0; i < sstruct.nStages; i++) 
		for (int i = 0; i < nstages; i++)
		{
			if (i == j) ss = matrix[j][i];
			if (i == (j+1)) dd = matrix[j][i];
		}
#if RS_CONTAIN
#if SEASONAL
		pSpecies->setSurv(hab,season,j,0,ss+dd);
		if ((ss+dd) > 0.0)
			pSpecies->setDev(hab,season,j,0,dd/(ss+dd));
		else
			pSpecies->setDev(hab,season,j,0,0.0);
#else
		pSpecies->setSurv(hab,j,0,ss+dd);
		if ((ss+dd) > 0.0)
			pSpecies->setDev(hab,j,0,dd/(ss+dd));
		else
			pSpecies->setDev(hab,j,0,0.0);
#endif // SEASONAL 
#else
#if SEASONAL
		pSpecies->setSurv(season,j,0,ss+dd);
		if ((ss+dd) > 0.0)
			pSpecies->setDev(season,j,0,dd/(ss+dd));
		else
			pSpecies->setDev(season,j,0,0.0);
#else
		pSpecies->setSurv(j,0,ss+dd);
		if ((ss+dd) > 0.0f)
			pSpecies->setDev(j,0,dd/(ss+dd));
		else
			pSpecies->setDev(j,0,0.0);
#endif // SEASONAL 
#endif // RS_CONTAIN 
	}
}
else { // complex sexual model
#if RSDEBUG
//DEBUGLOG << "Read_TransitionMatrix(): complex sexual model, sstruct.nStages = "
//	<< sstruct.nStages << endl;
#endif
	// create a temporary matrix
//	matrix = new float *[sstruct.nStages*2];
//	matrixsize = sstruct.nStages*2;
//	for (int j = 0; j < sstruct.nStages*2; j++)
//		matrix[j] = new float [sstruct.nStages*2-1];
	matrix = new float *[nstages*2];
	matrixsize = nstages*2;
	for (int j = 0; j < nstages*2; j++)
		matrix[j] = new float [nstages*2-1];

//	for (int i = 0; i < sstruct.nStages*2-1; i++) 
	for (int i = 0; i < nstages*2-1; i++)
	{ // i = row; j = coloumn
		tmfile >> header;
#if RSDEBUG
//DEBUGLOG << "Read_TransitionMatrix{}: i = " << i << " header = " << header << endl;
#endif
//		for (int j = 0; j < sstruct.nStages*2; j++) tmfile >> matrix[j][i];
		for (int j = 0; j < nstages*2; j++) tmfile >> matrix[j][i];
		if (i == 0) {
			tmfile >> minAge; pSpecies->setMinAge(i,0,minAge); pSpecies->setMinAge(i,1,minAge);
		}
		else {
			tmfile >> minAge;
			if (i%2) pSpecies->setMinAge((i+1)/2,1,minAge);	// odd lines  - males
			else 		 pSpecies->setMinAge(i/2,0,minAge);			// even lines - females
		}
	}
#if RSDEBUG
//	DEBUGLOG << endl;
//for (int ii = 0; ii < sstruct.nStages*2-1; ii++) { // row (0 = juvs, 1,2 = stage 1)
//	for (int jj = 0; jj < sstruct.nStages*2; jj++) { // column (m f m f)
//		DEBUGLOG << matrix[jj][ii] << " " ;   // matrix[column][row]
//	}
//	DEBUGLOG << endl;
//}
//	DEBUGLOG << endl;
#endif

	ii = 1;
//	for (int j = 2; j < sstruct.nStages*2; j++) 
	for (int j = 2; j < nstages*2; j++)
	{
#if RS_CONTAIN
#if SEASONAL
		if (j%2 == 0)
			pSpecies->setFec(hab,season,ii,1,matrix[j][0]);
		else {
			pSpecies->setFec(hab,season,ii,0,matrix[j][0]);
			ii++;
		}
#else
		if (j%2 == 0)
			pSpecies->setFec(hab,ii,1,matrix[j][0]);
		else {
			pSpecies->setFec(hab,ii,0,matrix[j][0]);
			ii++;
		}
#endif // SEASONAL 
#else
#if SEASONAL
		if (j%2 == 0)
			pSpecies->setFec(season,ii,1,matrix[j][0]);
		else {
			pSpecies->setFec(season,ii,0,matrix[j][0]);
			if (matrix[j][0] > 0.0) breeding = true;
			ii++;
		}
#else
		if (j%2 == 0)
			pSpecies->setFec(ii,1,matrix[j][0]);
		else {
			pSpecies->setFec(ii,0,matrix[j][0]);
			ii++;
		}
#endif // SEASONAL 
#endif // RS_CONTAIN 
	}
	// survival and development of male juveniles
#if RS_CONTAIN
#if SEASONAL
	pSpecies->setSurv(hab,season,0,1,(matrix[0][0]+matrix[0][1]));
	if ((matrix[0][0]+matrix[0][1]) > 0.0)
		pSpecies->setDev(hab,season,0,1,(matrix[0][1]/(matrix[0][0]+matrix[0][1])));
	else
		pSpecies->setDev(hab,season,0,1,0.0);
	// survival and development of female juveniles
	pSpecies->setSurv(hab,season,0,0,(matrix[1][0]+matrix[1][2]));
	if ((matrix[1][0]+matrix[1][2]) > 0.0)
		pSpecies->setDev(hab,season,0,0,(matrix[1][2]/(matrix[1][0]+matrix[1][2])));
	else
		pSpecies->setDev(hab,season,0,0,0.0);
#else
	pSpecies->setSurv(hab,0,1,(matrix[0][0]+matrix[0][1]));
	if ((matrix[0][0]+matrix[0][1]) > 0.0)
		pSpecies->setDev(hab,0,1,(matrix[0][1]/(matrix[0][0]+matrix[0][1])));
	else
		pSpecies->setDev(hab,0,1,0.0);
	// survival and development of female juveniles
	pSpecies->setSurv(hab,0,0,(matrix[1][0]+matrix[1][2]));
	if ((matrix[1][0]+matrix[1][2]) > 0.0)
		pSpecies->setDev(hab,0,0,(matrix[1][2]/(matrix[1][0]+matrix[1][2])));
	else
		pSpecies->setDev(hab,0,0,0.0);
#endif // SEASONAL 
#else
#if SEASONAL
	pSpecies->setSurv(season,0,1,(matrix[0][0]+matrix[0][1]));
	if ((matrix[0][0]+matrix[0][1]) > 0.0)
		pSpecies->setDev(season,0,1,(matrix[0][1]/(matrix[0][0]+matrix[0][1])));
	else
		pSpecies->setDev(season,0,1,0.0);
	// survival and development of female juveniles
	pSpecies->setSurv(season,0,0,(matrix[1][0]+matrix[1][2]));
	if ((matrix[1][0]+matrix[1][2]) > 0.0)
		pSpecies->setDev(season,0,0,(matrix[1][2]/(matrix[1][0]+matrix[1][2])));
	else
		pSpecies->setDev(season,0,0,0.0);
#else
	pSpecies->setSurv(0,1,(matrix[0][0]+matrix[0][1]));
	if ((matrix[0][0]+matrix[0][1]) > 0.0)
		pSpecies->setDev(0,1,(matrix[0][1]/(matrix[0][0]+matrix[0][1])));
	else
		pSpecies->setDev(0,1,0.0);
	// survival and development of female juveniles
	pSpecies->setSurv(0,0,(matrix[1][0]+matrix[1][2]));
	if ((matrix[1][0]+matrix[1][2]) > 0.0)
		pSpecies->setDev(0,0,(matrix[1][2]/(matrix[1][0]+matrix[1][2])));
	else
		pSpecies->setDev(0,0,0.0);
#endif // SEASONAL 
#endif // RS_CONTAIN 
	// survival and development of stages 1+
	ii = 1;
//	for (int j = 2; j < sstruct.nStages*2; j++) 
	for (int j = 2; j < nstages*2; j++)
	{
		ss = 0.0; dd = 0.0;
		if (j%2 == 0){ // males
//			for (int i = 0; i < sstruct.nStages*2-1; i++) 
			for (int i = 0; i < nstages*2-1; i++)
			{
				if (j == i+1) ss = matrix[j][i];
				if (j == i-1) dd = matrix[j][i];
			}
#if RS_CONTAIN
#if SEASONAL
			pSpecies->setSurv(hab,season,ii,1,(ss+dd));
			if ((ss+dd) > 0.0)
				pSpecies->setDev(hab,season,ii,1,dd/(ss+dd));
			else
				pSpecies->setDev(hab,season,ii,1,0.0);
#else
			pSpecies->setSurv(hab,ii,1,(ss+dd));
			if ((ss+dd) > 0.0)
				pSpecies->setDev(hab,ii,1,dd/(ss+dd));
			else
				pSpecies->setDev(hab,ii,1,0.0);
#endif // SEASONAL 
#else
#if SEASONAL
			pSpecies->setSurv(season,ii,1,(ss+dd));
			if ((ss+dd) > 0.0)
				pSpecies->setDev(season,ii,1,dd/(ss+dd));
			else
				pSpecies->setDev(season,ii,1,0.0);
#else
			pSpecies->setSurv(ii,1,(ss+dd));
			if ((ss+dd) > 0.0)
				pSpecies->setDev(ii,1,dd/(ss+dd));
			else
				pSpecies->setDev(ii,1,0.0);
#endif // SEASONAL 
#endif // RS_CONTAIN 
		}
		else{ // females
//			for (int i = 0; i < sstruct.nStages*2; i++) 
			for (int i = 0; i < nstages*2; i++) 
			{
				if (j == i+1) ss = matrix[j][i];
				if (j == i-1) dd = matrix[j][i];
			}
#if RS_CONTAIN
#if SEASONAL
			pSpecies->setSurv(hab,season,ii,0,(ss+dd));
			if ((ss+dd) > 0.0)
				pSpecies->setDev(hab,season,ii,0,dd/(ss+dd));
			else
				pSpecies->setDev(hab,season,ii,0,0.0);
#else
			pSpecies->setSurv(hab,ii,0,(ss+dd));
			if ((ss+dd) > 0.0)
				pSpecies->setDev(hab,ii,0,dd/(ss+dd));
			else
				pSpecies->setDev(hab,ii,0,0.0);
#endif // SEASONAL 
#else
#if SEASONAL
			pSpecies->setSurv(season,ii,0,(ss+dd));
			if ((ss+dd) > 0.0)
				pSpecies->setDev(season,ii,0,dd/(ss+dd));
			else
				pSpecies->setDev(season,ii,0,0.0);
#else
			pSpecies->setSurv(ii,0,(ss+dd));
			if ((ss+dd) > 0.0)
				pSpecies->setDev(ii,0,dd/(ss+dd));
			else
				pSpecies->setDev(ii,0,0.0);
#endif // SEASONAL 
#endif // RS_CONTAIN 
			ii++;
		}
	}
//	for (int j = 0; j < sstruct.nStages*2; j++) delete[] matrix[j];
//	delete[] matrix;
}
#if SEASONAL
pSpecies->setBreeding(season,breeding);
#endif // SEASONAL 

#if RSDEBUG
DEBUGLOG << "Read_TransitionMatrix(): matrix = " << matrix;
DEBUGLOG << " matrix[1][1] = " << matrix[1][1] << endl;
#endif

if (matrix != NULL) {
	for (int j = 0; j < matrixsize; j++) delete[] matrix[j];
	delete[] matrix;
	matrix = NULL; matrixsize = 0;
}

return 0;
}

//---------------------------------------------------------------------------
int ReadStageWeights(int option)
{
string header;
int i,j,n;
float f;
demogrParams dem = pSpecies->getDemogr();
stageParams sstruct = pSpecies->getStage();

if (dem.repType != 2) n = sstruct.nStages;
else n = sstruct.nStages*NSEXES;

#if RSDEBUG
DEBUGLOG << "Read_StageWeights(): option = " << option << " n = " << n << endl;
#endif
switch (option) {

case 1: { // fecundity
	// create stage weights matrix
	pSpecies->createDDwtFec(n);
#if RSDEBUG
DEBUGLOG << "Read_StageWeights(): completed fecundity matrix creation" << endl;
#endif
	for (i = 0; i < n+1; i++) fdfile >> header;
	// read coefficients
	for (i = 0; i < n; i++) {
		fdfile >> header;
		for (j = 0; j < n; j++) {
			fdfile >> f; pSpecies->setDDwtFec(j,i,f);
		}
	}
#if RSDEBUG
DEBUGLOG << "Read_StageWeights(): completed reading fecundity weights matrix " << endl;
#endif
	break;
}

case 2: { // development
	//create stage weights matrix
	pSpecies->createDDwtDev(n);
#if RSDEBUG
DEBUGLOG << "Read_StageWeights(): completed development matrix creation" << endl;
#endif
	for (i = 0; i < n+1; i++) ddfile >> header;
	//read coefficients
	for (i = 0; i < n; i++){
		ddfile >> header;
		for (j = 0; j < n; j++) {
			ddfile >> f; pSpecies->setDDwtDev(j,i,f);
		}
	}
#if RSDEBUG
DEBUGLOG << "Read_StageWeights(): completed reading development weights matrix " << endl;
#endif
	break;
}

case 3: { // sstruct.survival
	//create stage weights matrix
	pSpecies->createDDwtSurv(n);
#if RSDEBUG
DEBUGLOG << "Read_StageWeights(): completed survival matrix creation" << endl;
#endif
	for (i = 0; i < n+1; i++) sdfile >> header;
	//read coefficients
	for (i = 0; i < n; i++){
		sdfile >> header;
		for (j = 0; j < n; j++) {
			sdfile >> f; pSpecies->setDDwtSurv(j,i,f);
		}
	}
#if RSDEBUG
DEBUGLOG << "Read_StageWeights(): completed reading survival weights matrix " << endl;
#endif
	break;
}

}

return 0;
}

//---------------------------------------------------------------------------

#if SEASONAL

int ReadSeasonFile(const short option,const short nseasons) {
#if RSDEBUG
DEBUGLOG << "ReadSeasonFile(): option=" << option << " nseasons=" << nseasons << endl;
#endif
short season;        
string name;
//demogrParams dem = pSpecies->getDemogr();
stageParams sstruct = pSpecies->getStage();
string Inputs = paramsSim->getDir(1);

// read header line
string header;
int nheaders = 4;
//#if PARTMIGRN
//nheaders += 0;
//#endif // PARTMIGRN 
for (int i = 0; i < nheaders; i++) seasonfile >> header;

for (int i = 0; i < nseasons; i++) {
	seasonfile >> season;
	seasonfile >> name; // 'name' is TransMatrixFile
#if RSDEBUG
DEBUGLOG << "ReadSeasonFile(): i=" << i << " season=" << season << " name=" << name << endl;
#endif
	tmfile.open((Inputs+name).c_str());
	ReadTransitionMatrix(sstruct.nStages,sexesDem,0,season);    
	tmfile.close(); tmfile.clear();
	extrmevent e;
	seasonfile >> e.prob >> e.mort;
	pSpecies->setExtreme(season,e);
}

// ADD SURVIVAL SCHEDULE AND/OR STAGE WEIGHTS HERE?????

return 0;
}

int ReadExtremeFile(Landscape *pLandscape,const short nseasons) {
#if RSDEBUG
DEBUGLOG << "ReadExtremeFile(): nseasons=" << nseasons << endl;
#endif

// read header line
string header;
int nheaders;
if (patchmodel) nheaders = 4;	
else nheaders = 5;
for (int i = 0; i < nheaders; i++) extremefile >> header;

// Read data lines;
extEvent e;        
e.year = -98765;
extremefile >> e.year;
while (e.year != -98765) {
	extremefile >> e.season;
	if (patchmodel) {
		extremefile >> e.patchID; e.x = e.y = 0;
	}
	else {
		extremefile >> e.x >> e.y; e.patchID = 0;
	}
	extremefile >> e.probMort;
	pLandscape->addExtEvent(e);

	e.year = -98765;
	extremefile >> e.year;
	if (extremefile.eof()) e.year = -98765;
} // end of while loop

if (extremefile.is_open()) extremefile.close();
extremefile.clear();

return 0;
}

#endif // SEASONAL


//---------------------------------------------------------------------------
int ReadEmigration(int option)
{
int error = 0;

if (option == 0) { // open file and read header line
	emigFile.open(emigrationFile.c_str());
	string header;
#if GOBYMODEL
	for (int i = 0; i < 26; i++) emigFile >> header;
#else
#if GROUPDISP
	for (int i = 0; i < 27; i++) emigFile >> header;
#else
	for (int i = 0; i < 25; i++) emigFile >> header;
#endif // GROUPDISP
#endif // GOBYMODEL
	return 0;
}
if (option == 9) { // close file
	if (emigFile.is_open()) {
		emigFile.close(); emigFile.clear();
	}
	return 0;
}

int ffff,iiii,jjjj,kkkk,llll;
int Nlines,simulation,firstsimul = 0,stage,sex,emigstage;
#if GOBYMODEL
float asocD;
#endif
float	ep,d0,alpha,beta,epMean,epSD,d0Mean,d0SD,alphaMean,alphaSD,betaMean,betaSD;
float epScale,d0Scale,alphaScale,betaScale;
bool firstline = true;
demogrParams dem = pSpecies->getDemogr();
stageParams sstruct = pSpecies->getStage();
emigRules emig = pSpecies->getEmig();
emigTraits etraits{};
emigParams eparams{};

// set no.of lines assuming maximum stage- and sex-dependency
if (sstruct.nStages == 0) Nlines = sexesDisp;
else Nlines = sstruct.nStages * sexesDisp;

#if RSDEBUG
//DEBUGLOG << "ReadEmigration(): Nlines = " << Nlines << " dem.stageStruct = " << dem.stageStruct
//	<< " sstruct.nStages = " << sstruct.nStages << " sexesDisp = " << sexesDisp << endl;
#endif

for (int line = 0; line < Nlines; line++) {

	emigFile >> simulation >> iiii >> ffff >> jjjj >> kkkk >> llll >> emigstage;
	if (firstline) {
		firstsimul = simulation;
		if (iiii == 0) emig.densDep = false; else emig.densDep = true;
		if (jjjj == 0) emig.stgDep = false; else emig.stgDep = true;
		if (kkkk == 0) emig.sexDep = false; else emig.sexDep = true;
		if (llll == 0) emig.indVar = false; else emig.indVar = true;
		if (emigstage >= 0 && emigstage < sstruct.nStages) emig.emigStage = emigstage;
		else emig.emigStage = 0;
		// update no.of lines according to known stage- and sex-dependency
		if (emig.stgDep) {
			if (emig.sexDep) Nlines = sstruct.nStages * sexesDisp;
			else Nlines = sstruct.nStages;
		}
		else {
			if (emig.sexDep) Nlines = sexesDisp;
			else Nlines = 1;
		}
		if (ffff == 0) pSpecies->setFullKernel(false); else pSpecies->setFullKernel(true);
		pSpecies->setEmig(emig);
	}

#if RSDEBUG
//DEBUGLOG << "ReadEmigration(): Nlines=" << Nlines
//	<< " emig.densDep=" << emig.densDep
//	<< " emig.stgDep=" << emig.stgDep
//	<< " emig.sexDep=" << emig.sexDep
//	<< " emig.indVar=" << emig.indVar
//	<< endl;
#endif
	if (simulation != firstsimul) { // serious problem
		error = 300;
	}
	emigFile >> stage >> sex;

// ERROR MESSAGES SHOULD NEVER BE ACTIVATED ---------------------------------
#if RSDEBUG
//DEBUGLOG << "ReadEmigration(): line = " << line
//	<< " dem.stageStruct = " << dem.stageStruct
//	<< " emig.stgDep = " << emig.stgDep << endl;
#endif
if (dem.repType == 0) {
	if (emig.sexDep) error = 301;
}
if (dem.stageStruct) {
//	if (emig.indVar) error = 302;
}
else {
//	cout << endl << "***** pSpecies = " << pSpecies << endl << endl;
	if (emig.stgDep) error = 303;
}
//---------------------------------------------------------------------------

#if GROUPDISP
	int groupdisp;
	float groupmean;
	emigFile >> groupdisp >> groupmean;
	if (firstline) {
		if (groupdisp == 1) emig.groupdisp = true; else emig.groupdisp = false;
		emig.groupmean = groupmean;
		pSpecies->setEmig(emig);
	}
#endif
#if GOBYMODEL
	emigFile >> asocD;
	if (firstline) {
		emig.asocD = asocD;
		pSpecies->setEmig(emig);
	}
#endif
	emigFile >> ep >> d0 >> alpha >> beta >> epMean >> epSD >> d0Mean >> d0SD;
	emigFile >> alphaMean >> alphaSD >> betaMean >> betaSD;
	emigFile >> epScale >> d0Scale >> alphaScale >> betaScale;

	if (emig.sexDep) {
		if (emig.stgDep) {
			if (emig.densDep) {
				etraits.d0 = d0; etraits.alpha = alpha; etraits.beta = beta;
			}
			else {
				etraits.d0 = ep; etraits.alpha = etraits.beta = 0.0;
			}
			pSpecies->setEmigTraits(stage,sex,etraits);
		}
		else { // !emig.stgDep
			if (emig.indVar) {
				if (emig.densDep) {
					eparams.d0Mean = d0Mean; eparams.d0SD = d0SD;
					eparams.alphaMean = alphaMean; eparams.alphaSD = alphaSD;
					eparams.betaMean = betaMean; eparams.betaSD = betaSD;
				}
				else {
					eparams.d0Mean = epMean; eparams.d0SD = epSD;
					eparams.alphaMean = eparams.betaMean = 0.0;
					eparams.alphaSD = eparams.betaSD = 0.000001f;
				}
				pSpecies->setEmigParams(0,sex,eparams);
			}
			else { // !emig.indVar
				if (emig.densDep) {
					etraits.d0 = d0; etraits.alpha = alpha; etraits.beta = beta;
				}
				else {
					etraits.d0 = ep; etraits.alpha = etraits.beta = 0.0;
				}
				pSpecies->setEmigTraits(0,sex,etraits);
			}
		}
	}
	else { // !emig.sexDep
		if (emig.stgDep) {
			if (emig.densDep) {
				etraits.d0 = d0; etraits.alpha = alpha; etraits.beta = beta;
				pSpecies->setEmigTraits(stage,0,etraits);
			}
			else {
				etraits.d0 = ep; etraits.alpha = etraits.beta = 0.0;
				pSpecies->setEmigTraits(stage,0,etraits);
			}
		}
		else { // !emig.stgDep
			if (emig.densDep) {
				etraits.d0 = d0; etraits.alpha = alpha; etraits.beta = beta;
			}
			else {
				etraits.d0 = ep; etraits.alpha = etraits.beta = 0.0;
			}
			pSpecies->setEmigTraits(0,0,etraits);
#if RSDEBUG
//DEBUGLOG << "ReadEmigration(): case 0: emigP = " << ep << endl;
#endif
			if (emig.densDep) {
				eparams.d0Mean = d0Mean; eparams.d0SD = d0SD;
				eparams.alphaMean = alphaMean; eparams.alphaSD = alphaSD;
				eparams.betaMean = betaMean; eparams.betaSD = betaSD;
			}
			else {
				eparams.d0Mean = epMean; eparams.d0SD = epSD;
				eparams.alphaMean = eparams.betaMean = 0.0;
				eparams.alphaSD = eparams.betaSD = 0.000001f;
			}
			pSpecies->setEmigParams(0,0,eparams);
		}
	}

	if (stage == 0 && sex == 0) {
		emigScales scale = pSpecies->getEmigScales();
		if (emig.densDep) {
			scale.d0Scale = d0Scale; scale.alphaScale = alphaScale; scale.betaScale = betaScale;
		}
		else {
			scale.d0Scale = epScale; scale.alphaScale = scale.betaScale = 0.000001f;
		}
		pSpecies->setEmigScales(scale);
	}

	firstline = false;

} // end of Nlines for loop

return error;
}

//---------------------------------------------------------------------------
int ReadTransfer(int option, Landscape *pLandscape)
{
int iiii,jjjj,kkkk,Nlines,simulation,firstsimul = 0,stageDep,sexDep,stage,sex;
float tttt;
bool firstline = true;
int error = 0;
landParams paramsLand = pLandscape->getLandParams();
demogrParams dem = pSpecies->getDemogr();
stageParams sstruct = pSpecies->getStage();
trfrRules trfr = pSpecies->getTrfr();                   
#if RSDEBUG
DEBUGLOG << "ReadTransfer(): option=" << option
	<< " paramsLand.generated=" << paramsLand.generated
	<< " paramsLand.rasterType=" << paramsLand.rasterType
	<< " trfr.moveModel=" << trfr.moveModel
#if RS_CONTAIN
	<< " trfr.kernType=" << trfr.kernType
#else
	<< " trfr.twinKern=" << trfr.twinKern
#endif
	<< endl;
#endif

if (option == 0) { // open file and read header line
	transFile.open(transferFile.c_str());
	string header;
	int nheaders = 0;
	if (trfr.moveModel) {
#if RSDEBUG
DEBUGLOG << "ReadTransfer(): creating cost/mortality matrix, dimension=";
if (paramsLand.generated)
DEBUGLOG << paramsLand.nHab;
else
DEBUGLOG << paramsLand.nHabMax;
DEBUGLOG << endl;
#endif
		if (paramsLand.generated)
			pSpecies->createHabCostMort(paramsLand.nHab);
		else
			pSpecies->createHabCostMort(paramsLand.nHabMax);
		if (trfr.moveType == 1) { // SMS
//			int standardcols = 25;
			int standardcols = 23;
#if TEMPMORT
			standardcols += 1;
#endif
			if (paramsLand.generated) {
				nheaders = standardcols + 6; // artificial landscape
			}
			else { // real landscape
				if (paramsLand.rasterType == 0)
					nheaders = standardcols + 2 + 2*paramsLand.nHabMax; // habitat codes
				else nheaders = standardcols + 2; // habitat quality
			}
		}
		else { // CRW
			if (paramsLand.generated) {
				nheaders = 13;
			}
			else {
				if (paramsLand.rasterType == 0) nheaders = 13 + paramsLand.nHabMax;
				else nheaders = 13;
			}
		}
	}
	else { // dispersal kernel
#if RS_CONTAIN
		switch (trfr.kernType) {			
		case 0: // negative exponential kernel
		case 1: // double negative exponential kernel
			nheaders = 23; break;
		case 2: // 2Dt kernel
			nheaders = 10; break;
		case 3: // WALD kernel
			nheaders = 12 + sstruct.nStages; break;
		}
#else
		nheaders = 23;
#endif // RS_CONTAIN 
	}
#if RSDEBUG
DEBUGLOG << "ReadTransfer(): option=" << option << " nheaders=" << nheaders
	<< endl;
#endif
	for (int i = 0; i < nheaders; i++) transFile >> header;
	return 0;
}

if (option == 9) { // close file
	if (transFile.is_open()) {
		transFile.close(); transFile.clear();
	}
	return 0;
}

int TransferType; // new local variable to replace former global variable
#if RS_CONTAIN
if (trfr.moveModel) TransferType = trfr.moveType; 
else { 
	TransferType = 0; // negative exponential kernel
	if (trfr.kernType == 2) TransferType = 3; // 2Dt kernel
	if (trfr.kernType == 3) TransferType = 4; // WALD kernel
}
#else
if (trfr.moveModel) TransferType = trfr.moveType; else TransferType = 0;
#endif // RS_CONTAIN 

int sexKernels = 0;
trfrKernTraits k{};
trfrMovtTraits move{};
trfrKernParams kparams{};
trfrScales scale = pSpecies->getTrfrScales();
string CostsFile;
trfrSMSParams smsparams;
#if TEMPMORT
string MortFile;
#endif

switch (TransferType) {

case 0: // negative exponential dispersal kernel

	firstline = true;

	// set no.of lines assuming maximum stage- and sex-dependency
	if (sstruct.nStages == 0) Nlines = sexesDisp;
	else Nlines = sstruct.nStages * sexesDisp;

	for (int line = 0; line < Nlines; line++) {

#if RS_CONTAIN
		transFile >> simulation >> stageDep >> sexDep >> trfr.kernType >> jjjj >> kkkk;
#else
		transFile >> simulation >> stageDep >> sexDep >> iiii >> jjjj >> kkkk;
#endif // RS_CONTAIN 
		if (firstline) {
			firstsimul = simulation;
#if !RS_CONTAIN
			if (iiii == 0) trfr.twinKern = false; else trfr.twinKern = true;
#endif // RS_CONTAIN 
			if (jjjj == 0) trfr.distMort = false; else trfr.distMort = true;
			sexKernels = 2*stageDep + sexDep;
			if (kkkk == 0) trfr.indVar = false; else trfr.indVar = true;
			if (sexDep) trfr.sexDep = true; else trfr.sexDep = false;
			// update no.of lines according to known stage- and sex-dependency
			if (stageDep) {
				trfr.stgDep = true;
				if (sexDep) Nlines = sstruct.nStages * sexesDisp;
				else Nlines = sstruct.nStages;
			}
			else {
				trfr.stgDep = false;
				if (sexDep) Nlines = sexesDisp;
				else Nlines = 1;
			}
			pSpecies->setTrfr(trfr);
		}
		if (simulation != firstsimul) { // serious problem
			error = 400;
		}
		transFile >> stage >> sex;

		if (dem.repType == 0) {
			if (sexKernels == 1 || sexKernels == 3) error = 401;
		}
		if (dem.stageStruct) {
//			if (trfr.indVar) error = 402;
		}
		else{
			if (sexKernels == 2 || sexKernels == 3) error = 403;
		}
		if (sexKernels == 2 || sexKernels == 3) {
		}

		switch (sexKernels) {

		case 0: // no sex / stage dependence
			transFile >> k.meanDist1 >> k.meanDist2 >> k.probKern1;
			pSpecies->setKernTraits(0,0,k,(float)paramsLand.resol);
			if (trfr.indVar) {
				transFile >> kparams.dist1Mean >> kparams.dist1SD
					>> kparams.dist2Mean >> kparams.dist2SD
					>> kparams.PKern1Mean >> kparams.PKern1SD;
//				MAXDist = kparams.maxDist1;
				pSpecies->setKernParams(0,0,kparams,(float)paramsLand.resol);
			}
			else {
				for (int i = 0; i < 6; i++) transFile >> tttt;
			}
			break;

		case 1: // sex-dependent
			if (trfr.indVar)
			{
#if RS_CONTAIN
				if (trfr.kernType == 1) 
#else
				if (trfr.twinKern) 
#endif // RS_CONTAIN 
				{
					for (int i = 0; i < 3; i++) transFile >> tttt;
					transFile >> kparams.dist1Mean >> kparams.dist1SD
						>> kparams.dist2Mean >> kparams.dist2SD
						>> kparams.PKern1Mean >> kparams.PKern1SD;
				}
				else { // single kernel
					for (int i = 0; i < 3; i++) transFile >> tttt;
					transFile >> kparams.dist1Mean >> kparams.dist1SD;
					kparams.dist2Mean = kparams.dist1Mean; kparams.dist2SD = kparams.dist1SD;
					kparams.PKern1Mean = 0.999f; kparams.PKern1SD = 0.001f;
					for (int i = 0; i < 4; i++) transFile >> tttt;
				}
				pSpecies->setKernParams(0,sex,kparams,(float)paramsLand.resol);
			}
			else { // not varKernels
#if RS_CONTAIN
				if (trfr.kernType == 1) 
#else
				if (trfr.twinKern) 
#endif // RS_CONTAIN 
				{
					transFile >> k.meanDist1 >> k.meanDist2 >> k.probKern1;
					for (int i = 0; i < 6; i++) transFile >> tttt;
				}
				else {
					transFile >> k.meanDist1; k.meanDist2 = k.meanDist1; k.probKern1 = 1.0;
					for (int i = 0; i < 8; i++) transFile >> tttt;
				}
			pSpecies->setKernTraits(0,sex,k,(float)paramsLand.resol);
			}
			break;

		case 2: // stage-dependent
#if RS_CONTAIN
			if (trfr.kernType == 1) 
#else
			if (trfr.twinKern) 
#endif // RS_CONTAIN 
			{
				transFile >> k.meanDist1 >> k.meanDist2 >> k.probKern1;
				for (int i = 0; i < 6; i++) transFile >> tttt;
			}
			else {
				transFile >> k.meanDist1; k.meanDist2 = k.meanDist1; k.probKern1 = 1.0;
				for (int i = 0; i < 8; i++) transFile >> tttt;
			pSpecies->setKernTraits(stage,0,k,(float)paramsLand.resol);
			}
			break;

		case 3: // sex- & stage-dependent
#if RS_CONTAIN
			if (trfr.kernType == 1) 
#else
			if (trfr.twinKern) 
#endif // RS_CONTAIN 
			{
				transFile >> k.meanDist1 >> k.meanDist2 >> k.probKern1;
				for (int i = 0; i < 6; i++) transFile >> tttt;
			}
			else {
				transFile >> k.meanDist1; k.meanDist2 = k.meanDist1; k.probKern1 = 1.0;
				for (int i = 0; i < 8; i++) transFile >> tttt;
			pSpecies->setKernTraits(stage,sex,k,(float)paramsLand.resol);
			}
			break;
		} // end of switch (sexkernels)

		if (trfr.indVar) {
//			if (!trfr.indVar) error = 411;
//			if (dem.stageStruct) error = 412;
			if (stage == 0 && sex == 0) {
				transFile >> scale.dist1Scale >> scale.dist2Scale >> scale.PKern1Scale;
				pSpecies->setTrfrScales(scale);
			}
			else for (int i = 0; i < 3; i++) transFile >> tttt;
		}
		else for (int i = 0; i < 3; i++) transFile >> tttt;

		// mortality
		if (stage == 0 && sex == 0) {
			trfrMortParams mort{};
			transFile >> mort.fixedMort >> mort.mortAlpha >> mort.mortBeta;
			pSpecies->setMortParams(mort);
		}
		else for (int i = 0; i < 3; i++) transFile >> tttt;

		if (firstline) pSpecies->setTrfr(trfr);
		firstline = false;

	} // end of lines for loop

	break; // end of negative exponential dispersal kernel

case 1: // SMS

	transFile >> simulation >> iiii >> move.pr >> move.prMethod >> move.dp;
	if (iiii == 0) trfr.indVar = false; else trfr.indVar = true;
	transFile >> move.memSize >> move.gb >> move.goalType >> tttt >> iiii;
	if (move.goalType == 2) { // dispersal bias
		move.alphaDB = tttt; move.betaDB = iiii;
	}
#if RSDEBUG
DEBUGLOG << "ReadTransfer(): simulation=" << simulation << " indVar=" << trfr.indVar
	<< " PR=" << move.pr << " PRmethod=" << move.prMethod << endl;
DEBUGLOG << "ReadTransfer(): dp=" << move.dp << " MemSize=" << move.memSize
	<< " gb=" << move.gb << " goaltype=" << move.goalType << endl;
#endif
	smsparams = pSpecies->getSMSParams(0,0);
	transFile >> smsparams.dpMean >> smsparams.dpSD >> smsparams.gbMean >> smsparams.gbSD
		>> smsparams.alphaDBMean >> smsparams.alphaDBSD >> smsparams.betaDBMean >> smsparams.betaDBSD;
	transFile >> scale.dpScale >> scale.gbScale >> scale.alphaDBScale >> scale.betaDBScale;
	if (trfr.indVar) {
		pSpecies->setSMSParams(0,0,smsparams);
		pSpecies->setTrfrScales(scale);
	}
#if TEMPMORT
	transFile >> jjjj >> trfr.smType >> move.stepMort >> MortFile;  
#if RSDEBUG
DEBUGLOG << "ReadTransfer(): MortFile=" << MortFile << endl;
#endif
	if (trfr.smType == 2) {
		mortfilename = paramsSim->getDir(1) + MortFile;
		ReadMortalities(mortfilename);
	}
#else
	transFile >> jjjj >> iiii >> move.stepMort;
	if (iiii == 0) trfr.habMort = false; else trfr.habMort = true;
#endif // TEMPMORT 
	if (jjjj == 0) move.straigtenPath = false; else move.straigtenPath = true;

#if RSDEBUG
#if TEMPMORT
DEBUGLOG << "ReadTransfer(): SMtype=" << trfr.smType << " SMconst=" << move.stepMort << endl;
#else
DEBUGLOG << "ReadTransfer(): SMtype=" << trfr.habMort << " SMconst=" << move.stepMort << endl;
#endif // TEMPMORT 
#endif // RSDEBUG

	if (!paramsLand.generated) { // real landscape
		if (paramsLand.rasterType == 0) { // habitat codes
#if TEMPMORT
			if (trfr.smType == 1)
#else
			if (trfr.habMort)
#endif // TEMPMORT 
			{ // habitat-dependent step mortality
				for (int i = 0; i < paramsLand.nHabMax; i++)
				{
					transFile >> tttt;
					pSpecies->setHabMort(i,tttt);
#if RSDEBUG
//DEBUGLOG << "ReadTransfer(): nHabMax = " << paramsLand.nHabMax
//	<< " i = " << i << " mortality = " << tttt << endl;
#endif
				}
			}
			else { // constant step mortality
				for (int i = 0; i < paramsLand.nHabMax; i++) transFile >> tttt;
			}
		}
		else { // habitat quality
			// no columns to be read until CostMap
		}
	}
	else { // artificial landscape
#if TEMPMORT
		if (trfr.smType == 1)
#else
		if (trfr.habMort) 
#endif // TEMPMORT 
		{ // habitat-dependent step mortality
			// values are for habitat (hab=1) then for matrix (hab=0)
			transFile >> tttt; pSpecies->setHabMort(1,tttt);
			transFile >> tttt; pSpecies->setHabMort(0,tttt);
#if RSDEBUG
DEBUGLOG << "ReadTransfer(): nHab=" << paramsLand.nHab << endl;
DEBUGLOG << "ReadTransfer(): MortHabitat=" << pSpecies->getHabMort(1)
	<< " MortMatrix=" << pSpecies->getHabMort(0)
	<< endl;
#endif
		}
		else { // constant step mortality
			for (int i = 0; i < paramsLand.nHab; i++) transFile >> tttt;
		}
	}

	if (name_costfile != "NULL") trfr.costMap = true;
	else trfr.costMap = false;

	if (!paramsLand.generated) { // real landscape
		if (paramsLand.rasterType == 0) { // habitat codes
			if (trfr.costMap)
			{
				for (int i = 0; i < paramsLand.nHabMax; i++) transFile >> tttt;
			}
			else { // not costMap
				for (int i = 0; i < paramsLand.nHabMax; i++) {
					transFile >> tttt; iiii = (int)tttt;
					pSpecies->setHabCost(i,iiii);
#if RSDEBUG
DEBUGLOG << "ReadTransfer(): nHabMax=" << paramsLand.nHabMax << " i=" << i
	<< " tttt=" << tttt << " habCost[i]=" << pSpecies->getHabCost(i) << endl;
#endif
				}
			}
		}
		else { // habitat quality
			// no further columns to be read
		}
	}
	else { // artificial landscape
		if (trfr.costMap) // should not occur 
		{
			transFile >> tttt >> tttt;
		}
		else { // not costMap
			// costs are for habitat (hab=1) then for matrix (hab=0)
			transFile >> tttt; iiii = (int)tttt;
			pSpecies->setHabCost(1,iiii);
			transFile >> tttt; iiii = (int)tttt;
			pSpecies->setHabCost(0,iiii);
		}
	}
	pSpecies->setTrfr(trfr);
	pSpecies->setMovtTraits(move);

	break; // end of SMS

case 2: // CRW

{trfrCRWParams mparams{};

transFile >> simulation >> iiii;
if (iiii == 0) trfr.indVar = false; else trfr.indVar = true;

transFile >> move.stepLength >> move.rho
>> mparams.stepLgthMean >> mparams.stepLgthSD >> mparams.rhoMean >> mparams.rhoSD;
transFile >> scale.stepLScale >> scale.rhoScale;
#if TEMPMORT
transFile >> trfr.smType >> move.stepMort >> jjjj;
#else
transFile >> jjjj >> iiii >> move.stepMort;
if (iiii == 0) trfr.habMort = false; else trfr.habMort = true;
#endif // TEMPMORT 
if (jjjj == 0) move.straigtenPath = false; else move.straigtenPath = true;
pSpecies->setTrfrScales(scale);
#if RSDEBUG
DEBUGLOG << "ReadTransfer(): simulation=" << simulation
<< " paramsLand.rasterType=" << paramsLand.rasterType
<< " trfr.indVar=" << trfr.indVar
<< " move.stepLength=" << move.stepLength << " move.rho=" << move.rho
<< " mparams.stepLgthMean=" << mparams.stepLgthMean << " mparams.rhoMean=" << mparams.rhoMean
<< " move.straigtenPath=" << move.straigtenPath
<< endl;
#endif

//Habitat-dependent per step mortality
#if TEMPMORT
if (trfr.smType == 1 && paramsLand.rasterType != 0) error = 434;
#else
if (trfr.habMort && paramsLand.rasterType != 0) error = 434;
#endif // TEMPMORT 

if (!paramsLand.generated && paramsLand.rasterType == 0) { // real habitat codes landscape
#if TEMPMORT
	if (trfr.smType == 1)
#else
	if (trfr.habMort)
#endif // TEMPMORT 
	{ // habitat-dependent step mortality
		for (int i = 0; i < paramsLand.nHabMax; i++) {
			transFile >> tttt;
			pSpecies->setHabMort(i, tttt);
		}
	}
	else { // constant step mortality
		for (int i = 0; i < paramsLand.nHabMax; i++) transFile >> tttt;
	}
}
pSpecies->setTrfr(trfr);
pSpecies->setMovtTraits(move);
pSpecies->setCRWParams(0, 0, mparams);
}

break; // end of CRW


#if RS_CONTAIN

case 3: // 2Dt kernel

	trfr2Dt t2;

	transFile >> simulation >> iiii;
	if (iiii == 0) trfr.distMort = false; else trfr.distMort = true;

	transFile >> t2.u0Kernel1 >> t2.p0Kernel1 >> t2.u0Kernel2 >> t2.p0Kernel2 
		>> t2.propKernel1;
	pSpecies->setTrfr2Dt(t2);

	// mortality
	trfrMortParams mort;
	transFile >> mort.fixedMort >> mort.mortAlpha >> mort.mortBeta;
	pSpecies->setMortParams(mort);
	
#if RSDEBUG
DEBUGLOG << "ReadTransfer(): simulation=" << simulation
	<< " trfr.distMort=" << trfr.distMort
	<< " t2.u0Kernel1=" << t2.u0Kernel1 << " t2.p0Kernel1=" << t2.p0Kernel1
	<< " t2.u0Kernel2=" << t2.u0Kernel2 << " t2.p0Kernel2=" << t2.p0Kernel2
	<< " t2.propKernel1=" << t2.propKernel1 << " mort.fixedMort=" << mort.fixedMort
	<< endl;
#endif

	pSpecies->setTrfr(trfr);

	break; // end of 2Dt kernel

case 4: // WALD kernel

	trfrWald w;

	transFile >> simulation >> iiii;
	if (iiii == 0) trfr.distMort = false; else trfr.distMort = true;

	transFile >> w.meanU >> w.sigma_w >> w.hc;
	for (int i = 1; i < sstruct.nStages; i++) {
		float hr;
		transFile >> hr;
		pSpecies->setTrfrHr(hr,i);	
	}
	transFile >> w.vt >> w.kappa >> w.meanDirn >> w.sdDirn;
	pSpecies->setTrfrWald(w);

	// mortality
	trfrMortParams m;
	transFile >> m.fixedMort >> m.mortAlpha >> m.mortBeta;
	pSpecies->setMortParams(m);

#if RSDEBUG
DEBUGLOG << "ReadTransfer(): simulation=" << simulation
	<< " trfr.distMort=" << trfr.distMort
	<< " w.meanU=" << w.meanU << " w.sigma_w=" << w.sigma_w
	<< " w.hc=" << w.hc << " w.vt=" << w.vt
	<< " w.meanDirn=" << w.meanDirn << " w.sdDirn=" << w.sdDirn
	<< " m.fixedMort=" << m.fixedMort
	<< endl;
#endif
	
	pSpecies->setTrfr(trfr);

	break; // end of 2Dt kernel
	
#endif // RS_CONTAIN 

	
default:
	error = 440;
	break;
} // end of switch (TransferType)

return error;
}

#if TEMPMORT

//---------------------------------------------------------------------------
int ReadMortalities(string mortfile) {

//emigRules emig = pSpecies->getEmig();
//trfrRules trfr = pSpecies->getTrfr();
//settleType sett = pSpecies->getSettle();

int error = 0;

mortFile.open(mortfile.c_str());
pSpecies->clearMortalities();

int year;
float gradient;
#if RSDEBUG
//string msg = "No. of chromosomes set is " + Int2Str(nchromset);
//MessageDlg(msg.c_str(),mtWarning, TMsgDlgButtons() << mbOK,0);
#endif
year = -98765;
mortFile >> year;
while (year != -98765) {
	mortFile >> gradient;
#if RSDEBUG
DEBUGLOG << "ReadMortalities(): year=" << year
	<< " gradient=" << gradient 
	<< endl;
#endif
	pSpecies->addMortChange(year,gradient);
	year = -98765;
	mortFile >> year;
} ;

mortFile.close(); mortFile.clear();

return error;
}

#endif // TEMPMORT 

//---------------------------------------------------------------------------
// NOTE that stage- and sex-dependent settlement parameters are set for
// ALL stage/sex combinations, even if the species has stage- and/or
// sex-independent settlement rules
int ReadSettlement(int option)
{

int Nlines,simulation,firstsimul = 0,stageDep,sexDep,stage,sex;
bool firstline = true;
int error = 0;
demogrParams dem = pSpecies->getDemogr();
stageParams sstruct = pSpecies->getStage();
trfrRules trfr = pSpecies->getTrfr();
settleType sett = pSpecies->getSettle();
settleRules srules;
settleSteps ssteps{};
settleTraits settleDD{};
settParams sparams{};
int sexSettle = 0,settType = 0,densdep,indvar,findmate;
#if GOBYMODEL
float alphaSasoc,betaSasoc;
#endif

#if RSDEBUG
//DEBUGLOG << "ReadSettlement(): option=" << option << " transfer=" << transfer 
//	<< " trfr.moveModel=" << trfr.moveModel << endl;
#endif

if (option == 0) { // open file and read header line
	settFile.open(settleFile.c_str());
	string header;
	int nheaders = 0;
#if GOBYMODEL
	if (trfr.moveModel) nheaders = 25;
#else
	if (trfr.moveModel) nheaders = 23;
#endif
	else nheaders = 7;
	for (int i = 0; i < nheaders; i++) {
		settFile >> header;
#if RSDEBUG
//DEBUGLOG << "ReadSettlement(): i=" << i << " header=" << header << endl;
#endif
	}
	return 0;
}
if (option == 9) { // close file
	if (settFile.is_open()) {
		settFile.close(); settFile.clear();
	}
	return 0;
}

firstline = true;

// set no.of lines assuming maximum stage- and sex-dependency
if (sstruct.nStages == 0) Nlines = sexesDisp;
else Nlines = sstruct.nStages * sexesDisp;

for (int line = 0; line < Nlines; line++) {

	settFile >> simulation >> stageDep >> sexDep >> stage >> sex;
#if RS_CONTAIN
	if (transfer == 0 || transfer == 3 || transfer == 4) 
#else
	if (transfer == 0) 
#endif // RS_CONTAIN 
	{ // dispersal kernel
		settFile >> settType >> findmate;
	}
	else {
		settFile >> densdep >> indvar >> findmate;
	}
	if (firstline) {
		firstsimul = simulation;
		sexSettle = 2*stageDep + sexDep;
		if (stageDep == 1) sett.stgDep = true; else sett.stgDep = false;
		if (sexDep == 1) sett.sexDep = true; else sett.sexDep = false;
#if RS_CONTAIN
		if (transfer == 0 || transfer == 3 || transfer == 4) sett.indVar = false;  // dispersal kernel
#else
		if (transfer == 0) sett.indVar = false;  // dispersal kernel
#endif // RS_CONTAIN 
		else if (indvar == 1) sett.indVar = true; else sett.indVar = false;
		pSpecies->setSettle(sett);
		// update no.of lines according to known stage- and sex-dependency
		if (stageDep) {
			if (sexDep) Nlines = sstruct.nStages * sexesDisp;
			else Nlines = sstruct.nStages;
		}
		else {
			if (sexDep) Nlines = sexesDisp;
			else Nlines = 1;
		}
	}
	if (simulation != firstsimul) { // serious problem
		error = 500;
	}

	if (trfr.moveModel) {

		if (dem.repType == 0) {
			if (sexSettle == 1 || sexSettle == 3) error = 508;
		}
		if (!dem.stageStruct) {
			if (sexSettle == 2 || sexSettle == 3) error = 509;
		}

#if RSDEBUG
DEBUGLOG << "ReadSettlement(): sexSettle=" << sexSettle << endl;
#endif
		settFile >> ssteps.minSteps >> ssteps.maxSteps >> ssteps.maxStepsYr;
#if GOBYMODEL
		settFile >> alphaSasoc >> betaSasoc;
		if (firstline) {
			sett.alphaSasoc = alphaSasoc; sett.betaSasoc = betaSasoc;
			pSpecies->setSettle(sett);
		}
#endif
		settFile >> settleDD.s0 >> settleDD.alpha >> settleDD.beta;
		settFile >> sparams.s0Mean >> sparams.s0SD >> sparams.alphaSMean >> sparams.alphaSSD
			>> sparams.betaSMean >> sparams.betaSSD
			>> sparams.s0Scale >> sparams.alphaSScale >> sparams.betaSScale;

		switch (sexSettle) {

		case 0: // no sex- / stage-dependence
			srules = pSpecies->getSettRules(0,0);
			if (densdep == 1) srules.densDep = true; else srules.densDep = false;
			if (findmate == 1) srules.findMate = true; else srules.findMate = false;
			pSpecies->setSettRules(0,0,srules);
			pSpecies->setSteps(0,0,ssteps);
			if (srules.densDep) {
				if (sett.indVar) {
					pSpecies->setSettParams(0,0,sparams);
				}
				else {
					pSpecies->setSettTraits(0,0,settleDD);
				}
			}
			if (dem.stageStruct) { // model is structured - also set parameters for all stages
				for (int i = 1; i < sstruct.nStages; i++) {
					pSpecies->setSettRules(i,0,srules);
					pSpecies->setSteps(i,0,ssteps);
					pSpecies->setSettTraits(i,0,settleDD);
					if (dem.repType > 0) { // model is sexual - also set parameters for males
						pSpecies->setSettRules(i,1,srules);
						pSpecies->setSteps(i,1,ssteps);
						if (srules.densDep && !sett.indVar) pSpecies->setSettTraits(i,1,settleDD);
					}
				}
			}
			else {
				if (dem.repType > 0) { // model is sexual - also set parameters for males
					pSpecies->setSettRules(0,1,srules);
					pSpecies->setSteps(0,1,ssteps);
					if (srules.densDep) {
						if (sett.indVar) pSpecies->setSettParams(0,1,sparams);
						else pSpecies->setSettTraits(0,1,settleDD);
					}
				}
			}
			break;

		case 1: // sex-dependent
			srules = pSpecies->getSettRules(0,sex);
			if (densdep == 1) srules.densDep = true; else srules.densDep = false;
			if (findmate == 1) srules.findMate = true; else srules.findMate = false;
			pSpecies->setSettRules(0,sex,srules);
			pSpecies->setSteps(0,sex,ssteps);
#if RSDEBUG
DEBUGLOG << "ReadSettlement(): stage=" << stage << " sex=" << sex
	<< " sparams.s0Mean=" << sparams.s0Mean
	<< endl;
#endif
			if (srules.densDep) {
				if (sett.indVar) {
					pSpecies->setSettParams(0,sex,sparams);
				}
				else {
					pSpecies->setSettTraits(0,sex,settleDD);
				}
			}
			if (dem.stageStruct) { // model is structured - also set parameters for all stages
				for (int i = 1; i < sstruct.nStages; i++) {
					pSpecies->setSettRules(i,sex,srules);
					pSpecies->setSteps(i,sex,ssteps);
					if (srules.densDep && !sett.indVar) pSpecies->setSettTraits(i,sex,settleDD);
				}
			}
			break;

		case 2: // stage-dependent
			srules = pSpecies->getSettRules(stage,0);
			if (densdep == 1) srules.densDep = true; else srules.densDep = false;
			if (findmate == 1) srules.findMate = true; else srules.findMate = false;
			pSpecies->setSettRules(stage,0,srules);
			pSpecies->setSteps(stage,0,ssteps);
			if (srules.densDep) {
				if (sett.indVar) {
					if (stage == 0) pSpecies->setSettParams(0,0,sparams);
				}
				else {
					pSpecies->setSettTraits(stage,0,settleDD);
				}
			}
			if (dem.repType > 0) { // model is sexual - also set parameters for males
				pSpecies->setSettRules(stage,1,srules);
				pSpecies->setSteps(stage,1,ssteps);
				if (srules.densDep) {
					if (sett.indVar) {
						if (stage == 0) pSpecies->setSettParams(0,1,sparams);
					}
					else pSpecies->setSettTraits(stage,1,settleDD);
				}
			}
			break;

		case 3: // sex- & stage-dependent
			srules = pSpecies->getSettRules(stage,sex);
			if (densdep == 1) srules.densDep = true; else srules.densDep = false;
			if (findmate == 1) srules.findMate = true; else srules.findMate = false;
			pSpecies->setSettRules(stage,sex,srules);
			pSpecies->setSteps(stage,sex,ssteps);
			if (srules.densDep) {
				if (sett.indVar) {
					if (stage == 0) pSpecies->setSettParams(0,sex,sparams);
				}
				else {
					pSpecies->setSettTraits(stage,sex,settleDD);
				}
			}
			break;
		}

	} // end of movement model
	else { // dispersal kernel

		if (dem.repType == 0) {
			if (sett.sexDep) error = 501;
		}
		if (!dem.stageStruct) {
			if (sett.stgDep) error = 502;
		}

		switch (sexSettle) {

		case 0: //no sex / stage dependence
			if ((settType == 1 || settType == 3) && dem.stageStruct == false) error = 503;
			if (findmate == 1 && dem.repType == 0) error = 504;
			srules = pSpecies->getSettRules(0,0);
			switch (settType) {
			case 0:
				srules.wait = false; srules.go2nbrLocn = false;
				break;
			case 1:
				srules.wait = true; srules.go2nbrLocn = false;
				break;
			case 2:
				srules.wait = false; srules.go2nbrLocn = true;
				break;
			case 3:
				srules.wait = true; srules.go2nbrLocn = true;
				break;
			}
			if (findmate == 1) srules.findMate = true; else srules.findMate = false;
			pSpecies->setSettRules(0,0,srules);
			if (dem.stageStruct) { // model is structured - also set parameters for all stages
				for (int i = 0; i < sstruct.nStages; i++) {
					pSpecies->setSettRules(i,0,srules);
					if (dem.repType > 0) { // model is sexual - also set parameters for males
						pSpecies->setSettRules(i,1,srules);
					}
				}
			}
			else {
				if (dem.repType > 0) { // model is sexual - also set parameters for males
					pSpecies->setSettRules(0,1,srules);
				}
			}
			break;

		case 1: //sex dependent
			if ((settType == 1 || settType == 3) && dem.stageStruct == false) error = 505;
			srules = pSpecies->getSettRules(0,sex);
			switch (settType) {
			case 0:
				srules.wait = false; srules.go2nbrLocn = false;
				break;
			case 1:
				srules.wait = true; srules.go2nbrLocn = false;
				break;
			case 2:
				srules.wait = false; srules.go2nbrLocn = true;
				break;
			case 3:
				srules.wait = true; srules.go2nbrLocn = true;
				break;
			}
			if (findmate == 1) srules.findMate = true; else srules.findMate = false;
			pSpecies->setSettRules(0,sex,srules);
			if (dem.stageStruct) { // model is structured - also set parameters for all stages
				for (int i = 1; i < sstruct.nStages; i++) {
					pSpecies->setSettRules(i,sex,srules);
				}
			}
			break;

		case 2: //stage dependent
			if (findmate == 1 && dem.repType == 0) error = 507;
			srules = pSpecies->getSettRules(stage,0);
			switch (settType) {
			case 0:
				srules.wait = false; srules.go2nbrLocn = false;
				break;
			case 1:
				srules.wait = true; srules.go2nbrLocn = false;
				break;
			case 2:
				srules.wait = false; srules.go2nbrLocn = true;
				break;
			case 3:
				srules.wait = true; srules.go2nbrLocn = true;
				break;
			}
			if (findmate == 1) srules.findMate = true; else srules.findMate = false;
			pSpecies->setSettRules(stage,0,srules);
			if (dem.repType > 0) { // model is sexual - also set parameters for males
				pSpecies->setSettRules(stage,1,srules);
			}
			break;

		case 3: //sex & stage dependent
			srules = pSpecies->getSettRules(stage,sex);
			switch (settType) {
			case 0:
				srules.wait = false; srules.go2nbrLocn = false;
				break;
			case 1:
				srules.wait = true; srules.go2nbrLocn = false;
				break;
			case 2:
				srules.wait = false; srules.go2nbrLocn = true;
				break;
			case 3:
				srules.wait = true; srules.go2nbrLocn = true;
				break;
			}
			if (findmate == 1) srules.findMate = true; else srules.findMate = false;
			pSpecies->setSettRules(stage,sex,srules);
			break;

		} // end of switch (sexSettle)

	} // end of dispersal kernel

	firstline = false;

} // end of for line loop

return error;
}

//---------------------------------------------------------------------------
int ReadGenetics(int option)
{
#if RSDEBUG
DEBUGLOG << "ReadGenetics(): option=" << option  << " geneticsFile= " << geneticsFile
	<< endl;
#endif
emigRules emig = pSpecies->getEmig();
trfrRules trfr = pSpecies->getTrfr();
settleType sett = pSpecies->getSettle();

int simulation,arch;
string archfile;
int error = 0;
genomeData g{};
demogrParams dem = pSpecies->getDemogr();

if (option == 0) { // open file and read header line
	genFile.open(geneticsFile.c_str());
	string header;
	int nheaders = 8;
	for (int i = 0; i < nheaders; i++) genFile >> header;
	return 0;
}

if (option == 9) { // close file
	if (genFile.is_open()) {
		genFile.close(); genFile.clear();
	}
	return 0;
}

genFile >> simulation >> arch >> g.nLoci >> archfile >> g.probMutn >> g.probCrossover
	>> g.alleleSD >> g.mutationSD;
if (dem.repType == 0) g.diploid = false; else g.diploid = true;
#if RSDEBUG
DEBUGLOG << "ReadGenetics(): simulation=" << simulation << " arch=" << arch
	<< " g.nLoci=" << g.nLoci << " archfile=" << archfile
	<< " g.probMutn=" << g.probMutn << " g.probCrossover=" << g.probCrossover
	<< " g.alleleSD=" << g.alleleSD << " g.mutationSD=" << g.mutationSD
	<< endl;
#endif

g.neutralMarkers = false;
if (arch == 0) { // no architecture file
	g.trait1Chromosome = true;
	pSpecies->set1ChromPerTrait(g.nLoci);
}
else { // architecture file
	g.trait1Chromosome = false;
	g.nLoci = 0;
	if (!(emig.indVar || trfr.indVar || sett.indVar)) {
		g.neutralMarkers = true;
	}
	error = ReadArchFile(archfile);
}

pSpecies->setGenomeData(g);

return error;
}

//---------------------------------------------------------------------------
int ReadArchFile(string archfile) {

emigRules emig = pSpecies->getEmig();
trfrRules trfr = pSpecies->getTrfr();
settleType sett = pSpecies->getSettle();

int error = 0;
string Inputs = paramsSim->getDir(1);
genfilename = Inputs + archfile;

// in order to determine no. of defined traits in architecture file, parse it again,
// which sets the variable fileNtraits correctly
bArchFile.open(genfilename.c_str());
ParseArchFile();
bArchFile.close(); bArchFile.clear();

// re-open the file to read data and set up trait maps
archFile.open(genfilename.c_str());

string paramname;
int nchromosomes,nloci,traitnum,chrom,locus;
// set no. of chromosomes
archFile >> paramname >> nchromosomes;
pSpecies->setNChromosomes(nchromosomes);
int nchromset = pSpecies->getNChromosomes();
#if VCL
#if RSDEBUG
//string msg = "No. of chromosomes set is " + Int2Str(nchromset);
//MessageDlg(msg.c_str(),mtWarning, TMsgDlgButtons() << mbOK,0);
#endif
#endif
if (nchromset <= 0) error = 1;
if (emig.indVar || trfr.indVar || sett.indVar) {
	pSpecies->setTraitData(fileNtraits);
}
else { // neutral markers only
	pSpecies->setTraitData(0);
}
// set no. of loci for each chromosome
archFile >> paramname;
for (int i = 0; i < nchromosomes; i++) {
	archFile >> nloci;
	pSpecies->setNLoci(i,nloci);
}
if (emig.indVar || trfr.indVar || sett.indVar) {
	// set trait maps
	paramname = "XXXyyyZZZ";
	archFile >> paramname;
	while (paramname != "XXXyyyZZZ") {
		archFile >> traitnum >> paramname >> nloci;
		pSpecies->setTraitMap(traitnum,nloci);
		for (int allele = 0; allele < nloci; allele++) {
			chrom = locus = -999999;
			archFile >> chrom >> locus;
			pSpecies->setTraitAllele(traitnum,allele,chrom,locus);
		}
		paramname = "XXXyyyZZZ";
		archFile >> paramname;
	} ;
}

archFile.close(); archFile.clear();

// any loci not contributing to a trait are recorded as neutral
if (emig.indVar || trfr.indVar || sett.indVar) {
	pSpecies->setNeutralLoci(false);
}
else { // model has neutral markers only
	pSpecies->setNeutralLoci(true);
}

return error;
}

#if RS_CONTAIN

//---------------------------------------------------------------------------
int ReadManageFile(int option,Landscape *pLandscape)
{
#if RSDEBUG
DEBUGLOG << "ReadManageFile(): option=" << option  << " name_managefile= " << name_managefile
	<< endl;
#endif

int error = 0;
demogrParams dem = pSpecies->getDemogr();
stageParams sstruct = pSpecies->getStage();

if (option == 0) { // open file and read header line
	amfile.open(name_managefile.c_str());
	string header;
	int nheaders = 19;
//	if (dem.stageStruct) nheaders += sstruct.nStages;
	if (dem.stageStruct) nheaders += 2;
	for (int i = 0; i < nheaders; i++) amfile >> header;
	return 0;
}

if (option == 9) { // close file
	if (amfile.is_open()) {
		amfile.close(); amfile.clear();
	}
	return 0;
}

int simulation,Nlines,line,cullNstages,cullstage,iiii;
culldata c;
cullstagedata cstage;       
float distdecay,cullmaxrate,cullalpha,cullbeta;
damageparams d;

Nlines = line = 1; // assume one line is to be read until cullNstages is read
for (int i = 0; i < Nlines; i++) {

	amfile >> simulation;

	//	amfile >> c.method >> c.timing >> c.popnThreshold >> c.cullMaxRate >> c.maxNpatches;
	amfile >> c.method >> c.timing >> distdecay >> c.maxNpatches >> c.densThreshold >> c.countCV;
	amfile >> cullNstages >> cullstage >> c.cullRate >> cullmaxrate >> cullalpha >> cullbeta; 
	if (c.cullRate == 2) {
		cstage.cullMaxRate = cullmaxrate; cstage.cullAlpha = cullalpha; cstage.cullBeta = cullbeta;
	}
	else {
		c.cullMaxRate = cullmaxrate; c.cullAlpha = cullalpha; c.cullBeta = cullbeta;
	}      
//	amfile >> iiii;
//	if (iiii == 1) c.edgeBias = true; else c.edgeBias = false;
	if (line == 1) {
		if (dem.stageStruct) pCull->resetCullStage();
		pCull->setCullData(c);
		if (c.method > 1) {
			landParams ppLand = pLandscape->getLandParams();
			if (!ppLand.dmgLoaded) error = 701;
		}
		pLandscape->setAlpha(distdecay);	
	}
	pCull->setCullStage(cullstage,true);
	if (c.cullRate == 2) {
		pCull->setCullStageData(cstage,cullstage); 		
	}
#if RSDEBUG
//DEBUGLOG << "ReadManageFile(): simulation=" << simulation << " alpha=" << alpha
//	<< endl;
#endif

//	if (dem.stageStruct) {
//		pCull->resetCullStage();
//		int cull;
//		for (int i = 0; i < sstruct.nStages; i++) {
//			amfile >> cull;
//			if (cull == 0) pCull->setCullStage(i,false); else pCull->setCullStage(i,true);
//		}
//	}

	amfile >> d.timing >> d.type >> d.occOption >> d.stage >> d.alphaOccupancy >> d.betaOccupancy 
		>> d.alphaTraversal >> d.betaTraversal;
	if (line == 1) pDamageParams->setDamageParams(d);

	if (line == 1) {
		if (cullNstages > 1) Nlines = cullNstages;		
	}	 
	line++;
	
}

return error;
}

#endif // RS_CONTAIN 

//---------------------------------------------------------------------------
int ReadInitialisation(int option, Landscape *pLandscape)
{
landParams paramsLand = pLandscape->getLandParams();
demogrParams dem = pSpecies->getDemogr();
stageParams sstruct = pSpecies->getStage();
initParams init = paramsInit->getInit();
string Inputs = paramsSim->getDir(1);

int simulation,maxcells;
float check;
int error = 0;

if (option == 0) { // open file and read header line
	initFile.open(initialFile.c_str());
	string header;
	int nheaders = 17;
	if (dem.stageStruct) nheaders += 1 + sstruct.nStages-1;
	for (int i = 0; i < nheaders; i++) initFile >> header;
	return 0;
}

if (option == 9) { // close file
	if (initFile.is_open()) {
		initFile.close(); initFile.clear();
	}
	return 0;
}

initFile >> simulation >> init.seedType >> init.freeType >> init.spDistType;
if (init.seedType == 1 && paramsLand.spDist == false) error = 601;
if (paramsLand.patchModel) initFile >> init.initDens >> init.indsHa;
else initFile >> init.initDens >> init.indsCell;
initFile >> init.minSeedX >> init.maxSeedX >> init.minSeedY >> init.maxSeedY
	>> init.nSeedPatches >> init.nSpDistPatches
	>> init.initFrzYr >> init.restrictRows >> init.restrictFreq >> init.finalFrzYr
	>> init.indsFile;
#if RSDEBUG
DEBUGLOG << "ReadInitialisation(): simulation=" << simulation
	<< " seedType=" << init.seedType << " freeType=" << init.freeType
	<< " spDistType=" << init.spDistType
	<< " maxSeedX=" << init.maxSeedX << " maxSeedY=" << init.maxSeedY
	<< " initFrzYr=" << init.initFrzYr << " restrictRows=" << init.restrictRows
	<< " restrictFreq=" << init.restrictFreq << " finalFrzYr=" << init.finalFrzYr
	<< " indsFile=" << init.indsFile
	<< endl;
#endif
init.restrictRange = false;
if (init.seedType == 0 && init.restrictRows > 0) init.restrictRange = true;

if (dem.stageStruct) {
	float p;
	initFile >> init.initAge;
	check = 0.0;
	for (int i = 1; i < sstruct.nStages; i++) {
	 initFile >> p;
	 check += p;
	 paramsInit->setProp(i,p);
	}
	if (check < 1.0 || check > 1.0)
	{ // this condition should not occur - WHAT COULD BE DONE?? ABORT WITH ERROR CODE ...
#if RSDEBUG
DEBUGLOG << "ReadInitialisation(): check = " << check << endl;
#endif
	}
}

paramsInit->setInit(init);
switch (init.seedType) {
case 0: // // free initialisation
	if (init.minSeedX < 0) init.minSeedX = 0;
	if (init.minSeedY < 0) init.minSeedY = 0;
	if (init.maxSeedX < 0) init.maxSeedX = paramsLand.maxX;
	if (init.maxSeedY < 0) init.maxSeedY = paramsLand.maxY;
	if (init.minSeedY > init.maxSeedY || init.minSeedX > init.maxSeedX) {
#if RSDEBUG
DEBUGLOG << "ReadInitialisation(): maxSeedX=" << init.maxSeedX
	<< " paramsLand.maxX=" << paramsLand.maxX
	<< " maxSeedY=" << init.maxSeedY
	<< " paramsLand.maxY=" << paramsLand.maxY << endl;
#endif
	error = 603;
	}
	maxcells = (init.maxSeedY - init.minSeedY)*(init.maxSeedX - init.minSeedX);
	if (init.freeType == 0 && init.nSeedPatches > maxcells) error = 602;
	break;
case 1: // from species distribution
	break;
case 2: // from initial individuals file
	if (init.indsFile != prevInitialIndsFile) {
		// read and store the list of individuals to be initialised
		ReadInitIndsFile(0,pLandscape,(Inputs + init.indsFile));
		prevInitialIndsFile = init.indsFile;
	}
	break;
default:
	break;
	;
}

return error;
}

//---------------------------------------------------------------------------
int ReadInitIndsFile(int option,Landscape *pLandscape,string indsfile) {
string header;
landParams paramsLand = pLandscape->getLandParams();
demogrParams dem = pSpecies->getDemogr();
//stageParams sstruct = pSpecies->getStage();
initParams init = paramsInit->getInit();

if (option == 0) { // open file and read header line
	initIndsFile.open(indsfile.c_str());
	string header;
	int nheaders = 3;
	if (paramsLand.patchModel) nheaders++;
	else nheaders += 2;
	if (dem.repType > 0) nheaders++;
	if (dem.stageStruct) nheaders  += 2;
	for (int i = 0; i < nheaders; i++) initIndsFile >> header;
	paramsInit->resetInitInds();
//	return 0;
}

if (option == 9) { // close file
	if (initIndsFile.is_open()) {
		initIndsFile.close(); initIndsFile.clear();
	}
	return 0;
}

// Read data lines;
initInd iind{};
int ninds;
int totinds = 0;
iind.year = -98765;
initIndsFile >> iind.year;
while (iind.year != -98765) {
	initIndsFile >> iind.species;
	if (paramsLand.patchModel) {
		initIndsFile >> iind.patchID; iind.x = iind.y = 0;
	}
	else {
		initIndsFile >> iind.x >> iind.y; iind.patchID = 0;
	}
	initIndsFile >> ninds;
	if (dem.repType > 0) initIndsFile >> iind.sex;
	else iind.sex = 0;
	if (dem.stageStruct) {
		initIndsFile >> iind.age >> iind.stage;
	}
	else {
		iind.age = iind.stage = 0;
	}
	for (int i = 0; i < ninds; i++) {
		totinds++;
		paramsInit->addInitInd(iind);
	}

	iind.year = -98765;
	initIndsFile >> iind.year;
	if (initIndsFile.eof()) iind.year = -98765;
} // end of while loop

if (initIndsFile.is_open()) initIndsFile.close();
initIndsFile.clear();

return totinds;
}

#if VIRTUALECOLOGIST

//---------------------------------------------------------------------------
int ReadVirtEcol(int option)
{
#if RSDEBUG
DEBUGLOG << "ReadVirtEcol(): option=" << option  << " virtEcolFile= " << virtEcolFile
	<< endl;
#endif

int simulation;
string samplefile,patchfile;
int error = 0;

virtParams virt = paramsSim->getVirtParams();
sampleLimits lim;

if (option == 0) { // open file and read header line
	bVirtEcolFile.open(virtEcolFile.c_str());
	string header;
	int nheaders = 19;
	for (int i = 0; i < nheaders; i++) bVirtEcolFile >> header;
	return 0;
}

if (option == 9) { // close file
	if (bVirtEcolFile.is_open()) {
		bVirtEcolFile.close(); bVirtEcolFile.clear();
	}
	return 0;
}

int landgen,outgen,sampleloci,minX,maxX,minY,maxY;
bVirtEcolFile >> simulation >> virt.patchMethod >> virt.rowsFront
		>> virt.maxNPatches >> virt.maxPatchNum >> patchfile
		>> virt.minIndsPatch >> virt.maxIndsPatch >> virt.stgMethod
		>> lim.minX >> lim.maxX >> lim.minY >> lim.maxY
		>> landgen >> virt.outStart >> virt.outInt >> outgen >> sampleloci >> samplefile;
if (landgen == 1) virt.landscapeGenetics = true; else virt.landscapeGenetics = false;
if (outgen == 1) virt.outGenomes = true; else virt.outGenomes = false;
if (sampleloci == 1) pSpecies->setSampleAllLoci(false); else pSpecies->setSampleAllLoci(true);

#if RSDEBUG
DEBUGLOG << "ReadVirtEcol(): simulation=" << simulation << " patchMethod=" << virt.patchMethod
	<< " maxNPatches=" << virt.maxNPatches << " maxPatchNum=" << virt.maxPatchNum 
	<< " patchfile=" << patchfile << " maxIndsPatch=" << virt.maxIndsPatch
	<< " landscapeGenetics=" << (int)virt.landscapeGenetics << " samplefile=" << samplefile
	<< endl;
#endif
paramsSim->setVirtParams(virt);
if (virt.patchMethod == 1) paramsSim->setLimits(lim);

if (virt.patchMethod == 3) { // read sample patch numbers from file
	error = ReadPatchFile(patchfile);
}

if (virt.landscapeGenetics) {
	if (samplefile != "NULL") { // sample loci from file
		error = ReadSampleFile(samplefile);
	}
}

return error;
}

//---------------------------------------------------------------------------
int ReadSampleFile(string samplefile) {
#if RSDEBUG
DEBUGLOG << "ReadSampleFile(): samplefile=" << samplefile
	<< endl;
#endif

int error = 0;
string Inputs = paramsSim->getDir(1);
locfilename = Inputs + samplefile;

bSampleFile.open(locfilename.c_str());

string paramname;
int nsamples,chromosome,locus;
bSampleFile >> paramname >> nsamples;
#if RSDEBUG
DEBUGLOG << "ReadSampleFile(): paramname=" << paramname << " nsamples=" << nsamples
	<< endl;
#endif

pSpecies->resetSampleLoci();
bSampleFile >> paramname;
int nchromosomes = pSpecies->getNChromosomes();
int nloci;
bool sampleOK = true;
for (int i = 0; i < nsamples; i++) {
	chromosome = locus = -999;
	bSampleFile >> chromosome >> locus;
	if (chromosome < 0 || chromosome >= nchromosomes || locus < 0) sampleOK = false;
	else {
		// check specified locus exists
		nloci = pSpecies->getNLoci(chromosome);
		if (locus >= nloci) sampleOK = false;
	}
	if (sampleOK) {
		pSpecies->addSampleLocus(chromosome,locus);
	}
	else {
		error = 701;
		i = nsamples;
		break;
	}
}

bSampleFile.close(); bSampleFile.clear();
#if RSDEBUG
DEBUGLOG << "ReadSampleFile(): error=" << error
	<< endl;
#endif

return error;
}

//---------------------------------------------------------------------------
int ReadPatchFile(string patchfile) {
#if RSDEBUG
DEBUGLOG << "ReadPatchFile(): patchfile=" << patchfile
	<< endl;
#endif

int error = 0;
string Inputs = paramsSim->getDir(1);
patchfilename = Inputs + patchfile;

bPatchFile.open(patchfilename.c_str());

int nsamples,patchnum;

paramsSim->clearSamplePatches();
patchnum = -98765;
bPatchFile >> patchnum;
#if RSDEBUG
DEBUGLOG << "ReadPatchFile(): FIRST patchnum=" << patchnum 
	<< endl;
#endif
while (patchnum != -98765) {
	if (patchnum > 0) paramsSim->addSamplePatch(patchnum);
	else error++;
	patchnum = -98765;
	bPatchFile >> patchnum;
}

bPatchFile.close(); bPatchFile.clear();
#if RSDEBUG
DEBUGLOG << "ReadPatchFile(): error=" << error
	<< endl;
#endif

return error;
}

#endif // VIRTUALECOLOGIST

//---------------------------------------------------------------------------
void RunBatch(int nSimuls, int nLandscapes)
{
int land_nr;
int t0,t1,t00,t01;
int read_error;
bool params_ok;
simParams sim = paramsSim->getSim();

#if SEASONAL
demogrParams dem = pSpecies->getDemogr();
#endif

Landscape *pLandscape = NULL;  		// pointer to landscape

#if RSDEBUG
DEBUGLOG << endl;
DEBUGLOG << "RunBatch(): nSimuls=" << nSimuls << " nLandscapes=" << nLandscapes << endl;
DEBUGLOG << "RunBatch(): landtype=" << landtype << " maxNhab=" << maxNhab << endl;
#endif

t0 = (int)time(0);

//int batch_line = 0;

string name = paramsSim->getDir(2) + "Batch" + Int2Str(sim.batchNum) + "_RS_log.csv";
if (rsLog.is_open()) {
	rsLog.close(); rsLog.clear();
}
rsLog.open(name.c_str());
if (!rsLog.is_open()) {
	cout << endl << "Error - unable to open Batch" << sim.batchNum 
		<< "_RS_log.csv file - aborting batch run" << endl;
	return;
}
rsLog << "Event,Number,Reps,Years,Time" << endl;
#if RSDEBUG
rsLog << "WARNING,***** RSDEBUG mode is active *****,,," << endl;
#endif
#if !RS_EMBARCADERO || !LINUX_CLUSTER
#if GROUPDISP || RS_ABC
rsLog << "RANDOM SEED,,,," << RS_random_seed << endl;
#else
rsLog << "RANDOM SEED," << RS_random_seed << ",,," << endl;
#endif
#endif

// Open landscape batch file and read header record
if (ReadLandFile(0)) {
	cout << endl << "Error opening landFile - aborting batch run" << endl;
	return;
}

for (int j = 0; j < nLandscapes; j++) {
#if RSDEBUG
DEBUGLOG << endl;
#endif
	// create new landscape
	if (pLandscape != NULL) delete pLandscape;
	pLandscape = new Landscape;
	bool landOK = true;

	t00 = (int)time(0);
	land_nr = ReadLandFile(1,pLandscape);
	if (land_nr <= 0) { // error condition
		string msg = "Error code " + Int2Str(-land_nr)
			+ " returned from reading LandFile - aborting batch run";
		cout << endl << msg << endl;
		MemoLine(msg.c_str());
		ReadLandFile(9); // close the landscape file
		return;
	}

	MemoLine(("Starting landscape " + Int2Str(land_nr) + "...").c_str());

#if RSDEBUG
DEBUGLOG << endl << "RunBatch(): j=" << j << " land_nr=" << land_nr
	<< " landtype=" << landtype;
if (landtype != 9)
	DEBUGLOG << " name_landscape=" << name_landscape
		<< " name_patch=" << name_patch 
		<< " name_costfile=" << name_costfile 
		<< " name_sp_dist=" << name_sp_dist;
DEBUGLOG << endl;
#endif
	landParams paramsLand = pLandscape->getLandParams();
	paramsLand.patchModel = patchmodel;
	paramsLand.resol = resolution;
	paramsLand.rasterType = landtype;
	if (landtype == 9) {
		paramsLand.generated = true;
		paramsLand.nHab = 2;
	}
	else {
		paramsLand.generated = false;
		if (name_dynland == "NULL") paramsLand.dynamic = false;
		else paramsLand.dynamic = true;
	}
	paramsLand.nHabMax = maxNhab;
	paramsLand.spDist = speciesdist;
	paramsLand.spResol = distresolution;
	pLandscape->setLandParams(paramsLand,sim.batchMode);

	if (landtype != 9) { // imported landscape
		string hname = paramsSim->getDir(1) + name_landscape;
		int landcode;
		string cname;
		if (name_costfile == "NULL" || name_costfile == "none") cname = "NULL";
		else cname = paramsSim->getDir(1) + name_costfile;
#if RS_CONTAIN
		string dname;
		if (name_damagefile == "NULL") dname = name_damagefile;
		else dname = paramsSim->getDir(1) + name_damagefile;
#endif // RS_CONTAIN 
		if (paramsLand.patchModel) {
			string pname = paramsSim->getDir(1) + name_patch;
#if RS_CONTAIN
#if SEASONAL
			landcode = pLandscape->readLandscape(nseasons,0,hname,pname,cname,dname);
#else
			landcode = pLandscape->readLandscape(0,hname,pname,cname,dname);
#endif // SEASONAL 
#else
#if SEASONAL
			landcode = pLandscape->readLandscape(nseasons,0,hname,pname,cname);
#else
#if RSDEBUG
int t02a = time(0);
#endif
			landcode = pLandscape->readLandscape(0,hname,pname,cname);
#if RSDEBUG
int t02b = time(0);
DEBUGLOG << "RunBatch(): TIME for readLandscape() " << t02b-t02a << endl;
#endif
#endif // SEASONAL 
#endif // RS_CONTAIN 
		}
		else {
#if RS_CONTAIN
#if SEASONAL
			landcode = pLandscape->readLandscape(nseasons,0,hname," ",cname,dname);
#else
			landcode = pLandscape->readLandscape(0,hname," ",cname,dname);
#endif // SEASONAL 
#else
#if SEASONAL
			landcode = pLandscape->readLandscape(nseasons,0,hname," ",cname);
#else
			landcode = pLandscape->readLandscape(0,hname," ",cname);
#endif // SEASONAL 
#endif // RS_CONTAIN 
    }
		if (landcode != 0) {
			rsLog << "Landscape," << land_nr << ",ERROR,CODE," << landcode << endl;
			cout << endl << "Error reading landscape " << land_nr << " - aborting" << endl;
			landOK = false;
		}
		if (paramsLand.dynamic) {
#if RSDEBUG
int t03a = time(0);
#endif
			landcode = ReadDynLandFile(pLandscape);
#if RSDEBUG
int t03b = time(0);
DEBUGLOG << "RunBatch(): TIME for ReadDynLandFile() " << t03b-t03a << endl;
#endif
			if (landcode != 0) {
				rsLog << "Landscape," << land_nr << ",ERROR,CODE," << landcode << endl;
				cout << endl << "Error reading landscape " << land_nr << " - aborting" << endl;
				landOK = false;
			}
		}
		if (landtype == 0) {
			pLandscape->updateHabitatIndices();
		}
#if RSDEBUG
landParams tempLand = pLandscape->getLandParams();
DEBUGLOG << "RunBatch(): j=" << j
	<< " land_nr=" << land_nr
	<< " landcode=" << landcode
	<< " nHab=" << tempLand.nHab
	<< endl;
#endif

		// species distribution
															 
		if (paramsLand.spDist) { // read initial species distribution
			// WILL NEED TO BE CHANGED FOR MULTIPLE SPECIES ...
			string distname = paramsSim->getDir(1) + name_sp_dist;
			landcode = pLandscape->newDistribution(pSpecies,distname);
			if (landcode == 0) {
			}
			else {
				rsLog << "Landscape," << land_nr << ",ERROR,CODE," << landcode << endl;
				cout << endl << "Error reading initial distribution for landscape "
					<< land_nr << " - aborting" << endl;
				landOK = false;
			}
		}
		paramsSim->setSim(sim);
#if RSDEBUG
DEBUGLOG << "RunBatch(): j=" << j
	<< " spDist=" << paramsLand.spDist
	<< endl;
#endif

#if SPATIALMORT
		// optional spatial mortality files
		if (name_mortfile[0] != "NULL") {
			string mname[2];
			mname[0] = paramsSim->getDir(1) + name_mortfile[0];
			mname[1] = paramsSim->getDir(1) + name_mortfile[1];
			landcode = pLandscape->readMortalityFiles(mname[0],mname[1]);
			if (landcode != 0) {
				cout << endl << "Error reading mortality files for landscape "
					<< land_nr << " - aborting" << endl;
				rsLog << "Landscape," << land_nr << ",ERROR,CODE," << landcode << endl;
				landOK = false;
			}
		}
#endif

		if (landOK) {
			t01 = (int)time(0);
			rsLog << "Landscape," << land_nr << ",,," << t01-t00 << endl;

#if VCL
		// batch mode is being invoked by the v1.0 method from the VCL version of the model
		// set up colours for landscape in case any simulations require maps to be saved
		rgb colour;
		// COLOURS ARE COPIED FROM Formland.cpp, WHICH IS NOT THE MOST EFFICIENT METHOD...
		int RED[21]   = {  0,250,200,100,200,150,153,155,128,230,  0,  0,  0,  0,200,60,0,204,255,128,  0};
		int GREEN[21] = {200,200,200,250,150,150,128,100, 26,140,100,128,  0,180,200,60,0,179,255,102,  0};
		int BLUE[21]  = { 50,150,100,100,250,150,  0, 60,128,166,  0,115,255,190,200,60,0,  0,128,255,128};
		landParams paramsLand = pLandscape->getLandParams();
#if RSDEBUG
DEBUGLOG << "RunBatch(): land_nr = " << land_nr << " paramsLand.nHab = " << paramsLand.nHab
	<< endl;
#endif
		for (int i = 0; i < paramsLand.nHab; i++) {
			if (paramsLand.nHab < 21) {
				colour.r = RED[i]; colour.g = GREEN[i]; colour.b = BLUE[i];
			}
			else {
				colour.r = colour.g = colour.b = 0;
			}
			pLandscape->addColour(colour);
		}
		if (sim.saveMaps) {
			pLandscape->setLandMap();
			pLandscape->drawLandscape(0,0,0);
		}
#endif

		} // end of landOK condition

	} // end of imported landscape

	if (landOK) {

		// Open all other batch files and read header records
		if (ReadParameters(0,pLandscape)) {
			cout << endl << "Error opening ParameterFile - aborting batch run" << endl;
			return;
		}
#if RSDEBUG
//bool pppp = parameters.is_open();
//DEBUGLOG << "RunBatch(): parameterFile = " << parameterFile
//	<< " parameters.open() = " << pppp << endl;
#endif
		if (stagestruct) {
#if SEASONAL
			ReadStageStructure(0,pLandscape);
#else
			ReadStageStructure(0);
#endif  
		}
		ReadEmigration(0);
		ReadTransfer(0,pLandscape);
		ReadSettlement(0);
		if (geneticsFile != "NULL") ReadGenetics(0);
#if RS_CONTAIN
		if (name_managefile != "NULL") ReadManageFile(0,pLandscape);
#endif // RS_CONTAIN 
		ReadInitialisation(0,pLandscape);
#if VIRTUALECOLOGIST
		if (virtEcolFile != "NULL") ReadVirtEcol(0);
#endif

		// nSimuls is the total number of lines (simulations) in
		// the batch and is set in the control function
		string msgsim = "Simulation,";
		string msgerr = ",ERROR CODE,";
		string msgabt = ",simulation aborted";
		for (int i = 0; i < nSimuls; i++) {    
			t00 = (int)time(0);
			params_ok = true;
			read_error = ReadParameters(1,pLandscape);
			simParams sim = paramsSim->getSim();
			if (read_error) {
				rsLog << msgsim << sim.simulation << msgerr << read_error << msgabt << endl;
				params_ok = false;
			}
			if (stagestruct) {
#if SEASONAL
				ReadStageStructure(1,pLandscape);
#else
				ReadStageStructure(1);
#endif  
			}
			read_error = ReadEmigration(1);
			if (read_error) {
				rsLog << msgsim << sim.simulation << msgerr << read_error << msgabt << endl;
				params_ok = false;
			}
			read_error = ReadTransfer(1,pLandscape);
			if (read_error) {
				rsLog << msgsim << sim.simulation << msgerr << read_error << msgabt << endl;
				params_ok = false;
			}
			read_error = ReadSettlement(1);
			if (read_error) {
				rsLog << msgsim << sim.simulation << msgerr << read_error << msgabt << endl;
				params_ok = false;
			}
			if (params_ok) {
#if RSDEBUG
DebugGUI("RunBatch(): simulation i=" + Int2Str(i));
#endif
				pSpecies->setNChromosomes(0);
				pSpecies->setTraits();
			}
			if (geneticsFile == "NULL") {
				// use default genetics parameters
				// (by setting illegal values except for diploid)
				genomeData g{};
				g.nLoci = -1; g.probMutn = g.probCrossover = g.alleleSD = g.mutationSD = -1.0;
				if (reproductn == 0) g.diploid = false; else g.diploid = true;
				g.neutralMarkers = g.trait1Chromosome = false;

				pSpecies->setGenomeData(g);
			}
			else {
				read_error = ReadGenetics(1);
				if (read_error) {
					rsLog << msgsim << sim.simulation << msgerr << read_error << msgabt << endl;
					params_ok = false;
				}
			}
#if RS_CONTAIN
			if (name_managefile != "NULL") read_error = ReadManageFile(1,pLandscape);
			if (read_error) {
				rsLog << msgsim << sim.simulation << msgerr << read_error << msgabt << endl;
				params_ok = false;
			}
#endif // RS_CONTAIN 
			read_error = ReadInitialisation(1,pLandscape);
			if (read_error) {
				rsLog << msgsim << sim.simulation << msgerr << read_error << msgabt << endl;
				params_ok = false;
			}
#if VIRTUALECOLOGIST
			if (virtEcolFile == "NULL") {
				// no virtual ecologist
				paramsSim->setVirtEcol(false);
			}
			else {
				paramsSim->setVirtEcol(true);
				read_error = ReadVirtEcol(1);
				if (read_error) {
					rsLog << msgsim << sim.simulation << msgerr << read_error << msgabt << endl;
					params_ok = false;
				}
			}
#endif

			if (params_ok) {
				simParams sim = paramsSim->getSim();
#if RS_ABC

				ABCmain(pLandscape);

#else

#if RSDEBUG
DEBUGLOG << endl << "RunBatch(): i=" << i
	<< " simulation=" << sim.simulation << " landFile=" << landFile
	<< " outRange=" << sim.outRange << " outIntRange=" << sim.outIntRange
	<< endl;
#endif

				cout << endl << "Running simulation nr. " << Int2Str(sim.simulation)
					<< " on landscape no. " << Int2Str(land_nr) << endl;

				MemoLine(("Starting simulation " + Int2Str(sim.simulation) + "...").c_str());

				// for batch processing, include landscape number in parameter file name
				OutParameters(pLandscape);

				RunModel(pLandscape,i);
#if RSDEBUG
//DEBUGLOG << endl << "RunBatch(): real landscape, i = " << i
//	<< " simulation = " << sim.simulation << " landFile = " << landFile
//	<< endl;
#endif

#endif // RS_ABC

				t01 = (int)time(0);
				rsLog << msgsim << sim.simulation << "," << sim.reps
					<< "," << sim.years << "," << t01-t00 << endl;
			} // end of if (params_ok)
			else {
				cout << endl << "Error in reading parameter file(s)" << endl;
			}
#if VCL
			if (stopRun) break;
#endif
		} // end of nSimuls for loop

		// close input files
		ReadParameters(9,pLandscape);
#if SEASONAL
		if (stagestruct) ReadStageStructure(9,pLandscape);
#else
		if (stagestruct) ReadStageStructure(9);
#endif  
		ReadEmigration(9);
		ReadTransfer(9,pLandscape);
		ReadSettlement(9);
		if (geneticsFile != "NULL") ReadGenetics(9);
#if RS_CONTAIN
		if (name_managefile != "NULL") ReadManageFile(9,pLandscape);
#endif // RS_CONTAIN 
		ReadInitialisation(9,pLandscape);
#if VIRTUALECOLOGIST
		if (virtEcolFile != "NULL") ReadVirtEcol(9);
#endif

//		if (landtype != 9) 
		if (pLandscape != NULL) 
		{
			delete pLandscape; pLandscape = NULL;
		}

	} // end of landOK condition
#if VCL
	if (stopRun) break;
#endif

} // end of nLandscapes loop

ReadLandFile(9); // close the landFile

// Write performance data to log file
t1 = (int)time(0);
rsLog << endl << "Batch,,,," << t1-t0 << endl;

if (rsLog.is_open()) {
	rsLog.close(); rsLog.clear();
}

}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


#if RS_RCPP
void EOFerrorR(string filename)
{
	std::cout << "*** Did not read to EOF in " << filename << std::endl;
}
void StreamErrorR(string filename)
{
	std::cout << "*** Corrupted file stream in " << filename << std::endl << "You can try to use a different file encoding, like UTF-8." << std::endl;
#if RSDEBUG
	DEBUGLOG << "Corrupted file stream in " << filename << std::endl;
#endif
}
#endif // RS_RCPP
