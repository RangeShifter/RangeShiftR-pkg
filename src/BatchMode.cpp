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

#include "BatchMode.h"
//---------------------------------------------------------------------------

ifstream controlfile;
// Note - all batch files are prefixed 'b' here for reasons concerned with RS v1.0
ifstream bParamFile,bLandFile,bDynLandFile;
ifstream bSpDistFile,bStageStructFile,bTransMatrix;
ifstream bStageWeightsFile;
ifstream bEmigrationFile,bTransferFile,bSettlementFile;
ifstream bGeneticsFile,bArchFile,bInitFile,bInitIndsFile;

ofstream batchlog;

ofstream rsLog; // performance log for recording simulation times, etc.

// NOTE: THE STREAMS USED TO READ THE DATA AT RUN TIME COULD TAKE THE SAME NAMES AS
// USED DURING PARSING (ABOVE)
ifstream parameters;
ifstream ssfile,tmfile,fdfile,ddfile,sdfile;
ifstream emigFile,transFile,settFile,genFile,archFile,initFile,initIndsFile;
ifstream landfile,dynlandfile;

// global variables passed between parsing functions...
int batchnum;
int patchmodel,resolution,landtype,maxNhab,speciesdist,distresolution;
int reproductn;
int repseasons;
int stagestruct,stages,transfer;
int sexesDem;		// no. of explicit sexes for demographic model
int sexesDisp;	// no. of explicit sexes for dispersal model
int firstsimul = 0;
int fileNtraits; // no. of traits defined in genetic architecture file
//rasterdata landraster,patchraster,spdistraster;
rasterdata landraster;
// ...including names of the input files
string parameterFile;
string landFile;
string name_landscape,name_patch,name_dynland,name_sp_dist,name_costfile;
string stageStructFile,transMatrix;
string emigrationFile,transferFile,settleFile,geneticsFile,initialFile;
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
	if (reproductn < 0 || reproductn > 2) {
		BatchError(filetype,-999,2,"Reproduction"); errors++;
	}
	else {
		switch (reproductn) {
		case 0: { sexesDem = 1; sexesDisp = 1; break; }
		case 1: { sexesDem = 1; sexesDisp = 2; break; }
		case 2: { sexesDem = 2; sexesDisp = 2; break; }
		}
		b.reproductn = reproductn; b.sexesDem = sexesDem; b.sexesDisp = sexesDisp;
	}
}
else controlFormatError = true; // wrong control file format

controlfile >> paramname >> repseasons;
if (paramname == "RepSeasons") {
	if (repseasons < 1) {
		BatchError(filetype,-999,11,"RepSeasons"); errors++;
	}
	else b.repseasons = repseasons;
}
else controlFormatError = true; // wrong control file format

controlfile >> paramname >> stagestruct;
if (paramname == "StageStruct") {
	if (stagestruct < 0 || stagestruct > 1) {
		BatchError(filetype,-999,1,"StageStruct"); errors++;
	}
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
	if (transfer < 0 || transfer > 2) {
		BatchError(filetype,-999,2,"Transfer"); errors++;
	}
	else b.transfer = transfer;
}
else controlFormatError = true; // wrong control file format

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
		b.nSimuls = ParseParameterFile();
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

if (controlFormatError) {
	CtrlFormatError();
	b.ok = false;
}

if (controlfile.is_open()) 	{ controlfile.close(); controlfile.clear(); }
if (batchlog.is_open()) 		{ batchlog.close(); batchlog.clear(); }

// NOTE: THE FOLLOWING ELEMENTS COULD BE REMOVED FROM b ...
parameterFile = b.parameterFile;
landFile = b.landFile;
stageStructFile = b.stageStructFile;
emigrationFile = b.emigrationFile;
transferFile = b.transferFile;
settleFile = b.settleFile;
geneticsFile = b.geneticsFile;
initialFile = b.initFile;

return b;

}

//---------------------------------------------------------------------------
int ParseParameterFile(void)
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
bParamFile >> header; if (header != "ac" ) errors++;
bParamFile >> header; if (header != "std" ) errors++;
bParamFile >> header; if (header != "minR" ) errors++;
bParamFile >> header; if (header != "maxR" ) errors++;
bParamFile >> header; if (header != "minK" ) errors++;
bParamFile >> header; if (header != "maxK" ) errors++;
bParamFile >> header; if (header != "LocalExt" ) errors++;
bParamFile >> header; if (header != "LocalExtProb" ) errors++;
bParamFile >> header; if (header != "PropMales" ) errors++;
bParamFile >> header; if (header != "Harem" ) errors++;
bParamFile >> header; if (header != "bc" ) errors++;
bParamFile >> header; if (header != "Rmax" ) errors++;
for (i = 0; i < maxNhab; i++) {
	Kheader = "K" + Int2Str(i+1);
	bParamFile >> header; if (header != Kheader ) Kerrors++;
}
bParamFile >> header; if (header != "OutStartPop" ) errors++;
bParamFile >> header; if (header != "OutStartInd" ) errors++;
bParamFile >> header; if (header != "OutStartGenetic" ) errors++;
bParamFile >> header; if (header != "OutStartTraitCell" ) errors++;
bParamFile >> header; if (header != "OutStartTraitRow" ) errors++;
bParamFile >> header; if (header != "OutStartConn" ) errors++;
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
	bParamFile >> infloat;
	if (envstoch && (infloat < 0.0 || infloat >= 1.0)) { BatchError(filetype,line,20,"ac"); errors++; }
	bParamFile >> infloat;
	if (envstoch && (infloat <= 0.0 || infloat > 1.0)) { BatchError(filetype,line,20,"std"); errors++; }
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
//string fname,header,intext,ftype,costfile;
string fname,header,intext,ftype;
int j,inint,line;
float infloat;
rasterdata patchraster,spdistraster,costraster;
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

		// check cost map filename
		ftype = "CostMapFile";
		bLandFile >> name_costfile;
		if (name_costfile == "NULL") {
			if (transfer == 1) { // SMS
				if (landtype == 2) {
					BatchError(filetype,line,0," "); errors++;
					batchlog << ftype << " is required for a habitat quality landscape" << endl;
				}
			}
		}
		else {
			if (transfer == 1) { // SMS
				fname = indir + name_costfile;
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
				int something = ParseDynamicFile(indir,name_costfile);
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

	// check costs change filename
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
float infloat;
int errors = 0;
int simuls = 0;
int prevsimul;
bool checkfile;
vector <string> transfiles,wtsfiles;
string filetype = "StageStructFile";

// Parse header line;
bStageStructFile >> header; if (header != "Simulation" ) errors++;
bStageStructFile >> header; if (header != "PostDestructn" ) errors++;
bStageStructFile >> header; if (header != "PRep" ) errors++;
bStageStructFile >> header; if (header != "RepInterval" ) errors++;
bStageStructFile >> header; if (header != "MaxAge" ) errors++;
bStageStructFile >> header; if (header != "TransMatrixFile" ) errors++;
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
	bStageStructFile >> inint;
	if (inint < 0 || inint > 1) { BatchError(filetype,line,1,"PostDestructn"); errors++; }
	bStageStructFile >> infloat;
	if (infloat <= 0 || infloat > 1.0) { BatchError(filetype,line,20,"PRep"); errors++; }
	bStageStructFile >> inint;
	if (inint < 0) { BatchError(filetype,line,19,"RepInterval"); errors++; }
	bStageStructFile >> inint;
	if (inint < 2) { BatchError(filetype,line,12,"MaxAge"); errors++; }

	bStageStructFile >> filename;
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

	bStageStructFile >> inint;
	if (inint < 0 || inint > 2) { BatchError(filetype,line,2,"SurvSched"); errors++; }
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
if (totfecundity <= 0.0) {
	BatchError(filetype,line,10,"Total fecundity"); errors++;
}
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

//---------------------------------------------------------------------------
int ParseEmigFile(void)
{
string header;
int simul;
int densdep,usefullkern,stagedep,sexdep,indvar,emigstage,stage,sex;
bool densdepset,indvarset;
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
simCheck current, prev;
simul = -98765;
prev.simul = -999;
prev.simlines = prev.reqdsimlines = 0;
bEmigrationFile >> simul;
// first simulation number must match first one in parameterFile
if (simul != firstsimul) {
	BatchError(filetype,line,111,"Simulation"); errors++;
}
current.simul = 0; //dummy line to prevent warning message in VisualStudio 2019
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
float mortProb,slope,inflPoint;
float morthab,mortmatrix;
int costhab,costmatrix;
float SL,rho;
float StepLMean,StepLSD,RhoMean,RhoSD,StepLScale,RhoScale;

vector <string> costsfiles;

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
	bTransferFile >> header; if (header != "StraightenPath" ) errors++;
	bTransferFile >> header; if (header != "SMtype" ) errors++;
	bTransferFile >> header; if (header != "SMconst" ) errors++;
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
simCheck current, prev;
simul = -98765;
prev.simul = -999;
prev.simlines = prev.reqdsimlines = 0;
bTransferFile >> simul;
// first simulation number must match first one in parameterFile
if (simul != firstsimul) {
	BatchError(filetype,line,111,"Simulation"); errors++;
}
current.simul = 0; //dummy line to prevent warning message in VisualStudio 2019
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
		bTransferFile >> straightenPath	>> smtype >> smconst;
		if (straightenPath < 0 || straightenPath > 1) {
			BatchError(filetype,line,1,"StraightenPath"); errors++;
		}
		if (landtype == 2) // habitat quality landscape 
		{ // must have constant mortality
			if (smtype != 0) {
				BatchError(filetype,line,0," "); errors++;
				batchlog << "SMtype must be 0 for LandType 2" << endl;
			}
		}
		else {
			if (smtype < 0 || smtype > 1) {
				BatchError(filetype,line,1,"SMtype"); errors++;
			}
		}
		if (smtype == 0) 
		{
			if (smconst < 0.0 || smconst >= 1.0) {
				BatchError(filetype,line,20,"SMconst"); errors++;
			}
		}
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


//---------------------------------------------------------------------------
int ParseSettleFile(void)
{
string header;
int simul,stagedep,sexdep,stage,sex,settletype;
int densdep,indvar,findmate,minSteps,maxSteps,maxStepsYear;
float s0,alphaS,betaS;
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
if (transfer == 0) 
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
simCheck current, prev;
simul = -98765;
prev.simul = -999;
prev.simlines = prev.reqdsimlines = 0;
bSettlementFile >> simul;
// first simulation number must match first one in parameterFile
if (simul != firstsimul) {
	BatchError(filetype,line,111,"Simulation"); errors++;
}
current.simul = 0; //dummy line to prevent warning message in VisualStudio 2019
while (simul != -98765) {
	if (transfer == 0) 
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
simCheck current,prev;
simul = -98765;
prev.simul = -999;
prev.simlines = prev.reqdsimlines = 0;
bGeneticsFile >> simul;
// first simulation number must match first one in parameterFile
if (simul != firstsimul) {
	BatchError(filetype,line,111,"Simulation"); errors++;
}
current.simul = 0; //dummy line to prevent warning message in VisualStudio 2019
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
simCheck current, prev;
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
current.simul = 0; //dummy line to prevent warning message in VisualStudio 2019
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
simCheck current;
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
	if (ppGenLand.continuous) ppLand.nHab = 2;		
	else ppLand.nHab = 1;
#if RSDEBUG
DEBUGLOG << "ReadLandFile(): ppLand.landNum=" << ppLand.landNum
	<< " continuous=" << ppGenLand.continuous << " ppLand.nHab=" << ppLand.nHab
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
	if (landtype == 2) ppLand.nHab = 1; // habitat quality landscape has one habitat class
#if RSDEBUG
DEBUGLOG << "ReadLandFile(): ppLand.landNum=" << ppLand.landNum
	<< " name_landscape=" << name_landscape
	<< " name_patch=" << name_patch
	<< " name_costfile=" << name_costfile
	<< " name_dynland=" << name_dynland
	<< " name_sp_dist=" << name_sp_dist
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
// as from v1.1, there is just one pair of min & max values,
// which are attributes of the species
// ULTIMATELY, THE PARAMETER FILE SHOULD HAVE ONLY TWO COLUMNS ...
//parameters >> env.ac >> env.std >> env.minR >> env.maxR >> env.minK >> env.maxK;
float minR,maxR,minK,maxK; 
parameters >> env.ac >> env.std >> minR >> maxR >> minK >> maxK;
if (env.inK) {
	float minKK,maxKK;
	minKK = minK * (((float)paramsLand.resol * (float)paramsLand.resol) / 10000.0f);
	maxKK = maxK * (((float)paramsLand.resol * (float)paramsLand.resol) / 10000.0f);
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

parameters >> dem.propMales >> dem.harem >> dem.bc >> dem.lambda;
pSpecies->setDemogr(dem);

float k;

if (landtype == 9) { // artificial landscape
	// only one value of K is read, but it must be applied as the second habitat if the
	// landscape is discrete (the first is the matrix where K = 0) or as the first 
	// (only) habitat if the landscape is continuous
	genLandParams genland = pLandscape->getGenLandParams(); 
	int nhab;
	if (genland.continuous) nhab = 1;		
	else nhab = 2;  
	pSpecies->createHabK(nhab);
	parameters >> k;
	k *= (((float)paramsLand.resol*(float)paramsLand.resol))/10000.0f;
	if (genland.continuous) {
		pSpecies->setHabK(0,k);		
	}
	else {
		pSpecies->setHabK(0,0);
		pSpecies->setHabK(1,k);
	}
}
else {
	pSpecies->createHabK(paramsLand.nHabMax);
	for (int i = 0; i < paramsLand.nHabMax; i++) {
		parameters >> k;
		k *= ((float)paramsLand.resol*(float)paramsLand.resol)/10000.0f;
		pSpecies->setHabK(i,k);
	}
}

#if RSDEBUG
DEBUGLOG << "ReadParameters(): dem.lambda=" << dem.lambda
	<< " habK[0] = " << pSpecies->getHabK(0)
	<< " nHabMax = " << paramsLand.nHabMax
	<< endl;
#endif

parameters >> sim.outStartPop >> sim.outStartInd >> sim.outStartGenetic
	>> sim.outStartTraitCell >> sim.outStartTraitRow >> sim.outStartConn;
parameters >> sim.outIntRange >> sim.outIntOcc >> sim.outIntPop >> sim.outIntInd
	>> sim.outIntGenetic >> sim.outGenType >> jjjj
	>> sim.outIntTraitCell >> sim.outIntTraitRow >> sim.outIntConn;
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
int ReadStageStructure(int option)
{
string name;
int simulation,postDestructn;
stageParams sstruct = pSpecies->getStage();
string Inputs = paramsSim->getDir(1);

if (option == 0) { // open file and read header line
	ssfile.open(stageStructFile.c_str());
	string header;
	int nheaders = 18;
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
ssfile >> postDestructn >> sstruct.probRep >> sstruct.repInterval >> sstruct.maxAge;
if (postDestructn == 1) sstruct.disperseOnLoss = true;
else sstruct.disperseOnLoss = false;

ssfile >> name; 
// 'name' is TransMatrixFile
tmfile.open((Inputs+name).c_str());
ReadTransitionMatrix(sstruct.nStages,sexesDem,0,0);
tmfile.close(); tmfile.clear();
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

pSpecies->setStage(sstruct);

if (sstruct.devDens || sstruct.survDens) {
	pSpecies->setDensDep(devCoeff,survCoeff);
}

return 0;
}

//---------------------------------------------------------------------------
int ReadTransitionMatrix(short nstages,short nsexesDem,short hab,short season)  
{
//#if RS_CONTAIN
//int hab = 0; // TEMPORARY set suitable habitat to 0
//#endif // RS_CONTAIN 
int ii;
int minAge;
float ss,dd; 
string header;
demogrParams dem = pSpecies->getDemogr();
//stageParams sstruct = pSpecies->getStage();

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
		pSpecies->setFec(j,0,matrix[j][0]);
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
		pSpecies->setSurv(j,0,ss+dd);
		if ((ss+dd) > 0.0f)
			pSpecies->setDev(j,0,dd/(ss+dd));
		else
			pSpecies->setDev(j,0,0.0);
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
		if (j%2 == 0)
			pSpecies->setFec(ii,1,matrix[j][0]);
		else {
			pSpecies->setFec(ii,0,matrix[j][0]);
			ii++;
		}
	}
	// survival and development of male juveniles
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
			pSpecies->setSurv(ii,1,(ss+dd));
			if ((ss+dd) > 0.0)
				pSpecies->setDev(ii,1,dd/(ss+dd));
			else
				pSpecies->setDev(ii,1,0.0);
		}
		else{ // females
//			for (int i = 0; i < sstruct.nStages*2; i++) 
			for (int i = 0; i < nstages*2; i++) 
			{
				if (j == i+1) ss = matrix[j][i];
				if (j == i-1) dd = matrix[j][i];
			}
			pSpecies->setSurv(ii,0,(ss+dd));
			if ((ss+dd) > 0.0)
				pSpecies->setDev(ii,0,dd/(ss+dd));
			else
				pSpecies->setDev(ii,0,0.0);
			ii++;
		}
	}
//	for (int j = 0; j < sstruct.nStages*2; j++) delete[] matrix[j];
//	delete[] matrix;
}

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


//---------------------------------------------------------------------------
int ReadEmigration(int option)
{
int error = 0;

if (option == 0) { // open file and read header line
	emigFile.open(emigrationFile.c_str());
	string header;
	for (int i = 0; i < 25; i++) emigFile >> header;
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
float	ep,d0,alpha,beta,epMean,epSD,d0Mean,d0SD,alphaMean,alphaSD,betaMean,betaSD;
float epScale,d0Scale,alphaScale,betaScale;
bool firstline = true;
demogrParams dem = pSpecies->getDemogr();
stageParams sstruct = pSpecies->getStage();
emigRules emig = pSpecies->getEmig();
emigTraits etraits;
emigParams eparams;

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
	<< " trfr.twinKern=" << trfr.twinKern
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
		nheaders = 23;
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
if (trfr.moveModel) TransferType = trfr.moveType; else TransferType = 0;

int sexKernels = 0;
trfrKernTraits k;
trfrMovtTraits move;
trfrKernParams kparams;
trfrScales scale = pSpecies->getTrfrScales();
string CostsFile;
trfrSMSParams smsparams;

switch (TransferType) {

case 0: // negative exponential dispersal kernel

	firstline = true;

	// set no.of lines assuming maximum stage- and sex-dependency
	if (sstruct.nStages == 0) Nlines = sexesDisp;
	else Nlines = sstruct.nStages * sexesDisp;

	for (int line = 0; line < Nlines; line++) {

		transFile >> simulation >> stageDep >> sexDep >> iiii >> jjjj >> kkkk;
		if (firstline) {
			firstsimul = simulation;
			if (iiii == 0) trfr.twinKern = false; else trfr.twinKern = true;
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
			pSpecies->setKernTraits(0,0,k,paramsLand.resol);
			if (trfr.indVar) {
				transFile >> kparams.dist1Mean >> kparams.dist1SD
					>> kparams.dist2Mean >> kparams.dist2SD
					>> kparams.PKern1Mean >> kparams.PKern1SD;
//				MAXDist = kparams.maxDist1;
				pSpecies->setKernParams(0,0,kparams,(double)paramsLand.resol);
			}
			else {
				for (int i = 0; i < 6; i++) transFile >> tttt;
			}
			break;

		case 1: // sex-dependent
			if (trfr.indVar)
			{
				if (trfr.twinKern) 
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
				if (trfr.twinKern) 
				{
					transFile >> k.meanDist1 >> k.meanDist2 >> k.probKern1;
					for (int i = 0; i < 6; i++) transFile >> tttt;
				}
				else {
					transFile >> k.meanDist1; k.meanDist2 = k.meanDist1; k.probKern1 = 1.0;
					for (int i = 0; i < 8; i++) transFile >> tttt;
				}
			pSpecies->setKernTraits(0,sex,k,paramsLand.resol);
			}
			break;

		case 2: // stage-dependent
			if (trfr.twinKern) 
			{
				transFile >> k.meanDist1 >> k.meanDist2 >> k.probKern1;
				for (int i = 0; i < 6; i++) transFile >> tttt;
			}
			else {
				transFile >> k.meanDist1; k.meanDist2 = k.meanDist1; k.probKern1 = 1.0;
				for (int i = 0; i < 8; i++) transFile >> tttt;
			pSpecies->setKernTraits(stage,0,k,paramsLand.resol);
			}
			break;

		case 3: // sex- & stage-dependent
			if (trfr.twinKern) 
			{
				transFile >> k.meanDist1 >> k.meanDist2 >> k.probKern1;
				for (int i = 0; i < 6; i++) transFile >> tttt;
			}
			else {
				transFile >> k.meanDist1; k.meanDist2 = k.meanDist1; k.probKern1 = 1.0;
				for (int i = 0; i < 8; i++) transFile >> tttt;
			pSpecies->setKernTraits(stage,sex,k,paramsLand.resol);
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
			trfrMortParams mort;
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
	transFile >> jjjj >> iiii >> move.stepMort;
	if (iiii == 0) trfr.habMort = false; else trfr.habMort = true;
	if (jjjj == 0) move.straigtenPath = false; else move.straigtenPath = true;

#if RSDEBUG
DEBUGLOG << "ReadTransfer(): SMtype=" << trfr.habMort << " SMconst=" << move.stepMort << endl;
#endif // RSDEBUG

	if (!paramsLand.generated) { // real landscape
		if (paramsLand.rasterType == 0) { // habitat codes
			if (trfr.habMort)
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
		if (trfr.habMort) 
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

	{
	trfrCRWParams mparams;

	transFile >> simulation >> iiii;
	if (iiii == 0) trfr.indVar = false; else trfr.indVar = true;

	transFile >> move.stepLength >> move.rho
		>> mparams.stepLgthMean >> mparams.stepLgthSD >> mparams.rhoMean >> mparams.rhoSD;
	transFile >> scale.stepLScale >> scale.rhoScale;
	transFile >> jjjj >> iiii >> move.stepMort;
	if (iiii == 0) trfr.habMort = false; else trfr.habMort = true;
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
	if (trfr.habMort && paramsLand.rasterType != 0) error = 434;

	if (!paramsLand.generated && paramsLand.rasterType == 0) { // real habitat codes landscape
		if (trfr.habMort)
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


default:
	error = 440;
	break;
} // end of switch (TransferType)

return error;
}

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
settleSteps ssteps;
settleTraits settleDD;
settParams sparams;
int sexSettle = 0,settType = 0,densdep,indvar,findmate;

#if RSDEBUG
//DEBUGLOG << "ReadSettlement(): option=" << option << " transfer=" << transfer 
//	<< " trfr.moveModel=" << trfr.moveModel << endl;
#endif

if (option == 0) { // open file and read header line
	settFile.open(settleFile.c_str());
	string header;
	int nheaders = 0;
	if (trfr.moveModel) nheaders = 23;
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
	if (transfer == 0) 
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
		if (transfer == 0) sett.indVar = false;  // dispersal kernel
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
genomeData g;
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
initInd iind;
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

//---------------------------------------------------------------------------
void RunBatch(int nSimuls, int nLandscapes)
{
int land_nr;
int t0,t1,t00,t01;
int read_error;
bool params_ok;
simParams sim = paramsSim->getSim();

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
rsLog << "RANDOM SEED," << RS_random_seed << ",,," << endl;

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
		if (paramsLand.patchModel) {
			string pname = paramsSim->getDir(1) + name_patch;
#if RSDEBUG
int t02a = time(0);
#endif
			landcode = pLandscape->readLandscape(0,hname,pname,cname);
#if RSDEBUG
int t02b = time(0);
DEBUGLOG << "RunBatch(): TIME for readLandscape() " << t02b-t02a << endl;
#endif
		}
		else {
			landcode = pLandscape->readLandscape(0,hname," ",cname);
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

		if (landOK) {
			t01 = (int)time(0);
			rsLog << "Landscape," << land_nr << ",,," << t01-t00 << endl;

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
			ReadStageStructure(0);
		}
		ReadEmigration(0);
		ReadTransfer(0,pLandscape);
		ReadSettlement(0);
		if (geneticsFile != "NULL") ReadGenetics(0);
		ReadInitialisation(0,pLandscape);

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
				ReadStageStructure(1);
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
				genomeData g;
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
			read_error = ReadInitialisation(1,pLandscape);
			if (read_error) {
				rsLog << msgsim << sim.simulation << msgerr << read_error << msgabt << endl;
				params_ok = false;
			}

			if (params_ok) {
				simParams sim = paramsSim->getSim();

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

				t01 = (int)time(0);
				rsLog << msgsim << sim.simulation << "," << sim.reps
					<< "," << sim.years << "," << t01-t00 << endl;
			} // end of if (params_ok)
			else {
				cout << endl << "Error in reading parameter file(s)" << endl;
			}
		} // end of nSimuls for loop

		// close input files
		ReadParameters(9,pLandscape);
		if (stagestruct) ReadStageStructure(9);
		ReadEmigration(9);
		ReadTransfer(9,pLandscape);
		ReadSettlement(9);
		if (geneticsFile != "NULL") ReadGenetics(9);
		ReadInitialisation(9,pLandscape);

//		if (landtype != 9) 
		if (pLandscape != NULL) 
		{
			delete pLandscape; pLandscape = NULL;
		}

	} // end of landOK condition

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


