/*----------------------------------------------------------------------------
 *
 *	Copyright (C) 2020 Anne-Kathleen Malchow, Greta Bocedi, Stephen C.F. Palmer, Justin M.J. Travis, Damaris Zurell
 *
 *	This file is part of RangeShiftR.
 *
 *	RangeShiftR is free software: you can redistribute it and/or modify
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
 *	along with RangeShiftR. If not, see <https://www.gnu.org/licenses/>.
 *
 --------------------------------------------------------------------------*/


/*------------------------------------------------------------------------------

 RangeShifter v2.0 Main

 Entry level function for the R-package RangeshiftR.

 For compilation with GCC g++

 For full details of RangeShifter, please see:
 Bocedi G., Palmer S.C.F., Pe'er G., Heikkinen R.K., Matsinos Y.G., Watts K.
 and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
 eco-evolutionary dynamics and species' responses to environmental changes.
 Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

 Author: Anne-Kathleen Malchow, Humboldt University Berlin
 Jette Reeg, University of Potsdam
 large parts modified from 'Main.cpp' and 'BatchMode.cpp' created by
 Steve Palmer, Roslyn Henry and Théo Pannetier, University of Aberdeen

 ------------------------------------------------------------------------------*/

#include "Rinterface.h"

class msghdrs1;

string habmapname, patchmapname, distnmapname; // req'd for compilation, but not used
string costmapname, genfilename;               // ditto
vector<string> hfnames;                        // ditto

paramGrad* paramsGrad;   // pointer to environmental gradient parameters
paramStoch* paramsStoch; // pointer to environmental stochasticity parameters
paramInit* paramsInit;   // pointer to initialisation parameters
paramSim* paramsSim;     // pointer to simulation parameters

Species* pSpecies; // pointer to species
Community* pComm;  // pointer to community
RSrandom* pRandom; // pointer to random number routines
Management* pManagement; // pointer to management routines

#if RSDEBUG
ofstream DEBUGLOG;
ofstream MUTNLOG;
#endif
// ofstream batchlog;
ofstream rsLog; // performance log for recording simulation times, etc.


// global variables passed between parsing functions...
int batchnum;
int patchmodel, resolution, landtype, maxNhab, speciesdist, distresolution;
int reproductn;
int repseasons;
int stagestruct, stages, gTransferType;
int sexesDem;  // no. of explicit sexes for demographic model
int gNbSexesDisp; // no. of explicit sexes for dispersal model
int gFirstSimulNb;
int fileNtraits; // no. of traits defined in genetic architecture file
bool gHasGenetics; // instead of static bool anyIndVar;
bool gHasNeutralGenetics; // instead of static bool anyIndVar;
bool gHasGeneticLoad; // instead of static bool anyIndVar;
DispersalTraitInputOptions gDispTraitOpt;
//vector<int> gNbTraitFileRows;

int translocation;
rasterdata landraster,patchraster,spdistraster,costsraster;

string name_landscape, name_patch, name_sp_dist, name_costfile;

string msgnlines = "No. of lines for final Simulation ";
string msgshldbe = " should be ";
string msgresol0 = "*** Resolution of ";
string msgresol1 = " does not match set Resolution ";
string msgresol2 = " does not match set Distribution Resolution ";
string msghdrs0 = "*** Headers of ";
string msghdrs1 = " do not match headers of LandscapeFile";
string msgpatch = " is required for patch-based model";
string msgmatch = " must match the specification exactly";
string msgcase = " case-sensitive parameter names";



//---------------------------------------------------------------------------
// Returns input value less next highest power of 2 (for x > 2)
int power2check(int x) {
    if (x < 2) return 0;
    int r = x % 2;
    while (r == 0) {
        x /= 2; r = x % 2;
    }
    return x;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//------ R interface --------------------------------------------------------
//---------------------------------------------------------------------------


//[[Rcpp::export(name = "run_from_R")]]
Rcpp::List BatchMainR(std::string dirpath, Rcpp::S4 ParMaster)
{
    // int i,t0,t1,Nruns;
    int t0, t1;
    int nSimuls, nLandscapes; // no. of simulations and landscapes in batch

    t0 = time(0);
    //clear_outPop();

    // set up parameter objects
    paramsGrad = new paramGrad;
    paramsStoch = new paramStoch;
    paramsInit = new paramInit;
    paramsSim = new paramSim(dirpath);

    // set up working directory and log file names
    // paramsSim->setDir(dirpath);

#if RSDEBUG
    Rcpp::Rcout << endl << "Working directory: " << paramsSim->getDir(0) << endl;
#endif

    bool errorfolder = CheckDirectory(dirpath);
    // if(errorfolder) {
    //     Rcpp::Rcout << endl << "***** Invalid working directory: " << paramsSim->getDir(0) << endl << endl;
    //     Rcpp::Rcout << "***** Working directory must contain Inputs, Outputs and Output_Maps folders" << endl << endl;
    //     Rcpp::Rcout << "*****" << endl;
    //     Rcpp::Rcout << "***** Simulation ABORTED " << endl;
    //     Rcpp::Rcout << "*****" << endl;
    //     return Rcpp::List::create(Rcpp::Named("Errors") = 666);
    // }

#if RSDEBUG
    // set up debugging log file
    string logname = paramsSim->getDir(2) + "DebugLog.txt";
    DEBUGLOG.open(logname.c_str());
    logname = paramsSim->getDir(2) + "MutnLog.txt";
    MUTNLOG.open(logname.c_str());
    if(DEBUGLOG.is_open())
        Rcpp::Rcout << endl << "Main(): DEBUGLOG is open" << endl << endl;
    else
        Rcpp::Rcout << endl << "Main(): DEBUGLOG is NOT open" << endl << endl;
#else
#endif


    // set global variables
    Rcpp::S4 control("ControlParams");
    // ParMaster = global.get(ParMaster_name);
    control = Rcpp::as<Rcpp::S4>(ParMaster.slot("control"));
    setglobalvarsR(control);

    nSimuls = 1;
    nLandscapes = 1;

    //////////////////////// batchfiles ParseControlR(string dirpath) ///////////////////

    //string indir  = dirpath + "Inputs/";
    string outdir = dirpath + "Outputs/";

    // Doublecheck global variables
    Rcpp::Rcout << "Checking Control parameters " << endl;

    int errors = 0;
    string filetype = "Control parameters";

    if(batchnum < 0) {
        BatchErrorR(filetype, -999, 19, "BatchNum");
        errors++;
    }

    if(patchmodel < 0 || patchmodel > 1) {
        BatchErrorR(filetype, -999, 1, "PatchModel");
        errors++;
    }

    if(resolution < 1) {
        BatchErrorR(filetype, -999, 11, "Resolution");
        errors++;
    }

    if(landtype != 0 && landtype != 2 && landtype != 9) {
        BatchErrorR(filetype, -999, 0, "LandType");
        Rcpp::Rcout << "LandType must be 0, 2 or 9" << endl;
        errors++;
    } else {
        if(landtype == 9 && patchmodel) {
            BatchErrorR(filetype, -999, 0, "LandType");
            Rcpp::Rcout << "LandType may not be 9 for a patch-based model" << endl;
            errors++;
        }
    }

    if(landtype == 0) { // raster with unique habitat codes
        if(maxNhab < 2) {
            BatchErrorR(filetype, -999, 12, "MaxHabitats");
            errors++;
        }
    } else { // raster with habitat quality OR artificial landscape
        if(maxNhab != 1) {
            BatchErrorR(filetype, -999, 0, " ");
            errors++;
            Rcpp::Rcout << "MaxHabitats must be 1 for LandType = " << landtype << endl;
        }
        /*else
         {
         if (landtype == 9) // artificial landscape
         // although the user enters 1, the actual number of habitats is 2
         b.maxNhab = 2; // only changes maxNhab in batchfile b, the global variable remains maxNhab=1
         }*/
    }

    if(speciesdist < 0 || speciesdist > 1) {
        BatchErrorR(filetype, -999, 1, "SpeciesDist");
        errors++;
    } else {
        if(speciesdist != 0 && landtype == 9) {
            BatchErrorR(filetype, -999, 0, "SpeciesDist");
            Rcpp::Rcout << "SpeciesDist must be 0 for an artificial landscape" << endl;
            errors++;
        }
    }

    if(speciesdist == 1) { // distribution resolution is required
        if(distresolution < resolution) {
            BatchErrorR(filetype, -999, 0, "DistResolution");
            Rcpp::Rcout << "DistResolution may not be less than Resolution" << endl;
            errors++;
        } else {
            if(distresolution % resolution) {
                BatchErrorR(filetype, -999, 0, "DistResolution");
                Rcpp::Rcout << "DistResolution must be an integer multiple of Resolution" << endl;
                errors++;
            }
        }
    }

    if(reproductn < 0 || reproductn > 2) {
        BatchErrorR(filetype, -999, 2, "Reproduction");
        errors++;
    }
    else {
        switch(reproductn) {
        case 0: {
        sexesDem = 1;
        gNbSexesDisp = 1;
        break;
    }
        case 1: {
            sexesDem = 1;
            gNbSexesDisp = 2;
            break;
        }
        case 2: {
            sexesDem = 2;
            gNbSexesDisp = 2;
            break;
        }
        }
    }

    if(repseasons < 1) {
        BatchErrorR(filetype, -999, 11, "RepSeasons");
        errors++;
    }

    if(stagestruct < 0 || stagestruct > 1) {
        BatchErrorR(filetype, -999, 1, "StageStruct");
        errors++;
    }

    if(stagestruct) {
        if(stages < 2 || stages > gMaxNbStages) {
            BatchErrorR(filetype, -999, 0, " ");
            errors++;
            Rcpp::Rcout << "Stages must be between 2 and " << gMaxNbStages << endl;
        }
    } else { // non-stage-structured model must have 2 stages
        stages = 2;
    }

    if( gTransferType < 0 ||  gTransferType > 2) {
        BatchErrorR(filetype, -999, 2, "Transfer");
        errors++;
    }

    if(errors > 0) { // terminate batch error checking

        Rcpp::Rcout << endl << "*** Control parameters must be corrected before proceeding" << endl;
        return Rcpp::List::create(Rcpp::Named("Errors") = errors);
    }

    //////////////////////////////////////


    // set up species
    // FOR MULTI-SPECIES MODEL, THERE WILL BE AN ARRAY OF SPECIES POINTERS
    // OR A COMMUNITY CLASS TO HOLD THE SPECIES
    pSpecies = new Species;
    demogrParams dem = pSpecies->getDemogrParams();
    stageParams sstruct = pSpecies->getStageParams();
    transferRules trfr = pSpecies->getTransferRules();

    // create new Management
    pManagement = new Management;
    managementParams m = pManagement->getManagementParams();

    if(errors == 0) {
        //   nSimuls = b.nSimuls;
        //   nLandscapes = b.nLandscapes;
        dem.repType = reproductn;
        dem.repSeasons = repseasons;
        if(stagestruct == 0)
            dem.stageStruct = false;
        else
            dem.stageStruct = true;
        sstruct.nStages = stages;
        if( gTransferType == 0)
            trfr.usesMovtProc = false;
        else {
            trfr.usesMovtProc = true;
            trfr.moveType =  gTransferType;
        }
        if(translocation == 0)
            m.translocation = false;
        else
            m.translocation = true;

        pSpecies->setDemogr(dem);
        pSpecies->setStage(sstruct);
        pSpecies->setTrfrRules(trfr);
        pManagement->setManagementParams(m);

        simParams sim = paramsSim->getSim();
        sim.batchMode = true;
        sim.batchNum = batchnum;
        paramsSim->setSim(sim);
        Rcpp::Rcout << endl << "Control Parameters checked" << endl;
    } else {
        Rcpp::Rcout << endl << "Error in parsing parameters - " << endl;
    }

#if RSDEBUG
    DEBUGLOG << "BatchMainR(): dem.repType = " << dem.repType << endl;
#endif

    // set up random number stream
    std::int64_t seed = Rcpp::as<std::int64_t>(control.slot("seed"));

    if(seed == 0) { // don't set seed from R
#if RSDEBUG
        pRandom = new RSrandom(666);  // fixed debug seed
#else
        pRandom = new RSrandom(-1);  // random seed
#endif
    }
    else pRandom = new RSrandom(seed);

    //Rcpp::RNGScope rngScope;

    Rcpp::List list_outPop;
    if(errors == 0) {
        Rcpp::Rcout << endl << "Run Simulation(s)";
        if(seed > 0) {
            Rcpp::Rcout << " with seed " << seed;
        }else{
            Rcpp::Rcout << " with random seed";
        }
        Rcpp::Rcout << " ..." << endl;

        list_outPop = RunBatchR(nSimuls, nLandscapes, ParMaster);
    }

#if RSDEBUG
    if(DEBUGLOG.is_open()) {
        DEBUGLOG.close();
        DEBUGLOG.clear();
    }
    if(MUTNLOG.is_open()) {
        MUTNLOG.close();
        MUTNLOG.clear();
    }
#endif

    simParams sim = paramsSim->getSim();

    delete paramsGrad;
    delete paramsStoch;
    delete paramsInit;
    delete paramsSim;
    delete pSpecies;
    delete pManagement;
    delete pRandom;

    t1 = time(0);
    Rcpp::Rcout << endl << "***** Elapsed time: " << t1 - t0 << " seconds" << endl << endl;

    Rcpp::Rcout << "*****" << endl;
    Rcpp::Rcout << "***** Simulation completed " << endl;
    Rcpp::Rcout << "***** Outputs folder: " << outdir << endl;
    Rcpp::Rcout << "*****" << endl;

    if(sim.ReturnPopRaster && sim.outIntPop > 0) {
        // return Rcpp::List::create(Rcpp::Named("runs") = errors);
        return list_outPop;
    } else {
        return Rcpp::List::create(Rcpp::Named("Errors") = errors);
    }
}

//---------------------------------------------------------------------------
// TODO: Change return value to integer code: error codes as negative values!
bool ReadLandParamsR(Landscape* pLandscape, Rcpp::S4 ParMaster)
{
    landParams ppLand = pLandscape->getLandParams();
    genLandParams ppGenLand = pLandscape->getGenLandParams();

    Rcpp::S4 LandParamsR("LandParams");
    LandParamsR = Rcpp::as<Rcpp::S4>(ParMaster.slot("land"));

    int errors = 0;

    if(landtype == 9) { // artificial landscape
#if RSDEBUG
        DEBUGLOG << "ReadLandParamsR(): artificial: " << endl;
#endif
        ppLand.nHab = 2;
        ppLand.rasterType = 9;
        ppLand.landNum = Rcpp::as<int>(LandParamsR.slot("LandNum"));
        ppGenLand.fractal = Rcpp::as<bool>(LandParamsR.slot("fractal"));
        ppGenLand.continuous = Rcpp::as<bool>(LandParamsR.slot("continuous"));
        ppLand.dimX = Rcpp::as<int>(LandParamsR.slot("dimX"));
        ppLand.dimY = Rcpp::as<int>(LandParamsR.slot("dimY"));
        ppGenLand.minPct = Rcpp::as<float>(LandParamsR.slot("minPct"));
        ppGenLand.maxPct = Rcpp::as<float>(LandParamsR.slot("maxPct"));
        ppGenLand.propSuit = Rcpp::as<float>(LandParamsR.slot("propSuit"));
        ppGenLand.hurst = Rcpp::as<float>(LandParamsR.slot("hurst"));
        ppLand.maxX = ppLand.dimX - 1;
        ppLand.maxY = ppLand.dimY - 1;

        // check dimensions  -->  these shouldn't occur, as already checked for on R-level
        if(ppGenLand.fractal && ppLand.maxX > ppLand.maxY) {
            //return -901;
            // fix it by swapping X and Y:
            ppLand.maxX = ppLand.dimY - 1;
            ppLand.maxY = ppLand.dimX - 1;
            ppLand.dimX = ppLand.maxX + 1;
            ppLand.dimY = ppLand.maxY + 1;
        }
        if(ppGenLand.fractal) {
            if((ppLand.dimX < 3 || ppLand.dimX % 2 != 1) // why only check for uneven number, not ((power of 2)-1) ?
                   || (ppLand.dimY < 3 || ppLand.dimY % 2 != 1)) {
                //return -902;  // -> no action, should be covered on R-level
            }
        }
        // SCFP 26/9/13 - min and max habitat percentages need to be set for all types of
        // fractal landscape (including discrete), as they are passed to the fractal generator
        // NOTE that will not have been checked for a discrete landscape
        if(ppGenLand.fractal && !ppGenLand.continuous) {
            ppGenLand.minPct = 1;
            ppGenLand.maxPct = 100;
        }
    } else { // imported raster map

        ppLand.landNum = Rcpp::as<int>(LandParamsR.slot("LandNum"));
        ppLand.nHab = Rcpp::as<int>(LandParamsR.slot("Nhabitats")); // no longer necessary to read no. of habitats from landFile

        Rcpp::IntegerVector dynland_years;
        Rcpp::StringVector habitatmaps, patchmaps, costmaps;
        dynland_years = Rcpp::as<Rcpp::IntegerVector>(LandParamsR.slot("DynamicLandYears"));
        if(dynland_years.size() == 1 && dynland_years[0] == 0 ) ppLand.dynamic = false;
        else ppLand.dynamic = true;
        if(ppLand.dynamic) {
            habitatmaps = Rcpp::as<Rcpp::StringVector>(LandParamsR.slot("LandscapeFile"));
            patchmaps = Rcpp::as<Rcpp::StringVector>(LandParamsR.slot("PatchFile"));
            costmaps = Rcpp::as<Rcpp::StringVector>(LandParamsR.slot("CostsFile"));
            name_landscape = habitatmaps(0);
            name_patch = patchmaps(0);
            name_costfile = costmaps(0);
        } else {
            name_landscape = Rcpp::as<string>(LandParamsR.slot("LandscapeFile"));
            name_patch = Rcpp::as<string>(LandParamsR.slot("PatchFile"));
            name_costfile = Rcpp::as<string>(LandParamsR.slot("CostsFile"));
        }
        if(!patchmodel && name_patch != "NULL") Rcpp::Rcout << "PatchFile must be NULL in a cell-based model!" << endl;

        name_sp_dist = Rcpp::as<string>(LandParamsR.slot("SpDistFile"));

        if(landtype == 2)
            ppLand.nHab = 1; // habitat quality landscape has one habitat class

        // CHECK IMPORTED RASTER FILES
        string indir = paramsSim->getDir(1);
        string ftype, fname;
        string filetype = "LandFile";
        //rasterdata patchraster, spdistraster; //declared globally

        // check landscape filename
        ftype = "LandscapeFile";
        fname = indir + name_landscape;
        landraster = ParseRasterHead(fname);
        if(landraster.ok) {
            if(landraster.cellsize == resolution)
                Rcpp::Rcout << ftype << " headers OK: " << fname << endl;
            else {
                errors++;
                Rcpp::Rcout << msgresol0 << ftype << " " << fname << msgresol1 << endl;
            }
        } else {
            errors++;
            if(landraster.errors == -111)
                OpenErrorR(ftype, fname);
            else
                FormatErrorR(fname, landraster.errors);
        }

        // check patch map filename
        ftype = "PatchFile";
        if(name_patch == "NULL") {
            if(patchmodel) {
                BatchErrorR(filetype, -999, 0, " ");
                errors++;
                Rcpp::Rcout << ftype << msgpatch << endl;
            }
        } else {
            if(patchmodel) {
                fname = indir + name_patch;
                patchraster = ParseRasterHead(fname);
                if(patchraster.ok) {
                    // check resolutions match
                    if(patchraster.cellsize == resolution) {
                        if(!errors) {
                            if(patchraster.cellsize == landraster.cellsize) {
                                // check that extent matches landscape extent
                                if(patchraster.ncols == landraster.ncols && patchraster.nrows == landraster.nrows) {
                                    // check origins match
                                    if((int)patchraster.xllcorner == (int)landraster.xllcorner &&
                                       (int)patchraster.yllcorner == (int)landraster.yllcorner) {
                                        Rcpp::Rcout << ftype << " headers OK: " << fname << endl;
                                    } else {
                                        Rcpp::Rcout << "*** Origin co-ordinates of " << ftype << msghdrs1 << endl;
                                        errors++;
                                    }
                                } else {
                                    Rcpp::Rcout << "*** Extent of " << ftype << " " << fname << msghdrs1 << endl;
                                    errors++;
                                }
                            } else {
                                Rcpp::Rcout << msgresol0 << ftype << " " << fname << msghdrs1 << endl;
                                errors++;
                            }
                        }
                    } else {
                        Rcpp::Rcout << msgresol0 << ftype << " " << fname << msgresol1 << endl;
                        errors++;
                    }
                } else {
                    errors++;
                    if(patchraster.errors == -111)
                        OpenErrorR(ftype, fname);
                    else
                        FormatErrorR(fname, patchraster.errors);
                }
            }
        }

        // check cost map filename
        ftype = "CostMapFile";
        if (name_costfile == "NULL") {
            if ( gTransferType == 1) { // SMS
                if (landtype == 2) { // habitat quality
                    BatchErrorR(filetype, -999, 0, " ");
                    errors++;
                    Rcpp::Rcout << ftype << " is required for a habitat quality landscape" << endl;
                }
            }
        }
        else {
            if ( gTransferType == 1) { // SMS
                fname = indir + name_costfile;
                costsraster = ParseRasterHead(fname);
                if(costsraster.ok) {
                    // check resolutions match
                    if(costsraster.cellsize == resolution) {
                        if(!errors) {
                            if(costsraster.cellsize == landraster.cellsize) {
                                // check that extent matches landscape extent
                                if(costsraster.ncols == landraster.ncols && costsraster.nrows == landraster.nrows) {
                                    // check origins match
                                    if((int)costsraster.xllcorner == (int)landraster.xllcorner &&
                                       (int)costsraster.yllcorner == (int)landraster.yllcorner) {
                                        Rcpp::Rcout << ftype << " headers OK: " << fname << endl;
                                    } else {
                                        Rcpp::Rcout << "*** Origin co-ordinates of " << ftype << msghdrs1 << endl;
                                        errors++;
                                    }
                                } else {
                                    Rcpp::Rcout << "*** Extent of " << ftype << " " << fname << msghdrs1 << endl;
                                    errors++;
                                }
                            } else {
                                Rcpp::Rcout << msgresol0 << ftype << " " << fname << msghdrs1 << endl;
                                errors++;
                            }
                        }
                    } else {
                        Rcpp::Rcout << msgresol0 << ftype << " " << fname << msgresol1 << endl;
                        errors++;
                    }
                } else {
                    errors++;
                    if(costsraster.errors == -111)
                        OpenErrorR(ftype, fname);
                    else
                        FormatErrorR(fname, costsraster.errors);
                }
            }
            else {
                BatchErrorR(filetype, -999, 0, " ");
                errors++;
                Rcpp::Rcout << ftype << " must be NULL if transfer model is not SMS" << endl;
            }
        }

        // check dynamic landscape filename
        // ...most checks are done at reading time in ReadDynLandR()
        ftype = "Dynamic landscape";
        if(ppLand.dynamic) {
            // check valid years
            if(dynland_years[0]!=0) {
                errors++;
                Rcpp::Rcout << "First year in dynamic landscape must be 0." << endl;
            } else {
                for(int i=1; i<dynland_years.size(); i++ ) {
                    if(dynland_years[i-1] >= dynland_years[i]) {
                        errors++;
                        Rcpp::Rcout << "Year in dynamic landscape must strictly increase." << endl;
                    }
                }
            }
            if(dynland_years.size() != habitatmaps.size()) {
                errors++;
                Rcpp::Rcout << "Dynamic landscape: Years must have as many elements as habitat maps." << endl;
            }
            if(patchmodel) {
                if( dynland_years.size() != patchmaps.size() ||
                    habitatmaps.size()   != patchmaps.size() ) {
                    errors++;
                    Rcpp::Rcout << "Dynamic landscape: Patchmaps must have as many elements as Years and habitat maps." << endl;
                }
            }
            if (name_costfile != "NULL") {
                if( dynland_years.size() != costmaps.size() ||
                    habitatmaps.size()   != costmaps.size() ) {
                    errors++;
                    Rcpp::Rcout << "Dynamic landscape: Costmaps must have as many elements as Years and habitat maps." << endl;
                }
            }
            if(errors==0) {
                // store land changes
                string landchangefile,patchchangefile;
                landChange chg;
                for(int i=1; i<dynland_years.size(); i++ ) {
                    chg.chgnum = i;
                    chg.chgyear = dynland_years[i];
                    chg.habfile = indir + habitatmaps(i);
                    if(patchmodel) chg.pchfile = indir + patchmaps(i);
                    else chg.pchfile = "NULL";
                    if (name_costfile == "NULL") chg.costfile = "none";
                    else chg.costfile = indir + costmaps(i);
                    pLandscape->addLandChange(chg);
                }
            }
        }

        // check initial distribution map filename
        ftype = "Species Distribution map";
        if(name_sp_dist == "NULL") {
            if(speciesdist) {
                BatchErrorR(filetype, -999, 0, " ");
                errors++;
                Rcpp::Rcout << ftype << " is required as SpeciesDist is 1 in Control" << endl;
            }
        } else {
            if(speciesdist) {
                fname = indir + name_sp_dist;
                spdistraster = ParseRasterHead(fname);
                if(spdistraster.ok) {
                    if(spdistraster.cellsize == distresolution) {
                        if(!errors) {
                            // check origins match
                            if((int)spdistraster.xllcorner == (int)landraster.xllcorner && (int)spdistraster.yllcorner == (int)landraster.yllcorner) {
                                // check extents match
                                if(spdistraster.cellsize == landraster.cellsize) {
                                    // same resolution
                                    if(spdistraster.ncols == landraster.ncols && spdistraster.nrows == landraster.nrows) {
                                        Rcpp::Rcout << ftype << " headers OK: " << fname << endl;
                                    } else {
                                        Rcpp::Rcout << "*** Extents of " << ftype << msghdrs1 << endl;
                                        errors++;
                                    }
                                } else { // different resolution
                                    if((spdistraster.cellsize % landraster.cellsize)==0) {
                                        double coarse = spdistraster.cellsize/landraster.cellsize;
                                        double rightx = ceil(landraster.ncols/coarse);
                                        double righty = ceil(landraster.nrows/coarse);
                                        if( spdistraster.ncols == int(rightx) && spdistraster.nrows == int(righty) ) {
                                            Rcpp::Rcout << ftype << " headers OK: " << fname << endl;
                                        } else {
                                            Rcpp::Rcout << "*** Extents of " << ftype << msghdrs1 << endl;
                                            errors++;
                                        }
                                    } else {
                                        BatchErrorR(filetype, -999, 0, "DistResolution");
                                        Rcpp::Rcout << "Resolution of initial distribution must be an integer multiple of Landscape resolution" << endl;
                                        errors++;
                                    }
                                }
                            } else {
                                Rcpp::Rcout << "*** Origin co-ordinates of " << ftype << msghdrs1 << endl;
                                errors++;
                            }
                        }
                    } else {
                        Rcpp::Rcout << msgresol0 << ftype << " " << fname << msgresol2 << endl;
                        errors++;
                    }
                } else {
                    errors++;
                    if(spdistraster.errors == -111)
                        OpenErrorR(ftype, fname);
                    else
                        FormatErrorR(fname, spdistraster.errors);
                }
            }
        }

    }

    pLandscape->setLandParams(ppLand, true);
    pLandscape->setGenLandParams(ppGenLand);

#if RSDEBUG
    // DEBUGLOG << "ReadLandParamsR(): NHab=" << ppLand.nHab << endl;
    DEBUGLOG << "ReadLandParamsR(): ppLand.landNum=" << ppLand.landNum << endl;
#endif

    if(errors) return false; //=landOK
    else return true;
}

//---------------------------------------------------------------------------

int ReadDynLandR(Landscape *pLandscape, Rcpp::S4 LandParamsR)
{

#if RSDEBUG
    DEBUGLOG << "ReadDynLandR(): pLandscape=" << pLandscape << endl;
#endif

    Rcpp::StringVector habitatmaps, patchmaps, costmaps;
    habitatmaps = Rcpp::as<Rcpp::StringVector>(LandParamsR.slot("LandscapeFile"));
    if (patchmodel) {
        patchmaps = Rcpp::as<Rcpp::StringVector>(LandParamsR.slot("PatchFile"));
    }
    bool costs = false;
    costmaps = Rcpp::as<Rcpp::StringVector>(LandParamsR.slot("CostsFile"));
    if (costmaps(0) != "NULL") costs = true;


    //------------ int ParseDynamicFile(string indir) {

    string indir = paramsSim->getDir(1);
    string fname,ftype;
    wstring header;
    wifstream hfile,pfile,cfile;
    int errors,ncols,nrows,cellsize,habnodata,pchnodata,costnodata;
    double xllcorner,yllcorner;

    errors = ncols = nrows = cellsize = habnodata = pchnodata = costnodata = 0;
    xllcorner = yllcorner = 0.0;

    if (patchmodel) {
        pLandscape->createPatchChgMatrix();
    }
    if (costs) {
        pLandscape->createCostsChgMatrix();
    }

    for(int i=1; i < habitatmaps.size(); i++ ) {

        // Habitat change file
        fname = indir + habitatmaps(i);

        // open file
#if RSWIN64
        hfile.open(fname.c_str());
#else
        hfile.open(fname, std::ios::binary);
#endif
        if(!hfile.is_open()) {
            OpenErrorR("Dynamic landscape habitat map ",  fname);
#if RSDEBUG
            DEBUGLOG << "Dynamic landscape habitat map failed to open: " << fname << std::endl;
#endif
            hfile.clear();
            return -212;
        } else {
#if RSDEBUG
            DEBUGLOG << "Dynamic landscape habitat map #" << i << " open to read" << std::endl;
#endif
#if !RSWIN64
            // check BOM for UTF-16
            if(check_bom(fname) == "utf16")
                // apply BOM-sensitive UTF-16 facet
                hfile.imbue(std::locale(hfile.getloc(), new std::codecvt_utf16<wchar_t, 0x10ffff, std::consume_header>));
#endif
            // ASCII header
            hfile >> header;
            if (!hfile.good()) {
#if RSDEBUG
                DEBUGLOG << "ReadDynLandR(): failed to read landscape habitat map #" << i << ": " << fname << std::endl;
#endif
                errors = -1112;
                hfile.close();
                hfile.clear();
                return errors;
            }
            if (header != L"ncols" && header != L"NCOLS") errors++;
            hfile >> ncols;

            hfile >> header >> nrows;
            if (header != L"nrows" && header != L"NROWS") errors++;

            hfile >> header >> xllcorner;
            if (header != L"xllcorner" && header != L"XLLCORNER") errors++;

            hfile >> header >> yllcorner;
            if (header != L"yllcorner" && header != L"YLLCORNER") errors++;
            double tmpcellsize;
            hfile >> header >> tmpcellsize;
            cellsize = (int) tmpcellsize;
            if (header != L"cellsize" && header != L"CELLSIZE") errors++;

            hfile >> header >> habnodata;
            if (header != L"NODATA_value" && header != L"NODATA_VALUE") errors++;

            if (errors > 0)  {
                FormatErrorR(fname,errors);
#if RSDEBUG
                DEBUGLOG << "ReadDynLandR(): failed to read Raster header of landscape habitat map #" << i << ": " << fname << std::endl;
#endif
                hfile.close();
                hfile.clear();
            } else {
                // check resolution match
                if (cellsize == resolution) {
                    // check that extent matches landscape extent
                    if(ncols == landraster.ncols && nrows == landraster.nrows) {
                        // check origins match
                        if((int)xllcorner == (int)landraster.xllcorner &&
                           (int)yllcorner == (int)landraster.yllcorner) {
                            Rcpp::Rcout << "Dynamic landscape habitat map #" << i << ", headers OK: " << fname << endl;
                        } else {
                            Rcpp::Rcout << "*** Origin co-ordinates of " << fname << msghdrs1 << endl;
                            errors++;
                        }
                    } else {
                        Rcpp::Rcout << "*** Extent of " << fname << msghdrs1 << endl;
                        errors++;
                    }
                } else {
                    Rcpp::Rcout << "*** Resolution of " << fname << msghdrs1 << endl;
                    errors++;
                }
                if (errors > 0)  {
                    FormatErrorR(fname,errors);
#if RSDEBUG
                    DEBUGLOG << "ReadDynLandR(): Errors in raster header of landscape habitat map #" << i << ": " << fname << std::endl;
#endif
                    hfile.close();
                    hfile.clear();
                }
            } // end of reading ASCII header
        }

        // Do the same for corresponding patch map, if applicable
        if (patchmodel) {
            // Patch change file
            fname = indir + patchmaps(i);

            // open file
#if RSWIN64
            pfile.open(fname.c_str());
#else
            pfile.open(fname, std::ios::binary);
#endif
            if(!pfile.is_open()) {
                OpenErrorR("Dynamic landscape patch map ",  fname);
#if RSDEBUG
                DEBUGLOG << "Dynamic landscape patch map failed to open: " << fname << std::endl;
#endif
                pfile.clear();
                return -213;
            } else {
#if RSDEBUG
                DEBUGLOG << "Dynamic landscape patch map #" << i << " open to read" << std::endl;
#endif
#if !RSWIN64
                // check BOM for UTF-16
                if(check_bom(fname) == "utf16")
                    // apply BOM-sensitive UTF-16 facet
                    pfile.imbue(std::locale(pfile.getloc(), new std::codecvt_utf16<wchar_t, 0x10ffff, std::consume_header>));
#endif
                // ASCII header
                pfile >> header;
                if (!pfile.good()) {
#if RSDEBUG
                    DEBUGLOG << "ReadDynLandR(): failed to read landscape patch map #" << i << ": " << fname << std::endl;
#endif
                    errors = -1113;
                    pfile.close();
                    pfile.clear();
                    return errors;
                }
                if (header != L"ncols" && header != L"NCOLS") errors++;
                pfile >> ncols;

                pfile >> header >> nrows;
                if (header != L"nrows" && header != L"NROWS") errors++;

                pfile >> header >> xllcorner;
                if (header != L"xllcorner" && header != L"XLLCORNER") errors++;

                pfile >> header >> yllcorner;
                if (header != L"yllcorner" && header != L"YLLCORNER") errors++;

                double tmpcellsize;
                pfile >> header >> tmpcellsize;
                cellsize = (int) tmpcellsize;
                if (header != L"cellsize" && header != L"CELLSIZE") errors++;

                pfile >> header >> pchnodata;
                if (header != L"NODATA_value" && header != L"NODATA_VALUE") errors++;

                if (errors > 0)  {
                    FormatErrorR(fname,errors);
#if RSDEBUG
                    DEBUGLOG << "ReadDynLandR(): failed to read Raster header of landscape patch map #" << i << ": " << fname << std::endl;
#endif
                    pfile.close();
                    pfile.clear();
                } else {
                    // check resolution match
                    if (cellsize == resolution) {
                        // check that extent matches landscape extent
                        if(ncols == landraster.ncols && nrows == landraster.nrows) {
                            // check origins match
                            if((int)xllcorner == (int)landraster.xllcorner &&
                               (int)yllcorner == (int)landraster.yllcorner) {
                                Rcpp::Rcout << "Dynamic landscape patch map #" << i << ", headers OK: " << fname << endl;
                            } else {
                                Rcpp::Rcout << "*** Origin co-ordinates of " << fname << msghdrs1 << endl;
                                errors++;
                            }
                        } else {
                            Rcpp::Rcout << "*** Extent of " << fname << msghdrs1 << endl;
                            errors++;
                        }
                    } else {
                        Rcpp::Rcout << "*** Resolution of " << fname << msghdrs1 << endl;
                        errors++;
                    }
                    if (errors > 0)  {
                        FormatErrorR(fname,errors);
#if RSDEBUG
                        DEBUGLOG << "ReadDynLandR(): Errors in raster header of landscape patch map #" << i << ": " << fname << std::endl;
#endif
                        pfile.close();
                        pfile.clear();
                    }
                } // end of reading ASCII header
            }
        } // end of if(patchmodel)
        else {
            pfile.clear();
        }

        // Do the same for corresponding cost map, if applicable
        if (costs) {
            // Cost change file
            fname = indir + costmaps(i);

            // open file
#if RSWIN64
            cfile.open(fname.c_str());
#else
            cfile.open(fname, std::ios::binary);
#endif
            if(!cfile.is_open()) {
                OpenErrorR("Dynamic SMS cost map ",  fname);
#if RSDEBUG
                DEBUGLOG << "Dynamic SMS cost map failed to open: " << fname << std::endl;
#endif
                cfile.clear();
                return -214;
            } else {
#if RSDEBUG
                DEBUGLOG << "Dynamic SMS cost map #" << i << " open to read" << std::endl;
#endif
#if !RSWIN64
                // check BOM for UTF-16
                if(check_bom(fname) == "utf16")
                    // apply BOM-sensitive UTF-16 facet
                    cfile.imbue(std::locale(cfile.getloc(), new std::codecvt_utf16<wchar_t, 0x10ffff, std::consume_header>));
#endif
                // ASCII header
                cfile >> header;
                if (!cfile.good()) {
#if RSDEBUG
                    DEBUGLOG << "ReadDynLandR(): failed to read SMS cost map #" << i << ": " << fname << std::endl;
#endif
                    errors = -1113;
                    cfile.close();
                    cfile.clear();
                    return errors;
                }
                if (header != L"ncols" && header != L"NCOLS") errors++;
                cfile >> ncols;

                cfile >> header >> nrows;
                if (header != L"nrows" && header != L"NROWS") errors++;

                cfile >> header >> xllcorner;
                if (header != L"xllcorner" && header != L"XLLCORNER") errors++;

                cfile >> header >> yllcorner;
                if (header != L"yllcorner" && header != L"YLLCORNER") errors++;

                double tmpcellsize;
                cfile >> header >> tmpcellsize;
                cellsize = (int) tmpcellsize;
                if (header != L"cellsize" && header != L"CELLSIZE") errors++;

                cfile >> header >> costnodata;
                if (header != L"NODATA_value" && header != L"NODATA_VALUE") errors++;

                if (errors > 0)  {
                    FormatErrorR(fname,errors);
#if RSDEBUG
                    DEBUGLOG << "ReadDynLandR(): failed to read Raster header of SMS cost map #" << i << ": " << fname << std::endl;
#endif
                    cfile.close();
                    cfile.clear();
                } else {
                    // check resolution match
                    if (cellsize == resolution) {
                        // check that extent matches landscape extent
                        if(ncols == landraster.ncols && nrows == landraster.nrows) {
                            // check origins match
                            if((int)xllcorner == (int)landraster.xllcorner &&
                               (int)yllcorner == (int)landraster.yllcorner) {
                                Rcpp::Rcout << "Dynamic SMS cost map #" << i << ", headers OK: " << fname << endl;
                            } else {
                                Rcpp::Rcout << "*** Origin co-ordinates of " << fname << msghdrs1 << endl;
                                errors++;
                            }
                        } else {
                            Rcpp::Rcout << "*** Extent of " << fname << msghdrs1 << endl;
                            errors++;
                        }
                    } else {
                        Rcpp::Rcout << "*** Resolution of " << fname << msghdrs1 << endl;
                        errors++;
                    }
                    if (errors > 0)  {
                        FormatErrorR(fname,errors);
#if RSDEBUG
                        DEBUGLOG << "ReadDynLandR(): Errors in raster header of SMS cost map #" << i << ": " << fname << std::endl;
#endif
                        cfile.close();
                        cfile.clear();
                    }
                } // end of reading ASCII header
            }
        } // end of if(costs)
        else {
            cfile.clear();
        }

        // Now read raster data of Habitat and, if applicable, Patch and/or Cost maps:
        int imported = 0;
        if (errors == 0)  {
            imported = pLandscape->readLandChange(i-1, costs, hfile, pfile, cfile, habnodata, pchnodata, costnodata);
            if (imported != 0) {
                if(hfile.is_open()) hfile.close();
                hfile.clear();
                if(patchmodel) {
                    if(pfile.is_open()) pfile.close();
                    pfile.clear();
                }
                if (costs) {
                    if(cfile.is_open()) cfile.close();
                    cfile.clear();
                }
                return imported;
            }
            if (patchmodel) {
                pLandscape->recordPatchChanges(i);
            }
            if (costs) {
                pLandscape->recordCostChanges(i);
            }
        }

        // Close files
        if(hfile.is_open()) hfile.close();
        hfile.clear();
        if(patchmodel) {
            if(pfile.is_open()) pfile.close();
            pfile.clear();
        }
        if (costs) {
            if(cfile.is_open()) cfile.close();
            cfile.clear();
        }
    } // end of loop over landscape changes i

    if(patchmodel) {
        // record changes back to original landscape for multiple replicates
        pLandscape->recordPatchChanges(0);
        pLandscape->deletePatchChgMatrix();
    }
    if (costs) {
        pLandscape->recordCostChanges(0);
        pLandscape->deleteCostsChgMatrix();
    }
#if RSDEBUG
    DEBUGLOG << "ReadDynLandR(): finished" << endl;
#endif
    return 0;
}


//---------------------------------------------------------------------------

int ReadParametersR(Landscape* pLandscape, Rcpp::S4 ParMaster)
{
    Rcpp::S4 ParamParamsR("SimulationParams");
    ParamParamsR = Rcpp::as<Rcpp::S4>(ParMaster.slot("simul"));

    Rcpp::S4 DemogParamsR("DemogParams");
    DemogParamsR = Rcpp::as<Rcpp::S4>(ParMaster.slot("demog"));

    Rcpp::S4 LandParamsR("LandParams");
    LandParamsR = Rcpp::as<Rcpp::S4>(ParMaster.slot("land"));

    int error = 0;
    landParams paramsLand = pLandscape->getLandParams();
    envStochParams env = paramsStoch->getStoch();
    demogrParams dem = pSpecies->getDemogrParams();
    simParams sim = paramsSim->getSim();

    sim.simulation = Rcpp::as<int>(ParamParamsR.slot("Simulation"));
    sim.reps = Rcpp::as<int>(ParamParamsR.slot("Replicates"));
    sim.years = Rcpp::as<int>(ParamParamsR.slot("Years"));
#if RSDEBUG
    DEBUGLOG << "ReadParametersR(): paramsSim = " << paramsSim << endl;
    DEBUGLOG << "ReadParametersR(): simulation = " << sim.simulation << " reps = " << sim.reps
             << " years = " << sim.years << endl;
#endif

    int iiii, gradType, shift_begin, shift_stop;
    float grad_inc, opt_y, f, optEXT, shift_rate;
    bool shifting;

    // Landscape boundary and any 'no-data' areas are absorbing?:
    sim.absorbing = Rcpp::as<bool>(ParamParamsR.slot("Absorbing"));

    // Environmental gradient: 0 = none, 1  = carrying capacity (or 1/b), 2 = growth rate (or fecundity), 3 = local
    // extinction probability. Must be 0 for patch-based models. NOTE: change in code numbers from v1.0
    gradType = Rcpp::as<int>(ParamParamsR.slot("Gradient"));
    // following not required for Gradient = 0
    grad_inc = Rcpp::as<float>(ParamParamsR.slot("GradSteep"));
    opt_y = Rcpp::as<float>(ParamParamsR.slot("Optimum"));
    f = Rcpp::as<float>(ParamParamsR.slot("f"));
    optEXT = Rcpp::as<float>(ParamParamsR.slot("ExtinctOptim"));
    paramsGrad->setGradient(gradType, grad_inc, opt_y, f, optEXT);

    shifting = Rcpp::as<bool>(ParamParamsR.slot("Shifting"));
    shift_rate = Rcpp::as<float>(ParamParamsR.slot("ShiftRate"));
    shift_begin = Rcpp::as<int>(ParamParamsR.slot("ShiftStart"));
    shift_stop = Rcpp::as<int>(ParamParamsR.slot("ShiftEnd"));
    if(shifting)
        paramsGrad->setShifting(shift_rate, shift_begin, shift_stop);
    else
        paramsGrad->noShifting();

    // Environmental stochasticity: 0 = none, 1 = global, 2 = local:
    iiii = Rcpp::as<int>(ParamParamsR.slot("EnvStoch"));
    if(iiii == 0)
        env.stoch = false;
    else {
        env.stoch = true;
        if(iiii == 2)
            env.local = true;
        else
            env.local = false;
    }
    // For a patch-based model, EnvStoch = 2 is NOT allowed:
    if(paramsLand.patchModel && env.local)
        error = 101;

    // Environmental stochasticity type: 0 = in growth rate, 1 = in carrying capacity
    // Not required for EnvStoch = 0; stochasticity in carrying capacity is allowed for an artificial landscape only
    env.inK = Rcpp::as<bool>(ParamParamsR.slot("EnvStochType"));

    // as from v1.1, there is just one pair of min & max values,
    // which are attributes of the species
    // ULTIMATELY, THE PARAMETER FILE SHOULD HAVE ONLY TWO COLUMNS ...
    env.ac = Rcpp::as<float>(ParamParamsR.slot("ac")); // Temporal autocorrelation coefficient
    env.std = Rcpp::as<float>(ParamParamsR.slot(
        "std")); // Amplitude of stochastic fluctuations: standard deviation of a normal distribution (having mean 0)
    float minR, maxR, minK, maxK;
    minR = Rcpp::as<float>(ParamParamsR.slot("minR")); // not required for EnvStoch = 0 or EnvStochType = 1
    maxR = Rcpp::as<float>(ParamParamsR.slot("maxR")); // -"-
    minK = Rcpp::as<float>(ParamParamsR.slot("minK")); // not required for EnvStoch = 0 or EnvStochType = 0
    maxK = Rcpp::as<float>(ParamParamsR.slot("maxK")); // -"-
    if(env.inK) {
        float minKK, maxKK;
        minKK = minK * (((float)(pow(paramsLand.resol, double(2)))) / 10000.0);
        maxKK = maxK * (((float)(pow(paramsLand.resol, double(2)))) / 10000.0);
        pSpecies->setMinMax(minKK, maxKK);
    } else
        pSpecies->setMinMax(minR, maxR);

    env.localExt = Rcpp::as<bool>(ParamParamsR.slot("LocalExt"));
    // For a patch-based model no LocalExt allowed:
    if(paramsLand.patchModel && env.localExt)
        error = 102;
    env.locExtProb = Rcpp::as<float>(ParamParamsR.slot("LocalExtProb"));
    paramsStoch->setStoch(env);

    dem.propMales = Rcpp::as<float>(DemogParamsR.slot("PropMales"));
    dem.harem = Rcpp::as<float>(DemogParamsR.slot("Harem"));
    dem.bc = Rcpp::as<float>(DemogParamsR.slot("bc"));
    dem.lambda = Rcpp::as<float>(DemogParamsR.slot("Rmax"));

    pSpecies->setDemogr(dem);

    float k;
    if(landtype == 9) { // artificial landscape
        // only one value of K is read, but the first 'habitat' is the matrix where K = 0
        pSpecies->createHabK(2);
        k = Rcpp::as<float>(LandParamsR.slot("K_or_DensDep"));
        k *= ((double)(pow(paramsLand.resol, double(2)))) / 10000.0;
        pSpecies->setHabK(0,0);
        pSpecies->setHabK(1,k);
    } else {
        pSpecies->createHabK(paramsLand.nHabMax);
        Rcpp::NumericVector k_vec;
        k_vec = Rcpp::as<Rcpp::NumericVector>(LandParamsR.slot("K_or_DensDep"));
        for (int i = 0; i < paramsLand.nHabMax; i++) {
            k = k_vec[i] * ((double)(pow(paramsLand.resol, double(2)))) / 10000.0;
            pSpecies->setHabK(i,k);
        }
    }

#if RSDEBUG
    DEBUGLOG << "ReadParametersR(): dem.lambda = " << dem.lambda
             << " habK[0] = " << pSpecies->getHabK(0)
             << " nHabMax = " << paramsLand.nHabMax << endl;
#endif


    // Output start years
    sim.outStartPop = Rcpp::as<int>(ParamParamsR.slot("OutStartPop"));
    sim.outStartInd = Rcpp::as<int>(ParamParamsR.slot("OutStartInd"));
    sim.outStartTraitCell = Rcpp::as<int>(ParamParamsR.slot("OutStartTraitCell"));
    sim.outStartTraitRow = Rcpp::as<int>(ParamParamsR.slot("OutStartTraitRow"));
    sim.outStartConn = Rcpp::as<int>(ParamParamsR.slot("OutStartConn"));
    sim.outStartPaths = Rcpp::as<int>(ParamParamsR.slot("OutStartPaths"));
    // Output intervals
    sim.outIntRange = Rcpp::as<int>(ParamParamsR.slot("OutIntRange"));
    sim.outIntOcc = Rcpp::as<int>(ParamParamsR.slot("OutIntOcc"));
    sim.outIntPop = Rcpp::as<int>(ParamParamsR.slot("OutIntPop"));
    sim.outIntInd = Rcpp::as<int>(ParamParamsR.slot("OutIntInd"));

    if(sim.outIntRange > 0)
        sim.outRange = true;
    else
        sim.outRange = false;
    if(sim.outIntOcc > 0)
        sim.outOccup = true;
    else
        sim.outOccup = false;
    if(sim.outIntPop > 0)
        sim.outPop = true;
    else
        sim.outPop = false;
    if(sim.outIntInd > 0)
        sim.outInds = true;
    else
        sim.outInds = false;

    sim.outIntTraitCell = Rcpp::as<int>(ParamParamsR.slot("OutIntTraitCell"));
    sim.outIntTraitRow = Rcpp::as<int>(ParamParamsR.slot("OutIntTraitRow"));
    sim.outIntConn = Rcpp::as<int>(ParamParamsR.slot("OutIntConn"));
    sim.outIntPaths = Rcpp::as<int>(ParamParamsR.slot("OutIntPaths"));
    if(sim.outIntTraitCell > 0)
        sim.outTraitsCells = true;
    else
        sim.outTraitsCells = false;
    if(sim.outIntTraitRow > 0)
        sim.outTraitsRows = true;
    else
        sim.outTraitsRows = false;
    if(sim.outIntConn > 0)
        sim.outConnect = true;
    else
        sim.outConnect = false;
    if(sim.outIntPaths > 0)
        sim.outPaths = true;
    else
        sim.outPaths = false;

    if(sim.outOccup && sim.reps < 2)
        error = 103;
    if(paramsLand.patchModel) {
        if(sim.outTraitsRows)
            error = 104;
    } else {
        if(sim.outConnect)
            error = 105;
    }
#if RSDEBUG
    DEBUGLOG << "ReadParametersR(): outRange = " << sim.outRange << " outInt = " << sim.outIntRange << endl;
#endif

    sim.saveMaps = Rcpp::as<bool>(ParamParamsR.slot("SaveMaps"));
    sim.mapInt = Rcpp::as<int>(ParamParamsR.slot("MapsInterval"));
    sim.saveVisits = Rcpp::as<bool>(ParamParamsR.slot("SMSHeatMap"));
    sim.drawLoaded = Rcpp::as<bool>(ParamParamsR.slot("DrawLoadedSp"));
    // sim.saveInitMap = false;
#if RS_RCPP
    sim.ReturnPopRaster = Rcpp::as<bool>(ParamParamsR.slot("ReturnPopRaster"));
    sim.CreatePopFile = Rcpp::as<bool>(ParamParamsR.slot("CreatePopFile"));
#endif

    paramsSim->setSim(sim);

    return error;
}

//---------------------------------------------------------------------------


int ReadStageStructureR(Rcpp::S4 ParMaster)
{
    Rcpp::S4 DemogParamsR("DemogParams");
    DemogParamsR = Rcpp::as<Rcpp::S4>(ParMaster.slot("demog"));
    Rcpp::S4 StagesParamsR("StagesParams");
    StagesParamsR = Rcpp::as<Rcpp::S4>(DemogParamsR.slot("StageStruct"));

    demogrParams dem = pSpecies->getDemogrParams();
    stageParams sstruct = pSpecies->getStageParams();
    int matrixsize, i, j, stg;
    float ss, dd, devCoeff, survCoeff;
    Rcpp::NumericMatrix trmatrix, wtsmatrix;
    Rcpp::IntegerVector minAge;

#if RSDEBUG
    DEBUGLOG << "ReadStageStructureR(): sstruct.nStages = " << sstruct.nStages << endl;
#endif
    // int simulation;
    // simulation = Rcpp::as<int>(StagesParamsR.slot("Simulation")); //Must match simulation numbers in ParamParams
    sstruct.disperseOnLoss = Rcpp::as<bool>(StagesParamsR.slot("PostDestructn"));
    sstruct.probRep = Rcpp::as<float>(StagesParamsR.slot("PRep"));
    sstruct.repInterval = Rcpp::as<int>(StagesParamsR.slot("RepInterval"));
    sstruct.maxAge = Rcpp::as<int>(StagesParamsR.slot("MaxAge"));
    trmatrix = Rcpp::as<Rcpp::NumericMatrix>(StagesParamsR.slot("TransMatrix"));
    // parsing is done on R-level
    minAge = Rcpp::as<Rcpp::IntegerVector>(StagesParamsR.slot("MinAge"));

    // Store Transition matrix:
    if(dem.repType != 2) { // asexual or implicit sexual model
        matrixsize = sstruct.nStages;
        for(i = 0; i < matrixsize; i++) {
            // set minimum ages
            pSpecies->setMinAge(i, 0, minAge(i));
            // set fecundities
            pSpecies->setFec(i, 0, (float)trmatrix(0, i));
            // set survival and development probalities
            ss = (float)trmatrix(i, i); // survival prob
            if((i + 1) != matrixsize) {
                dd = (float)trmatrix(i + 1, i);
            } // development prob
            else {
                dd = 0.0;
            }
            pSpecies->setSurv(i, 0, ss + dd);
            pSpecies->setDev(i, 0, dd / (ss + dd));
        }
    } else {                              // complex sexual model
        matrixsize = sstruct.nStages * 2; // juv (i=0) sex-independent as ending stage -> columns = rows + 1
        stg = 1;
        for(i = 1; i < (matrixsize - 1); i++) { // loop over rows
            // set minAge and fecundities
            if(i % 2) {                         // odd columns  -> males (1)
                pSpecies->setMinAge(stg, 1, minAge(i));
                pSpecies->setFec(stg, 1, (float)trmatrix(0, i + 1));
            } else { // even columns -> females (0)
                pSpecies->setMinAge(stg, 0, minAge(i));
                pSpecies->setFec(stg, 0, (float)trmatrix(0, i + 1));
                stg++;
            }
        }
        stg = 0;
        for(i = 0; i < matrixsize; i++) { // loop over columns
            // survival and development
            if(i != 0)
                ss = (float)trmatrix(i - 1, i); // survival prob
            else
                ss = (float)trmatrix(i, i);
            if((i + 2)  != matrixsize && (i + 1) != matrixsize)
                dd = (float)trmatrix(i + 1, i); // development prob
            else
                dd = 0.0;
            if(i % 2 == 0) { // even rows -> males (1)
                pSpecies->setSurv(stg, 1, ss + dd);
                pSpecies->setDev(stg, 1, dd / (ss + dd));
            } else { // odd rows -> females (0)
                pSpecies->setSurv(stg, 0, ss + dd);
                pSpecies->setDev(stg, 0, dd / (ss + dd));
                stg++;
            }
        }
    }
#if RSDEBUG
    DEBUGLOG << "Read_Transition Matrix: matrix = " << trmatrix << endl;
#endif

    // Survival schedule: 0 = At reproduction; 1 = Between reproductive events; 2 = Annually
    sstruct.survival = Rcpp::as<int>(StagesParamsR.slot("SurvSched"));

    if(dem.repType != 2)
        matrixsize = sstruct.nStages;
    else
        matrixsize = sstruct.nStages * gMaxNbSexes;
    // Fecundity
    sstruct.fecDens = Rcpp::as<bool>(StagesParamsR.slot("FecDensDep")); // Density-dependence in reproduction
    sstruct.fecStageDens =
        Rcpp::as<bool>(StagesParamsR.slot("FecStageWts")); // stage-specific density dependence of fecundity
    if(sstruct.fecStageDens) {
        wtsmatrix = Rcpp::as<Rcpp::NumericMatrix>(StagesParamsR.slot("FecStageWtsMatrix"));
        pSpecies->createDDwtFec(matrixsize);
        for(i = 0; i < matrixsize; i++) {
            for(j = 0; j < matrixsize; j++) {
                pSpecies->setDDwtFec(i, j, wtsmatrix(i, j));
            }
        }
#if RSDEBUG
        DEBUGLOG << "Read_StageWeights(): completed reading fecundity weights matrix " << endl;
#endif
    }

    // Development
    sstruct.devDens = Rcpp::as<bool>(StagesParamsR.slot("DevDensDep")); // Density-dependence in development
    devCoeff = Rcpp::as<float>(StagesParamsR.slot("DevDensCoeff"));
    sstruct.devStageDens =
        Rcpp::as<bool>(StagesParamsR.slot("DevStageWts")); // stage-specific density dependence of development
    if(sstruct.devStageDens) {
        wtsmatrix = Rcpp::as<Rcpp::NumericMatrix>(StagesParamsR.slot("DevStageWtsMatrix"));
        pSpecies->createDDwtDev(matrixsize);
        for(i = 0; i < matrixsize; i++) {
            for(j = 0; j < matrixsize; j++) {
                pSpecies->setDDwtDev(i, j, wtsmatrix(i, j));
            }
        }
#if RSDEBUG
        DEBUGLOG << "Read_StageWeights(): completed reading development weights matrix " << endl;
#endif
    }

    // Survival
    sstruct.survDens = Rcpp::as<bool>(StagesParamsR.slot("SurvDensDep")); // Density-dependence in survival
    survCoeff = Rcpp::as<float>(StagesParamsR.slot("SurvDensCoeff"));
    sstruct.survStageDens =
        Rcpp::as<bool>(StagesParamsR.slot("SurvStageWts")); // stage-specific density dependence of survival
    if(sstruct.survStageDens) {
        wtsmatrix = Rcpp::as<Rcpp::NumericMatrix>(StagesParamsR.slot("SurvStageWtsMatrix"));
        pSpecies->createDDwtSurv(matrixsize);
        for(i = 0; i < matrixsize; i++) {
            for(j = 0; j < matrixsize; j++) {
                pSpecies->setDDwtSurv(i, j, wtsmatrix(i, j));
            }
        }
#if RSDEBUG
        DEBUGLOG << "Read_StageWeights(): completed reading survival weights matrix " << endl;
#endif
    }

    pSpecies->setStage(sstruct);
    if(sstruct.devDens || sstruct.survDens) {
        pSpecies->setDensDep(devCoeff, survCoeff);
    }

    return 0;
}

//---------------------------------------------------------------------------

int ReadEmigrationR(Rcpp::S4 ParMaster)
{

    Rcpp::S4 DispParamsR("DispersalParams");
    DispParamsR = Rcpp::as<Rcpp::S4>(ParMaster.slot("dispersal"));
    Rcpp::S4 EmigParamsR("EmigrationParams");
    EmigParamsR = Rcpp::as<Rcpp::S4>(DispParamsR.slot("Emigration"));

    int error = 0;
    int Nlines, emigstage, stage, sex, offset;
    Rcpp::NumericMatrix EmigMatrix;
    Rcpp::NumericVector EmigScalesVec;
    demogrParams dem = pSpecies->getDemogrParams();
    stageParams sstruct = pSpecies->getStageParams();
    emigRules emig = pSpecies->getEmigRules();
    emigTraits emigrationTraits;

    emig.densDep = Rcpp::as<bool>(EmigParamsR.slot("DensDep"));
    pSpecies->setFullKernel(Rcpp::as<bool>(EmigParamsR.slot("UseFullKern"))); // Emigration rate derived from kernel. Only for kernel-based transfer and DensDep = 0
    emig.stgDep = Rcpp::as<bool>(EmigParamsR.slot("StageDep")); // Stage-dependent emigration. Must be 0 if IndVar is 1
    emig.sexDep = Rcpp::as<bool>(EmigParamsR.slot("SexDep"));   // Sex-dependent emigration.
    emig.indVar =
        Rcpp::as<bool>(EmigParamsR.slot("IndVar")); // Inter-individual variability. Must be 0 if StageDep is 1
    if(stagestruct && emig.indVar){
        emigstage = Rcpp::as<int>(EmigParamsR.slot("EmigStage")); // Stage which emigrates. Required for stage-strucutred population having IndVar = 1
        if(emigstage >= 0 && emigstage < sstruct.nStages)
            emig.emigStage = emigstage;
        else
            emig.emigStage = 0;
    } else {
        emig.emigStage = 0;
    }

    pSpecies->setEmigRules(emig);

    if(!dem.repType && emig.sexDep)
        error = 301;
    if(!dem.stageStruct && emig.stgDep)
        error = 303;

    // no.of lines according to known stage- and sex-dependency and corresponding column offset
    if(emig.stgDep) {
        if(emig.sexDep) {
            Nlines = sstruct.nStages * gNbSexesDisp;
            offset = 2;
        } else {
            Nlines = sstruct.nStages;
            offset = 1;
        }
    } else {
        if(emig.sexDep) {
            Nlines = gNbSexesDisp;
            offset = 1;
        } else {
            Nlines = 1;
            offset = 0;
        }
    }

#if RSDEBUG
    DEBUGLOG << "ReadEmigrationR(): Nlines = " << Nlines << " emig.densDep = " << emig.densDep
             << " emig.indVar = " << emig.indVar << " sexesDisp = " << sexesDisp << endl;
#endif

    if (!emig.indVar){
        EmigMatrix = Rcpp::as<Rcpp::NumericMatrix>(EmigParamsR.slot("EmigProb"));

        for(int line = 0; line < Nlines; line++) {

        if(emig.stgDep) {
            if(emig.sexDep) {
                stage = (int)EmigMatrix(line, 0);
                sex = (int)EmigMatrix(line, 1);
            } else {
                stage = (int)EmigMatrix(line, 0);
                sex = 0;
            }
        } else {
            if(emig.sexDep) {
                stage = 0;
                sex = (int)EmigMatrix(line, 0);
            } else {
                stage = 0;
                sex = 0;
            }
        }

        if(emig.densDep) {
            emigrationTraits.d0 = (float)EmigMatrix(line, offset + 0);
            emigrationTraits.alpha = (float)EmigMatrix(line, offset + 1);
            emigrationTraits.beta = (float)EmigMatrix(line, offset + 2);

            pSpecies->setSpEmigTraits(stage, sex, emigrationTraits);
        } else {
            emigrationTraits.d0 = (float)EmigMatrix(line, offset + 0);
            emigrationTraits.alpha = emigrationTraits.beta = 0.0;

            pSpecies->setSpEmigTraits(stage, sex, emigrationTraits);
            }
        } // end of Nlines for loop

    } else{
        gHasGenetics = true;

        // set the emigration probability parameters to -9
        if(emig.stgDep) { // it should not go in here because it cannot be stage dependent for individual variability
            if(emig.sexDep) {
                for (int i = 0; i < sstruct.nStages; i++) {
                    for (int j = 0; j < gNbSexesDisp; j++) {
                        if(emig.densDep){
                            emigrationTraits.d0 = -9.0;
                            emigrationTraits.alpha = -9.0;
                            emigrationTraits.beta = -9.0;
                            pSpecies->setSpEmigTraits(i, j, emigrationTraits);
                        } else{
                            emigrationTraits.d0 = -9.0;
                            emigrationTraits.alpha = emigrationTraits.beta = 0.0;
                            pSpecies->setSpEmigTraits(i, j, emigrationTraits);
                        }
                    }
                }
            } else {
                for (int i = 0; i < sstruct.nStages; i++) {
                    if(emig.densDep){
                        emigrationTraits.d0 = -9.0;
                        emigrationTraits.alpha = -9.0;
                        emigrationTraits.beta = -9.0;
                        pSpecies->setSpEmigTraits(i, 0, emigrationTraits);
                    } else{
                        emigrationTraits.d0 = -9.0;
                        emigrationTraits.alpha = emigrationTraits.beta = 0.0;
                        pSpecies->setSpEmigTraits(i, 0, emigrationTraits);
                    }

                }
            }
        } else {
            if(emig.sexDep) {
                for (int j = 0; j < gNbSexesDisp; j++) {
                    if(emig.densDep){
                        emigrationTraits.d0 = -9.0;
                        emigrationTraits.alpha = -9.0;
                        emigrationTraits.beta = -9.0;
                        pSpecies->setSpEmigTraits(0, j, emigrationTraits);
                    } else{
                        emigrationTraits.d0 = -9.0;
                        emigrationTraits.alpha = emigrationTraits.beta = 0.0;
                        pSpecies->setSpEmigTraits(0, j, emigrationTraits);
                    }
                }
            } else {
                    if(emig.densDep){
                        emigrationTraits.d0 = -9.0;
                        emigrationTraits.alpha = -9.0;
                        emigrationTraits.beta = -9.0;
                        pSpecies->setSpEmigTraits(0, 0, emigrationTraits);
                    } else{
                        emigrationTraits.d0 = -9.0;
                        emigrationTraits.alpha = emigrationTraits.beta = 0.0;
                        pSpecies->setSpEmigTraits(0, 0, emigrationTraits);
                    }
            }
        }
    } // if indVar

    return error;
}

//---------------------------------------------------------------------------

int ReadTransferR(Landscape* pLandscape, Rcpp::S4 ParMaster)
{
    Rcpp::S4 DispParamsR("DispersalParams");
    DispParamsR = Rcpp::as<Rcpp::S4>(ParMaster.slot("dispersal"));

    int Nlines, offset, stage, sex;
    int error = 0;

    landParams paramsLand = pLandscape->getLandParams();
    demogrParams dem = pSpecies->getDemogrParams();
    stageParams sstruct = pSpecies->getStageParams();
    transferRules trfr = pSpecies->getTransferRules();

#if RSDEBUG
    DEBUGLOG << "ReadTransferR(): TransferType=" << trfr.moveModel
             << " paramsLand.generated=" << paramsLand.generated
             << " paramsLand.rasterType=" << paramsLand.rasterType
             << " trfr.moveModel=" << trfr.moveModel
             << " trfr.twinKern=" << trfr.twinKern
             << endl;
#endif

    // Create Costs vector of species
    if(trfr.usesMovtProc) {
#if RSDEBUG
        DEBUGLOG << "ReadTransferR(): creating cost/mortality matrix, dimension=";
        if(paramsLand.generated)
            DEBUGLOG << paramsLand.nHab;
        else
            DEBUGLOG << paramsLand.nHabMax;
        DEBUGLOG << endl;
#endif
        if(paramsLand.generated) {
            pSpecies->createHabCostMort(paramsLand.nHab);
        } else {
            pSpecies->createHabCostMort(paramsLand.nHabMax);
        }
    }

    int TransferType; // new local variable to replace former global variable
    if(trfr.usesMovtProc)
        TransferType = trfr.moveType;
    else
        TransferType = 0;

    // trfrKernTraits k; // andere Implementierung für individuelle Variabilität
    // trfrMovtTraits movt; // andere Implementierung für individuelle Variabilität
    trfrKernelParams kparams;
    trfrMortParams mort;
    // trfrScales scale; // andere Implementierung für individuelle Variabilität
    string CostsFile;
    trfrMovtParams smsparams;
    trfrMovtParams mparams;
    Rcpp::NumericVector ReadVec;

    switch(TransferType) {
    case 0: { // dispersal kernel

        Rcpp::S4 TransParamsR("DispersalKernel");
        TransParamsR = Rcpp::as<Rcpp::S4>(DispParamsR.slot("Transfer"));

        Rcpp::NumericMatrix DispMatrix;
        // Rcpp::NumericVector DispScalesVec;

        // simulation = Rcpp::as<int>(TransParamsR.slot("Simulation")); // REMOVED in R-interface //Must match
        // simulation numbers in ParamParams
        trfr.stgDep = Rcpp::as<bool>(TransParamsR.slot("StageDep")); // Stage-dependent transfer. Must be 0 if IndVar is 1
        trfr.sexDep = Rcpp::as<bool>(TransParamsR.slot("SexDep")); // Sex-dependent transfer.
        trfr.twinKern = Rcpp::as<bool>(TransParamsR.slot("DoubleKernel")); // 0 = negative exponential; 1 = double negative exponential
        trfr.distMort = Rcpp::as<bool>(TransParamsR.slot("DistMort")); // Distance-dependent mortality
        trfr.indVar = Rcpp::as<bool>(TransParamsR.slot("IndVar"));

        // some error checks
        if(dem.repType == 0) {
            if(trfr.sexDep) {
                error = 401;
            }
        }
        if(dem.stageStruct) {
        } // if (trfr.indVar) error = 402;
        else {
            if(trfr.stgDep) {
                error = 403;
            }
        }

        // set no. of lines according to known stage- and sex-dependency and corresponding column offset
        if(trfr.stgDep) {
            if(trfr.sexDep) {
                Nlines = sstruct.nStages * gNbSexesDisp;
                offset = 2;
            } else {
                Nlines = sstruct.nStages;
                offset = 1;
            }
        } else {
            if(trfr.sexDep) {
                Nlines = gNbSexesDisp;
                offset = 1;
            } else {
                Nlines = 1;
                offset = 0;
            }
        }

        if(trfr.indVar){
            // parameters for individual variability are set in ReadTraits()
            stage = 0;
            kparams.meanDist1 = -9;
            kparams.meanDist2 = -9;
            kparams.probKern1 = -9;
            if(trfr.sexDep) {
                for (sex=0; sex<gNbSexesDisp; sex++) {
                    pSpecies->setSpKernTraits(stage, sex, kparams, paramsLand.resol);
                }
            } else {
                pSpecies->setSpKernTraits(stage, 0, kparams, paramsLand.resol);
            }
        }
        else {
            DispMatrix = Rcpp::as<Rcpp::NumericMatrix>(TransParamsR.slot("Distances")); // only if not indVar

            for(int line = 0; line < Nlines; line++) {
                if(trfr.stgDep) {
                    if(trfr.sexDep) {
                        stage = (int)DispMatrix(line, 0);
                        sex = (int)DispMatrix(line, 1);
                    } else {
                        stage = (int)DispMatrix(line, 0);
                        sex = 0;
                    }
                } else {
                    if(trfr.sexDep) {
                        stage = 0;
                        sex = (int)DispMatrix(line, 0);
                    } else {
                        stage = 0;
                        sex = 0;
                    }
                }

                if(trfr.twinKern) {
                        kparams.meanDist1 = (float)DispMatrix(line, offset + 0);
                        kparams.meanDist2 = (float)DispMatrix(line, offset + 1);
                        kparams.probKern1 = (float)DispMatrix(line, offset + 2);

                        pSpecies->setSpKernTraits(stage, sex, kparams, paramsLand.resol);
                } else { // single kernel
                        kparams.meanDist1 = (float)DispMatrix(line, offset + 0);
                        kparams.meanDist2 = kparams.meanDist1;
                        kparams.probKern1 = 1.0;

                        pSpecies->setSpKernTraits(stage, sex, kparams, paramsLand.resol);
                }
            } // end of Nlines for-loop
        }

        // mortality
        mort.fixedMort = Rcpp::as<float>(TransParamsR.slot("MortProb"));
        mort.mortAlpha = Rcpp::as<float>(TransParamsR.slot("Slope"));
        mort.mortBeta  = Rcpp::as<float>(TransParamsR.slot("InflPoint"));
        pSpecies->setMortParams(mort);

        pSpecies->setTrfrRules(trfr);

    }
        break; // end of dispersal kernel

    case 1: { // SMS

        Rcpp::S4 TransParamsR("StochMove");
        TransParamsR = Rcpp::as<Rcpp::S4>(DispParamsR.slot("Transfer"));

        // simulation = Rcpp::as<int>(TransParamsR.slot("Simulation")); // REMOVED in R-interface // Must match simulation numbers in ParamParams
        smsparams.pr = Rcpp::as<short>(TransParamsR.slot("PR")); // Perceptual range (cells)
        smsparams.prMethod = Rcpp::as<short>(TransParamsR.slot("PRMethod")); // Perceptual range method: 1 = arithmetic mean; 2 = harmonic mean; 3 = weighted arithmtic mean
        smsparams.memSize = Rcpp::as<short>(TransParamsR.slot("MemSize"));    // No. of previous steps over which to calculate current direction to apply DP [1-14]
        smsparams.goalType = Rcpp::as<short>(TransParamsR.slot("GoalType"));   // Goal type: 0 (none) or 2 (dispersal bias)
        trfr.indVar = Rcpp::as<bool>(TransParamsR.slot("IndVar"));

        if(trfr.indVar) {
            smsparams.dp = -9;
        } else {
            ReadVec = Rcpp::as<Rcpp::NumericVector>(TransParamsR.slot("DP")); // Directional persistence, Must be >= 1.0
            if(ReadVec.size() == 1) {
                smsparams.dp = (float)ReadVec[0];
            } else {
                error = 436;
            }
        }
        if(trfr.indVar) {
            smsparams.gb = -9;
        } else {
            ReadVec = Rcpp::as<Rcpp::NumericVector>(TransParamsR.slot("GoalBias")); // Goal bias strength, Must be >= 1.0
            if(ReadVec.size() == 1) {
                smsparams.gb = (float)ReadVec[0];
            } else {
                error = 436;
            }
        }
        if(smsparams.goalType == 2) {	// dispersal bias
            if(trfr.indVar) {
                smsparams.alphaDB = -9;
            } else {
                ReadVec = Rcpp::as<Rcpp::NumericVector>(TransParamsR.slot("AlphaDB")); // Dispersal bias decay rate (> 0)
                if(ReadVec.size() == 1) {
                    smsparams.alphaDB = (float)ReadVec[0];
                } else {
                    error = 436;
                }
            }
            if(trfr.indVar) {
                smsparams.betaDB = -9;
            } else {
                ReadVec = Rcpp::as<Rcpp::NumericVector>(TransParamsR.slot("BetaDB")); // Dispersal bias decay inflection point (no. of steps) (> 0)
                if(ReadVec.size() == 1) {
                    smsparams.betaDB = (float)ReadVec[0];
                } else {
                    error = 436;
                }
            }
        }
#if RSDEBUG
        DEBUGLOG << "ReadTransferR(): SMS" << endl
                 << " indVar=" << trfr.indVar << " PR=" << movt.pr << " PRmethod=" << movt.prMethod << endl;
        DEBUGLOG << "ReadTransferR(): dp=" << movt.dp << " MemSize=" << movt.memSize << " gb=" << movt.gb
                 << " goaltype=" << movt.goalType << endl;
#endif

        smsparams.straightenPath = Rcpp::as<bool>(TransParamsR.slot("StraightenPath")); // Straighten path after decision not to settle?

        // Mortality
        Rcpp::NumericVector HabMortVec;
        HabMortVec = Rcpp::as<Rcpp::NumericVector>(TransParamsR.slot("StepMort"));
        if(HabMortVec.size() == 1) {
            trfr.habMort = false;                 // Per-step mortality type: 0 = constant
            smsparams.stepMort = (float)HabMortVec[0]; // Constant per-step mortality probability
        } else {
            trfr.habMort = true; // Per-step mortality type: 1 = habitat-dependent
            smsparams.stepMort = -9;
            if(paramsLand.generated) {
                // values are for habitat (hab=1) then for matrix (hab=0)
                pSpecies->setHabMort(1, (double)HabMortVec[1]);
                pSpecies->setHabMort(0, (double)HabMortVec[0]);
#if RSDEBUG
                DEBUGLOG << "ReadTransferR(): Generated Landscpae with MortHabitat=" << pSpecies->getHabMort(1)
                         << " MortMatrix=" << pSpecies->getHabMort(0) << endl;
#endif
            }
            if(paramsLand.rasterType == 0) {
#if RSDEBUG
                DEBUGLOG << "ReadTransferR(): nHabMax = " << paramsLand.nHabMax << endl;
#endif
                for(int i = 0; i < paramsLand.nHabMax; i++) {
                    pSpecies->setHabMort(i, (double)HabMortVec[i]);
#if RSDEBUG
                    DEBUGLOG << "ReadTransferR(): Habitat #" << i << ": mortality = " << pSpecies->getHabMort(i)
                             << endl;
#endif
                }
            }
        }


#if RSDEBUG
        DEBUGLOG << "ReadTransferR(): SMtype=" << trfr.habMort << " SMconst=" << movt.stepMort << endl;
#endif

        // Costs
        trfr.costMap = Rcpp::as<bool>(TransParamsR.slot("CostMap"));

        // read habitat costs for land types
        Rcpp::NumericVector HabCostVec;
        if(!trfr.costMap) {
            HabCostVec = Rcpp::as<Rcpp::NumericVector>(TransParamsR.slot("Costs"));
#if RSDEBUG
            DEBUGLOG << "ReadTransferR(): Read Habitat Cost vector of length " << HabCostVec.size() << endl;
#endif
        }

        if(!paramsLand.generated) {          // real landscape
            if(paramsLand.rasterType == 0) { // habitat codes
                if(!trfr.costMap) {          // habitat costs
                    for(int i = 0; i < paramsLand.nHabMax; i++) {
                        pSpecies->setHabCost(i, (int)HabCostVec[i]);
#if RSDEBUG
                        DEBUGLOG << "ReadTransferR(): Habitat #" << i << ": cost = " << pSpecies->getHabCost(i) << endl;
#endif
                    }
                }
            } else { // habitat quality
                // should have trfr.costMap = 1
            }
        } else {               // artificial landscape
            if(trfr.costMap) { // should not occur
                // should have trfr.costMap = 0
            } else { // habitat costs
                // costs are for habitat (hab=1) then for matrix (hab=0)
                pSpecies->setHabCost(1, (int)HabCostVec[1]);
                pSpecies->setHabCost(0, (int)HabCostVec[0]);
            }
        }
        pSpecies->setTrfrRules(trfr);
        pSpecies->setSpMovtTraits(smsparams);
    }
        break; // end of SMS

    case 2: { // CRW

        Rcpp::S4 TransParamsR("CorrRW");
        TransParamsR = Rcpp::as<Rcpp::S4>(DispParamsR.slot("Transfer"));

        // simulation = Rcpp::as<int>(TransParamsR.slot("Simulation"));// REMOVED in R-interface  //Must match
        // simulation numbers in ParamParams
        trfr.indVar = Rcpp::as<bool>(TransParamsR.slot("IndVar"));

        if(trfr.indVar) {
        	mparams.stepLength = -9; // Step length initial mean (m); Required for IndVar = 1
        } else {
            ReadVec = Rcpp::as<Rcpp::NumericVector>(TransParamsR.slot("StepLength")); // Step length params
            if(ReadVec.size() == 1) {
                mparams.stepLength = (float)ReadVec[0]; // Step length (m); Required for IndVar = 0; must be > 0
            } else {
                error = 436;
            }
        }

        if(trfr.indVar) {
        	mparams.rho = -9; // Step correlation coefficient initial mean (m); Required for IndVar = 1
        } else {
            ReadVec = Rcpp::as<Rcpp::NumericVector>(TransParamsR.slot("Rho")); // Step correlation coefficient params
            if(ReadVec.size() == 1) {
                mparams.rho =
                    (float)ReadVec[0]; // Step correlation coefficient; Required for IndVar = 0; must be > 0.0 and < 1.0
            } else {
                error = 436;
            }
        }

        mparams.straightenPath = Rcpp::as<bool>(TransParamsR.slot("StraightenPath")); // Straighten path after decision not to settle?
        // pSpecies->setTrfrScales(scale);

#if RSDEBUG
        DEBUGLOG << "ReadTransferR():"
                 << " paramsLand.rasterType=" << paramsLand.rasterType << " trfr.indVar=" << trfr.indVar
                 << " move.stepLength=" << movt.stepLength << " move.rho=" << movt.rho
                 << " mparams.stepLgthMean=" << mparams.stepLgthMean << " mparams.rhoMean=" << mparams.rhoMean
                 << " move.straigtenPath=" << movt.straigtenPath << endl;
#endif

        // Mortality
        Rcpp::NumericVector HabMortVec;
        HabMortVec = Rcpp::as<Rcpp::NumericVector>(TransParamsR.slot("StepMort"));

        if(HabMortVec.size() == 1) {
            trfr.habMort = 0;                     // Per-step mortality type: 0 = constant
            mparams.stepMort = (float)HabMortVec[0]; // Constant per-step mortality probability
        } else {
            trfr.habMort = 1; // Per-step mortality type: 1 = habitat-dependent
            mparams.stepMort = -9;
            if(paramsLand.generated) {
                // values are for habitat (hab=1) then for matrix (hab=0)
                pSpecies->setHabMort(1, (double)HabMortVec[1]);
                pSpecies->setHabMort(0, (double)HabMortVec[0]);
#if RSDEBUG
                DEBUGLOG << "ReadTransferR(): Generated Landscpae with MortHabitat=" << pSpecies->getHabMort(1)
                         << " MortMatrix=" << pSpecies->getHabMort(0) << endl;
#endif
            }
            if(paramsLand.rasterType == 0) {
#if RSDEBUG
                DEBUGLOG << "ReadTransferR(): nHabMax = " << paramsLand.nHabMax << endl;
#endif
                for(int i = 0; i < paramsLand.nHabMax; i++) {
                    pSpecies->setHabMort(i, (double)HabMortVec[i]);
#if RSDEBUG
                    DEBUGLOG << "ReadTransferR(): Habitat #" << i << ": mortality = " << pSpecies->getHabMort(i)
                             << endl;
#endif
                }
            }
        }
        if(trfr.habMort && paramsLand.rasterType != 0)
            error = 434; // habitat percentage landscape cant have habitat-dependent mortality

#if RSDEBUG
        DEBUGLOG << "ReadTransferR(): SMtype=" << trfr.habMort << " SMconst=" << movt.stepMort << endl;
#endif

        pSpecies->setTrfrRules(trfr);
        pSpecies->setSpMovtTraits(mparams);

    }
        break; // end of CRW

    default:
        error = 440;
    } // end of switch (TransferType)

    if(trfr.indVar) gHasGenetics = true;

    return error;
}

//---------------------------------------------------------------------------
// NOTE that stage- and sex-dependent settlement parameters are set for
// ALL stage/sex combinations, even if the species has stage- and/or
// sex-independent settlement rules
int ReadSettlementR(Rcpp::S4 ParMaster)
{

    Rcpp::S4 DispParamsR("DispersalParams");
    DispParamsR = Rcpp::as<Rcpp::S4>(ParMaster.slot("dispersal"));
    Rcpp::S4 SettleParamsR("SettlementParams");
    SettleParamsR = Rcpp::as<Rcpp::S4>(DispParamsR.slot("Settlement"));

    int Nlines, stage, sex, offset, sexSettle, settType;
    int error = 0;
    bool findmate, densdep;

    demogrParams dem = pSpecies->getDemogrParams();
    stageParams sstruct = pSpecies->getStageParams();
    transferRules trfr = pSpecies->getTransferRules();
    settleType sett = pSpecies->getSettle();
    settleRules srules;
    settleSteps ssteps = pSpecies->getSteps(0,0);
    settleTraits settleDD; // = pSpecies->getSettTraits(0,0);
    sett.stgDep = Rcpp::as<bool>(SettleParamsR.slot("StageDep")); // Stage-dependent settlement.
    sett.sexDep = Rcpp::as<bool>(SettleParamsR.slot("SexDep"));   // Sex-dependent settlement.
    if ( gTransferType == 0) {
        // dispersal kernel                                         // dispersal kernel
        sett.indVar = false;
        if(dem.repType == 0) {
            if(sett.sexDep)
                error = 501; // sex-dependent settlement is not possible with asexual models
        }
        if(!dem.stageStruct) {
            if(sett.stgDep)
                error = 502; // stage-dependent settlement is not possible without stage structure
        }
    } else { // movement process
        sett.indVar = Rcpp::as<bool>(SettleParamsR.slot("IndVar"));
        densdep = Rcpp::as<bool>(SettleParamsR.slot("DensDep"));
        if(dem.repType == 0) {
            if(sett.sexDep)
                error = 508; // sex-dependent settlement is not possible with asexual models
        }
        if(!dem.stageStruct) {
            if(sett.stgDep)
                error = 509; // stage-dependent settlement is not possible without stage structure
        }
    }
    pSpecies->setSettle(sett);

    // no.of lines according to known stage- and sex-dependency and corresponding column offset
    if(sett.stgDep) {
        if(sett.sexDep) {
            Nlines = sstruct.nStages * gNbSexesDisp;
            offset = 2;
        } else {
            Nlines = sstruct.nStages;
            offset = 1;
        }
    } else {
        if(sett.sexDep) {
            Nlines = gNbSexesDisp;
            offset = 1;
        } else {
            Nlines = 1;
            offset = 0;
        }
    }
#if RSDEBUG
    DEBUGLOG << "ReadSettlementR(): sett.stgDep = " << sett.stgDep << ", sett.sexDep = " << sett.sexDep
             << ", sett.indVar = " << sett.indVar << endl;
#endif
    // check if SettleParamsR.slots are single values or a matrix and assign them accordingly:
    Rcpp::NumericMatrix mat(1,1); // temporary storage for single value

    Rcpp::NumericMatrix FindMate;
    bool constFindMate = false;
    if(Rf_isMatrix(SettleParamsR.slot("FindMate"))){
        FindMate = Rcpp::as<Rcpp::NumericMatrix>(SettleParamsR.slot("FindMate"));
    } else {
        mat(0,0) = Rcpp::as<int>(SettleParamsR.slot("FindMate"));
        FindMate = mat;
        constFindMate = true;
    }

    Rcpp::NumericMatrix MinSteps;
    bool constMinSteps = false;
    if(Rf_isMatrix(SettleParamsR.slot("MinSteps"))){
        MinSteps = Rcpp::as<Rcpp::NumericMatrix>(SettleParamsR.slot("MinSteps"));
    } else {
        mat(0,0) = Rcpp::as<int>(SettleParamsR.slot("MinSteps"));
        MinSteps = mat;
        constMinSteps = true;
    }


    Rcpp::NumericMatrix MaxSteps;
    bool constMaxSteps = false;
    if(Rf_isMatrix(SettleParamsR.slot("MaxSteps"))){
        MaxSteps = Rcpp::as<Rcpp::NumericMatrix>(SettleParamsR.slot("MaxSteps"));
    } else {
        mat(0,0) = Rcpp::as<int>(SettleParamsR.slot("MaxSteps"));
        MaxSteps = mat;
        constMaxSteps = true;
    }

    Rcpp::NumericMatrix MaxStepsYr;
    bool constMaxStepsYr = false;
    if(Rf_isMatrix(SettleParamsR.slot("MaxStepsYear"))){
        MaxStepsYr = Rcpp::as<Rcpp::NumericMatrix>(SettleParamsR.slot("MaxStepsYear"));
    } else {
        mat(0,0) = Rcpp::as<int>(SettleParamsR.slot("MaxStepsYear"));
        MaxStepsYr = mat;
        constMaxStepsYr = true;
    }

    sexSettle = 2 * sett.stgDep + sett.sexDep;

    Rcpp::NumericMatrix SettleCondMatrix;

    if(Rf_isMatrix(SettleParamsR.slot("Settle"))){
        SettleCondMatrix = Rcpp::as<Rcpp::NumericMatrix>(SettleParamsR.slot("Settle"));
    } else {
        mat(0,0) = Rcpp::as<int>(SettleParamsR.slot("Settle"));
        SettleCondMatrix = mat;
    }

    for(int line = 0; line < Nlines; line++) {
        // FindMate
        // determine stage and sex of this line
        if(sett.stgDep) {
            if(sett.sexDep) {
                stage = (int)FindMate(line, 0);
                sex = (int)FindMate(line, 1);
            } else {
                stage = (int)FindMate(line, 0);
                sex = 0;
            }
        } else {
            if(sett.sexDep) {
                stage = 0;
                sex = (int)FindMate(line, 0);
            } else {
                stage = 0;
                sex = 0;
            }
        }

        if(trfr.usesMovtProc) { // ...movement process
            if(constFindMate) {
                findmate = (bool)FindMate(0, 0);
            } else {
                findmate = (bool)FindMate(line, offset);
            }
            if(findmate && dem.repType == 0)
                error = 504;

            switch(sexSettle) {

            case 0: { // no sex- / stage-dependence
                    srules = pSpecies->getSettRules(0, 0);
                    srules.densDep = densdep;
                    srules.findMate = findmate;
                    pSpecies->setSettRules(0, 0, srules);

                    if(dem.stageStruct) { // model is structured - also set parameters for all stages
                        for(int i = 1; i < sstruct.nStages; i++) {
                            pSpecies->setSettRules(i, 0, srules);
                            if(dem.repType > 0) {                        // model is sexual - also set parameters for males
                                pSpecies->setSettRules(i, 1, srules);
                            }
                        }
                    } else {                  // see comment above (at case label)
                        if(dem.repType > 0) { // model is sexual - also set parameters for males
                            pSpecies->setSettRules(0, 1, srules);
                        }
                    }
                }
                break;

            case 1: { // sex-dependent
                srules = pSpecies->getSettRules(0, sex);
                srules.densDep = densdep;
                srules.findMate = findmate;
                pSpecies->setSettRules(0, sex, srules);
                if(dem.stageStruct) { // model is structured - also set parameters for all stages
                    for(int i = 1; i < sstruct.nStages; i++) {
                        pSpecies->setSettRules(i, sex, srules);
                    }
                }
            }
                break;

            case 2: { // stage-dependent
                srules = pSpecies->getSettRules(stage, 0);
                srules.densDep = densdep;
                srules.findMate = findmate;
                pSpecies->setSettRules(stage, 0, srules);
                if(dem.repType > 0) { // model is sexual - also set parameters for males
                    pSpecies->setSettRules(stage, 1, srules);
                }
            }
                break;

            case 3: { // sex- & stage-dependent
                srules = pSpecies->getSettRules(stage, sex);
                srules.densDep = densdep;
                srules.findMate = findmate;
                pSpecies->setSettRules(stage, sex, srules);
            }
                break;
            } // end sexSettle for FindMate

        }
        // read find mate conditions for...
        else { // ...dispersal kernel
            if(constFindMate) {
                findmate = (bool)FindMate(0, 0);
            } // Mating requirements to settle, required for a sexual population only
            else {
                findmate = (bool)FindMate(line, offset);
            }
            if(findmate && dem.repType == 0)
                error = 504;

            switch(sexSettle) {
            case 0: { // no sex / stage dependence
                    if(findmate && dem.repType == 0)
                        error = 504;
                    srules = pSpecies->getSettRules(0, 0);
                    srules.findMate = findmate;
                    if(dem.stageStruct) { // model is structured - also set parameters for all stages
                        for(int i = 0; i < sstruct.nStages; i++) {
                            pSpecies->setSettRules(i, 0, srules);
                            if(dem.repType > 0) { // model is sexual - also set parameters for males
                                pSpecies->setSettRules(i, 1, srules);
                            }
                        }
                    } else {
                        pSpecies->setSettRules(0, 0, srules);
                        if(dem.repType > 0) { // model is sexual - also set parameters for males
                            pSpecies->setSettRules(0, 1, srules);
                        }
                    }
                }
                break;

            case 1: { // sex dependent
                srules = pSpecies->getSettRules(0, sex);
                srules.findMate = findmate;
                pSpecies->setSettRules(0, sex, srules);
                if(dem.stageStruct) { // model is structured - also set parameters for all stages
                    for(int i = 1; i < sstruct.nStages; i++) {
                        pSpecies->setSettRules(i, sex, srules);
                    }
                }
            }
                break;

            case 2: { // stage dependent
                if(findmate && dem.repType == 0)
                    error = 507;
                srules = pSpecies->getSettRules(stage, 0);
                srules.findMate = findmate;
                pSpecies->setSettRules(stage, 0, srules);
                if(dem.repType > 0) { // model is sexual - also set parameters for males
                    pSpecies->setSettRules(stage, 1, srules);
                }
            }
                break;

            case 3: { // sex & stage dependent
                srules = pSpecies->getSettRules(stage, sex);
                srules.findMate = findmate;
                pSpecies->setSettRules(stage, sex, srules);
            }
                break;

            } // end of switch (sexSettle)

        } // end of dispersal kernel

        // MinSteps
        // determine stage and sex of this line
        if(sett.stgDep) {
            if(sett.sexDep) {
                stage = (int)MinSteps(line, 0);
                sex = (int)MinSteps(line, 1);
            } else {
                stage = (int)MinSteps(line, 0);
                sex = 0;
            }
        } else {
            if(sett.sexDep) {
                stage = 0;
                sex = (int)MinSteps(line, 0);
            } else {
                stage = 0;
                sex = 0;
            }
        }

        if(trfr.usesMovtProc) { // ...movement process
            if(constMinSteps) {
                ssteps.minSteps = (int)MinSteps(0, 0);
            } else {
                ssteps.minSteps = (int)MinSteps(line, offset);
            }

            switch(sexSettle) {

            case 0: { // no sex- / stage-dependence
                srules = pSpecies->getSettRules(0, 0);
                pSpecies->setSteps(0, 0, ssteps);

                if(dem.stageStruct) { // model is structured - also set parameters for all stages
                    for(int i = 1; i < sstruct.nStages; i++) {
                        pSpecies->setSteps(i, 0, ssteps);
                        if(dem.repType > 0) {                        // model is sexual - also set parameters for males
                            pSpecies->setSteps(i, 1, ssteps);
                        }
                    }
                } else {                  // see comment above (at case label)
                    if(dem.repType > 0) { // model is sexual - also set parameters for males
                        pSpecies->setSteps(0, 1, ssteps);
                    }
                }
            }
                break;

            case 1: { // sex-dependent
                srules = pSpecies->getSettRules(0, sex);
                pSpecies->setSteps(0, sex, ssteps);
                if(dem.stageStruct) { // model is structured - also set parameters for all stages
                    for(int i = 1; i < sstruct.nStages; i++) {
                        pSpecies->setSteps(i, sex, ssteps);
                    }
                }


            }
                break;

            case 2: { // stage-dependent
                srules = pSpecies->getSettRules(stage, 0);
                pSpecies->setSteps(stage, 0, ssteps);
                if(dem.repType > 0) { // model is sexual - also set parameters for males
                    pSpecies->setSteps(stage, 1, ssteps);
                }
            }
                break;

            case 3: { // sex- & stage-dependent
                srules = pSpecies->getSettRules(stage, sex);
                pSpecies->setSteps(stage, sex, ssteps);
            }
                break;
            } // end sexSettle
        } // end if MovementModel for MinSteps

        // MaxSteps
        // determine stage and sex of this line
        if(sett.stgDep) {
            if(sett.sexDep) {
                stage = (int)MaxSteps(line, 0);
                sex = (int)MaxSteps(line, 1);
            } else {
                stage = (int)MaxSteps(line, 0);
                sex = 0;
            }
        } else {
            if(sett.sexDep) {
                stage = 0;
                sex = (int)MaxSteps(line, 0);
            } else {
                stage = 0;
                sex = 0;
            }
        }

        if(trfr.usesMovtProc) { // ...movement process
            if(constMaxSteps) {
                ssteps.maxSteps = (int)MaxSteps(0, 0);
            } else {
                ssteps.maxSteps = (int)MaxSteps(line, offset);
            }

            switch(sexSettle) {

            case 0: { // no sex- / stage-dependence
                srules = pSpecies->getSettRules(0, 0);
                pSpecies->setSteps(0, 0, ssteps);

                if(dem.stageStruct) { // model is structured - also set parameters for all stages
                    for(int i = 1; i < sstruct.nStages; i++) {
                        pSpecies->setSteps(i, 0, ssteps);
                        if(dem.repType > 0) {                        // model is sexual - also set parameters for males
                            pSpecies->setSteps(i, 1, ssteps);
                        }
                    }
                } else {                  // see comment above (at case label)
                    if(dem.repType > 0) { // model is sexual - also set parameters for males
                        pSpecies->setSteps(0, 1, ssteps);
                    }
                }
            }
                break;

            case 1: { // sex-dependent
                srules = pSpecies->getSettRules(0, sex);
                pSpecies->setSteps(0, sex, ssteps);
                if(dem.stageStruct) { // model is structured - also set parameters for all stages
                    for(int i = 1; i < sstruct.nStages; i++) {
                        pSpecies->setSteps(i, sex, ssteps);
                    }
                }


            }
                break;

            case 2: { // stage-dependent
                srules = pSpecies->getSettRules(stage, 0);
                pSpecies->setSteps(stage, 0, ssteps);
                if(dem.repType > 0) { // model is sexual - also set parameters for males
                    pSpecies->setSteps(stage, 1, ssteps);
                }
            }
                break;

            case 3: { // sex- & stage-dependent
                srules = pSpecies->getSettRules(stage, sex);
                pSpecies->setSteps(stage, sex, ssteps);
            }
                break;
            } // end sexSettle

        } // End Movement model for MaxSteps

        // MaxStepsYr
        // determine stage and sex of this line
        if(sett.stgDep) {
            if(sett.sexDep) {
                stage = (int)MaxStepsYr(line, 0);
                sex = (int)MaxStepsYr(line, 1);
            } else {
                stage = (int)MaxStepsYr(line, 0);
                sex = 0;
            }
        } else {
            if(sett.sexDep) {
                stage = 0;
                sex = (int)MaxStepsYr(line, 0);
            } else {
                stage = 0;
                sex = 0;
            }
        }

        if(trfr.usesMovtProc) { // ...movement process
            if(constMaxStepsYr) {
                ssteps.maxStepsYr = (int)MaxStepsYr(0, 0);
            } else {
                ssteps.maxStepsYr = (int)MaxStepsYr(line, offset);
            }

            switch(sexSettle) {

            case 0: { // no sex- / stage-dependence
                srules = pSpecies->getSettRules(0, 0);
                pSpecies->setSteps(0, 0, ssteps);

                if(dem.stageStruct) { // model is structured - also set parameters for all stages
                    for(int i = 1; i < sstruct.nStages; i++) {
                        pSpecies->setSteps(i, 0, ssteps);
                        if(dem.repType > 0) {                        // model is sexual - also set parameters for males
                            pSpecies->setSteps(i, 1, ssteps);
                        }
                    }
                } else {                  // see comment above (at case label)
                    if(dem.repType > 0) { // model is sexual - also set parameters for males
                        pSpecies->setSteps(0, 1, ssteps);
                    }
                }
            }
                break;

            case 1: { // sex-dependent
                srules = pSpecies->getSettRules(0, sex);
                pSpecies->setSteps(0, sex, ssteps);
                if(dem.stageStruct) { // model is structured - also set parameters for all stages
                    for(int i = 1; i < sstruct.nStages; i++) {
                        pSpecies->setSteps(i, sex, ssteps);
                    }
                }
            }
                break;

            case 2: { // stage-dependent
                srules = pSpecies->getSettRules(stage, 0);
                pSpecies->setSteps(stage, 0, ssteps);
                if(dem.repType > 0) { // model is sexual - also set parameters for males
                    pSpecies->setSteps(stage, 1, ssteps);
                }
            }
                break;

            case 3: { // sex- & stage-dependent
                srules = pSpecies->getSettRules(stage, sex);
                pSpecies->setSteps(stage, sex, ssteps);
            }
                break;
            } // end sexSettle

        } // End Movement model

        // Settle
        // determine stage and sex of this line
        if(!sett.indVar){
            if(sett.stgDep) {
                if(sett.sexDep) {
                    stage = (int)SettleCondMatrix(line, 0);
                    sex = (int)SettleCondMatrix(line, 1);
                } else {
                    stage = (int)SettleCondMatrix(line, 0);
                    sex = 0;
                }
            } else {
                if(sett.sexDep) {
                    stage = 0;
                    sex = (int)SettleCondMatrix(line, 0);
                } else {
                    stage = 0;
                    sex = 0;
                }
            }
        }

        if(trfr.usesMovtProc) { // ...movement process
            if(densdep) {
                if (sett.indVar) {
                    gHasGenetics = true;
                    settleDD.s0 = -9;
                    settleDD.alpha = -9;
                    settleDD.beta = -9;
                } else{
                    settleDD.s0 = (float)SettleCondMatrix(
                        line, offset + 0); // Max. settlement probability for density reaction norm. Required for
                    // DensDep = 1 and IndVar = 0; 0.0 < S0 <= 1.0
                    settleDD.alpha =
                        (float)SettleCondMatrix(line, offset + 1); // Required for DensDep = 1 and IndVar = 0
                    settleDD.beta =
                        (float)SettleCondMatrix(line, offset + 2); // Required for DensDep = 1 and IndVar = 0
                }
            }

            switch(sexSettle) {

            case 0: { // no sex- / stage-dependence
                srules = pSpecies->getSettRules(0, 0);
                // Loops vereinfachen indem ich direkt zu Anfang die densDep Bedingung überprüfe?
                if(srules.densDep) {
                    pSpecies->setSpSettTraits(0, 0, settleDD);
                    if(dem.stageStruct) { // model is structured - also set parameters for all stages
                        for(int i = 1; i < sstruct.nStages; i++) {
                            pSpecies->setSpSettTraits(i, 0, settleDD); //  /!\ different to ReadSettlement()
                            if(dem.repType > 0) {                        // model is sexual - also set parameters for males
                                pSpecies->setSpSettTraits(i, 1, settleDD);
                            }
                        }
                    } else {                  // see comment above (at case label)
                        if(dem.repType > 0) { // model is sexual - also set parameters for males
                            pSpecies->setSpSettTraits(0, 1, settleDD);
                        }
                    }
                }
            }
                break;

            case 1: { // sex-dependent

                if(sett.indVar){
                    for (int s = 0; s < gNbSexesDisp; s++){
                        pSpecies->setSpSettTraits(0, s, settleDD);
                        if(dem.stageStruct) { // model is structured - also set parameters for all stages
                            for(int i = 1; i < sstruct.nStages; i++) {
                                pSpecies->setSpSettTraits(i, s, settleDD);
                            }
                        }
                    }
                } else {
                    srules = pSpecies->getSettRules(0, sex);
                    // s. Kommentar oben
                    if(srules.densDep) {
                        pSpecies->setSpSettTraits(0, sex, settleDD);
                        if(dem.stageStruct) { // model is structured - also set parameters for all stages
                            for(int i = 1; i < sstruct.nStages; i++) {
                                pSpecies->setSpSettTraits(i, sex, settleDD);
                            }
                        }
                    }
                }

            }
                break;

            case 2: { // stage-dependent (cannot include individual variablility, so no need to account for that)
                srules = pSpecies->getSettRules(stage, 0);
                if(srules.densDep) {
                    pSpecies->setSpSettTraits(stage, 0, settleDD);
                    if(dem.repType > 0) { // model is sexual - also set parameters for males
                        pSpecies->setSpSettTraits(stage, 1, settleDD);
                    }
                }
            }
                break;

            case 3: { // sex- & stage-dependent (cannot include individual variablility, so no need to account for that)
                srules = pSpecies->getSettRules(stage, sex);
                if(srules.densDep) {
                    pSpecies->setSpSettTraits(stage, sex, settleDD);
                }
            }
                break;
            } // end sexSettle

        } // end of movement model for SettleCondMatrix

        // read settlement conditions for...
        else { // ...dispersal kernel

            settType = (int)SettleCondMatrix(line,
                        offset); // Settlement rule if the arrival cell/patch is unsuitable: 0 = die, 1 = wait, 2 = randomly
            // choose a suitable cell/patch or die, 3 = randomly choose a suitable cell/patch or wait.
            // Options 1 and 3 may be chosen for a stage-structured population only

            switch(sexSettle) {
            case 0: { // no sex / stage dependence
                if((settType == 1 || settType == 3) && !dem.stageStruct)
                    error = 503;
                srules = pSpecies->getSettRules(0, 0);
                switch(settType) {
                case 0:
                    srules.wait = false;
                    srules.go2nbrLocn = false;
                    break;
                case 1:
                    srules.wait = true;
                    srules.go2nbrLocn = false;
                    break;
                case 2:
                    srules.wait = false;
                    srules.go2nbrLocn = true;
                    break;
                case 3:
                    srules.wait = true;
                    srules.go2nbrLocn = true;
                    break;
                }
                if(dem.stageStruct) { // model is structured - also set parameters for all stages
                    for(int i = 0; i < sstruct.nStages; i++) {
                        pSpecies->setSettRules(i, 0, srules);
                        if(dem.repType > 0) { // model is sexual - also set parameters for males
                            pSpecies->setSettRules(i, 1, srules);
                        }
                    }
                } else {
                    pSpecies->setSettRules(0, 0, srules);
                    if(dem.repType > 0) { // model is sexual - also set parameters for males
                        pSpecies->setSettRules(0, 1, srules);
                    }
                }
            }
                break;

            case 1: { // sex dependent
                if((settType == 1 || settType == 3) && dem.stageStruct == false)
                    error = 505;
                srules = pSpecies->getSettRules(0, sex);
                switch(settType) {
                case 0:
                    srules.wait = false;
                    srules.go2nbrLocn = false;
                    break;
                case 1:
                    srules.wait = true;
                    srules.go2nbrLocn = false;
                    break;
                case 2:
                    srules.wait = false;
                    srules.go2nbrLocn = true;
                    break;
                case 3:
                    srules.wait = true;
                    srules.go2nbrLocn = true;
                    break;
                }
                pSpecies->setSettRules(0, sex, srules);
                if(dem.stageStruct) { // model is structured - also set parameters for all stages
                    for(int i = 1; i < sstruct.nStages; i++) {
                        pSpecies->setSettRules(i, sex, srules);
                    }
                }
            }
                break;

            case 2: { // stage dependent
                srules = pSpecies->getSettRules(stage, 0);
                switch(settType) {
                case 0:
                    srules.wait = false;
                    srules.go2nbrLocn = false;
                    break;
                case 1:
                    srules.wait = true;
                    srules.go2nbrLocn = false;
                    break;
                case 2:
                    srules.wait = false;
                    srules.go2nbrLocn = true;
                    break;
                case 3:
                    srules.wait = true;
                    srules.go2nbrLocn = true;
                    break;
                }
                srules.findMate = findmate;
                pSpecies->setSettRules(stage, 0, srules);
                if(dem.repType > 0) { // model is sexual - also set parameters for males
                    pSpecies->setSettRules(stage, 1, srules);
                }
            }
                break;

            case 3: { // sex & stage dependent
                srules = pSpecies->getSettRules(stage, sex);
                switch(settType) {
                case 0:
                    srules.wait = false;
                    srules.go2nbrLocn = false;
                    break;
                case 1:
                    srules.wait = true;
                    srules.go2nbrLocn = false;
                    break;
                case 2:
                    srules.wait = false;
                    srules.go2nbrLocn = true;
                    break;
                case 3:
                    srules.wait = true;
                    srules.go2nbrLocn = true;
                    break;
                }
                pSpecies->setSettRules(stage, sex, srules);
            }
                break;

            } // end of switch (sexSettle)

        } // end of dispersal kernel

    } // end of for line loop

    return error;
}


//---------------------------------------------------------------------------

int ReadInitialisationR(Landscape* pLandscape, Rcpp::S4 ParMaster)
{

    Rcpp::S4 InitParamsR("InitialisationParams");
    InitParamsR = Rcpp::as<Rcpp::S4>(ParMaster.slot("init"));

    Rcpp::NumericVector PropStages;

    landParams paramsLand = pLandscape->getLandParams();
    demogrParams dem = pSpecies->getDemogrParams();
    stageParams sstruct = pSpecies->getStageParams();
    initParams init = paramsInit->getInit();
    string Inputs = paramsSim->getDir(1);

    int maxcells; //,simulation
    float p, check;
    int error = 0;

    // simulation = Rcpp::as<int>(InitParamsR.slot("Simulation"));// REMOVED in R-interface
    init.seedType = Rcpp::as<int>(
        InitParamsR.slot("InitType")); // Initialisation type: 0 = free initialisation, 1 = from species distribution, 2
    // = from initial individuals file. Must be 0 for artificial landscapes.
    init.freeType = Rcpp::as<int>(InitParamsR.slot(
        "FreeType")); // If SeedType = 0: Type of free initialisation (i.e. not from specified distribution): 0 = Random
    // (given no. of cells/patches), 1 = All suitable cells/patches.
    init.spDistType = Rcpp::as<int>(
        InitParamsR.slot("SpType")); // If SeedType = 1: Type of initialisation from species distribution: 0 = All
    // suitable cells/patches within the distribution cells, 1 = Some randomly chosen
    // suitable cells/patches within distribution cells
    if(init.seedType == 1 && paramsLand.spDist == false)
        error = 601;

    if(paramsLand.patchModel) {
        init.initDens = Rcpp::as<int>(InitParamsR.slot(
            "InitDens")); // Initial density of individuals in each patch: 0 = at carrying capacity; 1
        // = at half carrying capacity; 2 = specified density. Patch carrying capacity
        // is determined by the carrying capacities of the cells forming the patch.
        init.indsHa =
            Rcpp::as<float>(InitParamsR.slot("IndsHaCell")); // If initDens = 2: Density of individuals to initialise
        // per hectare; required for initDens = 2.
    } else {
        init.initDens = Rcpp::as<int>(
            InitParamsR.slot("InitDens")); // Initial no. of individuals in each cell: 0 = at carrying capacity; 1 = at
        // half carrying capacity; 2 = specified no. of individuals
        init.indsCell =
            Rcpp::as<int>(InitParamsR.slot("IndsHaCell")); // If initDens = 2: No. of individuals to initialise in each
        // cell; required for initDens = 2, otherwise set to -9
    }

    init.minSeedX =
        Rcpp::as<int>(InitParamsR.slot("minX")); // Required for SeedType = 0 only; minX and minY may not be less than
    // 0, maxX and maxY may not be less than corresponding minimum
    init.maxSeedX = Rcpp::as<int>(InitParamsR.slot("maxX")); // -"-
    init.minSeedY = Rcpp::as<int>(InitParamsR.slot("minY")); // -"-
    init.maxSeedY = Rcpp::as<int>(InitParamsR.slot("maxY")); // -"-
    init.nSeedPatches = Rcpp::as<int>(
        InitParamsR.slot("NrCells")); // If SeedType = 0 and FreeType = 0: No. of cells/patches to initialise
    init.nSpDistPatches = Rcpp::as<int>(InitParamsR.slot(
        "NrCells")); // If  SeedType = 1 and SpType = 1: No. of species dist cells/patches to initialise
    init.initFrzYr =
        Rcpp::as<int>(InitParamsR.slot("InitFreezeYear")); // Year until which species is confined to initial range
    // limits. Must be >= 0 for SeedType = 0.
    init.restrictRows =
        Rcpp::as<int>(InitParamsR.slot("RestrictRows")); // No. of rows at northern front to restrict range. Must be >0
    // if applied for SeedType = 0, otherwise 0.
    init.restrictFreq =
        Rcpp::as<int>(InitParamsR.slot("RestrictFreq")); // Frequency (years) at which range is restricted to northern
    // front. Must be > 0 if RestrictRows is > 0.
    init.finalFrzYr =
        Rcpp::as<int>(InitParamsR.slot("FinalFreezeYear")); // Year after which species is confined to current range
    // limits. Must be > InitFreezeYear if applied, otherwise 0.
    init.indsFile = Rcpp::as<string>(InitParamsR.slot("InitIndsFile")); // Name of the initial individuals file (*.txt). Required for SeedType = 2, otherwise NULL.
#if RSDEBUG
    DEBUGLOG << "ReadInitialisationR():"
    //<< " simulation=" << simulation
      << " seedType=" << init.seedType
      << " freeType=" << init.freeType << " spDistType=" << init.spDistType << " maxSeedX=" << init.maxSeedX
      << " maxSeedY=" << init.maxSeedY << " initFrzYr=" << init.initFrzYr
      << " restrictRows=" << init.restrictRows << " restrictFreq=" << init.restrictFreq
      << " finalFrzYr=" << init.finalFrzYr << " indsFile=" << init.indsFile << endl;
#endif
    init.restrictRange = false;
    if(init.seedType == 0 && init.restrictRows > 0)
        init.restrictRange = true;

    if(dem.stageStruct) {
        init.initAge = Rcpp::as<int>(InitParamsR.slot("InitAge")); // Initial age distribution within each stage: 0 = lowest possible age, 1 = randomised, 2 =
        // quasi-equilibrium. Required for StageStruct = 1 only - otherwise OMIT COLUMNS.
        PropStages = Rcpp::as<Rcpp::NumericVector>(InitParamsR.slot("PropStages")); // Proportion of the initial individuals in stage class i>0 (No juveniles
        // are initialized). Required for StageStruct = 1 only (number of columns
        // is one less than number of stages) - otherwise OMIT COLUMNS
        if(init.seedType != 2) {
            check = 0.0;
            for(int i = 1; i < sstruct.nStages; i++) {
                p = (float)PropStages(i);
                check += p;
                paramsInit->setProp(i, p);
            }
            if(check < 1.0 || check > 1.0) {
                // this condition should not occur - WHAT COULD BE DONE?? ABORT WITH ERROR CODE ...
#if RSDEBUG
                DEBUGLOG << "ReadInitialisation(): check = " << check << endl;
#endif
                throw logic_error("The proportion of initial individuals in each stage doesn not sum to 1.");
            }
        }
    }
    paramsInit->setInit(init);
    switch(init.seedType) {
    case 0: { // free initialisation
        if(init.minSeedX < 0)
            init.minSeedX = 0;
        if(init.minSeedY < 0)
            init.minSeedY = 0;
        if(init.maxSeedX < 0 || init.maxSeedX > paramsLand.maxX)
            init.maxSeedX = paramsLand.maxX; // added upper boundary checks
        if(init.maxSeedY < 0 || init.maxSeedY > paramsLand.maxY)
            init.maxSeedY = paramsLand.maxY; // added upper boundary checks
        if(init.minSeedY > init.maxSeedY || init.minSeedX > init.maxSeedX) {
#if RSDEBUG
            DEBUGLOG << "ReadInitialisationR(): maxSeedX=" << init.maxSeedX << " paramsLand.maxX=" << paramsLand.maxX
                     << " maxSeedY=" << init.maxSeedY << " paramsLand.maxY=" << paramsLand.maxY << endl;
#endif
            error = 603;
        }
        maxcells = (init.maxSeedY - init.minSeedY) * (init.maxSeedX - init.minSeedX);
        if(init.freeType == 0 && init.nSeedPatches > maxcells)
            error = 602;
    }
        break;
    case 1: // from species distribution
        break;
    case 2: { // from initial individuals file
        // if (init.indsFile != prevInitialIndsFile) {
        // read and store the list of individuals to be initialised
        error = ReadInitIndsFileR(0, pLandscape); //open, parse, read header and lines, store in vector "initinds"
        // prevInitialIndsFile = init.indsFile;
        //}
    }
        break;
    default:
        ;
    }
    return error;
}

//---------------------------------------------------------------------------

int ReadGeneticsR(Rcpp::S4 GeneParamsR, Landscape* pLandscape)
{
    bool outputGeneValues, outputWeirCockerham, outputWeirHill; // what should be calculated in the output?
    int outputStartGenetics, outputGeneticInterval; // when should the output be calculated?
    set<int>patchList; // for which patches should the output be calculated?

    //not ideal to reset these in here
    pSpecies->resetGeneticParameters();

    int genomeSize = Rcpp::as<int>(GeneParamsR.slot("GenomeSize")); // how many loci are there?

    set<int> chrEnds;

    if(GeneParamsR.slot("ChromosomeEnds") != R_NilValue){
        Rcpp::IntegerVector ChromosomeEnds = Rcpp::as<Rcpp::IntegerVector>(GeneParamsR.slot("ChromosomeEnds")); // where do the chromosomes end?
        for (int i = 0; i < ChromosomeEnds.size(); i++) {
            chrEnds.insert(ChromosomeEnds[i]);
        }
    }
    else {
        chrEnds.insert(genomeSize - 1);
    }


    float recombinationRate = Rcpp::as<float>(GeneParamsR.slot("RecombinationRate"));

    outputGeneValues = Rcpp::as<bool>(GeneParamsR.slot("OutputGeneValues"));
    outputWeirCockerham = Rcpp::as<bool>(GeneParamsR.slot("OutputFstatsWeirCockerham"));
    outputWeirHill = Rcpp::as<bool>(GeneParamsR.slot("OutputFstatsWeirHill"));


    if(GeneParamsR.slot("OutputStartGenetics") != R_NilValue){
        outputStartGenetics = Rcpp::as<int>(GeneParamsR.slot("OutputStartGenetics"));
    } else {
        outputStartGenetics = -9;
    }
    if(GeneParamsR.slot("OutputInterval") != R_NilValue){
        outputGeneticInterval = Rcpp::as<int>(GeneParamsR.slot("OutputInterval"));
    } else {
        outputGeneticInterval = -9;
    }

    Rcpp::StringVector inPatches;
    string patchSamplingOption;
    int nPatchesToSample = 0;

    if(GeneParamsR.slot("PatchList") != R_NilValue){

        inPatches = Rcpp::as<Rcpp::StringVector>(GeneParamsR.slot("PatchList"));

        if (inPatches[0] != "all" && inPatches[0] != "random" && inPatches[0] != "random_occupied") {
            // then must be a list of indices
            patchSamplingOption = "list";
            for (int i = 0; i < inPatches.size(); i++) {
                patchList.insert(Rcpp::as<int>(inPatches[i]));
            }
            if (patchList.contains(0)) throw logic_error("Patch sampling: ID 0 is reserved for the matrix and should not be sampled.");
        }
        else  {
            patchSamplingOption = inPatches[0];
            if (inPatches[0] == "random" || inPatches[0] == "random_occupied"){
                if(GeneParamsR.slot("NbrPatchToSample") != R_NilValue)
                    nPatchesToSample = Rcpp::as<int>(GeneParamsR.slot("NbrPatchToSample"));
                else throw logic_error("You must provide the number of patches to sample if PatchList is random or random_occupied.");
                // patchList remains empty, filled when patches are sampled every gen
            }

        }
    }

    string NbInds;
    if (GeneParamsR.slot("nIndividualsToSample") != R_NilValue){
        if (Rf_isInteger(GeneParamsR.slot("nIndividualsToSample"))){
            int NbIndsInt = Rcpp::as<int>(GeneParamsR.slot("nIndividualsToSample"));
            NbInds = to_string(NbIndsInt);
        }
        else {
            NbInds = Rcpp::as<string>(GeneParamsR.slot("nIndividualsToSample"));
        }
    }
    const string strNbInds = NbInds;

    const int nbStages = pSpecies->getStageParams().nStages;
    set<int> stagesToSampleFrom;
    Rcpp::StringVector Stages;

    if (GeneParamsR.slot("Stages") != R_NilValue) {
        Stages = Rcpp::as<Rcpp::StringVector>(GeneParamsR.slot("Stages"));
        if (Stages[0] == "all") {
            for (int i = 0; i < nbStages; i++) {
                stagesToSampleFrom.insert(i);
            }
        }
        else {
            for (int i = 0; i < Stages.size(); i++) {
                stagesToSampleFrom.insert(Rcpp::as<int>(Stages[i]));
            }
        }
    }


    pSpecies->setGeneticParameters(chrEnds, genomeSize, recombinationRate,
                                   patchList, strNbInds, stagesToSampleFrom, nPatchesToSample);

    paramsSim->setGeneticSim(patchSamplingOption, outputGeneValues, outputWeirCockerham, outputWeirHill, outputStartGenetics, outputGeneticInterval);
    return 0;
}

//---------------------------------------------------------------------------

int ReadTraitsR(Rcpp::S4 TraitsParamsR)
{
    Rcpp::S4 GeneticLoadParamsR("GeneticLoadParams");
    Rcpp::S4 EmigrationTraitsParamsR("EmigrationTraitsParams");
    Rcpp::S4 KernelTraitsParamsR("KernelTraitsParams");
    Rcpp::S4 CorrRWTraitsParamsR("CorrRWTraitsParams");
    Rcpp::S4 SMSTraitsParamsR("SMSTraitsParams");
    Rcpp::S4 SettlementTraitsParamsR("SettlementTraitsParams");

    settleType sett = pSpecies->getSettle();
    emigRules emig = pSpecies->getEmigRules();
    transferRules trfr = pSpecies->getTransferRules();
    int TransferType = trfr.usesMovtProc ? trfr.moveType : 0;

    pSpecies->clearTraitTable();
    const int genomeSize = pSpecies->getGenomeSize();

    string TraitTypeStr;

    if(gHasNeutralGenetics){
        Rcpp::S4 NeutralTraitsParamsR("NeutralTraitsParams");
        NeutralTraitsParamsR = Rcpp::as<Rcpp::S4>(TraitsParamsR.slot("Neutral"));

        // TraitType
        string TraitTypeR = "neutral";// NEUTRAL; // each column corresponds to a TraitType

        // Positions and number of Positions
        set<int> positions;
        int NbOfPositionsR = -9;

        if(Rf_isString(NeutralTraitsParamsR.slot("Positions"))){
            string PositionsR = Rcpp::as<string>(NeutralTraitsParamsR.slot("Positions"));
                if(PositionsR == "random"){
                    NbOfPositionsR = Rcpp::as<int>(NeutralTraitsParamsR.slot("NbOfPositions"));
                    if(NbOfPositionsR > 0){
                        positions = selectRandomLociPositions(NbOfPositionsR, genomeSize);
                    }
                }
                else throw logic_error("NeutralTraits(): If positions are random you must provide the number of positions (>0).");
        }
        else {
            Rcpp::NumericVector PositionsR = Rcpp::as<Rcpp::NumericVector>(NeutralTraitsParamsR.slot("Positions")); // Positions are provided as numeric vector
            // use a for loop to insert the values at PositionsR into the set positions
            for (int i = 0; i < PositionsR.size(); i++) {
                if (static_cast<int>(PositionsR[i]) <= genomeSize) {
                    positions.insert(static_cast<int>(PositionsR[i]));
                }
                else throw logic_error("NeutralTraits(): Loci positions must be smaller than genome size");
            }
        }

        // Expression type
        string ExpressionTypeR = "#";

        // Initial distribution parameters
        string initDistR = Rcpp::as<string>(NeutralTraitsParamsR.slot("InitialDistribution"));
        if(initDistR != "uniform") initDistR == "#";

        Rcpp::NumericVector initParamsR = {0,Rcpp::as<int>(NeutralTraitsParamsR.slot("InitialParameters"))};

        // Initial dominance distribution parameters not applicable for neutral traits
        string initDomDistR = "#";
        Rcpp::NumericVector initDomParamsR = {0,0};

        // Dominance distribution parameters not applicable for neutral traits
        string DominanceDistR = "#";
        Rcpp::NumericVector DominanceParamsR = {0,0};

        // Mutation parameters
        bool isInherited = true;
        string MutationDistR = Rcpp::as<string>(NeutralTraitsParamsR.slot("MutationDistribution"));
        Rcpp::NumericVector MutationParamsR = Rcpp::as<Rcpp::NumericVector>(NeutralTraitsParamsR.slot("MutationParameters"));
        float MutationRateR =  Rcpp::as<float>(NeutralTraitsParamsR.slot("MutationRate"));

        // sex dependency
        int sexdep = 2; // NA = 2, FEM = 0 , MAL = 1; depending on row

        // Output values
        bool isOutputR = Rcpp::as<bool>(NeutralTraitsParamsR.slot("OutputValues"));

        setUpSpeciesTrait(TraitTypeR,
                          positions,
                          ExpressionTypeR,
                          initDistR,
                          initParamsR,
                          initDomDistR,
                          initDomParamsR,
                          DominanceDistR,
                          DominanceParamsR,
                          isInherited,
                          MutationDistR,
                          MutationParamsR,
                          MutationRateR,
                          sexdep,
                          isOutputR);
    }

    if(gHasGeneticLoad){
        Rcpp::S4 GeneticLoadParamsR("GeneticLoadParams");
        GeneticLoadParamsR = Rcpp::as<Rcpp::S4>(TraitsParamsR.slot("GeneticLoad"));

        // TraitType
        string TraitTypeR = "genetic_load";// each column corresponds to one genetic load

        // nb of loads
        int nbLoads = Rcpp::as<int>(GeneticLoadParamsR.slot("NbGeneticLoads"));

        // Positions and number of Positions
        Rcpp::List PositionsRList = Rcpp::as<Rcpp::List>(GeneticLoadParamsR.slot("Positions"));
        Rcpp::NumericVector NbOfPositionsRvec;
        if (GeneticLoadParamsR.slot("NbOfPositions")!= R_NilValue){
            NbOfPositionsRvec = Rcpp::as<Rcpp::NumericVector> (GeneticLoadParamsR.slot("NbOfPositions"));
        }

        // Expression type
        string ExpressionTypeR = "multiplicative"; // not applicable for genetic loads

        // Initial distribution parameters
        Rcpp::StringVector initDistRvec;
        if (GeneticLoadParamsR.slot("InitialDistribution")!= R_NilValue){
            initDistRvec =  Rcpp::as<Rcpp::StringVector>(GeneticLoadParamsR.slot("InitialDistribution"));
        } else {
            initDistRvec = {"#"};
        }
        Rcpp::Rcout << "initDistRvec read..." << endl;
        Rcpp::NumericMatrix initParamsRmat;
        if(GeneticLoadParamsR.slot("InitialParameters")!= R_NilValue){
            initParamsRmat = Rcpp::as<Rcpp::NumericMatrix>(GeneticLoadParamsR.slot("InitialParameters")); // as a matrix: columns for parameter, rows for load nb
        } else {
            initParamsRmat = Rcpp::NumericMatrix(0, 0);// empty matrix
        }
        Rcpp::Rcout << "initParamsRmat read..." << endl;

        // Initial dominance distribution parameters applicable for genetic loads
        Rcpp::StringVector initDomDistRvec;
        if (GeneticLoadParamsR.slot("InitialDomDistribution")!= R_NilValue){
            initDomDistRvec =  Rcpp::as<Rcpp::StringVector>(GeneticLoadParamsR.slot("InitialDomDistribution"));
        } else {
            initDomDistRvec = {"#"};
        }
        Rcpp::Rcout << "initDomDistRvec read..." << endl;

        Rcpp::NumericMatrix initDomParamsRmat;
        if(GeneticLoadParamsR.slot("InitialDomParameters")!= R_NilValue){
            initDomParamsRmat = Rcpp::as<Rcpp::NumericMatrix>(GeneticLoadParamsR.slot("InitialDomParameters")); // as a matrix: columns for parameter, rows for load nb
        } else {
            initDomParamsRmat = Rcpp::NumericMatrix(0, 0);
        }
        Rcpp::Rcout << "initDomParamsRmat read..." << endl;

        // Dominance distribution parameters
        Rcpp::StringVector DominanceDistRvec = Rcpp::as<Rcpp::StringVector>(GeneticLoadParamsR.slot("DominanceDistribution")); // check if values are NA -> then '#'
        Rcpp::NumericMatrix DominanceParamsRmat = Rcpp::as<Rcpp::NumericMatrix>(GeneticLoadParamsR.slot("DominanceParameters")); // as a matrix: columns for parameter, rows for load nb

        // Mutation parameters
        bool isInherited = true;
        Rcpp::StringVector MutationDistRvec = Rcpp::as<Rcpp::StringVector>(GeneticLoadParamsR.slot("MutationDistribution"));
        Rcpp::NumericMatrix MutationParamsRmat = Rcpp::as<Rcpp::NumericMatrix>(GeneticLoadParamsR.slot("MutationParameters"));
        Rcpp::NumericVector MutationRateRvec = Rcpp::as<Rcpp::NumericVector>(GeneticLoadParamsR.slot("MutationRate"));

        // sex dependency
        int sexdep = 2; // NA = 2, FEM = 0 , MAL = 1; depending on row

        // Output values
        Rcpp::LogicalVector isOutputRvec = Rcpp::as<Rcpp::LogicalVector>(GeneticLoadParamsR.slot("OutputValues"));

        for (int l = 0; l < nbLoads; l++){
            set<int> positions;
            int NbOfPositionsR = -9;
            // check if PositionsR[l] is a string
            if(Rf_isString(PositionsRList[l])){
                string pos = Rcpp::as<string>(PositionsRList[l]);
                if(pos == "random"){
                    NbOfPositionsR = (int)NbOfPositionsRvec[l]; // here is the error message
                    if(NbOfPositionsR > 0){
                        positions = selectRandomLociPositions(NbOfPositionsR, genomeSize);
                    }
                }
                else throw logic_error("GeneticLoad(): If positions are random you must provide the number of positions (>0).");
            }
            else {
                Rcpp::NumericVector PositionsR = Rcpp::as<Rcpp::NumericVector>(PositionsRList[l]); // Positions are provided as numeric vector
                // use a for loop to insert the values at PositionsR into the set positions
                for (int i = 0; i < PositionsR.size(); i++) {
                    if (static_cast<int>(PositionsR[i]) <= genomeSize) {
                        positions.insert(static_cast<int>(PositionsR[i]));
                    }
                    else throw logic_error("GeneticLoad(): Loci positions must be smaller than genome size");
                }
            }


            string initDistR = "#";
            if (initDistRvec[0] != "#") initDistR = (string)initDistRvec[l]; // it is checked beforehand whether the correct distributions are provided
            // else initDistR = "#";

            Rcpp::NumericVector initParamsR = {0,0};
            if (initParamsRmat.size() > 0) initParamsR = initParamsRmat.row(l);
            // else initParamsR = {0,0};

            string initDomDistR = "#";
            if (initDomDistRvec[0] != "#") initDomDistR = (string)initDomDistRvec[l]; // it is checked beforehand whether the correct distributions are provided
            // else initDomDistR = "#";

            Rcpp::NumericVector initDomParamsR = {0,0};
            if (initDomParamsRmat.size() >0) initDomParamsR = initDomParamsRmat.row(l);
            // else initDomParamsR = {0,0};

            string DominanceDistR = (string)DominanceDistRvec[l]; // it is checked beforehand whether the correct distributions are provided
            Rcpp::NumericVector DominanceParamsR = {0.2}; // DominanceParamsRmat.row(l);

            string MutationDistR =  (string)MutationDistRvec[l];
            Rcpp::NumericVector MutationParamsR = MutationParamsRmat.row(l);
            float MutationRateR = (float) MutationRateRvec[l];

            bool isOutputR = isOutputRvec[l] == 1;

            setUpSpeciesTrait(TraitTypeR,
                              positions,
                              ExpressionTypeR,
                          initDistR,
                          initParamsR,
                              initDomDistR,
                              initDomParamsR,
                          DominanceDistR,
                              DominanceParamsR,
                          isInherited,
                              MutationDistR,
                          MutationParamsR,
                          MutationRateR,
                          sexdep,
                          isOutputR);
        }
    }

    // Dispersal traits
    // emigration traits
    if(emig.indVar){
        EmigrationTraitsParamsR = Rcpp::as<Rcpp::S4>(TraitsParamsR.slot("EmigrationGenes"));

        // number of expected trait types
        int nbTraits;
        if (emig.densDep) nbTraits = 3; // emigration_d0, emigration_alpha, emigration_beta
        else nbTraits = 1; // emigration_d0

        if (emig.sexDep) nbTraits *= 2; // each sex has a different trait
        else nbTraits *= 1; // only one trait for both sexes

        // Positions and number of Positions
        Rcpp::List PositionsRList = Rcpp::as<Rcpp::List>(EmigrationTraitsParamsR.slot("Positions"));
        Rcpp::NumericVector NbOfPositionsRvec;
        if (EmigrationTraitsParamsR.slot("NbOfPositions")!= R_NilValue){
            NbOfPositionsRvec = Rcpp::as<Rcpp::NumericVector> (EmigrationTraitsParamsR.slot("NbOfPositions"));
        }

        // Expression type
        Rcpp::StringVector ExpressionTypeRvec = Rcpp::as<Rcpp::StringVector> (EmigrationTraitsParamsR.slot("ExpressionType")); // not applicable for genetic loads

        // Initial distribution
        Rcpp::StringVector InitialDistRvec = Rcpp::as<Rcpp::StringVector>(EmigrationTraitsParamsR.slot("InitialDistribution"));
        Rcpp::NumericMatrix InitialParamsRmat = Rcpp::as<Rcpp::NumericMatrix>(EmigrationTraitsParamsR.slot("InitialParameters")); // as a matrix: columns for parameter, rows for emigration trait

        Rcpp::LogicalVector isInheritedRvec = Rcpp::as<Rcpp::LogicalVector>(EmigrationTraitsParamsR.slot("IsInherited"));

        // Initial dominance distribution
        string initDomDistR = "#"; // not applicable for dispersal traits
        Rcpp::NumericVector initDomParamsR = {0,0}; // not applicable for dispersal traits

        // Dominance distribution
        string DominanceDistR = "#"; // not applicable for dispersal traits
        Rcpp::NumericVector DominanceParamsR = {0,0}; // not applicable for dispersal traits

        // Mutation parameters
        Rcpp::StringVector MutationDistRvec = Rcpp::as<Rcpp::StringVector>(EmigrationTraitsParamsR.slot("MutationDistribution"));
        Rcpp::NumericMatrix MutationParamsRmat = Rcpp::as<Rcpp::NumericMatrix>(EmigrationTraitsParamsR.slot("MutationParameters"));
        Rcpp::NumericVector MutationRateRvec = Rcpp::as<Rcpp::NumericVector>(EmigrationTraitsParamsR.slot("MutationRate"));

        // Output values
        Rcpp::LogicalVector isOutputRvec = Rcpp::as<Rcpp::LogicalVector>(EmigrationTraitsParamsR.slot("OutputValues"));

        string TraitTypeR;
        int sexdep;

        for (int l = 0; l < nbTraits; l++){
            if (emig.densDep){
                if (emig.sexDep){
                    // expecting emigration_d0, emigration_alpha, emigration_beta
                    if (l == 0 || l == 1) TraitTypeR = "emigration_d0";
                    else if (l == 2 || l == 3) TraitTypeR = "emigration_alpha";
                    else if (l == 4 || l == 5) TraitTypeR = "emigration_beta";
                    //if l is even -> female; sexdep=0
                    if(l % 2 == 0) sexdep = 0; // FEM
                    //if l is odd -> male; sexdep=1
                    else sexdep = 1; // MAL
                }
                else {
                    // expecting emigration_d0, emigration_alpha, emigration_beta
                    if (l == 0) TraitTypeR = "emigration_d0";
                    else if (l == 1) TraitTypeR = "emigration_alpha";
                    else if (l == 2) TraitTypeR = "emigration_beta";
                    sexdep=2;
                }

            }
            else {
                // expecting only emigration_d0
                TraitTypeR = "emigration_d0";
                if (emig.sexDep){
                    if(l % 2 == 0) sexdep = 0; // FEM
                    else sexdep = 1; // MAL
                }
                else {
                    sexdep=2;
                }
            }

            set<int> positions;
            int NbOfPositionsR = -9;
            // check if PositionsR[l] is a string
            if(Rf_isString(PositionsRList[l])){
                string pos = Rcpp::as<string>(PositionsRList[l]);
                if(pos == "random"){
                    NbOfPositionsR = (int)NbOfPositionsRvec[l]; // here is the error message
                    if(NbOfPositionsR > 0){
                        positions = selectRandomLociPositions(NbOfPositionsR, genomeSize);
                    }
                }
                else throw logic_error("EmigrationGenes(): If positions are random you must provide the number of positions (>0).");
            }
            else {
                Rcpp::NumericVector PositionsR = Rcpp::as<Rcpp::NumericVector>(PositionsRList[l]); // Positions are provided as numeric vector
                // use a for loop to insert the values at PositionsR into the set positions
                for (int i = 0; i < PositionsR.size(); i++) {
                    if (static_cast<int>(PositionsR[i]) <= genomeSize) {
                        positions.insert(static_cast<int>(PositionsR[i]));
                    }
                    else throw logic_error("EmigrationGenes(): Loci positions must be smaller than genome size");
                }
            }

            string ExpressionTypeR = (string)ExpressionTypeRvec[l];

            string initDistR = (string)InitialDistRvec[l]; // it is checked beforehand whether the correct distributions are provided
            Rcpp::NumericVector initParamsR = InitialParamsRmat.row(l);

            bool isInherited = isInheritedRvec[l] == 1;

            string MutationDistR =  (string)MutationDistRvec[l];
            Rcpp::NumericVector MutationParamsR = MutationParamsRmat.row(l);

            float MutationRateR = (float) MutationRateRvec[l];
            bool isOutputR = isOutputRvec[l] == 1;

            setUpSpeciesTrait(TraitTypeR,
                              positions,
                              ExpressionTypeR,
                              initDistR,
                              initParamsR,
                              initDomDistR,
                              initDomParamsR,
                              DominanceDistR,
                              DominanceParamsR,
                              isInherited,
                              MutationDistR,
                              MutationParamsR,
                              MutationRateR,
                              sexdep,
                              isOutputR);
        }

    }

    // settlement traits
    if(sett.indVar){
        SettlementTraitsParamsR = Rcpp::as<Rcpp::S4>(TraitsParamsR.slot("SettlementGenes"));

        // number of expected trait types
        int nbTraits = 3;

        if (sett.sexDep) nbTraits *= 2; // each sex has a different trait

        // Positions and number of Positions
        Rcpp::List PositionsRList = Rcpp::as<Rcpp::List>(SettlementTraitsParamsR.slot("Positions"));
        Rcpp::NumericVector NbOfPositionsRvec;
        if (SettlementTraitsParamsR.slot("NbOfPositions")!= R_NilValue){
            NbOfPositionsRvec = Rcpp::as<Rcpp::NumericVector> (SettlementTraitsParamsR.slot("NbOfPositions"));
        }

        // Expression type
        Rcpp::StringVector ExpressionTypeRvec = Rcpp::as<Rcpp::StringVector> (SettlementTraitsParamsR.slot("ExpressionType")); // not applicable for genetic loads

        // Initial distribution
        Rcpp::StringVector InitialDistRvec = Rcpp::as<Rcpp::StringVector>(SettlementTraitsParamsR.slot("InitialDistribution"));
        Rcpp::NumericMatrix InitialParamsRmat = Rcpp::as<Rcpp::NumericMatrix>(SettlementTraitsParamsR.slot("InitialParameters")); // as a matrix: columns for parameter, rows for Settlement trait

        Rcpp::LogicalVector isInheritedRvec = Rcpp::as<Rcpp::LogicalVector>(SettlementTraitsParamsR.slot("IsInherited"));

        // Initial dominance distribution
        string initDomDistR = "#"; // not applicable for dispersal traits
        Rcpp::NumericVector initDomParamsR = {0,0}; // not applicable for dispersal traits

        // Dominance distribution
        string DominanceDistR = "#"; // not applicable for dispersal traits
        Rcpp::NumericVector DominanceParamsR = {0,0}; // not applicable for dispersal traits

        // Mutation parameters
        Rcpp::StringVector MutationDistRvec = Rcpp::as<Rcpp::StringVector>(SettlementTraitsParamsR.slot("MutationDistribution"));
        Rcpp::NumericMatrix MutationParamsRmat = Rcpp::as<Rcpp::NumericMatrix>(SettlementTraitsParamsR.slot("MutationParameters"));
        Rcpp::NumericVector MutationRateRvec = Rcpp::as<Rcpp::NumericVector>(SettlementTraitsParamsR.slot("MutationRate"));

        // Output values
        Rcpp::LogicalVector isOutputRvec = Rcpp::as<Rcpp::LogicalVector>(SettlementTraitsParamsR.slot("OutputValues"));

        string TraitTypeR;
        int sexdep;

        for (int l = 0; l < nbTraits; l++){
                if (sett.sexDep){
                    // expecting emigration_d0, emigration_alpha, emigration_beta
                    if (l == 0 || l == 1) TraitTypeR = "settlement_s0";
                    else if (l == 2 || l == 3) TraitTypeR = "settlement_alpha";
                    else if (l == 4 || l == 5) TraitTypeR = "settlement_beta";
                    //if l is even -> female; sexdep=0
                    if(l % 2 == 0) sexdep = 0; // FEM
                    //if l is odd -> male; sexdep=1
                    else sexdep = 1; // MAL
                }
                else {
                    // expecting emigration_d0, emigration_alpha, emigration_beta
                    if (l == 0) TraitTypeR = "settlement_d0";
                    else if (l == 1) TraitTypeR = "settlement_alpha";
                    else if (l == 2) TraitTypeR = "settlement_beta";
                    sexdep=2;
                }

            set<int> positions;
            int NbOfPositionsR = -9;
            // check if PositionsR[l] is a string
            if(Rf_isString(PositionsRList[l])){
                string pos = Rcpp::as<string>(PositionsRList[l]);
                if(pos == "random"){
                    NbOfPositionsR = (int)NbOfPositionsRvec[l]; // here is the error message
                    if(NbOfPositionsR > 0){
                        positions = selectRandomLociPositions(NbOfPositionsR, genomeSize);
                    }
                }
                else throw logic_error("SettlementGenes(): If positions are random you must provide the number of positions (>0).");
            }
            else {
                Rcpp::NumericVector PositionsR = Rcpp::as<Rcpp::NumericVector>(PositionsRList[l]); // Positions are provided as numeric vector
                // use a for loop to insert the values at PositionsR into the set positions
                for (int i = 0; i < PositionsR.size(); i++) {
                    if (static_cast<int>(PositionsR[i]) <= genomeSize) {
                        positions.insert(static_cast<int>(PositionsR[i]));
                    }
                    else throw logic_error("SettlementGenes(): Loci positions must be smaller than genome size");
                }
            }

            string ExpressionTypeR = (string)ExpressionTypeRvec[l];

            string initDistR = (string)InitialDistRvec[l]; // it is checked beforehand whether the correct distributions are provided
            Rcpp::NumericVector initParamsR = InitialParamsRmat.row(l);

            bool isInherited = isInheritedRvec[l] == 1;

            string MutationDistR =  (string)MutationDistRvec[l];
            Rcpp::NumericVector MutationParamsR = MutationParamsRmat.row(l);

            float MutationRateR = (float) MutationRateRvec[l];
            bool isOutputR = isOutputRvec[l] == 1;

            setUpSpeciesTrait(TraitTypeR,
                              positions,
                              ExpressionTypeR,
                              initDistR,
                              initParamsR,
                              initDomDistR,
                              initDomParamsR,
                              DominanceDistR,
                              DominanceParamsR,
                              isInherited,
                              MutationDistR,
                              MutationParamsR,
                              MutationRateR,
                              sexdep,
                              isOutputR);
        }

    }

    switch(TransferType){
        case 0: { // kernel
            if(trfr.indVar){
                KernelTraitsParamsR = Rcpp::as<Rcpp::S4>(TraitsParamsR.slot("KernelGenes"));

                // number of expected traits
                int nbTraits;
                if(trfr.twinKern){
                    nbTraits = 3;
                }
                else {
                    nbTraits = 1;
                }
                if(trfr.sexDep){
                    nbTraits *= 2;
                }

                // Positions and number of Positions
                Rcpp::List PositionsRList = Rcpp::as<Rcpp::List>(KernelTraitsParamsR.slot("Positions"));
                Rcpp::NumericVector NbOfPositionsRvec;
                if (KernelTraitsParamsR.slot("NbOfPositions")!= R_NilValue){
                    NbOfPositionsRvec = Rcpp::as<Rcpp::NumericVector> (KernelTraitsParamsR.slot("NbOfPositions"));
                }

                // Expression type
                Rcpp::StringVector ExpressionTypeRvec = Rcpp::as<Rcpp::StringVector> (KernelTraitsParamsR.slot("ExpressionType")); // not applicable for genetic loads

                // Initial distribution
                Rcpp::StringVector InitialDistRvec = Rcpp::as<Rcpp::StringVector>(KernelTraitsParamsR.slot("InitialDistribution"));
                Rcpp::NumericMatrix InitialParamsRmat = Rcpp::as<Rcpp::NumericMatrix>(KernelTraitsParamsR.slot("InitialParameters")); // as a matrix: columns for parameter, rows for Settlement trait

                Rcpp::LogicalVector isInheritedRvec = Rcpp::as<Rcpp::LogicalVector>(KernelTraitsParamsR.slot("IsInherited"));

                // Initial dominance distribution
                string initDomDistR = "#"; // not applicable for dispersal traits
                Rcpp::NumericVector initDomParamsR = {0,0}; // not applicable for dispersal traits

                // Dominance distribution
                string DominanceDistR = "#"; // not applicable for dispersal traits
                Rcpp::NumericVector DominanceParamsR = {0,0}; // not applicable for dispersal traits

                // Mutation parameters
                Rcpp::StringVector MutationDistRvec = Rcpp::as<Rcpp::StringVector>(KernelTraitsParamsR.slot("MutationDistribution"));
                Rcpp::NumericMatrix MutationParamsRmat = Rcpp::as<Rcpp::NumericMatrix>(KernelTraitsParamsR.slot("MutationParameters"));
                Rcpp::NumericVector MutationRateRvec = Rcpp::as<Rcpp::NumericVector>(KernelTraitsParamsR.slot("MutationRate"));

                // Output values
                Rcpp::LogicalVector isOutputRvec = Rcpp::as<Rcpp::LogicalVector>(KernelTraitsParamsR.slot("OutputValues"));

                string TraitTypeR;
                int sexdep;

                for (int l = 0; l < nbTraits; l++){
                    // expecting crw_stepLength, crw_stepCorrelation
                    if(trfr.twinKern){
                        if(trfr.sexDep){
                            // should be two lines
                            if(l==0 || l == 1)    TraitTypeR = "kernel_meanDistance1";
                            if(l==2 || l == 3)    TraitTypeR = "kernel_meanDistance2";
                            if(l==4 || l == 5)    TraitTypeR = "kernel_probability";
                            if(l % 2 == 0) sexdep = 0; // FEML
                            else sexdep = 1; // MAL
                        }
                        else {
                            if(l==0)    TraitTypeR = "kernel_meanDistance1";
                            if(l==1)    TraitTypeR = "kernel_meanDistance2";
                            if(l==2)    TraitTypeR = "kernel_probability";
                            sexdep = 2;
                        }
                    }
                    else {
                        if(trfr.sexDep){
                            // should be two lines
                            TraitTypeR = "kernel_meanDistance1";
                            if(l % 2 == 0) sexdep = 0; // FEML
                            else sexdep = 1; // MAL
                        }
                        else {
                            //should only be one line
                            TraitTypeR = "kernel_meanDistance1";
                            sexdep = 2;
                        }
                    }


                    set<int> positions;
                    int NbOfPositionsR = -9;
                    // check if PositionsR[l] is a string
                    if(Rf_isString(PositionsRList[l])){
                        string pos = Rcpp::as<string>(PositionsRList[l]);
                        if(pos == "random"){
                            NbOfPositionsR = (int)NbOfPositionsRvec[l]; // here is the error message
                            if(NbOfPositionsR > 0){
                                positions = selectRandomLociPositions(NbOfPositionsR, genomeSize);
                            }
                        }
                        else throw logic_error("KernelGenes(): If positions are random you must provide the number of positions (>0).");
                    }
                    else {
                        Rcpp::NumericVector PositionsR = Rcpp::as<Rcpp::NumericVector>(PositionsRList[l]); // Positions are provided as numeric vector
                        // use a for loop to insert the values at PositionsR into the set positions
                        for (int i = 0; i < PositionsR.size(); i++) {
                            if (static_cast<int>(PositionsR[i]) <= genomeSize) {
                                positions.insert(static_cast<int>(PositionsR[i]));
                            }
                            else throw logic_error("KernelGenes(): Loci positions must be smaller than genome size");
                        }
                    }

                    string ExpressionTypeR = (string)ExpressionTypeRvec[l];

                    string initDistR = (string)InitialDistRvec[l]; // it is checked beforehand whether the correct distributions are provided
                    Rcpp::NumericVector initParamsR = InitialParamsRmat.row(l);

                    bool isInherited = isInheritedRvec[l] == 1;

                    string MutationDistR =  (string)MutationDistRvec[l];
                    Rcpp::NumericVector MutationParamsR = MutationParamsRmat.row(l);

                    float MutationRateR = (float) MutationRateRvec[l];
                    bool isOutputR = isOutputRvec[l] == 1;

                    setUpSpeciesTrait(TraitTypeR,
                                      positions,
                                      ExpressionTypeR,
                                      initDistR,
                                      initParamsR,
                                      initDomDistR,
                                      initDomParamsR,
                                      DominanceDistR,
                                      DominanceParamsR,
                                      isInherited,
                                      MutationDistR,
                                      MutationParamsR,
                                      MutationRateR,
                                      sexdep,
                                      isOutputR);
                }
            }
        }
        break;

        case 1: { // SMS
            // SMS traits
            trfrMovtParams smstraits = pSpecies->getSpMovtTraits();
            if(trfr.indVar){
                SMSTraitsParamsR = Rcpp::as<Rcpp::S4>(TraitsParamsR.slot("SMSGenes"));

                // number of expected traits
                int nbTraits = 1;
                if (smstraits.gb==2) nbTraits = 4;


                // Positions and number of Positions
                Rcpp::List PositionsRList = Rcpp::as<Rcpp::List>(SMSTraitsParamsR.slot("Positions"));
                Rcpp::NumericVector NbOfPositionsRvec;
                if (SMSTraitsParamsR.slot("NbOfPositions")!= R_NilValue){
                    NbOfPositionsRvec = Rcpp::as<Rcpp::NumericVector> (SMSTraitsParamsR.slot("NbOfPositions"));
                }

                // Expression type
                Rcpp::StringVector ExpressionTypeRvec = Rcpp::as<Rcpp::StringVector> (SMSTraitsParamsR.slot("ExpressionType")); // not applicable for genetic loads

                // Initial distribution
                Rcpp::StringVector InitialDistRvec = Rcpp::as<Rcpp::StringVector>(SMSTraitsParamsR.slot("InitialDistribution"));
                Rcpp::NumericMatrix InitialParamsRmat = Rcpp::as<Rcpp::NumericMatrix>(SMSTraitsParamsR.slot("InitialParameters")); // as a matrix: columns for parameter, rows for Settlement trait

                Rcpp::LogicalVector isInheritedRvec = Rcpp::as<Rcpp::LogicalVector>(SMSTraitsParamsR.slot("IsInherited"));

                // Initial dominance distribution
                string initDomDistR = "#"; // not applicable for dispersal traits
                Rcpp::NumericVector initDomParamsR = {0,0}; // not applicable for dispersal traits

                // Dominance distribution
                string DominanceDistR = "#"; // not applicable for dispersal traits
                Rcpp::NumericVector DominanceParamsR = {0,0}; // not applicable for dispersal traits

                // Mutation parameters
                Rcpp::StringVector MutationDistRvec = Rcpp::as<Rcpp::StringVector>(SMSTraitsParamsR.slot("MutationDistribution"));
                Rcpp::NumericMatrix MutationParamsRmat = Rcpp::as<Rcpp::NumericMatrix>(SMSTraitsParamsR.slot("MutationParameters"));
                Rcpp::NumericVector MutationRateRvec = Rcpp::as<Rcpp::NumericVector>(SMSTraitsParamsR.slot("MutationRate"));

                // Output values
                Rcpp::LogicalVector isOutputRvec = Rcpp::as<Rcpp::LogicalVector>(SMSTraitsParamsR.slot("OutputValues"));

                string TraitTypeR;
                int sexdep;

                for (int l = 0; l < nbTraits; l++){
                    // expecting crw_stepLength, crw_stepCorrelation
                    if (l == 0) {
                        TraitTypeR = "sms_directionalPersistence";
                    }
                    else if (l == 1) {
                        TraitTypeR = "sms_goalBias";
                    }
                    else if (l == 2) {
                        TraitTypeR = "sms_alphaDB";
                    }
                    else {
                        TraitTypeR = "sms_betaDB";
                    }

                    sexdep=2;

                    set<int> positions;
                    int NbOfPositionsR = -9;
                    // check if PositionsR[l] is a string
                    if(Rf_isString(PositionsRList[l])){
                        string pos = Rcpp::as<string>(PositionsRList[l]);
                        if(pos == "random"){
                            NbOfPositionsR = (int)NbOfPositionsRvec[l]; // here is the error message
                            if(NbOfPositionsR > 0){
                                positions = selectRandomLociPositions(NbOfPositionsR, genomeSize);
                            }
                        }
                        else throw logic_error("SMSGenes(): If positions are random you must provide the number of positions (>0).");
                    }
                    else {
                        Rcpp::NumericVector PositionsR = Rcpp::as<Rcpp::NumericVector>(PositionsRList[l]); // Positions are provided as numeric vector
                        // use a for loop to insert the values at PositionsR into the set positions
                        for (int i = 0; i < PositionsR.size(); i++) {
                            if (static_cast<int>(PositionsR[i]) <= genomeSize) {
                                positions.insert(static_cast<int>(PositionsR[i]));
                            }
                            else throw logic_error("SMSGenes(): Loci positions must be smaller than genome size");
                        }
                    }

                    string ExpressionTypeR = (string)ExpressionTypeRvec[l];

                    string initDistR = (string)InitialDistRvec[l]; // it is checked beforehand whether the correct distributions are provided
                    Rcpp::NumericVector initParamsR = InitialParamsRmat.row(l);

                    bool isInherited = isInheritedRvec[l] == 1;

                    string MutationDistR =  (string)MutationDistRvec[l];
                    Rcpp::NumericVector MutationParamsR = MutationParamsRmat.row(l);

                    float MutationRateR = (float) MutationRateRvec[l];
                    bool isOutputR = isOutputRvec[l] == 1;

                    setUpSpeciesTrait(TraitTypeR,
                                      positions,
                                      ExpressionTypeR,
                                      initDistR,
                                      initParamsR,
                                      initDomDistR,
                                      initDomParamsR,
                                      DominanceDistR,
                                      DominanceParamsR,
                                      isInherited,
                                      MutationDistR,
                                      MutationParamsR,
                                      MutationRateR,
                                      sexdep,
                                      isOutputR);
                }

            }
        }
        break;

        case 2: { //CorrRW
            // CorrRW traits
            if(trfr.indVar){
                CorrRWTraitsParamsR = Rcpp::as<Rcpp::S4>(TraitsParamsR.slot("CorrRWGenes"));

                // number of expected traits
                int nbTraits = 2;

                // Positions and number of Positions
                Rcpp::List PositionsRList = Rcpp::as<Rcpp::List>(CorrRWTraitsParamsR.slot("Positions"));
                Rcpp::NumericVector NbOfPositionsRvec;
                if (CorrRWTraitsParamsR.slot("NbOfPositions")!= R_NilValue){
                    NbOfPositionsRvec = Rcpp::as<Rcpp::NumericVector> (CorrRWTraitsParamsR.slot("NbOfPositions"));
                }

                // Expression type
                Rcpp::StringVector ExpressionTypeRvec = Rcpp::as<Rcpp::StringVector> (CorrRWTraitsParamsR.slot("ExpressionType")); // not applicable for genetic loads

                // Initial distribution
                Rcpp::StringVector InitialDistRvec = Rcpp::as<Rcpp::StringVector>(CorrRWTraitsParamsR.slot("InitialDistribution"));
                Rcpp::NumericMatrix InitialParamsRmat = Rcpp::as<Rcpp::NumericMatrix>(CorrRWTraitsParamsR.slot("InitialParameters")); // as a matrix: columns for parameter, rows for Settlement trait

                Rcpp::LogicalVector isInheritedRvec = Rcpp::as<Rcpp::LogicalVector>(CorrRWTraitsParamsR.slot("IsInherited"));

                // Initial dominance distribution
                string initDomDistR = "#"; // not applicable for dispersal traits
                Rcpp::NumericVector initDomParamsR = {0,0}; // not applicable for dispersal traits

                // Dominance distribution
                string DominanceDistR = "#"; // not applicable for dispersal traits
                Rcpp::NumericVector DominanceParamsR = {0,0}; // not applicable for dispersal traits

                // Mutation parameters
                Rcpp::StringVector MutationDistRvec = Rcpp::as<Rcpp::StringVector>(CorrRWTraitsParamsR.slot("MutationDistribution"));
                Rcpp::NumericMatrix MutationParamsRmat = Rcpp::as<Rcpp::NumericMatrix>(CorrRWTraitsParamsR.slot("MutationParameters"));
                Rcpp::NumericVector MutationRateRvec = Rcpp::as<Rcpp::NumericVector>(CorrRWTraitsParamsR.slot("MutationRate"));

                // Output values
                Rcpp::LogicalVector isOutputRvec = Rcpp::as<Rcpp::LogicalVector>(CorrRWTraitsParamsR.slot("OutputValues"));

                string TraitTypeR;
                int sexdep;

                for (int l = 0; l < nbTraits; l++){
                    // expecting crw_stepLength, crw_stepCorrelation
                    if (l == 0) {
                        TraitTypeR = "crw_stepLength";
                    }
                    else {
                        TraitTypeR = "crw_stepCorrelation";
                    }
                    sexdep=2;

                    set<int> positions;
                    int NbOfPositionsR = -9;
                    // check if PositionsR[l] is a string
                    if(Rf_isString(PositionsRList[l])){
                        string pos = Rcpp::as<string>(PositionsRList[l]);
                        if(pos == "random"){
                            NbOfPositionsR = (int)NbOfPositionsRvec[l]; // here is the error message
                            if(NbOfPositionsR > 0){
                                positions = selectRandomLociPositions(NbOfPositionsR, genomeSize);
                            }
                        }
                        else throw logic_error("CorrRWGenes(): If positions are random you must provide the number of positions (>0).");
                    }
                    else {
                        Rcpp::NumericVector PositionsR = Rcpp::as<Rcpp::NumericVector>(PositionsRList[l]); // Positions are provided as numeric vector
                        // use a for loop to insert the values at PositionsR into the set positions
                        for (int i = 0; i < PositionsR.size(); i++) {
                            if (static_cast<int>(PositionsR[i]) <= genomeSize) {
                                positions.insert(static_cast<int>(PositionsR[i]));
                            }
                            else throw logic_error("CorrRWGenes(): Loci positions must be smaller than genome size");
                        }
                    }

                    string ExpressionTypeR = (string)ExpressionTypeRvec[l];

                    string initDistR = (string)InitialDistRvec[l]; // it is checked beforehand whether the correct distributions are provided
                    Rcpp::NumericVector initParamsR = InitialParamsRmat.row(l);

                    bool isInherited = isInheritedRvec[l] == 1;

                    string MutationDistR =  (string)MutationDistRvec[l];
                    Rcpp::NumericVector MutationParamsR = MutationParamsRmat.row(l);

                    float MutationRateR = (float) MutationRateRvec[l];
                    bool isOutputR = isOutputRvec[l] == 1;

                    setUpSpeciesTrait(TraitTypeR,
                                      positions,
                                      ExpressionTypeR,
                                      initDistR,
                                      initParamsR,
                                      initDomDistR,
                                      initDomParamsR,
                                      DominanceDistR,
                                      DominanceParamsR,
                                      isInherited,
                                      MutationDistR,
                                      MutationParamsR,
                                      MutationRateR,
                                      sexdep,
                                      isOutputR);
                }

            }
        }
        break;
    }


    // SMS traits

    // Kernel traits

    return 0;

}
//---------------------------------------------------------------------------
// Convert integer to corresponding SexType value
const sex_t intToSex(const int& i) {
    if (i == 0) return FEM;
    else if (i == 1) return MAL;
    else return NA;

}
//---------------------------------------------------------------------------
// Convert string to corresponding TraitType value, if valid
TraitType stringToTraitType(const std::string& str) {
    // Non-dispersal traits
    if (str == "neutral") return NEUTRAL;
    else if (str == "genetic_load") return GENETIC_LOAD;
    // Sex-invariant dispersal traits
    else if (str == "emigration_d0") return E_D0; // EP uses d0 for trait data
    else if (str == "emigration_alpha") return E_ALPHA;
    else if (str == "emigration_beta") return E_BETA;
    else if (str == "settlement_s0") return S_S0;
    else if (str == "settlement_alpha") return S_ALPHA;
    else if (str == "settlement_beta") return S_BETA;
    else if (str == "kernel_meanDistance1") return KERNEL_MEANDIST_1;
    else if (str == "kernel_meanDistance2") return KERNEL_MEANDIST_2;
    else if (str == "kernel_probability") return KERNEL_PROBABILITY;
    else if (str == "crw_stepLength") return CRW_STEPLENGTH;
    else if (str == "crw_stepCorrelation") return CRW_STEPCORRELATION;
    else if (str == "sms_directionalPersistence") return SMS_DP;
    else if (str == "sms_goalBias") return SMS_GB;
    else if (str == "sms_alphaDB") return SMS_ALPHADB;
    else if (str == "sms_betaDB") return SMS_BETADB;
    else return INVALID_TRAIT;
}
//---------------------------------------------------------------------------
// Add sex to TraitType
TraitType addSexDepToTrait(const TraitType& t, const sex_t& sex) {
    if (sex == FEM) {
        if (t == E_D0) return E_D0_F; // EP uses d0 for trait data
        else if (t == E_ALPHA) return E_ALPHA_F;
        else if (t == E_BETA) return E_BETA_F;
        else if (t == S_S0) return S_S0_F;
        else if (t == S_ALPHA) return S_ALPHA_F;
        else if (t == S_BETA) return S_BETA_F;
        else if (t == KERNEL_MEANDIST_1) return KERNEL_MEANDIST_1_F;
        else if (t == KERNEL_MEANDIST_2) return KERNEL_MEANDIST_2_F;
        else if (t == KERNEL_PROBABILITY) return KERNEL_PROBABILITY_F;
        else return INVALID_TRAIT;
    }
    else if (sex == MAL) {
        if (t == E_D0) return E_D0_M; // EP uses d0 for trait data
        else if (t == E_ALPHA) return E_ALPHA_M;
        else if (t == E_BETA) return E_BETA_M;
        else if (t == S_S0) return S_S0_M;
        else if (t == S_ALPHA) return S_ALPHA_M;
        else if (t == S_BETA) return S_BETA_M;
        else if (t == KERNEL_MEANDIST_1) return KERNEL_MEANDIST_1_M;
        else if (t == KERNEL_MEANDIST_2) return KERNEL_MEANDIST_2_M;
        else if (t == KERNEL_PROBABILITY) return KERNEL_PROBABILITY_M;
        else return INVALID_TRAIT;
    }
    else return INVALID_TRAIT;
}

//---------------------------------------------------------------------------
GenParamType strToGenParamType(const string& str) {
    if (str == "mean")
        return MEAN;
    else if (str == "sd")
        return SD;
    else if (str == "min")
        return MIN;
    else if (str == "max")
        return MAX;
    else if (str == "shape")
        return SHAPE;
    else if (str == "scale")
        return SCALE;
    else return INVALID;
}
//---------------------------------------------------------------------------
map<GenParamType, float> NumericToParameterMap(string distributionString, Rcpp::NumericVector parameter) {
    // string is the distribution type; distribution type determines which parameter is in first or second column

    map<GenParamType, float> paramMap;
    if (distributionString != "#") {
        if (distributionString == "uniform"){
            paramMap.emplace(GenParamType::MIN, (float) parameter[0]);
            paramMap.emplace(GenParamType::MAX, (float) parameter[1]);
        }
        else if (distributionString == "normal"){
            paramMap.emplace(GenParamType::MEAN, (float) parameter[0]);
            paramMap.emplace(GenParamType::SD, (float) parameter[1]);
        }
        else if (distributionString == "gamma"){
            paramMap.emplace(GenParamType::SHAPE, (float) parameter[0]);
            paramMap.emplace(GenParamType::SCALE, (float) parameter[1]);
        }
        else if (distributionString == "scaled"){
            paramMap.emplace(GenParamType::MEAN, (float) parameter[0]);
        }
        else if (distributionString == "negExp"){
            paramMap.emplace(GenParamType::MEAN, (float) parameter[0]);
        }
        else if (distributionString == "KAM"){
            paramMap.emplace(GenParamType::MAX, (float) parameter[0]);
        }
        else if (distributionString == "SSM"){
            paramMap.emplace(GenParamType::MAX, (float) parameter[0]);
        }
    }
    return paramMap;
}

//---------------------------------------------------------------------------
// select random loci positions
set<int> selectRandomLociPositions(int nbLoci, const int& genomeSize) {
    set<int> positions;
    if (nbLoci > genomeSize) throw logic_error("Number of random loci exceeds genome size.");
    int rndLocus;
    for (int i = 0; i < nbLoci; ++i)
    {
        do {
            rndLocus = pRandom->IRandom(0, genomeSize - 1);
        } while (positions.contains(rndLocus));
        positions.insert(rndLocus);
    }
    return positions;
}
//---------------------------------------------------------------------------
// Convert string to corresponding ExpressionType value, if valid
ExpressionType stringToExpressionType(const std::string& str) {
    if (str == "average") return AVERAGE;
    else if (str == "additive") return ADDITIVE;
    else if (str == "multiplicative") return MULTIPLICATIVE;
    else if (str == "#") return NOTEXPR;
    else throw logic_error(str + " is not a valid gene expression type.");
}
//---------------------------------------------------------------------------
// Convert string to corresponding DistributionType value, if valid
DistributionType stringToDistributionType(const std::string& str) {
    if (str == "#") return NONE;
    else if (str == "uniform") return UNIFORM;
    else if (str == "normal") return NORMAL;
    else if (str == "gamma") return GAMMA;
    else if (str == "scaled") return SCALED;
    else if (str == "negExp") return NEGEXP;
    else if (str == "KAM") return KAM;
    else if (str == "SSM") return SSM;
    else throw logic_error(str + " is not a valid distribution type.");
}
//---------------------------------------------------------------------------
// set up genes
// Set up a trait from input parameters and add it Species
void setUpSpeciesTrait(string TraitTypeR,
                       set<int> positions,
                       string ExpressionTypeR,
                       string initDistR,
                       Rcpp::NumericVector initParamsR,
                       string initDominanceDistR,
                       Rcpp::NumericVector initDominanceParamsR,
                       string DominanceDistR,
                       Rcpp::NumericVector DominanceParamsR,
                       bool isInherited,
                       string MutationDistR,
                       Rcpp::NumericVector MutationParamsR,
                    float MutationRateR,
                       int sexdep,
                       bool isOutputR) {

    TraitType traitType = stringToTraitType(TraitTypeR);
    const sex_t sex = intToSex(sexdep);
    if (sex != NA) traitType = addSexDepToTrait(traitType, sex);

    const ExpressionType expressionType = stringToExpressionType(ExpressionTypeR); // NOTEXPR;

    const DistributionType initDist = stringToDistributionType(initDistR);
    const map<GenParamType, float> initParams = NumericToParameterMap(initDistR,initParamsR);

    const DistributionType initDomDist = stringToDistributionType(initDominanceDistR);
    const map<GenParamType, float> initDomParams = NumericToParameterMap(initDominanceDistR,initDominanceParamsR);

    const DistributionType dominanceDist = stringToDistributionType(DominanceDistR);
    const map<GenParamType, float> dominanceParams = NumericToParameterMap(DominanceDistR,DominanceParamsR);

    DistributionType mutationDistribution = isInherited ?
    stringToDistributionType(MutationDistR) :
        DistributionType::NONE;

    float mutationRate = isInherited ? MutationRateR : 0.0;

    map<GenParamType, float> mutationParameters;
    if (isInherited) {
        mutationParameters = NumericToParameterMap(MutationDistR, MutationParamsR);
    }

    int ploidy = gNbSexesDisp;

    const bool isOutput = isOutputR;

    // Create species trait
    unique_ptr<SpeciesTrait> trait(new SpeciesTrait(
        traitType, sex,
        positions, expressionType,
        initDist, initParams,
            initDomDist, initDomParams,
        isInherited, mutationRate,
        mutationDistribution, mutationParameters,
            dominanceDist, dominanceParams,
        ploidy,
        isOutput
    ));

    pSpecies->addTrait(traitType, *trait);
};
//---------------------------------------------------------------------------


int ReadTranslocationR(Landscape* pLandscape, Rcpp::S4 ParMaster)
{
    int error = 0;
    // create new Management
    if(pManagement != NULL)
        delete pManagement;
    pManagement = new Management;
    // get landscape parameter (to distinguish between patch and cell model)
    landParams paramsLand = pLandscape->getLandParams();
    // get demographic parameter (to distinguish between sexual and asexual reproduction, and stage structured or not)
    demogrParams dem = pSpecies->getDemogrParams();
    // get simulation parameter (to get the number of years)
    simParams sim = paramsSim->getSim();

    // get default values
    managementParams m = pManagement->getManagementParams();
    translocationParams t = pManagement->getTranslocationParams();


    Rcpp::S4 ManageParamsR("ManagementParams");
    ManageParamsR = Rcpp::as<Rcpp::S4>(ParMaster.slot("management"));
    Rcpp::S4 TranslocationParamsR("TranslocationParams");
    TranslocationParamsR = Rcpp::as<Rcpp::S4>(ManageParamsR.slot("Translocation"));

    double catching_rate_R = Rcpp::as<double>(TranslocationParamsR.slot("catching_rate"));

    Rcpp::IntegerVector translocation_years_R;
    translocation_years_R = Rcpp::as<Rcpp::IntegerVector>(TranslocationParamsR.slot("years"));

    Rcpp::IntegerMatrix translocation_matrix_R;
    translocation_matrix_R = Rcpp::as<Rcpp::IntegerMatrix>(TranslocationParamsR.slot("TransLocMat"));

    // activate translocation if translocation_years_R is not empty
    if(translocation_years_R[0] > 0) m.translocation = true;

    // set and update the management parameters
    pManagement->setManagementParams(m);

    if(m.translocation) {
        // check if catcing rate as valid value:
        if(catching_rate_R > 0 || catching_rate_R <= 1) {
            // parse catching_rate_R to catching_rate
            t.catching_rate = catching_rate_R;
        } else{
            Rcpp::Rcout << "ReadTranslocationR(): catching_rate must be between 0 and 1" << std::endl;
            error = 600;
            return error;
        }


        // parse translocation_years_R to translocation_years
        for (int i = 0; i < translocation_years_R.size(); i++) {
            // check if the year is within the simulated years
            if(translocation_years_R[i] > 0 || translocation_years_R[i] <= sim.years) {
                // check if year is not already in the vector
                if(std::find(t.translocation_years.begin(), t.translocation_years.end(), translocation_years_R[i]) == t.translocation_years.end()) {
                    t.translocation_years.push_back(translocation_years_R[i]);
                }else{
                    Rcpp::Rcout << "ReadTranslocationR(): translocation_years[" << i << "]=" << translocation_years_R[i] << " already exists. Make sure that you only include each year once." << std::endl;
                    error = 600;
                    return error;
                }
            } else{
                Rcpp::Rcout << "ReadTranslocationR(): translocation_years[" << i << "]=" << translocation_years_R[i] << " is not within the simulated years" << std::endl;
                error = 600;
                return error;
            }

        }

        for (int i = 0; i < t.translocation_years.size(); i++) {
            Rcpp::Rcout << "ReadTranslocationR(): t.translocation_years[" << i << "]=" << t.translocation_years[i] << std::endl;
        }



        Rcpp::Rcout << "ReadTranslocationR(): catching_rate: "<< t.catching_rate << std::endl;



        for (int i = 0; i < translocation_matrix_R.nrow(); ++i) {

            // extract the year in the first column and initialize the year for each map if needed
            int year = translocation_matrix_R(i,0);
            if (t.source.find(year) == t.source.end()) {
                // not found so add a new key
                t.source.insert(std::pair<int, std::vector<locn>>(year, std::vector<locn>()));
                t.target.insert(std::pair<int, std::vector<locn>>(year, std::vector<locn>()));
                t.nb.insert(std::pair<int, std::vector<int>>(year, std::vector<int>()));
                // if(dem.stageStruct) {
                t.min_age.insert(std::pair<int, std::vector<int>>(year, std::vector<int>()));
                t.max_age.insert(std::pair<int, std::vector<int>>(year, std::vector<int>()));
                t.stage.insert(std::pair<int, std::vector<int>>(year, std::vector<int>()));
                // }
                // if(dem.repType!=0) {
                t.sex.insert(std::pair<int, std::vector<int>>(year, std::vector<int>()));
                // }

            }

            // push_back the source to the source map
            locn s;
            if(paramsLand.patchModel){ // if patch model, the x is the patch ID
                // only if patch ID exists? otherwise exit?
                if(translocation_matrix_R(i,1) <= pLandscape->patchCount() && translocation_matrix_R(i,1) > 0){
                    s.x = translocation_matrix_R(i,1);
                    s.y = -9;
                } else{
                    Rcpp::Rcout << "ReadTranslocationR(): patch ID of source location " << translocation_matrix_R(i,1) << " is either 0 or greater than the number of overall patches." << std::endl;
                    error = 600;
                    return error;
                }

            } else {
                // if cell-based model
                // check if x and y values are inside boundaries
                bool data = false;
                data = pLandscape->checkDataCell(translocation_matrix_R(i,1), translocation_matrix_R(i,2));
                if(data == false){ // cell is out of boundary
                    Rcpp::Rcout << "ReadTranslocationR(): source location " << translocation_matrix_R(i,1) << " is outside the landscape." << std::endl;
                    error = 600;
                    return error;
                } else{ // cell is within landscape
                    s.x = translocation_matrix_R(i,1);
                    s.y = translocation_matrix_R(i,2);
                };
            };

            t.source[year].push_back(s);

            // push_back the target to the target map
            if(paramsLand.patchModel){ // if patch model, the x is the patch ID
                // only if patch ID exists? otherwise exit?
                if(translocation_matrix_R(i,2) <= pLandscape->patchCount() && translocation_matrix_R(i,2) > 0){
                    s.x = translocation_matrix_R(i,2);
                    s.y = -9;
                } else{
                    Rcpp::Rcout << "ReadTranslocationR(): patch ID of target location " << translocation_matrix_R(i,2) << " is either 0 or greater than the number of overall patches." << std::endl;
                    error = 600;
                    return error;
                }
            } else {
                // if cell-based model
                // check if x and y values are inside boundaries and habitat
                if(translocation_matrix_R(i,3) > 0 && translocation_matrix_R(i,3) <= paramsLand.maxX){
                    s.x = translocation_matrix_R(i,3);
                } else{
                    Rcpp::Rcout << "ReadTranslocationR(): x value of target location " << translocation_matrix_R(i,3) << " is either 0 or greater than the maximum x value." << std::endl;
                    error = 600;
                    return error;
                }
                if(translocation_matrix_R(i,4) > 0 && translocation_matrix_R(i,4) <= paramsLand.maxY){
                    s.y = translocation_matrix_R(i,4);
                } else{
                    Rcpp::Rcout << "ReadTranslocationR(): y value of target location " << translocation_matrix_R(i,4) << " is either 0 or greater than the maximum y value." << std::endl;
                    error = 600;
                    return error;
                }
            }

            t.target[year].push_back(s);

            // push_back the number of individuals to the nb map
            if(paramsLand.patchModel){
                if((int)translocation_matrix_R(i,3)>0){
                    t.nb[year].push_back((int)translocation_matrix_R(i,3));
                }else{
                    Rcpp::Rcout << "ReadTranslocationR(): number of individuals to be translocated is 0 or a negative value." << std::endl;
                    error = 600;
                    return error;
                }
            } else {
                // cell-based
                if((int)translocation_matrix_R(i,5)>0){
                    t.nb[year].push_back((int)translocation_matrix_R(i,5));
                }else{
                    Rcpp::Rcout << "ReadTranslocationR(): number of individuals to be translocated is 0 or a negative value." << std::endl;
                    error = 600;
                    return error;
                }
            }

            // only if stage structured
            if(dem.stageStruct) {
                // push_back the minimal age of the individuals to the min_age map
                // the maximal age of the individuals to the max_age map
                // and the stage of the individuals to the stage map
                if(paramsLand.patchModel){
                    if ((int)translocation_matrix_R(i,4)>=0 || (int)translocation_matrix_R(i,4)==-9){
                        t.min_age[year].push_back((int)translocation_matrix_R(i,4));
                    } else{
                        Rcpp::Rcout << "ReadTranslocationR(): minimal age of the individuals to be translocated is a negative value which is not -9." << std::endl;
                        error = 600;
                        return error;
                    }
                    if ((int)translocation_matrix_R(i,5)>0 || (int)translocation_matrix_R(i,5)==-9){
                        t.max_age[year].push_back((int)translocation_matrix_R(i,5));
                    } else{
                        Rcpp::Rcout << "ReadTranslocationR(): maximal age of the individuals to be translocated is 0 or a negative value which is not -9." << std::endl;
                        error = 600;
                        return error;
                    }
                    if ((int)translocation_matrix_R(i,6)>=0 || (int)translocation_matrix_R(i,6)==-9){
                        t.stage[year].push_back((int)translocation_matrix_R(i,6));
                    } else{
                        Rcpp::Rcout << "ReadTranslocationR(): stage of the individuals to be translocated is a negative value which is not -9." << std::endl;
                        error = 600;
                        return error;
                    }
                } else {
                    if ((int)translocation_matrix_R(i,6)>=0 | (int)translocation_matrix_R(i,6)==-9){
                        t.min_age[year].push_back((int)translocation_matrix_R(i,6));
                    } else{
                        Rcpp::Rcout << "ReadTranslocationR(): minimal age of the individuals to be translocated is a negative value which is not -9." << std::endl;
                        error = 600;
                        return error;
                    }
                    if ((int)translocation_matrix_R(i,7)>0 | (int)translocation_matrix_R(i,7)==-9){
                        t.max_age[year].push_back((int)translocation_matrix_R(i,7));
                    } else{
                        Rcpp::Rcout << "ReadTranslocationR(): maximal age of the individuals to be translocated is 0 or a negative value which is not -9." << std::endl;
                        error = 600;
                        return error;
                    }
                    if ((int)translocation_matrix_R(i,8)>=0 | (int)translocation_matrix_R(i,8)==-9){
                        t.stage[year].push_back((int)translocation_matrix_R(i,8));
                    } else{
                        Rcpp::Rcout << "ReadTranslocationR(): stage of the individuals to be translocated is a negative value which is not -9." << std::endl;
                        error = 600;
                        return error;
                    }
                }
            } else{
                t.min_age[year].push_back(-9);
                t.max_age[year].push_back(-9);
                t.stage[year].push_back(-9);
            }

            if(dem.repType!=0) {
                // push_back the sex of the individuals to the sex map
                if(paramsLand.patchModel){
                    if ((int)translocation_matrix_R(i,7) == 0 | (int)translocation_matrix_R(i,7) == 1 | (int)translocation_matrix_R(i,7) == -9){
                        t.sex[year].push_back((int)translocation_matrix_R(i,7));
                    } else{
                        Rcpp::Rcout << "ReadTranslocation(): sex of the individuals to be translocated is not 0, 1 or -9." << std::endl;
                        error = 600;
                        return error;
                    }
                } else {
                    // cell-based
                    if ((int)translocation_matrix_R(i,9) == 0 | (int)translocation_matrix_R(i,9) == 1 | (int)translocation_matrix_R(i,9) == -9){
                        t.sex[year].push_back((int)translocation_matrix_R(i,9));
                    } else{
                        Rcpp::Rcout << "ReadTranslocation(): sex of the individuals to be translocated is not 0, 1 or -9." << std::endl;
                        error = 600;
                        return error;
                    }
                }
            } else{
                t.sex[year].push_back(-9);
            }
        }

        // check input
        // loop over t.source map and print out the content
        for (std::map<int, std::vector<locn>>::iterator it = t.source.begin(); it != t.source.end(); ++it) {
            Rcpp::Rcout << "ReadTranslocationR(): t.source[" << it->first << "]: ";
            for (int i = 0; i < it->second.size(); i++) {
                Rcpp::Rcout << it->second[i].x << " " << it->second[i].y << " ";
            }
            Rcpp::Rcout << std::endl;
        }

        // check input
        // loop over t.target map and print out the content
        for (std::map<int, std::vector<locn>>::iterator it = t.target.begin(); it != t.target.end(); ++it) {
            Rcpp::Rcout << "ReadTranslocationR(): t.target[" << it->first << "]: ";
            for (int i = 0; i < it->second.size(); i++) {
                Rcpp::Rcout << it->second[i].x << " " << it->second[i].y << " ";
            }
            Rcpp::Rcout << std::endl;
        }

        // check input
        // loop over t.nb map and print out the content
        for (std::map<int, std::vector<int>>::iterator it = t.nb.begin(); it != t.nb.end(); ++it) {
            Rcpp::Rcout << "ReadTranslocationR(): t.nb[" << it->first << "]: ";
            for (int i = 0; i < it->second.size(); i++) {
                Rcpp::Rcout << it->second[i] << " ";
            }
            Rcpp::Rcout << std::endl;
        }

        if(dem.stageStruct) {
            // check input
            // loop over t.min_age map and print out the content
            for (std::map<int, std::vector<int>>::iterator it = t.min_age.begin(); it != t.min_age.end(); ++it) {
                Rcpp::Rcout << "ReadTranslocationR(): t.min_age[" << it->first << "]: ";
                for (int i = 0; i < it->second.size(); i++) {
                    Rcpp::Rcout << it->second[i] << " ";
                }
                Rcpp::Rcout << std::endl;
            }

            // check input
            // loop over t.max_age map and print out the content
            for (std::map<int, std::vector<int>>::iterator it = t.max_age.begin(); it != t.max_age.end(); ++it) {
                Rcpp::Rcout << "ReadTranslocationR(): t.max_age[" << it->first << "]: ";
                for (int i = 0; i < it->second.size(); i++) {
                    Rcpp::Rcout << it->second[i] << " ";
                }
                Rcpp::Rcout << std::endl;
            }

            // check input
            // loop over t.stage map and print out the content
            for (std::map<int, std::vector<int>>::iterator it = t.stage.begin(); it != t.stage.end(); ++it) {
                Rcpp::Rcout << "ReadTranslocationR(): t.stage[" << it->first << "]: ";
                for (int i = 0; i < it->second.size(); i++) {
                    Rcpp::Rcout << it->second[i] << " ";
                }
                Rcpp::Rcout << std::endl;
            }
        }

        // only if sexual reproduction
        if(dem.repType!=0) {

            // check input
            // loop over t.stage map and print out the content
            for (std::map<int, std::vector<int>>::iterator it = t.sex.begin(); it != t.sex.end(); ++it) {
                Rcpp::Rcout << "ReadTranslocationR(): t.sex[" << it->first << "]: ";
                for (int i = 0; i < it->second.size(); i++) {
                    Rcpp::Rcout << it->second[i] << " ";
                }
                Rcpp::Rcout << std::endl;
            }
        }

        // set and update the translocation parameters
        pManagement->setTranslocationParams(t);

    }
    return error;
}
//---------------------------------------------------------------------------

Rcpp::List RunBatchR(int nSimuls, int nLandscapes, Rcpp::S4 ParMaster)
{
    int land_nr;
    int t0, t1, t00, t01;
    int read_error;
    bool params_ok;
    simParams sim = paramsSim->getSim();

    Rcpp::List list_outPop;
    Landscape* pLandscape = NULL; // pointer to landscape

#if RSDEBUG
    DEBUGLOG << endl;
    DEBUGLOG << "RunBatchR(): nSimuls=" << nSimuls << " nLandscapes=" << nLandscapes << endl;
    DEBUGLOG << "RunBatchR(): landtype=" << landtype << " maxNhab=" << maxNhab << endl;
#endif

    t0 = time(0);

    // int batch_line = 0;

    string name = paramsSim->getDir(2) + "Batch" + to_string(sim.batchNum) + "_RS_log.csv";
    if(rsLog.is_open()) {
        rsLog.close();
        rsLog.clear();
    }
    rsLog.open(name.c_str());
    if(!rsLog.is_open()) {
        Rcpp::Rcout << endl
                    << "Error - unable to open Batch" << sim.batchNum << "_RS_log.csv file - aborting batch run"
                    << endl;
        return Rcpp::List::create(Rcpp::Named("Errors") = -678);
    }
    rsLog << "Event,Number,Reps,Years,Time" << endl;
#if RSDEBUG
    rsLog << "WARNING,***** RSDEBUG mode is active *****,,," << endl;
#endif
#if RS_RCPP
    rsLog << "RNG SEED,,,," << RS_random_seed << endl;
#endif

    // loop over landscapes

    for(int j = 0; j < nLandscapes; j++) {

#if RSDEBUG
        DEBUGLOG << endl;
#endif
        // create new landscape
        if(pLandscape != NULL)
            delete pLandscape;
        pLandscape = new Landscape;
        bool landOK = true;

        t00 = time(0);

        landOK = ReadLandParamsR(pLandscape, ParMaster);
        //land_nr = ReadLandParamsR(pLandscape, ParMaster);
        land_nr = j; // TODO: ReadLandParamsR() is supposed to return land_nr; this is a temporary replacement; should I adapt this?

        if(!landOK) {
            // to get more precise error codes: use return values in ReadLandParamsR() similar to batch version
            // then it would be
            // if(land_nr<= 0){
            //  string msg = "Error code " + to_string(-land_nr)
            // + " returned from reading LandFile - aborting batch run";
            // cout << endl << msg << endl;
            rsLog << "Error reading landscape ASCII haeders - aborting" << endl;
            Rcpp::Rcout << "Error reading landscape ASCII haeders - aborting" << endl;
        } else {


#if RSDEBUG
            DEBUGLOG << endl << "RunBatchR(): j=" << j << " land_nr=" << land_nr << " landtype=" << landtype;
            if(landtype != 9)
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
                if(landtype == 9) {
                    paramsLand.generated = true;
                    paramsLand.nHab = 2;
                } else {
                    paramsLand.generated = false;
                    /*if(name_dynland == "NULL")
                     paramsLand.dynamic = false;
                     else
                     paramsLand.dynamic = true;*/
                }
                paramsLand.nHabMax = maxNhab;
                paramsLand.spDist = speciesdist;
                paramsLand.spResol = distresolution;
                pLandscape->setLandParams(paramsLand, sim.batchMode);

                if(landtype != 9) { // imported landscape
                    string hname = paramsSim->getDir(1) + name_landscape;
                    int landcode;
                    string cname;
                    if (name_costfile == "NULL" || name_costfile == "none") cname = "NULL";
                    else cname = paramsSim->getDir(1) + name_costfile;
                    if(paramsLand.patchModel) {
                        string pname = paramsSim->getDir(1) + name_patch;
                        landcode = pLandscape->readLandscape(0, hname, pname, cname);
                    } else
                        landcode = pLandscape->readLandscape(0, hname, " ", cname);
                    if(landcode != 0) {
                        rsLog << "Landscape," << land_nr << ",ERROR,CODE," << landcode << endl;
                        Rcpp::Rcout << endl << "Error reading landscape " << land_nr << " - aborting" << endl;
                        landOK = false;
                    }
                    if(paramsLand.dynamic) {
                        Rcpp::S4 LandParamsR("LandParams");
                        LandParamsR = Rcpp::as<Rcpp::S4>(ParMaster.slot("land"));
                        landcode = ReadDynLandR(pLandscape, LandParamsR);
                        if(landcode != 0) {
                            rsLog << "Landscape," << land_nr << ",ERROR,CODE," << landcode << endl;
                            Rcpp::Rcout << endl << "Error reading landscape " << land_nr << " - aborting" << endl;
                            landOK = false;
                        }
                    }
                    if(landtype == 0) {
                        pLandscape->updateHabitatIndices();
                    }
#if RSDEBUG
                    landParams tempLand = pLandscape->getLandParams();
                    DEBUGLOG << "RunBatchR(): j=" << j << " land_nr=" << land_nr << " landcode=" << landcode
                             << " nHab=" << tempLand.nHab << endl;
#endif

                    // species distribution

                    if(paramsLand.spDist) { // read initial species distribution
                        // WILL NEED TO BE CHANGED FOR MULTIPLE SPECIES ...
                        string distname = paramsSim->getDir(1) + name_sp_dist;
                        landcode = pLandscape->newDistribution(pSpecies, distname);
                        if(landcode == 0) {
                        } else {
                            rsLog << "Landscape," << land_nr << ",ERROR,CODE," << landcode << endl;
                            Rcpp::Rcout << endl
                                        << "Error reading initial distribution for landscape " << land_nr << " - aborting"
                                        << endl;
                            landOK = false;
                        }
                    }
                    paramsSim->setSim(sim);
#if RSDEBUG
                    DEBUGLOG << "RunBatchR(): j=" << j << " spDist=" << paramsLand.spDist << endl;
#endif

                    if(landOK) {
                        t01 = time(0);
                        rsLog << "Landscape," << land_nr << ",,," << t01 - t00 << endl;

                    } // end of landOK condition

                } // end of imported landscape
        }
        if(landOK) {

            // Open all other batch files and read header records

            // nSimuls is the total number of lines (simulations) in
            // the batch and is set in the control function
            string msgsim = "Simulation,";
            string msgerr = ",ERROR CODE,";
            string msgabt = ",simulation aborted";
            for(int i = 0; i < nSimuls; i++) { // this loop is useless at the moment since nSimuls is set to one in R entry function BatchMainR()
                t00 = time(0);
                params_ok = true;
                read_error = ReadParametersR(pLandscape, ParMaster);
                simParams sim = paramsSim->getSim();
                if(read_error) {
                    rsLog << msgsim << sim.simulation << msgerr << read_error << msgabt << endl;
                    params_ok = false;
                }
                if(stagestruct) {
                    ReadStageStructureR(ParMaster);
                }
                read_error = ReadEmigrationR(ParMaster);
                if(read_error) {
                    rsLog << msgsim << sim.simulation << msgerr << read_error << msgabt << endl;
                    params_ok = false;
                }
                read_error = ReadTransferR(pLandscape, ParMaster);
                if(read_error) {
                    rsLog << msgsim << sim.simulation << msgerr << read_error << msgabt << endl;
                    params_ok = false;
                }
                read_error = ReadSettlementR(ParMaster);
                if(read_error) {
                    rsLog << msgsim << sim.simulation << msgerr << read_error << msgabt << endl;
                    params_ok = false;
                }
                if(translocation){
                    read_error = ReadTranslocationR(pLandscape, ParMaster);
                    if(read_error) {
                        rsLog << msgsim << sim.simulation << msgerr << read_error << msgabt << endl;
                        params_ok = false;
                    }
                }
                if(params_ok) {
#if RSDEBUG
                    DebugGUI("RunBatchR(): simulation i=" + Int2Str(i));
#endif
                    // veraltet
                    // pSpecies->setNChromosomes(0);
                    // pSpecies->setTraits();
                }
                // veraltet
                // Rcpp::S4 GeneParamsR("GeneticsParams");
                // GeneParamsR = Rcpp::as<Rcpp::S4>(ParMaster.slot("gene"));
                // if (anyIndVar || Rcpp::as<int>(GeneParamsR.slot("Architecture")) == 1) {
                // 	read_error = ReadGeneticsR(GeneParamsR);
                // 	if(read_error) {
                // 		rsLog << msgsim << sim.simulation << msgerr << read_error << msgabt << endl;
                // 		params_ok = false;
                // 	}
                // } else {
                // 	// use default genetics parameters
                // 	// (by setting illegal values except for diploid)
                // 	genomeData g;
                // 	g.nLoci = -1;
                // 	g.probMutn = g.probCrossover = g.alleleSD = g.mutationSD = -1.0;
                // 	if(reproductn == 0)
                // 		g.diploid = false;
                // 	else
                // 		g.diploid = true;
                // 	g.neutralMarkers = g.trait1Chromosome = false;
                //
                // 	pSpecies->setGenomeData(g);
                // }
                read_error = ReadInitialisationR(pLandscape, ParMaster);
                if(read_error) {
                    rsLog << msgsim << sim.simulation << msgerr << read_error << msgabt << endl;
                    params_ok = false;
                }

                Rcpp::Rcout << "ReadInitialisationR() done." << endl;

                if (gHasGenetics) { // genetics need to be set and traits file need to be provided
                    Rcpp::S4 GeneParamsR("GeneticsParams");
                    GeneParamsR = Rcpp::as<Rcpp::S4>(ParMaster.slot("gene"));

                    read_error = ReadGeneticsR(GeneParamsR, pLandscape);
#if RSDEBUG
                    DEBUGLOG << "ReadGeneticsFile()" << endl;
#endif
                    if (read_error) {
                        rsLog << msgsim << sim.simulation << msgerr << read_error << msgabt << endl;
                        params_ok = false;
                    }
                    Rcpp::S4 TraitsParamsR("TraitsParams");
                    TraitsParamsR = Rcpp::as<Rcpp::S4>(GeneParamsR.slot("Traits"));
                    read_error = ReadTraitsR(TraitsParamsR);
#if RSDEBUG
                    DEBUGLOG << "ReadTraitsFile()" << endl;
#endif
                    if (read_error) {
                        rsLog << msgsim << sim.simulation << msgerr << read_error << msgabt << endl;
                        params_ok = false;
                    }
                }


                if(params_ok) {
                    simParams sim = paramsSim->getSim();

#if RSDEBUG
                    DEBUGLOG << endl
                             << "RunBatchR(): i=" << i << " simulation=" << sim.simulation
                    //<< " landFile=" << landFile
                      << " outRange=" << sim.outRange << " outIntRange=" << sim.outIntRange << endl;
#endif

                    Rcpp::Rcout << endl
                                << "Running simulation nr. " << to_string(sim.simulation)
                    //<< " on landscape no. " << Int2Str(land_nr)
                      << endl;


                    // for batch processing, include landscape number in parameter file name
                    OutParameters(pLandscape);

                    // run the model
                    list_outPop = RunModel(pLandscape, i);
#if RSDEBUG
                    // DEBUGLOG << endl << "RunBatchR(): real landscape, i = " << i
                    //	<< " simulation = " << sim.simulation << " landFile = " << landFile
                    //	<< endl;
#endif

                    t01 = time(0);
                    rsLog << msgsim << sim.simulation << "," << sim.reps << "," << sim.years << "," << t01 - t00
                          << endl;
                } // end of if (params_ok)
                else {
                    Rcpp::Rcout << endl << "Error in reading parameter file(s)... see RS log." << endl;
                }
            } // end of nSimuls for loop

            // close input files

            //		if (landtype != 9) {
            if (pLandscape != NULL) {
                delete pLandscape;
                pLandscape = NULL;
            }
            if(pManagement != NULL) {
                delete pManagement;
                pManagement = NULL;
            }
        } // end of landOK condition

    } // end of nLandscapes loop

    // Write performance data to log file
    t1 = time(0);
    rsLog << endl << "Batch,,,," << t1 - t0 << endl;

    if(rsLog.is_open()) {
        rsLog.close();
        rsLog.clear();
    }

    return list_outPop;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void setglobalvarsR(Rcpp::S4 control)
{

    // nSimuls = 2;                                   //read from parameterfile
    // nLandscapes = 1;                               //read from landfile, but not yet used in R-class RSparams!

    batchnum = Rcpp::as<int>(control.slot("batchnum"));
    patchmodel = Rcpp::as<int>(control.slot("patchmodel"));
    resolution = Rcpp::as<int>(control.slot("resolution"));
    landtype = Rcpp::as<int>(control.slot("landtype"));
    maxNhab = Rcpp::as<int>(control.slot("maxNhab"));
    speciesdist = Rcpp::as<int>(control.slot("speciesdist"));
    distresolution = Rcpp::as<int>(control.slot("distresolution"));
    reproductn = Rcpp::as<int>(control.slot("reproductn"));
    repseasons = Rcpp::as<int>(control.slot("repseasons"));
    stagestruct = Rcpp::as<int>(control.slot("stagestruct"));
    stages = Rcpp::as<int>(control.slot("stages"));
    gTransferType = Rcpp::as<int>(control.slot("transfer"));
    translocation = Rcpp::as<int>(control.slot("translocation"));
    gHasNeutralGenetics = Rcpp::as<int>(control.slot("neutralgenetics"));
    gHasGeneticLoad = Rcpp::as<int>(control.slot("geneticload"));
    // gHasGenetics should be true if gHasNeutralGenetics or gHasGeneticLoads is true
    gHasGenetics = gHasNeutralGenetics || gHasGeneticLoad;

#if RSDEBUG
    /*
     Function slotNames("slotNames");
     CharacterVector snames = slotNames(obj);
     for (int i = 0; i < snames.size(); i++) {
     SEXP slot = obj.slot(Rcpp::as<std::string>(snames[i]));
     // do something with slot
     }
     */
#endif
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

#if !RSWIN64
// check for UTF-16 encoding
string check_bom(string file)
{
    /*
     const char *UTF_16_BE_BOM = "\xFE\xFF";
     const char *UTF_16_LE_BOM = "\xFF\xFE";
     const char *UTF_8_BOM     = "\xEF\xBB\xBF";
     const char *UTF_32_BE_BOM = "\x00\x00\xFE\xFF";
     const char *UTF_32_LE_BOM = "\xFF\xFE\x00\x00";
     */
    string enc = "undef";

    ifstream infile(file, ios::in | ios::binary);
    if(!infile) {
        enc = "error";
        infile.clear();
    } else {
        char buffer[20];
        infile.read(buffer, 4);
        if(buffer[0] == '\xFE' && buffer[1] == '\xFF') { // UTF-16 LE
            enc = "utf16";
        }
        if(buffer[0] == '\xFF' && buffer[1] == '\xFE') { // UTF-16 BE
            enc = "utf16";
        }
        if(buffer[0] == '\xEF' && buffer[1] == '\xBB' && buffer[2] == '\xBF') { // UTF-8
            enc = "utf8";
        }
        if(buffer[0] == '\x00' && buffer[1] == '\x00' && buffer[2] == '\xFE' && buffer[3] == '\xFF') { // UTF-32 BE
            enc = "utf32";
        }
        if(buffer[0] == '\xFF' && buffer[1] == '\xFE' && buffer[2] == '\x00' && buffer[3] == '\x00') { // UTF-32 BE
            enc = "utf32";
        }
        /*
         if (size >= 3) {
         if (memcmp(data, UTF_8_BOM, 3) == 0)
         return "UTF-8";
         }
         if (size >= 4) {
         if (memcmp(data, UTF_32_LE_BOM, 4) == 0)
         return "UTF-32-LE";
         if (memcmp(data, UTF_32_BE_BOM, 4) == 0)
         return "UTF-32-BE";
         }
         if (size >= 2) {
         if (memcmp(data, UTF_16_LE_BOM, 2) == 0)
         return "UTF-16-LE";
         if (memcmp(data, UTF_16_BE_BOM, 2) == 0)
         return "UTF-16-BE";
         }*/
        infile.close();
    }
    return enc;
}
#endif

//---------------------------------------------------------------------------

rasterdata ParseRasterHead(string file)
{
    wifstream infile;
    rasterdata r;
    wstring header;
    int inint;

    r.ok = true;
    r.errors = r.ncols = r.nrows = r.cellsize = 0;
    r.xllcorner = r.yllcorner = 0.0;
    r.utf = false;

    // open file
#if RSWIN64
    infile.open(file.c_str());
#else
    infile.open(file, std::ios::binary);
#endif
    if (infile.is_open()) {
#if RSDEBUG
        DEBUGLOG << "Parsing raster file " << file << std::endl;
#endif
#if !RSWIN64
        // check BOM for UTF-16
        if(check_bom(file) == "utf16") {
            // apply BOM-sensitive UTF-16 facet
            infile.imbue(std::locale(infile.getloc(), new std::codecvt_utf16<wchar_t, 0x10ffff, std::consume_header>));
            r.utf = true;
        } else {
            r.utf = false;
        }
#endif
        // parse ASCII header
        infile >> header;
        if (!infile.good()) {
#if RSDEBUG
            DEBUGLOG << "ParseRasterHead(): failed to read " << file << std::endl;
#endif
            r.ok = false;
            r.errors = -211;
            infile.close();
            infile.clear();
            return r;
        }
        if (header != L"ncols" && header != L"NCOLS") r.errors++;
        infile >> r.ncols;

        infile >> header >> r.nrows;
        if (header != L"nrows" && header != L"NROWS") r.errors++;

        infile >> header >> r.xllcorner;
        if (header != L"xllcorner" && header != L"XLLCORNER") r.errors++;

        infile >> header >> r.yllcorner;
        if (header != L"yllcorner" && header != L"YLLCORNER") r.errors++;

        double tmpcellsize;
        infile >> header >> tmpcellsize;
        if (header != L"cellsize" && header != L"CELLSIZE") r.errors++;
        r.cellsize = (int) tmpcellsize;

        infile >> header >> inint;
        if (header != L"NODATA_value" && header != L"NODATA_VALUE") r.errors++;


        if (r.errors > 0) r.ok = false;

    } else {
        r.ok = false;
        r.errors = -111;
        //OpenErrorR("Raster file ", file);
#if RSDEBUG
        DEBUGLOG << "Raster file failed to open: " << file << std::endl;
#endif
    }
    infile.close();
    infile.clear();
    return r;
}

//----------------------------------------------------------------------------------------------

int ReadInitIndsFileR(int option, Landscape* pLandscape)
{
    landParams paramsLand = pLandscape->getLandParams();
    demogrParams dem = pSpecies->getDemogrParams();
    //stageParams sstruct = pSpecies->getStage();
    initParams init = paramsInit->getInit();

    string indsfile = paramsSim->getDir(1) + init.indsFile;
    wifstream initIndsFile;

    wstring header;
    string filetype = "InitIndsFile";
    string filename, ftype2, fname;

    int line, prevyear;
    //int year, sex, species, patchID, x, y, ninds, age, stage;

    int errors = 0;

    if(option == 0) { // open file, parse and read header and lines
        // open file
#if RSWIN64
        initIndsFile.open(indsfile.c_str());
#else
        initIndsFile.open(indsfile, std::ios::binary);
#endif
        if(!initIndsFile.is_open()) {
            OpenErrorR("Initial individuals file", indsfile);
#if RSDEBUG
            DEBUGLOG << "Initial individuals file failed to open: " << indsfile << std::endl;
#endif
            return -21;
        } else {
#if RSDEBUG
            DEBUGLOG << "Initial individuals file open to read" << std::endl;
#endif
#if !RSWIN64
            // check BOM for UTF-16
            if(check_bom(indsfile) == "utf16")
                // apply BOM-sensitive UTF-16 facet
                initIndsFile.imbue(std::locale(initIndsFile.getloc(), new std::codecvt_utf16<wchar_t, 0x10ffff, std::consume_header>));
#endif

            // Check right headers format
            initIndsFile >> header;
            if(!initIndsFile.good()) {
                Rcpp::Rcout << "Initial individuals failed to read" << std::endl;
#if RSDEBUG
                DEBUGLOG << "Initial individuals failed to read" << std::endl;
#endif
                return -211;
            }
            if(header != L"Year")
                errors++;
            initIndsFile >> header;
            if(header != L"Species")
                errors++;
            if(patchmodel) {
                initIndsFile >> header;
                if(header != L"PatchID")
                    errors++;
            } else {
                initIndsFile >> header;
                if(header != L"X")
                    errors++;
                initIndsFile >> header;
                if(header != L"Y")
                    errors++;
            }
            initIndsFile >> header;
            if(header != L"Ninds")
                errors++;
            if(reproductn > 0) {
                initIndsFile >> header;
                if(header != L"Sex")
                    errors++;
            }
            if(stagestruct) {
                initIndsFile >> header;
                if(header != L"Age")
                    errors++;
                initIndsFile >> header;
                if(header != L"Stage")
                    errors++;
            }
            // Report any errors in headers, and if so, terminate validation
            if(errors > 0) {
                FormatErrorR(filetype, errors);
                return -111;
            }

            paramsInit->resetInitInds();

            // Read data lines
            initInd iind;
            int ninds;
            int totinds = 0;

            line = 1;
            iind.year = prevyear = -98765;
            initIndsFile >> iind.year;

            while(iind.year != -98765) {
                // Year
                if(iind.year < 0) {
                    BatchErrorR(filetype, line, 19, "Year");
                    errors++;
                } else {
                    if(iind.year < prevyear) {
                        BatchErrorR(filetype, line, 2, "Year", "previous Year");
                        errors++;
                    }
                }
                prevyear = iind.year;

                // Species
                initIndsFile >> iind.species;
                if(iind.species != 0) {
                    BatchErrorR(filetype, line, 0, " ");
                    errors++;
                    Rcpp::Rcout << "Species must be 0" << endl;
                }

                // Patch | Coordinates
                if(paramsLand.patchModel) {
                    initIndsFile >> iind.patchID;
                    if(iind.patchID < 1) {
                        BatchErrorR(filetype, line, 11, "PatchID");
                        errors++;
                        iind.x = iind.y = 0;
                    }
                } else {
                    initIndsFile >> iind.x >> iind.y;
                    if(iind.x < 0 || iind.y < 0) {
                        BatchErrorR(filetype, line, 19, "X and Y");
                        errors++;
                        iind.patchID = 0;
                    }
                }

                // No of individuals
                initIndsFile >> ninds;
                if(ninds < 1) {
                    BatchErrorR(filetype, line, 11, "Ninds");
                    errors++;
                }

                // Sex
                if(dem.repType > 0){
                    initIndsFile >> iind.sex;
                    if(iind.sex < 0 || iind.sex > 1) {
                        BatchErrorR(filetype, line, 1, "Sex");
                        errors++;
                    }
                }
                else iind.sex = 0;

                // Stage
                if(dem.stageStruct) {
                    initIndsFile >> iind.age >> iind.stage;
                    if(iind.age < 1) {
                        BatchErrorR(filetype, line, 11, "Age");
                        errors++;
                    }
                    if(iind.stage < 1) {
                        BatchErrorR(filetype, line, 11, "Stage");
                        errors++;
                    }
                    if(iind.stage >= stages) {
                        BatchErrorR(filetype, line, 4, "Stage", "no. of stages");
                        errors++;
                    }
                } else {
                    iind.age = iind.stage = 0;
                }

                for(int i = 0; i < ninds; i++) {
                    totinds++;
                    paramsInit->addInitInd(iind);
                }

                line++;
                iind.year = -98765;				// finished current line
                if(!errors){					// check for format errors
                    if(!initIndsFile.eof())		// check for EOF to end loop
                        initIndsFile >> iind.year;		// read next year
                }
            } // end of while loop over lines
            if(!initIndsFile.eof()) {
                EOFerrorR(filetype);
                errors++;
            }

        } // end of "file is open"

        if(initIndsFile.is_open())
            initIndsFile.close();
        initIndsFile.clear();

        Rcpp::Rcout << "Initial individuals file OK:" << indsfile << std::endl;

        return errors; //totinds;

    } // end of option 0

    if(option == 9) { // close file
        if(initIndsFile.is_open()) {
            initIndsFile.close();
        }
        initIndsFile.clear();
        return 0;
    }
    return -1;
}

//---------------------------------------------------------------------------

void BatchErrorR(string filename, int line, int option, string fieldname)
{
    if(line == -999) { // message does not cite line number
        Rcpp::Rcout << "*** Error in " << filename << ": ";
    } else {
        Rcpp::Rcout << "*** Error in " << filename << " at line " << line << ": ";
    }
    switch(option) {
    case 0:
        break;
    case 1:
        Rcpp::Rcout << fieldname << " must be 0 or 1";
        break;
    case 2:
        Rcpp::Rcout << fieldname << " must be 0, 1 or 2";
        break;
    case 3:
        Rcpp::Rcout << fieldname << " must be 0, 1, 2 or 3";
        break;
    case 4:
        Rcpp::Rcout << fieldname << " must be from 0 to 4";
        break;
    case 5:
        Rcpp::Rcout << fieldname << " must be from 0 to 5";
        break;
    case 6:
        Rcpp::Rcout << fieldname << " must be from 0 to 6";
        break;
    case 7:
        Rcpp::Rcout << fieldname << " must be from 0 to 7";
        break;
    case 10:
        Rcpp::Rcout << fieldname << " must be greater than zero";
        break;
    case 11:
        Rcpp::Rcout << fieldname << " must be 1 or more";
        break;
    case 12:
        Rcpp::Rcout << fieldname << " must be 2 or more";
        break;
    case 13:
        Rcpp::Rcout << fieldname << " must be 3 or more";
        break;
    case 18:
        Rcpp::Rcout << fieldname << " must be greater than 1.0";
        break;
    case 19:
        Rcpp::Rcout << fieldname << " must be 0 or more";
        break;
    case 20:
        Rcpp::Rcout << fieldname << " must be between 0 and 1";
        break;
    case 21:
        Rcpp::Rcout << fieldname << " must be greater than 1";
        break;
    case 33:
        Rcpp::Rcout << fieldname << " must be 1, 2 or 3";
        break;
    case 44:
        Rcpp::Rcout << fieldname << " must be from 1 to 4";
        break;
    case 55:
        Rcpp::Rcout << fieldname << " must be from 1 to 5";
        break;
    case 66:
        Rcpp::Rcout << fieldname << " must be from 1 to 6";
        break;
    case 100:
        Rcpp::Rcout << fieldname << " must be between 0 and 100";
        break;
    case 111:
        Rcpp::Rcout << fieldname << " must match the first Simulation in ParameterFile";
        break;
    case 222:
        Rcpp::Rcout << "Simulation numbers must be sequential integers";
        break;
    case 333:
        Rcpp::Rcout << "No. of " << fieldname << " columns must equal max. no. of habitats (" << maxNhab
                    << ") and be sequentially numbered starting from 1";
        break;
    case 444:
        Rcpp::Rcout << "No. of " << fieldname << " columns must be one fewer than no. of stages, i.e. " << stages - 1
                    << ", and be sequentially numbered starting from 1";
        break;
    case 555:
        Rcpp::Rcout << "No. of " << fieldname << " columns must equal no. of stages, i.e. " << stages
                    << ", and be sequentially numbered starting from 0";
        break;
    case 666:
        Rcpp::Rcout << fieldname << " must be a unique positive integer";
        break;
    default:
        Rcpp::Rcout << "*** Unspecified error regarding parameter " << fieldname;
    }
    if(option != 0)
        Rcpp::Rcout << endl;
}

void BatchErrorR(string filename,int line,int option,string fieldname,string fieldname2)
{
    if (line == -999) { // message does not cite line number
        Rcpp::Rcout << "*** Error in " << filename << ": ";
    } else {
        Rcpp::Rcout << "*** Error in " << filename << " at line " << line <<": ";
    }
    switch (option) {
    case 0:
        break;
    case 1:
        Rcpp::Rcout << fieldname << " must be greater than " << fieldname2;
        break;
    case 2:
        Rcpp::Rcout << fieldname << " must be greater than or equal to " << fieldname2;
        break;
    case 3:
        Rcpp::Rcout << fieldname << " must be less than or equal to " << fieldname2;
        break;
    case 4:
        Rcpp::Rcout << fieldname << " must be less than " << fieldname2;
        break;
    default:
        Rcpp::Rcout << "*** Unspecified error regarding parameters " << fieldname
                    << " and " << fieldname2;
        ;
    }
    if (option != 0) Rcpp::Rcout << endl;
}

void ArchFormatErrorR(void)
{
    //batchlog << "*** Format error in ArchFile: case-sensitive parameter names "
    //	<< "must match the specification exactly" << endl;
    Rcpp::Rcout << "*** Format error in ArchFile:" << msgcase << msgmatch << endl;
}

void FormatErrorR(string filename, int errors)
{
    Rcpp::Rcout << "*** Format error in header line of ";
    if(errors == 0) {
        Rcpp::Rcout << filename << endl;
    } else {
        Rcpp::Rcout << filename << ": " << errors << " error";
        if(errors > 1)
            Rcpp::Rcout << "s";
        Rcpp::Rcout << " detected" << endl;
    }
}

void OpenErrorR(string ftype, string fname)
{
    Rcpp::Rcout << "*** Unable to open " << ftype << " " << fname << std::endl;
}

void EOFerrorR(string filename)
{
    Rcpp::Rcout << "*** Did not read to EOF in " << filename << std::endl;
}

void StreamErrorR(string filename)
{
    Rcpp::Rcout << "*** Corrupted file stream in " << filename << std::endl << "Too few entries? Unsuppoerted file encoding? (You might try to use a different one, like UTF-8.)" << std::endl;
#if RSDEBUG
    DEBUGLOG << "Corrupted file stream in " << filename << std::endl;
#endif
}

//---------------------------------------------------------------------------
