#---------------------------------------------------------------------------
#
#	Copyright (C) 2020-2022 Anne-Kathleen Malchow, Greta Bocedi, Stephen C.F. Palmer, Justin M.J. Travis, Damaris Zurell
#
#	This file is part of RangeShiftR.
#
#	RangeShiftR is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
#
#	RangeShiftR is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with RangeShiftR. If not, see <https://www.gnu.org/licenses/>.
#
#----------------------------------------------------------------------------



#### CLASS LANDPARAMS

# from RS 'Land' file
# can take one of two forms: 'ArtificialLandscape' or 'ImportedLandscape'

LandParams <- setClass("LandParams", slots = c(LandNum = "integer_OR_numeric")
                                    , prototype = list(LandNum = 1L)
                        )
    # landscape number must be unique

setValidity("LandParams", function(object) {
    msg <- NULL
    if (anyNA(object@LandNum) || length(object@LandNum)!=1) {
        msg <- c(msg, "LandNum must be set and of length 1!")
    }
    else {
        if (object@LandNum < 1) {
            msg <- c(msg, "Landscape number must be positive.")
        }
    }
    if (is.null(msg)) TRUE else msg}
)
setMethod("initialize", "LandParams", function(.Object, ...) {
    .Object <- callNextMethod()
    .Object}
)



### CLASS ARTIFICIALLANDSCAPE

#' Create an Artificial Landscape
#'
#' An artificial landscape can be generated with a random or fractal spatial structure, and may be using binary or continuous habitat values to characterise each cell.
#'
#' @usage ArtificialLandscape(propSuit = 0.5, K_or_DensDep = 10, Resolution = 100, dimX = 65, dimY = 65,
#'                     fractal = FALSE, hurst, continuous = FALSE, minPct, maxPct)
#' @param propSuit Proportion of suitable habitat cells, defaults to \eqn{0.5}.
#' @param K_or_DensDep determines the demographic density dependence of the modelled species and is given in units of the number of individuals per hectare (defaults to \eqn{10}).
#' Its precise interpretation depends on the structure of the species' \code{\link[RangeShiftR]{Demography}}:\cr
#' If combined with a \code{\link[RangeShiftR]{StageStructure}}d model, \code{K_or_DensDep} will be used as the \emph{strength of demographic density dependence} \ifelse{html}{\out{b<sup>-1</sup>}}{\eqn{1/b}}.
#' If combined with a non-structured model, \code{K_or_DensDep} will be interpreted as \emph{limiting carrying capacity} \eqn{K}.\cr
#' @param Resolution Cell size in meters, defaults to \eqn{100}. (integer)
#' @param dimX,dimY Number of cells along the x- and y-axis, both default to \eqn{65}. (integer) \cr
#' If \code{fractal=TRUE}, \code{dimX} and \code{dimY} must be powers of \eqn{2} plus \eqn{1} (\eqn{2^n+1}) and \code{dimX} \eqn{\le} \code{dimY} .
#' @param fractal If \code{FALSE} (default), a random landscape is generated. Each cell has a certain probability of being a habitat.\cr
#' If \code{TRUE}, a fractal landscape is generated. These landscapes contain a greater structure than random landscapes but less than a completely deterministic one.
#' @param hurst Required if \code{fractal=TRUE}: Hurst exponent. Can be any number in the open interval \eqn{(0,1)}.
#' @param continuous Use continuous or binary habitat values to describe each cell?\cr If \code{FALSE} (default), the resulting landscape
#' is binary, with the two values \emph{Habitat} and \emph{Matrix} (i.e. Non-habitat).\cr If \code{TRUE}, each cell is given a continuous value
#' between \eqn{0} and \eqn{1}, which describes the percentage of habitat cover or habitat quality.
#' within a cell. The effective \code{K_or_DensDep} of that cell is calculated as the respective fraction of the value of \eqn{K_or_DensDep}.
#' @param minPct,maxPct Required if \code{continuous=TRUE}: Minimum and Maximum percentage of habitat cover within a cell. Can be any number in the open interval \eqn{(0,1)},
#' \code{maxPct} may be exactly \eqn{1}. \code{minPct} must be smaller than \code{maxPct}.
#' @details For theoretical studies which might be related to fundamental questions in eco-evolutionary dynamics or strategic questions concerning
#' conservation ecology, it is often desirable to use artificial landscapes.\cr
#' Fractal landscapes are characterised by possessing greater structure than a completely random landscape, but less than a
#' completely deterministic one \insertCite{with1997application}{RangeShiftR} – but note that the spatial structure of landscapes fragmented
#' by human activities is often not fractal in nature and, depending upon the research question, other landscape generators may be more
#' appropriate \insertCite{@see review by @pe2013simple}{RangeShiftR}.
#'
#' Internally generated artificial landscapes may not be dynamic.
#'
#' The fractal landscape generator implemented in \emph{RangeShiftR} uses the midpoint displacemeant algorithm \insertCite{saupe1988algorithms}{RangeShiftR}.
#' The Hurst exponent, often referred to as \eqn{H}, describes the degree of spatial autocorrelation of the landscape configuration.
#' Values close to \eqn{0} represent a low autocorrelation but the generated landscapes still aren't completely spatially independent.
#' On the contrary, values close to \eqn{1} represent high autocorrelation, i.e. high habitat aggregation.\cr
#'
#' Note that more complex algorithms are available for providing fractals where setting \eqn{H = 0.0} results in no spatial autocorrelation
#' \insertCite{@see @chipperfield2011updated}{RangeShiftR}. For applications where the embedded algorithm is not sufficient, we recommend
#' to import landscapes generated by these alternative algorithms.
#'
#'
#' @references
#'         \insertAllCited{}
#' @return A parameter object of class "ArtificialLandscape"
#' @author Anne-Kathleen Malchow
#' @name ArtificialLandscape
#' @export ArtificialLandscape
ArtificialLandscape <- setClass("ArtificialLandscape", slots = c(propSuit = "numeric",
                                                                 K_or_DensDep = "integer_OR_numeric",
                                                                 Resolution = "integer_OR_numeric",
                                                                 dimX = "integer_OR_numeric",
                                                                 dimY = "integer_OR_numeric",
                                                                 fractal = "logical",
                                                                 hurst = "numeric",
                                                                 continuous = "logical",
                                                                 minPct = "numeric",
                                                                 maxPct = "numeric")
                           , prototype = list(propSuit = 0.5,
                                              K_or_DensDep = 10L,
                                              Resolution = 100L,
                                              dimX = 65L,
                                              dimY = 65L,
                                              fractal = FALSE,
                                              #hurst,
                                              continuous = FALSE)
                                              #minPct,
                                              #maxPct,
                           , contains = "LandParams")

setValidity("ArtificialLandscape", function(object) {
    msg <- NULL
    if (anyNA(object@propSuit) || length(object@propSuit)!=1) {
        msg <- c(msg, "Proportion of suitable habitat must be set and of length 1!")
    }
    else {
        if (object@propSuit<0 || object@propSuit>1) {
            msg <- c(msg, "Proportion of suitable habitat must be in the interval [0,1].")
        }
    }
    if (anyNA(object@K_or_DensDep) || length(object@K_or_DensDep)==0) {
        msg <- c(msg, "K_or_DensDep must be set!")
    }
    else {
        if (object@K_or_DensDep<0) {
            msg <- c(msg, "K_or_DensDep must not be smaller than 0.")
        }
    }
    if (anyNA(object@Resolution) || length(object@Resolution)!=1) {
        msg <- c(msg, "Resolution of landscape must be set and of length 1!")
    }
    else {
        if (object@Resolution < 1) {
            msg <- c(msg, "Resolution of landscape must be positive.")
        }
    }
    if (anyNA(object@dimX) || length(object@dimX)!=1) {
        msg <- c(msg, "dimX must be set!")
    }
    if (anyNA(object@dimY) || length(object@dimY)!=1) {
        msg <- c(msg, "dimY must be set!")
    }
    if (anyNA(object@fractal) || length(object@fractal)!=1) {
        msg <- c(msg, "fractal must be set!")
    }
    if (anyNA(object@continuous) || length(object@continuous)!=1) {
        msg <- c(msg, "continuous must be set!")
    }
    if (is.null(msg)) {
        if (object@fractal) {
            if (object@dimX < 3) {
                msg <- c(msg, "Number of cells in any direction must be at least 3 for a fractal landscape.")
            }
            else {
                if (!isPowerOf2(object@dimX-1)) {
                    msg <- c(msg, "Number of cells in any direction must be a power of 2 plus 1 for a fractal landscape.")
                }
            }
            if (object@dimY < 3) {
                msg <- c(msg, "Number of cells in any direction must be at least 3 for a fractal landscape.")
            }
            else {
                if (!isPowerOf2(object@dimY-1)) {
                    msg <- c(msg, "Number of cells in any direction must be a power of 2 plus 1 for a fractal landscape.")
                }
            }
            if (object@dimY < object@dimX) {
                msg <- c(msg, "Y-dimension may not be less than X-dimension for a fractal landscape.")
            }
            if (anyNA(object@hurst) || length(object@hurst)!=1) {
                msg <- c(msg, "Hurst exponent must be set for a fractal landscape.")
            }
            else {
                if (object@hurst<0 || object@hurst>1) {
                    msg <- c(msg, "Hurst exponent must be in the open interval (0,1).")
                }
                if (object@hurst==0 || object@hurst==1) {
                    msg <- c(msg, "Hurst exponent must be in the open interval (0,1). It can't be exactly 0 or 1.")
                }
            }
        }
        else { # non-fractal landscape
            if (object@dimX < 1) {
                msg <- c(msg, "Number of cells in any direction must be at least 1 for a non-fractal landscape.")
            }
            if (object@dimY < 1) {
                msg <- c(msg, "Number of cells in any direction must be at least 1 for a non-fractal landscape.")
            }
        }
        if(object@continuous) {
            if (anyNA(object@minPct) || length(object@minPct)!=1) {
                minP.ok = FALSE
                msg <- c(msg, "Minimum habitat percentage must be set for a continuous landscape.")
            }
            else {
                if (object@minPct<=0 || object@minPct>=1) {
                    minP.ok = FALSE
                    msg <- c(msg, "Minimum habitat percentage must be in the open interval (0,1).")
                }
                else {minP.ok = TRUE}
            }
            if (anyNA(object@maxPct) || length(object@maxPct)!=1) {
                maxP.ok = FALSE
                msg <- c(msg, "Maximum habitat percentage must be set for a continuous landscape.")
            }
            else {
                if (object@maxPct<=0 || object@maxPct>1) {
                    maxP.ok = FALSE
                    msg <- c(msg, "Maximum habitat percentage must be in the half-open interval (0,1].")
                }
                else {maxP.ok = TRUE}
            }
            if (minP.ok && maxP.ok) {
                if (object@maxPct <= object@minPct) {
                    msg <- c(msg, "Maximum habitat percentage may not be less than Minimum habitat percentage.")
                }
            }
        }
    }
    if (is.null(msg)) TRUE else msg}
)
setMethod("initialize", "ArtificialLandscape", function(.Object,...) {
    this_func = "ArtificialLandscape(): "
    args <- list(...)
    .Object <- callNextMethod()
    if ( length(args) == 0 ) {
        validObject(.Object)
    }
    if (!.Object@fractal) {
        .Object@hurst = -9L
        if (!is.null(args$hurst)) {
            warning(this_func, "Hurst exponent", warn_msg_ignored, "for a non-fractal landscape.", call. = FALSE)
       }
    }
    if(!.Object@continuous) {
        .Object@minPct = -9L
        if (!is.null(args$minPct)) {
            warning(this_func, "Minimum habitat percentage", warn_msg_ignored, "for a discrete landscape.", call. = FALSE)
        }
        .Object@maxPct = -9L
        if (!is.null(args$maxPct)) {
            warning(this_func, "Maximum habitat percentage", warn_msg_ignored, "for a discrete landscape.", call. = FALSE)
        }
    }
    .Object}
)
setMethod("show", "ArtificialLandscape", function(object){
    cat(" Artificial landscape: ")
    #cat(" Artificial landscape #", object@LandNum, ": ")
    if(object@fractal) {
        cat("fractal structure, ")
    }
    else {
        cat("random structure, ")
    }
    if(object@continuous) {
        cat("continuous habitat cover/quality\n")
    }
    else {
        cat("binary habitat/matrix code\n")
    }
    cat ("   Size         :", paste(object@dimX), "x", paste(object@dimY), "cells\n")
    cat ("   Resolution   :      ", paste(object@Resolution) , "meters\n")
    cat ("   Proportion of suitable habitat:", object@propSuit,"\n")
    cat ("   K or 1/b     :", object@K_or_DensDep,"\n")
    if(object@fractal) {
        cat ("   Hurst exponent    : H =", object@hurst, "\n")
    }
    if(object@continuous) {
        cat ("   Minimum habitat percentage: minPct =", object@minPct, "\n")
        cat ("   Maximum habitat percentage: maxPct =", object@maxPct, "\n")
    }
})


### CLASS IMPORTEDLANDSCAPE

#' Import a Landscape from file
#'
#' @description Provide the filename(s) (\code{LandscapeFile}) or the matrices (\code{LandscapeMatrix}) of the map(s) to be imported, their
#' resolution, origin coordinates in case you use matrices and, if applicable, the number of habitat codes (\code{Nhabitats})
#' as well as their respective demographic density dependence (\code{K_or_DensDep}).
#'
#' For a dynamic landscape, the year in which each landscape is loaded has to be provided.
#'
#' Other, optional input maps are:\cr
#' - Patch map(s) to define habitat patches - as filenames or matrices,\cr
#' - SMS cost map(s) to define landscape resistance to - as filenames or matrices,\cr
#' - a distribution map to define an initial species distribution - as filenames or matrices.
#'
#' @usage ImportedLandscape(LandscapeFile,
#'                  LandscapeMatrix,
#'                  Resolution = 100,
#'                  OriginCoords = c(0,0),
#'                  HabPercent = FALSE,
#'                  Nhabitats, K_or_DensDep = 10,
#'                  PatchFile =  NULL,
#'                  PatchMatrix =  list(),
#'                  CostsFile =  NULL,
#'                  CostsMatrix = list(),
#'                  DynamicLandYears = 0,
#'                  SpDistFile = NULL,
#'                  SpDistMatrix = list(),
#'                  SpDistResolution,
#'                  demogScaleLayers, nrDemogScaleLayers)
#' @param LandscapeFile Filename(s) of the landscape habitat map(s) which shall be imported from the Inputs-folder. See the Details for information on the required format.
#' @param LandscapeMatrix List of matrices of the landscape habitat raster(s); contains each cells' habitat suitability or land cover index.
#' @param Resolution Cell size in meters, defaults to \eqn{100}. (integer)
#' @param OriginCoords X- and Y-coordinates of the map origin given in meters as a vector of length 2. Defaults to \eqn{0,0}.
#' @param HabPercent If \code{FALSE} (default), unique integer habitat codes are expected in the imported map to characterise the habitat of each cell. This requires to set \code{Nhabitats}. \cr
#' If \code{TRUE}, continuous values are expected, ranging from \eqn{0.0} to \eqn{100.0}, that represent percentages of habitat cover or quality.\cr
#' Make sure your imported landscape file uses the specified standard (see Details below).
#' @param Nhabitats Required if \code{HabPercent=FALSE}: Number of habitats in the imported landscape if unique integer habitat codes are used. (integer)
#' @param K_or_DensDep determines the demographic density dependence of the modelled species and is given in units of the number of individuals per hectare (defaults to \eqn{10}).
#' Its precise interpretation depends on the structure of the species' \code{\link[RangeShiftR]{Demography}}:\cr
#' If combined with a \code{\link[RangeShiftR]{StageStructure}}d model, \code{K_or_DensDep} will be used as the \emph{strength of demographic density dependence} \ifelse{html}{\out{b<sup>-1</sup>}}{\eqn{1/b}}.
#' If combined with a non-structured model, \code{K_or_DensDep} will be interpreted as \emph{limiting carrying capacity} \eqn{K}.\cr
#' The expected format:\cr
#' If \code{HabPercent=FALSE}: a vector of length \code{Nhabitats}, specifying the respective \code{K_or_DensDep} for every habitat code.\cr
#' If \code{HabPercent=TRUE}: \code{K_or_DensDep} is interpreted as the maximum \code{K_or_DensDep} reached in cells with \eqn{100}\% habitat. All other cells hold the respective fraction of \eqn{K_or_DensDep}.
#' @param PatchFile Filename(s) of the patch map(s) which shall be imported, default is \code{NULL}.
#' @param PatchMatrix List of matrices of the patch raster(s); contains each cells' patch index.
#' @param CostsFile Filename(s) of the SMS cost map(s) which shall be imported, default is \code{NULL}.
#' @param CostsMatrix List of matrices of the SMS cost raster(s); contains each cells' SMS cost.
#' @param DynamicLandYears Integer vector indicating the years of landscape changes. For a non-dynamic landscape its only entry is \eqn{0} (default).
#' For a dynamic landscape, \code{DynamicLandYears} lists the years in which the corresponding habitat maps in \code{LandscapeFile} and - if applicable - their respective patch and/or costs
#' maps (in \code{PatchFile},\code{CostsFile}) are loaded and used in the simulation. More details below.
#' @param SpDistFile Filename of the species initial distribution map which shall be imported (*.txt). Default is \code{NULL}.
#' @param SpDistMatrix List of one matrix containing the species' initial distribution raster.
#' @param SpDistResolution Required if \code{SpDistFile} is given: Cell size of the distribution map in meters. (integer) Must be an integer multiple of the landscape resolution.
#' @param demogScaleLayers List of arrays that describe additional landscape layers which can be used to locally scale certain demographic rates and thus allow them to vary spatially.
#' The list must contain equally sized 3D-arrays, one for each element in \code{DynamicLandYears}, which are interpreted as stacked layers.
#' The arrays' first two dimensions correspond to the x- and y-dimensions of the maps in \code{LandscapeMatrix} and must match in resolution and offset;
#' the third dimension indexes the various layers (of which there are \code{nrDemogScaleLayers}).
#' Can only be used in combination with habitat quality maps, i.e. when \code{HabPercent=TRUE}. Contain percentage values ranging from \eqn{0.0} to \eqn{100.0}.
#' @param nrDemogScaleLayers number of additional landscape layers for spatial demographic scaling
#'
#'
#' @details \emph{RangeShiftR} requires all landscape related inputs to be either text files in ESRI ASCII raster format or matrices.
#'
#' Text files must have the following six header lines:
#' \tabular{ll}{\code{ncols} \tab Number of columns \cr
#' \code{nrows} \tab Number of rows \cr
#' \code{xllcorner} \tab x-coordinate (longitude) of the lower-left corner \cr
#' \code{yllcorner} \tab y-coordinate (latitude) of the lower-left corner \cr
#' \code{cellsize} \tab Resolution (in meters) \cr
#' \code{NODATA_value} \tab Value for cells having missing data (usually -9999) }
#'
#' The rest of the file is a grid containing a value for each cell, one line per row.
#'
#' If you want to use matrices for landscape inputs, you must provide the origin coordinates \code{OriginCoords}.
#' \code{OriginCoords} is the map origin given as a vector of length 2 where the first entry is the x-coordinate (longitude) of the lower-left corner and
#' the second entry is the y-coordinate (latitude) of the lower-left corner.\cr
#'
#' \code{RangeShiftR} can use two different types of habitat maps:\cr
#' \itemize{
#'     \item \emph{Raster with habitat codes} (\code{HabPercent=FALSE})\cr In this option each habitat or land-cover type has a unique integer code. Each cell in the file contains a single habitat code and \eqn{100} percent coverage is assumed for the cell. The landscape is therefore composed of discrete habitat cells. The codes are required to be sequential integers starting from \eqn{1} and ranging to \code{Nhabitats}.\cr
#'     \item \emph{Raster with habitat quality} (\code{HabPercent=TRUE})\cr Each cell in the landscape is assigned a continuous percentage value between \eqn{0.0} and \eqn{100.0} of the maximum \code{K_or_DensDep}. There are no explicit habitat or land-cover types. This allows integrating different methods for calculating the habitat suitability for a given species. For example, qualities can result from
#'  different methods of suitability modelling, which incorporate multiple variables like habitat types, elevation, climate, etc. In the current version of the program, a straight-line relationship between the habitat quality and the actual density-dependence of the population dynamics (i.e. \eqn{K} or \eqn{1/b} in a non-structured or structured population, respectively) is assumed.
#'  Therefore, the quality should be scaled accordingly in case of a curvilinear relationship.
#'  }
#'
#' \emph{Patch map} \cr
#' The simulation can be run as a \emph{patch-based model} on the same habitat map described above. An additional file/matrix must be provided through \code{PatchFile} / \code{PatchMatrix}: a raster map of the same landscape, where
#' each cell contains the ID number of the patch to which it belongs. Each patch must have a unique positive integer ID. The ID of every cell that does not belong to a patch (i.e. non-habitat/matrix) must be zero.
#' Note that a single patch is the unit at which the density dependence in the population dynamics acts. Therefore, a patch can be discontinuous, i.e. it can contain cells that do not belong to the patch if they
#' are assumed not to affect the dynamics, or on the other hand, patch cells that are not physically contiguous to the rest of the patch cells.
#'
#' \emph{Costs layer} \cr
#' Only used if the simulation is run with \code{\link[RangeShiftR]{SMS}} as its transfer module and the landscapes resistance to movement is given via a costs raster map (see argument \code{Costs} in \code{SMS()}).
#' In this case, the specified raster \code{CostsFile} / \code{CostsMatrix} has to match the landscape raster in extent, coordinates and resolution, and each cell contains a cost value, with the minimal possible cost being \eqn{1}.
#' Importing a cost layer is the only option when the landscape comprises habitat coverage or quality.
#'
#' \emph{Initial distribution} \cr
#' A \emph{species distribution map} can be overlaid on top of the habitat map and can be used to define an initial distribution. The map or raster matrix is provided through \code{SpDistFile} / \code{SpDistMatrix} and must be aligned with the landscape map. The extent of the map does not have to be necessarily
#' the same as the landscape. The resolution can be the same or coarser, provided that it is a multiple of the landscape resolution. For example, if the landscape cell size is \eqn{250m}, the species distribution can be at the resolution of \eqn{250m}, \eqn{500m}, \eqn{750m}, \eqn{1000m} etc.
#' Each cell of the species distribution map must contain either \eqn{0} (species absent or not recorded) or \eqn{1} (species present).
#'
#'  \emph{Demographic scaling layers} \cr
#' A number of additional landscape layers can be provided in \code{demogScaleLayers} to locally scale certain demographic rates and thus allow them to vary spatially.
#' This can only be used in combination with habitat quality maps, i.e. when \code{HabPercent=TRUE}, and with a stage-structured population model.
#' The additional layers contain percentage values ranging from \eqn{0.0} to \eqn{100.0}. Chosen demographic rates for a specified stage and sex can be mapped to one of these scaling layers using the
#' parameters \code{FecLayer}, \code{DevLayer}, and \code{SurvLayer} in \code{\link[RangeShiftR]{StageStructure}}. If a demographic rate varies spatially, its value in the transition matrix \code{TransMatrix})
#' is now interpreted as a maximum value, while the realised local value in each cell or patch is determined as this maximum value scaled by the percentage given in the respective mapped scaling layer.
#' For a patch-based landscape, the scaling percentage of a patch is given by the average percentage of its constituent cells.
#' Potential density-dependence mediated by the strength 1/b still takes effect also for spatially-varying demographic rates. The respective base values φ_0, σ_0 or γ_0 are then replaced by their locally scaled values.
#'
#' \emph{Dynamic landscapes} \cr
#' An imported landscape may be dynamic, i.e. the attributes of cells (either habitat class or quality index) and its patch number (if the model is patch-based) may be changed at specified years during the course of
#' a simulation. Note that any landscape change occurs at the start of the year, i.e. before the first/only reproductive season. In a patch-based model, the shape of patches may change, patches may
#' be removed and new patches may be created where there was previously inter-patch matrix. Thus some populations may be extirpated (in a non-structured population, all individuals die; in a stage-structured population,
#' all individuals either die or have an immediate opportunity to disperse), and new populations may arise from colonisation of newly suitable areas.
#'
#' However, there are certain restrictions. Any part of the original landscape which was a ‘no-data’ region (e.g. the sea or land beyond a study area boundary) must remain in that state for the whole simulation.
#' The identity of patches is not cross-checked between changes, and care must therefore be taken to ensure consistency; otherwise, a patch (and its resident population) can jump to a distant location or be split into
#' two or more disjunct parts, with unpredictable and possibly weird consequences. It is legitimate for a patch to be split into two or more separate patches (e.g. by construction of a motorway or some other barrier),
#' but any existing population will remain with the part (if any) which retains the original patch number, and populations within the other parts (having a new patch number) must arise through colonisation.
#' Possible ways to work around this restriction include:
#' \itemize{
#'    \item Assign to all post-change parts of the original patch a new, unique patch number and specify that dispersal is allowed after population destruction (which is possible only for a structured population), in which
#' case some colonisation of the new patches should occur. Note that the connectivity matrix will be misleading in such cases, as every successful ‘disperser’ will appear to have moved from patch N to patch M (where M
#' is the new patch number).
#'    \item Instead of a single original patch, define two (or more) distinct but adjacent patches in the original landscape, so that they each retain their own populations when they become separated by the landscape change.
#' }
#'
#' A dynamic landscape can be specified using the slots \code{LandscapeFile} / \code{LandscapeMatrix} (, \code{PatchFile} / \code{PatchMatrix}, \code{CostsFile}) / \code{CostsMatrix} and \code{DynamicLandYears}.
#' \code{LandscapeFile} (and \code{PatchFile}, \code{CostsFile}) take a character vector with the file names of the maps to be loaded.
#' \code{LandscapeMatrix} (and \code{PatchMatrix}, \code{CostsMatrix}) take a list with the matrices representing the raster maps to be used.
#'
#' Please be aware, that all landscape related input must be either matrices or file names. You cannot mix both types of inputs.
#'
#' All provided maps must agree in resolution, extent and origin. \code{DynamicLandYears} is a number vector that contains the years, in which these landscapes shall be loaded; it must have the same ordering so
#' that years and maps can be matched. If a specific map is used multiple times, it must be listed each time nevertheless.
#'
#' @return A parameter object of class ImportedLandscape
#' @author Anne-Kathleen Malchow, Jette Reeg
#' @name ImportedLandscape
#' @export ImportedLandscape
ImportedLandscape <- setClass("ImportedLandscape", slots = c(LandscapeFile = "character",
                                                             LandscapeMatrix = "list",
                                                             OriginCoords = "ANY",
                                                             Resolution = "integer_OR_numeric",
                                                             HabPercent = "logical",
                                                             Nhabitats = "integer_OR_numeric", # not used in RS anymore. In R is used to define maxNhab in ControlParams
                                                             K_or_DensDep = "integer_OR_numeric",
                                                             PatchFile = "ANY",          # sets the patchmodel -switch in class ControlParams when added
                                                             PatchMatrix = "list",
                                                             CostsFile = "ANY",
                                                             CostsMatrix = "list",
                                                             SpDistFile = "ANY",         # sets the speciesdist -switch in class ControlParams when added
                                                             SpDistMatrix = "list",
                                                             SpDistResolution = "integer_OR_numeric",
                                                             DynamicLandYears = "integer_OR_numeric",
                                                             nrDemogScaleLayers = "integer",
                                                             demogScaleLayers = "list") #= "data.frame")
                              , prototype = list(#LandscapeFile,
                                                LandscapeMatrix = list(),
                                                 Resolution = 100L,
                                                 HabPercent = FALSE,
                                                 OriginCoords = NULL,
                                                 #Nhabitats,
                                                 K_or_DensDep = 10L,
                                                 PatchFile = NULL,
                                                 PatchMatrix = list(),
                                                 CostsFile = NULL,
                                                 CostsMatrix = list(),
                                                 SpDistFile = NULL,
                                                 SpDistMatrix = list(),
                                                 #SpDistResolution,
                                                 DynamicLandYears = 0L,
                                                 nrDemogScaleLayers = 0L,
                                                 demogScaleLayers= list()) #= data.frame())
                              , contains = "LandParams")

setValidity("ImportedLandscape", function(object) {
    msg <- NULL

    land_ncol <- 0
    land_nrow <- 0
    if ((anyNA(object@LandscapeFile) || length(object@LandscapeFile)==0) && # No LandscapeFile provided
        (anyNA(object@LandscapeMatrix) || length(object@LandscapeMatrix)==0)) { # No Landscape Matrix provided
        msg <- c(msg, "No landscape filename nor matrix was given.")
    }
    else {
        if (length(object@LandscapeMatrix) > 0) { # matrix input
            if(any(sapply(object@LandscapeMatrix, class)[1,] != "matrix")){
                msg <- c(msg, "All elements of LandscapeMatrix list must be of class matrix.")
            }
            else{
                land_ncol <- ncol(object@LandscapeMatrix[[1]])
                land_nrow <- nrow(object@LandscapeMatrix[[1]])
                if( (any(sapply(object@LandscapeMatrix, ncol) != land_ncol)) || (any(sapply(object@LandscapeMatrix, nrow) != land_nrow)) ){
                    msg <- c(msg, "All elements of LandscapeMatrix list must have the same ncol and nrow.")
                }
            }
        }
    }

    if (anyNA(object@Resolution) || length(object@Resolution)!=1) {
        msg <- c(msg, "Resolution of landscape must be given and of length 1!")
    }
    else {
        if (object@Resolution < 1) {
            msg <- c(msg, "Resolution of landscape must be positive.")
        }
    }
    if (anyNA(object@HabPercent) || length(object@HabPercent)!=1) {
        msg <- c(msg, "HabPercent must be set!")
    }
    if (length(object@LandscapeMatrix) > 0) {
        if (!is.numeric(object@OriginCoords) || anyNA(object@OriginCoords) || length(object@OriginCoords)!=2) {
            msg <- c(msg, "Origin coordinates must be set and of length 2!")
        } else{
            if ( any(object@OriginCoords < 0) ) {
                msg <- c(msg, "OriginCoords of landscape must be positive.")
            }
        }
    }else{
        if (!is.null(object@OriginCoords)) {
            msg <- c(msg, "OriginCoords must be NULL if you use files as landscape inputs!")
        }
    }
    if (anyNA(object@K_or_DensDep) || length(object@K_or_DensDep)==0) {
        msg <- c(msg, "K_or_DensDep must be set!")
    }
    else {
        if (any(object@K_or_DensDep<0)) {
            msg <- c(msg, "K_or_DensDep must not be smaller than zero for any habitat type.")
        }
    }
    if (!object@HabPercent) {
        if (anyNA(object@Nhabitats) || length(object@Nhabitats)!=1) {
            msg <- c(msg, "Number of habitat codes must be given for the imported landscape. If your landscape is specified by a habitat cover percentage, please set \"HabPercent = TRUE\". ")
        }
        else {
            if (object@Nhabitats < 2){
                msg <- c(msg, "Number of habitat codes must be at least 2 when habitat codes are used.")
            }
        }
    }
    if (is.null(msg)) {
        if (object@HabPercent) {
            if (length(object@K_or_DensDep) != 1) {
                msg <- c(msg, "K_or_DensDep must be of length 1, specifying the maximum K or 1/b for 100% habitat.")
            }
        }
        else {
            if (length(object@K_or_DensDep) != object@Nhabitats) {
                msg <- c(msg, "K_or_DensDep must be a vector of length 'Nhabitats', specifying the respective K or 1/b for every habitat code.")
            }
        }
    }
    # if LandscapeFile is not NULL (aka given), PatchMatrix must be an empty list and PatchFile must be either
    # ALL NULL (cell-based) or has the same length as LandscapeFile (patch-based, potentially dynamic)
    if(length(object@LandscapeFile) > 0) { #or is.null(object@LandscapeFile)?
        if(length(object@PatchMatrix) > 0) {
            msg <- c(msg, "If LandscapeFile is given, you must provide all landscape inputs as files. PatchMatrix must be an empty list!")
        }
        if(length(object@PatchFile) > 0) { # if patch-based
            if (anyNA(object@PatchFile) || length(object@PatchFile)!=length(object@LandscapeFile)) {
                msg <- c(msg, "If patch-based, PatchFile must be the same length as LandscapeFile!")
            }
        }
    } else{ # if length(LandscapeFile==0) -> LandscapeMatrix!
        if(length(object@PatchFile) > 0) {
            msg <- c(msg, "If LandscapeMatrix is given, you must provide all landscape inputs as matrices. PatchFile must be NULL!")
        }
        if(length(object@PatchMatrix) > 0) {
            if (anyNA(object@PatchMatrix) || length(object@PatchMatrix)!=length(object@LandscapeMatrix)) {
                msg <- c(msg, "If patch-based, PatchMatrix must be the same length as LandscapeMatrix!")
            }
            if(any(sapply(object@PatchMatrix, class)[1,] != "matrix")){
                msg <- c(msg, "All elements of the PatchMatrix list must be of class matrix.")
            }
            else{
                if( (any(sapply(object@PatchMatrix, ncol) != land_ncol)) || (any(sapply(object@PatchMatrix, nrow) != land_nrow)) ){
                    msg <- c(msg, "All elements of PatchMatrix list must have the same ncol and nrow as the LandscapeMatrix list")
                }
            }
        }
    }
    # if LandscapeFile is not NULL (aka given), CostsMatrix must be an empty list and CostsFile must be either
    # NULL (no SMS) or has the same length as LandscapeFile (SMS costs, potentially dynamic)
    if(length(object@LandscapeFile) > 0) { #or is.null(object@LandscapeFile)?
        if(length(object@CostsMatrix) > 0) {
            msg <- c(msg, "If LandscapeFile is given, you must provide all landscape inputs as files. CostsMatrix must be an empty list!")
        }
        if(length(object@CostsFile) > 0) { # if SMS costs are given
            if (anyNA(object@CostsFile) || length(object@CostsFile)!=length(object@LandscapeFile)) {
                msg <- c(msg, "If SMS costs are given, CostsFile must be the same length as LandscapeFile!")
            }
        }
    } else{ # if length(LandscapeFile==0) -> LandscapeMatrix!
        if(length(object@CostsFile) > 0) {
            msg <- c(msg, "If LandscapeMatrix is given, you must provide all landscape inputs as matrices. CostsFile must be NULL!")
        }
        if(length(object@CostsMatrix) > 0) {
            if (anyNA(object@CostsMatrix) || length(object@CostsMatrix)!=length(object@LandscapeMatrix)) {
                msg <- c(msg, "If SMS costs are given, CostsMatrix must be the same length as LandscapeMatrix!")
            } else{
                if(any(sapply(object@CostsMatrix, class)[1,] != "matrix")){
                    msg <- c(msg, "All elements of CostsMatrix list must be of class matrix.")
                }
                else{
                    if( (any(sapply(object@CostsMatrix, ncol) != land_ncol)) || (any(sapply(object@CostsMatrix, nrow) != land_nrow)) ){
                        msg <- c(msg, "All elements of CostsMatrix list must have the same ncol and nrow as the LandscapeMatrix list")
                    }
                }
            }
        }
    }

    if(length(object@LandscapeFile) > 0) { # File input
        if(length(object@SpDistMatrix) > 0) {
            msg <- c(msg, "If LandscapeFile is given, you must provide all landscape inputs as files. SpDistMatrix must be an empty list!")
        } else{ # SpDistFile can be used
            if(!is.null(object@SpDistFile) && length(object@SpDistFile) != 1){
                msg <- c(msg, "Only one Species distribution file may be given.")
            }
        }

    } else{ # Matrix input
        if(length(object@SpDistFile) > 0) {
            msg <- c(msg, "If LandscapeMatrix is given, you must provide all landscape inputs as matrices. SpDistFile must be NULL!")
        } else{ # SpDistMatrix can be used
            if(length(object@SpDistMatrix) > 0){
                if(length(object@SpDistMatrix) != 1){
                    msg <- c(msg, "Only one Species distribution matrix may be given.")
                } else{
                    if( class(object@SpDistMatrix[[1]])[1] != "matrix") {
                        msg <- c(msg, "Species distribution must be of class matrix.")
                    }
                    else{
                        coarse <- object@SpDistResolution / object@Resolution
                        if( (ncol(object@SpDistMatrix[[1]]) != ceiling(land_ncol/coarse)) || (nrow(object@SpDistMatrix[[1]]) != ceiling(land_nrow/coarse)) ){
                            msg <- c(msg, "Extent of SpDistMatrix must match that of the LandscapeMatrix.")
                        }
                    }
                }
            }
        }
    }

    if (!is.null(object@SpDistFile) || length(object@SpDistMatrix)>0) {
        if (anyNA(object@SpDistResolution) || length(object@SpDistResolution)!=1) {
            msg <- c(msg, "Resolution of Species distribution must be set and of length 1!")
        }
        else {
            if (object@SpDistResolution < 1) {
                msg <- c(msg, "Resolution of Species distribution must be positive.")
            }
            else {
                if (object@SpDistResolution < object@Resolution) {
                    msg <- c(msg, "Resolution of Species distribution may not be less than Landscape Resolution.")
                }
                else {
                    if (object@SpDistResolution %% object@Resolution) {
                        msg <- c(msg, "SpDistResolution must be an integer multiple of Resolution.")
                    }
                }
            }
        }
    }

    if(anyNA(object@DynamicLandYears) || length(object@DynamicLandYears)==0) {
        msg <- c(msg, "DynamicLandYears must be set!")
    }
    else {
        if(length(object@LandscapeFile) > 0){
            if(length(object@LandscapeFile) != length(object@DynamicLandYears)){
                msg <- c(msg, "LandscapeFile and DynamicLandYears must have the same number of entries!")
            }
            else{
                if(object@DynamicLandYears[1] != 0){
                    msg <- c(msg, "The first entry of DynamicLandYears must be 0!")
                }
                else{
                    if(!all(sort(object@DynamicLandYears) == object@DynamicLandYears)){
                        msg <- c(msg, "DynamicLandYears must contain subsequent years!")
                    }
                }
            }
        } else {
            if(length(object@LandscapeMatrix) != length(object@DynamicLandYears)){
                msg <- c(msg, "LandscapeMatrix and DynamicLandYears must have the same number of entries!")
            }
            else{
                if(object@DynamicLandYears[1] != 0){
                    msg <- c(msg, "The first entry of DynamicLandYears must be 0!")
                }
                else{
                    if(!all(sort(object@DynamicLandYears) == object@DynamicLandYears)){
                        msg <- c(msg, "DynamicLandYears must contain subsequent years!")
                    }
                }
            }

        }
    }

    # demographic spatial variation
    # only valid for LandscapeMatrix and HabPercent=TRUE
    if(length(object@LandscapeFile) == 0 && object@HabPercent){
        if(length(object@nrDemogScaleLayers) != 1){
            msg <- c(msg, "nrDemogScaleLayers must be of length 1.")
        }
        else {
            if( object@nrDemogScaleLayers < 0){
                msg <- c(msg, "nrDemogScaleLayers must be positive.")
            }
            else {
                if (length(object@demogScaleLayers) > 0){ # scaling layers are given
                    if( any( sapply(object@demogScaleLayers, class) != "array")){
                        msg <- c(msg, "demogScaleLayers must be a list of arrays.")
                    }
                    else{
                        if( length(object@demogScaleLayers) != length(object@DynamicLandYears) ){
                            msg <- c(msg, "demogScaleLayers must be a list that contains an array for each element in DynamicLandYears.")
                        }
                        else{
                            ds_dims <- sapply(object@demogScaleLayers, dim) # (size of dimensions per array)
                            if( !"matrix" %in% class(ds_dims) ){
                                msg <- c(msg, "the arrays in demogScaleLayers must have the same dimensionality.")
                            }
                            else{
                                if( dim(ds_dims)[1] != 3 ){
                                    msg <- c(msg, "demogScaleLayers must be a list that contains 3-dimensional arrays.")
                                }
                                else{
                                    if( !all(apply(ds_dims, 2, function(col) identical(col, ds_dims[,1]))) ){
                                        msg <- c(msg, "all arrays in demogScaleLayers must have the same size.")
                                    }
                                    else{
                                        if( object@nrDemogScaleLayers > 0 &&  ds_dims[3,1] != object@nrDemogScaleLayers ){
                                            msg <- c(msg, "nrDemogScaleLayers must give the number of layers contained in each element (array) of demogScaleLayers.")
                                        }
                                        else{
                                            if( ds_dims[1,1] != land_nrow || ds_dims[2,1] != land_ncol ){
                                                msg <- c(msg, "All elements of demogScaleLayers list must have the same ncol and nrow as the LandscapeFile list")
                                            }
                                            else{
                                                ds_vals <- c(unlist(object@demogScaleLayers))
                                                ds_vals <- ds_vals[!is.na(ds_vals)]
                                                if( any( ds_vals < 0) || any( ds_vals > 100 ) ){
                                                    msg <- c(msg, "All elements of the arrays in demogScaleLayers must be values between 0 and 100.")
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    } else{
        if(length(object@demogScaleLayers) > 0 || object@nrDemogScaleLayers > 0){
            msg <- c(msg, "Demographic scaling layers can only be used with habitat quality maps (HabPercent=TRUE) and Landscape matrix inputs.")
        }
    }
   if (is.null(msg)) TRUE else msg}
)
setMethod("initialize", "ImportedLandscape", function(.Object, ...) {
    this_func = "ImportedLandscape(): "
    args <- list(...)
    .Object <- callNextMethod()
    if ( length(args) == 0 ) {
        validObject(.Object)
    }
    if (.Object@HabPercent) {
        .Object@Nhabitats = 1L
        if (!is.null(args$Nhabitats)) {
            warning(this_func, "Nhabitats", warn_msg_ignored, "for continuous habitat percentage landscape.", call. = FALSE)
        }
    }
    if (is.null(.Object@SpDistFile) && length(.Object@SpDistMatrix) == 0) {
        .Object@SpDistResolution = -9
        if (!is.null(args$SpDistResolution)) {
            warning(this_func, "Resolution of Species distribution", warn_msg_ignored, "since no map file is given.", call. = FALSE)
        }
    }
    if (.Object@HabPercent && length(.Object@LandscapeMatrix) > 0) {
        if (is.null(args$nrDemogScaleLayers)) {
            if (length(.Object@demogScaleLayers) > 0) {
                .Object@nrDemogScaleLayers = dim(.Object@demogScaleLayers[[1]])[3]
            }
        }
        else {
            if (length(.Object@demogScaleLayers) == 0) {
                warning(this_func, "nrDemogScaleLayers", warn_msg_ignored, "since no demogScaleLayers maps are given.", call. = FALSE)
            }
        }
    }
    else {
        if (!is.null(args$demogScaleLayers) | !is.null(args$nrDemogScaleLayers)) {
            warning(this_func, "demogScaleLayers", warn_msg_ignored, "since they can only be used in combination with habitat quality maps and landscape matrices as inputs.", call. = FALSE)
        }
    }
    .Object}
)
setMethod("show", "ImportedLandscape", function(object){
    cat(" Landscape imported from file")
    if(length(object@DynamicLandYears)==1) {
        if(!is.null(object@LandscapeFile)) cat(":\n  ", paste(object@LandscapeFile))
        else cat(":\n  used a landscape matrix")
    }
    cat("\n")
    if(object@HabPercent) {
        cat("   with continuous habitat percentages,\n")
        cat("   at 100%: K or 1/b =", paste(object@K_or_DensDep), "[inds per ha].\n")
    }
    else {
        cat("   with", paste(object@Nhabitats), "unique integer habitat code(s)\n")
        cat("   K or 1/b        : ", paste(object@K_or_DensDep), "[inds per ha].\n")
    }
    if (!is.null(object@PatchFile[1])) {
        cat("   Patches imported from file \n")
        if(length(object@DynamicLandYears)==1) {
            cat(paste(object@PatchFile),"\n")
        }
    }
    if(length(object@PatchMatrix)>0) {
        cat("   Patches imported from matrix\n")
    }
    if (!is.null(object@CostsFile[1])) {
        cat("   SMS costs imported from file \n")
        if(length(object@DynamicLandYears)==1) {
            cat(paste(object@CostsFile),"\n")
        }
    }
    if(length(object@CostsMatrix)>0) {
        cat("   SMS costs imported from matrix\n")
    }
    cat ("   Resolution      :", paste(object@Resolution),"\n")
    if(length(object@DynamicLandYears)>1) {
        if(!is.null(object@LandscapeFile)){
            if (is.null(object@PatchFile[1])) {
                if (is.null(object@CostsFile[1])) {
                    cat("   Land changes in\n    Year   to  Habitat file:\n")
                }
                else {
                    cat("   Land changes in\n    Year   to  Habitat file     Costs file:\n")
                }
            }
            else {
                if (is.null(object@CostsFile[1])) {
                    cat("   Land changes in\n    Year   to  Habitat file     Patch file:\n")
                }
                else {
                    cat("   Land changes in\n    Year   to  Habitat file     Patch file     Costs file:\n")
                }
            }
        }
        else {
            if (length(object@PatchMatrix)>0) {
                if (length(object@CostsMatrix)>0) {
                    cat("   Land changes in\n    Year   to  Habitat matrix, patch matrix and costs matrix:\n")
                }
                else {
                    cat("   Land changes in\n    Year   to  Habitat matrix and patch matrix:\n")
                }
            }
            else {
                if (length(object@CostsMatrix)>0) {
                    cat("   Land changes in\n    Year   to  Habitat matrix and costs matrix:\n")
                }
                else {
                    cat("   Land changes in\n    Year   to  Habitat matrix:\n")
                }
            }
        }

        for (a in 1:length(object@DynamicLandYears)) {
            if(!is.null(object@LandscapeFile)){
                if (is.null(object@PatchFile[1])) {
                    if (is.null(object@CostsFile[1])) {
                        cat("   ",paste(object@DynamicLandYears[a]),"  ",paste(object@LandscapeFile[a]),"\n")
                    }
                    else {
                        cat("   ",paste(object@DynamicLandYears[a]),"  ",paste(object@LandscapeFile[a]),"  ",paste(object@CostsFile[a]),"\n")
                    }
                }
                else {
                    if (is.null(object@CostsFile[1] )) {
                        cat("   ",paste(object@DynamicLandYears[a]),"  ",paste(object@LandscapeFile[a]),"  ",paste(object@PatchFile[a]),"\n")
                    }
                    else {
                        cat("   ",paste(object@DynamicLandYears[a]),"  ",paste(object@LandscapeFile[a]),"  ",paste(object@PatchFile[a]),"  ",paste(object@CostsFile[a]),"\n")
                    }
                }
            }
            else {
                cat("   ", paste(object@DynamicLandYears[a]), "\n")
            }

        }
    }
    if(!is.null(object@SpDistFile)) {
        cat(" Initial Species Distribution imported from file:\n  ", paste(object@SpDistFile), "\n")
        cat ("   Resolution      :", paste(object@SpDistResolution),"\n")
    }
    if(length(object@SpDistMatrix)>0) {
        cat(" Initial Species Distribution imported from matrix\n")
    }
    #cat ("\n")
    }
)


### Helper functions

## Check if a number is a power of 2
isPowerOf2 <- function(x) {
    n1s <- sum(as.numeric(intToBits(x)))
    if (n1s == 1) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

# RS manual 3.1.1 (page 51) - 3.1.2 (page 54)
