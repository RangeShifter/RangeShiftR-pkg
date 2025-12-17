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


# -----
#
# Run a RangeShiftR simulation
#
# -----

#' Run the Simulation
#'
#' @param RSparams A parameter master object (class 'RSparams'), contains all parameters to specify the simulation settings.
#' Can be generated using the function \code{\link[RangeShiftR]{RSsim}}.
#' @param dirpath File path to RS working directory; must contain the folders 'Inputs', 'Outputs', 'Output_Maps'.\cr
#' If NULL, the current \code{R} working directory will be used.
#' @return returns an error code
#' @export
RunRS <- function(RSparams, dirpath = getwd()){
    if (missing(RSparams)) {
        stop("Missing parameter object")
    }
    else{
        if (class(RSparams)[1] != "RSparams") {
            stop("Parameter object must be of class RSparams")
        }
        else {
            validObject(RSparams)
        }
    }
    if (is.null(dirpath)) {
        dirpath = getwd()
    }

    out = run_from_R(dirpath, RSparams)

    if (class(out)=="list" && is.null(out$Errors)) {
        if ( length(out)>0 ) {
            resol = RSparams@control@resolution
		if (RSparams@control@threadsafe){
			 if(class(RSparams@land)=="ImportedLandscape") llcorner = RSparams@land@OriginCoords
				    else  llcorner = c(0,0)
				    raster_list <- lapply(out, function(x) {
				        r <- terra::rast(x)
				        ext(r) <- c(llcorner[1], ncol(out[[1]])*resol+llcorner[1], llcorner[2], nrow(out[[1]])*resol+llcorner[2])
				        return(r)
				    })
				    return(terra::rast(raster_list))
		} else {
    		    raster_list <- lapply(out, function(x) {
    		        r <- terra::rast(x)
    		        ext(r) <- c(0, ncol(out[[1]])*resol, 0, nrow(out[[1]])*resol)
    		        return(r)
    		    })
			    return(terra::rast(raster_list))
		}
        }
        else return(NULL)
    }
    else {
        return(out)
    }
}
