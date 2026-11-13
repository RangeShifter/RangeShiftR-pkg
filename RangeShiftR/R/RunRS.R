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
#' Runs a \emph{RangeShifter} simulation with the parameters given in the parameter master object.
#'
#' @param RSparams A parameter master object (class 'RSparams'), contains all parameters to specify the simulation settings.
#' Can be generated using the function \code{\link[RangeShiftR]{RSsim}}.
#' @param dirpath File path to RS working directory; must contain the folders 'Inputs', 'Outputs', 'Output_Maps'.\cr
#' If NULL, the current \code{R} working directory will be used.
#' @details
#' The figure shows the order of processes in a simulation. Each replicate starts with the initialisation of the landscape
#' and the initial individuals. Years are numbered from \eqn{0} to \code{Years}\eqn{-1}, and a replicate stops early
#' if the population goes extinct. At the start of each year the environment is updated (dynamic landscape,
#' environmental gradient, carrying capacities). This is followed by \code{RepSeasons} reproductive seasons. Management
#' actions (see \code{\link[RangeShiftR]{Management}}) are applied only in the first reproductive season of the
#' specified years. In each season, individuals reproduce, then disperse in the three phases emigration, transfer and settlement.
#'
#' The timing of survival and development depends on the survival schedule \code{SurvSched}
#' (see \code{\link[RangeShiftR]{StageStructure}}):
#' \itemize{
#'   \item \code{SurvSched = 0}: survival and development of adults directly after reproduction; of juveniles after dispersal.
#'   \item \code{SurvSched = 1}: survival and development of all stages after dispersal.
#'   \item \code{SurvSched = 2}: development after dispersal in each season; survival once at the end of the year.
#'   \item Non-structured models: all adults die after reproduction (non-overlapping generations); the offspring
#'         are subject to survival after dispersal.
#' }
#' At the end of each year individuals age, and individuals exceeding the maximum age die.
#' After the last year, the final outputs are written for \code{Year = Years}.
#' The timing of the individual outputs within a year is shown in \code{\link[RangeShiftR]{Simulation}}.
#'
#' \if{html}{\figure{simulation_flow.svg}{options: width="680" alt="Flow of a RangeShifter simulation"}}
#' \if{latex}{\figure{simulation_flow.pdf}{options: width=5in}}
#'
#' @return
#' \itemize{
#'   \item With \code{ReturnPopMatrix = TRUE} (see \code{\link[RangeShiftR]{Simulation}}): a named list of matrices,
#'         one per replicate, output year and reproductive season, named \code{rep<r>_year<y>_season<s>_NInd}; with \code{ReturnStages} additionally
#'         \code{rep<r>_year<y>_season<s>_NInd_stage<i>} (and \code{rep<r>_year<y>_season<s>_NJuv} for stage \eqn{0}).
#'   \item With \code{ReturnPopDataFrame = TRUE}: a data frame with the columns \code{Rep}, \code{Year}, \code{RepSeason}, \code{PatchID},
#'         \code{totalAbundance} and, with \code{ReturnStages}, \code{NInd_Stage<i>}.
#'   \item Otherwise \code{NULL}. All results are written to the folder \file{Outputs}.
#'   \item If an error occurs, a list with the element \code{Errors} containing the error code.
#' }
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
            if(RSparams@simul@ReturnPopMatrix){
                resol = RSparams@control@resolution
		        if (RSparams@control@threadsafe){
# 			        if(class(RSparams@land)=="ImportedLandscape") llcorner = RSparams@land@OriginCoords
# 				        else  llcorner = c(0,0)
#
# 			        raster_list <- vector("list", length(out))
#
			        # for(i in seq_along(out)) {
			        #
			        #     r <- terra::rast(out[[i]])
			        #     # ext(r) <- c(llcorner[1], ncol(out[[1]])*resol+llcorner[1], llcorner[2], nrow(out[[1]])*resol+llcorner[2])
			        #
			        #     raster_list[[i]] <- r
			        #
			        #     out[[i]] <- 0
			        # }
#
# 			        gc()
#
#
# 			#         raster_list <- lapply(out, function(x) {
# 			# 	        r <- terra::rast(x)
# 			# 	        ext(r) <- c(llcorner[1], ncol(out[[1]])*resol+llcorner[1], llcorner[2], nrow(out[[1]])*resol+llcorner[2])
# 			# 	        return(r)
# 			# 	    })
#
#                     out <- NULL
# 				    return(terra::rast(raster_list))
		            return(out)
		        } else {
# 		            raster_list <- vector("list", length(out))
#
# 		            for(i in seq_along(out)) {
#
# 		                r <- terra::rast(out[[i]])
# 		                ext(r) <- c(0, ncol(out[[i]]) * resol,
# 		                            0, nrow(out[[i]]) * resol)
#
# 		                raster_list[[i]] <- r
#
# 		                out[[i]] <- NULL
# 		            }
#
# 		            gc()
#
#     #     		    raster_list <- lapply(out, function(x) {
#     #     		        r <- terra::rast(x)
#     #     		        ext(r) <- c(0, ncol(out[[1]])*resol, 0, nrow(out[[1]])*resol)
#     #     		        return(r)
#     # 		        })
# 			        return(terra::rast(raster_list))
		            return(out)
		        }
            }
            if(RSparams@simul@ReturnPopDataFrame){
                stage_idx <- which(RSparams@simul@ReturnStages) - 1L # ReturnStages[1] is stage 0 (juveniles)
                df <- do.call(
                    rbind,
                    lapply(names(out), function(nm) {

                        mat <- out[[nm]]

                        # extract Rep, Year and RepSeason from list name
                        rep    <- as.integer(sub(".*rep([0-9]+).*", "\\1", nm))
                        year   <- as.integer(sub(".*year([0-9]+).*", "\\1", nm))
                        season <- as.integer(sub(".*season([0-9]+).*", "\\1", nm))

                        # matrix -> data.frame
                        d <- as.data.frame(mat, stringsAsFactors = FALSE)

                        # base column names
                        colnames(d)[1:2] <- c("PatchID", "totalAbundance")

                        # stage column names based on ReturnStages
                        if (ncol(d) > 2) {
                            colnames(d)[3:ncol(d)] <- paste0(
                                "NInd_Stage",
                                stage_idx[seq_len(ncol(d) - 2)]
                            )
                        }

                        # add Rep, Year and RepSeason
                        d$Rep       <- rep
                        d$Year      <- year
                        d$RepSeason <- season

                        # reorder columns
                        id_cols <- c("Rep", "Year", "RepSeason")
                        d[, c(id_cols, colnames(d)[!colnames(d) %in% id_cols])]
                    })
                )
                return(df)
            }

        }
        else return(NULL)
    }
    else {
        return(out)
    }
}
