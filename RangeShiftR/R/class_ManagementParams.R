#---------------------------------------------------------------------------
#
#	Copyright (C) 2024 Jette Reeg, Anne-Kathleen Malchow, Greta Bocedi, Stephen C.F. Palmer, Justin M.J. Travis, Damaris Zurell
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

# Sublass holding the 'Translocation 'parameter

#' Set Translocation Parameters
#'
#' Translocation transfers individuals from one location to another in a given year. You should set at least one year for a translocation event,
#' but you can also define multiple years. You can one or multiple locations (patches or cells) from which individuals should be catched and selected for translocation.
#' The number and characteristics of selected individuals may vary between years. Individuals are selected according to the given minimal and maximal age, stage,
#' sex nd dispersal status. Additionally, you can define a catching success rate to define the success of really catching an individual.
#'
#'
#' @usage Translocation(years = 1,
#'                      TranLocMat = 0,
#'                      catching_rate = 0.5
#'                       )
#'
#' @param years Vector of years in which translocation events should take place
#' @param TranLocMat Matrix of translocation events. Each row represents a translocation event. Columns represent:
#'
#'  - the year of the event,
#'  - the source location (\code{Patch_ID} in case of cell-based models or \code{X} and \code{Y} location in case of cell-based models),
#'  - the target location,
#'  - the number of individuals which are tried to be catched,
#'  - minimal age of each individual,
#'  - maximal age of the individual,
#'  - stage of the individual,
#'  - sex of the individual
#' @param catching_rate Catching success rate
#'
#' @details
#'
#' In each \code{year} of a translocation event, individuals are selected from a given source location and transferred to a new target location.
#' The success of catching an individual in the source location is defined by the \code{catching_rate}.
#' To select only individuals within a certain range of characteristics, a Translocation Matrix \code{TranLocMat} is defined in which the characteristics of the individuals are defined for each translocation event.
#' Each row of the translocation matrix \code{TranLocMat} holds a unique combination of source and target location as well as characteristics of catched individuals.
#' You may add more than one translocation event per year, i.e. there may be multiple source and target locations for each year of a translocation event or multiple individual characteristics for a certain pair of source and target location.
#'
#' In the columns of the \code{TranLocMat} the year, source and target location as well as the number of individuals and the characteristics are defined:
#'
#' - \code{year} of the translocation event
#' - \code{Patch_ID} or \code{X} and \code{Y} location of the source location from which individuals are catched
#' - \code{Patch_ID} or \code{X} and \code{Y} location of the target location to which individuals are transferred
#' - \code{nb_catch} how many individuals of the given set of characteristics are tried to be catched
#' - \code{min_age} minimal age of the individual: 0-MaxAge, Set to 0 to ignore lower boundary
#' - \code{max_age} maximal age of the individual: \code{min_age}-MaxAge, , Set 0 to ignore upper boundary
#' - \code{stage} of the individual: Only for non-stage structured models, otherwise set to 'NA'
#' - \code{sex} of the individual: Only for sexual models, otherwise set to 'NA'
#'
#' To avoid unsuccessful translocation events, you should set the range of allowed characteristics as broad as possible.
#'
#' @references
#'         \insertAllCited{}
#' @return a parameter object of class "ManagementParams"
#' @author Jette Reeg
#' @name Translocation
#' @export Translocation
Translocation <- setClass("TranslocationParams", slots = c(years = "numeric",
                                                           TranLocMat = "matrix",
                                                           catching_rate = "numeric")
                          , prototype = list(years = 0L,
                                             TranLocMat = matrix(0, nrow = 8, ncol = 8),
                                             catching_rate = 1L
                          )
)

setValidity("TranslocationParams", function(object) {
    msg <- NULL
    # Check if years is not NA, has at least one entry and is either a numeric or integer
    if (any(is.na(object@years)) || length(object@years)==0) {
        msg <- c(msg, "Years must be defined")
    }else{
        if (!is.integer(object@years) && !is.numeric(object@years)) {
            msg <- c(msg, "Years must be numeric or integer")
        }
    }

    # Check if TranLocMat is not NA, is a matrix and has at least as many rows and years
    if (any(is.na(object@TranLocMat)) && !is.matrix(object@TranLocMat)){
        msg <- c(msg, "TranLocMat must be defined and be a matrix")
    }else {
        if(nrow(object@TranLocMat) < length(object@years)) {
            msg <- c(msg, "TranLocMat must have at least as many rows as years of Transclocation events")
        } else {
            if(!all(sort(object@TranLocMat[,1]) == object@TranLocMat[,1])){
                msg <- c(msg, "Translocation matrix must contain subsequent years!")
            } else{
                if (ncol(object@TranLocMat) != 8) {
                    msg <- c(msg, "TranLocMat must have 8 columns: year, source location, target location, number of individuals, min age, max age, stage.")
                }
            }

        }
    }

    # Check if catching_rate is not NA and is a numeric
    if (is.na(object@catching_rate)) {
        msg <- c(msg, "Catching rate must be defined")
    }else{
        if (!is.numeric(object@catching_rate)) {
            msg <- c(msg, "Catching rate must be numeric")
        }
    }

    if (is.null(msg)) TRUE else msg}
)

setMethod("initialize", "TranslocationParams", function(.Object,...) {
    this_func = "Translocation(): "
    args <- list(...)
    .Object <- callNextMethod()
    if ( length(args) == 0 ) {
        validObject(.Object)
    }
    .Object}
)


setMethod("show", "TranslocationParams", function(object){
    cat(" Translocation:\n")
    cat("  Years of translocation events: ", object@years, "\n")
    cat("  Catching rate: ", object@catching_rate, "\n")
    cat("  Translocation Matrix: \n")
    print(object@TranLocMat)
})



#### MANAGEMENT ####

# Superclass holding the subclasses 'Translocation', but eventually also 'Hunting' and 'Poaching'

#' Set Management Parameters
#'
#' Set all management parameters you wish to apply on your simulations. This includes for now only the translocation of individuals.
#'
#'
#' @usage Management(Translocation = Translocation())
#'
#' @param Translocation Individuals are selected in one location and then transferred to a new location. You can define the year(s) of the translocation event(s),
#' the source and target location, the number of individuals to be translocated and catching rates to define the success of really catching an individual.
#' See \code{\link{Translocation}} for more details.
#'
#' @details TODO: Add general information on management options and translocations
#'
#' @references
#'         \insertAllCited{}
#' @return a parameter object of class "ManagementParams"
#' @author Jette Reeg
#' @name Management
#' @export Management
Management <- setClass("ManagementParams", slots = c(Translocation = "TranslocationParams")
                      , prototype = list(Translocation = Translocation())
)

setValidity("ManagementParams", function(object) {
    msg <- NULL
    validObject(object@Translocation)
    if (is.null(msg)) TRUE else msg}
)
setMethod("initialize", "ManagementParams", function(.Object,...) {
    this_func = "Management(): "
    args <- list(...)
    .Object <- callNextMethod()
    if ( length(args) == 0 ) {
        validObject(.Object)
    }
    .Object}
)
setMethod("show", "ManagementParams", function(object){
    cat(" Management: \n  Translocation:\n")
    print(object@Translocation)
})

