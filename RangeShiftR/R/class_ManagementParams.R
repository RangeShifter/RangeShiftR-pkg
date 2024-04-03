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

# Subclass holding the 'Translocation 'parameter

#' Set Translocation Parameters
#'
#' This function translocates individuals from one location to another in a given year.
#' You should set at least one year for a translocation event,
#' but you can also define multiple years. You can define one or multiple locations (patches or cells) from which individuals
#' should be selected and to which individuals are translocated.
#' The number and characteristics of selected individuals may vary between years.
#' Individuals may be selected according to their age, stage and sex, depending on the type of population.
#'
#' Each unique combination of source and target location
#' and criteria for individuals to be selected represent one translocation event.
#'
#' Additionally, you can define a catching success rate to define the success of catching a selected individual.
#'
#'
#' @usage Translocation(years = 1,
#'                      TransLocMat = 0,
#'                      catching_rate = 0.5
#'                       )
#'
#' @param years Vector of years in which translocation events should take place
#' @param TransLocMat Matrix of translocation events. Each row represents a translocation event. Columns represent:
#'
#'  - the year of the event,\cr
#'  - the source location (\code{Patch_ID} in case of cell-based models or \code{X} and \code{Y} location in case of cell-based models),\cr
#'  - the target location (\code{Patch_ID} in case of cell-based models or \code{X} and \code{Y} location in case of cell-based models),\cr
#'  - the number of individuals which are tried to be catched,\cr
#'  - minimal age of each individual,\cr
#'  - maximal age of the individual,\cr
#'  - stage of the individual,\cr
#'  - sex of the individual\cr
#' @param catching_rate Catching success rate
#'
#' @details
#'
#' To select individuals with certain criteria for translocation, a Translocation Matrix \code{TransLocMat} needs to be generated in which the characteristics of the individuals are defined for each translocation event.
#'
#' In the columns of the \code{TransLocMat} the year, source and target location as well as the number of individuals and the characteristics are defined in the following order:
#'
#' - \code{year} of the translocation event \cr
#' - \code{Patch_ID} or \code{X} and \code{Y} location of the source location from which individuals are catched. \cr
#' - \code{Patch_ID} or \code{X} and \code{Y} location of the target location to which individuals are transferred. \cr
#' - \code{nb_catch} how many individuals of the given set of characteristics are tried to be catched \cr
#' - \code{min_age} minimal age of the individual in the interval of 0-MaxAge. Set to -9 to ignore. \cr
#' - \code{max_age} maximal age of the individual in the interval of \code{min_age}-MaxAge. Set to -9 to ignore.\cr
#' - \code{stage} of the individual. Set to -9 to ignore. \cr
#' - \code{sex} of the individual: Only for sexual models, otherwise set to -9 to ignore  \cr
#'
#' \code{min_age}, \code{max_age} and \code{stage} can only be defined for stage structured models. For non stage structured models, or to ignore the characteristic, set the value to -9.
#' \code{sex} can only be defined for sexual models. For non sexual models, or to ignore the characteristic, set the value to -9.
#'
#' To avoid unsuccessful translocation events, you should set the range of allowed characteristics as broad as possible.
#'
#' Each row of the translocation matrix \code{TransLocMat} holds a unique combination of a source and a target location, the number of individuals to catch as well as characteristics of these individuals.
#' You may add more than one translocation event per year, i.e. there may be multiple source and target locations for each year of a translocation
#' event or multiple individual characteristics for a certain pair of source and target location.
#'
#' In each \code{year} of a translocation event, individuals matching the given criteria in the source location are collected and a given number of individuals \code{nb_catch}
#' are sampled from this subset. The success of catching one of these sampled individuals in the source location is defined by the \code{catching_rate}.
#' Successfully catched individuals are then translocated to the target location and the individual will be assigend the status 4 (successfully dispersed) to avoid further dispersal.
#'
#' The source as well as the target location are required to be habitat patches or cells of the landscape. Otherwise the translocation event will not be skipped.
#'
#' @references
#'         \insertAllCited{}
#' @return a parameter object of class "ManagementParams"
#' @author Jette Reeg
#' @name Translocation
#' @export Translocation
Translocation <- setClass("TranslocationParams", slots = c(years = "numeric",
                                                           TransLocMat = "matrix",
                                                           catching_rate = "numeric")
                          , prototype = list(years = -9L,
                                             TransLocMat = matrix(0, nrow = 8, ncol = 8),
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

    # Check if TransLocMat is not NA, is a matrix and has at least as many rows and years
    if (any(is.na(object@TransLocMat)) && !is.matrix(object@TransLocMat)){
        msg <- c(msg, "TransLocMat must be defined and be a matrix")
    }else {
        if(nrow(object@TransLocMat) < length(object@years)) {
            msg <- c(msg, "TransLocMat must have at least as many rows as years of Transclocation events")
        } else {
            if(!all(sort(object@TransLocMat[,1]) == object@TransLocMat[,1])){
                msg <- c(msg, "Translocation matrix must contain subsequent years!")
            } else{
                if (ncol(object@TransLocMat) != 8 && ncol(object@TransLocMat) != 10) { # 8 is only true for patch-based models; for cell based models it should be 10
                    msg <- c(msg, "TransLocMat must have 8 or 10 columns: year, source location (patch ID OR 2 columns X and Y), target location (patch ID OR 2 columns X and Y), number of individuals, min age, max age, stage.")
                }
            }

        }
    }

    # Check if TransLocMat has valid entries:
    # for patch based models:
    # Check if second column (source patch ID) is numeric and ID exists
    # check if third column (target patch ID) is numeric and ID exists
    # check if fourth column (number of individuals) is numeric and greater than 0
    # check if fifth column (min age) is numeric and equal or greater than 0 or -9
    # check if sixth column (max age) is numeric and between min age and MaxAge or -9
    # check if seventh column (stage) is numeric and either between 0 and number of stages or -9
    # check if eighth column (sex) is numeric and either 0, 1 or -9
    # for cell based models:
    # Check if second and fourth column (source and target X cells) is numeric and within X boundaries of the landscape
    # Check if third and fifth column (source and target Y cells) is numeric and within Y boundaries of the landscape
    # check if sixth column (number of individuals) is numeric and greater than 0
    # check if seventh column (min age) is numeric and equal or greater than 0 or -9
    # check if eighth column (max age) is numeric and between min age and MaxAge or -9
    # check if nineth column (stage) is numeric and either between 0 and number of stages or -9
    # check if tenth column (sex) is numeric and either 0, 1 or -9


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
    print(object@TransLocMat)
})



#### MANAGEMENT ####

# define this ClassUnion so that the 'Translocation' slot in the parameter master class 'RSparams' can be FALSE for not considering translocations
setClassUnion("TranslocationSlot", c("logical", "TranslocationParams"))

# Superclass holding the subclasses 'Translocation', but eventually also 'Hunting' and 'Poaching'

#' Set Management Parameters
#'
#' Set all management parameters you wish to apply on your simulations. This includes for now only the translocation of individuals.
#'
#'
#' @usage Management(Translocation = Translocation())
#'
#' @param Translocation Individuals are selected in one location and then transfered to a new location. You can define the year(s) of the translocation event(s),
#' the source and target location, the number of individuals to be translocated and characteristics of individuals to be sampled from. The catching rate defines the success catching an individual.
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
Management <- setClass("ManagementParams", slots = c(Translocation = "TranslocationSlot")
                      , prototype = list(Translocation = FALSE)
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

