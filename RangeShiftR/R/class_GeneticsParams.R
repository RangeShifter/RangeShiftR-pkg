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

# define ClassUnion for Positions: vector of integers or "random"
setClassUnion("character_OR_integer", c("character", "integer"))

### CLASS GENETICSPARAMS

# from RS 'Traits' file

### SUBCLASS NEUTRALTRAITS

#' Set genetic traits structure for neutral traits
#'
#' @param Positions Loci positions coding for trait within genome. Must be in the
#' range 0-(GenomeSize-1), specified in Genetics file. Positions can overlap across
#' traits - there will be no pleiotropy, but this will influence genetic linkage.,
#' @param NbOfPositions Only specify if above is set to ‘random’, else must be blank (NULL)
#' @param InitialDistribution Distribution from which to draw initial allele values from.
#' If \code{uniform}. Initialise with random characters between 0 – max. Can be left blank (NULL)
#' and assume every individual is the ‘wildtype’ (max value provided in  \code{InitialParameters} ) and new mutations alter that.
#' @param InitialParameters Maximal value for the uniform distribution.
#' @param MutationDistribution Distribution for mutations to draw from. Can be either 'KAM' or 'SSM'. KAM (k-alleles model) is randomly
#' draw a value between 0 and max (see MutationParameters).
#' SSM (single-step mutation) is to move in a stepwise manner, A to B, B to C.
#' @param MutationParameters Parameters for the above distribution: maximal value for KAM or SSM (cannot exceed 255)
#' @param MutationRate Mutation rate applicable to this type of loci. Must be between 0.0 and 1.0
#' @param OutputValues If OutputGeneValues in \code{\link[RangeShiftR]{Genetics}} is
#' enabled, should allele values for this gene be written to output? Ignored if OutputGeneValues is set to \code{FALSE}.
#'
#' @return a parameter object of class "NeutralTraitsParams"
#' @author Jette Reeg
#' @name NeutralTraits
#' @export NeutralTraits
NeutralTraits<- setClass("NeutralTraitsParams", slots = c(Positions = "ANY", # vector of numbers or "random"
                                                         NbOfPositions = "ANY", # random or list of integer values
                                                         InitialDistribution = "character", # uniform
                                                         InitialParameters = "integer_OR_numeric", # max value
                                                         MutationDistribution = "character", # KAM or SSM
                                                         MutationParameters = "integer_OR_numeric", # max
                                                         MutationRate = "integer_OR_numeric", # float
                                                         OutputValues = "logical")
                   , prototype = list(Positions = "random", # "random" or list of integer values
                                      NbOfPositions = NULL, # numeric, only of positions random
                                      InitialDistribution = NULL, # uniform (neutral + dispersal), normal (dispersal), NULL (genetic load)
                                      InitialParameters = 2, # neutral traits: only max value; dispersal: two values: either min/max oder mean+sd, not applicable for genetic load
                                      MutationDistribution = "KAM", # neutral: "KAM" or "SSM", genetic load: "gamma", "uniform", "normal", "negExp", dispersal: uniform or normal
                                      MutationParameters = 2, # single value or 2 values
                                      MutationRate = 0.0, # numeric
                                      OutputValues = FALSE
                   ))
setValidity("NeutralTraitsParams", function(object) {
    msg <- NULL

    # Check Position and NbOfPositions
    isNumeric <- class(object@Positions) == "integer"
    isCharacter <- class(object@Positions) == "character"
    if (is.null(object@Positions) || (!isNumeric && !isCharacter)) {
        msg <- c(msg, "Positions must be integer or character.")
    }
    if (!isNumeric){
        if(isCharacter && object@Positions != "random") {
            msg <- c(msg, "Positions in neutral genetics must be either a vector of integers, or random.")
        }
    }
    if (isCharacter && object@Positions == "random") {
        if (is.null(object@NbOfPositions) && object@NbOfPositions <= 0) {
            msg <- c(msg, "NbrOfPositions must be a strictly positive integrer.")
        }
    } else if (!is.null(object@NbOfPositions)) {
        msg <- c(msg, "If Positions is not random, NbrOfPositions must not be set (NULL).")
    }
    # Check InitialDistribution
    if (!is.null(object@InitialDistribution) && object@InitialDistribution != "uniform") {
        msg <- c(msg,"InitialDistribution must be either uniform or left blank (NuLL) for the neutral trait.")
    }
    # Check InitialParameters
    if (object@InitialDistribution == "uniform") {
        # should only be one value (max)
        if (length(object@InitialParameters)!=1) {
            msg <- c(msg, "For neutral traits with uniform initialisation, InitialParameters must only supply the maximal value")
            }
            else {
                if (object@InitialParameters > 255) {
                    msg <- c(msg, "For neutral trait with uniform initialisation, max parameter must be between 0 and 255.")
                }
            }
        }
    ## if not uniform then initDist must be blank, no params
    else if (is.null(object@InitialDistribution)) {
        msg <- c(msg, "For neutral trait with no initialisation, InitialParameters must not be set (NULL)")

    }
    # Check mutation rate
    if (!is.null(object@MutationRate) && (object@MutationRate < 0 || object@MutationRate > 1)) {
        msg <- c(msg, "MutationRate must be between 0 and 1.")
    }
    # Check MutationDistribution and MutationParameters
    if (object@MutationDistribution == "KAM" || object@MutationDistribution == "SSM") {
        if (length(object@MutationParameters) != 1) {
            msg <- c(msg, "For a neutral trait, mutationParams must have form max=int.")
        }
        else {
            if (object@MutationParameters > 255) {
                msg <- c(msg, "For the neutral trait mutation max parameter must be between 0 and 255.")
            }
        }
    }
    else if (!is.null(object@MutationDistribution)){
        msg <- c(msg, "For a neutral trait, mutationDistribution must be either KAM or SSM.")
    }

    if (object@OutputValues != TRUE && object@OutputValues != FALSE) {
        msg <- c(msg, "OutputValues for neutral genetics must be either TRUE or FALSE.")
    }
    if (is.null(msg)) TRUE else msg}
)
setMethod("initialize", "NeutralTraitsParams", function(.Object, ...) {
    this_func = "NeutralTraits(): "
    args <- list(...)
    .Object <- callNextMethod()
    if ( length(args) == 0 ) {
        validObject(.Object)
    }
    .Object
})
setMethod("show", "NeutralTraitsParams", function(object){
    cat("   Neutral Genetics: \n")
    if(is.numeric(object@Positions)) cat("     Loci positions coding for trait: ", object@Positions, "\n")
    if(!is.numeric(object@Positions) && object@Positions=="random") cat("    Loci positions coding for trait randomly chosen with ", object@NbOfPositions, " positions\n")
    cat("     Initial distribution: ", object@InitialDistribution, "\n")
    cat("     Initial parameters: ", object@InitialParameters, "\n")
    cat("     Mutation distribution: ", object@MutationDistribution, "\n")
    cat("     Mutation parameters: ", object@MutationParameters, "\n")
    cat("     Mutation rate: ", object@MutationRate, "\n")
    if(object@OutputValues) cat("     Allel values for gene is written to output")
})


### SUBCLASS GENETICLOADTRAITS

#' Set genetic structure for genetic fitness traits
#'
#' @usage GeneticLoadTraits(NbGeneticLoads = 1, Positions = list("random"), NbOfPositions = 10,
#' DominanceDistribution = "normal", DominanceParameters = matrix(c(0.5,0.1), nrow=1),
#' MutationDistribution = "normal", MutationParameters = matrix(c(0.5,0.2), nrow=1),
#' MutationRate = 0.0001, OutputValues = FALSE)
#'
#' @param NbGeneticLoads Number of genetic loads
#' @param Positions Loci positions coding for that trait within genome. Should be provided as a list of strings (if random) and/or vectors of integers (if not random)
#' @param NbOfPositions Only specify when the \code{Positions} of a genetic load trait are set to ‘random’, else must be blank (NULL)
#' @param DominanceDistribution Distribution of dominance values. Can be \code{gamma}, \code{uniform}, \code{normal}, \code{negExp}, \code{scaled}. Should be provided as a vector of strings if \code{NbGeneticLoads} > 1
#' @param DominanceParameters Parameters for the dominance distribution: You must provide two colums for \code{uniform}, \code{normal} and \code{gamma} distributions: min and max (uniform), mean and sd (normal) or shape and scale (gamma) or one column for \code{negExp}, \code{scaled}: mean
#' If genetic loads have different \code{DominanceDistribution} and one require two columns you need to set the second value to NA in case of \code{negExp} or \code{scaled} distribution.
#' Each row in the matrix corresponds to a genetic load trait.
#' @param MutationDistribution Distribution for mutations to draw from. Can be \code{gamma}, \code{uniform}, \code{normal}, \code{negExp}. Should be provided as a vector of strings if \code{NbGeneticLoads} > 1
#' @param MutationParameters Parameters for the mutation distribution: You must provide two colums for \code{uniform}, \code{normal} and \code{gamma} distributions: min and max (uniform), mean and sd (normal) or shape and scale (gamma) or one column for \code{negExp}: mean
#' If genetic loads have different \code{DominanceDistribution} and one require two columns you need to set the second value to NA in case of \code{negExp} distribution.
#' @param MutationRate Mutation rate applicable to this type of loci. Must be between 0.0 and 1.0. Should be provided as a vector if multiple genetic loads are specified.
#' @param OutputValues If OutputGeneValues in GeneticsFile is enabled, should allele values for this gene be written to output? Ignored if OutputGeneValues is set to FALSE. Should be provided as a vector if multiple genetic loads are specified.
#'
#' The expression type of genetic load traits is always multiplicative.
#'
#' @return a parameter object of class "GeneticLoadParams"
#' @author Jette Reeg
#' @name GeneticLoadTraits
#' @export GeneticLoadTraits
GeneticLoadTraits<- setClass("GeneticLoadParams", slots = c(
    NbGeneticLoads = "integer_OR_numeric", # number of genetic loads
    Positions = "list",# "random" or list of integer values
    NbOfPositions = "ANY", # numeric, only where positions are random; otherwise NA
    DominanceDistribution = "character", # ‘gamma’, ‘uniform’, ‘normal’, ‘negExp’, ‘scaled’
    DominanceParameters = "matrix", # 2 values for min/max, mean/sd, shape/scale or one value: mean
    MutationDistribution = "character", # ‘gamma’, ‘uniform’, ‘normal’,‘negExp’
    MutationParameters = "matrix", #  2 values for min/max, mean/sd, shape/scale or one value: mean
    MutationRate = "numeric", # float
    OutputValues = "logical")
    , prototype = list(
        NbGeneticLoads = 1L,
        Positions = list("random"),
        NbOfPositions = 2L,
        DominanceDistribution = "normal",
        DominanceParameters = matrix(c(0.5,0.1), nrow=1),
        MutationDistribution = "normal",
        MutationParameters = matrix(c(0.5,0.2), nrow=1),
        MutationRate = 0.001,
        OutputValues = FALSE
    ))
setValidity("GeneticLoadParams", function(object) {
    msg <- NULL
    # only 5 GeneticLoads are allowed
    if (object@NbGeneticLoads < 1 || object@NbGeneticLoads > 5) {
        msg <- c(msg, "Number of genetic loads must be between 1 and 5.")
    }
    # Check Position and NbOfPositions
    # Positions must be of type list?
    if(class(object@Positions) != "list") {
        msg <- c(msg, "GeneticLoad(): Positions must be provided as a list.")
    }
    # NbOfPositions must be either numeric, integer or NULL
    if (!is.null(object@NbOfPositions) && class(object@NbOfPositions) != "numeric" && class(object@NbOfPositions) != "integer") {
        msg <- c(msg, "GeneticLoad(): NbrOfPositions must be either NULL (if all positions are given) or numeric (if at least one genetic load has random positions).")
    }

    if (length(object@Positions) != object@NbGeneticLoads) {
        msg <- c(msg, "For each genetic load you must provide the positions.")
    } else if (all(object@Positions == "random")){ # if all positions are random
        if(length(object@NbOfPositions[!is.na(object@NbOfPositions)]) != object@NbGeneticLoads ) {
            msg <- c(msg, "GeneticLoad(): For each genetic load with random positions you must provide the number of positions.")
            }
        } else{ # if NOT all positions are random
            isNumeric <- sapply(Positions, is.numeric)
            if (!all(isNumeric)) { # if not all are numeric,
                if(object@Positions[isNumeric==FALSE] != "random"){ # then those not numeric must be random
                    msg <- c(msg, "GeneticLoad(): Positions in genetic loads must be either a vector of integers or random.")
                }
                if (any(is.na(object@NbOfPositions[object@Positions == "random"])) || any(object@NbOfPositions[object@Positions == "random"] <= 0)){ # if number of positions are NA or smaller than 0
                    msg <- c(msg, "GeneticLoad(): NbrOfPositions must be set to a strictly positive integer for random positions.")
                }
                if (any(!is.na(object@NbOfPositions[object@Positions != "random"]))) { # if there are NbOfPositions supplied for non-random positions
                    msg <- c(msg, "GeneticLoad(): if Positions is not random NbrOfPositions must be not be set (NA).")
                }
            }
            else { # if all positions are not random
               if (!is.null(object@NbOfPositions)) {
                    msg <- c(msg, "GeneticLoad(): If positions are not random, you must not specify the number of positions (NbOfPositions).")
                }
            }
        }

    # Check DominanceDistribution
    if (!is.null(object@DominanceDistribution)){
        if(length(object@DominanceDistribution) != object@NbGeneticLoads) {
            msg <- c(msg, "For each genetic load you must provide the DominanceDistribution.")
        } else if (nrow(object@DominanceParameters) != object@NbGeneticLoads) {
            msg <- c(msg, "If you have set DominanceDistributions you must provide the DominanceParameters for each genetic load. Use one row for each genetic load.")
        } else {
            if (any(object@DominanceDistribution == "normal")) { # if any distribution is normal
                # two values for mean and sd
                if (ncol(object@DominanceParameters) !=2 || # if DominanceParameters has not 2 columns OR
                    any(!is.numeric(object@DominanceParameters[object@DominanceDistribution=="normal"])) || # if entries are not numeric
                    any(is.na(object@DominanceParameters[object@DominanceDistribution=="normal"]))) { # if entries are NA
                    msg <- c(msg,"For a normal dominance distribution, DominanceParams must provide two values for mean (first column) and sd (second column)")
                }
            }
            if (any(object@DominanceDistribution == "gamma")) {
                # two values for shape and scale
                if (ncol(object@DominanceParameters) !=2 || # if DominanceParameters has not 2 columns OR
                    any(!is.numeric(object@DominanceParameters[object@DominanceDistribution=="gamma"])) || # if entries are not numeric
                    any(is.na(object@DominanceParameters[object@DominanceDistribution=="gamma"]))) { # if entries are NA
                    msg <- c(msg,"For a gamma dominance distribution, DominanceParams must provide two values for shape (first column) and scale (second column)")
                }
            }
            if (any(object@DominanceDistribution == "uniform")) {
                # two values for min and max
                if (ncol(object@DominanceParameters) !=2 || # if DominanceParameters has not 2 columns OR
                    any(!is.numeric(object@DominanceParameters[object@DominanceDistribution=="uniform"])) || # if entries are not numeric
                    any(is.na(object@DominanceParameters[object@DominanceDistribution=="uniform"]))) { # if entries are NA
                    msg <- c(msg,"For a uniform dominance distribution, DominanceParams must provide two values for min (first column) and max (second column)")
                }
            }
            if (all(object@DominanceDistribution == "negExp" || object@DominanceDistribution == "scaled")) { # if it is only negExp or scaled
                # one value for mean and one NA
                if (ncol(object@DominanceParameters) !=1 || # if DominanceParameters has more than 1 column
                    !is.numeric(object@DominanceParameters) ||
                    any(is.na(object@DominanceParameters))
                     ) {
                    msg <- c(msg,"For negative exponential and scaled dominance distribution, DominanceParams must provide only one column for mean.")
                }
            } else{
                if (any(object@DominanceDistribution == "scaled" || object@DominanceDistribution == "negExp")) { # if only some are scaled or negative exponential
                    # one value for mean and one NA if other distributions need 2 values
                    if (ncol(object@DominanceParameters) !=2 || # if DominanceParameters has not 2 columns OR
                        !is.numeric(object@DominanceParameters) || # if entries are not numeric
                        !all(is.na(object@DominanceParameters[object@DominanceDistribution=="scaled",2])) || # second column is not NA
                        !all(is.na(object@DominanceParameters[object@DominanceDistribution=="negExp",2])) || # second column is not NA
                        any(is.na(object@DominanceParameters[object@DominanceDistribution=="scaled",1])) || # first column is NA
                        any(is.na(object@DominanceParameters[object@DominanceDistribution=="negExp",1])) # first column is NA
                        ) {
                        msg <- c(msg,"For the scaled or negative exponential dominance distribution, DominanceParameters must provide only one value for mean (first column) and the second column need to be NA if other genetic loads use other dominance distributions.")
                    }
                }
            }

            if (any(object@DominanceDistribution != "normal" && object@DominanceDistribution != "gamma" &&
                    object@DominanceDistribution != "uniform" && object@DominanceDistribution != "negExp" && object@DominanceDistribution != "scaled")) {
                msg <- c(msg, "DominanceDistribution must be either normal, gamma, uniform, negExp or scaled for genetic load traits.")
            }
        }
    }

    # Check mutation rate
    if(!is.numeric(object@MutationRate) || length(object@MutationRate) != object@NbGeneticLoads){
        msg <- c(msg, "You must provide the mutation rate for each genetic load as a numeric vector.")
    } else {
        if (!is.numeric(object@MutationRate) ||  any(object@MutationRate < 0.0 || object@MutationRate > 1.0)) {
            msg <- c(msg, "MutationRate must be between 0.0 and 1.0.")
        }
    }

    # Check MutationDistribution and MutationParameters
    if (!is.null(object@MutationDistribution)){
        if (length(object@MutationDistribution) != object@NbGeneticLoads){
            msg <- c(msg, "For each genetic load you must provide the MutationDistribution.")
        } else if (nrow(object@MutationParameters) != object@NbGeneticLoads) {
            msg <- c(msg, "For each genetic load you must provide the MutationParameters.")
        } else {
            if (any(object@MutationDistribution == "uniform", object@MutationDistribution == "normal", object@MutationDistribution == "gamma")){
                if (ncol(object@MutationParameters) !=2){
                    msg <- c(msg,"MutationParams must provide two values for uniform/normal/gamma distribution: min/mean/shape (first column) and max/sd/scale (second column)")
                }
            }
            if (any(object@MutationDistribution == "uniform")) {
                # two values for min and max
                if (any(!is.numeric(object@MutationParameters[object@MutationDistribution=="uniform"])) ||
                    any(is.na(object@MutationParameters[object@MutationDistribution=="uniform"]))) {
                    msg <- c(msg,"For a uniform mutation distribution, MutationParams must provide two values for min (first column) and max (second column)")
                }
            }
            if (any(object@MutationDistribution == "normal")) {
                # two values for meand and sd
                if (any(!is.numeric(object@MutationParameters[object@MutationDistribution=="normal"])) ||
                    any(is.na(object@MutationParameters[object@MutationDistribution=="normal"]))) {
                    msg <- c(msg,"For a normal mutation distribution, MutationParams must provide two values for mean (first column) and sd (second column)")
                }
            }
            if (any(object@MutationDistribution == "gamma")) {
                # two values for shape and scale
                if (any(!is.numeric(object@MutationParameters[object@MutationDistribution=="gamma"])) ||
                    any(is.na(object@MutationParameters[object@MutationDistribution=="gamma"]))) {
                    msg <- c(msg,"For a gamma mutation distribution, MutationParams must provide two values for shape (first column) and scale (second column)")
                }
            }
            if (all(object@MutationDistribution == "negExp")) { # if it is only negExp
                # one value for mean
                if (ncol(object@MutationParameters) !=1 || # if MutationParameters has more than 1 column
                    !is.numeric(object@DominanceParameters) ||
                    any(is.na(object@DominanceParameters))
                ) {
                    msg <- c(msg,"For negative exponential mutation distribution, MutationParams must provide only one column for mean.")
                }
            } else{
                if (any(object@MutationDistribution == "negExp")) { # if only some are scaled or negative exponential
                    # one value for mean and one NA if other distributions need 2 values
                    if (ncol(object@MutationParameters) !=2 || # if DominanceParameters has not 2 columns OR
                        !is.numeric(object@MutationParameters) || # if entries are not numeric
                        !all(is.na(object@MutationParameters[object@MutationDistribution=="negExp",2])) || # second column is not NA
                        any(is.na(object@MutationParameters[object@MutationDistribution=="negExp",1])) # first column is NA
                    ) {
                        msg <- c(msg,"For the negative exponential mutation distribution, MutationParameters must provide only one value for mean (first column) and the second column need to be NA if other genetic loads use other mutation distributions.")
                    }
                }
            }

            if (any(object@MutationDistribution != "normal" && object@MutationDistribution != "gamma" &&
                    object@MutationDistribution != "uniform" && object@MutationDistribution != "negExp")) {
                msg <- c(msg, "MutationDistribution must be either normal, gamma, uniform or negExp for genetic load traits.")
            }

        }
    }

    if (object@OutputValues != TRUE && object@OutputValues != FALSE) {
        msg <- c(msg, "OutputValues for genetic loads must be either TRUE or FALSE.")
    }
    if (is.null(msg)) TRUE else msg
})
setMethod("initialize", "GeneticLoadParams", function(.Object, ...) {
    this_func = "GeneticLoad(): "
    args <- list(...)
    .Object <- callNextMethod()
    if ( length(args) == 0 ) {
        validObject(.Object)
    }
    .Object
})
setMethod("show", "GeneticLoadParams", function(object){
})

### SUBCLASS EMIGRATIONTRAITS

#' Set genetic traits structure for Emigration traits
#'
#'
#' @return a parameter object of class "EmigrationTraitsParams"
#' @author Jette Reeg
#' @name EmigrationTraits
#' @export EmigrationTraits
EmigrationTraits<- setClass("EmigrationTraitsParams", slots = c(Positions = "character_OR_integer", #
                                                                NbOfPositions = "character_OR_integer", # random or list of integer values
                                                                ExpressionType = "character", # additive or average
                                                                InitialDistribution = "character", # uniform or normal
                                                                InitialParameters = "integer_OR_numeric", # min and max value or mean and sd
                                                                IsInherited = "logical", # T/F
                                                                MutationDistribution = "character", # uniform or normal
                                                                MutationParameters = "integer_OR_numeric", # min mx or mean sd
                                                                MutationRate = "integer_OR_numeric", # float
                                                                OutputValues = "logical")
                                                                , prototype = list(
                                                                    # ExprSex = FALSE, # is TRUE as soon as Emigration is sexdependent
                                                                    # TraitType = NULL, # dispersal: "E_0", "E_alpha", "E_beta"; cannot be 0 is determined by emigration settings
                                                                    Positions = NULL, # "random" or list of integer values
                                                                    NbOfPositions = NULL, # numeric, only of positions random
                                                                    ExpressionType = NULL, # dispersal: "additive" or "average"; geneticload: "multiplicative"
                                                                    InitialDistribution = NULL, # uniform (neutral + dispersal), normal (dispersal), NULL (genetic load)
                                                                    InitialParameters = NULL, # neutral traits: only max value; dispersal: two values: either min/max oder mean+sd, not applicable for genetic load
                                                                    IsInherited = FALSE, # only for dispersal
                                                                    MutationDistribution = NULL, # neutral: "KAM" or "SSM", genetic load: "gamma", "uniform", "normal", "negExp", dispersal: uniform or normal
                                                                    MutationParameters = NULL, # single value or 2 values
                                                                    MutationRate = NULL, # numeric
                                                                    OutputValues = FALSE
                                                                ))
setValidity("EmigrationTraitsParams", function(object) {
    msg <- NULL
    # I would expect for each parameter a matrix depending on the parameter master -> most things need to be checked in RSparams, where I can extract SexDep, DensDep, etc.
    # Check Position and NbOfPositions
    patternPositions <-  "^\"?(([0-9]+-)?[0-9]+,)*([0-9]+-)?[0-9]+\"?$"

    isMatch <- grepl(patternPositions, object@Positions)
    if (!all(isMatch) && object@Positions[isMatch==FALSE] != "random") {
        msg <- c(msg, "In EmigrationTraits(): Positions must be either a comma-separated list of integer ranges, or random.")
    }
    if (any(is.na(object@NbOfPositions[object@Positions == "random"])) || object@NbOfPositions[object@Positions == "random"] <= 0){
        msg <- c(msg, "NbrOfPositions must be set to a strictly positive integrer.")
    } else if (!is.na(object@NbOfPositions[object@Positions != "random"])) {
        msg <- c(msg, "In EmigrationTraits(): if Positions is not random NbrOfPositions must be not be set (NA).")
    }

    # Check ExpressionType must be additive or average
    if (!all(object@ExpressionType %in% c("additive", "average"))) {
        msg <- c(msg, "In EmigrationTraits(): ExpressionType must be either additive or average.")
    }

    # Check InitialDistribution and InitialParameter: Distribution must be uniform or normal
    if (is.null(object@InitialDistribution) || !all(object@InitialDistribution %in% c("uniform", "normal"))) {
        msg <- c(msg,"In EmigrationTraits(): InitialDistribution must be either uniform or normal.")
    }
    # Check InitialParameters: must be numeric values
    if (any(!is.numeric(object@InitialParameters)) || is.null(object@InitialParameters)) {
        msg <- c(msg, "In EmigrationTraits(): InitialParameters must be provided.")
    }

    # Check IsInherited: must be TRUE or FALSE
    if (!is.logical(object@IsInherited) || any(is.na(object@IsInherited))) {
        msg <- c(msg, "In EmigrationTraits(): IsInherited must be either TRUE or FALSE." )
        }
    # Check MutationRate;
    if (object@MutationRate[object@IsInherited == "TRUE"] < 0.0 || object@MutationRate[object@IsInherited == "TRUE"] > 1.0) {
        msg <- c(msg, "In EmigrationTraits(): MutationRate must be between 0.0 and 1.0.")
    }   else if (any(!is.na(object@MutationRate[object@IsInherited == "FALSE"]))) {
            msg <- c(msg, "In EmigrationTraits(): If isInherited if off, mutationRate must be blank (NA).")
    }

    # Check MutationDistribution
    # in all places, where the isInherited matrix is TRUE, MutationDistribution must be uniform or normal
    if (any(object@MutationDistribution[object@IsInherited == "TRUE"] != "uniform" && object@MutationDistribution[object@IsInherited == "TRUE"] != "normal")) {
        msg <- c(msg, "In EmigrationTraits(): MutationDistribution must be either uniform or normal if they are inherited.")
    }
    # in all places where isInherited matrix is FALSE, MutationDistribution must be NA
    if (any(!is.na(object@MutationDistribution[object@IsInherited == "FALSE"]))) {
        msg <- c(msg, "In EmigrationTraits(): If isInherited is off, MutationDistribution must be blank (NA).")
    }

    # Check MutationParameters
    # in all places, where the isInherited matrix is TRUE, MutationParameters must be numeric
    if (any(!is.numeric(object@MutationParameters[object@IsInherited == "TRUE"])) || any(is.na(object@MutationParameters[object@IsInherited == "TRUE"]))){
        msg <- c(msg, "In EmigrationTraits(): MutationParameters must be numeric if they are inherited.")
    }
    # in all places where isInherited matrix is FALSE, MutationParameters must be NA
    if (any(!is.na(object@MutationParameters[object@IsInherited == "FALSE"]))) {
        msg <- c(msg, "In EmigrationTraits(): If isInherited is off, MutationParameters must be blank (NA).")
    }

    # Check OutputValues
    if (object@OutputValues != TRUE && object@OutputValues != FALSE) {
        msg <- c(msg, "In EmigrationTraits(): OutputValues for emigration traits must be either TRUE or FALSE.")
    }

    if (is.null(msg)) TRUE else msg
})
setMethod("initialize", "EmigrationTraitsParams", function(.Object, ...) {
    this_func = "EmigrationTraits(): "
    args <- list(...)
    .Object <- callNextMethod()
    if ( length(args) == 0 ) {
        validObject(.Object)
    }
    .Object
})
setMethod("show", "EmigrationTraitsParams", function(object){
})

### SUBCLASS SETTLEMENTTRAITS

#' Set genetic traits structure for neutral traits
#'
#' Can I have a list of different dispersaltraits() depending on the dispersal parameter? (E0. Ealpa, Ebeta,)
#'
#' @return a parameter object of class "SettlementTraitsParams"
#' @author Jette Reeg
#' @name SettlementTraits
#' @export SettlementTraits
SettlementTraits<- setClass("SettlementTraitsParams", slots = c(Positions = "character_OR_integer", #
                                                                 NbOfPositions = "character_OR_integer", # random or list of integer values
                                                                 ExpressionType = "character", # additive or average
                                                                 InitialDistribution = "character", # uniform or normal
                                                                 InitialParameters = "integer_OR_numeric", # min and max value or mean and sd
                                                                 IsInherited = "logical", # T/F
                                                                 MutationDistribution = "character", # uniform or normal
                                                                 MutationParameters = "integer_OR_numeric", # min mx or mean sd
                                                                 MutationRate = "integer_OR_numeric", # float
                                                                 OutputValues = "logical")
                            , prototype = list(Positions = NULL, # "random" or list of integer values
                                NbOfPositions = NULL, # numeric, only of positions random
                                ExpressionType = NULL, # dispersal: "additive" or "average"; geneticload: "multiplicative"
                                InitialDistribution = NULL, # uniform (neutral + dispersal), normal (dispersal), NULL (genetic load)
                                InitialParameters = NULL, # neutral traits: only max value; dispersal: two values: either min/max oder mean+sd, not applicable for genetic load
                                IsInherited = FALSE, # only for dispersal
                                MutationDistribution = NULL, # neutral: "KAM" or "SSM", genetic load: "gamma", "uniform", "normal", "negExp", dispersal: uniform or normal
                                MutationParameters = NULL, # single value or 2 values
                                MutationRate = NULL, # numeric
                                OutputValues = FALSE
                            ))
setValidity("SettlementTraitsParams", function(object) {
    msg <- NULL
    # I would expect for each parameter a matrix depending on the parameter master -> most things need to be checked in RSparams, where I can extract SexDep, DensDep, etc.
    # Check Position and NbOfPositions
    patternPositions <-  "^\"?(([0-9]+-)?[0-9]+,)*([0-9]+-)?[0-9]+\"?$"

    isMatch <- grepl(patternPositions, object@Positions)
    if (!all(isMatch) && object@Positions[isMatch==FALSE] != "random") {
        msg <- c(msg, "In SettlementTraits(): Positions must be either a comma-separated list of integer ranges, or random.")
    }
    if (any(is.na(object@NbOfPositions[object@Positions == "random"])) || object@NbOfPositions[object@Positions == "random"] <= 0){
        msg <- c(msg, "NbrOfPositions must be set to a strictly positive integrer.")
    } else if (!is.na(object@NbOfPositions[object@Positions != "random"])) {
        msg <- c(msg, "In SettlementTraits(): if Positions is not random NbrOfPositions must be not be set (NA).")
    }

    # Check ExpressionType must be additive or average
    if (!all(object@ExpressionType %in% c("additive", "average"))) {
        msg <- c(msg, "In SettlementTraits()ExpressionType must be either additive or average.")
    }

    # Check InitialDistribution and InitialParameter: Distribution must be uniform or normal
    if (is.null(object@InitialDistribution) || !all(object@InitialDistribution %in% c("uniform", "normal"))) {
        msg <- c(msg,"In SettlementTraits(): InitialDistribution must be either uniform or normal.")
    }
    # Check InitialParameters: must be numeric values
    if (any(!is.numeric(object@InitialParameters)) || is.null(object@InitialParameters)) {
        msg <- c(msg, "In SettlementTraits(): InitialParameters must be provided.")
     }

     # Check IsInherited: must be TRUE or FALSE
     if (!is.logical(object@IsInherited) || any(is.na(object@IsInherited))) {
         msg <- c(msg, "In SettlementTraits(): IsInherited must be either TRUE or FALSE." )
     }
     # Check MutationRate;
     if (object@MutationRate[object@IsInherited == "TRUE"] < 0.0 || object@MutationRate[object@IsInherited == "TRUE"] > 1.0) {
         msg <- c(msg, "In SettlementTraits(): MutationRate must be between 0.0 and 1.0.")
     }   else if (any(!is.na(object@MutationRate[object@IsInherited == "FALSE"]))) {
         msg <- c(msg, "In SettlementTraits(): If isInherited if off, mutationRate must be blank (NA).")
     }

     # Check MutationDistribution
     # in all places, where the isInherited matrix is TRUE, MutationDistribution must be uniform or normal
     if (any(object@MutationDistribution[object@IsInherited == "TRUE"] != "uniform" && object@MutationDistribution[object@IsInherited == "TRUE"] != "normal")) {
         msg <- c(msg, "In SettlementTraits(): MutationDistribution must be either uniform or normal if they are inherited.")
     }
     # in all places where isInherited matrix is FALSE, MutationDistribution must be NA
     if (any(!is.na(object@MutationDistribution[object@IsInherited == "FALSE"]))) {
         msg <- c(msg, "In SettlementTraits(): If isInherited is off, MutationDistribution must be blank (NA).")
     }

     # Check MutationParameters
     # in all places, where the isInherited matrix is TRUE, MutationParameters must be numeric
     if (any(!is.numeric(object@MutationParameters[object@IsInherited == "TRUE"])) || any(is.na(object@MutationParameters[object@IsInherited == "TRUE"]))){
         msg <- c(msg, "In SettlementTraits(): MutationParameters must be numeric if they are inherited.")
     }
     # in all places where isInherited matrix is FALSE, MutationParameters must be NA
     if (any(!is.na(object@MutationParameters[object@IsInherited == "FALSE"]))) {
         msg <- c(msg, "In SettlementTraits(): If isInherited is off, MutationParameters must be blank (NA).")
     }

     # Check OutputValues
     if (object@OutputValues != TRUE && object@OutputValues != FALSE) {
         msg <- c(msg, "In SettlementTraits(): OutputValues for emigration traits must be either TRUE or FALSE.")
     }

     if (is.null(msg)) TRUE else msg
})
setMethod("initialize", "SettlementTraitsParams", function(.Object, ...) {
    this_func = "SettlementTraits(): "
    args <- list(...)
    .Object <- callNextMethod()
    if ( length(args) == 0 ) {
        validObject(.Object)
    }
    .Object
})
setMethod("show", "SettlementTraitsParams", function(object){
})

### SUBCLASS CRWTRAITS

#' Set genetic traits structure for CRW traits
#'
#' Can I have a list of different CRWtraits() depending on the dispersal parameter? (E0. Ealpa, Ebeta,)
#'
#' @return a parameter object of class "CRWTraitsParams"
#' @author Jette Reeg
#' @name CRWTraits
#' @export CRWTraits
CRWTraits<- setClass("CRWTraitsParams", slots = c(Positions = "character_OR_integer", #
                                                                 NbOfPositions = "character_OR_integer", # random or list of integer values
                                                                 ExpressionType = "character", # additive or average
                                                                 InitialDistribution = "character", # uniform or normal
                                                                 InitialParameters = "integer_OR_numeric", # min and max value or mean and sd
                                                                 IsInherited = "logical", # T/F
                                                                 MutationDistribution = "character", # uniform or normal
                                                                 MutationParameters = "integer_OR_numeric", # min mx or mean sd
                                                                 MutationRate = "integer_OR_numeric", # float
                                                                 OutputValues = "logical")
                            , prototype = list(Positions = NULL, # "random" or list of integer values
                                NbOfPositions = NULL, # numeric, only of positions random
                                ExpressionType = NULL, # dispersal: "additive" or "average"; geneticload: "multiplicative"
                                InitialDistribution = NULL, # uniform (neutral + dispersal), normal (dispersal), NULL (genetic load)
                                InitialParameters = NULL, # neutral traits: only max value; dispersal: two values: either min/max oder mean+sd, not applicable for genetic load
                                IsInherited = FALSE, # only for dispersal
                                MutationDistribution = NULL, # neutral: "KAM" or "SSM", genetic load: "gamma", "uniform", "normal", "negExp", dispersal: uniform or normal
                                MutationParameters = NULL, # single value or 2 values
                                MutationRate = NULL, # numeric
                                OutputValues = FALSE
                            ))
setValidity("CRWTraitsParams", function(object) {
    msg <- NULL
    # I would expect for each parameter a matrix depending on the parameter master -> most things need to be checked in RSparams, where I can extract SexDep, DensDep, etc.
    # Check Position and NbOfPositions
    patternPositions <-  "^\"?(([0-9]+-)?[0-9]+,)*([0-9]+-)?[0-9]+\"?$"

    isMatch <- grepl(patternPositions, object@Positions)
    if (!all(isMatch) && object@Positions[isMatch==FALSE] != "random") {
        msg <- c(msg, "In CRWTraits(): Positions must be either a comma-separated list of integer ranges, or random.")
    }
    if (any(is.na(object@NbOfPositions[object@Positions == "random"])) || object@NbOfPositions[object@Positions == "random"] <= 0){
        msg <- c(msg, "NbrOfPositions must be set to a strictly positive integrer.")
    } else if (!is.na(object@NbOfPositions[object@Positions != "random"])) {
        msg <- c(msg, "In CRWTraits(): if Positions is not random NbrOfPositions must be not be set (NA).")
    }

    # Check ExpressionType must be additive or average
    if (!all(object@ExpressionType %in% c("additive", "average"))) {
        msg <- c(msg, "In CRWTraits(): ExpressionType must be either additive or average.")
    }

    # Check InitialDistribution and InitialParameter: Distribution must be uniform or normal
    if (is.null(object@InitialDistribution) || !all(object@InitialDistribution %in% c("uniform", "normal"))) {
        msg <- c(msg,"In CRWTraits(): InitialDistribution must be either uniform or normal.")
    }
    # Check InitialParameters: must be numeric values
    if (any(!is.numeric(object@InitialParameters)) || is.null(object@InitialParameters)) {
        msg <- c(msg, "In CRWTraits(): InitialParameters must be provided.")
    }

    # Check IsInherited: must be TRUE or FALSE
    if (!is.logical(object@IsInherited) || any(is.na(object@IsInherited))) {
        msg <- c(msg, "In CRWTraits(): IsInherited must be either TRUE or FALSE." )
    }
    # Check MutationRate;
    if (object@MutationRate[object@IsInherited == "TRUE"] < 0.0 || object@MutationRate[object@IsInherited == "TRUE"] > 1.0) {
        msg <- c(msg, "In CRWTraits(): MutationRate must be between 0.0 and 1.0.")
    }   else if (any(!is.na(object@MutationRate[object@IsInherited == "FALSE"]))) {
        msg <- c(msg, "In CRWTraits(): If isInherited if off, mutationRate must be blank (NA).")
    }

    # Check MutationDistribution
    # in all places, where the isInherited matrix is TRUE, MutationDistribution must be uniform or normal
    if (any(object@MutationDistribution[object@IsInherited == "TRUE"] != "uniform" && object@MutationDistribution[object@IsInherited == "TRUE"] != "normal")) {
        msg <- c(msg, "In CRWTraits(): MutationDistribution must be either uniform or normal if they are inherited.")
    }
    # in all places where isInherited matrix is FALSE, MutationDistribution must be NA
    if (any(!is.na(object@MutationDistribution[object@IsInherited == "FALSE"]))) {
        msg <- c(msg, "In CRWTraits(): If isInherited is off, MutationDistribution must be blank (NA).")
    }

    # Check MutationParameters
    # in all places, where the isInherited matrix is TRUE, MutationParameters must be numeric
    if (any(!is.numeric(object@MutationParameters[object@IsInherited == "TRUE"])) || any(is.na(object@MutationParameters[object@IsInherited == "TRUE"]))){
        msg <- c(msg, "In CRWTraits(): MutationParameters must be numeric if they are inherited.")
    }
    # in all places where isInherited matrix is FALSE, MutationParameters must be NA
    if (any(!is.na(object@MutationParameters[object@IsInherited == "FALSE"]))) {
        msg <- c(msg, "In CRWTraits(): If isInherited is off, MutationParameters must be blank (NA).")
    }

    # Check OutputValues
    if (object@OutputValues != TRUE && object@OutputValues != FALSE) {
        msg <- c(msg, "In CRWTraits(): OutputValues for emigration traits must be either TRUE or FALSE.")
    }

    if (is.null(msg)) TRUE else msg
})
setMethod("initialize", "CRWTraitsParams", function(.Object, ...) {
    this_func = "CRWTraits(): "
    args <- list(...)
    .Object <- callNextMethod()
    if ( length(args) == 0 ) {
        validObject(.Object)
    }
    .Object
})
setMethod("show", "CRWTraitsParams", function(object){
})

### SUBCLASS KERNELTRAITS

#' Set genetic traits structure for kernel traits
#'
#' Can I have a list of different kerneltraits() depending on the dispersal parameter? (E0. Ealpa, Ebeta,)
#'
#' @return a parameter object of class "KernelTraitsParams"
#' @author Jette Reeg
#' @name KernelTraits
#' @export KernelTraits
KernelTraits<- setClass("KernelTraitsParams", slots = c(Positions = "character_OR_integer", #
                                                   NbOfPositions = "character_OR_integer", # random or list of integer values
                                                   ExpressionType = "character", # additive or average
                                                   InitialDistribution = "character", # uniform or normal
                                                   InitialParameters = "integer_OR_numeric", # min and max value or mean and sd
                                                   IsInherited = "logical", # T/F
                                                   MutationDistribution = "character", # uniform or normal
                                                   MutationParameters = "integer_OR_numeric", # min mx or mean sd
                                                   MutationRate = "integer_OR_numeric", # float
                                                   OutputValues = "logical")
                     , prototype = list(Positions = NULL, # "random" or list of integer values
                         NbOfPositions = NULL, # numeric, only of positions random
                         ExpressionType = NULL, # dispersal: "additive" or "average"; geneticload: "multiplicative"
                         InitialDistribution = NULL, # uniform (neutral + dispersal), normal (dispersal), NULL (genetic load)
                         InitialParameters = NULL, # neutral traits: only max value; dispersal: two values: either min/max oder mean+sd, not applicable for genetic load
                         IsInherited = FALSE, # only for dispersal
                         MutationDistribution = NULL, # neutral: "KAM" or "SSM", genetic load: "gamma", "uniform", "normal", "negExp", dispersal: uniform or normal
                         MutationParameters = NULL, # single value or 2 values
                         MutationRate = NULL, # numeric
                         OutputValues = FALSE
                     ))
setValidity("KernelTraitsParams", function(object) {
    msg <- NULL
    # I would expect for each parameter a matrix depending on the parameter master -> most things need to be checked in RSparams, where I can extract SexDep, DensDep, etc.
    # Check Position and NbOfPositions
    patternPositions <-  "^\"?(([0-9]+-)?[0-9]+,)*([0-9]+-)?[0-9]+\"?$"

    isMatch <- grepl(patternPositions, object@Positions)
    if (!all(isMatch) && object@Positions[isMatch==FALSE] != "random") {
        msg <- c(msg, "In KernelTraits(): Positions must be either a comma-separated list of integer ranges, or random.")
    }
    if (any(is.na(object@NbOfPositions[object@Positions == "random"])) || object@NbOfPositions[object@Positions == "random"] <= 0){
        msg <- c(msg, "NbrOfPositions must be set to a strictly positive integrer.")
    } else if (!is.na(object@NbOfPositions[object@Positions != "random"])) {
        msg <- c(msg, "In KernelTraits(): if Positions is not random NbrOfPositions must be not be set (NA).")
    }

    # Check ExpressionType must be additive or average
    if (!all(object@ExpressionType %in% c("additive", "average"))) {
        msg <- c(msg, "In KernelTraits(): ExpressionType must be either additive or average.")
    }

    # Check InitialDistribution and InitialParameter: Distribution must be uniform or normal
    if (is.null(object@InitialDistribution) || !all(object@InitialDistribution %in% c("uniform", "normal"))) {
        msg <- c(msg,"In KernelTraits(): InitialDistribution must be either uniform or normal.")
    }
    # Check InitialParameters: must be numeric values
    if (any(!is.numeric(object@InitialParameters)) || is.null(object@InitialParameters)) {
        msg <- c(msg, "In KernelTraits(): InitialParameters must be provided.")
    }

    # Check IsInherited: must be TRUE or FALSE
    if (!is.logical(object@IsInherited) || any(is.na(object@IsInherited))) {
        msg <- c(msg, "In KernelTraits(): IsInherited must be either TRUE or FALSE." )
    }
    # Check MutationRate;
    if (object@MutationRate[object@IsInherited == "TRUE"] < 0.0 || object@MutationRate[object@IsInherited == "TRUE"] > 1.0) {
        msg <- c(msg, "In KernelTraits(): MutationRate must be between 0.0 and 1.0.")
    }   else if (any(!is.na(object@MutationRate[object@IsInherited == "FALSE"]))) {
        msg <- c(msg, "In KernelTraits(): If isInherited if off, mutationRate must be blank (NA).")
    }

    # Check MutationDistribution
    # in all places, where the isInherited matrix is TRUE, MutationDistribution must be uniform or normal
    if (any(object@MutationDistribution[object@IsInherited == "TRUE"] != "uniform" && object@MutationDistribution[object@IsInherited == "TRUE"] != "normal")) {
        msg <- c(msg, "In KernelTraits(): MutationDistribution must be either uniform or normal if they are inherited.")
    }
    # in all places where isInherited matrix is FALSE, MutationDistribution must be NA
    if (any(!is.na(object@MutationDistribution[object@IsInherited == "FALSE"]))) {
        msg <- c(msg, "In KernelTraits(): If isInherited is off, MutationDistribution must be blank (NA).")
    }

    # Check MutationParameters
    # in all places, where the isInherited matrix is TRUE, MutationParameters must be numeric
    if (any(!is.numeric(object@MutationParameters[object@IsInherited == "TRUE"])) || any(is.na(object@MutationParameters[object@IsInherited == "TRUE"]))){
        msg <- c(msg, "In KernelTraits(): MutationParameters must be numeric if they are inherited.")
    }
    # in all places where isInherited matrix is FALSE, MutationParameters must be NA
    if (any(!is.na(object@MutationParameters[object@IsInherited == "FALSE"]))) {
        msg <- c(msg, "In KernelTraits(): If isInherited is off, MutationParameters must be blank (NA).")
    }

    # Check OutputValues
    if (object@OutputValues != TRUE && object@OutputValues != FALSE) {
        msg <- c(msg, "In KernelTraits(): OutputValues for emigration traits must be either TRUE or FALSE.")
    }

    if (is.null(msg)) TRUE else msg
})
setMethod("initialize", "KernelTraitsParams", function(.Object, ...) {
    this_func = "KernelTraits(): "
    args <- list(...)
    .Object <- callNextMethod()
    if ( length(args) == 0 ) {
        validObject(.Object)
    }
    .Object
})
setMethod("show", "KernelTraitsParams", function(object){
})

### SUBCLASS SMSTRAITS

#' Set genetic traits structure for SMS traits
#'
#' Can I have a list of different SMStraits() depending on the dispersal parameter? (E0. Ealpa, Ebeta,)
#'
#' @return a parameter object of class "SMSTraitsParams"
#' @author Jette Reeg
#' @name SMSTraits
#' @export SMSTraits
SMSTraits<- setClass("SMSTraitsParams", slots = c(Positions = "character_OR_integer", #
                                                         NbOfPositions = "character_OR_integer", # random or list of integer values
                                                         ExpressionType = "character", # additive or average
                                                         InitialDistribution = "character", # uniform or normal
                                                         InitialParameters = "integer_OR_numeric", # min and max value or mean and sd
                                                         IsInherited = "logical", # T/F
                                                         MutationDistribution = "character", # uniform or normal
                                                         MutationParameters = "integer_OR_numeric", # min mx or mean sd
                                                         MutationRate = "integer_OR_numeric", # float
                                                         OutputValues = "logical")
                        , prototype = list(Positions = NULL, # "random" or list of integer values
                            NbOfPositions = NULL, # numeric, only of positions random
                            ExpressionType = NULL, # dispersal: "additive" or "average"; geneticload: "multiplicative"
                            InitialDistribution = NULL, # uniform (neutral + dispersal), normal (dispersal), NULL (genetic load)
                            InitialParameters = NULL, # neutral traits: only max value; dispersal: two values: either min/max oder mean+sd, not applicable for genetic load
                            IsInherited = FALSE, # only for dispersal
                            MutationDistribution = NULL, # neutral: "KAM" or "SSM", genetic load: "gamma", "uniform", "normal", "negExp", dispersal: uniform or normal
                            MutationParameters = NULL, # single value or 2 values
                            MutationRate = NULL, # numeric
                            OutputValues = FALSE
                        ))
setValidity("SMSTraitsParams", function(object) {
    msg <- NULL
    # I would expect for each parameter a matrix depending on the parameter master -> most things need to be checked in RSparams, where I can extract SexDep, DensDep, etc.
    # Check Position and NbOfPositions
    patternPositions <-  "^\"?(([0-9]+-)?[0-9]+,)*([0-9]+-)?[0-9]+\"?$"

    isMatch <- grepl(patternPositions, object@Positions)
    if (!all(isMatch) && object@Positions[isMatch==FALSE] != "random") {
        msg <- c(msg, "In SMSTraits(): Positions must be either a comma-separated list of integer ranges, or random.")
    }
    if (any(is.na(object@NbOfPositions[object@Positions == "random"])) || object@NbOfPositions[object@Positions == "random"] <= 0){
        msg <- c(msg, "NbrOfPositions must be set to a strictly positive integrer.")
    } else if (!is.na(object@NbOfPositions[object@Positions != "random"])) {
        msg <- c(msg, "In SMSTraits(): if Positions is not random NbrOfPositions must be not be set (NA).")
    }

    # Check ExpressionType must be additive or average
    if (!all(object@ExpressionType %in% c("additive", "average"))) {
        msg <- c(msg, "In SMSTraits(): ExpressionType must be either additive or average.")
    }

    # Check InitialDistribution and InitialParameter: Distribution must be uniform or normal
    if (is.null(object@InitialDistribution) || !all(object@InitialDistribution %in% c("uniform", "normal"))) {
        msg <- c(msg,"In SMSTraits(): InitialDistribution must be either uniform or normal.")
    }
    # Check InitialParameters: must be numeric values
    if (any(!is.numeric(object@InitialParameters)) || is.null(object@InitialParameters)) {
        msg <- c(msg, "In SMSTraits(): InitialParameters must be provided.")
    }

    # Check IsInherited: must be TRUE or FALSE
    if (!is.logical(object@IsInherited) || any(is.na(object@IsInherited))) {
        msg <- c(msg, "In SMSTraits(): IsInherited must be either TRUE or FALSE." )
    }
    # Check MutationRate;
    if (object@MutationRate[object@IsInherited == "TRUE"] < 0.0 || object@MutationRate[object@IsInherited == "TRUE"] > 1.0) {
        msg <- c(msg, "In SMSTraits(): MutationRate must be between 0.0 and 1.0.")
    }   else if (any(!is.na(object@MutationRate[object@IsInherited == "FALSE"]))) {
        msg <- c(msg, "In SMSTraits(): If isInherited if off, mutationRate must be blank (NA).")
    }

    # Check MutationDistribution
    # in all places, where the isInherited matrix is TRUE, MutationDistribution must be uniform or normal
    if (any(object@MutationDistribution[object@IsInherited == "TRUE"] != "uniform" && object@MutationDistribution[object@IsInherited == "TRUE"] != "normal")) {
        msg <- c(msg, "In SMSTraits(): MutationDistribution must be either uniform or normal if they are inherited.")
    }
    # in all places where isInherited matrix is FALSE, MutationDistribution must be NA
    if (any(!is.na(object@MutationDistribution[object@IsInherited == "FALSE"]))) {
        msg <- c(msg, "In SMSTraits(): If isInherited is off, MutationDistribution must be blank (NA).")
    }

    # Check MutationParameters
    # in all places, where the isInherited matrix is TRUE, MutationParameters must be numeric
    if (any(!is.numeric(object@MutationParameters[object@IsInherited == "TRUE"])) || any(is.na(object@MutationParameters[object@IsInherited == "TRUE"]))){
        msg <- c(msg, "In SMSTraits(): MutationParameters must be numeric if they are inherited.")
    }
    # in all places where isInherited matrix is FALSE, MutationParameters must be NA
    if (any(!is.na(object@MutationParameters[object@IsInherited == "FALSE"]))) {
        msg <- c(msg, "In SMSTraits(): If isInherited is off, MutationParameters must be blank (NA).")
    }

    # Check OutputValues
    if (object@OutputValues != TRUE && object@OutputValues != FALSE) {
        msg <- c(msg, "In SMSTraits(): OutputValues for emigration traits must be either TRUE or FALSE.")
    }

    if (is.null(msg)) TRUE else msg
})
setMethod("initialize", "SMSTraitsParams", function(.Object, ...) {
    this_func = "SMSTraits(): "
    args <- list(...)
    .Object <- callNextMethod()
    if ( length(args) == 0 ) {
        validObject(.Object)
    }
    .Object
})
setMethod("show", "SMSTraitsParams", function(object){
})

### SUBCLASS TRAITSPARAMS

# define this ClassUnion so that the 'Neutral' slot in the genetic class 'GeneticParams' can be FALSE if neutral genetics should not be modelled
setClassUnion("NeutralSlot", c("logical", "NeutralTraitsParams"))
setClassUnion("GeneticLoadSlot", c("logical", "GeneticLoadParams"))
setClassUnion("EmigrationTraitsSlot", c("logical", "EmigrationTraitsParams"))
setClassUnion("SettlementTraitsSlot", c("logical", "SettlementTraitsParams"))
setClassUnion("CRWTraitsSlot", c("logical", "CRWTraitsParams"))
setClassUnion("SMSTraitsSlot", c("logical", "SMSTraitsParams"))
setClassUnion("KernelTraitsSlot", c("logical", "KernelTraitsParams"))


#' Set genetic traits structure
#'
#' @return a parameter object of class "TraitsParams"
#' @author Jette Reeg
#' @name Traits
#' @export Traits
Traits <- setClass("TraitsParams", slots = c(Neutral = "NeutralSlot",
                                             GeneticLoad = "GeneticLoadSlot",
                                             EmigrationGenes = "EmigrationTraitsSlot",
                                             SettlementGenes = "SettlementTraitsSlot",
                                             CRWGenes = "CRWTraitsSlot",
                                             SMSGenes = "SMSTraitsSlot",
                                             KernelGenes = "KernelTraitsSlot"
                                            )
                       , prototype = list(Neutral = FALSE, # NeutralTraits(),
                                          GeneticLoad = FALSE, # GeneticLoadTraits(), # could this also be a list of multiple GeneticLoads?
                                          EmigrationGenes = FALSE, # EmigrationTraits(), # NULL E_D0, E_Alpha, E_Beta
                                          SettlementGenes = FALSE, # SettlementTraits(), # NULL, S_S0, S_Alpha, S_Beta
                                          CRWGenes = FALSE, # CRWTraits(), # NULL, # CRW_STEPLENGTH, CRW_STEPCORRELATION
                                          SMSGenes = FALSE, # SMSTraits(), # NULL, # SMS_BETADB
                                          KernelGenes = FALSE # KernelTraits() # NULL # KERNEL_MEANDIST_1, KERNEL_MEANDIST_2, KERNEL_PROBABILITY
                       ))
setValidity("TraitsParams", function(object) {
    msg <- NULL

    # Check Neutral
    if (!is.logical(object@Neutral) && !is(object@Neutral,"NeutralTraitsParams")) {
        msg <- c(msg, "In Traits(): Neutral must be of class NeutralTraitsParams.")
    }
    # Check GeneticLoad
    if (!is.logical(object@GeneticLoad) && !is(object@GeneticLoad, "GeneticLoadParams")) {
        msg <- c(msg, "In Traits(): GeneticLoad must be of class GeneticLoadParams.")
    }
    # Check EmigrationGenes
    if (!is.logical(object@EmigrationGenes) && !is(object@EmigrationGenes, "EmigrationTraitsParams")) {
        msg <- c(msg, "In Traits(): EmigrationGenes must be of class EmigrationTraitsParams.")
    }
    # Check SettlementGenes
    if (!is.logical(object@SettlementGenes) && !is(object@SettlementGenes, "SettlementTraitsParams")) {
        msg <- c(msg, "In Traits(): SettlementGenes must be of class SettlementTraitsParams.")
    }
    # Check CRWGenes
    if (!is.logical(object@CRWGenes) && !is(object@CRWGenes, "CRWTraitsParams")) {
        msg <- c(msg, "In Traits(): CRWGenes must be of class CRWTraitsParams.")
    }
    # Check SMSGenes
    if (!is.logical(object@SMSGenes) && !is(object@SMSGenes, "SMSTraitsParams")) {
        msg <- c(msg, "In Traits(): SMSGenes must be of class SMSTraitsParams.")
    }
    # Check KernelGenes
    if (!is.logical(object@KernelGenes) && !is(object@KernelGenes, "KernelTraitsParams")) {
        msg <- c(msg, "In Traits(): KernelGenes must be of class KernelTraitsParams.")
    }
    if (is.null(msg)) TRUE else msg
})
setMethod("initialize", "TraitsParams", function(.Object, ...) {
    this_func = "Traits(): "
    args <- list(...)
    .Object <- callNextMethod()
    if ( length(args) == 0 ) {
        validObject(.Object)
    }
    .Object
})
setMethod("show", "TraitsParams", function(object){
    if(class(object@Neutral) == "NeutralTraitsParams") print(object@Neutral)
    if(class(object@GeneticLoad) == "GeneticLoadParams") print(object@GeneticLoad)
    if(class(object@EmigrationGenes) == "EmigrationTraitsParams") print(object@EmigrationGenes)
    if(class(object@SettlementGenes) == "SettlementTraitsParams") print(object@SettlementGenes)
    if(class(object@CRWGenes) == "CRWTraitsParams") print(object@CRWGenes)
    if(class(object@SMSGenes) == "SMSTraitsParams") print(object@SMSGenes)
    if(class(object@KernelGenes) == "KernelTraitsParams") print(object@KernelGenes)
})


#' Set Genetics parameters
#'
#' @description Set genetics parameters\cr
#'
#' Controls heritability and evolution of traits (if inter-individual variability is enabled (\code{IndVar=TRUE}) for at least one (dispersal) trait).
#' Provides control over the genome size, the number of loci, the recombination rate (if the species is diploid) and a
#' flexible mapping of traits to chromosomes, allowing linkage, pleiotropy and neutral alleles to be incorporated. It is also possible to model
#' neutral alleles when no adaptive traits are present.
#'
#' @usage Genetics(GenomeSize = 10,
#'          ChromosomeEnds = 1, RecombinationRate = 0.0,
#'          OutputGeneValues = FALSE, OutputNeutralStatistics = FALSE,
#'          OutputFstatsWeirCockerham = FALSE, OutputFstatsWeirHill = FALSE,
#'          OutputStartGenetics = NULL, OutputInterval = NULL,
#'          PatchList = NULL, NbrPatchToSample = NULL,
#'          nIndividualsToSample = NULL, Stages = NULL, Traits = Traits()
#'          )
#'
#' @param GenomeSize Maximum size of genome (number of loci)
#' @param ChromosomeEnds Where the genome is split into chromosomes, if empty
#' assumed one chromosome is equal to GenomeSize. These areas recombine with a
#' probability of 0.5.
#' @param RecombinationRate Recombination rate (through chromosomal crossover)
#' across the whole genome (in addition to the chromosomeEnds above).
#' @param OutputGeneValues Output the values of all alleles for all genes of all
#' sampled individuals. Does not output the resulting trait values: mean and SD
#' of dispersal and genetic fitness traits are output in the TraitsXPatch,
#' TraitsXCell and/or TraitsXrow output files. Enables the geneValues output
#' files.
#' @param OutputFstatsWeirCockerham Calculate F-statistics (including global and
#' per-locus estimates) according to Weir & Cockerham (1984)'s method-of-moments
#' approach. Enables the neutralGenetics and perLocusNeutralGenetics output files.
#' @param OutputFstatsWeirHill Calculate F-statistics calculated according to the
#' estimators of Weir & Hill (2002), including global estimates corrected for
#' unequal sample sizes, population- (i.e. patch-) specific estimates, and pairwise
#' estimates. Enables the neutralGenetics and pairwisePatchNeutralGenetics output files.
#' @param OutputStartGenetics Which year should RangeShifter start to produce the output files listed above?
#' @param OutputInterval How frequently to output genetic output, including gene values and neutral statistics.
#' @param PatchList Which patches are to be sampled for output.  Patches can be
#' specified according to their patch number, as per the patch layer in a patch-based
#' model. Or sampled randomly or all patches can be chosen. In a cell-based landscape
#' random is the only option with number of patches (=cells) specified.
#' @param NbrPatchToSample If PatchList=random or random_occupied then this specifies
#' the number of patches to sample randomly. Random: The chosen sample patches remain
#' the same throughout the simulation, i.e. do not vary between years or replicates
#' unless artificially generated landscape that is generated afresh between replicates.
#' Random_occupied: patches are re-sampled every generation among all patches containing at least 1 individual.
#' @param nIndividualsToSample The number of individuals to sample in a patch. If nInds < nIndividualsToSample then sampled individuals = nInds
#' @param Stages The age stages to sample from.
#' @param Traits The genetic traits to be modelled.

#' @details TBD

#'
#'
# #' @references \insertAllCited{}
#' @return a parameter object of class "GeneticsParams"
#' @author Jette Reeg
#' @name Genetics
#' @export Genetics
Genetics <- setClass("GeneticsParams", slots = c(GenomeSize = "integer_OR_numeric",
                                                 ChromosomeEnds = "integer_OR_numeric", # NULL or vector
                                                 RecombinationRate = "integer_OR_numeric", # NULL or numeric
                                                 OutputGeneValues = "logical",
                                                 OutputFstatsWeirCockerham = "logical",
                                                 OutputFstatsWeirHill = "logical",
                                                 OutputStartGenetics = "integer_OR_numeric", # positive integer if any output is TRUE or NULL
                                                 OutputInterval = "integer_OR_numeric",
                                                 PatchList = "character_OR_integer", # vector of integers or a string
                                                 NbrPatchToSample = "integer_OR_numeric", # NULL or integer
                                                 nIndividualsToSample = "character_OR_integer", # character or integer
                                                 Stages = "character_OR_integer", # vector
                                                 Traits = "TraitsParams")
                                     , prototype = list(GenomeSize = 0L,
                                                        ChromosomeEnds = 0L, # NULL or vector
                                                        RecombinationRate = 0.0, # NULL or numeric
                                                        OutputGeneValues = FALSE,
                                                        OutputFstatsWeirCockerham = FALSE,
                                                        OutputFstatsWeirHill = FALSE,
                                                        OutputStartGenetics = 0L, # positive integer if any output is TRUE or NULL
                                                        OutputInterval = 0L,
                                                        PatchList = "all", # vector or string
                                                        NbrPatchToSample = 0L, # NULL or integer
                                                        nIndividualsToSample = "all", # NULL or integer
                                                        Stages = "all", # vector
                                                        Traits = Traits())

)
setValidity('GeneticsParams', function(object){
    msg <- NULL
    if (is.null(msg)) TRUE else msg
})
setMethod('initialize', 'GeneticsParams', function(.Object, ...) {
    this_func = "Genetics(): "
    args <- list(...)
    .Object <- callNextMethod()
    if ( length(args) == 0 ) {
        validObject(.Object)
    }
    .Object
    })
setMethod("show", "GeneticsParams", function(object){
    cat(" Genetics: \n")
    cat("   Genome size =", object@GenomeSize, "\n ")
    if(length(object@ChromosomeEnds) > 0){
        cat("   Genome is slitted into chromosomes at =", object@ChromosomeEnds, "\n ")
    } else cat( "   Genome is not slitted into chromosomes \n")

    cat("   Recombination rate: ", object@RecombinationRate, "\n")
    cat("   Output genetic values: ", object@OutputGeneValues, "\n")
    cat("   Output Fstats after Weir Cockerham: ", object@OutputFstatsWeirCockerham, "\n")
    cat("   Output Fstats after Weir Hill: ", object@OutputFstatsWeirHill, "\n")
    cat("   Start genetic output at year: ", object@OutputStartGenetics, "and output every ",object@OutputInterval ," year \n")
    cat("   Patches to sample: ", object@PatchList, "\n")
    if(object@PatchList=="random" || object@PatchList=="random_occupied"){
        cat("   Number of patches to sample: ", object@NbrPatchToSample, "\n")
    }
    cat("   Number of individuals to sample: ", object@nIndividualsToSample, "\n")
    cat("   Stages to sample: ", object@Stages, "\n")
    print(object@Traits)
})
