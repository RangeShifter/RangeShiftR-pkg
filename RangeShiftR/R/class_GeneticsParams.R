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

### CLASS GENETICSPARAMS

# from RS 'Traits' file

### SUBCLASS NEUTRALTRAITS

#' Set up structure of neutral traits
#'
#' A method to specify neutral traits in the genetic module.
#'
#' #' @usage NeutralTraits(Positions = "random", NbOfPositions = 10,
#' Positions = "random", # "random" or list of integer values
#' NbOfPositions = 10, # numeric, only of positions random
#' InitialDistribution = NULL, # uniform (neutral + dispersal), normal (dispersal)
#' InitialParameters = 2, # neutral traits: only max value; dispersal: two values: either min/max oder mean+sd, not applicable for genetic load
#' MutationDistribution = "KAM", # neutral: "KAM" or "SSM", genetic load: "gamma", "uniform", "normal", "negExp", dispersal: uniform or normal
#' MutationParameters = 2, # single value or 2 values
#' MutationRate = 0.0, # numeric
#' OutputValues = FALSE)
#'
#' @param Positions Loci positions coding for trait within genome. Must be in the
#' range 0-(\code{GenomeSize}-1), specified in \code{\link[RangeShiftR]{Genetics}}. Positions can overlap across
#' traits - there will be no pleiotropy, but this will influence genetic linkage.,
#' @param NbOfPositions Only specify if above is set to \code{"random"}, else must be blank (\code{NULL})
#' @param InitialDistribution Distribution from which to draw initial allele values from.
#' If \code{uniform}. Initialise with random characters between 0 – \code{max}. Note that possible values start at 0, so \code{max=0}
#' specifies a monomorphic initial population.
#' @param InitialParameters Maximal value for the uniform distribution.
#' @param MutationDistribution Distribution for mutations to draw from. Can be either \code{"KAM"} or \code{"SSM"}. \code{KAM} (k-alleles model) is randomly
#' draw a value between 0 and \code{max} (see \code{MutationParameters}).
#' \code{SSM} (single-step mutation) is to move in a stepwise manner, A to B, B to C.
#' @param MutationParameters Parameters for the above distribution: maximal value for \code{KAM} or \code{SSM} (cannot exceed 255)
#' @param MutationRate Mutation rate applicable to this type of loci. Must be between 0.0 and 1.0
#' @param OutputValues If OutputGeneValues in \code{\link[RangeShiftR]{Genetics}} is
#' enabled, should allele values for this gene be written to output? Ignored if OutputGeneValues is set to \code{FALSE}.
#'
#' @details
#'
#' Neutral trait does not have any phenotypic effect during the simulation. It is used to compute F-statistics and other measures of neutral variation.
#' For neutral trait you must specify:
#'
#' -	The number and positions of genes controlling the trait: \code{Positions} and \code{NbOfPositions} \cr
#' -	A distribution to sample initial values from: \code{InitialDistribution} and \code{InitialParameters} \cr
#' -	A mutation rate for all genes controlling the trait. \code{MutationRate} \cr
#' -	A distribution to sample mutations from. \code{MutationDistribution} and \code{MutationParameters} \cr
#'
#' The user specifies the number of possible alleles for neutral loci (up to 256), via the maximum parameter of the mutation distribution.
#' Initial values are either identical for all sites (equal to the max value) or sampled in a uniform distribution (between 0 and the maximum value).
#' Mutations are either sampled in a uniform distribution between 0 and the max parameter (k-allele model, KAM, \insertCite{peng2012}{RangeShiftR}) or
#' added as increments (random -1 or +1 changes) of the previous value (symmetric stepwise model, SSM, \insertCite{peng2012}{RangeShiftR}).
#'
#' Dominance values and inheritance are not applicable for neutral traits.
#'
#'@references
#'         \insertAllCited{}
#' @return a parameter object of class "NeutralTraitsParams"
#' @author Jette Reeg
#' @name NeutralTraits
#' @export NeutralTraits
NeutralTraits<- setClass("NeutralTraitsParams", slots = c(Positions = "ANY", # vector of numbers or "random"
                                                         NbOfPositions = "ANY", # integer value or NULL
                                                         InitialDistribution = "character", # uniform
                                                         InitialParameters = "integer_OR_numeric", # max value
                                                         MutationDistribution = "character", # KAM or SSM
                                                         MutationParameters = "integer_OR_numeric", # max
                                                         MutationRate = "integer_OR_numeric", # float
                                                         OutputValues = "logical")
                   , prototype = list(Positions = "random", # "random" or list of integer values
                                      NbOfPositions = 10, # numeric, only if positions random
                                      InitialDistribution = NULL, # uniform
                                      InitialParameters = 2, # neutral traits: only max value;
                                      MutationDistribution = "KAM", # neutral: "KAM" or "SSM"
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
    if (object@InitialDistribution != "uniform") {
        msg <- c(msg,"InitialDistribution must be uniform for the neutral trait.")
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
        msg <- c(msg, "For neutral trait InitialDistribution must be uniform.")

    }
    # Check mutation rate
    if (!is.null(object@MutationRate) && (any(object@MutationRate < 0) || any(object@MutationRate > 1))) {
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
    if(object@OutputValues) cat("     Allel values for gene is written to output.  \n")
})


### SUBCLASS GENETICLOADTRAITS

#' Set genetic structure for genetic fitness traits
#'
#' @description
#' A method to specify genetic fitness traits in the genetic module. You can specify up to 5 genetic loads.
#'
#'
#' @usage GeneticLoadTraits(NbGeneticLoads = 1, Positions = list("random"), NbOfPositions = 10,
#' InitialDistribution = "normal", InitialParameters = matrix(c(0.5,0.1), nrow=1),
#' InitialDomDistribution = "normal", InitialDomParameters = matrix(c(0.5,0.1), nrow=1),
#' DominanceDistribution = "normal", DominanceParameters = matrix(c(0.5,0.1), nrow=1),
#' MutationDistribution = "normal", MutationParameters = matrix(c(0.5,0.2), nrow=1),
#' MutationRate = 0.0001, OutputValues = FALSE)
#'
#' @param NbGeneticLoads Number of genetic loads
#' @param Positions Loci positions coding for that trait within genome. Should be provided as a list of strings (if random) and/or vectors of integers (if not random)
#' @param NbOfPositions Only specify when the \code{Positions} of a genetic load trait are set to \code{"random"}, else must be blank (\code{NULL})
#' @param InitialDistribution Distribution from which to draw initial allele values from. Can be \code{gamma}, \code{uniform}, \code{normal}, \code{negExp}. Should be provided as a vector of strings if \code{NbGeneticLoads > 1}
#' @param InitialParameters Parameters for the initial distribution: You must provide two colums for \code{uniform}, \code{normal} and \code{gamma} distributions:
#' min and max (\code{uniform}), mean and sd (\code{normal}) or shape and scale (\code{gamma}) or one column for \code{negExp}: mean
#' If genetic loads have different \code{InitialDistribution} and one requires two columns you need to set the second value to \code{NA} in case of \code{negExp} distribution.
#' Each row in the matrix corresponds to a genetic load trait.
#' @param InitialDomDistribution Distribution from which to draw initial dominance values. Can be \code{gamma}, \code{uniform}, \code{normal}, \code{negExp}, \code{scaled}. Should be provided as a vector of strings if \code{NbGeneticLoads > 1}
#' @param InitialDomParameters Parameters for the initial dominance distribution: You must provide two colums for \code{uniform}, \code{normal} and \code{gamma} distributions:
#' min and max (\code{uniform}), mean and sd (\code{normal}) or shape and scale (\code{gamma}) or one column for \code{negExp}: mean
#' If genetic loads have different \code{InitialDomParameters} and one requires two columns you need to set the second value to \code{NA} in case of \code{negExp} distribution.
#' Each row in the matrix corresponds to a genetic load trait.
#' @param DominanceDistribution Distribution of dominance values. Can be \code{gamma}, \code{uniform}, \code{normal}, \code{negExp}, \code{scaled}. Should be provided as a vector of strings if \code{NbGeneticLoads > 1}
#' @param DominanceParameters Parameters for the dominance distribution: You must provide two colums for \code{uniform}, \code{normal} and \code{gamma} distributions:
#'  min and max (\code{uniform}), mean and sd (\code{normal}) or shape and scale (\code{gamma}) or one column for \code{negExp}, \code{scaled}: mean
#' If genetic loads have different \code{DominanceDistribution} and one requires two columns you need to set the second value to \code{NA} in case of \code{negExp} or \code{scaled} distribution.
#' Each row in the matrix corresponds to a genetic load trait.
#' @param MutationDistribution Distribution for mutations to draw from. Can be \code{gamma}, \code{uniform}, \code{normal}, \code{negExp}. Should be provided as a vector of strings if \code{NbGeneticLoads > 1}
#' @param MutationParameters Parameters for the mutation distribution: You must provide two colums for \code{uniform}, \code{normal} and \code{gamma} distributions: min and max (\code{uniform}), mean and sd (\code{normal})
#'  or shape and scale (\code{gamma}) or one column for \code{negExp}: mean
#' If genetic loads have different \code{DominanceDistribution} and one require two columns you need to set the second value to \code{NA} in case of \code{negExp} distribution.
#' @param MutationRate Mutation rate applicable to this type of loci. Must be between \code{0.0} and \code{1.0}. Should be provided as a vector if multiple genetic loads are specified.
#' @param OutputValues If OutputGeneValues in GeneticsFile is enabled, should allele values for this gene be written to output? Ignored if OutputGeneValues is set to \code{FALSE}.
#' Should be provided as a vector if multiple genetic loads are specified.
#'
#' @details
#' The alleles of genetic fitness traits represent deleterious mutations which combined expression reduce the viability of juveniles,
#' i.e. the genetic load. Immediately after birth, all newly born individuals are checked for viability via a Bernoulli trial.
#' The probability of an individual passing the test and surviving birth is equal to its genetic fitness.
#' Genetic fitness is 1 by default, but every allele reduces this value by an amount that depends on
#' its selection coefficient \eqn{s}, and (for diploid systems) the value of its dominance coefficient \eqn{h} relative
#' to that of the other allele it is paired with.
#'
#' More precisely, the genetic fitness \eqn{W} of an individual is the product of the contributions \eqn{w} of each genetic load locus.
#' The contribution \ifelse{html}{\out{w<sub>i</sub>}}{\eqn{w_i}}  of locus \eqn{i} with alleles \eqn{A} and \eqn{B} is:
#'
#' \ifelse{html}{\out{&emsp;&emsp; w<sub>i</sub> = 1 - h<sub>i</sub>s<sub>A</sub> - (1 - h<sub>i</sub>)s<sub>B</sub>}}{\deqn{w_i = 1 - h_i s_A - (1 - h_i) s_b}}
#'
#' \ifelse{html}{\out{&emsp;&emsp; h<sub>i</sub> = h<sub>A</sub> / (h<sub>A</sub> + h<sub>B</sub>)}}{\deqn{h_i = h_A / ( h_A + h_B )}}
#'
#' Selection and dominance coefficients for new alleles (mutations) are drawn from distributions specified by the user:
#' either a uniform, normal, gamma, or negative exponential, and are not additive.
#' Dominance coefficients can additionally be sampled from a scaled uniform distribution between zero and a maximum value
#' that depends on the selection coefficient. This maximum value is equal to \ifelse{html}{\out{e<sup>-ks<sub>i</sub></sup>}}{\deqn{exp(-ks_i)}} where \ifelse{html}{\out{s<sub>i</sub>}}{\eqn{s_i}}is the selection coefficient
#' for the locus and
#'
#' \ifelse{html}{\out{&emsp;&emsp; k = -log (2h<sub>D</sub>)/s<sub>D</sub>;&ast;h<sub>D</sub>}}{\deqn{k = -log(2h_D) / s_D}}
#'
#' \ifelse{html}{\out{h<sub>D</sub>}}{\eqn{h_D}} is the desired mean dominance coefficient; \ifelse{html}{\out{s<sub>D</sub>}}{\eqn{s_D}} is the mean of the selection coefficient distribution,
#' calculated after the parameters entered for the (selection coefficients) mutation distribution.
#'
#' While selection coefficients should typically be zero or positive, to represent the effect of deleterious mutations,
#' negative values up to -1 are allowed and may arise if the mutation distribution specified by the user allows it.
#' In this case, negative values would represent (universally) beneficial mutations that counteract the effect of genetic load.
#' The total genetic fitness is however bounded to 1.
#'
#' To allow more flexibility, there can be up to 5 such genetic load traits, each with a potentially different distribution
#' of mutations and/or dominance distribution and associated parameters.
#'
#' The expression type of genetic load traits is always multiplicative.
#'
#' Initial values and inheritance is not applicable for genetic load traits.
#'
#' @return a parameter object of class "GeneticLoadParams"
#' @author Jette Reeg
#' @name GeneticLoadTraits
#' @export GeneticLoadTraits
GeneticLoadTraits<- setClass("GeneticLoadParams", slots = c(
    NbGeneticLoads = "integer_OR_numeric", # number of genetic loads
    Positions = "list",# "random" or list of integer values
    NbOfPositions = "ANY", # numeric, only where positions are random; otherwise NA
    InitialDistribution = "ANY", #‘gamma’, ‘uniform’, ‘normal’,‘negExp’
    InitialParameters = "ANY", # 2 values for min/max, mean/sd, shape/scale or one value: mean
    InitialDomDistribution = "ANY", #  ‘gamma’, ‘uniform’, ‘normal’, ‘negExp’, ‘scaled’
    InitialDomParameters = "ANY", # 2 values for min/max, mean/sd, shape/scale or one value: mean
    DominanceDistribution = "character", # ‘gamma’, ‘uniform’, ‘normal’, ‘negExp’, ‘scaled’ # character vector
    DominanceParameters = "matrix", # 2 values for min/max, mean/sd, shape/scale or one value: mean
    MutationDistribution = "character", # ‘gamma’, ‘uniform’, ‘normal’,‘negExp’
    MutationParameters = "matrix", #  2 values for min/max, mean/sd, shape/scale or one value: mean
    MutationRate = "numeric", # float
    OutputValues = "logical")
    , prototype = list(
        NbGeneticLoads = 1L,
        Positions = list("random"),
        NbOfPositions = 2L,
        InitialDistribution = NULL, # "normal",
        InitialParameters = NULL, # matrix(c(0.5,0.1), nrow=1),
        InitialDomDistribution = NULL, # "normal",
        InitialDomParameters = NULL, # matrix(c(0.5,0.1), nrow=1),
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
            isNumeric <- sapply(object@Positions, is.numeric)
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

    # Check InitialDistribution given or NA
    if (!is.null(object@InitialDistribution)){
        # if it is given, it should be character
        if (class(object@InitialDistribution) != "character") {
            msg <- c(msg, "InitialDistribution must be a character vector.")
        }
        # and also InitialParameters should be numeric
        if (class(object@InitialParameters)[1] != "matrix") {
            msg <- c(msg, "InitialParameters must be a matrix.")
        }
        if(length(object@InitialDistribution) != object@NbGeneticLoads) {
            msg <- c(msg, "If you want to have initial allel distribution for one genetic load, you must provide it for all genetic load. If you don't want to set the initial allel value for some genetic loads, set it to NA for those.")
        } else if (nrow(object@InitialParameters) != object@NbGeneticLoads) {
            msg <- c(msg, "If you have set InitialDistributions for at least one genetic load you must provide the InitialParameters for each genetic load. Use one row for each genetic load.")
        } else {
            if (any(object@InitialDistribution == "normal")) { # if any distribution is normal
                # two values for mean and sd
                if (ncol(object@InitialParameters) !=2 || # if DominanceParameters has not 2 columns OR
                    any(!is.numeric(object@InitialParameters[object@InitialDistribution=="normal"])) || # if entries are not numeric
                    any(is.na(object@InitialParameters[object@InitialDistribution=="normal"]))) { # if entries are NA
                    msg <- c(msg,"For a normal initial distribution, InitialParams must provide two values for mean (first column) and sd (second column)")
                }
            }
            if (any(object@InitialDistribution == "gamma")) {
                # two values for shape and scale
                if (ncol(object@InitialParameters) !=2 || # if DominanceParameters has not 2 columns OR
                    any(!is.numeric(object@InitialParameters[object@InitialDistribution=="gamma"])) || # if entries are not numeric
                    any(is.na(object@InitialParameters[object@InitialDistribution=="gamma"]))) { # if entries are NA
                    msg <- c(msg,"For a gamma initial distribution, InitialParams must provide two values for shape (first column) and scale (second column)")
                }
            }
            if (any(object@InitialDistribution == "uniform")) {
                # two values for min and max
                if (ncol(object@InitialParameters) !=2 || # if DominanceParameters has not 2 columns OR
                    any(!is.numeric(object@InitialParameters[object@InitialDistribution=="uniform"])) || # if entries are not numeric
                    any(is.na(object@InitialParameters[object@InitialDistribution=="uniform"]))) { # if entries are NA
                    msg <- c(msg,"For a uniform initial distribution, InitialParams must provide two values for min (first column) and max (second column)")
                }
            }
            if (all(object@InitialDistribution == "negExp")) { # if it is only negExp or scaled
                # one value for mean and one NA
                if (ncol(object@InitialParameters) !=1 || # if DominanceParameters has more than 1 column
                    !is.numeric(object@InitialParameters) ||
                    any(is.na(object@InitialParameters))
                ) {
                    msg <- c(msg,"For negative exponential and scaled initial distribution, InitialParams must provide only one column for mean.")
                }
            } else{
                if (any(object@InitialDistribution == "negExp")) { # if only some are scaled or negative exponential
                    # one value for mean and one NA if other distributions need 2 values
                    if (ncol(object@InitialParameters) !=2 || # if DominanceParameters has not 2 columns OR
                        !is.numeric(object@InitialParameters) || # if entries are not numeric
                        !all(is.na(object@InitialParameters[object@InitialDistribution=="negExp",2])) || # second column is not NA
                        any(is.na(object@InitialParameters[object@InitialDistribution=="negExp",1])) # first column is NA
                    ) {
                        msg <- c(msg,"For the negative exponential initial distribution, InitialParameters must provide only one value for mean (first column) and the second column need to be NA if other genetic loads use other initial distributions.")
                    }
                }
            }

            # if any InitialDistribution is not 'normal', 'gamma', 'uniform', 'negExp' or 'scaled' OR NA
            if (!all(object@InitialDistribution %in% c("uniform", "normal", "gamma", "negExp", NA))) {
                msg <- c(msg, "InitialDistribution must be either normal, gamma, uniform, negExp, scaled or NA (if not initialized for a genetic load) for genetic load traits.")
            }
        }
    }

    # Check InitialDomDistribution
    if (!is.null(object@InitialDomDistribution)){
        if (class(object@InitialDomDistribution) != "character") {
            msg <- c(msg, "InitialDistribution must be a character vector.")
        }
        # and also InitialParameters should be numeric
        if (class(object@InitialDomParameters)[1] != "matrix") {
            msg <- c(msg, "InitialParameters must be a matrix.")
        }
        if(length(object@InitialDomDistribution) != object@NbGeneticLoads) {
            msg <- c(msg, "If you want to have initial dominance distribution for one genetic load, you must provide it for all genetic load. If you don't want to set the initial dominance value for some genetic loads, set it to NA for those.")
        } else if (nrow(object@InitialDomParameters) != object@NbGeneticLoads) {
            msg <- c(msg, "If you have set InitialDomDistributions for at least one genetic load you must provide the InitialDomParameters for each genetic load. Use one row for each genetic load.")
        } else {
            if (any(object@InitialDomDistribution == "normal")) { # if any distribution is normal
                # two values for mean and sd
                if (ncol(object@InitialDomParameters) !=2 || # if DominanceParameters has not 2 columns OR
                    any(!is.numeric(object@InitialDomParameters[object@InitialDomDistribution=="normal"])) || # if entries are not numeric
                    any(is.na(object@InitialDomParameters[object@InitialDomDistribution=="normal"]))) { # if entries are NA
                    msg <- c(msg,"For a normal initial dominance distribution, InitialDomParameters must provide two values for mean (first column) and sd (second column)")
                }
            }
            if (any(object@InitialDomDistribution == "gamma")) {
                # two values for shape and scale
                if (ncol(object@InitialDomParameters) !=2 || # if DominanceParameters has not 2 columns OR
                    any(!is.numeric(object@InitialDomParameters[object@InitialDomDistribution=="gamma"])) || # if entries are not numeric
                    any(is.na(object@InitialDomParameters[object@InitialDomDistribution=="gamma"]))) { # if entries are NA
                    msg <- c(msg,"For a gamma initial dominance distribution, InitialDomParameters must provide two values for shape (first column) and scale (second column)")
                }
            }
            if (any(object@InitialDomDistribution == "uniform")) {
                # two values for min and max
                if (ncol(object@InitialDomParameters) !=2 || # if DominanceParameters has not 2 columns OR
                    any(!is.numeric(object@InitialDomParameters[object@InitialDomDistribution=="uniform"])) || # if entries are not numeric
                    any(is.na(object@InitialDomParameters[object@InitialDomDistribution=="uniform"]))) { # if entries are NA
                    msg <- c(msg,"For a uniform initial dominance distribution, InitialDomParameters must provide two values for min (first column) and max (second column)")
                }
            }
            if (all(object@InitialDomDistribution == "negExp" || object@InitialDomDistribution == "scaled")) { # if it is only negExp or scaled
                # one value for mean and one NA
                if (ncol(object@InitialDomParameters) !=1 || # if DominanceParameters has more than 1 column
                    !is.numeric(object@InitialDomParameters) ||
                    any(is.na(object@InitialDomParameters))
                ) {
                    msg <- c(msg,"For negative exponential and scaled initial dominance distribution, InitialParameters must provide only one column for mean.")
                }
            } else{
                if (any(object@InitialDomDistribution == "scaled" || object@InitialDomDistribution == "negExp")) { # if only some are scaled or negative exponential
                    # one value for mean and one NA if other distributions need 2 values
                    if (ncol(object@InitialDomParameters) !=2 || # if DominanceParameters has not 2 columns OR
                        !is.numeric(object@InitialDomParameters) || # if entries are not numeric
                        !all(is.na(object@InitialDomParameters[object@InitialDomDistribution=="scaled",2])) || # second column is not NA
                        !all(is.na(object@InitialDomParameters[object@InitialDomDistribution=="negExp",2])) || # second column is not NA
                        any(is.na(object@InitialDomParameters[object@InitialDomDistribution=="scaled",1])) || # first column is NA
                        any(is.na(object@InitialDomParameters[object@InitialDomDistribution=="negExp",1])) # first column is NA
                    ) {
                        msg <- c(msg,"For the scaled or negative exponential initial dominance distribution, InitialParameters must provide only one value for mean (first column) and the second column need to be NA if other genetic loads use other initial dominance distributions.")
                    }
                }
            }

            # if any InitialDomDistribution is not 'normal', 'gamma', 'uniform', 'negExp' or 'scaled' OR NA
            if (!all(object@InitialDomDistribution %in% c("uniform", "normal", "gamma", "negExp", "scaled", NA))) {
                msg <- c(msg, "InitialDomDistribution must be either normal, gamma, uniform, negExp, scaled or NA (if not initialized for a genetic load) for genetic load traits.")
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
            # if any DominanceDistribution is not 'normal', 'gamma', 'uniform', 'negExp' or 'scaled' OR NA
            if (!all(object@DominanceDistribution %in% c("uniform", "normal", "gamma", "negExp", "scaled"))) {
                msg <- c(msg, "DominanceDistribution must be either normal, gamma, uniform, negExp, or scaled for genetic load traits.")
            }
        }
    } else{ # cannot be NULL
        msg <- c(msg, "You must provide the DominanceDistribution for each genetic load.")
    }

    # Check mutation rate
    if(!is.numeric(object@MutationRate) || length(object@MutationRate) != object@NbGeneticLoads){
        msg <- c(msg, "You must provide the mutation rate for each genetic load as a numeric vector.")
    } else {
        if (!is.numeric(object@MutationRate) ||  (any(object@MutationRate < 0.0) || any(object@MutationRate > 1.0))) {
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

    if (!all(object@OutputValues %in% c(TRUE, FALSE))) {
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
    cat("   Genetic Loads: \n")
    cat("     Number of genetic loads: ", object@NbGeneticLoads, "\n")
    for (i in 1:object@NbGeneticLoads){
        cat("     Configuration of genetic load ", i, ": \n")
        if(is.numeric(object@Positions[i])) cat("     Loci positions coding for trait: ", object@Positions[i], "\n")
        if(!is.numeric(object@Positions[i]) && object@Positions[i]=="random") cat("    Loci positions coding for trait randomly chosen with ", object@NbOfPositions[i], " positions\n")
        cat("       Initial allel distribution: ", object@InitialDistribution[i], "\n")
        cat("       Initial allel distribution parameter: ", object@InitialParameters[i,], "\n")
        cat("       Initial dominance distribution: ", object@InitialDomDistribution[i], "\n")
        cat("       Initial dominance parameter: ", object@InitialDomParameters[i,], "\n")
        cat("       Dominance distribution: ", object@DominanceDistribution[i], "\n")
        cat("       Dominance parameter: ", object@DominanceParameters[i,], "\n")
        cat("       Mutation distribution: ", object@MutationDistribution[i], "\n")
        cat("       Mutation parameters: ", object@MutationParameters[i,], "\n")
        cat("       Mutation rate: ", object@MutationRate, "\n")
        if(object@OutputValues) cat("       Allel values for gene is written to output. \n")
    }

})

### SUBCLASS EMIGRATIONTRAITS

#' Set genetic traits structure for emigration traits
#'
#' @description
#'
#' Depending of the settings in \code{\link[RangeShiftR]{Emigration}}, up to three emigration traits evolve:
#'
#' - The probability of emigration \ifelse{html}{\out{D<sub>0</sub>}}{\eqn{D_0}} / \eqn{d} \cr
#' - The slope of the density dependent emigration function (if \code{DensDep=TRUE}) \ifelse{html}{\out{&alpha;}}{\eqn{α}} \cr
#' - The density threshold of the density dependent emigration function (if \code{DensDep=TRUE}) \ifelse{html}{\out{&beta;}}{\eqn{β}} \cr
#'
#' If emigration is sex dependent, \code{SexDep=TRUE}, you must provide details for both sexes.
#'
#' This results in following number of emigration traits that need to be specifed:
#'
#' - 1 entry if emigration is neither density dependent nor sex dependent (\code{DensDep=FALSE} and \code{SexDep=FALSE}): emgiration probability \eqn{d} \cr
#' - 2 entries if emigration is sex dependent (\code{SexDep=TRUE}): emigration probability of females \eqn{d(f)} and males \eqn{d(m)} \cr
#' - 3 entries if emigration probability is only density dependent \code{DensDep=TRUE}:  \ifelse{html}{\out{D<sub>0</sub>}}{\eqn{D_0}}, \ifelse{html}{\out{&alpha;}}{\eqn{α}} and \ifelse{html}{\out{&beta;}}{\eqn{β}} \cr
#' - 6 entries if emigration probability is both density dependent (\code{DensDep=TRUE}) and sex dependent (\code{SexDep=TRUE}): \ifelse{html}{\out{D<sub>0</sub>(f), D<sub>0</sub>(m)}}{\eqn{D_0(f),D_0(m)}}, \ifelse{html}{\out{&alpha;(f), &alpha;(m)}}{\eqn{α(f), α(m)}} and \ifelse{html}{\out{&beta;(f), &beta;(m)}}{\eqn{β(f), β(m)}} \cr
#'
#' The entries of the trait parameters must be provided in the same order as the emigration traits are listed above. If parameters expect a matrix, the rows must match the order of kernel traits listed above.
#' @details
#'
#' Traits set to evolve cannot simultaneously be stage-dependent.
#'
#' The alleles of each trait can be expressed according to an additive model (allele values across all loci are summed) or be averaged.
#'
#' Mutations are additive and can be sampled in either a \code{normal} or \code{uniform} distribution.
#'
#' Initial allele values are sampled in a \code{normal} or \code{uniform} distribution.
#'
#' Emigration traits can also be **not** inherited, that is, allele values are resampled from the initial distribution for every new individual.
#'
#' Dominance values are not applicable for emigration traits.
#'
#' @usage EmigrationTraits(Positions = list("random"), NbOfPositions = 10,
#' ExpressionType = "additive",
#' InitialDistribution = "normal", InitialParameters = matrix(c(0.5,0.1), nrow=1),
#' IsInherited = FALSE,
#' MutationDistribution = "normal", MutationParameters = matrix(c(0.5,0.2), nrow=1),
#' MutationRate = 0.0001, OutputValues = FALSE)
#'
#' @param Positions Loci positions coding for the trait within genome. Should be provided as a list. Entries can either be a string (\code{"random"}) and/or vectors of integers.
#' The length must be equal to the number of required emigration traits (see above) and the sequence must match the sequence of the emigration traits listed above.
#' @param NbOfPositions Only specify when the \code{Positions} of the emigration trait is set to \code{"random"}, else must be blank (\code{NULL}).
#' The length must be equal to the number of required emigration traits (see above) and the sequence must match the sequence of the emigration traits listed above.
#' @param ExpressionType Type of expression for the emigration trait. Can be either \code{additive} or \code{average}.
#' The length must be equal to the number of required emigration traits (see above) and the sequence must match the sequence of the emigration traits listed above.
#' @param InitialDistribution Distribution of the initial values. Can be \code{uniform} or \code{normal}. Should be provided as a vector of strings.
#' The length must be equal to the number of required emigration traits (see above) and the sequence must match the sequence of the emigration traits listed above.
#' @param InitialParameters Parameters for the initial distribution: You must provide two colums min and max  for \code{uniform} distribution and mean and sd for \code{normal} distribution.
#' Each row in the matrix corresponds to an emigration trait. The number of rows must be equal to the number of required emigration traits (see above) and the sequence must match the sequence of the emigration traits listed above.
#' @param IsInherited Should the emigration trait be inherited? Can be either \code{TRUE} or \code{FALSE}.
#' The length must be equal to the number of required emigration traits (see above) and the sequence must match the sequence of the emigration traits listed above.
#' @param MutationDistribution Distribution for mutations to draw from. Can be \code{uniform} or \code{normal}.
#' The length must be equal to the number of required emigration traits (see above) and the sequence must match the sequence of the emigration traits listed above.
#' @param MutationParameters Parameters for the mutation distribution: You must provide two colums: min and max for \code{uniform} distribution and mean and sd for \code{normal} distribution.
#' Each row in the matrix corresponds to an emigration trait. The number of rows must be equal to the number of required emigration traits (see above) and the sequence must match the sequence of the emigration traits listed above.
#' @param MutationRate Mutation rate applicable to this type of loci. Must be between \code{0.0} and \code{1.0}.
#' The length must be equal to the number of required emigration traits (see above) and the sequence must match the sequence of the emigration traits listed above.
#' @param OutputValues If OutputGeneValues in GeneticsFile is enabled, should allele values for this gene be written to output? Ignored if OutputGeneValues is set to \code{FALSE}.
#' The length must be equal to the number of required emigration traits (see above) and the sequence must match the sequence of the emigration traits listed above.
#'

#'
#' @return a parameter object of class "EmigrationTraitsParams"
#' @author Jette Reeg
#' @name EmigrationTraits
#' @export EmigrationTraits
EmigrationTraits<- setClass("EmigrationTraitsParams", slots = c(Positions = "list", #
                                                                NbOfPositions = "ANY", # random or list of integer values
                                                                ExpressionType = "character", # additive or average
                                                                InitialDistribution = "character", # uniform or normal
                                                                InitialParameters = "matrix", # min and max value or mean and sd
                                                                IsInherited = "logical", # T/F
                                                                MutationDistribution = "character", # uniform or normal
                                                                MutationParameters = "matrix", # min mx or mean sd
                                                                MutationRate = "numeric", # float
                                                                OutputValues = "logical")
                                                                , prototype = list(
                                                                    # ExprSex = FALSE, # is TRUE as soon as Emigration is sexdependent
                                                                    # TraitType = NULL, # dispersal: "E_0", "E_alpha", "E_beta"; cannot be 0 is determined by emigration settings
                                                                    Positions = list("random"), # "random" or list of integer values
                                                                    NbOfPositions = 2L, # numeric, only of positions random
                                                                    ExpressionType = "average", # dispersal: "additive" or "average"
                                                                    InitialDistribution = "uniform", # uniform , normal (dispersal)
                                                                    InitialParameters = matrix(c(0.5,0.1), nrow=1), # dispersal: two values: either min/max oder mean+sd
                                                                    IsInherited = FALSE, # only for dispersal
                                                                    MutationDistribution = "uniform", # dispersal: uniform or normal
                                                                    MutationParameters = matrix(c(0.5,0.1), nrow=1), # single value or 2 values
                                                                    MutationRate = c(0.001), # numeric
                                                                    OutputValues = FALSE
                                                                ))
setValidity("EmigrationTraitsParams", function(object) {
    msg <- NULL
    # Check Position and NbOfPositions
    # Positions must be of type list?
    if(class(object@Positions) != "list") {
        msg <- c(msg, "In EmigrationTraits(): Positions must be provided as a list.")
    }
    # NbOfPositions must be either numeric, integer or NULL
    if (!is.null(object@NbOfPositions) && class(object@NbOfPositions) != "numeric" && class(object@NbOfPositions) != "integer") {
        msg <- c(msg, "In EmigrationTraits(): NbrOfPositions must be either NULL (if all positions are given) or numeric (if at least one emigration trait has random positions).")
    }

    if (all(object@Positions == "random")){ # if all positions are random
        if(length(object@NbOfPositions[!is.na(object@NbOfPositions)]) != length(object@Positions) ) {
            msg <- c(msg, "In EmigrationTraits(): For each emigration trait with random positions you must provide the number of positions.")
        }
    } else{ # if NOT all positions are random
        isNumeric <- sapply(object@Positions, is.numeric)
        if (!all(isNumeric)) { # if not all positions are numeric,
            if(object@Positions[isNumeric==FALSE] != "random"){ # then those not numeric must be random
                msg <- c(msg, "In EmigrationTraits(): Positions in emigration traits must be either a vector of integers or random.")
            }
            if (any(is.na(object@NbOfPositions[object@Positions == "random"])) || any(object@NbOfPositions[object@Positions == "random"] <= 0)){ # if number of positions are NA or smaller than 0
                msg <- c(msg, "In EmigrationTraits(): NbrOfPositions must be set to a strictly positive integer for random positions.")
            }
            if (any(!is.na(object@NbOfPositions[object@Positions != "random"]))) { # if there are NbOfPositions supplied for non-random positions
                msg <- c(msg, "In EmigrationTraits(): if Positions is not random NbrOfPositions must be not be set (NA).")
            }
        }
        else { # if all positions are not random
            if (!is.null(object@NbOfPositions)) {
                msg <- c(msg, "In EmigrationTraits(): If positions are not random, you must not specify the number of positions (NbOfPositions).")
            }
        }
    }

    # Check ExpressionType must be additive or average
    if(length(object@ExpressionType) != length(object@Positions)){
        msg <- c(msg, "In EmigrationTraits(): You must provide the ExpressionType for each emigration trait.")
    }
    if (!all(object@ExpressionType %in% c("additive", "average"))) {
        msg <- c(msg, "In EmigrationTraits(): ExpressionType must be either additive or average.")
    }

    # Check InitialDistribution and InitialParameter: Distribution must be uniform or normal and the length of the list must be the same as the number of positions
    if (!all(object@InitialDistribution %in% c("uniform", "normal"))) {
        msg <- c(msg, "In EmigrationTraits(): InitialDistribution must be either normal, or uniform.")
    }

    if(length(object@InitialDistribution) != length(object@Positions)) {
        msg <- c(msg, "In EmigrationTraits(): For each emigration parameter you must provide the InitialDistribution.")
    } else if (nrow(object@InitialParameters) != length(object@Positions)) {
        msg <- c(msg, "In EmigrationTraits(): For each emigration parameter you must provide the InitialParameters. Use one row for each emigration parameter. Check the R help file ?EmigrationTraits for the structure of the matrix.")
    } else {
        # two columns are necessary for mean and sd or min and max
        if (ncol(object@InitialParameters) !=2 || # if DominanceParameters has not 2 columns OR
            any(!is.numeric(object@InitialParameters)) || # if entries are not numeric
            any(is.na(object@InitialParameters))) { # if entries are NA
            msg <- c(msg,"In EmigrationTraits(): For the initial distributions, InitialParams must provide two values for mean (normal) or min (uniform) (first column) and sd (normal) or max (uniform) (second column)")
        }
    }

    # Check IsInherited: must be TRUE or FALSE
    if (length(object@IsInherited) != length(object@Positions)){
        msg <- c(msg, "In EmigrationTraits(): You must provide IsInherited for each emigration trait.")
    }
    if (!all(is.logical(object@IsInherited))) {
        msg <- c(msg, "In EmigrationTraits(): IsInherited must be either TRUE or FALSE." )
    }

    # Check mutation rate
    if(!is.numeric(object@MutationRate) || length(object@MutationRate) != length(object@Positions)){
        msg <- c(msg, "In EmigrationTraits(): You must provide the mutation rate for each emigration trait as a numeric vector.")
    } else {
        if (!is.numeric(object@MutationRate) ||  (any(object@MutationRate < 0.0) || any(object@MutationRate > 1.0))) {
            msg <- c(msg, "In EmigrationTraits(): MutationRate must be between 0.0 and 1.0.")
        }
    }

    # Check MutationDistribution and MutationParameters
    if (!is.null(object@MutationDistribution)){
        if (length(object@MutationDistribution) != length(object@Positions)){
            msg <- c(msg, "In EmigrationTraits(): For each emigration trait you must provide the MutationDistribution.")
        } else if (nrow(object@MutationParameters) != length(object@Positions)) {
            msg <- c(msg, "In EmigrationTraits(): For each emigration trait you must provide the MutationParameters.")
        } else if(!all(object@MutationDistribution %in% c("uniform", "normal"))){
                msg <- c(msg, "In EmigrationTraits(): MutationDistribution must be either normal or uniform for emigration traits.")
        }  else {
            if (ncol(object@MutationParameters) !=2){
                msg <- c(msg,"In EmigrationTraits(): MutationParams must provide two values for uniform and normal distribution: min/mean (first column) and max/sd (second column)")
            }

            if (any(!is.numeric(object@MutationParameters)) ||
                any(is.na(object@MutationParameters))) {
                msg <- c(msg,"In EmigrationTraits(): For a uniform or normal mutation distribution, MutationParams must provide two values for min (first column) and max (second column)")
            }
        }
    }


    # Check OutputValues
    if (length(object@OutputValues) != length(object@Positions) && !all(object@OutputValues  %in% c(TRUE, FALSE))) {
        msg <- c(msg, "In EmigrationTraits(): OutputValues must be provided for all emigration traits and must be either TRUE or FALSE.")
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
    cat("  Emigration traits: \n")
    if (length(object@Positions) == 1) {
        cat("     Trait type: (1) D0 \n")
    }
    if (length(object@Positions) == 2) {
        cat("     Trait types: (1) D0 female, (2) D0 male  \n")
    }
    if (length(object@Positions) == 3) {
        cat("     Trait types: (1) D0, (2) Alpha, (3) Beta  \n")
    }
    if (length(object@Positions) == 6) {
        cat("     Trait types: (1) D0 female, (2) D0 male, (3) Alpha female, (4) Alpha male, (5) Beta female, (6) Beta male \n")
    }
    for (i in 1:length(object@Positions)){
        cat("     Configuration of emigration trait ", i, ": \n")
        if(is.numeric(object@Positions[i])) cat("     Loci positions coding for trait: ", object@Positions[i], "\n")
        if(!is.numeric(object@Positions[i]) && object@Positions[i]=="random") cat("    Loci positions coding for trait randomly chosen with ", object@NbOfPositions[i], " positions\n")
        cat("       Expression type: ", object@ExpressionType[i], "\n")
        cat("       Initial distribution: ", object@InitialDistribution[i], "\n")
        cat("       Initial parameter: ", object@InitialParameters[i,], "\n")
        cat("       IsInherited: ", object@IsInherited[i], "\n")
        cat("       Mutation distribution: ", object@MutationDistribution[i], "\n")
        cat("       Mutation parameters: ", object@MutationParameters[i,], "\n")
        cat("       Mutation rate: ", object@MutationRate, "\n")
        if(object@OutputValues[i]) cat("       Allel values for gene is written to output. \n")
    }
})

### SUBCLASS SETTLEMENTTRAITS

#' Set genetic traits structure for settlement traits
#'
#' @description
#' Only if settlement is density-dependent (\code{DensDep = TRUE}) in \code{\link[RangeShiftR]{Settlement}}, settlement traits can be evolvable.
#'
#' Three settlement traits can evolve:
#'
#'  - The probability of settlement \ifelse{html}{\out{S<sub>0</sub>}}{\eqn{S_0}} \cr
#'  - The slope of the density-dependent settlement function \ifelse{html}{\out{&alpha;}}{\eqn{α}} \cr
#'  - The density threshold of the density-dependent settlement function \ifelse{html}{\out{&beta;}}{\eqn{β}} \cr
#'
#' If settlement is sex dependent, \code{SexDep=TRUE}, you must provide details for both sexes.
#'
#' This results in following number of settlement traits that need to be specifed:
#'
#' - 3 entries/rows if settlement probability is not sex dependent \code{SexDep=FALSE}:  \ifelse{html}{\out{S<sub>0</sub>}}{\eqn{S_0}}, \eqn{α} and \eqn{β} \cr
#' - 6 entries/rows if settlement probability is sex dependent (\code{SexDep=TRUE}): \ifelse{html}{\out{S<sub>0</sub>(f), S<sub>0</sub>(m)}}{\eqn{S_0(f),S_0(m)}}, \ifelse{html}{\out{&alpha;(f), &alpha;(m)}}{\eqn{α(f), α(m)}} and \ifelse{html}{\out{&beta;(f), &beta;(m)}}{\eqn{β(f), β(m)}} \cr
#'
#' The entries of the trait parameters must be provided in the same order as the kernel traits are listed above. If parameters expect a matrix, the rows must match the order of kernel traits listed above.
#'
#' @details
#'
#' Traits set to evolve cannot simultaneously be stage-dependent.
#'
#' The alleles of each trait can be expressed according to an additive model (allele values across all loci are summed) or be averaged.
#'
#' Mutations are additive and can be sampled in either a \code{normal} or \code{uniform} distribution.
#'
#' Initial allele values are sampled in a \code{normal} or \code{uniform} distribution.
#'
#' Settlement traits can also be **not** inherited, that is, allele values are resampled from the initial distribution for every new individual.
#'
#' Dominance values are not applicable for settlement traits.
#'
#'
#' @usage SettlementTraits(Positions = list("random","random","random"), NbOfPositions = c(10, 10, 10),
#' ExpressionType = rep("additive",3),
#' InitialDistribution = rep("normal",3), InitialParameters = matrix(c(rep(0.5,3),(rep(0.1,3), nrow=3),
#' IsInherited = rep(FALSE,3),
#' MutationDistribution = rep("normal",3), MutationParameters = matrix(c(rep(0.5,3),(rep(0.2,3), nrow=3),
#' MutationRate = rep(0.0001,3), OutputValues = rep(FALSE,3))
#'
#' @param Positions Loci positions coding for the trait within genome. Should be provided as a list. Entries can either be a string (\code{"random"}) and/or vectors of integers.
#' The length must be equal to the number of required settlement traits (see above) and the sequence must match the sequence of the settlement traits listed above.
#' @param NbOfPositions Only specify when the \code{Positions} of the settlement trait is set to \code{"random"}, else must be blank (\code{NULL}).
#' The length must be equal to the number of required settlement traits (see above) and the sequence must match the sequence of the settlement traits listed above.
#' @param ExpressionType Type of expression for the settlement trait. Can be either \code{additive} or \code{average}.
#' The length must be equal to the number of required settlement traits (see above) and the sequence must match the sequence of the settlement traits listed above.
#' @param InitialDistribution Distribution of the initial values. Can be \code{uniform} or \code{normal}. Should be provided as a vector of strings.
#' The length must be equal to the number of required settlement traits (see above) and the sequence must match the sequence of the settlement traits listed above.
#' @param InitialParameters Parameters for the initial distribution: You must provide two colums min and max  for \code{uniform} distribution and mean and sd for \code{normal} distribution.
#' Each row in the matrix corresponds to an settlement trait. The number of rows must be equal to the number of required settlement traits (see above) and the sequence must match the sequence of the settlement traits listed above.
#' @param IsInherited Should the settlement trait be inherited? Can be either TRUE or FALSE.
#' The length must be equal to the number of required settlement traits (see above) and the sequence must match the sequence of the settlement traits listed above.
#' @param MutationDistribution Distribution for mutations to draw from. Can be \code{uniform} or \code{normal}.
#' The length must be equal to the number of required settlement traits (see above) and the sequence must match the sequence of the settlement traits listed above.
#' @param MutationParameters Parameters for the mutation distribution: You must provide two colums: min and max for \code{uniform} distribution and mean and sd for \code{normal} distribution.
#' Each row in the matrix corresponds to an settlement trait. The number of rows must be equal to the number of required settlement traits (see above) and the sequence must match the sequence of the settlement traits listed above.
#' @param MutationRate Mutation rate applicable to this type of loci. Must be between 0.0 and 1.0.
#' The length must be equal to the number of required settlement traits (see above) and the sequence must match the sequence of the settlement traits listed above.
#' @param OutputValues If OutputGeneValues in GeneticsFile is enabled, should allele values for this gene be written to output? Ignored if OutputGeneValues is set to \code{FALSE}.
#' The length must be equal to the number of required settlement traits (see above) and the sequence must match the sequence of the settlement traits listed above.
#'
#'
#' @return a parameter object of class "SettlementTraitsParams"
#' @author Jette Reeg
#' @name SettlementTraits
#' @export SettlementTraits
SettlementTraits<- setClass("SettlementTraitsParams", slots = c(Positions = "list", #
                                                                NbOfPositions = "ANY", # random or list of integer values
                                                                ExpressionType = "character", # additive or average
                                                                InitialDistribution = "character", # uniform or normal
                                                                InitialParameters = "matrix", # min and max value or mean and sd
                                                                IsInherited = "logical", # T/F
                                                                MutationDistribution = "character", # uniform or normal
                                                                MutationParameters = "matrix", # min mx or mean sd
                                                                MutationRate = "numeric", # float
                                                                OutputValues = "logical")
                            , prototype = list(
                                Positions = list("random", "random", "random"), # "random" or list of integer values
                                NbOfPositions = c(2, 2, 2), # numeric, only of positions random
                                ExpressionType = rep("additive",3), # dispersal: "additive" or "average"
                                InitialDistribution = rep("uniform",3), # uniform , normal (dispersal)
                                InitialParameters = matrix(c(0.5,0.5,0.5,0.1,0.1,0.1), nrow=3), # dispersal: two values: either min/max oder mean+sd
                                IsInherited = rep(FALSE, 3), # only for dispersal
                                MutationDistribution = rep("uniform",3), # dispersal: uniform or normal
                                MutationParameters = matrix(c(0.5,0.5,0.5,0.1,0.1,0.1), nrow=3), # single value or 2 values
                                MutationRate = rep(0.001,3), # numeric
                                OutputValues = rep(FALSE,3)
                            ))
setValidity("SettlementTraitsParams", function(object) {
    msg <- NULL
    # Check Position and NbOfPositions
    # Positions must be of type list?
    if(class(object@Positions) != "list") {
        msg <- c(msg, "In SettlementTraits(): Positions must be provided as a list.")
    }
    # NbOfPositions must be either numeric, integer or NULL
    if (!is.null(object@NbOfPositions) && class(object@NbOfPositions) != "numeric" && class(object@NbOfPositions) != "integer") {
        msg <- c(msg, "In SettlementTraits(): NbrOfPositions must be either NULL (if all positions are given) or numeric (if at least one settlement trait has random positions).")
    }

    if (all(object@Positions == "random")){ # if all positions are random
        if(length(object@NbOfPositions[!is.na(object@NbOfPositions)]) != length(object@Positions) ) {
            msg <- c(msg, "In SettlementTraits(): For each settlement trait with random positions you must provide the number of positions.")
        }
    } else{ # if NOT all positions are random
        isNumeric <- sapply(object@Positions, is.numeric)
        if (!all(isNumeric)) { # if not all positions are numeric,
            if(object@Positions[isNumeric==FALSE] != "random"){ # then those not numeric must be random
                msg <- c(msg, "In SettlementTraits(): Positions in settlement traits must be either a vector of integers or random.")
            }
            if (any(is.na(object@NbOfPositions[object@Positions == "random"])) || any(object@NbOfPositions[object@Positions == "random"] <= 0)){ # if number of positions are NA or smaller than 0
                msg <- c(msg, "In SettlementTraits(): NbrOfPositions must be set to a strictly positive integer for random positions.")
            }
            if (any(!is.na(object@NbOfPositions[object@Positions != "random"]))) { # if there are NbOfPositions supplied for non-random positions
                msg <- c(msg, "In SettlementTraits(): if Positions is not random NbrOfPositions must be not be set (NA).")
            }
        }
        else { # if all positions are not random
            if (!is.null(object@NbOfPositions)) {
                msg <- c(msg, "In SettlementTraits(): If positions are not random, you must not specify the number of positions (NbOfPositions).")
            }
        }
    }

    # Check ExpressionType must be additive or average
    if(length(object@ExpressionType) != length(object@Positions)){
        msg <- c(msg, "In SettlementTraits(): You must provide the ExpressionType for each settlement trait.")
    }
    if (!all(object@ExpressionType %in% c("additive", "average"))) {
        msg <- c(msg, "In SettlementTraits(): ExpressionType must be either additive or average.")
    }

    # Check InitialDistribution and InitialParameter: Distribution must be uniform or normal and the length of the list must be the same as the number of positions
    if (!all(object@InitialDistribution %in% c("uniform", "normal"))) {
        msg <- c(msg, "In SettlementTraits(): InitialDistribution must be either normal, or uniform.")
    }

    if(length(object@InitialDistribution) != length(object@Positions)) {
        msg <- c(msg, "In SettlementTraits(): For each settlement parameter you must provide the InitialDistribution.")
    } else if (nrow(object@InitialParameters) != length(object@Positions)) {
        msg <- c(msg, "In SettlementTraits(): For each settlement parameter you must provide the InitialParameters. Use one row for each settlement parameter. Check the R help file ?SettlementTraits for the structure of the matrix.")
    } else {
        # two columns are necessary for mean and sd or min and max
        if (ncol(object@InitialParameters) !=2 || # if DominanceParameters has not 2 columns OR
            any(!is.numeric(object@InitialParameters)) || # if entries are not numeric
            any(is.na(object@InitialParameters))) { # if entries are NA
            msg <- c(msg,"In SettlementTraits(): For the initial distributions, InitialParams must provide two values for mean (normal) or min (uniform) (first column) and sd (normal) or max (uniform) (second column)")
        }
    }

    # Check IsInherited: must be TRUE or FALSE
    if (length(object@IsInherited) != length(object@Positions)){
        msg <- c(msg, "In SettlementTraits(): You must provide IsInherited for each settlement trait.")
    }
    if (!all(is.logical(object@IsInherited))) {
        msg <- c(msg, "In SettlementTraits(): IsInherited must be either TRUE or FALSE." )
    }

    # Check mutation rate
    if(!is.numeric(object@MutationRate) || length(object@MutationRate) != length(object@Positions)){
        msg <- c(msg, "In SettlementTraits(): You must provide the mutation rate for each settlement trait as a numeric vector.")
    } else {
        if (!is.numeric(object@MutationRate) ||  (any(object@MutationRate < 0.0) || any(object@MutationRate > 1.0))) {
            msg <- c(msg, "In SettlementTraits(): MutationRate must be between 0.0 and 1.0.")
        }
    }

    # Check MutationDistribution and MutationParameters
    if (!is.null(object@MutationDistribution)){
        if (length(object@MutationDistribution) != length(object@Positions)){
            msg <- c(msg, "In SettlementTraits(): For each settlement trait you must provide the MutationDistribution.")
        } else if (nrow(object@MutationParameters) != length(object@Positions)) {
            msg <- c(msg, "In SettlementTraits(): For each settlement trait you must provide the MutationParameters.")
        } else if(!all(object@MutationDistribution %in% c("uniform", "normal"))){
            msg <- c(msg, "In SettlementTraits(): MutationDistribution must be either normal or uniform for settlement traits.")
        }  else {
            if (ncol(object@MutationParameters) !=2){
                msg <- c(msg,"In SettlementTraits(): MutationParams must provide two values for uniform and normal distribution: min/mean (first column) and max/sd (second column)")
            }

            if (any(!is.numeric(object@MutationParameters)) ||
                any(is.na(object@MutationParameters))) {
                msg <- c(msg,"In SettlementTraits(): For a uniform or normal mutation distribution, MutationParams must provide two values for min (first column) and max (second column)")
            }
        }
    }


    # Check OutputValues
    if (length(object@OutputValues) != length(object@Positions) && !all(object@OutputValues  %in% c(TRUE, FALSE))) {
        msg <- c(msg, "In SettlementTraits(): OutputValues must be provided for all settlement traits and must be either TRUE or FALSE.")
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
    cat("  Settlement traits: \n")
    if (length(object@Positions) == 3) {
        cat("     Trait types: (1) S0, (2) Alpha, (3) Beta  \n")
    }
    if (length(object@Positions) == 6) {
        cat("     Trait types: (1) S0 female, (2) S0 male, (3) Alpha female, (4) Alpha male, (5) Beta female, (6) Beta male \n")
    }
    for (i in 1:length(object@Positions)){
        cat("     Configuration of settlement trait ", i, ": \n")
        if(is.numeric(object@Positions[i])) cat("     Loci positions coding for trait: ", object@Positions[i], "\n")
        if(!is.numeric(object@Positions[i]) && object@Positions[i]=="random") cat("    Loci positions coding for trait randomly chosen with ", object@NbOfPositions[i], " positions\n")
        cat("       Expression type: ", object@ExpressionType[i], "\n")
        cat("       Initial distribution: ", object@InitialDistribution[i], "\n")
        cat("       Initial parameter: ", object@InitialParameters[i,], "\n")
        cat("       IsInherited: ", object@IsInherited[i], "\n")
        cat("       Mutation distribution: ", object@MutationDistribution[i], "\n")
        cat("       Mutation parameters: ", object@MutationParameters[i,], "\n")
        cat("       Mutation rate: ", object@MutationRate, "\n")
        if(object@OutputValues[i]) cat("       Allel values for gene is written to output. \n")
    }
})

### SUBCLASS KERNELTRAITS

#' Set genetic traits structure for kernel traits
#'
#' @description
#'
#' Depending on the settings of \code{\link[RangeShiftR]{DispersalKernel}}, the following traits can be evolvable:
#'
#' - The mean parameters(s) of the (two) dispersel kernel(s) \cr
#' - The probability of using the first kernel, if a double kernel is used (\code{DoubleKernel=TRUE}). \cr
#'
#' If dispersal kernels are sex dependent, \code{SexDep=TRUE}, you must provide details for both sexes.
#'
#' This results in the following number of dispersal kernel traits that need to be specified:
#'
#' - 1 if \code{DoubleKernel=FALSE} and \code{SexDep=FALSE}: \ifelse{html}{\out{&delta;}}{\eqn{δ}} \cr
#' - 2 if \code{DoubleKernel=FALSE} and \code{SexDep=TRUE}: \ifelse{html}{\out{&delta;(f)}}{\eqn{δ}(f)}, \ifelse{html}{\out{&delta;(m)}}{\eqn{δ}(m)}\cr
#' - 3 if \code{DoubleKernel=TRUE} and \code{SexDep=FALSE}: \ifelse{html}{\out{&delta;<sub>1</sub>, &delta;<sub>2</sub>, p<sub>I</sub>}}{\eqn{δ_1, δ_2, p_I}}\cr
#' - 6 if \code{DoubleKernel=TRUE} and \code{SexDep=TRUE}: \ifelse{html}{\out{&delta;<sub>1</sub>(f), &delta;<sub>1</sub>(m), &delta;<sub>2</sub>(f), &delta;<sub>2</sub>(m), p<sub>I</sub>(f), p<sub>I</sub>(m)}}{\eqn{δ_1(f), δ_1(m), δ_2(f), δ_2(m), p_I(f), p_I(m)}} \cr
#'
#' The entries of the trait parameters must be provided in the same order as the kernel traits are listed above. If parameters expect a matrix, the rows must match the order of kernel traits listed above.
#'
#' @details
#'
#' Traits set to evolve cannot simultaneously be stage-dependent.
#'
#' The alleles of each trait can be expressed according to an additive model (allele values across all loci are summed) or be averaged.
#'
#' Mutations are additive and can be sampled in either a \code{normal} or \code{uniform} distribution.
#'
#' Initial allele values are sampled in a \code{normal} or \code{uniform} distribution.
#'
#' Dispersal kernel traits can also be **not** inherited, that is, allele values are resampled from the initial distribution for every new individual.
#'
#' Dominance values are not applicable for dispersal kernel traits.
#'
#' @usage KernelTraits(Positions = list("random"), NbOfPositions = c(10),
#' ExpressionType = rep("additive",1),
#' InitialDistribution = rep("normal",1), InitialParameters = matrix(c(rep(0.5,1),(rep(0.1,1), nrow=1),
#' IsInherited = rep(FALSE,1),
#' MutationDistribution = rep("normal",1), MutationParameters = matrix(c(rep(0.5,1),(rep(0.2,1), nrow=1),
#' MutationRate = rep(0.0001,1), OutputValues = rep(FALSE,1))
#'
#' @param Positions Loci positions coding for the trait within genome. Should be provided as a list. Entries can either be a string (\code{"random"}) and/or vectors of integers.
#' The length must be equal to the number of required kernel traits (see above) and the sequence must match the sequence of the kernel traits listed above.
#' @param NbOfPositions Only specify when the \code{Positions} of the kernel trait is set to ‘random’, else must be blank (\code{NULL}).
#' The length must be equal to the number of required kernel traits (see above) and the sequence must match the sequence of the kernel traits listed above.
#' @param ExpressionType Type of expression for the emigration trait. Can be either \code{additive} or \code{average}.
#' The length must be equal to the number of required kernel traits (see above) and the sequence must match the sequence of the kernel traits listed above.
#' @param InitialDistribution Distribution of the initial values. Can be \code{uniform} or \code{normal}. Should be provided as a vector of strings.
#' The length must be equal to the number of required kernel traits (see above) and the sequence must match the sequence of the kernel traits listed above.
#' @param InitialParameters Parameters for the initial distribution: You must provide two colums min and max  for \code{uniform} distribution and mean and sd for \code{normal} distribution.
#' Each row in the matrix corresponds to an emigration trait. The number of rows must be equal to the number of required emigration traits (see above) and the sequence must match the sequence of the emigration traits listed above.
#' @param IsInherited Should the emigration trait be inherited? Can be either \code{TRUE} or \code{FALSE}.
#' The length must be equal to the number of required kernel traits (see above) and the sequence must match the sequence of the kernel traits listed above.
#' @param MutationDistribution Distribution for mutations to draw from. Can be \code{uniform} or \code{normal}.
#' The length must be equal to the number of required kernel traits (see above) and the sequence must match the sequence of the kernel traits listed above.
#' @param MutationParameters Parameters for the mutation distribution: You must provide two colums: min and max for \code{uniform} distribution and mean and sd for \code{normal} distribution.
#' Each row in the matrix corresponds to an emigration trait. The number of rows must be equal to the number of required emigration traits (see above) and the sequence must match the sequence of the emigration traits listed above.
#' @param MutationRate Mutation rate applicable to this type of loci. Must be between 0.0 and 1.0.
#' The length must be equal to the number of required kernel traits (see above) and the sequence must match the sequence of the kernel traits listed above.
#' @param OutputValues If OutputGeneValues in GeneticsFile is enabled, should allele values for this gene be written to output? Ignored if OutputGeneValues is set to \code{FALSE}.
#' The length must be equal to the number of required kernel traits (see above) and the sequence must match the sequence of the kernel traits listed above.
#'
#' @return a parameter object of class "KernelTraitsParams"
#' @author Jette Reeg
#' @name KernelTraits
#' @export KernelTraits
KernelTraits<- setClass("KernelTraitsParams", slots = c(Positions = "list", #
                                                        NbOfPositions = "ANY", # random or list of integer values
                                                        ExpressionType = "character", # additive or average
                                                        InitialDistribution = "character", # uniform or normal
                                                        InitialParameters = "matrix", # min and max value or mean and sd
                                                        IsInherited = "logical", # T/F
                                                        MutationDistribution = "character", # uniform or normal
                                                        MutationParameters = "matrix", # min mx or mean sd
                                                        MutationRate = "numeric", # float
                                                        OutputValues = "logical")
                        , prototype = list(
                            Positions = list("random"), # "random" or list of integer values
                            NbOfPositions = 2L, # numeric, only of positions random
                            ExpressionType = "additive", # dispersal: "additive" or "average"
                            InitialDistribution = "uniform", # uniform , normal (dispersal)
                            InitialParameters = matrix(c(0.5,0.1), nrow=1), # dispersal: two values: either min/max oder mean+sd
                            IsInherited = FALSE, # only for dispersal
                            MutationDistribution = "uniform", # dispersal: uniform or normal
                            MutationParameters = matrix(c(0.5,0.1), nrow=1), # single value or 2 values
                            MutationRate = c(0.001), # numeric
                            OutputValues = FALSE
                        ))
setValidity("KernelTraitsParams", function(object) {
    msg <- NULL
    # Check Position and NbOfPositions
    # Positions must be of type list?
    if(class(object@Positions) != "list") {
        msg <- c(msg, "In KernelTraits(): Positions must be provided as a list.")
    }
    # NbOfPositions must be either numeric, integer or NULL
    if (!is.null(object@NbOfPositions) && class(object@NbOfPositions) != "numeric" && class(object@NbOfPositions) != "integer") {
        msg <- c(msg, "In KernelTraits(): NbrOfPositions must be either NULL (if all positions are given) or numeric (if at least one kernel trait has random positions).")
    }
    if(!(length(object@Positions) %in% c(1,2,3,6))){
        msg <- c(msg, "In KernelTraits(): You must provide 1 (DoubleKernel=FALSE, SexDep=FALSE), 2 (SexDep=TRUE), 3 (DoubleKernel=TRUE, SexDep=FALSE) or (DoubleKernel=TRUE, SexDep=TRUE) positions for the kernel traits.")
    }
    if (all(object@Positions == "random")){ # if all positions are random
        if(length(object@NbOfPositions[!is.na(object@NbOfPositions)]) != length(object@Positions) ) {
            msg <- c(msg, "In KernelTraits(): For each kernel trait with random positions you must provide the number of positions.")
        }
    } else{ # if NOT all positions are random
        isNumeric <- sapply(object@Positions, is.numeric)
        if (!all(isNumeric)) { # if not all positions are numeric,
            if(object@Positions[isNumeric==FALSE] != "random"){ # then those not numeric must be random
                msg <- c(msg, "In KernelTraits(): Positions in kernel traits must be either a vector of integers or random.")
            }
            if (any(is.na(object@NbOfPositions[object@Positions == "random"])) || any(object@NbOfPositions[object@Positions == "random"] <= 0)){ # if number of positions are NA or smaller than 0
                msg <- c(msg, "In KernelTraits(): NbrOfPositions must be set to a strictly positive integer for random positions.")
            }
            if (any(!is.na(object@NbOfPositions[object@Positions != "random"]))) { # if there are NbOfPositions supplied for non-random positions
                msg <- c(msg, "In KernelTraits(): if Positions is not random NbrOfPositions must be not be set (NA).")
            }
        }
        else { # if all positions are not random
            if (!is.null(object@NbOfPositions)) {
                msg <- c(msg, "In KernelTraits(): If positions are not random, you must not specify the number of positions (NbOfPositions).")
            }
        }
    }

    # Check ExpressionType must be additive or average
    if(!(length(object@ExpressionType) %in% c(1,2,3,6))){
        msg <- c(msg, "In KernelTraits(): You must provide the ExpressionType for each kernel trait.")
    }
    if (!all(object@ExpressionType %in% c("additive", "average"))) {
        msg <- c(msg, "In KernelTraits(): ExpressionType must be either additive or average.")
    }

    # Check InitialDistribution and InitialParameter: Distribution must be uniform or normal and the length of the list must be the same as the number of positions
    if (!(length(object@InitialDistribution) %in% c(1,2,3,6))){
        msg <- c(msg, "In KernelTraits(): You must provide the InitialDistribution for each kernel trait.")
    }
    if (!all(object@InitialDistribution %in% c("uniform", "normal"))) {
        msg <- c(msg, "In KernelTraits(): InitialDistribution must be either normal, or uniform.")
    }

    if (!(nrow(object@InitialParameters) %in% c(1,2,3,6))) {
        msg <- c(msg, "In KernelTraits(): For each kernel parameter you must provide the InitialParameters. Use one row for each kernel parameter.")
    } else {
        # two columns are necessary for mean and sd or min and max
        if (ncol(object@InitialParameters) !=2 || # if DominanceParameters has not 2 columns OR
            any(!is.numeric(object@InitialParameters)) || # if entries are not numeric
            any(is.na(object@InitialParameters))) { # if entries are NA
            msg <- c(msg,"In KernelTraits(): For the initial distributions, InitialParams must provide two values for mean (normal) or min (uniform) (first column) and sd (normal) or max (uniform) (second column)")
        }
    }

    # Check IsInherited: must be TRUE or FALSE
    if (!(length(object@IsInherited) %in% c(1,2,3,6))){
        msg <- c(msg, "In KernelTraits(): You must provide IsInherited for each kernel trait.")
    }
    if (!all(is.logical(object@IsInherited))) {
        msg <- c(msg, "In KernelTraits(): IsInherited must be either TRUE or FALSE." )
    }

    # Check mutation rate
    if(!is.numeric(object@MutationRate) || !(length(object@MutationRate) %in% c(1,4))){
        msg <- c(msg, "In KernelTraits(): You must provide the mutation rate for each kernel trait as a numeric vector.")
    } else {
        if (!is.numeric(object@MutationRate) ||  (any(object@MutationRate < 0.0) || any(object@MutationRate > 1.0))) {
            msg <- c(msg, "In KernelTraits(): MutationRate must be between 0.0 and 1.0.")
        }
    }

    # Check MutationDistribution and MutationParameters
    if (!is.null(object@MutationDistribution)){
        if (!(length(object@MutationDistribution) %in% c(1,2,3,6))){
            msg <- c(msg, "In KernelTraits(): For each kernel trait you must provide the MutationDistribution.")
        } else if (!(nrow(object@MutationParameters) %in% c(1,2,3,6))) {
            msg <- c(msg, "In KernelTraits(): For each kernel trait you must provide the MutationParameters.")
        } else if(!all(object@MutationDistribution %in% c("uniform", "normal"))){
            msg <- c(msg, "In KernelTraits(): MutationDistribution must be either normal or uniform for kernel traits.")
        }  else {
            if (ncol(object@MutationParameters) !=2){
                msg <- c(msg,"In KernelTraits(): MutationParams must provide two values for uniform and normal distribution: min/mean (first column) and max/sd (second column)")
            }

            if (any(!is.numeric(object@MutationParameters)) ||
                any(is.na(object@MutationParameters))) {
                msg <- c(msg,"In KernelTraits(): For a uniform or normal mutation distribution, MutationParams must provide two values for min (first column) and max (second column)")
            }
        }
    }

    # Check OutputValues
    if (!(length(object@OutputValues) %in% c(1,2,3,6)) && !all(object@OutputValues  %in% c(TRUE, FALSE))) {
        msg <- c(msg, "In KernelTraits(): OutputValues must be provided for all kernel traits and must be either TRUE or FALSE.")
    }

    if (is.null(msg)) TRUE else msg
})
setMethod("initialize", "KernelTraitsParams", function(.Object, ...) {
    this_func = "CorrRWTraits(): "
    args <- list(...)
    .Object <- callNextMethod()
    if ( length(args) == 0 ) {
        validObject(.Object)
    }
    .Object
})
setMethod("show", "KernelTraitsParams", function(object){
    cat("  Kernel traits: \n")
    if (length(object@Positions) == 1) {
        cat("     Trait types: (1) Delta  \n")
    }
    if (length(object@Positions) == 2) {
        cat("     Trait types: (1) Delta(f) (2) Delta(m)  \n")
    }
    if (length(object@Positions) == 3) {
        cat("     Trait types: (1) Delta_1 (2) Delta_2 (3) P_I  \n")
    }
    if (length(object@Positions) == 6) {
        cat("     Trait types: (1) Delta_1(f) (2) Delta_1(m) (3) Delta_2(f) (4) Delta_2(m) (5) P_I(f) (6) P_I(m)  \n")
    }

    for (i in 1:length(object@Positions)){
        cat("     Configuration of CorrRW trait ", i, ": \n")
        if(is.numeric(object@Positions[i])) cat("     Loci positions coding for trait: ", object@Positions[i], "\n")
        if(!is.numeric(object@Positions[i]) && object@Positions[i]=="random") cat("    Loci positions coding for trait randomly chosen with ", object@NbOfPositions[i], " positions\n")
        cat("       Expression type: ", object@ExpressionType[i], "\n")
        cat("       Initial distribution: ", object@InitialDistribution[i], "\n")
        cat("       Initial parameter: ", object@InitialParameters[i,], "\n")
        cat("       IsInherited: ", object@IsInherited[i], "\n")
        cat("       Mutation distribution: ", object@MutationDistribution[i], "\n")
        cat("       Mutation parameters: ", object@MutationParameters[i,], "\n")
        cat("       Mutation rate: ", object@MutationRate, "\n")
        if(object@OutputValues[i]) cat("       Allel values for gene is written to output. \n")
    }
})

### SUBCLASS SMSTRAITS


#' Set genetic traits structure for SMS traits
#'
#' @description
#'
#' Depending on the settings of \code{\link[RangeShiftR]{SMS}}, the following traits can be evolvable:
#'
#' - The directional persistence of the individual \code{DP} \cr
#' - The \code{GoalBias}, that is the tendency for an individual to move away from its natal patch \cr
#' - The slope of the distance dependent decay of the goal bias \code{AlphaDB} \cr
#' - The distance threshold of the distance dependent devay of the goal bias \code{BetaDB} \cr
#'
#' If dispersal kernels are sex dependent, \code{SexDep=TRUE}, you must provide details for both sexes.
#'
#' This results in the following number of SMS traits that need to be specified:
#'
#' - 1 if \code{GoalType = 0}: directional persistence \code{DP} \cr
#' - 4 if \code{GoalType = 2}: \code{DP}, \code{GoalBias}, \code{AlphaDB}, \code{BetaDB} \cr
#'
#' The entries of the trait parameters must be provided in the same order as the SMS traits are listed above. If parameters expect a matrix, the rows must match the order of SMS traits listed above.
#'
#' @details
#'
#' Traits set to evolve cannot simultaneously be stage-dependent.
#'
#' The alleles of each trait can be expressed according to an additive model (allele values across all loci are summed) or be averaged.
#'
#' Mutations are additive and can be sampled in either a \code{normal} or \code{uniform} distribution.
#'
#' Initial allele values are sampled in a \code{normal} or \code{uniform} distribution.
#'
#' SMS traits can also be **not** inherited, that is, allele values are resampled from the initial distribution for every new individual.
#'
#' Dominance values are not applicable for SMS traits.
#'
#' @usage SMSTraits(Positions = list("random"), NbOfPositions = c(10),
#' ExpressionType = rep("additive",1),
#' InitialDistribution = rep("normal",1), InitialParameters = matrix(c(rep(0.5,1),(rep(0.1,1), nrow=1),
#' IsInherited = rep(FALSE,1),
#' MutationDistribution = rep("normal",1), MutationParameters = matrix(c(rep(0.5,1),(rep(0.2,1), nrow=1),
#' MutationRate = rep(0.0001,1), OutputValues = rep(FALSE,1))
#'
#' @param Positions Loci positions coding for the trait within genome. Should be provided as a list. Entries can either be a string (\code{"random"}) and/or vectors of integers.
#' The length must be 1 (\code{GoalType=0}) or 4 (\code{GoalType=2}) for the required SMS traits \code{DP} and, if \code{GoalType=2} \code{GoalBias}, , \code{AlphaDB} and , \code{BetaDB}.
#' @param NbOfPositions Only specify when the \code{Positions} of the SMS trait is set to ‘random’, else must be blank (\code{NULL}).
#' The length must be 1 (\code{GoalType=0}) or 4 (\code{GoalType=2}) for the required SMS traits \code{DP} and, if \code{GoalType=2} \code{GoalBias}, , \code{AlphaDB} and , \code{BetaDB}.
#' @param ExpressionType Type of expression for the emigration trait. Can be either \code{additive} or \code{average}.
#' The length must be 1 (\code{GoalType=0}) or 4 (\code{GoalType=2}) for the required SMS traits \code{DP} and, if \code{GoalType=2} \code{GoalBias}, , \code{AlphaDB} and , \code{BetaDB}.
#' @param InitialDistribution Distribution of the initial values. Can be \code{uniform} or \code{normal}. Should be provided as a vector of strings.
#' The length must be 1 (\code{GoalType=0}) or 4 (\code{GoalType=2}) for the required SMS traits \code{DP} and, if \code{GoalType=2} \code{GoalBias}, , \code{AlphaDB} and , \code{BetaDB}.
#' @param InitialParameters Parameters for the initial distribution: You must provide two colums min and max  for \code{uniform} distribution and mean and sd for \code{normal} distribution.
#' Each row in the matrix corresponds to an SMS trait.
#' The number of rows must be 1 (\code{GoalType=0}) or 4 (\code{GoalType=2}) for the required SMS traits \code{DP} and, if \code{GoalType=2} \code{GoalBias}, , \code{AlphaDB} and , \code{BetaDB}.
#' @param IsInherited Should the emigration trait be inherited? Can be either \code{TRUE} or \code{FALSE}.
#' The length must be 1 (\code{GoalType=0}) or 4 (\code{GoalType=2}) for the required SMS traits \code{DP} and, if \code{GoalType=2} \code{GoalBias}, , \code{AlphaDB} and , \code{BetaDB}.
#' @param MutationDistribution Distribution for mutations to draw from. Can be \code{uniform} or \code{normal}.
#' The length must be 1 (\code{GoalType=0}) or 4 (\code{GoalType=2}) for the required SMS traits \code{DP} and, if \code{GoalType=2} \code{GoalBias}, , \code{AlphaDB} and , \code{BetaDB}.
#' @param MutationParameters Parameters for the mutation distribution: You must provide two colums: min and max for \code{uniform} distribution and mean and sd for \code{normal} distribution.
#' Each row in the matrix corresponds to an SMS trait.
#' The number of rows must be 1 (\code{GoalType=0}) or 4 (\code{GoalType=2}) for the required SMS traits \code{DP} and, if \code{GoalType=2} \code{GoalBias}, , \code{AlphaDB} and , \code{BetaDB}.
#' @param MutationRate Mutation rate applicable to this type of loci. Must be between 0.0 and 1.0.
#' The length must be 1 (\code{GoalType=0}) or 4 (\code{GoalType=2}) for the required SMS traits \code{DP} and, if \code{GoalType=2} \code{GoalBias}, , \code{AlphaDB} and , \code{BetaDB}.
#' @param OutputValues If OutputGeneValues in GeneticsFile is enabled, should allele values for this gene be written to output? Ignored if OutputGeneValues is set to \code{FALSE}.
#' The length must be 1 (\code{GoalType=0}) or 4 (\code{GoalType=2}) for the required SMS traits \code{DP} and, if \code{GoalType=2} \code{GoalBias}, , \code{AlphaDB} and , \code{BetaDB}.
#'
#'
#' @return a parameter object of class "SMSTraitsParams"
#' @author Jette Reeg
#' @name SMSTraits
#' @export SMSTraits
SMSTraits<- setClass("SMSTraitsParams", slots = c(Positions = "list", #
                                                        NbOfPositions = "ANY", # random or list of integer values
                                                        ExpressionType = "character", # additive or average
                                                        InitialDistribution = "character", # uniform or normal
                                                        InitialParameters = "matrix", # min and max value or mean and sd
                                                        IsInherited = "logical", # T/F
                                                        MutationDistribution = "character", # uniform or normal
                                                        MutationParameters = "matrix", # min mx or mean sd
                                                        MutationRate = "numeric", # float
                                                        OutputValues = "logical")
                        , prototype = list(
                            Positions = list("random"), # "random" or list of integer values
                            NbOfPositions = 2L, # numeric, only of positions random
                            ExpressionType = "additive", # dispersal: "additive" or "average"
                            InitialDistribution = "uniform", # uniform , normal (dispersal)
                            InitialParameters = matrix(c(0.5,0.1), nrow=1), # dispersal: two values: either min/max oder mean+sd
                            IsInherited = FALSE, # only for dispersal
                            MutationDistribution = "uniform", # dispersal: uniform or normal
                            MutationParameters = matrix(c(0.5,0.1), nrow=1), # single value or 2 values
                            MutationRate = c(0.001), # numeric
                            OutputValues = FALSE
                        ))
setValidity("SMSTraitsParams", function(object) {
    msg <- NULL
    # Check Position and NbOfPositions
    # Positions must be of type list?
    if(class(object@Positions) != "list") {
        msg <- c(msg, "In SMSTraits(): Positions must be provided as a list.")
    }
    # NbOfPositions must be either numeric, integer or NULL
    if (!is.null(object@NbOfPositions) && class(object@NbOfPositions) != "numeric" && class(object@NbOfPositions) != "integer") {
        msg <- c(msg, "In SMSTraits(): NbrOfPositions must be either NULL (if all positions are given) or numeric (if at least one SMS trait has random positions).")
    }
    if(!(length(object@Positions) %in% c(1,4))){
        msg <- c(msg, "In SMSTraits(): You must provide 1 (GoalType=0) or 4 (GoalType=2) positions for the SMS traits DP and, if GoalType=2 GoalBias, AlphaDB and BetaDB.")
    }
    if (all(object@Positions == "random")){ # if all positions are random
        if(length(object@NbOfPositions[!is.na(object@NbOfPositions)]) != length(object@Positions) ) {
            msg <- c(msg, "In SMSTraits(): For each SMS trait with random positions you must provide the number of positions.")
        }
    } else{ # if NOT all positions are random
        isNumeric <- sapply(object@Positions, is.numeric)
        if (!all(isNumeric)) { # if not all positions are numeric,
            if(object@Positions[isNumeric==FALSE] != "random"){ # then those not numeric must be random
                msg <- c(msg, "In SMSTraits(): Positions in SMS traits must be either a vector of integers or random.")
            }
            if (any(is.na(object@NbOfPositions[object@Positions == "random"])) || any(object@NbOfPositions[object@Positions == "random"] <= 0)){ # if number of positions are NA or smaller than 0
                msg <- c(msg, "In SMSTraits(): NbrOfPositions must be set to a strictly positive integer for random positions.")
            }
            if (any(!is.na(object@NbOfPositions[object@Positions != "random"]))) { # if there are NbOfPositions supplied for non-random positions
                msg <- c(msg, "In SMSTraits(): if Positions is not random NbrOfPositions must be not be set (NA).")
            }
        }
        else { # if all positions are not random
            if (!is.null(object@NbOfPositions)) {
                msg <- c(msg, "In SMSTraits(): If positions are not random, you must not specify the number of positions (NbOfPositions).")
            }
        }
    }

    # Check ExpressionType must be additive or average
    if(!(length(object@ExpressionType) %in% c(1,4))){
        msg <- c(msg, "In SMSTraits(): You must provide the ExpressionType for each SMS trait.")
    }
    if (!all(object@ExpressionType %in% c("additive", "average"))) {
        msg <- c(msg, "In SMSTraits(): ExpressionType must be either additive or average.")
    }

    # Check InitialDistribution and InitialParameter: Distribution must be uniform or normal and the length of the list must be the same as the number of positions
    if (!(length(object@InitialDistribution) %in% c(1,4))){
        msg <- c(msg, "In SMSTraits(): You must provide the InitialDistribution for each SMS trait.")
    }
    if (!all(object@InitialDistribution %in% c("uniform", "normal"))) {
        msg <- c(msg, "In SMSTraits(): InitialDistribution must be either normal, or uniform.")
    }

    if (!(nrow(object@InitialParameters) %in% c(1,4))) {
        msg <- c(msg, "In SMSTraits(): For each SMS parameter you must provide the InitialParameters. Use one row for each SMS parameter.")
    } else {
        # two columns are necessary for mean and sd or min and max
        if (ncol(object@InitialParameters) !=2 || # if DominanceParameters has not 2 columns OR
            any(!is.numeric(object@InitialParameters)) || # if entries are not numeric
            any(is.na(object@InitialParameters))) { # if entries are NA
            msg <- c(msg,"In SMSTraits(): For the initial distributions, InitialParams must provide two values for mean (normal) or min (uniform) (first column) and sd (normal) or max (uniform) (second column)")
        }
    }

    # Check IsInherited: must be TRUE or FALSE
    if (!(length(object@IsInherited) %in% c(1,4))){
        msg <- c(msg, "In SMSTraits(): You must provide IsInherited for each SMS trait.")
    }
    if (!all(is.logical(object@IsInherited))) {
        msg <- c(msg, "In SMSTraits(): IsInherited must be either TRUE or FALSE." )
    }

    # Check mutation rate
    if(!is.numeric(object@MutationRate) || !(length(object@MutationRate) %in% c(1,4))){
        msg <- c(msg, "In SMSTraits(): You must provide the mutation rate for each SMS trait as a numeric vector.")
    } else {
        if (!is.numeric(object@MutationRate) ||  (any(object@MutationRate < 0.0) || any(object@MutationRate > 1.0))) {
            msg <- c(msg, "In SMSTraits(): MutationRate must be between 0.0 and 1.0.")
        }
    }

    # Check MutationDistribution and MutationParameters
    if (!is.null(object@MutationDistribution)){
        if (!(length(object@MutationDistribution) %in% c(1,4))){
            msg <- c(msg, "In SMSTraits(): For each SMS trait you must provide the MutationDistribution.")
        } else if (!(nrow(object@MutationParameters) %in% c(1,4))) {
            msg <- c(msg, "In SMSTraits(): For each SMS trait you must provide the MutationParameters.")
        } else if(!all(object@MutationDistribution %in% c("uniform", "normal"))){
            msg <- c(msg, "In SMSTraits(): MutationDistribution must be either normal or uniform for SMS traits.")
        }  else {
            if (ncol(object@MutationParameters) !=2){
                msg <- c(msg,"In SMSTraits(): MutationParams must provide two values for uniform and normal distribution: min/mean (first column) and max/sd (second column)")
            }

            if (any(!is.numeric(object@MutationParameters)) ||
                any(is.na(object@MutationParameters))) {
                msg <- c(msg,"In SMSTraits(): For a uniform or normal mutation distribution, MutationParams must provide two values for min (first column) and max (second column)")
            }
        }
    }


    # Check OutputValues
    if (!(length(object@OutputValues) %in% c(1,4)) && !all(object@OutputValues  %in% c(TRUE, FALSE))) {
        msg <- c(msg, "In SMSTraits(): OutputValues must be provided for all SMS traits and must be either TRUE or FALSE.")
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
    cat("  SMS traits: \n")
    if(length(object@Positions) == 1) cat("     Trait types: (1) DP  \n")
    if(length(object@Positions) == 4) cat("     Trait types: (1) DP, (2) GoalBias, (3) AlphaDB, (4) BetaDB  \n")

    for (i in 1:length(object@Positions)){
        cat("     Configuration of SMS trait ", i, ": \n")
        if(is.numeric(object@Positions[i])) cat("     Loci positions coding for trait: ", object@Positions[i], "\n")
        if(!is.numeric(object@Positions[i]) && object@Positions[i]=="random") cat("    Loci positions coding for trait randomly chosen with ", object@NbOfPositions[i], " positions\n")
        cat("       Expression type: ", object@ExpressionType[i], "\n")
        cat("       Initial distribution: ", object@InitialDistribution[i], "\n")
        cat("       Initial parameter: ", object@InitialParameters[i,], "\n")
        cat("       IsInherited: ", object@IsInherited[i], "\n")
        cat("       Mutation distribution: ", object@MutationDistribution[i], "\n")
        cat("       Mutation parameters: ", object@MutationParameters[i,], "\n")
        cat("       Mutation rate: ", object@MutationRate, "\n")
        if(object@OutputValues[i]) cat("       Allel values for gene is written to output. \n")
    }
})

### SUBCLASS CorrRWTRAITS

#' Set genetic traits structure for CorrRW traits
#'
#' @description
#'
#' Depending on the settings of \code{\link[RangeShiftR]{CorrRW}}, the following traits can be evolvable:
#'
#' - The length of a single transfer step \code{Steplength} \cr
#' - The angle correlation parameter of the random walk \code{Rho} \cr
#'
#' The entries of the trait parameters must be provided in the same order as the CorrRW traits are listed above (first \code{StepLength} and second \code{Rho}). If parameters expect a matrix, the rows must match the order of CorrRW traits listed above.
#'
#' @details
#'
#' Traits set to evolve cannot simultaneously be stage-dependent.
#'
#' The alleles of each trait can be expressed according to an additive model (allele values across all loci are summed) or be averaged.
#'
#' Mutations are additive and can be sampled in either a \code{normal} or \code{uniform} distribution.
#'
#' Initial allele values are sampled in a \code{normal} or \code{uniform} distribution.
#'
#' CorrRW traits can also be **not** inherited, that is, allele values are resampled from the initial distribution for every new individual.
#'
#' Dominance values are not applicable for CorrRW traits.
#'
#' @usage CorrRWTraits(Positions = list("random","random"), NbOfPositions = c(10, 10),
#' ExpressionType = rep("additive",2),
#' InitialDistribution = rep("normal",2), InitialParameters = matrix(c(rep(0.5,2),(rep(0.1,2), nrow=2),
#' IsInherited = rep(FALSE,2),
#' MutationDistribution = rep("normal",2), MutationParameters = matrix(c(rep(0.5,2),(rep(0.2,2), nrow=2),
#' MutationRate = rep(0.0001,2), OutputValues = rep(FALSE,2))
#'
#' @param Positions Loci positions coding for the trait within genome. Should be provided as a list. Entries can either be a string (\code{"random"}) and/or vectors of integers.
#' The length must be 2 for the required CorrRW traits \code{Steplength} and \code{Rho}.
#' @param NbOfPositions Only specify when the \code{Positions} of the CorrRW trait is set to ‘random’, else must be blank (NULL).
#' The length must be 2 for the required CorrRW traits \code{Steplength} and \code{Rho}.
#' @param ExpressionType Type of expression for the emigration trait. Can be either \code{additive} or \code{average}.
#' The length must be 2 for the required CorrRW traits \code{Steplength} and \code{Rho}.
#' @param InitialDistribution Distribution of the initial values. Can be \code{uniform} or \code{normal}. Should be provided as a vector of strings.
#' The length must be 2 for the required CorrRW traits \code{Steplength} and \code{Rho}.
#' @param InitialParameters Parameters for the initial distribution: You must provide two colums min and max  for \code{uniform} distribution and mean and sd for \code{normal} distribution.
#' Each row in the matrix corresponds to an CorrRW trait. The number of rows must be 2 for the required CorrRW traits \code{Steplength} and \code{Rho}.
#' @param IsInherited Should the emigration trait be inherited? Can be either \code{TRUE} or \code{FALSE}.
#' The length must be 2 for the required CorrRW traits \code{Steplength} and \code{Rho}.
#' @param MutationDistribution Distribution for mutations to draw from. Can be \code{uniform} or \code{normal}.
#' The length must be 2 for the required CorrRW traits \code{Steplength} and \code{Rho}.
#' @param MutationParameters Parameters for the mutation distribution: You must provide two colums: min and max for \code{uniform} distribution and mean and sd for \code{normal} distribution.
#' Each row in the matrix corresponds to an CorrRW trait. The number of rows must be 2 for the required CorrRW traits \code{Steplength} and \code{Rho}.
#' @param MutationRate Mutation rate applicable to this type of loci. Must be between 0.0 and 1.0.
#' The length must be 2 for the required CorrRW traits \code{Steplength} and \code{Rho}.
#' @param OutputValues If OutputGeneValues in GeneticsFile is enabled, should allele values for this gene be written to output? Ignored if OutputGeneValues is set to \code{FALSE}.
#' The length must be 2 for the required CorrRW traits \code{Steplength} and \code{Rho}.
#'
#'
#' @return a parameter object of class "CorrRWTraitsParams"
#' @author Jette Reeg
#' @name CorrRWTraits
#' @export CorrRWTraits
CorrRWTraits<- setClass("CorrRWTraitsParams", slots = c(Positions = "list", #
                                                        NbOfPositions = "ANY", # random or list of integer values
                                                        ExpressionType = "character", # additive or average
                                                        InitialDistribution = "character", # uniform or normal
                                                        InitialParameters = "matrix", # min and max value or mean and sd
                                                        IsInherited = "logical", # T/F
                                                        MutationDistribution = "character", # uniform or normal
                                                        MutationParameters = "matrix", # min mx or mean sd
                                                        MutationRate = "numeric", # float
                                                        OutputValues = "logical")
                        , prototype = list(
                            Positions = list("random", "random"), # "random" or list of integer values
                            NbOfPositions = rep(2,2), # numeric, only of positions random
                            ExpressionType = rep("additive",2), # dispersal: "additive" or "average"
                            InitialDistribution = rep("uniform",2), # uniform , normal (dispersal)
                            InitialParameters = matrix(c(0.5,0.5,0.1,0.1), nrow=2), # dispersal: two values: either min/max oder mean+sd
                            IsInherited = rep(FALSE,2), # only for dispersal
                            MutationDistribution = rep("uniform",2), # dispersal: uniform or normal
                            MutationParameters = matrix(c(0.5,0.5,0.1,0.1), nrow=2), # single value or 2 values
                            MutationRate = rep(0.001,2), # numeric
                            OutputValues = rep(FALSE,2)
                        ))
setValidity("CorrRWTraitsParams", function(object) {
    msg <- NULL
    # Check Position and NbOfPositions
    # Positions must be of type list?
    if(class(object@Positions) != "list") {
        msg <- c(msg, "In CorrRWTraits(): Positions must be provided as a list.")
    }
    # NbOfPositions must be either numeric, integer or NULL
    if (!is.null(object@NbOfPositions) && class(object@NbOfPositions) != "numeric" && class(object@NbOfPositions) != "integer") {
        msg <- c(msg, "In CorrRWTraits(): NbrOfPositions must be either NULL (if all positions are given) or numeric (if at least one CorrRW trait has random positions).")
    }
    if(length(object@Positions) != 2){
        msg <- c(msg, "In CorrRWTraits(): You must provide two positions for the CorrRW traits Steplength and Rho.")
    }
    if (all(object@Positions == "random")){ # if all positions are random
        if(length(object@NbOfPositions[!is.na(object@NbOfPositions)]) != length(object@Positions) ) {
            msg <- c(msg, "In CorrRWTraits(): For each CorrRW trait with random positions you must provide the number of positions.")
        }
    } else{ # if NOT all positions are random
        isNumeric <- sapply(object@Positions, is.numeric)
        if (!all(isNumeric)) { # if not all positions are numeric,
            if(object@Positions[isNumeric==FALSE] != "random"){ # then those not numeric must be random
                msg <- c(msg, "In CorrRWTraits(): Positions in CorrRW traits must be either a vector of integers or random.")
            }
            if (any(is.na(object@NbOfPositions[object@Positions == "random"])) || any(object@NbOfPositions[object@Positions == "random"] <= 0)){ # if number of positions are NA or smaller than 0
                msg <- c(msg, "In CorrRWTraits(): NbrOfPositions must be set to a strictly positive integer for random positions.")
            }
            if (any(!is.na(object@NbOfPositions[object@Positions != "random"]))) { # if there are NbOfPositions supplied for non-random positions
                msg <- c(msg, "In CorrRWTraits(): if Positions is not random NbrOfPositions must be not be set (NA).")
            }
        }
        else { # if all positions are not random
            if (!is.null(object@NbOfPositions)) {
                msg <- c(msg, "In CorrRWTraits(): If positions are not random, you must not specify the number of positions (NbOfPositions).")
            }
        }
    }

    # Check ExpressionType must be additive or average
    if(length(object@ExpressionType) != 2){
        msg <- c(msg, "In CorrRWTraits(): You must provide the ExpressionType for each CorrRW trait.")
    }
    if (!all(object@ExpressionType %in% c("additive", "average"))) {
        msg <- c(msg, "In CorrRWTraits(): ExpressionType must be either additive or average.")
    }

    # Check InitialDistribution and InitialParameter: Distribution must be uniform or normal and the length of the list must be the same as the number of positions
    if (length(object@InitialDistribution) != 2){
        msg <- c(msg, "In CorrRWTraits(): You must provide the InitialDistribution for each CorrRW trait.")
    }
    if (!all(object@InitialDistribution %in% c("uniform", "normal"))) {
        msg <- c(msg, "In CorrRWTraits(): InitialDistribution must be either normal, or uniform.")
    }

    if (nrow(object@InitialParameters) != 2) {
        msg <- c(msg, "In CorrRWTraits(): For each CorrRW parameter you must provide the InitialParameters. Use one row for each CorrRW parameter.")
    } else {
        # two columns are necessary for mean and sd or min and max
        if (ncol(object@InitialParameters) !=2 || # if DominanceParameters has not 2 columns OR
            any(!is.numeric(object@InitialParameters)) || # if entries are not numeric
            any(is.na(object@InitialParameters))) { # if entries are NA
            msg <- c(msg,"In CorrRWTraits(): For the initial distributions, InitialParams must provide two values for mean (normal) or min (uniform) (first column) and sd (normal) or max (uniform) (second column)")
        }
    }

    # Check IsInherited: must be TRUE or FALSE
    if (length(object@IsInherited) != 2){
        msg <- c(msg, "In CorrRWTraits(): You must provide IsInherited for each CorrRW trait.")
    }
    if (!all(is.logical(object@IsInherited))) {
        msg <- c(msg, "In CorrRWTraits(): IsInherited must be either TRUE or FALSE." )
    }

    # Check mutation rate
    if(!is.numeric(object@MutationRate) || length(object@MutationRate) != 2){
        msg <- c(msg, "In CorrRWTraits(): You must provide the mutation rate for each CorrRW trait as a numeric vector.")
    } else {
        if (!is.numeric(object@MutationRate) ||  (any(object@MutationRate < 0.0) || any(object@MutationRate > 1.0))) {
            msg <- c(msg, "In CorrRWTraits(): MutationRate must be between 0.0 and 1.0.")
        }
    }

    # Check MutationDistribution and MutationParameters
    if (!is.null(object@MutationDistribution)){
        if (length(object@MutationDistribution) != 2){
            msg <- c(msg, "In CorrRWTraits(): For each CorrRW trait you must provide the MutationDistribution.")
        } else if (nrow(object@MutationParameters) != 2) {
            msg <- c(msg, "In CorrRWTraits(): For each CorrRW trait you must provide the MutationParameters.")
        } else if(!all(object@MutationDistribution %in% c("uniform", "normal"))){
            msg <- c(msg, "In CorrRWTraits(): MutationDistribution must be either normal or uniform for CorrRW traits.")
        }  else {
            if (ncol(object@MutationParameters) !=2){
                msg <- c(msg,"In CorrRWTraits(): MutationParams must provide two values for uniform and normal distribution: min/mean (first column) and max/sd (second column)")
            }

            if (any(!is.numeric(object@MutationParameters)) ||
                any(is.na(object@MutationParameters))) {
                msg <- c(msg,"In CorrRWTraits(): For a uniform or normal mutation distribution, MutationParams must provide two values for min (first column) and max (second column)")
            }
        }
    }


    # Check OutputValues
    if (length(object@OutputValues) != 2 && !all(object@OutputValues  %in% c(TRUE, FALSE))) {
        msg <- c(msg, "In CorrRWTraits(): OutputValues must be provided for all CorrRW traits and must be either TRUE or FALSE.")
    }

    if (is.null(msg)) TRUE else msg
})
setMethod("initialize", "CorrRWTraitsParams", function(.Object, ...) {
    this_func = "CorrRWTraits(): "
    args <- list(...)
    .Object <- callNextMethod()
    if ( length(args) == 0 ) {
        validObject(.Object)
    }
    .Object
})
setMethod("show", "CorrRWTraitsParams", function(object){
    cat("  CorrRW traits: \n")
    cat("     Trait types: (1) StepLength, (2) Rho  \n")

    for (i in 1:length(object@Positions)){
        cat("     Configuration of CorrRW trait ", i, ": \n")
        if(is.numeric(object@Positions[i])) cat("     Loci positions coding for trait: ", object@Positions[i], "\n")
        if(!is.numeric(object@Positions[i]) && object@Positions[i]=="random") cat("    Loci positions coding for trait randomly chosen with ", object@NbOfPositions[i], " positions\n")
        cat("       Expression type: ", object@ExpressionType[i], "\n")
        cat("       Initial distribution: ", object@InitialDistribution[i], "\n")
        cat("       Initial parameter: ", object@InitialParameters[i,], "\n")
        cat("       IsInherited: ", object@IsInherited[i], "\n")
        cat("       Mutation distribution: ", object@MutationDistribution[i], "\n")
        cat("       Mutation parameters: ", object@MutationParameters[i,], "\n")
        cat("       Mutation rate: ", object@MutationRate, "\n")
        if(object@OutputValues[i]) cat("       Allel values for gene is written to output. \n")
    }
})

### SUBCLASS TRAITSPARAMS

# define this ClassUnion so that the 'Neutral' slot in the genetic class 'GeneticParams' can be FALSE if neutral genetics should not be modelled, etc.
setClassUnion("NeutralSlot", c("logical", "NeutralTraitsParams"))
setClassUnion("GeneticLoadSlot", c("logical", "GeneticLoadParams"))
setClassUnion("EmigrationTraitsSlot", c("logical", "EmigrationTraitsParams"))
setClassUnion("SettlementTraitsSlot", c("logical", "SettlementTraitsParams"))
setClassUnion("CorrRWTraitsSlot", c("logical", "CorrRWTraitsParams"))
setClassUnion("SMSTraitsSlot", c("logical", "SMSTraitsParams"))
setClassUnion("KernelTraitsSlot", c("logical", "KernelTraitsParams"))


#' Set genetic traits structure
#'
#' @description
#' Three types of traits can be made evolvable and parameterised with their own genetic architecture:
#'
#' - Dispersal traits correspond to the main parameters controlling each phase of dispersal (emigration, transfer and settlement). \cr
#' - Genetic fitness traits represent genetic load, the accumulation of deleteriousmutations and their effect on the viability of newborn offspring. \cr
#' - Neutral trait does not have any phenotypic effect during the simulation. It is used to compute F-statistics and other measures of neutral variation. \cr
#'
#'
#' @usage Traits(Neutral = FALSE,
#' GeneticLoad = FALSE,
#' EmigrationGenes = FALSE,
#' SettlementGenes = FALSE,
#' CorrRWGenes = FALSE,
#' SMSGenes = FALSE,
#' KernelGenes = FALSE)
#'
#' @param Neutral If neutral traits should be modelled, define \code{\link[RangeShiftR]{NeutralTraits}}. If \code{FALSE} (default), the neutral traits are not modelled.
#' @param GeneticLoad If genetic load should be modelled, define \code{\link[RangeShiftR]{GeneticLoadTraits}}. There can be up to 5 genetic loads. If \code{FALSE} (default), the genetic load is not modelled.
#' @param EmigrationGenes If evolvable emigration traits should be modelled (\code{IndVar = TRUE} in \code{\link[RangeShiftR]{Emigration}}), define \code{\link[RangeShiftR]{EmigrationTraits}}. If \code{FALSE} (default), the emigration traits are not modelled/evolvable.
#' @param SettlementGenes If evolvable settlement traits should be modelled (\code{IndVar = TRUE} in \code{\link[RangeShiftR]{Settlement}}), define \code{\link[RangeShiftR]{SettlementTraits}}. If \code{FALSE} (default), the settlement traits are not modelled/evolvable.
#' @param CorrRWGenes If evolvable CorrRW traits should be modelled (\code{IndVar = TRUE} in \code{\link[RangeShiftR]{CorrRW}}), define \code{\link[RangeShiftR]{CorrRWTraits}}. If \code{FALSE} (default), the CorrRW traits are not modelled/evolvable.
#' @param SMSGenes If evolvable SMS traits should be modelled (\code{IndVar = TRUE} in \code{\link[RangeShiftR]{SMS}}), define \code{\link[RangeShiftR]{SMSTraits}}. If \code{FALSE} (default), the SMS traits are not modelled/evolvable.
#' @param KernelGenes If evolvable Kernel traits should be modelled (\code{IndVar = TRUE} in \code{\link[RangeShiftR]{Kernel}}), define \code{\link[RangeShiftR]{KernelTraits}}. If \code{FALSE} (default), the Kernel traits are not modelled/evolvable.
#'
#' @details
#' The parameterisation structure for each type of trait is defined in the corresponding subclasses:
#'
#' - Dispersal traits: \code{\link[RangeShiftR]{EmigrationTraits}}, \code{\link[RangeShiftR]{SettlementTraits}}, \code{\link[RangeShiftR]{CorrRWTraits}}, \code{\link[RangeShiftR]{SMSTraits}} and \code{\link[RangeShiftR]{KernelTraits}} \cr
#' - Genetic fitness traits: \code{\link[RangeShiftR]{GeneticLoadTraits}} \cr
#' - Neutral traits: \code{\link[RangeShiftR]{NeutralTraits}} \cr
#'
#' They follow a similar structure, but some parameters and their specific options are specific to the type of trait.
#'
#' -	The number and positions of genes controlling the trait. \cr
#' -	A rule for expression (restricted to dispersal traits) \cr
#' -	A mutation rate for all genes controlling the trait. \cr
#' -	A distribution to sample mutations from. \cr
#' -	A distribution to sample initial values from. \cr
#' -	A distribution to sample dominance coefficients from (restricted to genetic fitness traits). \cr
#'
#'
#'
#' @return a parameter object of class "TraitsParams"
#' @author Jette Reeg
#' @name Traits
#' @export Traits
Traits <- setClass("TraitsParams", slots = c(Neutral = "NeutralSlot",
                                             GeneticLoad = "GeneticLoadSlot",
                                             EmigrationGenes = "EmigrationTraitsSlot",
                                             SettlementGenes = "SettlementTraitsSlot",
                                             CorrRWGenes = "CorrRWTraitsSlot",
                                             SMSGenes = "SMSTraitsSlot",
                                             KernelGenes = "KernelTraitsSlot"
                                            )
                       , prototype = list(Neutral = FALSE, # NeutralTraits(),
                                          GeneticLoad = FALSE, # GeneticLoadTraits(), # could this also be a list of multiple GeneticLoads?
                                          EmigrationGenes = FALSE, # EmigrationTraits(), # NULL E_D0, E_Alpha, E_Beta
                                          SettlementGenes = FALSE, # SettlementTraits(), # NULL, S_S0, S_Alpha, S_Beta
                                          CorrRWGenes = FALSE, # CorrRWTraits(), # NULL, # CRW_STEPLENGTH, CRW_STEPCORRELATION
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
    # Check CorrRWGenes
    if (!is.logical(object@CorrRWGenes) && !is(object@CorrRWGenes, "CorrRWTraitsParams")) {
        msg <- c(msg, "In Traits(): CorrRWGenes must be of class CorrRWTraitsParams.")
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
    if(class(object@CorrRWGenes) == "CorrRWTraitsParams") print(object@CorrRWGenes)
    if(class(object@SMSGenes) == "SMSTraitsParams") print(object@SMSGenes)
    if(class(object@KernelGenes) == "KernelTraitsParams") print(object@KernelGenes)
})


#' Set Genetics parameters
#'
#' @description
#'
#' Controls heritability and evolution of traits (in case of dispersal traits only if inter-individual variability is enabled (\code{IndVar=TRUE})).
#' Provides control over the genome size, the number of loci, the recombination rate, output options and the genetic traits to be modelled
#'
#' @usage Genetics(GenomeSize = 10,
#'          ChromosomeEnds = 1, RecombinationRate = 0.0,
#'          OutputGeneValues = FALSE,
#'          OutputFstatsWeirCockerham = FALSE, OutputFstatsWeirHill = FALSE,
#'          OutputStartGenetics = NULL, OutputInterval = NULL,
#'          PatchList = NULL, NbrPatchesToSample = NULL,
#'          nIndividualsToSample = NULL, Stages = NULL, Traits = Traits()
#'          )
#'
#' @param GenomeSize Maximum size of genome (number of loci)
#' @param ChromosomeEnds Where the genome is split into chromosomes, if empty
#' assumed one chromosome is equal to GenomeSize. These areas recombine with a
#' probability of 0.5. Disabled for haploid organisms (no recombination).
#' @param RecombinationRate Recombination rate (through chromosomal crossover)
#' across the whole genome (in addition to the chromosomeEnds above). Disabled for haploid organisms (no recombination).
#' @param OutputGeneValues Output the values of all alleles for all genes of all sampled individuals.
#' Does not output the resulting trait values: mean and SD of dispersal and genetic fitness traits are
#' output in the TraitsXPatch, TraitsXCell and/or TraitsXrow output files. Enables the geneValues output files.
#' @param OutputFstatsWeirCockerham Calculate F-statistics. Enables the neutralGenetics and
#' perLocusNeutralGenetics output files.
#' @param OutputFstatsWeirHill Calculate F-statistics. Enables the neutralGenetics and pairwisePatchNeutralGenetics output files.
#' @param OutputStartGenetics Which year should RangeShifter start to produce the output files listed above?
#' @param OutputInterval How frequently to output genetic output, including gene values and neutral statistics.
#' @param PatchList Which patches are to be sampled for output.  Patches can be
#' specified according to their patch number, as per the patch layer in a patch-based
#' model. Or sampled randomly or all patches can be chosen. In a cell-based landscape
#' random is the only option with number of patches (=cells) specified.
#' @param NbrPatchesToSample If PatchList=random or random_occupied then this specifies
#' the number of patches to sample randomly. Random: The chosen sample patches remain
#' the same throughout the simulation, i.e. do not vary between years or replicates
#' unless artificially generated landscape that is generated afresh between replicates.
#' Random_occupied: patches are re-sampled every generation among all patches containing at least 1 individual.
#' @param nIndividualsToSample The number of individuals to sample in a patch. If nInds < nIndividualsToSample then sampled individuals = nInds
#' @param Stages The age stages to sample from.
#' @param Traits The genetic traits to be modelled. See \code{\link[RangeShiftR]{Traits}} for more information.

#' @details
#' The genome itself is not modelled explicitly, and is instead represented by a genome size,
#' the set of gene positions and a set of positions for the chromosome breaks (hence determining the number of chromosomes).
#'
#' The rate of recombination applies to the entire genome. Genetic linkage does occur, based on the distances between genes controlling different traits.
#' Multiple traits can use the same position, making the alleles of such genes completely linked, but mutations and expression of different traits are resolved independently,
#' such that pleiotropy is not possible.
#' If either option for sexual reproduction is selected, all individuals are diploid. Otherwise, asexual individuals are always haploid (and do not recombine).
#'
#' \emph{Output} \cr
#'
#' \emph{Sampling} \cr
#' The volume of genetic output grows quickly with the number of individuals in the simulation, and it is therefore crucial to first constitute an appropriate sample of the community. \cr
#' First, a set of patches is sampled from either a random subset of all patches in the landscape, or a pre-specified list of patches.
#' In either case, the same patches will be sampled through the simulation, and no check is conducted to verify if each patch does contain a population or is empty.
#' There is an additional option to sample patches from a random subset of occupied patches, which will change from one generation to the next. \cr
#' Second, a number of individuals are sampled from the population of each sampled patch. In stage-structured populations,
#' it is possible to select only certain stages to be sampled.\cr
#'
#' \emph{Allele values}\cr
#'
#' This output, if enabled (\code{OutputGeneValues = TRUE}), will write all alleles values of the sampled individuals for the selected trait(s) to \emph{Sim<sim_nb>_Land<landscape_nb>_Rep<replicate_nb>_geneValues},
#' along with some contextual information.\cr
#'
#' Each row corresponds to a single gene:\cr
#'
#'     1.	Year\cr
#'     2.	Generation\cr
#'     3.	Individual ID\cr
#'     4.	Trait type (e.g. “kernel_meanDist1”, matching the value given as input)\cr
#'     5.	Position in the genome\cr
#'     6.	Value of the allele on the first chromosome\cr
#'     7.	 Dominance coefficient for the allele on the first chromosome (this will be 0 except for genetic fitness traits)\cr
#'     8.	(If diploid) Value of the allele on the second chromosome\cr
#'     9.	 (If diploid) Dominance coefficient for the allele on the second chromosome\cr
#'
#' \emph{Neutral genetics} \cr
#'
#' The standard neutral genetics output, Sim<sim_nb>_Land<landscape_nb>_neutralGenetics, writes the following entries (one row per generation): \cr
#' 1.	Replicate number \cr
#' 2.	Year \cr
#' 3.	Generation \cr
#' 4.	Number of sampled patches with a non-zero population \cr
#' 5.	Total number of sampled individuals \cr
#' 6.	Standard Fst  (Cockerham’s θ) \cr
#' 7.	Standard Fis (Cockerham’s f) \cr
#' 8.	Standard Fit (Cockerham’s F) \cr
#' 9.	Global allelic diversity, calculated as the mean number of neutral alleles per locus for the entire sample. \cr
#' 10.	Local allelic diversity, calculated as the mean number of neutral alleles per locus for each sampled patch, then averaged over patches. \cr
#' 11.	Number of globally fixed alleles. \cr
#' 12.	Mean number of fixed alleles per patch. Note that this may differ from the number of globally fixed alleles, for example if one allele is fixed in a given patch but polymorphism exists in other patches. \cr
#' 13.	Observed heterozygosity Ho, calculated as the mean number of heterozygous loci per individual per locus. \cr
#'
#' RangeShifter estimates standard F-statistics as \insertCite{cockerham1969}{RangeShiftR}’s \ifelse{html}{\out{&theta;}}{\eqn{θ}} statistics, using 1) the classic method-of-moments estimator of \insertCite{weir1984}{RangeShiftR} (\code{OutputFstatsWeirCockerham=TRUE}) and/or
#' 2) the unequal sample-size generalization and extensions from \insertCite{weir2002}{RangeShiftR} (\code{OutputFstatsWeirHill=TRUE}). \cr
#'
#' In short, \ifelse{html}{\out{&theta;}}{\eqn{θ}} (Fst) measures the correlation between alleles within sub-populations (patches) relative to the complete sampled population,
#' f (Fis) the correlation between alleles within individuals relative to the sub-population and F (Fit) the correlation between alleles
#' within individuals relative to the complete sampled population (see \insertCite{holsinger2009}{RangeShiftR} for an introduction). \cr
#'
#' See the RangeShifter manual for more details on the calculation of these statistics. \cr
#'
#' \emph{Per-locus neutral genetics} \cr
#'
#' If the Weir and Cockerham method is enabled (\code{OutputFstatsWeirCockerham=TRUE}), RangeShifter outputs an additional file \emph{Sim<sim_nb>_Land<landscape_nb>_perLocusNeutralGenetics} with one entry for each neutral locus:\cr
#'
#' 1.	Year \cr
#' 2.	RepSeason \cr
#' 3.	Locus, the ID of locus \eqn{l} \cr
#' 4.	\ifelse{html}{\out{F<sub>st</sub>}}{\eqn{F_st}}, the value of \ifelse{html}{\out{F<sub>st,l</sub>}}{\eqn{F_st,l}} for locus \eqn{l} \cr
#' 5.	\ifelse{html}{\out{F<sub>is</sub>}}{\eqn{F_is}}, the value of \ifelse{html}{\out{F<sub>is</sub>}}{\eqn{F_is,l}} for locus \eqn{l} \cr
#' 6.	\ifelse{html}{\out{F<sub>it</sub>}}{\eqn{F_it}}, the value of \ifelse{html}{\out{F<sub>it</sub>}}{\eqn{F_it,l}} for locus \eqn{l} \cr
#' 7.	Het, the sample-level observed heterozygosity (\ifelse{html}{\out{H<sub>o</sub>}}{\eqn{H_o}}) for locus \eqn{l} \cr
#' 8.	One column \emph{patch_<i>_het} for each patch \emph{i} in the sample, indicating \ifelse{html}{\out{H<sub>o</sub>}}{\eqn{H_o}} for patch \emph{i} and locus \eqn{l} \cr
#'
#' \emph{Pairwise patch neutral genetics}\cr
#'
#' If the Weir and Hill method is enabled (\code{OutputFstatsWeirHill=TRUE}), a pairwise \ifelse{html}{\out{F<sub>st</sub>}}{\eqn{F_st}} matrix is constituted and filled with the corresponding
#' values of \ifelse{html}{\out{&beta;<sub>ii’</sub>}}{\eqn{β_ii’}} for each pair of patches in the sample. Values of ifelse{html}{\out{&beta;<sub>i</sub>}}{\eqn{β_i}} are also computed along the diagonal. \cr
#'
#' In this case, there is one row of output for each pair: \cr
#'     1.	Year\cr
#'     2.	Generation\cr
#'     3.	Patch ID of the first patch\cr
#'     4.	Patch ID of the second patch (same as 3. along the diagonal of the matrix) \cr
#'     5.	Pairwise \ifelse{html}{\out{F<sub>st</sub>}}{\eqn{F_st}} \cr
#'
#' \emph{Traits} \cr
#'
#' In the case of inter-individual variability and evolution of the dispersal traits,
#' it is possible to output the mean traits of the population. There are two types of traits output:\cr
#'
#' 1.	\emph{Mean traits by cell/patch (Sim0_TraitsXcell.txt or Sim0_TraitsXpatch.txt}). This file reports mean and standard deviation of the varying traits for each cell/patch,
#' for each replicate and reproductive season at the set year interval. \cr
#' 2.	\emph{Mean traits by row (Sim0_TraitsXrow.txt)}. The mean and standard deviation of the varying traits are computed at the row (\emph{y}) level,
#' pulling together all the populations occupying cells in \emph{y}. Values are reported for each replicate and reproductive season at the specified yearly interval.
#' This is particularly useful for analysing the structuring of traits along latitudinal gradients. It is possible to compute this output only for cell-based models. \cr
#'
#' Data for these outputs are collected at the same time as for the range and population outputs, i.e. before reproduction at each reproductive season at the set year interval and at the end of the simulation. \cr
#'
#' @references
#'         \insertAllCited{}
#'
#' @return a parameter object of class "GeneticsParams"
#' @author Jette Reeg
#' @name Genetics
#' @export Genetics
Genetics <- setClass("GeneticsParams", slots = c(GenomeSize = "integer_OR_numeric",
                                                 ChromosomeEnds = "ANY", # NULL or vector
                                                 RecombinationRate = "integer_OR_numeric", # NULL or numeric
                                                 OutputGeneValues = "logical",
                                                 OutputFstatsWeirCockerham = "logical",
                                                 OutputFstatsWeirHill = "logical",
                                                 OutputStartGenetics = "ANY", # positive integer if any output is TRUE or NULL
                                                 OutputInterval = "ANY",
                                                 PatchList = "ANY", # vector of integers or a string
                                                 NbrPatchesToSample = "ANY", # NULL or integer
                                                 nIndividualsToSample = "ANY", # "character_OR_integer", # character or integer
                                                 Stages = "ANY", # "character_OR_integer", # vector
                                                 Traits = "TraitsParams")
                                     , prototype = list(GenomeSize = 1L,
                                                        ChromosomeEnds = 0L, # NULL or vector
                                                        RecombinationRate = 0.0, # NULL or numeric
                                                        OutputGeneValues = FALSE,
                                                        OutputFstatsWeirCockerham = FALSE,
                                                        OutputFstatsWeirHill = FALSE,
                                                        OutputStartGenetics = NULL, # positive integer if any output is TRUE or NULL
                                                        OutputInterval = NULL,
                                                        PatchList = NULL, #"all", # vector or string
                                                        NbrPatchesToSample = NULL, #0L, # NULL or integer
                                                        nIndividualsToSample = NULL, #"all", # NULL or integer
                                                        Stages = NULL, #"all", # vector
                                                        Traits = Traits())

)
setValidity('GeneticsParams', function(object){
    msg <- NULL
    # Genome Size
    # must be numeric and >0
    if (!is.numeric(object@GenomeSize) || object@GenomeSize <= 0) {
        msg <- c(msg, "In Genetics(): GenomeSize must be a positive integer.")
    }
    # Chromosome Ends
    # must be NULL or numeric vector
    if (!is.null(object@ChromosomeEnds) && !is.numeric(object@ChromosomeEnds)) {
        msg <- c(msg, "In Genetics(): ChromosomeEnds must be a numeric vector or NULL.")
    }

    # Recombination Rate
    # should be checked in RSparams: should not be set for asexual models
    # should be between 0 and 0.5
    if (!is.numeric(object@RecombinationRate) || object@RecombinationRate < 0.0 || object@RecombinationRate > 0.5) {
        msg <- c(msg, "In Genetics(): RecombinationRate must be a float between 0.0 and 0.5.")
    }

    # OutputGeneValues
    # must be a boolean
    if (!is.logical(object@OutputGeneValues)) {
        msg <- c(msg, "In Genetics(): OutputGeneValues must be true or false.")
    }

    # OutputFstatsWeirCockerham
    # must be a boolean
    if (!is.logical(object@OutputFstatsWeirCockerham)) {
        msg <- c(msg, "In Genetics(): OutputFstatsWeirCockerham must be true or false.")
    }

    # OutputFstatsWeirHill
    # must be a boolean
    if (!is.logical(object@OutputFstatsWeirHill)) {
        msg <- c(msg, "In Genetics(): OutputFstatsWeirHill must be true or false.")
    }

    anyNeutral = object@OutputFstatsWeirCockerham || object@OutputFstatsWeirHill

    anyGeneticsOutput = object@OutputGeneValues == "TRUE" || anyNeutral

    if (anyGeneticsOutput) {
        if (is.null(object@OutputStartGenetics) || object@OutputStartGenetics < 0) {
            msg <- c(msg, "OutStartGenetics be greater than 0 if any genetic output option is TRUE.")
        }

        if (is.null(object@OutputInterval) || object@OutputInterval <= 0) { # check whether is.null()
            msg <- c(msg,"OutputInterval must be at least 1 if any genetic output option is TRUE.")
        }
    }
    else if (!is.null(object@OutputInterval) || !is.null(object@OutputStartGenetics)){
        msg <- c(msg, "OutStartGenetics and OutputInterval should be NULL if all genetic output options are FALSE.")
    }

    # Check PatchList
    if (anyGeneticsOutput) {
        if(!is.null(object@PatchList)){
        if (is.character(object@PatchList) || is.numeric(object@PatchList)) {
            if (!is.numeric(object@PatchList) && object@PatchList != "random" && object@PatchList != "all" && object@PatchList != "random_occupied") {
                msg <- c(msg,"PatchList must be either a vector of integers, or \"random\", \"random_occupied\" or \"all\".")
            }
        }
        } else{
            msg <- c(msg, "PatchList cannot be NULL if any genetic output option is TRUE.")
        }

    }
    else if (!is.null(object@PatchList)) {
        msg <- c(msg, "PatchList should be NULL if all genetic output options are FALSE.")
    }

    # Check NbrPatchesToSample
    if (anyGeneticsOutput) {
        if (is.character(object@PatchList)){
            if (object@PatchList == "random" || object@PatchList == "random_occupied") {
                if (is.null(object@NbrPatchesToSample) || object@NbrPatchesToSample == 0) {
                    msg <- c(msg, "NbrPatchesToSample cannot be NULL or 0 if PatchList is \"random\" or \"random_occupied\".")
                }
                else {
                    if (object@NbrPatchesToSample <= 0) {
                        msg <- c(msg, "NbrPatchesToSample must be a positive integer if PatchList is \"random\" or \"random_occupied\".")
                    }
                }
            }
        }
        else if (!is.null(object@NbrPatchesToSample) && object@NbrPatchesToSample != 0) {
            msg <- c(msg, "NbrPatchesToSample must be NULL or zero if PatchList is not \"random\" or \"random_occupied\".")
        }
    }

    # Check IndividualsToSample
    if (anyGeneticsOutput) {
        if (is.null(object@nIndividualsToSample)  || (is.numeric(object@nIndividualsToSample) &&  object@nIndividualsToSample == 0)) {
            msg <- c(msg, "nIndividualsToSample cannot be NULL or zero if any genetics output option is TRUE.")
        }
        else if ((is.character(object@nIndividualsToSample) &&  object@nIndividualsToSample != "all") || (is.numeric(object@nIndividualsToSample) && object@nIndividualsToSample <= 0)) {
            msg <- c(msg, "nIndividualsToSample must be either a positive integer or \"all\".")
        }
    }
    else if (!is.null(object@nIndividualsToSample)) {
        msg <- c(msg, "numberIndividualsToSample must be NULL if all genetics output options are FALSE.")
    }

    # Check Stages
    if (anyGeneticsOutput) {
        if (is.null(object@Stages)) {
            msg <- c(msg, "Stages cannot be NULL if any genetic output option is TRUE.")
        }
        else {
            if (!is.numeric(object@Stages) && (is.character(object@Stages) &&  object@Stages != "all")) {
                msg <- c(msg, "Stages must be either a vector of integers, or \"all\".")
            }
        }
    }
    else if (!is.null(object@Stages)) {
        msg <- c(msg, "Stages must be NULL if all genetic output options are FALSE.")
    }
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
        cat("     Genome is slitted into chromosomes at =", object@ChromosomeEnds, "\n ")
    } else cat( "     Genome is not slitted into chromosomes \n")

    cat("   Recombination rate: ", object@RecombinationRate, "\n")
    cat("   Output genetic values: ", object@OutputGeneValues, "\n")
    cat("   Output Fstats after Weir Cockerham: ", object@OutputFstatsWeirCockerham, "\n")
    cat("   Output Fstats after Weir Hill: ", object@OutputFstatsWeirHill, "\n")
    if(any(object@OutputGeneValues || object@OutputFstatsWeirCockerham || object@OutputFstatsWeirHill)){
        cat("   Start genetic output at year: ", object@OutputStartGenetics, "and output every ",object@OutputInterval ," year \n")
        cat("     Patches to sample: ", object@PatchList, "\n")
        if(object@PatchList=="random" || object@PatchList=="random_occupied"){
            cat("     Number of patches to sample: ", object@NbrPatchesToSample, "\n")
        }
        cat("     Number of individuals to sample: ", object@nIndividualsToSample, "\n")
        cat("     Stages to sample: ", object@Stages, "\n")
    } else{
        cat("   No genetic output. \n")
    }
    print(object@Traits)
})
