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

#' Set genetic traits structure for neutral traits
#'
#' @return a parameter object of class "NeutralTraitsParams"
#' @author Jette Reeg
#' @name NeutralTraits
#' @export NeutralTraits
NeutralTraits<- setClass("NeutralTraitsParams", slots = c(Positions = "integer_OR_numeric",
                                                         NbOfPositions = "integer_OR_numeric", # random or list of integer values
                                                         InitialDistribution = "character", # uniform
                                                         InitialParameters = "character", # max value
                                                         MutationDistribution = "character", # KAM or SSM
                                                         MutationParameters = "character", # max
                                                         MutationRate = "numeric", # float
                                                         OutputValues = "logical")
                   , prototype = list(Positions = NULL, # "random" or list of integer values
                                      NbOfPositions = NULL, # numeric, only of positions random
                                      InitialDistribution = NULL, # uniform (neutral + dispersal), normal (dispersal), NULL (genetic load)
                                      InitialParameters = NULL, # neutral traits: only max value; dispersal: two values: either min/max oder mean+sd, not applicable for genetic load
                                      MutationDistribution = NULL, # neutral: "KAM" or "SSM", genetic load: "gamma", "uniform", "normal", "negExp", dispersal: uniform or normal
                                      MutationParameters = NULL, # single value or 2 values
                                      MutationRate = NULL, # numeric
                                      OutputValues = FALSE
                   ))
setValidity("NeutralTraitsParams", function(object) {
    msg <- NULL
    patternPositions <-  "^\"?(([0-9]+-)?[0-9]+,)*([0-9]+-)?[0-9]+\"?$"
    # Check Position and NbOfPositions
    isMatch <- grepl(patternPositions, object@Positions)
    if (!isMatch && object@Positions != "random") {
        msg <- c(msg, "Positions in neutral genetics must be either a comma-separated list of integer ranges, or random.")
    }
    if (object@Positions == "random") {
        if (object@NbPositions <= 0) {
            msg <- c(msg, "NbrOfPositions must be a strictly positive integrer.")
        }
    } else if (!is.null(object@NbPositions)) {
        msg <- c(msg, "If Positions is not random NbrOfPositions must be not be set (NA).")
    }
    # Check InitialDistribution
    if (!is.null(object@InitialDistribution) && object@InitialDistribution != "uniform") {
        msg <- c(msg,"InitialDistribution must be either uniform or left blank (NA) for the neutral trait.")
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
    .Object
})
setMethod("show", "NeutralTraitsParams", function(object){

})


### SUBCLASS GENETICLOADTRAITS

#' Set genetic traits structure for genetic fitness
#'
#' @return a parameter object of class "GeneticLoadParams"
#' @author Jette Reeg
#' @name GeneticLoadTraits
#' @export GeneticLoadTraits
GeneticLoadTraits<- setClass("GeneticLoadParams", slots = c(
    NbGeneticLoads = "integer", # number of genetic loads
    Positions = "integer_OR_numeric",# "random" or list of integer values
    NbOfPositions = "integer_OR_numeric", # numeric, only of positions random
    ExpressionType = "integer_OR_numeric",# "multiplicative"
    DominanceDistribution = "character", # ‘gamma’, ‘uniform’, ‘normal’, ‘negExp’, ‘scaled’
    DominanceParameters = "character", # 2 values for min/max, mean/sd, shape/scale or one value: mean
    MutationDistribution = "character", # ‘gamma’, ‘uniform’, ‘normal’,‘negExp’
    MutationParameters = "character", #  2 values for min/max, mean/sd, shape/scale or one value: mean
    MutationRate = "numeric", # float
    OutputValues = "logical")
    , prototype = list(
        NbGeneticLoads = 1L,
        Positions = NULL,
        NbOfPositions = NULL,
        ExpressionType = "multiplicative",
        DominanceDistribution = NULL,
        DominanceParameters = NULL,
        MutationDistribution = NULL,
        MutationParameters = NULL,
        MutationRate = NULL,
        OutputValues = FALSE
    ))
setValidity("GeneticLoadParams", function(object) {
    msg <- NULL
    patternPositions <-  "^\"?(([0-9]+-)?[0-9]+,)*([0-9]+-)?[0-9]+\"?$"
    # only 5 GeneticLoads are allowed
    if (object@NbGeneticLoads < 1 || object@NbGeneticLoads > 5) {
        msg <- c(msg, "Number of genetic loads must be between 1 and 5.")
    }
    # Check Position and NbOfPositions
    if (length(object@Positions) != object@NbGeneticLoads) {
        msg <- c(msg, "For each genetic load you must provide the positions.")
    } else if (all(object@Positions == "random") && length(object@NbOfPositions) != object@NbGeneticLoads ) {
        msg <- c(msg, "For each genetic load you must provide the number of positions.")
        } else{
            isMatch <- grepl(patternPositions, object@Positions)
            if (!all(isMatch) && object@Positions[isMatch==FALSE] != "random") {
                msg <- c(msg, "Positions in genetic loads must be either a comma-separated list of integer ranges, or random.")
            }
            if (any(is.na(object@NbPositions[object@Positions == "random"])) || object@NbPositions[object@Positions == "random"] <= 0){
                msg <- c(msg, "NbrOfPositions must be set to a strictly positive integrer.")
            } else if (!is.na(object@NbOfPositions[object@Positions != "random"])) {
                msg <- c(msg, "In Genetic loads: if Positions is not random NbrOfPositions must be not be set (NA).")
            }
        }

    # Check ExpressionType
    if (length(object@ExpressionType)!=object@NbGeneticLoads){
        msg <- c(msg, "For each genetic load you must provide the ExpressionType.")
    } else {
        if (any(object@ExpressionType != "multiplicative")) {
        msg <- c(msg, "ExpressionType must be \"multiplicative\" for genetic load traits.")
        }
    }

    # Check DominanceDistribution
    if (!is.null(object@DominanceDistribution)){
        if(length(object@DepressionDistribution) != object@NbGeneticLoads) {
            msg <- c(msg, "For each genetic load you must provide the DominanceDistribution.")
        } else if (length(object@DominanceParameters) != object@NbGeneticLoads) {
            msg <- c(msg, "If you have set DominanceDistributions you must provide the DominanceParameters for each genetic load.")
        } else {
            if (any(object@DominanceDistribution == "normal")) { # if any distribution is normal
                # two values for mean and sd
                if (ncol(object@DominanceParameters[object@DominanceDistribution=="normal"]) !=2 || # if DominanceParameters has not 2 columns OR
                    any(!is.numeric(object@DominanceParameters[object@DominanceDistribution=="normal"])) || # if entries are not numeric
                    any(is.na(object@DominanceParameters[object@DominanceDistribution=="normal"]))) { # if entries are NA
                    msg <- c(msg,"For a normal dominance distribution, DominanceParams must provide two values for mean (first column) and sd (second column)")
                }
            }
            if (any(object@DominanceDistribution == "gamma")) {
                # two values for shape and scale
                if (ncol(object@DominanceParameters[object@DominanceDistribution=="gamma"]) !=2 || # if DominanceParameters has not 2 columns OR
                    any(!is.numeric(object@DominanceParameters[object@DominanceDistribution=="gamma"])) || # if entries are not numeric
                    any(is.na(object@DominanceParameters[object@DominanceDistribution=="gamma"]))) { # if entries are NA
                    msg <- c(msg,"For a gamma dominance distribution, DominanceParams must provide two values for shape (first column) and scale (second column)")
                }
            }
            if (any(object@DominanceDistribution == "uniform")) {
                # two values for min and max
                if (ncol(object@DominanceParameters[object@DominanceDistribution=="uniform"]) !=2 || # if DominanceParameters has not 2 columns OR
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
                        any(is.na(object@DominanceParameters[object@DominanceDistribution=="negExp",1]))
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
        if (nrow(object@MutationDistribution) != object@NbGeneticLoads){
            msg <- c(msg, "For each genetic load you must provide the MutationDistribution.")
        } else if (nrow(object@MutationParameters) != object@NbGeneticLoads) {
            msg <- c(msg, "For each genetic load you must provide the MutationParameters.")
        } else {
            if (any(object@MutationDistribution == "uniform")) {
                # two values for min and max
                if (ncol(object@MutationParameters[object@MutationDistribution=="uniform"]) !=2 ||
                    any(!is.numeric(object@MutationParameters[object@MutationDistribution=="uniform"])) ||
                    any(is.na(object@MutationParameters[object@MutationDistribution=="uniform"]))) {
                    msg <- c(msg,"For a uniform mutation distribution, MutationParams must provide two values for min (first column) and max (second column)")
                }
            }
            if (any(object@MutationDistribution == "normal")) {
                # two values for meand and sd
                if (ncol(object@MutationParameters[object@MutationDistribution=="normal"]) !=2 ||
                    any(!is.numeric(object@MutationParameters[object@MutationDistribution=="normal"])) ||
                    any(is.na(object@MutationParameters[object@MutationDistribution=="normal"]))) {
                    msg <- c(msg,"For a normal mutation distribution, MutationParams must provide two values for mean (first column) and sd (second column)")
                }
            }
            if (any(object@MutationDistribution == "gamma")) {
                # two values for shape and scale
                if (ncol(object@MutationParameters[object@MutationDistribution=="gamma"]) !=2 ||
                    any(!is.numeric(object@MutationParameters[object@MutationDistribution=="gamma"])) ||
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
EmigrationTraits<- setClass("EmigrationTraitsParams", slots = c(Positions = "integer_OR_numeric", #
                                                                NbOfPositions = "integer_OR_numeric", # random or list of integer values
                                                                ExpressionType = "integer_OR_numeric", # additive or average
                                                                InitialDistribution = "character", # uniform or normal
                                                                InitialParameters = "character", # min and max value or mean and sd
                                                                IsInherited = "logical", # T/F
                                                                MutationDistribution = "character", # uniform or normal
                                                                MutationParameters = "character", # min mx or mean sd
                                                                MutationRate = "numeric", # float
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
    if (any(is.na(object@NbPositions[object@Positions == "random"])) || object@NbPositions[object@Positions == "random"] <= 0){
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
SettlementTraits<- setClass("SettlementTraitsParams", slots = c(Positions = "integer_OR_numeric", #
                                                                 NbOfPositions = "integer_OR_numeric", # random or list of integer values
                                                                 ExpressionType = "integer_OR_numeric", # additive or average
                                                                 InitialDistribution = "character", # uniform or normal
                                                                 InitialParameters = "character", # min and max value or mean and sd
                                                                 IsInherited = "logical", # T/F
                                                                 MutationDistribution = "character", # uniform or normal
                                                                 MutationParameters = "character", # min mx or mean sd
                                                                 MutationRate = "numeric", # float
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
    if (any(is.na(object@NbPositions[object@Positions == "random"])) || object@NbPositions[object@Positions == "random"] <= 0){
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
CRWTraits<- setClass("CRWTraitsParams", slots = c(Positions = "integer_OR_numeric", #
                                                                 NbOfPositions = "integer_OR_numeric", # random or list of integer values
                                                                 ExpressionType = "integer_OR_numeric", # additive or average
                                                                 InitialDistribution = "character", # uniform or normal
                                                                 InitialParameters = "character", # min and max value or mean and sd
                                                                 IsInherited = "logical", # T/F
                                                                 MutationDistribution = "character", # uniform or normal
                                                                 MutationParameters = "character", # min mx or mean sd
                                                                 MutationRate = "numeric", # float
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
    if (any(is.na(object@NbPositions[object@Positions == "random"])) || object@NbPositions[object@Positions == "random"] <= 0){
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
KernelTraits<- setClass("KernelTraitsParams", slots = c(Positions = "integer_OR_numeric", #
                                                   NbOfPositions = "integer_OR_numeric", # random or list of integer values
                                                   ExpressionType = "integer_OR_numeric", # additive or average
                                                   InitialDistribution = "character", # uniform or normal
                                                   InitialParameters = "character", # min and max value or mean and sd
                                                   IsInherited = "logical", # T/F
                                                   MutationDistribution = "character", # uniform or normal
                                                   MutationParameters = "character", # min mx or mean sd
                                                   MutationRate = "numeric", # float
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
    if (any(is.na(object@NbPositions[object@Positions == "random"])) || object@NbPositions[object@Positions == "random"] <= 0){
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
SMSTraits<- setClass("SMSTraitsParams", slots = c(Positions = "integer_OR_numeric", #
                                                         NbOfPositions = "integer_OR_numeric", # random or list of integer values
                                                         ExpressionType = "integer_OR_numeric", # additive or average
                                                         InitialDistribution = "character", # uniform or normal
                                                         InitialParameters = "character", # min and max value or mean and sd
                                                         IsInherited = "logical", # T/F
                                                         MutationDistribution = "character", # uniform or normal
                                                         MutationParameters = "character", # min mx or mean sd
                                                         MutationRate = "numeric", # float
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
    if (any(is.na(object@NbPositions[object@Positions == "random"])) || object@NbPositions[object@Positions == "random"] <= 0){
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
    if (!is.logical(object@Neutral) && !is(object@GeneticLoad, "GeneticLoadParams")) {
        msg <- c(msg, "In Traits(): GeneticLoad must be of class GeneticLoadParams.")
    }
    # Check EmigrationGenes
    if (!is.logical(object@Neutral) && !is(object@EmigrationGenes, "EmigrationTraitsParams")) {
        msg <- c(msg, "In Traits(): EmigrationGenes must be of class EmigrationTraitsParams.")
    }
    # Check SettlementGenes
    if (!is.logical(object@Neutral) && !is(object@SettlementGenes, "SettlementTraitsParams")) {
        msg <- c(msg, "In Traits(): SettlementGenes must be of class SettlementTraitsParams.")
    }
    # Check CRWGenes
    if (!is.logical(object@Neutral) && !is(object@CRWGenes, "CRWTraitsParams")) {
        msg <- c(msg, "In Traits(): CRWGenes must be of class CRWTraitsParams.")
    }
    # Check SMSGenes
    if (!is.logical(object@Neutral) && !is(object@SMSGenes, "SMSTraitsParams")) {
        msg <- c(msg, "In Traits(): SMSGenes must be of class SMSTraitsParams.")
    }
    # Check KernelGenes
    if (!is.logical(object@Neutral) && !is(object@KernelGenes, "KernelTraitsParams")) {
        msg <- c(msg, "In Traits(): KernelGenes must be of class KernelTraitsParams.")
    }
    if (is.null(msg)) TRUE else msg
})
setMethod("initialize", "TraitsParams", function(.Object, ...) {
    .Object
})
setMethod("show", "TraitsParams", function(object){
})


#' Set Genetics parameters
#'
#' @description Set genetics parameters and architecture.\cr
#'
#' Controls heritability and evolution of traits (if inter-individual variability is enabled (\code{IndVar=TRUE}) for at least one (dispersal) trait).
#' Provides control over the number of chromosomes, the number of loci on each, the recombination rate (if the species is diploid) and a
#' flexible mapping of traits to chromosomes, allowing linkage, pleiotropy and neutral alleles to be incorporated. It is also possible to model
#' neutral alleles when no adaptive traits are present.
#'
#' @usage Genetics(Architecture = 0,
#'          NLoci = 1, ArchFile = "NULL",
#'          ProbMutn = 0.0, MutationSD = 0.1,
#'          ProbCross = 0.0,
#'          AlleleSD = 0.1)
#' @param Architecture Genetic Architecture: \cr 0 = One chromosome per trait (default), \cr 1 = Read from file (set \code{ArchFile}).
#' @param NLoci Required if \code{Architecture=0}: Number of loci per chromosome, defaults to \eqn{1} (integer).
#' @param ArchFile Required if \code{Architecture=1}: Name of the genetic architecture file.
#' @param ProbMutn Probability of mutation of each individual allele at meiosis, defaults to \eqn{0.0}.\cr
#' @param MutationSD Standard deviation of mutation magnitude, defaults to \eqn{0.1}.\cr Must be \eqn{> 0}.
#' Must be \eqn{0 \le}\code{ProbMutn} \eqn{\le 1}.
#' @param ProbCross Probability of crossover at each individual locus at meiosis, defaults to \eqn{0.0}.\cr
#' Must be \eqn{0 \le}\code{ProbCross} \eqn{\le 1}.
#' @param AlleleSD Standard deviation of initial allelic values around phenotypic value, defaults to \eqn{0.1}.\cr Must be \eqn{> 0}.
#' @details When inter-individual variability in at least one trait (at present limited to dispersal traits)
#' is enabled (\code{IndVar=TRUE}), each individual carries a genome coding for the varying trait values.\cr
#' If the reproductive model is \emph{asexual/female-only} (\code{ReproductionType} \eqn{=0}), the species is assumed to be haploid and chromosomes
#' hold a single allele at each locus. In this case, changes in the genotype, and hence also in the phenotype, can occur only through mutation (see below). \cr
#' In the case of \emph{sexual} models(\code{ReproductionType} \eqn{={1,2}}), the species is assumed to be diploid and chromosomes
#' hold two alleles at each locus. New-born offspring inherit one set of chromosomes from each of its parents.
#' Note that we use the term \emph{chromosome} here to represent both the single chromosomes of a haploid species and
#' the chromosome pair of a diploid species.\cr
#'
#' \emph{Simple genetic architecture}\cr
#' A simple type of architecture may be used by choosing the \emph{one chromosome per trait}
#' option (\code{Architecture} \eqn{=0}) and setting the number of loci per chromosome;
#' all chromosomes will carry the same number of loci.\cr
#'
#' In the 1.x versions of \emph{RangeShifter}, the trait value was held by the individual directly
#' as a ‘pseudo-gene’ (for diploid species, the mean of two such alleles controlled the phenotype),
#' but this did not allow for representation of more complex genetic phenomena, such as linkage
#' between traits. This simple type of implementation may still be represented approximately by
#' choosing the ‘one chromosome per trait’ option, and setting one locus per chromosome.
#'
#' \emph{Flexible genetic architecture}\cr
#' For a more realistic representation of heritable traits, an explicit genetic architecture
#' may be defined, which must be read from a text file (\code{Architecture} \eqn{=1}).
#' The file specifies how many chromosomes
#' each individual will have, the number of loci on each chromosome and which loci contribute
#' to the phenotypic value of each trait. Thus, for example, it is possible to model a species
#' exhibiting a single variable trait (e.g. the mean of a negative exponential dispersal kernel)
#' dependent on loci spreads across three chromosomes or a species exhibiting three trait (e.g.
#' density-dependent emigration) all of which are governed by loci located on a single chromosome.
#' In practice, most genetic architectures are likely to fall somewhere between these extremes,
#' i.e. there will be several chromosomes and traits will be mapped across them. In contrast to
#' \emph{RangeShifter} v1, whenever there are variable traits in a model, evolution is assumed, although
#' it can if desired be effectively eliminated for a haploid species by setting a very low mutation
#' probability (\code{ProbMutn}).\cr
#'
#' Care must be taken to specify the architecture file in the required format. This is an example
#' architecture file content:\cr\cr
#' NChromosomes	4\cr
#' NLoci	5 6 10 2\cr
#' Trait 0 NLoci 4 0 0 0 1 0 2 1 5\cr
#' Trait 1 NLoci 8 0 4 2 0 2 1 2 2 2 3 2 4 2 9 3 0\cr
#' Trait 2 NLoci 2 1 4 3 1 \cr
#'
#' The first line contains the keyword \code{NChromosomes} followed by an integer \eqn{n_C} that
#' specifies the number of chromosomes.\cr
#' The second line contains the keyword \code{NLoci} followed by \eqn{n_C} integers which specify
#' the number of loci on each chromosome.\cr
#' Then follows one line for each heritable trait, in the order they are defined in the model
#' (emigration / movement / settlement).
#' They begin with the keyword \code{Trait} followed by an consecutive (starting with zero) integer
#' trait ID. On the same line follows the keyword \code{NLoci} and an integer \eqn{n_Ti} that
#' specifies the number of loci the corresponding trait is mapped to. After this follow \eqn{n_Ti}
#' integer pairs specifying the respective chromosome and locus. (Note that numbering starts at \eqn{0}.)\cr
#'
#' Any loci which do not contribute to at least one trait are treated as neutral loci (see below).
#' Neutral loci may also be specified for a model in which there are no adaptive traits:
#' Specifying a genetic architecture file for a model having no adaptive traits will result in neutral
#' genetics being set up for the species.\cr
#'
#' \emph{From genotype to phenotype}\cr
#' All alleles are represented by integer values (positive or negative),
#' and the sum of alleles at all loci contributing to a trait (both alleles at each locus of a diploid
#' species) controls the phenotype. However, as phenotypic traits exist on several widely different
#' scales (e.g. emigration probability between \eqn{0} and \eqn{1}, dispersal kernel mean typically
#' many hundreds or thousands of metres), it is necessary to specify how the allelic scale relates
#' to the phenotypic scale. A scaling factor is specified for each trait (the parameter
#' \code{TraitScaleFactor} in each dispersal sub-module), which governs how large a
#' change of \eqn{100} units (which is the ‘integer base’ of the genome) on the allele scale will be on the
#' phenotypic scale. For example, if the scaling factor for density-independent emigration probability
#' is \eqn{0.1}, and a juvenile’s sum of all alleles contributing to that trait is \eqn{150} less than
#' its parent’s equivalent sum, then the juvenile’s emigration probability phenotype will be \eqn{0.15}
#' lower than its parent’s phenotype (but subject in this case to the constraint that it may not be
#' less than zero).\cr
#'
#' \emph{Mutation and Linkage}\cr
#' Mutation can occur when an individual allele is inherited by a new-born individual from its parent.
#' It is governed by two genome-level parameters: the mutation probability \code{ProbMutn}
#' that a mutation takes place, and standard deviation \code{MutationSD} which sets the typical magnitude.
#' Both are applied in a standard way at the level of the individual locus (unlike in version 1, in
#' which separate probabilities and magnitudes were applied for separate traits). When a mutation occurs
#' at a locus, a random number is drawn from a normal distribution having zero mean and the mutation standard
#' deviation. This number is multiplied by the integer base to yield an integer value which is added to the
#' allele before it is copied to the juvenile’s chromosome. The default mutation standard deviation of \eqn{0.1}
#' will therefore give mutations which mostly range between \eqn{-30} and \eqn{+30} at the allele scale.\cr
#'
#' In a diploid species, traits may be linked if their coding loci are on the same chromosomes; or they can
#' evolve independently, if their coding loci are on different chromosomes.
#' The degree of linkage, however, also depends on the crossover probability \code{ProbCross} specified for
#' the genome. If it is high, the degree of linkage is reduced, and the linked traits tend to evolve more
#' independently. \cr
#'
#' The crossover probability denoted the probability that,
#' when copying a multi-locus chromosome from a parent’s to the offspring’s genome, there will
#' be a change at the current locus from copying the parent’s maternally-inherited alleles to the
#' parent’s paternally-inherited alleles or vice versa.\cr
#' The crossover probability parameter is applied at the scale of the locus, i.e. during meiosis.
#' As a parent’s chromosomes are being inherited by its offspring, a crossover occurs at each locus with the
#' specified probability \code{ProbCross}. If the crossover probability is high, the degree of linkage is reduced, and two linked
#' traits tend to evolve more independently.\cr
#'
#' \emph{Pleiotropy and Neutral loci}\cr
#' It is possible that a particular locus can be specified more than once for a particular trait; this
#' increases its weighting relative to other loci contributing to the trait. A locus may also code for
#' more than one trait, i.e. it is pleiotropic. Large positive allele values at that locus will tend to
#' lead to larger than average values of both traits (but subject to any other loci coding separately
#' for the two traits). Thus the two traits are forced (to some extent) to be positively correlated. It
#' is not possible to specify inverse correlation between two traits in this way.\cr
#'
#' A locus which does not code for any trait is neutral, and therefore not subject directly to
#' selection, although the distribution across the population of allele values at that locus may vary
#' over time owing to genetic drift and/or linkage to loci which are under selection. Such neutral
#' loci may be used as markers for estimating relatedness at the individual or population level.\cr
#'
#' \emph{Genome initialisation}\cr
#' At model initialisation, individuals are assigned phenotypic traits drawn from a normal
#' distribution controlled by a specified mean phenotype and standard deviation, but also subject
#' to any phenotypic constraints (e.g. a probability must lie between \eqn{0} and \eqn{1}, step length
#' must be positive, etc.). The values for the traits initial mean and standard deviation are set in the appropriate
#' sub-modules, where they replace the constant values used when \code{IndVar = FALSE}.
#' The standard deviation for a trait may not be greater than the corresponding
#' scaling factor (\code{TraitScaleFactor}), but it may be substantially less if an initial population which is highly homogeneous
#' for a particular trait is required. \cr
#'
#' The genome actually controls individual variation relative to the
#' initial population mean value of a trait; thus if the sum of an individual’s alleles is negative, its
#' phenotype will be less than the initial mean value, and if its sum is positive, its phenotype will be
#' greater than the mean value. However, to prevent all alleles for a multi-loci trait being identical in
#' initial individuals, further random variation is applied at the allele scale (i.e. common to all traits),
#' which is drawn from a normal distribution having zero mean and a standard deviation specified in \code{AlleleSD}.
#' Note that therefore the observed variance in a trait value across the initial population may not match
#' exactly the specified variance for the trait.
#'
#' Side note: This differs from the implementation in version 1.x, which used a uniform distribution; hence
#' an initial population cannot be set up in v2.0 to have exactly the same properties as in v1, but a similar
#' equilibrium population should arise after a period of sufficiently strong selection for one or more traits.
#'
#'
# #' @references \insertAllCited{}
#' @return a parameter object of class "GeneticsParams"
#' @author Anne-Kathleen Malchow
#' @name Genetics
#' @export Genetics
Genetics <- setClass("GeneticsParams", slots = c(GenomeSize = "integer_OR_numeric",
                                                 ChromosomeEnds = "integer_OR_numeric", # NULL or vector
                                                 RecombinationRate = "integer_OR_numeric", # NULL or numeric
                                                 OutputGeneValues = "logical",
                                                 OutputNeutralStatistics = "logical",
                                                 OutputFstatsWeirCockerham = "logical",
                                                 OutputFstatsWeirHill = "logical",
                                                 OutputStartGenetics = "integer_OR_numeric", # positive integer if any output is TRUE or NULL
                                                 OutputInterval = "integer_OR_numeric",
                                                 PatchList = "integer_OR_numeric", # vector or string
                                                 NbrPatchToSample = "integer_OR_numeric", # NULL or integer
                                                 nIndividualsToSample = "character", # NULL or integer
                                                 Stages = "integer_OR_numeric", # vector
                                                 Traits = "TraitsParams")
                                     , prototype = list(GenomeSize = 0L,
                                                        ChromosomeEnds = 0L, # NULL or vector
                                                        RecombinationRate = NULL, # NULL or numeric
                                                        OutputGeneValues = FALSE,
                                                        OutputNeutralStatistics = FALSE,
                                                        OutputFstatsWeirCockerham = FALSE,
                                                        OutputFstatsWeirHill = FALSE,
                                                        OutputStartGenetics = NULL, # positive integer if any output is TRUE or NULL
                                                        OutputInterval = NULL,
                                                        PatchList = NULL, # vector or string
                                                        NbrPatchToSample = NULL, # NULL or integer
                                                        nIndividualsToSample = NULL, # NULL or integer
                                                        Stages = NULL, # vector
                                                        Traits = Traits())

)
setValidity('GeneticsParams', function(object){
    # msg <- NULL
    # if(anyNA(object@Architecture) || length(object@Architecture)!=1) {
    #     msg <- c(msg, "Architecture must be set!")
    # }
    # else {
    #     if(object@Architecture %in% c(0,1)) {
    #         if(object@Architecture == 0) {
    #             if(anyNA(object@NLoci) || length(object@NLoci)!=1) {
    #                 msg <- c(msg, "NLoci must be set if Architecture = 0 (one chromosome per trait)!")
    #             }
    #             else {
    #                 if(object@NLoci < 1) {
    #                     msg <- c(msg, "NLoci must be greater than 0!")
    #                 }
    #             }
    #         }
    #         if(object@Architecture == 1) {
    #             if (object@ArchFile == "NULL"){
    #                 msg <- c(msg, 'ArchFile is required if Architecture = 1 (load architecture file).')
    #             }
    #         }
    #     }else{
    #         msg <- c(msg, "Architecture must be either 0 or 1!")
    #     }
    # }
    # if(anyNA(object@ProbMutn) || length(object@ProbMutn)!=1) {
    #     msg <- c(msg, "ProbMutn must be set!")
    # }
    # else {
    #     if(object@ProbMutn < 0.0 || object@ProbMutn > 1.0) {
    #         msg <- c(msg, "ProbMutn must be within the closed interval [0,1]!")
    #     }
    # }
    # if(anyNA(object@ProbCross) || length(object@ProbCross)!=1) {
    #     msg <- c(msg, "ProbCross must be set!")
    # }
    # else {
    #     if(object@ProbCross < 0.0 || object@ProbCross > 1.0) {
    #         msg <- c(msg, "ProbCross must be within the closed interval [0,1]!")
    #     }
    # }
    # if(anyNA(object@AlleleSD) || length(object@AlleleSD)!=1) {
    #     msg <- c(msg, "AlleleSD must be set!")
    # }
    # else {
    #     if(object@AlleleSD <= 0.0) {
    #         msg <- c(msg, "AlleleSD must be strictly positive!")
    #     }
    # }
    # if(anyNA(object@MutationSD) || length(object@MutationSD)!=1) {
    #     msg <- c(msg, "MutationSD must be set!")
    # }
    # else {
    #     if(object@MutationSD <= 0.0) {
    #         msg <- c(msg, "MutationSD must be strictly positive!")
    #     }
    # }
    # if (is.null(msg)) TRUE else msg
    }
)
setMethod('initialize', 'GeneticsParams', function(.Object, ...) {
    # this_func = "Genetics(): "
    # args <- list(...)
    # .Object <- callNextMethod()
    # if(.Object@Architecture == 0) {
    #     if (!is.null(args$ArchFile)) {
    #         warning(this_func, "ArchFile", warn_msg_ignored, "since Architecture = 0 (one chromosome per trait).", call. = FALSE)
    #     }
    # }
    # if(.Object@Architecture == 1) {
    #     if (!is.null(args$NLoci)) {
    #         warning(this_func, "NLoci", warn_msg_ignored, "since Architecture = 1 (load architecture file).", call. = FALSE)
    #     }
    # }
    # .Object
    }
)
setMethod("show", "GeneticsParams", function(object){
    # cat(" Genetics: \n")
    # cat("   Architecture =", object@Architecture, ": ")
    # if (object@Architecture == 0) {
    #     cat("One chromosome per trait \n")
    #     cat("   with", object@NLoci)
    #     if (object@NLoci==1) cat(" locus")
    #     else cat(" loci")
    #     cat(" per chromosome.")
    # }
    # if (object@Architecture == 1) {
    #     cat("loaded from architecture file:\n")
    #     cat("   ", object@ArchFile)
    # }
    # cat("\n")
    # cat("   ProbMutn   =", object@ProbMutn, "\n")
    # cat("   ProbCross  =", object@ProbCross, "\n")
    # cat("   AlleleSD   =", object@AlleleSD, "\n")
    # cat("   MutationSD =", object@MutationSD, "\n")
})
