#' @author Anne-Kathleen Malchow, \email{malchowa@geo.hu-berlin.de}
#' @details
#' Bocedi, Greta, et al. "RangeShifter: a platform for modelling spatial eco‐evolutionary dynamics and species' responses to environmental changes." Methods in Ecology and Evolution 5.4 (2014): 388-396.
#'
#' @section Warning:
#' You must not call this function unless ...
#'
#' \subsection{Exceptions}{
#'    Apart from the following special cases...
#' }
#' @references \url{https://en.wikipedia.org/wiki/List_of_Crayola_crayon_colors}
#' @keywords internal
"_PACKAGE"


#' @useDynLib RangeShiftR
#' @importFrom Rcpp sourceCpp
NULL


# #' @references \url{https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12162%4010.1111/%28ISSN%292041-210X.Biogeographys}


.onAttach <- function(libname, pkgname) {
    packageStartupMessage("\nRangeshiftR version 1.0.0 (01.10.2020)\n",
                          "Copyright (C) 2020 Greta Bocedi, Justin Travis, Steve Palmer, Anne-Kathleen Malchow, Damaris Zurell\n\n",
                          "This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.\n",
                          "You are welcome to redistribute it and/or modify it under certain conditions;\n",
                          "type 'RangeShiftR_license()' for details.\n")
}

#' Display RangeShiftR license information
#'
#' @export
RangeShiftR_license <- function ()
{
    cat("\nRangeshiftR version 1.0.0 (01.10.2020)\n")
    cat("Copyright (C) 2020 Greta Bocedi, Justin Travis, Steve Palmer, Anne-Kathleen Malchow, Damaris Zurell\n\n")
    cat("This program is free software: you can redistribute it and/or modify\n")
    cat("it under the terms of the GNU General Public License as published by\n")
    cat("the Free Software Foundation, either version 3 of the License, or\n")
    cat("(at your option) any later version.\n\n")
    cat("This software is distributed in the hope that it will be useful,\n")
    cat("but WITHOUT ANY WARRANTY; without even the implied warranty of\n")
    cat("MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the\n")
    cat("GNU General Public License for more details.\n\n")
    cat("You should have received a copy of the GNU General Public License\n")
    cat("along with this software. Copies of the license can also be found at\n")
    cat("http://www.gnu.org/licenses/\n")
    cat("\n")
    cat("'Share and Enjoy.'\n\n")
}

#'\ Display citation
#'
#' @export
RangeShiftR_citation <- function ()
{
    citation(package = "RangeShiftR", lib.loc = NULL, auto = NULL)
}


##------- TODO

# complete package RangeShiftR.R , meta data

# write setters and getters

# think about hierarchy of sub-classes and "+"-methods

# add references to documentation

# new class EnvGradient to move these params into from SimulationParams? Maybe as sub-class to Land-Params?

# write new output (plotting) functions

# clean C++ code to address warnings




##------ define ClassUnions

#for Integer and Numerical so that 'Integer'-slots can also accept 'Numerical' input
setClassUnion("integer_OR_numeric", c("integer", "numeric"))

#for Matrix and Numerical so that 'Matrix'-slots can also accept 'Numerical' input when 1x1-Matrix is expected
setClassUnion("matrix_OR_numeric", c("matrix", "numeric"))

#for cost layer to accept habitat codes or raster map file
setClassUnion("numeric_OR_character", c("numeric", "character"))


## define error and warning messages
warn_msg_ignored = " will be ignored "


# -----
# Plotting parameterised probabilities functions
# -----

#' Plot parameterised probabilities
#'
#' Visualise the dependencies or distributions of probabilities of the various processes that are defined by the different parameter modules,
#' like e.g. fecundity or mortality.
#'
#' @param x a RangeShiftR parameter object
#' @param xmax,ymax upper plotting limits (lower limits are fixed to \eqn{0})
#' @param stage,sex stage(s) and sexe(s) to plot, default: all
#' @param ... various options, depend on the given parameter module type
#' @details
#' Available methods and their options:
#' \itemize{
#'   \item \code{\link[RangeShiftR]{Emigration}}: plot emigration probability
#'   \item \code{\link[RangeShiftR]{DispersalKernel}}:
#'   \itemize{
#'     \item \code{mortality=FALSE} - plot dispersal distance probability density  (default)
#'     \item \code{mortality= TRUE} - plot mortality probability
#'   }
#'   If a mixed kernel was defined (i.e. \code{DoubleKernel=TRUE}), plot the resulting dispersal probability by...
#'   \itemize{
#'     \item \code{combinekernels=FALSE} - ...plotting both kernels separately (default)
#'     \item \code{combinekernels= TRUE} - ...combining both kernels, i.e. \ifelse{html}{ \out{p(d; &delta;<sub>1</sub>,&delta;<sub>2</sub>) = p<sub>I</sub> p(d;&delta;<sub>1</sub>) + (1-p<sub>I</sub>) p(d;&delta;<sub>1</sub>) } }{\deqn{ p(d; δ_1,δ_2) = p_I p(d;δ_1) + (1-p_I) p(d;δ_2)} }
#'   }
#'   \item \code{\link[RangeShiftR]{StageStructure}}: plot fecundity as well as survival and development probabilities
#' }
#' @export
setGeneric("plotProbs", function(x, ...) standardGeneric("plotProbs") )

## density dependence function
densdep <- function(x, A0 = 1.0, alpha = 1.0, beta = 0.0) {
    A0/(1+exp(alpha*(beta-x)))
}

