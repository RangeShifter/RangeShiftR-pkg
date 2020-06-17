#' @author Anne-Kathleen Malchow, \email{malchowa@geo.hu-berlin.de}
#' @details The only function you're likely to need from roxygen2 is [roxygenize()].
#' Otherwise refer to the vignettes to see how to format the documentation.
#'
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


#' @useDynLib RangeshiftR
#' @importFrom Rcpp sourceCpp
NULL


# #' @references \url{https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12162%4010.1111/%28ISSN%292041-210X.Biogeographys}




##------- TODO

# add citation()

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
#'   \item \code{\link[RangeshiftR]{Emigration}}: plot emigration probability
#'   \item \code{\link[RangeshiftR]{DispersalKernel}}:
#'   \itemize{
#'     \item \code{mortality=FALSE} - plot dispersal distance probability density  (default)
#'     \item \code{mortality= TRUE} - plot mortality probability
#'   }
#'   If a mixed kernel was defined (i.e. \code{DoubleKernel=TRUE}), plot the resulting dispersal probability by...
#'   \itemize{
#'     \item \code{combinekernels=FALSE} - ...plotting both kernels separately (default)
#'     \item \code{combinekernels= TRUE} - ...combining both kernels, i.e. \ifelse{html}{ \out{p(d; &delta;<sub>1</sub>,&delta;<sub>2</sub>) = p<sub>I</sub> p(d;&delta;<sub>1</sub>) + (1-p<sub>I</sub>) p(d;&delta;<sub>1</sub>) } }{\deqn{ p(d; δ_1,δ_2) = p_I p(d;δ_1) + (1-p_I) p(d;δ_2)} }
#'   }
#'   \item \code{\link[RangeshiftR]{StageStructure}}: plot fecundity as well as survival and development probabilities
#' }
#' @export
setGeneric("plotProbs", function(x, ...) standardGeneric("plotProbs") )

## density dependence function
densdep <- function(x, A0 = 1.0, alpha = 1.0, beta = 0.0) {
    A0/(1+exp(alpha*(beta-x)))
}

