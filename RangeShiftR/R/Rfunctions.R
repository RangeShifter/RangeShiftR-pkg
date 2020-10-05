
# -----
#
# R-level functions
#
# -----


#-- define ClassUnions
#----------------------

#for Integer and Numerical so that 'Integer'-slots can also accept 'Numerical' input
setClassUnion("integer_OR_numeric", c("integer", "numeric"))

#for Matrix and Numerical so that 'Matrix'-slots can also accept 'Numerical' input when 1x1-Matrix is expected
setClassUnion("matrix_OR_numeric", c("matrix", "numeric"))

#for cost layer to accept habitat codes or raster map file
setClassUnion("numeric_OR_character", c("numeric", "character"))


#-- define error and warning messages
#----------------------
warn_msg_ignored = " will be ignored "


#-- density dependence function
#----------------------
densdep <- function(x, A0 = 1.0, alpha = 1.0, beta = 0.0) {
    A0/(1+exp(alpha*(beta-x)))
}


#-- validation function
#----------------------

#' Validate a given RS parameter object
#'
#' @export
validateRSparams <- function(x){validObject(x)}
