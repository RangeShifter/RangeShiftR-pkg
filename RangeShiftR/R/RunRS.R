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
            return(raster::stack(lapply(X = out, FUN = raster::raster, xmn=0, xmx=ncol(out[[1]])*resol, ymn=0, ymx=nrow(out[[1]])*resol)))
        }
        else return(NULL)
    }
    else {
        return(out)
    }
}
