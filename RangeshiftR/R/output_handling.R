
#-------------------
# Output handling and plotting functions
#-------------------


### READING OUTPUT FILES


#' Read 'range' file
#'
#' Read the Rangeshifter output file 'range' into a data.frame, if it was generated.
#' @param s RSmaster parameter object
#' @param dirpath RS directory path
#' @return a data.frame
#' @export
setGeneric("readRange", function(s,dirpath) standardGeneric("readRange") )

setMethod("readRange", c(s="RSparams", dirpath="character"), function(s,dirpath) {
    path <- paste0(dirpath, "Outputs/Batch", s@control@batchnum, "_Sim", s@simul@Simulation, "_Land", s@land@LandNum, "_Range.txt")
    if(file.exists(path)){
        return(read.table(path, h = T, sep = "\t"))
    }else {
        warning("The 'range' output file for this simulation does not exist.", call. = FALSE)
    }
})


#' Read 'pop' file
#'
#' Read the Rangeshifter output file 'pop' into a data.frame, if it was generated.
#' @param s RSmaster parameter object
#' @param dirpath RS directory path
#' @details The x- and y- coordinates (in a cell-based model) differ from those in the text file in that they will be shifted from the
#' lower left corner to the center of their respective cell.
#' @return a data.frame
#' @export
setGeneric("readPop", function(s,dirpath) standardGeneric("readPop") )

setMethod("readPop", c(s="RSparams", dirpath="character"), function(s,dirpath) {
    path <- paste0(dirpath, "Outputs/Batch", s@control@batchnum, "_Sim", s@simul@Simulation, "_Land", s@land@LandNum, "_Pop.txt")
    if(file.exists(path)){
        pop_table <- read.table(path, h = T, sep = "\t")
        if (sum(grepl('x',names(pop_table)))){
            pop_table$x <- (pop_table$x+0.5)*s@land@Resolution
            pop_table$y <- (pop_table$y+0.5)*s@land@Resolution
        }
        return(pop_table)
    }else {
        warning("The 'pop' output file for this simulation does not exist.", call. = FALSE)
    }
})



### PLOTTING


#' Plot Abundance
#'
#' Uses the Rangeshifter output data 'range' to generate abundance time series.
#' Plots the mean abundance over all replicates, and optionally the standard deviation and/or the single replicates.
#' @param s RSparams object or a data.frame in the 'range' file format
#' @param dirpath RS directory path; required if \code{s} is a \code{RSparams}
#' @param sd plot standard deviation? (default is \code{FALSE})
#' @param replicates plot the replicates? (default is \code{TRUE})
#' @param ylim upper limit to the y-axis
#' @export
setGeneric("plotAbundance", function(s,...) standardGeneric("plotAbundance") )

setMethod("plotAbundance", "data.frame", function(s, sd = FALSE, replicates = TRUE, ylim=NULL, ...) {
    # Calculate means
    rep_means <- aggregate(NInds~Year, data = s, FUN = "mean")
    # Calculate standard deviation
    if (sd) {
        rep_sd <- aggregate(NInds~Year, data = s, FUN = "sd")
        rep_sd[is.na(rep_sd)] <- 0
    }
    # Set y-limits
    if (is.null(ylim)) {
        if (replicates) {
            ylim <- c(min(s$NInds),max(s$NInds))
        }else {
            if (sd) {
                ylim <- c(max(0,min(rep_means$NInds-rep_sd$NInds)), max(rep_means$NInds+rep_sd$NInds))
            }else{
                ylim <- c(min(rep_means$NInds),max(rep_means$NInds))
            }
        }
    }
    # New plot
    plot(NULL, type = "n", ylab = "Abundance", xlab = "Year", xlim=c(0, max(s$Year)), ylim=ylim, ...)
    # Plot standard deviation
    if (sd) {
        polygon(c(rep_sd$Year,rev(rep_sd$Year)), c(rep_means$NInds+rep_sd$NInds, rev(pmax(0,rep_means$NInds-rep_sd$NInds))), border=NA, col='grey80')
    }
    # Plot abundance
    lines(rep_means$Year, rep_means$NInds, type = "l", lwd = 3, col = "red")
    # Plot replicates
    if (replicates) {
        for (i in 0:max(s$Rep)) {
            lines(s$Year[s$Rep==i], s$NInds[s$Rep==i], type = "l", lwd = 0.5)
        }
    }
})
setMethod("plotAbundance", "RSparams", function(s, dirpath, ...) {
    range_table <- try(readRange(s, dirpath))
    if ( class(range_table) == "data.frame") {
        plotAbundance(range_table, ...)
    }
})


#' Plot Occupancy
#'
#' Uses the Rangeshifter output data 'range' to generate occupancy time series.
#' Plots the mean occupancy over all replicates, and optionally the standard deviation and/or the single replicates.
#' @param s RSparams object or a data.frame in the 'range' file format
#' @param dirpath RS directory path; required if \code{s} is a \code{RSparams}
#' @param sd plot standard deviation? (default is \code{FALSE})
#' @param replicates plot the replicates? (default is \code{TRUE})
#' @param ylim upper limit to the y-axis
#' @export
setGeneric("plotOccupancy", function(s,...) standardGeneric("plotOccupancy") )

setMethod("plotOccupancy", "data.frame", function(s, sd = FALSE, replicates = TRUE, ylim=NULL, ...) {
    names(s)[grep('NOccup',names(s))] <- 'NOccup'
    # Calculate means
    rep_means <- aggregate(NOccup~Year, data = s, FUN = "mean")
    # Calculate standard deviation
    if (sd) {
        rep_sd <- aggregate(NOccup~Year, data = s, FUN = "sd")
        rep_sd[is.na(rep_sd)] <- 0
    }
    # Set y-limits
    if (is.null(ylim)) {
        if (replicates) {
            ylim <- c(min(s$NOccup),max(s$NOccup))
        } else {
            if (sd) {
                ylim <- c(max(0,min(rep_means$NOccup-rep_sd$NOccup)), max(rep_means$NOccup+rep_sd$NOccup))
            } else {
                ylim <- c(min(rep_means$NOccup),max(rep_means$NOccup))
            }
        }
    }
    # New plot
    plot(NULL, type = "n", ylab = "Occupancy", xlab = "Year", xlim=c(0, max(s$Year)), ylim=ylim, ...)
    # Plot standard deviation
    if (sd) {
        polygon(c(rep_sd$Year,rev(rep_sd$Year)), c(rep_means$NOccup+rep_sd$NOccup, rev(pmax(0,rep_means$NOccup-rep_sd$NOccup))), border=NA, col='grey80')
    }
    # Plot occupancy
    lines(rep_means$Year, rep_means$NOccup, type = "l", lwd = 3, col = "red")
    # plot replicates
    if (replicates) {
        for (i in 0:max(s$Rep)) {
            lines(s$Year[s$Rep==i], s$NOccup[s$Rep==i], type = "l", lwd = 0.5)
        }
    }
})
setMethod("plotOccupancy", "RSparams", function(s, dirpath, ...) {
    range_table <- try(readRange(s, dirpath))
    if ( class(range_table) == "data.frame") {
        plotOccupancy(range_table, ...)
    }
})
