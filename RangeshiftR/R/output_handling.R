
#-------------------
# Output handling and plotting functions
#-------------------


#---------------------------------------------------------

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


#---------------------------------------------------------

### PROCESS OUTPUT


#' ColonisationStats
#'
#' This function produces patch statistics and maps on occupancy probability and mean time to colonisation.
#'
#' It uses the Rangeshifter 'population' output data and, as of now, is only implemented for patch-based models.
#' @usage ColonisationStats(x, y = NULL, years = numeric(0), maps = FALSE)
#' ColonisationStats(x, y = NULL, years = numeric(0))
#' @param x,y Either the parameter master (\code{x}) and RS directory path (\code{y}) of the simulation or, alternatively,\cr
#' the population output as a dataframe (\code{x}). In this case \code{y} is an optional parameter taking the patch map(s) as a raster layer (stack),
#' which will then be used to produce maps of the patch statistics. For more info on this, see the Details.
#' @param years Years at which the probabilty of occupancy over all replicates will be calculated.
#' @param maps Only in the parameter master (=\code{x}) notation: For each given \code{year}, uses the current patch map (in case of a dynamic landscape, otherwise the static patch map)
#' to produce a raster of occupancy probabilties. For time to colonisation, uses the current patch map of last year of population record.
#' @details In the population dataframe (=\code{x}) notation, there are 3 options on how many maps to produce:
#' (1) y = NULL: no raster maps produced.
#' (2) y is a RasterLayer: All statistics projected onto the same patch map.
#' (3) y is a RasterStack with 2 layers: The first will be used for probabilty of occupancy for all years, the second will be used for time to colonisation.
#' (4) y is a RasterStack with length(years)+1 layers: The first ones will be used for probabilty of occupancy on each year, the last for time to colonisation.
#' @return a list with the elements
#' \code{occ_prob} (If one year is given: a named numeric with probabilty of occupancy at the given year over all replicates.
#' If a vector of years is given: a numeric array with rows named by years and columns by patch-IDs.),
#' \code{col_time} (a dataframe with first recorded year of colonisation for each patch; with replicates in rows and columns named by patch-IDs),
#' Optional:\cr
#' \code{patch_occ_prob} (a raster (stack) with the data from \code{occ_prob} stored in the respective patch cells per given year),
#' \code{patch_col_time} (a raster with the data from \code{col_time_mean} stored in their respective patch cells)
#' @export
setGeneric("ColonisationStats", function(x, ...) standardGeneric("ColonisationStats") )

setMethod("ColonisationStats", "data.frame", function(x, y = NULL, years = numeric(0)) { # y="BasicRaster"
#test <- function(x, y = NULL, years = numeric(0)) {
    if("PatchID" %in% names(x)){    # include cell-based models?

        if("PatchID" %in% names(x)) x <- x[x$PatchID!=0,]   # added this, ok?

        if(class(years) %in% c("integer","numeric") ){
            # check that requested years are in dataset
            if(!all(years %in% unique(x$Year))){
                years <- sort(years[years %in% unique(x$Year)])
                warning("ColonisationStats(): Some given years are not included in the population dataframe.", call. = FALSE)
            }
            if(length(years)==0) years==max(x$Year)
        }
        else{
            warning("ColonisationStats(): Years must be of class numeric or integer.", call. = FALSE)
            return(NULL)
        }

        patches <- sort(unique(x$PatchID))
        # Occupancy probability in a given year over all replicates
        occ_prob <- sapply(patches,
                           FUN=function(patch,n=length(unique(x$Rep))){
                               sapply(years,
                                      FUN=function(year){
                                          sum(subset(x,PatchID==patch & Year==year)$NInd>0)/n
                                          }
                                   )
                               }
                           )
        if(length(years)>1){
            colnames(occ_prob) <- patches
            rownames(occ_prob) <- years
        }else  names(occ_prob) <- patches

        # Time to colonisation:
        col_time <- data.frame(sapply(patches,
                                      FUN=function(p){sapply(sort(unique(x$Rep)),
                                                             FUN=function(r){
                                                                 ifelse(nrow(subset(x,PatchID==p & Rep==r & NInd>0)),
                                                                        min(subset(x,PatchID==p & Rep==r & NInd>0)$Year),
                                                                        NA)
                                                                 }
                                                             )
                                          }
                                      )
                               )
        names(col_time) <- patches
        col_time_mean <- colMeans(col_time, na.rm = TRUE)

        # create maps
        if(!is.null(y)){
            require('raster')

            # non-dynamic landscape
            if(class(y) == "RasterLayer" || (class(y) == "RasterStack" && length(y@layers)==1) ){

                # initialise output rasters
                if(class(y) == "RasterStack") y <- y[[1]]
                patch_occ_prob <- patch_col_time <- y
                # denote matrix with NA
                values(patch_occ_prob)[values(y)==0] <- values(patch_col_time)[values(y)==0] <- NA
                # all habitat patches to address those that never had a population
                values(patch_occ_prob)[values(y) >0] <- 0
                values(patch_col_time)[values(y) >0] <- -9

                # fill output rasters
                if(length(years)>1){
                    patch_outstack <- stack()
                    for (j in 1:length(years)){
                        patch_outstack <- addLayer(patch_outstack, patch_occ_prob)
                        for (i in patches){
                            values(patch_outstack[[j]])[values(y)==i] <- occ_prob[j,paste(i)]
                        }
                    }
                }
                else{
                    for (i in patches){
                        values(patch_occ_prob)[values(y)==i] <- occ_prob[paste(i)]
                    }
                    patch_outstack <- patch_occ_prob
                }

                for (i in patches){
                    values(patch_col_time)[values(y)==i] <- ifelse(is.na(col_time_mean[paste(i)]),-9,col_time_mean[paste(i)])
                }

                return(list(occ_prob=occ_prob, col_time=col_time, map_occ_prob=patch_outstack, map_col_time=patch_col_time))
            }

            # dynamic landscape
            if(class(y) == "RasterStack" && length(y@layers)>1 ){
                N_layers <- length(years)+1
                if( length(y@layers) != N_layers ){
                    warning("ColonisationStats(): Number of raster layers must be either 1 or number of years plus 1.", call. = FALSE)
                }
                else{
                    # initialise output rasters
                    patch_outstack <- y
                    # denote matrix with NA
                    values(patch_outstack)[values(y)==0] <- NA
                    # all habitat patches to address those that never had a population
                    values(patch_outstack)[values(y)>0] <- 0.0
                    values(patch_outstack[[N_layers]])[values(y[[N_layers]])>0] <- -9

                    # fill output rasters
                    for (i in patches){
                        if(N_layers==2){
                            values(patch_outstack[[1]])[values(y[[1]])==i] <- occ_prob[paste(i)]
                        }else {
                            for (j in 1:length(years)){
                                values(patch_outstack[[j]])[values(y[[j]])==i] <- occ_prob[j,paste(i)]
                            }
                        }
                        values(patch_outstack[[N_layers]])[values(y[[N_layers]])==i] <- ifelse(is.na(col_time_mean[paste(i)]),-9,col_time_mean[paste(i)])
                    }
                    return(list(occ_prob=occ_prob, col_time=col_time, map_occ_prob=patch_outstack[[-N_layers]], map_col_time=patch_outstack[[N_layers]]))
                }
            }
        }
        return(list(occ_prob=occ_prob, col_time=col_time))

    }else{warning("ColonisationStats(): This function is implemented for patch-based models only.", call. = FALSE)}
})

setMethod("ColonisationStats", "RSparams", function(x, y = NULL, years = numeric(0), maps = FALSE) {
    if(s@control@patchmodel){
        if(s@simul@OutIntPop>0){
            if(!is.null(y) & class(y)=="character" ){
                if(class(years) %in% c("integer","numeric") ){

                    # read population output
                    pop_df <- try(readPop(x, y))
                    if ( class(pop_df) == "try-error" ) warning("ColonisationStats(): Couldn't read population output for this simulation.", call. = FALSE)
                    if(length(years)==0) years <- max(pop_df$Year)

                    # read patch rasters if needed
                    if(maps){
                        require('raster')
                        #non-dynamic landscape
                        if(length(s@land@LandscapeFile)==1){
                            patch_r <- try(raster(paste0(dirpath, "Inputs/", s@land@PatchFile)))
                            if ( class(patch_r) == "try-error" ) warning("ColonisationStats(): Couldn't read patch raster file nr ", current , " for this simulation.", call. = FALSE)
                        }
                        #dynamic landscape
                        else{
                            patch_r <- stack()
                            # rasters for occ_prob output
                            for(year in years){
                                current <- which(s@land@DynamicLandYears == max(s@land@DynamicLandYears[s@land@DynamicLandYears<=year]) )
                                patch_curr <- try(raster(paste0(dirpath, "Inputs/", s@land@PatchFile[current])))
                                if ( class(patch_curr) == "try-error" ) warning("ColonisationStats(): Couldn't read patch raster file nr ", current , " for this simulation.", call. = FALSE)
                                else patch_r <- addLayer(patch_r ,patch_curr)
                            }
                            # rasters for col_time output
                            year <- max(pop_df$Year)
                            current <- which(s@land@DynamicLandYears == max(s@land@DynamicLandYears[s@land@DynamicLandYears<=year]) )
                            patch_curr <- try(raster(paste0(dirpath, "Inputs/", s@land@PatchFile[current])))
                            patch_r <- addLayer(patch_r ,patch_curr)
                        }
                        if(class(pop_df) != "try-error" & length(patch_r@layers)==(length(years)+1) ) ColonisationStats(pop_df,patch_r,years)
                    }else {
                        if(class(pop_df) != "try-error") ColonisationStats(pop_df,NULL,years)
                    }
                }else {
                    warning("ColonisationStats(): Years must be of class numeric or integer.", call. = FALSE)}
            }else {
                warning("ColonisationStats(): Dirpath must be set", call. = FALSE)}
        }else{
            warning("ColonisationStats(): This simulation has population output turned off (OutIntPop=0), but that is needed to calculate the colonisation statistics.", call. = FALSE)}
    }else {
        warning("ColonisationStats(): This function is implemented for patch-based models only.", call. = FALSE)}
})





#if(length(s@land@LandscapeFile)==1 && all(s@land@DynamicLandYears==0) ){
#}else {
#    warning("ColonisationStats(): This function is implemented for non-dynamic landscape models only. However, you can specify a patch raster and use the method for signature c(data.frame,RasterLayer),", call. = FALSE)}


#---------------------------------------------------------

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
    if (class(dirpath)=="character"){
        range_table <- try(readRange(s, dirpath))
        if ( class(range_table) == "data.frame") {
            plotAbundance(range_table, ...)
        }
    }else{
        warning("plotAbundance(): dirpath must be of type character.", call. = TRUE)
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
    lines(rep_means$Year, rep_means$NOccup, type = "l", lwd = 3, col = "blue")
    # plot replicates
    if (replicates) {
        for (i in 0:max(s$Rep)) {
            lines(s$Year[s$Rep==i], s$NOccup[s$Rep==i], type = "l", lwd = 0.5)
        }
    }
})
setMethod("plotOccupancy", "RSparams", function(s, dirpath, ...) {
    if (class(dirpath)=="character"){
        range_table <- try(readRange(s, dirpath))
        if ( class(range_table) == "data.frame") {
            plotOccupancy(range_table, ...)
        }
    }else{
        warning("plotOccupancy(): dirpath must be of type character.", call. = TRUE)
    }
})
