
#-------------------
# R-level functions
#-------------------


#' Define a RangeShiftR parameter master object
#'
#' Set up a parameter master that can be used in \code{\link[RangeshiftR]{RunRS}}() to run a simulation.\cr
#' All parameter modules can be added to an existing parameter master via the "+"-functions. However, note that the entire respective module will be overwritten.\cr
#'
#' @usage RSsim(batchnum = 1L,
#'       simul = Simulation(),
#'       land = ArtificialLandscape(),
#'       demog = Demography(Rmax = 1.5),
#'       dispersal = Dispersal(),
#'       gene = Genetics(),
#'       init = Initialise())
#' @include class_RSparams.R
#' @param batchnum Batch ID is part of output files names and can be used to prevent overwriting.
#' @param simul Set \code{\link[RangeshiftR]{Simulation}} parameters
#' @param land Set landscape parameters. Can be either \code{\link[RangeshiftR]{ArtificialLandscape}} or \code{\link[RangeshiftR]{ImportedLandscape}}.
#' @param demog Set \code{\link[RangeshiftR]{Demography}} parameters
#' @param dispersal Set \code{\link[RangeshiftR]{Dispersal}} parameters
#' @param gene Set \code{\link[RangeshiftR]{Genetics}} parameters
#' @param init Set \code{\link[RangeshiftR]{Initialise}} parameters
#' @return returns a RangeShiftR parameter master object (class 'RSparams')
#' @details The cell size (resolution) is specified by the user in meters. It is important to note an essential difference in spatial scale between
#' the cell-based and the patch-based version. In the cell-based model, the cell resolution represents the spatial scale at which the two fundamental
#' processes of population dynamics and dispersal happen. This means that all the density-dependencies in the model (reproduction, survival,
#' emigration, settlement, etc...) act at the cell scale and the same scale is used as a single step unit for discrete movement models. In
#' the patch-based version, two spatial scales are simultaneously present: the cell scale, which in this case is used just for the transfer phase
#' of dispersal (movements) and the patch scale, at which the density-dependences are acting. The choice of type of model and cell resolution (as
#' well as the definition/scale of patches) is of fundamental importance because, depending on the system and on the question being tackled, it can
#' systematically bias the outcomes of the model.
#'
#' The user also defines the temporal scales. There are three distinct temporal scales. The highest-level one has years as units and represents the
#' scale at which variations in the abiotic environment are modelled (RangeShiftR does not explicitly model within-year variability in conditions).
#' The intermediate scale is the species’ reproductive season. The model can be used to simulate the case where there is only one reproductive
#' season per year but it is also possible to simulate situations where there more than one per year or only one every \eqn{N} years. A single
#' reproductive event is always followed by dispersal. Finally, the smallest time scale is represented by the number of steps that emigrants take
#' during the movement phase of dispersal. This can be determined by a maximum number of steps, per-step mortality or both.
#'
#' Demographic stochasticity is fundamentally important for the dynamics of populations that are naturally small or have declined to low abundances owing to
#' anthropogenic pressures. Additionally, inter-individual variability within populations can have a major influence on dynamics. Modelling stochastic events
#' that happen to individuals is crucial for avoiding systematic overestimation of population viability or rate of spread (Clark et al. 2001; Kendall & Fox 2003;
#' Robert et al. 2003; Grimm & Railsback 2005; Jongejans et al. 2008; Travis et al. 2011). Thus, population dynamics in RangeShifter were constructed to be
#' fully individual-based and stochastic. Each reproductive individual produces a discrete number of offspring sampled from a Poisson distribution with a mean
#' that is influenced by the species’ demographic parameters and the local population density. As RangeShifter has been designed for modelling a variety of
#' species with different life-history traits, a range of different population models can be chosen, depending on the species being modelled and on the
#' available information. In all cases demographic stochasticity is implemented.
#'
#' \emph{Cell-based vs. patch-based model} RangeShiftR can be run as a cell-based or patch-based model (Bian 2003). It should be noted
#' that the selection between cell-based or patch-based model is of fundamental importance for population dynamics calculations because
#' it influences the spatial extent at which density dependence operates. In both cases, the landscape is represented as a grid with cells
#' belonging to a particular habitat type, holding proportions of different habitats or being assigned a habitat quality index. However,
#' when RangeShiftR is run using the cell-based setting, the cell is the scale at which processes such as population dynamics and dispersal
#' act. The individuals present in a cell define a distinct population, and density-dependencies for reproduction, emigration and settlement
#' all operate at this scale. Even in the case where two habitat cells are adjacent, they still hold separate populations. In contrast, in
#' the patch-based model, population dynamics happen at the patch level, a patch being an assemblage of landscape cells of potentially
#' different habitat types. Patches are not defined automatically by RangeShiftR. Rather, the user is required to define which cells belong
#' to which patch, taking into account the ecological understanding of the study species. Density-dependencies regarding reproduction,
#' development, survival, emigration and settlement will depend on the density of individuals in a patch. However, discrete step-wise
#' movements during the transfer phase will always use the cell as the resolution at which steps occur, thus retaining important information
#' about the landscape heterogeneity.
#'
#' The choice between cell- and patch-based modelling can be of crucial importance. While a cell-based model provides an excellent abstraction
#' of space for many theoretical studies, for some applied studies it may be insufficient. This is because the misrepresentation of population
#' dynamics and dispersal (in terms of the scale at which they operate) can lead to substantial biases in projections regarding, for example,
#' rate of range expansion and population persistence (Bocedi et al. 2012b). Ideally, the scales at which population dynamics and dispersal
#' processes are modelled (by choosing the cell resolution or by defining the patches) should be those that are relevant for the species.
#' Importantly, the patch-based implementation allows separating the scales used for population dynamics and movements. In this case, the
#' landscape can be modelled at very fine resolution in order to capture the features that are likely to influence movements (e.g. narrow linear
#' features) without constraining the local population dynamics to operate at too small a scale.
#'
#' @export
RSsim <- function(batchnum = 1L,
                  simul = NULL,
                  land = NULL,
                  demog = NULL,
                  dispersal = NULL,
                  gene = NULL,
                  init = NULL){
    args <- as.list(match.call())
    # filter for names in ... that are also in slots(ControlParams) and pass them on
    s <- RSparams(control = ControlParams(batchnum = batchnum),
                  simul = Simulation(),
                  land = ArtificialLandscape(),
                  demog = Demography(Rmax = 1.5),
                  dispersal = Dispersal(),
                  gene = Genetics(),
                  init = Initialise())
    if (!is.null(args$land))  s <- s + land
    if (!is.null(args$simul)) s <- s + simul
    if (!is.null(args$demog)) s <- s + demog
    if (!is.null(args$dispersal))  s <- s + dispersal
    if (!is.null(args$gene))  s <- s + gene
    if (!is.null(args$init))  s <- s + init
    # check validity
    validObject(s)
    return(s)
}


#' Validate a given RS parameter object
#'
#' @export
validateRSparams <- function(x){validObject(x)}

#' Run the Simulation
#'
#' @param RSparams A parameter master object (class 'RSparams'), contains all parameters to specify the simulation settings.
#' @param dirpath File path to RS working directory; must contain the folders 'Inputs', 'Outputs', 'Output_Maps'.\cr
#' If NULL, the current \code{R} working directory will be used.
#' @return returns an error code or, if \code{RSparams@simul@ReturnPopRaster = TRUE}, a data frame with population data output
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
    #ParMaster_name <- deparse(substitute(RSparams))
    #run_from_R(RSparams, dirpath)

    if(RSparams@control@seed>0){
        #print(paste0("set seed ",RSparams@control@seed))
        #set.seed(RSparams@control@seed)
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


### Define Addition operations

#' @export
setGeneric("+")

setMethod("+", signature(e1 = "RSparams", e2 = "SimulationParams"), function(e1, e2) {
    validObject(e2)
    e1@simul <- e2
    return(e1)}
)

setMethod("+", signature(e1 = "RSparams", e2 = "LandParams"), function(e1, e2) {
    validObject(e2)
    if (class(e2)[1] == "ImportedLandscape") {
        if (e2@PatchFile=="NULL") {
            e1@control@patchmodel = FALSE
        }
        else {
            e1@control@patchmodel = TRUE
        }
        e1@control@resolution = e2@Resolution
        if (e2@HabitatQuality) {
            e1@control@landtype = 2L
            e1@control@maxNhab = 1L
        }
        else { # habitat codes
            e1@control@landtype = 0L
            e1@control@maxNhab = e2@Nhabitats
        }
        if (e2@SpDistFile=="NULL") {
            e1@control@speciesdist = FALSE
            e1@control@distresolution = -9L
        }
        else {
            e1@control@speciesdist = TRUE
            e1@control@distresolution = e2@SpDistResolution
        }
    }
    if (class(e2)[1] == "ArtificialLandscape") {
        e1@control@patchmodel = FALSE
        e1@control@resolution = e2@Resolution
        e1@control@landtype = 9L
        e1@control@maxNhab = 1L
        e1@control@speciesdist = FALSE
        e1@control@distresolution = 0L
    }
    e1@land = e2
    return(e1)}
)

setMethod("+", signature(e1 = "RSparams", e2 = "DemogParams"), function(e1, e2) {
    validObject(e2)
    e1@control@reproductn = e2@ReproductionType
    if (class(e2@StageStruct)[1] == "StagesParams") {
        e1@control@repseasons = e2@StageStruct@RepSeasons
        e1@control@stagestruct = TRUE
        e1@control@stages = e2@StageStruct@Stages
        if (length(e2@StageStruct@MinAge)==1 && e2@StageStruct@MinAge==0) {
            if (e2@ReproductionType==2) {
                e2@StageStruct@MinAge <- rep(0L, 2*e2@StageStruct@Stages-1)
            }
            else {
                e2@StageStruct@MinAge <- rep(0L, e2@StageStruct@Stages)
            }
        }
    }
    else {
        e1@control@repseasons = 1L
        e1@control@stagestruct = FALSE
        e1@control@stages = 1L
    }
    e1@demog <- e2
    return(e1)}
)

setMethod("+", signature(e1 = "DemogParams", e2 = "StagesParams"), function(e1, e2) {
    validObject(e2)
    e1@StageStruct <- e2
    e1@Rmax <- -9L
    e1@bc <- -9L
    return(e1)}
)

setMethod("+", signature(e1 = "RSparams", e2 = "DispersalParams"), function(e1, e2) {
    validObject(e2)
    if (class(e2@Transfer)[1] == "DispersalKernel") {
        e1@control@transfer = 0
    }
    else {
        if (class(e2@Transfer)[1] == "StochMove") {
            e1@control@transfer = 1
        }
        else {
            if (class(e2@Transfer)[1] == "CorrRW") {
                e1@control@transfer = 2
            }
            else {
                e1@control@transfer = NA_integer_
            }
        }
    }
    e1@dispersal <- e2
    return(e1)}
)

setMethod("+", signature(e1 = "DispersalParams", e2 = "EmigrationParams"), function(e1, e2) {
    validObject(e2)
    e1@Emigration <- e2
    validObject(e1)
    return(e1)}
)

setMethod("+", signature(e1 = "DispersalParams", e2 = "TransferParams"), function(e1, e2) {
    validObject(e2)
    e1@Transfer <- e2
    validObject(e1)
    return(e1)}
)

setMethod("+", signature(e1 = "DispersalParams", e2 = "SettlementParams"), function(e1, e2) {
    validObject(e2)
    e1@Settlement <- e2
    validObject(e1)
    return(e1)}
)

setMethod("+", signature(e1 = "RSparams", e2 = "GeneticsParams"), function(e1, e2) {
    validObject(e2)
    e1@gene <- e2
    return(e1)}
)

setMethod("+", signature(e1 = "RSparams", e2 = "InitialisationParams"), function(e1, e2) {
    validObject(e2)
    e1@init <- e2
    return(e1)}
)