
#### CLASS CONTROLPARAMS

# from RS 'CONTROL.txt' file

ControlParams <- setClass("ControlParams", slots = c(
                                   #nSimuls = "integer_OR_numeric",         # Not yet used by R version, in C++ version its read from parameterfile
                                   #nLandscapes = "integer_OR_numeric",     # Not yet used by R version, in C++ version its read from landfile
                                   batchnum = "integer_OR_numeric",         # only variable to set from RSsim(), optional
                                   patchmodel = "logical",                  # set via +Land
                                   resolution = "integer_OR_numeric",       # set via +Land
                                   landtype = "integer_OR_numeric",         # set via +Land
                                   maxNhab = "integer_OR_numeric",          # set via +Land
                                   speciesdist = "logical",                 # set via +Land
                                   distresolution = "integer_OR_numeric",   # set via +Land
                                   reproductn = "integer_OR_numeric",       # set via +Demography
                                   repseasons = "integer_OR_numeric",       # set via +Demography if (StageStruct is Type "StagesParams") {repseasons = 1} else {demo@StageStruct@repseasons}
                                   stagestruct = "logical",                 # set via +Demography
                                   stages = "integer_OR_numeric",           # set via +Demography@StageStruct
                                   transfer = "integer_OR_numeric",         # set via +Dispersal     Transfer method: 0 = dispersal kernels, 1 = SMS, 2 = CRW)
                                   seed = "integer_OR_numeric")
                         ,prototype = list(#nSimuls = 1L,
                                  #nLandscapes = 1L,
                                  batchnum = 0L,
                                  patchmodel = FALSE,
                                  resolution = 100L,
                                  landtype = 9L,
                                  maxNhab = 1L,
                                  speciesdist = FALSE,
                                  distresolution = NA_integer_,
                                  reproductn = 0L,
                                  repseasons = 1L,
                                  stagestruct = FALSE,
                                  stages = NA_integer_,
                                  transfer = 0L,
                                  seed = 0L)
)
setValidity("ControlParams", function(object){
    msg <- NULL
    if(object@stagestruct) {
        stg <- as.integer(object@stages)
        if(stg < 2 || stg > 10) msg <- c(msg, "Number of Stages must be in the interval [2; 10]")
                                                          #paste("Unequal x,y lengths: ", length(object@x), ", ", length(object@y), sep="")
    }
    if (is.null(msg)) TRUE else msg}
)
setMethod("show", "ControlParams", function(object){cat(" Batch #", object@batchnum, "\n")})
