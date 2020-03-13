
### CLASS GENETICSPARAMS

# from RS 'Genetics' file

#' Set genetics parameters
#'
#' @description Set genetics parameters and control genetics types.\cr
#' Furthermore, optionally define a genetics.
#'
# #' @author Anne-Kathleen Malchow
#' @usage Genetics(Architecture = 0, NLoci = 1, ArchFile = "NULL",
#'          ProbMutn = 0.0, ProbCross = 0.0,
#'          AlleleSD = 0.1, MutationSD = 0.1)
#' @param Architecture Genetic Architecture: \cr 0 = One chromosome per trait (default), \cr 1 = Read from file (set \code{ArchFile}).
#' @param NLoci Required if \code{Architecture=0}: Number of loci per chromosome, defaults to \eqn{1} (integer).
#' @param ArchFile Required if \code{Architecture=1}: Name of the genetic architecture file.
#' @param ProbMutn Probability of mutation of an individual allele at meiosis, defaults to \eqn{0.0}.\cr
#' Must be \eqn{0 \le}\code{ProbMutn} \eqn{\le 1}.
#' @param ProbCross Probability of crossover at an individual locus at meiosis, defaults to \eqn{0.0}.\cr
#' Must be \eqn{0 \le}\code{ProbCross} \eqn{\le 1}.
#' @param AlleleSD Standard deviation of initial allelic values around phenotypic value, defaults to \eqn{0.1}.\cr Must be \eqn{> 0}.
#' @param MutationSD Standard deviation of mutation magnitude, defaults to \eqn{0.1}.\cr Must be \eqn{> 0}.
#' @details In RangeShifter (v2.0), any heritable variable trait (at present limited to dispersal traits) is
#' controlled by a separate genetics module. In previous versions, the trait value (e.g. emigration probability,
#' mean step length of a CRW) was held by the individual directly as a ‘pseudo-gene’ for diploid species,
#' the mean of two such alleles controlled the phenotype), but this did not allow for representation of
#' more complex genetic phenomena, such as linkage between traits. This simple type of implementation
#' may still be represented approximately by choosing the ‘one chromosome per trait’ option, and setting
#' one locus per chromosome. Note that we use the term ‘chromosome’ here to represent both the single
#' chromosomes of a haploid species and the chromosome pair of a diploid species.
#'
#' \emph{Flexible genetic architecture}\cr
#' However, for a more realistic representation of heritable traits, an explicit genetic architecture
#' may be defined, which must be read from a text file. The file specifies how many chromosomes
#' each individual will have, the number of loci on each chromosome and which loci contribute
#' to the phenotypic value of each trait. Thus, for example, it is possible to model a species
#' exhibiting a single variable trait (e.g. the mean of a negative exponential dispersal kernel)
#' dependent on loci spreads across three chromosomes or a species exhibiting three trait (e.g.
#' density-dependent emigration) all of which are governed by loci located on a singlechromosome.
#' In practice, most genetic architectures are likely to fall somewhere between these extremes,
#' i.e. there will be several chromosomes and traits will be mapped across them. In contrast to
#' RangeShifter v1, whenever there are variable traits in a model, evolution is assumed, although
#' it can if desired be effectively eliminated for a haploid species by setting a very low mutation
#' probability.
#'
#' Traits are specified in the order they are defined in the model (emigration / movement / settlement),
#' and each trait is mapped to one or more loci. Any loci which do not contribute to at least one trait
#' are treated as neutral loci (see below). Neutral loci may also be specified for a model in which
#' there are no adaptive traits. All alleles are represent by integer values (positive or negative),
#' and the sum of alleles at all loci contributing to a trait (both alleles at each locus of a diploid
#' species) controls the phenotype. However, as phenotypic traits exist on several widely different
#' scales (e.g. emigration probability between \eqn{0} and \eqn{1}, dispersal kernel mean typically
#' many hundreds or thousands of metres), it is necessary to specify how the allelic scale relates
#' to the phenotypic scale. A scaling factor is specified for each trait, which governs how large a
#' change of \eqn{100} units (the ‘integer base’ of the genome) on the allele scale will be on the
#' phenotypic scale. For example, if the scaling factor for density-independent emigration probability
#' is \eqn{0.1}, and a juvenile’s sum of all alleles contributing to that trait is \eqn{150} less than
#' its parent’s equivalent sum, then the juvenile’s emigration probability phenotype will be \eqn{0.15}
#' lower than its parent’s phenotype (but subject in this case to the constraint that it may not be
#' less than zero).
#'
#' \emph{Genome initialisation}\cr
#' At model initialisation, individuals are assigned phenotypic traits drawn from a normal
#' distribution controlled by a specified mean phenotype and standard deviation, but also subject
#' to any phenotypic constraints (e.g. a probability must lie between \eqn{0} and \eqn{1}, step length
#' must be positive, etc.). The standard deviation for a trait may not be greater than the corresponding
#' scaling factor, but it may be substantially less if an initial population which is highly homogeneous
#' for a particular trait is required. The genome actually controls individual variation relative to the
#' initial population mean value of a trait; thus if the sum of an individual’s alleles is negative, its
#' phenotype will be less than the initial mean value, and if its sum is positive, its phenotype will be
#' greater than the mean value. However, to prevent all alleles for a multi-loci trait being identical in
#' initial individuals, further random variation is applied at the allele scale (i.e. common to all traits),
#' which is drawn from a normal distribution having zero mean and a specified standard deviation. Note
#' that therefore the observed variance in a trait value across the initial population may not match exactly
#' the specified variance for the trait.
#'
#' Side note: This differs from the implementation in version 1, which used a uniform distribution; hence
#' an initial population cannot be set up in v2.0 to have exactly the same properties as in v1, but a similar
#' equilibrium population should arise after a period of sufficiently strong selection for one or more traits.
#'
#' \emph{Pleiotropy, neutral loci and mutation}\cr
#' It is possible that a particular locus can be specified more than once for a particular trait; this
#' increases its weighting relative to other loci contributing to the trait. A locus may also code for
#' more than one trait, i.e. it is pleiotropic. Large positive allele values at that locus will tend to
#' lead to larger than average values of both traits (but subject to any other loci coding separately
#' for the two traits). Thus the two traits are forced (to some extent) to be positively correlated. It
#' is not possible to specify inverse correlation between two traits in this way.
#'
#' A locus which does not code for any trait is neutral, and therefore not subject directly to
#' selection, although the distribution across the population of allele values at that locus may vary
#' over time owing to genetic drift and/or linkage to loci which are under selection. Such neutral
#' loci may be used as markers for estimating relatedness at the individual or population level.
#' Mutation is governed by two genome-level parameters, the mutation probability and standard deviation,
#' which are applied in a standard way at the level of the individual locus (unlike in version 1, in
#' which separate probabilities and magnitudes were applied for separate traits). When a mutation occurs
#' at a locus, a random number is drawn from a normal distribution having zero mean and the mutation standard
#' deviation. This number is multiplied by the integer base to yield an integer value which is added to the
#' allele before it is copied to the juvenile’s chromosome. The default mutation standard deviation of \eqn{0.1}
#' will therefore give mutations which mostly range between \eqn{-30} and \eqn{+30} at the allele scale.
#' @references ???
#' @name Genetics
#' @export Genetics
Genetics <- setClass("GeneticsParams", slots = c(Architecture = "integer_OR_numeric",
                                                 NLoci = "integer_OR_numeric",
                                                 ArchFile = "character",
                                                 ProbMutn = "numeric",
                                                 ProbCross = "numeric",
                                                 AlleleSD = "numeric",
                                                 MutationSD = "numeric")
                                     , prototype = list(Architecture = 0L,
                                                        NLoci = 1L,
                                                        ArchFile = "NULL",
                                                        ProbMutn = 0.0,
                                                        ProbCross = 0.0,
                                                        AlleleSD = 0.1,
                                                        MutationSD = 0.1)
)
setValidity('GeneticsParams', function(object){
    msg <- NULL
    if(is.na(object@Architecture) || length(object@Architecture)!=1) {
        msg <- c(msg, "Architecture must be set!")
    }
    else {
        if(object@Architecture %in% c(0,1)) {
            if(object@Architecture == 0) {
                if(is.na(object@NLoci) || length(object@NLoci)!=1) {
                    msg <- c(msg, "NLoci must be set if Architecture = 0 (one chromosome per trait)!")
                }
                else {
                    if(object@NLoci < 1) {
                        msg <- c(msg, "NLoci must be greater than 0!")
                    }
                }
            }
            if(object@Architecture == 1) {
                if (object@ArchFile == "NULL"){
                    msg <- c(msg, 'ArchFile is required if Architecture = 1 (load architecture file).')
                }
            }
        }else{
            msg <- c(msg, "Architecture must be either 0 or 1!")
        }
    }
    if(is.na(object@ProbMutn) || length(object@ProbMutn)!=1) {
        msg <- c(msg, "ProbMutn must be set!")
    }
    else {
        if(object@ProbMutn < 0.0 || object@ProbMutn > 1.0) {
            msg <- c(msg, "ProbMutn must be within the closed interval [0,1]!")
        }
    }
    if(is.na(object@ProbCross) || length(object@ProbCross)!=1) {
        msg <- c(msg, "ProbCross must be set!")
    }
    else {
        if(object@ProbCross < 0.0 || object@ProbCross > 1.0) {
            msg <- c(msg, "ProbCross must be within the closed interval [0,1]!")
        }
    }
    if(is.na(object@AlleleSD) || length(object@AlleleSD)!=1) {
        msg <- c(msg, "AlleleSD must be set!")
    }
    else {
        if(object@AlleleSD <= 0.0) {
            msg <- c(msg, "AlleleSD must be strictly positive!")
        }
    }
    if(is.na(object@MutationSD) || length(object@MutationSD)!=1) {
        msg <- c(msg, "MutationSD must be set!")
    }
    else {
        if(object@MutationSD <= 0.0) {
            msg <- c(msg, "MutationSD must be strictly positive!")
        }
    }
    if (is.null(msg)) TRUE else msg}
)
setMethod('initialize', 'GeneticsParams', function(.Object, ...) {
    this_func = "Genetics(): "
    args <- list(...)
    .Object <- callNextMethod()
    if(.Object@Architecture == 0) {
        if (!is.null(args$ArchFile)) {
            warning(this_func, "ArchFile", warn_msg_ignored, "since Architecture = 0 (one chromosome per trait).", call. = FALSE)
        }
    }
    if(.Object@Architecture == 1) {
        if (!is.null(args$NLoci)) {
            warning(this_func, "NLoci", warn_msg_ignored, "since Architecture = 1 (load architecture file).", call. = FALSE)
        }
    }
    .Object}
)
setMethod("show", "GeneticsParams", function(object){
    cat(" Genetics: \n")
    cat("   Architecture =", object@Architecture, ": ")
    if (object@Architecture == 0) {
        cat("One chromosome per trait \n")
        cat("   with", object@NLoci)
        if (object@NLoci==1) cat(" locus")
        else cat(" loci")
        cat(" per chromosome.")
    }
    if (object@Architecture == 1) {
        cat("loaded from architecture file:\n")
        cat("   ", object@ArchFile)
    }
    cat("\n")
    cat("   ProbMutn =", object@ProbMutn, "\n")
    cat("   ProbCross =", object@ProbCross, "\n")
    cat("   AlleleSD =", object@AlleleSD, "\n")
    cat("   MutationSD =", object@MutationSD, "\n")
})