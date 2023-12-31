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

# from RS 'Genetics' file

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
    if(anyNA(object@Architecture) || length(object@Architecture)!=1) {
        msg <- c(msg, "Architecture must be set!")
    }
    else {
        if(object@Architecture %in% c(0,1)) {
            if(object@Architecture == 0) {
                if(anyNA(object@NLoci) || length(object@NLoci)!=1) {
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
    if(anyNA(object@ProbMutn) || length(object@ProbMutn)!=1) {
        msg <- c(msg, "ProbMutn must be set!")
    }
    else {
        if(object@ProbMutn < 0.0 || object@ProbMutn > 1.0) {
            msg <- c(msg, "ProbMutn must be within the closed interval [0,1]!")
        }
    }
    if(anyNA(object@ProbCross) || length(object@ProbCross)!=1) {
        msg <- c(msg, "ProbCross must be set!")
    }
    else {
        if(object@ProbCross < 0.0 || object@ProbCross > 1.0) {
            msg <- c(msg, "ProbCross must be within the closed interval [0,1]!")
        }
    }
    if(anyNA(object@AlleleSD) || length(object@AlleleSD)!=1) {
        msg <- c(msg, "AlleleSD must be set!")
    }
    else {
        if(object@AlleleSD <= 0.0) {
            msg <- c(msg, "AlleleSD must be strictly positive!")
        }
    }
    if(anyNA(object@MutationSD) || length(object@MutationSD)!=1) {
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
    cat("   ProbMutn   =", object@ProbMutn, "\n")
    cat("   ProbCross  =", object@ProbCross, "\n")
    cat("   AlleleSD   =", object@AlleleSD, "\n")
    cat("   MutationSD =", object@MutationSD, "\n")
})
