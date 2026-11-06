#---------------------------------------------------------------------------
#
#	Copyright (C) 2020-2026 Anne-Kathleen Malchow, Jette Wolff, Théo Pannetier, Roslyn Henry, Greta Bocedi, Stephen C.F. Palmer, Justin M.J. Travis, Damaris Zurell
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


#' @details
#' \code{RangeShiftR} provides flexible and fast simulations of spatial eco-evolutionary dynamics.
#' It implements the individual-based simulation software \emph{RangeShifter},
#' \insertCite{bocedi2014,bocedi2020}{RangeShiftR},
#' making it available for all machines that can build and run \emph{R} \insertCite{malchow2021}{RangeShiftR}.
#'
#' \code{RangeShiftR} models the processes of demography, dispersal and evolution in an
#' inter-dependent way, offering substantial complexity in the corresponding modelling choices.
#' It is entirely open-source and aims to facilitate the application
#' of individual-based and mechanistic modelling to eco-evolutionary questions.
#'
#' The \emph{RangeShifter} project is presented on our website
#' \url{https://rangeshifter.github.io/}.
#'
#' To get started with the package, please find an overview and a range of tutorials here:
#' \url{https://rangeshifter.github.io/RangeshiftR-tutorials/}.
#'
#' \subsection{Workflow}{
#' The figure shows the typical workflow. The parameters of a simulation are set in parameter modules, each created by its own function:
#' \code{\link[RangeShiftR]{Simulation}}, \code{\link[RangeShiftR]{ImportedLandscape}} or \code{\link[RangeShiftR]{ArtificialLandscape}},
#' \code{\link[RangeShiftR]{Demography}}, \code{\link[RangeShiftR]{Dispersal}}, \code{\link[RangeShiftR]{Genetics}},
#' \code{\link[RangeShiftR]{Management}} and \code{\link[RangeShiftR]{Initialise}}.
#' Some modules contain sub-modules: \code{\link[RangeShiftR]{StageStructure}} in \code{Demography};
#' \code{\link[RangeShiftR]{Emigration}}, \code{\link[RangeShiftR]{DispersalKernel}}, \code{\link[RangeShiftR]{SMS}} or
#' \code{\link[RangeShiftR]{CorrRW}}, and \code{\link[RangeShiftR]{Settlement}} in \code{Dispersal};
#' \code{\link[RangeShiftR]{Traits}} with \code{\link[RangeShiftR]{NeutralTraits}}, \code{\link[RangeShiftR]{GeneticLoadTraits}}
#' and the dispersal traits in \code{Genetics}; \code{\link[RangeShiftR]{Translocation}} and \code{\link[RangeShiftR]{Harvesting}} in \code{Management}.
#'
#' \code{\link[RangeShiftR]{RSsim}} combines the modules into the parameter master object and validates it.
#' \code{\link[RangeShiftR]{RunRS}} runs the simulation (see there for the order of processes). The results are written to the folder
#' \file{Outputs}; population data can also be returned to \emph{R} (see \code{ReturnPopMatrix} and \code{ReturnPopDataFrame}
#' in \code{\link[RangeShiftR]{Simulation}}). They can be analysed with \code{\link[RangeShiftR]{readPop}},
#' \code{\link[RangeShiftR]{readRange}}, \code{\link[RangeShiftR]{plotAbundance}}, \code{\link[RangeShiftR]{plotOccupancy}},
#' \code{\link[RangeShiftR]{ColonisationStats}} and \code{\link[RangeShiftR]{SMSpathLengths}}.
#'
#' \if{html}{\figure{package_workflow.svg}{options: width="680" alt="RangeShiftR workflow"}}
#' \if{latex}{\figure{package_workflow.pdf}{options: width=5in}}
#' }
#'
#' @references
#'         \insertAllCited{}
#' @keywords internal
"_PACKAGE"


#' @useDynLib RangeShiftR
#' @importFrom Rcpp sourceCpp
#' @importFrom Rdpack reprompt
NULL


# Show start-up message upon loading the package
.onAttach <- function(libname, pkgname) {
    packageStartupMessage("RangeshiftR version 3.0.2 (21.10.2026)\n",
                          "Copyright (C) 2020-2026 Anne-Kathleen Malchow, Greta Bocedi, Stephen C.F. Palmer, Justin M.J. Travis, Jette Wolff, Théo Pannetier, Roslyn Henry, Damaris Zurell\n\n",
                          "This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.\n",
                          "You are welcome to redistribute it and/or modify it under certain conditions;\n",
                          "type 'RangeShiftR_license()' for details.\n")
}

# Unload dynamic library when unloading the package
.onUnload <- function (libpath) {
    library.dynam.unload("RangeShiftR", libpath)
}


#' Display RangeShiftR license information
#'
#' @export
RangeShiftR_license <- function ()
{
    cat("\nRangeshiftR version 3.0.2 (21.10.2026)\n")
    cat("Copyright (C) 2020-2026 Anne-Kathleen Malchow, Greta Bocedi, Stephen C.F. Palmer, Justin M.J. Travis, Jette Wolff, Théo Pannetier, Roslyn Henry, Damaris Zurell\n\n")
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


#' Display citation
#'
#' @export
RangeShiftR_citation <- function ()
{
    citation(package = "RangeShiftR", lib.loc = NULL, auto = NULL)
}
