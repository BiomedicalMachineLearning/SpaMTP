#' Cleaned ChEBI `(Chemical entities of biological interest)` reference dataset
#'
#' A dataset containing information on chemical entities, including their
#' formulas, exact masses, isomers, InChIKeys, and names.
#'
#' @format A data frame with 6 rows and 37 variables:
#' \describe{
#'   \item{formula}{Chemical formula of the entity (character)}
#'   \item{exactmass}{Exact mass of the entity (numeric)}
#'   \item{isomers}{Presence of isomers (character)}
#'   \item{isomers_inchikey}{InChIKeys for isomers (character)}
#'   \item{isomers_names}{Names of isomers (character)}
#'   \item{2M-H}{... (numeric)}  # Various other adducts
#'   ...
#' }
#'
#' @source ChEBI database (https://www.ebi.ac.uk/chebi/)
#'
#' @usage data(Chebi_1_names)
#'
#' @examples
#' head(Chebi_1_names)
#'
#' # Get information on a specific entity
#' filter(Chebi_1_names, formula == "C34H38Cl2N2O5")
Chebi_1_names <- load("data/Chebi_1_names.rda")
