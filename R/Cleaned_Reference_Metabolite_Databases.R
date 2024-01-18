#### THIS FILE ONLY CONTAINS DOCUMENTATION FOR THE ATTACHED CLEANED METABOLITE DATASETS ####


#' @title Lipidmaps_1_names: A cleaned version of the lipid database from LIPID MAPS
#'
#' @description This object contains a collection of lipids from the LIPID MAPS Structure Database.
#'
#' @format A data frame with 9493 rows and 37 variables:
#' \describe{
#'   \item{formula}{Chemical formula of the entity (character)}
#'   \item{exactmass}{Exact mass of the entity (numeric)}
#'   \item{isomers}{Presence of isomers (character)}
#'   \item{isomers_inchikey}{InChIKeys for isomers (character)}
#'   \item{isomers_names}{Names of isomers (character)}
#'   \item{2M-H}{Mass of the entity as a 2M-H ion (double)}
#'   ... (other variables related to mass as different ions)
#' }
#'
#' @source \url{https://www.lipidmaps.org/}
#'
#' @examples
#' # Access the molecular formula of the first lipid
#' Lipidmaps_1_names$formula[1]
#'
#' # Access the exact mass of the third lipid
#' Lipidmaps_1_names$exactmass[3]
#'
#' # View all information for a specific lipid
#' Lipidmaps_1_names[LIPIDMAPS_db$isomers == "LMSP0601GB03", ]
Lipidmaps_1_names <- load("data/Lipidmaps_1_names.rda")
