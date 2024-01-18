#### THIS FILE ONLY CONTAINS DOCUMENTATION FOR THE ATTACHED CLEANED METABOLITE DATASETS ####


#' @title HMDB_1_names: A cleaned version of the reference metabolomics dataset from the Human Metabolome Database (HMDB)
#'
#' @description A dataset containing metabolite information from the HMDB,
#' including chemical formulas, exact masses, isomers, InChIKeys, and names.
#'
#' @format A data frame with 26190 rows and 37 variables:
#' \describe{
#'   \item{formula}{Character string representing the molecular formula.}
#'   \item{exactmass}{Numeric value indicating the exact mass.}
#'   \item{isomers}{Character string listing possible isomers (if applicable).}
#'   \item{isomers_inchikey}{Character string listing InChIKeys for isomers.}
#'   \item{isomers_names}{Character string listing names of isomers.}
#'   ... (additional columns for various adduct forms)
#' }
#'
#' @source \url{https://hmdb.ca/}
#'
#' @references Wishart DS, et al. (2018). HMDB 4.0: the human metabolome database for 2018. Nucleic Acids Res. 46(D1):D608-D617.
#'
#' @examples
#' head(HMDB_data)
#' dim(HMDB_data)
HMDB_1_names <- load("data/HMDB_1_names.rda")


#' @title Chebi_1_names: Cleaned ChEBI `(Chemical entities of biological interest)` reference dataset
#'
#' @description
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




#' @title Lipidmaps_1_names: A cleaned version of the lipid database from LIPID MAPS
#'
#' @description
#' This object contains a collection of lipids from the LIPID MAPS Structure Database.
#'
#' @format A data frame with 9493 rows and 37 variables:
#' \describe{
#'   \item{formula}{Molecular formula of the lipid (character)}
#'   \item{exactmass}{Exact mass of the lipid (double)}
#'   \item{isomers}{Comma-separated list of LIPID MAPS identifiers for isomers (character)}
#'   \item{isomers_inchikey}{Comma-separated list of InChIKeys for isomers (character)}
#'   \item{isomers_names}{Comma-separated list of names for isomers (character)}
#'   \item{2M-H}{Mass of the lipid as a 2M-H ion (double)}
#'   \item{2M+ACN+H}{Mass of the lipid as a 2M+ACN+H ion (double)}
#'   \item{...}{Additional columns for other adducts (double)}
#' }
#'
#' @source \url{https://www.lipidmaps.org/}
#'
#' @examples
#' # Access the molecular formula of the first lipid
#' LIPIDMAPS_db$formula`[1]`

#' # Access the exact mass of the third lipid
#' LIPIDMAPS_db$exactmass`[3]`

#' # View all information for a specific lipid
#' LIPIDMAPS_db`[LIPIDMAPS_db$isomers == "LMSP0601GB03", ]`
Lipidmaps_1_names <- load("data/Lipidmaps_1_names.rda")



#' @title GNPS_db: A cleaned database of metabolites from GNPS
#'
#' @description
#' This object contains a collection of metabolites from the Global Natural Products Social Molecular Networking (GNPS) platform.
#'
#' @format A data frame with 489 rows and 37 variables:
#' \describe{
#'   \item{formula}{Molecular formula of the metabolite (character)}
#'   \item{exactmass}{Exact mass of the metabolite (double)}
#'   \item{isomers}{Comma-separated list of GNPS identifiers for isomers (character)}
#'   \item{isomers_inchikey}{Comma-separated list of InChIKeys for isomers (character)}
#'   \item{isomers_names}{Comma-separated list of names for isomers (character)}
#'   \item{2M-H}{Mass of the metabolite as a 2M-H ion (double)}
#'   \item{2M+ACN+H}{Mass of the metabolite as a 2M+ACN+H ion (double)}
#'   \item{...}{Additional columns for other adducts (double)}
#' }
#'
#' @source \url{https://gnps.ucsd.edu/: https://gnps.ucsd.edu/}
#'
#' @examples
#' # Access the molecular formula of the first metabolite
#' GNPS_db$formula`[1]`
#'
#' # Access the exact mass of the third metabolite
#' GNPS_db$exactmass`[3]`
#'
#' # View all information for a specific metabolite
#' GNPS_db`[GNPS_db$isomers == "CMP17912", ]`
GNPS_1_names <- load("data/GNPS_1_names.rda")





