#### THIS FILE ONLY CONTAINS DOCUMENTATION FOR THE ATTACHED CLEANED METABOLITE DATASETS ####


#' @title HMDB_db: A cleaned version of the reference metabolomics dataset from the Human Metabolome Database (HMDB)
#'
#' @description A dataset containing metabolite information from the HMDB,
#' including chemical formulas, exact masses, isomers, InChIKeys, and names.
#'
#' @format ## A data frame with 26190 rows and 37 variables:
#' \describe{
#'   \item{formula}{Chemical formula of the entity (character)}
#'   \item{exactmass}{Exact mass of the entity (numeric)}
#'   \item{isomers}{Presence of isomers (character)}
#'   \item{isomers_inchikey}{InChIKeys for isomers (character)}
#'   \item{isomers_names}{Names of isomers (character)}
#'   \item{2M-H}{Mass of the entity as a 2M-H ion (double)}
#'   \item{2M+ACN+H}{Mass of the entity as a 2M+ACN+H ion (double)}
#'   \item{2M+ACN+Na}{Mass of the entity as a 2M+ACN+Na ion (double)}
#'   \item{2M+FA-H}{Mass of the entity as a 2M+FA-H ion (double)}
#'   \item{2M+H}{Mass of the entity as a 2M+H ion (double)}
#'   \item{2M+Hac-H}{Mass of the entity as a 2M+Hac-H ion (double)}
#'   \item{2M+K}{Mass of the entity as a 2M+K ion (double)}
#'   \item{2M+NH4}{Mass of the entity as a 2M+NH4 ion (double)}
#'   \item{2M+Na}{Mass of the entity as a 2M+Na ion (double)}
#'   \item{3M-H}{Mass of the entity as a 3M-H ion (double)}
#'   \item{M+2ACN+H}{Mass of the entity as a M+2ACN+H ion (double)}
#'   \item{M+2K+H}{Mass of the entity as a M+2K+H ion (double)}
#'   \item{M+2Na-H}{Mass of the entity as a M+2Na-H ion (double)}
#'   \item{M+ACN+H}{Mass of the entity as a M+ACN+H ion (double)}
#'   \item{M+ACN+Na}{Mass of the entity as a M+ACN+Na ion (double)}
#'   \item{M+Br}{Mass of the entity as a M+Br ion (double)}
#'   \item{M+CH3OH+H}{Mass of the entity as a M+CH3OH+H ion (double)}
#'   \item{M+Cl}{Mass of the entity as a M+Cl ion (double)}
#'   \item{M+DMSO+H}{Mass of the entity as a M+DMSO+H ion (double)}
#'   \item{M+FA-H}{Mass of the entity as a M+FA-H ion (double)}
#'   \item{M+H}{Mass of the entity as a M+H ion (double)}
#'   \item{M+Hac-H}{Mass of the entity as a M+Hac-H ion (double)}
#'   \item{M+IsoProp+H}{Mass of the entity as a M+IsoProp+H ion (double)}
#'   \item{M+IsoProp+Na+H}{Mass of the entity as a M+IsoProp+Na+H ion (double)}
#'   \item{M+K}{Mass of the entity as a M+K ion (double)}
#'   \item{M+K-2H}{Mass of the entity as a M+K-2H ion (double)}
#'   \item{M+NH4}{Mass of the entity as a M+NH4 ion (double)}
#'   \item{M+Na}{Mass of the entity as a M+Na ion (double)}
#'   \item{M+Na-2H}{Mass of the entity as a M+Na-2H ion (double)}
#'   \item{M+TFA-H}{Mass of the entity as a M+TFA-H ion (double)}
#'   \item{M-H}{Mass of the entity as a M-H ion (double)}
#'   \item{M-H2O-H}{Mass of the entity as a M-H2O-H ion (double)}
#'   ...
#' }
#'
#' @source <https://hmdb.ca/>
#'
#' @references Wishart DS, et al. (2018). HMDB 4.0: the human metabolome database for 2018. Nucleic Acids Res. 46(D1):D608-D617.
HMDB_db <- load("data/HMDB_db.rda")



#' @title Lipidmaps_db: A cleaned version of the lipid database from LIPID MAPS
#'
#' @description This object contains a collection of lipids from the LIPID MAPS Structure Database.
#'
#' @format ## A data frame with 9493 rows and 37 variables:
#' \describe{
#'   \item{formula}{Chemical formula of the entity (character)}
#'   \item{exactmass}{Exact mass of the entity (numeric)}
#'   \item{isomers}{Presence of isomers (character)}
#'   \item{isomers_inchikey}{InChIKeys for isomers (character)}
#'   \item{isomers_names}{Names of isomers (character)}
#'   \item{2M-H}{Mass of the entity as a 2M-H ion (double)}
#'   \item{2M+ACN+H}{Mass of the entity as a 2M+ACN+H ion (double)}
#'   \item{2M+ACN+Na}{Mass of the entity as a 2M+ACN+Na ion (double)}
#'   \item{2M+FA-H}{Mass of the entity as a 2M+FA-H ion (double)}
#'   \item{2M+H}{Mass of the entity as a 2M+H ion (double)}
#'   \item{2M+Hac-H}{Mass of the entity as a 2M+Hac-H ion (double)}
#'   \item{2M+K}{Mass of the entity as a 2M+K ion (double)}
#'   \item{2M+NH4}{Mass of the entity as a 2M+NH4 ion (double)}
#'   \item{2M+Na}{Mass of the entity as a 2M+Na ion (double)}
#'   \item{3M-H}{Mass of the entity as a 3M-H ion (double)}
#'   \item{M+2ACN+H}{Mass of the entity as a M+2ACN+H ion (double)}
#'   \item{M+2K+H}{Mass of the entity as a M+2K+H ion (double)}
#'   \item{M+2Na-H}{Mass of the entity as a M+2Na-H ion (double)}
#'   \item{M+ACN+H}{Mass of the entity as a M+ACN+H ion (double)}
#'   \item{M+ACN+Na}{Mass of the entity as a M+ACN+Na ion (double)}
#'   \item{M+Br}{Mass of the entity as a M+Br ion (double)}
#'   \item{M+CH3OH+H}{Mass of the entity as a M+CH3OH+H ion (double)}
#'   \item{M+Cl}{Mass of the entity as a M+Cl ion (double)}
#'   \item{M+DMSO+H}{Mass of the entity as a M+DMSO+H ion (double)}
#'   \item{M+FA-H}{Mass of the entity as a M+FA-H ion (double)}
#'   \item{M+H}{Mass of the entity as a M+H ion (double)}
#'   \item{M+Hac-H}{Mass of the entity as a M+Hac-H ion (double)}
#'   \item{M+IsoProp+H}{Mass of the entity as a M+IsoProp+H ion (double)}
#'   \item{M+IsoProp+Na+H}{Mass of the entity as a M+IsoProp+Na+H ion (double)}
#'   \item{M+K}{Mass of the entity as a M+K ion (double)}
#'   \item{M+K-2H}{Mass of the entity as a M+K-2H ion (double)}
#'   \item{M+NH4}{Mass of the entity as a M+NH4 ion (double)}
#'   \item{M+Na}{Mass of the entity as a M+Na ion (double)}
#'   \item{M+Na-2H}{Mass of the entity as a M+Na-2H ion (double)}
#'   \item{M+TFA-H}{Mass of the entity as a M+TFA-H ion (double)}
#'   \item{M-H}{Mass of the entity as a M-H ion (double)}
#'   \item{M-H2O-H}{Mass of the entity as a M-H2O-H ion (double)}
#'   ...
#' }
#'
#' @source <https://www.lipidmaps.org/>
Lipidmaps_db <- load("data/Lipidmaps_db.rda")




#' @title Chebi_db: Cleaned ChEBI `(Chemical entities of biological interest)` reference dataset
#'
#' @description A dataset containing information on chemical entities, including their formulas, exact masses, isomers, InChIKeys, and names.
#'
#' @format ## A data frame with 46297 rows and 37 variables:
#' \describe{
#'   \item{formula}{Chemical formula of the entity (character)}
#'   \item{exactmass}{Exact mass of the entity (numeric)}
#'   \item{isomers}{Presence of isomers (character)}
#'   \item{isomers_inchikey}{InChIKeys for isomers (character)}
#'   \item{isomers_names}{Names of isomers (character)}
#'   \item{2M-H}{Mass of the entity as a 2M-H ion (double)}
#'   \item{2M+ACN+H}{Mass of the entity as a 2M+ACN+H ion (double)}
#'   \item{2M+ACN+Na}{Mass of the entity as a 2M+ACN+Na ion (double)}
#'   \item{2M+FA-H}{Mass of the entity as a 2M+FA-H ion (double)}
#'   \item{2M+H}{Mass of the entity as a 2M+H ion (double)}
#'   \item{2M+Hac-H}{Mass of the entity as a 2M+Hac-H ion (double)}
#'   \item{2M+K}{Mass of the entity as a 2M+K ion (double)}
#'   \item{2M+NH4}{Mass of the entity as a 2M+NH4 ion (double)}
#'   \item{2M+Na}{Mass of the entity as a 2M+Na ion (double)}
#'   \item{3M-H}{Mass of the entity as a 3M-H ion (double)}
#'   \item{M+2ACN+H}{Mass of the entity as a M+2ACN+H ion (double)}
#'   \item{M+2K+H}{Mass of the entity as a M+2K+H ion (double)}
#'   \item{M+2Na-H}{Mass of the entity as a M+2Na-H ion (double)}
#'   \item{M+ACN+H}{Mass of the entity as a M+ACN+H ion (double)}
#'   \item{M+ACN+Na}{Mass of the entity as a M+ACN+Na ion (double)}
#'   \item{M+Br}{Mass of the entity as a M+Br ion (double)}
#'   \item{M+CH3OH+H}{Mass of the entity as a M+CH3OH+H ion (double)}
#'   \item{M+Cl}{Mass of the entity as a M+Cl ion (double)}
#'   \item{M+DMSO+H}{Mass of the entity as a M+DMSO+H ion (double)}
#'   \item{M+FA-H}{Mass of the entity as a M+FA-H ion (double)}
#'   \item{M+H}{Mass of the entity as a M+H ion (double)}
#'   \item{M+Hac-H}{Mass of the entity as a M+Hac-H ion (double)}
#'   \item{M+IsoProp+H}{Mass of the entity as a M+IsoProp+H ion (double)}
#'   \item{M+IsoProp+Na+H}{Mass of the entity as a M+IsoProp+Na+H ion (double)}
#'   \item{M+K}{Mass of the entity as a M+K ion (double)}
#'   \item{M+K-2H}{Mass of the entity as a M+K-2H ion (double)}
#'   \item{M+NH4}{Mass of the entity as a M+NH4 ion (double)}
#'   \item{M+Na}{Mass of the entity as a M+Na ion (double)}
#'   \item{M+Na-2H}{Mass of the entity as a M+Na-2H ion (double)}
#'   \item{M+TFA-H}{Mass of the entity as a M+TFA-H ion (double)}
#'   \item{M-H}{Mass of the entity as a M-H ion (double)}
#'   \item{M-H2O-H}{Mass of the entity as a M-H2O-H ion (double)}
#'   ...
#' }
#'
#' @source <https://www.ebi.ac.uk/chebi/>
Chebi_db <- load("data/Chebi_db.rda")


#' @title GNPS_db: A cleaned database of metabolites from GNPS
#'
#' @description This object contains a collection of metabolites from the Global Natural Products Social Molecular Networking (GNPS) platform.
#'
#' @format ## A data frame with 489 rows and 37 variables:
#' \describe{
#'   \item{formula}{Chemical formula of the entity (character)}
#'   \item{exactmass}{Exact mass of the entity (numeric)}
#'   \item{isomers}{Presence of isomers (character)}
#'   \item{isomers_inchikey}{InChIKeys for isomers (character)}
#'   \item{isomers_names}{Names of isomers (character)}
#'   \item{2M-H}{Mass of the entity as a 2M-H ion (double)}
#'   \item{2M+ACN+H}{Mass of the entity as a 2M+ACN+H ion (double)}
#'   \item{2M+ACN+Na}{Mass of the entity as a 2M+ACN+Na ion (double)}
#'   \item{2M+FA-H}{Mass of the entity as a 2M+FA-H ion (double)}
#'   \item{2M+H}{Mass of the entity as a 2M+H ion (double)}
#'   \item{2M+Hac-H}{Mass of the entity as a 2M+Hac-H ion (double)}
#'   \item{2M+K}{Mass of the entity as a 2M+K ion (double)}
#'   \item{2M+NH4}{Mass of the entity as a 2M+NH4 ion (double)}
#'   \item{2M+Na}{Mass of the entity as a 2M+Na ion (double)}
#'   \item{3M-H}{Mass of the entity as a 3M-H ion (double)}
#'   \item{M+2ACN+H}{Mass of the entity as a M+2ACN+H ion (double)}
#'   \item{M+2K+H}{Mass of the entity as a M+2K+H ion (double)}
#'   \item{M+2Na-H}{Mass of the entity as a M+2Na-H ion (double)}
#'   \item{M+ACN+H}{Mass of the entity as a M+ACN+H ion (double)}
#'   \item{M+ACN+Na}{Mass of the entity as a M+ACN+Na ion (double)}
#'   \item{M+Br}{Mass of the entity as a M+Br ion (double)}
#'   \item{M+CH3OH+H}{Mass of the entity as a M+CH3OH+H ion (double)}
#'   \item{M+Cl}{Mass of the entity as a M+Cl ion (double)}
#'   \item{M+DMSO+H}{Mass of the entity as a M+DMSO+H ion (double)}
#'   \item{M+FA-H}{Mass of the entity as a M+FA-H ion (double)}
#'   \item{M+H}{Mass of the entity as a M+H ion (double)}
#'   \item{M+Hac-H}{Mass of the entity as a M+Hac-H ion (double)}
#'   \item{M+IsoProp+H}{Mass of the entity as a M+IsoProp+H ion (double)}
#'   \item{M+IsoProp+Na+H}{Mass of the entity as a M+IsoProp+Na+H ion (double)}
#'   \item{M+K}{Mass of the entity as a M+K ion (double)}
#'   \item{M+K-2H}{Mass of the entity as a M+K-2H ion (double)}
#'   \item{M+NH4}{Mass of the entity as a M+NH4 ion (double)}
#'   \item{M+Na}{Mass of the entity as a M+Na ion (double)}
#'   \item{M+Na-2H}{Mass of the entity as a M+Na-2H ion (double)}
#'   \item{M+TFA-H}{Mass of the entity as a M+TFA-H ion (double)}
#'   \item{M-H}{Mass of the entity as a M-H ion (double)}
#'   \item{M-H2O-H}{Mass of the entity as a M-H2O-H ion (double)}
#'   ...
#' }
#'
#' @source \url{https://gnps.ucsd.edu/: https://gnps.ucsd.edu/}
GNPS_db <- load("data/GNPS_db.rda")

