% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpatialMultiomicAlignment.R
\name{generate_new_SM_counts}
\alias{generate_new_SM_counts}
\title{Generate new MALDI counts matrix for equivalent Visium spots
\itemize{
\item This function takes an original Spatial Metabolomics counts matrix, along with a metadata table (\code{obs_x}) containing information about the correspondence between MALDI and Visium spots. It generates a new counts matrix where each MALDI spot has an associated aggregated count based on its corresponding Visium spots.
}}
\usage{
generate_new_SM_counts(original_SM, obs_x, assay, slots, verbose = TRUE)
}
\arguments{
\item{original_SM}{A Seurat Spatial Metabolomics object containing the original counts matrix}

\item{obs_x}{A metadata table with information about the correspondence between MALDI and Visium spots. It should have columns 'Visium_spot' and 'MALDI_barcodes'.}

\item{assay}{Character string defining the Seurat assay that contains the annotated counts and metadata corresponding to the m/z values.}

\item{slots}{Vector of character strings describing which slots to pull the relative intensity values from (default = c("counts", "data")).}

\item{verbose}{Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).}
}
\value{
A data frame representing the new counts matrix for equivalent Visium spots,where each row corresponds to a Visium spot and columns correspond to m/z features.
}
\description{
Generate new MALDI counts matrix for equivalent Visium spots
\itemize{
\item This function takes an original Spatial Metabolomics counts matrix, along with a metadata table (\code{obs_x}) containing information about the correspondence between MALDI and Visium spots. It generates a new counts matrix where each MALDI spot has an associated aggregated count based on its corresponding Visium spots.
}
}
\examples{
## Generate new MALDI counts matrix for equivalent Visium spots
# new_counts <- generate_new_SM_counts(SeuratObj, obs_x, assay = "Spatial")
}
