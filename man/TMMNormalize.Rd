% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SeuratPreprocessing.R
\name{TMMNormalize}
\alias{TMMNormalize}
\title{Performs TMM normalization between categories based on a specified ident}
\usage{
TMMNormalize(
  combined.obj,
  ident,
  refIdent,
  normalisation.type = "CPM",
  CPM.scale.factor = 1e+06,
  assay = "Spatial",
  slot = "counts"
)
}
\arguments{
\item{combined.obj}{Seurat object that contains groups being normalized.}

\item{ident}{Character string defining the column name or ident group to normalize between.}

\item{refIdent}{Character string specifying one class/group type to use as a reference for TMM normalisation.}

\item{normalisation.type}{Character string defining the normalization method to run. Options are either c("CPM", "TIC", "LogNormalize") which represent counts per million (CPM), Total Ion Current (TIC) normalization or Log Normalization, respectively (default = "CPM").}

\item{CPM.scale.factor}{Numeric value that sets the scale factor for pixel/spot level normalization. Following normalization the total intensity value across each pixel will equal this value (default = 1e6).}

\item{assay}{Character string defining the name of the Seurat Object assay to pull the corresponding intensity data from (default = "Spatial").}

\item{slot}{Character string defining the name of the slot within the Seurat Object assay to pull the corresponding intensity data from (default = "counts").}
}
\value{
Seurat object with count values normalised and corrected for between categories
}
\description{
Performs TMM normalization between categories based on a specified ident
}
\examples{
# norm.data <- TMMNormalize(SeuratObj, ident = "samples", refIdent = "sample1", normalisation.type = "CPM")
}
