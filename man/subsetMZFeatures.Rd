% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MZAnnotation.R
\name{subsetMZFeatures}
\alias{subsetMZFeatures}
\title{Subset a Seurat Spatial Metabolomic object by list of m/z's}
\usage{
subsetMZFeatures(data, features, assay = "Spatial")
}
\arguments{
\item{data}{A Seurat Spatial Metabolomic Object for subsetting.}

\item{features}{A list of character strings defining the features/mz values to subset against.}

\item{assay}{A character string identifying the Seurat assay which contains the count data being subset.}
}
\value{
A subset Seurat object containing only m/z values that were specified
}
\description{
Subset a Seurat Spatial Metabolomic object by list of m/z's
}
\examples{
# subsetMZFeatures(SeuratObj, c("mz-160","mz-170","mz-180"))
}
