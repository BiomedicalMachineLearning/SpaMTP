% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LoadData.R
\name{ReadSM_mtx}
\alias{ReadSM_mtx}
\title{Read a Spatial Metabolomics image matrix file
NOTE: This file must be in a format similar to the one below:}
\usage{
ReadSM_mtx(
  mtx.file,
  assay = "Spatial",
  verbose = TRUE,
  feature.start.column = 1,
  mz.prefix = NULL
)
}
\arguments{
\item{mtx.file}{Character string defining the path of the spatial metabolomic image matrix .csv file}

\item{assay}{Character string of the Seurat object assay name to store the relative intensity data (default = "Spatial").}

\item{verbose}{Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).}

\item{feature.start.column}{Numeric value defining the start index containing the x, y and m/z value columns within the table (default = 1).}

\item{mz.prefix}{character string matching the prefix string in front of each m/z name (deafult = NULL).}
}
\value{
A SpaMTP Seurat class object containing the intensity values in the counts slot of the designated assay
}
\description{
\if{html}{\out{<div class="sourceCode">}}\preformatted{     A data.frame: 5 × 5
   x	   y    mz1   mz2   mz3
 <int> <int> <dbl> <dbl> <dbl>
}\if{html}{\out{</div>}}

1	 0	   1	   0	   0	   11
2	 0	   2	   0	   0	   0
3	 0	   3	   0	   0	   0
4	 0	   4	   20	   0	   0
5	 0	   5	   0	   0	   0
}
\details{
Where:
- The first 2 columns are labeled x and y containing the respective x/y spatial coordinates
- The next columns are then the respective m/z values and their intensities for each spatial pixel
}
\examples{
# msi_data <- ReadSM_mtx("~/Documents/msi_mtx.csv")
}
