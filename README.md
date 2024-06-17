
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SpaMTP <img src="man/figures/logo.png" align="right" height="150" alt="" />

<br>
<!-- badges: start -->
x
<!-- badges: end -->

## New *R*-based User-Frendly Spatial Metabolomic, Transcriptomic and Proteomic Data Analysis Tool

<br>
<br>

The goal of SpaMTP is to integrate tools previously designed for
analysing different spatial metabolomic, transcriptomic and proteomic
data together into a user-freindly package. This will allow users to
seamlessly transition, integrate and co-analyse multimodality datasets
in particularly focusing on metabolomic/transcriptomics integration.
This package includes various functions to preform pre-processing,
spatial visualisation, various down-stream biological centered analyses,
data integration and data export of Spatial Metabolomic data. In
addition, this package has the ability to use both Cardinal (Spatial
Metabolomic based software) and Seurat (Spatial Transcriptomic based
software) functions.

<img src="inst/figures/SpaMTPSumFig.png" height="400" alt="" style="background-color: white;" />

## Another development version of the github is below
https://github.com/BCRL-tylu/SpaMTP

## Installation

You can install the development version of SpaMTP from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("agc888/SpaMTP")
```

Below is an example of how to set up SpaMTP using a conda environment

``` console
conda create -n SpaMTP -c conda-forge r-base=4.3.3 r-essentials r-devtools
conda activate SpaMTP
```

Within the SpaMTP environment open *R* and install ***SpaMTP***
packages

``` r

if (!require("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("BiomedicalMachineLearning/SpaMTP")
```

### Possible Installiation Errors

#### Cairo failed to install

``` console
checking if R was compiled with the RConn patch... no
checking for cairo.h... no
configure: error: Cannot find cairo.h! Please install cairo (http://www.cairographics.org/) and/or set CAIRO_CFLAGS/LIBS correspondingly.
ERROR: configuration failed for package ‘Cairo’
```

check these libraries exist:

``` console
conda install conda-forge::r-cairo
```

This should resolve any issues. Rerun:

``` r
BiocManager::install("scater")
```

#### Cardinal Failed to install

``` console
ERROR: dependency ‘EBImage’ is not available for package ‘Cardinal’

BiocManager::install("EBImage")

fftwtools.c:28:9: fatal error: fftw3.h: No such file or directory
   28 | #include<fftw3.h>
      |         ^~~~~~~~~
```

If EBImage failed to installed it ist most likely due to an issue with the installation of *fftwtools*

Try in *R*:

```r

if (!require("BiocManager", quietly = TRUE)) #Check if BiocManager is installed
    install.packages("BiocManager")

BiocManager::install("fftwtools")

```

Else try:

``` console
conda install bioconda::r-fftwtools
conda install conda-forge::fftw
```

If *install.packages(“fftwtools”)* still fails then try this:

``` console
conda install bioconda::bioconductor-cardinal  
```

For other issue please flag on github under *Issues*

## SpaMTP Vignette

To come:

``` r
## To Come ...
```
