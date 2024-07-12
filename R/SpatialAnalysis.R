
#### SpaMTP Spatial Analysis ##############################################################################################################################################################################
## Code is manipulated from Seurat

#' Compute the row variances for each m/z value
#'
#' @param x Matrix containing count values used for correlation
#'
#' @return Vector contating the row variances for each m/z value
#'
#' @examples
#' #HELPER FUNCTION
RowVar <- function(x) {
  .Call('_Seurat_RowVar', PACKAGE = 'Seurat', x)
}


#' Finds metabolites that display strong spatial patterns using MoransI
#'
#' @param object SpaMTP Seurat class object contating the intensity values for each m/z
#' @param assay Character string indicating which Seurat object assay to pull data form (default = "SPM").
#' @param slot Character string indicating the assay slot to use to pull expression values form (default = "counts").
#' @param image Character string defining the image to extract the tissue coordinates from (defualt = "slice1").
#' @param nfeatures Numeric values defining the top number of features to mark as the top spatially variable (default = 2000).
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#'
#' @return Returns SpaMTP object containing the MoransI pvalue and rank stored in the assay feature meta.data
#' @export
#'
#' @examples
#' # SpaMTP.obj <- FindSpatiallyVariableMetabolites(SpaMTP)
FindSpatiallyVariableMetabolites <- function(object, assay = "SPM", slot = "counts",image = "slice1", nfeatures = 2000, verbose = TRUE){

  features <- rownames(x = object[[assay]])
  spatial.location <- GetTissueCoordinates(object = object[[image]])
  data <- GetAssayData(object = object, slot = slot)
  data <- as.matrix(x = data[features, ])
  data <- data[RowVar(x = data) > 0, ]
  svf.info <- RunMoransI(data = data, pos = spatial.location, verbose = verbose)
  colnames(x = svf.info) <- paste0("MoransI_", colnames(x = svf.info))
  var.name <- paste0("moransi", ".spatially.variable")
  var.name.rank <- paste0(var.name, ".rank")
  svf.info[[var.name]] <- FALSE
  svf.info <- svf.info[order(svf.info[, 2], -abs(svf.info[, 1])), , drop = FALSE]
  svf.info[[var.name]][1:(min(nrow(x = svf.info), nfeatures))] <- TRUE
  svf.info[[var.name.rank]] <- 1:nrow(x = svf.info)
  object[[assay]][[names(x = svf.info)]] <- svf.info

  return(object)
}

#' Gets the top n number of spatially variable features
#'
#' @param object SpaMTP Seurat class object containing the intensity values for each m/z
#' @param assay Character string indicating which Seurat object assay to pull data form (default = "SPM").
#' @param n Numeric value defining the top number of metabolites to return (default = 10).
#'
#' @return Vector containing m/z feature names corresponding to the top n number of spatially variable metabolites
#' @export
#'
#' @examples
#' # features <- GetSpatiallyVariableMetabolites(SpaMTP, n = 6)
GetSpatiallyVariableMetabolites <- function(object, assay = "SPM", n = 10){

  return(rownames(object[[assay]][["moransi.spatially.variable.rank"]]%>%arrange(moransi.spatially.variable.rank))[1:n])
}
