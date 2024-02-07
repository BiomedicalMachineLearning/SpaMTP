library(Seurat)


#' Normalizes m/z intensity data stored in a Seurat Object
#'
#' @param data Seurat Object to be normalized.
#' @param normalisation.type Character string defining the normalization method to run. Options are either c("TIC", "LogNormalize") which represent Total Ion Current (TIC) normalization or Log Normalization, respectively (default = "TIC").
#' @param scale.factor Numeric value that sets the scale factor for pixel/spot level normalization. Following normalization the total intensity value across each pixel will equal this value. If scale.factor = NULL, TIC normalization will use a scale factor = number of m/z and Log Normalisation will use a scale factor = 10000 (default = NULL).
#' @param assay Character string defining the name of the Seurat Object assay to pull the corresponding intensity data from (default = "Spatial").
#' @param slot  Character string defining the name of the slot within the Seurat Object assay to pull the corresponding intensity data from (default = "counts").
#'
#' @return A Seurat Object with intensity values normalized. Normalized data is stored in the $data slot of the specified assay
#' @export
#'
#' @examples
#' # normalised_data <- NormalizeSeuratData(SeuratObject)
NormalizeSeuratData <- function(data, normalisation.type = 'TIC', scale.factor = NULL, assay = "Spatial", slot = "counts") {

  if (is.null(normalisation.type)) {
    stop("Error: no normalisation.type is select. Please enter either 'LogNormalize' or 'TIC'")
  }

  if (!(normalisation.type == "TIC" | normalisation.type == "LogNormalize")) {
    stop("Error: incorrect normalisation.type is select. Please enter either 'LogNormalize' or 'TIC'")
  }

  if (is.null(scale.factor)) {
    scale.factor <- 10000
  }

  if (normalisation.type == 'TIC') {
    scale.factor <- length(rownames(data))
  }

  normalised.data <- Seurat::NormalizeData(data, scale.factor = scale.factor)


  if (normalisation.type == 'TIC') {
    normalised.data[[assay]]$data <- exp(normalised.data[[assay]][slot]) -1
    normalised.data[[assay]]$data <- as(normalised.data[[assay]]$data, "sparseMatrix")
  }


  return(normalised.data)

}
