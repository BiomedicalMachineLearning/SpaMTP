library(Cardinal)
library(Seurat)


#### SpaMTP Cardinal to Seurat Functions ###############################################################################################################################################################################

#' Converts a Cardinal Object into a Seurat Object
#'
#' @param data A Cardinal Object that is being converted into a Seurat Object.
#' @param run_name A character string defining the run name of the Cardinal data to be converted to a Seurat Object
#' @param seurat.coord A Data.Frame containing two columns titled 'X_new' and 'Y_new' specifying the pixel coordinates of each data point. This is only required mapping Spatial Metabolomic data with a H&E image if the Cardinal Object coordinates are not the same as the H&E image coordinates. Else, set to NULL(default).
#'
#' @returns A Seurat Object containing the mz count data of the supplied Cardinal Object
#' @export
#'
#' @examples
#' # CardinalToSeurat(CardinalObj, run_name = "run_1", seurat.coord = NULL)
CardinalToSeurat <- function(data,run_name, seurat.coord = NULL){

  message("Convering Cardinal object to Seurat object .... ")
  run_data <- Cardinal::subsetPixels(data, Cardinal::run(data) == paste0(run_name))
  sparse_matrix <- Cardinal::spectra(run_data)
  pixel_data <- Cardinal::pixelData(run_data)

  if (!(is.null(seurat.coord))){
    message("Convering Cardinal Coordinates to Seurat Visium Coordinates specified in the seurat.coord file .... ")
    pixel_data[["x_coord",]] <- seurat.coord$X_new # changes coordinates to matched Visium Object
    pixel_data[["y_coord",]] <- seurat.coord$Y_new
    Cardinal::pixelData(run_data) <- pixel_data
  } else{
    if ("x" %in% colnames(data.frame(pixel_data)) & "y" %in% colnames(data.frame(pixel_data))){
      pixel_data_df <- data.frame(pixel_data)
      pixel_data[["x_coord",]] <- pixel_data_df$x
      pixel_data[["y_coord",]] <- pixel_data_df$y
      Cardinal::pixelData(run_data) <- pixel_data
    } else {
      warning("There is no column called 'x' and 'y' in pixelData(CardinalObject)")
      stop("x and y pixel columns do not exist")
    }
  }


  message("Generating Seurat Barcode Labels from Pixel Coordinates .... ")
  spot_name <- c()

  for(idx in seq(1,length(Cardinal::pixelData(run_data)[[1]]))){
    x_coord <- Cardinal::pixelData(run_data)[["x_coord",]][idx]
    y_coord <- Cardinal::pixelData(run_data)[["y_coord",]][idx]
    name <- paste0(x_coord,"_",y_coord)
    spot_name <- c(spot_name, name)
  }


  colnames(sparse_matrix)<- spot_name
  rownames(sparse_matrix)<- paste("mz-", Cardinal::featureData(run_data)@mz, sep = "")

  message("Constructing Seurat Object ....")
  mat <- as.matrix(sparse_matrix)


  seuratobj <- Seurat::CreateSeuratObject(mat, assay = "Spatial")

  message("Adding Pixel Metadata ....")
  seuratobj <- Seurat::AddMetaData(seuratobj,col.name = "sample", metadata = Cardinal::run(run_data))

  for (name in names(Cardinal::pixelData(run_data))){
    seuratobj <- Seurat::AddMetaData(seuratobj,col.name = name, metadata = Cardinal::pixelData(run_data)[[name,]])
  }

  message("Creating Centroids for Spatial Seurat Object ....")
  ## Add spatial data
  cents <- SeuratObject::CreateCentroids(data.frame(x = c(Cardinal::pixelData(run_data)[["x_coord",]]), y = c(Cardinal::pixelData(run_data)[["y_coord",]]), cell = c(spot_name)))


  segmentations.data <- list(
    "centroids" = cents,
    "segmentation" = NULL
  )

  coords <- SeuratObject::CreateFOV(
    coords = segmentations.data,
    type = c("segmentation", "centroids"),
    molecules = NULL,
    assay = "Spatial"
  )

  seuratobj[["fov"]] <- coords


  metadata <- data.frame("raw_mz" = sapply(strsplit(rownames(seuratobj), "-"), function(x) as.numeric(x[[2]])))
  rownames(metadata) <- rownames(seuratobj)


  seuratobj[["Spatial"]] <- Seurat::AddMetaData(object = seuratobj[["Spatial"]],
                                        metadata = metadata,
                                        col.name = 'raw_mz')


  seuratobj[["Spatial"]]@meta.data$mz_names <- rownames(seuratobj)

  return(seuratobj)
}
########################################################################################################################################################################################################################
