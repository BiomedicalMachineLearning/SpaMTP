library(Seurat)
library(RANN)
library(stats)
library(utils)
library(graphics)


#### SpaMTP MALDI to Visium Spot Merging Functions  ####################################################################################################################################################################


#' Splits a pixel into multiple spots
#'     - This function takes a central point and a pixel radius and splits the pixel into multiple smaller spots based on the specified configuration.
#'
#' @param center_point A numeric vector of length 2 representing the (x, y) coordinates of the center point of the pixel.
#' @param pixel_radius The radius of the original pixel.
#' @param pseudo_n The number of smaller spots to split the pixel into. Must be either 4 or 9 for square pixels and 4 for a circular pixel.
#' @param pixel_shape A character string indicating the shape of the original pixel. Options are "square" (default) or "circle".
#' @param show_split_diagram A logical indicating whether to plot a diagram showing the original pixel and the split spots (default = FALSE).
#'
#' @returns A matrix with two columns representing the (x, y) coordinates of the centers of the smaller spots after splitting the pixel.
#' @export
#'
#' @examples
#' # Split a square pixel into 4 spots
#' split_pixel(center_point = c(0, 0), pixel_radius = 1, pseudo_n = 4)
#'
#' # Split a circular pixel into 9 spots and display the split diagram
#' split_pixel(center_point = c(0, 0), pixel_radius = 1, pseudo_n = 9, pixel_shape = "circle", show_split_diagram = TRUE)
#'
#' # Split a square pixel into 16 spots and display the split diagram
#' split_pixel(center_point = c(0, 0), pixel_radius = 1, pseudo_n = 16, show_split_diagram = TRUE)
split_pixel <- function(center_point, pixel_radius, pseudo_n = 4, pixel_shape = "square", show_split_diagram = FALSE) {

  if (pseudo_n == 4 || ((pseudo_n == 9 || pseudo_n == 16) && pixel_shape == "square")) {

    x <- center_point[1]
    y <- center_point[2]

    if (pixel_shape == "circle") {

      R <- pixel_radius
      # Radius of each smaller circle
      r <- R * (sqrt(pseudo_n) - 1)
      dis <- r * sqrt(pseudo_n)
      # Coordinates of the centers of the smaller circles
      centers <- matrix(c(R - dis + x, R - dis + y, R - dis + x, -(R - dis) + y,
                          -(R - dis) + x, R - dis + y, -(R - dis) + x, -(R - dis) + y), ncol = 2, byrow = TRUE)

    } else {
      # print("Assuming Pixel is in the shape of a square, if not set 'pixel_shape = 'circle'")
      r <- (pixel_radius * 2) / (sqrt(pseudo_n) * 2)

      if (pseudo_n == 4) {
        centers <- matrix(c(x - r, y - r, x + r, y - r, x + r, y + r, x - r, y + r), ncol = 2, byrow = TRUE)
      } else if (pseudo_n == 9) {
        centers <- matrix(c(x - r * 0, y + r * 0, x - r * 0, y - r * 2, x - r * 2, y + r * 0,
                            x - r * 2, y + r * 2, x - r * 2, y - r * 2,
                            x + r * 2, y + r * 0, x + r * 0, y + r * 2,
                            x + r * 2, y - r * 2, x + r * 2, y + r * 2), ncol = 2, byrow = TRUE)
      } else {
        centers <- matrix(c(x - r * 1, y + r * 1, x + r * 1, y - r * 1, x + r * 1, y + r * 1,
                            x - r * 1, y - r * 1, x - r * 1, y - r * 3,
                            x - r * 1, y + r * 3, x - r * 3, y + r * 1,
                            x - r * 3, y - r * 1, x - r * 3, y - r * 3,
                            x - r * 3, y + r * 3, x + r * 1, y + r * 3,
                            x + r * 1, y - r * 3, x + r * 3, y + r * 1,
                            x + r * 3, y - r * 1, x + r * 3, y - r * 3,
                            x + r * 3, y + r * 3), ncol = 2, byrow = TRUE)
      }
    }



    if (show_split_diagram) {
      if (pixel_shape == "circle") {
        plot(c(x - pixel_radius * 2, x + pixel_radius * 2), c(y - pixel_radius * 2, y + pixel_radius * 2), type = "n", xlab = "X-axis", ylab = "Y-axis")
        graphics::symbols(centers[, 1], centers[, 2], circles = rep(r, pseudo_n), add = TRUE, inches = 1, col = "red", bg = "white")
        graphics::points(x, y, pch = 19, col = "blue")
        graphics::title(main = "Four Circles Inside a Larger Circle")

      } else {
        graphics::plot(c(x - pixel_radius * 2, x + pixel_radius * 2), c(y - pixel_radius * 2, y + pixel_radius * 2), type = "n", xlab = "X-axis", ylab = "Y-axis")
        graphics::rect(x - pixel_radius, y - pixel_radius, x + pixel_radius, y + pixel_radius, border = "blue", lty = 1, col = NA)
        graphics::symbols(centers[, 1], centers[, 2], circles = rep(r, pseudo_n), add = TRUE, inches = 1, col = "red", bg = "white")
        graphics::title(main = "Four Circles Inside a Square")
      }
    }

    return(centers)

  } else {
    warning("ERROR: pseudo_n must be either 4 or 9 for square pixels and 4 for circle")
    stop("n invalid integer")
  }
}



#' Increase the resolution of Spatial Metabolomics data for merging with Spatial Transcriptomic (Visium) data
#'    - This function takes a Spatial Metabolomics dataset and increases its resolution by generating pseudo-high-resolution spots. This is useful for merging the data with Spatial Transcriptomic (Visium) data that has a different resolution.
#'    - This function uses the split_pixel() function
#'
#' @param MALDI_adata A Seurat Spatial Metabolomics object containing spatial metabolomic information. It should have 'x_coord' and 'y_coord' columns representing spot coordinates.
#' @param res_increase An integer specifying the factor by which the resolution should be increased. res_increase = 4, will generates 4 times as many spots in each dimension (default = 4).
#'
#' @returns A data frame with increased resolution metadata, containing additional columns 'new_x_coord', 'new_y_coord', 'new_barcode', and 'old_barcode'. The 'new_x_coord' and 'new_y_coord' columns represent the coordinates of the pseudo-high-resolution spots. The 'new_barcode' column is a unique identifier for each pseudo-high-resolution spot, and 'old_barcode' retains the original spot identifier.
#' @export
#'
#' @examples
#' ## Increase resolution of MALDI dataset by a factor of 4
#' # increased_data <- increase_MALDI_res(SeuratObj, res_increase = 4)
increase_MALDI_res <- function(MALDI_adata, res_increase = 4) {
  message("Increasing the resolution of MALDI Pixel Data ...\n")

  # Assuming MALDI_adata$obs has columns 'x_coord' and 'y_coord'

  # Randomly sample 1000 indices
  sampled_indices <- sample(1:(nrow(MALDI_adata@meta.data) - 1), size = 1000, replace = TRUE)

  # Initialize a list to store the distances
  distances <- c()

  # Calculate distances for each pair
  for (i in sampled_indices) {
    point1 <- c(MALDI_adata@meta.data$x_coord[i], MALDI_adata@meta.data$y_coord[i])
    point2 <- c(MALDI_adata@meta.data$x_coord[i + 1], MALDI_adata@meta.data$y_coord[i + 1])
    distances <- c(distances, sqrt(sum((point1 - point2)^2)) / 2)
  }

  # Find the median distance
  median_distance <- stats::median(distances)

  message("Median Distance Between MALDI Spots: ", median_distance, "\n")

  # Generates new obs matrix with higher resolution
  new_meta_data <- data.frame(matrix(NA, nrow = nrow(MALDI_adata@meta.data) * 4, ncol = ncol(MALDI_adata@meta.data)))
  colnames(new_meta_data) <- colnames(MALDI_adata@meta.data)

  total_spots = nrow(MALDI_adata@meta.data)

  message("Generating psuedo-highres MALDI data: ")

  pb <- utils::txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = total_spots, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")


  for (spot_idx in 1:total_spots) {

    sub_spot_idx <- (spot_idx * 4)-3
    pseudo_centers <- split_pixel(c(MALDI_adata@meta.data$x_coord[spot_idx], MALDI_adata@meta.data$y_coord[spot_idx]), median_distance, pseudo_n = 4, pixel_shape = "square", show_split_diagram = FALSE)

    for (idx in seq(0, (dim(pseudo_centers)[1])-1)) {
      new_meta_data[sub_spot_idx + idx, ] <- MALDI_adata@meta.data[spot_idx, ]
      new_meta_data[sub_spot_idx + idx, "old_barcode"] <- rownames(MALDI_adata@meta.data[spot_idx,])
      new_meta_data[sub_spot_idx + idx, c("new_x_coord", "new_y_coord")] <- pseudo_centers[idx+1, ]
      new_meta_data[sub_spot_idx + idx, "new_barcode"] <- paste0(pseudo_centers[idx+1, 1], "_", pseudo_centers[idx+1, 2])
    }
    utils::setTxtProgressBar(pb, spot_idx)
  }
  utils::close(pb)
  return(new_meta_data)
}



#' Generate new MALDI counts matrix for equivalent Visium spots
#'    - This function takes an original Spatial Metabolomics counts matrix, along with a metadata table (`obs_x`) containing information about the correspondence between MALDI and Visium spots. It generates a new counts matrix where each MALDI spot has an associated aggregated count based on its corresponding Visium spots.
#'
#' @param original_MALDI A Seurat Spatial Metabolomics object containing the original counts matrix
#' @param obs_x A metadata table with information about the correspondence between MALDI and Visium spots. It should have columns 'Visium_spot' and 'MALDI_barcodes'.
#' @param assay Character string defining the Seurat assay that contains the annotated counts and metadata corresponding to the m/z values.
#'
#' @return A data frame representing the new counts matrix for equivalent Visium spots,where each row corresponds to a Visium spot and columns correspond to m/z features.
#' @export
#'
#' @examples
#' ## Generate new MALDI counts matrix for equivalent Visium spots
#' # new_counts <- generate_new_MALDI_counts(SeuratObj, obs_x, assay = "Spatial")
generate_new_MALDI_counts <- function(original_MALDI, obs_x, assay) {

  message("Merging MALDI counts ... ")

  counts_df <- t(as.data.frame(original_MALDI[[assay]]$counts))

  new_count_data <- data.frame(matrix(NA, nrow = 0, ncol = ncol(counts_df)))
  colnames(new_count_data) <- colnames(counts_df)


  total_spots = nrow(obs_x)
  pb <- utils::txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = total_spots, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")

  for (spot_idx in 1:total_spots) {
    spot_barcode <- obs_x$Visium_spot[spot_idx]
    barcode_list <- unlist(strsplit(obs_x$MALDI_barcodes[spot_idx], ", "))

    if (all(barcode_list[1] == barcode_list)) {
      new_count_data[spot_barcode, ] <- counts_df[barcode_list[1], ]
    } else {
      mini_raw_count_dfs <- lapply(barcode_list, function(MALDI_barcodes) counts_df[MALDI_barcodes, ])
      combined_df <- do.call(rbind, mini_raw_count_dfs)
      new_count_data[spot_barcode, ] <- colMeans(combined_df)
    }
    utils::setTxtProgressBar(pb, spot_idx)
  }
  utils::close(pb)

  return(new_count_data)
}


#' Converts and aggregates Spatial Metabolomic (MALDI) data to corresponding Spatial Transcriptomics (Visium) spots.
#'    - This function uses generate_new_MALDI_counts() and increase_MALDI_res()
#'
#' @param visium_adata A Seurat object representing the Visium Spatial Transcriptomics data.
#' @param MALDI_adata A Seurat object representing the Spatial Metabolomics data.
#' @param img_res Character string defining the image resolution associated with the Visium image pixel data (default = "hires").
#' @param new_library_id Character string specifying the library ID for the new MALDI data in the output Seurat object.
#' @param res_increase Integer value defining the factor by which the resolution of MALDI spots should be increased before assignment. It should be either 4 or 9, see increase_MALDI_res() documentation for specifics (Default = NULL).
#' @param annotations Boolean value indicating if the Spatial Metabolomics (MALDI) Seurat object contains annotations assigned to m/z values (default = FALSE).
#' @param assay Character string defining the Seurat assay that contains the annotated counts and metadata corresponding to the m/z values (default = "Spatial").
#' @param slice Character string of the image slice name in the Visium object (default = "slice1").
#'
#' @return A Seurat object with the Spatial Metabolomic data assigned to equivalent Spatial Transcripomics (Visium) spots.
#' @export
#'
#' @examples
#' ## Convert MALDI data to equivalent Visium spots
#' # convert_MALDI_to_visium_like_adata(VisiumObj, SeuratObj, img_res = "hires", new_library_id = "MALDI", res_increase = NULL)
convert_MALDI_to_visium_like_adata <- function(visium_adata, MALDI_adata, img_res = "hires", new_library_id = "MALDI", res_increase = NULL, annotations = FALSE, assay = "Spatial", slice = "slice1") {

  new_MALDI_obs <- MALDI_adata@meta.data

  if (!is.null(res_increase)) {
    if (res_increase == 4 || res_increase == 9) {
      new_MALDI_obs <- increase_MALDI_res(MALDI_adata, res_increase = 4)
    } else {
      message("Error: res_increase must be either 4 or 9\n")
      return(NULL)
    }
  }

  message("Assigning MALDI to Visium Spots ... \n")

  img <- visium_adata@images[[slice]]

  image_data <- visium_adata@images[[slice]]@coordinates
  image_data$imagerow_sf <- image_data$imagerow * visium_adata@images[[slice]]@scale.factors[[img_res]]
  image_data$imagecol_sf <- image_data$imagecol * visium_adata@images[[slice]]@scale.factors[[img_res]]

  dis <- abs((stats::lm(image_data$imagerow_sf ~image_data$col))$coefficients[[2]])
  radius <- 2 * dis / 100 * 55 / 2

  new_coords <- as.matrix(new_MALDI_obs[, c("new_y_coord", "new_x_coord")])
  query_coords <- as.matrix(image_data[, c("imagerow_sf", "imagecol_sf")])
  v_points <- RANN::nn2(new_coords,query_coords, treetype = "kd",searchtype = "radius", radius = radius)$nn.idx

  obs_ <- data.frame(
    index = rownames(new_MALDI_obs),
    nFeature_Spatial = new_MALDI_obs$nFeature_Spatial,
    old_barcode = new_MALDI_obs$old_barcode,
    Visium_spot = "Not_assigned")

  for (i in 1:nrow(v_points)) {
    obs_[v_points[i,], "Visium_spot"] <- rownames(visium_adata@meta.data)[i]
  }

  obs_x <- obs_[, c("Visium_spot", "nFeature_Spatial", "old_barcode")]

  obs_x <- obs_x %>%
    dplyr::group_by(Visium_spot) %>%
    dplyr::summarize(MALDI_barcodes = toString(unique(old_barcode)))

  obs_x <- data.frame(obs_x)
  counts_x <- generate_new_MALDI_counts(MALDI_adata, obs_x, assay = assay)

  message("Generating new MALDI Anndata Object ... ")

  seuratobj <- Seurat::CreateSeuratObject(counts = t(counts_x), assay = "Spatial")

  rownames(obs_x) <- rownames(seuratobj@assays$Spatial@cells)
  seuratobj@meta.data <- obs_x

  vis_subset <- subset(visium_adata, cells = rownames(seuratobj@meta.data))
  seuratobj <- subset(seuratobj, cells = rownames(vis_subset@meta.data))
  seuratobj@images[["slice1"]] <- vis_subset@images$slice1

  if (annotations){
    seuratobj[[assay]]@meta.data <- MALDI_adata[[assay]]@meta.data
  }

  return(seuratobj)

}


########################################################################################################################################################################################################################


