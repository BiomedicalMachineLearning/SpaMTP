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
#' split_pixel(center_point = c(0, 0), pixel_radius = 1, pseudo_n = 4, pixel_shape = "circle", show_split_diagram = TRUE)
#'
#' # Split a circular pixel into 9 spots and display the split diagram
#' split_pixel(center_point = c(0, 0), pixel_radius = 1, pseudo_n = 9, pixel_shape = "square", show_split_diagram = TRUE)
#'
#' # Split a square pixel into 16 spots and display the split diagram
#' split_pixel(center_point = c(0, 0), pixel_radius = 1, pseudo_n = 16, pixel_shape = "square", show_split_diagram = TRUE)
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
#' @param SM.data A Seurat Spatial Metabolomics object containing spatial metabolomic information. It should have 'x_coord' and 'y_coord' columns representing spot coordinates.
#' @param res_increase An integer specifying the factor by which the resolution should be increased. res_increase = 4, will generates 4 times as many spots in each dimension (default = 4).
#'
#' @returns A data frame with increased resolution metadata, containing additional columns 'new_x_coord', 'new_y_coord', 'new_barcode', and 'old_barcode'. The 'new_x_coord' and 'new_y_coord' columns represent the coordinates of the pseudo-high-resolution spots. The 'new_barcode' column is a unique identifier for each pseudo-high-resolution spot, and 'old_barcode' retains the original spot identifier.
#' @export
#'
#' @examples
#' ## Increase resolution of MALDI dataset by a factor of 4
#' # increased_data <- increase_MALDI_res(SeuratObj, res_increase = 4)
increase_MALDI_res <- function(SM.data, res_increase = 4) {
  message("Increasing the resolution of MALDI Pixel Data ...\n")

  # Assuming SM.data$obs has columns 'x_coord' and 'y_coord'

  # Randomly sample 1000 indices
  sampled_indices <- sample(1:(nrow(SM.data@meta.data) - 1), size = 1000, replace = TRUE)

  # Initialize a list to store the distances
  distances <- c()

  # Calculate distances for each pair
  for (i in sampled_indices) {
    point1 <- c(SM.data@meta.data$x_coord[i], SM.data@meta.data$y_coord[i])
    point2 <- c(SM.data@meta.data$x_coord[i + 1], SM.data@meta.data$y_coord[i + 1])
    distances <- c(distances, sqrt(sum((point1 - point2)^2)) / 2)
  }

  # Find the median distance
  median_distance <- stats::median(distances)

  message("Median Distance Between MALDI Spots: ", median_distance, "\n")

  # Generates new obs matrix with higher resolution
  new_meta_data <- data.frame(matrix(NA, nrow = nrow(SM.data@meta.data) * 4, ncol = ncol(SM.data@meta.data)))
  colnames(new_meta_data) <- colnames(SM.data@meta.data)

  total_spots = nrow(SM.data@meta.data)

  message("Generating psuedo-highres MALDI data: ")

  pb <- utils::txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = total_spots, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")


  for (spot_idx in 1:total_spots) {

    sub_spot_idx <- (spot_idx * 4)-3
    pseudo_centers <- split_pixel(c(SM.data@meta.data$x_coord[spot_idx], SM.data@meta.data$y_coord[spot_idx]), median_distance, pseudo_n = 4, pixel_shape = "square", show_split_diagram = FALSE)

    for (idx in seq(0, (dim(pseudo_centers)[1])-1)) {
      new_meta_data[sub_spot_idx + idx, ] <- SM.data@meta.data[spot_idx, ]
      new_meta_data[sub_spot_idx + idx, "old_barcode"] <- rownames(SM.data@meta.data[spot_idx,])
      new_meta_data[sub_spot_idx + idx, c("new_x_coord", "new_y_coord")] <- pseudo_centers[idx+1, ]
      new_meta_data[sub_spot_idx + idx, "new_barcode"] <- paste0(pseudo_centers[idx+1, 1], "_", pseudo_centers[idx+1, 2])
    }
    utils::setTxtProgressBar(pb, spot_idx)
  }
  close(pb)
  return(new_meta_data)
}



#' Generate new MALDI counts matrix for equivalent Visium spots
#'    - This function takes an original Spatial Metabolomics counts matrix, along with a metadata table (`obs_x`) containing information about the correspondence between MALDI and Visium spots. It generates a new counts matrix where each MALDI spot has an associated aggregated count based on its corresponding Visium spots.
#'
#' @param original_MALDI A Seurat Spatial Metabolomics object containing the original counts matrix
#' @param obs_x A metadata table with information about the correspondence between MALDI and Visium spots. It should have columns 'Visium_spot' and 'MALDI_barcodes'.
#' @param assay Character string defining the Seurat assay that contains the annotated counts and metadata corresponding to the m/z values.
#' @param slots Vector of character strings describing which slots to pull the relative intensity values from (default = c("counts", "data")).
#'
#' @return A data frame representing the new counts matrix for equivalent Visium spots,where each row corresponds to a Visium spot and columns correspond to m/z features.
#' @export
#'
#' @examples
#' ## Generate new MALDI counts matrix for equivalent Visium spots
#' # new_counts <- generate_new_MALDI_counts(SeuratObj, obs_x, assay = "Spatial")
generate_new_MALDI_counts <- function(original_MALDI, obs_x, assay, slots) {

  message("Merging MALDI counts ... ")

  data_list <- list()
  new_data_list <- list()

  for (slot in slots) {
    data_list[[slot]] <- t(as.data.frame(original_MALDI[[assay]][slot]))
    new_data_list[[slot]] <- data.frame(matrix(NA, nrow = 0, ncol = ncol(data_list[[slot]])))
    colnames(new_data_list[[slot]]) <- colnames(data_list[[slot]])
  }


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
      for (data_slot in names(data_list)){
        new_data_list[[data_slot]][spot_barcode, ] <- data_list[[data_slot]][barcode_list[1], ]
      }
    } else {
      for (data_slot in names(data_list)){
        mini_raw_count_dfs <- lapply(barcode_list, function(MALDI_barcodes) data_list[[data_slot]][MALDI_barcodes, ])
        combined_df <- do.call(rbind, mini_raw_count_dfs)
        new_data_list[[data_slot]][spot_barcode, ] <- colMeans(combined_df)
      }
    }
    utils::setTxtProgressBar(pb, spot_idx)
  }
  close(pb)

  return(new_data_list)
}



#' Converts and aggregates Spatial Metabolomic (MALDI) data to corresponding Spatial Transcriptomics (Visium) spots.
#'    - This function uses generate_new_MALDI_counts() and increase_MALDI_res()
#'
#' @param SM.data A Seurat object representing the Spatial Metabolomics data.
#' @param ST.data A Seurat object representing the Visium Spatial Transcriptomics data.
#' @param img_res Character string defining the image resolution associated with the Visium image pixel data (default = "hires").
#' @param res_increase Integer value defining the factor by which the resolution of MALDI spots should be increased before assignment. It should be either 4 or 9, see increase_MALDI_res() documentation for specifics (Default = NULL).
#' @param annotations Boolean value indicating if the Spatial Metabolomics (MALDI) Seurat object contains annotations assigned to m/z values (default = FALSE).
#' @param assay Character string defining the Seurat assay that contains the annotated counts and metadata corresponding to the m/z values (default = "Spatial").
#' @param slice Character string of the image slice name in the Visium object (default = "slice1").
#' @param slots Vector of character strings describing which slots from the Spatial Metabolomic Seurat object to adjust for the new overlayed object (default = c("counts", "data")).
#' @param new_SpM.assay Character string defining the assay name of the new overlayed Seurat object containing all updated metabolomic data (default = "SPM").
#' @param add.ST Boolean indicating whether to add the Spatial Metabolomic data from the Seurat Visium Object to the new updated Seruat Object (default = TRUE).
#' @param ST.assay Character string specifying the current assay to use to extract transcriptional data from (default = "Spatial").
#' @param ST.layers Vector of character strings defining the relative slots to extract from the initial Visium object, to add to the new Seurat Object if ST.assay == TRUE (default = c("counts", "data")).
#' @param new_SpT.assay Character string defining the assay name of the new overlayed Seurat object containing all updated transcriptomics data (default = "SPT").
#' @param verbose Boolean value indicating whether to print pregression update messages and progress bar (default = TRUE).
#'
#' @return A Seurat object with the Spatial Metabolomic data assigned to equivalent Spatial Transcripomics (Visium) spots.
#' @export
#'
#' @examples
#' ## Convert MALDI data to equivalent Visium spots
#' # convert_MALDI_to_visium_like_adata(VisiumObj, SeuratObj, img_res = "hires", new_library_id = "MALDI", res_increase = NULL)
AlignSpatialOmics <- function(SM.data, ST.data, res_increase = NULL, annotations = FALSE, assay = "Spatial", slots = c("counts"), img_res = "hires", slice = "slice1", new_SpM.assay = "SPM", add.ST = TRUE, ST.assay = "Spatial", ST.layers = c("counts"), new_SpT.assay = "SPT", verbose = FALSE) {


  new_MALDI_obs <- SM.data@meta.data

  if (!is.null(res_increase)) {
    if (res_increase == 4 || res_increase == 9) {
      new_MALDI_obs <- increase_MALDI_res(SM.data, res_increase = 4)
    } else {
      message("Error: res_increase must be either 4 or 9\n")
      return(NULL)
    }
  }

  message("Assigning MALDI to Visium Spots ... \n")

  img <- ST.data@images[[slice]]

  image_data <- ST.data@images[[slice]]@coordinates
  image_data$imagerow_sf <- image_data$imagerow * ST.data@images[[slice]]@scale.factors[[img_res]]
  image_data$imagecol_sf <- image_data$imagecol * ST.data@images[[slice]]@scale.factors[[img_res]]

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
    obs_[v_points[i,], "Visium_spot"] <- rownames(ST.data@meta.data)[i]
  }

  obs_x <- obs_[, c("Visium_spot", "nFeature_Spatial", "old_barcode")]

  obs_x <- obs_x %>%
    dplyr::group_by(Visium_spot) %>%
    dplyr::summarize(MALDI_barcodes = toString(unique(old_barcode)))

  obs_x <- data.frame(obs_x)
  counts_x <- generate_new_MALDI_counts(SM.data, obs_x, assay = assay, slots = slots)

  message("Generating new MALDI Seurat Object ... ")
  seuratobj <- Seurat::CreateSeuratObject(counts = t(counts_x[[names(counts_x)[1]]]), assay = new_SpM.assay) #creates a new assay with the spatial metabolomics

  if (names(counts_x)[1] != "counts"){
    seuratobj[[new_SpM.assay]][names(counts_x)[1]] <- seuratobj[[new_SpM.assay]]
    seuratobj[[new_SpM.assay]]$counts <- NULL
  }

  if (length(names(counts_x)) > 1){
    for (slot_name in names(counts_x)[2:length(names(counts_x))]){
      seuratobj[[new_SpM.assay]][slot_name] <- t(counts_x[[slot_name]])
    }
  }


  rownames(obs_x) <- rownames(seuratobj[[new_SpM.assay]]@cells)
  seuratobj@meta.data[colnames(obs_x)] <- obs_x



  vis_subset <- subset(ST.data, cells = rownames(seuratobj@meta.data))
  seuratobj <- subset(seuratobj, cells = rownames(vis_subset@meta.data))
  seuratobj@images[["slice1"]] <- vis_subset@images[[slice]]



  if (annotations){
    seuratobj[[new_SpM.assay]]@meta.data <- SM.data[[assay]]@meta.data
  }

  if (add.ST){
    visium_sub <- subset(ST.data, cells = intersect(colnames(seuratobj), colnames(ST.data)))

    # add Transcriptomic assay to combined object
    counts_layer <- NULL
    data_layer <- NULL

    for (layer in ST.layers) {
      if ( layer == "counts" ) {
        counts_layer <- visium_sub[[ST.assay]]$counts
      }  else if (layer == "data") {
        data_layer <- visium_sub[[ST.assay]]$data
      } else {
        stop("ST.layers given does not match requirments. Must be 'counts' and/or 'data'")
      }
    }
    st_assay <- CreateAssay5Object(counts = counts_layer, data = data_layer)
    seuratobj[[new_SpT.assay]] <- st_assay


  }
  return(seuratobj)

}


########################################################################################################################################################################################################################


