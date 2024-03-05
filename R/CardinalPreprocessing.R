library(Cardinal)


#### SpaMTP Cardinal Preprocessing Functions ############################################################################################################################################################################



#' Calculates new coordinates for merged Cardinal Objects
#'    - When merging objects this will shift the plotting coordinates appropriately to ensure that images do not overlap when displayed on the same plot
#'
#' @param index An integer defining the index of the Cardinal Object in the list of Cardinal Objects to merge.
#' @param ncol An integer defining the number of columns present when plotting the merged Cardinal Object.
#' @param padding An integer defining the pixel padding between plots
#'
#' @returns A list containing new x and y coordinate values
#' @export
#'
#' @examples
#' calculate_coordinates(index = 1, ncol = 2, padding = 200)
calculate_coordinates <- function(index, ncol, padding) {
  row <- ceiling(index / ncol)
  col <- (index - 1) %% ncol + 1
  x <- (col - 1) * padding
  y <- (row - 1) * padding
  return(list(x = x, y = y))
}


#' Sets the centroid value of a Cardinal Object = TRUE
#'    - This is required to merge Cardinal objects together
#'
#' @param data A Cardinal Object required for merging.
#'
#' @returns A Cardinal Object with centroid = TRUE
#' @export
#'
#' @examples
#' # set_centroid_to_true(cardinalObj)
set_centroid_to_true <- function(data){
  Cardinal::centroided(data) <- TRUE
  return(data)
}



#' Merges Cardinal Objects together
#'    - This is required for combined analysis and plotting
#'
#' @param data_list List of Cardinal Objects being merged together.
#' @param shift.image Boolean value describing if to shift the image coordinates for plotting to prevent overlay (default = FALSE).
#' @param ncols An integer defining the number of columns present when plotting the merged Cardinal Object (default = 2).
#' @param padding An integer defining the pixel padding between plots (default = 200).
#'
#' @returns A Cardinal Object with values merged form each individual Cardinal Object givin.
#' @export
#'
#' @examples
#' # cardinal_object_list <- list(cardinalObj1, cardinalObj2, cardinalObj3, cardinalObj4)
#' # MergeCardinalData(cardinal_object_list)
MergeCardinalData <- function(data_list, shift.image = FALSE, ncols = 2, padding = 200){

  ### NOTE: mass.range and resolution of each sample in the data list must be the same to merge

  message("Setting Centroids as TRUE ...... ")

  data_list <- lapply(data_list, set_centroid_to_true)

  if (shift.image){
    ### Changes pixel coordinates so that they do not overlap on the same plot
    message("Shifting Pixel Coordinates to Combining Samples ...... ")

    n <- 1
    while ( n <= length(data_list)){

      data_copy <- data_list[[n]]
      data_copy_pixel_data <- Cardinal::pixelData(data_copy)
      coordnate_list <- calculate_coordinates(n,ncols, padding)
      data_copy_pixel_data@coord$x <- data_copy_pixel_data@coord$x + as.numeric(coordnate_list$x) # Adds values to pixels so they dont overlap on the plots
      data_copy_pixel_data@coord$x <- data_copy_pixel_data@coord$y + as.numeric(coordnate_list$y)
      Cardinal::pixelData(data_list[[n]]) <- data_copy_pixel_data
      n <- n+1
    }
  }

  return(Cardinal::combine(data_list))

}


########################################################################################################################################################################################################################
