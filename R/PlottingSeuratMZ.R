library(Seurat)
library(ggplot2)


#### SpaMTP Seurat Plotting Functions #################################################################################################################################################################################

#' Finds the nearest m/z peak to a given value in the specified Seurat Object
#'
#' @param data Seurat Spatial Metabolomic object containing mz values
#' @param target_mz Numeric value defining the target m/z peak
#'
#' @returns String of the closest m/z value within the given dataset
#' @export
#'
#' @examples
#' # find_nearest(SeuratObj, target_mz = 400.01)
find_nearest <- function(data, target_mz){

  numbers <- as.numeric(gsub("mz-", "", SeuratObject::Features(data)))
  closest_number <- numbers[which.min(abs(numbers - target_mz))]
  return(paste0("mz-",closest_number))
}



#' Identifies all mz peaks within a plus-minus range of the target_mz
#'    - This function uses find_nearest()
#'
#' @param data Seurat Spatial Metabolomic object containing mz values
#' @param target_mz Numeric value defining the target m/z peak
#' @param plus_minus Numeric value defining the range/threshold either side of the target peak to also be identified
#'
#' @returns Vector of m/z values within plus-minus range of the target mz value
#'
#' @examples
#' # plusminus(SeuratObj, target_mz = 400.01, plus_minus = 0.005)
plusminus <- function(data, target_mz, plus_minus){
  feature_list <- c()
  center <- find_nearest(data, target_mz)

  center_value <- as.numeric(gsub("mz-", "", center))
  up_value <- center_value + as.numeric(plus_minus)
  low_value <- center_value - as.numeric(plus_minus)

  upper <- find_nearest(data, up_value)
  lower <- find_nearest(data, low_value)

  feature_list <- c(feature_list, center)
  if (!(upper %in% feature_list)){
    feature_list <- c(feature_list, upper)
  }
  if (!(lower %in% feature_list)){
    feature_list <- c(feature_list, lower)
  }

  if (length(feature_list) == 1){
    warning("No other m/z peaks in plusminus range -> increase to include more peaks")
  }
  return(feature_list)
}



#' Bins multiple m/z values into one.
#'
#' @param data Seurat Spatial Metabolomic object containing mz values
#' @param mz_list Vector of m/z names to be binned into one value
#'
#' @returns Vector of binned mz counts for each spot/pixel
#' @export
#'
#' @examples
#' # mz_values <- plusminus(SeuratObj, 448.2, 0.05)
#' # bin.mz(SeuratObj, mz_values)
bin.mz <- function(data, mz_list){
  data_copy <- data
  assay_counts <- data_copy@assays$Spatial$counts
  selected_genes <- assay_counts[mz_list, , drop = FALSE]

  if (is.null(selected_genes)) {
    stop("One or more genes not found in the assay counts.")
  }


  return(colSums(selected_genes))
}



#' Visualise expression of m/z values in a spatial context
#'      - This is for plotting Seurat Spatial Metabolomic data without an image(i.e. H&E image)
#'      - This function inherits off Seurat::ImageFeaturePlot(). Look here for more detailed documentation about inputs.
#'
#' @param object Seurat Spatial Metabolomic Object to Visualise.
#' @param mzs Vector of numeric m/z values to plot (e.g. c(400.1578, 300.1)). The function find_nearest() is used to automatically find the nearest m/z value to the ones given.
#' @param plusminus Numeric value defining the range/threshold either side of the target peak/peaks to be binned together for plotting (default = NULL).
#' @param fov Character string of name of FOV to plot (default = NULL).
#' @param boundaries A vector of segmentation boundaries per image to plot (default = NULL).
#' @param cols Vector of character strings defining colours used for plotting (default = c("lightgrey", "firebrick1")).
#' @param size Numeric value defining the point size for cells/spots when plotting (default = 0.5).
#' @param min.cutoff Vector of numeric value describing the minimum cutoff values for each m/z feature (default = NA).
#' @param max.cutoff Vector of numeric value describing the maximum cutoff values for each m/z feature (default = NA).
#' @param split.by Character string defining a factor in the Seurat Object metadata to split the feature plot by (default = NULL).
#' @param molecules Vector of character strings describing molecules to plot (default = NULL).
#' @param mols.size Numeric value for the point size of molecules to plot (default = 0.1).
#' @param mols.cols A vector of colours for the molecules to plot (default = NULL).
#' @param nmols Integer of the max number of each molecule specified in 'molecules' to be plot (default = 1000).
#' @param alpha Numeric value between 0 and 1 defining the image alpha (default = 1).
#' @param border.color Character string specifying the colour of each spot/cell border (default = white).
#' @param border.size Numeric value for the thickness of the cell segmentation border (default = NULL).
#' @param dark.background Boolean value indicating if the plot background is coloured black (default = FALSE).
#' @param blend Boolean value indicating whether to scale and blend expression values to visualize coexpression of two features (default = FALSE).
#' @param blend.threshold Numeric value defining the color cutoff from weak signal to strong signal; ranges from 0 to 1 (default = 0.5).
#' @param crop Boolean value of whether to crop the plots to area with cells only (default = FALSE).
#' @param cells Vector of character strings defining a group of cells to plot (default = NUll; plots all cells).
#' @param scale Set color scaling across multiple plots; c("features", "all", "none").
#' @param overlap Overlay boundaries from a single image to create a single plot (default = FALSE).
#' @param axes Boolean defining if to keep axes and panel background (default = FALSE).
#' @param combine Boolean value stating if to combine plots into a single patchworked ggplot object (default = TRUE).
#' @param coord.fixed Boolean value of it to plot cartesian coordinates with fixed aspect ratio (default = TURE).
#'
#' @returns A ggplot showing the Spatial representation of expression data of specified m/z values
#' @export
#'
#' @examples
#' # ImageMZPlot(SeuratObj, mzs = c(400.678, 300.1))
#' # ImageMZPlot(SeuratObj, mzs = c(400.678, 300.1), plusminus = 0.05)
ImageMZPlot <- function(object,
                        mzs,
                        plusminus = NULL,
                        fov = NULL,
                        boundaries = NULL,
                        cols = if (isTRUE(x = blend)) {
                          c("lightgrey", "#ff0000", "#00ff00")
                        } else {
                          c("lightgrey", "firebrick1")
                        },
                        size = 0.5,
                        min.cutoff = NA,
                        max.cutoff = NA,
                        split.by = NULL,
                        molecules = NULL,
                        mols.size = 0.1,
                        mols.cols = NULL,
                        nmols = 1000,
                        alpha = 1,
                        border.color = "white",
                        border.size = NULL,
                        dark.background = TRUE,
                        blend = FALSE,
                        blend.threshold = 0.5,
                        crop = FALSE,
                        cells = NULL,
                        scale = c("feature", "all", "none"),
                        overlap = FALSE,
                        axes = FALSE,
                        combine = TRUE,
                        coord.fixed = TRUE
){

  if (is.null(mzs)){
    stop("No mz values have been supplied")
  } else{

    mz_list <- c()
    for (target_mz in mzs){
      mz_string <- find_nearest(object, target_mz)
      mz_list <- c(mz_list, mz_string)
    }
  }

  if (!(is.null(plusminus))){

    data_copy <- object
    col_names_to_plot <- c()
    plot_titles <- c()

    for (target_mz in mz_list){
      mz_integer <- as.numeric(strsplit(target_mz, "-")[[1]][2])

      meta_col <- bin.mz(data_copy, plusminus(data_copy, mz_integer, plusminus))

      col_name <- paste0(target_mz,"_plusminus_", plusminus)
      plot_name <- paste0("mz: ", round(mz_integer, 3)," \\u00b1 ", plusminus)

      col_names_to_plot <- c(col_names_to_plot,col_name)
      plot_titles <- c(plot_titles,plot_name)
      data_copy[[col_name]] = meta_col
    }

    plot <- Seurat::ImageFeaturePlot(object = data_copy,
                             features = col_names_to_plot,
                             fov = fov,
                             boundaries = boundaries,
                             cols = cols,
                             size = size,
                             min.cutoff = min.cutoff,
                             max.cutoff = max.cutoff,
                             split.by = split.by,
                             molecules = molecules,
                             mols.size = mols.size,
                             mols.cols = mols.cols,
                             nmols = nmols,
                             alpha = alpha,
                             border.color = border.color,
                             border.size = border.size,
                             dark.background = dark.background,
                             crop = crop,
                             cells = cells,
                             scale = scale,
                             overlap = overlap,
                             axes = axes,
                             combine = combine,
                             coord.fixed = coord.fixed
    )

    for (plot_idx in seq(1, length(plot_titles))){
      plot[[plot_idx]] <- plot[[plot_idx]] +
        ggplot2::ggtitle(plot_titles[[plot_idx]])+
        ggplot2::labs(fill = plot_titles[[plot_idx]])
    }

  } else {

    plot <- Seurat::ImageFeaturePlot(object = object,
                             features = mz_list,
                             fov = fov,
                             boundaries = boundaries,
                             cols = cols,
                             size = size,
                             min.cutoff = min.cutoff,
                             max.cutoff = max.cutoff,
                             split.by = split.by,
                             molecules = molecules,
                             mols.size = mols.size,
                             mols.cols = mols.cols,
                             nmols = nmols,
                             alpha = alpha,
                             border.color = border.color,
                             border.size = border.size,
                             dark.background = dark.background,
                             crop = crop,
                             cells = cells,
                             scale = scale,
                             overlap = overlap,
                             axes = axes,
                             combine = combine,
                             coord.fixed = coord.fixed
    )

    for (plot_idx in seq(1, length(mz_list))){
      mz_integer <- as.numeric(strsplit(mz_list[plot_idx], "-")[[1]][2])
      plot[[plot_idx]] <- plot[[plot_idx]] +
        ggplot2::ggtitle(paste0("mz: ",round(mz_integer,3)))+
        ggplot2::labs(fill = paste0("mz: ",round(mz_integer,3)))
    }
  }


  return(plot)

}



#' Visualise expression of metabolites in a spatial context
#'      - This is for plotting Seurat Spatial Metabolomic data without an image(i.e. H&E image) that has been annotated using annotate.SeuratMALDI()
#'      - This function inherits off Seurat::ImageFeaturePlot(). Look here for more detailed documentation about inputs.
#'
#' @param object Seurat Spatial Metabolomic Object to Visualise.
#' @param metabolites Vector of metabolite names to plot (e.g. c("Glucose", "Glutamine")). The Seurat Object provided must contain annotations in the respective assay metadata.
#' @param plusminus Numeric value defining the range/threshold either side of the target peak/peaks to be binned together for plotting (default = NULL).
#' @param fov Character string of name of FOV to plot (default = NULL).
#' @param boundaries A vector of segmentation boundaries per image to plot (default = NULL).
#' @param cols Vector of character strings defining colours used for plotting (default = c("lightgrey", "firebrick1")).
#' @param size Numeric value defining the point size for cells/spots when plotting (default = 0.5).
#' @param min.cutoff Vector of numeric value describing the minimum cutoff values for each m/z feature (default = NA).
#' @param max.cutoff Vector of numeric value describing the maximum cutoff values for each m/z feature (default = NA).
#' @param split.by Character string defining a factor in the Seurat Object metadata to split the feature plot by (default = NULL).
#' @param molecules Vector of character strings describing molecules to plot (default = NULL).
#' @param mols.size Numeric value for the point size of molecules to plot (default = 0.1).
#' @param mols.cols A vector of colours for the molecules to plot (default = NULL).
#' @param nmols Integer of the max number of each molecule specified in 'molecules' to be plot (default = 1000).
#' @param alpha Numeric value between 0 and 1 defining the spot alpha (default = 1).
#' @param border.color Character string specifying the colour of each spot/cell border (default = white).
#' @param border.size Numeric value for the thickness of the cell segmentation border (default = NULL).
#' @param dark.background Boolean value indicating if the plot background is coloured black (default = FALSE).
#' @param blend Boolean value indicating whether to scale and blend expression values to visualize coexpression of two features (default = FALSE).
#' @param blend.threshold Numeric value defining the color cutoff from weak signal to strong signal; ranges from 0 to 1 (default = 0.5).
#' @param crop Boolean value of whether to crop the plots to area with cells only (default = FALSE).
#' @param cells Vector of character strings defining a group of cells to plot (default = NUll; plots all cells).
#' @param scale Set color scaling across multiple plots; c("features", "all", "none").
#' @param overlap Overlay boundaries from a single image to create a single plot (default = FALSE).
#' @param axes Boolean defining if to keep axes and panel background (default = FALSE).
#' @param combine Boolean value stating if to combine plots into a single patchworked ggplot object (default = TRUE).
#' @param coord.fixed Boolean value of it to plot cartesian coordinates with fixed aspect ratio (default = TURE).
#' @param assay Character string defining the Seurat assay that contains the annotated metadata corresponding to the m/z values.
#' @param column.name Character string defining the column name where the annotations are stored in the slot meta.data (default = "all_IsomerNames").
#' @param plot.exact Boolean value describing if to only plot exact matches to the metabolite search terms, else will plot all metabolites which contain serach word in name (default = TRUE).
#'
#' @returns A ggplot showing the Spatial representation of expression data of specified metabolites.
#' @export
#'
#' @examples
#' # ImageMZPlot(SeuratObj, mzs = c("Glucose", "Glutamine"))
#' # ImageMZPlot(SeuratObj, mzs = c("Glucose", "Glutamine"), plusminus = 0.05)
ImageMZAnnotationPlot <- function(object,
                                  metabolites,
                                  plusminus = NULL,
                                  fov = NULL,
                                  boundaries = NULL,
                                  cols = if (isTRUE(x = blend)) {
                                    c("lightgrey", "#ff0000", "#00ff00")
                                  } else {
                                    c("lightgrey", "firebrick1")
                                  },
                                  size = 0.5,
                                  min.cutoff = NA,
                                  max.cutoff = NA,
                                  split.by = NULL,
                                  molecules = NULL,
                                  mols.size = 0.1,
                                  mols.cols = NULL,
                                  nmols = 1000,
                                  alpha = 1,
                                  border.color = "white",
                                  border.size = NULL,
                                  dark.background = TRUE,
                                  blend = FALSE,
                                  blend.threshold = 0.5,
                                  crop = FALSE,
                                  cells = NULL,
                                  scale = c("feature", "all", "none"),
                                  overlap = FALSE,
                                  axes = FALSE,
                                  combine = TRUE,
                                  coord.fixed = TRUE,
                                  assay = "Spatial",
                                  column.name = "all_IsomerNames",
                                  plot.exact = TRUE

){
  mzs <- c()
  for (metabolite in metabolites){
    met.row <- SearchAnnotations(object,metabolite, assay = assay, column.name = column.name,search.exact = plot.exact)

    if (dim(met.row)[1] != 1){
      warning(paste("There are either none or multiple entries for the metabolite: ", metabolite,
                    "\n please check using SearchAnnotations() and FindDuplicateAnnotations"))
      stop("n entries != 1")
    }

    all_annots <- unlist(strsplit(met.row[[column.name]], "; "))

    if (length(all_annots) != 1){
      warning(paste("There are another ", (length(all_annots)-1),
                    "annotations assocated with this metabolite ( ",metabolite," )... These being: ",
                    "\n", paste(list(all_annots)), "\n Use SearchAnnotations() to see ..."))
    }

    mz <- met.row$raw_mz
    mzs <- c(mzs,mz)
  }

  plot <- ImageMZPlot(object,
                      mzs,
                      plusminus = plusminus,
                      fov = fov,
                      boundaries = boundaries,
                      cols = cols,
                      size = size,
                      min.cutoff = min.cutoff,
                      max.cutoff = max.cutoff,
                      split.by = split.by,
                      molecules = molecules,
                      mols.size = mols.size,
                      mols.cols = mols.cols,
                      nmols = nmols,
                      alpha = alpha,
                      border.color = border.color,
                      border.size = border.size,
                      dark.background = dark.background,
                      blend = blend,
                      blend.threshold = blend.threshold,
                      crop = crop,
                      cells = cells,
                      scale = scale,
                      overlap = overlap,
                      axes = axes,
                      combine = combine,
                      coord.fixed = coord.fixed
  )

  for (i in 1:length(plot)){

    plusmin_str <- ""

    if (!(is.null(plusminus))){
      plusmin_str <- paste0(" \\u00b1 ", plusminus)
    }

    plot[[i]] <- plot[[i]] +
      ggplot2::ggtitle(paste0(metabolites[i], plusmin_str))+
      ggplot2::labs(fill = paste0(metabolites[i], plusmin_str))

  }

  return(plot)

}



#' Visualise expression of m/z values in a spatial context for Spatial Seurat Objects with H&E images.
#'      - This is for plotting Seurat Spatial Metabolomic data without an image(i.e. H&E image)
#'      - This function inherits off Seurat::SpatialFeaturePlot(). Look here for more detailed documentation about inputs.
#'
#' @param object Seurat Spatial Metabolomic Object to Visualise.
#' @param mzs Vector of numeric m/z values to plot (e.g. c(400.1578, 300.1)). The function find_nearest() is used to automatically find the nearest m/z value to the ones given.
#' @param plusminus Numeric value defining the range/threshold either side of the target peak/peaks to be binned together for plotting (default = NULL).
#' @param images Character string of the name of the image to plot (default = NULL).
#' @param crop Boolean value indicating if to crop the plot to focus on only points being plotted (default = TRUE).
#' @param slot Character string indicating the assay slot to use to pull expression values form (default = "counts").
#' @param keep.scale Character string describing how to handle the color scale across multiple plots. Check Seurat::SpatialFeaturePlot() for all options (default = "feature").
#' @param min.cutoff Vector of numeric value describing the minimum cutoff values for each m/z feature (default = NA).
#' @param max.cutoff Vector of numeric value describing the maximum cutoff values for each m/z feature (default = NA).
#' @param ncol Integer defining the number of columns if plotting multiple plots (default = NULL).
#' @param combine Boolean value stating if to combine plots into a single patchworked ggplot object (default = TRUE).
#' @param pt.size.factor Numeric value defining the point size for spots when plotting (default = 1.6).
#' @param alpha Numeric value between 0 and 1 defining the spot alpha (default = 1).
#' @param image.alpha Numeric value between 0 and 1 defining the image alpha (default = 1).
#' @param stroke Numeric value describing the width of the border around the spot (default. =0.25).
#' @param interactive Boolean value of if to launch an interactive SpatialDimPlot or SpatialFeaturePlot session, see Seurat::ISpatialDimPlot() or Seurat::ISpatialFeaturePlot() for more details (default = FALSE).
#' @param information An optional dataframe or matrix of extra infomation to be displayed on hover (default = NULL).
#'
#' @returns A ggplot showing the Spatial representation of expression data of specified m/z values for Spatial data with H&E Image
#' @export
#'
#' @examples
#' # SpatialMZPlot(SeuratObj, mzs = c(400.678, 300.1))
#' # SpatialMZPlot(SeuratObj, mzs = c(400.678, 300.1), plusminus = 0.05)
SpatialMZPlot <- function(object,
                        mzs,
                        plusminus = NULL,
                        images = NULL,
                        crop = TRUE,
                        slot = "counts",
                        keep.scale = "feature",
                        min.cutoff = NA,
                        max.cutoff = NA,
                        ncol = NULL,
                        combine = TRUE,
                        pt.size.factor = 1.6,
                        alpha = c(1, 1),
                        image.alpha = 1,
                        stroke = 0.25,
                        interactive = FALSE,
                        information = NULL
                      ){

  if (is.null(mzs)){
    stop("No mz values have been supplied")
  } else{

    mz_list <- c()
    for (target_mz in mzs){
      mz_string <- find_nearest(object, target_mz)
      mz_list <- c(mz_list, mz_string)
    }
  }

  if (!(is.null(plusminus))){

    data_copy <- object
    col_names_to_plot <- c()
    plot_titles <- c()

    for (target_mz in mz_list){
      mz_integer <- as.numeric(strsplit(target_mz, "-")[[1]][2])

      meta_col <- bin.mz(data_copy, plusminus(data_copy, mz_integer, plusminus))

      col_name <- paste0(target_mz,"_plusminus_", plusminus)
      plot_name <- paste0("mz: ", round(mz_integer, 3)," \\u00b1 ", plusminus)

      col_names_to_plot <- c(col_names_to_plot,col_name)
      plot_titles <- c(plot_titles,plot_name)
      data_copy[[col_name]] = meta_col
    }

    plot <- Seurat::SpatialFeaturePlot(object = data_copy,
                             features = col_names_to_plot,
                             images = images,
                             crop = crop,
                             slot = slot,
                             keep.scale = keep.scale,
                            min.cutoff = min.cutoff,
                            max.cutoff = max.cutoff,
                            ncol = ncol,
                            combine = combine,
                            pt.size.factor = pt.size.factor,
                            alpha = alpha,
                            image.alpha = image.alpha,
                            stroke = stroke,
                            interactive = interactive,
                            information = information
    )

    for (plot_idx in seq(1, length(plot_titles))){
      plot[[plot_idx]] <- plot[[plot_idx]] +
        ggplot2::labs(fill = plot_titles[[plot_idx]])
    }

  } else {

    plot <- Seurat::SpatialFeaturePlot(object = object,
                             features = mz_list,
                             images = images,
                             crop = crop,
                             slot = slot,
                             keep.scale = keep.scale,
                            min.cutoff = min.cutoff,
                            max.cutoff = max.cutoff,
                            ncol = ncol,
                            combine = combine,
                            pt.size.factor = pt.size.factor,
                            alpha = alpha,
                            image.alpha = image.alpha,
                            stroke = stroke,
                            interactive = interactive,
                            information = information
    )

    for (plot_idx in seq(1, length(mz_list))){
      mz_integer <- as.numeric(strsplit(mz_list[plot_idx], "-")[[1]][2])
      plot[[plot_idx]] <- plot[[plot_idx]] +
        ggplot2::labs(fill = paste0("mz: ",round(mz_integer,3)))
    }
  }


  return(plot)

}



#' Visualise expression of metabolites in a spatial context for Spatial Seurat Objects with H&E images.
#'      - This is for plotting Seurat Spatial Metabolomic data without an image(i.e. H&E image)
#'      - This function inherits off Seurat::ImageFeaturePlot(). Look here for more detailed documentation about inputs.
#'
#' @param object Seurat Spatial Metabolomic Object to Visualise.
#' @param metabolites Vector of metabolite names to plot (e.g. c("Glucose", "Glutamine")). The Seurat Object provided must contain annotations in the respective assay metadata.
#' @param plusminus Numeric value defining the range/threshold either side of the target peak/peaks to be binned together for plotting (default = NULL).
#' @param images Character string of the name of the image to plot (default = NULL).
#' @param assay Character string indicating which Seurat object assay to pull data form (default = "Spatial").
#' @param crop Boolean value indicating if to crop the plot to focus on only points being plotted (default = TRUE).
#' @param slot Character string indicating the assay slot to use to pull expression values form (default = "counts").
#' @param keep.scale Character string describing how to handle the color scale across multiple plots. Check Seurat::SpatialFeaturePlot() for all options (default = "feature").
#' @param min.cutoff Vector of numeric value describing the minimum cutoff values for each m/z feature (default = NA).
#' @param max.cutoff Vector of numeric value describing the maximum cutoff values for each m/z feature (default = NA).
#' @param ncol Integer defining the number of columns if plotting multiple plots (default = NULL).
#' @param combine Boolean value stating if to combine plots into a single patchworked ggplot object (default = TRUE).
#' @param pt.size.factor Numeric value defining the point size for spots when plotting (default = 1.6).
#' @param alpha Numeric value between 0 and 1 defining the spot alpha (default = 1).
#' @param image.alpha Numeric value between 0 and 1 defining the image alpha (default = 1).
#' @param stroke Numeric value describing the width of the border around the spot (default. =0.25).
#' @param interactive Boolean value of if to launch an interactive SpatialDimPlot or SpatialFeaturePlot session, see Seurat::ISpatialDimPlot() or Seurat::ISpatialFeaturePlot() for more details (default = FALSE).
#' @param information An optional dataframe or matrix of extra infomation to be displayed on hover (default = NULL).
#' @param column.name Character string defining the column name where the annotations are stored in the slot meta.data (default = "all_IsomerNames").
#' @param plot.exact Boolean value describing if to only plot exact matches to the metabolite search terms, else will plot all metabolites which contain serach word in name (default = TRUE).
#'
#' @returns A ggplot showing the Spatial representation of expression data of specified m/z values for Spatial data with H&E Image
#' @export
#'
#' @examples
#' # SpatialMZAnnotationPlot(SeuratObj, mzs = c("Glucose", "Glutamine"))
#' # SpatialMZAnnotationPlot(SeuratObj, mzs = c("Glucose", "Glutamine"), plusminus = 0.05)
SpatialMZAnnotationPlot <- function(object,
                                    metabolites,
                                    plusminus = NULL,
                                    images = NULL,
                                    crop = TRUE,
                                    assay = "Spatial",
                                    slot = "counts",
                                    keep.scale = "feature",
                                    min.cutoff = NA,
                                    max.cutoff = NA,
                                    ncol = NULL,
                                    combine = TRUE,
                                    pt.size.factor = 1.6,
                                    alpha = c(1, 1),
                                    image.alpha = 1,
                                    stroke = 0.25,
                                    interactive = FALSE,
                                    information = NULL,
                                    column.name = "all_IsomerNames",
                                    plot.exact = TRUE

){
  mzs <- c()
  for (metabolite in metabolites){
    met.row <- SearchAnnotations(object,metabolite, assay = assay, column.name = column.name,search.exact = plot.exact)

    if (dim(met.row)[1] != 1){
      warning(paste("There are either none or multiple entries for the metabolite: ", metabolite,
                    "\n please check using SearchAnnotations() and FindDuplicateAnnotations"))
      stop("n entries != 1")
    }

    all_annots <- unlist(strsplit(met.row[[column.name]], "; "))

    if (length(all_annots) != 1){
      warning(paste("There are another ", (length(all_annots)-1),
                    "annotations assocated with this metabolite ( ",metabolite," )... These being: ",
                    "\n", paste(list(all_annots)), "\n Use SearchAnnotations() to see ..."))
    }

    mz <- met.row$raw_mz
    mzs <- c(mzs,mz)
  }

  plot <- SpatialMZPlot(object,
                        mzs,
                        plusminus = plusminus,
                        images = images,
                        crop = crop,
                        slot = slot,
                        keep.scale = keep.scale,
                        min.cutoff = min.cutoff,
                        max.cutoff = max.cutoff,
                        ncol = ncol,
                        combine = combine,
                        pt.size.factor = pt.size.factor,
                        alpha = alpha,
                        image.alpha = image.alpha,
                        stroke = stroke,
                        interactive = interactive,
                        information = information
  )

  for (i in 1:length(plot)){

    plusmin_str <- ""

    if (!(is.null(plusminus))){
      plusmin_str <- paste0(" \\u00b1 ", plusminus)
    }

    plot[[i]] <- plot[[i]]+
      ggplot2::labs(fill = paste0(metabolites[i], plusmin_str))

  }

  return(plot)

}







########################################################################################################################################################################################################################

