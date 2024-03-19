library(Seurat)
library(ggridges)


#' Normalizes m/z intensity data stored in a Seurat Object
#'
#' @param data Seurat Object to be normalized.
#' @param normalisation.type Character string defining the normalization method to run. Options are either c("TIC", "LogNormalize", "RC") which represent Total Ion Current (TIC) normalization, Log Normalization or counts per million (RC), respectively (default = "TIC").
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

  if (!(normalisation.type == "TIC" | normalisation.type == "LogNormalize"| normalisation.type == "RC")) {
    stop("Error: incorrect normalisation.type is select. Please enter either 'LogNormalize', 'TIC' or 'RC'")
  }

  if (is.null(scale.factor)) {
    scale.factor <- 10000
  }


  if (normalisation.type == 'TIC') {
    scale.factor <- length(rownames(data))
    normalisation.type <- "RC"
  }

  normalised.data <- Seurat::NormalizeData(data, normalization.method = normalisation.type, scale.factor = scale.factor)

  return(normalised.data)

}

TMMNormalize <- function(combined.obj, ident, refIdent = NULL, normalisation.type = "CPM", CPM.scale.factor = 1e6, assay = "Spatial", slot = "counts") {

  data_list <- list()
  Idents(combined.obj) <- ident
  if (length(unique(Idents(combined.obj))) <= 1){
    stop("Specified ident has 0 or 1 catagory in seurat object. Length of Idents must be > 1 for TMM (between sample) normalisation factors to be calculated")

  }

  if (is.null(refIdent)){
    stop("A refIdent must be specified")
  }

  if (!(is.null(refIdent))) {
    if (!(refIdent %in% unique(Idents(combined.obj)))){
      stop("The refIdent supplied is not present in the ident column. Please specify a group that is found within the ident column specififed")
    }
  }

  if (normalisation.type == "CPM") {
    normalisation.type <- "RC"
  }


  for (name in unique(Idents(combined.obj))){
    suppressWarnings({
      suppressMessages({
        sub <- subset_opt(combined.obj, idents = name)
      })
    })

    data_list[[name]] <- sub
  }

  df <- data.frame(mz = rownames(data_list[[1]]))
  rownames(df) <- df$mz

  for (dataset in names(data_list)){
    rowsum <- rowSums(data_list[[dataset]][[assay]][slot])
    df[dataset] <- rowsum
  }
  mtx <- Matrix::as.matrix(df[names(data_list)])

  factors <- edgeR::calcNormFactors(mtx, method = "TMM", refColumn = refIdent)

  norm_data_list <- list()
  for (name in names(data_list)){
    norm.data <- NormalizeSeuratDataX(data_list[[name]], normalisation.type = normalisation.type, scale.factor = (CPM.scale.factor * factors[name]), assay = assay, slot = slot)
    norm_data_list[[name]] <- norm.data
  }

  merged.data <- JoinLayers(merge(norm_data_list[[1]], y = norm_data_list[2: length(names(norm_data_list))]), merge.data = TRUE)
return(merged.data)
}



########### Pre-processing plots #############
statPlot <- function (seurat.obj, group.by = NULL, assay = "Spatial", slot = "counts", title = "VlnPlot", x.lab = "var", y.lab = "intensity", bottom.cutoff = NULL, top.cutoff = NULL, log.data = FALSE){

  data_list <- list()
  if (!(is.null(group.by))){
    Idents(seurat.obj) <- group.by
    for (ident in unique(Idents(seurat.obj))){
      suppressWarnings({
        suppressMessages({
          sub <- subset_opt(seurat.obj, idents = ident)
        })
      })

      data_list[[ident]] <- sub
    }
  } else {
    data_list[["data"]] <- seurat.obj
  }

  df <- data.frame(mz = rownames(data_list[[1]]))
  rownames(df) <- df$mz

  for (dataset in names(data_list)){
    rowsum <- rowSums(data_list[[dataset]][[assay]][slot])
    df[dataset] <- rowsum
  }

  df2 <- tidyr::pivot_longer(df, cols =  names(data_list), names_to = "var", values_to = "x")

  df2 <- df2 %>% dplyr::arrange(x)

  if (!(is.null(bottom.cutoff))){
    df2 <- df2 %>%
      dplyr::group_by(var) %>%
      dplyr::mutate(bottom_cutoff = quantile(x, bottom.cutoff))

    df2 <- df2 %>%  dplyr::group_by(var) %>%dplyr::filter(x >= bottom_cutoff)
    message(paste0("Removing bottom ", bottom.cutoff*100, "% of datapoints"))
  }

  if (!(is.null(top.cutoff))){
    df2 <- df2 %>%
      dplyr::group_by(var) %>%
      dplyr::mutate(top_cutoff = quantile(x, 1- top.cutoff))

    df2 <- df2 %>%  dplyr::group_by(var) %>% dplyr::filter(x <= top_cutoff)
    message(paste0("Removing top ", top.cutoff*100, "% of datapoints"))

  }

  if (log.data) {
    df2$x <- log10(df2$x + 1)
  }

  return(df2)
}

MZRidgePlot <- function (seurat.obj, group.by = NULL, assay = "Spatial", slot = "counts", title = "VlnPlot", x.lab = "var", y.lab = "intensity", bottom.cutoff = NULL, top.cutoff = NULL, bins = 1000,log.data = FALSE, cols = NULL){
  data <- statPlot(seurat.obj = seurat.obj,
                   group.by = group.by,
                   assay = assay,
                   slot = slot,
                   title = title,
                   x.lab = x.lab,
                   y.lab = y.lab,
                   bottom.cutoff = bottom.cutoff,
                   top.cutoff = top.cutoff,
                   log.data = log.data)

  ridge_plot <- ggplot(data, aes(y=var, x=x,  fill=var)) +
    ggridges::geom_density_ridges(alpha=0.6, stat="binline", bins=bins) +
    ggridges::theme_ridges() +  labs(title = title, x = x.lab, y = y.lab)

  return(ridge_plot)
}

MZVlnPlot <- function (seurat.obj, group.by = NULL, assay = "Spatial", slot = "counts", title = "VlnPlot", x.lab = "var", y.lab = "intensity", show.points = TRUE, bottom.cutoff = NULL, top.cutoff = NULL,log.data = FALSE, cols = NULL){

  data <- statPlot(seurat.obj = seurat.obj,
                   group.by = group.by,
                   assay = assay,
                   slot = slot,
                   title = title,
                   x.lab = x.lab,
                   y.lab = y.lab,
                   bottom.cutoff = bottom.cutoff,
                   top.cutoff = top.cutoff,
                   log.data = log.data)

  violin_plot <- ggplot(data, aes(x = var , y = x, fill = var)) +
    geom_violin() +
    theme_classic() +
    labs(title = title, x = x.lab, y = y.lab)

  if (show.points){
    violin_plot <- violin_plot + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)
  }


  return(violin_plot)
}

MZBoxPlot <- function (seurat.obj, group.by = NULL, assay = "Spatial", slot = "counts", title = "VlnPlot", x.lab = "var", y.lab = "intensity", show.points = TRUE, bottom.cutoff = NULL, top.cutoff = NULL,log.data = FALSE, cols = NULL){
  data <- statPlot(seurat.obj = seurat.obj,
                   group.by = group.by,
                   assay = assay,
                   slot = slot,
                   title = title,
                   x.lab = x.lab,
                   y.lab = y.lab,
                   bottom.cutoff = bottom.cutoff,
                   top.cutoff = top.cutoff,
                   log.data = log.data)

  box_plot <- ggplot(data, aes(x = var , y = x, fill = var)) +
    geom_boxplot() +
    theme_classic() +
    labs(title = title, x = x.lab, y = y.lab)

  if (show.points){
    box_plot <- box_plot + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)
  }

  if (!(is.null(cols))){
    box_plot <- box_plot + scale_fill_manual(values = cols)
  }

  return(box_plot)
}













