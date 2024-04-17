#' Helper function for suppressing function progress messages
#'
#' @param message_text Character string containing the message being shown
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the messsage will be suppressed (default = TRUE).
#'
#'
#' @examples
#' verbose_message("Finished!", TRUE)
verbose_message <- function(message_text, verbose) {
  if (verbose) {
    message(message_text)
  }
}



#cite: https://github.com/alikhuseynov/add-on_R/blob/develop/R/subset_obj_seurat_v2.R

#'@importFrom magrittr %>% %<>%
NULL

#' Intermediate solution to \code{subset()}:
#' subset FOVs/centroids if selected cells are NOT found in each FOV
#' NOTE: some code parts and args are taken from SeuratObject

#' Function params/args:
#' @param object An S4 object or A \code{FOV} object
#' @param subset Logical expression indicating features/variables to keep
#' @param cells A vector of cells to keep; if \code{NULL}, defaults to all cells
#' @param idents A vector of identity classes to keep
#' @param features A vector of feature names or indices to keep
#' @param Update.slots If to update slots of an object
#' @param Update.object If to update final object, default to TRUE.
#' @param ... Arguments passed to \code{subset()} and other methods
#'
#' @return A subset Seurat object
#' @export
#'
#' @examples
#' # sub <- subset_obt(seurat.obj, idents = "Sample1")
subset_SPM <- function(
    object = NULL,
    subset = NULL,
    cells = NULL,
    idents = NULL,
    features = NULL,
    Update.slots = TRUE,
    Update.object = TRUE,
    ...)
{

  if (Update.slots) {
    message("Updating object slots..")
    object %<>% UpdateSlots()
  }

  message("Cloing object..")
  obj_subset <- object

  # sanity check - use only cell ids (no indices)
  if (all(is.integer(cells))) {
    cells <- Cells(obj_subset)[cells]
  }

  if (!missing(subset) || !is.null(idents)) {
    message("Extracting cells matched to `subset` and/or `idents`")
  }

  if (class(obj_subset) == "FOV") {
    message("object class is `FOV` ")
    cells <- Cells(obj_subset)
  } else if (!class(obj_subset) == "FOV" && !missing(subset)) {
    subset <- enquo(arg = subset)
    # cells to keep in the object
    cells <-
      WhichCells(object = obj_subset,
                 cells = cells,
                 idents = idents,
                 expression = subset,
                 return.null = TRUE, ...)
  } else if (!class(obj_subset) == "FOV" && !is.null(idents)) {
    cells <-
      WhichCells(object = obj_subset,
                 cells = cells,
                 idents = idents,
                 return.null = TRUE, ...)
  } else if (is.null(cells)) {
    cells <- Cells(obj_subset)
  }

  # added support for object class `FOV`
  if (class(obj_subset) == "FOV") {
    message("Matching cells for object class `FOV`..")
    cells_check <- any(obj_subset %>% Cells %in% cells)
  } else {
    # check if cells are present in all FOV
    message("Matching cells in FOVs..")
    cells_check <-
      lapply(Images(obj_subset) %>% seq,
             function(i) {
               any(obj_subset[[Images(obj_subset)[i]]][["centroids"]] %>% Cells %in% cells)
             }) %>% unlist
  }

  if (all(cells_check)) {
    message("Cell subsets are found in all FOVs!", "\n",
            "Subsetting object..")
    obj_subset %<>% base::subset(cells = cells,
                                 idents = idents,
                                 features = features,
                                 ...)
    # subset FOVs
    message("Subsetting FOVs..")
    fovs <-
      lapply(Images(obj_subset) %>% seq, function(i) {
        base::subset(x = obj_subset[[Images(obj_subset)[i]]],
                     cells = cells,
                     idents = idents,
                     features = features,
                     ...)
      })
    # replace subsetted FOVs
    for (i in fovs %>% seq) { obj_subset[[Images(object)[i]]] <- fovs[[i]] }

  } else {
    # if cells are present only in one or several FOVs:
    # subset FOVs
    fovs <-
      lapply(Images(obj_subset) %>% seq, function(i) {
        if (any(obj_subset[[Images(obj_subset)[i]]][["centroids"]] %>% Cells %in% cells)) {
          message("Cell subsets are found only in FOV: ", "\n", Images(obj_subset)[i])
          message("Subsetting Centroids..")
          base::subset(x = obj_subset[[Images(obj_subset)[i]]],
                       cells = cells,
                       idents = idents,
                       features = features,
                       ...)
        }
      })
    # remove FOVs with no matching cells
    message("Removing FOVs where cells are NOT found: ", "\n",
            paste0(Images(object)[which(!cells_check == TRUE)], "\n"))
    # replace subsetted FOVs
    for (i in fovs %>% seq) { obj_subset[[Images(object)[i]]] <- fovs[[i]] }

    # subset final object
    message("..subset final object")
    obj_subset %<>%
      base::subset(cells = cells,
                   idents = idents,
                   features = features,
                   ...)
  }

  if (Update.object && !class(obj_subset) == "FOV") {
    message("Updating object..")
    obj_subset %<>% UpdateSeuratObject() }

  message("Object is ready!")
  return(obj_subset)

}



########################################################################
