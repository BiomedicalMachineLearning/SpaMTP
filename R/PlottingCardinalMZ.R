library(Cardinal)
library(graphics)

#' Plots mean mass spectrometry intensity values with m/z peak annotation labels
#' @param data Cardinal object which contains annotations in the feature metadata slot
#' @param annotations Vector of character strings describing the annotations labels to show on the plot (default = NULL).
#' @param main Character string defining the main title of the plot (default = "mass intensity plot").
#' @param annotation.column Character string defining the column name of the feature metadata slot where the annotation names are stored (default = "feature_metadata.all_IsomerNames").
#'
#' @return Mass intensity plot which has the specified m/z values annotated
#' @export
#'
#' @examples
#' # AnnotatedIntensityPlot(CardinalObj, annotations = "Glutamine")
AnnotatedIntensityPlot <- function (data, annotations = NULL, main = "mass intensity plot", annotation.column = "feature_metadata.all_IsomerNames"){

  data.plot <- Cardinal::plot(data)

  if (!(is.null(annotations))){
    all_annotations <- Cardinal::featureData(data)[[annotation.column]]
    annotations_to_plot <- all_annotations
    matching <- annotations_to_plot %in% annotations
    annotations_to_plot[!matching] <- NA
    annotations <- annotations_to_plot
  }

  print(Cardinal::plot(data, main  = main))
  graphics::text(x =data.plot$facets[[1]][[1]]$x,
       y = data.plot$facets[[1]][[1]]$y,
       labels = annotations,
       cex = 0.6, pos = 4, col = "red")
  final.plot <- recordPlot()
  return(final.plot)
}
