library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(SingleCellExperiment)
library(scater)
library(edgeR)
library(pheatmap)
library(limma)
library(utils)
library(stats)

#### SpaMTP Differential Peaks Analysis Functions ########################################################################################################################################################################################


#' Runs pooling of a merged Seurat Dataset to generate psudo-replicates for each sample
#'       - This function is used by run_edgeR_annotations()
#'
#' @param data.filt A Seurat Object containing count values for pooling.
#' @param idents A character string defining the idents column to pool the data against.
#' @param n An integer defining the amount of psudo-replicates to generate for each sample (default = 3).
#'
#' @returns A SinglCellExpereiment object which contains pooled (n)-psudo-replicate counts data based on the Seurat Object input
#' @export
#'
#' @examples
#' # run_pooling <- list(seuratObj, idents = "sample", n = 3)
run_pooling <- function(data.filt, idents, n) {

  cell_metadata <- data.filt@meta.data
  samples <- unique(cell_metadata[[idents]])
  message("Pooling one sample into 3 replicates...")
  nrg <- n
  for(i in c(1:length(samples))){
    set.seed(i)
    wo<-which(cell_metadata[[idents]]== samples[i])
    cell_metadata[wo,'orig.ident2']<-paste(samples[i],sample(c(1:n),length(wo)
                                                             ,replace=T,prob=rep(1/nrg,nrg)),sep='_')
  }
  gene_data <- row.names(data.filt)
  filtered.sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = data.filt@assays$Spatial$counts),
                                       colData = cell_metadata)
  #dim(filtered.sce)

  tempf=strsplit(filtered.sce@colData[["orig.ident2"]],'_')
  pid=NULL
  for(i in 1:length(tempf)){
    pidone=tempf[[i]]
    if(length(pidone)!=3){
      pidone=c(pidone[1],'yes',pidone[2])
    }
    pid=rbind(pid,pidone)
  }

  filtered.sce@colData$type=pid[,2]

  summed <- scater::aggregateAcrossCells(filtered.sce,
                                 id=SingleCellExperiment::colData(filtered.sce)[,'orig.ident2'])

  return(summed)
}




#' Runs EdgeR analysis for pooled data
#'       - This function is used by run_edgeR_annotations()
#'
#' @param pooled_data A SingleCellExperiment object which contains the pooled pseudo-replicate data.
#' @param seurat_data A Seurat object containing the merged Xenium data being analysed (this is subset).
#' @param ident A character string defining the ident column to perform differential expression analysis against.
#' @param output_dir A character string defining the ident column to perform differential expression analysis against.
#' @param run_name A character string defining the title of this DE analysis (will be used when saving DEPs to .csv file).
#' @param n An integer that defines the number of pseudo-replicates per sample (default = 3).
#' @param logFC_threshold A numeric value indicating the logFC threshold to use for defining significant genes (default = 1.2).
#' @param annotations A Boolean value indicating if the Seurat Object contains annotations. This requires annotate.SeuratMALDI() to be run (default = FALSE).
#' @param assay A character string defining the assay where the mz count data and annotations are stored (default = "Spatial").
#' @param intercept A Boolean value indicating whether to control for difference between samples in the alternative group (default = FALSE).
#'
#' @returns A list which contains the relative output requested by the "to_return" variable
#' @export
#'
#' @examples
#' # pooled_obj <- run_pooling(SeuratObj, "sample", n = 3)
#' # run_DE(pooled_obj, SeuratObj, "sample", "~/Documents/DE_output/", "run_1", n = 3, logFC_threshold = 1.2, annotations = TRUE, assay = "Spatial")
run_DE <- function(pooled_data, seurat_data, ident, output_dir, run_name, n, logFC_threshold, annotations, assay, intercept = FALSE){

  message(paste("Running edgeR DE Analysis for ", run_name, " -> with samples [", paste(unique(unlist(seurat_data@meta.data$sample)), collapse = ", "), "]"))

  annotation_result <- list()
  for (condition in unique(seurat_data@meta.data[[ident]])){
    message(paste0("Starting condition: ",condition))
    groups <- SingleCellExperiment::colData(pooled_data)[[ident]]
    groups <- gsub(condition, "Comp_A", groups)
    groups <- ifelse(groups != "Comp_A", "Comp_B", groups)



    y <- edgeR::DGEList(SingleCellExperiment::counts(pooled_data), samples=SingleCellExperiment::colData(pooled_data)$orig.ident2, group = groups)


    y$samples$condition <- groups
    y$samples$ident <- sub("_(.*)", "", y$samples$samples)

    keep <- edgeR::filterByExpr(y, group = groups, min.count = 2, min.total.count = 10)
    y <- y[keep,]
    y <- edgeR::calcNormFactors(y)

    intercept_category <- as.character(SingleCellExperiment::colData(pooled_data)[[ident]])


    if (intercept==TRUE){
      message("This code has not been fixed .... DO NOT RUN INTERCEPT = TRUE")
      sample_idx <- which(as.character(unique(SingleCellExperiment::colData(pooled_data)[[ident]])) == condition)
      #print(sample_idx)
      design <- stats::model.matrix(~groups+intercept_category)

      #print(design)
      to_remove <- -(1 + as.numeric(sample_idx))


      if (sample_idx == 1){
        design <- design
      } else{
        design <- design[,to_remove]
      }

    }

    if(intercept==FALSE) {
      design <- stats::model.matrix(~groups)
    }


    y <- edgeR::estimateDisp(y, design, robust=TRUE)
    fit <- edgeR::glmQLFit(y, design, robust=TRUE)
    res <- edgeR::glmTreat(fit, coef=ncol(fit$design), lfc=log2(logFC_threshold))
    summary(limma::decideTests(res))
    res$table$regulate <- dplyr::recode(as.character(limma::decideTests(res)),"0"="Normal","1"="Up","-1"="Down")
    de_group_edgeR <- res$table[order(res$table$PValue),]
    table(limma::decideTests(res))
    res <- edgeR::topTags(res,n = nrow(y))
    res$table$regulate <- "Normal"
    res$table$regulate[res$table$logFC>0 & res$table$FDR<0.05] <- "Up"
    res$table$regulate[res$table$logFC<0 & res$table$FDR<0.05] <- "Down"
    de_group_edgeR <- res$table[order(res$table$FDR),]
    table(de_group_edgeR$regulate)

    if (annotations){

      annotation.data <- seurat_data[[assay]]@meta.data
      if (!("mz_names" %in% colnames(annotation.data))){
        stop("Warning: There is no annotations in seurat_data[[assay]]@meta.data")
      } else {
        rownames(annotation.data) <- annotation.data$mz_names
        annotation.data_subset <- annotation.data[rownames(de_group_edgeR),]
        de_group_edgeR$all_IsomerNames <- annotation.data_subset$all_IsomerNames
      }
    }

    if (!(is.null(output_dir))){
      utils::write.csv(de_group_edgeR, paste0(output_dir,condition,"_",run_name, ".csv"))
    }

    y$DEPs <- de_group_edgeR
    annotation_result[[condition]] <- y

  }
  return(annotation_result)
}


#' Finds all differentially expressed m/z values between comparison groups specified
#'       - This function uses run_pooling() and run_DE() to pool and run EdgeR analysis
#'
#' @param data A Seurat object containing mz values for differential expression analysis.
#' @param ident A character string defining the metadata column or groups to compare mz values between.
#' @param n An integer that defines the number of psudo-replicates (pools) per sample (default = 3).
#' @param logFC_threshold A numeric value indicating the logFC threshold to use for defining significant genes (default = 1.2).
#' @param DE_output_dir A character string defining the directory path for all output files to be stored. This path must a new directory. Else, set to NULL as default.
#' @param run_name A character string defining the title of this DE analysis that will be used when saving DEPs to .csv file (default = 'FindAllDEPs').
#' @param annotations A Boolean value indicating if the Seurat Object contains annotations. This requires annotate.SeuratMALDI() to be run (default = FALSE).
#' @param assay A character string defining the assay where the mz count data and annotations are stored (default = "Spatial").
#'
#' @returns Returns an list() contains the EdgeR DE results. Pseudo-bulk counts are stored in $counts and DEPs are in $DEPs.
#' @export
#'
#' @examples
#' # FindAllDEPs(SeuratObj, "sample",DE_output_dir = "~/Documents/DE_output/", annotations = TRUE)
FindAllDEPs <- function(data, ident, n = 3, logFC_threshold = 1.2, DE_output_dir = NULL, run_name = "FindAllDEPs", annotations = FALSE, assay = "Spatial"){

  if (!(is.null(DE_output_dir))){
    if (dir.exists(DE_output_dir)){
      warning("Please supply a directory path that doesn't already exist")
      stop("dir.exists(DE_output_dir) = TRUE")
    } else{
      dir.create(DE_output_dir)
    }
  }

  #Step 1: Run Pooling to split each unique ident into 'n' number of pseudo-replicate pools
  pooled_data <- run_pooling(data,ident, n = n)

  #Step 2: Run EdgeR to calculate differentially expressed m/z peaks
  DEP_results <- run_DE(pooled_data, data, ident = ident, output_dir = DE_output_dir, run_name = run_name, n=n, logFC_threshold=logFC_threshold, annotations = annotations, assay = assay)

  # Returns an EDGEr object which contains the pseudo-bulk counts in $counts and DEPs in $DEPs
  return(DEP_results)

}




#' Generates a Heatmap of DEPs generated from edgeR analysis run using FindAllDEPs().
#'       - this function uses pheatmap() to plot data
#'
#' @param edgeR_output A list containing outputs from edgeR analysis (from FindAllDEPs()). This includes pseudo-bulked counts and DEPs.
#' @param n A numeric integer that defines the number of UP and DOWN regulated peaks to plot (default = 25).
#' @param FDR_threshold An integer that defines the FDR_threshold to use for defining most significant results (default = 0.05).
#' @param scale A character string indicating if the values should be centered and scaled in either the row direction or the column direction, or none. Corresponding values are "row", "column" and "none"
#' @param color A vector of colors used in heatmap (default = ggplot2::colorRampPalette(c("navy", "white", "red"))(50)).
#' @param cluster_cols A boolean values determining if columns should be clustered or hclust object (default = F).
#' @param cluster_rows A boolean values determining if rows should be clustered or hclust object (default = T).
#' @param fontsize_row A numeric value defining the fontsize of rownames (default = 15).
#' @param fontsize_col A numeric value defining the fontsize of colnames (default = 15).
#' @param cutree_cols A numeric value defining the number of clusters the columns are divided into, based on the hierarchical clustering(using cutree), if cols are not clustered, the argument is ignored (default = 9).
#' @param silent A boolean value indicating if the plot should not be draw (default = TRUE).
#' @param plot_annotations A boolean value indicating if metabolite annotation names should be plot rather then m/z values. Annotations = TRUE must be used in FindAllDEPs() for edgeR output to include annotations (default = FALSE).
#'
#' @returns A heatmap plot of significantly differentially expressed peaks defined in the edgeR ouput object.
#' @export
#'
#' @examples
#' # DEPs <- FindAllDEPs(SeuratObj, "sample")
#'
#' # DEPsHeatmap(DEPs)
DEPsHeatmap <- function(edgeR_output, n = 25, FDR_threshold = 0.05, scale="row", color= ggplot2::colorRampPalette(c("navy", "white", "red"))(50),cluster_cols = F,cluster_rows = T, fontsize_row = 15, fontsize_col = 15, cutree_cols = 9, silent = TRUE, plot_annotations = FALSE){

  message("Generating Heatmap .......")

  degs <- edgeR_output$DEPs
  degs$gene <- rownames(degs)
  degs <- subset(degs, FDR < FDR_threshold)
  degs <- degs[order(degs$logFC),]
  top <- utils::head(degs, n=n)
  bot <- utils::tail(degs, n=n)
  degs <- merge(x = top, y = bot, all = TRUE)
  degs <- subset(degs, FDR < FDR_threshold)
  rownames(degs) <- degs$gene

  col_annot <- data.frame(sample = edgeR_output$samples$ident)
  row.names(col_annot) <- colnames(as.data.frame(edgeR::cpm(edgeR_output,log=TRUE)))

  mtx <- as.matrix(as.data.frame(edgeR::cpm(edgeR_output,log=TRUE))[rownames(degs),])
  if (plot_annotations){
    if (is.null(edgeR_output$DEPs)){
      warning("There are no annotations present in the edgeR_output object. Run 'annotate.SeuratMALDI()' prior to 'FindAllDEPs' and set annotations = TRUE .....\n Heatmap will plot default m/z values ... ")
    } else{
      DEPs.subset <- edgeR_output$DEPs[rownames(edgeR_output$DEPs) %in% rownames(mtx),]
      rownames(mtx) <- DEPs.subset$all_IsomerNames
    }
  }

  p <- pheatmap::pheatmap(mtx,scale=scale,color=color,cluster_cols = cluster_cols, annotation_col=col_annot, cluster_rows = cluster_rows,
                fontsize_row = fontsize_row, fontsize_col = fontsize_col, cutree_cols = cutree_cols, silent = silent)
  return(p)
}


########################################################################################################################################################################################################################
