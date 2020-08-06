
#' Plot marker genes on a UMAP
#'
#' @param umap a data frame with rownames = cell id, and 2 columns ("UMAP1", "UMAP2")
#' @param count_df a long format data frame with 3 columns ("cell", "gene" and "count")
#' @param geneID ID of gene to plot in the UMAP
#' @param title title of the plot
#'
#' @return ggplot object
#'
#' @importFrom ggplot2 geom_point theme_grey ggtitle labs scale_fill_gradient2
#' @export
#'
#' @examples
#'
plotMarker <- function(umap, count_df, geneID, title) {
  cm <- dplyr::filter(count_df, gene == geneID)
  umap$count <- cm[match(rownames(umap), cm$cell), "count"]$count

  ggplot(umap, aes(UMAP1, UMAP2, col = as.factor(louvain), fill = count)) +
    geom_point(shape = 21) + theme_grey(base_size = 16) +
    labs(col = "Louvain", fill = "Counts", title = title) +
    scale_fill_gradient2(low = "grey90", mid = "grey60", high = "black", na.value = "white")
}


#' Regions detected as a function of counts
#'
#' @param matrix A region * cell matrix of counts
#'
#' @return ggplot object
#' @export
#'
#' @examples
#'
#'
regionDetectionPlot <- function(mat) {
  tc <- colSums(mat)
  tcy <- colSums(mat > 0)
  plot(tc, tcy, log = "xy", xlab = "log(total Counts)", ylab = "log(No. of regions with non-zero counts)")
}


