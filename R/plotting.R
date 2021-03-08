
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


#' Plot no. of regions detected as a function of counts
#'
#' @param matrix A region * cell matrix of counts
#'
#' @return ggplot object
#' @export
#'
#' @examples
#'
#'
plotRegionDetection <- function(mat) {
  tc <- colSums(mat)
  tcy <- colSums(mat > 0)
  plot(tc, tcy, log = "xy", xlab = "log(total Counts)", ylab = "log(No. of regions with non-zero counts)")
}


#' Plot variance explained for a PC matrix
#'
#' @param matrix N*PC matrix
#'
#' @return ggplot object
#' @export
#'
#' @examples
#'
#'
plotVarExplained <- function(mat) {
  v <- data.frame(apply(mat, 2, var)/sum(apply(mat, 2, var)))
  colnames(v) <- "VarExplained"
  if(is.null(colnames(mat))) {
    names <- paste0("PC_", 1:nrow(v))
  } else {
    names <- rownames(v)
  }
  v$PC <- factor(names, levels = names)

  ggplot(v, aes(PC, VarExplained, group=1)) + geom_point(size=3) + geom_line()

}


#' Plot numeric information on a 384-well plate layout
#'
#' @param barcodeFile File with barcodes arranged in order: 1-384
#' @param bcInfo Data frame with row.names=barcodes (matching barcodeFile) and columns with
#'               values to plot (eg. total counts)
#' @param varToPlot Which column of bcInfo to plot
#' @param groupInfo (optional) list which highlights certain wells. This can be used to provide
#'                  well grouping information. For example, if last 4 columns of a plate contain
#'                  control samples, or if last 8 corner wells are left empty etc. The list should
#'                  contain 3 values row (top to bottom, numeric), col (left to right, numeric), and
#'                  group (character/factor).
#'
#' @return Plate plot (ggplot object)
#' @export
#'
#' @examples
#'
#'

plotPlateLayout <- function(barcodeFile, bcInfo, varToPlot="count", groupInfo=NA) {

  bc <- read.delim(barcodeFile, header= F, stringsAsFactors = F)$V1
  ## prep a data frame with plate layout and matching barcodes
  g <- expand.grid(x = LETTERS[1:16], y = 1:24)
  g$n <- rep(0:15, 24)
  g$idx <- 24*g$n + g$y
  g$bc <- bc[g$idx]
  g$y <- as.factor(g$y)
  g$x <- as.factor(g$x)
  g$x <- with(g, factor(x, levels = rev(levels(x))))

  ## plot
  warning("Not yet fully implemented!")
  ggplot(g, aes(y, x, fill = idx)) +
    geom_tile() + ggtitle(label = y) +
    scale_fill_viridis_c() +
    theme_bw(base_size = 16) +
    labs(x="", y="", fill=varToPlot)

  ## map the groupInfo into this plot
#  if(!(is.na(groupInfo))) {
#    vals <- t(vapply(gr$rows, function(x) (24*(x-1))+gr$cols, numeric(3L)))
#    frame <- g[g$idx %in% vals, ]
#  }

  }


