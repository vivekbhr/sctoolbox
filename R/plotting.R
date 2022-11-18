
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


#' Plot multiple UMAPs in the same plot where the points (cells) are connectec
#'
#' For example tChIC data with ChIC + RNA UMAP
#'
#' @param umaps names list of data.frames with columns (UMAP1, UMAP2) and cell IDs in rownames.
#' @param color_by Column ID to color the UMAP by
#' @param center_grp Name of the list element to use as "center" of the joint UMAP
#' @param highlight_subgroup Name of the subgroup in the "color_by" column to highlight. All other points will be colored grey.
#' @param match_df In case cellIDs are not matched with themselves, provide a dataframe with source->target cell matchings.
#'                 NOTE: this only works for 2 modalities, columns of the data. frame must be named "source" and "target"
#'
#'
#' @return ggplot2 object
#' @export
#'
#' @examples
#'
#'
plotMultiUMAP <- function(umaps,
                          color_by='annotation',
                          center_grp= 'rna',
                          color_lines=TRUE,
                          match_df=NA,
                          highlight_subgroup=NA) {

  umap_merge <- plyr::ldply(umaps, function(x) as.data.frame(x) %>% tibble::rownames_to_column("cellID"))

  ## update the position of other plots w.r.t center
  crange <- range(umap_merge[umap_merge$.id == center_grp, 'UMAP1'])
  padding = max(abs(crange)) * (40/100) # 40% padding
  modes = setdiff(unique(umap_merge$.id), center_grp)

  for(i in 1:length(modes)) {
    m = i+1
    if(i %% 2 == 0) {
      # 2,4,6.. : add coords to the left
      to_subtract = (m*min(crange)) - padding# -ve number
      umap_merge[umap_merge$.id == modes[i], "UMAP1"] %<>% add(to_subtract)

    } else {
      # 1,3,5.. :add coords to the right
      to_add = (m*max(crange)) + padding
      umap_merge[umap_merge$.id == modes[i], "UMAP1"] %<>% add(to_add)
    }
    i = m
  }

  # plot
  subgrp_colors <- get_colors(umap_merge[[color_by]])

  if(!is.na(highlight_subgroup)) {
    umap_merge$subgroup <- umap_merge[[color_by]] == highlight_subgroup
    umap_merge[is.na(umap_merge$subgroup), "subgroup"] <- FALSE
    umap_merge_subgroup <- umap_merge[umap_merge$subgroup, ] #subset (only for subgroup)
  }

  p <- ggplot(umap_merge, aes_string('UMAP1', 'UMAP2', fill=color_by)) +
    theme_minimal(base_size = 14)

  if(!(is.na(match_df))) {
    rownames(umap_merge) <- umap_merge$cellID
    umap_merge[match_df$source, 'cellID'] <- match_df$target
  }

 if(isTRUE(color_lines) & is.na(highlight_subgroup)) {
     p <- p %+% geom_line(aes_string(group = 'cellID', color=color_by), size=0.5, alpha=0.05)
   } else {
     if(!is.na(highlight_subgroup)) {
       p <- p %+% geom_line(data = umap_merge[umap_merge$subgroup, ],
                            aes_string(group = 'cellID'),
                            size=0.5, alpha=0.05,
                            color='black')
     } else {
       p <- p %+% geom_line(aes_string(group = 'cellID'), size=0.5, alpha=0.05, color='grey')
     }
   }

  p <- p %+%  theme(legend.position = "top", axis.text.x = element_blank(),
              axis.ticks = element_line(linetype = "solid" )) +
              guides(fill = guide_legend(title.position = "top"), color = "none")


  if(!is.na(highlight_subgroup)) {
    p <- p %+% geom_point(data=umap_merge_subgroup, size=1.5, pch=21) +
               scale_fill_manual(values = subgrp_colors) +
               geom_point(data=umap_merge[!umap_merge$subgroup, ], size=1.5, pch=21,
                          alpha = 0.01, fill="grey60")
  } else {
    p <- p %+% geom_point(size=1.5, pch=21) + scale_fill_manual(values = subgrp_colors)
  }

  p + NULL
}
