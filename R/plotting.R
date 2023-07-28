
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
#' @param point_size numeric, size of the points
#' @param point_stroke numeric, stroke size of the points
#' @param line_size numeric, size of the connecting lines
#' @param line_alpha numeric, transparancy of the connecting lines
#' @param highlight_subgroup Name of the subgroup in the "color_by" column to highlight. All other points will be colored grey.
#' @param subgroup_colors vector A named vactor of colors, where names = categories present in "color_by" column, values are colors
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
                          padding = 0.4,
                          color_lines=TRUE,
                          point_size=1.5,
                          point_stroke=0.05,
                          line_size=0.5,
                          line_alpha=0.05,
                          highlight_subgroup=NA,
                          subgroup_colors=NULL,
                          match_df=NA) {

  umap_merge <- plyr::ldply(umaps, function(x) as.data.frame(x) %>% tibble::rownames_to_column("cellID"))
  # rescale individual UMAPs
  minmax <- function(x) ((x - min(x)) / (max(x) - min(x)))
  umap_merge %<>% dplyr::group_by(.id) %>% dplyr::mutate(UMAP1 = minmax(UMAP1), UMAP2 = minmax(UMAP2))

  ## update the position of other plots w.r.t center
  modes = setdiff(unique(umap_merge$.id), center_grp)
  crange <- c(0, 1)
#  range(umap_merge[umap_merge$.id == center_grp, 'UMAP1'])
  m = 0
  for(i in 1:length(modes)) {
    if(i %% 2 != 0) {
      to_add <- crange[2] + padding + m
      # 1,3,5.. :add coords to the right
      umap_merge[umap_merge$.id == modes[i], "UMAP1"] %<>% add(to_add)
      crange[2] %<>% add(to_add)

      } else {
      # 2,4,6.. : add coords to the left
      to_subtract = crange[1] - padding - m# -ve number
      umap_merge[umap_merge$.id == modes[i], "UMAP1"] %<>% add(to_subtract)
      crange[1] %<>% add(to_subtract)
    }
    m %<>% add(1)
  }

  # plot
  if(is.null(subgroup_colors)) {
    # fetch colors automatically
    subgrp_colors <- get_colors(umap_merge[[color_by]])
  } else {
    # TODO: check if enough colors are supplied
    subgrp_colors <- subgroup_colors
  }


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
     p <- p %+% geom_line(aes_string(group = 'cellID', color=color_by), size=line_size, alpha=line_alpha)
   } else {
     if(!is.na(highlight_subgroup)) {
       p <- p %+% geom_line(data = umap_merge[umap_merge$subgroup, ],
                            aes_string(group = 'cellID'),
                            size=line_size, alpha=line_alpha,
                            color='black')
     } else {
       p <- p %+% geom_line(aes_string(group = 'cellID'), size=line_size, alpha=line_alpha, color='grey')
     }
   }

  p <- p %+%  theme(legend.position = "top",
                    axis.text = element_blank(),
                    axis.ticks = element_blank(), #element_line(linetype = "solid"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank() ) +
              guides(fill = guide_legend(title.position = "top",
                                         override.aes = list(size=5)), color = "none")


  if(!is.na(highlight_subgroup)) {
    p <- p %+% geom_point(data=umap_merge_subgroup, size=point_size, pch=21, stroke=point_stroke) +
               scale_fill_manual(values = subgrp_colors) + scale_color_manual(values = subgrp_colors) +
               geom_point(data=umap_merge[!umap_merge$subgroup, ], size=point_size, pch=21,
                          alpha = 0.01, fill="grey60")
  } else {
    p <- p %+% geom_point(size=point_size, pch=21, stroke=point_stroke) +
               scale_fill_manual(values = subgrp_colors) +
               scale_color_manual(values = subgrp_colors)
  }

  p + NULL
}





#' Make a DimPlot (UMAP/PCA) from a SingleCellExperiment Object,
#' where cluster labels are numeric, and the annotation is separated
#'
#' @param object A SingleCellExperiment object (with reducedDim slot),
#'               or a data.frame with with colnames: "Dim1", "Dim2", that are to be plotted.
#'               There must be additional columns present in the input data frame or in colData(sce),
#'               that would be used for coloring/labelling the plot.
#' @param dimtoplot slot from reducedDimNames(sce) to plot
#' @param vartoplot column from colData(sce) to plot (must be categorical)
#' @param point_size numeric, size of the points
#' @param point_stroke numeric, stroke size of the points
#' @param numeric_annotation set it to TRUE to assign numeric IDs to each annotation (this avoids label overlap)
#' @param separate_legend if TRUE, legend is plotted in a new plot. By default
#'                        legend is stitched to the right of the plot.
#' @param manual_colors vector A named vactor of colors, where names = categories present in "vartoplot" column, values are colors
#' @param return_objects boolean, whether to return plot objects as list, instead of plotting
#'
#' @return ggplot object
#' @export
#'
#' @examples
#'
complexDimPlot <- function(sce,
                           dimtoplot = "X_umap",
                           vartoplot = "hour",
                           point_size = 1.5,
                           point_stroke = 0.05,
                           numeric_annotation = TRUE,
                           separate_legend = FALSE,
                           manual_colors = NULL,
                           return_objects = FALSE
                           ) {
  # Get 2D embedding
  if(class(sce) == "SingleCellExperiment") {
    umap <- reducedDim(sce, dimtoplot)
    colnames(umap) <- c("Dim1", "Dim2")
    umap %<>% cbind(colData(sce) %>% as.data.frame())
  } else if(class(sce) == "data.frame") {
    umap <- sce
  } else {
    stop("Please provide a sce object or a data.frame!")
  }
  if(vartoplot == "annotation") {
    message("Renaming vartoplot to annotation_orig")
    umap$annotation_orig <- umap[,vartoplot]
    vartoplot <- "annotation_orig"
  }
  umap[,vartoplot] %<>% factor() # only catagoricals are allowed
  # Annotation
  df <- umap %>% dplyr::group_by( !!as.name(vartoplot) ) %>%
    dplyr::summarise(Dim1 = mean(Dim1),
                     Dim2 = mean(Dim2),
                     annotation=unique( !!as.name(vartoplot) )) %>%
    dplyr::mutate(annot_num = factor(1:length(annotation)),
                  annotation = paste(annot_num, annotation, sep=": ")) %>%
    as.data.frame()
  df$annotation <- factor(df$annotation, levels=df$annotation[order(df$annot_num)], ordered=TRUE)

  # colors
  if(is.null(manual_colors)){
    manual_colors <- sctoolbox::get_colors(umap[,vartoplot])
  }
  manual_num_colors <- manual_colors
  names(manual_num_colors) <- df[match(names(manual_num_colors),
                                       df[,vartoplot]), "annotation"]

  plt <- ggplot(umap, aes_string("Dim1", "Dim2", fill=vartoplot, color=vartoplot)) +
    geom_point(size=point_size, pch=21, stroke=point_stroke) +
    scale_fill_manual(values = manual_colors) +
    scale_color_manual(values = manual_colors) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none",
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
          )


  if(isTRUE(numeric_annotation)) {
    ## plot annotation separately in another step
    # P1: with umap + numeric labels
     plt <- plt + geom_point(data=df, aes(Dim1, Dim2), size = 7, shape = 21, fill = "white") +
            geom_text(data=df, aes(Dim1, Dim2, label=annot_num), fontface="bold")

    # P2: only label + annotation
    for_num_legend <- ggplot(df, aes(Dim1, Dim2, color=annotation)) + geom_point() +
      scale_color_manual(values = manual_num_colors) +
      theme(legend.title=element_blank(),
            legend.text = element_text(face = "bold"))
    num_legend <- cowplot::get_legend(for_num_legend)
    out <- list(plot = plt,
                legend = num_legend,
                label_df = df)
    # Stitch P1+P2
    if(isTRUE(separate_legend)) {
        grid::grid.draw(plt)
        grid::grid.newpage()
        grid::grid.draw(num_legend)
        } else {
          cowplot::plot_grid(plt,
                           cowplot::plot_grid(num_legend, nrow = 1),
                           ncol = 2,
                           rel_widths = c(4,1),
                           rel_heights = c(1,1) )
        }
    } else {
      ## Plot annotation directly on UMAP if numeric_annotation=FALSE
      df$annotation <- df[,vartoplot]
      if(isTRUE(separate_legend)) message("Legends are plotted directly on the Dimplot")

      out <- plt + geom_label_repel(data=df, aes(Dim1, Dim2, label=annotation), fontface="bold") + NULL
    }

  if(isTRUE(return_objects)){
    return(out)
  }
}
