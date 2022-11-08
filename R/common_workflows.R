

#' Monocle ATAC/DNA workflow
#'
#' @param counts matrix: either wide (regions*cells) or tidy (regions, cell, count) format
#' @param isTidy bool: are the counts in tidy format?
#'
#' @return cellDataSet object
#' @export
#'
#' @examples
#'
#'

monocleDNAworkflow <- function(counts, nDim=20, binary=FALSE) {
  if(class(counts) == "cell_data_set") {
    cds <- counts
  } else {
    ## make CDS
    cds <- cicero::make_atac_cds(counts, binarize = binary)
  }
  ## cluster
  set.seed(2017)
  cds %<>% monocle3::detect_genes()
  cds %<>% monocle3::estimate_size_factors()
  cds %<>% monocle3::preprocess_cds(method = "LSI", num_dim = nDim)
  cds %<>% monocle3::reduce_dimension(reduction_method = 'UMAP', preprocess_method = "LSI")

  ## PAGA and graph learning
  cds %<>% monocle3::cluster_cells(reduction_method = "UMAP")
  cds %<>% monocle3::learn_graph()

  return(cds)
}


#' Seurat workflow
#'
#' @param counts count matrix (rows = regions, columns = cells)
#'
#' @return seurat object
#' @export
#'
#' @examples
#'

seuratWorkflow <- function(counts) {

  seu <- Seurat::CreateSeuratObject(counts)
  seu %<>% Seurat::FindVariableFeatures(selection.method = "vst", nfeatures = 3000)
  seu %<>% Seurat::ScaleData(features = rownames(seu))
  seu %<>% Seurat::RunPCA(features = VariableFeatures(seu))
  seu %<>% Seurat::FindNeighbors(dims = 1:20)
  seu %<>% Seurat::FindClusters(resolution = 0.5)
  seu %<>% Seurat::RunUMAP(dims = 1:20)
  mat <- seu@meta.data
  mat %<>% cbind(seu@reductions$umap@cell.embeddings)

  ggplot2::ggplot(mat, aes(UMAP_1, UMAP_2, fill = log10(nCount_DNA))) +
    ggplot2::geom_point(pch=21) + ggplot2::scale_fill_viridis_c() +
    ggplot2::ggtitle("Total Counts") + ggplot2::theme_grey(base_size = 14)

  return(seu)
}


#' vivekbhr's workflow for LSA and annotation
#'
#' @param binCounts_filt Filtered counts in genomic bins (rows = regions, columns = cells)
#' @param genes_gr GenomicRanges with gene names and their locations.
#' @param qc_df (optional) DataFrame with per-cell quality control (or other) information
#' @param gene_symbol data frame with ensembl gene ID and corresponding symbol (from biomart)
#'
#' @return Annotated LSA object (as list)
#' @export
#'
#' @examples
#'

lsaWorkflowVB <- function(binCounts_filt, genes_gr, qc_df=NA, gene_symbols=NA) {

  ## Run LSA
  filt_lsa <- lsa_wrapper(binCounts_filt, nPC = 30, nK = 20, outPdf = NA)
  filt_lsa$umap$totalCount <- colSums(binCounts_filt)
  if(!(is.na(qc_df))) {
    filt_lsa$umap <- qc_df %>% mutate(id = paste(sample, bc, sep ="_")) %>% merge(filt_lsa$umap, ., by.x=0, by.y="id")
  } else {
    filt_lsa$umap$Row.names <- rownames(filt_lsa$umap)
  }

  # plot
  ggplot(filt_lsa$umap, aes(UMAP1, UMAP2, color = factor(louvain))) + geom_point() +
    scale_color_brewer(palette = "Paired") + labs(color = "Cluster")

  ## Get top genes per cluster
  topGenes <- topGenesByCluster(data.frame(cells = filt_lsa$umap$Row.names,
                                           cluster = filt_lsa$umap$louvain),
                                binCounts_filt, n = 50)
  topGenes %<>% lapply(charToGRanges) %>% lapply(function(x) subsetByOverlaps(genes_gr, x))
  topGenes_names <- topGenes %>% lapply(function(x) return(x$gene_id)) %>% lapply(sanitizeGeneID) %>%
    lapply(function(x)  gene_symbols[match(x,  gene_symbols$ensembl_gene_id), "external_gene_name"])

  ## return the list object
  filt_lsa$topGenes <- topGenes
  filt_lsa$topGenes_names <- topGenes_names
  return(filt_lsa)
}


# Wrapper for SnapATAC workflow up until PCA step.
# Args:
#   snap_file (string): path to snap file
#   promooter.df (dataframe): dataframe with promoter definitions as shown in SnapATAC tutorials (dataframe not file name; see examples in markdown)
#   blacklist.df (dataframe): dataframe with blacklist region definitions as shown in SnapATAC tutorials (dataframe not file name; see examples in markdown)
#   fragment_number_threshold (int): threshold for number of unique fragments per cell
#   promoter_ratio_range (c(float, float)): vector with lower and upper bound of acceptable fraction of reads in promoter regions for cell filtering as used in SnapATAC tutorial.
#   window_z_range (c(float, float)): vector with lower and upper bound of acceptable window z-scores for non-zero entries for site filtering as used in SnapATAC tutorial.
#   sample_name (string): sample_name provided to SnapATAC
#   pc.num (int): total PCs to compute
# Returns:
#   snap object: SnapATAC object
# Notes:
#   This uses single core because multithreaded implementation interferes with Knitr. In running this code, any do.par=FALSE and num.cores=1 could be changed as needed.
snapatac_workflow = function(region_cell_matrix, snap_file = NA, promoter.df=NULL, blacklist.df=NULL,
                             fragment_number_threshold=500, promoter_ratio_range=c(0.2, 0.8),
                             window_z_range=c(-1.5, 1.5), sample_name='default', pc.num=50) {
  if (!is.na(snap_file)) {
    x.sp = createSnap(
      file=snap_file,
      sample=sample_name,
      do.par=FALSE,
      num.cores=1)

    plotBarcode(
      obj=x.sp,
      pdf.file.name=NULL,
      pdf.width=7,
      pdf.height=7,
      col="grey",
      border="grey",
      breaks=50
    )

    x.sp = filterCells(
      obj=x.sp,
      subset.names=c("fragment.num", "UMI"),
      low.thresholds=c(fragment_number_threshold, fragment_number_threshold),
      high.thresholds=c(Inf, Inf)
    )
    x.sp = addBmatToSnap(x.sp, bin.size=5000, num.cores=1)

  } else {
    x.sp <- newSnap()
    x.sp@bmat <- t(region_cell_matrix)
    x.sp@barcode <- colnames(region_cell_matrix)
    x.sp@feature <- GRanges(seqnames = rownames(x.sp) %>% gsub("(chr.*):([0-9]*)-([0-9]*)", "\\1", .),
                            IRanges(start = rownames(x.sp) %>% gsub("(chr.*):([0-9]*)-([0-9]*)", "\\2", .) %>%
                                      as.numeric(),
                                    end = rownames(x.sp) %>% gsub("(chr.*):([0-9]*)-([0-9]*)", "\\2", .) %>%
                                      as.numeric()))
  }


  # Optionally filter cells based on ratio of reads in promoters
  if (!is.null(promoter.df)) {
    promoter.gr = GRanges(promoter.df[,1], IRanges(promoter.df[,2], promoter.df[,3]))
    ov = findOverlaps(x.sp@feature, promoter.gr)
    idy = queryHits(ov)
    promoter_ratio = SnapATAC::rowSums(x.sp[,idy, mat="bmat"], mat="bmat") / SnapATAC::rowSums(x.sp, mat="bmat")
    plot(log(SnapATAC::rowSums(x.sp, mat="bmat") + 1,10), promoter_ratio, cex=0.5, col="grey", xlab="log(count)", ylab="FIP Ratio", ylim=c(0,1 ))
    idx = which(promoter_ratio > promoter_ratio_range[1] & promoter_ratio < promoter_ratio_range[2])
    x.sp = x.sp[idx,]
  }

  x.sp = makeBinary(x.sp, mat="bmat");

  # Filter out non-standard contigs if present
  idy2 = grep("chrM|random", x.sp@feature)

  if (!is.null(blacklist.df)) {
    black_list.gr = GRanges(
      blacklist.df[,1],
      IRanges(blacklist.df[,2], blacklist.df[,3])
    )
    idy1 = queryHits(findOverlaps(x.sp@feature, black_list.gr))

  } else {
    # No blacklist provided, so just ignore
    idy1 = c()
  }

  idy = unique(c(idy1, idy2))
  x.sp = x.sp[,-idy, mat="bmat"]

  # Filter based on frequency
  x.sp = filterBins(
    x.sp,
    low.threshold=window_z_range[1],
    high.threshold=window_z_range[2],
    mat="bmat"
  )

  #plotBinCoverage(
  #  x.sp,
  #  pdf.file.name=NULL,
  #  col="grey",
  #  border="grey",
  #  breaks=10,
  #  xlim=c(-6,6)
  #)

  x.sp = runJaccard(
    obj = x.sp,
    tmp.folder=tempdir(),
    mat = "bmat",
    max.var=2000,
    seed.use=10
  )

  x.sp = runNormJaccard(
    obj = x.sp,
    tmp.folder=tempdir(),
    ncell.chunk=1000,
    method="zscore",
    row.center=TRUE,
    row.scale=TRUE,
    low.threshold=-5,
    high.threshold=5,
    do.par=FALSE,
    num.cores=1,
    seed.use=10
  )

  x.sp = runDimReduct(
    x.sp,
    pc.num=pc.num,
    input.mat="jmat",
    method="svd",
    center=TRUE,
    scale=FALSE,
    seed.use=10
  )

  rownames(x.sp@bmat) = x.sp@barcode
  colnames(x.sp@bmat) = as.character(1:ncol(x.sp@bmat))

  return(x.sp)
}

# Reperform Jaccard + PCA on SnapATAC object (used for redoing after modifying matrix)
# Args:
#   snao_obj (snap object): snap object
#   pc.num (int): total PCs to compute
# Returns:
#   snap object: SnapATAC object
# Notes:
#   This uses single core because multithreaded implementation interferes with Knitr. In running this code, any do.par=FALSE and num.cores=1 could be changed as needed.
snapatac_rerun_jaccard = function(snap_obj, pc.num=50) {

  snap_obj = runJaccard(
    obj = snap_obj,
    tmp.folder=tempdir(),
    mat = "bmat",
    max.var=2000,
    ncell.chunk=1000,
    do.par=FALSE,
    num.cores=1,
    seed.use=10
  )

  snap_obj = runNormJaccard(
    obj = snap_obj,
    tmp.folder=tempdir(),
    ncell.chunk=1000,
    method="normOVE",
    row.center=TRUE,
    row.scale=TRUE,
    low.threshold=-5,
    high.threshold=5,
    do.par=FALSE,
    num.cores=1,
    seed.use=10
  )

  snap_obj = runDimReduct(
    snap_obj,
    pc.num=pc.num,
    input.mat="jmat",
    method="svd",
    center=TRUE,
    scale=FALSE,
    seed.use=10
  )
}

