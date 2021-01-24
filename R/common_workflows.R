

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

monocleDNAworkflow <- function(counts, nDim=20) {
  ## make CDS
  cds <- cicero::make_atac_cds(counts, binarize = TRUE)
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
#'
#' @return Annotated LSA object (as list)
#' @export
#'
#' @examples
#'

lsaWorkflowVB <- function(binCounts_filt, genes_gr, qc_df=NA) {

  ## Run LSA
  filt_lsa <- lsa_wrapper(binCounts_filt, nPC = 30, nK = 20, outPdf = NA)
  filt_lsa$umap$totalCount <- colSums(binCounts_filt)
  if(!(is.na(qc_df))) {
    filt_lsa$umap <- qc_df %>% mutate(id = paste(sample, bc, sep ="_")) %>% merge(filt_lsa$umap, ., by.x=0, by.y="id")
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
    lapply(function(x) gene_to_symbol[match(x, gene_to_symbol$ensembl_gene_id), "external_gene_name"])

  ## return the list object
  filt_lsa$topGenes <- topGenes
  filt_lsa$topGenes_names <- topGenes_names
  return(filt_lsa)
}

