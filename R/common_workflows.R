

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

monocleDNAworkflow <- function(counts) {
  ## make CDS
  cds <- cicero::make_atac_cds(counts, binarize = TRUE)
  ## cluster
  set.seed(2017)
  cds %<>% monocle3::detect_genes()
  cds %<>% monocle3::estimate_size_factors()
  cds %<>% monocle3::preprocess_cds(method = "LSI", num_dim = 20)
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
