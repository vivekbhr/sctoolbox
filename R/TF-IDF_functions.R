
safe_tfidf_multiply = function(tf, idf) {
  tf = Matrix::t(tf)
  tf@x <- tf@x * rep.int(idf, diff(tf@p))
  tf = Matrix::t(tf)
  return(tf)
}


#' Perform TF-IDF transformation of a matrix
#'
#' @param bmat sparse matrix of class dgcMatrix (cells = columns, genes/regions = rows)
#' @param frequencies logical, divide bmat by colSums (if FALSE simply use bmat for TF matrix)
#' @param log_scale_tf logical, take log1p of the term-frequencies?
#' @param scale_factor size factor to multiply with (if converting to frequencies)
#' @param threshold_for_idf count threshold to convert data to binary. Use for IDF
#'                          (IDF is not be used with non-binary data)
#'
#' @return tf-idf (sparse) matrix
#' @export
#'
#' @examples
#'
tfidf = function(bmat, frequencies=TRUE, log_scale_tf=TRUE, scale_factor=100000, threshold_for_idf = 1) {
  # Use either raw counts or divide by total   counts in each cell
  if (frequencies) {
    # "term frequency" method
    tf = Matrix::t(Matrix::t(bmat) / Matrix::colSums(bmat))
  } else {
    # "raw count" method
    tf = bmat
  }
  # check/convert the tf to sparse matrix
  if(class(tf) == "matrix" | attr(class(tf), "package") != "Matrix") {
    tf <- Matrix::Matrix(tf, sparse = TRUE)
  }
  # Either TF method can optionally be log scaled
  if (log_scale_tf) {
    if (frequencies) {
      tf@x = log1p(tf@x * scale_factor)
    } else {
      tf@x = log1p(tf@x * 1)
    }
  }

  # IDF w/ "inverse document frequency smooth" method
  ## in order to make this useful to non-binarized matrix, I have to reduce this to binarized matrix
  if(is.null(threshold_for_idf)) {
    idf = log(1 + ncol(bmat) / Matrix::rowSums(bmat))
  } else {
    idf = log(1 + ncol(bmat) / Matrix::rowSums(bmat >= threshold_for_idf))
  }

  # TF-IDF
  tf_idf_counts = safe_tfidf_multiply(tf, idf)
  tf_idf_counts@x[!is.finite(tf_idf_counts@x)] <- 0
  rownames(tf_idf_counts) = rownames(bmat)
  colnames(tf_idf_counts) = colnames(bmat)
  return(tf_idf_counts)
}

#' Run fast PCA (irlba) on a feature*sample matrix
#' @param mat feature*sample matrix (eg. rows= genes/regions, cols=cells), counts must already be normalized and scaled
#' @param dims how many dims to calculate?
#'
#' @return list of feature*PC and cell\*PC output matrices
#' @export
#'
#' @examples
#'
do_pca <- function(mat, dims=50) {
  mat@x[!is.finite(mat@x)] <- 0
  pca.results = irlba::irlba(Matrix::t(mat), nv=dims)
  final_result = pca.results$u %*% diag(pca.results$d)
  rownames(final_result) = colnames(mat)
  colnames(final_result) = paste0('PC_', 1:dims)
  return(final_result)
}

#' Quick UMAP and louvain cluster of a feature*PC data frame
#'
#' @param x data frame with rownames=features/samples, colnames=PCs
#' @param PCs vector of which principal components (columns) to select
#' @param nK No. of neighbors for louvain clustering
#' @param spread UMAP spread
#' @param mindist UMAP min_dist
#'
#' @return data.frame with UMAP1/2 and louvain clusters
#' @export
#'
#' @examples
#'
quick_umap_louvain <- function(x, PCs=2:50, nK=30, spread=5, mindist=0.1){

  lsa_umap <- x[,PCs] %>%
    uwot::umap(spread = spread, min_dist = mindist,
               n_neighbors = nK, scale = "none") %>%
    as.data.frame() %>%
    magrittr::set_colnames(c("UMAP1", "UMAP2")) %>%
    magrittr::set_rownames(rownames(x))

  lsa_umap$louvain <- x[,PCs] %>% cosine_similarity() %>%
    makeKNNgraph(nk = nK, replace_dist = FALSE, method = "RANN") %>%
    igraph::cluster_louvain() %>% .$membership
  return(lsa_umap)
}


#' Wrapper for LSA
#'
#' @param counts a count matrix of genes/regions (rows) * cells(columns), or a tibble of "gene", "cell" and "count"
#' @param nPC No. of Dims to reduce the data to
#' @param nK No. of neighbours louvain clustering and UMAP
#' @param outPdf File name to save the output UMAP
#'
#' @return list object with pca and umap output for both cells and regions.
#' @export
#'
#' @examples
#'

lsa_wrapper <- function(counts, nPC, nK, umap_spread=5, umap_mindist=0.1, outPdf=NA) {

  if(tibble::is_tibble(counts)) {
    counts <- reshape2::dcast(counts, gene ~ cell, fun.aggregate = sum) %>%
              column_to_rownames("gene") %>% as.matrix() %>%
              Matrix::Matrix(sparse = TRUE)
  }

  counts_tfidf <- tfidf(counts)
  pca.results <- irlba::irlba(Matrix::t(counts_tfidf), nv=nPC)
  region_pca <- pca.results$v %*% diag(pca.results$d)
  cells_pca <- pca.results$u %*% diag(pca.results$d)
  rownames(region_pca) <- rownames(counts_tfidf)
  rownames(cells_pca) <- colnames(counts)
  colnames(region_pca) <- colnames(cells_pca) <- paste0('PC_', 1:nPC)

  pca_list <- list(region = region_pca, cells = cells_pca)
  umap_out <- quick_umap_louvain(pca_list$cells, PCs=2:nPC, nK=nK, spread=umap_spread, mindist=umap_mindist)

  p <- ggplot(umap_out, aes(UMAP1, UMAP2, col = as.factor(louvain))) +
       geom_point() + theme_grey(base_size = 16) +
       labs(col = "Louvain")

  if(!is.na(outPdf)) ggsave(outPdf, p)
  lsa_out <- list(pca = pca_list, umap = umap_out)
  return(lsa_out)
}

