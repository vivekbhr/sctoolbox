
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
#' @return tf-idf matrix
#' @export
#'
#' @examples
#'
tfidf = function(bmat, frequencies=TRUE, log_scale_tf=TRUE, scale_factor=100000, threshold_for_idf = 1) {
  # Use either raw counts or divide by total counts in each cell
  if (frequencies) {
    # "term frequency" method
    tf = Matrix::t(Matrix::t(bmat) / Matrix::colSums(bmat))
  } else {
    # "raw count" method
    tf = bmat
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
  rownames(tf_idf_counts) = rownames(bmat)
  colnames(tf_idf_counts) = colnames(bmat)
  return(tf_idf_counts)
}

#' Run fast PCA (irlba) on TF-IDF matrix
#'
#' @param mat matrix
#' @param dims how many dims to calculate?
#'
#' @return PCA output matrix
#' @export
#'
#' @examples
#'
do_pca = function(mat, dims=50) {
  pca.results = irlba::irlba(Matrix::t(mat), nv=dims)
  final_result = pca.results$u %*% diag(pca.results$d)
  rownames(final_result) = colnames(mat)
  colnames(final_result) = paste0('PC_', 1:dims)
  return(final_result)
}
