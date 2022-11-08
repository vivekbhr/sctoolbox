## Functions to do logistic regression on transcript counts

#' fit a linear model and do a likelihood ratio test
#'
#' @param counts vector of counts
#' @param cell_group vector with cell groups
#' @param cutoff count cutoff
#'
#' @return pvalue of LR test
#' @export
#'
#' @examples
#'
#'
LR <- function(counts, cell_group, cutoff = 0.9 * nrow(counts))

  {
  if(is.vector(counts))
  {
    if(sum(counts  == 0) > cutoff)
    {
      return(NA)
    }
  }
  else{
    zeros <- apply(counts, 2, function(x) sum(x ==0) > cutoff)
    if(sum(zeros) == ncol(counts))
    {
      return(NA)
    }
    counts <- counts[,!zeros]
  }
  table <- data.frame(counts, cell_group)
  full <-  glm(cell_group~., data=table, family=binomial)
  null <-  glm(cell_group~1, data=table, family=binomial)
  lrtest <- lmtest::lrtest(full, null)
  lrtest$Pr[2]
}

#' Gene-wise likelihood ratio test
#'
#' @param counts tibble of counts (gene, cell, counts)
#' @param genes vector with gene names (for testing)
#' @param cell_group vector with cell groups (cluster labels)
#' @param filter_threshold lower count threshold
#'
#' @return data frame with pvalue per gene
#' @export
#'
#' @examples
#'
#'
LR_fit <- function(counts, genes, cell_group, filter_threshold)
{
  table <- data.frame(genes = genes, index = 1:length(genes))
  indices <- dplyr::summarise(dplyr::group_by(table, genes), indices = list(index))
  cutoff <- filter_threshold * nrow(counts)
  pvalues <- sapply(indices$indices, function(x) {
    LR(counts[,x], cell_group, cutoff)
  })
  n_transcripts <- sapply(indices$indices, length)
  data.frame(pvalues, genes =  indices$genes, n_transcripts)
}



