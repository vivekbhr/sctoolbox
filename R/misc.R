#' Quick look into a large matrix/df
#'
#' @param m matrix/df
#'
#' @return matrix with 10 rows and cols
#' @export
#'
#' @examples
#'
look <- function(m) {
  m[1:10, 1:10]
}


#' Return n distinct color as hex code
#'
#' @param n
#'
#' @return
#' @export
#'
#' @examples
#'
#'
load_colorList <- function(n) {
    p1 <- c(
      '#e0346c',
      '#29a5db',
      '#682a4d',
      '#2a23e5',
      '#584880',
      '#69e186',
      '#a4d352',
      '#93be65',
      '#c80b54',
      '#695c20',
      '#8d7e8c',
      '#8a2eca',
      '#ff4752',
      '#ab3b94',
      '#7b7a78')
    return(p1[1:n])
  }

#' Load commonly used packages
#'
#' @param category choose from "general", "topicmodels", "stats",
#'                 "granges", "rna" or "all"
#'
#' @return Nothing. Just loads packages
#' @export
#'
#' @examples
#'
loadPackages <- function(category = "general") {
  if (category %in% c("general", "all")) {
    # general
    library(magrittr) # for %$% operator
    library(ggplot2)
    library(dplyr)
  }

  if (category %in% c("topicmodels", "all")) {
    # LDA
    library(tidytext)
    library(ldatuning)
    library(topicmodels)
  }

  if (category %in% c("stats", "all")) {
    # clustering and vis
    library(Matrix)
    library(irlba)
    library(uwot)
    library(fnn)
    library(matrixcalc)
    library(igraph)
  }

  if (category %in% c("granges", "all")) {
    # GRanges
    library(GenomicAlignments)
    library(GenomicFeatures)
    library(Rsamtools)
    library(BiocParallel)
    library(biomaRt)
  }

  if (category %in% c("rna", "all")) {
  # RNAseq
  library(Seurat)
  }

}


#' Remove version name from ensembl gene IDs
#'
#' @param geneIDs char vector: ensembl gene IDs
#'
#' @return char vecter of santized gene IDs
#' @export
#'
#' @examples
#'

sanitizeGeneID <- function(geneIDs) {
  g <- gsub("([a-z|A-Z][0-9]*)\\.[0-9]*", "\\1" , geneIDs)
  return(g)
}

