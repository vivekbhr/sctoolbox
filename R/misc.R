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


#' convert char to GRanges
#'
#' @param char "chr_start_end"
#'
#' @return
#' @export
#'
#' @examples
#'
charToGRanges <- function(char) {
  chrom <- gsub("(.*)_([0-9]*)_([0-9]*)", "\\1", char)
  start <- gsub("(.*)_([0-9]*)_([0-9]*)", "\\2", char)
  end <- gsub("(.*)_([0-9]*)_([0-9]*)", "\\3", char)
  GenomicRanges::GRanges(chrom, IRanges::IRanges(as.numeric(start), as.numeric(end)))
}


#' prepare .txt file to create pseudobulk bigwigs from clustered UMAP (using snakemake workflow)
#'
#' @param umap umap data frame, from standard LSA/LDA wrapper
#' @param prefix file name prefix (match it with bam file prefix to run snakemake afterwards)
#' @param suffix file name suffix (match it with bam file suffix to run snakemake afterwards)
#'
#' @return
#' @export
#'
#' @examples
#'
#'
prepForPseudobulk <- function(umap, prefix, suffix = ".txt") {
  umap$plate <- gsub("(.*)_([ATGC]*)", "\\2", rownames(umap))

  lapply(unique(umap$plate), function(x){
    df <- umap %>% rownames_to_column("cell") %>% as_tibble() %>% filter(plate == x) %>%
      mutate(cell = gsub("(.*)_([ATGC]*)", "\\2", cell)) %>% dplyr::select(cell, louvain)
    write.table(df, file = paste0(prefix, x, suffix), sep = "\t",
                quote = F, col.names = F, row.names = F)

    })

}
