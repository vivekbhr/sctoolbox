#' Quick look into a large matrix/df
#'
#' @param m matrix/ dataframe
#'
#' @return matrix with 10 rows and cols
#' @export
#'
#' @examples
#'
look <- function(m) {
  m[1:10, 1:10]
}

make_char <- function(...){
  as.character(match.call()[-1])
}


#' Return hex codes for distinct colors mapped to input names
#'
#' @param names vector of names (each unique name is mapped to a hex code)
#'
#' @return clist vector of colors which can be used in ggplot2::scale_color_manual
#' @export
#'
#' @examples
#'
#'
get_colors <- function(names, color_space=NA) {

  scheme <- rwantshue::iwanthue(seed = 0, force_init = TRUE) # recreate with a seed
  if(is.na(color_space)) {
    color_space <- list(
      c(0, 360),	# hue range [0,360]
      c(0, 80),		# chroma range [0,100]
      c(0, 100))		# lightness range [0,100]
  }

  alis <- unique(names)
  clist <- scheme$hex(length(alis), force_mode = FALSE, quality = 50, color_space = color_space)
  names(clist) <- alis

  return(clist)
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

#' Convert counts to TPM (or TP<scale_factor>)
#'
#' @param counts Count matrix of genes*samples
#' @param effLen vector of effective gene lengths per gene
#'
#' @return matrix of TPM normalized counts
#' @export
#'
#' @examples
#'
countToTpm <- function(counts, effLen, scale_factor=1e6) {
  rate <- log(counts+1) - log(effLen+1)# counts/length
  denom <- log(sum(exp(rate), na.rm = T))# sum(counts)
  exp(rate - denom + log(scale_factor))#log([counts/length] / [sum(counts)/scale_factor])
}

#' Convert counts to FPKM
#'
#' @param counts Count matrix of genes*samples
#' @param effLen vector of effective gene lengths per gene
#'
#' @return matrix of RPKM normalized counts
#' @export
#'
#' @examples
#'
countToFpkm <- function(counts, effLen) {
  N <- sum(counts)
  exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}


#' Convert FPKM to TPM
#'
#' @param fpkm vector of FPKM values
#'
#' @return
#' @export
#'
#' @examples
#'
fpkmToTpm <- function(fpkm) {
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

#' Convert counts to Effective counts
#'
#' @param counts vector of counts per gene
#' @param len actual length per gene
#' @param effLen effective length per gene
#'
#' @return
#' @export
#'
#' @examples
#'
countToEffCounts <- function(counts, len, effLen) {
  counts * (len / effLen)
}
