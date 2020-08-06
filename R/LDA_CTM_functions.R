#' Load params for LDA or CTM
#'
#' @param out which param to load ('LDA' or 'CTM')
#'
#' @return list of params
#' @export
#'
#' @examples
#'
loadLDAparams <- function(out){
  control_list_gibbs <- list(
    burnin = 250, #2500,
    iter = 500, #5000,
    seed = 0:4,
    nstart = 5,
    best = TRUE
  )

  control_list_ctm <- list(
    seed = 0:4,
    nstart = 5,
    best = TRUE
  )
  if(out == "lda") {
    return(control_list_gibbs)
  } else {
    return(control_list_ctm)
  }

}


#' Get the topic to document/cell (gamma) distribution from a CTM output
#'
#' @param CTM_object output of topicmodels::CTM function
#'
#' @return tibble with gamma values
#' @export
#'
#' @examples
#'
#'
tidy_ctm_gamma  <- function(CTM_object){
  CTM_object %>%
    slot("gamma")  %>%
    as_data_frame()  %>%
    mutate (document = row_number()) %>%
    gather(topic, gamma, -document) %>%
    mutate(topic = strtoi(stringr::str_sub(topic,2)))
}

#' Get the word/region to topic (beta) distribution from a CTM output
#'
#' @param CTM_object output of topicmodels::CTM function
#'
#' @return tibble with beta values
#' @export
#'
#' @examples
#'
tidy_ctm_beta  <- function(CTM_object){
  Terms  <- CTM_object %>%
    slot("terms")
  CTM_object %>%
    slot("beta")  %>%
    as_data_frame() %>%
    setNames(Terms) %>%
    mutate (topic = row_number()) %>%
    gather(term, beta, -topic) %>%
    mutate(beta = exp(beta))
}

#' Plot the topic to document/cell (gamma) distribution from LDA/CTM output
#'
#' @param tidymodel tidy output of LDA/CTM
#'
#' @return ggplot object
#' @export
#'
#' @examples
#'
plotGamma <- function(tidymodel) {
  tidymodel %>%
    group_by(document) %>%
    arrange(desc(gamma)) %>%
    slice(1) %>%
    ungroup() %>%
    ggplot(aes(x=gamma)) +
    geom_histogram(bins = 20) +
    xlab("maximum gamma per document") +
    geom_vline(aes(xintercept = 1/20),
               color="darkred")

}

## helper function to sort
## correlation matrix
reorder_cor <- function(x){
  ord <- corrplot::corrMatOrder(x)
  x[ord,ord]
}

## helper function to extract
## lower triangle of matrix
get_lower_tri<-function(x){
  x[upper.tri(x)] <- NA
  return(x)
}

#' Extract topic correlations from a CTM model output
#'
#' @param ctm_gamma gamma matrix output of CTM
#' @param ntopics number of topics
#' @param plot logical. Whether to plot topic*topic correlation across documents.
#'
#' @return data frame of correlations (+plot object if plot= TRUE)
#' @export
#'
#' @examples
#'
getTopicCorrelations <- function(ctm_gamma, ntopics, plot = TRUE) {
  ## get corr matrix (cor between topics)
  out <- ctm_gamma %>%
    pivot_wider(names_from = topic, values_from = gamma) %>%
    select(-document) %>%
    cor() %>%
    reorder_cor() %>%
    get_lower_tri() %>%
    as_data_frame() %>%
    mutate(topic1 = forcats::as_factor(paste(names(.)))) %>%
    gather(topic2, correlation, - topic1) %>%
    mutate(topic2 = factor(topic2, levels=levels(topic1)))

  ## plot corr matrix
  if(isTRUE(plot)) {
    p <- ggplot(out, aes(as.numeric(topic1), as.numeric(topic2), fill=correlation)) +
      geom_tile(color="white") +
      scale_fill_gradient2(low = "firebrick4", high = "dodgerblue4",
                           midpoint = 0, limit = c(-1,1), space = "Lab",
                           na.value = "white",
                           name="Pearson\nCorrelation") +
      scale_x_continuous(
        breaks=1:ntopics, labels = levels(cor_data$topic1), name="")+
      scale_y_continuous(
        breaks=1:ntopics, labels = levels(cor_data$topic2), name="")+
      theme_minimal()+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
      coord_fixed()
  }
  if (isTRUE(plot)) {
    return(list(out, p))
  } else {
    return(out)
  }
}


#' Wrapper for LDA, UMAP and clustering
#'
#' @param geneCounts geneCounts in long format (gene, cell, counts)
#' @param nTopic No. of topics for LDA
#' @param nK No. of nearest neighbours for louvain and UMAP
#'
#' @return list with LDA and umap output
#' @export
#'
#' @examples
#'
#'
ldaWrapper <- function(geneCounts, binarize = FALSE, nTopic = 20, nK = 30) {

  ### LDA
  if(isTRUE(binarize)) {
    geneCounts[geneCounts$count > 1,]$count <- 1
  }
  counts_dtm <- geneCounts %>% tidytext::cast_dtm(document = cell, term = gene, value = count)
  ldaparams <- loadLDAparams(out = "lda")
  ldaOut <- list()
  ldaOut$dtm <- counts_dtm
  ldaOut$lda <- topicmodels::LDA(k = nTopic, x=counts_dtm, method="Gibbs", control=ldaparams)
  ldaOut$gamma <- tidytext::tidy(ldaOut$lda, matrix = "gamma")
  ldaOut$beta <- tidytext::tidy(ldaOut$lda, matrix = "beta")

  ## umap
  lda_gamma <- tidyr::spread(ldaOut$gamma, topic, gamma)
  ldaOut$umap <- lda_gamma[2:(nTopic+1)] %>% as.matrix() %>%
    uwot::umap(spread = 5, min_dist = 0.1,
               n_neighbors = nK, scale = "none") %>%
    as.data.frame() %>%
    magrittr::set_colnames(c("UMAP1", "UMAP2")) %>%
    magrittr::set_rownames(lda_gamma$document)
  ## clustering
  ldaOut$umap$louvain <- lda_gamma %>% as.matrix() %>%
    set_rownames(.[,1]) %>% .[,2:(nTopic+1)] %>%
    makeKNNgraph(nk = nK, replace_dist = TRUE, method = "RANN") %>%
    igraph::cluster_louvain() %>% .$membership

  ldaOut$umap$hdbscan <- dbscan::hdbscan(ldaOut$umap[1:2], 10)$cluster
  ## plot
  ggplot(ldaOut$umap, aes(UMAP1, UMAP2, col = as.factor(louvain))) +
    geom_point() + theme_grey(base_size = 16) +
    labs(col = "Louvain") + scale_color_brewer(palette = "Set1")

  ## return
  return(ldaOut)
}




#' Count filters and plots
#'
#' @param counts long format count matrix (with columns: cell, gene, count)
#' @param minGeneCount min Counts for a gene to be kept
#' @param maxGeneCount max Counts for a gene to be kept (remove outliers)
#' @param minCellCount min Counts for a Cell to be kept
#' @param minCellsForGenes to filter genes, min number of cells with nonzero counts
#' @param minGenesForCells to filter cells, min number of genes with nonzero counts
#'
#' @return a list, with filtered geneCounts in tibble (cell,gene, count), and wide matrix (genes*cells)
#' @export
#'
#' @examples
#'
countFilter <- function(counts, minGeneCount = 0, maxGeneCount = Inf,
                        minCellCount = 100, minCellsForGenes = 0, minGenesForCells = 0) {

  ## filter
  cell_sums <- counts %>% group_by(cell) %>% dplyr::summarise(counts = sum(count))
  gene_sums <- counts %>% group_by(gene) %>% dplyr::summarise(counts = sum(count))
  genes_tokeep <- gene_sums[gene_sums$counts >= minGeneCount & gene_sums$counts <= maxGeneCount, ]$gene
  cells_tokeep <- counts %>% dplyr::group_by(cell) %>% dplyr::summarise(counts = sum(count)) %>%
    filter(counts >= minCellCount) %>% .$cell

  print(paste0("Cells kept at minCellCount = ", minCellCount, " : ", length(cells_tokeep)))
  print(paste0("Genes kept at minGeneCount = ", minGeneCount, " and maxGeneCount = ",
               maxGeneCount, " : ", length(genes_tokeep)))

  # make a histogram with filter thresholds
  pl1 <- ggplot(gene_sums, aes(counts)) + geom_histogram() + scale_x_log10() +
    labs(x = "Total counts across cells", y = "nGenes")
  pl2 <- ggplot(cell_sums, aes(counts)) + geom_histogram() + scale_x_log10() +
         labs(x = "Total counts across regions", y = "nCells") #+ geom_vline(minCellCount)
  print(gridExtra::grid.arrange(pl1, pl2))

  counts %<>% dplyr::filter(gene %in% genes_tokeep)
  counts %<>% dplyr::filter(cell %in% cells_tokeep)

  ## Filter
  if(minCellsForGenes != 0 | minGenesForCells != 0) {
    p1 <- counts %>% group_by(cell) %>% summarise(counts = sum(count > 0))
    p2 <- counts %>% group_by(gene) %>% summarise(counts = sum(count > 0))
    keep <- list(cells = p1[p1$counts >= minGenesForCells, "cell"],
                 genes = p2[p2$counts >= minCellsForGenes, "gene"])
    print(paste0("Cells kept at minGenesForCells = ", minGenesForCells, " : ", nrow(keep$cells)))
    print(paste0("Genes kept at minCellsForGenes = ", minCellsForGenes, " : ", nrow(keep$genes)))

    # make a histogram with filter thresholds
    pl1 <- ggplot(p1, aes(counts)) + geom_histogram() + scale_x_log10() +
      labs(x = "Total No. of regions detected", y = "nCells") #+
      #geom_vline(minGenesForCells)
    pl2 <- ggplot(p2, aes(counts)) + geom_histogram() + scale_x_log10() +
      labs(x = "Total No. of cells detected", y = "nGenes") #+
      #geom_vline(minCellsForGenes)
    print(gridExtra::grid.arrange(pl1, pl2))

    counts %<>% dplyr::filter(cell %in% keep$cells$cell, gene %in% keep$genes$gene)
  }

  ## convert to wide
  counts_wide <- counts %>% tidyr::spread(cell, count, fill = 0)
  names <- counts_wide$gene
  counts_wide <- as.matrix(counts_wide[,2:ncol(counts_wide)])
  rownames(counts_wide) <- names

  return(list(counts_tibble = counts, counts_mat = counts_wide))

  }
