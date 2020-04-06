
## calculate similarity matrix (using cosine similarity)

#' calculate cell*cell cosine similarity
#'
#' @param matrix input matrix (rows = cells, column = genes)
#'
#' @return cells*cells matrix
#' @export
#'
#' @examples
#'

cosine_similarity <- function(matrix) {
  sim <- matrix / sqrt(rowSums(matrix * matrix))
  sim <- sim %*% t(sim)
  return(sim)
}


#' Convert a cellxPC or a cellxTopic matrix to a KNN matrix
#'
#' @param mat Input matrix after dim reduction (cell x PC or cell x Topic)
#' @param nk no. of K
#' @param replace_dist add the distance back to the KNN graph, if TRUE, the distance
#'                     to the K-nearest neighbours are retained in the output graph (weighted)
#'                     otherwise the distances are removed, returning unweighted graph
#' @param method KNN building method (RANN or fnn)
#'
#' @return igraph objet
#' @export
#'
#' @examples
#'
#'

makeKNNgraph <- function(mat, nk, replace_dist = FALSE, method = "RANN") {

  if(method == "fnn") {
    knn.info<- FNN::get.knn(mat, k = nk) ## doesn't make much diff
    knn <- knn.info$nn.index
  } else {
    knn.info<- RANN::nn2(mat, k = nk)
    knn <- knn.info$nn.idx
  }

  ## convert to adjacancy matrix
  adj <- matrix(0, nrow(mat), nrow(mat))
  rownames(adj) <- colnames(adj) <- rownames(mat)

  if(isTRUE(replace_dist)) {
    for(i in seq_len(nrow(mat))) {
      adj[i,rownames(mat)[knn[i,]]] <- mat[knn[i,]]
    }
    graph_obj <- igraph::graph.adjacency(adj, weighted = TRUE, mode = "undirected")
  } else {
    for(i in seq_len(nrow(mat))) {
      adj[i,rownames(mat)[knn[i,]]] <- 1
    }
    graph_obj <- igraph::graph.adjacency(adj, mode = "undirected")
  }

  return(graph_obj)
}



#' Get top-genes per cluster based on T-test
#'
#' @param cluster_df Data frame with 2 columns "cell ID", "cluster no."
#' @param count_matrix Wide format count matrix with rownames = genes, colnames = cells
#' @param n No. of top genes to extract (NA = return df with all significant genes at padj < 0.05)
#'
#' @return vector of top genes, or (if n = NA) data frame with t-test output
#' @export
#'
#' @examples
#'

topGenesByCluster <- function(cluster_df, count_matrix, n = NA) {
  colnames(cluster_df) <- c("cells", "cluster")
  cl <- unique(cluster_df$cluster)
  lapply(cl, function(x){
    cells_test <- cluster_df[cluster_df$cluster == x, "cells"]
    cells_control <- cluster_df$cells[!(cluster_df$cells %in% cells_test)]
    out <- apply(count_matrix, 1, function(c) t.test(c[cells_test], c[cells_control], alternative = "greater")$p.value)
    outdf <- data.frame(pvalue = out)
    outdf$padj <- p.adjust(outdf$pvalue, method = "bonferroni")
    outdf <- outdf[order(outdf$padj), ]
    return(outdf)
    }) -> final
  names(final) <- paste0("cluster_", cl)

  if (!is.na(n)) {
    final <- lapply(final, function(x) return(rownames(x[1:n,])) )
  } else {
    final <- lapply(final, function(x) return(x[x$padj < 0.05,]))
  }

  return(final)
}




#' Plot the counts for gene lists to compare bw clusters
#'
#' @param geneList list with names = clusterIDs and values = geneIDs
#' @param counts count tibble (gene, cell, count)
#' @param umap umap output (rownames=  cells)
#' @param cl_col which column contain cluster info in the umap?
#'
#' @return ggplot objects
#' @export
#'
#' @examples
#'
#'
plotGenesPerCluster <- function(geneList, counts, umap, cl_col = "louvain") {
  lapply(names(geneList), function(cl) {
    df2 <- dplyr::filter(counts, gene %in% geneList[[cl]])
    df2$cluster <- umap[match(df2$cell, rownames(umap)), cl_col]
    ggplot(df2, aes(as.factor(cluster), count, group = cluster)) +
      geom_boxplot(position = "dodge") + geom_jitter(alpha = 0.2) +
      scale_y_log10() + labs(title = cl)
  })
}



#' PLot topic probability on UMAP
#'
#' @param ldaOut output of lda Wrapper
#' @param nTopic nTopics to use
#' @param outPdf output pdf file
#'
#' @return pdf file with ggplot of topic probs
#' @export
#'
#' @examples
#'
plotTopicProbs <- function(ldaOut, nTopic = 20, outPdf) {
  lda_umap <- ldaOut$umap[,1:2]
  lda_gamma <- tidyr::spread(ldaOut$gamma, topic, gamma)
  lda_umap %<>% merge(lda_gamma, by.x = 0, by.y = 1)
  colnames(lda_umap)[4:(nTopic+3)] <- paste0("topic_", 1:nTopic)

  pdf(outPdf)
  lapply(paste0("topic_", 1:nTopic), function(x){
  print(ggplot(lda_umap, aes_string("UMAP1", "UMAP2", fill = x)) +
          geom_point(shape = 21) + theme_grey(base_size = 16) +
          ggtitle(x) +
          labs(fill = "Probability") +
          scale_fill_gradient(low = "white", high = "black"))
  })
  dev.off()
}

