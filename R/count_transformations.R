
#' Normalize a single-cell count matrix based on DEseq sizefactor
#'
#' @param m matrix (cells in cols, genes in rows)
#'
#' @return matrix (normalized)
#' @export
#'
#' @examples
#'
deseq_norm <- function(m)
{
  ncols <- dim(m)[2]
  condition <- c(0, rep(1, ncols-1))
  m <- round(m)
  colData <- data.frame(row.names=colnames(m), condition = factor(condition))
  dds <- DESeq2::DESeqDataSetFromMatrix(countData=m, colData=colData, design = ~condition)
  dds <- DESeq2::estimateSizeFactors(dds)
  sizeFactors <- dds$sizeFactor
  m <- sweep(m, 2, sizeFactors, '/')
  m
}




#' Obtain average CPM per group from a SingleCellExperiment object
#'
#' @param sce SingleCellExperiment object
#' @param group_var variable in colData(sce) that can be used to group cells
#' @param assays slot from assayNames(sce) which we want to transform
#' @param aggtype type of aggregation (sum/mean)
#'
#' @return dataframe with genes/regions in rows and groups in columns
#' @export
#'
#' @examples
#'
cpm_bygroup <- function(sce, group_var, assays, aggtype='mean') {

  group_names <- unique(colData(sce)[[group_var]])

  if(aggtype=="mean"){
    group_counts <- plyr::ldply(group_names, function(name){
      sce_subset <- sce[,colData(sce)[[group_var]]==name]

      count_list <- lapply(assays, function(aname){
        assay(sce_subset, "counts") <- assay(sce_subset, aname)
        counts <- edgeR::aveLogCPM.SummarizedExperiment(sce_subset)
        return(counts)
      })
      names(count_list) <- assays

      out <- data.frame(var = rownames(sce_subset),
                        group = name,
                        count_list)

      return(out)
    })

  } else if(aggtype=="sum") {

    group_counts <- plyr::ldply(assays, function(aname){

      count_list <- plyr::ldply(group_names, function(name){
        sce_subset <- sce[,colData(sce)[[group_var]]==name]
        sums <- Matrix::rowSums(assay(sce_subset, aname))
        return(sums)
      })# returns groups=rows, genes=columns
      rownames(count_list) <- group_names

      dge <- edgeR::DGEList(counts = t(count_list), group = 1:nrow(count_list))
      out <- edgeR::cpm(dge)#log2()?
      out <- reshape2::melt(out)
      out$assay <- aname
      return(out)
    })
    group_counts %<>% tidyr::pivot_wider(names_from = assay, values_from = value)
    colnames(group_counts)[1:2] <- c("var", "group")
  } else {
    message("Aggtype not supported")
  }

  return(group_counts)
}



#' Get transcript per 10k counts per group
#'
#' @param sce SingleCellExperiment object, with cell groups and lengths, counts should be raw
#' @param group_var column in colData(sce) to make cell groups
#' @param length_var column in rowData(sce) which indicate gene lengths
#' @param assays which assays (from assayNames(sce)) to use
#' @param aggtype either "sum" (sum by group then take tpt), or "mean" (take tpt per cell then mean by group)
#'
#' @return data.frame (rows=genes, col = groups)
#' @export
#'
#' @examples
#'
tpt_bygroup <- function(sce, group_var, length_var, assays, aggtype="sum") {

  group_names <- unique(colData(sce)[[group_var]])

  if(aggtype=="sum") {
    # sum across cells per group, then take tpt
    group_counts <- plyr::ldply(assays, function(aname){
      count_list <- plyr::ldply(group_names, function(name){
        sce_subset <- sce[,colData(sce)[[group_var]]==name]
        sums <- Matrix::rowSums(assay(sce_subset, aname))
        return(sums)
      })# returns groups=rows, genes=columns
      rownames(count_list) <- group_names
      return(t(count_list))
    })

    tpt <- apply(group_counts, 2, countToTpm, effLen = rowData(sce)[,length_var], scale_factor=1e04)
    rownames(tpt) <- rownames(sce)

  } else if(aggtype=="mean") {
    # take tpt first, then mean them across cells per group
    tpt <- plyr::ldply(assays, function(aname) {
      sce_tpt <- apply(assay(sce, aname), 2, countToTpm, effLen = rowData(sce)[,length_var], scale_factor=1e04)
      # avg tpt per group
      count_list <- plyr::ldply(group_names, function(name){
        tpt_subset <- sce_tpt[,colData(sce)[[group_var]]==name]
        avg <- rowMeans(tpt_subset)
        return(avg)
      })# returns groups=rows, genes=columns

      rownames(count_list) <- group_names
      return(t(count_list))
    } )
    rownames(tpt) <- rownames(sce)
  }

  return(tpt)
}




#' Collate binCounts at gene level based on overlap
#'
#' @param mat Count matrix (bin-level) with regions as rownames (format: chr_start_end)
#' @param genes_gr Genes as GRanges object, column "gene_id" should indicate gene names/IDs
#' @param subset List of genes to subset from genes_gr for overlap
#' @param overlapGap Max distance between gene and bin for overlap
#'
#' @return Sparse matrix (dgcMatrix) with collates counts (genes=rows, cells=columns)
#' @export
#'
#' @examples
#'

getGeneActivity <- function(mat, genes_gr, subset=NULL, overlapGap=10000){

  gr <- sctoolbox::charToGRanges(rownames(mat))
  if(!is.null(subset)){
    gr_gene <- genes_gr[subset]
  } else {
    gr_gene <- genes_gr
  }

  gr_gene <- gr_gene[!duplicated(gr_gene$gene_id)]
  if(length(gr_gene)==0){
    warning("No gene found in GRanges")
  }
  # get peaks overlaping genes +- overlapGap
  ol <- findOverlaps(gr_gene, gr, maxgap = overlapGap)
  qhits <- unique(queryHits(ol))
  ol %<>% as.data.frame()
  gr_gene <- gr_gene[qhits]

  scores <- lapply(qhits, function(x){
    s_hits <- unique(ol[ol$queryHits == x, "subjectHits"])
    if(length(s_hits) > 1){
      return(colSums(mat[s_hits,]))
    } else {
      return(mat[s_hits, ])
    }
  }) %>% do.call(rbind, .)

  if(is.null(dim(scores)) | length(gr_gene)==0) {
    warning("scores are empty!")
    out <- 0
  } else {
    out <- Matrix::Matrix(scores, sparse = TRUE)
    rownames(out) <- gr_gene$gene_id
  }
  return(out)

}
