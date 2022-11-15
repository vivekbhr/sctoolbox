#' filter Reads based on nucleotide content at 5-prime
#'
#' @param bamdf bam file as data frame
#' @param pos_strand_expected expected Nucleoties in read and genome for alignments in + strand c("A", "AT")
#' @param neg_strand_expected expected Nucleoties in read and genome for alignments in - strand c("A", "AT")
#' @param plotFile file to save the plot of nucleotide content
#' @param bsgenome a BSgenome object
#'
#' @return vector with T/F with filtering result for each read
#' @export
#'
#' @examples
#'
#'
#'

filterReadMotifs <- function(bamdf,
                        pos_strand_expected,
                        neg_strand_expected,
                        plotFile = NA,
                        bsgenome) {

  bamdf$pos <- as.numeric(bamdf$pos)
  bamdf$ref_width <- cigarWidthAlongReferenceSpace(bamdf$cigar)
  bamdf$ref_pos <- ifelse(bamdf$strand == "+", bamdf$pos,
                          bamdf$pos + bamdf$ref_width)
  bamdf$firstbase <- substr(bamdf$seq, 1, 1)
  bamdf$lastbase <- substr(bamdf$seq, bamdf$qwidth, bamdf$qwidth)
  bamdf$query_base <- ifelse(bamdf$strand == "+", bamdf$firstbase,
                             bamdf$lastbase)
  gr <- GRanges(bamdf$rname, IRanges(bamdf$ref_pos-1, bamdf$ref_pos),
                "+")#bamdf$strand
  bamdf$ref_bases <- as.character(
    suppressWarnings(getSeq(bsgenome, gr)))
  keep <- ifelse(bamdf$strand == "+",
                 ifelse(bamdf$query_base == pos_strand_expected[1] &
                          bamdf$ref_bases %in% pos_strand_expected[2:length(pos_strand_expected)], TRUE, FALSE),
                 ifelse(bamdf$query_base == neg_strand_expected[1] &
                          bamdf$ref_bases %in% neg_strand_expected[2:length(neg_strand_expected)], TRUE, FALSE)
  )
  ### make plot
  if(!(is.na(plotFile))) {

    split(bamdf, bamdf$strand) %>%
      lapply(function(x) table(x$ref_bases) %>% as.data.frame) %>%
      plyr::ldply(data.frame) -> refpos
    split(bamdf, bamdf$strand) %>%
      lapply(function(x) table(x$query_base) %>% as.data.frame) %>%
      plyr::ldply(data.frame) -> querypos

    p1 <- ggplot(querypos, aes(Var1, Freq)) +
      geom_bar(stat = "identity", position = "dodge") + facet_grid(.id~.) +
      labs(x = "Base", y = "Reads (from 1Mil)", title = "1st base in R1") +
      theme_gray(base_size=14) + theme(axis.text.x=element_text(face = "bold"))

    p2 <- ggplot(refpos, aes(Var1, Freq)) +
      geom_bar(stat = "identity", position = "dodge") + facet_grid(.id~.) +
      labs(x = "Base", y = "Reads (from 1Mil)", title = "Genomic context")+
      theme_gray(base_size=14) + theme(axis.text.x=element_text(face = "bold"))

    final_p <- p1 + p2
    ggsave(plot = final_p, filename=plotFile, height=6, width=12)

  }

  return(keep)
}


#' Get counts for given motif at 5'-cut-site of MNAse
#'
#' @param plus_strand_expected vector. expected motif at forward-oriented reads and their genomic location
#' @param minus_strand_expected vector. expected motif at forward-oriented reads and their genomic location
#' @param plot fileName to plot
#' @param bsGenome bsGenome object
#' @param bamfile bam file object (using Rsamtools::bamFile)
#'
#' @return numeric with no. of reads passing the filter
#' @export
#'
#' @examples
#'
#'
#'

get5pNucCounts <- function(bamfile,
                           plus_strand_expected = c('A', 'TA'),
                           minus_strand_expected = c('T', 'TA'),
                           plot = NA, bsGenome = NA) {

  bcCounts <- parallel::mclapply(barcodes$V1, function(bc){
    bfreads <- Rsamtools::scanBam(bamfile,
                                  param=ScanBamParam(what=c("rname", "strand", "pos", "qwidth", "cigar", "seq"),
                                                     tag = "BC", tagFilter = list("BC" = bc)) )
    df <- as.data.frame(bfreads[[1]])
    filt <- filterReadMotifs(df[df$rname %in% standardChroms, ],
                        pos_strand_expected = plus_strand_expected,
                        neg_strand_expected = minus_strand_expected,
                        plotFile = plot,
                        bsgenome = bsGenome)
    return(sum(filt))
  }, mc.cores=20)

  return(bcCounts)
}

## get fragment length by looking at +/- strand
getFragLen <- function(x, bam, regions = NA){
  bfreads <- readGAlignments(bam,param=ScanBamParam(tag = "BC", tagFilter = list("BC" = x)) )
  if(!is.na(regions)) {
    bfreads %<>% subsetByOverlaps(ctcf)
  }
  bfreads %<>% split(strand(bfreads))
  bfreads %<>% lapply(GRanges)
  bfreads$`+` %<>% resize(1)
  bfreads$`-` %<>% resize(1, fix="end")
  distanceToNearest(bfreads$`+`, bfreads$`-`, ignore.strand = TRUE) %>%
    as.data.frame() -> d
  return(d[d$distance <= 1000, "distance"])

}

## bin the fragment lengths vector
binFragLen <- function(x, bc){
  names(x) <- bc
  x %<>% plyr::ldply(data.frame)
  colnames(x) <- c("barcodes", "hist")
  ## define breaks
  breaks = seq(0, 1000, 5)
  x$group_tags <- cut(x$hist,
                      breaks=breaks,
                      include.lowest=TRUE,
                      right=FALSE)

  x %<>% dplyr::group_by(barcodes, group_tags) %>% dplyr::summarise(n = n())
  return(x)
}


#' Get distance between R1-cuts on opposite strand per cell
#'
#' @param bam BAM file
#' @param barcodes barcodes (chr vector)
#' @param subsetRegions GRanges to subset to (if needed)
#' @param nCores Int, No. of cores
#' @param linePlot chr, name of line plot
#' @param heatmap chr, name of heatmap
#'
#' @return data.frame with barcode and distance between cuts
#' @export
#'
#' @examples
#'
getDistanceBetweenCuts <- function(bam, barcodes, subsetRegions = NA, nCores = 2,
                                   linePlot = "fragSizes_from_scCuts.png",
                                   heatmap = "fragSizes_from_scCuts_heatmap.png") {

  read_dists <- parallel::mclapply(barcodes, getFragLen, bam = bam, regions = subsetRegions, mc.cores = nCores)
  read_dists_df <- binFragLen(read_dists, barcodes)
  if(!is.na(linePlot)) {
    p1 <- ggplot(read_dists_df$i1, aes(group_tags, n, col = barcodes, group = barcodes)) + geom_line() + theme(legend.position = "none")
    ggsave(filename = linePlot, plot = p1)
  }

  read_dists_df %<>% as.data.frame()
  ## heatmap
  hmdf <- tidyr::spread(read_dists_df, group_tags, n) %>% as.data.frame() %>% dplyr::select(-barcodes) %>% as.matrix()
  hmdf[is.na(hmdf)] <- 0
  colnames(hmdf) %<>% gsub("\\[|\\|\\])", "", .) %>% gsub(",", "-", .) %>% paste0("bin_", .)

  if(!is.na(heatmap)) {
    pheatmap::pheatmap(hmdf[rev(order(rowMeans(hmdf))),], cluster_cols = FALSE,
                       cluster_rows = FALSE, scale = "row", color = rainbow(50), filename = heatmap)
  }
  return(read_dists_df)

}
