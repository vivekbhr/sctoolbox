#' filter Reads based on nucleotide content at 5-prime
#'
#' @param bamdf bam file as data frame
#' @param pos_strand_expected
#' @param neg_strand_expected
#' @param plotFile
#' @param bsgenome
#'
#' @return
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
#' @return
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
