

#' Get total counts per cell
#'
#' @param barcodelist vector with barcodes
#' @param bamFiles bamfiles (path)
#'
#' @return data frame with counts per cell
#' @export
#'
#' @examples
#'
#'
totalCounts_perCell <- function(barcodelist, bamFiles) {
  tpc <- BiocParallel::bplapply(barcodeList, function(bc){
    o <- countBam(bamFiles, ignore.strand = TRUE,
                  param = ScanBamParam(tag = "BC", tagFilter = list("BC" = bc)),
                  BPPARAM = MulticoreParam(workers = 6))
    return(o)
  }, BPPARAM = BiocParallel::MulticoreParam(workers = 5))

  tpc <- plyr::ldply(totalCounts_perCell, function(x) return(x$records))
  tpc$barcodes <- barcodeList
  return(tpc)
}



#' get counts in a given GRanges per cell
#'
#' @param regions GRanges of regions to count
#' @param barcodelist vector, barcode list
#' @param bamFiles bam files path
#'
#' @return counts per region per cell
#' @export
#'
#' @examples
#'
countsInRegions_perCell <- function(regions, barcodelist, bamFiles) {
    BiocParallel::bplapply(barcodelist, function(bc){
      print(bc)
      o <- GenomicAlignments::summarizeOverlaps(regions, bamFiles, ignore.strand = TRUE,
                             param = Rsamtools::ScanBamParam(tag = "BC", tagFilter = list("BC" = bc)),
                             mode = "Union",
                             BPPARAM = MulticoreParam(workers = 5))
      return(o)
    }, BPPARAM = BiocParallel::MulticoreParam(workers = 6)) -> counts

  return(counts)
}

