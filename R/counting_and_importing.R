

#' Get total counts per barcode from BAM files
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



#' get counts in a given region (GRanges) per cell from BAM files
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



#' Import a loom file as SingleCellExperiment using loomR API
#'
#' @description This function is only needed in case the loom file contains
#'   "obsm" and "varm" entries, which causes loomExperiment to error. Otherwise
#'   simply use LoomExperiment::import(). In this case the "obsm" and "varm"
#'   entries are added to rowData and colData.
#'
#' @param loomfile Path to loom file
#' @param obs_names Attribute that stores observation IDs (eg. cell IDs) in loom file
#' @param var_names Attribute that stores variable IDs (eg. gene IDs) in loom file
#' @param reducedDimNames PCA/UMAP attribute names (from adata.obsm.keys())
#'
#' @return SingleCellExperiment object
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom magrittr %>%
#' @importFrom S4Vectors DataFrame
#' @export
#'
#' @examples
#'
#'
loom_to_sce <- function(loomfile,
                        obs_names="obs_names",
                        var_names="var_names",
                        reducedDimNames=c("X_pca", "X_umap")) {

  adata = loomR::connect(loomfile, skip.validate = T)

  message("Storing reducedDims")
  rDIMs <- lapply(reducedDimNames, function(n){
    d <- adata$col.attrs[[n]]$dims
    if(length(d) > 1) {
      out <- t(adata$col.attrs[[n]][,])
      colnames(out) <- paste0(toupper(n), 1:ncol(out))
      return(out)
    }
  })
  names(rDIMs) <- reducedDimNames
  rDIMs <- rDIMs[!sapply(rDIMs, is.null)]

  message("Storing all cell info to col_attrs")
  sapply(names(adata$col.attrs), function(n){
    d <- adata$col.attrs[[n]]$dims
    if(length(d) == 1) {
      return(adata$col.attrs[[n]][])
    } else if(!(n %in% reducedDimNames)) {
      out <- t(adata$col.attrs[[n]][,])
      colnames(out) <- paste0(toupper(n), 1:ncol(out))
      return(out)
    }
  }, simplify = FALSE) %>% do.call(cbind, .) %>% DataFrame() -> col_attrs
  rownames(col_attrs) <- col_attrs[[obs_names]]

  message("Storing all var info as row_attrs (also includes adata.varm)")
  sapply(names(adata$row.attrs), function(n){
    d <- adata$row.attrs[[n]]$dims
    if(length(d) == 1) {
      return(adata$row.attrs[[n]][])
    } else {
      out <- t(adata$row.attrs[[n]][,])
      colnames(out) <- paste0(toupper(n), 1:ncol(out))
      return(out)
    }
  }, simplify = FALSE) %>% do.call(cbind, .) %>% DataFrame() -> row_attrs
  rownames(row_attrs) <- row_attrs[[var_names]]

  message("Storing all layers as sparse matrices")
  sapply(names(adata[['layers']]), function(layer){
    print(layer)
    name = paste0("layers/", layer)
    mat <- adata[[name]][,]
    out <- Matrix::t(Matrix::Matrix(mat, sparse = T))
    rownames(out) <- rownames(row_attrs)
    colnames(out) <- rownames(col_attrs)
    return(out)
  }) -> layers

  message("Converting to SingleCellExperiment")
  sce <- SingleCellExperiment(layers,
                              rowData = row_attrs,
                              colData = col_attrs,
                              reducedDims = rDIMs)

  return(sce)

}
