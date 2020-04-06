#' Regions detected as a function of counts
#'
#' @param matrix A region * cell matrix of counts
#'
#' @return ggplot object
#' @export
#'
#' @examples
#'
#'
regionDetectionPlot <- function(mat) {
  tc <- colSums(mat)
  tcy <- colSums(mat > 0)
  plot(tc, tcy, log = "xy", xlab = "log(total Counts)", ylab = "log(No. of regions with non-zero counts)")
}
