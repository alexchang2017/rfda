#' Adjust the bandwidth selected by gcv_locPoly1d
#'
#' @param bw_opt The bandwidth valued selected by gcv_locPoly1d.
#' @param sparsity The sparsity of data. Check it in \code{\link{checkSparsity}}.
#' @param kernel A string. It could be 'gauss', 'gaussvar', 'epan' or 'quar'.
#' @param drv An integer, the order of derivative.
#' @return A adjusted bandwidth.
#' @examples
#' data("regularExData", package = 'rfda')
#' sparsity <- checkSparsity(regularExData, "sampleID", "t")
#' degree <- 1
#' regBwCand <- bwCandChooser(regularExData, "sampleID", "t", sparsity, "gauss", degree)
#' bw_gcv <- gcv_locPoly1d(regBwCand, regularExData$t, regularExData$y,
#'                         rep(1, nrow(regularExData)), "gauss", degree-1, degree)
#' adjGcvBw1d(bw_gcv, sparsity, "gauss", degree-1)
#' @export
adjGcvBw1d <- function(bw_opt, sparsity, kernel, drv){
  if (kernel == "gauss"){
    bwAdjFac <- ifelse(sparsity == 2, c(1.1, 0.8, 0.8), c(1.1, 1.2, 2))
  } else if (kernel == "epan"){
    bwAdjFac <- ifelse(sparsity == 2, c(1.1, 1.0, 1.1), c(1.1, 1.2, 1.5))
  }
  facTake <- ifelse(drv > 2, 2L, ifelse(drv >= 0, as.integer(drv)+1, 0L))
  return(bw_opt * bwAdjFac[facTake])
}

# adjGcvBw2d <- function(){
#   return(0.0)
# }

