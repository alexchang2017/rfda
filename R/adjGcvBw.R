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
