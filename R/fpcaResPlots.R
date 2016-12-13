#' The design plot of the functional data
#'
#' @param x An object used to select a method.
#'   Default is an data.frame or data.table containing the variables in model.
#' @param ... Other arguments passed into \code{\link{levelplot}}.
#' @return A \code{\link{levelplot}} object.
#' @rdname designPlot
#' @importFrom lattice levelplot
#' @importFrom grDevices heat.colors
#' @export designPlot
designPlot <- function(x, ...) UseMethod("designPlot")

#' @param id.var An character to indicate the subject id of data.
#' @param time.var An character to indicate the time points of data.
#' @param diagonal A logical value. Plot diagonal points or not.
#' @param col.regions Default is \code{heat.colors(25)}.
#'   The color vector to be used. Please refer to \code{\link{levelplot}}.
#' @examples
#' data("irregularExData", package = 'rfda')
#' designPlot(irregularExData, "sampleID", "t")
#' @rdname designPlot
#' @method designPlot default
#' @export
designPlot.default <- function(x, id.var, time.var, diagonal = FALSE, col.regions = heat.colors(25), ...) {
  setDT(x)
  setnames(x, time.var, "t")
  graphDT <- x[ , .(t1 = rep(t, length(t)), t2 = rep(t, each = length(t))), by = id.var] %>>%
    `[`(j = .(cnt = .N), by = .(t1, t2))
  if (!diagonal)
    graphDT <- graphDT[t1 != t2]
  return(levelplot(cnt ~ t1 + t2, graphDT, col.regions = col.regions, ...))
}

#' @examples
#' \dontrun{
#' data("irregularExData", package = 'rfda')
#' fpcaResIrregular <- FPCA(y ~ t, "sampleID", irregularExData)
#' designPlot(fpcaResIrregular)
#' }
#' @rdname designPlot
#' @method designPlot fpcaRes
#' @export
designPlot.fpcaRes <- function(x, ...){
  return(designPlot(x$data, x$id.var, x$time.var, ...))
}

#' screeplot for
#'
#' @param x An object used to select a method.
#'   Default is a numeric vector. Generally, it is the eigenvalues of covariance.
#' @param ... Other arguments passed into \code{\link{barchart}}.
#' @return A \code{\link{barchart}} object.
#' @rdname screePlot
#' @importFrom data.table data.table
#' @importFrom lattice barchart panel.barchart panel.lines panel.abline
#' @exportClass fpcaRes
#' @export screePlot
screePlot <- function(x, ...) UseMethod("screePlot")

#' @param n The number of components to take.
#' @param ylim The plot limit of y-axis.
#' @param yAxisLabels The location of y-axis Labels.
#' @param col.bar The colour of bar.
#' @param col.line The colour of line.
#' @param pch.line The type of points.
#' @param lty.line The line type of line.
#' @examples
#' set.seed(100)
#' randCov <- matrix(rnorm(100), 10, 10)
#' randCov <- (randCov + t(randCov)) / 2
#' diag(randCov) <- diag(randCov) + 2
#' eigRes <- eigen(randCov, TRUE)
#' idx <- eigRes$values > 0
#' screePlot(eigRes$values[idx], sum(idx))
#' screePlot(eigRes$values[idx], 4, ylim = c(0, 0.8),
#'           col.line = "red", pch.line = 17, lty.line = 4)
#' @rdname screePlot
#' @method screePlot default
#' @export
screePlot.default <-  function(x, n, ylim = c(0, 1.02), yAxisLabels = seq(0, 1, by = 0.2),
                       col.bar = "gray", col.line = "blue", pch.line = 16, lty.line = 2, ...) {
  graphDT <- data.table(comp = paste0("Comp. ", 1L:n), FVE = x[1:n] / sum(x)) %>>% `[`(j = cumFVE := cumsum(FVE))
  return(barchart(FVE ~ comp, graphDT, xlab = "", ylab = "Fraction of Variation Explained",
                  ylim = ylim, col = col.bar,
                  panel = function(...) {
                    panel.barchart(...)
                    panel.lines(1L:length(graphDT$comp), graphDT$cumFVE,
                                col = col.line, type = "b", lty = lty.line, pch = pch.line)
                    panel.abline(h = seq(0, 1, by = 0.1), lty = 2, col = "gray")
                  }, scales = list(y = list(at = yAxisLabels,
                                            labels = sprintf("%2.0f%%", 100 * yAxisLabels))),
                  key = list(corner = c(0.95, 0.55), lines = list(lty = lty.line, col = col.line),
                             text = list("Cum. FVE")), ...))
}


#' @examples
#' \dontrun{
#' data("irregularExData", package = 'rfda')
#' fpcaResIrregular <- FPCA(y ~ t, "sampleID", irregularExData)
#' screePlot(fpcaResIrregular)
#' }
#' @rdname screePlot
#' @method screePlot fpcaRes
#' @export
screePlot.fpcaRes <- function(x, ...){
  return(screePlot(x$FVE, x$numFPC, ...))
}
