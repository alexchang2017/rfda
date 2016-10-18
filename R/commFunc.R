#' trapz: trapezoidal rule to approximate the integral values
#'
#' Returns approximation of integral.
#'
#' @param x A vector with n elements, \code{x[i]} is a support, \code{i = 1, ..., n}.
#' If y is NULL, support is taken as \code{seq(1, length(x), by = 1)}.
#' @param y \code{y[i, j]} is jth values on corresponding value of \code{x[i], i = 1, ..., n}.
#' If y is vector, the length of y must be equal to the lenght of x.
#' If y is matrix, the number of rows must be equal to the lenght of x.
#' @return A value, the approximation of integral.
#' @section Reference:
#' Kai Habel, trapz, Octave.
#' @examples
#' # case 1
#' x <- c(1, 4, 9, 16, 25)
#' trapz(x) # 42
#'
#' # case 2
#' x <- matrix(c(1,4,9,16,25,1,8,27,64,125), 5)
#' trapz(x) # 42 162
#' @rdname trapz
#' @export
trapz <- function(x, y = NULL){
  if (is.null(y))
  {
    y <- x
    if (is.vector(x))
      x <- 1:length(x)
    else
      x <- 1:nrow(x)
  }
  return(as.vector(trapz_cpp(x, as.matrix(y))))
}

# function to perform data binning
#' @importFrom data.table setnames
#' @importFrom utils head
binData <- function(data, numBins){
  # check data
  assert_that(!is.na(numBins), is.finite(numBins), numBins > 0, numBins - floor(numBins) < 1e-6)

  # find the boudaries to split data
  boundaries <- seq(min(data$timePnt), max(data$timePnt), length.out = numBins + 1)
  # find the middle points to stand time points of binned data
  newTimePnts <- head(boundaries, numBins) + diff(boundaries) / 2
  # average the data in the interval for data binning
  newDataDT <- data %>>% `[`(j = idx_agg := findInterval(timePnt, boundaries, TRUE), by = .(subId,variable)) %>>%
    `[`(j = .(value = mean(value), timePnt = newTimePnts[idx_agg]), by = .(subId,variable,idx_agg)) %>>%
    `[`(j = idx_agg := NULL)
  return(newDataDT)
}

# sub-function for bwCandChooser
#' @importFrom utils head tail
find_max_diff_f <- function(t, lag_n){
  assert_that(!is.na(lag_n), is.finite(lag_n), lag_n > 0, lag_n - floor(lag_n) < 1e-6)

  sort_t <- sort(t)
  n <- length(t)
  if (n < lag_n)
    return(NA)
  if (lag_n > 1)
    return(max(tail(sort_t, n - lag_n + 1) - head(sort_t, n - lag_n + 1)))
  else
    return(max(diff(sort_t))/2)
}

#' Find the candidates of bandwidths for locPoly1d
#'
#' @param data An data.frame or data.table containing the variables in model.
#' @param id.var A string. The variable name of subject ID.
#' @param timeVarName A string. The variable name of time points.
#' @param sparsity The sparsity of data which is tested by \code{\link{checkSparsity}}.
#' @param kernel The kernel function used in locPoly1d.
#' @param degree The degree of polyminal function used in locPoly1d.
#' @return The candidates of bandwidths
#' @examples
#' data("regularExData", package = 'rfda')
#' bwCandChooser(regularExData, "sampleID", "t", 2, "gauss", 1)
#' @rdname bwCandChooser
#' @export
bwCandChooser <- function(data, id.var, timeVarName, sparsity, kernel, degree){
  # check data
  assert_that(is.data.frame(data), is.character(id.var), is.character(timeVarName),
              !is.na(sparsity), is.finite(sparsity), sparsity %in% c(0, 1, 2),
              kernel %in% c('gauss','epan','gaussvar','quar'),
              !is.na(degree), is.finite(degree), degree > 0, degree - floor(degree) < 1e-6)

  # get the range of time points
  r <- diff(range(data[[timeVarName]]))
  # get the minimum bandwidth given sparsity of data
  if (sparsity == 0){
    dstar <- find_max_diff_f(data[[timeVarName]], degree + 2)
    if (!is.na(dstar)){
      if (dstar > r / 4){
        dstar <- dstar * 0.75
        message(sprintf("The min bandwidth choice is too big, reduce to %.6f", dstar))
      }
      minBW <- 2.5 * dstar
    } else {
      minBW <- NA
    }
  } else if (sparsity == 1){
    minBW <- 2.0 * find_max_diff_f(data[[timeVarName]], degree + 1)
  } else if (sparsity == 2){
    minBW <- 1.5 * find_max_diff_f(data[[timeVarName]], degree + 1)
  }

  # use range / 2 if kernel is gaussian and minimum is not found
  if ((is.na(minBW) || minBW < 0) && kernel == "gauss"){
    message("Data is too sparse, use the range / 2 as minimum bandwidth.")
    minBW <- 0.5 * r;
  } else if ((is.na(minBW) || minBW < 0) && kernel != "gauss"){
    stop("Data is too sparse, no suitable bandwidth can be found! Try Gaussian kernel instead!\n");
  }

  # find the candidates
  minBW <- min(minBW, r)
  q <- (r / minBW / 4)^(1/9)
  return(sort(q^(0:9) * minBW, decreasing = FALSE))
}

#' Find the candidates of bandwidths for locLinear2d
#'
#' @param dataAllGrid An data.table containing the id of subjects nameing \code{sampleID} and
#'   all timepoints combinations by variables naming \code{t1} and \code{t2} in model. (see examples.)
#' @return The candidates of bandwidths
#' @examples
#' require(data.table)
#' require(pipeR)
#'
#' data("sparseExData", package = 'rfda')
#' sparsity <- checkSparsity(sparseExData, "sampleID", "t")
#' sparseExData %>>% data.table %>>% `[`( , .(t1 = rep(t, length(t)),
#'   t2 = rep(t, each=length(t))), by = .(sampleID)) %>>%
#'   bwCandChooser2("sampleID", c("t1", "t2"), sparsity, "gauss", 1)
#'
#' data("regularExData", package = 'rfda')
#' sparsity <- checkSparsity(regularExData, "sampleID", "t")
#' regularExData %>>% data.table %>>% `[`( , .(t1 = rep(t, length(t)),
#'   t2 = rep(t, each=length(t))), by = .(sampleID)) %>>%
#'   bwCandChooser2("sampleID", c("t1", "t2"), sparsity, "gauss", 1)
#' @rdname bwCandChooser
#' @importFrom data.table data.table setorder is.data.table
#' @export
bwCandChooser2 <- function(dataAllGrid, id.var, timeVarName, sparsity, kernel, degree){
  # cehck data
  assert_that(is.data.table(dataAllGrid), is.character(id.var), all(is.character(timeVarName)),
              length(timeVarName) == 2L, !is.na(sparsity), is.finite(sparsity), sparsity %in% c(0, 1, 2),
              kernel %in% c('gauss','epan','gaussvar','quar'),
              !is.na(degree), is.finite(degree), degree > 0, degree - floor(degree) < 1e-6)

  # get output grid
  xout <- unique(dataAllGrid$t1)
  # get range of time points
  r <- diff(range(dataAllGrid$t1))
  # get the minimum bandwidth given sparsity of data
  if (sparsity == 0){
    outGrid <- data.table(expand.grid(t1 = range(xout), t2 = xout))
    b <- dataAllGrid[t1 != t2, .(t1, t2)] %>>% rbind(outGrid) %>>% unique %>>% setorder(t2, t1) %>>% `$`(t1)
    minBW <- max(find_max_diff_f(xout, degree + 2), max(diff(b)) / 2.0)
  } else if (sparsity == 1){
    minBW <- 2.0 * find_max_diff_f(xout, degree + 1)
  } else if (sparsity == 2){
    minBW <- 1.5 * find_max_diff_f(xout, degree + 2)
  }

  # shrink the minimum bandwidth if kernel is gaussian
  if (kernel == "gauss") {
    if (is.na(minBW) || minBW < 0){
      message("Data is too sparse, use the max(t) as minimum bandwidth.")
      minBW <- max(xout)
    }
    minBW <- 0.2 * minBW;
  } else if ((is.na(minBW) || minBW < 0) && kernel != "gauss") {
    stop("Data is too sparse, no suitable bandwidth can be found! Try Gaussian kernel instead!\n");
  }

  # find the candidates
  minBW <- min(minBW, r / 4)
  q <- (r / minBW / 4)^(1/9)
  bwMat <- matrix(rep(q^(0:9) * minBW, 2), 10)
  return(bwMat[order(bwMat[,1], bwMat[,2], decreasing = FALSE), ])
}

#' Find the candidates of bandwidths for locLinear2d
#'
#' @param data A data.frame containing sample id, observed time points and correponding observed values.
#' @param subid The column name of the id of subjects.
#' @param sparsity A numeric vector between 0 and 1. The proportion of data will be extracted.
#'   The length of sparsity must 1 or the number of observation.
#'   (The number of observation does not mean the number of rows of data, see examples.)
#' @return A data.frame after sparsifying.
#' @examples
#' require(ggplot2)
#' tp <- seq(1, 10, len = 101)
#' DT <- funcDataGen(3, tp, function(x) sin(x), function(x) rep(1, length(x)), "BesselJ")
#' sparseDT <- sparsify(DT, "sampleID", 0.85)
#' ggplot(sparseDT, aes(x = t, y = y, color = factor(sampleID))) +
#'   geom_line() + geom_point() + labs(color = "Sample ID")
#'
#' message("The number of observation is ", no <- length(unique(DT$sampleID)))
#' sparseDT2 <- sparsify(DT, "sampleID", runif(no))
#' ggplot(sparseDT2, aes(x = t, y = y, color = factor(sampleID))) +
#'   geom_line() + geom_point() + labs(color = "Sample ID")
#' @importFrom data.table between setDT
#' @export
sparsify <- function(data, subid, sparsity){
  # cehck data
  assert_that(is.data.frame(data), subid %in% names(data), all(between(sparsity, 0, 1, FALSE)))

  setDT(data)
  uniSubIds <- unique(data[[subid]])
  if (length(sparsity) != length(uniSubIds) && length(sparsity) != 1)
    stop("The length of sparsity must 1 or the number of observation.")
  if (length(sparsity) == 1)
    sparsity <- rep(sparsity, length(uniSubIds))
  if (length(uniSubIds))

  sparseDT <- mapply(function(dt, p) dt[sample(nrow(dt), round(nrow(dt)*p))],
                     split(data, data[[subid]]), 1-sparsity, SIMPLIFY = FALSE) %>>% rbindlist
  return(sparseDT)
}
