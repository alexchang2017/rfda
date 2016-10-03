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
  return(trapz_cpp(x, as.matrix(y)))
}

#' @importFrom data.table setnames
binData <- function(data, numBins){
  if (numBins <= 0 || abs(numBins - floor(numBins)) > 0 || is.na(numBins) || is.infinite(numBins))
    stop("numBins is not positive integer.")
  boundaries <- seq(min(data$timePnt), max(data$timePnt), length.out = numBins + 1)
  newTimePnts <- head(boundaries, numBins) + diff(boundaries) / 2
  newDataDT <- data %>>% `[`( , idx_agg := findInterval(timePnt, boundaries, rightmost.closed = TRUE),
                              by = "subId,variable") %>>%
    `[`( , .(new_value = mean(value), new_timePnt = newTimePnts[idx_agg]), by = "subId,variable,idx_agg") %>>%
    setnames(c("new_value", "new_timePnt"), c("value", "timePnt")) %>>% `[`( , idx_agg := NULL)
  return(newDataDT)
}

# sub-function for bwCandChooser
find_max_diff_f <- function(t, lag_n){
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
  r <- diff(range(data[[timeVarName]]))

  if (sparsity == 0){
    dstar <- find_max_diff_f(data[[timeVarName]], degree + 2)
    if (dstar > r / 4){
      dstar <- dstar * 0.75
      message(sprintf("The min bandwidth choice is too big, reduce to %.6f", minBW))
    }
    minBW <- 2.5 * dstar
  } else if (sparsity == 1){
    minBW <- 2.0 * find_max_diff_f(data[[timeVarName]], degree + 1)
  } else if (sparsity == 2){
    minBW <- 1.5 * find_max_diff_f(data[[timeVarName]], degree + 1)
  } else {
    stop("sparsity must be 0, 1 or 2.\n");
  }

  if (is.na(minBW) && kernel == "gauss"){
    minBW <- 0.5 * r;
  } else if (minBW < 0 && kernel != "gauss"){
    stop("The data is too sparse, no suitable bandwidth can be found! Try Gaussian kernel instead!\n");
  }

  minBW <- min(minBW, r)
  q <- (r / minBW / 4)^(1/9)
  return(q^(0:9) * minBW)
}
