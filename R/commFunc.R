#' trapz: trapezoidal rule to approximate the integral values
#'
#' Returns approximation of integral.
#'
#' @param x A vector with \code{n} elements, \code{x[i]} is a support, \code{i = 1, ..., n}.
#' If \code{y} is \code{NULL}, support is taken as \code{seq(1, length(x), by = 1)}.
#' @param y \code{y[i, j]} is jth values on corresponding value of \code{x[i], i = 1, ..., n}.
#' If \code{y} is vector, the length of \code{y} must be equal to the lenght of \code{x}.
#' If \code{y} is matrix, the number of rows must be equal to the lenght of \code{x}.
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
#' @export
trapz <- function(x, y = NULL){
  if (is.null(y)) {
    y <- x
    x <- switch(is.matrix(x) + 1, {1:length(x)}, {1:nrow(x)})
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

#' Find the candidates of bandwidths for locPoly1d and locLinear2d
#'
#' The difference between \code{bwCandChooser2} and \code{bwCandChooser3} is whether the
#' candidates of bandwidths are the same on the x-axis and y-axis.
#' In our application, \code{bwCandChooser2} is used in finding the candidates of bandwidth of covariance
#' surface and \code{bwCandChooser3} is used in finding the candidates of bandwidth of cross-covariance surface.
#'
#' @param data An data.frame or data.table containing the variables in model.
#' @param id.var A string. The variable name of subject ID.
#' @param timeVarName A string. The variable name of time points.
#' @param sparsity The sparsity of data which is tested by \code{\link{checkSparsity}}.
#' @param kernel A string. It could be 'gauss', 'gaussvar', 'epan' or 'quar'.
#' @param degree An integer, the degree of polynomial.
#' @return The candidates of bandwidths
#' @examples
#' ## examples for bwCandChooser
#' data("regularExData", package = 'rfda')
#' bwCandChooser(regularExData, "sampleID", "t", 2, "gauss", 1)
#' @rdname bwCandChooser
#' @export
bwCandChooser <- function(data, id.var, timeVarName, sparsity, kernel, degree = 1){
  # check data
  assert_that(is.data.frame(data), is.character(id.var), is.character(timeVarName),
              !is.na(sparsity), is.finite(sparsity), sparsity %in% c(0, 1, 2),
              kernel %in% c('gauss','epan','gaussvar','quar'),
              !is.na(degree), is.finite(degree), degree > 0, degree - floor(degree) < 1e-6)

  # get the range of time points
  r <- diff(range(data[[timeVarName]]))
  # get the minimum bandwidth given sparsity of data
  if (sparsity == 0) {
    dstar <- find_max_diff_f(data[[timeVarName]], degree + 2)
    minBW <- ifelse(!is.na(dstar), ifelse(dstar > r/4, 0.75, 1) * 2.5 * dstar, NA)
  } else if (sparsity == 1) {
    minBW <- 2.0 * find_max_diff_f(data[[timeVarName]], degree + 1)
  } else if (sparsity == 2) {
    minBW <- 1.5 * find_max_diff_f(data[[timeVarName]], degree + 1)
  }

  # use range / 2 if kernel is gaussian and minimum is not found
  if ((is.na(minBW) || minBW < 0) && kernel == "gauss") {
    message("Data is too sparse, use the range / 2 as minimum bandwidth.")
    minBW <- 0.5 * r
  } else if ((is.na(minBW) || minBW < 0) && kernel != "gauss") {
    stop("Data is too sparse, no suitable bandwidth can be found! Try Gaussian kernel instead!\n");
  }

  # find the candidates
  minBW <- min(minBW, r)
  q <- (r / minBW / 4)^(1/9)
  return(sort(q^(0:9) * minBW, decreasing = FALSE))
}

#' @param dataAllGrid An data.table containing the grid of timepoints with
#'   naming \code{t1} and \code{t2} in model. (see examples.)
#' @examples
#'
#' ## examples for bwCandChooser2
#' require(data.table)
#' require(pipeR)
#'
#' data("sparseExData", package = 'rfda')
#' sparsity <- checkSparsity(sparseExData, "sampleID", "t")
#' sparseExData %>>% data.table %>>% `[`( , .(t1 = rep(t, length(t)),
#'   t2 = rep(t, each = length(t))), by = .(sampleID)) %>>%
#'   bwCandChooser2(sparsity, "gauss", 1)
#'
#' data("regularExData", package = 'rfda')
#' sparsity <- checkSparsity(regularExData, "sampleID", "t")
#' regularExData %>>% data.table %>>% `[`( , .(t1 = rep(t, length(t)),
#'   t2 = rep(t, each = length(t))), by = .(sampleID)) %>>%
#'   bwCandChooser2(sparsity, "gauss", 1)
#' @rdname bwCandChooser
#' @importFrom data.table data.table setorder is.data.table
#' @export
bwCandChooser2 <- function(dataAllGrid, sparsity, kernel, degree = 1){
  # cehck data
  assert_that(is.data.table(dataAllGrid), !is.na(sparsity), is.finite(sparsity), sparsity %in% c(0, 1, 2),
              kernel %in% c("gauss", "epan", "gaussvar", "quar"), !is.na(degree), is.finite(degree),
              degree > 0, degree - floor(degree) < 1e-6)

  # get output grid
  xout <- unique(dataAllGrid$t1)
  # get range of time points
  r <- diff(range(dataAllGrid$t1))
  # get the minimum bandwidth given sparsity of data
  if (sparsity == 0) {
    outGrid <- data.table(expand.grid(t1 = range(xout), t2 = xout))
    b <- dataAllGrid[t1 != t2, .(t1, t2)] %>>% rbind(outGrid) %>>% unique %>>% setorder(t2, t1) %>>% `$`(t1)
    minBW <- max(find_max_diff_f(xout, degree + 2), max(diff(b)) / 2.0)
  } else if (sparsity == 1) {
    minBW <- 2.0 * find_max_diff_f(xout, degree + 1)
  } else if (sparsity == 2) {
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

#' @examples
#'
#' ## examples for bwCandChooser3
#' # These examples are demo cases, we does not use this to find the candidates of
#' # bandwidths for smoothing covariance surface.
#' require(data.table)
#' require(pipeR)
#'
#' data("sparseExData", package = 'rfda')
#' sparsity <- checkSparsity(sparseExData, "sampleID", "t")
#' sparseExData %>>% data.table %>>% `[`( , .(t1 = rep(t, length(t)),
#'   t2 = rep(t, each = length(t))), by = .(sampleID)) %>>%
#'   bwCandChooser3(sparsity, "gauss", 1)
#'
#' data("regularExData", package = 'rfda')
#' sparsity <- checkSparsity(regularExData, "sampleID", "t")
#' regularExData %>>% data.table %>>% `[`( , .(t1 = rep(t, length(t)),
#'   t2 = rep(t, each = length(t))), by = .(sampleID)) %>>%
#'   bwCandChooser3(sparsity, "gauss", 1)
#' @rdname bwCandChooser
#' @importFrom data.table data.table setorder is.data.table
#' @export
bwCandChooser3 <- function(dataAllGrid, sparsity, kernel, degree = 1){
  # cehck data
  assert_that(is.data.table(dataAllGrid), !is.na(sparsity), is.finite(sparsity), sparsity %in% c(0, 1, 2),
              kernel %in% c("gauss", "epan", "gaussvar", "quar"), !is.na(degree), is.finite(degree),
              degree > 0, degree - floor(degree) < 1e-6)

  # get output grid
  xout <- unique(dataAllGrid$t1)
  # get range of time points
  r <- diff(range(dataAllGrid$t1))
  # get the minimum bandwidth given sparsity of data
  if (sparsity == 0) {
    outGrid <- data.table(expand.grid(t1 = range(xout), t2 = xout))
    b <- dataAllGrid[t1 != t2, .(t1, t2)] %>>% rbind(outGrid) %>>% unique %>>% setorder(t2, t1) %>>% `$`(t1)
    minBW <- max(find_max_diff_f(xout, degree + 2), max(diff(b)) / 2.0)
  } else if (sparsity == 1) {
    minBW <- 2.0 * find_max_diff_f(xout, degree + 1)
  } else if (sparsity == 2) {
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
  q <- (r / minBW / 4)^(1/4)
  bwMat <- (q^(0:4) * minBW) %>>% expand.grid(.) %>>% as.matrix %>>% `colnames<-`(NULL)
  return(bwMat[order(bwMat[,1], bwMat[,2], decreasing = FALSE), ])
}

#' Adjustment of optimal bandwidth choosed by gcv
#'
#' The usage of this function can be found in the examples of \code{\link{gcvLocPoly1d}} and
#' \code{\link{gcvLocLinear2d}}.
#'
#' @param bwOpt A numeric. The optimal bandwidth choosed by gcv.
#' @param sparsity The sparsity of data which is tested by \code{\link{checkSparsity}}.
#' @param kernel A string. It could be 'gauss', 'gaussvar', 'epan' or 'quar'.
#' @param drv An integer, the order of derivative.
#' @return An adjusted bandwidth.
#' @export
adjGcvBw <- function(bwOpt, sparsity, kernel, drv = 0){
  if (kernel == "gauss") {
    bwAdjFac <- switch(as.integer(sparsity == 2) + 1, c(1.1, 1.2, 2), c(1.1, 0.8, 0.8))
  } else if (kernel == "epan") {
    bwAdjFac <- switch(as.integer(sparsity == 2) + 1, c(1.1, 1.2, 1.5), c(1.1, 1.0, 1.1))
  }
  facTake <- ifelse(drv > 2, 2L, ifelse(drv >= 0, as.integer(drv) + 1, 0L))
  return(bwOpt * bwAdjFac[facTake])
}

#' Find the candidates of bandwidths for locLinear2d
#'
#' @param data A data.frame containing sample id, observed time points and correponding observed values.
#' @param subid The column name of the id of subjects.
#' @param sparsityRate A numeric vector between \code{0} and \code{1}. The proportion of data will be extracted.
#'   The length of sparsity must \code{1} or the number of curves (\code{n} in \code{\link{get_FPCA_opts}}).
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
#' @importFrom data.table between
#' @export
sparsify <- function(data, subid, sparsityRate){
  # cehck data
  assert_that(is.data.frame(data), subid %in% names(data), all(between(sparsityRate, 0, 1, FALSE)))

  # convert data to data.table with copy (not change the data)
  data <- data.table(data)
  #  get the unique subject id
  uniSubIds <- unique(data[[subid]])
  # check the length of sparsityRate
  if (length(sparsityRate) != length(uniSubIds) && length(sparsityRate) != 1)
    stop("The length of sparsityRate must 1 or the number of observation.")
  if (length(sparsityRate) == 1)
    sparsityRate <- rep(sparsityRate, length(uniSubIds))
  # sparsify data
  sparseDT <- mapply(function(dt, p) dt[sample(nrow(dt), round(nrow(dt)*p))],
                     split(data, data[[subid]]), 1 - sparsityRate, SIMPLIFY = FALSE) %>>% rbindlist
  return(sparseDT)
}

#' Find the candidates of bandwidths for locLinear2d
#'
#' @param DT A data.table containing list or vector in the cell.
#'   The cells in each row must have the same number of elements.
#' @param unnestCols The column names to unnest.
#' @return A unnested data.table.
#' @examples
#' require(data.table)
#' DT <- unnest(data.table(V1 = list(c(1,3,5), c(1,7)), V2 = list(c(2,5,3), c(4,6)), V3 = 1:2))
#'
#' require(jsonlite)
#' jsonDataFile <- system.file("extdata", "funcdata.json", package = "rfda")
#' # Following line may have a parse error with message "premature EOF has occured".
#' \dontrun{
#'   DT <- unnest(data.table(fromJSON(jsonDataFile)))
#' }
#' @importFrom data.table .SD
#' @importFrom plyr laply
#' @export
unnest <- function(DT, unnestCols = NULL){
  # check the columns to unnest
  if (is.null(unnestCols)) {
    unnestCols <- names(DT)[laply(DT, function(x) any(class(x) %in% "list"))]
    message("Automatically recognize the nested columns: ", paste0(unnestCols, collapse = ", "))
  }
  # check unnestCols is in the DT
  if (any(!unnestCols %in% names(DT)))
    stop(sprintf("The columns, %s, does not in the DT.",
                 paste0(unnestCols[!unnestCols %in% names(DT)], collapse = ", ")))
  # get the group by variable
  groupbyVar <- setdiff(names(DT), unnestCols)
  # generate the expression to remove group by variable
  chkExpr <- paste0(groupbyVar, "=NULL", collapse = ",") %>>% (paste0("`:=`(", ., ")"))
  # check the lengths of each cell in list-column are all the same
  chkLenAllEqual <- DT[ , lapply(.SD, function(x) laply(x, length)), by = groupbyVar] %>>%
    `[`(j = eval(parse(text = chkExpr))) %>>% as.matrix %>>% apply(1, diff) %>>% `==`(0) %>>% all
  if (!chkLenAllEqual)
    stop("The length in each cell is not equal.")

  # generate unnest expression
  expr <- unnestCols %>>% (paste0(., "=unlist(",  ., ")")) %>>%
    paste0(collapse = ",") %>>% (paste0(".(", ., ")"))
  # return unnested data.table
  return(DT[ , eval(parse(text = expr)), by = groupbyVar])
}
