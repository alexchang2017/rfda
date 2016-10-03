#' Test whether data is regular
#'
#' If the time points of each subject are equal, then return 2.
#' If the ratio of the number of pooled time points and
#' the number of unique pooled time points multiply the number of subjects is
#' less than 75%, we say that data is regular with missing values and return 1.
#' Otherwise, data is sparse or irregular and return 0.
#'
#' @param data An data.frame or data.table containing the variables in model.
#' @param id.var A string. The variable name of subject ID.
#' @param timeVarName A string. The variable name of time points.
#' @return An integer. 0 for sparse and irregular data. 1 for regular data with missing values.
#'   2 for completely balanced data.
#' @examples
#' # sparse case
#' data("sparseExData", package = 'rfda')
#' checkSparsity(sparseExData, "sampleID", "t") # [1] 0
#' # regular case with missing
#' data("irregularExData", package = 'rfda')
#' checkSparsity(irregularExData, "sampleID", "t") # [1] 1
#' # regular case
#' data("regularExData", package = 'rfda')
#' checkSparsity(regularExData, "sampleID", "t") # [1] 2
#' @export
checkSparsity <- function(data, id.var, timeVarName){
  propNonSparse <- nrow(data) / length(unique(data[[timeVarName]])) / length(unique(data[[id.var]]))
  return(ifelse(propNonSparse == 1, 2, ifelse(propNonSparse > 0.75, 1, 0)))
}

#' Calculating the functional principal components
#'
#' Return mean, covariance, eigen functions and fpc scores.
#'
#' @param formula An object of class "formula", a description of FPCA model.
#' The RHS of `~` is the responses and LHS is time variables.
#' Notice that the names of variables must be a-z, A-Z or _.
#' @param id.var An character to indicate the subject id of data.
#' @param data An data.frame or data.table containing the variables in model.
#' @param options An list containing the options to fit FPCA model.
#' @return An list containing mean, covariance, eigen functions,
#'   functional principal components scores, etc. Please see "Details" for more details.
#' @seealso \code{\link{get_FPCA_opts}} for input options.
#' @examples
#' \dontrun{
#' # sparse case
#' data("sparseExData", package = 'rfda')
#'
#' # regular with missing values case
#' data("irregularExData", package = 'rfda')
#'
#' # regular case
#' data("regularExData", package = 'rfda')
#' }
#'
#' @importFrom plyr is.formula
#' @importFrom RcppParallel setThreadOptions
#' @importFrom data.table data.table melt.data.table setnames
#' @importFrom stats median
#' @importFrom utils modifyList
#' @export
FPCA <- function(formula, id.var, data, options = list()){
  # formula = as.formula("y ~ t")
  # formula = as.formula("y + y2 ~ t")
  # id.var = "sampleID"
  # data("irregularExData", package = "rfda")
  # data <- irregularExData %>>% data.table %>>% `[`( , y2 := y*2 + rnorm(nrow(.)))
  # options <- list()
  assert_that(is.formula(formula), is.character(id.var),
              length(id.var) == 1, is.data.frame(data))

  # check the formula
  chkFmLHS <- as.character(formula[[2]]) %>>% (grepl("[+a-zA-z_]", .)) %>>% all
  chkFmRHS <- as.character(formula[[3]]) %>>% (grepl("[+a-zA-z_]", .)) %>>% length %>>% `==`(1)
  chkFormla <- chkFmLHS || chkFmRHS
  message("Checking the formula...")
  assert_that(chkFormla, msg = "Check failed")

  # set the number of thread be used
  if (FPCA_opts$ncpus != 0)
    setThreadOptions(FPCA_opts$ncpus)

  # find the names of variables and name of variable indicating time points
  varName <- setdiff(all.vars(formula), as.character(formula[[3]]))
  timeVarName <- as.character(formula[[3]])

  # get the full options of FPCA and check
  default_FPCA_opts <- get_FPCA_opts(length(varName))
  optNamesUsed <- names(options) %in% default_FPCA_opts
  FPCA_opts <- modifyList(default_FPCA_opts, options[optNamesUsed]) %>>%
    chk_FPCA_opts(nrow(data))
  if (any(!optNamesUsed))
    paste(names(options)[!optNamesUsed], collapse = ", ") %>>%
      sprintf(fmt = "Ignoring the non-found options %s.") %>>% message

  # get the sparsity of data
  sparsity <- checkSparsity(data, id.var, timeVarName)

  # melt table to get a data.table to get a simple data.table and remove the NA, NaN and Inf.
  # additionally, give names for id.var and timeVarName
  dataDT <- melt.data.table(data.table(data), id.vars = c(id.var, timeVarName),
                            measure.vars = varName, variable.factor = FALSE) %>>%
    `[`(!is.na(value) & is.finite(value)) %>>%
    setnames(c(id.var, timeVarName) , c("subId", "timePnt"))

  # find the number of observations for each observed function
  subIdInsuffSize <- dataDT[, .(numObv = .N) ,by = c("subId", "variable")][numObv <= 1] %>>%
    `$`(subId) %>>% unique
  # remove the
  dataDT <- dataDT[!subId %in% subIdInsuffSize]

  # binning data
  if (FPCA_opts$numBins == -1 || FPCA_opts$numBins >= 10)
  {
    if (FPCA_opts$numBins == -1){
      numObv <- nrow(dataDT[variable == dataDT$variable[1]])
      numObvEach <- dataDT[variable == dataDT$variable[1] , .(numObv = .N), by = "timePnt"] %>>%
        `$`(numObv)
      criterionNum <- ifelse(sparsity == 0, stats::median(numObvEach), max(numObvEach))
      if (criterionNum > 400) {
        FPCA_opts$numBin <- 400
      } else if (numObv > 5000 && criterionNum > 20) {
        FPCA_opts$numBin <- min(criterionNum, ceiling(max(20, (5000-numObv)*19/2250+400)))
      }
    }

    if (FPCA_opts$numBins > 0) {
      dataDT <- binData(dataDT, FPCA_opts$numBins)
      sparsity <- checkSparsity(dataDT[variable == dataDT$variable[1]], "subId", "timePnt")
    }
  } else if (FPCA_opts$numBins < 10 && FPCA_opts$numBins != 0)
  {
    warning('The number of bins must be at least 10! No binning will be performed!')
  }

  # get weight
  if (FPCA_opts$weighted)
  {
    byVars <- switch(as.character(sparsity), "0" = c("variable", "subId"),
                     "1" = c("variable", "timePnt"), "2" = c("variable", "timePnt"))
    dataDT <- merge(dataDT, dataDT[ , .(weight = 1/.N), by = byVars], by = byVars)
  } else
  {
    dataDT[ , weight := 1]
  }

  # Initializing allTimePnts is based on the unique time points of the pooled data + the unique
  # time points of "newdata", the output time grid.
  allTimePnts <- sort(unique(c(dataDT$timePnt, FPCA_opts$newdata)))
  # Initializing sampledTimePnts is based on the number of grid to be chosen in the range of
  # all time span.
  sampledTimePnts <- seq(min(allTimePnts), max(allTimePnts), length.out = FPCA_opts$numGrid)

  # get mean functions of variables
  dataList <- split(dataDT, dataDT$variable)

  validMFList <- length(FPCA_opts$userMeanFuncList) == length(dataList) &&
    all(sapply(FPCA_opts$userMeanFuncList, length) == length(sampledTimePnts)) &&
    all(sapply(FPCA_opts$userMeanFuncList, function(mf) all(!is.na(mf)) && all(!is.infinite(mf))))
  if (validMFList) {
    MFList <- FPCA_opts$userMeanFuncList
    MDFList <- lapply(FPCA_opts$userMeanFuncList, function(mf){
      as.vector(interp1(sampledTimePnts, as.matrix(mf), allTimePnts, 'spline'))
    })
    bwOptLocPoly1d <- NA
  } else {
    MFRes <- lapply(dataList, function(dat){
      bwCand <- bwCandChooser(dat, "subId", "timePnt", sparsity, FPCA_opts$kernel, 1)
      bwOptLocPoly1d <- gcv_locPoly1d(bwCand, dat$timePnt, dat$value,
                                      dat$weight, FPCA_opts$kernel, 0, 1)
      bwOptLocPoly1d <- adjGcvBw1d(bwOptLocPoly1d, sparsity, FPCA_opts$kernel, 0)
      if (FPCA_opts$bwMean == -1)
        bwOptLocPoly1d <- sqrt(find_max_diff_f(dat[["timePnt"]], 2) * bwOptLocPoly1d)
      return(list(locPoly1d(bwOptLocPoly1d, dat$timePnt, dat$value, dat$weight,
                sampledTimePnts, FPCA_opts$kernel, 0, 1) %>>% as.vector,
                locPoly1d(bwOptLocPoly1d, dat$timePnt, dat$value, dat$weight,
                          allTimePnts, FPCA_opts$kernel, 0, 1) %>>% as.vector))
    })
    MFList <- lapply(MFRes, `[[`, 1)
    MDFList <- lapply(MFRes, `[[`, 2)
    rm(MFRes)
  }

  if (sapply(MFList, function(mf) all(is.na(mf))) %>>% any)
    stop(paste0("The bandwidth of mean function is not appropriately!\n",
                "If it is chosen automatically"))

  if (FPCA_opts$normMethod == "variance"){
    # get rawCov

    # get smoothed covariance

    # normalization

  } else if (FPCA_opts$normMethod == "IQR") {
    # get IQR

    # normalization

  }

  # find cross-covariance

  # decide the number of FPC

  # calculation of FPC scores

  return(1)
}
