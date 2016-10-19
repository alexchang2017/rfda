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

#' @importFrom data.table setorder rbindlist .N
getRawCrCov <- function(demeanDataDT){
  # geneerate the all combinations of t1,t2 and varaibles
  baseDT <- demeanDataDT[ , .(t1 = rep(timePnt, length(timePnt)), t2 = rep(timePnt, each=length(timePnt)),
                              value.var1 = rep(value, length(timePnt))), by = .(variable, subId)]
  # calculation of raw cross-covariance
  rawCrCovDT <- split(demeanDataDT, demeanDataDT$variable) %>>%
    lapply(function(dt) merge(baseDT[variable >= dt$variable[1]], dt, suffixes = c("1", "2"),
                              by.x = c("subId", "t2"), by.y = c("subId", "timePnt"))) %>>%
    rbindlist %>>% setnames("value", "value.var2") %>>%
    `[`(j = .(sse = sum(value.var1 * value.var2), cnt = .N), by = .(variable1, variable2, t1, t2)) %>>%
    setorder(variable1, variable2, t1, t2) %>>% `[`(j = weight := 1)
  return(rawCrCovDT)
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
#' @importFrom plyr is.formula
#' @importFrom RcppParallel setThreadOptions
#' @importFrom data.table data.table melt.data.table setnames is.data.table rbindlist .N .SD
#' @importFrom stats median
#' @importFrom utils modifyList
#' @export
FPCA <- function(formula, id.var, data, options = list()){
  # library(plyr); library(data.table); library(pipeR); library(assertthat)
  # formula = as.formula("y ~ t")
  # formula = as.formula("y + y2 ~ t")
  # formula = as.formula("y + y2 + y3 ~ t")
  # id.var = "sampleID"
  # data("irregularExData", package = "rfda")
  # data <- irregularExData %>>% data.table %>>% `[`( , `:=`(y2 = y*0.5 + rnorm(nrow(.)), y3 = y*rnorm(nrow(.))))
  # # data("sparseExData", package = "rfda")
  # # data <- sparseExData %>>% data.table %>>% `[`( , `:=`(y2 = y*0.5 + rnorm(nrow(.)), y3 = y*rnorm(nrow(.))))
  # options <- list()
  assert_that(is.formula(formula), is.character(id.var), length(id.var) == 1, is.data.frame(data))

  # check the formula
  chkFmLHS <- as.character(formula[[2]]) %>>% (grepl("[+a-zA-z_]", .)) %>>% all
  chkFmRHS <- as.character(formula[[3]]) %>>% (grepl("[+a-zA-z_]", .)) %>>% length %>>% `==`(1)
  chkFormla <- chkFmLHS || chkFmRHS
  message("Checking the formula...", appendLF = FALSE)
  message(ifelse(chkFormla, " Pass...", " Failed... Stop Now!"))

  # find the names of variables and name of variable indicating time points
  varName <- setdiff(all.vars(formula), as.character(formula[[3]]))
  timeVarName <- as.character(formula[[3]])

  # get the full options of FPCA and check
  default_FPCA_opts <- get_FPCA_opts(length(varName))
  optNamesUsed <- names(options) %in% default_FPCA_opts
  FPCA_opts <- modifyList(default_FPCA_opts, options[optNamesUsed]) %>>%
    chk_FPCA_opts(nrow(data))
  message(ifelse(all(optNamesUsed), "All options are checked...",
                 paste(names(options)[!optNamesUsed], collapse = ", ") %>>%
                   sprintf(fmt = "Ignoring the non-found options %s.")))

  # set the number of thread be used
  if (FPCA_opts$ncpus != 0)
    setThreadOptions(FPCA_opts$ncpus)

  # get the sparsity of data
  message("Checking and transforming data...")
  sparsity <- checkSparsity(data, id.var, timeVarName)

  # melt table to get a data.table to get a simple data.table and remove the NA, NaN and Inf.
  # additionally, give names for id.var and timeVarName
  dataDT <- melt.data.table(data.table(data), id.vars = c(id.var, timeVarName),
                            measure.vars = varName, variable.factor = TRUE) %>>%
    `[`(j = variable := as.integer(variable)) %>>% `[`(!is.na(value) & is.finite(value)) %>>%
    setnames(c(id.var, timeVarName) , c("subId", "timePnt"))

  # find the number of observations for each observed function
  subIdInsuffSize <- dataDT[, .(numObv = .N) ,by = .(subId, variable)][numObv <= 1] %>>%
    `$`(subId) %>>% unique
  # remove the
  dataDT <- dataDT[!subId %in% subIdInsuffSize]

  # binning data
  if (FPCA_opts$numBins == -1 || FPCA_opts$numBins >= 10)
  {
    # find the number of bins
    if (FPCA_opts$numBins == -1){
      numObv <- nrow(dataDT[variable == 1])
      numObvEach <- dataDT[variable == 1 , .(numObv = .N), by = timePnt] %>>% `$`(numObv)
      # find the number to decide whether to implement data binning
      decNum <- ifelse(sparsity == 0, stats::median(numObvEach), max(numObvEach))
      if (decNum > 400) {
        FPCA_opts$numBin <- 400
      } else if (numObv > 5000 && decNum > 20) {
        FPCA_opts$numBin <- min(decNum, ceiling(max(20, (5000-numObv)*19/2250+400)))
      }
    }

    # binning data and re-find the sparsity
    if (FPCA_opts$numBins > 0) {
      message("Start to implement data binning...")
      dataDT <- binData(dataDT, FPCA_opts$numBins)
      sparsity <- checkSparsity(dataDT[variable == 1], "subId", "timePnt")
    }
  }

  # get weight
  if (FPCA_opts$weight)
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

  # validate the list of user-specified mean functions
  validMFList <- !is.null(FPCA_opts$userMeanFuncList) && is.data.table(FPCA_opts$userMeanFuncList) &&
    all(c("timePnt", "value", "variable") %in% names(FPCA_opts$userMeanFuncList)) &&
    all(!unlist(FPCA_opts$userMeanFuncList[ , lapply(.SD, function(x) any(is.na(x) || is.infinite(x)))]))
  # get smoothing mean functions
  if (validMFList) {
    message("Use the user-specified mean functions...")
    MFRes <- lapply(list(sampledTimePnts, allTimePnts), function(v){
      FPCA_opts$userMeanFuncList %>>%
        `[`(j = .(timePnt = v, value = interp1(timePnt, value, v, "spline")),
             by = .(variable))
    })
    bwOptLocPoly1d <- NA
  } else {
    message("Get the smoothed mean functions...")
    # use gcv to get mean functions
    MFRes <- lapply(dataList, function(dat){
      # get the candidates of bandwidths
      bwCand <- bwCandChooser(dat, "subId", "timePnt", sparsity, FPCA_opts$bwKernel, 1)
      # get the optimal bandwidth with gcv
      bwOptLocPoly1d <- gcvLocPoly1d(bwCand, dat$timePnt, dat$value, dat$weight, FPCA_opts$bwKernel, 0, 1)
      # adjust the bandwidth
      bwOptLocPoly1d <- adjGcvBw1d(bwOptLocPoly1d, sparsity, FPCA_opts$bwKernel, 0)

      # Geometric mean of the minimum bandwidth and the GCV bandwidth
      if (FPCA_opts$bwMean == -1)
        bwOptLocPoly1d <- sqrt(find_max_diff_f(dat[["timePnt"]], 2) * bwOptLocPoly1d)

      meanFunc <- locPoly1d(bwOptLocPoly1d, dat$timePnt, dat$value, dat$weight,
                            sampledTimePnts, FPCA_opts$bwKernel, 0, 1)
      meanFuncDense <- locPoly1d(bwOptLocPoly1d, dat$timePnt, dat$value, dat$weight,
                                 allTimePnts, FPCA_opts$bwKernel, 0, 1)
      return(list(data.table(variable = dat$variable[1], bwOpt = bwOptLocPoly1d),
                  data.table(timePnt = sampledTimePnts, value = meanFunc, variable = unique(dat$variable)),
                  data.table(timePnt = allTimePnts, value = meanFuncDense, variable = unique(dat$variable))))
    }) %>>% rbindTransList
  }

  if (sapply(MFRes[2:3], function(dt) all(is.na(dt$value))) %>>% any)
    stop(paste0("The bandwidth of mean function is not appropriately!\n",
                "If it is chosen automatically, please provide your own mean functions."))

  # calculation of demeaned data
  demeanDataDT <- merge(dataDT, MFRes[[3]], by = c("timePnt", "variable"), suffixes = c(".ori", ".mean")) %>>%
    `[`(j = value := (value.ori - value.mean))

  # get raw covariance
  rawCov <- getRawCrCov(demeanDataDT)
  if (FPCA_opts$weight){
    if (sparsity == 0){
      rawCov <- rawCov[ , weight := NULL] %>>%
        merge(dataDT[, .(weight = weight[which.max(subId)]), by = .(variable,timePnt)],
              by.x = c("variable1", "t1"), by.y = c("variable", "timePnt"))
    } else {
      rawCov[ , weight := 1/cnt]
    }
  }

  if (FPCA_opts$methodNorm == "variance"){
    message("Start to normalize data with smoothed variances...")
    # get rawCov

    # get smoothed covariance

    # normalization

  } else if (FPCA_opts$methodNorm == "IQR") {
    message("Start to normalize data with smoothed IQRs...")
    # get IQR
    ## diff(qnorm(probs))
    # normalization

  } else {
    message("Not perform the normalization...")
  }

  # decide the number of FPC

  # calculation of FPC scores

  return(1)
}
