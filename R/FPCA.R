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

#' Get the raw cross-covariance surface
#'
#' The usage of this function can be found in the example of \code{\link{gcvLocLinear2d}}.
#'
#' @param demeanDataDT A data.table with five columns \code{timePnt}, \code{value}, \code{variable} and
#'   \code{subId}. \code{value} is the observation minus smoothed mean.
#'   \code{timePnt} is corresponding time points. \code{variable} is the name of observed variable.
#'   \code{subId} is the id of subject.
#' @return A expand grid of cross-covariance surface.
#' @export
#' @importFrom data.table setorder .N setDT setnames
#' @importFrom plyr ddply
getRawCrCov <- function(demeanDataDT){
  # geneerate the all combinations of t1,t2 and varaibles
  baseDT <- demeanDataDT[ , .(t1 = rep(timePnt, length(timePnt)), t2 = rep(timePnt, each=length(timePnt)),
                              value.var1 = rep(value, length(timePnt))), by = .(variable, subId)]
  # calculation of raw cross-covariance
  rawCrCovDT <- do.call("ddply", list(demeanDataDT, "variable", function(df){
    merge(baseDT[variable >= df$variable[1]], df, suffixes = c("1", "2"),
          by.x = c("subId", "t2"), by.y = c("subId", "timePnt"))
  })) %>>% setDT %>>% `[`(j = variable := NULL) %>>% setnames("value", "value.var2") %>>%
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
#' @importFrom plyr is.formula llply laply dlply
#' @importFrom RcppParallel setThreadOptions
#' @importFrom data.table data.table melt.data.table setnames is.data.table .N .SD setDT
#' @importFrom stats median
#' @importFrom utils modifyList combn
#' @export
FPCA <- function(formula, id.var, data, options = list()){
  # library(plyr); library(data.table); library(pipeR); library(assertthat)
  # formula = as.formula("y ~ t")
  # formula = as.formula("y + y2 ~ t")
  # formula = as.formula("y + y2 + y3 ~ t")
  # id.var = "sampleID"
  # data("irregularExData", package = "rfda")
  # data <- irregularExData %>>% data.table %>>% `[`( , `:=`(y2 = y*sin(t), y3 = y*cos(t)))
  # # data("sparseExData", package = "rfda")
  # # data <- sparseExData %>>% data.table %>>% `[`( , `:=`(y2 = y*sin(t), y3 = y*cos(t)))
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
  optNamesUsed <- names(options) %in% names(default_FPCA_opts)
  FPCA_opts <- modifyList(default_FPCA_opts, options[optNamesUsed]) %>>% chk_FPCA_opts
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
    (~ varNameMapping <- levels(.$variable)) %>>%
    `[`(j = variable := as.integer(variable)) %>>% `[`(!is.na(value) & is.finite(value)) %>>%
    setnames(c(id.var, timeVarName) , c("subId", "timePnt"))

  # find the number of observations for each observed function
  subIdInsuffSize <- dataDT[, .(numObv = .N) ,by = .(subId, variable)][numObv <= 1] %>>%
    `$`(subId) %>>% unique
  # remove the
  dataDT <- dataDT[!subId %in% subIdInsuffSize]

  # binning data
  if (FPCA_opts$numBins == -1 || FPCA_opts$numBins >= 10) {
    # find the number of bins
    if (FPCA_opts$numBins == -1) {
      numObv <- nrow(dataDT[variable == 1])
      numObvEach <- dataDT[variable == 1 , .(numObv = .N), by = timePnt] %>>% `$`(numObv)
      # find the number to decide whether to implement data binning
      criBinNum <- ifelse(sparsity == 0, median(numObvEach), max(numObvEach))
      # minimun bin number
      minBinNum <- min(criBinNum, ceiling(max(20, (5000 - numObv)*19/2250+400)))
      FPCA_opts$numBin <- ifelse(criBinNum > 400, 400, ifelse(numObv > 5000 && criBinNum > 20, minBinNum, -1))
    }

    # binning data and re-find the sparsity
    if (FPCA_opts$numBins > 0) {
      message("Start to implement data binning...")
      dataDT <- binData(dataDT, FPCA_opts$numBins)
      sparsity <- checkSparsity(dataDT[variable == 1], "subId", "timePnt")
    }
  }

  # get weight
  if (FPCA_opts$weight) {
    byVars <- switch(as.character(sparsity), "0" = c("variable", "subId"),
                     "1" = c("variable", "timePnt"), "2" = c("variable", "timePnt"))
    dataDT <- merge(dataDT, dataDT[ , .(weight = 1/.N), by = byVars], by = byVars)
  } else {
    dataDT[ , weight := 1]
  }

  # Initializing allTimePnts is based on the unique time points of the pooled data + the unique
  # time points of "newdata", the output time grid.
  allTimePnts <- sort(unique(c(dataDT$timePnt, FPCA_opts$newdata)))
  # Initializing sampledTimePnts is based on the number of grid to be chosen in the range of
  # all time span.
  sampledTimePnts <- seq(min(allTimePnts), max(allTimePnts), length.out = FPCA_opts$numGrid)

  # validate the list of user-specified mean functions
  validMFList <- !is.null(FPCA_opts$userMeanFunc) && is.data.table(FPCA_opts$userMeanFunc) &&
    all(c("timePnt", "value", "variable") %in% names(FPCA_opts$userMeanFunc)) &&
    all(!unlist(FPCA_opts$userMeanFunc[ , lapply(.SD, function(x) any(is.na(x) || is.infinite(x)))])) &&
    all(unique(FPCA_opts$userMeanFunc$variable) %in% varNameMapping)
  # get smoothing mean functions
  if (validMFList) {
    message("Use the user-specified mean functions...")
    FPCA_opts$userMeanFunc[ , variable := match(variable, varNameMapping)]
    MFRes <- c(list(data.table(variable = 1:length(varName), value = NA_real_)),
      llply(list(sampledTimePnts, allTimePnts), function(v){
      FPCA_opts$userMeanFunc %>>%
        `[`(j = .(timePnt = v, value = interp1(timePnt, value, v, "spline")),
             by = .(variable))
    }))
  } else {
    message("Get the smoothed mean functions...")
    if (is.null(FPCA_opts$bwMean)) {
      # use default
      FPCA_opts$bwMean <- data.table(variable = 1:length(varName), value = -1)
    } else {
      assert_that(is.data.frame(FPCA_opts$bwMean), all(c("variable", "value") %in% names(FPCA_opts$bwMean)),
                  all(FPCA_opts$bwMean$value > 0 || FPCA_opts$bwMean$value %in% c(-1, -2)),
                  is.numeric(FPCA_opts$bwMean$value), all(!is.na(FPCA_opts$bwMean$value)),
                  all(is.finite(FPCA_opts$bwMean$value)))
      FPCA_opts$bwMean$variable <- as.character(FPCA_opts$bwMean$variable)
      setDT(FPCA_opts$bwMean)
      FPCA_opts$bwMean <- data.table(variable = varNameMapping) %>>%
        merge(FPCA_opts$bwMean, by = "variable", all.x = TRUE) %>>%
        `[`(j = `:=`(variable = match(variable, varNameMapping),
                     value = ifelse(is.na(value), -1, value))) %>>%
        setorder(variable)
    }
    message("The setting of bandwidth used in smoothing mean functions is ")
    message("    ", paste0(varNameMapping, ": ", FPCA_opts$bwMean$value) %>>% paste0(collapse = ", "))

    bwValues <- split(FPCA_opts$bwMean, FPCA_opts$bwMean$variable)
    # use gcv to get mean functions
    MFRes <- do.call("dlply", list(dataDT, "variable", function(dat){
      bwOpt <- FPCA_opts$bwMean[variable == dat$variable[1]]$value
      if (bwOpt > 0) {
        bwOptLocPoly1d <- bwOpt
      } else {
        # get the candidates of bandwidths
        bwCand <- bwCandChooser(dat, "subId", "timePnt", sparsity, FPCA_opts$bwKernel, 1)
        # get the optimal bandwidth with gcv
        bwOptLocPoly1d <- gcvLocPoly1d(bwCand, dat$timePnt, dat$value, dat$weight, FPCA_opts$bwKernel, 0, 1)
        # adjust the bandwidth
        bwOptLocPoly1d <- adjGcvBw(bwOptLocPoly1d, sparsity, FPCA_opts$bwKernel, 0)
      }

      # Geometric mean of the minimum bandwidth and the GCV bandwidth
      if (bwOpt == -1)
        bwOptLocPoly1d <- sqrt(find_max_diff_f(dat[["timePnt"]], 2) * bwOptLocPoly1d)

      meanFunc <- locPoly1d(bwOptLocPoly1d, dat$timePnt, dat$value, dat$weight,
                            sampledTimePnts, FPCA_opts$bwKernel, 0, 1)
      meanFuncDense <- locPoly1d(bwOptLocPoly1d, dat$timePnt, dat$value, dat$weight,
                                 allTimePnts, FPCA_opts$bwKernel, 0, 1)
      return(list(data.table(variable = dat$variable[1], value = bwOptLocPoly1d),
                  data.table(timePnt = sampledTimePnts, value = meanFunc, variable = unique(dat$variable)),
                  data.table(timePnt = allTimePnts, value = meanFuncDense, variable = unique(dat$variable))))
    })) %>>% rbindTransList
    message("The bandwidth of mean functions is \n    ",
            paste0(varNameMapping, ": ", sprintf("%.6f", MFRes[[1]]$value)) %>>% paste0(collapse = ", "))
  }

  if (laply(MFRes[2:3], function(dt) all(is.na(dt$value))) %>>% any)
    stop(paste0("The bandwidth of mean function is not appropriately!\n",
                "If it is chosen automatically, please provide your own mean functions."))

  # calculation of demeaned data
  demeanDataDT <- merge(dataDT, MFRes[[3]], by = c("timePnt", "variable"), suffixes = c(".ori", ".mean")) %>>%
    `[`(j = value := (value.ori - value.mean))

  # get raw covariance
  rawCov <- getRawCrCov(demeanDataDT)
  if (FPCA_opts$weight) {
    if (sparsity == 0) {
      rawCov <- rawCov[ , weight := NULL] %>>%
        merge(dataDT[, .(weight = weight[which.max(subId)]), by = .(variable,timePnt)],
              by.x = c("variable1", "t1"), by.y = c("variable", "timePnt"))
    } else {
      rawCov[ , weight := 1/cnt]
    }
  }

  # validate the list of user-specified covariance surface
  # if (length(varName) == 1) {
  #   cfListNames <- paste(varName, varName, sep = "-")
  # } else {
  #   cfListNames <- combn(varName, 2) %>>% (paste(.[1,], .[2,], sep = "-"))
  # }
  # validCFList <- !is.null(FPCA_opts$userCovFunc) && is.list(FPCA_opts$userCovFunc) &&
  #   !is.null(names(FPCA_opts$userCovFunc)) && all(names(FPCA_opts$userCovFunc) %in% cfListNames) &&
  #   all(laply(FPCA_opts$userCovFunc, is.matrix)) && all(laply(FPCA_opts$userCovFunc, dim) == 2) &&
  #   all(laply(FPCA_opts$userCovFunc, function(x) nrow(x) == ncol(x))) &&
  #   all(!is.null(laply(FPCA_opts$userCovFunc, colnames))) &&
  #   all(!is.null(laply(FPCA_opts$userCovFunc, rownames))) &&
  #   all(laply(FPCA_opts$userCovFunc, rownames) == laply(FPCA_opts$userCovFunc, colnames))

  # get smoothing covariance surface
  if (validCFList) {
    # message("Use the user-specified covariance surface...")
    # userGrid <- llply(FPCA_opts$userCovFunc, function(m) as.numeric(rownames(m)))
    # outBwDT <- expand.grid(1:length(varName), 1:length(varName)) %>>% data.table %>>%
    #   setnames(paste0("variable", 1:2)) %>>% `[`(j = value := NA_real_)
    # CSRes <- c(list(outBwDT), mapply(function(m, grid){
    #   interp2(grid, grid, m, sampledTimePnts, sampledTimePnts, 'spline')
    #   }, FPCA_opts$userCovFunc, userGrid, SIMPLIFY = FALSE))
  } else {
    # split raw covariance and smooth separately
    # splitFactor <- paste(rawCov$variable1, rawCov$variable2, sep = "-")
    # dataList2 <- split(rawCov, splitFactor)
    # get smoothing mean functions
  }

  # if (any(!laply(CSRes[[2]], is.na)))
  #   stop(paste0("The bandwidth of covariance function is not appropriate!\n",
  #               "If it is chosen automatically, please provide your own covariance surface."))

  if (FPCA_opts$methodNorm == "smoothCov") {
    message("Start to normalize data with smoothed variances...")
    # get rawCov

    # get smoothed covariance

    # normalization

  } else if (FPCA_opts$methodNorm == "quantile") {
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
