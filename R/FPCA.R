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
#' Return mean, cross-covariance, eigen functions and fpc scores.
#'
#' @param formula An object of class "formula", a description of FPCA model.
#' The RHS of `~` is the responses and LHS is time variables.
#' Notice that the names of variables must be a-z, A-Z or _.
#' @param id.var An character to indicate the subject id of data.
#' @param data An data.frame or data.table containing the variables in model.
#' @param options An list containing the options to fit FPCA model.
#' @return An list containing mean, cross-covariance, eigen functions,
#'   functional principal components scores, etc. Please see "Details" for more details.
#' @seealso \code{\link{get_FPCA_opts}} for input options.
#' @examples
#' \dontrun{
#' # sparse case
#' data("sparseExData", package = 'rfda')
#' fpcaResSparse <- FPCA(y ~ t, "sampleID", sparseExData)
#'
#' # regular with missing values case
#' data("irregularExData", package = 'rfda')
#' fpcaResIrregular <- FPCA(y ~ t, "sampleID", irregularExData)
#'
#' # regular case
#' data("regularExData", package = 'rfda')
#' fpcaResRegular <- FPCA(y ~ t, "sampleID", regularExData)
#' }
#' @importFrom plyr is.formula llply mapvalues
#' @importFrom RcppParallel setThreadOptions
#' @importFrom data.table data.table is.data.table as.data.table
#' @importFrom data.table melt.data.table dcast.data.table
#' @importFrom data.table setnames setDT setkey setkeyv
#' @importFrom data.table CJ tstrsplit .N .SD
#' @importFrom stats median qnorm cov cov2cor quantile
#' @importFrom utils modifyList combn type.convert
#' @export
FPCA <- function(formula, id.var, data, options = list()){
  # library(plyr); library(data.table); library(pipeR); library(assertthat); library(testthat)
  # formula = as.formula("y ~ t")
  # formula = as.formula("y + y2 ~ t")
  # formula = as.formula("y + y2 + y3 ~ t")
  # id.var = "sampleID"
  # data("irregularExData", package = "rfda")
  # data("sparseExData", package = "rfda")
  # data("regularExData", package = "rfda")
  # # data <- irregularExData %>>% data.table %>>% `[`( , `:=`(y2 = y*sin(t), y3 = y*cos(t)))
  # # data <- sparseExData %>>% data.table %>>% `[`( , `:=`(y2 = y*sin(t), y3 = y*cos(t)))
  # data <- regularExData %>>% data.table %>>% `[`( , `:=`(y2 = y*sin(t), y3 = y*cos(t)))
  # options <- list()

  #### check inputs ####
  assert_that(is.formula(formula), is.character(id.var), length(id.var) == 1, is.data.frame(data))
  chkFmLHS <- as.character(formula[[2]]) %>>% (grepl("[+a-zA-z_]", .)) %>>% all
  chkFmRHS <- as.character(formula[[3]]) %>>% (grepl("[+a-zA-z_]", .)) %>>% length %>>% `==`(1)
  chkFormla <- chkFmLHS || chkFmRHS
  message("Checking the formula...", appendLF = FALSE)
  message(ifelse(chkFormla, " Pass...", " Failed... Stop Now!"))
  numCurves <- length(unique(data[[id.var]]))
  if (numCurves < 2)
    stop("The number of curves cannot be less than 2!")

  # find the names of variables and name of variable indicating time points
  timeVarName <- as.character(formula[[3]])
  varName <- setdiff(all.vars(formula), timeVarName)

  #### get options ####
  # get the full options of FPCA and check
  default_FPCA_opts <- get_FPCA_opts(length(varName), length(unique(data[[id.var]])))
  optNamesUsed <- names(options) %in% names(default_FPCA_opts)
  FPCA_opts <- modifyList(default_FPCA_opts, options[optNamesUsed]) %>>% chk_FPCA_opts
  message(ifelse(all(optNamesUsed), "All options are checked...",
                 paste(names(options)[!optNamesUsed], collapse = ", ") %>>%
                   sprintf(fmt = "Ignoring the non-found options %s.")))
  timePntsRange <- range(data[[timeVarName]])
  if (!all(FPCA_opts$newdata >= timePntsRange[1] & FPCA_opts$newdata <= timePntsRange[2]))
    stop("The value of newdata (FPCA options) must in the range of input time points.")
  rm(default_FPCA_opts, optNamesUsed)

  #### check parallel ####
  # set the number of thread be used
  if (FPCA_opts$ncpus != 0)
    setThreadOptions(FPCA_opts$ncpus)

  #### check sparsity ####
  # get the sparsity of data
  message("Checking and transforming data...")
  sparsity <- checkSparsity(data, id.var, timeVarName)
  if (FPCA_opts$methodFPCS == "IN" && sparsity != 2)
    stop("The case methodFPCS = 'IN' is only avaiable for regular data, please use other methodFPCS.")

  #### transform data ####
  # melt table to get a data.table to get a simple data.table and remove the NA, NaN and Inf.
  # additionally, give names for id.var and timeVarName
  dataDT <- melt.data.table(data.table(data), id.vars = c(id.var, timeVarName),
                            measure.vars = varName, variable.factor = TRUE) %>>%
    (~ varNameMapping <- levels(.$variable)) %>>%
    `[`(j = variable := as.integer(variable)) %>>% `[`(!is.na(value) & is.finite(value)) %>>%
    setnames(c(id.var, timeVarName) , c("subId", "timePnt"))
  # function to map variable names to integers
  mapVarNames <- function(x) as.integer(mapvalues(x, varNameMapping, 1L:length(varNameMapping),
                                                  warn_missing = FALSE))

  # find the number of observations for each observed function
  numObvDT <- dataDT[, .(numObv = .N) ,by = .(subId, variable)]
  if (any(numObvDT$numObv <= 1)) {
    message("The number of points of some cuvres is less than 2, now remove them!")
    # remove the observations with insufficient data size
    subIdInsuffSize <- numObvDT[numObv <= 1, .(subId)] %>>% unlist %>>% unique
    warning("Remove the observation with ", id.var, " = ", paste0(subIdInsuffSize, collapse = ", "), ".")
    dataDT <- dataDT[!subId %in% subIdInsuffSize]
    numCurves <- length(unique(dataDT$subId))
    if (numCurves < 2)
      stop("The number of curves cannot be less than 2!")
    FPCA_opts$maxNumFPC <- min(FPCA_opts$maxNumFPC, numCurves - 2)
    rm(subIdInsuffSize)
  }

  #### binning data ####
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
      rm(numObv, numObvEach, criBinNum, minBinNum)
    }

    # binning data and re-find the sparsity
    if (FPCA_opts$numBins > 0) {
      message("Start to implement data binning...")
      dataDT <- binData(dataDT, FPCA_opts$numBins)
      sparsity <- checkSparsity(dataDT[variable == 1], "subId", "timePnt")
    }
  }

  #### get weight ####
  dataDtKeys <- c("variable", "timePnt")
  if (FPCA_opts$weight) {
    byVars <- switch(as.character(sparsity), "0" = c("variable", "subId"),
                     "1" = dataDtKeys, "2" = dataDtKeys)
    setkeyv(dataDT, byVars)
    dataDT <- merge(dataDT, dataDT[ , .(weight = 1/.N), by = byVars], by = byVars)
  } else {
    dataDT[ , weight := 1]
  }
  # change the keys of dataDT
  setkeyv(dataDT, dataDtKeys)

  #### get time points ####
  # Initializing allTimePnts is based on the unique time points of the pooled data + the unique
  # time points of "newdata", the output time grid.
  allTimePnts <- sort(unique(c(dataDT$timePnt, FPCA_opts$newdata)))
  # Initializing workTimePnts is based on the number of grid to be chosen in the range of
  # all time span.
  workTimePnts <- seq(timePntsRange[1], timePntsRange[2], length.out = FPCA_opts$numGrid)

  #### get mean functions ####
  # validate the list of user-specified mean functions
  validMFList <- !is.null(FPCA_opts$userMeanFunc) && is.data.table(FPCA_opts$userMeanFunc) &&
    all(c("timePnt", "value", "variable") %in% names(FPCA_opts$userMeanFunc)) &&
    all(!unlist(FPCA_opts$userMeanFunc[ , lapply(.SD, function(x) any(is.na(x) || is.infinite(x)))])) &&
    all(unique(FPCA_opts$userMeanFunc$variable) %in% varNameMapping)
  # get smoothing mean functions
  if (validMFList) {
    #### use input mean functions ####
    message("Use the user-specified mean functions...")
    FPCA_opts$userMeanFunc[ , variable := mapVarNames(variable)]
    MFRes <- c(list(data.table(variable = 1:length(varName), value = NA_real_)),
               llply(list(workTimePnts, allTimePnts), function(v){
                 FPCA_opts$userMeanFunc %>>%
                   `[`(j = .(timePnt = v, value = interp1(timePnt, value, v, "spline")),
                       by = .(variable))
               }))
  } else {
    #### calculate mean functions ####
    message("Get the smoothed mean functions...")
    if (is.null(FPCA_opts$bwMean)) {
      # use default
      FPCA_opts$bwMean <- data.table(variable = 1:length(varName), value = -1)
    } else {
      assert_that(is.data.frame(FPCA_opts$bwMean), all(c("variable", "value") %in% names(FPCA_opts$bwMean)),
                  all(FPCA_opts$bwMean$value > 0 || FPCA_opts$bwMean$value %in% c(-1, -2)),
                  is.numeric(FPCA_opts$bwMean$value), all(!is.na(FPCA_opts$bwMean$value)),
                  all(is.finite(FPCA_opts$bwMean$value)))
      setDT(FPCA_opts$bwMean)
      FPCA_opts$bwMean <- FPCA_opts$bwMean[ , variable := mapVarNames(as.character(variable))] %>>%
        setkeyv("variable") %>>% merge(data.table(variable = 1:length(varName), key = "variable"),
                                       by = "variable", all.y = TRUE) %>>%
        `[`(j = `:=`(variable = variable, value = ifelse(is.na(value), -1, value))) %>>% setorder(variable)
    }
    message("The setting of bandwidth used in smoothing mean functions is ")
    message("    ", paste0(varNameMapping, ": ", FPCA_opts$bwMean$value) %>>% paste0(collapse = ", "))

    # use gcv to get smoothed mean functions
    MFRes <- split(dataDT, by = "variable") %>>% llply(function(dat){
      meanBwOpt <- FPCA_opts$bwMean[variable == dat$variable[1]]$value
      if (meanBwOpt > 0) {
        bwOptLocPoly1d <- meanBwOpt
      } else {
        # get the candidates of bandwidths
        bwCand <- bwCandChooser(dat, "subId", "timePnt", sparsity, FPCA_opts$bwKernel, 1)
        # get the optimal bandwidth with gcv
        bwOptLocPoly1d <- with(dat, gcvLocPoly1d(bwCand, timePnt, value, weight, FPCA_opts$bwKernel, 0, 1))
        # adjust the bandwidth
        if (meanBwOpt ==  -1)
          bwOptLocPoly1d <- adjGcvBw(bwOptLocPoly1d, sparsity, FPCA_opts$bwKernel, 0)
      }

      meanFunc <- with(dat, locPoly1d(bwOptLocPoly1d, timePnt, value, weight,
                                      workTimePnts, FPCA_opts$bwKernel, 0, 1))
      meanFuncDense <- with(dat, locPoly1d(bwOptLocPoly1d, timePnt, value, weight,
                                           allTimePnts, FPCA_opts$bwKernel, 0, 1))
      return(list(data.table(variable = dat$variable[1], value = bwOptLocPoly1d),
                  data.table(timePnt = workTimePnts, value = meanFunc,
                             variable = unique(dat$variable), key = dataDtKeys),
                  data.table(timePnt = allTimePnts, value = meanFuncDense,
                             variable = unique(dat$variable), key = dataDtKeys)))
    }) %>>% rbindTransList
    message("The bandwidth of mean functions is \n    ",
            paste0(varNameMapping, ": ", sprintf("%.6f", MFRes[[1]]$value)) %>>% paste0(collapse = ", "))
  }

  #### check mean functions ####
  stopifnot(all(sapply(MFRes[2:3], function(dt) all(!is.na(dt$value)))))

  #### get demeaned values ####
  dataDT <- merge(dataDT, MFRes[[3]], by = dataDtKeys, suffixes = c("", ".mean")) %>>%
    `[`(j = value.demean := (value - value.mean))

  #### get raw cross-covariance ####
  rawCov <- getRawCrCov(dataDT)
  if (FPCA_opts$weight) {
    if (sparsity == 0) {
      rawCov <- rawCov[ , weight := NULL] %>>%
        merge(dataDT[, .(weight = weight[which.max(subId)]), by = dataDtKeys],
              by.x = c("variable1", "t1"), by.y = dataDtKeys)
    } else {
      rawCov[ , weight := 1/cnt]
    }
  }
  rawCovKeys <- paste0("variable", 1:2)
  setkeyv(rawCov, rawCovKeys)

  #### get cross-covariance functions ####
  # validate the list of user-specified cross-covariance functions
  cfListNames <- switch(as.integer(length(varName) == 1) + 1,
                        {rbind(data.table(variable1 = varName, variable2 = varName),
                               data.table(t(combn(varName, 2))) %>>% setnames(rawCovKeys)) %>>%
                            setorder(variable1, variable2) %>>% with(paste(variable1, variable2, sep = "-"))},
                        {paste(varName, varName, sep = "-")})

  validCFList <- !is.null(FPCA_opts$userCovFunc) && is.list(FPCA_opts$userCovFunc) &&
    !is.null(names(FPCA_opts$userCovFunc)) && all(cfListNames %in% names(FPCA_opts$userCovFunc)) &&
    all(sapply(FPCA_opts$userCovFunc, is.matrix)) &&
    all(!sapply(FPCA_opts$userCovFunc, function(x) any(is.na(x) || is.infinite(x)))) &&
    all(sapply(FPCA_opts$userCovFunc, function(x) nrow(x) == ncol(x))) &&
    all(!is.null(sapply(FPCA_opts$userCovFunc, colnames))) &&
    all(!is.null(sapply(FPCA_opts$userCovFunc, rownames))) &&
    all(sapply(FPCA_opts$userCovFunc, function(m) colnames(m) == rownames(m))) &&
    all(sapply(FPCA_opts$userCovFunc, function(x) {
      class(type.convert(colnames(x))) %in% c("integer", "numeric")
    }))
  rm(cfListNames)

  # get smoothing cross-covariance functions
  if (validCFList) {
    #### use input cross-covariance functions ####
    message("Use the user-specified cross-covariance functions...")
    message("Check the bwCov is given correctly...")
    assert_that(is.data.frame(FPCA_opts$bwCov),
                all(c(rawCovKeys, paste0("value", 1:2)) %in% names(FPCA_opts$bwCov)),
                is.numeric(FPCA_opts$bwCov$value1), all(!is.na(FPCA_opts$bwCov$value1)),
                all(is.finite(FPCA_opts$bwCov$value1)), is.numeric(FPCA_opts$bwCov$value2),
                all(!is.na(FPCA_opts$bwCov$value2)), all(is.finite(FPCA_opts$bwCov$value2)),
                all(FPCA_opts$bwCov$value1 > 0), all(FPCA_opts$bwCov$value2 > 0))
    setDT(FPCA_opts$bwCov, key = rawCovKeys)
    FPCA_opts$bwCov <- FPCA_opts$bwCov[variable1 == variable2] %>>%
      `[`(j = `:=`(variable1 = mapVarNames(as.character(variable1)),
                   variable2 = mapVarNames(as.character(variable2))))
    if (nrow(FPCA_opts$bwCov) != length(varName))
      stop("The bandwidths of covariance function (variable1 = variable2) must be ",
           "assigned with real positive values.")
    CFRes <- c(list(FPCA_opts$bwCov), list(llply(FPCA_opts$userCovFunc, function(m){
      grid <- as.numeric(rownames(m))
      interp2(grid, grid, m, workTimePnts, workTimePnts, 'spline') %>>%
        `colnames<-`(workTimePnts) %>>% `rownames<-`(workTimePnts)
    })))
  } else {
    #### calculate cross-covariance functions ####
    message("Get the smoothed cross-covariance functions...")
    if (is.null(FPCA_opts$bwCov)) {
      # use default
      FPCA_opts$bwCov <- unique(rawCov, by = rawCovKeys) %>>%
        `[`(j = value1 := ifelse(variable1 == variable2, -2, -3)) %>>%
        `[`(j = value2 := value1) %>>% `[`(j = .(variable1, variable2, value1, value2))
    } else {
      assert_that(is.data.frame(FPCA_opts$bwCov),
                  all(c(rawCovKeys, paste0("value", 1:2)) %in% names(FPCA_opts$bwCov)),
                  is.numeric(FPCA_opts$bwCov$value1), all(!is.na(FPCA_opts$bwCov$value1)),
                  all(is.finite(FPCA_opts$bwCov$value1)), is.numeric(FPCA_opts$bwCov$value2),
                  all(!is.na(FPCA_opts$bwCov$value2)), all(is.finite(FPCA_opts$bwCov$value2)),
                  all(FPCA_opts$bwCov$value1 > 0 || FPCA_opts$bwCov$value1 %in% c(-1, -2, -3)),
                  all(FPCA_opts$bwCov$value2 > 0 || FPCA_opts$bwCov$value2 %in% c(-1, -2, -3)))
      setDT(FPCA_opts$bwCov, key = rawCovKeys)

      FPCA_opts$bwCov <- FPCA_opts$bwCov[ , `:=`(variable1 = mapVarNames(as.character(variable1)),
                                                 variable2 = mapVarNames(as.character(variable2)))]
      if (any(FPCA_opts$bwCov$variable1 < FPCA_opts$bwCov$variable2))
        FPCA_opts$bwCov[ , tmp := variable1][ , `:=`(variable1 = variable2, variable2 = tmp)][ , tmp := NULL]

      FPCA_opts$bwCov <- merge(with(rawCov, CJ(variable1 = variable1, variable2 = variable2, unique = TRUE)),
                               FPCA_opts$bwCov, by = rawCovKeys, all.x = TRUE) %>>%
        `[`(j = value1 := ifelse(is.na(value1), ifelse(variable1 == variable2, -2, -3), value1)) %>>%
        `[`(j = value2 := ifelse(is.na(value2), ifelse(variable1 == variable2, -2, -3), value2)) %>>%
        setorder(variable1, variable2)
    }
    CFSettings <- with(FPCA_opts$bwCov, paste0(varNameMapping[variable1], "-", varNameMapping[variable2],
                                               ": (", value1, ",", value2, ")"))
    message("The setting of bandwidth used in smoothing cross-covariance functions is ")
    message("    ", paste0(CFSettings, collapse = ", "))
    rm(CFSettings)

    # smoothing diagonal cross-covariance matrix
    CFRes1 <- split(rawCov[variable1 == variable2], by = rawCovKeys) %>>% llply(function(dat){
      dat <- dat[t1 != t2]
      covBwOpt <- FPCA_opts$bwCov[variable1 == dat$variable1[1] & variable2 == dat$variable2[1],
                                  .(value1, value2)] %>>% unlist
      if (all(covBwOpt > 0)) {
        bwOptLocLinear2d <- unlist(covBwOpt)
      } else {
        # get the candidates of bandwidths
        bwCand <- bwCandChooser2(dat, sparsity, FPCA_opts$bwKernel, 1)
        # get the optimal bandwidth with gcv
        bwOptLocLinear2d <- with(dat, gcvLocLinear2d(bwCand, cbind(t1, t2), sse, weight, cnt,
                                                     FPCA_opts$bwKernel, FPCA_opts$bwNumGrid))
        # adjust the bandwidth
        if (all(covBwOpt ==  -1))
          bwOptLocLinear2d <- adjGcvBw(bwOptLocLinear2d, sparsity, FPCA_opts$bwKernel, 0)
      }

      list(data.table(variable1 = dat$variable1[1], variable2 = dat$variable2[1],
                      value1 = bwOptLocLinear2d[1], value2 = bwOptLocLinear2d[2]),
           with(dat, locLinear2d(bwOptLocLinear2d, cbind(t1, t2), sse, weight, cnt,
                                 workTimePnts, workTimePnts, FPCA_opts$bwKernel)))
    })

    # smoothing cross-covariance matrix
    CFRes2 <- split(rawCov[variable1 != variable2], by = rawCovKeys) %>>% llply(function(dat){
      covBwOpt <- FPCA_opts$bwCov[variable1 == dat$variable1[1] & variable2 == dat$variable2[1],
                                  .(value1, value2)] %>>% unlist
      if (all(covBwOpt > 0)) {
        bwOptLocLinear2d <- unlist(covBwOpt)
      } else if (all(covBwOpt == -3)) {
        bwOptLocLinear2d <- c(CFRes1[[dat$variable1[1]]][[1]]$value1, CFRes1[[dat$variable2[1]]][[1]]$value2)
      } else {
        # get the candidates of bandwidths
        bwCand <- bwCandChooser3(dat, sparsity, FPCA_opts$bwKernel, 1)
        # get the optimal bandwidth with gcv
        bwOptLocLinear2d <- with(dat, gcvLocLinear2d(bwCand, cbind(t1, t2), sse, weight, cnt,
                                                     FPCA_opts$bwKernel, FPCA_opts$bwNumGrid))
        # adjust the bandwidth
        if (all(covBwOpt ==  -1))
          bwOptLocLinear2d <- adjGcvBw(bwOptLocLinear2d, sparsity, FPCA_opts$bwKernel, 0)
      }

      list(data.table(variable1 = dat$variable1[1], variable2 = dat$variable2[1],
                      value1 = bwOptLocLinear2d[1], value2 = bwOptLocLinear2d[2]),
           with(dat, locLinear2d(bwOptLocLinear2d, cbind(t1, t2), sse, weight, cnt,
                                 workTimePnts, workTimePnts, FPCA_opts$bwKernel)))
    })

    # combine the smoothing results
    CFResNames <- c(names(CFRes1), names(CFRes2)) %>>% strsplit("\\.") %>>%
      sapply(function(v) paste0(varNameMapping[as.integer(v)], collapse = "-"))
    CFRes <- transCFRes(c(CFRes1, CFRes2), workTimePnts, CFResNames)
    rm(CFRes1, CFRes2)

    CFResBwValues <- with(CFRes[[1]], paste(sprintf("%.4f", value1), sprintf("%.4f", value2), sep = ","))
    message("The bandwidth of cross-covariance functions is \n    ",
            paste0(CFResNames, ": (", CFResBwValues, ")") %>>% paste0(collapse = ", "))
    rm(CFResBwValues)
  }

  #### check cross-covariance functions ####
  stopifnot(all(sapply(CFRes[[2]], function(x) all(!is.na(x)))))

  #### estimate measurement error ####
  measErrVarDT <- getMeasErr(rawCov, CFRes[[1]], sparsity, FPCA_opts) %>>%
    (ifelse(is.na(.), FPCA_opts$minMeasErr, .)) %>>%
    (data.table(variable = as.integer(names(.)), value = .))
  message("Get the variance of measurement error: ",
          paste0(sprintf("%s: %.4f", varNameMapping[measErrVarDT$variable], measErrVarDT$value),
                 collapse = ", "))

  #### Leaving out the data in the boundary ####
  if (FPCA_opts$outPercent > 0) {
    message("Leave out ", sprintf("%.2f%%", FPCA_opts$outPercent * 100), " data in the boundary...")
    # get the range to keep data
    userRange <- quantile(c(dataDT$timePnt, FPCA_opts$newdata), FPCA_opts$outPercent %>>%
                            (c(./2, 1-./2)))
    # leave outPersent data
    dataDT <- dataDT[timePnt >= userRange[1] & timePnt <= userRange[2]]
    # get the new time points
    allTimePntsOld <- allTimePnts
    workTimePntsOld <- workTimePnts
    allTimePnts <- sort(unique(c(dataDT$timePnt, FPCA_opts$newdata)))
    workTimePnts <- seq(min(allTimePnts), max(allTimePnts), length.out = FPCA_opts$numGrid)
    # get the new mean functions
    MFRes[[2]] <- MFRes[[2]][ , .(value = interp1(workTimePntsOld, value, workTimePnts, "linear"),
                                  timePnt = workTimePnts), by = .(variable)]
    MFRes[[3]] <- MFRes[[3]][ , .(value = interp1(allTimePntsOld, value, allTimePnts, "linear"),
                                  timePnt = allTimePnts), by = .(variable)]
    # get the new cross-covariance functions
    CFRes[[2]] <- llply(CFRes[[2]], function(mat){
      interp2(workTimePntsOld, workTimePntsOld, mat, workTimePnts, workTimePnts, "linear") %>>%
        ((. + t(.)) / 2) %>>% `dimnames<-`(list(workTimePnts, workTimePnts))
    })
  }

  #### get the cross-covariance functions ####
  message("Get the cross-covariance functions...")
  # initialize cross-covariance function matrix
  CFMatNames <- rep(1:length(varName), each = FPCA_opts$numGrid)
  CFMat <- matrix(0, FPCA_opts$numGrid * length(varName), FPCA_opts$numGrid * length(varName),
                  dimnames = list(CFMatNames, CFMatNames))
  # find the corresponding variables in the list of cross-covariance functions
  orderCFResDT <- tstrsplit(names(CFRes[[2]]), "-") %>>% as.data.table %>>%
    setnames(rawCovKeys) %>>% `[`(j = lapply(.SD, mapVarNames))
  # fill the cross-covariance function matrix
  for (i in 1:nrow(orderCFResDT)) {
    idxList <- llply(orderCFResDT[i, ], function(x){
      seq((x-1)*FPCA_opts$numGrid+1, x*FPCA_opts$numGrid)
    })
    CFMat[idxList[[1]], idxList[[2]]] <- CFRes[[2]][[i]]
    if (diff(unlist(orderCFResDT[i, ])) != 0)
      CFMat[idxList[[2]], idxList[[1]]] <- t(CFRes[[2]][[i]])
  }
  rm(orderCFResDT, idxList)

  #### get normalized data if needed ####
  # start to normalize data and get the working cross-correlation matrix for find the FPC scores
  if (FPCA_opts$methodNorm == "smoothCov") {
    message("Normalize data with smoothed variances...")
    # acquire names of smoothed variance
    smoothVarNames <- names(CFRes[[2]][1:length(varName)]) %>>% strsplit("-") %>>% sapply(`[`, 1)
    # acquire smoothed variance
    smoothVarMat <- sapply(CFRes[[2]][1:length(varName)], diag)
    # let negative variance be the minimum positive smoothed variance
    smoothVarDT <- smoothVarMat %>>% `colnames<-`(smoothVarNames) %>>% data.table(keep.rownames = TRUE) %>>%
      melt.data.table("rn") %>>% `[`(j = `:=`(variable = mapVarNames(variable), timePnt = as.numeric(rn))) %>>%
      `[`(j = value2 := ifelse(value < 0, min(value[value > 0]), value), by = .(variable)) %>>%
      setkeyv(dataDtKeys) %>>% setorder(variable, timePnt)

    # get the working cross-correlation matrix to find the FPC scores
    message("Get the working cross-correlation functions...")
    diag(CFMat) <- smoothVarDT$value2
    varDT <- smoothVarDT[ , .(timePnt, variable, value = value2)] %>>% setkey(variable, timePnt)
    CFMat2 <- cov2cor(CFMat)

    # get new dataDT with normalized values
    smoothVarDT2 <- smoothVarDT[ , .(value = interp1(timePnt, value, allTimePnts, "spline"),
                                     timePnt = allTimePnts), by = .(variable)] %>>%
      `[`(j = value := ifelse(value < 0, min(value[value > 0]), value), by = .(variable)) %>>%
      setkeyv(dataDtKeys) %>>% setorder(variable, timePnt)
    dataDT <- merge(dataDT, smoothVarDT2, by = dataDtKeys, suffixes = c("", ".var"),
                    all.x = TRUE) %>>% `[`(j = value := (value.demean / sqrt(value.var)))
    rm(smoothVarNames, smoothVarMat, smoothVarDT, smoothVarDT2)
  } else if (FPCA_opts$methodNorm == "quantile") {
    message("Normalize data with smoothed IQRs...")
    # acquire quantiles and convert it to variance estimator
    quantRes <- split(dataDT, by = "variable") %>>% llply(function(dat){
      bwOpt <- MFRes[[1]][variable == dat$variable[1]]$value
      quantFuncs <- with(dat, locQuantPoly1d(bwOpt, FPCA_opts$quantileProbs, timePnt, value,
                                             weight, workTimePnts, FPCA_opts$bwKernel, 0, 1))
      quantFuncsDense <- with(dat, locQuantPoly1d(bwOpt, FPCA_opts$quantileProbs, timePnt, value,
                                                  weight, allTimePnts, FPCA_opts$bwKernel, 0, 1))
      iqr <- (quantFuncs[ , 2] - quantFuncs[ , 1]) / diff(qnorm(FPCA_opts$quantileProbs))
      iqrDense <- (quantFuncsDense[ , 2] - quantFuncsDense[ , 1]) / diff(qnorm(FPCA_opts$quantileProbs))
      return(list(data.table(timePnt = workTimePnts, value = iqr**2, variable = unique(dat$variable),
                             key = dataDtKeys),
                  data.table(timePnt = allTimePnts, value = iqrDense**2,
                             variable = unique(dat$variable), key = dataDtKeys)))
    }) %>>% rbindTransList %>>% llply(function(DT) setorder(DT, variable, timePnt))

    # get the working cross-correlation matrix to find the FPC scores
    message("Get the working cross-correlation functions...")
    diag(CFMat) <- quantRes[[1]]$value
    varDT <- quantRes[[1]] %>>% setkey(variable, timePnt)
    CFMat2 <- cov2cor(CFMat)

    # get new dataDT with normalized values
    dataDT <- merge(dataDT, quantRes[[2]], by = dataDtKeys, suffixes = c("", ".var"),
                    all.x = TRUE) %>>% `[`(j = value := (value.demean / sqrt(value.var)))
    rm(quantRes)
  } else if (FPCA_opts$methodNorm == "no") {
    message("Not perform the normalization...")
    dataDT[ , value.var := 1]
    message("Use original cross-covariance functions...")
    varDT <- CJ(timePnt = workTimePnts, variable = 1:length(varName), value = 1) %>>%
      setkey(variable, timePnt)
    CFMat2 <- CFMat
  }

  #### find the eigen functions and the corresponding eigenvalues ####
  message("Get the eigen functions and the corresponding eigenvalues...")
  eigRes <- getEigRes(CFMat2, as.integer(rownames(CFMat2)), workTimePnts, MFRes[[2]]$value, allTimePnts)

  #### get fitted cross-covariance/cross-correlation functions ####
  message("Get the fitted cross-covariance functions...")
  fittedCFMat2 <- eigRes$eigFuncsWork %*% diag(eigRes$eigVals) %*% t(eigRes$eigFuncsWork)
  fittedCFMat <- sweep(fittedCFMat2, 1, sqrt(varDT$value), "*") %>>% sweep(2, sqrt(varDT$value), "*")

  #### calculate FPC scores ####
  message("Get the functional principal component scores...")
  if (FPCA_opts$methodFPCS == "CE") {
    if (FPCA_opts$rho != "no") {
      if (numCurves > 2048) {
        # getRho
      } else {
        # getRho
      }
    }
    # getFpcScoresCE
  } else if (FPCA_opts$methodFPCS == "IN") {
    valueCol <- ifelse(FPCA_opts$methodNorm == "no", "value.demean", "value")
    splitVar <- rep(1:length(varName), each = length(allTimePnts))
    yMat <- dcast.data.table(dataDT, variable + timePnt ~ subId, sum, value.var = valueCol) %>>%
      setkeyv(dataDtKeys) %>>% `[`(j = `:=`(eval(dataDtKeys), list(NULL, NULL))) %>>% as.matrix
    fpcScores <- getFpcScoresIN(allTimePnts, splitVar, yMat, eigRes$eigFuncs, FPCA_opts$shrink,
                                eigRes$eigVals, measErrVarDT$value)
  } else if (FPCA_opts$methodFPCS == "LS") {
    # getFpcScoresLS
  } else if (FPCA_opts$methodFPCS == "WLS") {
    # getFpcScoresWLS
  }

  #### calculate FVE ####
  FVE <- eigRes$eigVals %>>% (cumsum(.) / sum(.))
  # find the number of FPC with FVE criterion
  numFPC_FVE <- Position(function(x) x > FPCA_opts$FVE_threshold, FVE)
  message("The number of FPC decided by FVE is ", numFPC_FVE, " with threshold ",
          sprintf("%.2f%%", FPCA_opts$FVE_threshold * 100), "...")

  #### decide the number of FPC ####
  if (is.character(FPCA_opts$numFPC))
    message("The criteria to decide the number of FPC is ", FPCA_opts$numFPC, "...")
  if (is.numeric(FPCA_opts$numFPC))
    message("The the number of FPC is assigned by user with value ", FPCA_opts$numFPC, "...")
  if (is.character(FPCA_opts$numFPC) && FPCA_opts$numFPC != "AIC_R") {
    # find the number of FPC with AIC or BIC
    if (FPCA_opts$numFPC %in% c("AIC", "BIC")) {
      # resNumFPC_IC <- getNumFpcByIC()
      # assign(FPCA_opts$numFPC, resNumFPC_IC$criterion)
      # numFPC <- resNumFPC_IC$numFPC
    }

    # numFPC <- numFPC_FVE
    # reset numFPC if it is NA
    # if (is.na(numFPC)) {
    #   FPCA_opts$numFPC <- "FVE"
    #   numFPC <- numFPC_FVE
    # }
  } else if (is.numeric(FPCA_opts$numFPC) || FPCA_opts$numFPC == "AIC_R") {
    numFPC <- ifelse(FPCA_opts$numFPC == "AIC_R", (FPCA_opts$numGrid - 2L) * length(varName), FPCA_opts$numFPC)
    if (numFPC > length(eigRes$eigVals))
      warning("There are only ", length(eigRes$eigVals), " avaiable eigenfunctions, it is less than ",
              numFPC, ", so reset numFPC to ", length(eigRes$eigVals), ".")
    numFPC <- min(length(eigRes$eigVals), numFPC)
  }
  return(1)
}
