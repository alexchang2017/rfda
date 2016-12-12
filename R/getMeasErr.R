#' Estimate the variance of measurement error
#' @param rawCov A data.table with columns variable1, variable2, t1, t2, sse, cnt and weight.
#'   The sse and cnt are the sum of squared error and count at (t1, t2) for the
#'   cross-covariance surface of variable1 and variable2 and weight is corresponding weight.
#' @param smoothCfBwDT A data.table with columns: variable1, variable2, value1 and value2.
#'   The columns, value1 and value2, are the bandwidth values for the cross-covariance surface
#'   of variable1 and variable2.
#' @param sparsity An integer. The sparsity of data.
#' @param FPCA_opts The options for FPCA.
#' @noRd
getMeasErr <- function(rawCov, smoothCfBwDT, sparsity, FPCA_opts){
  timeRange <- range(rawCov$t1)
  workTimePnts <- seq(timeRange[1], timeRange[2], length.out = FPCA_opts$bwNumGrid)
  cutTimePnts <- quantile(workTimePnts, FPCA_opts$measErrOut %>>% (c(./2, 1-./2)))
  idx <- workTimePnts >= cutTimePnts[1] & workTimePnts <= cutTimePnts[2]

  return(mapply(function(dt1, dt2){
    varFunc <- with(dt1, locPoly1d(dt2$value1, t1, sse / cnt, weight, workTimePnts, FPCA_opts$bwKernel, 0, 1))
    rotVarFunc <- with(dt1, locLinearRotate2d(c(dt2$value1, dt2$value2), cbind(t1, t2), sse, weight, cnt,
                                              cbind(workTimePnts, workTimePnts), FPCA_opts$bwKernel))
    if (is.na(varFunc) || is.na(rotVarFunc))
      return(NA)
    return(max(trapz(workTimePnts[idx], rotVarFunc[idx] - varFunc[idx]) / diff(timeRange), FPCA_opts$minMeasErr))
  }, split(rawCov[t1 != t2 & variable1 == variable2], by = "variable1"),
  split(smoothCfBwDT[variable1 == variable2], by = "variable1")))
}
