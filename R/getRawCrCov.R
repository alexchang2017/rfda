#' Get the raw cross-covariance surface
#'
#' The usage of this function can be found in the example of \code{\link{gcvLocLinear2d}}.
#'
#' @param dataDT A data.table with five columns \code{timePnt}, \code{value.demean}, \code{variable} and
#'   \code{subId}. \code{value.demean} is the observation minus smoothed mean.
#'   \code{timePnt} is corresponding time points. \code{variable} is the name of observed variable.
#'   \code{subId} is the id of subject.
#' @return A expand grid of cross-covariance surface.
#' @export
#' @importFrom data.table setorder setnames setkey .N
#' @importFrom plyr llply
getRawCrCov <- function(dataDT){
  # geneerate the all combinations of t1,t2 and varaibles
  baseDT <- dataDT[ , .(t1 = rep(timePnt, length(timePnt)), t2 = rep(timePnt, each=length(timePnt)),
                        value.var1 = rep(value.demean, length(timePnt))), by = .(variable, subId)]
  # set the keys of data.table
  setkey(baseDT, subId, t2)
  setkey(dataDT, subId, timePnt)
  # calculation of raw cross-covariance
  rawCrCovDT <- split(dataDT, by = "variable") %>>% llply(function(df){
    merge(baseDT[variable >= df$variable[1]], df, suffixes = c("1", "2"),
          by.x = c("subId", "t2"), by.y = c("subId", "timePnt"))
  }) %>>% rbindlist %>>% setnames("value.demean", "value.var2") %>>%
    `[`(j = .(sse = sum(value.var1 * value.var2), cnt = .N), by = .(variable1, variable2, t1, t2)) %>>%
    setorder(variable1, variable2, t1, t2) %>>% `[`(j = weight := 1)
  return(rawCrCovDT)
}
