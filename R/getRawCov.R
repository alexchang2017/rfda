#' @importFrom data.table setorder
getRawCov <- function(dataDT, MDF_dt){
  demeanDataDT <- merge(dataDT, MDF_dt, by = c("timePnt", "variable"), suffixes = c(".ori", ".mean")) %>>%
    `[`( , value := (value.ori - value.mean)) %>>% `[`( , .(timePnt, variable, subId, value))

  rawCovDT <- demeanDataDT[ , .(t1 = rep(timePnt, length(timePnt)), t2 = rep(timePnt, each=length(timePnt)),
                    value.t1 = rep(value, length(timePnt)), value.t2 = rep(value, ach=length(timePnt))),
                by = .(variable, subId)] %>>%
    `[`( , .(sse = sum(value.t1 * value.t2), cnt = .N), by = .(variable,t1,t2)) %>>%
    setorder(variable, t1, t2) %>>% `[`( , weight := 1)
  return(rawCovDT)
}
