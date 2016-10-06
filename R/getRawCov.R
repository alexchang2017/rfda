getRawCov <- function(allTimePnts, dataDT, MDF_dt, sparsity){
  if (sparsity == 0){
    # tbd
  } else {
    rawCov <- expand.grid(t1 = allTimePnts, t2 = allTimePnts, subId = unique(dataDT$subId),
                          variable = unique(dataDT$variable)) %>>% data.table

    demeanDataDT <- merge(dataDT, MDF_dt, by = c("timePnt", "variable"), suffixes = c(".ori", ".mean")) %>>%
      `[`( , value := (value.ori - value.mean)) %>>% `[`( , .(timePnt, variable, subId, value))

    rawCovDT <- merge(rawCov, demeanDataDT, by.x = c("t1", "variable", "subId"),
          by.y = c("timePnt", "variable", "subId")) %>>% setnames("value", "value.t1") %>>%
      merge(demeanDataDT, by.x = c("t2", "variable", "subId"),
            by.y = c("timePnt", "variable", "subId")) %>>% setnames("value", "value.t2") %>>%
      `[`( , .(sse = sum(value.t1 * value.t2), cnt = .N), by = .(variable,t1,t2)) %>>%
      setorder(variable, t1, t2) %>>% `[`( , weight := 1)
  }
  return(rawCovDT)
}
