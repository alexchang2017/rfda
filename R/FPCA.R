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
#' @param FPCA_opts An list containing the options to fit FPCA model.
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
#' @importFrom stringr str_detect
#' @importFrom data.table data.table melt.data.table
#' @importFrom stats quantile
#' @export
FPCA <- function(formula, id.var, data, FPCA_opts = get_FPCA_opts()){
  # formula = as.formula("y ~ t")
  # id.var = "sampleID"
  # data = irregularExData %>>% data.table %>>% `[`( , y2 := y*2 + rnorm(nrow(.)))
  assert_that(is.formula(formula), is.character(id.var),
              length(id.var) == 1, is.data.frame(data))

  chkFmLHS <- as.character(formula[[2]]) %>>% str_detect("[+a-zA-z_]") %>>% all
  chkFmRHS <- as.character(formula[[3]]) %>>% str_detect("[+a-zA-z_]") %>>% length %>>% `==`(1)
  chkFormla <- chkFmLHS || chkFmRHS
  message("Checking the formula...")
  assert_that(chkFormla, msg = "Check failed")

  if (!FPCA_opts$parallel)
    RcppParallel::setThreadOptions(1)

  varName <- setdiff(all.vars(formula), as.character(formula[[3]]))
  timeVarName <- as.character(formula[[3]])
  sparsity <- checkSparsity(data, id.var, timeVarName)

  dataDT <- melt.data.table(data.table(data), id.vars = c(id.var, timeVarName),
                            measure.vars = varName, variable.factor = FALSE)
  if (length(unique(dataDT$variable)) > 1){

  }
  #  data[ , .N, by = timeVarName]




  return(1)
}
