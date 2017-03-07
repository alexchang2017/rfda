#' FPCA options settings
#'
#' Allow the users to set and examine a variety of FPCA options which affect
#' the way in which FPCA calculate its results.
#'
#' For convinient explianation, we denote the number of variables (functions) as \code{p},
#' the number of curves of a variable (function) as \code{n}, the number of time points as \code{nt} and
#' the number of observations is denoted by \code{N (= n * nt)}.
#' The options of FPCA:
#' \itemize{
#' \item \code{bwMean}: A data.frame or data.table with two columns naming \code{variable} and \code{value}.
#'   Default value is \code{NULL}. If \code{bwMean} is \code{NULL},
#'   \code{bwMean = -1} will be used in all variables.
#'   First column is the names of variables. Second column is the value of \code{bwMean}.
#'   The \code{bwMean} for variables can be different. The value of \code{bwMean} can be following:
#'   \itemize{
#'     \item Any positive number: user-specified bandwidth.
#'     \item \code{-1}: Geometric mean of the minimum bandwidth and the GCV bandwidth. [default]
#'     \item \code{-2}: Use generalized cross-validation to choose bandwidth automatically.
#'   }
#'   For Gaussian/Epanechnikov kernel, the optimal bandwidth from GCV would be adjusted.
#' \item \code{bwCov}: A data.frame or data.table with four columns naming \code{variable1}, \code{variable2},
#'   \code{value1} and \code{value2}. Default value is \code{NULL}. If \code{bwCov} is \code{NULL},
#'   \code{bwCov = c(-1, -1)} will be used in all variables.
#'   First two columns is the names of variables. Last two column is the value of \code{bwCov}.
#'   The \code{bwCov} for variables can be different. The value of \code{bwCov} can be following:
#'   \itemize{
#'     \item Any positive number: user-specified bandwidth.
#'     \item \code{-1, -1}: Geometric mean of the minimum bandwidth and the GCV bandwidth.
#'     \item \code{-2, -2}: Use generalized cross-validation to choose bandwidth automatically.
#'       [default for the cases \code{variable1} is equal to \code{variable2}]
#'     \item \code{-3, -3}: Only used in choosing the bandwidth of Cross-Covariance. Use the bandwidths in
#'       smoothing covariance surface. [default for the cases \code{variable1} is not equal to \code{variable2}]
#'   }
#'   If a row is \code{-1, -2}, we only take \code{-1}.
#'   For Gaussian/Epanechnikov kernel, the optimal bandwidth from GCV would be adjusted.
#' \item \code{bwNumGrid}: Default is \code{30}.
#'    An integer is the number of support points of covariance surface. (for GCV)
#'   A smaller \code{bwNumGrid} accelerate process at less accuracy.
#' \item \code{bwKernel}: A character string to define the kernel to be used in the
#'   smoothing procedures of mean and covariance surface.
#'   \itemize{
#'     \item \code{epan}: Epanechnikov kernel, \code{f(x) = 0.75*(1-x^2), -1 <= x <= 1}.
#'        [Default for dense designs with n_i >= 20]
#'     \item \code{gauss}: Gaussian kernel, \code{f(x) = exp(-x^2/2)/sqrt(2*pi), -4 <= x <= 4}.
#'       [Default for sparse designs, regular designs with missings, dense designs for \code{n_i < 20}]
#'     \item \code{gaussvar}: A variant of Gaussian kernel,
#'       \code{f(x) = exp(-x^2/2)/sqrt(2*pi)*(1.25-x^2/4), -4 <= x <= 4}.
#'     \item \code{quar}: quartic kernel, \code{f(x) = (1-x^2)^2*15/16, -1 <= x <= 1}.
#'   }
#'   Note 1: The Gaussian kernel is overall best for sparse designs but is slower than the other kernels
#'     and if computational speed is of concern then one may wish to use the Epanechnikov kernel
#'     also for the case of sparse designs.
#' \item \code{numBins}: The number of bins to implement data binning.
#'   \itemize{
#'     \item Any positive integer (\code{>=10}): Use \code{numBins} to implement data binning.
#'     \item \code{-1}: Activate an automatic process to decide whether to implement data binning.
#'     \item \code{0}: Not implement data binning. [Default]
#'   }
#' \item \code{errTerm}: To indicate whether there is the measurement error in the model.
#'   \itemize{
#'     \item \code{FALSE}: no measurement error given in the model.
#'     \item \code{TRUE}: measurement error is given in the model. [Default]
#'   }
#' \item \code{numGrid}: Default is \code{51}. An integer, number of support points in each direction of
#'   covariance surface when performing functional principal component analysis.
#'   \code{numGrid} must be greater than the number of functional principal components (\code{numFPC}).
#' \item \code{weight}: sample weight used in the local weighted least-square.
#'   \itemize{
#'     \item \code{FALSE}: weights of all samples are \code{1}. [Default]
#'     \item \code{TRUE}: weights of all samples are the inverse of number of observation for each subject,
#'       ie, \code{1/ni}. If \code{weight} is used, the criterion \code{'AIC'}, \code{'BIC'} in \code{numFPC}
#'       will be calculated with \code{weight}.
#'   }
#' \item \code{numFPC}: A string or any real positive integer.
#'   \code{'AIC'}, \code{'BIC'}, \code{'FVE'} or \code{'AIC_R'} can be selected for selecting
#'   by one of these criterion. A real positive integer is also avaiable.
#'   \itemize{
#'     \item \code{AIC}: To use AIC with pseudo-likelihood of measurements (marginal likelihood).
#'     \item \code{BIC}: To use BIC with pseudo-likelihood of measurements (marginal likelihood).
#'     \item \code{FVE}: The fraction of variance explained. Use scree plot approach to select number of
#'       functional principal components. The threshold is set by \code{FVE_threshold}.
#'       Please see \code{FVE_threshold} below. [Default]
#'     \item \code{AIC_R}: It is set to be \code{numGrid}. (This option is default for the functional regression.)
#'     \item \code{Any positive integer}: A user-specified number of functional principal components.
#'   }
#'   Note: \code{'BIC'} and \code{'FVE'} produce the most parsimonious models.
#' \item \code{FVE_threshold}: A positive number is between 0 and 1. It is the fraction of variance explained.
#'   This option is only used in the case \code{numFPC} = \code{'FVE'} to select the number of functional
#'   principal components that explain at least \code{FVE_threshold} of total variation.
#'   [Default is 0.85.]
#' \item \code{maxNumFPC}: Default is \code{min(20, n-1)}. An integer, the maximum number of functional
#'   principal components. This option is only used in the case  \code{numFPC} = \code{'AIC'} or \code{'BIC'}.
#'   Otherwise, \code{maxNumFPC * numVar} will be used in the multivariate functional data.
#'   Note: The number of fpc would use \code{min(x, k-1)},
#'     where \code{x} is the number of fpc selected by AIC, BIC, FVE or AIC_R or user-defined,
#'     \code{k} is the number of the real potitive eigen values.
#' \item \code{outPercent}: Default is \code{0}. A positive number is between \code{0} and \code{0.5}.
#'   This indicates that we will leave out \code{outPercent} data in the boundary.
#'   When performing local linear smoothing for mean functions and cross-covariance surface,
#'   all the data are used, but the output grid will be restricted within the reduced range.
#'   This option is used to alleviate the boundary effect.
#' \item \code{minMeasErr}: Default is \code{1e-6}. A numeric, the minimum of the measurement error variance.
#' \item \code{measErrOut}: Default is \code{0.5}. A positive number is between \code{0} and \code{0.5}.
#'   This indicates that we will leave out \code{measErrOut} data in the boundary for
#'   estimating the measurement error variance.
#' \item \code{methodFPCS}: A string, \code{'CE'}, \code{'IN'}, \code{'LS'} or \code{'WLS'}.
#'   The method to estimate fpc scores.
#'   \itemize{
#'     \item \code{CE}: The conditional expectation method is applied.
#'     \item \code{IN}: The integration method is applied.
#'     \item \code{LS}: The least-square method is applied.
#'     \item \code{WLS}: The weighted least-square method is applied. [Default]
#'   }
#'   Note: \code{'CE'}, \code{'LS'} or \code{'WLS'} can be applied for sparse data, irregular or regular data,
#'   but \code{'IN'} only be applied for regular data.
#' \item \code{rho}: A truncation threshold for the measurement error variance.
#'   This option is only used in the case \code{methodFPCS} = \code{'CE'}.
#'   \itemize{
#'     \item \code{cv}: To use a fixed rule to split the training data and validation data to
#'       find the optimal value of truncation threshold with cross-validation.
#'       (The results can be reproduced.) [Default]
#'     \item \code{cv-random}: To randomly split the training data and validation data to
#'       find the optimal value of truncation threshold with cross-validation.
#'       (The results can not be reproduced.)
#'     \item \code{no}: The truncation threshold will not be used.
#'     \item Any positive number: user-specified the measurement error variance.
#'       The length of rho must be equal to the number of variables.
#'   }
#' \item \code{shrink}: Whether to apply shrinkage method to estimate the FPC scores.
#'   (only for regular data.)
#'   \itemize{
#'     \item \code{TRUE}: The shrinkage method only used in the case \code{methodFPCS} = \code{'IN'} and
#'       \code{errTerm} = \code{TRUE}. Otherwise, this will be re-set to \code{FALSE}.
#'     \item \code{FALSE}: Not perform shrinkage method on eastimation of the FPC scores. [Default]
#'   }
#' \item \code{methodNorm}: The method to normalize the data.
#'   The default for univariate functional data is \code{'no'}.
#'   The default for multivariate functional data is \code{'quantile'}.
#'   \itemize{
#'     \item \code{no}: The normalization on data is not applied. Default value for univariate functional data.
#'     \item \code{quantile}: The normalization performed by substracting smoothed mean function and
#'       deviding by \code{0.75*IQR} of all curves. Default value for multivariate functional data.
#'     \item \code{smoothCov}: The normalization performed by substracting smoothed mean function and
#'       deviding by the variance of smoothed covariance surface. (Take the minimum positive value as threshold.).
#'       Please find the reference, Chiou, Chen and Yang. (2014), in \code{\link{rfda}}.
#'   }
#' \item \code{quantileProbs}: Default is \code{(0.25, 0.75)}.
#'   A positive numeric vector is between \code{0} and \code{1}.
#'   The probabilities of quantiles to approximate the variances.
#'   This option is only used in the case \code{methodNorm} = \code{'quantile'}.
#' \item \code{ncpus}: The number of threads used in computation.
#'   \itemize{
#'     \item \code{0}: To use all threads in computation. [Default]
#'     \item Any positive integer: To use user-specified number of threads in computation.
#'   }
#' \item \code{userMeanFunc}: Default is \code{NULL}. A data.table with three columns \code{timePnt}, \code{value}
#'   and \code{variable}. \code{value} is the mean function at specific time points.
#'   \code{timePnt} is the corresponding time points. \code{variable} is the name of observed variable.
#'   Please see the examples of \code{\link{FPCA}}.
#' \item \code{userCovFunc}: Default is \code{NULL}. A list of matrices with names \code{variable1}-\code{variable2}.
#'   The matrix is the cross-covariance surface of \code{variable1} and \code{variable2}.
#'   The colnames and rownames are the grid values of cross-covariance surface.
#'   If \code{userCovFunc} is given, \code{bwCov} must be given with real positive values of
#'   \code{value1} and \code{value2}.
#'   Please see the examples of \code{\link{FPCA}}.
#' }
#'
#' @param numVar The number of functional data to fit.
#' @param numCurves The number of curves of a variable.
#' @return An list of options to fit FPCA model used in \code{\link{FPCA}}.
#' @export
get_FPCA_opts <- function(numVar, numCurves){
  assert_that(length(numVar) == 1, is.numeric(numVar), is.finite(numVar), numVar > 0)
  assert_that(length(numCurves) == 1, is.numeric(numCurves), is.finite(numCurves), numCurves > 0)
  return(list(
    bwMean = NULL, bwCov = NULL, bwNumGrid = 30L, bwKernel = "gauss", numBins = 0L,
    errTerm = TRUE, numGrid = 51L, weight = FALSE, numFPC = "FVE", FVE_threshold = 0.85,
    maxNumFPC = min(20L, numCurves - 1L), outPercent = 0, minMeasErr = 1e-6, measErrOut = 0.5,
    methodFPCS = "WLS", rho = "cv", shrink = FALSE, methodNorm = ifelse(numVar == 1L, "no", "quantile"),
    quantileProbs = c(0.25, 0.75), ncpus = 0L, userMeanFunc = NULL, userCovFunc = NULL
  ))
}

#' @importFrom utils combn
chk_FPCA_opts <- function(optns, numVar){
  # check bwMean
  assert_that(is.null(optns$bwMean) || is.data.frame(optns$bwMean))

  # check bwCov
  assert_that(is.null(optns$bwCov) || is.data.frame(optns$bwCov))

  # check bwNumGrid
  assert_that(length(optns$bwNumGrid) == 1, is.numeric(optns$bwNumGrid), is.finite(optns$bwNumGrid),
              optns$bwNumGrid > 0)
  # check bwKernel
  assert_that(is.character(optns$bwKernel), optns$bwKernel %in% c("gauss", "gaussvar", "epan", "quar"))

  # check numBines
  if (!is.null(optns$numBins))
    assert_that(length(optns$numBins) == 1, is.numeric(optns$numBins), is.finite(optns$numBins),
                optns$numBins >= 0 || optns$numBins == -1)
  # check errTerm
  assert_that(length(optns$errTerm) == 1, is.logical(optns$errTerm), !is.na(optns$errTerm))
  # check numGrid
  assert_that(length(optns$numGrid) == 1, is.numeric(optns$numGrid), is.finite(optns$numGrid),
              optns$numGrid > 0)
  # check weight
  assert_that(length(optns$weight) == 1, is.logical(optns$weight), !is.na(optns$weight))

  # check numFPC
  assert_that(length(optns$numFPC) == 1)
  if (is.numeric(optns$numFPC)) {
    assert_that(is.finite(optns$numFPC), optns$numFPC > 0, optns$numFPC - floor(optns$numFPC) < 1e-6)
    if (optns$numFPC > optns$numGrid - 2) {
      warning(paste0("numFPC can only be less than or equal to numGrid-2!",
                     " Reset it to be numGrid-2 now!"))
      optns$numFPC <- optns$numGrid - 2
    }
  } else {
    assert_that(is.character(optns$numFPC), optns$numFPC %in% c("AIC", "BIC", "FVE", "AIC_R"))
    if (optns$numFPC == "FVE") {
      assert_that(length(optns$FVE_threshold) == 1, is.numeric(optns$FVE_threshold),
                  optns$FVE_threshold > 0, optns$FVE_threshold <= 1)
    } else if (optns$numFPC %in% c("AIC", "BIC")) {
      assert_that(length(optns$maxNumFPC) == 1, optns$maxNumFPC > 0,
                  is.numeric(optns$maxNumFPC), is.finite(optns$maxNumFPC))
      if (optns$maxNumFPC > optns$numGrid - 2) {
        warning(paste0("maxNumFPC can only be less than or equal to numGrid-2!",
                       " Reset it to be numGrid-2 now!"))
        optns$maxNumFPC<- optns$numGrid - 2
      }
    }
  }

  # check minMeasErr
  assert_that(length(optns$minMeasErr) == 1, is.finite(optns$minMeasErr), !is.na(optns$minMeasErr),
              is.numeric(optns$minMeasErr), optns$minMeasErr >= 0)
  # check measErrOut
  assert_that(length(optns$measErrOut) == 1, is.finite(optns$measErrOut), !is.na(optns$measErrOut),
              is.numeric(optns$measErrOut), optns$measErrOut >= 0, optns$measErrOut <= 0.5)

  # check methodFPCS
  assert_that(is.character(optns$methodFPCS), optns$methodFPCS %in% c("CE", "IN", "LS", "WLS"))

  # check shrink
  assert_that(length(optns$shrink) == 1, is.logical(optns$shrink), !is.na(optns$shrink))
  if (optns$shrink && (!optns$errTerm || optns$methodFPCS != "IN")) {
    warning("The shrinkage method only had effects when methodFPCS = 'IN' ",
            "and errTerm = TRUE! Reset to shrink = FALSE now!")
    optns$shrink <- FALSE
  }

  # check outPercent
  assert_that(length(optns$outPercent) == 1, is.finite(optns$outPercent), !is.na(optns$outPercent),
              is.numeric(optns$outPercent), optns$outPercent >= 0, optns$outPercent <= 0.5)

  # check rho
  if (is.character(optns$rho)) {
    assert_that(length(optns$rho) == 1, optns$rho %in% c("cv", "cv-random", "no"))
  } else if (is.numeric(optns$rho)) {
    assert_that(length(optns$rho) == numVar, all(is.finite(optns$rho)), all(!is.na(optns$rho)), all(optns$rho > 0))
  }

  # check outPercent
  assert_that(length(optns$outPercent) == 1, is.finite(optns$outPercent), !is.na(optns$outPercent),
              is.numeric(optns$outPercent), optns$outPercent >= 0, optns$outPercent <= 0.5)

  # check methodNorm
  assert_that(length(optns$methodNorm) == 1, optns$methodNorm %in% c("no", "quantile", "smoothCov"))

  # check quantileProbs
  assert_that(length(optns$quantileProbs) == 2, all(is.finite(optns$quantileProbs)),
              all(!is.na(optns$quantileProbs)), is.numeric(optns$quantileProbs),
              all(optns$quantileProbs >= 0), all(optns$quantileProbs <= 1),
              optns$quantileProbs[2] > optns$quantileProbs[1])

  # check ncpus
  assert_that(length(optns$ncpus) == 1, is.finite(optns$ncpus), !is.na(optns$ncpus),
              optns$ncpus - floor(optns$ncpus) < 1e-6, optns$ncpus >= 0)
  return(optns)
}
