#' example data: regular with missing
#'
#' This data is simulated. It is come from a random design, 200 subjects with 1 to 20
#' measurements on each subject. Time interval is [0,10]. The number of measurements
#' for each subject is uniformly distributed on {1,2,3, ..., 20}. The timepoints for
#' each subject are equally spaced values on [0,10]. The missing data is randomly drawn in each subject.
#' The number of true principal components is 2. One is N(0,9) and another one follows N(0,4).
#' The true mean function: t+sin(t). The two eigenfunctions are -sqrt(0.2)*cos(2*t*pi/10) and
#' sqrt0.2)*sin(2*t*pi/10), respectively. The measurement error follows standard normal distribution.
#'
#' @format A data.frame with 2956 observations on 3 variables.
#'
#' \itemize{
#'   \item \code{response}: the observed measurements for all subjects
#'   \item \code{timePoint}: the time points of samples for which corresponding measurements \code{response}.
#'   \item \code{sampleID}: the id of samples
#' }
#'
#' @source The example modified from example.m in PACE 2.17. You can find the generating script at
#'   \url{https://github.com/ChingChuan-Chen/rfda/blob/master/testData/exDataGeneration.m}.
#'
#' @docType data
#' @name irregularExData
#' @usage irregularExData
NULL

