#' Calculating the functional principal components
#'
#' Return mean, covariance, eigen functions and fpc scores.
#'
#' @param formula an object of class "formula", a symbolic description of the model to be fitted.
#' The RHS of `~` is the response and LHS is time variables.
#'
#' @param data
#'
#' @param y A list with n elements, \code{y[[i]]} is the vector of observed measurements for
#'   the ith subject, \code{i=1,...,n}.
#' @param t A list with n elements, \code{t[[i]]} is the vector of time points for the
#'   ith subject for which corresponding measurements \code{y[[i]]} are available,
#'   \code{i=1,...,n}.
#' @param FPCA_opts An R6 object. Please see \code{\link{get_FPCA_opts}} for more details.
#'   The default is NULL.
#' @param ... Extra arguments for \code{\link{get_FPCA_opts}}.
#' @return An R6 object containing mean, covariance, eigen functions,
#'   functional principal components scores, etc. Please see "Details" for more details.
#' @seealso \code{\link{get_FPCA_opts}} for input options.
#'   \code{\link{FPCA_options}} for option class.
#'   \code{\link{FPCA_results}} for output class.
#' @examples
#' library(magrittr)
#' \dontrun{
#' ## Following three examples may take long time.
#' library(dplyr)
#' # Taiwan freeway 5 data
#' data("TaiwanFreeway5vdData", package = 'rpace')
#' df <- TaiwanFreeway5vdData %>% filter(vdid == 'nfbVD-5S-9.063') %>%
#'   mutate(sampleID = sprintf("%04i_%02i_%02i", year, month, day), time = hour + (minute+1)/60)
#' splitData_TWFW5 <- df2FPCA_list_f("sampleID", c("volume", "time"), df)
#' res_FPCA_TWFW5 <- splitData_TWFW5 %$% FPCA(volume, time)
#'
#' # sparse case
#' data("exampleData_sparse", package = 'rpace')
#' splitData_s <- df2FPCA_list_f("sampleID", c("response", "timePoint"), exampleData_sparse)
#' res_FPCA_s <- splitData_s %$% FPCA(response, timePoint)
#'
#' # regular with missing values case
#' data("exampleData_regular_with_missing", package = 'rpace')
#' splitData_rm <- df2FPCA_list_f("sampleID", c("response", "timePoint"),
#'   exampleData_regular_with_missing)
#' res_FPCA_rm <- splitData_rm %$% FPCA(response, timePoint)
#' }
#'
#' # regular case
#' data("exampleData_regular", package = 'rpace')
#' splitData_r <- df2FPCA_list_f("sampleID", c("response", "timePoint"), exampleData_regular)
#' res_FPCA_r <- splitData_r %$% FPCA(response, timePoint)
#'
#' # convert FPCA_results to a list
#' res_FPCA_r_list <- res_FPCA_r$conv2list()
#'
#' # use the eval function to get the prediction at different time points
#' new_y_pred_1 <- res_FPCA_r$eval(c(1, 2, 6), seq(0.1, 0.5, 0.1))
#' new_y_pred_2 <- res_FPCA_r$eval(c(1, 2, 6), list(c(0.1, 0.2), c(0.3, 0.4, 0.5), c(0.6, 0.8, 1)))
#'
#' # use the predict function to get the prediction with different measurements
#' # at different time points
#' new_y_pred_3 <- res_FPCA_r$predict(list(c(2, 5, 6), c(3, 7, 9), c(4, 8, 11)), c(1, 3, 5))
#' new_y_pred_4 <- res_FPCA_r$predict(list(c(2, 5, 6), c(3, 7, 9), c(4, 8, 11)),
#'                                    list(c(1, 2, 3), c(3, 4, 5), c(6, 8, 10)))
#' @importFrom purrr map
#' @importFrom purrr map_lgl
#' @importFrom purrr map_int
#' @importFrom purrr map2
#' @importFrom purrr map2_lgl
#' @importFrom plyr llply
#' @importFrom stats quantile
#' @export
FPCA <- function(y, t, FPCA_opts = NULL, ...){
  if (is.null(FPCA_opts))
    FPCA_opts <- get_FPCA_opts()
  FPCA_opts$set(...)
  # check data
  assert_that(
    all(purrr::map_lgl(y, is.vector, mode = 'numeric')),
    all(purrr::map_lgl(y, ~!any(is.na(.)))),
    all(purrr::map_int(y, length) > 0),
    all(purrr::map_lgl(t, is.vector, mode = 'numeric')),
    all(purrr::map_lgl(t, ~!any(is.na(.)))),
    all(purrr::map_int(t, length) > 0),
    all(purrr::map_lgl(t, ~!any(duplicated(.)))),
    length(y) == length(t),
    all(purrr::map2_lgl(t, y, ~ length(.x) == length(.y))),
    all(class(FPCA_opts) %in% c("R6", "FPCA_options")))
  if (all(purrr::map_int(t, length) == 1))
    stop('FPCA is aborted because the data do not contain repeated measurements!')
  regular <- isregular(t) # the regularity of data
  FPCA_opts$check(length(t), regular) # check options
  if (is.null(FPCA_opts$get('numBins')))
  {
    resBinning <- binData(y, t, regular, FPCA_opts$get('verbose'))
  } else if (FPCA_opts$get('numBins') >= 10)
  {
    resBinning <- binData(y, t, regular, FPCA_opts$get('verbose'), FPCA_opts$get('numBins'))
  } else if (FPCA_opts$get('numBins') == 0)
  {
    resBinning <- list(newy = NULL, newt = NULL)
  } else if (FPCA_opts$get('numBins') < 10)
  {
    warning('The number of bins must be at least 10! No binning will be performed!')
    resBinning <- list(newy = NULL, newt = NULL)
  }
  if (!is.null(resBinning$newy))
  {
    y <- resBinning$newy
    t <- resBinning$newt
    if (regular == 0)
      regular <- 1
  }
  rm(resBinning)

  # pool all the subjects and their corresponding time points into 1 x N vectors
  # tt: 1 x N vector to hold observed time points from all subjects
  # yy: 1 x N vector to hold the observed measurements from all subjects
  yy <- unlist(y)
  tt <- unlist(t)
  # Initial out1 is based on the unique time points of the pooled data + the unique
  # time points of "newdata", the output time grid. When newdata = NULL, output
  # "out1" is equivalent to be the unique sorted pooled time points; otherwise, it
  # corresponds to the unique "newdata".
  out1 <- sort(unique(c(tt, FPCA_opts$get('newdata'))))
  out21 <- seq(min(out1), max(out1), length.out = FPCA_opts$get('numGrid'))
  res_FPCA <- FPCA_results$new()
  res_FPCA$set(regular = regular, out1 = out1, out21 = out21,
               y = y , t = t, opts_input = FPCA_opts)

  if (FPCA_opts$get('weight'))
  {
    if (regular == 0)
    {
      w <- purrr::map(t, ~rep(1/length(.), length(.)))
    } else
    {
      tmp_t <- sort(unique(tt))
      t_loc <- purrr::map(t, ~match(., tmp_t))
      count_t <- table(unlist(t_loc))
      w <- purrr::map(t_loc, ~1/count_t[.])
    }
  } else
  {
    w <- purrr::map(t, ~rep(1, length(.)))
  }
  ww <- w %>% do.call(c, .)

  if (FPCA_opts$get('verbose') == 'on')
    message('++++ I: Obtain smoothed mean curve. ++++')

  mu <- FPCA_opts$get('meanFunc')
  if (!is.null(mu) && (length(mu) == length(out1)))
  {
    muDense <- interp1(out1, as.matrix(mu), out21, 'spline') %>% as.vector
    bw_mu <- FPCA_opts$get('bwMean')
  } else
  {
    if (!is.null(FPCA_opts$get('meanFunc')))
      warning('The length of input meanFunc is not the same as the length of out1.')
    if (FPCA_opts$get('bwMean') == 0)
    {
      sampleID <- purrr::map(1:length(y), ~rep(., length(y[[.]]))) %>% do.call(c, .)
      bw_mu <- cvfda_lwls(yy, tt, ww, FPCA_opts$get('bwKernel'), 1, 1, 0,
                          regular, FPCA_opts$get('verbose'), sampleID, 0, FPCA_opts$get('parallel'))
    } else if (FPCA_opts$get('bwMean') %in% c(-1, -2))
    {
      bw_mu <- gcv_lwls(yy, tt, ww, FPCA_opts$get('bwKernel'), 1, 1, 0,
                        regular, FPCA_opts$get('verbose'), 0, FPCA_opts$get('parallel'))
      if (is.na(bw_mu))
        stop('FPCA is aborted because the observed data is too sparse to estimate the mean function!')
      if (FPCA_opts$get('bwMean') == -1)
      {
        bw_mu <- sqrt(find_mini_bw_f(tt, 2L) * bw_mu)
      }
      if (FPCA_opts$get('verbose') == 'on')
        message(sprintf('Adjusted GCV bandwidth choice for mean function: %0.6f', bw_mu))
    } else if (FPCA_opts$get('bwMean') > 0)
    {
      bw_mu <- FPCA_opts$get('bwMean')
    }
    mu <- lwls(bw_mu, FPCA_opts$get('bwKernel'), 1, 1, 0, tt, yy, ww, out1, 0.0,
               FPCA_opts$get('parallel')) %>% as.vector
    muDense <- lwls(bw_mu, FPCA_opts$get('bwKernel'), 1, 1, 0, tt, yy, ww, out21, 0.0,
                    FPCA_opts$get('parallel')) %>% as.vector
  }
  if (all(is.na(mu)))
    stop("The bandwidth of mean function is not appropriate! You may change the method of choosing bandwidth!")

  if (FPCA_opts$get('verbose') == 'on')
    message('++++ II: Choose bandwidth of smoothing covariance surface. ++++')

  rcov <- getRawCov(y, t, w, mu, out1, regular, 0)
  rcov %<>% plyr::llply(function(x){
    if (ncol(x) == 1)
      as.vector(x)
    else
      x
  })

  xcov <- FPCA_opts$get('covFunc')
  if (!is.null(xcov) && (length(xcov) == FPCA_opts$get('numGrid')**2))
  {
    bw_cov <- FPCA_opts$get('bwCov')
  } else
  {
    if (!is.null(FPCA_opts$get('covFunc')))
      warning('The size of input covFunc is not numGrid x numGrid.')
    if (all(FPCA_opts$get('bwCov') == 0))
    {
      bw_cov <- cv_mullwlsn(y, t, w, mu, FPCA_opts$get('bwcvSize'), regular, FPCA_opts$get('errorTerm'),
                            FPCA_opts$get('bwKernel'), rcov, FPCA_opts$get('verbose'),
                            FPCA_opts$get('parallel')) %>% as.vector
    } else if (all(FPCA_opts$get('bwCov') %in% c(-1, -2)))
    {
      bw_cov <- gcv_mullwlsn(t, FPCA_opts$get('bwNumGrid'), regular, FPCA_opts$get('errorTerm'),
                             FPCA_opts$get('bwKernel'), rcov, FPCA_opts$get('verbose'), 0,
                             FPCA_opts$get('parallel')) %>% as.vector
      if (any(is.na(bw_cov)))
        stop('FPCA is aborted because the observed data is too sparse to estimate the covariance function!')
      if (all(FPCA_opts$get('bwCov') == -2))
      {
        minbwcov <- find_mini_bw_f2(t, out1, regular)
        bw_cov <- sqrt(minbwcov*bw_cov)
      }
      if (FPCA_opts$get('verbose') == 'on')
        message(sprintf('Adjusted GCV bandwidth choice for COV function : %0.6f, %0.6f.', bw_cov[1], bw_cov[2]))
    } else if (all(FPCA_opts$get('bwCov') > 0))
    {
      bw_cov <- FPCA_opts$get('bwCov')
    }

    rcov1 <- rcov
    if (FPCA_opts$get('errorTerm'))
    {
      tneq <- rcov1$tpairn[ ,1] != rcov1$tpairn[ ,2]
      rcov1$tpairn <- rcov1$tpairn[tneq, ]
      rcov1$cxxn <- rcov1$cyy[tneq]
      rcov1$win <- rcov1$win[tneq]
      rcov1$indx <- rcov1$indx[tneq]
      rcov1$count <- rcov1$count[tneq]
    }
    if (regular == 2)
      rcov1$count <- rep(1, length(rcov1$count))
    # smooth raw covariance
    xcov <- mullwlsk(bw_cov, FPCA_opts$get('bwKernel'), rcov1$tpairn, rcov1$cxxn, rcov1$win,
                     out21, out21, rcov1$count, FPCA_opts$get('parallel'))
  }
  if (all(is.na(xcov)))
    stop("The bandwidth of covariance function is not appropriate! You may change the method of choosing bandwidth!")
  # transform the smoothed covariance matrix to guarantee it is a symmetric matrix.
  xcov <- (xcov + t(xcov)) / 2

  if (FPCA_opts$get('positiveVar'))
  {
    tmp_var <- diag(xcov)
    if (all(tmp_var < 0))
      stop("The 'positiveVar' option can't work out because the diagnal of xcov are all negative!")
    if (any(tmp_var < 0))
    {
      message(sprintf('Warning: Setting %d negative variance to minumum variance.', sum(tmp_var < 0)))
      tmp_var[tmp_var < 0] <- min(tmp_var[tmp_var > 0])
      xcov[cbind(which(tmp_var < 0), which(tmp_var < 0))] <- min(tmp_var[tmp_var > 0])
    }
    xvar <- sqrt(tmp_var %*% t(tmp_var))
    if (any(xvar < abs(xcov)))
    {
      loc <- xvar < abs(xcov)
      xcov[loc] <- sign(xcov[loc]) * xvar[loc]
    }
  }

  return(1)
}
