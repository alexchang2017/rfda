#' Perform cubic spline on data for interpolation
#'
#' Using cubic spline interpolation to find the values of a cubic function at
#' the values of correspondant values. The extrapolation is used, please be caution
#' in using the values which \code{xi} is larger than \code{max(x)} and smaller than \code{min(x)}.
#'
#' @param x A vector with n elements, \code{x[i]} is a support, \code{i = 1, ..., n}.
#'   If \code{x} is not sorted, it will be sorted. If \code{x} is not unique, the corresponding \code{y} values
#'   will be averaged.
#' @param y \code{y[i, j]} is jth values on corresponding value of \code{x[i]}, i = 1, ..., n.
#'   If \code{y} is vector, the length of \code{y} must be equal to the lenght of \code{x}.
#'   If \code{y} is matrix, the number of rows or the number of columns must be equal to the lenght of \code{x}.
#' @param xi A vector with m elements, \code{xi[k]} is the point which you want to interpolate,
#'   \code{k = 1, ..., m}.
#' @return A vector or matrix (depends on \code{y}) with the interpolated values corresponding to \code{xi}.
#' @section Reference:
#' Cleve Moler, Numerical Computing with MATLAB, chapter 3,
#'   \url{http://www.mathworks.com/moler/index_ncm.html}. \cr
#' Kai Habel, David Bateman, spline, Octave.
#' @examples
#' library(ggplot2)
#' plot_res <- function(x, y, xx, yy){
#'   ggplot(data.frame(x, y) , aes(x=x, y=y)) + geom_point() +
#'     geom_line(aes(x=x, y=y, colour = "spline"), data = data.frame(x=xx, y=yy)) +
#'     labs(title='Results of Interpolation', x='', y='')
#' }
#' x = 0:10
#' y = sin(x)
#' xx = seq(0, 10, 0.2)
#' yy = spline_f(x, as.matrix(y), xx)
#' plot_res(x, y, xx, yy)
#'
#' x <- c(0.8, 0.3, 0.1, 0.6, 0.9, 0.5, 0.2, 0.0, 0.7, 1.0, 0.4)
#' y <- matrix(c(x**2-0.6*x+1, 0.5*x**3-2*x**2+2*x+1), length(x))
#' xx <- seq(0, 1, len=81)
#' yy <- spline_f(x, y, xx)
#' plot_res(x, y[,1], xx, yy[,1])
#' plot_res(x, y[,2], xx, yy[,2])
#'
#' # example in spline function of MatLab
#' x <- seq(0, 2, 0.5) * pi
#' y <- matrix(c(0,1,0,-1,0,1,0,1,0,1,0,-1,0,1), 7)
#' yy <- spline_f(x, y, seq(0,2*pi,len=101))
#' ggplot(data.frame(x = y[2:5,1], y = y[2:5,2]) , aes(x=x, y=y)) + geom_point() +
#'   geom_path(aes(x=x, y=y), data = data.frame(x=yy[,1], y=yy[,2]), colour = "blue")
#' @export
spline_f <- function(x, y, xi){
  z <- spline_cpp(x, as.matrix(y), xi)
  if (ncol(z) == 1)
    z <- as.vector(z)
  return(z)
}

#' 1-D data interpolation.
#'
#' Returns interpolated values of a 1-D function at specific query points using linear interpolation.
#' The extrapolation is used, please be caution in using the values which \code{xi} is larger than
#' \code{max(x)} and smaller than \code{min(x)}.
#'
#' @param x A vector with n elements, \code{x[i]} is a support, \code{i = 1, ..., n}.
#'   If \code{x} is not sorted, it will be sorted. If \code{x} is not unique, the corresponding \code{y} values
#'   will be averaged.
#' @param y If \code{y} is vector, the length of \code{y} must be equal to the lenght of \code{x}.
#'   If \code{y} is matrix, the number of rows or the number of columns must be equal to the lenght of \code{x}.
#'   If the number of rows is equal to the lenght of \code{x, y[i, j]} is jth values on corresponding
#'   value of \code{x[i], i = 1, ..., n}.
#' @param xi A vector with m elements, \code{xi[k]} is the point which you want to interpolate,
#'   \code{k = 1, ..., m}.
#' @param method A string "linear" or "spline", the method of interpolation.
#' @return A vector or matrix (depends on \code{y}) with the interpolated values corresponding to
#'   \code{xi}.
#' @section Reference:
#' Cleve Moler, Numerical Computing with MATLAB, chapter 3,
#'   \url{http://www.mathworks.com/moler/index_ncm.html}. \cr
#' Nir Krakauer, Paul Kienzle, VZLU Prague, interp1, Octave.
#' @examples
#' library(ggplot2)
#' plot_res <- function(x, y, xi, yl, ys){
#'   ggplot(data.frame(x, y) , aes(x=x, y=y)) + geom_point() +
#'     geom_line(aes(x=x, y=y, colour = "linear"), data = data.frame(x=xi, y=yl)) +
#'     geom_line(aes(x=x, y=y, colour = "spline"), data = data.frame(x=xi, y=ys)) +
#'     scale_colour_manual(values = c("linear"="red", "spline"="blue")) +
#'     labs(title='Results of Interpolation', x='', y='', colour = 'Interpolation')
#' }
#' x <- c(0.8, 0.3, 0.1, 0.6, 0.9, 0.5, 0.2, 0.0, 0.7, 1.0, 0.4)
#' y <- matrix(c(x**2 - 0.6*x, 0.2*x**3 - 0.6*x**2 + 0.5*x), length(x))
#' xi <- seq(0, 1, len=81)
#' yl <- interp1(x, y, xi, 'linear')
#' ys <- interp1(x, y, xi, 'spline')
#' plot_res(x, y[,1], xi, yl[,1], ys[,1])
#' plot_res(x, y[,2], xi, yl[,2], ys[,2])
#'
#' x <- seq(0, 2*pi, pi/4)
#' y <- sin(x)
#' xi <- seq(0, 2*pi, pi/16)
#' yl <- interp1(x, as.matrix(y), xi, 'linear')
#' ys <- interp1(x, as.matrix(y), xi, 'spline')
#' plot_res(x, y, xi, yl, ys)
#' @export
interp1 <- function(x, y, xi, method){
  z <- interp1_cpp(x, as.matrix(y), xi, method)
  if (ncol(z) == 1)
    z <- as.vector(z)
  return(z)
}

#' One-dimensional kernel local polynominal smoother
#'
#' Perform one-dimensional kernel local polynominal smoother for data \code{(x,y)} with weight \code{w} on \code{xout}.
#'
#' @param bandwidth A single numerical value. The kernel smoothing parameter.
#' @param x A vector, the variable of of x-axis.
#' @param y A vector, the variable of of y-axis. \code{y[i]} is corresponding value of \code{x[i]}.
#' @param w A vector, the weight of data. \code{w[i]} is corresponding value of \code{x[i]}.
#' @param xout A vector, vector of output time points. It should be a sorted vecotr.
#' @param kernel A string. It could be 'gauss', 'gaussvar', 'epan' or 'quar'.
#' @param drv An integer, the order of derivative.
#' @param degree An integer, the degree of polynomial.
#' @return A estimated value on \code{xout} by one-dimensional kernel local polynominal smoother.
#' @examples
#' N <- 100
#' x <- runif(N, 0, 10)
#' y <- rnorm(N)
#' xout <- sort(runif(200, 0, 10))
#' est <- locPoly1d(1.2, x, y, rep(1, N), xout, 'gauss', 0, 1)
#' require(ggplot2)
#' ggplot(data.frame(x,y), aes(x,y)) + geom_point() +
#'   geom_line(aes(xout, est), data = data.frame(xout, est))
#' @export
locPoly1d <- function(bandwidth, x, y, w, xout, kernel, drv, degree){
  return(as.vector(locPoly1d_cpp(bandwidth, x, y, w, xout, kernel, drv, degree)))
}

locLinearRotate2d <- function(bandwidth, x, y, w, count, outMat, kernel){
  return(as.vector(locLinearRotate2d_cpp(bandwidth, x, y, w, count, outMat, kernel)))
}
