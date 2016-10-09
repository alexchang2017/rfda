#' Perform cubic spline on data for interpolation
#'
#' Using cubic spline interpolation to find the values of a cubic function at
#' the values of correspondant values. The extrapolation is used, please be caution
#' in using the values which xi is larger than max(x) and smaller than min(x).
#'
#' @param x A vector with n elements, \code{x[i]} is a support, i = 1, ..., n.
#'   If x is not sorted, it will be sorted. If x is not unique, the corresponding y values
#'   will be averaged.
#' @param y \code{y[i, j]} is jth values on corresponding value of \code{x[i]}, i = 1, ..., n.
#'   If y is vector, the length of y must be equal to the lenght of x.
#'   If y is matrix, the number of rows or the number of columns must be equal to the lenght of x.
#' @param xi A vector with m elements, \code{xi[k]} is the point which you want to interpolate,
#'   k = 1, ..., m.
#' @return A vector or matrix (depends on y) with the interpolated values corresponding to \code{xi}.
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
#' The extrapolation is used, please be caution in using the values which xi is larger than
#' max(x) and smaller than min(x).
#'
#' @param x A vector with n elements, x[i] is a support, i = 1, ..., n.
#'   If x is not sorted, it will be sorted. If x is not unique, the corresponding y values
#'   will be averaged.
#' @param y If y is vector, the length of y must be equal to the lenght of x.
#'   If y is matrix, the number of rows or the number of columns must be equal to the lenght of x.
#'   If the number of rows is equal to the lenght of x, y[i, j] is jth values on corresponding
#'   value of x[i], i = 1, ..., n.
#' @param xi A vector with m elements, xi[k] is the point which you want to interpolate,
#'   k = 1, ..., m.
#' @param method A string "linear" or "spline", the method of interpolation.
#' @return A vector or matrix (depends on y) with the interpolated values corresponding to
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
