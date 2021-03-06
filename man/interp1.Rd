% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interp.R
\name{interp1}
\alias{interp1}
\title{1-D data interpolation.}
\usage{
interp1(x, y, xi, method)
}
\arguments{
\item{x}{A vector with n elements, \code{x[i]} is a support, \code{i = 1, ..., n}.
If \code{x} is not sorted, it will be sorted. If \code{x} is not unique, the corresponding \code{y} values
will be averaged.}

\item{y}{If \code{y} is vector, the length of \code{y} must be equal to the lenght of \code{x}.
If \code{y} is matrix, the number of rows or the number of columns must be equal to the lenght of \code{x}.
If the number of rows is equal to the lenght of \code{x, y[i, j]} is jth values on corresponding
value of \code{x[i], i = 1, ..., n}.}

\item{xi}{A vector with m elements, \code{xi[k]} is the point which you want to interpolate,
\code{k = 1, ..., m}.}

\item{method}{A string "linear" or "spline", the method of interpolation.}
}
\value{
A vector or matrix (depends on \code{y}) with the interpolated values corresponding to
  \code{xi}.
}
\description{
Returns interpolated values of a 1-D function at specific query points using linear interpolation.
The extrapolation is used, please be caution in using the values which \code{xi} is larger than
\code{max(x)} and smaller than \code{min(x)}.
}
\section{Reference}{

Cleve Moler, Numerical Computing with MATLAB, chapter 3,
  \url{http://www.mathworks.com/moler/index_ncm.html}. \cr
Nir Krakauer, Paul Kienzle, VZLU Prague, interp1, Octave.
}
\examples{
library(ggplot2)
plot_res <- function(x, y, xi, yl, ys){
  ggplot(data.frame(x, y) , aes(x=x, y=y)) + geom_point() +
    geom_line(aes(x=x, y=y, colour = "linear"), data = data.frame(x=xi, y=yl)) +
    geom_line(aes(x=x, y=y, colour = "spline"), data = data.frame(x=xi, y=ys)) +
    scale_colour_manual(values = c("linear"="red", "spline"="blue")) +
    labs(title='Results of Interpolation', x='', y='', colour = 'Interpolation')
}
x <- c(0.8, 0.3, 0.1, 0.6, 0.9, 0.5, 0.2, 0.0, 0.7, 1.0, 0.4)
y <- matrix(c(x**2 - 0.6*x, 0.2*x**3 - 0.6*x**2 + 0.5*x), length(x))
xi <- seq(0, 1, len=81)
yl <- interp1(x, y, xi, 'linear')
ys <- interp1(x, y, xi, 'spline')
plot_res(x, y[,1], xi, yl[,1], ys[,1])
plot_res(x, y[,2], xi, yl[,2], ys[,2])

x <- seq(0, 2*pi, pi/4)
y <- sin(x)
xi <- seq(0, 2*pi, pi/16)
yl <- interp1(x, as.matrix(y), xi, 'linear')
ys <- interp1(x, as.matrix(y), xi, 'spline')
plot_res(x, y, xi, yl, ys)
}

