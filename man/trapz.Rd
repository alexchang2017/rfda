% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/commFunc.R
\name{trapz}
\alias{trapz}
\title{trapz: trapezoidal rule to approximate the integral values}
\usage{
trapz(x, y = NULL)
}
\arguments{
\item{x}{A vector with \code{n} elements, \code{x[i]} is a support, \code{i = 1, ..., n}.
If \code{y} is \code{NULL}, support is taken as \code{seq(1, length(x), by = 1)}.}

\item{y}{\code{y[i, j]} is jth values on corresponding value of \code{x[i], i = 1, ..., n}.
If \code{y} is vector, the length of \code{y} must be equal to the lenght of \code{x}.
If \code{y} is matrix, the number of rows must be equal to the lenght of \code{x}.}
}
\value{
A value, the approximation of integral.
}
\description{
Returns approximation of integral.
}
\section{Reference}{

Kai Habel, trapz, Octave.
}
\examples{
# case 1
x <- c(1, 4, 9, 16, 25)
trapz(x) # 42

# case 2
x <- matrix(c(1,4,9,16,25,1,8,27,64,125), 5)
trapz(x) # 42 162
}

