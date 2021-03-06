% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{splitMat}
\alias{splitMat}
\title{trapz: trapezoidal rule to approximate the integral values}
\usage{
splitMat(m, margin, f)
}
\arguments{
\item{m}{A numeric matrix to be divided into list of matrices.}

\item{margin}{The margin of the matrix to split.}

\item{f}{A integer vector to split the matrix.}
}
\value{
A list of matrices.
}
\description{
Returns approximation of integral.
}
\examples{
x <- matrix(rnorm(30), 6, 8)
splitMat(x, 1, rep(1:3, each = 2))
splitMat(x, 2, rep(1:4, each = 2))
}

