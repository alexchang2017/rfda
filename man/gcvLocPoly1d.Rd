% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{gcvLocPoly1d}
\alias{gcvLocPoly1d}
\title{Find the optimal bandwidth for one-dimensional kernel local polynominal smoother}
\usage{
gcvLocPoly1d(bwCand, x, y, w, kernel, drv, degree)
}
\arguments{
\item{bwCand}{A numerical vector for the candidates of bandwidth.}

\item{x}{A vector, the variable of of x-axis.}

\item{y}{A vector, the variable of of y-axis. \code{y[i]} is corresponding value of \code{x[i]}.}

\item{w}{A vector, the weight of data. \code{w[i]} is corresponding value of \code{x[i]}.}

\item{kernel}{A string. It could be 'gauss', 'gaussvar', 'epan' or 'quar'.}

\item{drv}{An integer, the order of derivative.}

\item{degree}{An integer, the degree of polynomial.}
}
\value{
A optimal bandwidth selected by minimizing gcv scores.
}
\description{
Find the optimal bandwidth used in \code{\link{locPoly1d}}.
}
\examples{
data("regularExData", package = 'rfda')
bwCand <- bwCandChooser(regularExData, "sampleID", "t", 2, "gauss", 1)
w <- rep(1, nrow(regularExData))
bwOpt <- gcvLocPoly1d(bwCand, regularExData$t, regularExData$y, w, "gauss", 0, 1)
bwOpt <- adjGcvBw(bwOpt, 2, "gauss", 0)
}

